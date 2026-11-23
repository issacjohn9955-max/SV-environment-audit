"""Build a read-only audit of genotype inventory, previous registry and climate evidence."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd
from PIL import Image

VARIABLES = ('bio3', 'bio15', 'bio18', 'wind', 'srad')
POINT = re.compile(r'^(-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)$')


def table(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, encoding='utf-8-sig')


def check_ids(table, name):
    if table[name].eq('').any() or not table[name].is_unique:
        raise ValueError(f'Blank or duplicate ID in {name}')


def write(table, path):
    table.to_csv(path, sep='\t', index=False, float_format='%.9g')


def sha(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def origin_table(path):
    section = path.read_text().split('## 附录 A', 1)[1].split('## 附录 B', 1)[0]
    result = []
    for line in section.splitlines():
        if not re.match(r'^\| .* \| \d+ \|', line):
            continue
        breed, size, rows, legacy, tier, region, coords = [s.strip() for s in line.strip('|').split(' | ')]
        m = POINT.fullmatch(coords)
        lon, lat = map(float, m.groups()) if m else (None, None)
        result.append(dict(breed=breed, reviewed_registry_n=int(size), legacy_rows=rows,
                           legacy_n=legacy, origin_tier=tier, origin_basis=region,
                           listed_lon=lon, listed_lat=lat))
    frame = pd.DataFrame(result)
    check_ids(frame, 'breed')
    if set(frame.origin_tier) - set('ABCDEU') or len(frame) == 0:
        raise ValueError('Unexpected origin evidence codes.')
    return frame


def climate(points, raster_dir):
    subset = points[points.origin_tier.isin(('A', 'B'))].dropna(subset=['listed_lon', 'listed_lat']).copy()
    row = np.floor((90 - subset.listed_lat.to_numpy())*12).astype(int)
    col = np.floor((180 + subset.listed_lon.to_numpy())*12).astype(int)
    if ((row < 0) | (row >= 2160) | (col < 0) | (col >= 4320)).any():
        raise ValueError('Representative coordinate out of raster bounds.')

    def sample(layer):
        with Image.open(raster_dir/f'wc2.1_5m_{layer}.tif') as image:
            scale, tie = image.tag_v2[33550], image.tag_v2[33922]
            if image.size != (4320, 2160) or not np.isclose(scale[:2], 1/12).all() or tie[3:5] != (-180, 90):
                raise ValueError(f'Unexpected WorldClim grid for {layer}.')
            values = np.asarray(image)[row, col].astype(float)
        values[(~np.isfinite(values)) | (values <= -1e30) | (values == -9999)] = np.nan
        return values

    for name in ('3', '15', '18'):
        subset[f'bio{name}'] = sample(f'bio_{name}')
    for name, conversion in [('wind', 1), ('srad', 86.4)]:
        subset[name] = np.mean([sample(f'{name}_{month:02d}') for month in range(1, 13)], axis=0)/conversion
    subset['grid_row_0based'] = row
    subset['grid_column_0based'] = col
    subset['climate_complete'] = subset[list(VARIABLES)].notna().all(axis=1)
    return subset


def header_check(record):
    expected, path = record
    if not path:
        return 'vcf_path_missing', ''
    try:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8') as file:
            for line in file:
                if line.startswith('#CHROM\t'):
                    ids = line.rstrip('\n\r').split('\t')[9:]
                    status = ('exact_single_sample' if ids == [expected] else
                              'exact_multisample' if expected in ids else 'expected_id_absent')
                    return status, '|'.join(ids)
                if not line.startswith('##'):
                    break
    except (OSError, EOFError, UnicodeError):
        return 'vcf_header_read_error', ''
    return 'vcf_header_missing', ''


def main(args):
    env_dir = args.output/'environment'
    audit_dir = args.output/'review'
    env_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    latest = table(args.inventory)
    old = pd.read_excel(args.registry, sheet_name='个体明细', dtype=str, keep_default_na=False)
    old = old.rename(columns={'个体ID': 'id', '品种/类群': 'previous_registry_breed'})
    check_ids(latest, 'id')
    check_ids(old, 'id')
    headers = latest[['vcf_sample_id', 'vcf']].itertuples(index=False, name=None)
    with ThreadPoolExecutor(max_workers=8) as pool:
        status, observed = zip(*pool.map(header_check, headers))
    current = latest[['id', 'vcf_sample_id', 'breed', 'population_code', 'sources', 'genotype_type']].copy()
    current['vcf_header_status'] = status
    current['observed_vcf_header_ids'] = observed
    joined = current.merge(old[['id', 'previous_registry_breed']], on='id', how='outer',
                           indicator='inventory_vs_previous', validate='one_to_one')
    joined['breed_label_exact_match'] = (
        joined.breed.eq(joined.previous_registry_breed) & joined.inventory_vs_previous.eq('both')
    )

    points = origin_table(args.origin_review)
    extracted = climate(points, args.rasters)
    points['coordinate_review_status'] = np.where(points.origin_tier.isin(('A', 'B')),
                                                  'representative_point_review_needed_per_batch',
                                                  'hold_or_exclude')
    aliases = table(args.breed_crosswalk)
    check_ids(aliases, 'raw_breed')
    aliases = aliases.set_index('raw_breed')
    if not aliases.origin_review_breed.isin(points.breed).all():
        raise ValueError('Breed crosswalk refers to unknown origin review labels.')
    joined['origin_review_breed'] = joined.breed.map(aliases.origin_review_breed).fillna(joined.breed)
    joined['breed_crosswalk_evidence'] = joined.breed.map(aliases.evidence)
    lookup = points[['breed', 'origin_tier', 'coordinate_review_status']].rename(
        columns={'breed': 'origin_review_breed'})
    joined = joined.merge(lookup, on='origin_review_breed', how='left', validate='many_to_one')
    kept = extracted[['breed', 'listed_lon', 'listed_lat', 'climate_complete', *VARIABLES]].rename(
        columns={'breed': 'origin_review_breed'})
    joined = joined.merge(kept, on='origin_review_breed', how='left', validate='many_to_one')
    joined['coordinate_review_status'] = joined.coordinate_review_status.fillna('no_documented_breed_link')
    joined['individual_environment_review'] = 'pending_batch_or_origin_verification'
    exact_environment = np.isfinite(joined[list(VARIABLES)].to_numpy(float)).all(axis=1)
    joined.loc[~exact_environment, 'individual_environment_review'] = 'no_complete_candidate_climate'
    joined.loc[joined.origin_tier.eq('E'), 'individual_environment_review'] = 'commercial_exclude'
    identity_hold = joined.origin_tier.eq('U') | joined.id.isin([f'CRR{i}' for i in range(2303436, 2303440)])
    library_hold = joined.id.isin([f'ERR{i}' for i in range(5633123, 5633139)])
    joined.loc[identity_hold, 'individual_environment_review'] = 'identity_conflict_hold'
    joined.loc[library_hold, 'individual_environment_review'] = 'library_type_hold'
    joined.loc[identity_hold | library_hold, [*VARIABLES, 'listed_lon', 'listed_lat']] = np.nan
    joined.loc[identity_hold | library_hold, 'climate_complete'] = False
    exact_environment = np.isfinite(joined[list(VARIABLES)].to_numpy(float)).all(axis=1)

    historic = table(args.historical_samples)
    check_ids(historic, 'sample')
    ids = subprocess.check_output([args.bcftools, 'query', '-l', str(args.historical_bcf)], text=True).splitlines()
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate ID in historical BCF sample header.')
    env_columns = ['sample', 'canonical_sample', 'breed', 'group', 'environment_tier', 'lon', 'lat', *VARIABLES]
    previous = historic[env_columns].rename(columns={'sample': 'bcf_sample', 'breed': 'historical_breed'})
    previous['in_historical_bcf_header'] = previous.bcf_sample.isin(ids)
    previous['exact_canonical_id_in_latest'] = previous.canonical_sample.isin(current.id)
    previous['historical_use'] = 'previous_405_cohort_only'
    write(previous, audit_dir/'historical_main405_alignment.tsv')
    write(points, env_dir/'origin_evidence_1726_snapshot.tsv')
    write(extracted, env_dir/'candidate_climate_by_breed.tsv')
    write(joined, audit_dir/'individual_alignment.tsv')
    issues = joined.loc[
        ~joined.inventory_vs_previous.eq('both') | ~joined.breed_label_exact_match |
        ~joined.vcf_header_status.eq('exact_single_sample') |
        ~joined.individual_environment_review.eq('pending_batch_or_origin_verification')
    ]
    write(issues, audit_dir/'review_queue.tsv')
    payload = {
        'latest_inventory_unique': len(latest), 'prior_registry_unique': len(old),
        'inventory_vs_previous': dict(Counter(joined.inventory_vs_previous)),
        'raw_breed_label_differences_among_shared_ids': int(
            (joined.inventory_vs_previous.eq('both') & ~joined.breed_label_exact_match).sum()),
        'vcf_header_status': dict(Counter(status)),
        'origin_evidence_by_tier': dict(Counter(points.origin_tier)),
        'representative_points_with_complete_climate': int(extracted.climate_complete.sum()),
        'candidate_individuals_with_complete_climate': int(exact_environment.sum()),
        'individual_environment_review': dict(Counter(joined.individual_environment_review)),
        'historical_bcf_header_n': len(ids), 'historical_environment_n': len(previous),
        'historical_env_ids_absent_from_bcf': previous.loc[
            ~previous.in_historical_bcf_header, 'bcf_sample'].tolist(),
        'historical_environment_groups': historic.group.nunique(),
        'source_sha256': {name: sha(getattr(args, name)) for name in
                          ('inventory', 'registry', 'origin_review', 'breed_crosswalk', 'historical_samples')},
        'scope': 'VCF header identity and preliminary environment review; cohort approval pending.'
    }
    report = json.dumps(payload, ensure_ascii=False, indent=2) + '\n'
    (audit_dir/'audit_summary.json').write_text(report)
    print(report, end='')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for parameter in ('inventory', 'registry', 'origin-review', 'breed-crosswalk', 'rasters', 'historical-samples', 'historical-bcf', 'output'):
        parser.add_argument(f'--{parameter}', type=Path, required=True)
    parser.add_argument('--bcftools', default='bcftools')
    main(parser.parse_args())
