"""Test coordinate-level real SV frequencies against individual climate variables."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from partial_rank import partial_rank


def run(frequency, design, settings, pc_count, out):
    groups = pd.read_csv(design, sep='\t', index_col='id', dtype={'id': str})
    y = pd.read_csv(frequency, sep='\t', index_col='id', dtype={'id': str})
    if not groups.index.is_unique or not y.index.equals(groups.index) or not y.columns.is_unique:
        raise ValueError('Frequency and group design must have the same ordered unique group IDs and unique loci.')
    variables = settings['environments']
    covariates = [f'PC{i}' for i in range(1, pc_count + 1)]
    if settings['adjust_geography']:
        covariates.extend(('lon', 'lat'))
    rho, p, n, detail = partial_rank(
        y.to_numpy(float), groups[variables].to_numpy(float), groups[covariates].to_numpy(float),
        min_groups=8, return_diagnostics=True
    )
    files = []
    for j, environment in enumerate(variables):
        valid = np.isfinite(p[:, j])
        q = np.full(len(valid), np.nan)
        if valid.any():
            q[valid] = stats.false_discovery_control(p[valid, j], method='bh')
        files.append(pd.DataFrame({
            'locus': y.columns, 'environment': environment, 'pc_count': pc_count,
            'rho': rho[:, j], 'p': p[:, j], 'q': q, 'n_groups': n, 'df': detail['df'],
            'test_status': detail['reason'][:, j]
        }))
    result = pd.concat(files, ignore_index=True)
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out/'associations.tsv.gz', sep='\t', index=False)
    result[result.q < settings['alpha']].to_csv(out/'discoveries.tsv', sep='\t', index=False)
    result.groupby(['environment', 'test_status']).size().rename('loci').to_csv(out/'testability.tsv', sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pc_count', type=int, choices=(3, 5, 10))
    parser.add_argument('prepared', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--frequency', type=Path, help='Optional matching group-frequency TSV including partially observed loci.')
    args = parser.parse_args()
    run(args.frequency or args.prepared/'rda_frequency.tsv.gz', args.prepared/'group_design.tsv',
        json.loads((args.prepared/'settings.json').read_text()), args.pc_count, args.out)
