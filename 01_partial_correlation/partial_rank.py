"""Partial rank kernel from the 2026-11-18 observed-SV workflow."""
import numpy as np
from scipy import stats

def partial_rank(y,x,cov,min_groups=8,relative_tolerance=1e-10,return_diagnostics=False):
 """Partial rank correlation with explicit non-estimable states and relative residual checks.

 Ranks of response/environment are centred before projection. Continuous covariates
 are centred/scaled for a stable QR (same column space, no covariate ranking).
 Missing patterns use their own observed rows. No untestable p enters BH downstream.
 """
 y=np.asarray(y,float);x=np.asarray(x,float);cov=np.asarray(cov,float)
 assert y.ndim==x.ndim==cov.ndim==2 and len(y)==len(x)==len(cov)
 shape=(y.shape[1],x.shape[1]);rhos=np.full(shape,np.nan);p=rhos.copy();ns=np.zeros(y.shape[1],int)
 reasons=np.full(shape,'insufficient_observed_groups',dtype='U40');dfs=np.full(y.shape[1],np.nan);ranks=np.zeros(y.shape[1],int)
 masks=np.isfinite(y);_,patterns=np.unique(np.packbits(masks.T,axis=1),axis=0,return_inverse=True)
 for pattern in np.unique(patterns):
  cols=np.flatnonzero(patterns==pattern);valid=masks[:,cols[0]]&np.isfinite(x).all(1)&np.isfinite(cov).all(1);n=int(valid.sum());ns[cols]=n
  if n<min_groups:continue
  c=cov[valid].copy()
  if c.shape[1]:
   c-=c.mean(0);scale=np.sqrt((c*c).mean(0));c=np.divide(c,scale,out=np.zeros_like(c),where=scale>0)
  design=np.column_stack([np.ones(n),c]);rank=np.linalg.matrix_rank(design);ranks[cols]=rank;df=n-rank-1;dfs[cols]=df
  if rank!=design.shape[1]:reasons[cols]='rank_deficient_covariates';continue
  if n<rank+3:reasons[cols]='insufficient_residual_df';continue
  q,_=np.linalg.qr(design,mode='reduced')
  xr=stats.rankdata(x[valid],axis=0);xr-=xr.mean(0);xnorm=np.linalg.norm(xr,axis=0);xx=xr-q@(q.T@xr);xres=np.linalg.norm(xx,axis=0)
  xconstant=np.ptp(x[valid],axis=0)==0;xgood=(~xconstant)&(xres>relative_tolerance*xnorm)&(xnorm>0)
  for start in range(0,len(cols),10000):
   cc=cols[start:start+10000];raw=y[np.ix_(valid,cc)];yconstant=np.ptp(raw,axis=0)==0
   yr=stats.rankdata(raw,axis=0);yr-=yr.mean(0);ynorm=np.linalg.norm(yr,axis=0);yy=yr-q@(q.T@yr);yres=np.linalg.norm(yy,axis=0)
   ygood=(~yconstant)&(ynorm>0)&(yres>relative_tolerance*ynorm)
   reasons[cc]='degenerate_response';reasons[cc[yconstant]]='constant_response'
   reasons[np.ix_(cc[ygood],np.flatnonzero(xconstant))]='constant_environment'
   reasons[np.ix_(cc[ygood],np.flatnonzero((~xconstant)&(~xgood)))]='degenerate_environment'
   yc=np.flatnonzero(ygood);xc=np.flatnonzero(xgood)
   if len(yc) and len(xc):
    rr=(yy[:,yc]/yres[yc]).T@(xx[:,xc]/xres[xc]);rr=np.clip(rr,-1,1)
    with np.errstate(divide='ignore',invalid='ignore'):tt=rr*np.sqrt(df/(1-rr*rr))
    rhos[np.ix_(cc[yc],xc)]=rr;p[np.ix_(cc[yc],xc)]=2*stats.t.sf(abs(tt),df);reasons[np.ix_(cc[yc],xc)]='tested'
 diagnostics={'reason':reasons,'df':dfs,'design_rank':ranks,'relative_tolerance':relative_tolerance}
 return (rhos,p,ns,diagnostics) if return_diagnostics else (rhos,p,ns)
