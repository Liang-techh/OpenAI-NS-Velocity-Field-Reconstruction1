"""Unconstrained numerical profile-drift minimum in the two-patch basis.

This diagnoses basis capacity only: momentum, moment, cone and shape gates
are absent. The output must not be adopted as a physical candidate.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def run():
    path=ROOT/'localized_drift_constraint.npz'
    with np.load(path,allow_pickle=False) as cache:
        G,g,b=cache['gram'],cache['linear'],float(cache['constant'])
    scales=np.sqrt(np.maximum(np.diag(G),0.))
    active=scales>0
    H=G[np.ix_(active,active)]/scales[active,None]/scales[None,active]
    H=(H+H.T)/2
    rhs=g[active]/scales[active]
    eigenvalues,Q=np.linalg.eigh(H)
    rows=[]
    for rcond in (1e-8,1e-10,1e-12):
        keep=eigenvalues>rcond*eigenvalues.max()
        y=-Q[:,keep]@((Q[:,keep].T@rhs)/eigenvalues[keep])
        c=np.zeros(len(g));c[active]=y/scales[active]
        square=float(c@G@c+2*g@c+b)
        rows.append(dict(relative_eigenvalue_cutoff=rcond,rank=int(keep.sum()),
                         drift_squared=square,drift=float(np.sqrt(max(square,0.))),
                         physical_coefficients=c.tolist()))
    report=dict(status='completed',accepted=False,pde_validated=False,
                scale_recursion_established=False,
                scope='Unconstrained numerical least-squares capacity of velocity columns only; NOT a feasible joint candidate or rigorous continuum lower bound.',
                cache_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                parent_drift=float(np.sqrt(b)),active_velocity_columns=int(active.sum()),
                zero_response_columns=int((~active).sum()),
                minimum_normalized_eigenvalue=float(eigenvalues.min()),rows=rows)
    (ROOT/'localized_drift_floor.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({**{k:v for k,v in report.items() if k!='rows'},
                      'rows':[{k:v for k,v in row.items() if k!='physical_coefficients'} for row in rows]}),flush=True)


if __name__=='__main__':run()
