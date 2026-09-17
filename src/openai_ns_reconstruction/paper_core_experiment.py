"""Reproducible off-grid radial-degree study; no full-NS acceptance claim."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .paper_core_series import PaperCoreSeries
from .paper_core_reference import PaperCoreReference
from .paper_profile_residual import residual


def run(output, degrees=(2,4,6,8), nodes=(129,257), step=5e-5, sigma=.1, radii=(.001,.01,.1,.3), eta_count=None):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    results=[]
    for n in nodes:
        for degree in degrees:
            start=time.monotonic()
            model=PaperCoreSeries(PaperCoreReference(sigma=sigma),maxdegree=degree,eta_nodes=n)
            rows=[]
            for X in radii:
                for eta in ([-.5,-.1,0,.1,.5] if eta_count is None else np.linspace(-.95,.95,eta_count)):
                    row=dict(X=X,eta=eta)
                    try:
                        r=residual(model.profile,X,eta,model.reference.h,step)
                        row.update(residual=r.tolist(),F=float(model.F(X,eta)),U=float(model.U(X,eta)))
                        if not np.all(np.isfinite(r)):
                            row['error']='nonfinite residual'
                    except (ValueError,ArithmeticError) as e:
                        row['error']=str(e)
                    rows.append(row)
            by_radius=[]
            for X in radii:
                valid=[v['residual'] for v in rows if v['X']==X and 'error' not in v]
                by_radius.append(dict(X=X,valid=len(valid),total=5 if eta_count is None else eta_count,
                    maximum=np.max(np.abs(valid),axis=0).tolist() if valid else None))
            entry=dict(degree=degree,eta_nodes=n,step=step,seconds=time.monotonic()-start,
                       by_radius=by_radius,rows=rows,metadata=model.metadata())
            results.append(entry)
            (output/'convergence.json').write_text(json.dumps(dict(
                status='development samples; leading equations only; no convergence certificate',
                residual_order=['angular_4.13','axial_4.13','pressure_4.7'],runs=results),indent=2)+'\n')
            print(json.dumps({k:entry[k] for k in ['degree','eta_nodes','seconds','by_radius']}),flush=True)
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',default='artifacts/function_first/core_series')
    p.add_argument('--degrees',type=int,nargs='+',default=[2,4,6,8])
    p.add_argument('--nodes',type=int,nargs='+',default=[129,257])
    p.add_argument('--step',type=float,default=5e-5)
    p.add_argument('--sigma',type=float,default=.1)
    p.add_argument('--radii',type=float,nargs='+',default=[.001,.01,.1,.3])
    p.add_argument('--eta-count',type=int)
    args=p.parse_args()
    run(args.output,args.degrees,args.nodes,args.step,args.sigma,args.radii,args.eta_count)
