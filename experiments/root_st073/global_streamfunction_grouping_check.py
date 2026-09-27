"""Compare grouped streamfunction evaluation with the frozen scalar implementation."""
import hashlib
import json
import subprocess
import time
import types
from pathlib import Path
import numpy as np
import global_axial_extension as current

ROOT=Path(__file__).resolve().parent
BASE='08f7b6ee'


def run():
    source=subprocess.check_output(['git','show',BASE+':experiments/root_st073/global_axial_extension.py'],cwd=ROOT,text=True)
    previous=types.ModuleType('scalar_global_extension')
    previous.__file__=str(ROOT/'global_axial_extension.py')
    exec(compile(source,previous.__file__,'exec'),previous.__dict__)
    _,localized,_,_,snapshot,_=current.build_candidate()
    rows=[]
    for dk in (0.,.001):
        tau=snapshot['inputs']['mean']['tau']*2**(-dk)
        X=localized.join_X*np.array([.2,.9,1.2,localized.ratio**2*1.1])
        xx,ee,aa=np.meshgrid(X,[-.38,0.,.38],np.arange(8)*np.pi/4,indexing='ij')
        points=localized.inner.from_similarity(xx.ravel(),ee.ravel(),tau,angle=aa.ravel())
        t=time.perf_counter();old=previous._joined_streamfunction(localized.joined,points,tau);scalar_seconds=time.perf_counter()-t
        t=time.perf_counter();new=current._joined_streamfunction(localized.joined,points,tau);grouped_seconds=time.perf_counter()-t
        error=float(np.max(np.abs(new-old)))
        if not np.allclose(new,old,rtol=1e-13,atol=1e-20):
            raise AssertionError('Grouped primitive differs from scalar construction')
        rows.append(dict(delta_k=dk,point_count=len(points),max_abs_difference=error,
                         scalar_seconds=scalar_seconds,grouped_seconds=grouped_seconds,
                         speedup=scalar_seconds/grouped_seconds))
    report=dict(status='completed',reference_commit=BASE,
        current_source_sha256=hashlib.sha256((ROOT/'global_axial_extension.py').read_bytes()).hexdigest(),
        scope='Numerical equivalence and timing for the joined streamfunction helper; no PDE acceptance.',rows=rows)
    (ROOT/'global_streamfunction_grouping_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    run()
