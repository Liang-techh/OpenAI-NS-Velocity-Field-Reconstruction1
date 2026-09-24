from pathlib import Path
import time,json
import numpy as np
from full_radial import *
ROOT=Path(__file__).resolve().parent
path=ROOT/'evidence/full_order_calibration.json'
if path.exists():raise FileExistsError(path)
rows=[]
for N in [2,4,6,8,10,12]:
 f=FullRadialField(Parameters(order=N));t=time.monotonic();mx=0.;div=0.;axis=[]
 for k in [0,3,6]:
  xx,ee=np.meshgrid([.003,.009,.015625],[-.5,-.25,0,.25,.5]);d=f.evaluate_similarity(xx,ee,.5*2**(-k))
  vals=np.linalg.norm(d['residual'],axis=-1);mx=max(mx,float(vals.max()));div=max(div,float(abs(d['divergence']).max()))
  axis.append(dict(k=k,max=float(vals.max())))
 row=dict(order=N,max=mx,divergence=div,scales=axis,seconds=time.monotonic()-t);rows.append(row)
 path.write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
