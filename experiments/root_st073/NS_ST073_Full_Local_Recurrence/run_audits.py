from pathlib import Path
import json,time
from audit import holdout,physical,fd_audit,save,ROOT
for seed in [9237391,9237392]:holdout(seed)
for k in [0,3,6]:
 for n in [8,12,18]:
  dst=ROOT/'evidence'/f'physical_k{k}_n{n}.json'
  if dst.exists():continue
  t=time.monotonic();row=physical(k,n);row['seconds']=time.monotonic()-t;save(dst,row);print('PHYS',k,n,row['ST073-F'],flush=True)
dst=ROOT/'evidence/independent_fd.json'
if not dst.exists():save(dst,fd_audit())
print('ALL AUDITS COMPLETE',flush=True)
