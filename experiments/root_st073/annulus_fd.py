"""Independent Cartesian derivative check at new annulus points."""
from pathlib import Path
import sys,json,numpy as np
from radial_continuation import ROOT,FullRadialField
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from local_field import independent_fd

def run():
 f=FullRadialField.load(ROOT/'radial_continuation/candidate.json');rows=[]
 for k in (.4,2.7,5.5):
  tau=.5*2**(-k);points=f.from_similarity(np.array([.026,.044]),np.array([-.22,.37]),tau,angle=[.3,.8]);reference=f.evaluate(points,tau)['residual']
  for hs,ht in ((.002,.0005),(.001,.0005),(.001,.00025)):
   residual,div=independent_fd(f,points,tau,hs*np.sqrt(f.nu*tau),ht*tau)
   rows.append(dict(k=k,space_factor=hs,time_factor=ht,maximum_momentum_residual=float(np.max(np.linalg.norm(residual,axis=1))),reference_difference=float(np.max(np.linalg.norm(residual-reference,axis=1))),divergence_max=float(np.max(np.abs(div)))))
   print(k,hs,ht,rows[-1]['maximum_momentum_residual'],flush=True)
 report=dict(rows=rows,scope='Six off-calibration physical points inside new annulus, independent Cartesian/time FD, separate spatial and time step variations. Not a uniform bound or global/exterior validation.',global_field_ready=False)
 (ROOT/'radial_continuation/independent_fd.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
