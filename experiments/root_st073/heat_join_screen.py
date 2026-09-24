"""Paper heat-exterior direct-join screen against actual corrected ST073 traces."""
import sys,json,numpy as np
from radial_continuation import ROOT,FullRadialField
from interface_stress import gradient
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from heat_exterior import physical

def run():
 f=FullRadialField.load(ROOT/'radial_continuation/candidate.json');x=3/64
 p=f.from_similarity([x],[0.],.5);u=f.evaluate(p,.5)['velocity'][0,1]
 coefficient=float(u/physical(p,.5)['velocity'][0,1]);rows=[]
 for k in (0,3,6):
  tau=.5*2**(-k);eta=np.array([-.4,0,.4]);pts=f.from_similarity(np.full(3,x),eta,tau);inner=f.evaluate(pts,tau);outer=physical(pts,tau,c=coefficient)
  for i,e in enumerate(eta):
   jac=gradient(f,x,e,tau);core_shear=f.nu*(jac[1,0]+jac[0,1]);r=pts[i,0];h=r*1e-4
   plus=pts[i:i+1].copy();minus=plus.copy();plus[0,0]+=h;minus[0,0]-=h
   derivative=(physical(plus,tau,c=coefficient)['velocity'][0,1]-physical(minus,tau,c=coefficient)['velocity'][0,1])/(2*h)
   outer_shear=f.nu*(derivative-outer['velocity'][i,1]/r)
   rows.append(dict(k=k,eta=float(e),velocity_jump=(outer['velocity'][i]-inner['velocity'][i]).tolist(),inner_rtheta_stress=float(core_shear),outer_rtheta_stress=float(outer_shear),rtheta_stress_jump=float(outer_shear-core_shear)))
 report=dict(heat_amplitude=coefficient,fit='One amplitude matched to midplane swirl at k0 only, then frozen',rows=rows,source='OpenAI Navier-Stokes PDF Appendix A heat exterior; existing shipped evaluator reused',scope='Direct interface compatibility screen, not moment repair, wave realization or finite total energy. Pure swirl exterior carries no axial/radial flow.',global_field_ready=False)
 out=ROOT/'heat_join';out.mkdir(exist_ok=True);(out/'screen.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows[6:],indent=2))
if __name__=='__main__':run()
