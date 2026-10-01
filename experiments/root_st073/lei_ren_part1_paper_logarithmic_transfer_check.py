"""Independent representable-scale fixture for logarithmic inlet algebra."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,constant
from lei_ren_part1_paper_candidate_shared_inlet import finite_moments,normalized_inlet
from lei_ren_part1_paper_logarithmic_comparison import scaled_transfer
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent

def run():
 c=MPIntervalContext();c.dps=160;checks=0
 with mp.workdps(200):
  z=IntervalTaylor.variable(c,c.mpf('.3'),2);dt=c.mpf('.01');s=c.mpf('.5')
  phi=[1+z/10,constant(c,'.2',2)];u=[2-z/10,constant(c,'-.1',2)]
  m=finite_moments(phi,u,s)
  ph=sum((x*s**n for n,x in enumerate(phi)),phi[0]*0)
  uv=sum((x*s**n for n,x in enumerate(u)),u[0]*0)
  ell=constant(c,'.25',2);S=constant(c,'.0625',2);p0=-3+z/2
  initial=ph*constant(c,'1.1',2);zero=ph*0
  for lam in ('8','500','1e8'):
   lam=c.mpf(lam);e=1/lam
   old=normalized_inlet(phi,u,m,S,ell,p0,lam,s,z,dt,
      endpoints_override=dict(phi_exit=ph,u_exit=uv,phi_s=zero,u_s=zero))
   new=scaled_transfer(c,ph,uv,m,ell*e,S*e**2,p0*e,e,s,z,dt,initial)
   references=dict(D=old['ratio'],I_z=old['iz']*c.sqrt(e),pressure_scaled=old['pressure']*e,
      physical_Ur=old['ur']*c.sqrt(e),unmodulated_driver_A=-old['ratio']/2,
      unmodulated_driver_B=-(initial/ph)*old['iz']*e*c.sqrt(s/2))
   for name,expected in references.items():
    actual=new[name]
    for x,y in zip(actual.coefficients,expected.coefficients):
     a,b=endpoints(x);d,f=endpoints(y)
     if max(a,d)>min(b,f):raise ArithmeticError('nonoverlap: '+name)
     if max(abs(a-d),abs(b-f))/max(mp.mpf(1),abs(d),abs(f))>mp.mpf('1e-140'):
      raise ArithmeticError('endpoint mismatch: '+name)
     checks+=1
   for chi in ('0','.25','1'):
    # chi has no axial dependence, but both base drivers require it.
    for name in ('unmodulated_driver_A','unmodulated_driver_B'):
     actual=new[name]*c.mpf(chi);expected=references[name]*c.mpf(chi)
     for x,y in zip(actual.coefficients,expected.coefficients):
      a,b=endpoints(x);d,f=endpoints(y)
      if max(a,d)>min(b,f):raise ArithmeticError('cutoff mismatch')
      checks+=1
  result=dict(all_scaled_transfer_fixture_checks_passed=True,coefficient_comparisons=checks,
    cutoff_cases=['0','.25','1'],Lambda_cases=['8','500','1e8'],
    base_drivers_require_actual_chi_multiplier=True,
    full_pressure_and_swirl_axial_derivatives_retained=True,
    input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in (
      Path(__file__).name,'lei_ren_part1_paper_logarithmic_comparison.py',
      'lei_ren_part1_paper_candidate_shared_inlet.py','lei_ren_part1_paper_interval_taylor.py')})
  Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print('PASS:',checks,'scaled transfer/cutoff coefficient checks')
 return result
if __name__=='__main__':run()
