"""Representable h fixtures for the exact shared phase change of variables.

The fixture h values test algebra only; they are not paper parameter
admission or replacements for the implicit production parameter.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
HERE=Path(__file__).parent

def run():
 c=MPIntervalContext();c.dps=160;checks=0
 with mp.workdps(200):
  for power in (3,10,1000):
   h=c.mpf(2)**(-power)
   for xi_string in ('0','.25','.5','.75','1','1.25','1.5','1.75','2'):
    xi=c.mpf(xi_string);y=h*xi
    for actual,expected in (
      (alpha_box(c,xi,c.mpf(1)),alpha_box(c,y,h)),
      (h+(1-h)*alpha_box(c,xi+1,c.mpf(1)),h+(1-h)*alpha_box(c,y+h,h))):
     if endpoints(actual)!=endpoints(expected):raise ArithmeticError('phase cutoff differs')
     checks+=1
   # A known constant ODE integrated in phase must retain physical dy=h dxi.
   n=32;step=h/n;acc=c.mpf(0)
   for _ in range(n):acc+=step
   if not endpoints(acc)[0]<=endpoints(h)[0]<=endpoints(acc)[1]:raise ArithmeticError('dy lost')
   checks+=1
  names=('lei_ren_part1_paper_shared_comparison.json','lei_ren_part1_paper_shared_exit_bridge.json',
    'lei_ren_part1_paper_shared_exit_continuation.json','lei_ren_part1_paper_shared_exit_switch.json',
    'lei_ren_part1_paper_shared_R110_cone.json')
  receipts=[json.loads((HERE/n).read_bytes()) for n in names]
  sha=receipts[0]['parameter_family_sha256']
  for r in receipts:
   if r['parameter_family_sha256']!=sha or not r['h_b_equals_exit_epsilon_by_definition']:
    raise ValueError('mixed coupled parameter family')
   if r['full_Section9_parameter_admission']:raise ValueError('unproved admission promoted')
   for n,d in r['input_hashes'].items():
    if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('stale dependency: '+n)
  result=dict(phase_cutoff_and_dy_checks=checks,fixture_h_powers_of_two=[-3,-10,-1000],
    production_family_sha256=sha,all_five_receipts_share_same_h_parameter=True,
    fixture_is_algebra_only=True,full_Section9_parameter_admission=False,
    input_hashes={**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
      Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'lei_ren_part1_paper_interval_comparison_enclosure.py':hashlib.sha256((HERE/'lei_ren_part1_paper_interval_comparison_enclosure.py').read_bytes()).hexdigest()})
  Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print('PASS:',checks,'phase/dy checks; five receipts share one parameter family')
 return result
if __name__=='__main__':run()
