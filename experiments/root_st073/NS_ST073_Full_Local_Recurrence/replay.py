"""Recompute LOCAL frozen momentum only. Exit0 is not original global NS admission."""
from pathlib import Path
import argparse,json,hashlib
from audit import physical,ROOT,check,save
from audit_phase2 import physical as vphysical,binding

def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--candidate',choices=['ST073-F','ST073-V'],default='ST073-V');p.add_argument('--k',type=float,default=6);p.add_argument('--order',type=int,default=18);p.add_argument('--out',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
 if not 0<=a.k<=6:p.error('Only the registered finite k[0,6] window is supported')
 if not 4<=a.order<=48:p.error('Quadrature order must be between4and48')
 fr=binding() if a.candidate=='ST073-V' else check()
 sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['full_radial.py','audit.py','audit_phase2.py','upstream/source_coordinates.py','upstream/local_field.py','upstream/core_series.py','upstream/general_core.py','replay.py']}
 identity=dict(candidate=a.candidate,k=a.k,quadrature_order=a.order,model_sha256=fr['model_sha256'],code_sha256=sources)
 if a.out.exists():
  if not a.resume:raise FileExistsError('Output exists; use --resume for identity verification')
  envelope=json.loads(a.out.read_text());payload=envelope['payload']
  if payload['identity']!=identity or digest(payload)!=envelope['payload_sha256']:raise ValueError('Existing report identity or digest mismatch')
  print('Verified existing LOCAL report; no computation or fitting. Global field remains unconstructed.')
  return int(not payload['local_sampled_momentum_pass'])
 r=vphysical(a.k,a.order) if a.candidate=='ST073-V' else physical(a.k,a.order)
 vals=r['full'] if a.candidate=='ST073-V' else r['ST073-F'];maximum=vals['max'] if a.candidate=='ST073-V' else vals['sampled_max'];norm=vals['L2'] if a.candidate=='ST073-V' else vals['full_local_L2'];passed=maximum<.001 and norm<.001
 payload=dict(identity=identity,report=r,local_sampled_momentum_pass=bool(passed),pde_validated=False,global_field_ready=False,original_global_benchmark='NOT EVALUATED / NOT CONSTRUCTED',scope='Finite small local source-coordinate domain, unforced; no global energy or matching',optimizer_executed=False)
 save(a.out,dict(payload=payload,payload_sha256=digest(payload)))
 print('LOCAL sampled momentum:', 'PASS' if passed else 'FAIL','| original global candidate: NOT CONSTRUCTED')
 return int(not passed)
if __name__=='__main__':raise SystemExit(main())
