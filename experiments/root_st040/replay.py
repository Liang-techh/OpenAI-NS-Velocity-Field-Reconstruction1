"""Reload small frozen recipes; no training or hidden source coefficients required."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import bootstrap
from temporal_fit import TemporalModel,TemporalBiasModel
from structure_audit import audit
from validate import validate

HERE=Path(__file__).resolve().parent

def load(ident='ST042'):
 data=json.loads((HERE/'recipes.json').read_text())
 if ident not in data:raise ValueError('Unknown frozen recipe')
 d=data[ident]
 if d['pde_validated'] is not False or d['source_correspondence_verified'] is not False:raise ValueError('Unsupported scientific promotion')
 if hashlib.sha256(bootstrap.PARENT.read_bytes()).hexdigest()!=d['parent_sha256']:raise ValueError('Parent field changed')
 m=(TemporalBiasModel if d['bias'] else TemporalModel)(bootstrap.PARENT)
 c=np.asarray(d['modifiers'],float)
 if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Bad modifier vector')
 raw=m.candidate(c)
 u,p=m.f.fields(raw,np.asarray(d['reference_points']),np.asarray(d['reference_times']))
 if not (np.allclose(u,d['reference_velocity'],rtol=0,atol=d['reference_absolute_tolerance']) and np.allclose(p,d['reference_pressure'],rtol=0,atol=d['reference_absolute_tolerance'])):raise ValueError('Recipe replay disagrees with frozen evaluations')
 return m.f,raw

def main():
 p=argparse.ArgumentParser();p.add_argument('--id',default='ST042',choices=['ST041','ST042']);p.add_argument('--out',type=Path,required=True);p.add_argument('--validate',action='store_true');p.add_argument('--seed',type=int,default=9174101);p.add_argument('--audit',action='store_true');a=p.parse_args()
 if a.out.exists() and any(a.out.iterdir()):p.error('Output directory must be empty')
 a.out.mkdir(parents=True,exist_ok=True);f,raw=load(a.id);candidate=a.out/'candidate.json'
 f.save(raw,candidate,dict(recipe=a.id,scope='Numerically reconstructed recipe; not original metadata bytes',pde_validated=False))
 if a.audit:(a.out/'structure.json').write_text(json.dumps(audit(f,raw,97),indent=2)+'\n')
 if a.validate:
  r=validate(candidate,a.out/'validation.json',a.seed)
  if not r['all_numeric_gates_pass']:raise SystemExit(1)
 print('Replayed '+a.id+'; pde_validated=false')
if __name__=='__main__':main()
