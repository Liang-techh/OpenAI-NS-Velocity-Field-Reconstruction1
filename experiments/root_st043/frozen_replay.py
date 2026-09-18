"""Replay frozen joint-correction recipes, with no fitting or external field data."""
from __future__ import annotations
import argparse,json,tempfile
from pathlib import Path
import numpy as np
import bootstrap
from replay import load as load_parent
from localized_model import LocalizedModel
from coupled_fit import CoupledModel
from structural_report import report
from validate import validate

HERE=Path(__file__).resolve().parent

def construct(d):
 if d.get('schema')!='coupled_ST042_correction_v1':raise ValueError('Unknown correction schema')
 if d.get('pde_validated') is not False or d.get('source_correspondence_verified') is not False:raise ValueError('Unsupported scientific claim')
 parent_records=json.loads((HERE.parent/'root_st040/recipes.json').read_text())
 if d.get('parent_recipe')!='ST042' or d.get('parent_original_sha256')!=parent_records['ST042']['frozen_raw_file_sha256']:raise ValueError('Parent identity changed')
 f,raw=load_parent('ST042');opts=d['model'];constructor=LocalizedModel if opts['localized'] else CoupledModel
 with tempfile.TemporaryDirectory() as folder:
  path=Path(folder)/'parent.json';f.save(raw,path,{'scope':'replayed parent metadata'})
  m=constructor(path,radial=opts['radial'],axial=opts['axial'],time_degree=opts['time_degree'])
 c=np.asarray(d['modifiers'],float)
 if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Invalid modifier array')
 limits=np.asarray(m.bounds)
 if np.any(c<limits[:,0]-1e-12) or np.any(c>limits[:,1]+1e-12):raise ValueError('Modifier bounds violated')
 raw=m.candidate(c);u,p=m.f.fields(raw,np.asarray(d['reference_points']),np.asarray(d['reference_times']))
 if not np.allclose(u,d['reference_velocity'],atol=d['reference_atol'],rtol=0) or not np.allclose(p,d['reference_pressure'],atol=d['reference_atol'],rtol=0):raise ValueError('Frozen numerical reference mismatch')
 return m.f,raw

def load(ident='ST045-H'):
 entries=json.loads((HERE/'recipes.json').read_text())
 if ident not in entries:raise ValueError('Unknown frozen candidate')
 return construct(entries[ident])

def main():
 p=argparse.ArgumentParser();p.add_argument('--id',default='ST045-H',choices=['ST045-G','ST045-H','ST044-F']);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9174301);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
 if a.out.exists() and any(a.out.iterdir()):p.error('Output directory must be empty')
 a.out.mkdir(parents=True,exist_ok=True);f,raw=load(a.id);candidate=a.out/'candidate.json';f.save(raw,candidate,{'recipe':a.id,'scope':'frozen mathematical field replay, not original metadata bytes','pde_validated':False})
 if a.structure:
  fp,rp=load_parent('ST042');parent=a.out/'parent_ST042.json';fp.save(rp,parent)
  (a.out/'structure.json').write_text(json.dumps(report(candidate,parent),indent=2)+'\n')
 if a.validate:
  r=validate(candidate,a.out/'validation.json',a.seed)
  if not r['all_numeric_gates_pass']:raise SystemExit(1)
 print('Frozen recipe reconstructed; pde_validated=false')
if __name__=='__main__':main()
