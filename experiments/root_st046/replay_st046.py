"""Reconstruct a frozen modifier recipe and run independent validation. No fit."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,tempfile
from pathlib import Path
import numpy as np
import acceleration_fit
from localized_model import LocalizedModel
from frozen_replay import load as previous_load
from validate import validate
from mechanism_audit import audit

HERE=Path(__file__).resolve().parent

def reconstruct(ident,records=None):
    records=json.loads((HERE/'recipes.json').read_text()) if records is None else records
    rec=records[ident]
    if rec.get('pde_validated') is not False or rec.get('parent_id')!='ST045-H':raise ValueError('Unsupported claim/parent')
    data=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(data).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Modifier checksum')
    c=np.load(io.BytesIO(data),allow_pickle=False)
    f,raw=previous_load('ST045-H')
    with tempfile.TemporaryDirectory() as t:
        path=Path(t)/'parent.json';f.save(raw,path)
        m=LocalizedModel(path,rec['radial'],rec['axial'],rec['time_degree'])
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Malformed modifiers')
    lo,hi=np.array(m.bounds).T
    if np.any((c<lo-1e-12)|(c>hi+1e-12)):raise ValueError('Modifier bound violation')
    out=m.candidate(c);m.f.coefficients(out)
    ref=rec['reference'];u,p=m.f.fields(out,np.asarray(ref['points']),np.asarray(ref['times']))
    if not (np.allclose(u,ref['velocity'],rtol=1e-8,atol=1e-8) and np.allclose(p,ref['pressure'],rtol=1e-8,atol=1e-8)):
        raise ValueError('Reconstructed field differs from frozen numerical references')
    return m.f,out

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9174601);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Output must be empty')
    a.out.mkdir(parents=True,exist_ok=True);f,x=reconstruct(a.id);path=a.out/'candidate.json'
    f.save(x,path,dict(recipe_id=a.id,scope='Mathematical-field replay; metadata and file hash differ from original raw JSON'))
    if a.structure:(a.out/'structure.json').write_text(json.dumps(audit(path),indent=2)+'\n')
    if a.validate:
        r=validate(path,a.out/'validation.json',a.seed)
        failed=[k for k,v in r['gates'].items() if not v]
        print('FAIL: '+', '.join(failed) if failed else 'All reported numeric gates pass; not a continuum proof')
        return 1 if failed else 0
    return 0
if __name__=='__main__':raise SystemExit(main())
