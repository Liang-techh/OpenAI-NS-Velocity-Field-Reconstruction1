"""Frozen modifier replay from the immutable mathematical ST047-E parent."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,tempfile
from pathlib import Path
import numpy as np
import boundary_shear
from localized_model import LocalizedModel
from spacetime import Family
from validate import validate
from structure_and_particles import diagnostics
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
PARENT_SHA='dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0'

def parent_field():
    path=ROOT/'artifacts/research/ST047-E/candidate.json'
    if path.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest()!=PARENT_SHA:raise ValueError('Parent raw checksum mismatch')
        return Family.load(path)
    from replay_st047 import reconstruct
    return reconstruct('ST047-E')

def reconstruct(ident,records=None):
    records=json.loads((HERE/'recipes.json').read_text()) if records is None else records
    rec=records[ident]
    if rec.get('parent_id')!='ST047-E' or rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Wrong parent or unsupported scientific acceptance')
    data=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(data).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Modifier checksum mismatch')
    c=np.load(io.BytesIO(data),allow_pickle=False);f,raw=parent_field()
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'parent.json';f.save(raw,p)
        m=LocalizedModel(p,rec['radial'],rec['axial'],rec['time_degree'])
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Invalid modifier vector')
    lo,hi=np.asarray(m.bounds).T
    if np.any((c<lo-1e-10)|(c>hi+1e-10)):raise ValueError('Modifier bounds violated')
    result=m.candidate(c);ref=rec['reference'];u,p=m.f.fields(result,np.array(ref['points']),np.array(ref['times']))
    if not(np.allclose(u,ref['velocity'],atol=1e-8,rtol=1e-8) and np.allclose(p,ref['pressure'],atol=1e-8,rtol=1e-8)):raise ValueError('Frozen numerical reference mismatch')
    return m.f,result

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9174801);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Refuse to overwrite an existing output directory')
    a.out.mkdir(parents=True,exist_ok=True);f,raw=reconstruct(a.id);path=a.out/'candidate.json'
    f.save(raw,path,dict(recipe=a.id,scope='Frozen mathematical field replay; raw JSON metadata/hash can differ'))
    if a.structure:diagnostics(path,a.out)
    if a.validate:
        r=validate(path,a.out/'validation.json',a.seed);failed=[k for k,v in r['gates'].items() if not v]
        print('FAIL: '+', '.join(failed) if failed else 'All numerical gates passed; not a continuum proof')
        return int(bool(failed))
    return 0
if __name__=='__main__':raise SystemExit(main())
