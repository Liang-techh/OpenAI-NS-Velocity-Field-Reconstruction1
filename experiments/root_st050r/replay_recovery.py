"""Recover immutable mathematical fields, without rerunning optimization."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,sys,tempfile
from pathlib import Path
import numpy as np
import pressure_morph as pm
from spacetime import Family
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
PARENT_HASHES={'13127759debcf40b3f3d6c3627747ace67f1471b825b2b0e8c53ae79a9911344','6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09'}

def parent_field():
    p=ROOT/'artifacts/research/ST048-S/candidate.json'
    if p.exists():
        if hashlib.sha256(p.read_bytes()).hexdigest() not in PARENT_HASHES:raise ValueError('Parent checksum mismatch')
        return Family.load(p)
    sys.path.insert(0,str(ROOT/'experiments/root_st048'))
    from replay_st048 import reconstruct
    return reconstruct('ST048-S')

def reconstruct(ident,records=None):
    records=json.loads((HERE/'recipes.json').read_text()) if records is None else records
    rec=records[ident]
    if rec.get('parent_id')!='ST048-S' or rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Wrong parent or unsupported claim')
    data=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(data).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Modifier checksum mismatch')
    c=np.load(io.BytesIO(data),allow_pickle=False);f,raw=parent_field()
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'parent.json';f.save(raw,p);m=pm.LocalizedModel(p,rec['radial'],rec['axial'],rec['td'])
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Malformed modifiers')
    lo,hi=np.array(m.bounds).T
    if np.any((c<lo-1e-10)|(c>hi+1e-10)):raise ValueError('Modifier bound violation')
    raw=m.candidate(c);r=rec['reference'];u,p=m.f.fields(raw,np.array(r['points']),np.array(r['times']))
    if not (np.allclose(u,r['velocity'],rtol=1e-8,atol=1e-8) and np.allclose(p,r['pressure'],rtol=1e-8,atol=1e-8)):raise ValueError('Numerical reference mismatch')
    return m.f,raw

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9175091);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Refuse nonempty output directory')
    a.out.mkdir(parents=True,exist_ok=True);f,raw=reconstruct(a.id);path=a.out/'candidate.json';f.save(raw,path,dict(recipe=a.id,scope='Exact mathematical recipe; raw JSON metadata and hash differ from training artifact'))
    if a.structure:
        from diagnostics import report
        fb,rb=parent_field();parent=a.out/'parent.json';fb.save(rb,parent);report(path,parent,a.out/'structure.json')
    if a.validate:
        from validate import validate
        r=validate(path,a.out/'validation.json',a.seed);failed=[k for k,v in r['gates'].items() if not v]
        print('FAIL: '+', '.join(failed) if failed else 'Reported numerical gates pass; not a continuum certificate')
        return int(bool(failed))
    return 0
if __name__=='__main__':raise SystemExit(main())
