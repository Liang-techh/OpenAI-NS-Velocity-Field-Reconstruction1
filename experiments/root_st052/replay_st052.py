"""Checksum-checked reconstruction of the frozen field; never fits coefficients."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,tempfile
from pathlib import Path
import numpy as np
from minimax_exchange import ROOT,EdgeModel,PARENT_SHA
from spacetime import Family
HERE=Path(__file__).resolve().parent

def parent_field():
    path=ROOT/'artifacts/research/ST051-B/candidate.json'
    if path.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest()!=PARENT_SHA:raise ValueError('Wrong parent raw bytes')
        return Family.load(path)
    from replay_st051 import reconstruct
    return reconstruct('ST051-B')

def reconstruct(rec=None):
    if rec is None:rec=json.loads((HERE/'recipe.json').read_text())
    if rec.get('parent_sha256')!=PARENT_SHA or rec.get('parent_id')!='ST051-B':raise ValueError('Parent identity mismatch')
    if rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Unsupported scientific claim')
    payload=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(payload).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Modifier checksum mismatch')
    c=np.load(io.BytesIO(payload),allow_pickle=False);f,raw=parent_field()
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'parent.json';f.save(raw,p);m=EdgeModel(p,edge_modes=True)
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Invalid modifier array')
    lo,hi=np.array(m.bounds).T
    if np.any((c<lo-1e-10)|(c>hi+1e-10)):raise ValueError('Modifier bounds violated')
    raw=m.candidate(c);ref=rec['reference'];u,p=m.f.fields(raw,np.array(ref['points']),np.array(ref['times']))
    if not (np.allclose(u,ref['velocity'],atol=1e-8,rtol=1e-8) and np.allclose(p,ref['pressure'],atol=1e-8,rtol=1e-8)):raise ValueError('Frozen u/p reference mismatch')
    return m.f,raw

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9175291);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Output must be empty')
    a.out.mkdir(parents=True,exist_ok=True);f,raw=reconstruct();fn=a.out/'candidate.json';f.save(raw,fn,dict(recipe='ST052-M',scope='Frozen mathematical field; regenerated metadata can change file SHA'))
    if a.structure:
        from audit_st051 import audit
        from replay_st051 import parent_field as old_parent
        f0,r0=parent_field();parent=a.out/'ST051-B.json';f0.save(r0,parent)
        fs,rs=old_parent('ST048-S');sf=a.out/'ST048-S.json';fs.save(rs,sf)
        audit(fn,parent,sf,a.out/'structure.json',seed=9175293)
    if a.validate:
        from validate import validate
        result=validate(fn,a.out/'validation.json',a.seed);bad=[k for k,v in result['gates'].items() if not v]
        print('FAIL: '+', '.join(bad) if bad else 'Reported gates pass; no continuum or source-identity certificate')
        return int(bool(bad))
    return 0
if __name__=='__main__':raise SystemExit(main())
