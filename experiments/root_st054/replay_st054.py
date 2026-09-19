"""Frozen full-pressure completion replay; no fitting, no free force."""
from __future__ import annotations
import argparse,hashlib,json,sys
from functools import lru_cache
from pathlib import Path
import numpy as np
from pressure_completion import ROOT,Family
HERE=Path(__file__).resolve().parent
PARENT_HASHES={'ST052-M':{'2f815d641e576b24d9896506199efc971b70dd09e33efdefe35470fd3aafc0e5','e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da'},'ST053-Q':{'07a41ca9b45481673b32624a4655f3733f8487b0c6b4b75aa8d8c51161e65a8c','cb7a85c95d13fffcd8265c441fe156c26c3b53722807ea5f4374adc566b2eb29'}}
@lru_cache(maxsize=2)
def parent_field(ident):
    if ident not in PARENT_HASHES:raise ValueError('Unknown parent identity')
    p=ROOT/f'artifacts/research/{ident}/candidate.json'
    if p.is_file():
        if hashlib.sha256(p.read_bytes()).hexdigest() not in PARENT_HASHES[ident]:raise ValueError('Parent checksum mismatch')
        f,r=Family.load(p)
    else:
        sys.path.insert(0,str(ROOT/('experiments/root_st052' if ident=='ST052-M' else 'experiments/root_st053')))
        if ident=='ST052-M':
            from replay_st052 import reconstruct as prior
        else:
            from replay_st053 import reconstruct as prior
        f,r=prior()
    r.setflags(write=False)
    return f,r

def reconstruct(ident,records=None):
    records=json.loads((HERE/'recipes.json').read_text()) if records is None else records
    rec=records[ident]
    expected={'ST054-M3':'ST052-M','ST054-Q2':'ST053-Q'}
    if rec.get('parent_id')!=expected.get(ident):raise ValueError('Wrong parent identity')
    if rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Unsupported scientific claim')
    p=(HERE/rec['delta_file']).resolve()
    if not p.is_relative_to((HERE/'pressure_deltas').resolve()):raise ValueError('Pressure path escape')
    payload=p.read_bytes()
    if hashlib.sha256(payload).hexdigest()!=rec['delta_sha256']:raise ValueError('Pressure checksum mismatch')
    import io
    dq=np.load(io.BytesIO(payload),allow_pickle=False);f,base=parent_field(rec['parent_id'])
    if dq.shape!=(f.n,) or not np.isfinite(dq).all():raise ValueError('Malformed pressure coefficients')
    raw=base.copy();raw[2*f.n:3*f.n]+=dq
    if np.max(abs(raw[2*f.n:3*f.n]))>100:raise ValueError('Pressure bound violation')
    ref=rec['reference'];u,p=f.fields(raw,np.array(ref['points']),np.array(ref['times']))
    if not (np.allclose(u,ref['velocity'],atol=1e-8,rtol=1e-8) and np.allclose(p,ref['pressure'],atol=1e-8,rtol=1e-8)):raise ValueError('Reference mismatch')
    np.testing.assert_array_equal(raw[:2*f.n],base[:2*f.n]);np.testing.assert_array_equal(raw[-2:],base[-2:])
    return f,raw

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9175491);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Output directory is not empty')
    a.out.mkdir(parents=True,exist_ok=True);f,raw=reconstruct(a.id);path=a.out/'candidate.json';f.save(raw,path,dict(id=a.id,scope='Frozen mathematical field; regenerated JSON metadata may differ'))
    if a.structure:
        from audit_pressure import audit
        rec=json.loads((HERE/'recipes.json').read_text())[a.id];fb,b=parent_field(rec['parent_id']);parent=a.out/'parent.json';fb.save(b,parent);audit(parent,path,a.out/'structure.json')
    if a.validate:
        from validate import validate
        r=validate(path,a.out/'validation.json',a.seed);bad=[k for k,v in r['gates'].items() if not v]
        print('FAIL: '+', '.join(bad) if bad else 'Sampled numerical gates pass; no continuum certificate')
        return int(bool(bad))
    return 0
if __name__=='__main__':raise SystemExit(main())
