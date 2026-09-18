"""Replay only: checked parent plus frozen modifiers, not new fitting."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,sys,tempfile
from pathlib import Path
import numpy as np
import continuation
from localized_model import LocalizedModel
from spacetime import Family
from validate import validate
from structure_and_particles import diagnostics
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def reconstruct(ident,records=None):
    records=json.loads((HERE/'recipes.json').read_text()) if records is None else records
    rec=records[ident]
    if rec.get('parent_id')!='ST046-A' or rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:
        raise ValueError('Unsupported parent or scientific claim')
    data=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(data).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Modifier checksum')
    c=np.load(io.BytesIO(data),allow_pickle=False)
    local=ROOT/'artifacts/research/ST046-A/candidate.json'
    if local.exists():
        if hashlib.sha256(local.read_bytes()).hexdigest()!=rec['parent_local_reconstruction_sha256']:raise ValueError('Local reconstructed parent checksum')
        f,raw=Family.load(local)
    else:
        # On GitHub use the existing #366 reconstruction; do not duplicate it.
        sys.path.insert(0,str(ROOT/'experiments/root_st046'))
        from replay_st046 import reconstruct as parent_reconstruct
        f,raw=parent_reconstruct('ST046-A')
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'parent.json';f.save(raw,path)
        m=LocalizedModel(path,rec['radial'],rec['axial'],rec['time_degree'])
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Malformed modifiers')
    lo,hi=np.asarray(m.bounds).T
    if np.any((c<lo-1e-10)|(c>hi+1e-10)):raise ValueError('Modifier bounds')
    raw=m.candidate(c);r=rec['reference'];u,p=m.f.fields(raw,np.array(r['points']),np.array(r['times']))
    if not(np.allclose(u,r['velocity'],atol=1e-8,rtol=1e-8) and np.allclose(p,r['pressure'],atol=1e-8,rtol=1e-8)):
        raise ValueError('Frozen field reference mismatch')
    return m.f,raw

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9174701);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Output must be empty')
    a.out.mkdir(parents=True,exist_ok=True);f,r=reconstruct(a.id);path=a.out/'candidate.json'
    f.save(r,path,{'recipe':a.id,'scope':'Mathematical replay; raw JSON metadata/hash differ from original frozen file'})
    if a.structure:diagnostics(path,a.out)
    if a.validate:
        result=validate(path,a.out/'validation.json',a.seed);failed=[k for k,v in result['gates'].items() if not v]
        print('FAIL: '+', '.join(failed) if failed else 'All reported numerical gates pass; not a continuum proof')
        return int(bool(failed))
    return 0
if __name__=='__main__':raise SystemExit(main())
