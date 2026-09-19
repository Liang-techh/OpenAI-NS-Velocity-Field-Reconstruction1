"""Frozen ST053 mathematical field replay; no coefficient fitting."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,sys,tempfile
from pathlib import Path
import numpy as np
from moment_step import ROOT
from aligned_continuation import EdgeModel,Family
HERE=Path(__file__).resolve().parent
PARENT_ORIGINAL='e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da'
PARENT_LOCAL='482471204a2c2217a95dde23017a8a09ff35f4ddd8c5b31b5ae48141edaa37ae'
def parent_field():
    p=ROOT/'artifacts/research/ST052-M/candidate.json'
    if p.is_file():
        if hashlib.sha256(p.read_bytes()).hexdigest() not in (PARENT_ORIGINAL,PARENT_LOCAL):raise ValueError('Wrong parent bytes')
        return Family.load(p)
    sys.path.insert(0,str(ROOT/'experiments/root_st052'))
    from replay_st052 import reconstruct
    return reconstruct()
def reconstruct(rec=None):
    rec=json.loads((HERE/'recipe.json').read_text()) if rec is None else rec
    if rec.get('parent_original_sha256')!=PARENT_ORIGINAL or rec.get('parent_id')!='ST052-M':raise ValueError('Wrong parent identity')
    if rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Unsupported scientific claim')
    payload=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(payload).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Wrong modifier hash')
    c=np.load(io.BytesIO(payload),allow_pickle=False);f,raw=parent_field()
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'parent.json';f.save(raw,p);m=EdgeModel(p,edge_modes=True)
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Malformed modifiers')
    lo,hi=np.array(m.bounds).T
    if np.any((c<lo-1e-10)|(c>hi+1e-10)):raise ValueError('Modifier outside bounds')
    raw=m.candidate(c);ref=rec['reference'];u,p=m.f.fields(raw,np.array(ref['points']),np.array(ref['times']))
    if not(np.allclose(u,ref['velocity'],atol=1e-8,rtol=1e-8) and np.allclose(p,ref['pressure'],atol=1e-8,rtol=1e-8)):raise ValueError('Reference mismatch')
    return m.f,raw

def main():
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);a.add_argument('--seed',type=int,default=9175391)
    a.add_argument('--validate',action='store_true');a.add_argument('--structure',action='store_true');args=a.parse_args()
    if args.out.exists() and any(args.out.iterdir()):a.error('Output not empty')
    args.out.mkdir(parents=True,exist_ok=True);f,raw=reconstruct();candidate=args.out/'candidate.json'
    f.save(raw,candidate,dict(recipe='ST053-Q',scope='Frozen mathematical field; JSON metadata may differ'))
    if args.structure:
        from audit_st051 import audit
        from replay_st051 import parent_field as old_parent
        fp,rp=parent_field();par=args.out/'parent.json';fp.save(rp,par)
        fs,rs=old_parent('ST048-S');ref=args.out/'ST048-S.json';fs.save(rs,ref)
        audit(candidate,par,ref,args.out/'structure.json',seed=9175393)
    if args.validate:
        from validate import validate
        r=validate(candidate,args.out/'validation.json',args.seed);bad=[k for k,v in r['gates'].items() if not v]
        print('FAIL: '+', '.join(bad) if bad else 'Reported numerical gates pass; no continuum certificate')
        return int(bool(bad))
    return 0
if __name__=='__main__':raise SystemExit(main())
