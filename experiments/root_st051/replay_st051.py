"""Replay frozen modifiers only. No optimizer or external data fetch is executed."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,tempfile
from pathlib import Path
import numpy as np
from aligned_continuation import EdgeModel,Family,ROOT
HERE=Path(__file__).resolve().parent
PARENT_HASH='f5b5a2869991387e503a522b599af86ea814ba02bdae2f8b56fd16674f7dab19'
S_HASHES={'13127759debcf40b3f3d6c3627747ace67f1471b825b2b0e8c53ae79a9911344','6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09'}

def parent_field(ident='ST050R-P'):
    if ident not in ('ST050R-P','ST048-S'):raise ValueError('Unsupported parent')
    path=ROOT/f'artifacts/research/{ident}/candidate.json'
    if path.is_file():
        expected={PARENT_HASH} if ident=='ST050R-P' else S_HASHES
        if hashlib.sha256(path.read_bytes()).hexdigest() not in expected:raise ValueError('Parent hash mismatch')
        return Family.load(path)
    from replay_recovery import reconstruct,parent_field as recover_S
    return recover_S() if ident=='ST048-S' else reconstruct(ident)

def reconstruct(ident,records=None):
    rec=json.loads((HERE/'recipes'/f'{ident}.json').read_text()) if records is None else records[ident]
    if rec.get('parent_id')!='ST050R-P' or rec.get('parent_original_raw_sha256')!=PARENT_HASH:raise ValueError('Parent identity mismatch')
    if rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Unsupported scientific claim')
    payload=base64.b64decode(rec['modifiers_npy_b64'],validate=True)
    if hashlib.sha256(payload).hexdigest()!=rec['modifiers_sha256']:raise ValueError('Modifier hash mismatch')
    c=np.load(io.BytesIO(payload),allow_pickle=False);f,raw=parent_field()
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/'parent.json';f.save(raw,path);m=EdgeModel(path,edge_modes=rec['edge_modes'])
    if c.shape!=(m.dim,) or not np.isfinite(c).all():raise ValueError('Invalid modifiers')
    lo,hi=np.array(m.bounds).T
    if np.any((c<lo-1e-10)|(c>hi+1e-10)):raise ValueError('Modifier bound violation')
    raw=m.candidate(c);ref=rec['reference'];u,p=m.f.fields(raw,np.array(ref['points']),np.array(ref['times']))
    if not (np.allclose(u,ref['velocity'],atol=1e-8,rtol=1e-8) and np.allclose(p,ref['pressure'],atol=1e-8,rtol=1e-8)):raise ValueError('Frozen field reference mismatch')
    return m.f,raw

def main():
    a=argparse.ArgumentParser();a.add_argument('--id',required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--seed',type=int,default=9175191);a.add_argument('--validate',action='store_true');a.add_argument('--structure',action='store_true');args=a.parse_args()
    out=args.out
    if out.exists() and any(out.iterdir()):a.error('Refusing nonempty output directory')
    out.mkdir(parents=True,exist_ok=True);f,r=reconstruct(args.id);path=out/'candidate.json'
    f.save(r,path,dict(recipe=args.id,scope='Reconstructed mathematical field; JSON metadata differs from original local raw candidate'))
    if args.structure:
        from audit_st051 import audit
        paths=[]
        for ident in ('ST050R-P','ST048-S'):
            f0,r0=parent_field(ident);p=out/f'{ident}_parent.json';f0.save(r0,p);paths.append(p)
        audit(path,*paths,out/'structure.json')
    if args.validate:
        from validate import validate
        rep=validate(path,out/'validation.json',args.seed)
        failed=[k for k,v in rep['gates'].items() if not v]
        print('FAIL: '+', '.join(failed) if failed else 'Reported numerical gates pass, not a continuum certificate')
        return int(bool(failed))
    return 0
if __name__=='__main__':raise SystemExit(main())
