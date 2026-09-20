"""Reconstruct frozen fields from checked increments, without optimization."""
from __future__ import annotations
import argparse,base64,hashlib,io,json,sys,tempfile,zlib
from pathlib import Path
import numpy as np
from joint_step import Family,JointModel,ROOT
HERE=Path(__file__).resolve().parent
PARENT_SHA='d772949621e0f7703a9d6b36f28828532ba3e5ec678dec14db5ce14eb8b851d5'

def parent_field():
    p=ROOT/'artifacts/research/ST054-Q2/candidate.json'
    if p.is_file():
        if hashlib.sha256(p.read_bytes()).hexdigest()!=PARENT_SHA:raise ValueError('Parent checksum mismatch')
        return Family.load(p)
    sys.path.insert(0,str(ROOT/'experiments/root_st054'))
    from replay_st054 import reconstruct
    return reconstruct('ST054-Q2')

def reconstruct(ident,record=None):
    if ident not in ('ST056-J2','ST056-M1'):raise ValueError('Unsupported frozen field')
    rec=json.loads((HERE/f'{ident}.json').read_text()) if record is None else record
    if rec.get('id')!=ident or rec.get('parent')!='ST054-Q2' or rec.get('parent_original_sha256')!=PARENT_SHA:raise ValueError('Parent or field identity mismatch')
    if rec.get('pde_validated') is not False or rec.get('source_correspondence_verified') is not False:raise ValueError('Unsupported scientific claim')
    if rec.get('delta_codec')!='base64(zlib(npy))':raise ValueError('Unsupported codec')
    packed=base64.b64decode(rec['delta'],validate=True)
    decomp=zlib.decompressobj();payload=decomp.decompress(packed,20000)
    if not decomp.eof or decomp.unused_data or decomp.unconsumed_tail:raise ValueError('Invalid or oversized delta')
    if hashlib.sha256(payload).hexdigest()!=rec['delta_npy_sha256']:raise ValueError('Delta checksum mismatch')
    c=np.load(io.BytesIO(payload),allow_pickle=False)
    if c.shape!=(1620,) or not np.isfinite(c).all():raise ValueError('Invalid delta array')
    f,raw=parent_field();gamma=rec['initial_swirl_multiplier']
    expected=1. if ident=='ST056-J2' else .958132801471448
    if gamma!=expected:raise ValueError('Initial-data multiplier mismatch')
    raw=raw.copy();raw[f.n:2*f.n]*=gamma
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'parent.json';f.save(raw,p);m=JointModel(p)
    child=m.candidate(c);refs=rec['references'];u,p=m.f.fields(child,np.array(refs['points']),np.array(refs['times']))
    if not np.allclose(u,refs['velocity'],rtol=1e-8,atol=1e-8) or not np.allclose(p,refs['pressure'],rtol=1e-8,atol=1e-8):raise ValueError('Frozen reference mismatch')
    return m.f,child

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',default='ST056-J2',choices=['ST056-J2','ST056-M1']);p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=int,default=9205691);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):p.error('Refusing nonempty output directory')
    a.out.mkdir(parents=True,exist_ok=True);f,raw=reconstruct(a.id);path=a.out/'candidate.json'
    f.save(raw,path,dict(id=a.id,scope='Frozen mathematical reconstruction; original JSON metadata/hash differ; no scientific acceptance'))
    if a.structure:
        from audit import structure
        f0,r0=parent_field();parent=a.out/'parent.json';f0.save(r0,parent);structure(path,parent,a.out/'structure.json')
    if a.validate:
        from validate import validate
        report=validate(path,a.out/'validation.json',a.seed);failed=[k for k,v in report['gates'].items() if not v]
        print('FAIL: '+', '.join(failed) if failed else 'Sampled gates pass, not a continuum certificate')
        return int(bool(failed))
    return 0
if __name__=='__main__':raise SystemExit(main())
