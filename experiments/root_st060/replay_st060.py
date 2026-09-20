"""Reproduce a frozen ST060 field and its diagnostics without optimization."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st059r'))
from checkpoint import atomic_json,atomic_bytes,filehash
from validate_task import execute

def main():
    p=argparse.ArgumentParser();p.add_argument('--id',default='ST060-Q');p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seed',type=int,default=9206191);p.add_argument('--validate',action='store_true');p.add_argument('--structure',action='store_true');p.add_argument('--resume',action='store_true');a=p.parse_args()
    freeze=json.loads((ROOT/'evidence_st060/freeze.json').read_text())
    if a.id not in freeze['candidates']:p.error('Candidate not in the frozen verification set')
    item=freeze['candidates'][a.id];original=ROOT/item['path']
    if filehash(original)!=item['sha256']:raise ValueError('Frozen candidate bytes changed')
    a.out.mkdir(parents=True,exist_ok=True);copy=a.out/'candidate.json'
    if copy.exists():
        if not a.resume:raise FileExistsError('Existing output; use --resume to verify it')
        if filehash(copy)!=item['sha256']:raise ValueError('Output candidate changed')
    else:atomic_bytes(copy,original.read_bytes())
    if a.structure:
        from verify_round import audit
        audit(copy,ROOT/freeze['candidates']['ST059R-V']['path'],a.out/'diagnostics')
    code=0
    if a.validate:
        r=execute(copy,a.out/'validation.json',a.seed,a.resume);code=r['scientific_exit_code']
    atomic_json(a.out/'replay_receipt.json',dict(candidate_id=a.id,candidate_sha256=item['sha256'],seed=a.seed,scientific_exit_code=code if a.validate else None,optimizer_executed=False,pde_validated=False))
    print('Frozen fields replayed; see original scientific gate result, not a continuum certificate.')
    return code
if __name__=='__main__':raise SystemExit(main())
