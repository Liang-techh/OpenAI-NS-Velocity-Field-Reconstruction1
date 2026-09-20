"""ST059: original independent Cartesian validator with honest scientific exit status."""
import sys,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/root_st030'))
from validate import validate
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True,type=Path);p.add_argument('--seed',required=True,type=int);a=p.parse_args()
 if a.out.exists():p.error('Refusing to overwrite independent evidence')
 a.out.parent.mkdir(parents=True,exist_ok=True);r=validate(a.candidate,a.out,a.seed);failed=[k for k,v in r['gates'].items() if not v];print('FAIL: '+', '.join(failed) if failed else 'Sampled gates pass, no continuum certificate')
 raise SystemExit(int(bool(failed)))
