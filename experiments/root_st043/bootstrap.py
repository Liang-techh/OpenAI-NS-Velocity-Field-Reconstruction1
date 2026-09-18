from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
for name in ('root_st030','root_st040'):
    sys.path.insert(0,str(ROOT/'experiments'/name))
PARENT=ROOT/'artifacts/research/ST006/candidate.json'
