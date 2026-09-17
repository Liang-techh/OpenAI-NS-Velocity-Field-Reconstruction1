"""Non-destructive publication layout check and live checkout inventory."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research_baseline import verify_integrity

def main():
    verify_integrity()
    required=['README.md','docs/README.md','docs/REPOSITORY_GUIDE.md','docs/RESEARCH_STATUS.md','docs/BRANCH_AND_PR_GUIDE.md','artifacts/research/README.md','artifacts/research/experiment_index.json','research_baseline/__init__.py','experiments/root_st030/validate.py','tests/test_published_baseline.py']
    missing=[s for s in required if not (ROOT/s).is_file()]
    if missing:raise ValueError('Missing publication paths: '+', '.join(missing))
    print(json.dumps({'layout':'ok','required_paths':len(required),'top_level_file_counts':{p.name:sum(f.is_file() for f in p.rglob('*') if '.git' not in f.parts and '__pycache__' not in f.parts) for p in ROOT.iterdir() if p.is_dir() and not p.name.startswith('.')},'pde_validated':False},indent=2))
if __name__=='__main__':main()
