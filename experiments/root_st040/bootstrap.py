"""Locate immutable published parent and reuse existing research evaluator."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
LEGACY = ROOT / 'experiments/root_st030'
if not LEGACY.is_dir():
    raise RuntimeError('Run from a complete repository checkout or the supplied archive')
sys.path.insert(0, str(LEGACY))
PARENT = ROOT / 'artifacts/research/ST006/candidate.json'
