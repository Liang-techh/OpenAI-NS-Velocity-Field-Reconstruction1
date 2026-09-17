"""Expose the isolated full-momentum research regressions to the constrained CI lane."""
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1] / 'experiments' / 'root_st001'
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('root_st001_regressions', ROOT / 'test_spacetime.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
for name in dir(module):
    if name.startswith('test_'):
        globals()[name] = getattr(module, name)
