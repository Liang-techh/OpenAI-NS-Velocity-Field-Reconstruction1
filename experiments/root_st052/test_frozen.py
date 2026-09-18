"""Check actual frozen field, and reject false claims/corrupt recipes."""
from pathlib import Path
import copy,json
import numpy as np
import pytest
from replay_st052 import reconstruct

def test_frozen_recipe_and_mutations():
    p=Path(__file__).with_name('recipe.json');rec=json.loads(p.read_text())
    f,r=reconstruct(rec)
    assert np.isfinite(r).all()
    bad=copy.deepcopy(rec);bad['pde_validated']=True
    with pytest.raises(ValueError,match='scientific claim'):reconstruct(bad)
    bad=copy.deepcopy(rec);bad['modifiers_sha256']='0'*64
    with pytest.raises(ValueError,match='checksum'):reconstruct(bad)
    bad=copy.deepcopy(rec);bad['parent_id']='invented'
    with pytest.raises(ValueError,match='Parent identity'):reconstruct(bad)
