from pathlib import Path
import numpy as np
import pytest
from openai_ns_reconstruction.constrained_localized_swirl import LocalizedSwirlCandidate
from openai_ns_reconstruction.constrained_axial_swirl import AxialSwirlCandidate
from openai_ns_reconstruction.constrained_momentum_budget import angular_moment


def test_nested_basis_and_roundtrip(tmp_path):
    root=Path(__file__).resolve().parents[1]
    old=LocalizedSwirlCandidate.load(root/'artifacts/constrained/localized_swirl/candidate.json')
    c=AxialSwirlCandidate(old.parent,old.coefficients+(0.,)*18)
    x=np.random.default_rng(7).uniform(-1.5,1.5,(32,3));t=np.linspace(.25,.75,32)
    np.testing.assert_allclose(c.velocity(x,t),old.velocity(x,t),atol=1e-14)
    full=AxialSwirlCandidate(old.parent,tuple(np.linspace(-.5,.5,36)))
    np.testing.assert_allclose(full.velocity(x,.25),old.parent.velocity(x,.25),atol=1e-14)
    assert abs(angular_moment(lambda p,t:full.velocity(p,t)-old.parent.velocity(p,t),.7,96))<1e-10
    path=tmp_path/'candidate.json';full.save(path)
    np.testing.assert_allclose(AxialSwirlCandidate.load(path).velocity(x,t),full.velocity(x,t),atol=1e-14)
    with pytest.raises(ValueError):AxialSwirlCandidate(old.parent,(2.,)*36)
