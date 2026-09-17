from pathlib import Path
import numpy as np
from openai_ns_reconstruction.constrained_axial_swirl import AxialSwirlCandidate
from openai_ns_reconstruction.constrained_quintic_swirl import QuinticSwirlCandidate,elevate_cubic


def test_degree_elevation_preserves_velocity_and_time_derivative(tmp_path):
    root=Path(__file__).resolve().parents[1]
    old=AxialSwirlCandidate.load(root/'artifacts/constrained/axial_swirl_grid48/candidate.json')
    new=QuinticSwirlCandidate(old.parent,elevate_cubic(old.coefficients))
    rng=np.random.default_rng(57);x=rng.uniform(-1.7,1.7,(64,3));t=rng.uniform(.252,.748,64)
    np.testing.assert_allclose(new.velocity(x,t),old.velocity(x,t),atol=1e-14)
    h=1e-5
    np.testing.assert_allclose((new.velocity(x,t+h)-new.velocity(x,t-h))/(2*h),(old.velocity(x,t+h)-old.velocity(x,t-h))/(2*h),atol=1e-10)
    path=tmp_path/'candidate.json';new.save(path)
    np.testing.assert_allclose(QuinticSwirlCandidate.load(path).velocity(x,t),new.velocity(x,t),atol=1e-14)
    assert max(abs(v) for v in new.coefficients)<=1
