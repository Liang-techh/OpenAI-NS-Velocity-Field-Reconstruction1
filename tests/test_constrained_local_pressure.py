from pathlib import Path
import numpy as np
from openai_ns_reconstruction.constrained_inner_swirl import InnerSwirlCandidate
from openai_ns_reconstruction.constrained_local_pressure import LocalPressureCandidate


def test_local_pressure_preserves_velocity_and_compact_support(tmp_path):
    root=Path(__file__).resolve().parents[1]
    base=InnerSwirlCandidate.load(root/'artifacts/constrained/inner_swirl_pressure/candidate.json')
    c=LocalPressureCandidate(base,tuple(np.linspace(-.5,.5,27)))
    p=np.array([[.4,.2,.5],[0.,0.,0.],[2.1,0.,0.],[0.,0.,2.1]])
    np.testing.assert_array_equal(c.velocity(p,.6),base.velocity(p,.6))
    assert np.max(np.abs(c.pressure_basis(p[2:],.6)))==0
    assert np.linalg.norm(c.pressure(p[:2],.6)-base.pressure(p[:2],.6))>0
    np.testing.assert_allclose(c.pressure(p,.6),c.pressure(p*np.array([1,1,-1]),.6),atol=1e-14)
    path=tmp_path/'candidate.json';c.save(path)
    np.testing.assert_allclose(LocalPressureCandidate.load(path).pressure(p,.6),c.pressure(p,.6),atol=1e-14)
