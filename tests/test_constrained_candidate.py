import numpy as np
import pytest

from openai_ns_reconstruction.constrained_candidate import CompactCandidate


def test_divergence_axis_support_and_core_signs():
    c = CompactCandidate(radial_shape=0.13, axial_shape=-0.17).normalized()
    points=np.random.default_rng(25).uniform(-1.5,1.5,(80,3))
    h=1e-5
    div=sum((c.velocity(points+np.eye(3)[j]*h,0.4)[:,j]
             -c.velocity(points-np.eye(3)[j]*h,0.4)[:,j])/(2*h) for j in range(3))
    assert np.max(np.abs(div)) < 2e-7
    np.testing.assert_array_equal(c.velocity([[2,0,0],[0,0,2]],0.25),0)
    np.testing.assert_array_equal(c.pressure([[2,0,0],[0,0,2]],0.25),0)
    assert np.isfinite(c.velocity([0,0,0.2],0.5)).all()
    for t in (0.25,0.5,0.75):
        u=c.velocity([0.1*np.sqrt(1-t),0,0.1*(1-t)**0.495],t)
        assert u[0]<0<u[1] and u[2]>0


def test_energy_normalization_scaling_and_artifact(tmp_path):
    c=CompactCandidate().normalized()
    assert abs(c.energy(order=128)-1)<1e-8
    profiles=[]
    for t in (0.25,0.5,0.75):
        tau=1-t
        u=c.velocity([0.1*np.sqrt(tau),0,0.1*tau**0.495],t)
        profiles.append(u*np.array([np.sqrt(tau),tau**0.505,tau**0.505]))
    np.testing.assert_allclose(profiles,[profiles[0]]*3,atol=1e-13)
    path=tmp_path/'candidate.json'
    c.save(path)
    assert CompactCandidate.load(path)==c


def test_invalid_parameters_and_time_rejected():
    with pytest.raises(ValueError):
        CompactCandidate(swirl_ratio=0)
    with pytest.raises(ValueError):
        CompactCandidate(radial_width=2)
    with pytest.raises(ValueError):
        CompactCandidate().velocity([0,0,0],1)
