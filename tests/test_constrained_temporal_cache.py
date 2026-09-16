"""The optimization cache must reproduce the original nonlinear FD residual."""
from pathlib import Path
import numpy as np
from openai_ns_reconstruction.constrained_outer_momentum import AngularMomentumCandidate
from openai_ns_reconstruction.constrained_temporal_swirl import TemporalSwirlCandidate
from openai_ns_reconstruction.constrained_temporal_cache import cached_residual
from openai_ns_reconstruction.constrained_optimize import training_residual


def test_cached_quadratic_residual_matches_direct():
    root=Path(__file__).resolve().parents[1]
    parent=AngularMomentumCandidate.load(root/'artifacts/constrained/continued_pressure/candidate.json')
    zero=TemporalSwirlCandidate(parent)
    rng=np.random.default_rng(82);x=rng.uniform(-1.7,1.7,(48,3));t=rng.uniform(.26,.74,48)
    evaluate=cached_residual(zero,x,t)
    for coefficients in (np.zeros(15),rng.uniform(-.9,.9,15)):
        trial=TemporalSwirlCandidate(parent,tuple(coefficients))
        expected=training_residual(trial,parent.force,x,t,.01)
        np.testing.assert_allclose(evaluate(coefficients),expected,atol=1e-9,rtol=1e-9)

def test_cache_parameter_derivatives():
    from openai_ns_reconstruction.constrained_temporal_cache import cached_energy
    root=Path(__file__).resolve().parents[1]
    parent=AngularMomentumCandidate.load(root/'artifacts/constrained/continued_pressure/candidate.json')
    zero=TemporalSwirlCandidate(parent);rng=np.random.default_rng(29)
    x=rng.uniform(-1.7,1.7,(32,3));t=rng.uniform(.26,.74,32)
    a=rng.uniform(-.5,.5,15);direction=rng.uniform(-1,1,15);h=1e-6
    for evaluate in (cached_residual(zero,x,t),cached_energy(zero,[.25,.5,.75])):
        numerical=(evaluate(a+h*direction)-evaluate(a-h*direction))/(2*h)
        np.testing.assert_allclose(evaluate.jacobian(a)@direction,numerical,atol=1e-8,rtol=1e-7)
    energy=cached_energy(zero,[.25,.5,.75])
    trial=TemporalSwirlCandidate(parent,tuple(a))
    np.testing.assert_allclose(energy(a),[trial.energy(t,24) for t in [.25,.5,.75]],atol=1e-12)
