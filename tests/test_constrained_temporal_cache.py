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
