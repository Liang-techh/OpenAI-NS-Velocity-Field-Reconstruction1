import numpy as np
from openai_ns_reconstruction.constrained_candidate import CompactCandidate
from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_optimize import training_residual
from openai_ns_reconstruction.constrained_variable_projection import fit_linear_coefficients


def test_bounded_linear_fit_improves_fixed_velocity_residual():
    rng=np.random.default_rng(291)
    x=rng.uniform(-1.5,1.5,(160,3));t=rng.uniform(.3,.7,len(x))
    c=CompactCandidate().normalized()
    before=training_residual(c,RestrictedForce(),x,t,.01)
    fitted,force=fit_linear_coefficients(c,x,t,.01)
    after=training_residual(fitted,force,x,t,.01)
    assert np.linalg.norm(after)<np.linalg.norm(before)
    np.testing.assert_array_equal(fitted.velocity(x,t),c.velocity(x,t))
    assert 0<=force.a<=10 and 0<=force.c<=10
    assert all(abs(getattr(fitted,k))<=100 for k in
               ('pressure_constant','pressure_radial','pressure_axial'))
