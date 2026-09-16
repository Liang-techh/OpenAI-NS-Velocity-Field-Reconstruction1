import numpy as np
from openai_ns_reconstruction.constrained_tensor_candidate import TensorCandidate
from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_outer_momentum import AngularMomentumCandidate


def test_outer_swirl_preserves_core_and_support_and_roundtrip(tmp_path):
    base=TensorCandidate().normalized();c=AngularMomentumCandidate(base,RestrictedForce(),outer_shape=(-1.2,.4),pressure_coefficients=(.2,)*18)
    p=np.array([[.01,0,.01],[0,0,.2],[.1,.1,-.1]])
    np.testing.assert_array_equal(c.velocity(p,.5),base.velocity(p,.5))
    np.testing.assert_array_equal(c.velocity([[2,0,0],[0,0,2]],.5),0)
    out=tmp_path/'outer.json';c.save(out);loaded=AngularMomentumCandidate.load(out)
    np.testing.assert_allclose(c.velocity([.8,.1,.3],.5),loaded.velocity([.8,.1,.3],.5))

    np.testing.assert_allclose(c.pressure([.8,.1,.3],.5),loaded.pressure([.8,.1,.3],.5))
    np.testing.assert_array_equal(c.pressure([[2,0,0],[0,0,2]],.5),0)
