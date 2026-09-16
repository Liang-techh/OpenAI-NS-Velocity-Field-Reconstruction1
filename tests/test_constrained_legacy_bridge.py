import numpy as np
from openai_ns_reconstruction.constrained_candidate import CompactCandidate
from openai_ns_reconstruction.local_field import LocalField, sum_vector_fields


def test_legacy_curl_path_and_composition_reproduce_candidate():
    c=CompactCandidate(swirl_radial_time=1.2,swirl_axial_time=-.7,
                       radial_shape=.1,axial_shape=-.1,
                       poloidal_radial_time=1.1,poloidal_axial_time=-.9).normalized()
    legacy=c.as_legacy_local_field()
    points=np.random.default_rng(45).uniform(-1.2,1.2,(20,3))
    points=np.vstack((points,[0,0,.3],[2,0,0]))
    actual=np.array([legacy.velocity(*p,.5,eps=1e-5) for p in points])
    np.testing.assert_allclose(actual,c.velocity(points,.5),atol=2e-8,rtol=2e-6)
    # Legacy correction assembly: splitting a potential into halves must
    # preserve its curl, including the axis. This exercises the real old API.
    half=lambda x,y,z,t:0.5*c.vector_potential([x,y,z],t)
    composed=LocalField(sum_vector_fields([half,half]),legacy.azimuthal_scalar)
    np.testing.assert_allclose(composed.velocity(.2,.1,.3,.5),
                               c.velocity([.2,.1,.3],.5),atol=2e-8)
