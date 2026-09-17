import numpy as np
from openai_ns_reconstruction.eq45_supported_delivery import default_field
from openai_ns_reconstruction.constrained_eq45_bipolar_seed import bipolar_seed
from openai_ns_reconstruction.constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


def test_central_motion_parity_and_serialization(tmp_path):
    base=default_field().candidate;c=bipolar_seed(base)
    for t in [.25,.5,.75]:
        up=c.velocity_xyz(.1,0,.1,t);down=c.velocity_xyz(.1,0,-.1,t)
        assert up[0]<0 and up[1]>0 and up[2]>0 and down[2]<0
        np.testing.assert_allclose(down,up*[1,1,-1],atol=1e-14)
        assert c.velocity_xyz(.1,0,0,t)[2]==0
        np.testing.assert_array_equal(c.velocity_xyz(2,0,0,t),[0,0,0])
    assert c.parent.profile_basis.swirl_coefficients==base.parent.profile_basis.swirl_coefficients
    c.save_json(tmp_path/'candidate.json')
    loaded=Eq45SupportedVelocityCandidate.load_json(tmp_path/'candidate.json')
    assert loaded.sha256==c.sha256
    assert not loaded.to_dict()['truth_boundary']['pde_validated']
