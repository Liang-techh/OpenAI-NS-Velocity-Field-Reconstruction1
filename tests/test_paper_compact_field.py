import numpy as np
from openai_ns_reconstruction.paper_compact_field import PaperCompactField


def test_preserves_inner_and_exact_zero_exterior():
    f=PaperCompactField()
    for p in [[.1,0,.1],[0,0,.2],[-.2,.1,-.1]]:
        np.testing.assert_allclose(f.velocity(*p,.5),f.core.velocity(*p,.5),rtol=1e-12,atol=1e-12)
    for p in [[2,0,0],[0,0,-2],[1.5,0,0],[1e200,1e200,1e200]]:
        np.testing.assert_array_equal(f.velocity(*p,.5),[0,0,0])


def test_divergence_in_radial_and_axial_transition():
    f=PaperCompactField();step=2e-5
    for p in [[.56,.05,.1],[.4,.3,.7],[.6,0,1.5],[0,0,.8]]:
        p=np.array(p);div=0.
        for i in range(3):
            d=np.eye(3)[i]*step
            div+=(f.velocity(*(p+d),.5)[i]-f.velocity(*(p-d),.5)[i])/(2*step)
        assert abs(div)<1e-5
