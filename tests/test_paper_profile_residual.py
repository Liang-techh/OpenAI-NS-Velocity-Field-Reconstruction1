import math
import numpy as np
from scipy.integrate import quad
from scipy.special import hyp1f1
from openai_ns_reconstruction.profiles import LeadingProfile
from openai_ns_reconstruction.paper_profile_residual import residual


def test_exact_equatorial_swirl_equation():
    # At eta=0, U=0, the angular equation is the Kummer ODE:
    # 2X F''+(4-X)F'-(1+h)F=0. This is an independent
    # nonzero check of signs, radial factors and pressure differentiation.
    h=.005
    F=lambda X,e: float(hyp1f1(1+h,2,X/2))
    p=LeadingProfile(E=lambda X,e: math.sqrt(2*X)*F(X,e),
                     U=lambda X,e:0.,dU_deta=lambda X,e:0.,F=F,
                     Pi=lambda X,e:quad(lambda s:F(s,e)**2,0,X)[0])
    for X in [.03,.1,.3]:
        np.testing.assert_allclose(residual(p,X,0,h,1e-4),0,atol=2e-7)

