import numpy as np
from openai_ns_reconstruction.constrained_pressure_circulation import loop_circulation


def test_polynomial_circulation():
    def velocity(p,t):return np.stack((p[:,2],np.zeros(len(p)),p[:,0]**2),axis=-1)
    def force(p,t):return np.zeros_like(p)
    # Acceleration=(r²,0,2rz), so its meridional curl is2z.
    value=loop_circulation(velocity,force,[.2,.8,.1,.7],time=.5,order=8,step=.001)
    expected=(.8-.2)*(.7**2-.1**2)
    np.testing.assert_allclose(value,expected,atol=1e-8)
