import math
import numpy as np
import pytest
from openai_ns_reconstruction.paper_core_reference import PaperCoreReference
from openai_ns_reconstruction.coordinates import similarity_coordinates


def test_source_series_axis_and_coordinate_identity():
    f = PaperCoreReference()
    X,e = .03,.2
    H=(.5-f.h)*e+(1-e*e)*(4*e+f.j)
    s=f.Lambda*X*H*H/(H*H+f.sigma**2)
    series=sum((-s/2)**n/(math.factorial(n)*math.factorial(n+1)) for n in range(20))
    assert f.F(X,e)/f.amplitude(e) == pytest.approx(series, abs=1e-14)
    chart=similarity_coordinates(.1,.2,.5,f.h)
    assert chart.q-.2**2*chart.q**(2*f.h) == pytest.approx(.5,abs=1e-12)
    axis=f.velocity(0,0,.2,.5)
    c=similarity_coordinates(0,.2,.5,f.h)
    np.testing.assert_allclose(axis,[0,0,c.q**(-c.A)*(4*c.eta+f.j)])
    with pytest.raises(ValueError): f.velocity(3,0,0,.5)


def test_derivative_and_incompressibility():
    f=PaperCoreReference()
    step=1e-5
    for e in [-.6,0,.3]:
        derivative=(f.U(.03,e+step)-f.U(.03,e-step))/(2*step)
        assert f.dU(.03,e) == pytest.approx(derivative,rel=1e-8)
    for point in [[.1,.02,.1],[.2,-.1,-.2],[0,0,0]]:
        p=np.array(point,dtype=float)
        div=0.
        for i in range(3):
            d=np.eye(3)[i]*step
            div+=(f.velocity(*(p+d),.5)[i]-f.velocity(*(p-d),.5)[i])/(2*step)
        assert abs(div)<2e-7
