import numpy as np
import pytest
from openai_ns_reconstruction.velocity_components import VelocityField,velocity,u,v,w


def test_components_axis_exterior_and_broadcast():
    field=VelocityField();points=np.array([[.1,0,.1],[0,0,.2],[.3,-.2,.4]])
    actual=field.at_points(points,.5)
    np.testing.assert_allclose(actual,field.candidate.velocity(points,.5),atol=1e-14)
    assert all(isinstance(a,float) for a in velocity(.1,0,.1,.5))
    np.testing.assert_allclose((u(.1,0,.1,.5),v(.1,0,.1,.5),w(.1,0,.1,.5)),actual[0])
    assert velocity(1e200,0,0,.5)==(0.,0.,0.)
    assert field.grid([0,.1],[0],[.1],[.25,.5]).shape==(2,2,1,1,3)
    for time in [.24,.76,float('nan')]:
        with pytest.raises(ValueError):field.components(0,0,0,time)


def test_chunking_keeps_point_order():
    field=VelocityField();x=np.linspace(-1,1,4200)
    actual=field.at_points(np.column_stack((x,np.zeros(len(x)),np.full(len(x),.1))),.5)
    indices=[0,100,4095,4096,4199]
    for i in indices:np.testing.assert_allclose(actual[i],field.components(x[i],0,.1,.5),atol=1e-14)
