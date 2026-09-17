from pathlib import Path
import numpy as np
from openai_ns_reconstruction.constrained_quintic_swirl import QuinticSwirlCandidate
from openai_ns_reconstruction.constrained_inner_swirl import InnerSwirlCandidate
from openai_ns_reconstruction.constrained_momentum_budget import angular_moment


def test_inner_modes_reach_collar_but_preserve_core_and_moment(tmp_path):
    root=Path(__file__).resolve().parents[1]
    old=QuinticSwirlCandidate.load(root/'artifacts/constrained/quintic_swirl/candidate.json')
    embedded=InnerSwirlCandidate(old.parent,old.coefficients+(0.,)*30)
    x=np.array([[.32,0,.8857],[.18,.12,.4]])
    np.testing.assert_allclose(embedded.velocity(x,.7),old.velocity(x,.7),atol=1e-14)
    c=InnerSwirlCandidate(old.parent,old.coefficients+(.4,)*30)
    assert np.linalg.norm(c.velocity(x,.7)-old.velocity(x,.7))>1e-3
    times=np.linspace(.25,.75,21);tau=1-times
    core=np.column_stack((.1*np.sqrt(tau),np.zeros(21),.1*tau**.495))
    np.testing.assert_allclose(c.velocity(core,times),old.velocity(core,times),atol=1e-14)
    np.testing.assert_allclose(c.velocity(x,.25),old.velocity(x,.25),atol=1e-14)
    assert abs(angular_moment(lambda p,t:c.velocity(p,t)-old.velocity(p,t),.7,96))<1e-10
    h=1e-5;div=np.zeros(len(x))
    for j in range(3):
        d=np.eye(3)[j]*h;div+=(c.velocity(x+d,.7)[:,j]-c.velocity(x-d,.7)[:,j])/(2*h)
    np.testing.assert_allclose(div,0,atol=1e-7)
    p=tmp_path/'candidate.json';c.save(p)
    np.testing.assert_allclose(InnerSwirlCandidate.load(p).velocity(x,.7),c.velocity(x,.7),atol=1e-14)
