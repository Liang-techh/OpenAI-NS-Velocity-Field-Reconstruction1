from pathlib import Path
import numpy as np
import pytest
from openai_ns_reconstruction.constrained_local_pressure import LocalPressureCandidate
from openai_ns_reconstruction.constrained_poloidal import PoloidalCandidate


@pytest.mark.parametrize("anchor_core",[False,True])
def test_streamfunction_divergence_core_and_initial(tmp_path,anchor_core):
    root=Path(__file__).resolve().parents[1];old=LocalPressureCandidate.load(root/'artifacts/constrained/local_pressure/candidate.json')
    c=PoloidalCandidate(old.base,old.pressure_coefficients,tuple(np.linspace(-.1,.1,27)),anchor_core)
    x=np.array([[.4,.2,.5],[.7,-.3,.8],[.15,.05,.4]])
    np.testing.assert_allclose(c.velocity(x,.25),old.velocity(x,.25),atol=1e-14)
    t=np.linspace(.25,.75,21);tau=1-t;core=np.column_stack((.1*np.sqrt(tau),np.zeros(21),.1*tau**.495))
    np.testing.assert_allclose(c.velocity(core,t),old.velocity(core,t),atol=1e-14)
    assert np.linalg.norm(c.velocity(x,.6)-old.velocity(x,.6))>.01
    h=1e-5;div=np.zeros(len(x))
    for j in range(3):
        d=np.eye(3)[j]*h;div+=(c.velocity(x+d,.6)[:,j]-c.velocity(x-d,.6)[:,j])/(2*h)
    np.testing.assert_allclose(div,0,atol=1e-7)
    p=tmp_path/'candidate.json';c.save(p)
    np.testing.assert_allclose(PoloidalCandidate.load(p).velocity(x,.6),c.velocity(x,.6),atol=1e-14)
