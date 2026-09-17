import numpy as np
from annular_curl_modes import CurlModes

def test_complete_curl_jets_and_support():
    m=CurlModes(ms=(1,2),radial_degrees=(0,),axial_degrees=(0,1));rng=np.random.default_rng(9172905)
    x=rng.uniform(-1.1,1.1,(15,3));u,j,L,p=m.jets(x)
    assert np.max(abs(np.trace(j,axis1=2,axis2=3)))<3e-13
    h=1e-4;lap=np.zeros_like(u)
    for k in range(3):
        e=np.eye(3)[k]*h;up=m.jets(x+e)[0];um=m.jets(x-e)[0]
        assert np.max(abs((up-um)/(2*h)-j[:,:,:,k]))<2e-5
        lap+=(up-2*u+um)/h**2
    assert np.max(abs(lap-L))<3e-4
    probe=np.array([[.1,0,.1],[0,0,0],[1.9,0,0],[0,1.9,0],[.5,0,1.9]])
    for a in m.jets(probe):assert np.max(abs(a))==0.
