import numpy as np
from openai_ns_reconstruction.constrained_optimize import training_residual


def test_training_operator_with_nonzero_convection():
    class Polynomial:
        def velocity(self,x,t):
            return np.asarray(t)[...,None]*x[:,[1,2,0]]**2
        def pressure(self,x,t):
            return x[:,0]*x[:,2]
    def force(x,t):
        a,b,c=x.T
        ut=x[:,[1,2,0]]**2
        conv=2*t[:,None]**2*np.column_stack((b*c*c,c*a*a,a*b*b))
        gp=np.column_stack((c,np.zeros(len(x)),a))
        return ut+conv+gp-0.02*t[:,None]
    x=np.random.default_rng(145).uniform(-1,1,(50,3))
    t=np.linspace(.3,.7,len(x))
    r=training_residual(Polynomial(),force,x,t,.01)
    assert np.max(np.abs(r))<1e-8
