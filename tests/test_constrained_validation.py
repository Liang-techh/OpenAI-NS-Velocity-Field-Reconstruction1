import numpy as np
import pytest
from openai_ns_reconstruction.constrained_validation import residual

@pytest.mark.parametrize('t',[0.25,0.5,0.75])
def test_polynomial_manufactured_solution_and_force_sign_detection(t):
    # Divergence-free u=(t*y^2,0,0); p=x*z. Independent exact force.
    def u(x,t): return np.column_stack((t*x[:,1]**2,np.zeros((len(x),2))))
    def p(x,t): return x[:,0]*x[:,2]
    def f(x,t): return np.column_stack((x[:,1]**2+x[:,2]-0.02*t,np.zeros(len(x)),x[:,0]))
    x=np.random.default_rng(9).uniform(-1,1,(50,3))
    r=residual(u,p,f,x,t)
    assert np.max(np.abs(r['momentum']))<1e-10
    assert np.max(np.abs(r['divergence']))<1e-10
    wrong=residual(u,p,lambda x,t:-f(x,t),x,t)
    assert np.max(np.abs(wrong['momentum']))>1

def test_structural_gate_rejects_energy_collapse():
    import json
    from pathlib import Path
    from openai_ns_reconstruction.constrained_candidate import CompactCandidate
    from openai_ns_reconstruction.constrained_validation import structure_metrics
    cfg=json.loads(Path('configs/constraints.json').read_text())
    c=CompactCandidate().normalized()
    assert structure_metrics(c,cfg)['sampled_constraints_pass']
    class CollapsedEnergy:
        velocity=c.velocity
        def energy(self,time,order):return 0.09
    assert not structure_metrics(CollapsedEnergy(),cfg)['sampled_constraints_pass']
