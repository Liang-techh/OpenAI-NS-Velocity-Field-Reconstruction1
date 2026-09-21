import copy, hashlib, json
import numpy as np
import pytest
import openai_ns_reconstruction.kokuno_a4_current_cartesian_exterior_to_xr_divergence_independent_audit as a4

class LinearSolenoidal:
    D=1.0; X_h=1.0; X_R=10.0
    def velocity(self,x,y,z,t):
        x,y,z,t=np.broadcast_arrays(*[np.asarray(v,float) for v in (x,y,z,t)])
        return np.stack((1.2*x+.3*y+.1*t,-.4*x-.7*y+.2*z,-.5*z+.2*x),axis=-1)

class FakeField(LinearSolenoidal):
    X_h=2.0;X_R=2.0*np.exp(5.0)

def _seal(r):
    q=copy.deepcopy(r);q.pop('receipt_sha256',None)
    r['receipt_sha256']=hashlib.sha256(a4._canon(q).encode()).hexdigest();return r

def synthetic_receipt():
    base={'step':0.005,'sampled_max':2e-7,'pooled_weighted_rms':1e-7,'pooled_volume_l2_estimate':1e-7,'normalized_sampled_max':2e-7,'normalized_weighted_rms':1e-7,'speed_rms':1.0,'seam':{'sampled_max':2e-7},'axis':{'sampled_max':2e-7},'per_zone':[],'worst':{}}
    r={'schema':a4.SCHEMA,'upstream':{'head':a4.UPSTREAM_HEAD,'source_blob':a4.UPSTREAM_SOURCE_BLOB},'protocol':{'seed':a4.SEED,'steps':list(a4.STEPS),'times':list(a4.TIMES),'eta_interval':list(a4.ETA),'zones':[list(x) for x in a4.ZONES],'points_per_zone_time':a4.N_PER_ZONE_TIME,'seam_offsets':list(a4.SEAM_OFFSETS),'divergence_gate':a4.DIV_GATE,'stability_factor':a4.STABILITY_FACTOR,'stability_floor':a4.STABILITY_FLOOR,'speed_floor':a4.SPEED_FLOOR,'mutation_epsilon':a4.MUT_EPS,'mutation_detection_floor':a4.MUT_DETECT},'resolutions':[copy.deepcopy(base),copy.deepcopy(base),copy.deepcopy(base)],'checks':{'save_reload_semantic_identity':True,'config_parameter_mutation_changes_semantic_identity':True,'velocity_mutation_detected':True,'offgrid_order_invariant':True,'post_XR_fail_closed':True},'truth_boundary':copy.deepcopy(a4.TRUTH)}
    r['audit_pass']=a4._pass(r);return _seal(r)

def test_fd2_manufactured_solenoidal():
    f=LinearSolenoidal();p=np.array([[.2,-.3,.4],[.5,.1,-.2],[-.4,.7,.3]]);t=np.array([.31,.47,.71])
    for h in a4.STEPS:
        J=a4.fd2_jacobian(f,p,t,h);assert np.max(np.abs(np.trace(J,axis1=1,axis2=2)))<2e-13

def test_sampler_is_deterministic_fresh_exterior_and_weighted():
    f=FakeField();p1=a4.make_probes(f);p2=a4.make_probes(f)
    np.testing.assert_array_equal(p1['p'],p2['p']);np.testing.assert_array_equal(p1['w'],p2['w'])
    assert len(p1['p'])==len(a4.TIMES)*len(a4.ZONES)*a4.N_PER_ZONE_TIME
    assert np.all(p1['w']>0);assert len(p1['seam'])==len(a4.TIMES)*len(a4.SEAM_OFFSETS)

def test_manufactured_metrics_and_mutation_detection():
    f=FakeField();P=a4.make_probes(f);m=a4.metrics(f,P,a4.STEPS[-1]);mm=a4.metrics(a4.Mut(f),P,a4.STEPS[-1])
    assert m['sampled_max']<2e-12;assert mm['sampled_max']>=9.9e-4

def test_frozen_gate_and_laundering_rejection():
    r=synthetic_receipt();a4.enforce_receipt(r)
    bad=copy.deepcopy(r);bad['resolutions'][-1]['sampled_max']=2e-5;bad['audit_pass']=True;_seal(bad)
    with pytest.raises(ValueError,match='laundering'):a4.enforce_receipt(bad)

def test_protocol_truth_and_checksum_mutations_fail_closed():
    r=synthetic_receipt();bad=copy.deepcopy(r);bad['protocol']['divergence_gate']=2e-5;_seal(bad)
    with pytest.raises(ValueError,match='protocol drift'):a4.enforce_receipt(bad)
    bad=copy.deepcopy(r);bad['truth_boundary']['pde_validated']=True;_seal(bad)
    with pytest.raises(ValueError,match='truth drift'):a4.enforce_receipt(bad)
    bad=copy.deepcopy(r);bad['resolutions'][-1]['sampled_max']*=1.01
    with pytest.raises(ValueError,match='checksum'):a4.enforce_receipt(bad)
