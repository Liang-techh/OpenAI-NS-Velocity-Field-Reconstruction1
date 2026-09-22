import copy, importlib.util, json, tempfile
from pathlib import Path
import numpy as np, pytest
from openai_ns_reconstruction import kokuno_current_main_pulse_leading_oscillatory_identity as m

def _has_a1():
 try: return importlib.util.find_spec(m.A1_MODULE) is not None
 except ModuleNotFoundError: return False

def test_public_contract_and_truth_boundary():
 c=m.public_contract(); assert c["public_inputs"]==["x","y","z","t"]; assert c["forbidden_velocity_inputs_present"]==[]
 for k in ("full_source_xi_11_current_cartesian_materialized","source_exact_main_pulse_amplitude_materialized","velocity_export_ready","paper_exact","pde_validated"): assert c[k] is False
 t=m._truth(); assert t["current_main_pulse_leading_plus_frozen_complete_curl_oscillation_materialized"] is True; assert t["agent3_correction_velocity_composed"] is False

def test_enforcer_rejects_truth_and_domain_promotion():
 r={"schema":m.RECEIPT_SCHEMA,"task":m.TASK,"source_provenance":m.SOURCE,"truth_boundary":m._truth(),"domain":{"xi_materializable_max":9.0,"full_source_xi_11_current_cartesian_materialized":False},"checks":{"inner_additive_max_abs":0.0,"pulse_additive_max_abs":0.0,"inner_oscillatory_vector_rms":1e-4,"axis_oscillatory_max_abs":0.0,"outside_support_oscillatory_max_abs":0.0,"save_load_velocity_replay_max_abs":0.0,"post_materializable_X_boundary_fail_closed":True}}
 m.enforce_receipt(r)
 q=copy.deepcopy(r); q["truth_boundary"]["pde_validated"]=True
 with pytest.raises(AssertionError): m.enforce_receipt(q)
 q=copy.deepcopy(r); q["domain"]["xi_materializable_max"]=11.0
 with pytest.raises(AssertionError): m.enforce_receipt(q)

@pytest.mark.skipif(not _has_a1(), reason="exact A1 #1100 sibling supplied by dedicated CI")
def test_exact_main_pulse_composition_receipt():
 f=m.default_field(); assert 0<f.xi_materializable_max<11
 x,y,z,t,_,_=m._pulse_points(f); e=m._comp._evaluate_with_backend(f.leading_backend,x,y,z,t)
 np.testing.assert_allclose(f.velocity(x,y,z,t),e.leading+e.oscillatory,rtol=0,atol=m.COMP_ATOL)
 r=m.materialize_receipt(); m.enforce_receipt(r); assert r["checks"]["post_materializable_X_boundary_fail_closed"] is True

@pytest.mark.skipif(not _has_a1(), reason="exact A1 #1100 sibling supplied by dedicated CI")
def test_save_load_rejects_rehashed_runtime_and_source_mutations():
 f=m.default_field()
 with tempfile.TemporaryDirectory() as d:
  p=Path(d)/"c.json"; f.save_candidate(p); raw=json.loads(p.read_text())
  q=copy.deepcopy(raw); o=q["configuration"]["oscillatory_runtime"]; o["payload"]["mutation_probe"]="x"; o["payload_sha256"]=m._id._sha256(o["payload"]); q["semantic_sha256"]=m._id._sha256(q["configuration"]); p.write_text(json.dumps(q))
  with pytest.raises(ValueError,match="oscillatory runtime"): m.CurrentMainPulseLeadingOscillatoryField.load_candidate(p)
  q=copy.deepcopy(raw); q["configuration"]["source_provenance"]["corrected_release_date"]="2099-01-01"; q["semantic_sha256"]=m._id._sha256(q["configuration"]); p.write_text(json.dumps(q))
  with pytest.raises(ValueError,match="configuration/provenance/truth"): m.CurrentMainPulseLeadingOscillatoryField.load_candidate(p)
