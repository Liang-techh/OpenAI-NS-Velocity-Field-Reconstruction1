"""A2 identity-bound main-pulse-prefix leading + frozen complete-curl oscillation.

Consumes exact A1 #1100 without changing its leading field.  The A1 field is only
its finite-float64-X main-pulse prefix and uses a repository-autonomous principal
pulse amplitude.  The A2 oscillation is the already-frozen bounded complete-curl
runtime.  Nothing here claims full xi=11, terminal/global velocity, paper exactness,
or Navier--Stokes validation.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from functools import lru_cache
import importlib, inspect, json, math, tempfile
from pathlib import Path
from typing import Any, Mapping
import numpy as np
from . import kokuno_current_partial_composite_identity_save_load as _id
from . import kokuno_current_partial_leading_oscillatory_velocity as _comp
from . import kokuno_source_shell_edge_localized_harmonic as _parent

TASK="K2-OSC-100"
SCHEMA="kokuno-a2-current-main-pulse-leading-oscillatory-v1"
SAVED_SCHEMA=SCHEMA+"-saved"
RECEIPT_SCHEMA=SCHEMA+"-receipt"
PARENT_PR=1099; PARENT_HEAD="3f0988ff016f41126ec546681cba696c5fa8001d"
PARENT_BLOB="a16fba5d2c77fe952d39ae0dae9be2c04b2031fd"
A1_PR=1100; A1_HEAD="2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156"
A1_MODULE="openai_ns_reconstruction.kokuno_pa16_current_cartesian_main_pulse"
A1_CLASS="KokunoPA16CurrentCartesianMainPulse"
A1_BLOB="efcce36612ffc1cb3a74f5eb9c9ae6c86748f60d"
SOURCE={"repository":"KokunoYumeto/yang-mills-interacting-workbench","commit":"143f6773feb424ad9ed3a8d116653200f20346b7","path":"navier-stokes/navier_stokes_workbench.tex","blob_sha1":"205a99807302e21a51c5eaf223390c0dfc42bcd0","corrected_release_date":"2026-09-09","classification":"structural_math_provenance_only","paper_exact_claim":False}
COMP_ATOL=1e-13; REPLAY_ATOL=2e-12; SIGNAL_FLOOR=1e-12
INNER=np.array([[.32,.11,-.20],[.41,-.17,-.08],[-.36,.24,.05],[-.52,-.16,.16],[.58,.21,.22],[.47,-.31,-.14]])
IT=np.array([.31,.39,.47,.55,.63,.71])
AX=np.array([[0,0,-.10],[0,0,0],[0,0,.12]]); AT=np.array([.37,.51,.67])
OUT=np.array([[3.5,0,0],[0,-3.5,.4],[2.8,2.8,-.5]]); OT=np.array([.41,.53,.69])
PF=np.array([.15,.50,.85]); PE=np.array([-.29,.07,.35]); PT=np.array([.43,.55,.67]); PTH=np.array([.27,.93,1.51])

A1_TRUTH={
 "public_reconstruction_source":True,"exact_current_pulse_entry_parent_consumed":True,
 "source_main_pulse_kernel_formula_consumed":True,"repository_autonomous_principal_amplitude_consumed":True,
 "current_cartesian_main_pulse_prefix_materialized":True,"current_pulse_M_and_M_eta_history_integrated":True,
 "current_v0_rebuilt_from_full_primitive_history":True,"velocity_interface_vectorized":True,
 "configuration_serializable":True,"source_exact_amplitude_root_materialized":False,
 "source_hidden_parameters_recovered":False,"source_pulse_end_MJ_corrections_materialized":False,
 "full_source_xi_11_current_cartesian_materialized":False,"source_terminal_tail_schedule_bound_into_current_velocity":False,
 "source_exterior_heat_replacement_materialized":False,"outer_global_leading_velocity_materialized":False,
 "unified_global_cartesian_velocity_export_ready":False,"matched_global_pressure_materialized":False,
 "restricted_forcing_materialized":False,"heldout_ns_residual_assessed":False,
 "same_protocol_comparable_to_st006":False,"pde_validated":False,"paper_exact":False,
 "openai_field_identified":False,"blowup_proved":False}

def _truth():
 return {"current_main_pulse_finite_X_leading_consumed":True,
 "current_main_pulse_leading_plus_frozen_complete_curl_oscillation_materialized":True,
 "full_concrete_oscillatory_runtime_digest_bound":True,"identity_preserving_save_load_available":True,
 "field_summands_changed_by_this_increment":False,"source_exact_main_pulse_amplitude_materialized":False,
 "full_source_xi_11_current_cartesian_materialized":False,"source_pulse_end_MJ_corrections_materialized":False,
 "terminal_global_leading_velocity_materialized":False,"agent3_correction_velocity_composed":False,
 "matched_global_pressure_materialized":False,"restricted_forcing_materialized":False,
 "complete_ns_residual_assessed":False,"same_protocol_comparable_to_st006":False,
 "residual_reduction_claimed":False,"velocity_export_ready":False,"paper_exact":False,"pde_validated":False}

def _parent_blob():
 p=getattr(_parent,"__file__",None)
 if not p or _comp._git_blob_sha1(p)!=PARENT_BLOB: raise RuntimeError("A2 #1099 parent source drifted")
 return PARENT_BLOB

def _a1_identity(b):
 c=type(b)
 if c.__module__!=A1_MODULE or c.__name__!=A1_CLASS: raise RuntimeError("not exact A1 #1100 backend")
 p=inspect.getsourcefile(c)
 if not p or _comp._git_blob_sha1(p)!=A1_BLOB: raise RuntimeError("A1 #1100 source drifted")
 tr=getattr(b,"truth_boundary",None)
 if not isinstance(tr,Mapping): raise RuntimeError("A1 #1100 truth missing")
 for k,v in A1_TRUTH.items():
  if tr.get(k) is not v: raise RuntimeError(f"A1 #1100 truth drifted at {k}")
 cfg=b.configuration(); sem=str(b.semantic_sha256)
 if not isinstance(cfg,Mapping) or len(sem)!=64: raise RuntimeError("A1 #1100 identity malformed")
 D=float(b.D); xp=float(b.X_p); xe=float(b.X_materializable_end); xm=float(b.xi_materializable_max)
 if not(D>0 and xp>0 and xe>xp and 0<xm<11 and all(np.isfinite([D,xp,xe,xm]))): raise RuntimeError("A1 #1100 prefix invalid")
 return {"semantic_sha256":sem,"configuration_sha256":_id._sha256(cfg),"D":D,"X_p":xp,"X_end":xe,"xi_max":xm}

def _load(cfg=None):
 try: m=importlib.import_module(A1_MODULE)
 except ModuleNotFoundError as e: raise RuntimeError("exact A1 #1100 runtime unavailable") from e
 C=getattr(m,A1_CLASS,None)
 if C is None: raise RuntimeError("exact A1 #1100 class unavailable")
 b=C() if cfg is None else C.from_configuration(cfg)
 _a1_identity(b); return b

def _pulse_points(f):
 xi=PF*f.xi_materializable_max
 X=np.array([f.leading_backend.X_from_xi(float(v)) for v in xi])
 q=(1-PT)/(1-PE*PE); z=q**f.D*PE; r=np.sqrt(2*q*X)
 return r*np.cos(PTH),r*np.sin(PTH),z,PT,X,xi

@dataclass(frozen=True)
class CurrentMainPulseLeadingOscillatoryField:
 leading_backend: Any
 _li:dict=field(init=False,repr=False,compare=False); _op:dict=field(init=False,repr=False,compare=False); _os:str=field(init=False,repr=False,compare=False)
 def __post_init__(self):
  _parent_blob(); li=_a1_identity(self.leading_backend); op=_id._oscillatory_payload()
  object.__setattr__(self,"_li",li); object.__setattr__(self,"_op",op); object.__setattr__(self,"_os",_id._sha256(op))
 @property
 def D(self): return self._li["D"]
 @property
 def X_p(self): return self._li["X_p"]
 @property
 def X_materializable_end(self): return self._li["X_end"]
 @property
 def xi_materializable_max(self): return self._li["xi_max"]
 @property
 def oscillatory_runtime_sha256(self): return self._os
 @property
 def truth_boundary(self): return _truth()
 def velocity(self,x,y,z,t): return _comp._evaluate_with_backend(self.leading_backend,x,y,z,t).velocity
 def configuration(self):
  cfg=self.leading_backend.configuration()
  return {"schema":SCHEMA,"task":TASK,"parent_agent2":{"pr":PARENT_PR,"head":PARENT_HEAD,"source_blob_sha1":PARENT_BLOB},
   "agent1":{"pr":A1_PR,"head":A1_HEAD,"module":A1_MODULE,"class":A1_CLASS,"source_blob_sha1":A1_BLOB,"configuration":json.loads(_id._canonical_json(cfg)),"configuration_sha256":self._li["configuration_sha256"],"semantic_sha256":self._li["semantic_sha256"]},
   "oscillatory_runtime":{"payload":self._op,"payload_sha256":self._os,"public_z_source_blob_sha1":_id.PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1},
   "source_provenance":SOURCE,"composition":"u_main_pulse_prefix = u_lead_A1_1100 + u_osc_frozen_complete_curl","truth_boundary":_truth()}
 @property
 def semantic_sha256(self): return _id._sha256(self.configuration())
 def save_candidate(self,path):
  p={"schema":SAVED_SCHEMA,"configuration":self.configuration(),"semantic_sha256":self.semantic_sha256}; Path(path).write_text(json.dumps(p,indent=2,sort_keys=True)+"\n"); return p
 @classmethod
 def from_configuration(cls,p):
  if not isinstance(p,Mapping) or p.get("schema")!=SCHEMA or p.get("source_provenance")!=SOURCE or p.get("truth_boundary")!=_truth(): raise ValueError("current main-pulse configuration/provenance/truth drifted")
  pa=p.get("parent_agent2",{}); a=p.get("agent1",{}); o=p.get("oscillatory_runtime",{})
  if pa!={"pr":PARENT_PR,"head":PARENT_HEAD,"source_blob_sha1":PARENT_BLOB}: raise ValueError("A2 parent identity drifted")
  if any([a.get("pr")!=A1_PR,a.get("head")!=A1_HEAD,a.get("module")!=A1_MODULE,a.get("class")!=A1_CLASS,a.get("source_blob_sha1")!=A1_BLOB]): raise ValueError("A1 #1100 provenance drifted")
  cfg=a.get("configuration"); cs=a.get("configuration_sha256"); ss=a.get("semantic_sha256")
  if not isinstance(cfg,Mapping) or _id._sha256(cfg)!=cs: raise ValueError("A1 configuration digest mismatch")
  b=_load(cfg); li=_a1_identity(b)
  if li["configuration_sha256"]!=cs or li["semantic_sha256"]!=ss: raise ValueError("A1 semantic identity drifted")
  op=o.get("payload"); os=o.get("payload_sha256"); cur=_id._oscillatory_payload()
  if not isinstance(op,Mapping) or _id._sha256(op)!=os or _id._sha256(cur)!=os or _id._canonical_json(op)!=_id._canonical_json(cur): raise ValueError("oscillatory runtime drifted")
  obj=cls(b)
  if obj.configuration()!=p: raise ValueError("reconstructed configuration drifted")
  return obj
 @classmethod
 def load_candidate(cls,path):
  raw=json.loads(Path(path).read_text()); cfg=raw.get("configuration")
  if raw.get("schema")!=SAVED_SCHEMA or not isinstance(cfg,Mapping) or _id._sha256(cfg)!=raw.get("semantic_sha256"): raise ValueError("saved current main-pulse identity invalid")
  obj=cls.from_configuration(cfg)
  if obj.semantic_sha256!=raw["semantic_sha256"]: raise ValueError("reloaded semantic identity drifted")
  return obj

@lru_cache(maxsize=1)
def default_field(): return CurrentMainPulseLeadingOscillatoryField(_load())
def velocity(x,y,z,t): return default_field().velocity(x,y,z,t)

def public_contract():
 bad={"amplitude","phase","scale","orientation","support","residual","forcing","pressure","viscosity","gain","threshold","mean","stress","correction","tolerance"}
 ps=set(inspect.signature(velocity).parameters)
 return {"public_inputs":list(inspect.signature(velocity).parameters),"forbidden_velocity_inputs_present":sorted(ps&bad),"full_source_xi_11_current_cartesian_materialized":False,"source_exact_main_pulse_amplitude_materialized":False,"velocity_export_ready":False,"paper_exact":False,"pde_validated":False}

def materialize_receipt():
 f=default_field(); inner=_comp._evaluate_with_backend(f.leading_backend,INNER[:,0],INNER[:,1],INNER[:,2],IT); axis=_comp._evaluate_with_backend(f.leading_backend,AX[:,0],AX[:,1],AX[:,2],AT)
 px,py,pz,pt,pX,pxi=_pulse_points(f); pulse=_comp._evaluate_with_backend(f.leading_backend,px,py,pz,pt); out=np.asarray(_comp.velocity_osc_batch(OUT,OT))
 with tempfile.TemporaryDirectory() as d:
  p=Path(d)/"c.json"; f.save_candidate(p); g=CurrentMainPulseLeadingOscillatoryField.load_candidate(p); replay=max(np.max(abs(g.velocity(INNER[:,0],INNER[:,1],INNER[:,2],IT)-inner.velocity)),np.max(abs(g.velocity(px,py,pz,pt)-pulse.velocity)))
 eta=.17; t=.53; th=.71; X=f.X_materializable_end*(1+1e-10); q=(1-t)/(1-eta*eta); z=q**f.D*eta; r=math.sqrt(2*q*X); closed=False
 try: f.velocity(r*math.cos(th),r*math.sin(th),z,t)
 except ValueError: closed=True
 return {"schema":RECEIPT_SCHEMA,"task":TASK,"parent_head":PARENT_HEAD,"a1_head":A1_HEAD,"source_provenance":SOURCE,"semantic_sha256":f.semantic_sha256,"domain":{"X_p":f.X_p,"X_materializable_end":f.X_materializable_end,"xi_materializable_max":f.xi_materializable_max,"full_source_xi_11_current_cartesian_materialized":False},"probes":{"pulse_X":pX.tolist(),"pulse_xi":pxi.tolist()},"checks":{"inner_additive_max_abs":float(np.max(abs(inner.velocity-inner.leading-inner.oscillatory))),"pulse_additive_max_abs":float(np.max(abs(pulse.velocity-pulse.leading-pulse.oscillatory))),"inner_oscillatory_vector_rms":float(np.sqrt(np.mean(np.sum(inner.oscillatory**2,axis=-1)))),"pulse_oscillatory_vector_rms_observation":float(np.sqrt(np.mean(np.sum(pulse.oscillatory**2,axis=-1)))),"axis_oscillatory_max_abs":float(np.max(abs(axis.oscillatory))),"outside_support_oscillatory_max_abs":float(np.max(abs(out))),"save_load_velocity_replay_max_abs":float(replay),"post_materializable_X_boundary_fail_closed":closed},"truth_boundary":_truth()}

def enforce_receipt(r):
 if r.get("schema")!=RECEIPT_SCHEMA or r.get("task")!=TASK or r.get("truth_boundary")!=_truth() or r.get("source_provenance")!=SOURCE: raise AssertionError("receipt identity/truth drifted")
 c=r.get("checks",{}); d=r.get("domain",{}); xm=float(d.get("xi_materializable_max",np.nan))
 if float(c.get("inner_additive_max_abs",np.inf))>COMP_ATOL or float(c.get("pulse_additive_max_abs",np.inf))>COMP_ATOL: raise AssertionError("additive closure failed")
 if float(c.get("inner_oscillatory_vector_rms",0))<SIGNAL_FLOOR: raise AssertionError("inner oscillation vacuous")
 if float(c.get("axis_oscillatory_max_abs",np.inf))!=0 or float(c.get("outside_support_oscillatory_max_abs",np.inf))!=0: raise AssertionError("oscillatory support/axis contract failed")
 if float(c.get("save_load_velocity_replay_max_abs",np.inf))>REPLAY_ATOL or c.get("post_materializable_X_boundary_fail_closed") is not True: raise AssertionError("replay/boundary gate failed")
 if not(np.isfinite(xm) and 0<xm<11) or d.get("full_source_xi_11_current_cartesian_materialized") is not False: raise AssertionError("finite-X prefix boundary drifted")
