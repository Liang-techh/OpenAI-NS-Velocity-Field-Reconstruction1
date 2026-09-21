"""Independent A4 audit of A1 #947 current Xi prefix moments.

Reference numerics use only save/reloaded public ``profile_values(X,eta)`` and
piecewise composite Simpson quadrature.  A1's Gauss--Legendre moment path is
only the quantity under audit.  This is source-coordinate evidence, not a
Cartesian/full-NS validation; the final 1e-3 momentum and 1e-5 divergence gates
remain unchanged and unevaluated here.
"""
from __future__ import annotations

import argparse, copy, hashlib, json, math, tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_pa10_actual_xi_prefix_moments import KokunoPA10ActualXiPrefixMoments

TASK = "K4-VAL-091"
SCHEMA = "kokuno-a4-actual-xi-prefix-moments-independent-audit-v1"
A1_PR = 947
A1_HEAD = "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b"
A1_SOURCE_BLOB = "92dce65c9a8793e06497a7e7347c832861032238"
MOMENT_NAMES = ("M", "I", "J", "S", "C_p")
FROZEN_FINAL_GATES = {
    "normalized_momentum_max": 1e-3,
    "normalized_momentum_volume_l2": 1e-3,
    "divergence_max": 1e-5,
    "divergence_volume_l2": 1e-5,
}

@dataclass(frozen=True)
class FrozenXiMomentAuditProtocol:
    seed: int = 9173631
    random_eta_count: int = 10
    simpson_panels_per_segment: tuple[int, int, int] = (48, 96, 192)
    relative_rms_gate: float = 5e-3
    relative_sampled_max_gate: float = 2e-2
    refinement_ratio_gate: float = 6.0
    refinement_floor: float = 5e-9
    nontrivial_discrepancy_l2_floor: float = 1e-14
PROTOCOL = FrozenXiMomentAuditProtocol()

TRUTH_BOUNDARY = {
    "current_candidate_side_Xi_prefix_moments_materialized": True,
    "current_candidate_side_PA15_discrepancy_materialized": True,
    "current_Xi_prefix_moments_independently_audited": False,
    "current_Xi_PA15_discrepancy_independently_audited": False,
    "source_prepared_upstream_five_moment_discrepancy_materialized": False,
    "source_T_sh_lower_bound_verified": False,
    "pa16_repair_applied": False,
    "inner_to_outer_join_completed": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "correction_velocity_materialized": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "heldout_complete_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
}

def _rms(x: Any) -> float:
    a = np.asarray(x, float)
    return float(np.sqrt(np.mean(a*a)))

def _finite(x: Any, name: str) -> np.ndarray:
    a = np.asarray(x, float)
    if np.any(~np.isfinite(a)):
        raise ValueError(f"{name} must be finite")
    return a

def _canonical(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)

def public_contract() -> dict[str, Any]:
    return {
        "task": TASK, "schema": SCHEMA, "agent1_pr": A1_PR,
        "agent1_exact_head": A1_HEAD, "agent1_source_blob": A1_SOURCE_BLOB,
        "reference": "save/reloaded public profile_values + piecewise composite Simpson",
        "production_under_audit": ["physical_prefix_moments_at_Xi", "actual_PA15_scaled_moments", "incoming_PA15_discrepancy"],
        "protocol": asdict(PROTOCOL), "moment_names": list(MOMENT_NAMES),
        "frozen_final_gates": dict(FROZEN_FINAL_GATES), "truth_boundary": dict(TRUTH_BOUNDARY),
    }

def _require(field: Any) -> None:
    for n in ("profile_values", "physical_prefix_moments_at_Xi", "actual_PA15_scaled_moments",
              "incoming_PA15_discrepancy", "save_configuration", "semantic_sha256", "eta_interval",
              "X_0", "X_1", "X_b_ref", "X_i", "log_X_R", "log_x_i", "log_P_star", "C"):
        if not hasattr(field, n): raise TypeError(f"missing public surface {n}")
    if float(field.X_b_ref) != 100.0 or float(field.X_i) != 110.0:
        raise ValueError("audit frozen to X=100 -> Xi=110")

def _densities(field: Any, X: np.ndarray, eta: float) -> np.ndarray:
    v = field.profile_values(X, np.full_like(X, eta))
    F, U = _finite(v["F_actual_prefix"], "F"), _finite(v["U_actual_prefix"], "U")
    H = 2*X*F
    return np.stack((U, H, U*H, U*U-X*F*F, F*F), axis=1)

def _simpson_segment(field: Any, eta: float, a: float, b: float, panels: int) -> np.ndarray:
    if panels < 2 or panels % 2: raise ValueError("panels must be positive even")
    X = np.linspace(a, b, panels+1)
    w = np.ones(panels+1); w[1:-1:2] = 4; w[2:-1:2] = 2
    return ((b-a)/panels/3.0) * (w @ _densities(field, X, eta))

def independent_physical_prefix_moments(field: Any, eta: float, panels: int) -> np.ndarray:
    _require(field)
    lo, hi = map(float, field.eta_interval)
    if not lo <= eta <= hi: raise ValueError("eta outside executable interval")
    segs = ((0.0,float(field.X_0)), (float(field.X_0),float(field.X_1)),
            (float(field.X_1),100.0), (100.0,110.0))
    return sum((_simpson_segment(field, eta, a, b, panels) for a,b in segs), np.zeros(5))

def _scaled(field: Any, p: np.ndarray) -> np.ndarray:
    s1, s15 = math.exp(-float(field.log_X_R)), math.exp(-1.5*float(field.log_X_R))
    return np.array([p[0]*s1,p[1]*s15,p[2]*s15,p[3]*s1,p[4]], float)

def independent_ideal_PA15_scaled_moments(field: Any, eta: float) -> np.ndarray:
    f, lx, lp = 1/(1+eta*eta), float(field.log_x_i), float(field.log_P_star)
    x = math.exp(lx)
    I = (5*math.sqrt(2)/8)*f*math.exp(lp+1.6*lx)
    Eterm = (5/12)*f*f*math.exp(2*lp+1.2*lx)
    Cp = 2.5*f*f/(float(field.C)**2)
    return np.array([4*eta*x,I,4*eta*I,16*eta*eta*x-Eterm,Cp], float)

def _row_metrics(prod: np.ndarray, ref: np.ndarray) -> dict[str, dict[str,float]]:
    out = {}
    for j,n in enumerate(MOMENT_NAMES):
        e = ref[:,j]-prod[:,j]
        out[n] = {
            "absolute_rms": _rms(e), "absolute_max": float(np.max(np.abs(e))),
            "relative_rms": _rms(e)/max(_rms(ref[:,j]),_rms(prod[:,j]),1e-30),
            "relative_sampled_max": float(np.max(np.abs(e)))/max(float(np.max(np.abs(ref[:,j]))),float(np.max(np.abs(prod[:,j]))),1e-30),
        }
    return out

def _metric_gate(m: Mapping[str,Mapping[str,float]]) -> bool:
    return all(r["relative_rms"] <= PROTOCOL.relative_rms_gate and r["relative_sampled_max"] <= PROTOCOL.relative_sampled_max_gate for r in m.values())

def _stability(levels: list[np.ndarray]) -> dict[str,dict[str,float|bool]]:
    out = {}
    for j,n in enumerate(MOMENT_NAMES):
        c,m,f = levels[0][:,j], levels[1][:,j], levels[2][:,j]
        scale = max(_rms(m),_rms(f),1e-30); cm, mf = _rms(c-m)/scale, _rms(m-f)/scale
        ratio = 1e300 if mf == 0 else cm/mf
        out[n] = {"coarse_to_medium_scale_normalized_rms":cm,"medium_to_fine_scale_normalized_rms":mf,
                  "ratio":ratio,"numerical_floor_hit":mf<=PROTOCOL.refinement_floor,
                  "passed":bool(mf<=PROTOCOL.refinement_floor or ratio>=PROTOCOL.refinement_ratio_gate)}
    return out

def _sample_eta(field: Any) -> np.ndarray:
    lo,hi=map(float,field.eta_interval); margin=max(.08*(hi-lo),1e-3)
    r=np.random.default_rng(PROTOCOL.seed).uniform(lo+margin,hi-margin,PROTOCOL.random_eta_count)
    c=np.array([0.,-1e-8,1e-8,-1e-6,1e-6])
    if np.any((c<lo)|(c>hi)): raise ValueError("eta interval excludes frozen center probes")
    return np.concatenate((r,c))

def _checksum_valid_C_drift_rejected(config: Mapping[str,Any]) -> bool:
    p=copy.deepcopy(dict(config)); params=p["outer_schedule"]["parameters"]; params["C"]=999.0
    p["outer_schedule"]=KokunoOuterReservedPatchSchedule(**params).to_payload()
    try: KokunoPA10ActualXiPrefixMoments.from_configuration(p)
    except (TypeError,ValueError): return True
    return False

def _worst(eta: np.ndarray, prod: np.ndarray, ref: np.ndarray) -> dict[str,dict[str,float]]:
    out={}
    for j,n in enumerate(MOMENT_NAMES):
        e=np.abs(ref[:,j]-prod[:,j]); i=int(np.argmax(e))
        out[n]={"eta":float(eta[i]),"production":float(prod[i,j]),"independent":float(ref[i,j]),"absolute_error":float(e[i])}
    return out

def materialize_actual_xi_prefix_moment_audit(field: Any) -> dict[str,Any]:
    _require(field); original_sha=str(field.semantic_sha256)
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-xi-moment-") as d:
        path=Path(d)/"candidate.json"; config=field.save_configuration(path); rebound=KokunoPA10ActualXiPrefixMoments.load_configuration(path)
    if rebound.semantic_sha256 != original_sha: raise ValueError("save/reload semantic identity changed")
    eta=_sample_eta(rebound)
    prod_p=np.stack([rebound.physical_prefix_moments_at_Xi(float(e)) for e in eta])
    prod_s=np.stack([rebound.actual_PA15_scaled_moments(float(e)) for e in eta])
    prod_d=_finite(rebound.incoming_PA15_discrepancy(eta),"production discrepancy")
    plev=[]; slev=[]; dlev=[]
    for panels in PROTOCOL.simpson_panels_per_segment:
        p=np.stack([independent_physical_prefix_moments(rebound,float(e),panels) for e in eta])
        s=np.stack([_scaled(rebound,row) for row in p])
        ideal=np.stack([independent_ideal_PA15_scaled_moments(rebound,float(e)) for e in eta])
        plev.append(p); slev.append(s); dlev.append(s-ideal)
    pm=[_row_metrics(prod_p,x) for x in plev]; sm=[_row_metrics(prod_s,x) for x in slev]; dm=[_row_metrics(prod_d,x) for x in dlev]
    ps,ss,ds=_stability(plev),_stability(slev),_stability(dlev)
    dl2=np.sqrt(np.sum(dlev[-1]**2,axis=1))
    p_mut=prod_p.copy(); p_mut[:,0]*=.99
    row=int(np.argmax([_rms(dlev[-1][:,j]) for j in range(5)])); d_mut=prod_d.copy(); d_mut[:,row]*=-1
    mutations={"production_M_times_0p99_detected":not _metric_gate(_row_metrics(p_mut,plev[-1])),
               "sign_flipped_discrepancy_row_detected":not _metric_gate(_row_metrics(d_mut,dlev[-1])),
               "sign_flipped_discrepancy_row":MOMENT_NAMES[row]}
    gates={"physical_moment_match":_metric_gate(pm[-1]),"actual_PA15_scaled_match":_metric_gate(sm[-1]),
           "incoming_discrepancy_match":_metric_gate(dm[-1]),"physical_resolution_stable":all(v["passed"] for v in ps.values()),
           "scaled_resolution_stable":all(v["passed"] for v in ss.values()),"discrepancy_resolution_stable":all(v["passed"] for v in ds.values()),
           "discrepancy_nontrivial":float(np.max(dl2))>PROTOCOL.nontrivial_discrepancy_l2_floor,
           "save_reload_semantic_identity":rebound.semantic_sha256==original_sha,
           "checksum_valid_C_drift_rejected":_checksum_valid_C_drift_rejected(config),
           "mutation_firewall":bool(mutations["production_M_times_0p99_detected"] and mutations["sign_flipped_discrepancy_row_detected"])}
    passed=bool(all(gates.values())); truth=dict(TRUTH_BOUNDARY)
    truth["current_Xi_prefix_moments_independently_audited"]=passed; truth["current_Xi_PA15_discrepancy_independently_audited"]=passed
    receipt={"task":TASK,"schema":SCHEMA,"agent1":{"pr":A1_PR,"exact_head":A1_HEAD,"source_blob":A1_SOURCE_BLOB},
             "artifact_semantic_sha256":original_sha,"protocol":asdict(PROTOCOL),"eta_samples":eta.tolist(),
             "physical_resolution_metrics":pm,"scaled_resolution_metrics":sm,"discrepancy_resolution_metrics":dm,
             "physical_resolution_stability":ps,"scaled_resolution_stability":ss,"discrepancy_resolution_stability":ds,
             "physical_worst_witness":_worst(eta,prod_p,plev[-1]),"discrepancy_worst_witness":_worst(eta,prod_d,dlev[-1]),
             "max_independent_discrepancy_l2":float(np.max(dl2)),"mutation_checks":mutations,"gates":gates,"passed":passed,
             "frozen_final_gates":dict(FROZEN_FINAL_GATES),"truth_boundary":truth}
    receipt["receipt_sha256"]=hashlib.sha256(_canonical(receipt).encode()).hexdigest(); return receipt

def enforce_audit(receipt: Mapping[str,Any]) -> None:
    if receipt.get("passed") is not True: raise RuntimeError("frozen A4 Xi-prefix-moment audit failed:\n"+json.dumps({"gates":receipt.get("gates"),"physical_resolution_metrics":receipt.get("physical_resolution_metrics"),"discrepancy_resolution_metrics":receipt.get("discrepancy_resolution_metrics"),"mutation_checks":receipt.get("mutation_checks")},indent=2,sort_keys=True))
    t=receipt["truth_boundary"]
    for k in ("outer_global_leading_velocity_materialized","matched_global_pressure_materialized","restricted_forcing_materialized","correction_velocity_materialized","heldout_complete_ns_residual_assessed","same_protocol_comparable_to_st006","pde_validated"):
        if t.get(k) is not False: raise RuntimeError(f"scoped audit illegally promoted {k}")

def _main() -> int:
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--output",type=Path,required=True); ap.add_argument("--no-enforce",action="store_true"); a=ap.parse_args()
    r=materialize_actual_xi_prefix_moment_audit(KokunoPA10ActualXiPrefixMoments()); a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"passed":r["passed"],"artifact_semantic_sha256":r["artifact_semantic_sha256"],"receipt_sha256":r["receipt_sha256"],"max_independent_discrepancy_l2":r["max_independent_discrepancy_l2"],"gates":r["gates"]},indent=2,sort_keys=True))
    if not a.no_enforce: enforce_audit(r)
    return 0

if __name__ == "__main__": raise SystemExit(_main())
