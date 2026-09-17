"""Target-free F(2,0) radial-swirl capacity audit for the capped bipolar field."""
from __future__ import annotations
from dataclasses import replace
from typing import Any, Iterable
import numpy as np
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .eq45_supported_delivery import Eq45SupportedDeliveryField

TASK_ID = "CR003-BIPOLAR-F20-RADIAL-SWIRL-CAPACITY-035"
F10, F20 = (1, 0), (2, 0)
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "force_or_pressure_fitted": False,
    "held_out_pde_residual_evaluated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _source_field() -> Eq45SupportedVelocityCandidate:
    return Eq45SupportedDeliveryField.load_candidate(
        "artifacts/bipolar_joint_capped/candidate.json"
    ).candidate


def _extend_radial2(basis: Eq45CompactProfileBasis) -> Eq45CompactProfileBasis:
    modes = tuple((i, j) for i in range(3) for j in range(basis.eta_degree + 1))
    index = {m: k for k, m in enumerate(modes)}
    phi = np.zeros(len(modes)); swirl = np.zeros(len(modes))
    for k, mode in enumerate(basis.mode_indices):
        phi[index[mode]] = basis.phi_coefficients[k]
        swirl[index[mode]] = basis.swirl_coefficients[k]
    return Eq45CompactProfileBasis(
        radial_degree=2, eta_degree=basis.eta_degree,
        phi_coefficients=tuple(phi), swirl_coefficients=tuple(swirl),
        x_cut=basis.x_cut, eta_cut=basis.eta_cut,
        cutoff_power=basis.cutoff_power, coefficient_limit=basis.coefficient_limit,
    )


def _embedded_field() -> Eq45SupportedVelocityCandidate:
    base = _source_field()
    basis = _extend_radial2(base.parent.profile_basis)
    return replace(base, parent=replace(base.parent, profile_basis=basis))


def _with_swirl_delta(field, mode, delta):
    basis = field.parent.profile_basis
    k = basis.mode_indices.index(mode)
    c = list(basis.swirl_coefficients); c[k] += float(delta)
    if abs(c[k]) > basis.coefficient_limit:
        raise ValueError("perturbation exceeds coefficient_limit")
    return replace(field, parent=replace(field.parent,
        profile_basis=replace(basis, swirl_coefficients=tuple(c))))


def _ring(r, z, n=16):
    th = 2*np.pi*np.arange(n)/n
    return np.column_stack((r*np.cos(th), r*np.sin(th), np.full(n, z)))


def _probes():
    inner = np.concatenate([_ring(r,z) for r in (.30,.50,.70) for z in (-.4,0,.4)])
    outer = np.concatenate([_ring(r,z) for r in (1.05,1.25,1.45) for z in (-.4,0,.4)])
    return inner, outer, np.concatenate((inner, outer))


def _response(field, mode, step, points, times):
    plus, minus = _with_swirl_delta(field, mode, step), _with_swirl_delta(field, mode, -step)
    return np.concatenate([((plus.at_points(points,t)-minus.at_points(points,t))/(2*step)).reshape(-1) for t in times])


def _region(field, mode, step, points, time=.5):
    plus, minus = _with_swirl_delta(field, mode, step), _with_swirl_delta(field, mode, -step)
    d = (plus.at_points(points,time)-minus.at_points(points,time))/(2*step)
    x,y = points[:,0], points[:,1]; r=np.hypot(x,y)
    ur=(x*d[:,0]+y*d[:,1])/r; ut=(-y*d[:,0]+x*d[:,1])/r; uz=d[:,2]
    rms=lambda a: float(np.sqrt(np.mean(np.asarray(a)**2)))
    return {"swirl_rms":rms(ut), "poloidal_rms":float(np.sqrt(np.mean(ur*ur+uz*uz)))}


def audit_bipolar_f20_capacity(*, times: Iterable[float]=(.3125,.5,.75), coefficient_steps: Iterable[float]=(.02,.01)) -> dict[str, Any]:
    times=tuple(map(float,times)); steps=tuple(map(float,coefficient_steps))
    if len(steps)!=2 or not steps[0]>steps[1]>0: raise ValueError("need two decreasing positive steps")
    base, field = _source_field(), _embedded_field(); inner,outer,points=_probes()
    embedding=max(float(np.max(np.abs(base.at_points(points,t)-field.at_points(points,t)))) for t in times)
    mats=[np.column_stack([_response(field,F10,h,points,times),_response(field,F20,h,points,times)]) for h in steps]
    a=mats[-1]; s=np.linalg.svd(a,compute_uv=False); tol=max(a.shape)*np.finfo(float).eps*s[0]
    rank=int(np.sum(s>tol)); cond=float(s[0]/s[-1]) if rank==2 else float("inf")
    refine=float(np.linalg.norm(mats[0]-a)/max(np.linalg.norm(a),1e-300))
    f10,f20=a[:,0],a[:,1]; proj=f10*(np.dot(f10,f20)/max(np.dot(f10,f10),1e-300))
    novelty=float(np.linalg.norm(f20-proj)/max(np.linalg.norm(f20),1e-300))
    loc={}
    for mode,name in ((F10,"F10"),(F20,"F20")):
        i,o=_region(field,mode,steps[-1],inner),_region(field,mode,steps[-1],outer)
        loc[name]={"inner":i,"outer":o,"outer_to_inner_swirl_response":o["swirl_rms"]/max(i["swirl_rms"],1e-300),
                   "swirl_fraction":(i["swirl_rms"]+o["swirl_rms"])/max(i["swirl_rms"]+o["swirl_rms"]+i["poloidal_rms"]+o["poloidal_rms"],1e-300)}
    return {"task_id":TASK_ID,"source_candidate_sha256":base.sha256,"embedding_max_abs_velocity_error":embedding,
            "response_singular_values":[float(v) for v in s],"response_rank":rank,"response_condition_number":cond,
            "finite_difference_refinement_relative_change":refine,"F20_novelty_outside_F10_span":novelty,
            "localization_at_t_050":loc,"F20_vs_F10_outer_selectivity_ratio":loc["F20"]["outer_to_inner_swirl_response"]/max(loc["F10"]["outer_to_inner_swirl_response"],1e-300),
            "parameter_growth":{"spatial_modes_added":1,"scalar_coefficients_added_if_selected":1},"truth_boundary":dict(TRUTH_BOUNDARY)}
