"""Whole-axis analytic candidate core moment enclosures at s=4.

Uses the accepted fixed-point Xh norm and prescribed axis jets, not point
interpolation. These bounds are deliberately conservative functional data.
They do not constitute a finite whole-axis core atlas or transition matching.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_pressure_function import load_datum, pressure_jets
from lei_ren_part1_paper_factored_core_positivity import squared_axis_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets

HERE = Path(__file__).resolve().parent
NORM = "lei_ren_part1_paper_candidate_core_tail_budget_Lambda120.json"
CENTER = "lei_ren_part1_paper_candidate_shared_inlet_Lambda120.json"
LAMBDA = "1e120"
J = "1e-14"
DELTA = "1e-200"
LOGC = "5e151"


def restore(ctx, record):
    return ctx.mpf([mp.make_mpf(tuple(record["lower_exact_mpf_tuple"])),
                    mp.make_mpf(tuple(record["upper_exact_mpf_tuple"]))])


def symmetric(ctx, value):
    hi = endpoints(value)[1]
    return ctx.mpf([-hi, hi])


def norm_coefficient_bound(ctx, norm, h, smax, k):
    """Ordinary axial coefficient bound with the identically zero n=0 row.

    Xh gives |f[n,k]| <= norm*binom(n+k,k)/(20^n*h^k*(n+1)^2*(k+1)^2).
    Drop only (n+1)^(-2) and sum n>=1 exactly by its generating function.
    """
    r = smax/20
    return norm*((1-r)**(-k-1)-1)/(h**k*(k+1)**2)


def global_core_moments(ctx, centers, scaled_exit="4", order=3):
    if order != 3:
        raise ValueError("This receipt currently supplies ordinary axial order3")
    z = ctx.mpf(centers)
    if endpoints(z)[0] < -1 or endpoints(z)[1] > 1:
        raise ValueError("Require real centers in [-1,1]")
    s = ctx.mpf(scaled_exit)
    if endpoints(s)[0] <= 0 or endpoints(s)[1] > mp.mpf("4.1"):
        raise ValueError("Require 0<s<=4.1")
    norm = json.loads((HERE/NORM).read_text(encoding="utf-8"))
    if mp.mpf(norm["Lambda"]) != mp.mpf(LAMBDA):
        raise ValueError("Wrong analytic candidate Lambda")
    for name, digest in norm["input_hashes"].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError("Analytic norm dependency changed: " + name)
    datum, pressure_hashes = load_datum()
    if datum["accepted_schedule_sha256"] != norm["accepted_schedule_sha256"]:
        raise ValueError("Core and functional pressure datum disagree")
    h = restore(ctx, norm["Xh_parameter"])
    phinorm = restore(ctx, norm["analytic_Phi_Xh_norm_upper"])
    psinorm = restore(ctx, norm["analytic_Psi_Xh_norm_upper"])
    phi0 = ctx.mpf([endpoints(restore(ctx, norm["Phi_real_lower"]))[0],
                    endpoints(restore(ctx, norm["Phi_real_upper"]))[1]])
    phi = IntervalTaylor(ctx, [phi0]+[
        symmetric(ctx, norm_coefficient_bound(ctx, phinorm, h, s, k)) for k in range(1, 4)])
    psi = IntervalTaylor(ctx, [symmetric(ctx, norm_coefficient_bound(ctx, psinorm, h, s, k))
                              for k in range(4)])
    eps = 1/ctx.mpf(LAMBDA)
    axis = uniform_axis_jets(ctx, radius=1, j=J, Lambda=LAMBDA, logC=LOGC,
                             delta=DELTA, length=4, axial_interval=z)
    F0 = IntervalTaylor(ctx, axis["F0"])
    S = IntervalTaylor(ctx, squared_axis_rows(ctx, axis["gradient_coefficients"],
                                             ctx.mpf(LAMBDA), axis["F0_interval"], 4))
    U0 = IntervalTaylor(ctx, axis["U0"])
    U = U0+psi*eps
    U2 = U*U
    U2 = IntervalTaylor(ctx, [U[0]**2]+list(U2.coefficients[1:]))
    phi2 = phi*phi
    # The jets enclose the field at every integration radius, so positive
    # radial weights can be integrated without quadrature or interpolation.
    normalized = dict(theta=phi*s**2, z=U*s,
                      theta_z=(phi*U)*s**2,
                      z_theta=U2*s-S*phi2*(eps*s**2/2), p=phi2*s)
    physical = dict(theta=normalized["theta"]*F0*eps**2,
                    z=normalized["z"]*eps,
                    theta_z=normalized["theta_z"]*F0*eps**2,
                    z_theta=normalized["z_theta"]*eps,
                    p=normalized["p"]*S*eps)
    p0 = IntervalTaylor(ctx, pressure_jets(ctx, z, 3, datum)["physical_pressure_coefficients"])
    return dict(center_interval=z, scaled_exit=s, physical_exit=s*eps,
                ordinary_axial_Taylor_order=3,
                normalized_moment_jets=normalized, physical_moment_jets=physical,
                pressure_axis_jets=p0, recovered_pressure_exit=p0+physical["p"],
                uniform_field_coefficient_enclosures=dict(Phi=phi, Psi=psi, Uz=U),
                F0_axis_jet_enclosure=F0,
                accepted_schedule_sha256=norm["accepted_schedule_sha256"],
                proof_hashes={NORM: hashlib.sha256((HERE/NORM).read_bytes()).hexdigest(), **pressure_hashes},
                analytic_functional_core_moments_enclosed=True,
                finite_whole_axis_core_atlas_generated=False,
                transition_moments_generated=False,
                terminal_functional_five_moment_closure=False)


def packed(packet):
    result = dict(packet)
    for key in ("normalized_moment_jets", "physical_moment_jets", "uniform_field_coefficient_enclosures"):
        result[key] = {name: list(jet.coefficients) for name, jet in packet[key].items()}
    for key in ("pressure_axis_jets", "recovered_pressure_exit", "F0_axis_jet_enclosure"):
        result[key] = list(packet[key].coefficients)
    return result


def run():
    ctx = MPIntervalContext()
    ctx.dps = 260
    with mp.workdps(300):
        packet = global_core_moments(ctx, ["-1", "1"])
        center = json.loads((HERE/CENTER).read_text(encoding="utf-8"))
        if center["Lambda"] != LAMBDA or center["accepted_schedule_sha256"] != packet["accepted_schedule_sha256"]:
            raise ValueError("Independent center data has different construction parameters")
        count = 0
        for name, rows in center["analytic_normalized_moment_enclosures"].items():
            for k, row in enumerate(rows):
                lo, hi = endpoints(restore(ctx, row))
                gl, gh = endpoints(packet["normalized_moment_jets"][name][k])
                if not gl <= lo <= hi <= gh:
                    raise AssertionError(("Global analytic moment envelope missed center data", name, k))
                count += 1
        report = dict(**packed(packet),
                      input_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                    for name in (NORM, CENTER, "lei_ren_part1_paper_uniform_axis_jets.py",
                                                 "lei_ren_part1_paper_candidate_pressure_function.py",
                                                 "lei_ren_part1_paper_factored_core_positivity.py",
                                                 "lei_ren_part1_paper_interval_taylor.py")},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      independent_center_moment_containment_count=count,
                      positive_radial_weight_integrals_bounded_without_quadrature=True,
                      zero_axis_rows_Phi_minus1_and_Psi_used=True,
                      field_bound_formula="norm*((1-s/20)^(-k-1)-1)/(h^k*(k+1)^2)",
                      amplitude_fitted_or_reset=False,
                      pressure_fit_or_anchor_change=False,
                      original_parameter_errors_enclosed=False)
        Path(__file__).with_suffix(".json").write_text(json.dumps(encode(report), indent=2)+"\n", encoding="utf-8")
        print("Whole-axis analytic core moment envelopes generated; center containment checks", count,
              "actual transition and terminal matching remain open", flush=True)
        return report


if __name__ == "__main__":
    run()
