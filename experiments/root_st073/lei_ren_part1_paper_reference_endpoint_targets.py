"""Directed functional reference moment targets at the physical endpoint Rh.

Implements Lei--Ren (9.3) using the accepted outer schedule. Radial factors
remain separate in logarithmic form; s=Lambda*R is never used as Rh.
The targets are analytic reference data, not achieved transition moments.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).resolve().parent
SOURCE = "lei_ren_part1_paper_coherent_pressure_source_alignment.json"
ACCEPTED_SHA = "736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4"
NAMES = ("theta", "z", "theta_z", "z_theta", "p")


def accepted_parameters():
    record = json.loads((HERE/SOURCE).read_text(encoding="utf-8"))
    accepted = record["accepted_schedule"]
    digest = hashlib.sha256(json.dumps(accepted["inputs"], sort_keys=True).encode()).hexdigest()
    if digest != ACCEPTED_SHA or accepted["sha256"] != digest:
        raise ValueError("Accepted endpoint parameters changed")
    return accepted["inputs"]


def reference_targets(ctx, Z, order=3, parameters=None):
    """Return ordinary axial jets and physical radial scale logarithms.

    For an interval Z, coefficients enclose derivatives at every real center
    in that interval. They do not represent a single Taylor expansion about
    an interval. Multiply each jet by exp(moment_log_scales[name]) to obtain
    physical moment jets, or retain the factorization for comparisons.
    """
    if not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    if parameters is None:
        parameters = accepted_parameters()
    elif hashlib.sha256(json.dumps(parameters, sort_keys=True).encode()).hexdigest() != ACCEPTED_SHA:
        raise ValueError("Reference targets require the authoritative accepted schedule")
    center = ctx.mpf(Z)
    lo, hi = endpoints(center)
    if lo < -1 or hi > 1:
        raise ValueError("Require real axial centers in [-1,1]")
    logRh = ctx.mpf(parameters["logRref"])-5
    logPstar = ctx.mpf(parameters["logPstar"])
    delta = ctx.mpf(parameters["delta"])
    # q=1+Z^2. Keep the exact real square dependency in coefficient0,
    # especially when the center interval crosses zero.
    qcoeff = [1+center**2]
    if order >= 1:
        qcoeff.append(2*center)
    if order >= 2:
        qcoeff.append(ctx.mpf(1))
    qcoeff.extend(ctx.mpf(0) for _ in range(max(0, order+1-len(qcoeff))))
    q = IntervalTaylor(ctx, qcoeff)
    z = IntervalTaylor.variable(ctx, center, order)
    z2 = z*z
    z2 = IntervalTaylor(ctx, [center**2]+list(z2.coefficients[1:]))
    uh = q.reciprocal()*ctx.exp(logPstar-ctx.mpf(1)/2)
    theta = uh*(ctx.mpf(5)/8)
    jets = dict(theta=theta, z=z*4, theta_z=z*theta*4,
                z_theta=z2*16-uh*uh*(ctx.mpf(5)/12),
                p=uh*uh*(ctx.mpf(5)/2))
    angular_scale = logRh*(ctx.mpf(3)/2)+ctx.log(ctx.mpf(2))/2
    return dict(center_Z=center, ordinary_axial_Taylor_order=order,
                delta=delta, logPstar=logPstar, logRh=logRh,
                reference_Utheta_at_Rh=uh,
                factored_moment_jets=jets,
                moment_log_scales=dict(theta=angular_scale, z=logRh,
                                       theta_z=angular_scale, z_theta=logRh, p=ctx.mpf(0)),
                moment_scale_definitions=dict(theta="Rh*sqrt(2Rh)", z="Rh",
                                              theta_z="Rh*sqrt(2Rh)", z_theta="Rh", p="1"),
                radial_coordinate="physical R", relative_log_endpoint_y="-5",
                core_scaled_exit_substituted_for_Rh=False,
                actual_transition_moments_generated=False,
                reference_endpoint_matching_achieved=False,
                terminal_functional_five_moment_closure=False)


def physical_moment_jets(ctx, target):
    return {name: jet*ctx.exp(target["moment_log_scales"][name])
            for name, jet in target["factored_moment_jets"].items()}


def run():
    params = accepted_parameters()
    ctx = MPIntervalContext()
    ctx.dps = 260
    with mp.workdps(340):
        samples = {}
        count = 0
        for text in ("-1", "0", ".3", "1"):
            target = reference_targets(ctx, text, 3, params)
            c = mp.mpf(text)
            amplitude = mp.exp(mp.mpf(params["logPstar"])-mp.mpf(1)/2)
            def reference(z):
                uh = amplitude/(1+z*z)
                theta = mp.mpf(5)/8*uh
                return dict(theta=theta, z=4*z, theta_z=4*z*theta,
                            z_theta=16*z*z-mp.mpf(5)/12*uh*uh,
                            p=mp.mpf(5)/2*uh*uh)
            for name in NAMES:
                for k, value in enumerate(target["factored_moment_jets"][name].coefficients):
                    expected = mp.diff(lambda z: reference(z)[name], c, k)/mp.factorial(k)
                    if c == 0 and ((name in ("theta", "z_theta", "p") and k % 2)
                                   or (name in ("z", "theta_z") and not k % 2)):
                        # Exact parity at zero; numerical differentiation of
                        # an even/odd function can leave an underflow-scale artifact.
                        expected = mp.mpf(0)
                    lo, hi = endpoints(value)
                    if not lo <= expected <= hi:
                        raise AssertionError(("Independent reference derivative mismatch", text, name, k))
                    count += 1
            samples[text] = target
        whole = reference_targets(ctx, ["-1", "1"], 3, params)
        for target in samples.values():
            for name in NAMES:
                for k, value in enumerate(target["factored_moment_jets"][name].coefficients):
                    lo, hi = endpoints(value)
                    wl, wh = endpoints(whole["factored_moment_jets"][name][k])
                    if not wl <= lo <= hi <= wh:
                        raise AssertionError("Whole-axis reference enclosure missed a point jet")
        # Physical scales are materialized only for this independent
        # pressure-free relation check; the public receipt keeps log factors.
        physical = physical_moment_jets(ctx, samples[".3"])
        ratio = physical["theta"][0]/samples[".3"]["factored_moment_jets"]["theta"][0]
        scale = ctx.exp(samples[".3"]["moment_log_scales"]["theta"])
        if max(endpoints(ratio)[0], endpoints(scale)[0]) > min(endpoints(ratio)[1], endpoints(scale)[1]):
            raise AssertionError("Physical angular moment scale mismatch")
        report = dict(accepted_schedule_sha256=ACCEPTED_SHA,
                      input_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                    for name in (SOURCE, "lei_ren_part1_paper_interval_taylor.py")},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      accepted_parameters=params, paper_equation="(9.3)",
                      independent_derivative_checks=count,
                      whole_axis_point_jet_containment_passed=True,
                      physical_scale_relation_passed=True,
                      point_targets={label: {**target,
                          "reference_Utheta_at_Rh": list(target["reference_Utheta_at_Rh"].coefficients),
                          "factored_moment_jets": {name: list(jet.coefficients)
                                                   for name, jet in target["factored_moment_jets"].items()}}
                                     for label, target in samples.items()},
                      whole_axis_target={**whole,
                          "reference_Utheta_at_Rh": list(whole["reference_Utheta_at_Rh"].coefficients),
                          "factored_moment_jets": {name: list(jet.coefficients)
                                                   for name, jet in whole["factored_moment_jets"].items()}},
                      full_reference_profile_reconstructed=False,
                      actual_transition_moments_generated=False,
                      reference_endpoint_matching_achieved=False,
                      terminal_functional_five_moment_closure=False)
        Path(__file__).with_suffix(".json").write_text(json.dumps(encode(report), indent=2)+"\n", encoding="utf-8")
        print("Reference endpoint targets: independent derivatives", count,
              "whole-axis enclosures passed; actual matching remains open", flush=True)
        return report


if __name__ == "__main__":
    run()
