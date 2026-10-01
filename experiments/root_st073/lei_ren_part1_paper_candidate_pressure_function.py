"""Functional accepted pressure datum: directed jets at real axial centers.

The same fourteen accepted radial masses are retained. Only the axial center
changes. Variable-beta flatten pressure is enclosed by a uniform Cauchy bound,
not fitted and not discarded. This does not certify a new candidate collar.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_pressure_axis_jets import _interval_from_exact
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).resolve().parent
DATUM = "lei_ren_part1_paper_candidate_pressure_axis_jets_refined.json"
ACCEPTED_SHA = "736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4"


def load_datum(path=HERE/DATUM):
    """Check the complete transitive receipt dependencies before using masses."""
    path = Path(path)
    checked = {}

    def visit(filename):
        filename = filename.resolve()
        if filename in checked:
            return
        raw = filename.read_bytes()
        checked[filename] = hashlib.sha256(raw).hexdigest()
        if filename.suffix == ".json":
            record = json.loads(raw)
            for name, digest in record.get("input_hashes", {}).items():
                dependency = filename.parent/name
                actual = hashlib.sha256(dependency.read_bytes()).hexdigest()
                if actual != digest:
                    raise ValueError("Pressure datum dependency changed: " + name)
                visit(dependency)

    visit(path)
    record = json.loads(path.read_bytes())
    if record["accepted_schedule_sha256"] != ACCEPTED_SHA:
        raise ValueError("Unexpected accepted pressure datum")
    if record["stage_count_total"] != 14 or not record["all_14_true_pressure_stages_included"]:
        raise ValueError("Incomplete accepted pressure datum")
    if record["pressure_units"] != "physical_P":
        raise ValueError("Pressure units are not physical_P")
    if record["pressure_parameter_order_truncated"] or record["finite_pressure_quadrature_approximation_used"]:
        raise ValueError("A truncated or quadrature-fit pressure datum is not accepted")
    return record, {str(p.relative_to(HERE)): digest for p, digest in checked.items()}


def q_jets(ctx, center, order):
    # Center is real, possibly an interval crossing zero. Use the exact square
    # dependency to keep 1+center^2 bounded away from zero on the whole axis.
    a = 1 + center**2
    b = 2*center
    coefficients = [1/a**2]
    previous = ctx.mpf(0)
    for n in range(order):
        current = coefficients[n]
        coefficients.append(-(b*(n+2)*current+(n+3)*previous)/(a*(n+1)))
        previous = current
    return coefficients


def pressure_jets(ctx, Z, order=127, datum=None):
    """Ordinary Taylor coefficients at every real center in Z, through order.

    If Z is an interval, each coefficient encloses the corresponding derivative
    at every center in that interval. It is not a Taylor expansion about an
    interval viewed as one point.
    """
    if not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    z = ctx.mpf(Z)
    zl, zh = endpoints(z)
    if zl < -1 or zh > 1:
        raise ValueError("Require real axial centers in [-1,1]")
    if datum is None:
        datum, _ = load_datum()
    read = lambda name: _interval_from_exact(ctx, datum[name])
    mass0 = read("fixed_beta0_mass_upper_normalized")
    mass2 = read("fixed_beta2_mass_upper_normalized")
    flatten = read("flatten_true_mass_upper_normalized")
    scale = read("Pstar_squared")
    if any(endpoints(m)[0] < 0 for m in (mass0, mass2, flatten)):
        raise ValueError("Pressure masses must be nonnegative")
    rho = ctx.mpf(".25")
    # For alpha=beta/2 in [0,1], factor q=1+alpha*z^2 at its two
    # imaginary roots. On any radius-rho disk about a real center c,
    # |q| >= (sqrt(1+alpha*c^2)-sqrt(alpha)*rho)^2 >= (1-rho)^2.
    q_lower = (1-rho)**2
    qcoef = q_jets(ctx, z, order)
    rows = []
    for n, coefficient in enumerate(qcoef):
        if n == 0:
            remainder = flatten  # q(real)^(-2) <= 1; the measure is positive.
        else:
            bound = endpoints(flatten/q_lower**2/rho**n)[1]
            remainder = ctx.mpf([-bound, bound])
        fixed = mass2*coefficient + (mass0 if n == 0 else 0)
        rows.append(-scale*(fixed+remainder))
    return dict(center_Z=z, ordinary_Taylor_order=order,
                physical_pressure_coefficients=rows,
                q_inverse_squared_coefficients=qcoef,
                flatten_Cauchy_radius=rho, flatten_q_modulus_lower=q_lower,
                accepted_schedule_sha256=ACCEPTED_SHA,
                pressure_units="physical_P", all_14_stages_included=True,
                terminal_functional_five_moment_closure=False,
                candidate_collar_pressure_compatibility_certified=False)


def run():
    datum, hashes = load_datum()
    ctx = MPIntervalContext()
    ctx.dps = 260
    with mp.workdps(300):
        at_center = pressure_jets(ctx, ".3", 127, datum)
        zero = q_jets(ctx, ctx.mpf(0), 16)
        for n, coefficient in enumerate(zero):
            exact = mp.mpf(0) if n % 2 else mp.mpf((-1)**(n//2)*(n//2+1))
            lo, hi = endpoints(coefficient)
            if not lo <= exact <= hi:
                raise AssertionError("Exact q-series fixture failed at order " + str(n))
        for n, coefficient in enumerate(at_center["physical_pressure_coefficients"]):
            previous = _interval_from_exact(ctx, datum["Taylor_rows"][n]["physical_pressure_coefficient"])
            lo, hi = endpoints(coefficient)
            pl, ph = endpoints(previous)
            if max(lo, pl) > min(hi, ph):
                raise AssertionError("Original center pressure jets disagree at order " + str(n))
        report = dict(input_hashes=hashes,
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      accepted_schedule_sha256=ACCEPTED_SHA,
                      samples={label: pressure_jets(ctx, center, 3, datum)
                               for label, center in (("left", "-1"), ("zero", "0"),
                                                     ("center", ".3"), ("right", "1"),
                                                     ("whole_axis", ["-1", "1"]))},
                      real_center_interval_supported=True,
                      exact_zero_center_q_fixture_count=17,
                      original_center_pressure_overlap_count=128,
                      all_original_pressure_stages_retained=True,
                      functional_pressure_datum_evaluable=True,
                      candidate_core_whole_axis_generated=False,
                      candidate_collar_pressure_compatibility_certified=False,
                      terminal_functional_five_moment_closure=False,
                      original_parameter_errors_enclosed=False)
        Path(__file__).with_suffix(".json").write_text(json.dumps(encode(report), indent=2)+"\n", encoding="utf-8")
        print("Functional accepted pressure jets: whole real axis; center overlaps128; exact q fixtures17", flush=True)
        return report


if __name__ == "__main__":
    run()
