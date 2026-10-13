"""Focused independent checker for the original Rin-to-infinity pressure values.

The producer composes the accepted r_minus-to-Rp affine source models with
the accepted normalized r_minus-to-Rc models.  This checker verifies that
composition from the persisted accepted archives, retains the positive scale
ratio as an uncertainty, and checks public evaluation without rebuilding any
ancestor owner.
"""

import gzip
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_external_pressure_values as current


def _read_raw():
    return json.loads(gzip.decompress(
        (current.HERE / current.NAME).read_bytes()))


def _iv(c, value):
    """Rebox an interval using its exact directed endpoints."""
    if hasattr(value, "_mpi_"):
        return c.make_mpf(value._mpi_)
    if isinstance(value, dict):
        if value.get("exact_zero", False):
            return c.mpf(0)
        if "lower_exact_mpf_tuple" in value:
            return c.make_mpf((tuple(value["lower_exact_mpf_tuple"]),
                               tuple(value["upper_exact_mpf_tuple"])))
        if "lower" in value and "upper" in value:
            return c.mpf([c.mpf(str(value["lower"])),
                          c.mpf(str(value["upper"]))])
    return c.mpf(value)


def _bits(value):
    if hasattr(value, "_mpi_"):
        return value._mpi_
    return _iv(MPIntervalContext(), value)._mpi_


def _ends(value):
    if hasattr(value, "_mpi_"):
        return (mp.make_mpf(tuple(value._mpi_[0])),
                mp.make_mpf(tuple(value._mpi_[1])))
    if isinstance(value, dict):
        if value.get("exact_zero", False):
            return mp.mpf(0), mp.mpf(0)
        if "lower_exact_mpf_tuple" in value:
            return (mp.make_mpf(tuple(value["lower_exact_mpf_tuple"])),
                    mp.make_mpf(tuple(value["upper_exact_mpf_tuple"])))
        if "lower" in value and "upper" in value:
            return mp.mpf(str(value["lower"])), mp.mpf(str(value["upper"]))
    scalar = mp.mpf(str(value))
    return scalar, scalar


def _upper(value):
    return _ends(value)[1]


def _lower(value):
    return _ends(value)[0]


def _contains(outer, inner):
    olo, ohi = _ends(outer)
    ilo, ihi = _ends(inner)
    return olo <= ilo and ihi <= ohi


def _absolute_upper(c, value):
    lo, hi = current.ep(value)
    return c.mpf(max(abs(lo), abs(hi)))


def _sym_equal(left, right):
    return sy.simplify(sy.expand_log(sy.together(left - right), force=True)) == 0


def independent_ratio_identity():
    """Prove the source scale identity symbolically and numerically."""
    lh, lp, lup, ll, lrc, lmu, lm, lphi, mu, power = sy.symbols(
        "lh lp lup ll lrc lmu lm lphi mu power", real=True)
    logE = (-1 / mu - 3 * sy.exp(-sy.Rational(2, 3) *
                                    (sy.log(4) + lmu))
            + sy.Rational(2, 3) * (sy.log(4) + lmu)
            - sy.log(6) / 2 - 2 * lmu)
    logA = lp + lup + logE - sy.Rational(13, 2) / mu - 13
    logRv = lrc - 2 - 60 * lmu + 13 / mu
    direct = (lh + 2 * lp - power * ll + lm - lphi - logA
              - power * logRv)
    produced = (lh + lp - power * (ll + lrc) - lup + lm - lphi
                + 13 + 2 * power + (60 * power + sy.Rational(4, 3)) * lmu
                + (sy.Rational(15, 2) - 13 * power) / mu
                + 3 * sy.exp(-sy.Rational(2, 3) * (sy.log(4) + lmu))
                - sy.Rational(2, 3) * sy.log(4) + sy.log(6) / 2)
    assert _sym_equal(direct, produced)

    c = MPIntervalContext()
    c.dps = 180
    phi = c.mpf("10.99") * c.sqrt(2 * c.pi) * (c.exp(5) - c.exp(3))
    values = dict(loghb=c.mpf(".37"), logP=c.mpf(".41"),
                  logLambda=c.mpf("1.23"), logRc=c.mpf(".29"),
                  logUp=c.mpf("-.62"), M=c.mpf([".73", ".77"]))
    errors = {}
    for p in (1, 2):
        got = current.ratio_log(c, p, values["loghb"], values["logP"],
                                values["logLambda"], values["logRc"],
                                values["logUp"], values["M"])
        lmu_i = -4 * values["logP"] - c.ln(1000)
        mu_i = c.exp(lmu_i)
        logE_i = (-1 / mu_i - 3 * c.exp(-c.mpf(2) / 3 *
                                          (c.ln(4) + lmu_i))
                   + c.mpf(2) / 3 * (c.ln(4) + lmu_i)
                   - c.ln(6) / 2 - 2 * lmu_i)
        logA_i = values["logP"] + values["logUp"] + logE_i \
            - c.mpf(13) / (2 * mu_i) - 13
        logRv_i = values["logRc"] - 2 - 60 * lmu_i + 13 / mu_i
        direct_i = (values["loghb"] + 2 * values["logP"]
                    - p * values["logLambda"] + c.ln(values["M"])
                    - c.ln(phi) - logA_i - p * logRv_i)
        assert _contains(got, direct_i) or _contains(direct_i, got)
        got_center = sum(_ends(got)) / 2
        direct_center = sum(_ends(direct_i)) / 2
        errors[str(p)] = mp.nstr(abs(got_center - direct_center), 30)
        assert abs(got_center - direct_center) < mp.mpf("1e-150")
    return dict(symbolic_log_ratio_identity=True,
                numerical_nonunit_scale_identity=True,
                powers=(1, 2), center_errors=errors,
                keeps_logE_inverse_mu=True,
                keeps_logRp_as_logRc_minus_two_minus_sixty_logmu=True)


def _assert_source_hashes(raw, owner):
    assert raw["source_family"] == owner.family
    assert raw["input_hashes"] == owner.hashes
    for name, digest in raw["input_hashes"].items():
        assert current.sha(name) == digest, name
    assert raw[current.GATE] is False
    assert all(not raw.get(name, False) for name in current.OPEN)
    assert raw["radial_scope"] == "Rin to infinity"
    assert raw["inner_cut_pressure_moments_not_included"] is True
    assert raw["ancestor_math_suites_not_rerun"] is True
    return True


def _assert_raw_live(raw, owner, samples):
    assert raw["samples"] == current.record(samples)
    assert raw["models"] == current.record(owner.models)
    assert raw["source_scales"] == current.record(owner.scale_data)
    assert raw["original_physical_units"] == owner.upper["original_physical_units"]
    assert raw["source_bindings"] == current.record(owner.bindings)
    return True


def _assert_geometry(owner):
    scale = owner.scale_data
    geometry = owner.upper["source_data"]["geometry"]
    assert scale["common_rminus"] == dict(
        lower_phase="sc/2", upper_phase="sc/2", Rin_phase="sc/4",
        same_original_family=True, shared_endpoint_counted_once=True)
    assert geometry["original_log_r_minus_formula"] == \
        "log(4)-4*logP-1000+hbB*sc/2"
    assert geometry["original_contiguous_17_chart_offsets"] is True
    assert geometry["source_radius_and_phase_Z_independent"] is True
    assert geometry["radius_caps_not_subtracted"] is True
    assert geometry["native_jacobian_positive_on_each_source_cell"] is True
    identities = geometry["exact_symbolic_contiguous_offset_identities"]
    assert len(identities) == 16
    assert all(row["residual"] == "0" for row in identities)
    windows = geometry["exact_17_chart_source_windows"]
    assert len(windows) == 17
    assert len({row["chart"] for row in windows}) == 17
    assert owner.upper["source_data"]["scope"] == dict(
        lower="original r_minus", upper="original Rp",
        active="17 charts r_minus to Rc", compact="Rc to 2Rc",
        quiet="2Rc to Rp", original_shared_Rc_cancelled_before_boxing=True,
        bounded_preRp_is_nonzero_remainder_not_signed_center=True)
    return dict(shared_Rin_phase="sc/4", shared_rminus_phase="sc/2",
                exact_contiguous_chart_offsets=16, active_source_charts=17,
                positive_native_jacobians=True,
                radius_caps_not_subtracted=True)


def _assert_units_and_scale(owner):
    units = owner.upper["original_physical_units"]
    assert units["source_A"] == \
        "Pstar*U0*exp(logE-13/(2*mu)-13)"
    assert units["source_Rv_formula"] == "Rp*exp(13/mu)"
    assert units["source_logRp_formula"] == "logRc-2-60*log(mu)"
    assert units["m13"] == dict(A_power=1, Rv_power=1,
                                inverse_mu_numerator=11,
                                inverse_mu_denominator=2)
    assert units["m15"] == dict(A_power=1, Rv_power=2,
                                inverse_mu_numerator=37,
                                inverse_mu_denominator=2)
    assert units["integer_powers_combined_before_boxing"] is True
    assert units["no_physical_scale_materialization"] is True
    assert "26/mu" not in json.dumps(units, sort_keys=True)

    rows = owner.scale_data["fields"]
    checks = {}
    for key, power in (("m13", 1), ("m15", 2)):
        row = rows[key]
        assert row["original_source_exponents"] == [1, 0, 1, -power]
        M = _iv(owner.c, row["coefficient_norm"])
        assert _lower(M) > 0
        ratio_log = _iv(owner.c, row["exact_positive_ratio_log"])
        assert _upper(ratio_log) < owner.c.mpf(-1000)
        ratio = _iv(owner.c, row["positive_ratio_enclosure"])
        ratio_lo, _ = current.ep(ratio)
        expected = owner.c.exp(-1000)
        assert ratio_lo == 0
        assert ratio._mpi_[1] == expected._mpi_[1]
        assert row["positive_ratio_is_bound_not_selected_zero"] is True
        checks[key] = dict(source_exponents=row["original_source_exponents"],
                           ratio_log_upper=str(_upper(ratio_log)),
                           positive_uncertain_ratio=True)
    assert scale_names(owner) == ['hb', 'F0_anchor', 'Pstar_squared', 'Lambda']
    return dict(units=units, fields=checks,
                common_exact_log_scale_formula=owner.scale_data["exact_log_scale_formula"],
                inverse_mu_coefficient=owner.scale_data["inverse_mu_coefficient"],
                no_scale_materialization=True)


def scale_names(owner):
    return owner.scale_data["original_basis_names"]


def _independent_composed_fields(owner, index):
    """Reproduce compose_models without calling the production composer."""
    c = owner.c
    local = owner.lower["cells"][index]
    parent = owner.upper["models"][index // 2]
    center = _iv(c, local["exact_center"])
    left = _iv(c, local["exact_Z_cell"][0])
    right = _iv(c, local["exact_Z_cell"][1])
    half = (right - left) / 2
    shift = center - _iv(c, parent["center"])
    fields = {}
    for key in ("m13", "m15"):
        old = parent["fields"][key]
        local_field = local["models"][key]
        ratio = _iv(c, owner.scale_data["fields"][key]["positive_ratio_enclosure"])
        a, b = [_iv(c, value)
                 for value in local_field["affine_coefficients_Z0_Z1"]]
        A = (_iv(c, old["center_intercept"])
             + _iv(c, old["center_slope"]) * shift + ratio * a)
        B = _iv(c, old["center_slope"]) + ratio * b
        lower_remainder = ratio * _iv(c, local_field["Taylor_remainder_abs_upper"])
        remainder = _iv(c, old["linear_Taylor_remainder"]) + lower_remainder
        source_error = _iv(c, old["full_nonzero_source_error"])
        alo, ahi = current.ep(A)
        blo, bhi = current.ep(B)
        ambiguity = (ahi - alo) / 2 + (bhi - blo) * half / 2
        error = current.absolute_upper(c, remainder + source_error + ambiguity)
        fields[key] = dict(center_intercept=A, center_slope=B,
                           inherited_outer_Taylor_remainder=
                           _iv(c, old["linear_Taylor_remainder"]),
                           lower_Taylor_remainder=lower_remainder,
                           linear_Taylor_remainder=remainder,
                           coefficient_ambiguity=ambiguity,
                           full_nonzero_source_error=source_error,
                           absolute_error=error)
    return left, right, center, half, fields


def _assert_models(owner):
    c = owner.c
    models = owner.models
    lower_cells = owner.lower["cells"]
    upper_models = owner.upper["models"]
    assert len(models) == 512
    assert len(lower_cells) == 512
    assert len(upper_models) == 256
    previous = None
    signed_count = 0
    positive_error_count = 0
    max_error = c.mpf(0)
    for index, model in enumerate(models):
        local = lower_cells[index]
        left, right, center, half, expected = _independent_composed_fields(owner, index)
        assert _bits(model["left"]) == _bits(left)
        assert _bits(model["right"]) == _bits(right)
        assert _bits(model["center"]) == _bits(center)
        assert _bits(model["half_width"]) == _bits(half)
        assert model["lower_model_index"] == index
        assert model["upper_model_index"] == index // 2
        assert previous is None or _bits(model["left"]) == _bits(previous)
        previous = model["right"]
        assert model["approximation"] == \
            "interval affine enclosure; does not select the physical field"
        for key in ("m13", "m15"):
            got = model["fields"][key]
            want = expected[key]
            for name in ("center_intercept", "center_slope",
                         "inherited_outer_Taylor_remainder",
                         "lower_Taylor_remainder", "linear_Taylor_remainder",
                         "coefficient_ambiguity", "full_nonzero_source_error",
                         "absolute_error"):
                assert _bits(got[name]) == _bits(want[name]), (index, key, name)
            assert got["lower_signed_coefficients_retained"] is True
            assert got["full_pressure_derivative_accuracy_not_claimed"] is True
            assert got["absolute_accuracy_certified"] is True
            assert _upper(got["absolute_error"]) < mp.mpf(".001")
            assert _upper(got["full_nonzero_source_error"]) > 0
            assert (_lower(got["center_intercept"]) < 0 or
                    _upper(got["center_intercept"]) > 0 or
                    _lower(got["center_slope"]) < 0 or
                    _upper(got["center_slope"]) > 0)
            signed_count += 1
            positive_error_count += 1
            max_error = max(max_error, _upper(got["absolute_error"]))
    assert _bits(previous) == _bits(c.mpf(1))
    return dict(models=512, parent_models=256, split_factor=2,
                continuous_full_Z_cover=True, signed_coefficients_retained=signed_count,
                positive_source_errors=positive_error_count,
                maximum_absolute_error_upper=mp.nstr(max_error, 30),
                absolute_error_threshold=".001",
                separate_lower_remainder_and_upper_source_error=True)


def _upper_archive_value_hull(owner, Z):
    """Evaluate accepted outer archive intervals independently of compose_models."""
    c = owner.c
    z = c.mpf(Z)
    values = {"m13": [], "m15": []}
    for model in owner.upper["models"]:
        left, right = _iv(c, model["left"]), _iv(c, model["right"])
        zlo, zhi = current.ep(z)
        if zhi < current.ep(left)[0] or zlo > current.ep(right)[1]:
            continue
        lo = max(zlo, current.ep(left)[0]); hi = min(zhi, current.ep(right)[1])
        offset = c.mpf([lo, hi]) - _iv(c, model["center"])
        for key, field in model["fields"].items():
            error = _absolute_upper(c, _iv(c, field["linear_Taylor_remainder"])
                                    + _iv(c, field["full_nonzero_source_error"])
                                    + _iv(c, field["coefficient_ambiguity"]))
            values[key].append(_iv(c, field["center_intercept"])
                               + _iv(c, field["center_slope"]) * offset
                               + c.mpf([-error, error]))
    assert all(values[key] for key in values)
    return {key: c.mpf([min(current.ep(v)[0] for v in rows),
                        max(current.ep(v)[1] for v in rows)])
            for key, rows in values.items()}


def _independent_sum_view(owner, Z):
    c = owner.c
    outer = _upper_archive_value_hull(owner, Z)
    lower_view = current.lower.evaluate_models(c, owner.lower, Z, True)
    result = {}
    for key in ("m13", "m15"):
        ratio = _iv(c, owner.scale_data["fields"][key]["positive_ratio_enclosure"])
        result[key] = outer[key] + ratio * _iv(c, lower_view["normalized_values"][key])
    return result


def _assert_cross_grid_views(owner):
    points = ("-1", "-.5", "0", ".123", ".371", ".5", "1")
    results = {}
    for point in points:
        live = owner._make_view(point)
        independent = _independent_sum_view(owner, point)
        for key in ("m13", "m15"):
            assert _contains(independent[key], live["full_pressure_values"][key]), \
                (point, key, _ends(independent[key]), _ends(live["full_pressure_values"][key]))
        assert live["source_family"] == owner.family
        assert live["radial_scope"] == "Rin to infinity"
        assert live["inner_cut_pressure_moments_not_included"] is True
        assert live["positive_lower_scale_ratio_retained"] is True
        results[point] = dict(active_models=live["active_models"],
                              independent_sum_contains_composed=True)
    return dict(points=results, end_seam_interior_queries=True,
                independent_outer_plus_normalized_lower=True)


def _public_default_precision(owner):
    before = mp.mp.dps
    with mp.workdps(15):
        public = current.CurrentOriginalN1ExternalPressureValues(require_checked=False)
        packet = public.evaluate(".371")
        view = public.report(packet)
        assert mp.mp.dps == 15
        assert view["source_family"] == owner.family
        assert view["radial_scope"] == "Rin to infinity"
        assert view["inner_cut_pressure_moments_not_included"] is True
        assert view["positive_lower_scale_ratio_retained"] is True
    assert mp.mp.dps == before
    return dict(default_precision=15, precision_restored=True,
                public_evaluate_report=True)


@current.source_precision
def run(owner=None, raw=None):
    began = time.monotonic()
    if owner is None:
        raise ValueError("Hot external-pressure owner required; checker does not rebuild ancestors")
    passed_raw = raw
    saved_raw = _read_raw()
    if passed_raw is not None:
        assert saved_raw == current.record(passed_raw)
    raw = saved_raw
    owner.assert_graph()
    state_before = owner.state()
    samples = {z: owner._make_view(z)
               for z in ("-1", "-.639", "0", ".211", ".371", ".639", "1")}
    owner.assert_graph()
    assert owner.state() == state_before
    _assert_source_hashes(raw, owner)
    _assert_raw_live(raw, owner, samples)
    identity = independent_ratio_identity()
    geometry = _assert_geometry(owner)
    units = _assert_units_and_scale(owner)
    models = _assert_models(owner)
    cross_grid = _assert_cross_grid_views(owner)
    public = _public_default_precision(owner)

    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(current.__file__).name] = current.sha(Path(current.__file__).name)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(
            independent_ratio_identity=identity,
            geometry_and_shared_Rin=geometry,
            physical_units_and_positive_ratio=units,
            all_512_composed_models=models,
            cross_grid_end_seam_interior=cross_grid,
            saved_raw_and_source_bindings=True,
            default_precision_public_query=public,
            inner_cut_pressure_moments_explicitly_excluded=True,
            full_derivative_accuracy_repair_global_stress_and_NS_not_claimed=True),
        input_hashes=hashes,
        **{current.GATE: True},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(receipt, separators=(",", ":"), ensure_ascii=False) + "\n",
        encoding="utf8", newline="\n")
    print("PASS_N1_EXTERNAL_PRESSURE_VALUES", len(owner.models), flush=True)
    return receipt


if __name__ == "__main__":
    run()
