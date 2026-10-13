"""Independent checker for the complete original n1 five-defect values.

The producer adds the signed original cut contribution to the already checked
external pressure models.  This checker keeps that composition explicit: it
recomputes each cut-sector sign enclosure, each model error addition, and the
external-plus-inner query sum.  It never selects an affine midpoint as a
physical value and leaves derivative, repair, and global gates open.
"""

import gzip
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as current


CELLS = current.CELLS
VELOCITY = current.VELOCITY
PRESSURE = current.PRESSURE
MOMENTS = ("m11", "m12", "m13", "m14", "m15")


def _read_raw():
    return json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))


def _iv(c, value):
    """Rebox a live or serialized interval without decimal round-tripping."""
    if hasattr(value, "_mpi_"):
        return c.make_mpf(value._mpi_)
    if isinstance(value, dict):
        if value.get("exact_zero", False):
            return c.mpf(0)
        if "lower_exact_mpf_tuple" in value:
            return c.make_mpf((tuple(value["lower_exact_mpf_tuple"]),
                               tuple(value["upper_exact_mpf_tuple"])))
        if "exact_mpf_tuple" in value:
            t = tuple(value["exact_mpf_tuple"])
            return c.make_mpf((t, t))
        if "lower" in value and "upper" in value:
            return c.mpf([c.mpf(str(value["lower"])),
                          c.mpf(str(value["upper"]))])
        if "value" in value and value.get("value") is not value:
            return _iv(c, value["value"])
    return c.mpf(value)


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
        if "exact_mpf_tuple" in value:
            t = mp.make_mpf(tuple(value["exact_mpf_tuple"]))
            return t, t
        if "lower" in value and "upper" in value:
            return mp.mpf(str(value["lower"])), mp.mpf(str(value["upper"]))
        if "value" in value and value.get("value") is not value:
            return _ends(value["value"])
    scalar = mp.mpf(str(value))
    return scalar, scalar


def _lo(value):
    return _ends(value)[0]


def _hi(value):
    return _ends(value)[1]


def _bits(value):
    if hasattr(value, "_mpi_"):
        return value._mpi_
    return _iv(MPIntervalContext(), value)._mpi_


def _contains(outer, inner):
    olo, ohi = _ends(outer)
    ilo, ihi = _ends(inner)
    return olo <= ilo and ihi <= ohi


def _absolute_upper(c, value):
    lo, hi = current.ep(value)
    return c.mpf(max(abs(lo), abs(hi)))


def _hull(c, values):
    return c.mpf([min(current.ep(value)[0] for value in values),
                  max(current.ep(value)[1] for value in values)])


def _assert_unit_identity():
    """Independently expand the common Phi*A*Rv^p logarithm."""
    lp, lu, lc, lm, p = sy.symbols("lp lu lc lm p", real=True)
    inv_mu = sy.exp(-lm)
    log_e = (-inv_mu - 3 * sy.exp(-sy.Rational(2, 3) *
                                    (sy.log(4) + lm))
             + sy.Rational(2, 3) * (sy.log(4) + lm)
             - sy.log(6) / 2 - 2 * lm)
    direct = (lp + lu + log_e - sy.Rational(13, 2) * inv_mu - 13
              + p * (lc - 2 - 60 * lm + 13 * inv_mu))
    phi = sy.log(sy.Symbol("Phi", positive=True))
    produced = (lp + lu + p * lc - 13 - 2 * p
                + (13 * p - sy.Rational(15, 2)) * sy.exp(-lm)
                - 3 * sy.exp(-sy.Rational(2, 3) * (sy.log(4) + lm))
                - (60 * p + sy.Rational(4, 3)) * lm
                + sy.Rational(2, 3) * sy.log(4) - sy.log(6) / 2)
    assert sy.simplify(sy.expand(direct - produced)) == 0

    c = MPIntervalContext()
    c.dps = 180
    scales = {"logP": c.mpf(".41"), "logUp": c.mpf("-.62"),
              "logRc": c.mpf(".29")}
    log_mu = -4 * scales["logP"] - c.ln(1000)
    mu_i = c.exp(log_mu)
    log_e_i = (-1 / mu_i - 3 * c.exp(-c.mpf(2) / 3 *
                                       (c.ln(4) + log_mu))
               + c.mpf(2) / 3 * (c.ln(4) + log_mu)
               - c.ln(6) / 2 - 2 * log_mu)
    phi_i = c.mpf("10.99") * c.sqrt(2 * c.pi) * (c.exp(5) - c.exp(3))
    errors = {}
    for power in (1, 2):
        got = current.pressure_unit_log(c, scales, power)
        direct_i = (c.ln(phi_i) + scales["logP"] + scales["logUp"]
                    + log_e_i - c.mpf(13) / (2 * mu_i) - 13
                    + power * (scales["logRc"] - 2 - 60 * log_mu
                               + 13 / mu_i))
        assert _contains(got, direct_i) or _contains(direct_i, got)
        got_center = sum(_ends(got)) / 2
        direct_center = sum(_ends(direct_i)) / 2
        error = abs(got_center - direct_center)
        assert error < mp.mpf("1e-150")
        errors[str(power)] = mp.nstr(error, 25)
    return dict(symbolic=True, nonunit_scales=True, powers=(1, 2),
                center_errors=errors, inverse_mu_terms_combined=True,
                no_radius_or_N_materialized=True,
                phi_definition="10.99*sqrt(2*pi)*(exp(5)-exp(3))")


def _term_expected(c, term, unit):
    source = term["source_term"]
    assert term["source_exponents"] == source["source_exponents"]
    coefficient = current.interval(c, source["ordinary_Z_derivatives"][0])
    source_log = current.interval(c, source["combined_source_log"])
    magnitude = _absolute_upper(c, coefficient)
    if current.ep(magnitude)[1] == 0:
        return dict(exact=True, value=c.mpf(0), absolute=c.mpf(0), log=None,
                    coefficient=coefficient)
    log_cap = source_log - unit + c.ln(magnitude)
    assert current.ep(log_cap)[1] < -1000
    eps = c.exp(-1000)
    eps_hi = current.ep(eps)[1]
    lo, hi = current.ep(coefficient)
    if lo >= 0:
        value = c.mpf([0, eps_hi])
        sign = "nonnegative"
    elif hi <= 0:
        value = c.mpf([-eps_hi, 0])
        sign = "nonpositive"
    else:
        value = c.mpf([-eps_hi, eps_hi])
        sign = "straddling"
    return dict(exact=False, value=value, absolute=eps, log=log_cap,
                coefficient=coefficient, sign=sign)


def _assert_inner_term(c, cell, tails):
    """Recompute every signed inner cut sector from its source record."""
    scales = tails.scale_data
    result = {}
    for key, power in (("m13", 1), ("m15", 2)):
        actual = cell["view"]["physical_cutoff_moment_values_at_Rin"][key]
        assert actual["ordinary_Z_derivative_order"] >= 1
        unit = current.pressure_unit_log(c, scales, power)
        terms = []
        for got in cell["inner_pressure"][key]["terms"]:
            want = _term_expected(c, got, unit)
            assert got["source_exponents"] == got["source_term"]["source_exponents"]
            assert got["exact_coefficient_zero"] is want["exact"]
            assert _bits(got["normalized_value"]) == _bits(want["value"])
            assert _bits(got["absolute_upper"]) == _bits(want["absolute"])
            if want["exact"]:
                assert "normalized_log_absolute_upper" not in got
            else:
                assert _bits(got["normalized_log_absolute_upper"]) == _bits(want["log"])
                assert got["actual_source_coefficient_sign_retained"] is True
                assert got["nonzero_bound_not_selected_zero"] is True
            terms.append(want)
        value = sum((row["value"] for row in terms), c.mpf(0))
        absolute = sum((row["absolute"] for row in terms), c.mpf(0))
        record = cell["inner_pressure"][key]
        assert _bits(record["unit_log"]) == _bits(unit)
        assert _bits(record["normalized_signed_value"]) == _bits(value)
        assert _bits(record["absolute_upper"]) == _bits(absolute)
        assert record["original_signed_Z1_sectors"] == actual
        assert record["ordinary_inner_Z1_source_preserved"] is True
        assert record["full_pressure_derivative_accuracy_not_claimed"] is True
        result[key] = dict(terms=len(terms), exact_zero=sum(row["exact"] for row in terms),
                           nonzero=len(terms) - sum(row["exact"] for row in terms),
                           absolute_upper=mp.nstr(_hi(record["absolute_upper"]), 20),
                           source_unit=mp.nstr(_lo(unit), 20))
    return result


def _assert_source_cells(owner, raw):
    assert len(owner.source_cells) == 4
    assert len(raw["source_cells"]) == 4
    assert owner.cutoff.family == owner.tails.family == owner.family
    assert all(name in owner.hashes for name in owner.cutoff.hashes)
    assert all(name in owner.hashes for name in owner.tails.hashes)
    family = owner.family
    assert all(name in family for name in
               ("actual_five_defect_family_sha256", "implicit_source_sha256",
                "datum_enclosure_sha256"))
    summaries = []
    for index, cell in enumerate(owner.source_cells):
        assert cell["ends"] == list(CELLS[index])
        assert raw["source_cells"][index] == current.record(cell)
        view = cell["view"]
        assert view["source_family"] == family
        assert view["query"]["chart"] == "first"
        assert view["query"].get("selected_fraction") in ("1", 1, None)
        radial = view["radial_cutoff"]
        assert radial["keep_fraction"] == "1/2"
        assert radial["zero_fraction"] == "3/4"
        assert radial["Z_independent"] is True
        assert radial["zero_radius_is_geometric_mean_of_Rkeep_Rin"] is True
        assert radial["zero_before_arithmetic_mean_by_AM_GM"] is True
        assert len(radial["value"]) == 2
        assert all(_iv(MPIntervalContext(), row["value"])._mpi_ ==
                   MPIntervalContext().mpf(0)._mpi_ for row in radial["value"])
        for key in VELOCITY:
            support = view["full_velocity_support_moments"][key]
            budget = view["full_velocity_support_moment_budgets"][key]
            assert support["ordinary_Z_derivative_order"] >= 1
            assert isinstance(support["terms"], list) and support["terms"]
            assert set(budget) >= {"Z0", "Z1"}
            assert isinstance(budget["Z0"], dict) and isinstance(budget["Z1"], dict)
        for key in PRESSURE:
            actual = view["physical_cutoff_moment_values_at_Rin"][key]
            assert actual["ordinary_Z_derivative_order"] >= 1
            assert actual["terms"]
        assert set(view["physical_cutoff_moment_values_at_Rin"]) == set(MOMENTS)
        assert view["physical_own_M_Z2"]["ordinary_Z_derivative_order"] >= 2
        assert set(view["physical_own_mean_K_V"]) == {"mean1", "K1", "V1"}
        assert view["physical_recovered_pressure_on_common_inner"]["ordinary_Z_derivative_order"] >= 1
        for flag in ("full_velocity_support_fixed_before_original_inner_exit",
                     "mean_V_pressure_recovered_not_multiplied_by_cutoff",
                     "original_hb_inverse_and_inverse_squared_cutoff_derivatives_retained",
                     "reconstructed_V_satisfies_own_n1_radial_recovery",
                     "untouched_core_and_keep_branch_preserved",
                     "leading_support_tail_two_moments_repair_n2_and_stress_still_open"):
            assert view[flag] is True
        assert view["unresolved_leading_pressure_forcing_tails"]["both_tails_uncomputed_and_not_set_zero"] is True
        summaries.append(dict(cell=cell["ends"], velocity_moments=list(VELOCITY),
                              pressure_moments=list(PRESSURE),
                              signed_inner=_assert_inner_term(owner.c, cell, owner.tails)))
    return dict(cells=4, exact_Z_cells=[list(x) for x in CELLS], summaries=summaries,
                cut_before_own_M_V_P=True, common_family=True,
                independent_P0_family_hashes=True)


def _assert_models(owner):
    c = owner.c
    assert len(owner.models) == 512
    assert len(owner.tails.models) == 512
    previous = None
    max_error = c.mpf(0)
    positive_errors = 0
    signed_fields = 0
    for index, model in enumerate(owner.models):
        parent = owner.tails.models[index]
        cell = owner.source_cells[index // 128]
        assert model["inner_cut_source_cell"] == cell["ends"]
        assert model["all_inner_and_external_pressure_parts_included"] is True
        # The composed model is serialized, while the accepted tail model
        # remains live IV data.  Compare their exact directed MPI endpoints.
        for geometry_key in ("left", "right", "center", "half_width"):
            assert _bits(model[geometry_key]) == _bits(parent[geometry_key])
        assert model["lower_model_index"] == parent["lower_model_index"]
        assert model["upper_model_index"] == parent["upper_model_index"]
        assert previous is None or _bits(model["left"]) == _bits(previous)
        previous = model["right"]
        assert model["approximation"] == parent["approximation"]
        for key in PRESSURE:
            old = parent["fields"][key]
            got = model["fields"][key]
            extra = current.interval(c, cell["inner_pressure"][key]["absolute_upper"])
            inherited = current.interval(c, old["full_nonzero_source_error"])
            expected_source = inherited + extra
            expected_error = current.interval(c, old["absolute_error"]) + extra
            expected_sign = current.interval(c, cell["inner_pressure"][key]["normalized_signed_value"])
            assert _bits(got["center_intercept"]) == _bits(old["center_intercept"])
            assert _bits(got["center_slope"]) == _bits(old["center_slope"])
            assert _bits(got["inherited_external_source_error"]) == _bits(inherited)
            assert _bits(got["inner_cut_nonzero_absolute_error"]) == _bits(extra)
            assert _bits(got["full_nonzero_source_error"]) == _bits(expected_source)
            assert _bits(got["absolute_error"]) == _bits(expected_error)
            assert _bits(got["inner_cut_sign_enclosure"]) == _bits(expected_sign)
            assert got["affine_coefficients_describe_external_part_only"] is True
            assert got["lower_signed_coefficients_retained"] is True
            assert got["full_pressure_derivative_accuracy_not_claimed"] is True
            assert got["absolute_accuracy_certified"] is (current.ep(expected_error)[1] < c.mpf(".001"))
            got_error = _iv(c, got["absolute_error"])
            got_source_error = _iv(c, got["full_nonzero_source_error"])
            got_sign = _iv(c, got["inner_cut_sign_enclosure"])
            assert current.ep(got_error)[1] < c.mpf(".001")
            assert current.ep(got_source_error)[1] > 0
            max_error = max(max_error, current.ep(got_error)[1])
            positive_errors += 1
            if current.ep(got_sign)[0] < 0 or current.ep(got_sign)[1] > 0:
                signed_fields += 1
    assert previous is not None and current.ep(_iv(c, previous))[1] == 1
    return dict(models=512, exact_external_error_addition=True,
                continuous_full_Z_cover=True, positive_source_errors=positive_errors,
                signed_inner_fields=signed_fields,
                maximum_absolute_error_upper=mp.nstr(max_error, 25),
                absolute_error_threshold=".001")


def _independent_view(owner, Z):
    """Compose accepted external values with inner cell values independently."""
    c = owner.c
    z = c.mpf(Z)
    tail = owner.tails._make_view(Z)
    zlo, zhi = current.ep(z)
    active = []
    for cell in owner.source_cells:
        left, right = c.mpf(cell["ends"][0]), c.mpf(cell["ends"][1])
        if zhi < current.ep(left)[0] or zlo > current.ep(right)[1]:
            continue
        active.append(cell)
    assert active
    inner = {}
    for key in PRESSURE:
        inner[key] = _hull(c, [current.interval(c, cell["inner_pressure"][key]["normalized_signed_value"])
                               for cell in active])
    return {key: tail["full_pressure_values"][key] + inner[key] for key in PRESSURE}, active


def _assert_query_composition(owner, raw):
    points = ("-1", "-.5", "0", ".123", ".371", ".5", "1")
    seams = 0
    for point in points:
        live = owner._make_view(point)
        expected, active = _independent_view(owner, point)
        if len(active) == 2:
            seams += 1
        assert live["source_family"] == owner.family
        assert live["moment_order"] == list(MOMENTS)
        assert live["active_source_cells"] == [cell["ends"] for cell in active]
        assert live["full_radial_support_included"] is True
        assert live["original_own_M_V_P_and_independent_P0_preserved"] is True
        assert live["original_m15_external_minus_sign_applied_once"] is True
        assert live["velocity_values_remain_factored_not_zeroed"] is True
        assert live["full_pressure_derivative_accuracy_not_claimed"] is True
        for key in PRESSURE:
            assert _bits(live["inner_cut_pressure_values"][key]) == _bits(
                _hull(owner.c, [current.interval(owner.c, cell["inner_pressure"][key]["normalized_signed_value"])
                                for cell in active]))
            assert _bits(live["external_pressure_values"][key]) == _bits(owner.tails._make_view(point)["full_pressure_values"][key])
            assert _bits(live["full_pressure_moment_values"][key]) == _bits(expected[key])
        for key in VELOCITY:
            assert len(live["full_velocity_moment_source_values"][key]) == len(active)
            for row, cell in zip(live["full_velocity_moment_source_values"][key], active):
                assert row["exact_Z_cell"] == cell["ends"]
                assert row["signed_physical_sectors"] == cell["view"]["physical_cutoff_moment_values_at_Rin"][key]
                assert row["signed_source_budgets"] == cell["view"]["full_velocity_support_moment_budgets"][key]
    assert raw["samples"] == current.record({z: owner._make_view(z)
                                               for z in ("-1", "-.639", "0", ".211", ".371", ".639", "1")})
    return dict(points=list(points), seam_queries=seams,
                independent_external_plus_inner=True,
                m15_external_minus_sign_applied_once=True,
                velocity_sectors_retained=True)


def _assert_hashes_and_scope(owner, raw):
    assert raw["source_family"] == owner.family
    assert raw["input_hashes"] == owner.hashes
    for name, digest in raw["input_hashes"].items():
        assert current.sha(name) == digest, name
    assert raw["source_bindings"] == current.record(owner.bindings)
    assert raw["full_radial_support_included"] is True
    assert raw["three_velocity_moments_remain_signed_physical_sectors"] is True
    assert raw["both_pressure_cut_and_external_parts_included"] is True
    assert raw["inner_cut_first_derivatives_not_full_pressure_accuracy"] is True
    assert raw[current.GATE] is False
    for name in current.OPEN:
        assert raw.get(name, False) is False
    for name in ("coefficient_F1_Uz1_P1_solved", "positive_order_five_moment_repair_closed",
                 "common_holomorphic_extension_assumption14_1_certified", "temporal_recursion",
                 "regional_admissible_stress_certified", "flat_remainder_certified",
                 "full_corrected_Navier_Stokes_solution", "n1_full_inner_interval_solution_certified",
                 "n1_radial_Taylor_remainder_certified", "completed_leading_supports_outside_Rin_certified",
                 "n1_Volterra_iteration_prefix_computed", "complete_signed_pressure_tails_certified",
                 "n1_second_and_third_repair_closed", "global_physical_pressure_defect_installed",
                 "uniform_full_pressure_derivative_accuracy_certified", "global_stress_certified",
                 "global_NS_solution_certified"):
        assert raw.get(name, False) is False
    return dict(source_hashes=len(raw["input_hashes"]), same_family=True,
                cut_before_own_M_V_P=True, independent_P0=True,
                all_later_gates_false=True)


def _public_default_precision(owner):
    before = mp.mp.dps
    state = owner.state()
    with mp.workdps(15):
        packet = owner.evaluate(".371")
        view = owner.report(packet)
        assert mp.mp.dps == 15
        assert view["source_family"] == owner.family
        assert set(view["full_pressure_moment_values"]) == set(PRESSURE)
        assert view["original_m15_external_minus_sign_applied_once"] is True
        assert owner.state() == state
    assert mp.mp.dps == before
    return dict(default_precision=15, precision_restored=True,
                public_evaluate_report=True, owner_state_restored=True,
                no_owner_reload=True)


@current.source_precision
def run(owner=None, raw=None):
    """Run against the hot checked parents; root supplies owner/raw."""
    began = time.monotonic()
    if owner is None:
        raise ValueError("Hot complete-five-defect owner required")
    saved = _read_raw()
    if raw is not None:
        assert saved == current.record(raw)
    raw = saved
    owner.assert_graph()
    state = owner.state()
    assert raw["source_cells"] == current.record(owner.source_cells)
    assert raw["pressure_models"] == current.record(owner.models)
    assert raw["samples"] == current.record({z: owner._make_view(z)
                                               for z in ("-1", "-.639", "0", ".211", ".371", ".639", "1")})
    owner.assert_graph()
    assert owner.state() == state
    scope = _assert_hashes_and_scope(owner, raw)
    cells = _assert_source_cells(owner, raw)
    models = _assert_models(owner)
    query = _assert_query_composition(owner, raw)
    identity = _assert_unit_identity()
    public = _public_default_precision(owner)
    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(current.__file__).name] = current.sha(Path(current.__file__).name)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(all_passed=True, source_family=owner.family,
                   checks=dict(independent_unit_identity=identity,
                               source_cells_and_signed_sectors=cells,
                               complete_512_pressure_models=models,
                               independent_query_composition=query,
                               saved_raw_source_bindings=scope,
                               default_precision_public_query=public,
                               derivative_repair_global_and_NS_not_claimed=True),
                   input_hashes=hashes, **{current.GATE: True},
                   **dict.fromkeys(current.OPEN, False),
                   execution_seconds=time.monotonic() - began)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(receipt, separators=(",", ":"), ensure_ascii=False) + "\n",
        encoding="utf8", newline="\n")
    print("PASS_N1_FIVE_DEFECT_VALUES", len(owner.source_cells), len(owner.models), flush=True)
    return receipt


if __name__ == "__main__":
    run()
