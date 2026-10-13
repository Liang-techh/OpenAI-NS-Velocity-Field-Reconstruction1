"""Independent checker for the original r_minus-to-infinity pressure source.

The producer supplies magnitude enclosures on the 17 actual source charts and
retains the compact-band/quiet archives.  This checker verifies those records
against the live owner, reconstructs the radial-free identities on an
independent smooth fixture, and keeps every omitted contribution as a positive
source error.  It does not certify full tails, repair, or global equations.
"""

import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source as current
import lei_ren_part1_paper_compliant_current_original_n1_preRc_pressure_endpoint as endpoint


def _read_raw():
    return json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))


def _interval_endpoints(value):
    if hasattr(value, "coefficients"):
        value = value[0]
    if hasattr(value, "_mpi_"):
        # Rebox the exact directed endpoints.  ``.a``/``.b`` can be ordinary
        # mpf views at the caller's ambient precision and would discard bits
        # when a checked source is inspected from the default 15-digit state.
        return (mp.make_mpf(tuple(value._mpi_[0])),
                mp.make_mpf(tuple(value._mpi_[1])))
    if isinstance(value, dict):
        if value.get("exact_zero", False):
            return mp.mpf(0), mp.mpf(0)
        if "log_absolute_upper" in value:
            return _interval_endpoints(value["log_absolute_upper"])
        if "lower_exact_mpf_tuple" in value:
            return (mp.make_mpf(tuple(value["lower_exact_mpf_tuple"])),
                    mp.make_mpf(tuple(value["upper_exact_mpf_tuple"])))
        if "lower" in value and "upper" in value:
            return mp.mpf(str(value["lower"])), mp.mpf(str(value["upper"]))
    scalar = mp.mpf(str(value))
    return scalar, scalar


def _upper(value):
    return _interval_endpoints(value)[1]


def _absolute_upper(value):
    if hasattr(value, "log"):
        return mp.mpf(0) if value.log is None else mp.exp(_upper(value.log))
    if isinstance(value, dict) and "log_absolute_upper" in value:
        return (mp.mpf(0) if value.get("exact_zero", False)
                else mp.exp(_upper(value["log_absolute_upper"])))
    return _upper(value)


def _interval_bits(value):
    if hasattr(value, "_mpi_"):
        return value._mpi_
    if isinstance(value, dict) and "lower_exact_mpf_tuple" in value:
        return (tuple(value["lower_exact_mpf_tuple"]),
                tuple(value["upper_exact_mpf_tuple"]))
    return None


def _canon(value):
    return sy.factor(sy.cancel(sy.together(sy.expand_log(sy.sympify(value), force=True))))


def _equal(left, right):
    return _canon(left - right) == 0


def _direct_fields(mean, radius, z, delta):
    d = 1 - z**2
    L = 1 - delta * z**2
    U = sy.cancel(mean + radius * sy.diff(mean, radius))
    q = sy.cancel((2 * z * U - (1 - delta) * z * mean - d * sy.diff(mean, z)) / L)
    qz = sy.diff(q, z)
    qr = sy.diff(q, radius)
    qy = radius * qr
    S = sy.cancel((q + (1 - delta) * z * qz / 2 + qy) / L
                  + q * (q / 2 + qy)
                  + U * (-2 * z * q + d * qz - 2 * z * qy) / L)
    n1p = sy.cancel(-S / 2 + 2 * qr + radius * sy.diff(q, radius, 2))
    return dict(U=U, q=q, qz=qz, qr=qr, qy=qy, S=S, n1p=n1p)


def independent_pressure_fixture():
    """Verify the pressure IBP identities and shared-Rc unit algebra independently."""
    R, Z = sy.symbols("R Z", positive=True, real=True)
    delta = sy.Rational(2, 13)
    mean = (sy.Rational(3, 11) * (Z + 2 * Z**2 + Z**3) / R
            + sy.Rational(2, 17) * Z**2 * (R - 1)**2 / R
            + sy.Rational(1, 19) * Z**3 * (R - sy.Rational(3, 2)) / R)
    fields = _direct_fields(mean, R, Z, delta)

    def boundary(value):
        return sy.cancel(value.subs(R, 2) - value.subs(R, 1))

    i0 = {}
    i1 = {}
    integrands = dict(
        qz=fields["qz"], Uqz=fields["U"] * fields["qz"],
        U2=fields["U"]**2, UUz=fields["U"] * sy.diff(fields["U"], Z),
        q=fields["q"], q2=fields["q"]**2)
    for name, value in integrands.items():
        i0[name] = sy.integrate(value, (R, 1, 2))
        i1[name] = sy.integrate(R * value, (R, 1, 2))

    L = 1 - delta * Z**2
    d = 1 - Z**2
    i0S = (boundary(R * fields["q"]) / L
           + (1 - delta) * Z * i0["qz"] / (2 * L)
           + boundary(R * fields["q"]**2) / 2
           + d * i0["Uqz"] / L
           - 2 * Z**2 * boundary(R * fields["U"]**2) / L**2
           - 2 * delta * Z**2 * i0["U2"] / L**2
           + 2 * Z * d * i0["UUz"] / L**2)
    i1S = (boundary(R**2 * fields["q"]) / L - i1["q"] / L
           + (1 - delta) * Z * i1["qz"] / (2 * L)
           + boundary(R**2 * fields["q"]**2) / 2 - i1["q2"] / 2
           + d * i1["Uqz"] / L
           - 2 * Z**2 * boundary(R**2 * fields["U"]**2) / L**2
           + 2 * (1 - delta) * Z**2 * i1["U2"] / L**2
           + 2 * Z * d * i1["UUz"] / L**2)
    direct_s0 = sy.integrate(fields["S"], (R, 1, 2))
    direct_s1 = sy.integrate(R * fields["S"], (R, 1, 2))
    direct0 = sy.integrate(fields["n1p"], (R, 1, 2))
    direct1 = sy.integrate(R * fields["n1p"], (R, 1, 2))
    assert _equal(i0S, direct_s0)
    assert _equal(i1S, direct_s1)
    assert _equal(-direct_s0 / 2 + boundary(fields["q"] + R * fields["qr"]), direct0)
    assert _equal(-direct_s1 / 2 + boundary(R**2 * fields["qr"]), direct1)

    # Shared Rc cancellation in Rc^p/(Phi*Aphys*Rv^p), with
    # log(Rv)=log(Rc)-2-60 log(mu)+13/mu.
    LP, LUp, LE, LPhi, LRc, Lmu, mu = sy.symbols(
        "logP logUp logE logPhi logRc logmu mu", real=True, positive=True)
    logA = LP + LUp + LE - sy.Rational(13, 2) / mu - 13
    logRv = LRc - 2 - 60 * Lmu + 13 / mu
    units = {}
    for power in (1, 2):
        lhs = power * LRc - (LPhi + logA + power * logRv)
        rhs = (-LP - LUp - LE - LPhi + 13
               + power * (2 + 60 * Lmu)
               + (sy.Rational(13, 2) - 13 * power) / mu)
        assert _equal(lhs, rhs)
        units[power] = sy.sstr(rhs)

    assert fields["U"] != 0 and fields["qz"] != 0
    assert boundary(fields["q"] + R * fields["qr"]) != 0
    assert boundary(R**2 * fields["qr"]) != 0
    return dict(
        direct_reduced_S_I0_and_I1=True,
        direct_endpoint_IBP=True,
        nonlinear_U_qZ_and_U_UZ_retained=True,
        viscous_qR_and_qRR_retained=True,
        finite_endpoint_terms_nonzero=True,
        shared_Rc_unit_ratio_exact=True,
        unit_log_formula_by_power=units,
        fixture_delta=str(delta),
    )


def _assert_source_hashes(raw, owner):
    assert raw["source_family"] == owner.family
    assert raw[current.GATE] is False
    assert all(not raw.get(name, False) for name in current.OPEN)
    assert raw["input_hashes"] == owner.hashes
    for name, digest in raw["input_hashes"].items():
        assert current.sha(name) == digest, name
    return True


def _expected_bindings():
    return {name: current.pulse.end.coefficients.function_binding(fn) for name, fn in (
        ("active_reduced_pressure_caps", current.preRc_coefficients),
        ("shared_Rc_unit_ratio", current.Rc_unit_log),
        ("complete_preRp_source_data", current.build_data),
        ("checked_chart_views", current.checked_views),
        ("original_endpoint_caps", endpoint.endpoint_caps),
        ("original_geometry", endpoint.geometry_proof),
        ("original_exact_beyond_Rv_support", current.support.source_proof),
        ("original_reduced_pressure", current.reduced.reduced_S_integrals))}


def _assert_raw_live(raw, owner, samples, models):
    assert raw["samples"] == current.record(samples)
    assert raw["models"] == current.record(models)
    assert raw["source_data"] == current.record(owner.data)
    assert raw["source_bindings"] == current.record(_expected_bindings())
    return True


def _assert_chart_views(owner):
    data_views = owner.data["actual_chart_views"]
    expected = {(tuple(ends), chart) for ends in current.pulse.CELLS
                for chart in owner.fields.charts}
    actual = set()
    source_cells = 0
    source_cells_per_Z = {tuple(ends): 0 for ends in current.pulse.CELLS}
    for item in data_views:
        ends = tuple(item["ends"])
        chart = item["chart"]
        view = item["view"]
        key = (ends, chart)
        assert key not in actual
        actual.add(key)
        assert view["source_family"] == owner.family
        assert view["whole_cell_range_not_endpoint_value"] is True
        assert view["exact_original_source_cell_union"] is True
        assert view["physical_common_scale_Pstar_applied_once"] is True
        partition = owner.fields.partitions[chart]
        assert current.record(view["exact_native_cell"]) == current.record(
            [partition[0], partition[-1]])
        count = len(partition) - 1
        source_cells += count
        source_cells_per_Z[ends] += count
    assert actual == expected and len(actual) == 68
    assert source_cells == 4 * 57
    assert all(count == 57 for count in source_cells_per_Z.values())
    assert len(owner.views) == 68
    return dict(chart_views=68, source_partition_cells=source_cells,
                source_partition_cells_total=source_cells,
                source_partition_cells_per_Z={str(ends): count
                                              for ends, count in source_cells_per_Z.items()},
                four_Z_cells=len(current.pulse.CELLS), all_chart_coverage=True)


def _assert_log_record(value):
    if isinstance(value, dict):
        if "log_absolute_upper" in value:
            if value.get("exact_zero", False):
                return
            assert value["log_absolute_upper"] is not None
        for item in value.values():
            _assert_log_record(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _assert_log_record(item)


def _independent_majorant_coefficients():
    """Return the reduced-S coefficients proved from the source inequalities."""
    q, qz, u, uz = sy.symbols("q q_Z U U_Z", nonnegative=True)
    common = qz + 2 * u * qz + 8 * u * uz
    s0 = sy.expand(4 * q + q**2 + 20 * u**2 + common)
    s1 = sy.expand(6 * q + sy.Rational(3, 2) * q**2 + 24 * u**2 + common)
    assert sy.expand(s0 - (4 * q + q**2 + qz + 2 * u * qz + 20 * u**2 + 8 * u * uz)) == 0
    assert sy.expand(s1 - (6 * q + sy.Rational(3, 2) * q**2 + qz + 2 * u * qz + 24 * u**2 + 8 * u * uz)) == 0
    return dict(
        S0={"q": "4", "q_squared": "1", "q_Z": "1", "U_q_Z": "2",
            "U_squared": "20", "U_U_Z": "8"},
        S1={"q": "6", "q_squared": "3/2", "q_Z": "1", "U_q_Z": "2",
            "U_squared": "24", "U_U_Z": "8"},
        assumptions=dict(abs_Z_le_one=True, d_le_one=True, L_lower="1/2",
                         delta_upper="1/2", endpoint_caps_added_after_half_S=True,
                         Rc_lower="1"),
    )


def _independent_rc_unit_log(c, pulse_owner, power):
    """Rebuild the exact cancelled Rc^p unit from its source factors."""
    end = pulse_owner.end
    logE = current.pulse.rebox(c, pulse_owner.coefficients.selected.amplitude.log_end_scale)
    return (-end.logP - end.logUp - logE + 13 - c.ln(pulse_owner.reference)
            + power * (2 + 60 * c.ln(pulse_owner.mu))
            + (c.mpf("6.5") - 13 * power) / pulse_owner.mu)


def _independent_active_caps(owner, cell):
    """Recompute active reduced-S and endpoint caps without preRc_coefficients."""
    c = owner.c
    bd = current.rg.Bounds(c)
    ends = tuple(cell["ends"])
    views = [item["view"] for item in owner.data["actual_chart_views"]
             if tuple(item["ends"]) == ends]
    assert len(views) == len(owner.fields.charts) == 17

    def maximum_field(name, row):
        return current.maximum(c, [current.log_cap(c, view[name][row]) for view in views])

    U = maximum_field("whole_cell_physical_U_Z3_magnitude", 0)
    UZ = maximum_field("whole_cell_physical_U_Z3_magnitude", 1)
    q = maximum_field("whole_cell_physical_q_Z2_magnitude", 0)
    qZ = maximum_field("whole_cell_physical_q_Z2_magnitude", 1)
    common = bd.sum(qZ, bd.scale(U * qZ, 2), bd.scale(U * UZ, 8))
    S0 = bd.sum(bd.scale(q, 4), q * q, bd.scale(U * U, 20), common)
    S1 = bd.sum(bd.scale(q, 6), bd.scale(q * q, c.mpf("1.5")),
                bd.scale(U * U, 24), common)
    endpoint_records = cell["endpoints"]
    lower = endpoint_records["lower"]["records"]
    upper = endpoint_records["upper"]["records"]
    visc0 = bd.sum(current.log_cap(c, lower["q_plus_q_y"]),
                   current.log_cap(c, upper["q_plus_q_y"]))
    visc1 = bd.sum(current.log_cap(c, lower["q_y"]),
                   current.log_cap(c, upper["q_y"]))
    m13 = bd.sum(bd.scale(S0, c.mpf(".5")), visc0)
    m15 = bd.sum(bd.scale(S1, c.mpf(".5")), visc1)
    return dict(fields=dict(U=U, U_Z=UZ, q=q, q_Z=qZ),
                reduced_S=dict(unweighted=S0, weighted=S1),
                viscous=dict(m13=visc0, m15=visc1), m13=m13, m15=m15,
                endpoint_caps_nonzero=(visc0.log is not None and visc1.log is not None))


def _assert_active_caps_and_threshold(owner):
    """Check all normalized source logs and their independent exp(-1000) budget."""
    coefficients = _independent_majorant_coefficients()
    c = owner.c
    threshold = c.mpf(-1000)
    band_raw, _ = current.checked_archive(current.band, current.band.GATE, owner.family)
    quiet_raw, _ = current.checked_archive(current.quiet, current.quiet.GATE, owner.family)
    band_views = {tuple(item["view"]["exact_Z_cell"]): item["view"]
                  for item in band_raw["samples"].values()}
    quiet_view = quiet_raw["actual_quiet_source"]
    normalized_count = 0
    cells = {}
    def box_text(value):
        lo, hi = _interval_endpoints(value)
        return {"lower": str(lo), "upper": str(hi)}
    for cell in owner.data["cells"]:
        ends = tuple(cell["ends"])
        active = _independent_active_caps(owner, cell)
        assert active["endpoint_caps_nonzero"]
        for key, power in (("m13", 1), ("m15", 2)):
            saved_active = current.log_cap(c, cell["active"][key])
            assert _interval_bits(saved_active.log) == _interval_bits(active[key].log)
            unit = _independent_rc_unit_log(c, owner.pulse, power)
            assert _interval_bits(unit) == _interval_bits(
                current.Rc_unit_log(c, owner.pulse, power))
            band_field = band_views[ends]["signed_pressure_" + key]
            # Compact-band signed records always expose the physical-radius
            # coefficient under log_absolute_coefficient.  The quiet source
            # retains the m13 log_coefficient name because its m13 record is
            # normalized before the shared Rc factor is restored.
            band_log = current.interval(c, band_field["log_absolute_coefficient"])
            qkey = "log_coefficient" if key == "m13" else "log_absolute_coefficient"
            quiet_log = current.interval(c, quiet_view[key][qkey])
            expected = dict(
                active=active[key].log + unit,
                compact_band=band_log + unit,
                quiet=quiet_log + unit)
            saved = cell["fields"][key]["normalized_log_caps"]
            assert set(saved) == set(expected)
            for name, value in expected.items():
                actual = current.interval(c, saved[name])
                assert _interval_bits(actual) == _interval_bits(value)
                assert _upper(value) < threshold
                normalized_count += 1
        cells[str(ends)] = dict(endpoint_caps_added_after_half_S=True,
                                 normalized_upper_below_minus_1000=True,
                                 m13_active_log=box_text(active["m13"].log),
                                 m15_active_log=box_text(active["m15"].log))
    assert _interval_endpoints(owner.pulse.end.logRc)[0] > 0
    return dict(independent_reduced_S_coefficients=coefficients,
                exact_cancelled_Rc_unit_formula=True, cells=cells,
                normalized_log_bounds=normalized_count,
                normalized_upper_below_minus_1000=True,
                exact_nonzero_error_cap="3*exp(-1000)",
                Rc_ge_one=True, finite_endpoint_caps_after_half_S=True)


def _assert_geometry_and_endpoints(owner):
    geometry = owner.data["geometry"]
    assert geometry["native_jacobian_positive_on_each_source_cell"] is True
    assert len(geometry["exact_symbolic_contiguous_offset_identities"]) == 16
    assert geometry["exact_symbolic_geometry_jacobians"]["actual_patch"]["formula"] == "1/x"
    assert geometry["exact_symbolic_geometry_jacobians"]["O2_axial"]["formula"] == "40*exp(40*x)"
    g = geometry["geometry"]
    assert g["exact_logRc_recipe_match"] is True
    assert g["exact_logRp_recipe_match"] is True
    assert g["Tw_gt_two_plus_log_two"] is True
    assert g["source_logRc_lower_positive"] is True
    saved = {tuple(item["ends"]): item["endpoints"] for item in owner.data["cells"]}
    assert set(saved) == set(current.pulse.CELLS)
    for ends in current.pulse.CELLS:
        live = endpoint.endpoint_caps(owner.fields, ends)
        assert current.record(live) == saved[ends]
        assert live["lower"]["proof"]["flat_driver"]["source_branch"] == \
            "source_owned_flat_singleton"
        assert live["upper"]["proof"]["q_y_equals_minus_q"] is True
        _assert_log_record(live["lower"]["records"])
        _assert_log_record(live["upper"]["records"])
    return dict(exact_geometry=True, exact_endpoint_cells=4,
                endpoint_log_bounds_checked=True, lower_singleton_flat_proof=True,
                upper_qy_minus_q_proof=True)


def _assert_archives(owner):
    band_raw, band_hashes = current.checked_archive(current.band, current.band.GATE, owner.family)
    quiet_raw, quiet_hashes = current.checked_archive(current.quiet, current.quiet.GATE, owner.family)
    support_raw, support_hashes = current.checked_archive(
        current.support, current.support.GATE, owner.family)
    for name, digest in {**band_hashes, **quiet_hashes, **support_hashes}.items():
        assert owner.hashes[name] == digest
    assert band_raw["source_family"] == owner.family
    assert quiet_raw["source_family"] == owner.family
    assert quiet_raw["actual_quiet_source"]["quiet_segment"] == "actual_2Rc_to_Rp"
    assert support_raw["source_family"] == owner.family
    assert set(support_raw["future_charts"]) == set(current.support.FUTURE)
    assert support_raw["endpoint"]["definition"] == \
        "Rv=exp(original_logRp+13/original_mu)"
    assert support_raw["endpoint"]["angular_and_pressure_not_zeroed"] is True
    for view in support_raw["future_charts"].values():
        assert view["exact_zero_unweighted_and_R_weighted_forcing_integrals"] == [0, 0]
    return dict(compact_band_archive_checked=True, quiet_archive_checked=True,
                beyond_Rv_support_archive_checked=True,
                beyond_Rv_exact_zero_integrals=True,
                beyond_Rv_meridional_only_scope=True,
                archive_source_family=True)


def _assert_models(owner, models):
    assert len(models) == 256
    pulse_raw = json.loads(gzip.decompress((current.HERE / current.pulse.NAME).read_bytes()))
    original_models = pulse_raw["models"]
    assert len(original_models) == 256
    previous = None
    count = 0
    unchanged = 0
    for index, model in enumerate(models):
        original = original_models[index]
        left = owner.c.mpf(-1) + owner.c.mpf(2) * index / 256
        right = owner.c.mpf(-1) + owner.c.mpf(2) * (index + 1) / 256
        assert _interval_bits(model["left"]) == left._mpi_
        assert _interval_bits(model["right"]) == right._mpi_
        for field_name in ("left", "right", "center", "half_width", "approximation",
                           "source_coefficient_cell"):
            assert model[field_name] == original[field_name]
        if previous is not None:
            assert _interval_bits(model["left"]) == _interval_bits(previous)
        previous = model["right"]
        assert model["radial_scope"] == "r_minus to infinity"
        assert tuple(model["source_coefficient_cell"]) in current.pulse.CELLS
        for key, field in model["fields"].items():
            assert key in ("m13", "m15")
            original_field = original["fields"][key]
            for field_name in ("center_intercept", "center_slope",
                               "linear_Z2_half_range", "linear_Taylor_remainder",
                               "coefficient_ambiguity"):
                assert field[field_name] == original_field[field_name]
            pre = field["preRp_nonzero_source_error"]
            full = field["full_nonzero_source_error"]
            absolute = field["absolute_error"]
            assert _interval_bits(pre) == _interval_bits(owner.preRp_error[key])
            expected_full = current.interval(
                owner.c, original_field["full_nonzero_source_error"]
            ) + owner.preRp_error[key]
            expected_absolute = current.interval(
                owner.c, original_field["absolute_error"]
            ) + owner.preRp_error[key]
            assert _interval_bits(full) == _interval_bits(expected_full)
            assert _interval_bits(absolute) == _interval_bits(expected_absolute)
            assert _absolute_upper(pre) > 0 and _upper(full) > 0 and _upper(absolute) > 0
            threshold = owner.c.mpf(".001")
            assert _upper(absolute) < _upper(threshold)
            assert field["absolute_accuracy_certified"] is True
            count += 1
            unchanged += 1
    assert _interval_bits(previous) == owner.c.mpf(1)._mpi_
    return dict(models=256, model_fields=count, adjacency=True,
                full_Z_cover=True, positive_nonzero_preRp_errors=True,
                baseline_model_fields_unchanged=(unchanged == 512),
                error_provenance_checked=True, absolute_error_upper_below='0.001')


def _assert_scopes(owner, raw, samples):
    scope = owner.data["scope"]
    assert scope["lower"] == "original r_minus"
    assert scope["upper"] == "original Rp"
    assert scope["active"] == "17 charts r_minus to Rc"
    assert scope["compact"] == "Rc to 2Rc"
    assert scope["quiet"] == "2Rc to Rp"
    assert scope["original_shared_Rc_cancelled_before_boxing"] is True
    assert scope["bounded_preRp_is_nonzero_remainder_not_signed_center"] is True
    for view in samples.values():
        assert view["radial_scope"] == "r_minus to infinity"
        assert view["lower_Rin_to_rminus_not_included"] is True
        assert view["preRp_bound_is_not_signed_zero"] is True
        assert view["pressure_derivative_accuracy_not_claimed"] is True
        assert view[current.GATE] is False
        assert all(not view[name] for name in current.OPEN)
        for value in view["preRp_nonzero_error"].values():
            assert _absolute_upper(value) > 0
    for cell in owner.data["cells"]:
        for key in ("m13", "m15"):
            row = cell["fields"][key]
            assert _interval_endpoints(row["each_nonzero_absolute_upper"]) == \
                _interval_endpoints(owner.c.exp(-1000))
            assert _interval_endpoints(row["total_nonzero_absolute_upper"]) == \
                _interval_endpoints(owner.c.exp(-1000) * 3)
            assert row["no_signed_zero_assigned"] is True
            assert row["all_pieces_counted_once"] is True
            for value in row["normalized_log_caps"].values():
                _assert_log_record(value)
    return dict(scope_and_open_gates=True, nonzero_errors=True,
                omitted_lower_Rin_segment_explicit=True)


def _public_default_precision(owner):
    before = mp.mp.dps
    with mp.workdps(15):
        public = current.CurrentOriginalN1PreRpPressureSource(
            owner.fields, owner.pulse, views=owner.views, require_checked=False)
        packet = public.evaluate(".371")
        view = public.report(packet)
        assert mp.mp.dps == 15
        assert view["source_family"] == owner.family
        assert view["radial_scope"] == "r_minus to infinity"
        assert view["preRp_bound_is_not_signed_zero"] is True
    assert mp.mp.dps == before
    return True


@current.source_precision
def run(owner=None, raw=None, samples=None, models=None):
    began = time.monotonic()
    if owner is None:
        raise ValueError("Hot accepted pre-Rp owner required; checker does not rebuild ancestors")
    # Always validate the persisted artifact itself.  A hot producer may hand
    # us a live raw object whose tuple endpoints have not yet passed through
    # JSON (the saved artifact necessarily contains lists).  Canonicalize that
    # object only for the equality check, then use the freshly read artifact
    # for every saved/live assertion below.
    passed_raw = raw
    saved_raw = _read_raw()
    if passed_raw is not None:
        assert saved_raw == current.record(passed_raw)
    raw = saved_raw
    if samples is None:
        # The producer's pure view builder is intentionally used once per
        # batch.  Repeated evaluate/report calls would rerun the complete
        # parent graph fourteen times for the seven source samples.
        owner.assert_graph()
        state_before = owner._state
        samples = {z: owner._make_view(z) for z in
                   ("-1", "-.639", "0", ".211", ".371", ".639", "1")}
        owner.assert_graph()
        assert owner._state == state_before
    if models is None:
        models = owner.affine_models()
    owner.assert_graph()
    _assert_source_hashes(raw, owner)
    _assert_raw_live(raw, owner, samples, models)
    fixture = independent_pressure_fixture()
    chart_evidence = _assert_chart_views(owner)
    endpoint_evidence = _assert_geometry_and_endpoints(owner)
    archive_evidence = _assert_archives(owner)
    model_evidence = _assert_models(owner, models)
    active_cap_evidence = _assert_active_caps_and_threshold(owner)
    scope_evidence = _assert_scopes(owner, raw, samples)
    public = _public_default_precision(owner)

    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(current.__file__).name] = current.sha(Path(current.__file__).name)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(
            independent_reduced_pressure_and_unit_fixture=fixture,
            all_17_chart_caps_and_57_source_cells=chart_evidence,
            four_endpoint_caps_and_geometry=endpoint_evidence,
            accepted_compact_and_quiet_archives=archive_evidence,
            saved_live_raw_and_source_bindings=True,
            all_256_preRp_models=model_evidence,
            independent_active_caps_and_threshold=active_cap_evidence,
            scope_and_nonzero_error_guards=scope_evidence,
            default_precision_public_query=public,
            full_tails_repair_derivative_accuracy_and_global_NS_not_claimed=True),
        input_hashes=hashes,
        **{current.GATE: True},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(receipt, separators=(",", ":"), ensure_ascii=False) + "\n",
        encoding="utf8", newline="\n")
    print("PASS_N1_PRE_RP_PRESSURE_SOURCE", len(models), len(samples), flush=True)
    return receipt


if __name__ == "__main__":
    run()
