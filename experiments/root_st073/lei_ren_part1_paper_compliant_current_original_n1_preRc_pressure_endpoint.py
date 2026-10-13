"""Endpoint source caps for the original pre-:math:`R_c` pressure piece.

The accepted all-chart physical owner exports whole-cell magnitude ranges, so
it deliberately does not expose endpoint values.  This adapter obtains the
two endpoint pieces needed by the reduced pressure identity from the same
source graph:

* at ``('sc_half',)`` the original finite-``N`` correction has zero incoming
  history and zero local flat-collar driver;
* at the right edge of ``O3_power`` the actual leading velocity and its local
  derivatives, together with the finite-``N`` velocity driver, are exactly
  zero, while the own mean memory is retained in the accepted physical
  whole-chart ``q`` cap.

All returned magnitudes remain ``LogUpper`` objects.  No radius, ``N``, or
large exponential is materialized.  This file does not install a receipt or
claim a complete pressure tail.
"""

import copy
import gzip
import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sy

from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_partial_mean as all_chart
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection as geometry_source


HERE = physical.HERE
PREFIX = physical.PREFIX
sha = physical.sha
ep = physical.ep
rg = physical.rg
current = physical.current
source = physical.source
source_precision = physical.source_precision
fingerprint = physical.fingerprint

NAME = PREFIX + "current_original_n1_preRc_pressure_endpoint.json"
GATE = "current_original_n1_preRc_pressure_endpoint_caps_installed"
OPEN = (
    "current_finite_N_whole_Z_correction_prefix_through_Rc",
    "whole_Z_five_moment_terminal_repair_closed",
    "global_completed_tensor_admissibility",
    "exact_heat_exterior_installed",
    "temporal_recursion",
    "full_corrected_Navier_Stokes_solution",
)

BRIDGE_RAW = PREFIX + "current_original_whole_Z_bridge_signed_functions.json.gz"
CONTROL_RAW = PREFIX + "current_original_whole_Z_all_N_control_connection.json.gz"
FRAME_RAW = PREFIX + "current_original_C3_Rp_native_frame.json.gz"
FRAME_RECEIPT = PREFIX + "current_original_C3_Rp_native_frame_check.json"
PARAMETERS = PREFIX + "compliant_pressure_source.json"


def _as_tuple(value):
    return tuple(_as_tuple(item) for item in value) if isinstance(value, list) else value


def _read_json(name):
    path = HERE / name
    if name.endswith(".gz"):
        return json.loads(gzip.open(path, "rt", encoding="utf8").read())
    return json.loads(path.read_bytes())


def _record(value):
    if hasattr(value, "record"):
        return value.record()
    return rg.record(value)


def _exact_zero_record(value):
    if isinstance(value, dict):
        if value.get("exact_zero", False):
            return True
        if "log_absolute_upper" in value:
            return value["log_absolute_upper"] is None
    return bool(getattr(value, "zero", False)) or getattr(value, "log", object()) is None


def _log_cap(c, value, label):
    """Decode a source ``record`` into a live directed ``LogUpper``."""
    rec = _record(value)
    if not isinstance(rec, dict) or "log_absolute_upper" not in rec:
        raise ValueError("Source logarithmic record missing for " + label)
    if rec.get("exact_zero", False) or rec["log_absolute_upper"] is None:
        return rg.Bounds(c).zero
    # ``packets.interval`` preserves the exact directed tuple in the saved
    # source record.  Do not use the coefficient interval as a substitute:
    # source scales (including Pstar) live in log_absolute_upper.
    return rg.LogUpper(c, current.packets.interval(c, rec["log_absolute_upper"]))


def _cap_records(values):
    return {name: value.record() for name, value in values.items()}


def _rows(value):
    return tuple(getattr(value, name) for name in ("value", "Z", "ZZ", "ZZZ"))


def _is_zero_bound(value):
    return bool(getattr(value, "zero", False)) or getattr(value, "log", object()) is None


def _source_jet(c, values, scale, label):
    if len(values) < 4:
        raise ValueError("Four ordinary source rows required for " + label)
    caps = tuple(rg.Bounds(c).constant(math.factorial(k)) * scale *
                 _log_cap(c, value, "%s[%d]" % (label, k))
                 for k, value in enumerate(values[:4]))
    return rg.JetBound(*caps)


def _linear_cap(c, ends, delta, U, Uy, UZ):
    """Cap the value of q+q_y without introducing m or m_Z."""
    z = IntervalTaylor.variable(c, c.mpf(ends), 2)
    L = 1 - z * z * delta
    d = 1 - z * z
    coefficients = (z * (1 + delta) / L, 2 * z / L, -d / L)
    bd = rg.Bounds(c)
    values = (U.value, Uy.value, UZ.value)
    pieces = []
    for coefficient, value in zip(coefficients, values):
        pieces.append(bd.constant(ep(abs(coefficient[0]))[1]) * value)
    return bd.sum(*pieces)


def _flat_driver_proof(fields, ends, chart):
    """Check the accepted source-owned flat branch without differentiating caps."""
    if chart == "first_micro":
        # The first partition cell reaches well beyond the source-owned flat
        # collar.  Query the genuine singleton at the left endpoint instead
        # of treating that whole cell as flat.
        packet = fields.modified.target.phase.leading_source(
            tuple(ends), chart, ("sc_half",), ("sc_half",))
        if packet.get("source_identity") != fields.family:
            raise ValueError("Flat singleton source family changed")
        if (_as_tuple(packet.get("exact_left")) != ("sc_half",) or
                _as_tuple(packet.get("exact_right")) != ("sc_half",)):
            raise ValueError("First-micro flat query was not the exact singleton")
        if packet.get("actual_phase_Z_exact_zero") is not True:
            raise ValueError("Flat singleton phase must be Z independent")
        if packet.get("source_owned_modification_support") != "unchanged_left":
            raise ValueError("First-micro singleton is not the unchanged-left source branch")
        if not packet.get("source_owned_collar_q_C0_Z_exact_zero", False):
            raise ValueError("Source-owned flat singleton q theorem missing")
        initial = packet.get("actual_source_owned_initial_condition", {})
        if not initial.get("q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar", False):
            raise ValueError("Source-owned two-sided singleton flat theorem missing")
        qrows = packet.get("original_q_C0_Z")
        if not isinstance(qrows, dict) or not qrows or not all(_is_zero_bound(row) for row in qrows.values()):
            raise ValueError("Flat singleton source q rows are not exact zero")
        return dict(
            chart=chart,
            exact_source_cell=[["sc_half"], ["sc_half"]],
            source_branch="source_owned_flat_singleton",
            source_owned_modification_support=packet["source_owned_modification_support"],
            source_owned_collar_flat=True,
            source_owned_initial_condition=copy.deepcopy(initial),
            all_local_B_rows_exact_zero=True,
            singleton_q_rows_exact_zero=True,
            caps_are_not_differentiated=True,
        )
    partition = tuple(fields.partitions[chart])
    if len(partition) < 2:
        raise ValueError("Original chart partition must contain a source cell")
    cell = fields.modified.ranges.cell(tuple(ends), chart, partition[0], partition[1])
    primitive = cell["primitive"]
    B = primitive["B"]
    bzero = all(_is_zero_bound(row) for row in rg.rows(B))
    proof = primitive.get("proof", {})
    branch = proof.get("branch")
    if not bzero or branch != "exact_original_flat":
        raise ValueError("Accepted exact flat C3 driver required at " + chart)
    return dict(
        chart=chart,
        exact_source_cell=[partition[0], partition[1]],
        source_branch=branch,
        source_owned_collar_flat=bool(proof.get("source_owned_collar_flat", False)),
        all_local_B_rows_exact_zero=bzero,
        local_B_rows=rg.record(B),
        caps_are_not_differentiated=True,
    )


def _zero_first_incoming_proof(fields, ends):
    transport, stored = fields.parent._transport(tuple(ends), "first_micro")
    rows = stored["actual_N_scaled_incoming_C3"]["m"]
    if len(rows) != 4 or not all(_exact_zero_record(row) for row in rows):
        raise ValueError("First original micro incoming own-mean history is not exact zero")
    return dict(
        exact_Z_cell=list(ends),
        saved_transport_source_identity=transport["source_identity"],
        all_four_incoming_N_scaled_mean_rows_exact_zero=True,
        incoming_rows=[copy.deepcopy(row) for row in rows],
        predecessor_history_retained=True,
    )


def _first_micro_caps(fields, ends):
    """Build physical lower-end caps from the actual source rows."""
    ends = tuple(ends)
    fields.assert_graph()
    source_view = physical.source_rows(fields, ends, "first_micro", ("sc_half",))
    packet = source_view["source_packet"]
    generic = packet.get("original_generic_source")
    background = packet.get("actual_background_source")
    if not isinstance(generic, dict) or not isinstance(background, dict):
        raise ValueError("Actual first-micro generic/background source packet required")
    histories = generic.get("common_own_five_histories_axial5")
    V = generic.get("common_velocity_V_axial5")
    Vy = background.get("V_y")
    if not isinstance(histories, dict) or not isinstance(V, (list, tuple)):
        raise ValueError("Actual first-micro V/m source rows missing")
    M = histories.get("m")
    if not isinstance(M, (list, tuple)) or not isinstance(Vy, (list, tuple)):
        raise ValueError("Actual first-micro m/V_y source rows missing")
    if not packet.get("source_owned_collar_q_C0_Z_exact_zero", False):
        raise ValueError("Source-owned two-sided left-collar zero theorem missing")

    c = fields.c
    bd = rg.Bounds(c)
    scale = _log_cap(c, source_view["Pstar"], "Pstar")
    U = _source_jet(c, V, scale, "first_micro V")
    m = _source_jet(c, M, scale, "first_micro m")
    Uy = rg.JetBound(scale * _log_cap(c, Vy[0], "first_micro V_y[0]"),
                     bd.zero, bd.zero, bd.zero,)
    # The ordinary source row 1 is the first Z derivative.  This is a source
    # row, not a derivative of a magnitude cap.
    UZ = rg.JetBound(U.Z, U.ZZ, U.ZZZ, bd.zero)
    q_rows = all_chart.partial.coefficient_Q_bounds(c, ends, fields.delta, U, m)
    q = q_rows[0]
    qZ = q_rows[1]
    qplus = _linear_cap(c, ends, fields.delta, U, Uy, UZ)
    qy = bd.sum(qplus, q)
    incoming = _zero_first_incoming_proof(fields, ends)
    flat = _flat_driver_proof(fields, ends, "first_micro")
    proof = dict(
        source_family=copy.deepcopy(fields.family),
        exact_query=("sc_half",),
        original_source_query=True,
        source_owned_collar_q_C0_Z_exact_zero=True,
        finite_N_initial_history_zero=True,
        finite_N_local_driver_exact_zero=True,
        finite_N_mean_and_delta_V_at_inlet_exact_zero=True,
        physical_common_Pstar_applied_once=True,
        source_quantity_paths=dict(
            V="original_generic_source.common_velocity_V_axial5",
            m='original_generic_source.common_own_five_histories_axial5["m"]',
            V_y="actual_background_source.V_y",
            Pstar="outer.owner.owner(ends).Pstar",
        ),
        incoming_history=incoming,
        flat_driver=flat,
        original_phase_Z_exact_zero=packet.get("actual_phase_Z_exact_zero") is True,
        no_radius_or_N_materialized=True,
    )
    if not proof["original_phase_Z_exact_zero"]:
        raise ValueError("Original phase must be Z independent at lower endpoint")
    caps = dict(U=U.value, U_y=Uy.value, U_Z=UZ.value, m=m.value,
                q=q, q_Z=qZ, q_plus_q_y=qplus, q_y=qy,
                m_y=bd.sum(U.value, m.value))
    return dict(chart="first_micro", query=("sc_half",), caps=caps,
                records=_cap_records(caps), proof=proof)


def _chart_view(fields, ends, chart):
    for entry in fields._packets.values():
        state = entry[1]
        view = state.get("view") if isinstance(state, dict) else None
        if (isinstance(view, dict) and tuple(view.get("exact_Z_cell", ())) == tuple(ends)
                and view.get("chart") == chart):
            return view
    packet = fields.range_chart(tuple(ends), chart)
    return fields.report(packet)


def _upper_caps(fields, ends):
    ends = tuple(ends)
    fields.assert_graph()
    source_packet = fields.modified.target.phase.leading_source(ends, "O3_power", (2, 1))
    generic = source_packet.get("original_generic_source", {})
    V = generic.get("common_velocity_V_axial5")
    if not isinstance(V, (list, tuple)) or len(V) < 4:
        raise ValueError("Actual O3 power V source rows missing")
    if not all(_exact_zero_record(value) for value in V[:4]):
        raise ValueError("Actual O3 power leading V and local Z rows must be exact zero")
    flat = _flat_driver_proof(fields, ends, "O3_power")
    view = _chart_view(fields, ends, "O3_power")
    qrec = view["whole_cell_physical_q_Z2_magnitude"][0]
    qcap = _log_cap(fields.c, qrec, "O3_power physical q whole-chart cap")
    urows = view["whole_cell_physical_U_Z3_magnitude"]
    if not all(_exact_zero_record(row) for row in urows):
        raise ValueError("Accepted O3 physical U rows are not exact zero")
    zero = rg.Bounds(fields.c).zero
    caps = dict(U=zero, U_y=zero, U_Z=zero, q=qcap, q_Z=_log_cap(
        fields.c, view["whole_cell_physical_q_Z2_magnitude"][1],
        "O3_power physical q_Z whole-chart cap"),
        q_plus_q_y=zero, q_y=qcap)
    proof = dict(
        source_family=copy.deepcopy(fields.family), chart="O3_power",
        query=(2, 1), exact_Rc_endpoint=True,
        original_O3_power_V_rows_exact_zero=True,
        finite_N_local_velocity_driver_exact_zero=True,
        physical_U_Uy_UZ_exact_zero=True,
        physical_mean_ode="m_y=U-m=-m",
        q_y_equals_minus_q=True,
        q_plus_q_y_exact_zero=True,
        q_cap_includes_finite_N_own_memory=True,
        q_cap_scope="accepted physical whole O3_power chart, value row",
        source_quantity_paths=dict(
            V="original_generic_source.common_velocity_V_axial5",
            q="accepted all-chart physical whole_cell_physical_q_Z2_magnitude",
            Pstar="outer.owner.owner(ends).Pstar",
        ),
        flat_driver=flat,
        original_phase_Z_exact_zero=source_packet.get("actual_phase_Z_exact_zero") is True,
        no_radius_or_N_materialized=True,
    )
    if not proof["original_phase_Z_exact_zero"]:
        raise ValueError("Original O3 phase must be Z independent")
    return dict(chart="O3_power", query=(2, 1), caps=caps,
                records=_cap_records(caps), proof=proof)


def symbolic_check():
    """Independent algebra check for the endpoint cancellation identity."""
    Z, delta, U, Uy, UZ, m, mZ = sy.symbols(
        "Z delta U U_y U_Z m m_Z", real=True)
    a = 1 - delta
    d = 1 - Z**2
    L = 1 - delta * Z**2
    q = (2 * Z * U - a * Z * m - d * mZ) / L
    qy = (2 * Z * Uy - a * Z * (U - m) - d * (UZ - mZ)) / L
    qplus = (Z * (1 + delta) * U + 2 * Z * Uy - d * UZ) / L
    residual = sy.factor(sy.together(q + qy - qplus))
    if residual != 0:
        raise AssertionError("q+q_y endpoint identity failed: " + str(residual))
    upper = sy.factor(sy.together(qy.subs({U: 0, Uy: 0, UZ: 0}) + q.subs(U, 0)))
    if upper != 0:
        raise AssertionError("upper q_y=-q identity failed: " + str(upper))
    return dict(
        passed=True,
        q_plus_q_y_formula=sy.sstr(qplus),
        q_plus_q_y_residual=sy.sstr(residual),
        upper_q_y_plus_q_residual=sy.sstr(upper),
        no_q_ZZ_used=True,
    )


@source_precision
def endpoint_caps(fields, ends):
    """Return lower and upper endpoint caps for one accepted Z cell."""
    if type(fields) is not physical.CurrentOriginalN1AllChartPhysicalPressureFields:
        raise ValueError("Checked all-chart physical pressure owner required")
    if not fields.acceptance_loaded:
        raise ValueError("Accepted all-chart physical pressure owner required")
    fields.assert_graph()
    ends = tuple(ends)
    lower = _first_micro_caps(fields, ends)
    upper = _upper_caps(fields, ends)
    proof = dict(
        source_family=copy.deepcopy(fields.family), exact_Z_cell=list(ends),
        original_first_micro_query=("sc_half",), original_O3_query=(2, 1),
        symbolic=symbolic_check(),
        common_physical_Pstar_once=True,
        endpoint_caps_are_LogUpper=True,
        no_signed_centers_or_radius_materialization=True,
        q_ZZ_not_required=True,
        full_pressure_tail_repair_and_global_claims_false=True,
    )
    return dict(source_family=copy.deepcopy(fields.family), exact_Z_cell=ends,
                lower=lower, upper=upper, lower_caps=lower["caps"],
                upper_caps=upper["caps"], proof=proof)


def _node_id(value):
    return getattr(value, "node", None)


def _directed_box(value):
    lo, hi = ep(value)
    return dict(lower=str(lo), upper=str(hi))


def _endpoint_bits(value):
    lo, hi = ep(value)
    return (lo._mpf_, hi._mpf_)


def _symbolic_geometry(pulse_owner, c):
    rho, Tw = sy.symbols("rho Tw", real=True)
    rc = rho
    rw = rc - 2
    rp = rw + Tw
    residual = sy.simplify((rp - rc) - (Tw - 2))
    if residual != 0:
        raise AssertionError("Rc/Rp geometry identity failed")
    end = getattr(pulse_owner, "end", None)
    if end is None:
        raise ValueError("Pulse owner must retain its accepted end-pressure owner")
    # The all-chart field owner and the accepted end owner can carry distinct
    # directed interval contexts (the former is 540 digits while the archived
    # end source may be 500 digits).  Recipe equality must be checked in the
    # context that created logRc/logRp; changing only ambient mp.dps does not
    # rebox those intervals.
    source_c = getattr(end, "c", None)
    if source_c is None:
        raise ValueError("Accepted end owner interval context is required")
    mu = getattr(end, "mu", None)
    log_rc = getattr(end, "logRc", None)
    log_rp = getattr(end, "logRp", None)
    if mu is None or log_rc is None or log_rp is None:
        raise ValueError("Accepted pulse mu/logRc/logRp source fields are required")
    radius_binding = getattr(end, "radius_binding", None)
    coefficients_owner = getattr(end, "coefficients", None)
    basis = getattr(coefficients_owner, "basis", None)
    selected_log_c = radius_binding.get("selected_logCstar") if isinstance(radius_binding, dict) else None
    if not isinstance(selected_log_c, dict):
        raise ValueError("Accepted native radius selected_logCstar source binding is required")
    # ``basis`` is a typed owner, while ``read_interval`` is the module-level
    # decoder used by the accepted end-pressure constructor.  Decode the
    # archived exact endpoint tuple directly so this proof does not depend on
    # an instance method that the basis owner does not expose.
    try:
        log_c = source_c.make_mpf((tuple(selected_log_c["lower_exact_mpf_tuple"]),
                                   tuple(selected_log_c["upper_exact_mpf_tuple"])))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Accepted selected_logCstar exact interval is required") from exc
    mu_box = mu
    if ep(mu_box)[0] <= 0:
        raise ValueError("Accepted pulse mu must be strictly positive")
    expected_log_rc = source_c.ln(110) + 10 * log_c + 11 * source_c.exp(40) + 124
    expected_log_rp = expected_log_rc - 2 - 60 * source_c.ln(mu_box)
    if _endpoint_bits(expected_log_rc) != _endpoint_bits(log_rc):
        raise ValueError("Accepted logRc does not equal the native selected_logCstar recipe")
    if _endpoint_bits(expected_log_rp) != _endpoint_bits(log_rp):
        raise ValueError("Accepted logRp does not equal the native logRc/mu recipe")
    log_rc_box = log_rc
    tw_from_mu = -source_c.mpf(60) * source_c.ln(mu_box)
    threshold = source_c.mpf(2) + source_c.ln(source_c.mpf(2))
    if ep(tw_from_mu)[0] <= ep(threshold)[1]:
        raise ValueError("Actual pulse Tw=-60 log(mu) does not exceed 2+log(2)")
    if ep(log_rc_box)[0] <= 0:
        raise ValueError("Accepted absolute logRc lower bound does not prove Rc >= 1")
    return dict(
        passed=True,
        source_mu="pulse_owner.end.mu",
        source_mu_interval=_directed_box(mu),
        source_Tw="-60*log(mu)",
        source_Tw_interval=_directed_box(tw_from_mu),
        source_geometry_context_dps=getattr(source_c, "dps", None),
        source_logRc_recipe="log(110)+10*selected_logCstar+11*exp(40)+124",
        source_logRp_recipe="logRc-2-60*log(mu)",
        selected_logCstar_interval=_directed_box(log_c),
        selected_logCstar_exact_singleton=(ep(log_c)[0] == ep(log_c)[1]),
        source_logRc="pulse_owner.end.logRc",
        source_logRc_interval=_directed_box(log_rc),
        source_logRp="pulse_owner.end.logRp",
        source_logRp_interval=_directed_box(log_rp),
        native_radius_binding_functions=copy.deepcopy(radius_binding.get("functions", {})),
        native_radius_expected_nodes=copy.deepcopy(radius_binding.get("expected_nodes", {})),
        native_frame_receipt=FRAME_RECEIPT,
        native_frame_receipt_sha256=sha(FRAME_RECEIPT),
        exact_logRc_recipe_match=True,
        exact_logRp_recipe_match=True,
        source_logRc_lower_positive=True,
        Tw_threshold="2+log(2)",
        Tw_threshold_interval=_directed_box(threshold),
        Rc_to_Rp_log_difference="Tw-2",
        Rc_to_Rp_ratio="exp(Tw-2)>2",
        Tw_gt_two_plus_log_two=True,
        exact_residual=sy.sstr(residual),
    )


@source_precision
def geometry_proof(fields, pulse_owner):
    """Validate the accepted source-bound 17-chart geometry and Rc/Rp order."""
    if type(fields) is not physical.CurrentOriginalN1AllChartPhysicalPressureFields:
        raise ValueError("Checked all-chart physical pressure owner required")
    if not fields.acceptance_loaded:
        raise ValueError("Accepted all-chart physical pressure owner required")
    fields.assert_graph()
    if pulse_owner is None:
        raise ValueError("Accepted same-family pulse owner required")
    pulse_family = getattr(pulse_owner, "family", getattr(pulse_owner, "identity", None))
    if pulse_family is not None and pulse_family != fields.family:
        raise ValueError("Pulse owner source family differs from all-chart physical source")

    control = _read_json(CONTROL_RAW)
    frame = _read_json(FRAME_RAW)
    frame_receipt = _read_json(FRAME_RECEIPT)
    if frame_receipt.get("all_passed") is not True or not frame_receipt.get(
            "current_original_Rp_native_exact_radius_identity_installed", False):
        raise ValueError("Accepted native Rc/Rp frame receipt required")
    if control.get("source_family") != fields.family or frame.get("source_family") != fields.family:
        raise ValueError("Geometry source family changed")
    windows = control.get("exact_current_integral_and_control_graph", {}).get("windows", [])
    charts = tuple(fields.charts)
    if len(charts) != 17 or tuple(w.get("chart") for w in windows) != charts:
        raise ValueError("All 17 original chart windows required in source order")
    offsets, domains, native_symbol = geometry_source.exact_geometry()
    if tuple(geometry_source.CHARTS) != charts:
        raise ValueError("Exact geometry chart order differs from all-chart owner")
    contiguous = []
    for left_chart, right_chart in zip(charts, charts[1:]):
        left_end = sy.simplify(offsets[left_chart].subs(native_symbol, domains[left_chart][1]))
        right_start = sy.simplify(offsets[right_chart].subs(native_symbol, domains[right_chart][0]))
        identity = sy.simplify(left_end - right_start)
        if identity != 0:
            raise ValueError("Exact symbolic chart-offset seam is not contiguous")
        contiguous.append(dict(left_chart=left_chart, right_chart=right_chart,
                               residual=sy.sstr(identity)))
    rows = []
    jacobians = {}

    def coordinate(value):
        return geometry_source.symbolic_coordinate(_as_tuple(value))

    for window, chart in zip(windows, charts):
        partition = tuple(fields.partitions[chart])
        source_partition = tuple(geometry_source.native_partitions()[chart])
        if tuple(_as_tuple(item) for item in source_partition) != tuple(_as_tuple(item) for item in partition):
            raise ValueError("Exact source partition changed for " + chart)
        domain_left, domain_right = domains[chart]
        if sy.simplify(coordinate(partition[0]) - domain_left) != 0:
            raise ValueError("Source partition left endpoint changed for " + chart)
        if sy.simplify(coordinate(partition[-1]) - domain_right) != 0:
            raise ValueError("Source partition right endpoint changed for " + chart)
        cells = window.get("actual_current_source_cells", [])
        if len(cells) != len(partition) - 1:
            raise ValueError("Original contiguous source partition changed for " + chart)
        for cell, left, right in zip(cells, partition, partition[1:]):
            if (_as_tuple(cell.get("exact_left")) != left or
                    _as_tuple(cell.get("exact_right")) != right):
                raise ValueError("Source cell endpoint changed for " + chart)
        jacobian = sy.simplify(sy.diff(offsets[chart], native_symbol))
        if chart == "actual_patch":
            if sy.simplify(jacobian - 1 / native_symbol) != 0:
                raise ValueError("Actual-patch Jacobian source form changed")
            jacobian_kind = "1/x > 0 on x in [1, exp(1)]"
        elif chart == "O2_axial":
            if sy.simplify(jacobian - 40 * sy.exp(40 * native_symbol)) != 0:
                raise ValueError("O2-axial Jacobian source form changed")
            jacobian_kind = "40*exp(40*x) > 0 on the original domain"
        else:
            if sy.simplify(sy.diff(jacobian, native_symbol)) != 0:
                raise ValueError("Ordinary source Jacobian is not constant for " + chart)
            jacobian_kind = "constant exact source slope; positivity from source mass/native width"
        jacobians[chart] = dict(formula=sy.sstr(jacobian), kind=jacobian_kind,
                                constant=chart not in ("actual_patch", "O2_axial"))
        cell_rows = []
        for cell, left, right in zip(cells, partition, partition[1:]):
            width = sy.simplify(coordinate(right) - coordinate(left))
            if width == 0:
                raise ValueError("Exact source cell has zero native width")
            # weights() enforces the original ordered native interval and
            # supplies the directed positive mass for that exact cell.
            weight = fields.modified.ranges.weights(("0", ".5"), chart, left, right, Fraction(0))
            mass = weight["original_mass"]
            mass_rec = _record(mass)
            if mass_rec.get("exact_zero", False):
                raise ValueError("Accepted native Jacobian/source mass is not positive for " + chart)
            coefficient = getattr(mass, "coefficient", None)
            if coefficient is None:
                coefficient = current.packets.interval(fields.c, mass_rec["coefficient_interval"])
            lo, hi = ep(coefficient)
            if lo <= 0:
                raise ValueError("Accepted native Jacobian/source mass is not positive for " + chart)
            cell_rows.append(dict(left=left, right=right,
                                  native_width_formula=sy.sstr(width),
                                  positive_native_width_by_source_order=True,
                                  positive_original_mass_record=mass_rec))
        rows.append(dict(chart=chart,
                         offset_node=window.get("offset"),
                         left_offset_node=window.get("left_offset"),
                         right_offset_node=window.get("right_offset"),
                         Jacobian_node=window.get("Jacobian"),
                         Jacobian_source=jacobians[chart],
                         cells=cell_rows,
                         source_partition_cell_count=len(cells)))

    funcs = frame.get("actual_Rp_native_frame_functions", {})
    if funcs.get("original_Rc_phase_offset") != rows[-1]["right_offset_node"]:
        raise ValueError("O3 power right offset is not the accepted Rc phase node")
    if funcs.get("original_logRc") is None or funcs.get("original_logRp") is None:
        raise ValueError("Accepted absolute Rc/Rp graph nodes missing")
    geometry = _symbolic_geometry(pulse_owner, fields.c)
    return dict(
        source_family=copy.deepcopy(fields.family), pulse_owner_type=type(pulse_owner).__name__,
        exact_17_chart_source_windows=rows,
        source_control_graph_sha256=sha(CONTROL_RAW),
        accepted_native_frame_sha256=sha(FRAME_RAW),
        accepted_native_frame_receipt_sha256=sha(FRAME_RECEIPT),
        original_Rc_phase_offset_node=funcs["original_Rc_phase_offset"],
        original_logRc_node=funcs["original_logRc"],
        original_logRp_node=funcs["original_logRp"],
        original_log_r_minus_formula="log(4)-4*logP-1000+hbB*sc/2",
        absolute_Rc_formula="original_O3_power_right_offset+log_r_minus",
        absolute_Rw_formula="logRc-2",
        absolute_Rp_formula="logRw+Tw",
        exact_symbolic_contiguous_offset_identities=contiguous,
        exact_symbolic_geometry_jacobians=jacobians,
        native_jacobian_positive_on_each_source_cell=True,
        original_contiguous_17_chart_offsets=True,
        geometry=geometry,
        source_radius_and_phase_Z_independent=True,
        radius_caps_not_subtracted=True,
    )


__all__ = [
    "NAME", "GATE", "OPEN", "endpoint_caps", "geometry_proof", "symbolic_check",
]


if __name__ == "__main__":
    print(json.dumps(symbolic_check(), sort_keys=True))
