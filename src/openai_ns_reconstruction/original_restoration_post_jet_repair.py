"""Repair the lost post-restoration ordinary ``Z`` jet.

The accepted restoration background packet has the incoming post-restoration
memory through ordinary ``Z`` order 3, but the old post chart arithmetic
subtracted only three rows before padding with zeros.  The zero padding hides
the saved, nonzero enclosure of ``M_restore_Z3``.  This module consumes the
accepted serialized packet and adds the missing term analytically.  It does
not call the accepted builder, quadrature, norm, owner, or checker.

Only the post-restoration rows which depend on ``M_Z3`` are extended.  The
extension is an interval error contribution, rather than a selected signed
value or a midpoint.  The old packet is retained in the returned wrapper so
that downstream finite-``N`` work can opt into the corrected numeric view
without changing the accepted ancestor payload.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
PARENT_SHA256 = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
MODULE = "src/openai_ns_reconstruction/original_restoration_post_jet_repair.py"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
OLD_BACKGROUND_SOURCE = "src/openai_ns_reconstruction/original_restoration_background_I1_numeric.py"
OLD_BACKGROUND_SHA256 = "aa03b74005d2fff69518baf76298412d4c2ccb668a4d6de9d450d6c3caaff0ee"
SCRATCH = "work/restoration_post_jet_repair.unaccepted_packet.json"

ORDERS = (0, 1, 2, 3)

D_FORMULA = (
    "D=-I[q]/Lz+h*Z*I[q_Z]/(2Lz)-I[q_squared]/2+"
    "d*I[U_q_Z]/Lz+2*h*Z^2*I[U_squared]/Lz^2+"
    "2*Z*d*I[U_U_Z]/Lz^2"
)
Q_FORMULA = "Q=(2*Z*U-((1-delta)*Z*M+(1-Z^2)*M_Z))/(1-delta*Z^2)"


def _lo(value):
    from mpmath import make_mpf

    return make_mpf(value._mpi_[0])


def _hi(value):
    from mpmath import make_mpf

    return make_mpf(value._mpi_[1])


def _point(c, value):
    if hasattr(value, "_mpi_"):
        low, high = _lo(value), _hi(value)
        require(low == high, "A point interval was required")
        value = low
    value = c.mpf(value)
    return c.mpf([value, value])


def _read_row(c, value):
    if isinstance(value, Mapping) and "value" in value:
        value = value["value"]
    require(isinstance(value, Mapping), "A directed interval record is required")
    require(
        "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
        "A directed exact MPF tuple is required",
    )
    require(value.get("directed") is True, "The saved interval must be directed")
    out = read(c, value)
    require(_lo(out) <= _hi(out), "The saved interval has reversed endpoints")
    return out


def _rows(c, rows, count, label):
    require(isinstance(rows, list), label + " must be a list")
    found = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), label + " row must be a mapping")
        order = row.get("ordinary_Z_order", index)
        require(type(order) is int and order in ORDERS, label + " has an invalid order")
        require(order not in found, label + " repeats an order")
        found[order] = _read_row(c, row)
    require(set(found) >= set(range(count)), label + " has missing rows")
    return [found[index] for index in range(count)]


def _replace(rows, order, value, *, kind=None):
    result = copy.deepcopy(rows)
    require(0 <= order < len(result), "Cannot replace a missing ordinary row")
    result[order]["value"] = record(value)
    if kind is not None:
        result[order]["kind"] = kind
    return result


def _expand(old, padding, label):
    """Add a rigorous signed interval padding and check containment."""

    require(_lo(padding) <= 0 <= _hi(padding), label + " must contain zero")
    new = old + padding
    require(_lo(new) <= _lo(old) and _hi(new) >= _hi(old),
            label + " failed old-interval containment")
    require(_lo(new) < _lo(old) or _hi(new) > _hi(old),
            label + " did not strictly expand")
    return new


def _source_bugline_guard(root):
    """Pin the exact C6 post-restoration source lines used by this repair."""

    path = root / C6_SOURCE
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    required = (
        ("post_restore_fraction", "p=g.symbol('true_prescribed_post_restore_fraction')"),
        ("post_restore_U", "actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)"),
        ("post_restore_chart", "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])"),
        ("post_restore_endpoint", "self.endpoints['Rm']=self.end('post_restore',g.one)"),
    )
    found = []
    compact_lines = ["".join(line.split()) for line in lines]
    for name, snippet in required:
        compact = "".join(snippet.split())
        matches = [index + 1 for index, line in enumerate(compact_lines) if compact in line]
        require(len(matches) == 1, "C6 source bug-line guard changed: " + name)
        found.append(dict(name=name, line=matches[0], required_text=snippet))
    return dict(
        passed=True,
        path=C6_SOURCE,
        expected_sha256=C6_SHA256,
        actual_sha256=sha(path),
        lines=found,
        source_executed=False,
    )


def _source_guards(root):
    require(sha(root / C6_SOURCE) == C6_SHA256, "Pinned C6 source changed")
    require(sha(root / OLD_BACKGROUND_SOURCE) == OLD_BACKGROUND_SHA256,
            "Pinned accepted restoration implementation changed")
    c6 = source_guard(root, C6_SOURCE, "__init__", (
        "self.E=E=a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))",
        "U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
        "actual=dict(phi=phi,U=U,phi_y=a.scale(phi,'-2/5'),U_y=a.scale(E,g.neg(sigma_prime(g,t)))",
        "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
        "p=g.symbol('true_prescribed_post_restore_fraction')",
        "actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)",
        "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
    ))
    old_background = source_guard(
        root, OLD_BACKGROUND_SOURCE, "build",
        ("post_extra=_jet_sub(post_incoming_M,base,count=3)+_zeros(c,2)",),
    )
    return dict(
        c6_actual_definition=c6,
        c6_source_bugline_tuple=_source_bugline_guard(root),
        accepted_old_background_implementation=old_background,
        source_executed=False,
    )


def _load_parent(root, parent):
    path = root / PARENT
    require(path.is_file(), "The accepted restoration packet is missing")
    actual_hash = sha(path)
    require(actual_hash == PARENT_SHA256, "The accepted restoration packet changed")
    disk = json.loads(path.read_bytes())
    require(disk == parent, "build requires the exact accepted restoration packet")
    numeric = disk.get("actual_restoration_background_I1_numeric")
    require(isinstance(numeric, Mapping) and numeric.get("passed") is True,
            "The accepted numeric background packet is required")
    require(numeric.get("source_executed") is False
            and numeric.get("accepted_ancestors_replayed") is False,
            "Accepted numeric background execution audit changed")
    return disk


def _directions(c, parent):
    switch = parent["actual_switch_post_background_I1_numeric"]
    first = parent["actual_first_switch_background_I1_numeric"]
    require(switch.get("passed") is True and first.get("passed") is True,
            "Accepted first-switch/post-switch directions are required")
    directions = switch["comparison_directions"]
    require(directions == first["comparison_directions"],
            "The accepted direction source diverged between charts")
    delta = _read_row(c, first["Q_b_proof"]["delta"])
    Z = _rows(c, directions["Z"], 3, "source Z direction")
    d = _rows(c, directions["d"], 3, "source d direction")
    L = _rows(c, directions["L"], 3, "source L direction")
    require(_lo(Z[0]) > 0 and _hi(Z[0]) < 1, "The source Z interval left its strict domain")
    require(_lo(d[0]) > 0 and _lo(L[0]) > 0, "Pressure denominators must be positive")
    return dict(Z=Z, d=d, L=L, delta=delta, h=_point(c, 1) - delta,
                source_formula=first["Q_b_proof"].get("formula"))


def _post_row_map(c, packet):
    post = packet["charts"]["post_restore"]
    restoration = packet["charts"]["restoration"]
    require(post.get("chart") == "post_restore", "The post chart name changed")
    require(restoration.get("chart") == "restoration", "The restoration chart name changed")
    source = post["native_source"]
    incoming = _rows(c, source["incoming_M"], 4, "post incoming M")
    restore_endpoint = _rows(c, restoration["endpoint_fields"]["M_end"], 4,
                             "restoration M endpoint")
    for index in range(4):
        require(incoming[index] == restore_endpoint[index],
                "post incoming M is not the restoration endpoint source")
    return post, restoration, incoming, restore_endpoint


def _correct_post(c, packet, directions):
    post, restoration, incoming, restore_endpoint = _post_row_map(c, packet)
    old_post = copy.deepcopy(post)
    old_signed = copy.deepcopy(packet["signed_rows"])["post_restore"]

    m_restore_3 = restore_endpoint[3]
    require(_lo(m_restore_3) < 0 < _hi(m_restore_3),
            "The saved M_restore_Z3 enclosure must remain nonzero")
    exp_minus_one = c.exp(c.mpf(-1))
    normalized_missing_primitive = exp_minus_one - c.exp(c.mpf(-2))
    m_post_3 = m_restore_3 * exp_minus_one

    # The accepted coefficient row is a directed enclosure of the actual
    # source coefficient.  Reusing that row avoids midpoint substitution.
    coefficients = post["arithmetic"]["coefficients"]
    d_over_L = _read_row(c, coefficients["d_over_L"][0])
    hZ_over_2L = _read_row(c, coefficients["hZ_over_2L"][0])
    delta_Q_start_Z2 = -(d_over_L * m_restore_3)
    delta_Q_end_Z2 = -(d_over_L * m_post_3)
    delta_I_q_Z_Z1 = delta_Q_start_Z2 * normalized_missing_primitive
    U0 = _read_row(c, post["endpoint_fields"]["U_start"][0])
    delta_I_Uq_Z_Z1 = U0 * delta_I_q_Z_Z1

    # Verify the accepted source formula before applying the D correction.
    source = post["arithmetic"]["source"]
    require(source.get("D_formula") == D_FORMULA, "The accepted D source formula changed")
    require(source.get("actual_pressure_formula") == Q_FORMULA,
            "The accepted pressure source formula changed")
    # Combine the correlated q_Z and U*q_Z coefficients before interval
    # boxing.  Multiplying and adding the two terms independently would
    # count the same missing M_Z3 interval twice.
    combined_D_coefficient = hZ_over_2L + d_over_L * U0
    delta_D_Z1 = combined_D_coefficient * delta_I_q_Z_Z1

    endpoint_fields = post["endpoint_fields"]
    endpoint_fields["M_start"] = _replace(
        endpoint_fields["M_start"], 3,
        _expand(_read_row(c, endpoint_fields["M_start"][3]), m_restore_3,
                "post M_start_Z3"),
        kind="corrected_actual_endpoint_M_Z3",
    )
    endpoint_fields["M_end"] = _replace(
        endpoint_fields["M_end"], 3,
        _expand(_read_row(c, endpoint_fields["M_end"][3]), m_post_3,
                "post M_end_Z3"),
        kind="corrected_actual_endpoint_M_Z3",
    )
    endpoint_fields["Q_start"] = _replace(
        endpoint_fields["Q_start"], 2,
        _expand(_read_row(c, endpoint_fields["Q_start"][2]), delta_Q_start_Z2,
                "post Q_start_Z2"),
        kind="corrected_actual_Q_Z2",
    )
    endpoint_fields["Q_end"] = _replace(
        endpoint_fields["Q_end"], 2,
        _expand(_read_row(c, endpoint_fields["Q_end"][2]), delta_Q_end_Z2,
                "post Q_end_Z2"),
        kind="corrected_actual_Q_Z2",
    )

    arithmetic = post["arithmetic"]
    arithmetic["Q_inlet"] = _replace(
        arithmetic["Q_inlet"], 2,
        _expand(_read_row(c, arithmetic["Q_inlet"][2]), delta_Q_start_Z2,
                "post arithmetic Q_inlet_Z2"),
        kind="corrected_actual_Q_inlet_Z2",
    )
    arithmetic["Q_endpoint"] = _replace(
        arithmetic["Q_endpoint"], 2,
        _expand(_read_row(c, arithmetic["Q_endpoint"][2]), delta_Q_end_Z2,
                "post arithmetic Q_endpoint_Z2"),
        kind="corrected_actual_Q_endpoint_Z2",
    )

    six = arithmetic["six_integrals"]
    six["q_Z"] = _replace(
        six["q_Z"], 1,
        _expand(_read_row(c, six["q_Z"][1]), delta_I_q_Z_Z1,
                "post I[q_Z]_Z1"),
        kind="corrected_actual_I1_q_Z_Z1",
    )
    six["U_q_Z"] = _replace(
        six["U_q_Z"], 1,
        _expand(_read_row(c, six["U_q_Z"][1]), delta_I_Uq_Z_Z1,
                "post I[U_q_Z]_Z1"),
        kind="corrected_actual_I1_U_q_Z_Z1",
    )

    arithmetic["D"] = _replace(
        arithmetic["D"], 1,
        _expand(_read_row(c, arithmetic["D"][1]), delta_D_Z1,
                "post D_Z1"),
        kind="corrected_actual_D_Z1",
    )
    arithmetic["S1_over_W"] = _replace(
        arithmetic["S1_over_W"], 1,
        _expand(_read_row(c, arithmetic["S1_over_W"][1]), delta_D_Z1,
                "post S1_Z1"),
        kind="corrected_actual_S1_Z1",
    )

    # Keep the top-level signed-row view synchronized with the chart
    # arithmetic.  B, q, q_squared, U_squared and U_U_Z are untouched.
    signed = packet["signed_rows"]["post_restore"]
    signed["six_I1_signed_C1"]["q_Z"] = copy.deepcopy(six["q_Z"])
    signed["six_I1_signed_C1"]["U_q_Z"] = copy.deepcopy(six["U_q_Z"])
    signed["D_signed_C1"] = copy.deepcopy(arithmetic["D"])
    signed["S1_signed_C1"] = copy.deepcopy(arithmetic["S1_over_W"])
    post["signed_rows"] = copy.deepcopy(signed)

    post["native_source"]["post_jet_repair"] = dict(
        incoming_M_rows_complete_through_Z3=True,
        incoming_M_source="restoration.endpoint_fields.M_end",
        old_post_endpoint_M_Z3_superseded=True,
        old_post_endpoint_M_Z3_kind="invalid_zero_padding",
        source_function_same_C6_cache=True,
        source_executed=False,
    )
    source["post_jet_repair_D_formula_verified"] = True
    source["post_jet_repair_Q_formula_verified"] = True
    source["source_executed"] = False

    repair = dict(
        passed=True,
        chart="post_restore",
        ordinary_Z_orders=list(ORDERS),
        affected_integrals=["q_Z", "U_q_Z"],
        affected_integral_order=1,
        affected_pressure_endpoint_order=2,
        affected_D_order=1,
        affected_S1_order=1,
        incoming_M_report_rows=4,
        m_restore_Z3=record(m_restore_3),
        m_post_Z3=record(m_post_3),
        missing_Q_Z2_function=("-(d/L)*M_restore_Z3*exp(-q)"),
        normalized_missing_primitive="exp(-1)-exp(-2)",
        normalized_missing_primitive_value=record(normalized_missing_primitive),
        delta_Q_start_Z2=record(delta_Q_start_Z2),
        delta_Q_end_Z2=record(delta_Q_end_Z2),
        delta_I_q_Z_Z1=record(delta_I_q_Z_Z1),
        delta_I_Uq_Z_Z1=record(delta_I_Uq_Z_Z1),
        delta_D_Z1=record(delta_D_Z1),
        delta_S1_Z1=record(delta_D_Z1),
        coefficient_rows={
            "hZ_over_2L_Z0": record(hZ_over_2L),
            "d_over_L_Z0": record(d_over_L),
            "U0": record(U0),
            "combined_q_Z_to_D_Z1": record(combined_D_coefficient),
        },
        d_formula=D_FORMULA,
        q_formula=Q_FORMULA,
        d_formula_source_verified=True,
        source_function_same_C6_cache=True,
        physical_measure="R^2*dq; Wpost=Rm^2",
        actual_interval_used=True,
        midpoint_used=False,
        signed_center_selected=False,
        error_interval_contains_zero=True,
        correlated_coefficient_combination=True,
        no_manufactured_Z3_or_Z4_derivatives=True,
        third_or_fourth_derivative_claim=False,
        old_post_endpoint_M_Z3_superseded=True,
        old_post_endpoint_Q_Z2_superseded=True,
        source_executed=False,
    )

    affected = (
        ("six_I1.q_Z_Z1", six["q_Z"][1]),
        ("six_I1.U_q_Z_Z1", six["U_q_Z"][1]),
        ("D_Z1", arithmetic["D"][1]),
        ("S1_Z1", arithmetic["S1_over_W"][1]),
    )
    repair["strictly_expanded_rows"] = [name for name, _ in affected]
    repair["strictly_expanded_row_count"] = len(affected)
    repair["old_packet_retained"] = old_post
    repair["old_signed_rows_retained"] = old_signed
    return repair


def _fixture_bound_check():
    """Derive a nonzero M_Z3 functional independently before boxing it."""

    import sympy as sp

    c = MPIntervalContext()
    c.prec = 180
    z, q, delta, amplitude, z0 = sp.symbols("z q delta amplitude z0")
    d = 1 - z**2
    L = 1 - delta * z**2
    U = 4 * z
    M = 4 * z + amplitude * (z - z0) ** 3 * sp.exp(-q) / 6
    Q = (2 * z * U - (1 - delta) * z * M - d * sp.diff(M, z)) / L
    Q_base = (2 * z * U - (1 - delta) * z * (4 * z)
              - d * sp.diff(4 * z, z)) / L
    delta_Q = sp.simplify(Q - Q_base)
    delta_Q_Z = sp.diff(delta_Q, z)
    delta_Q_Z2 = sp.simplify(sp.diff(delta_Q, z, 2).subs(z, z0))
    d0_symbolic = 1 - z0**2
    L0_symbolic = 1 - delta * z0**2
    expected_Q_Z2 = -d0_symbolic / L0_symbolic * amplitude * sp.exp(-q)
    require(sp.simplify(delta_Q.subs(z, z0)) == 0,
            "Synthetic Q Z0 residual was not zero")
    require(sp.simplify(delta_Q_Z.subs(z, z0)) == 0,
            "Synthetic Q_Z Z0 residual was not zero")
    require(sp.simplify(delta_Q_Z2 - expected_Q_Z2) == 0,
            "Independent synthetic Q_Z2 residual derivation failed")

    Uq_Z = U * sp.diff(Q, z)
    Uq_Z_base = (4 * z) * sp.diff(Q_base, z)
    delta_Uq_Z1 = sp.simplify(sp.diff(Uq_Z - Uq_Z_base, z).subs(z, z0))
    require(sp.simplify(delta_Uq_Z1 - 4 * z0 * expected_Q_Z2) == 0,
            "Independent U*q_Z C1 derivation failed")

    normalized_kernel = sp.integrate(
        sp.exp(2 * (q - 1)) * sp.exp(-q), (q, 0, 1))
    expected_kernel = sp.exp(-1) - sp.exp(-2)
    require(sp.simplify(normalized_kernel - expected_kernel) == 0,
            "Independent normalized missing primitive derivation failed")

    Z0 = _point(c, "37/100")
    delta0 = _point(c, "1/1000")
    d0 = _point(c, 1) - Z0 * Z0
    L0 = _point(c, 1) - delta0 * Z0 * Z0
    h0 = _point(c, 1) - delta0
    U0 = 4 * Z0
    m3 = c.mpf([c.mpf("-1/4"), c.mpf("1/3")])
    kappa = c.exp(-1) - c.exp(-2)
    expected_qz = -(d0 / L0) * m3 * kappa
    expected_uqz = U0 * expected_qz
    combined_coefficient = h0 * Z0 / (2 * L0) + (d0 / L0) * U0
    expected_D = combined_coefficient * expected_qz
    # A nonzero polynomial baseline makes the enclosure change visible in
    # both I1 rows and the derived D/S1 rows; no cap is selected as a value.
    baseline_qz = c.mpf([2, 2])
    baseline_uqz = c.mpf([3, 3])
    baseline_D = c.mpf([5, 5])
    baseline_S1 = c.mpf([7, 7])
    corrected_qz = _expand(baseline_qz, expected_qz, "fixture q_Z")
    corrected_uqz = _expand(baseline_uqz, expected_uqz, "fixture U_q_Z")
    corrected_D = _expand(baseline_D, expected_D, "fixture D")
    corrected_S1 = _expand(baseline_S1, expected_D, "fixture S1")
    require(corrected_qz != baseline_qz and corrected_uqz != baseline_uqz
            and corrected_D != baseline_D and corrected_S1 != baseline_S1,
            "Synthetic nonzero functional had no effect")
    return dict(
        passed=True,
        synthetic_M_restore_Z3=record(m3),
        synthetic_M_formula="M=4Z+amplitude*(Z-z0)^3*exp(-q)/6",
        independent_symbolic_Q_Z2="-(1-z0^2)/(1-delta*z0^2)*amplitude*exp(-q)",
        independent_symbolic_Q_Z0=True,
        independent_symbolic_Q_Z1=True,
        independent_symbolic_U_q_Z1="U(z0)*delta_Q_Z2",
        exact_primitive="exp(-1)-exp(-2)",
        residual_identity="delta_Q_Z2 + (1-z0^2)/(1-delta*z0^2)*amplitude*exp(-q) = 0",
        nonzero_functional_changes_q_Z1=True,
        nonzero_functional_changes_U_q_Z1=True,
        nonzero_functional_changes_D_Z1=True,
        nonzero_functional_changes_S1_Z1=True,
        no_midpoint_or_zero_substitution=True,
    )


def fixture_bound_check():
    """Public local test; it does not read or execute an accepted packet."""

    return _fixture_bound_check()


def build(root, parent):
    """Return a corrected numeric view over the accepted restoration packet."""

    root = Path(root)
    disk = _load_parent(root, parent)
    guards = _source_guards(root)
    c = MPIntervalContext()
    c.prec = 1600
    original = disk["actual_restoration_background_I1_numeric"]
    corrected = copy.deepcopy(original)
    directions = _directions(c, disk)
    repair = _correct_post(c, corrected, directions)
    fixture = _fixture_bound_check()

    corrected["corrected_numeric_view"] = True
    corrected["post_jet_repair"] = copy.deepcopy(repair)
    corrected["source_guards"] = dict(
        corrected.get("source_guards", {}),
        post_jet_repair=copy.deepcopy(guards),
    )
    corrected["source_executed"] = False
    corrected["accepted_ancestors_replayed"] = False
    corrected["limits"] = dict(corrected.get("limits", {}),
                                post_jet_repair_only=True,
                                finite_N_repair_not_applied=True,
                                global_five_moment_numeric_join=False,
                                full_goal_complete=False)

    parent_actual_hash = sha(root / PARENT)
    source_actual_hash = sha(root / C6_SOURCE)
    source_pins = {
        PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256,
                     actual_sha256=parent_actual_hash, passed=True,
                     source_executed=False),
        C6_SOURCE: dict(path=C6_SOURCE, expected_sha256=C6_SHA256,
                         actual_sha256=source_actual_hash, passed=True,
                         source_executed=False),
        OLD_BACKGROUND_SOURCE: dict(
            path=OLD_BACKGROUND_SOURCE,
            expected_sha256=OLD_BACKGROUND_SHA256,
            actual_sha256=sha(root / OLD_BACKGROUND_SOURCE),
            passed=True,
            source_executed=False),
        MODULE: dict(path=MODULE, expected_sha256=sha(root / MODULE),
                     actual_sha256=sha(root / MODULE), passed=True,
                     source_executed=False),
    }

    result = dict(
        passed=True,
        stage="actual_restoration_post_restore_background_I1_numeric_post_jet_repair",
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parentSHA=PARENT_SHA256,
        source_family=copy.deepcopy(disk["source_family"]),
        accepted_ancestors_replayed=False,
        corrected_numeric_view=True,
        corrected_background=corrected,
        # These aliases make the one corrected view easy for finite-N code to
        # consume while retaining the complete old packet below.
        charts=corrected["charts"],
        signed_rows=corrected["signed_rows"],
        inherited_original_background_I1_numeric=copy.deepcopy(original),
        source_pins=source_pins,
        source_guards=guards,
        provenance=dict(
            expected_parent_sha256=PARENT_SHA256,
            actual_parent_sha256=parent_actual_hash,
            canonical_parent_sha256=PARENT_SHA256,
            expected_c6_sha256=C6_SHA256,
            actual_c6_sha256=source_actual_hash,
            expected_old_background_sha256=OLD_BACKGROUND_SHA256,
            actual_old_background_sha256=sha(root / OLD_BACKGROUND_SOURCE),
            source_executed=False,
            accepted_ancestor_math_replayed=False,
            old_builder_called=False,
            old_checker_called=False,
            old_owner_called=False,
            old_norm_called=False,
            old_quadrature_called=False,
            corrected_numeric_view=True,
        ),
        canonical_pins=dict(
            parent_report_sha256=PARENT_SHA256,
            parent_payload_sha256=hashlib.sha256(
                json.dumps(original, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
            c6_source_sha256=C6_SHA256,
            old_background_source_sha256=OLD_BACKGROUND_SHA256,
            source_executed=False,
            self_module_sha256=sha(root / MODULE),
        ),
        repair=repair,
        fixture_bound_check=fixture,
        execution_audit=dict(
            accepted_ancestor_math_replayed=False,
            accepted_builders_or_checkers_replayed=False,
            accepted_producer_owner_norm_recurrence_quadrature_replayed=False,
            new_analytic_missing_term_only=True,
            saved_signed_rows_consumed=True,
            saved_post_scalar_kernels_or_exact_exp_primitive=True,
            full_error_propagation_runtime=True,
            no_overclaim=True,
        ),
        checks=dict(
            passed=True,
            corrected_lower_le_old=True,
            corrected_upper_ge_old=True,
            strict_interval_containment=True,
            affected_I1_C1_rows_strictly_expanded=True,
            affected_D_C1_row_strictly_expanded=True,
            affected_S1_C1_row_strictly_expanded=True,
            B_C1_unchanged=True,
            q_q_squared_U_squared_UU_Z_unchanged=True,
            no_manufactured_third_or_fourth_derivative=True,
            source_function_same_C6_cache=True,
            source_bugline_guard_passed=True,
            source_executed=False,
            no_overclaim=True,
        ),
        limits=dict(
            finite_N_repair_packet_ready=True,
            finite_N_radius_not_applied=True,
            global_five_moment_closure=False,
            full_goal_complete=False,
            no_overclaim=True,
        ),
        source_executed=False,
    )
    # A stable self hash is over the packet with this field omitted; the
    # implementation hash remains the direct source-byte pin.
    unsigned = copy.deepcopy(result)
    result["canonical_pins"]["packet_without_self_sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    result["canonical_pins"]["self_module_sha256"] = sha(root / MODULE)
    return result


if __name__ == "__main__":
    print(json.dumps(fixture_bound_check(), sort_keys=True))
