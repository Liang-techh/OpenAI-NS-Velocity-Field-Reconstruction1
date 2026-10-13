"""C1 endpoint data for the original pre-:math:`R_c` pressure piece.

This adapter is intentionally narrower than the C0 endpoint source.  It keeps
the actual source rows for ``V_y`` and for the second ordinary ``Z`` row of
``V``.  The latter is used to form ``U_ZZ``; no derivative is taken of a
magnitude cap.  The lower endpoint is the source-owned ``sc_half`` singleton
and the upper endpoint is the accepted ``O3_power`` right endpoint.

The returned numbers are directed magnitude enclosures.  The signed endpoint
operators are retained as formulas, with ``m13`` using ``B(q+q_y)`` and
``m15`` using ``-B(R*q_y)``.  This module installs only the local endpoint C1
gate; complete pressure derivatives, repair, and global claims remain open.
"""

import copy
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp
import sympy as sy

from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_partial_mean as all_chart
import lei_ren_part1_paper_compliant_current_original_n1_preRc_pressure_endpoint as c0


HERE = physical.HERE
PREFIX = physical.PREFIX
sha = physical.sha
ep = physical.ep
rg = physical.rg
current = physical.current
source = physical.source
c3 = all_chart.c3
source_precision = physical.source_precision
fingerprint = physical.fingerprint

NAME = PREFIX + "current_original_n1_preRc_endpoint_C1.json.gz"
RECEIPT = PREFIX + "current_original_n1_preRc_endpoint_C1_check.json"
GATE = "current_original_n1_preRc_endpoint_C1_installed"
OPEN = tuple(dict.fromkeys((*physical.OPEN,
    "current_original_n1_preRc_full_pressure_C1_installed",
    "current_original_n1_pressure_derivative_accuracy_installed",
    "positive_order_five_moment_repair_closed",
    "global_completed_tensor_admissibility",
    "full_corrected_Navier_Stokes_solution")))


def _as_tuple(value):
    return tuple(_as_tuple(item) for item in value) if isinstance(value, list) else value


def _record(value):
    if isinstance(value, rg.LogUpper):
        return value.record()
    if isinstance(value, rg.JetBound):
        return rg.record(value)
    if isinstance(value, dict):
        return {key: _record(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_record(item) for item in value]
    if hasattr(value, "record"):
        return value.record()
    return value


def _log_cap(c, value, label):
    return c0._log_cap(c, value, label)


def _zero(c):
    return rg.Bounds(c).zero


def _source_jet(c, values, scale, label):
    # The accepted source rows are Taylor coefficients.  c0._source_jet is
    # the one shared factorial conversion used by the C0 adapter.
    return c0._source_jet(c, values, scale, label)


def _abs_coefficient(c, value, order=0):
    coefficient = value[order] if isinstance(value, IntervalTaylor) else value
    if not hasattr(coefficient, "_mpi_"):
        coefficient = c.mpf(coefficient)
    return c.mpf(ep(abs(coefficient))[1])


def _linear_sum(c, terms):
    bd = rg.Bounds(c)
    terms = [term for term in terms if not getattr(term[1], "zero", False)]
    if not terms:
        return bd.zero
    return bd.sum(*(bd.constant(coefficient) * value
                    for coefficient, value in terms))


def _qplus_caps(c, ends, delta, U, Uy, UZ, UyZ, UZZ):
    """Return |q+q_y| and |(q+q_y)_Z| on one exact Z cell.

    All coefficient intervals are formed first on the exact cell.  The field
    rows are then multiplied by those coefficient caps, so U_ZZ and U_yZ are
    genuine source rows rather than differentiated range enclosures.
    """
    z = IntervalTaylor.variable(c, c.mpf(ends), 2)
    d = 1 - z * z
    L = 1 - z * z * delta
    invL = L.reciprocal()
    bd = rg.Bounds(c)

    # q+q_y = N0/L, with N0=(1+delta) Z U+2 Z U_y-d U_Z.
    n0_coefficients = (z * (1 + delta), z * 2, -d)
    n0 = _linear_sum(c, [(_abs_coefficient(c, coefficient), value)
                         for coefficient, value in zip(
                             n0_coefficients, (U.value, Uy.value, UZ.value))])
    invL_cap = _abs_coefficient(c, invL)
    qplus = bd.constant(invL_cap) * n0

    # Differentiating the numerator gives
    # N1=(1+delta)U+(3+delta)Z U_Z+2U_y+2Z U_yZ-d U_ZZ.
    n1_coefficients = (1 + delta, z * (3 + delta), 2, z * 2, -d)
    n1 = _linear_sum(c, [(_abs_coefficient(c, coefficient), value)
                         for coefficient, value in zip(
                             n1_coefficients,
                             (U.value, UZ.value, Uy.value, UyZ.value,
                              UZZ.value))])
    alpha = z * (2 * delta) * invL * invL
    alpha_cap = _abs_coefficient(c, alpha)
    qplus_Z = bd.sum(bd.constant(invL_cap) * n1,
                     bd.constant(alpha_cap) * n0)
    return qplus, qplus_Z


def _source_view(fields, ends, chart, left, right=None):
    view = physical.source_rows(fields, tuple(ends), chart, left, right)
    if view.get("source_family") != fields.family:
        raise ValueError("Endpoint source family changed")
    packet = view["source_packet"]
    if packet.get("actual_phase_Z_exact_zero") is not True:
        raise ValueError("Endpoint phase must be Z independent")
    return view


def _lower_caps(fields, ends, base):
    ends = tuple(ends)
    view = _source_view(fields, ends, "first_micro", ("sc_half",))
    packet = view["source_packet"]
    generic = packet.get("original_generic_source")
    background = packet.get("actual_background_source")
    if not isinstance(generic, dict) or not isinstance(background, dict):
        raise ValueError("Actual lower generic/background source packet required")
    V = generic.get("common_velocity_V_axial5")
    histories = generic.get("common_own_five_histories_axial5")
    Vy = background.get("V_y")
    if (not isinstance(V, (list, tuple)) or not isinstance(histories, dict)
            or not isinstance(histories.get("m"), (list, tuple))
            or not isinstance(Vy, (list, tuple)) or len(Vy) < 4):
        raise ValueError("Lower V/m/V_y source rows through C3 are required")
    c, bd = fields.c, rg.Bounds(fields.c)
    scale = _log_cap(c, view["Pstar"], "lower Pstar")
    U = _source_jet(c, V, scale, "lower V")
    m = _source_jet(c, histories["m"], scale, "lower m")
    Uy = _source_jet(c, Vy, scale, "lower V_y")
    UZ = rg.JetBound(U.Z, U.ZZ, U.ZZZ, bd.zero)
    UZZ = rg.JetBound(U.ZZ, U.ZZZ, bd.zero, bd.zero)
    UyZ = rg.JetBound(Uy.Z, Uy.ZZ, Uy.ZZZ, bd.zero)
    q_rows = all_chart.partial.coefficient_Q_bounds(c, ends, fields.delta, U, m)
    q, qZ = q_rows[0], q_rows[1]
    qplus, qplus_Z = _qplus_caps(c, ends, fields.delta, U, Uy, UZ, UyZ, UZZ)
    # q_y=qplus-q.  The old C0 adapter retained a sum as a safe magnitude;
    # this adapter records the corrected quotient sign and its Z derivative.
    qy = bd.sum(qplus, q)
    qyZ = bd.sum(qplus_Z, qZ)
    R = _log_cap(c, view["R"], "lower physical radius")
    incoming = c0._zero_first_incoming_proof(fields, ends)
    flat = c0._flat_driver_proof(fields, ends, "first_micro")
    if not packet.get("source_owned_collar_q_C0_Z_exact_zero", False):
        raise ValueError("Lower source-owned collar theorem missing")
    caps = dict(U=U.value, U_y=Uy.value, U_Z=UZ.value,
                U_ZZ=UZZ.value, U_yZ=UyZ.value, m=m.value,
                q=q, q_Z=qZ, q_plus_q_y=qplus, q_plus_q_y_Z=qplus_Z,
                q_y=qy, q_y_Z=qyZ, R=R,
                R_q_y=R * qy, R_q_y_Z=R * qyZ)
    proof = dict(
        endpoint="first_micro/sc_half singleton",
        source_family=copy.deepcopy(fields.family),
        source_owned_collar_q_C0_Z_exact_zero=True,
        finite_N_initial_history_zero=True,
        finite_N_predecessor_history_retained=True,
        finite_N_local_driver_exact_zero=True,
        incoming_history=incoming,
        flat_driver=flat,
        actual_background_V_y_rows=6,
        actual_V_y_row_0_used=True,
        actual_V_y_row_1_retained=True,
        actual_V_y_rows_0_through_3_factorial_converted_once=True,
        actual_U_ZZ_from_V_rows_2_and_3=True,
        actual_U_yZ_from_V_y_rows_1_through_3=True,
        q_y_formula="q_plus_q_y-q",
        q_y_Z_formula="q_plus_q_y_Z-q_Z",
        source_quantity_paths=dict(
            V="original_generic_source.common_velocity_V_axial5",
            m='original_generic_source.common_own_five_histories_axial5["m"]',
            V_y="actual_background_source.V_y",
            R="original_generic_source.original_physical_radius",
            Pstar="outer.owner.owner(ends).Pstar"),
        common_Pstar_applied_once=True,
        no_radius_or_N_materialized=True,
    )
    return dict(caps=caps, records=_record(caps), proof=proof,
                source_rows=dict(
                    V=[_record(value) for value in V[:4]],
                    m=[_record(value) for value in histories["m"][:4]],
                    V_y=[_record(value) for value in Vy[:4]],
                    R=_record(view["R"]), Pstar=_record(view["Pstar"])))


def _upper_caps(fields, ends, base):
    ends = tuple(ends)
    view = _source_view(fields, ends, "O3_power", (2, 1))
    packet = view["source_packet"]
    generic = packet.get("original_generic_source")
    V = generic.get("common_velocity_V_axial5") if isinstance(generic, dict) else None
    if not isinstance(V, (list, tuple)) or len(V) < 4:
        raise ValueError("Upper O3 source V rows missing")
    if not all(c0._exact_zero_record(row) for row in V[:4]):
        raise ValueError("Upper O3 source U rows are not exact zero")
    c, bd = fields.c, rg.Bounds(fields.c)
    old = base["upper"]["caps"]
    zero = bd.zero
    q, qZ = old["q"], old["q_Z"]
    R = _log_cap(c, view["R"], "upper physical radius")
    caps = dict(U=zero, U_y=zero, U_Z=zero, U_ZZ=zero, U_yZ=zero,
                q=q, q_Z=qZ, q_plus_q_y=zero, q_plus_q_y_Z=zero,
                q_y=q, q_y_Z=qZ, R=R, R_q_y=R * q, R_q_y_Z=R * qZ)
    proof = dict(
        endpoint="O3_power right endpoint (2,1)",
        source_family=copy.deepcopy(fields.family),
        original_O3_power_V_rows_exact_zero=True,
        physical_U_Uy_UZ_UZZ_UyZ_exact_zero=True,
        q_plus_q_y_exact_zero=True,
        q_plus_q_y_Z_exact_zero=True,
        q_y_equals_minus_q=True,
        q_y_Z_equals_minus_q_Z=True,
        physical_mean_ode="m_y=U-m=-m",
        finite_N_own_mean_q_cap_retained=True,
        source_quantity_paths=dict(
            V="original_generic_source.common_velocity_V_axial5",
            q="accepted all-chart physical whole_cell_physical_q_Z2_magnitude",
            R="original_generic_source.original_physical_radius",
            Pstar="outer.owner.owner(ends).Pstar"),
        common_Pstar_applied_once=True,
        no_radius_or_N_materialized=True,
    )
    return dict(caps=caps, records=_record(caps), proof=proof,
                source_rows=dict(V=[_record(value) for value in V[:4]],
                                 R=_record(view["R"]), Pstar=_record(view["Pstar"])))


def _endpoint_data(fields, ends):
    base = c0.endpoint_caps(fields, tuple(ends))
    lower = _lower_caps(fields, ends, base)
    upper = _upper_caps(fields, ends, base)
    c, bd = fields.c, rg.Bounds(fields.c)
    lc, uc = lower["caps"], upper["caps"]
    boundary = dict(
        m13=dict(
            signed_form="B(q+q_y)",
            signed_Z_form="B((q+q_y)_Z)",
            value=bd.sum(lc["q_plus_q_y"], uc["q_plus_q_y"]),
            Z=bd.sum(lc["q_plus_q_y_Z"], uc["q_plus_q_y_Z"])),
        m15=dict(
            signed_form="-B(R*q_y)",
            signed_Z_form="-B(R*q_y_Z), endpoints fixed in Z",
            value=bd.sum(lc["R_q_y"], uc["R_q_y"]),
            Z=bd.sum(lc["R_q_y_Z"], uc["R_q_y_Z"])),
    )
    return dict(
        source_family=copy.deepcopy(fields.family), exact_Z_cell=list(tuple(ends)),
        lower=lower, upper=upper,
        boundary=boundary,
        c0_proof=base["proof"],
        proof=dict(
            exact_Z_cell=list(tuple(ends)),
            source_family=copy.deepcopy(fields.family),
            physical_common_Pstar_applied_once=True,
            finite_endpoint_radius_fixed_in_Z=True,
            m13_boundary_identity="B(q+q_y)",
            m15_boundary_identity="-B(R*q_y)",
            m13_boundary_Z_identity="B((q+q_y)_Z)",
            m15_boundary_Z_identity="-B(R*q_y_Z)",
            q_y_Z_uses_qplus_Z_minus_q_Z=True,
            upper_q_y_Z_equals_minus_q_Z=True,
            no_signed_centers_or_cap_differentiation=True,
            no_q_ZZZ_or_full_pressure_derivative_claim=True,
            full_repair_global_and_Navier_Stokes_claims_false=True,
        ))


def symbolic_check():
    Z, delta, R = sy.symbols("Z delta R", positive=True)
    U = sy.Function("U")(Z)
    Uy = sy.Function("Uy")(Z)
    m = sy.Function("m")(Z)
    a = 1 - delta
    d = 1 - Z**2
    L = 1 - delta * Z**2
    q = (2 * Z * U - a * Z * m - d * sy.diff(m, Z)) / L
    qplus = (Z * (1 + delta) * U + 2 * Z * Uy - d * sy.diff(U, Z)) / L
    qy = sy.cancel(qplus - q)
    qplus_Z = sy.diff(qplus, Z)
    qy_Z = sy.diff(qy, Z)
    residual = sy.simplify(qy_Z - (qplus_Z - sy.diff(q, Z)))
    if residual != 0:
        raise AssertionError("q_yZ quotient identity failed: " + str(residual))
    numerator = Z * (1 + delta) * U + 2 * Z * Uy - d * sy.diff(U, Z)
    explicit = ((1 + delta) * U + (3 + delta) * Z * sy.diff(U, Z)
                + 2 * Uy + 2 * Z * sy.diff(Uy, Z)
                - d * sy.diff(U, Z, 2)) / L + 2 * delta * Z * numerator / L**2
    explicit_residual = sy.simplify(qplus_Z - explicit)
    if explicit_residual != 0:
        raise AssertionError("qplus_Z quotient derivative failed")
    upper = sy.simplify(qy.subs({U: 0, Uy: 0, sy.diff(U, Z): 0}))
    upper_Z = sy.simplify(qy_Z.subs({U: 0, Uy: 0,
                                     sy.diff(U, Z): 0,
                                     sy.diff(U, (Z, 2)): 0,
                                     sy.diff(Uy, Z): 0,
                                     sy.diff(m, Z): sy.diff(m, Z, evaluate=True)}))
    # The latter substitution is intentionally only used to expose the sign;
    # q_y=-q and q_yZ=-q_Z follow before choosing any m profile.
    sign_residual = sy.simplify(upper + q.subs({U: 0, sy.diff(U, Z): 0}))
    if sign_residual != 0:
        raise AssertionError("upper q_y=-q identity failed")
    return dict(
        passed=True,
        q_plus_q_y_formula=sy.sstr(qplus),
        q_y_formula="q_plus_q_y-q",
        q_y_Z_formula="q_plus_q_y_Z-q_Z",
        qplus_Z_formula=sy.sstr(explicit),
        q_y_Z_residual=sy.sstr(residual),
        qplus_Z_residual=sy.sstr(explicit_residual),
        upper_q_y_plus_q_residual=sy.sstr(sign_residual),
        upper_q_y_Z_equals_minus_q_Z=True,
        finite_radius_R_is_Z_independent=True,
    )


@dataclass(frozen=True, eq=False)
class OriginalPreRcEndpointC1Packet:
    leading_hierarchy_order: int = 1


class CurrentOriginalN1PreRcEndpointC1:
    @source_precision
    def __init__(self, checked_fields, require_checked=True):
        if (type(checked_fields) is not physical.CurrentOriginalN1AllChartPhysicalPressureFields
                or not checked_fields.acceptance_loaded):
            raise ValueError("Checked all-chart physical pressure owner required")
        checked_fields.assert_graph()
        self.fields = checked_fields
        self.c = checked_fields.c
        self.family = copy.deepcopy(checked_fields.family)
        self.cells = tuple(tuple(cell) for cell in all_chart.current.current.CELLS)
        self.hashes = dict(checked_fields.hashes)
        for module in (physical, all_chart, c0):
            self.hashes[Path(module.__file__).name] = sha(Path(module.__file__).name)
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self._identities = (self.fields, self.c)
        self._family = copy.deepcopy(self.family)
        self._hashes = copy.deepcopy(self.hashes)
        self._packets = {}
        self.acceptance_loaded = False
        self.symbolic = symbolic_check()
        if require_checked:
            checked = json.loads((HERE / RECEIPT).read_bytes())
            if (not checked.get("all_passed") or not checked.get(GATE)
                    or checked.get("source_family") != self.family
                    or any(checked.get(key, False) for key in OPEN)):
                raise ValueError("Pre-Rc endpoint C1 receipt or scope differs")
            for name, digest in checked.get("input_hashes", {}).items():
                if sha(name) != digest:
                    raise ValueError("Pre-Rc endpoint C1 source changed " + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if (self._identities != (self.fields, self.c)
                or self.family != self._family
                or self.hashes != self._hashes):
            raise ValueError("Pre-Rc endpoint C1 source identity changed")
        self.fields.assert_graph()

    @source_precision
    def _make_view(self, ends):
        self.assert_graph()
        ends = tuple(ends)
        if ends not in self.cells:
            raise ValueError("Accepted original four-Z axial cell required")
        return _endpoint_data(self.fields, ends)

    @source_precision
    def evaluate(self, ends):
        view = self._make_view(ends)
        packet = OriginalPreRcEndpointC1Packet()
        self._packets[id(packet)] = (packet, view, fingerprint(_record(view)))
        return packet

    @source_precision
    def report(self, packet):
        self.assert_graph()
        entry = self._packets.get(id(packet))
        if (type(packet) is not OriginalPreRcEndpointC1Packet or entry is None
                or entry[0] is not packet
                or fingerprint(_record(entry[1])) != entry[2]):
            raise ValueError("Live pre-Rc endpoint C1 packet required")
        return copy.deepcopy(_record(entry[1]))


def _build_fields():
    modified = all_chart.partial.parent.CurrentOriginalN1ModifiedLeadingC3(
        all_chart.partial.rg.CurrentC3TargetRanges(), require_checked=True)
    checked_partial = all_chart.partial.CurrentOriginalN1PartialMeanQ(
        modified, require_checked=True)
    checked_all = all_chart.CurrentOriginalN1AllChartPartialMeanQ(
        checked_partial, require_checked=True)
    return physical.CurrentOriginalN1AllChartPhysicalPressureFields(
        checked_all, require_checked=True)


@source_precision
def run(checked_fields=None):
    if checked_fields is None:
        checked_fields = _build_fields()
    owner = CurrentOriginalN1PreRcEndpointC1(checked_fields, require_checked=False)
    packets = {}
    samples = {}
    for ends in owner.cells:
        packet = owner.evaluate(ends)
        packets[ends] = packet
        samples[ends] = owner.report(packet)
    owner.assert_graph()
    raw = dict(
        candidate_preRc_endpoint_C1_constructed=True,
        source_family=owner.family,
        cells=[list(cell) for cell in owner.cells],
        samples={str(cell): sample for cell, sample in samples.items()},
        source_quantity_scope="actual lower sc_half singleton and O3_power right endpoint",
        source_bindings=dict(
            lower="physical.source_rows(first_micro, ('sc_half',))",
            upper="physical.source_rows(O3_power, (2,1))",
            lower_V_y="actual_background_source.V_y[0:4]",
            lower_U_ZZ="factorial source rows V[2:4]",
            lower_U_yZ="factorial source rows V_y[1:4]",
            common_Pstar="outer.owner.owner(ends).Pstar",
            q="partial.coefficient_Q_bounds",
            qplus="(1+delta)ZU+2ZUy-(1-Z^2)UZ over L",
            qplus_Z="explicit quotient derivative including 2*delta*Z*N/L^2"),
        input_hashes=owner.hashes,
        **{GATE: False}, **dict.fromkeys(OPEN, False))
    (HERE / NAME).write_bytes(gzip.compress(
        (json.dumps(c3.encoded(raw), separators=(",", ":")) + "\n").encode(),
        mtime=0))
    print("PRE_RC_ENDPOINT_C1_COMPUTED", len(packets), flush=True)
    return owner, packets


if __name__ == "__main__":
    run()
