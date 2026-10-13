"""Checked same-source inputs for the pre-:math:`R_p` pressure C1 route.

This adapter keeps three pieces separate:

* live 4 by 17 ordinary baseline ``U`` and ``m`` rows from the original
  physical source provider;
* the correlated ``R_c`` pressure-unit logarithms from the checked pulse
  owner; and
* a differentiated finite-N perturbation bound on the accepted ``R_c`` to
  ``2 R_c`` C3 repair band.

The compact-band bound differentiates the radial-free pressure identity from
the finite source rows.  It never differentiates an old magnitude cap and it
never materializes ``R_c``, ``N``, or any of the large physical exponentials.
The quiet ``2 R_c`` to ``R_p`` source is retained as an exact Z-independent
record.  This file installs only an input gate for a later pre-``R_p`` C1
adapter; full pressure derivative, repair, and global gates remain open.
"""
import copy
import gzip
import json
import math
import time
from functools import wraps
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_C3_repair_band as c3_repair
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q as partial
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source as preRp
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_source as pulse
import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source as quiet


HERE = physical.HERE
PREFIX = physical.PREFIX
sha = physical.sha
ep = physical.ep
rg = physical.rg
source_precision = physical.source_precision
SOURCE_DPS = 540

NAME = PREFIX + "current_original_n1_preRp_C1_inputs.json.gz"
RECEIPT = PREFIX + "current_original_n1_preRp_C1_inputs_check.json"
GATE = "current_original_n1_preRp_C1_inputs_installed"
OPEN = (
    "current_original_n1_preRp_C1_pressure_values_installed",
    "complete_signed_pressure_tails_certified",
    "n1_second_and_third_repair_closed",
    "global_physical_pressure_defect_installed",
    "temporal_recursion",
    "full_corrected_Navier_Stokes_solution",
)
CELLS = tuple(tuple(item) for item in pulse.CELLS)
SCRATCH_NAME = "tmp/n1_baseline_preRp_C3_20261013.json"


def _record(value):
    if isinstance(value, rg.LogUpper):
        # ``LogUpper.record`` retains its directed interval as an ``ivmpf``;
        # recurse through the accepted pulse serializer so JSON artifacts pin
        # exact endpoint tuples rather than live interval objects.
        return _record(value.record())
    if isinstance(value, rg.JetBound):
        return _record(rg.record(value))
    if isinstance(value, dict):
        return {key: _record(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_record(item) for item in value]
    if hasattr(value, "_mpi_"):
        return pulse.record(value)
    if isinstance(value, mp.mpf):
        return mp.nstr(value, n=max(80, mp.mp.dps))
    if hasattr(value, "record"):
        return _record(value.record())
    return value


def _read_json_gz(path):
    return json.loads(gzip.decompress(Path(path).read_bytes()))


def _interval(c, value):
    if hasattr(value, "_mpi_"):
        return c.make_mpf(value._mpi_)
    return c.make_mpf((tuple(value["lower_exact_mpf_tuple"]),
                       tuple(value["upper_exact_mpf_tuple"])))


def _raw_bound(c, value):
    if value.get("exact_zero", False):
        return rg.Bounds(c).zero
    return rg.LogUpper(c, _interval(c, value["log_absolute_upper"]))


def _raw_jet(c, values, label):
    if len(values) < 4:
        raise ValueError("Four ordinary source rows required for " + label)
    return rg.JetBound(*(_raw_bound(c, item) for item in values[:4]))


def _archive(module, gates, family):
    raw = _read_json_gz(module.HERE / module.NAME)
    receipt = json.loads((module.HERE / module.RECEIPT).read_bytes())
    if not receipt.get("all_passed") or raw.get("source_family") != family:
        raise ValueError("Checked same-family archive required: " + module.NAME)
    if receipt.get("source_family") != family:
        raise ValueError("Archive receipt family changed: " + module.NAME)
    for gate in gates:
        if not receipt.get(gate):
            raise ValueError("Required archive gate is absent: " + gate)
    hashes = {}
    for name, digest in receipt.get("input_hashes", {}).items():
        if sha(name) != digest:
            raise ValueError("Changed archive input " + name)
        hashes[name] = digest
    for name in (module.NAME, module.RECEIPT, Path(module.__file__).name):
        hashes[name] = sha(name)
    return raw, hashes


def _raw_ordinary_jet(c, bd, rows, label):
    if len(rows) < 4:
        raise ValueError("Four source Taylor rows required for " + label)
    return rg.JetBound(*(bd.row(rows[k] * math.factorial(k))
                         for k in range(4)))


def _same_source_family(fields, whole):
    if fields.family != whole.family:
        raise ValueError("Physical and pulse owners have different families")
    return copy.deepcopy(fields.family)


def _baseline(fields):
    """Build live baseline rows; no scratch value is used as a source."""
    fields.assert_graph()
    c = fields.c
    bd = rg.Bounds(c)
    result = {}
    for ends in CELLS:
        for chart in fields.charts:
            partition = tuple(fields.partitions[chart])
            source_view = physical.source_rows(
                fields, ends, chart, partition[0], partition[-1])
            op = fields.modified.ranges.phase.outer.owner.owner(tuple(ends))
            if source_view["Pstar"] is not op.Pstar:
                raise ValueError("Pstar source object changed")
            packet = source_view["source_packet"]
            if packet.get("actual_phase_Z_exact_zero") is not True:
                raise ValueError("Original phase is not Z independent")
            V = _raw_ordinary_jet(c, bd, source_view["V"], "baseline V")
            m = _raw_ordinary_jet(c, bd, source_view["m"], "baseline m")
            S = bd.row(source_view["Pstar"])
            U = bd.scaled(V, S)
            M = bd.scaled(m, S)
            key = "|".join(str(item) for item in ends) + "|" + chart
            result[key] = dict(
                exact_Z_cell=list(ends), chart=chart,
                exact_native_endpoints=[partition[0], partition[-1]],
                source_cell_count=len(partition) - 1,
                physical_U_Z3_magnitude=[item.record() for item in rg.rows(U)],
                physical_m_Z3_magnitude=[item.record() for item in rg.rows(M)],
                Pstar_record=S.record(),
                original_phase_Z_exact_zero=True,
                source_identity=fields.family,
                original_source_quantity_paths=dict(
                    V="original_generic_source.common_velocity_V_axial5",
                    m='original_generic_source.common_own_five_histories_axial5["m"]',
                    Pstar="outer.owner.owner(ends).Pstar",
                    R="original_generic_source.original_physical_radius"),
                baseline_rows_are_source_Taylor_coefficients=True,
                factorial_conversion_once=True,
                physical_scale_applied_once=True,
                z_independent_physical_scale=True,
                range_enclosure_not_function_value=True,
            )
    fields.assert_graph()
    return result


def _coefficient_envelopes(c):
    """Directed global coefficient and first-Z derivative envelopes."""
    z = c.mpf([-1, 1])
    delta = c.mpf([0, ".5"])
    a = 1 - delta
    d = 1 - z * z
    L = 1 - delta * z * z
    Lz = -2 * delta * z
    dz = -2 * z
    coeff = dict(
        invL=1 / L,
        az_over_2L=a * z / (2 * L),
        d_over_L=d / L,
        two_z2_over_L2=2 * z * z / (L * L),
        two_delta_z2_over_L2=2 * delta * z * z / (L * L),
        two_zd_over_L2=2 * z * d / (L * L),
    )
    deriv = dict(
        invL=-Lz / (L * L),
        az_over_2L=a * (L - z * Lz) / (2 * L * L),
        d_over_L=(dz * L - d * Lz) / (L * L),
        two_z2_over_L2=4 * z / (L * L) - 4 * z * z * Lz / (L * L * L),
        two_delta_z2_over_L2=4 * delta * z / (L * L)
        - 4 * delta * z * z * Lz / (L * L * L),
        two_zd_over_L2=2 * (d + z * dz) / (L * L)
        - 4 * z * d * Lz / (L * L * L),
    )
    return {key: max(abs(ep(coeff[key])[0]), abs(ep(coeff[key])[1]),
                     abs(ep(deriv[key])[0]), abs(ep(deriv[key])[1]))
            for key in coeff}


def _term_derivative(bd, envelopes, name, value, derivative, measure):
    coeff = bd.constant(envelopes[name])
    # ``value`` and ``derivative`` are scalar ``LogUpper`` records here;
    # ``Bounds.add`` is reserved for ``JetBound`` objects.  Use the scalar
    # directed sum explicitly so this path cannot dispatch to JetBound.add.
    return bd.scale(bd.sum(value, derivative), measure * envelopes[name])


def _compact_c1_cell(c, ends, raw, pstar, delta, k_log, log_rc):
    """Differentiate the finite-N compact perturbation in source units."""
    bd = rg.Bounds(c)
    profiles = raw["actual_selected_N_profile_y0_y1_y2_C3_bounds"]
    V = _raw_jet(c, profiles[0]["V"], "compact V")
    Vy = _raw_jet(c, profiles[1]["V"], "compact V_y")
    dm = _raw_jet(c, raw["actual_partial_correction_C3_bounds"]["m"],
                  "compact finite mean")
    U = bd.scaled(V, pstar)
    Uy = bd.scaled(Vy, pstar)
    m = bd.scaled(dm, pstar)
    my = bd.scaled(bd.add(V, dm), pstar)
    q = partial.coefficient_Q_bounds(c, ends, delta, U, m)
    qy = partial.coefficient_Q_bounds(c, ends, delta, Uy, my)
    dq, dqz, dqzz = q
    dqy, dqyz, _ = qy
    U0, U1, U2 = U.value, U.Z, U.ZZ
    uqz = U0 * dqz
    uqz_z = bd.sum(U1 * dqz, U0 * dqzz)
    u2 = U0 * U0
    u2_z = bd.scale(U0 * U1, 2)
    uuz = U0 * U1
    uuz_z = bd.sum(U1 * U1, U0 * U2)
    k = rg.LogUpper(c, k_log)
    q2_z = bd.sum(bd.scale(k * dqz, 2), bd.scale(dq * dqz, 2))
    env = _coefficient_envelopes(c)

    # x=R/Rc is in [1,2].  The factors 3, 5 and 2 below are the directed
    # B(x f), B(x^2 f), and integral[x dx] majorants.  These are applied to
    # differentiated quantities, never to an old C0 cap.
    s0 = bd.sum(
        _term_derivative(bd, env, "invL", dq, dqz, 3),
        _term_derivative(bd, env, "az_over_2L", dqz, dqzz, 1),
        bd.scale(q2_z, 3),
        _term_derivative(bd, env, "d_over_L", uqz, uqz_z, 1),
        _term_derivative(bd, env, "two_z2_over_L2", u2, u2_z, 5),
        _term_derivative(bd, env, "two_delta_z2_over_L2", u2, u2_z, 2),
        _term_derivative(bd, env, "two_zd_over_L2", uuz, uuz_z, 2),
    )
    s1 = bd.sum(
        _term_derivative(bd, env, "invL", dq, dqz, 5),
        _term_derivative(bd, env, "invL", dq, dqz, 2),
        _term_derivative(bd, env, "az_over_2L", dqz, dqzz, 2),
        bd.scale(q2_z, 5),
        bd.scale(q2_z, 2),
        _term_derivative(bd, env, "d_over_L", uqz, uqz_z, 2),
        _term_derivative(bd, env, "two_z2_over_L2", u2, u2_z, 5),
        _term_derivative(bd, env, "two_z2_over_L2", u2, u2_z, 2),
        _term_derivative(bd, env, "two_zd_over_L2", uuz, uuz_z, 2),
    )
    inv_rc = rg.LogUpper(c, c.mpf(-ep(log_rc)[0]))
    endpoint0 = inv_rc * bd.scale(bd.sum(dqz, dqyz), 3)
    endpoint1 = inv_rc * bd.scale(dqyz, 8)
    e13 = bd.sum(bd.scale(s0, c.mpf(".5")), endpoint0)
    e15 = bd.sum(bd.scale(s1, c.mpf(".5")), endpoint1)
    return dict(
        exact_Z_cell=list(ends), band_R_over_Rc=[1, 2],
        baseline_first_Z_derivative_exact_zero=True,
        finite_rows=dict(U=_record(U), U_y=_record(Uy), m=_record(m),
                         m_y=_record(my), q=_record(q), q_y=_record(qy)),
        differentiated_products=dict(
            q_squared_Z=q2_z.record(), U_q_Z_Z=uqz_z.record(),
            U_squared_Z=u2_z.record(), U_U_Z_Z=uuz_z.record()),
        endpoint_derivative_bounds=dict(
            q_plus_Rq_R_Z=endpoint0.record(),
            R_squared_q_R_Z=endpoint1.record(),
            inverse_Rc=inv_rc.record(),
            q_Z=dqz.record(), q_y_Z=dqyz.record()),
        reduced_S_Z_bounds=dict(unweighted=s0.record(), weighted=s1.record()),
        normalized_first_Z_upper=dict(m13=e13.record(), m15=e15.record()),
        coefficients={key: mp.nstr(value, n=90) for key, value in env.items()},
        physical_common_scale_Pstar_applied_once=True,
        selected_N_rows_source_bound=True,
        selected_N_division_not_repeated=True,
        finite_endpoint_q_and_q_y_derivatives_retained=True,
        no_differentiation_of_old_C0_caps=True,
    )


def _native_radius_source(c, family):
    # The accepted absolute native-frame binding is exported by the n1
    # repair-band pressure module retained by preRp, while its endpoint
    # adapter only provides endpoint caps/geometry.
    binding = preRp.band._accepted_native_radius_binding(family)
    log_c = physical.current.packets.interval(c, binding["selected_logCstar"])
    log_rc = c.ln(110) + 10 * log_c + 11 * c.exp(40) + 124
    return binding, log_rc


@source_precision
def _build_quiet_and_scale(c, whole, family):
    quiet_raw, quiet_hashes = _archive(
        quiet, (quiet.GATE,), family)
    qsrc = quiet_raw["actual_quiet_source"]
    if not qsrc.get("source_scaled_Z1_exact_zero"):
        raise ValueError("Quiet source Z1 identity is not exact")
    scales = {
        "m13": preRp.Rc_unit_log(c, whole, 1),
        "m15": preRp.Rc_unit_log(c, whole, 2),
    }
    quiet_record = dict(
        source_family=family,
        quiet_segment=qsrc["quiet_segment"],
        source_identity=qsrc["source_identity"],
        m13=qsrc["m13"], m15=qsrc["m15"], K=qsrc["K"],
        original_quiet_length=qsrc["original_quiet_length"],
        source_scaled_Z1_exact_zero=True,
        signed_coefficients_not_magnitude_caps=True,
        physical_radius_scale_retained_as_original_function=True,
        source_name=quiet.NAME, receipt_name=quiet.RECEIPT,
    )
    return quiet_record, scales, quiet_hashes


class OriginalPreRpC1InputsPacket:
    pass


class CurrentOriginalN1PreRpC1Inputs:
    @source_precision
    def __init__(self, fields, whole, repair_owner=None, require_checked=True):
        if type(fields) is not physical.CurrentOriginalN1AllChartPhysicalPressureFields:
            raise ValueError("Checked all-chart physical fields owner required")
        if type(whole) is not pulse.CurrentOriginalN1PulsePressureSource:
            raise ValueError("Checked pulse pressure owner required")
        if not fields.acceptance_loaded or not whole.acceptance_loaded:
            raise ValueError("Accepted physical and pulse owners required")
        fields.assert_graph(); whole.assert_graph()
        self.fields, self.pulse, self.c = fields, whole, fields.c
        self.family = _same_source_family(fields, whole)
        if repair_owner is not None:
            if type(repair_owner) is not c3_repair.CurrentC3RepairBand:
                raise ValueError("Checked C3 repair-band owner required")
            if not repair_owner.acceptance_loaded or repair_owner.identity != self.family:
                raise ValueError("C3 repair-band source family changed")
        self.repair_owner = repair_owner
        self.c3_raw, c3_hashes = _archive(
            c3_repair, c3_repair.GATES, self.family)
        self.baseline = _baseline(fields)
        self.quiet_source, self.scale_ratios, quiet_hashes = _build_quiet_and_scale(
            self.c, whole, self.family)
        self.radius_binding, self.log_Rc = _native_radius_source(self.c, self.family)
        self.compact = {}
        rows = self.c3_raw["actual_four_Z_C3_band_ranges"]
        if {tuple(row["exact_Z_cell"]) for row in rows} != set(CELLS):
            raise ValueError("C3 repair-band archive does not cover four cells")
        for row in rows:
            ends = tuple(row["exact_Z_cell"])
            pstar = rg.Bounds(self.c).row(
                physical.physical_source_scale(fields, ends))
            key = "|".join(str(item) for item in ends)
            self.compact[key] = _compact_c1_cell(
                self.c, ends, row, pstar, fields.delta,
                quiet.read_interval(self.c,
                                    self.quiet_source["K"]["log_coefficient"]),
                self.log_Rc)
        self.hashes = dict(fields.hashes)
        for source_hashes in (whole.hashes, c3_hashes, quiet_hashes):
            for name, digest in source_hashes.items():
                if name in self.hashes and self.hashes[name] != digest:
                    raise ValueError("Conflicting pre-Rp C1 source " + name)
                self.hashes[name] = digest
        for module in (physical, pulse, partial, preRp, quiet, c3_repair):
            for name in (Path(module.__file__).name,
                         getattr(module, "NAME", None),
                         getattr(module, "RECEIPT", None)):
                if name:
                    self.hashes[name] = sha(name)
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self.acceptance_loaded = False
        if require_checked:
            raw = _read_json_gz(HERE / NAME)
            receipt = json.loads((HERE / RECEIPT).read_bytes())
            if (not receipt.get("all_passed") or not receipt.get(GATE)
                    or receipt.get("source_family") != self.family
                    or any(receipt.get(key, False) for key in OPEN)):
                raise ValueError("Checked pre-Rp C1 inputs receipt required")
            for name, digest in receipt.get("input_hashes", {}).items():
                if sha(name) != digest:
                    raise ValueError("Changed pre-Rp C1 input " + name)
                self.hashes[name] = digest
            if raw.get("source_family") != self.family:
                raise ValueError("Pre-Rp C1 raw source family changed")
            self.acceptance_loaded = True
        self._state = _record((self.family, self.baseline, self.compact,
                               self.scale_ratios, self.quiet_source,
                               self.hashes))
        self._hashes = copy.deepcopy(self.hashes)
        self._packets = {}

    @source_precision
    def assert_graph(self):
        self.fields.assert_graph(); self.pulse.assert_graph()
        if self.hashes != self._hashes:
            raise ValueError("Pre-Rp C1 input source hashes changed")
        if self._state != _record((self.family, self.baseline, self.compact,
                                   self.scale_ratios, self.quiet_source,
                                   self.hashes)):
            raise ValueError("Pre-Rp C1 input source state changed")

    @source_precision
    def evaluate(self, ends):
        self.assert_graph()
        key = "|".join(str(item) for item in tuple(ends))
        if key not in self.compact:
            raise ValueError("Accepted original axial cell required")
        packet = OriginalPreRpC1InputsPacket()
        view = dict(
            source_family=self.family, exact_Z_cell=list(tuple(ends)),
            baseline={k: v for k, v in self.baseline.items()
                      if k.startswith(key + "|")},
            compact=self.compact[key], scale_ratios=_record(self.scale_ratios),
            quiet_source=copy.deepcopy(self.quiet_source),
            absolute_Rc_source=dict(
                native_frame_nodes=self.radius_binding["expected_nodes"],
                log_Rc=_record(self.log_Rc), radius_not_materialized=True),
            C1_inputs_only=True, no_signed_preRp_C1_center=True,
            complete_pressure_C1_not_installed=True,
            **{GATE: self.acceptance_loaded}, **dict.fromkeys(OPEN, False))
        self._packets[id(packet)] = (packet, view, _record(view))
        return packet

    @source_precision
    def evaluate_all(self):
        return {"|".join(key): self.evaluate(key) for key in CELLS}

    @source_precision
    def report(self, packet):
        self.assert_graph()
        item = self._packets.get(id(packet))
        if item is None or item[0] is not packet or item[2] != _record(item[1]):
            raise ValueError("Live pre-Rp C1 input packet required")
        return copy.deepcopy(item[1])


@source_precision
def run(fields, whole, repair_owner=None):
    began = time.monotonic()
    owner = CurrentOriginalN1PreRpC1Inputs(fields, whole, repair_owner,
                                            require_checked=False)
    packets = owner.evaluate_all()
    raw = dict(
        source_family=owner.family,
        baseline=_record(owner.baseline),
        compact=_record(owner.compact),
        scale_ratios=_record(owner.scale_ratios),
        quiet_source=owner.quiet_source,
        absolute_Rc_source=dict(
            native_frame_nodes=owner.radius_binding["expected_nodes"],
            log_Rc=_record(owner.log_Rc), radius_not_materialized=True),
        exact_Z_cells=[list(key) for key in CELLS],
        baseline_scope="4 axial cells x 17 exact whole-chart unions",
        source_bindings=dict(
            baseline_physical_source="physical.source_rows(fields,ends,chart,partition[0],partition[-1])",
            finite_repair_rows="accepted C3 repair-band actual_four_Z_C3_band_ranges",
            q_Z2="partial.coefficient_Q_bounds(c,ends,delta,U,m)",
            Rc_unit_log="preRp.Rc_unit_log(c,whole,power)",
            quiet_source="accepted quiet actual_quiet_source; source_scaled_Z1_exact_zero"),
        input_hashes=owner.hashes,
        **{GATE: False}, **dict.fromkeys(OPEN, False),
        execution_seconds=time.monotonic() - began)
    (HERE / NAME).write_bytes(gzip.compress(
        (json.dumps(raw, separators=(",", ":")) + "\n").encode(), mtime=0))
    print("ORIGINAL_N1_PRE_RP_C1_INPUTS_COMPUTED", len(owner.baseline),
          len(owner.compact), flush=True)
    return owner, packets


if __name__ == "__main__":
    raise SystemExit("run(fields, whole, repair_owner) requires retained owners")
