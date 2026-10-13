"""Signed n1 pressure moments on the accepted ``Rc``--``2 Rc`` band.

The C3 repair-band source already contains the actual finite-N velocity and
own-mean histories.  This adapter binds those histories to the physical
common unit once and uses the radial-free identities to make a signed value
enclosure.  The leading compact-band source is the accepted quiet coefficient
``k=K/Rc``; its exact values are ``k**2/8`` and ``-k**2*log(2)/4`` after the
``Rc`` and ``Rc**2`` factors are removed.  The finite-N part is carried as a
directed perturbation bound, including the two finite endpoint ``q_R`` terms.

This file deliberately stops at the compact C3 band.  It does not assert an
outer support endpoint, a higher-n repair, or a global pressure defect.
"""
import copy
from functools import wraps
import gzip
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_C3_repair_band as repair
import lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp as power
import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source as quiet
import lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q as partial
import lei_ren_part1_paper_compliant_current_original_n1_radial_free_pressure_integrals as radial_free


target = repair.target
current = repair.current
rg = repair.ranges
source = partial.source
HERE, PREFIX, sha, ep = repair.HERE, repair.PREFIX, repair.sha, repair.ep

NAME = PREFIX + 'current_original_n1_repair_band_pressure_source.json.gz'
RECEIPT = PREFIX + 'current_original_n1_repair_band_pressure_source_check.json'
GATE = 'current_original_n1_repair_band_pressure_source_installed'
OPEN = ('complete_signed_pressure_tails_certified',
        'n1_second_and_third_repair_closed',
        'outer_support_endpoint_installed',
        'global_physical_pressure_defect_installed')
NATIVE_FRAME_NAME = PREFIX + 'current_original_C3_Rp_native_frame.json.gz'
NATIVE_FRAME_RECEIPT = PREFIX + 'current_original_C3_Rp_native_frame_check.json'
NORM_NAME = PREFIX + 'physical_norm_family.json'
SOURCE_DPS = 540


def _source_precision_540(fn):
    """Run public source operations in the accepted directed precision."""
    @wraps(fn)
    def call(*args, **kwargs):
        with mp.workdps(SOURCE_DPS):
            return fn(*args, **kwargs)
    return call


def _as_ref(g, value):
    """Bind a source reference to the graph used by the pressure owner."""
    if isinstance(value, source.FunctionRef):
        return source.FunctionRef(g, value.node)
    if isinstance(value, int):
        return source.FunctionRef(g, value)
    if hasattr(value, 'node'):
        return source.FunctionRef(g, value.node)
    return value


def _clone_c3(g, value):
    """Copy a C3 function's node references into the power graph."""
    if isinstance(value, target.C3Function):
        return target.C3Function(*(_as_ref(g, q) for q in target.rows(value)))
    return value


def _clone_dict(g, value):
    if isinstance(value, target.C3Function):
        return _clone_c3(g, value)
    if isinstance(value, dict):
        return {key: _clone_dict(g, item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return type(value)(_clone_dict(g, item) for item in value)
    return value


def _read_json_gz(path):
    return json.loads(gzip.decompress(Path(path).read_bytes()))


def _accepted_native_radius_binding(family):
    """Load the accepted absolute ``Rc`` source identity without rebuilding it."""
    frame = _read_json_gz(HERE / NATIVE_FRAME_NAME)
    receipt = json.loads((HERE / NATIVE_FRAME_RECEIPT).read_bytes())
    gate = 'current_original_Rp_native_exact_radius_identity_installed'
    if (not receipt.get('all_passed') or not receipt.get(gate) or
            frame.get('source_family') != family or
            receipt.get('source_family') != family):
        raise ValueError('Accepted original absolute Rc/native-frame source required')
    for name, digest in receipt['input_hashes'].items():
        if sha(name) != digest:
            raise ValueError('Changed accepted native radius source ' + name)
    functions = frame['actual_Rp_native_frame_functions']
    expected = dict(original_Rc_phase_offset=4893,
                    original_log_r_minus=26684, original_logRc=26685)
    if any(functions.get(name) != node for name, node in expected.items()):
        raise ValueError('Accepted absolute Rc source node binding changed')
    norm = json.loads((HERE / NORM_NAME).read_bytes())
    selected = norm.get('selected_logCstar')
    if not isinstance(selected, dict) or selected.get('lower_exact_mpf_tuple') != selected.get('upper_exact_mpf_tuple'):
        raise ValueError('Original selected singleton logCstar required for absolute Rc')
    return dict(frame=frame, receipt=receipt, functions=functions,
                selected_logCstar=selected, expected_nodes=expected)


def _cell_key(ends):
    return tuple(str(item) for item in ends)


def _source_rows_equal(actual, expected):
    return current.current.previous.equivalent_rows(actual, expected)


def _raw_bound(c, record):
    """Decode one archived logarithmic magnitude row without materializing it."""
    if record.get('exact_zero'):
        return rg.LogUpper(c)
    return rg.LogUpper(c, current.packets.interval(c, record['log_absolute_upper']))


def _raw_jet(c, records):
    if len(records) < 4:
        raise ValueError('Four ordinary source rows required')
    return rg.JetBound(*(_raw_bound(c, item) for item in records[:4]))


def _record(value):
    if isinstance(value, rg.JetBound):
        return [item.record() for item in rg.rows(value)]
    if isinstance(value, rg.LogUpper):
        return value.record()
    if isinstance(value, dict):
        return {key: _record(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_record(item) for item in value]
    if hasattr(value, 'record'):
        return value.record()
    return value


def _log_sum(bd, *values):
    return bd.sum(*values)


def _scaled(bd, value, scalar):
    return bd.scale(value, scalar)


def _repair_rows(raw, key):
    return raw['actual_selected_N_profile_y0_y1_y2_C3_bounds'][0][key]


def _inverse_radius_upper(c, log_rc):
    """Directed upper bound for ``1/Rc`` from the absolute source log."""
    lower_log_rc = ep(log_rc)[0]
    if not mp.isfinite(lower_log_rc):
        raise ValueError('Finite absolute log Rc lower bound required')
    return rg.LogUpper(c, c.mpf(-lower_log_rc))


def _normalized_endpoint_errors(bd, dq, dqy, log_rc):
    """Normalize both finite endpoint terms by their physical Rc powers.

    With x=R/Rc, B(q+R q_R)/Rc and B(R^2 q_R)/Rc^2 each retain one
    factor 1/Rc.  The absolute source log is kept logarithmically.
    """
    inv_rc = _inverse_radius_upper(bd.c, log_rc)
    endpoint0 = inv_rc * _scaled(bd, _log_sum(bd, dq, dqy), 2)
    endpoint1 = inv_rc * _scaled(bd, dqy, 8)
    return endpoint0, endpoint1


def _uuz_bound(U):
    """Magnitude cap for the identity's ``U U_Z`` term."""
    return U.value * U.Z


def _delta_u2_bound(c, bd, u2, delta_value):
    """Retain a nonzero directed delta coefficient without float conversion."""
    delta_upper = c.mpf(ep(delta_value)[1])
    return _scaled(bd, u2, 8 * delta_upper)


def _finite_perturbation(c, ends, raw, delta, pstar, log_rc):
    """Return normalized magnitude errors for the two pressure moments.

    The input rows are source enclosures.  ``q`` is split into the exact
    compact-band profile ``q0=-k/x`` and a finite-N increment.  Every term in
    the radial-free formulas is bounded by a positive logarithmic enclosure;
    the q-linear weighted terms in the R*m15 identity are retained explicitly.
    """
    bd = rg.Bounds(c)
    S = pstar
    V = _raw_jet(c, _repair_rows(raw, 'V'))
    Vy = _raw_jet(c, raw['actual_selected_N_profile_y0_y1_y2_C3_bounds'][1]['V'])
    dm = _raw_jet(c, raw['actual_partial_correction_C3_bounds']['m'])
    U = bd.scaled(V, S)
    Uy = bd.scaled(Vy, S)
    m = bd.scaled(dm, S)
    # The finite part obeys m_y=V-m.  Using the separate partial m history
    # keeps the compact leading cancellation correlated; the archived
    # complete-history-y magnitude is deliberately not differentiated here.
    k = delta['k']
    zcap = bd.constant(1)
    m_y = bd.scaled(bd.add(V, dm), S)
    q = partial.coefficient_Q_bounds(c, ends, delta['value'], U, m)
    qy = partial.coefficient_Q_bounds(c, ends, delta['value'], Uy, m_y)
    dq, dqz = q[0], q[1]
    dqy = qy[0]
    # The source-independent coefficients below use |z|<=1, 1/L<=2 and
    # x=R/Rc in [1,2].  They are intentionally simple directed majorants;
    # no cap is differentiated or interpreted as a point value.
    delta_upper = ep(delta['value'])[1]
    if not delta_upper < c.mpf('.5'):
        raise ValueError('Original delta must satisfy delta < 1/2 for 1/L <= 2')
    kcap = k
    q2 = _log_sum(bd, _scaled(bd, kcap * dq, 2), dq * dq)
    u2 = U.value * U.value
    uqz = U.value * dqz
    # U U_Z uses the Z row of U itself.  Uy.Z would be U_{yZ} and is a
    # different derivative that does not occur in the radial-free identity.
    uuz = _uuz_bound(U)
    # S0 correction from the unweighted identity.
    s0 = _log_sum(bd,
                  _scaled(bd, dq, 6),                 # B(R q)/L, 2*3
                  dqz,                                  # a z I0(q_Z)/(2L)
                  _scaled(bd, q2, 2),                   # B(R q^2)/2
                  _scaled(bd, uqz, 2),                  # d I0(U q_Z)/L
                  _scaled(bd, u2, 32),                 # U^2 boundary
                  _delta_u2_bound(c, bd, u2, delta['value']),
                  _scaled(bd, uuz, 8))                  # U U_Z
    endpoint0, endpoint1 = _normalized_endpoint_errors(bd, dq, dqy, log_rc)
    e13 = _log_sum(bd, _scaled(bd, s0, c.mpf('.5')), endpoint0)
    # S1 includes both B(R^2 q)/L and -I1(q)/L.  Keeping both terms is
    # essential: they cancel for q0=-k/x and leave only finite increments.
    linear1 = _log_sum(bd, _scaled(bd, dq, 10), _scaled(bd, dq, 4))
    s1 = _log_sum(bd,
                  linear1,
                  _scaled(bd, dqz, 2),
                  _scaled(bd, q2, 4),
                  _scaled(bd, q2, 2),
                  _scaled(bd, uqz, 4),
                  _scaled(bd, u2, 64),
                  _scaled(bd, u2, 16),
                  _scaled(bd, uuz, 16))
    e15 = _log_sum(bd, _scaled(bd, s1, c.mpf('.5')), endpoint1)
    return dict(U=U, U_y=Uy, m=m, m_y=m_y, q_increment=q,
                q_y_increment=qy, error_m13=e13, error_m15=e15,
                endpoint_q_R_m13=endpoint0, endpoint_q_R_m15=endpoint1,
                source_scale_applied_once=True,
                finite_endpoint_q_R_retained=True,
                original_delta_upper_lt_half=True,
                inverse_L_majorant='1/L <= 2 from 0 <= delta < 1/2',
                radial_free_identity='current_original_n1_radial_free_pressure_integrals')


def _signed_record(c, log_interval, sign, error):
    """Combine a signed baseline interval with a positive error cap."""
    blo, bhi = ep(log_interval)
    if error.log is None:
        lo, hi = blo, bhi
        certified = True
    else:
        err = ep(error.log)[1]
        certified = err < blo
        if certified:
            # Do not ask mpmath to materialize exp(-10^400...).  Once the
            # logarithmic separation exceeds 100, a directed 1e-30 log slack
            # dominates the omitted correction by many orders of magnitude.
            if err < blo - 100:
                # Keep the padding inside the owner precision context.  The
                # endpoints returned by ``ep`` are ordinary mpf values, so
                # subtracting a 1e-30 slack from them directly would use
                # mpmath's global (usually 15 digit) precision and silently
                # round the padding away at the 540 digit source scale.
                # Re-enter the directed context before combining endpoints.
                slack = c.mpf('1e-30')
                blo_c = c.mpf(blo)
                bhi_c = c.mpf(bhi)
                lo = ep(blo_c - slack)[0]
                hi = ep(bhi_c + slack)[1]
            else:
                lo = c.log(c.exp(blo) - c.exp(err))
                hi = c.log(c.exp(bhi) + c.exp(err))
        else:
            lo, hi = c.ninf, c.inf
    interval = c.mpf([lo, hi]) if certified else c.mpf([c.mpf('-1e100'), c.mpf('1e100')])
    return dict(sign=sign, log_absolute_coefficient=quiet.record(interval),
                baseline_log=quiet.record(log_interval),
                perturbation_log_upper=error.record(),
                perturbation_smaller_than_baseline=certified,
                signed_value_enclosure=certified)


def baseline_pressure_coefficients(k):
    """Exact compact-band coefficients after factoring ``Rc`` powers."""
    log2 = sy.log(2) if hasattr(k, 'free_symbols') else mp.log(2)
    return dict(m13=k * k / 8, m15=-k * k * log2 / 4)


def _physical_scale_function(g, owner):
    node = g.node('current_original_repair_band_physical_source_scale',
                  source_family=owner.family, name='Pstar',
                  raw_quantity_path='outer.owner.owner(ends).Pstar',
                  source_projection='one common physical unit Pstar',
                  source_rows_provider=current.ast_binding(physical_source_scale),
                  Z_independent=True, physical_common_unit=True,
                  defining_quantity_not_a_range_value=True,
                  same_live_basis_context_and_ledger=True)
    return target.C3Function(node, g.zero, g.zero, g.zero)


def physical_source_scale(owner, ends):
    if not isinstance(owner, CurrentOriginalN1RepairBandPressureSource):
        raise ValueError('Retained repair-band pressure owner required')
    op = owner.repair.control.ranges.phase.outer.owner.owner(tuple(ends))
    return op.Pstar


def _original_log_r_minus(g, owner):
    """Restore the accepted microscopic origin before exponentiating Rc."""
    params = owner.power.control_base['parameters']
    ref = lambda name: _as_ref(g, params[name])
    return g.add(g.unary('log', g.constant(4)),
                 g.neg(g.mul(g.constant(4), ref('logP'))),
                 g.constant(-1000),
                 g.mul(ref('hbB'), ref('sc'), g.constant('1/2')))


def _absolute_log_rc(g, owner):
    """Accepted native-frame identity: log Rc = phase offset + log r_minus."""
    phase = _as_ref(g, owner.power.control_base['original_Rc_offset'])
    return g.add(phase, _original_log_r_minus(g, owner))


def _repair_band_radius(g, owner, x):
    phase = _as_ref(g, owner.power.control_base['original_Rc_offset'])
    origin = _original_log_r_minus(g, owner)
    log_rc = g.add(phase, origin)
    Rc = g.unary('exp', log_rc)
    return g.node('current_original_repair_band_physical_radius',
                  Rc=Rc.node, x=x.node, original_radius='Rc*x',
                  original_Rc_source='accepted Rp native frame: node 26685 = node 4893 + node 26684',
                  original_Rc_phase_offset_node=phase.node,
                  original_log_r_minus_source_node=origin.node,
                  original_logRc_source_node=log_rc.node,
                  original_Rc_is_absolute=True, radius_Z_independent=True,
                  source_family=owner.family,
                  defining_source=current.ast_binding(_repair_band_radius),
                  range_enclosure_is_not_selected_function_value=True)


def _outer_coordinate_substitution(g, body, owner, x, variable):
    return g.node('substitute_repair_band_coordinate_at_pressure_native',
                  body=body.node, free_query=x.node, native_variable=variable,
                  free_query_is_outer_pressure_endpoint=True,
                  inner_repair_band_native_variable='current_C3_repair_band_x',
                  inner_Duhamel_native_variables_preserved=True,
                  source_family=owner.family,
                  defining_source=current.ast_binding(_outer_coordinate_substitution),
                  same_original_source_and_global_phase=True)


def _signed_integrals(owner, fields):
    """Build signed value/Z1 integral handles over x in [1,2]."""
    g = owner.g
    x = _as_ref(g, owner.repair.functions['original_band_variable'])
    variable = 'repair_band_pressure_native'
    one, two = g.constant(1), g.constant(2)
    radius = _repair_band_radius(g, owner, x)
    Rc = g.unary('exp', _absolute_log_rc(g, owner))
    out = {}
    for name, value in fields['integrands'].items():
        # I0 fields integrate against dR=Rc*dx; I1 fields integrate against
        # R dR=(Rc*x)*Rc*dx.  The field names already encode the desired
        # moment and must not receive an extra radius factor.
        power_r = g.one if name.startswith('I0_') else radius
        rows = []
        for order, row in enumerate(physical.z1_rows(value)):
            body = g.mul(power_r, g.mul(_as_ref(g, row), Rc))
            body = _outer_coordinate_substitution(g, body, owner, x, variable)
            # The compact band has a fixed exact x union.  Keep this as a
            # source-bound graph quantity rather than expanding the enormous
            # generic integral AST; the radial-free identity is the numerical
            # evaluator and this node preserves its real body and measure.
            integral = g.node(
                'current_original_repair_band_pressure_signed_integral',
                body=body.node, free_query=x.node, native_variable=variable,
                lower=one.node, upper=two.node,
                measure=('physical dR=Rc*dx exactly once' if name.startswith('I0_')
                         else 'physical R*dR=R*Rc*dx exactly once'),
                ordinary_slow_Z_derivative_order=order,
                free_query_substituted_at_native=True,
                original_physical_radius_source='Rc*original_band_x',
                original_physical_Jacobian_applied_once=True,
                endpoint_independent_of_Z=True,
                compact_repair_band_only=True,
                defining_source=current.ast_binding(_signed_integrals),
                radial_free_identity='current_original_n1_radial_free_pressure_integrals')
            rows.append(integral)
        out[name] = physical.Z1Function(*rows)
    return out


@dataclass(frozen=True, eq=False)
class OriginalN1RepairBandPressurePacket:
    leading_hierarchy_order: int = 0


class CurrentOriginalN1RepairBandPressureSource:
    @_source_precision_540
    def __init__(self, checked_repair_owner=None, checked_power_owner=None,
                 checked_quiet_owner=None, require_checked=True):
        self.repair = checked_repair_owner if checked_repair_owner is not None else repair.CurrentC3RepairBand(require_checked=True)
        self.power = checked_power_owner if checked_power_owner is not None else power.CurrentC3PowerToRp(require_checked=True)
        self.quiet = checked_quiet_owner if checked_quiet_owner is not None else quiet.CurrentOriginalN1QuietPressureSource(require_checked=True)
        if not (self.repair.acceptance_loaded and self.power.acceptance_loaded and self.quiet.acceptance_loaded):
            raise ValueError('Accepted repair-band, power-to-Rp, and quiet pressure sources required')
        if self.repair.identity != self.power.identity or self.repair.identity != self.quiet.family:
            raise ValueError('Same original source family required')
        self.family = copy.deepcopy(self.repair.identity)
        self.native_radius_binding = _accepted_native_radius_binding(self.family)
        self.g = self.power.graph
        self.c = self.power.c
        self.logC = current.packets.interval(
            self.c, self.native_radius_binding['selected_logCstar'])
        # The accepted native frame proves
        #   log Rc = log(110) + 10 logCstar + 11 exp(40) + 124.
        # Keep this as a logarithmic source interval; never exponentiate it.
        self.log_Rc = (self.c.ln(110) + 10 * self.logC +
                       11 * self.c.exp(40) + 124)
        logP = _as_ref(self.g, self.power.control_base['parameters']['logP'])
        self.delta_node = self.g.unary(
            'exp', self.g.neg(self.g.add(
                self.g.mul(self.g.constant(4), logP),
                self.g.constant(30))))
        self.delta = self.c.exp(-4 * (self.c.exp(40) + 11) - 30)
        self.N = self.power.control_base['N']
        self.hashes = dict(self.repair.hashes)
        self.hashes.update(self.power.hashes)
        self.hashes.update(self.quiet.hashes)
        for name in (repair.NAME, repair.RECEIPT, power.NAME, power.RECEIPT,
                     quiet.NAME, quiet.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        for name, digest in self.native_radius_binding['receipt']['input_hashes'].items():
            if name in self.hashes and self.hashes[name] != digest:
                raise ValueError('Conflicting native radius dependency ' + name)
            self.hashes[name] = digest
        for name in (NATIVE_FRAME_NAME, NATIVE_FRAME_RECEIPT, NORM_NAME):
            self.hashes[name] = sha(name)
        self._quiet_raw = _read_json_gz(quiet.HERE / quiet.NAME)
        krec = self._quiet_raw['actual_quiet_source']['K']['log_coefficient']
        self.k_log = quiet.read_interval(self.c, krec)
        self.k = rg.LogUpper(self.c, self.k_log)
        self.functions = self._build_functions()
        self._raw_repair = _read_json_gz(repair.HERE / repair.NAME)
        self.rows = {_cell_key(row['exact_Z_cell']): row
                     for row in self._raw_repair['actual_four_Z_C3_band_ranges']}
        self._packets = {}
        self.acceptance_loaded = False
        if require_checked:
            receipt = json.loads((HERE / RECEIPT).read_bytes())
            if (not receipt.get('all_passed') or not receipt.get(GATE) or
                    receipt.get('source_family') != self.family or
                    any(receipt.get(key, False) for key in OPEN)):
                raise ValueError('Checked repair-band pressure source required')
            for name, digest in receipt['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('Changed repair-band pressure source ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._graph_length = len(self.g.nodes)
        self._delta = self.delta._mpi_
        self._hashes = copy.deepcopy(self.hashes)

    def _build_functions(self):
        g = self.g
        alg = target.C3Algebra(g)
        # Only the actual band rows used by pressure are copied.  Recursing
        # through the archived C3 proof dictionary would duplicate tens of
        # thousands of source nodes and is unnecessary for this adapter.
        source_functions = self.repair.functions
        profile = [{key: _clone_c3(g, value)
                    for key, value in row.items()}
                   for row in source_functions['profiles_at_band_x']['radial_y_rows']]
        complete_m = _clone_c3(g, source_functions['complete_histories']['m'])
        complete_y_m = _clone_c3(g, source_functions['complete_history_y']['m'])
        original_x = _as_ref(g, source_functions['original_band_variable'])
        U = alg.mul(_physical_scale_function(g, self), profile[0]['V'])
        m = alg.mul(_physical_scale_function(g, self), complete_m)
        q = partial.recover_Q_increment(g, U, m, self.delta_node)
        Uy = alg.mul(_physical_scale_function(g, self), profile[1]['V'])
        my = alg.mul(_physical_scale_function(g, self), complete_y_m)
        qy = partial.recover_Q_increment(g, Uy, my, self.delta_node)
        integrands = physical.nonradial_integrands(g, U, q)
        defining = _signed_integrals(self, dict(U=U, m=m, q=q, q_y=qy,
                                               integrands=integrands))
        return dict(U=U, m=m, m_y=alg.add(U, alg.neg(m)), q=q,
                    U_y=Uy, m_y_source=my, q_y=qy, integrands=integrands,
                    defining_integrals=defining,
                    source_band_functions=dict(profiles_at_band_x=profile,
                                               complete_history_m=complete_m,
                                               complete_history_y_m=complete_y_m),
                    original_band_x=original_x,
                    pressure_native_variable='repair_band_pressure_native',
                    physical_radius_scale='Rc',
                    finite_endpoint_q_R_source='q_y/R from actual band y rows')

    def assert_graph(self):
        if self.hashes != self._hashes or self.delta._mpi_ != self._delta:
            raise ValueError('Repair-band pressure source state changed')
        if not self.repair.acceptance_loaded or not self.power.acceptance_loaded or not self.quiet.acceptance_loaded:
            raise ValueError('Accepted source owners must remain live')

    @_source_precision_540
    def evaluate(self, ends):
        self.assert_graph()
        key = _cell_key(ends)
        if key not in self.rows:
            raise ValueError('Accepted original axial cell required')
        c = self.c
        op = self.repair.control.ranges.phase.outer.owner.owner(tuple(ends))
        pstar = rg.LogUpper(c, op.Pstar.record()['log_absolute_upper'])
        raw = self.rows[key]
        perturb = _finite_perturbation(
            c, tuple(ends), raw, dict(value=self.delta, k=self.k), pstar,
            self.log_Rc)
        base13 = 2 * self.k_log - c.ln(8)
        base15 = 2 * self.k_log + c.ln(c.ln(2) / 4)
        m13 = _signed_record(c, base13, 'positive', perturb['error_m13'])
        m15 = _signed_record(c, base15, 'negative', perturb['error_m15'])
        view = dict(
            source_family=self.family, exact_Z_cell=tuple(ends), band_R_over_Rc=[1, 2],
            signed_pressure_m13=m13, signed_pressure_m15=m15,
            normalized_baseline=dict(m13_over_Rc='k^2/8', m15_over_Rc2='-k^2*log(2)/4',
                                     k=quiet.record(self.k_log)),
            absolute_Rc_source=dict(native_frame_nodes=self.native_radius_binding['expected_nodes'],
                                    log_Rc=quiet.record(self.log_Rc),
                                    log_Rc_formula='log(110)+10*selected_logCstar+11*exp(40)+124',
                                    radius_not_materialized=True),
            perturbation_bounds=dict(error_m13=perturb['error_m13'].record(),
                                     error_m15=perturb['error_m15'].record(),
                                     endpoint_q_R_m13=perturb['endpoint_q_R_m13'].record(),
                                     endpoint_q_R_m15=perturb['endpoint_q_R_m15'].record()),
            finite_N_velocity_and_mean_retained=True,
            original_delta_upper_lt_half=True,
            inverse_L_majorant='1/L <= 2 from 0 <= delta < 1/2',
            physical_common_scale_Pstar_applied_once=True,
            original_Rc_scale_retained=True,
            radial_free_pressure_identity='current_original_n1_radial_free_pressure_integrals',
            finite_endpoint_q_R_retained=True,
            signed_functions_source='functions.U,m,q,q_y,defining_integrals',
            numerical_evidence_signed_value_enclosure=True,
            no_radius_or_N_materialization=True,
            no_outer_support_or_full_tail_claim=True,
            source_raw_parent_scope='accepted CurrentC3RepairBand actual_four_Z_C3_band_ranges',
            **{GATE: self.acceptance_loaded}, **dict.fromkeys(OPEN, False))
        packet = OriginalN1RepairBandPressurePacket()
        self._packets[id(packet)] = (packet, view, perturb)
        return packet

    def evaluate_all(self):
        return {key: self.evaluate(key) for key in self.rows}

    @_source_precision_540
    def report(self, packet):
        self.assert_graph()
        item = self._packets.get(id(packet))
        if item is None or item[0] is not packet:
            raise ValueError('Live repair-band pressure packet required')
        return copy.deepcopy(item[1])


def run(checked_repair_owner=None, checked_power_owner=None,
        checked_quiet_owner=None):
    began = time.monotonic()
    with mp.workdps(540):
        owner = CurrentOriginalN1RepairBandPressureSource(
            checked_repair_owner, checked_power_owner, checked_quiet_owner,
            require_checked=False)
        packets = owner.evaluate_all()
        # JSON object keys must be scalar.  Keep the exact cell in each value
        # while using a stable string key for the archived source packet.
        samples = {
            '|'.join(key): dict(exact_Z_cell=list(key), view=owner.report(packet))
            for key, packet in packets.items()
        }
        raw = dict(candidate_actual_repair_band_pressure_source_constructed=True,
                   source_family=owner.family, exact_band=[1, 2],
                   exact_Z_cells=[list(key) for key in owner.rows], samples=samples,
                   signed_function_source_nodes={name: _record(value)
                                                 for name, value in owner.functions.items()
                                                 if name in ('U', 'm', 'q', 'q_y')},
                   source_identity=dict(repair_receipt=repair.RECEIPT,
                                        power_receipt=power.RECEIPT,
                                        quiet_receipt=quiet.RECEIPT,
                                        same_original_family=True),
                   input_hashes=owner.hashes, **{GATE: False},
                   **dict.fromkeys(OPEN, False),
                   execution_seconds=time.monotonic() - began)
        (HERE / NAME).write_bytes(gzip.compress(
            (json.dumps(repair.control.encoded(raw), separators=(',', ':')) + '\n').encode(),
            mtime=0))
        print('REPAIR_BAND_SIGNED_PRESSURE_SOURCE_COMPUTED', len(packets), flush=True)
        return owner, packets


if __name__ == '__main__':
    run()
