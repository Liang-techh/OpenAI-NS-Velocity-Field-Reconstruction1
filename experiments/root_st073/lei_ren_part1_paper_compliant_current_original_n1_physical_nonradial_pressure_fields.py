"""Physical leading fields and nonradial pressure moment integrands.

This is a deliberately narrow continuation of the checked finite-N partial
mean/Q owner.  It retains the original signed leading V source, binds the
original own m history, and multiplies both by the one live common unit
``S=Pstar``.  The finite-N increment is the retained partial-mean graph, so
the physical fields are

    U = S (V + delta_V),  m = S (m_0 + delta_m),

with ``m_y=U-m`` and

    q = (2*Z*U-(1-delta)*Z*m-(1-Z**2)*m_Z)/(1-delta*Z**2).

Only first_micro and second_micro whole-cell magnitude bounds are exposed.
Signed graph functions remain source expressions; magnitude rows are never
used as point values.  This module stops before support, finite-endpoint
q_R, repair, higher n, or global pressure closure.
"""
import copy
import gzip
import json
import math
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q as partial


c3, current, source = partial.c3, partial.current, partial.source
rg = partial.rg
HERE, PREFIX, sha, ep = partial.HERE, partial.PREFIX, partial.sha, partial.ep
source_precision = partial.source_precision
fingerprint = partial.fingerprint
OPEN = tuple(dict.fromkeys(partial.OPEN))

NAME = PREFIX + 'current_original_n1_physical_nonradial_pressure_fields.json.gz'
RECEIPT = PREFIX + 'current_original_n1_physical_nonradial_pressure_fields_check.json'
GATE = 'current_original_n1_physical_nonradial_pressure_fields_installed'

CHARTS = ('first_micro', 'second_micro')
LIMITS = {'first_micro': (Fraction(1, 4), Fraction(1)),
          'second_micro': (Fraction(1), Fraction(2))}
I0_NAMES = ('q_Z', 'U_q_Z', 'U_squared', 'U_U_Z')
I1_NAMES = ('q', 'q_Z', 'q_squared', 'U_q_Z', 'U_squared', 'U_U_Z')


def rows(value):
    return (value.value, value.Z, value.ZZ, value.ZZZ)


@dataclass(frozen=True)
class Z1Function:
    """A signed function with exactly the two rows this pressure route uses."""

    value: object
    Z: object


@dataclass(frozen=True)
class Z1Bound:
    """A magnitude product bound with value and Z1 rows only."""

    value: object
    Z: object


def z1_rows(value):
    return (value.value, value.Z)


def _z1_encode(value):
    if isinstance(value, Z1Function):
        return [value.value.node, value.Z.node]
    if isinstance(value, Z1Bound):
        return [value.value.record(), value.Z.record()]
    if isinstance(value, dict):
        return {key if isinstance(key, str) else repr(key): _z1_encode(item)
                for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_z1_encode(item) for item in value]
    return c3.encoded(value)


def _bound_record(value):
    if isinstance(value, Z1Bound):
        return [value.value.record(), value.Z.record()]
    if isinstance(value, rg.JetBound):
        return rg.record(value)
    if isinstance(value, dict):
        return {key: _bound_record(item) for key, item in value.items()
                if key not in ('source_packet', 'source_family')}
    if isinstance(value, (list, tuple)):
        return [_bound_record(item) for item in value]
    if hasattr(value, 'record'):
        return value.record()
    return value


def _as_ref(g, value):
    return value if isinstance(value, source.FunctionRef) else source.FunctionRef(g, value)


def _equivalent_rows(actual, expected):
    return partial.current.current.previous.equivalent_rows(actual, expected)


def _same_source(value, c, flow, label):
    """Guard one live ordinary source object without evaluating its scale."""
    if value.ctx is not c or value.scale.bases is not flow.logs or value.ledger is not flow.ledger:
        raise ValueError('Same live source context/basis/ledger required for ' + label)


def source_rows(owner, ends, chart, left, right=None):
    """Read the actual unprojected leading source rows and physical radius.

    The row arrays are Taylor coefficients.  Consumers apply factorial(k)
    once when making ordinary Z derivatives.  No row endpoint is selected as
    a function value here.
    """
    if type(owner) is not CurrentOriginalN1PhysicalNonradialPressureFields:
        raise ValueError('This source provider requires the retained physical-field owner')
    ends = tuple(ends)
    if chart not in CHARTS:
        raise ValueError('Physical source rows are limited to the two actual micro charts')
    packet = owner.parent.parent.target.phase.leading_source(ends, chart, left, right)
    if (packet.get('source_identity') != owner.family or
            packet.get('exact_common_P0_axial5') is not
            owner.parent.parent.ranges.phase.outer.owner.owner(ends).P0 or
            packet.get('actual_phase_Z_exact_zero') is not True):
        raise ValueError('Same physical source family/P0 and Z-independent phase required')
    generic = packet.get('original_generic_source')
    if not isinstance(generic, dict):
        raise ValueError('Original generic source packet missing')
    required = ('common_velocity_V_axial5', 'common_own_five_histories_axial5',
                'original_physical_radius')
    if any(name not in generic for name in required):
        raise ValueError('Required original physical source quantity absent')
    histories = generic['common_own_five_histories_axial5']
    if 'm' not in histories or len(generic['common_velocity_V_axial5']) < 4 or len(histories['m']) < 4:
        raise ValueError('Original V/m Taylor rows through order 3 required')
    op = owner.parent.parent.ranges.phase.outer.owner.owner(ends)
    flow = op.flow
    if (len(op.P0) < 4 or
            not _equivalent_rows(generic.get('common_original_P0_axial5', ()), op.P0) or
            any(value.ctx is not owner.c or value.scale.bases is not flow.logs or
                value.ledger is not flow.ledger for value in op.P0)):
        raise ValueError('Equivalent original independent P0 Taylor rows required')
    for k, value in enumerate(generic['common_velocity_V_axial5']):
        _same_source(value, owner.c, flow, 'V Taylor row %d' % k)
    for k, value in enumerate(histories['m']):
        _same_source(value, owner.c, flow, 'm Taylor row %d' % k)
    _same_source(generic['original_physical_radius'], owner.c, flow, 'original physical radius')
    _same_source(op.Pstar, owner.c, flow, 'Pstar')
    return dict(V=generic['common_velocity_V_axial5'],
                m=histories['m'],
                R=generic['original_physical_radius'],
                Pstar=op.Pstar,
                source_packet=packet,
                exact_common_P0=op.P0,
                source_family=owner.family)


def physical_source_scale(owner, ends):
    """Return the actual common physical unit, without materializing it."""
    if type(owner) is not CurrentOriginalN1PhysicalNonradialPressureFields:
        raise ValueError('Physical scale provider requires retained owner')
    op = owner.parent.parent.ranges.phase.outer.owner.owner(tuple(ends))
    _same_source(op.Pstar, owner.c, op.flow, 'Pstar')
    return op.Pstar


def _query_substitution(g, body, chart, query):
    """Substitute the free physical query by the native integration variable."""
    pressure_variable = 'pressure_moment_native_' + chart
    return g.node('substitute_free_query_at_native_coordinate', body=body.node,
                   free_query=query.node, native_variable=pressure_variable,
                   bound_partial_mean_native_variable='native_' + chart,
                   free_query_is_outer_pressure_endpoint=True,
                   defining_source=current.ast_binding(_query_substitution),
                   source_query_substitution='free physical query -> pressure moment native coordinate',
                   inner_Duhamel_native_variable_preserved=True,
                   same_original_source_and_global_phase=True)


def _physical_radius(g, log_radius, chart, query):
    """Bind R to the original physical log-radius geometry."""
    return g.node('current_original_physical_radius', log_radius=log_radius.node,
                  free_query=query.node, native_chart=chart,
                  raw_quantity_path='original_generic_source.original_physical_radius',
                  radius_expression='exp(original physical log-radius)',
                  original_log_radius_is_absolute=True,
                  radius_Z_independent=True,
                  defining_source=current.ast_binding(_physical_radius),
                  range_enclosure_is_not_selected_function_value=True)


def _at_query(g, value, chart, query):
    """Evaluate an original native source body at the free physical query."""
    return c3.C3Function(*(partial.at_native(g, _as_ref(g, row), chart, query)
                           for row in rows(value)))


def _m_source_function(g, owner, chart):
    """Bind m rows to the actual retained generic-source history.

    The graph nodes carry explicit paths and the real source-row provider.
    Their ordinary derivative convention records the sole factorial(k)
    conversion from source Taylor coefficients.
    """
    template = owner.parent.parent.target.primitives[chart]['source_roots']['V']
    template_attrs = dict(g.nodes[template.value.node])
    coordinate = template_attrs.get('coordinate')
    Z_variable = template_attrs.get('Z_variable')
    if Z_variable is None:
        raise ValueError('Original leading source Z variable metadata missing')
    out = []
    for k in range(4):
        attrs = dict(source_family=owner.family, native_chart=chart,
                     quantity='m', Z_order=k, coordinate=coordinate,
                     Z_variable=Z_variable,
                     phase_independent_leading_source=True,
                     defining_quantity_not_a_range_value=True,
                     ordinary_slow_Z_derivative_order=k,
                     Taylor_coefficient_factorial=math.factorial(k),
                     source_projection='factorial(%d)*original_generic_source.common_own_five_histories_axial5["m"][%d]' % (k, k),
                     raw_quantity_path='original_generic_source.common_own_five_histories_axial5["m"][%d]' % k,
                     source_rows_provider=current.ast_binding(source_rows),
                     source_projection_binding=current.ast_binding(source_rows),
                     source_row_is_Taylor_coefficient=True,
                     factorial_conversion_applied_once=True,
                     original_leading_source_recipe=dict(g.nodes[template.value.node]).get('recipe'),
                     same_live_basis_context_and_ledger=True,
                     range_enclosure_is_not_selected_function_value=True)
        out.append(g.node('current_original_physical_source_quantity', **attrs))
    return c3.C3Function(*out)


def _physical_scale_function(g, owner):
    node = g.node('current_original_physical_source_scale',
                  source_family=owner.family, name='Pstar',
                  raw_quantity_path='outer.owner.owner(ends).Pstar',
                  source_rows_provider=current.ast_binding(physical_source_scale),
                  source_projection_binding=current.ast_binding(physical_source_scale),
                  Z_independent=True, physical_common_unit=True,
                  defining_quantity_not_a_range_value=True,
                  same_live_basis_context_and_ledger=True)
    return c3.C3Function(node, g.zero, g.zero, g.zero)


def _z1_shift(value):
    return Z1Function(value.Z, value.ZZ)


def _z1_mul(g, left, right):
    return Z1Function(g.mul(left.value, right.value),
                      g.add(g.mul(left.Z, right.value),
                            g.mul(left.value, right.Z)))


def _z1_add(g, *values):
    return Z1Function(g.add(*(value.value for value in values)),
                      g.add(*(value.Z for value in values)))


def nonradial_integrands(g, U, q):
    """The ten physical integrands, exposing genuine value/Z1 rows only."""
    u = Z1Function(U.value, U.Z)
    uz, qz = _z1_shift(U), _z1_shift(q)
    qr = Z1Function(q.value, q.Z)
    uqz, uu, uuz = (_z1_mul(g, u, qz), _z1_mul(g, u, u),
                    _z1_mul(g, u, uz))
    return dict(I0_q_Z=qz, I0_U_q_Z=uqz, I0_U_squared=uu, I0_U_U_Z=uuz,
                I1_q=qr, I1_q_Z=qz, I1_q_squared=_z1_mul(g, qr, qr),
                I1_U_q_Z=uqz, I1_U_squared=uu, I1_U_U_Z=uuz)


def _interval_key(left, right):
    def canonical(value):
        try:
            exact = rational(tuple(value))
            return [exact.numerator, exact.denominator]
        except ValueError:
            return list(value)
    return json.dumps([canonical(left), canonical(right)], separators=(',', ':'))


def _cell_key(cell):
    return _interval_key(cell['exact_left'], cell['exact_right'])


def _query_cells(base, left, right):
    target_left, target_right = rational(left), rational(right)
    selected = []
    for cell in base['actual_current_source_cells']:
        try:
            cell_left, cell_right = rational(tuple(cell['exact_left'])), rational(tuple(cell['exact_right']))
        except ValueError:
            continue
        if target_left <= cell_left and cell_right <= target_right:
            selected.append((cell_left, cell_right, cell))
    selected.sort(key=lambda item: item[0])
    if (not selected or selected[0][0] != target_left or
            selected[-1][1] != target_right or
            any(a[1] != b[0] for a, b in zip(selected, selected[1:]))):
        raise ValueError('Requested physical pressure cell is not an exact source-cell union')
    return [cell for _, _, cell in selected]


def _selected_defining_integrals(owner, chart, left, right):
    functions = owner.functions[chart]
    key = _interval_key(left, right)
    result = {}
    for name in functions['defining_integrals_by_cell']:
        try:
            result[name] = functions['defining_integrals_by_cell'][name][key]
        except KeyError as exc:
            raise ValueError('Requested physical pressure cell is not an exact source-cell union') from exc
    return result


def _original_log_r_minus(g, owner):
    """Restore the source microscopic origin before exponentiating R."""
    params = owner.parent.parent.target.phase.built['parameters']
    return g.add(g.unary('log', g.constant(4)),
                 g.neg(g.mul(g.constant(4), params['logP'])),
                 g.constant(-1000),
                 g.mul(params['hbB'], params['sc'], g.constant('1/2')))


def _field_functions(owner, chart):
    g = owner.g
    alg = c3.C3Algebra(g)
    query = owner.parent.functions[chart]['query']
    base = owner.parent.parent.target.phase.windows[chart]
    leading = owner.parent.parent.target.primitives[chart]['source_roots']
    scale = _physical_scale_function(g, owner)
    baseline_V = _at_query(g, leading['V'], chart, query)
    baseline_m = _at_query(g, _m_source_function(g, owner, chart), chart, query)
    delta_V = owner.parent.functions[chart]['delta_V']
    delta_m = owner.parent.functions[chart]['mean']
    U = alg.mul(scale, alg.add(baseline_V, delta_V))
    m = alg.mul(scale, alg.add(baseline_m, delta_m))
    q = partial.recover_Q_increment(g, U, m, owner.delta_node)
    q_Z = _z1_shift(q)
    U_Z = _z1_shift(U)
    m_y = alg.add(U, alg.neg(m))
    integrands = nonradial_integrands(g, U, q)
    offset = _as_ref(g, base['offset'])
    relative_log_radius = partial.at_native(g, offset, chart, query)
    absolute_log_radius = g.add(_original_log_r_minus(g, owner), relative_log_radius)
    radius = _physical_radius(g, absolute_log_radius, chart, query)
    jacobian = _as_ref(g, base['Jacobian'])
    defining_integrals_by_cell = {name: {} for name in integrands}
    for name, value in integrands.items():
        power = 1 if name.startswith('I0_') else 2
        # Build one defining integral per retained original native source cell
        # and per ordinary Z row.  The free query is substituted at evaluation.
        for cell in base['actual_current_source_cells']:
            lower = _as_ref(g, cell['lower'])
            upper = _as_ref(g, cell['upper'])
            rows_integral = []
            for j, row in enumerate(z1_rows(value)):
                body_q = g.mul(radius if power == 1 else g.mul(radius, radius),
                               g.mul(_as_ref(g, row), jacobian))
                body_native = _query_substitution(g, body_q, chart, query)
                rows_integral.append(current.controls.integral(
                    g, body_native, 'pressure_moment_native_' + chart, lower, upper,
                    measure=('physical dR=R*Jacobian*dnative exactly once' if power == 1
                             else 'physical R*dR=R^2*Jacobian*dnative exactly once'),
                     ordinary_slow_Z_derivative_order=j,
                    free_query_endpoint=query.node,
                    free_query_substituted_at_native=True,
                    original_physical_radius_source='original_generic_source.original_physical_radius',
                    original_physical_Jacobian_applied_once=True,
                    endpoint_independent_of_Z=True))
            defining_integrals_by_cell[name][_cell_key(cell)] = Z1Function(*rows_integral)
    ordered_cells = []
    for cell in base['actual_current_source_cells']:
        try:
            rational(tuple(cell['exact_left']))
            rational(tuple(cell['exact_right']))
        except ValueError:
            continue
        ordered_cells.append(cell)
    ordered_cells.sort(key=lambda cell: rational(tuple(cell['exact_left'])))
    for name, cell_integrals in defining_integrals_by_cell.items():
        unions = {}
        for first in range(len(ordered_cells)):
            total = None
            for last in range(first, len(ordered_cells)):
                cell = ordered_cells[last]
                total = (cell_integrals[_cell_key(cell)] if total is None
                         else _z1_add(g, total, cell_integrals[_cell_key(cell)]))
                unions[_interval_key(ordered_cells[first]['exact_left'],
                                     cell['exact_right'])] = total
        defining_integrals_by_cell[name] = unions
    return dict(query=query,scale=scale,baseline_V=baseline_V,baseline_m=baseline_m,
                delta_V=delta_V,delta_m=delta_m,U=U,m=m,m_y=m_y,q=q,q_Z=q_Z,U_Z=U_Z,
                radius=radius,integrands=integrands,
                defining_integrals_by_cell=defining_integrals_by_cell,
                source_geometry=dict(raw_quantity_path='original_generic_source.original_physical_radius',
                                     original_log_r_minus='log(4)-4*logP-1000+hbB*sc/2',
                                     relative_log_radius_offset=base['offset'],
                                     Jacobian=base['Jacobian'],
                                     pressure_moment_native_variable='pressure_moment_native_' + chart,
                                     bound_partial_mean_native_variable='native_' + chart,
                                     physical_measure='R^power*Jacobian*dnative with power 1 or 2'))


def _raw_ordinary_jet(c, bd, values, label):
    if len(values) < 4:
        raise ValueError('Four source Taylor rows required for ' + label)
    return rg.JetBound(*(bd.row(values[k] * math.factorial(k)) for k in range(4)))


def _physical_q_bounds(c, ends, delta, U, m):
    """Bound q's ordinary rows through Z2 from physical U/m rows."""
    z = IntervalTaylor.variable(c, c.mpf(ends), 2)
    L = 1 - z * z * delta
    d = 1 - z * z
    coefficients = (2 * z / L, -z * (1 - delta) / L, -d / L)
    bd = rg.Bounds(c)
    ur, mr = rg.rows(U), rg.rows(m)
    out = []
    for k in range(3):
        pieces = []
        for j in range(k + 1):
            for coefficient, row in zip(coefficients, (ur, mr, mr[1:])):
                magnitude = c.mpf(ep(abs(coefficient[j]) * math.factorial(j))[1])
                pieces.append(bd.constant(math.comb(k, j) * magnitude) * row[k - j])
        out.append(bd.sum(*pieces))
    return out


def _bound_z1_product(bd, left, right):
    return Z1Bound(left.value * right.value,
                   bd.sum(left.Z * right.value, left.value * right.Z))


def _bound_z1_shift(value):
    return Z1Bound(value.Z, value.ZZ)


def _bound_z1_scale(bd, value, factor):
    return Z1Bound(value.value * factor, value.Z * factor)


def _source_and_physical_bounds(owner, ends, chart, left, right, partial_packet):
    c = owner.c
    bd = rg.Bounds(c)
    source_view = source_rows(owner, ends, chart, left, right)
    op = owner.parent.parent.ranges.phase.outer.owner.owner(tuple(ends))
    f = op.flow
    if source_view['Pstar'] is not op.Pstar:
        raise ValueError('Pstar source object changed')
    # The checked parent retains the local finite-N state; no cap is treated
    # as a signed field value.
    entry = owner.parent._packets.get(id(partial_packet))
    if entry is None or entry[0] is not partial_packet:
        raise ValueError('Retained live partial mean/Q packet required')
    state = entry[2]
    partial_mean, partial_mean_y, partial_Q, scaled_mean, delta_V = state
    if not (isinstance(partial_mean, rg.JetBound) and isinstance(delta_V, rg.JetBound)):
        raise ValueError('Parent partial state must retain JetBound rows')
    baseline_V = _raw_ordinary_jet(c, bd, source_view['V'], 'baseline V')
    baseline_m = _raw_ordinary_jet(c, bd, source_view['m'], 'baseline m')
    S = bd.row(source_view['Pstar'])
    U = bd.scaled(bd.add(baseline_V, delta_V), S)
    m = bd.scaled(bd.add(baseline_m, partial_mean), S)
    q = _physical_q_bounds(c, ends, owner.delta, U, m)
    qjet = Z1Bound(q[0], q[1])
    qz = Z1Bound(q[1], q[2])
    uz = _bound_z1_shift(U)
    u1 = Z1Bound(U.value, U.Z)
    products = {
        'I0_q_Z': qz,
        'I0_U_q_Z': _bound_z1_product(bd, u1, qz),
        'I0_U_squared': _bound_z1_product(bd, u1, u1),
        'I0_U_U_Z': _bound_z1_product(bd, u1, uz),
        'I1_q': qjet,
        'I1_q_Z': qz,
        'I1_q_squared': _bound_z1_product(bd, qjet, qjet),
        'I1_U_q_Z': _bound_z1_product(bd, u1, qz),
        'I1_U_squared': _bound_z1_product(bd, u1, u1),
        'I1_U_U_Z': _bound_z1_product(bd, u1, uz),
    }
    weight = owner.parent.parent.ranges.weights(tuple(ends), chart, left, right, Fraction(0))
    width = bd.row(weight['original_mass'])
    R = bd.row(source_view['R'])
    measure0 = R * width
    measure1 = R * R * width
    integrated = {}
    for name, value in products.items():
        integrated[name] = _bound_z1_scale(
            bd, value, measure0 if name.startswith('I0_') else measure1)
    return dict(source_view=source_view,partial_state=state,baseline_V=baseline_V,
                baseline_m=baseline_m,S=S,U=U,m=m,q=q,qjet=qjet,q_Z=qz,U_Z=uz,
                original_mass=width,original_R=R,physical_dR_measure=measure0,
                physical_R_dR_measure=measure1,products=products,integrated=integrated,
                no_extra_Jacobian=True,ordinary_baseline_factorials_once=True,
                physical_common_scale_applied_once=True)


def rational(value):
    if (not isinstance(value, tuple) or len(value) != 2 or
            any(type(v) is not int for v in value) or value[1] <= 0):
        raise ValueError('Exact rational native tuple required')
    return Fraction(*value)


@dataclass(frozen=True, eq=False)
class OriginalPhysicalNonradialPressurePacket:
    leading_hierarchy_order: int = 0


class CurrentOriginalN1PhysicalNonradialPressureFields:
    """Checked physical fields composed on one retained partial-mean owner."""

    @source_precision
    def __init__(self, accepted_partial_owner, require_checked=True):
        if (type(accepted_partial_owner) is not partial.CurrentOriginalN1PartialMeanQ or
                not accepted_partial_owner.acceptance_loaded):
            raise ValueError('Accepted checked CurrentOriginalN1PartialMeanQ owner required')
        self.parent = accepted_partial_owner
        self.parent.assert_graph()
        self.g = accepted_partial_owner.g
        self.c = accepted_partial_owner.c
        self.family = copy.deepcopy(accepted_partial_owner.family)
        self.hashes = dict(accepted_partial_owner.hashes)
        self.N = accepted_partial_owner.N
        self.delta = accepted_partial_owner.delta
        self.delta_node = accepted_partial_owner.delta_node
        for name in (partial.NAME, partial.RECEIPT, Path(partial.__file__).name,
                     Path(__file__).name):
            self.hashes[name] = sha(name)
        self._delta = self.delta._mpi_
        # Bind only the two requested numerical charts.  Existing source
        # functions remain the original live graph prefix.
        self.functions = {chart: _field_functions(self, chart) for chart in CHARTS}
        self._identities = (self.parent, self.g, self.c, self.N, self.delta_node)
        self._graph_length = len(self.g.nodes)
        self._graph = fingerprint(self.g.nodes)
        self._functions = fingerprint(_z1_encode(self.functions))
        self._family = copy.deepcopy(self.family)
        self._hashes = copy.deepcopy(self.hashes)
        self._packets = {}
        self.acceptance_loaded = False
        if require_checked:
            checked = json.loads((HERE / RECEIPT).read_bytes())
            if (not checked.get('all_passed') or not checked.get(GATE) or
                    checked.get('source_family') != self.family or
                    any(checked.get(k, False) for k in OPEN)):
                raise ValueError('Physical nonradial pressure field receipt/scope differs')
            for name, digest in checked['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('Physical field source changed ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if ((self.parent, self.g, self.c, self.N, self.delta_node) != self._identities or
                self.g is not self.parent.g or not self.parent.acceptance_loaded or
                self.family != self._family or self.family != self.parent.family or
                self.hashes != self._hashes or len(self.g.nodes) < self._graph_length or
                fingerprint(self.g.nodes[:self._graph_length]) != self._graph or
                fingerprint(_z1_encode(self.functions)) != self._functions or
                self.delta._mpi_ != self._delta):
            raise ValueError('Physical nonradial pressure graph changed')
        self.parent.assert_graph()

    @source_precision
    def range_cell(self, ends, chart, left, right):
        self.assert_graph()
        began = time.monotonic()
        ends = tuple(ends)
        l, r = rational(left), rational(right)
        if chart not in CHARTS or not (LIMITS[chart][0] <= l < r <= LIMITS[chart][1]):
            raise ValueError('Supported ordered actual micro cell required')
        partial_packet = self.parent.range_cell(ends, chart, left, right)
        bounds = _source_and_physical_bounds(self, ends, chart, left, right, partial_packet)
        functions = self.functions[chart]
        defining_integrals = _selected_defining_integrals(self, chart, left, right)
        view = dict(
            source_family=self.family, exact_Z_cell=ends, chart=chart,
            exact_native_cell=[left, right],
            signed_physical_field_functions=_z1_encode(functions),
            signed_nonradial_integrands=_z1_encode(functions['integrands']),
            signed_defining_integrals=_z1_encode(defining_integrals),
            original_source_quantity_paths=dict(
                V='original_generic_source.common_velocity_V_axial5',
                m='original_generic_source.common_own_five_histories_axial5["m"]',
                R='original_generic_source.original_physical_radius',
                S='outer.owner.owner(ends).Pstar'),
            whole_cell_physical_U_Z3_magnitude=[v.record() for v in rg.rows(bounds['U'])],
            whole_cell_physical_m_Z3_magnitude=[v.record() for v in rg.rows(bounds['m'])],
            whole_cell_physical_q_Z2_magnitude=[v.record() for v in bounds['q']],
            whole_cell_value_Z1_integrand_magnitudes={
                k: [v.record() for v in (value.value, value.Z)]
                for k, value in bounds['products'].items()},
            whole_cell_integrated_value_Z1_magnitudes={
                k: [v.record() for v in (value.value, value.Z)]
                for k, value in bounds['integrated'].items()},
            original_positive_dlogR_width=bounds['original_mass'].record(),
            original_R_cap=bounds['original_R'].record(),
            physical_dR_measure=bounds['physical_dR_measure'].record(),
            physical_R_dR_measure=bounds['physical_R_dR_measure'].record(),
            retained_parent_partial_packet=True,
            baseline_rows_are_source_rows_not_range_values=True,
            ordinary_baseline_factorial_conversion_once=True,
            physical_common_scale_Pstar_applied_once=True,
            physical_radius_original_geometry=True,
            whole_cell_range_not_endpoint_value=True,
            numerical_evidence_absolute_whole_cell_magnitude_only=True,
            no_signed_numerical_centers_or_relative_width=True,
            no_finite_endpoint_q_R_or_support_claim=True,
            complete_pressure_tail_support_repair_n2_and_global_claims_not_installed=True,
            execution_seconds=time.monotonic() - began,
            **{GATE: self.acceptance_loaded}, **dict.fromkeys(OPEN, False))
        packet = OriginalPhysicalNonradialPressurePacket()
        state = (bounds, functions, partial_packet)
        self._packets[id(packet)] = (packet, view, state,
                                      fingerprint((view, _bound_record(bounds),
                                                   _z1_encode(functions))))
        return packet

    @source_precision
    def report(self, packet):
        self.assert_graph()
        entry = self._packets.get(id(packet))
        if (type(packet) is not OriginalPhysicalNonradialPressurePacket or
                entry is None or entry[0] is not packet):
            raise ValueError('Live physical nonradial pressure packet required')
        self.parent.report(entry[2][2])
        if fingerprint((entry[1], _bound_record(entry[2][0]),
                        _z1_encode(entry[2][1]))) != entry[3]:
            raise ValueError('Physical pressure packet or source ranges changed')
        return {**copy.deepcopy(entry[1]), GATE: self.acceptance_loaded}


def run(accepted_partial_owner=None):
    began = time.monotonic()
    with mp.workdps(540):
        if accepted_partial_owner is None:
            accepted_modified = partial.parent.CurrentOriginalN1ModifiedLeadingC3(
                partial.rg.CurrentC3TargetRanges(), require_checked=True)
            accepted_partial_owner = partial.CurrentOriginalN1PartialMeanQ(
                accepted_modified, require_checked=True)
        owner = CurrentOriginalN1PhysicalNonradialPressureFields(
            accepted_partial_owner, require_checked=False)
        packets = {chart: owner.range_cell(('0', '.5'), chart,
                                            (1, 2) if chart == 'first_micro' else (1, 1),
                                            (3, 4) if chart == 'first_micro' else (3, 2))
                   for chart in CHARTS}
        raw = dict(candidate_physical_nonradial_pressure_fields_constructed=True,
                   source_family=owner.family,
                   samples={key: owner.report(value) for key, value in packets.items()},
                   signed_physical_functions=_z1_encode(owner.functions),
                   original_graph_prefix_length=owner._graph_length,
                   input_hashes=owner.hashes, **{GATE: False},
                   **dict.fromkeys(OPEN, False), execution_seconds=time.monotonic() - began)
        (HERE / NAME).write_bytes(gzip.compress(
            (json.dumps(c3.encoded(raw), separators=(',', ':')) + '\n').encode(), mtime=0))
        print('PHYSICAL_NONRADIAL_PRESSURE_FIELDS_COMPUTED', len(packets), flush=True)
        return owner, packets


if __name__ == '__main__':
    run()
