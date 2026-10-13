"""Signed original n1 pressure C1 functions on r_minus <= R <= Rp.

This consumes genuine finite-N field graphs, not pressure cover centers.
All imported bodies retain their real integrands, coordinates and histories.
The missing inner/collar/pulse adapters are deliberately not replaced by zero.
"""
import copy
from functools import reduce
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator as linear
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as active
import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives as quiet
import lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1 as endpoint
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1 as accuracy

source, controls = linear.source, linear.controls
HERE, PREFIX, sha = linear.HERE, linear.PREFIX, linear.sha
NAME = PREFIX + 'current_original_n1_preRp_moment_functions.json.gz'
RECEIPT = PREFIX + 'current_original_n1_preRp_moment_functions_check.json'
GATE = 'current_original_n1_rminus_to_Rp_signed_pressure_C1_functions_installed'
OPEN = ('n1_complete_signed_defect_C1_function_adapter_installed',
        'n1_selected_linear_coefficient_functions_installed',
        'n1_repaired_own_M_V_P_installed', 'positive_order_five_moment_repair_closed',
        'temporal_recursion', 'full_corrected_Navier_Stokes_solution')


def encode(value):
    if type(value) is source.C1Function:
        return [value.value.node, value.Z.node]
    if type(value) is source.FunctionRef:
        return value.node
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value

# Only expression operands are relocated. Derivative orders, rational
# numerators, source-file node indices and certificate metadata are not.
REFERENCE_KEYS = {
    'exact_rational': (), 'bound_variable': (), 'mathematical_pi': (),
    'current_original_source_parameter': (), 'original_source_parameter': (),
    'exact_positive_integer_expression': (), 'shared_positive_integer': (),
    'sum': ('arguments',), 'product': ('arguments',),
    'negative': ('argument',), 'analytic_unary': ('argument',),
    'positive_quotient': ('numerator', 'denominator'),
    'exact_real_comparison': ('left', 'right'), 'logical_or': ('arguments',),
    'compact_raw_beta': ('argument',),
    'exact_original_source_flat_sigma': ('argument',),
    'current_original_leading_function_recipe': ('coordinate', 'Z_variable'),
    'current_original_physical_source_scale': (),
    'current_original_physical_source_quantity': ('coordinate', 'Z_variable'),
    'current_original_physical_radius': ('log_radius', 'free_query'),
    'substitute_original_native_coordinate': ('body', 'query_coordinate'),
    'substitute_free_query_at_native_coordinate': ('body', 'free_query'),
    'substitute_original_inverse_angle': ('body', 'inverse_angle'),
    'function_substitution': ('expression', 'variable', 'value'),
    'definite_integral': ('integrand', 'lower', 'upper'),
    'original_lazy_flat_branch': ('active_body', 'Delta', 'eta',
        'flat_predicate', 'original_collar_coordinate', 'original_collar_upper'),
    'original_monotone_phase_inverse': ('phase_function', 'target_phase',
        'angle_lower', 'angle_upper', 'Delta', 'eta'),
}


class GraphImport:
    """Iterative exact DAG import, with an explicit operand schema."""

    def __init__(self, graph, original):
        self.g, self.original = graph, original
        common = 0
        while common < min(len(graph.nodes), len(original.nodes)):
            if graph.nodes[common] != original.nodes[common]:
                break
            common += 1
        self.mapping = {i: source.FunctionRef(graph, i) for i in range(common)}
        self.common_prefix = common

    def _operands(self, node):
        row = self.original.nodes[node]
        if row['operation'] not in REFERENCE_KEYS:
            raise ValueError('Unspecified graph operand schema: ' + row['operation'])
        for key in REFERENCE_KEYS[row['operation']]:
            if key not in row or row[key] is None:
                continue
            items = row[key] if isinstance(row[key], list) else (row[key],)
            for item in items:
                if type(item) is not int or not 0 <= item < len(self.original.nodes):
                    raise ValueError('Invalid original graph operand ' + key)
                yield item

    def __call__(self, value):
        if type(value) is not source.FunctionRef or value.graph is not self.original:
            raise TypeError('A genuine source-graph FunctionRef is required')
        if not 0 <= value.node < len(self.original.nodes):
            raise ValueError('Source graph node absent')
        stack = [(value.node, False)]
        visiting = set()
        while stack:
            node, ready = stack.pop()
            if node in self.mapping:
                continue
            if not ready:
                if node in visiting:
                    raise ValueError('Source graph has a cyclic expression')
                visiting.add(node)
                stack.append((node, True))
                stack.extend((item, False) for item in reversed(list(self._operands(node)))
                             if item not in self.mapping)
                continue
            row = copy.deepcopy(self.original.nodes[node])
            op = row.pop('operation')
            for key in REFERENCE_KEYS[op]:
                if key in row and row[key] is not None:
                    row[key] = ([self.mapping[i].node for i in row[key]]
                                if isinstance(row[key], list) else self.mapping[row[key]].node)
            self.mapping[node] = self.g.node(op, **row)
            visiting.remove(node)
        return self.mapping[value.node]

    def c1(self, pair):
        return source.C1Function(self(pair.value), self(pair.Z))


class C1Algebra:
    def __init__(self, graph):
        self.g = graph

    def fixed(self, value):
        ref = value if type(value) is source.FunctionRef else self.g.constant(value)
        self.g.ids([ref])
        return source.C1Function(ref, self.g.zero)

    def add(self, *pairs):
        return reduce(self.g.c1add, pairs, self.fixed(0))

    def neg(self, pair):
        return self.g.c1scale(self.g.constant(-1), pair)

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def mul(self, *pairs):
        return reduce(self.g.c1mul, pairs, self.fixed(1))

    def scale(self, scalar, pair):
        return self.g.c1scale(self.fixed(scalar).value, pair)

    def div(self, left, right, certificate='original positive L=1-delta*Z^2'):
        g = self.g
        return source.C1Function(g.quotient(left.value, right.value, certificate),
            g.quotient(g.sub(g.mul(left.Z, right.value), g.mul(left.value, right.Z)),
                       g.mul(right.value, right.value), certificate + ' squared'))

    def at(self, pair, variable, value):
        g = self.g
        g.ids([variable, value, pair.value, pair.Z])
        return source.C1Function(*(g.node('function_substitution',
            expression=ref.node, variable=variable.node, value=value.node,
            Z_independent_substitution=True) for ref in (pair.value, pair.Z)))


def pressure_functions(g, delta, integrals, lower, upper):
    """Reduced-S identity, including both genuine finite viscous endpoints."""
    a = C1Algebra(g)
    z = source.C1Function(g.symbol('Z'), g.one)
    one = a.fixed(1)
    d = a.sub(one, a.mul(z, z))
    c = a.fixed(delta)
    h = a.sub(one, c)
    L = a.sub(one, a.mul(c, z, z))
    L2 = a.mul(L, L)
    half = lambda p: a.scale('1/2', p)
    boundary = lambda name: a.sub(upper[name], lower[name])
    i = integrals
    S0 = a.add(a.div(boundary('R_q'), L),
        half(a.div(a.mul(h, z, i['I0_q_Z']), L)),
        half(boundary('R_q_squared')), a.div(a.mul(d, i['I0_U_q_Z']), L),
        a.scale(-2, a.div(a.mul(z, z, boundary('R_U_squared')), L2)),
        a.scale(-2, a.div(a.mul(c, z, z, i['I0_U_squared']), L2)),
        a.scale(2, a.div(a.mul(z, d, i['I0_U_U_Z']), L2)))
    S1 = a.add(a.div(boundary('R_squared_q'), L),
        a.neg(a.div(i['I1_q'], L)), half(a.div(a.mul(h, z, i['I1_q_Z']), L)),
        half(boundary('R_squared_q_squared')), a.neg(half(i['I1_q_squared'])),
        a.div(a.mul(d, i['I1_U_q_Z']), L),
        a.scale(-2, a.div(a.mul(z, z, boundary('R_squared_U_squared')), L2)),
        a.scale(2, a.div(a.mul(h, z, z, i['I1_U_squared']), L2)),
        a.scale(2, a.div(a.mul(z, d, i['I1_U_U_Z']), L2)))
    return dict(m13=a.add(a.neg(half(S0)), boundary('q_plus_q_y')),
                m15=a.sub(half(S1), boundary('R_q_y')))


def endpoint_functions(g, U, q, qy, radius):
    a = C1Algebra(g)
    R = a.fixed(radius)
    return dict(R_q=a.mul(R, q), R_squared_q=a.mul(R, R, q),
        R_q_squared=a.mul(R, q, q), R_squared_q_squared=a.mul(R, R, q, q),
        R_U_squared=a.mul(R, U, U), R_squared_U_squared=a.mul(R, R, U, U),
        q_plus_q_y=a.add(q, qy), R_q_y=a.mul(R, qy))


def q_plus_q_y(g, delta, U, Uy):
    a = C1Algebra(g)
    z = source.C1Function(g.symbol('Z'), g.one)
    d = a.sub(a.fixed(1), a.mul(z, z))
    L = a.sub(a.fixed(1), a.mul(a.fixed(delta), z, z))
    Uz = source.C1Function(U.Z, U.ZZ)
    u = source.C1Function(U.value, U.Z)
    uy = source.C1Function(Uy.value, Uy.Z)
    return a.div(a.add(a.mul(a.add(a.fixed(1), a.fixed(delta)), z, u),
        a.scale(2, a.mul(z, uy)), a.neg(a.mul(d, Uz))), L)


def lower_velocity_y_source(fields, ends):
    """Actual unprojected leading V_y ordinary rows at the finite lower end."""
    view = endpoint._source_view(fields, tuple(ends), 'first_micro', ('sc_half',))
    packet = view['source_packet']
    if not packet.get('source_owned_collar_q_C0_Z_exact_zero'):
        raise ValueError('Original lower collar source theorem required')
    endpoint.c0._zero_first_incoming_proof(fields, tuple(ends))
    endpoint.c0._flat_driver_proof(fields, tuple(ends), 'first_micro')
    values = packet['actual_background_source']['V_y']
    flow = fields.modified.ranges.phase.outer.owner.owner(tuple(ends)).flow
    for k, value in enumerate(values[:4]):
        active._same_source(value, fields.c, flow, 'true lower V_y row ' + str(k))
    return tuple(value * math.factorial(k) for k, value in enumerate(values[:4]))


def active_functions(operator, fields):
    g = operator.graph
    importer = GraphImport(g, fields.g)
    ref = lambda i: source.FunctionRef(fields.g, i)
    algebra = C1Algebra(g)
    integrals = {}
    for chart in fields.charts:
        f = fields.functions[chart]
        key = active._interval_key(fields.partitions[chart][0], fields.partitions[chart][-1])
        for name, unions in f['defining_integrals_by_cell'].items():
            pair = importer.c1(unions[key])
            integrals[name] = algebra.add(integrals.get(name, algebra.fixed(0)), pair)
    delta = importer(fields.delta_node)
    c3 = active.c3
    pair = lambda p: importer.c1(p)
    endpoints = []
    for chart, side in (('first_micro', 0), ('O3_power', -1)):
        f = fields.functions[chart]
        cell = fields.modified.target.phase.windows[chart]['actual_current_source_cells'][side]
        coordinate = importer(ref(cell['lower' if side == 0 else 'upper']))
        query = importer(f['query'])
        at = lambda p: algebra.at(pair(p), query, coordinate)
        U, q = at(f['U']), at(f['q'])
        R = algebra.at(algebra.fixed(importer(f['radius'])), query, coordinate).value
        if side == 0:
            template = fields.modified.target.primitives[chart]['source_roots']['V']
            attrs = copy.deepcopy(fields.g.nodes[template.value.node])
            # The actual same-source provider supplies V_y, not a cap or a
            # coefficient selected from an interval. Retain its defining
            # leading recipe and live ordinary-row projection explicitly.
            rows = []
            for order in range(2):
                rows.append(g.node('current_original_physical_source_quantity',
                    source_family=fields.family, native_chart=chart, quantity='V_y',
                    coordinate=coordinate.node,
                    Z_variable=importer(ref(attrs['Z_variable'])).node,
                    Z_order=order, ordinary_slow_Z_derivative_order=order,
                    Taylor_coefficient_factorial=math.factorial(order),
                    raw_quantity_path='actual_background_source.V_y[%d]' % order,
                    source_projection_binding=active.current.ast_binding(lower_velocity_y_source),
                    original_leading_source_recipe=attrs['recipe'],
                    defining_quantity_not_a_range_value=True,
                    finite_N_lower_initial_history_and_local_driver_exact_zero=True))
            S = importer(f['scale'].value)
            Uy = source.C1Function(*(g.mul(S, v) for v in rows))
            Uc3 = c3.C3Function(at(f['U']).value, at(f['U']).Z,
                algebra.at(source.C1Function(importer(f['U'].ZZ), importer(f['U'].ZZZ)),
                           query, coordinate).value, g.zero)
            qp = q_plus_q_y(g, delta, Uc3, Uy)
            qy = algebra.sub(qp, q)
        else:
            # Checked endpoint theorem U=Uy=UZ=0; own nonlocal q is retained.
            qy = algebra.neg(q)
        endpoints.append(endpoint_functions(g, U, q, qy, R))
    return pressure_functions(g, delta, integrals, *endpoints), dict(
        exact_signed_integrals=encode(integrals),
        finite_endpoints=encode(endpoints), import_mapping=importer.mapping,
        original_graph_common_prefix_length=importer.common_prefix,
        original_all_chart_graph_node_count=len(fields.g.nodes))


def compact_functions(operator, scale, delta):
    g = operator.graph
    raw = operator.reports[linear.quiet.band.NAME]
    f = raw['actual_C3_repair_band_functions']
    c3 = active.c3
    ref = lambda i: source.FunctionRef(g, i)
    jet = lambda ids: c3.C3Function(*map(ref, ids))
    alg = c3.C3Algebra(g)
    S = c3.C3Function(scale, g.zero, g.zero, g.zero)
    profile = f['profiles_at_band_x']['radial_y_rows']
    U = alg.mul(S, jet(profile[0]['V']))
    m = alg.mul(S, jet(f['complete_histories']['m']))
    Uy = alg.mul(S, jet(profile[1]['V']))
    my = alg.mul(S, jet(f['complete_history_y']['m']))
    q = active.partial.recover_Q_increment(g, U, m, delta)
    qy = active.partial.recover_Q_increment(g, Uy, my, delta)
    products = active.physical.nonradial_integrands(g, U, q)
    x = ref(f['original_band_variable'])
    Rc = g.unary('exp', ref(operator.background['logRc']))
    R = g.mul(Rc, x)
    integrals = {}
    for name, p in products.items():
        factor = Rc if name.startswith('I0_') else g.mul(R, Rc)
        integrals[name] = source.C1Function(*(controls.integral(g,
            g.mul(factor, row), g.nodes[x.node]['name'], g.one, g.constant(2),
            measure='physical dR=Rc*dx' if name.startswith('I0_') else 'physical R*dR=Rc^2*x*dx',
            endpoint_independent_of_Z=True, ordinary_slow_Z_derivative_order=k)
            for k, row in enumerate((p.value, p.Z))))
    a = C1Algebra(g)
    boundary = []
    for edge in (g.one, g.constant(2)):
        at = lambda p: a.at(source.C1Function(p.value, p.Z), x, edge)
        boundary.append(endpoint_functions(g, at(U), at(q), at(qy), g.mul(Rc, edge)))
    return pressure_functions(g, delta, integrals, *boundary), dict(
        physical_integrals=encode(integrals), finite_endpoints=encode(boundary),
        actual_C3_original_band_variable=x.node, common_Pstar_applied_once=True)


def quiet_functions(operator, scale, quiet_owner):
    g = operator.graph
    raw = quiet_owner.frame_raw
    original = linear.realization.restore_graph(raw['exact_graph_nodes'])
    importer = GraphImport(g, original)
    f = raw['actual_Rp_native_frame_functions']
    ref = lambda i: importer(source.FunctionRef(original, i))
    Rc = g.unary('exp', source.FunctionRef(g, operator.background['logRc']))
    Rp = g.unary('exp', ref(f['original_logRp']))
    # The accepted exact quiet identity is m_common=kappa*Z. Its first
    # ordinary derivative is the genuine constant kappa, including Z=0.
    kappa = ref(f['actual_identified_Rp_histories_C3']['m'][1])
    K = g.mul(Rp, scale, kappa)
    lower = g.mul(g.constant(2), Rc)
    m13 = g.mul(g.constant('1/4'), K, K,
        g.sub(g.quotient(g.one, lower, 'positive physical 2Rc'),
              g.quotient(g.one, Rp, 'positive original Rp')))
    logratio = g.sub(ref(f['original_logRp']), g.add(g.unary('log', g.constant(2)),
                               source.FunctionRef(g, operator.background['logRc'])))
    m15 = g.mul(g.constant('-1/4'), K, K, logratio)
    return dict(m13=source.C1Function(m13, g.zero),
                m15=source.C1Function(m15, g.zero)), dict(
        K=K.node, kappa=kappa.node, lower=lower.node, upper=Rp.node,
        constant_Z_identity=quiet_owner.actual_source_identity(),
        original_native_frame_import_mapping=importer.mapping)


class CurrentOriginalN1PreRpMomentFunctions:
    def __init__(self, fields, operator=None):
        if type(fields) is not active.CurrentOriginalN1AllChartPhysicalPressureFields or not fields.acceptance_loaded:
            raise TypeError('Accepted live original all-chart fields required')
        fields.assert_graph()
        self.operator = operator or linear.CurrentOriginalN1ReservedLinearOperator()
        self.operator.assert_graph()
        if not self.operator.acceptance_loaded or self.operator.family != fields.family:
            raise ValueError('Same original checked linear operator required')
        selected = fields.g.nodes[fields.N.node]
        expected = self.operator.reports[linear.quiet.band.NAME]['actual_selected_repair_integer']
        if (selected['operation'] != 'exact_positive_integer_expression'
                or any(selected[key] != expected[key] for key in ('type', 'J', 'definition'))):
            raise ValueError('Active fields and compact controls selected N differ')
        self.fields, self.g, self.family = fields, self.operator.graph, copy.deepcopy(fields.family)
        self.quiet = quiet.CurrentOriginalN1QuietPressurePrimitives()
        if self.quiet.family != self.family:
            raise ValueError('Quiet source family differs')
        # Read the original provider once on every accepted full axial cell.
        for ends in accuracy.CELLS:
            lower_velocity_y_source(fields, ends)
            upper = endpoint._source_view(fields, ends, 'O3_power', (2, 1))
            values = upper['source_packet']['original_generic_source']['common_velocity_V_axial5']
            if not all(endpoint.c0._exact_zero_record(v) for v in values[:4]):
                raise ValueError('Actual upper O3 U and Z rows must vanish')
        self.pieces, self.bindings = {}, {}
        self.pieces['active'], self.bindings['active'] = active_functions(self.operator, fields)
        importer = GraphImport(self.g, fields.g)
        scale = importer(fields.functions['first_micro']['scale'].value)
        delta = importer(fields.delta_node)
        self.pieces['compact'], self.bindings['compact'] = compact_functions(self.operator, scale, delta)
        self.pieces['quiet'], self.bindings['quiet'] = quiet_functions(self.operator, scale, self.quiet)
        a = C1Algebra(self.g)
        self.functions = {key: a.add(*(piece[key] for piece in self.pieces.values()))
                          for key in ('m13', 'm15')}
        self.hashes = dict(self.operator.hashes)
        for incoming in (fields.hashes, self.quiet.hashes):
            for name, digest in incoming.items():
                if name in self.hashes and self.hashes[name] != digest:
                    raise ValueError('Conflicting same-function source ' + name)
                self.hashes[name] = digest
        for module in (active, quiet, endpoint, accuracy):
            for name in (Path(module.__file__).name, module.NAME, module.RECEIPT):
                self.hashes[name] = sha(name)
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self.prefix = copy.deepcopy(self.g.nodes)
        self.encoded_functions = encode(self.functions)

    def assert_graph(self):
        self.fields.assert_graph(); self.operator.assert_graph(); self.quiet.assert_source()
        if self.g.nodes[:len(self.prefix)] != self.prefix or encode(self.functions) != self.encoded_functions:
            raise ValueError('Signed pressure functions changed')


class SavedCurrentOriginalN1PreRpMomentFunctions:
    """Checked archive-only C1 handles in the executable linear operator graph.

    Restoring this adapter does not rebuild the leading-field ancestors or
    replace a true source expression by a signed numerical cover center.
    """

    def __init__(self):
        self.raw = linear.read(NAME)
        receipt = linear.read(RECEIPT)
        self.operator = linear.CurrentOriginalN1ReservedLinearOperator()
        if (not receipt.get('all_passed') or not receipt.get(GATE)
                or receipt['source_family'] != self.operator.family
                or self.raw['source_family'] != self.operator.family
                or any(receipt.get(key) for key in OPEN)):
            raise ValueError('Checked same-family signed pressure graph required')
        for name, digest in receipt['input_hashes'].items():
            if sha(name) != digest:
                raise ValueError('Changed signed pressure function source ' + name)
        self.g, self.family = self.operator.graph, copy.deepcopy(self.operator.family)
        nodes = self.raw['exact_graph_nodes']
        if nodes[:len(self.g.nodes)] != self.g.nodes:
            raise ValueError('Accepted linear operator prefix differs')
        for index in range(len(self.g.nodes), len(nodes)):
            row = copy.deepcopy(nodes[index])
            op = row.pop('operation')
            if self.g.node(op, **row).node != index:
                raise ValueError('Archived pressure function expression identity differs')
        ref = lambda i: source.FunctionRef(self.g, i)
        restore = lambda pair: source.C1Function(*map(ref, pair))
        self.functions = {key: restore(pair)
                          for key, pair in self.raw['signed_pressure_C1_functions'].items()}
        self.pieces = {label: {key: restore(pair) for key, pair in fields.items()}
                       for label, fields in self.raw['signed_pieces'].items()}
        self.hashes = copy.deepcopy(receipt['input_hashes'])
        self._prefix = copy.deepcopy(self.g.nodes)
        self._functions = encode(self.functions)
        self.acceptance_loaded = True

    def assert_graph(self):
        self.operator.assert_graph()
        if (self.g.nodes[:len(self._prefix)] != self._prefix
                or encode(self.functions) != self._functions):
            raise ValueError('Saved signed original pressure functions changed')


def run(fields):
    began = time.monotonic()
    owner = CurrentOriginalN1PreRpMomentFunctions(fields)
    bindings = copy.deepcopy(owner.bindings)
    for value in bindings.values():
        for key in ('import_mapping', 'original_native_frame_import_mapping'):
            if key in value:
                value[key] = {str(i): r.node for i, r in value[key].items()}
    raw = dict(source_family=owner.family, exact_graph_nodes=owner.g.nodes,
        signed_pressure_C1_functions=encode(owner.functions),
        signed_pieces=encode(owner.pieces), defining_bindings=bindings,
        radial_scope='original r_minus to Rp; active, compact, quiet counted once',
        inner_lower_collar_and_outer_pulse_not_included=True,
        numerical_enclosures_not_selected_function_coefficients=True,
        complete_five_defects_not_installed=True, input_hashes=owner.hashes,
        **dict.fromkeys((GATE,) + OPEN, False), execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',', ':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_PRE_RP_SIGNED_C1_FUNCTIONS_DEFINED', len(owner.g.nodes), flush=True)
    return owner, raw


if __name__ == '__main__':
    raise SystemExit('run(fields) consumes a retained checked all-chart source owner')
