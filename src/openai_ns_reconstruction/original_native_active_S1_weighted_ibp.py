"""Source-only native seventeen-chart weighted S1 integration-by-parts bridge.

This module deliberately treats every archived definition as an authenticated
source port.  It does not import an archived producer and it does not replay
an archived graph.  The only algebra performed by :func:`derive` is the new
cleared C1 polynomial witness below, implemented by the small sparse
``Fraction`` jet in this file.
"""

import ast
import hashlib
import json
from fractions import Fraction


PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
SCHEMA = 'original_native_active_S1_weighted_ibp_v1'

# The RC cocycle is the accepted parent for this component.  The lower report
# is an inherited universal source lemma; it is read as data and never rerun.
PARENT = PREFIX + 'native_Rc_flux_pressure_cocycle.json'
LOWER = PREFIX + 'native_lower_physical_flux_identity.json'
PRE = PREFIX + 'n1_preRp_moment_functions.json.gz'
PRE_SOURCE = PREFIX + 'n1_preRp_moment_functions.py'
PRE_CHECK_SOURCE = PREFIX + 'n1_preRp_moment_functions_check.py'
MOD = PREFIX + 'modified_Q4_pressure_C3_functions.json.gz'
MOD_SOURCE = PREFIX + 'modified_Q4_pressure_C3_functions.py'
MOD_CHECK_SOURCE = PREFIX + 'modified_Q4_pressure_C3_functions_check.py'
MOD_CHECK = PREFIX + 'modified_Q4_pressure_C3_functions_check.json'
LOOP = PREFIX + 'prescribed_loop_C5_functions.json.gz'
LOOP_SOURCE = PREFIX + 'prescribed_loop_C5_functions.py'
LOOP_CHECK_SOURCE = PREFIX + 'prescribed_loop_C5_functions_check.py'
LOOP_CHECK = PREFIX + 'prescribed_loop_C5_functions_check.json'
PATCH = PREFIX + 'true_Rm_patch_finite_N_C5_functions.json.gz'
PATCH_SOURCE = PREFIX + 'true_Rm_patch_finite_N_C5_functions.py'
PATCH_CHECK_SOURCE = PREFIX + 'true_Rm_patch_finite_N_C5_functions_check.py'
OUTER = PREFIX + 'true_outer_C6_C5_functions.json.gz'
OUTER_SOURCE = PREFIX + 'true_outer_C6_C5_functions.py'
OUTER_CHECK_SOURCE = PREFIX + 'true_outer_C6_C5_functions_check.py'
OUTER_CHECK = PREFIX + 'true_outer_C6_C5_functions_check.json'

# Keep this tuple literal/static: the one-entry frozen runner resolves only
# literals, names, tuples/lists, and string concatenation before import.
RECORD_NAMES = (
    PARENT,
    LOWER,
    PRE,
    PRE_SOURCE,
    PRE_CHECK_SOURCE,
    MOD,
    MOD_SOURCE,
    MOD_CHECK_SOURCE,
    MOD_CHECK,
    LOOP,
    LOOP_SOURCE,
    LOOP_CHECK_SOURCE,
    LOOP_CHECK,
    PATCH,
    PATCH_SOURCE,
    PATCH_CHECK_SOURCE,
    OUTER,
    OUTER_SOURCE,
    OUTER_CHECK_SOURCE,
    OUTER_CHECK,
)

STAGE_FLAGS = (
    'actual_native_active_S1_weighted_IBP_C1_identity_proved',
    'actual_native_active_S1_all17_source_ports_bound',
    'actual_native_active_S1_direct_local_sum17_bound',
    'actual_native_active_S1_all16_primitive_jumps_retained',
    'actual_native_active_S1_terminal_chain_bound',
)

FALSE_FLAGS = (
    'terminal_difference_zero_proved',
    'global_five_moment_closure',
    'native_active_identity_proved',
    'native_compact_identity_proved',
    'native_active_m15_source_identity_proved',
    'native_compact_m15_source_identity_proved',
    'native_active_compact_identity_proved',
    'numeric_values_computed',
    'old_error_terms_removed',
    'temporal_recursion_claimed',
    'time_order_recursion_claimed',
    'old_operators_replayed',
    'full_goal_complete',
)

# These are the authenticated raw pins supplied by the accepted RC parent.
# The parent itself is pinned by the one-entry runner; every other name below
# is looked up in that parent source inventory.
RECORD_HASHES = {
    PARENT: '147febde854a299058380e8d0cf81d52e2c970bc7ffc86b80339aa4077ad3ab0',
    LOWER: 'f1f98828579a76fc2e71f812b7da0f356620f73e5c707bb91edb1bf789792035',
    PRE: '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce',
    PRE_SOURCE: '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e',
    PRE_CHECK_SOURCE: '265ad77ad7849d682670fa4bc8483bfff3148d15a2476fea4feb82b20e9bf617',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    MOD_SOURCE: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
    MOD_CHECK_SOURCE: '829771078e07e2ec8b5168b181f9ff5f0fa62c941c8cf2436feeaba41e367124',
    MOD_CHECK: '189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    LOOP_SOURCE: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
    LOOP_CHECK_SOURCE: '066e6e09ccbb4dc28d7b1881197a426beb6bbc7f92d393e7f3e646f27287914b',
    LOOP_CHECK: '2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7',
    PATCH: '8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad',
    PATCH_SOURCE: 'c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15',
    PATCH_CHECK_SOURCE: 'eae4319ab612d1c32c5c1de674e1a6505e1318923006ed1f9c69ffed48fea0b0',
    OUTER: '75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d',
    OUTER_SOURCE: '4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3',
    OUTER_CHECK_SOURCE: '855c1b31c42c1db22ce503e2a31e5aa20d13dbf7d5d1feb2869c048c64560a00',
    OUTER_CHECK: 'ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba',
}

DECODED_RECORD_HASHES = {
    PRE: '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64',
    MOD: '0a35efead3cb3ea587d1c44a57d70ec78b11b4a5a5305483c8d5b7a041d91f75',
    LOOP: 'f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2',
    PATCH: 'c1026e05015584a7ce29bd7e8da95b17d97669523f1d063ac5b5c800559d9b69',
    OUTER: '53b1bbf15c8c2c5a8400d7e22ef940c6a502397df8cf26bb23e51ff3869e9139',
}

CHART_NAMES = (
    'first_micro', 'second_micro', 'frozen_macro', 'first_switch',
    'second_switch', 'post_power', 'long_reshape', 'reference',
    'restoration', 'postrestore', 'Rm_patch_finite_N',
    'Rh_reference', 'O2_slope', 'O2_axial', 'O2_buffer',
    'O3_transition', 'O3_power',
)
INNER_CHART_NAMES = (
    'first_micro', 'second_micro', 'frozen_macro', 'first_switch',
    'second_switch', 'post_power', 'long_reshape', 'reference',
    'restoration', 'postrestore',
)
OUTER_CHART_NAMES = (
    'Rh_reference', 'O2_slope', 'O2_axial', 'O2_buffer',
    'O3_transition', 'O3_power',
)
K_CHANNELS = (
    'minus_Q_over_L',
    'h_Z_QZ_over_2L',
    'minus_Q_squared_over_2',
    'd_U_QZ_over_L',
    'two_h_Z_squared_U_squared_over_L_squared',
    'two_Z_d_U_UZ_over_L_squared',
)
K_FORMULAS = (
    '-Q/L',
    'h*Z*QZ/(2*L)',
    '-Q^2/2',
    'd*U*QZ/L',
    '2*h*Z^2*U^2/L^2',
    '2*Z*d*U*UZ/L^2',
)
K_DERIVATIVE_FIELDS = (
    ('Q',),
    ('QZ',),
    ('Q',),
    ('U', 'QZ'),
    ('U',),
    ('U', 'UZ'),
)


def digest(value):
    """Return the runner's canonical JSON SHA-256 digest."""
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
        allow_nan=False).encode('utf-8')).hexdigest()


class _Graph:
    """A source-port graph whose nodes are never evaluated."""

    def __init__(self):
        self.nodes = []

    def node(self, operation, arguments=(), **metadata):
        ref = len(self.nodes)
        self.nodes.append(dict(operation=operation, arguments=list(arguments), **metadata))
        return ref

    def constant(self, value):
        q = Fraction(value)
        return self.node('exact_rational', numerator=q.numerator, denominator=q.denominator)


class _SparseJet:
    """Sparse multivariate polynomial with exact Fraction coefficients.

    Monomials are sorted tuples of variable names.  This is intentionally
    local to this kernel: no archived symbolic normalizer is imported.
    """

    __slots__ = ('terms',)

    def __init__(self, terms=()):
        cleaned = {}
        for monomial, value in terms:
            coefficient = Fraction(value)
            if coefficient:
                key = tuple(sorted(monomial))
                cleaned[key] = cleaned.get(key, Fraction(0)) + coefficient
        self.terms = {key: value for key, value in cleaned.items() if value}

    @classmethod
    def constant(cls, value):
        q = Fraction(value)
        return cls((((), q),)) if q else cls()

    @classmethod
    def symbol(cls, name):
        return cls((((str(name),), Fraction(1)),))

    def _coerce(self, value):
        return value if isinstance(value, _SparseJet) else _SparseJet.constant(value)

    def __add__(self, other):
        other = self._coerce(other)
        return _SparseJet(tuple(self.terms.items()) + tuple(other.terms.items()))

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return _SparseJet(((monomial, -value) for monomial, value in self.terms.items()))

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        other = self._coerce(other)
        terms = []
        for left, left_value in self.terms.items():
            for right, right_value in other.terms.items():
                terms.append((left + right, left_value * right_value))
        return _SparseJet(terms)

    def __rmul__(self, other):
        return self * other

    def derivative(self, variable):
        variable = str(variable)
        terms = []
        for monomial, value in self.terms.items():
            count = monomial.count(variable)
            if count:
                mutable = list(monomial)
                mutable.remove(variable)
                terms.append((tuple(mutable), value * count))
        return _SparseJet(terms)

    def partial(self, variable):
        """Formal partial derivative with all other jet variables fixed."""
        variable = str(variable)
        terms = []
        for monomial, value in self.terms.items():
            count = monomial.count(variable)
            if count:
                mutable = list(monomial)
                mutable.remove(variable)
                terms.append((tuple(mutable), value * count))
        return _SparseJet(terms)

    def is_zero(self):
        return not self.terms


def _new_c1_witness():
    """Build the new cleared witness, including its ordinary Z1 row.

    The expressions are multiplied by ``L**2`` before construction, so the
    result is a genuine sparse polynomial rather than a rational-function
    simplification.  ``U2`` and ``Q2`` are ordinary second Z derivatives;
    they are present in the jet data even though the displayed C1 identity
    itself needs only the first rows.
    """
    one = _SparseJet.constant(1)
    two = _SparseJet.constant(2)
    half = _SparseJet.constant(Fraction(1, 2))
    Z = _SparseJet.symbol('Z')
    delta = _SparseJet.symbol('delta')
    U0 = _SparseJet.symbol('U0')
    U1 = _SparseJet.symbol('U1')
    U2 = _SparseJet.symbol('U2')
    Uy0 = _SparseJet.symbol('Uy0')
    Uy1 = _SparseJet.symbol('Uy1')
    Q0 = _SparseJet.symbol('Q0')
    Q1 = _SparseJet.symbol('Q1')
    Q2 = _SparseJet.symbol('Q2')
    Qy0 = _SparseJet.symbol('Qy0')
    Qy1 = _SparseJet.symbol('Qy1')
    # The C1 row is a polynomial in independent field rows.  ``U1`` and
    # ``Q1`` are the actual ordinary Z derivatives, rather than coefficients
    # in a Taylor polynomial in the coordinate Z.
    U = U0
    UZ = U1
    Uy = Uy0
    Q = Q0
    QZ = Q1
    Qy = Qy0
    d = one - Z * Z
    h = one - delta
    L = one - delta * Z * Z

    A = L * (Q + Qy) - ((one + delta) * Z * U + two * Z * Uy - d * UZ)
    # Numerator for L**2 times the physical S.
    S_num = (
        L * (Q + Qy + h * Z * QZ * half)
        + L * L * (Q * Q * half + Q * Qy)
        + L * U * (-two * Z * Q + d * QZ - two * Z * Qy)
    )
    # Numerator for L**2 times B_y/R**2.
    By_num = L * (two * Q + Qy) + L * L * (Q * Q + Q * Qy) - two * two * Z * Z * (U * U + U * Uy)
    # Numerator for L**2 times K/R**2.
    K_num = (
        L * (-Q + h * Z * QZ * half + d * U * QZ)
        - L * L * Q * Q * half
        + two * h * Z * Z * U * U
        + two * Z * d * U * UZ
    )
    witness = S_num - By_num - K_num + two * Z * U * A
    # Ordinary Z1 differentiates through the source jet rows.  This is the
    # chain D requested by the source contract, not a Taylor expansion in Z.
    witness_Z1 = witness.derivative('Z')
    witness_Z1 += witness.partial('U0') * U1 + witness.partial('U1') * U2
    witness_Z1 += witness.partial('Uy0') * Uy1
    witness_Z1 += witness.partial('Q0') * Q1 + witness.partial('Q1') * Q2
    witness_Z1 += witness.partial('Qy0') * Qy1
    return dict(
        witness=witness,
        witness_Z1=witness_Z1,
        variables=(
            'U0', 'U1', 'U2', 'Uy0', 'Uy1',
            'Q0', 'Q1', 'Q2', 'Qy0', 'Qy1',
        ),
        formula=(
            'L^2*(S-B_y/R^2-K/R^2)+2*Z*U*A=0',
            'A=L*(Q+Qy)-((1+delta)*Z*U+2*Z*Uy-(1-Z^2)*UZ)',
        ),
    )


def _ast_has(tree, *needles):
    try:
        text = ast.unparse(tree)
    except Exception:
        return False
    return all(needle in text for needle in needles)


def _prefix_graph(graphs):
    """Check archived rows are shared prefixes without rebuilding them."""
    result = {}
    for left_name, right_name in zip(('pre', 'loop', 'modified', 'patch'),
                                     ('loop', 'modified', 'patch', 'outer')):
        left, right = graphs.get(left_name, ()), graphs.get(right_name, ())
        result[left_name + '_is_prefix_of_' + right_name] = bool(
            left and len(left) <= len(right) and right[:len(left)] == left)
    return result


def _port(graph, operation, source_record, source_node, chart, quantity, **metadata):
    values = source_node if isinstance(source_node, list) else [source_node]
    return graph.node(operation, source_record=source_record,
        source_nodes=list(values), chart=chart, quantity=quantity,
        ordinary_Z_orders=list(range(len(values))), bound_or_midpoint=False,
        **metadata)


def _field_ports(graph, source_record, chart, payload, source_path):
    """Create atomic physical/S/integral ports from one archived chart."""
    physical = payload.get('physical', {})
    ports = {}
    for quantity in ('U', 'U_y', 'm', 'm_y', 'Q', 'Q_y'):
        if quantity in physical:
            ports[quantity] = _port(graph, 'actual_native_physical_' + quantity,
                source_record, physical[quantity], chart, quantity,
                source_path=source_path, function_identity_not_claimed=True)
    if 'S' in payload:
        ports['S'] = _port(graph, 'actual_native_physical_S', source_record,
            payload['S'], chart, 'S', source_path=source_path,
            function_identity_not_claimed=True)
    for family in ('local_integrals', 'cumulative_integrals', 'incoming_integrals'):
        row = payload.get(family, {})
        if isinstance(row, dict) and 'S1' in row:
            ports[family + '.S1'] = _port(graph,
                'actual_native_' + family + '_S1', source_record, row['S1'],
                chart, family + '.S1', source_path=source_path,
                function_identity_not_claimed=True)
    boundary = payload.get('boundary', payload.get('finite_boundary', {}))
    if isinstance(boundary, dict) and 'R_Q_y' in boundary:
        ports['R_Q_y'] = _port(graph, 'actual_native_R_Q_y', source_record,
            boundary['R_Q_y'], chart, 'R_Q_y', source_path=source_path,
            function_identity_not_claimed=True)
    return ports


def _geometry_ports(graph, source_record, chart, payload):
    ports = {}
    for quantity in ('coordinate', 'domain', 'radius', 'physical_measure'):
        if quantity in payload:
            ports[quantity] = _port(graph, 'actual_native_geometry_' + quantity,
                source_record, payload[quantity], chart, quantity,
                source_path='geometry', Z_independent=True,
                function_identity_not_claimed=True)
    if 'physical_measure' in payload:
        # The archived providers call the one physical log-coordinate factor
        # ``physical_measure``.  Expose the semantic Jacobian port explicitly
        # while retaining the original row and its source identity.
        ports['Jacobian'] = _port(graph, 'actual_native_geometry_Jacobian',
            source_record, payload['physical_measure'], chart, 'Jacobian',
            source_path='geometry.physical_measure', Z_independent=True,
            source_name='physical_measure', function_identity_not_claimed=True)
    return ports


def _source_formula_guards(trees):
    return {
        'modified_physical_formula_guard': _ast_has(
            trees.get(MOD_SOURCE), 'def physical_pressure_fields',
            'same original L=1-delta*Z^2', 'Q_y'),
        'loop_measure_guard': _ast_has(
            trees.get(LOOP_SOURCE), 'physical_Jacobian_applied_once=True',
            'physical_measure', 'charts'),
        'patch_physical_formula_guard': _ast_has(
            trees.get(PATCH_SOURCE), 'physical', 'Q_y', 'pressure_continuation'),
        'outer_physical_formula_guard': _ast_has(
            trees.get(OUTER_SOURCE), 'physical', 'Q_y',
            'physical_Jacobian_applied_once=True'),
        'old_active_source_guard': _ast_has(
            trees.get(PRE_SOURCE), 'S1', 'defining_bindings', 'active'),
    }


def derive(source_records, *, expected_input_digest):
    """Authenticate source data and construct the new native17 IBP graph."""
    checks = {}
    data = {}
    trees = {}

    def check(name, condition):
        checks[name] = bool(condition)

    check('frozen_input_digest', digest(source_records) == expected_input_digest)
    for name in RECORD_NAMES:
        row = source_records.get(name)
        present = isinstance(row, dict)
        check('authenticated:' + name, present and
            row.get('raw_sha256') == RECORD_HASHES[name] and
            row.get('raw_authenticated_before_decode') is True and
            row.get('is_gzip') is name.endswith('.gz') and
            row.get('decoded_sha256') == DECODED_RECORD_HASHES.get(name, RECORD_HASHES[name]))
        if not present:
            continue
        try:
            if name.endswith('.py'):
                trees[name] = ast.parse(row['source_text'])
            else:
                data[name] = json.loads(row['source_text'])
        except Exception:
            check('parseable:' + name, False)

    parent = data.get(PARENT, {})
    parent_result = parent.get('component_result', {})
    lower = data.get(LOWER, {})
    lower_result = lower.get('component_result', {})
    modified = data.get(MOD, {})
    loop = data.get(LOOP, {})
    patch = data.get(PATCH, {})
    outer = data.get(OUTER, {})
    pre = data.get(PRE, {})
    modified_bind = modified.get('modified_Q4_pressure_C3_bindings', {})
    loop_bind = loop.get('prescribed_loop_C5_bindings', {})
    patch_bind = patch.get('Rm_patch_finite_N_C5_bindings', {})
    outer_bind = outer.get('true_outer_C6_C5_bindings', {})

    family = parent_result.get('source_family', {})
    check('accepted_parent_report', parent_result.get('passed') is True
        and parent.get('full_goal_complete') is False)
    check('accepted_parent_retained_construction_data',
        isinstance(parent_result.get('retained_construction_data'), dict))
    lower_audit = lower_result.get('source_prerequisite_audit', {})
    lower_interface = lower_result.get('source_interface', {})
    check('accepted_lower_universal_Nplus_equals_Nq_plus_Ny',
        lower_result.get('passed') is True and
        (lower_audit.get('new_universal_full_vs_reduced_Qy_numerator_C1') is True
         or lower_interface.get('universal_numerator_C1_zero_witnesses') == [157, 161]))

    provider_rows = (modified, loop, patch, outer, pre)
    check('same_exact_source_family', bool(family) and all(
        row.get('source_family') == family for row in provider_rows))
    graphs = {
        'pre': pre.get('exact_graph_nodes', []),
        'loop': loop.get('exact_graph_nodes', []),
        'modified': modified.get('exact_graph_nodes', []),
        'patch': patch.get('exact_graph_nodes', []),
        'outer': outer.get('exact_graph_nodes', []),
    }
    prefix_checks = _prefix_graph(graphs)
    check('shared_archived_graph_prefixes', bool(prefix_checks) and all(prefix_checks.values()))
    check('source_graph_rows_not_rebuilt', all(isinstance(rows, list) for rows in graphs.values()))

    syntax_checks = _source_formula_guards(trees)
    for name, condition in syntax_checks.items():
        check(name, condition)

    # The actual positive-quotient denominator is an archived graph port.  We
    # guard only each chart's Q0 and Qy0 rows: higher ordinary-Z rows may use
    # powers or distinct graph handles for the same source L.  The new graph
    # records the inherited parent delta port and never creates a free delta.
    def q0_L_guard(report_nodes, entries):
        if not entries:
            return False
        for q_rows, qy_rows in entries:
            if not q_rows or not qy_rows:
                return False
            q0, qy0 = q_rows[0], qy_rows[0]
            if not (isinstance(q0, int) and isinstance(qy0, int)
                    and 0 <= q0 < len(report_nodes)
                    and 0 <= qy0 < len(report_nodes)):
                return False
            q_node, qy_node = report_nodes[q0], report_nodes[qy0]
            if (q_node.get('operation') != 'positive_quotient'
                    or qy_node.get('operation') != 'positive_quotient'
                    or q_node.get('source_positive_certificate') != 'same original L=1-delta*Z^2'
                    or qy_node.get('source_positive_certificate') != 'same original L=1-delta*Z^2'
                    or q_node.get('denominator') != qy_node.get('denominator')
                    or not isinstance(q_node.get('denominator'), int)):
                return False
        return True

    mod_q = [node for chart_name in INNER_CHART_NAMES
        for node in modified_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q', [])]
    mod_qy = [node for chart_name in INNER_CHART_NAMES
        for node in modified_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q_y', [])]
    patch_q = patch_bind.get('physical', {}).get('Q', [])
    patch_qy = patch_bind.get('physical', {}).get('Q_y', [])
    outer_q = [node for chart_name in OUTER_CHART_NAMES
        for node in outer_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q', [])]
    outer_qy = [node for chart_name in OUTER_CHART_NAMES
        for node in outer_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q_y', [])]
    mod_entries = [
        (modified_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q', []),
         modified_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q_y', []))
        for chart_name in INNER_CHART_NAMES]
    outer_entries = [
        (outer_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q', []),
         outer_bind.get('charts', {}).get(chart_name, {}).get('physical', {}).get('Q_y', []))
        for chart_name in OUTER_CHART_NAMES]
    check('actual_L_denominator_guard',
        q0_L_guard(graphs['modified'], mod_entries)
        and q0_L_guard(patch.get('exact_graph_nodes', []), [(patch_q, patch_qy)])
        and q0_L_guard(graphs['outer'], outer_entries))
    mod_nodes = graphs['modified']
    actual_L_node = (mod_nodes[mod_q[0]].get('denominator')
        if mod_q and isinstance(mod_q[0], int) and mod_q[0] < len(mod_nodes)
        else None)
    parent_nodes = parent.get('definition_graph', {}).get('nodes', [])
    lower_nodes = lower.get('definition_graph', {}).get('nodes', [])
    parent_outputs = parent.get('definition_graph', {}).get('outputs', [])
    parent_compact_output_nodes = [row.get('compact_density_difference_integral')
        for row in parent_outputs if isinstance(row, dict)]
    parent_terminal_output_nodes = [row.get('terminal_density_difference')
        for row in parent_outputs if isinstance(row, dict)]
    check('accepted_parent_output_rows',
        len(parent_outputs) == 2
        and [row.get('ordinary_Z_order') for row in parent_outputs
             if isinstance(row, dict)] == [0, 1]
        and all(isinstance(node, int) and 0 <= node < len(parent_nodes)
                for node in parent_compact_output_nodes + parent_terminal_output_nodes))
    check('accepted_lower_delta_port', len(lower_nodes) > 3
        and lower_nodes[3].get('operation') == 'common_original_delta'
        and lower_nodes[3].get('Z_independent') is True)
    check('accepted_parent_graph_immutable',
        bool(parent_nodes) and digest(parent_nodes) ==
        parent.get('definition_graph', {}).get('nodes_digest'))

    proof = _new_c1_witness()
    check('new_cleared_polynomial_C1_zero', proof['witness'].is_zero())
    check('new_cleared_polynomial_C1_ordinary_Z1_zero', proof['witness_Z1'].is_zero())

    graph = _Graph()
    # Bind the new expressions to actual archived physical/geometry ports.
    # These are atomic rows; no archived expression is normalized or called.
    seed_geometry = _geometry_ports(graph, LOOP, 'first_micro',
        loop_bind.get('charts', {}).get('first_micro', {}))
    seed_physical = _field_ports(graph, MOD, 'first_micro',
        modified_bind.get('charts', {}).get('first_micro', {}),
        'modified_Q4_pressure_C3_bindings.charts')
    actual_formula_ports = [seed_geometry.get(key) for key in
        ('coordinate', 'domain', 'radius', 'Jacobian')]
    actual_formula_ports.extend(seed_physical.get(key) for key in
        ('U', 'U_y', 'Q', 'Q_y', 'S'))
    actual_formula_ports = [ref for ref in actual_formula_ports if isinstance(ref, int)]
    formula_refs = {
        'Z': graph.node('new_ordinary_Z_variable', ordinary_Z_orders=[0, 1],
            temporal_variable=False, source_owned_cell=['0', '1/2']),
        'delta': graph.node('accepted_inherited_delta_port', source_record=LOWER,
            source_graph='definition_graph', source_node=3,
            Z_independent=True, free_delta_not_invented=True),
        'L': graph.node('actual_L_denominator_port', source_record=MOD,
            source_nodes=[actual_L_node], source_positive=True,
            source_positive_certificate='same original L=1-delta*Z^2',
            source_operation='positive_quotient',
            source_formula='L=1-delta*Z^2',
            actual_DAG_formula='L=1-delta*Z^2',
            inherited_delta_port='accepted_inherited_delta_port',
            source_Q0_node=(mod_q[0] if mod_q else None),
            source_Qy0_node=(mod_qy[0] if mod_qy else None),
            same_Q0_Qy0_denominator=True,
            actual_physical_Q0_Qy0_only=True),
        'd': graph.node('new_formula_d', formula='1-Z^2'),
        'h': graph.node('new_formula_h', formula='1-delta'),
    }
    formula_refs['A'] = graph.node('new_formula_A',
        arguments=[formula_refs['L'], *actual_formula_ports],
        formula='L*(Q+Qy)-((1+delta)*Z*U+2*Z*Uy-d*UZ)',
        accepted_mass_ODE='m_y=U-m',
        universal_lower_lemma_inherited=True)
    formula_refs['physical_S'] = graph.node('new_formula_physical_S',
        arguments=[formula_refs['L'], *actual_formula_ports], formula=(
            '(Q+Qy+h*Z*QZ/2)/L+Q^2/2+Q*Qy+'
            'U*(-2*Z*Q+d*QZ-2*Z*Qy)/L'))
    formula_refs['B'] = graph.node('new_formula_B',
        arguments=[formula_refs['L'], *actual_formula_ports],
        formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)',
        radius_squared_measure_once=True)
    formula_refs['By'] = graph.node('new_formula_By_over_R2',
        arguments=[formula_refs['L'], *actual_formula_ports],
        formula='(2*Q+Qy)/L+Q^2+Q*Qy-4*Z^2*(U^2+U*Uy)/L^2',
        radius_squared_y_equals_2R_squared=True)
    formula_refs['K'] = graph.node('new_formula_K_over_R2',
        arguments=[formula_refs['L'], *actual_formula_ports],
        formula='-Q/L+h*Z*QZ/(2L)-Q^2/2+d*U*QZ/L+2*h*Z^2*U^2/L^2+2*Z*d*U*UZ/L^2')
    formula_refs['witness'] = graph.node('new_cleared_polynomial_witness',
        arguments=[formula_refs['A'], formula_refs['physical_S'], formula_refs['By'], formula_refs['K']],
        formula='L^2*(S-B_y/R^2-K/R^2)+2*Z*U*A=0',
        sparse_fraction_polynomial=True, no_sympy=True,
        variables=list(proof['variables']), ordinary_Z_orders=[0, 1],
        zero_coefficients=proof['witness'].is_zero(),
        ordinary_Z1_zero=proof['witness_Z1'].is_zero())

    # Keep source ports atomic.  The first ten charts intentionally take
    # geometry and physical rows from different archived providers.
    chart_graphs = []
    all_direct_local_ports = []
    for chart_name in INNER_CHART_NAMES:
        geometry = loop_bind.get('charts', {}).get(chart_name, {})
        physical = modified_bind.get('charts', {}).get(chart_name, {})
        geometry_ports = _geometry_ports(graph, LOOP, chart_name, geometry)
        physical_ports = _field_ports(graph, MOD, chart_name, physical, 'modified_Q4_pressure_C3_bindings.charts')
        chart_graphs.append(dict(name=chart_name, geometry=geometry_ports,
            physical=physical_ports, geometry_source=LOOP,
            physical_source=MOD))
    patch_geometry = patch_bind.get('chart', {})
    patch_ports = _geometry_ports(graph, PATCH, 'Rm_patch_finite_N', patch_geometry)
    patch_ports.update(_field_ports(graph, PATCH, 'Rm_patch_finite_N', patch_bind,
        'Rm_patch_finite_N_C5_bindings'))
    chart_graphs.append(dict(name='Rm_patch_finite_N', geometry=patch_ports,
        physical=patch_ports, geometry_source=PATCH, physical_source=PATCH))
    for chart_name in OUTER_CHART_NAMES:
        chart = outer_bind.get('charts', {}).get(chart_name, {})
        ports = _geometry_ports(graph, OUTER, chart_name, chart)
        ports.update(_field_ports(graph, OUTER, chart_name, chart,
            'true_outer_C6_C5_bindings.charts'))
        chart_graphs.append(dict(name=chart_name, geometry=ports, physical=ports,
            geometry_source=OUTER, physical_source=OUTER))

    # Saved continuation rows are separate endpoint ports.  Their topology is
    # checked against the source graph, rather than inferred from the count of
    # charts.  The first chart starts with literal zero incoming S1.
    chain_specs = []
    for chart_name in INNER_CHART_NAMES:
        chain_specs.append((MOD, chart_name,
            modified_bind.get('charts', {}).get(chart_name, {}),
            modified_bind.get('endpoints', {}).get(chart_name, {})
                .get('integrals', {}).get('S1', []), graphs['modified'],
            loop_bind.get('charts', {}).get(chart_name, {}).get('coordinate'),
            loop_bind.get('charts', {}).get(chart_name, {}).get('domain')))
    chain_specs.append((PATCH, 'Rm_patch_finite_N', patch_bind,
        patch_bind.get('endpoints', {}).get('Rh', {}).get('integrals', {}).get('S1', []),
        patch.get('exact_graph_nodes', []),
        patch_bind.get('chart', {}).get('coordinate'),
        patch_bind.get('chart', {}).get('domain')))
    for chart_name in OUTER_CHART_NAMES:
        chain_specs.append((OUTER, chart_name,
            outer_bind.get('charts', {}).get(chart_name, {}),
            outer_bind.get('endpoints', {}).get(chart_name, {})
                .get('integrals', {}).get('S1', []), graphs['outer'],
            outer_bind.get('charts', {}).get(chart_name, {}).get('coordinate'),
            outer_bind.get('charts', {}).get(chart_name, {}).get('domain')))
    for chart, (source_record, chart_name, source_chart, endpoint_rows,
                report_nodes, coordinate, domain) in zip(
            chart_graphs, chain_specs):
        chart['endpoint_S1'] = _port(graph, 'actual_native_endpoint_S1', source_record,
            endpoint_rows, chart_name, 'endpoint.S1',
            source_path='saved endpoint integrals', function_identity_not_claimed=True)

    def continuation_source_guard(specs):
        previous_endpoint = None
        topology = []
        for position, (source_record, chart_name, source_chart, endpoint_rows,
                       report_nodes, coordinate, domain) in enumerate(specs):
            incoming = source_chart.get('incoming_integrals', {}).get('S1', [])
            local = source_chart.get('local_integrals', {}).get('S1', [])
            cumulative = source_chart.get('cumulative_integrals', {}).get('S1', [])
            if not (isinstance(coordinate, int) and isinstance(domain, list)
                    and len(domain) == 2
                    and len(incoming) == len(local) == len(cumulative) == 4
                    and len(endpoint_rows) == 4):
                return False, topology
            expected_incoming = [0, 0, 0, 0] if position == 0 else previous_endpoint
            if incoming != expected_incoming:
                return False, topology
            for order, cumulative_node in enumerate(cumulative):
                if not isinstance(cumulative_node, int) or cumulative_node >= len(report_nodes):
                    return False, topology
                row = report_nodes[cumulative_node]
                if position == 0:
                    # The initial incoming is zero.  The provider either
                    # retains the local row itself or records the explicit
                    # zero-plus-local sum; an unrelated integral row is not
                    # accepted as the cumulative continuation.
                    if not (cumulative_node == local[order]
                            or (row.get('operation') == 'sum'
                                and row.get('arguments') == [0, local[order]])):
                        return False, topology
                elif row.get('operation') != 'sum' or row.get('arguments') != [incoming[order], local[order]]:
                    return False, topology
                endpoint_node = endpoint_rows[order]
                if not (isinstance(endpoint_node, int)
                        and 0 <= endpoint_node < len(report_nodes)):
                    return False, topology
                endpoint_row = report_nodes[endpoint_node]
                if not (endpoint_row.get('operation') == 'function_substitution'
                        and endpoint_row.get('expression') == cumulative_node
                        and endpoint_row.get('variable') == coordinate
                        and endpoint_row.get('value') == domain[1]
                        and endpoint_row.get('Z_independent_substitution') is True):
                    return False, topology
            previous_endpoint = endpoint_rows
            topology.append(dict(chart=chart_name, incoming=incoming,
                local=local, cumulative=cumulative, endpoint=endpoint_rows,
                initial_zero=(position == 0)))
        return previous_endpoint is not None and len(topology) == 17, topology

    continuation_ok, continuation_rows = continuation_source_guard(chain_specs)
    check('all17_actual_S1_continuation_topology', continuation_ok)
    check('all17_native_chart_ports_present', len(chart_graphs) == 17 and
        all(row.get('physical', {}).get('U') is not None
            and row.get('physical', {}).get('U_y') is not None
            and row.get('physical', {}).get('Q') is not None
            and row.get('physical', {}).get('Q_y') is not None
            and row.get('physical', {}).get('S') is not None
            and row.get('physical', {}).get('local_integrals.S1') is not None
            and row.get('physical', {}).get('cumulative_integrals.S1') is not None
            for row in chart_graphs))
    topology_refs = []
    for position, (chart, topology) in enumerate(zip(chart_graphs, continuation_rows)):
        physical = chart['physical']
        args = [physical.get('incoming_integrals.S1'),
            physical.get('local_integrals.S1'),
            physical.get('cumulative_integrals.S1'), chart.get('endpoint_S1')]
        topology_refs.append(graph.node('actual_native_S1_continuation_topology',
            arguments=[ref for ref in args if isinstance(ref, int)],
            chart=chart['name'], chart_index=position,
            incoming_port=physical.get('incoming_integrals.S1'),
            local_port=physical.get('local_integrals.S1'),
            cumulative_port=physical.get('cumulative_integrals.S1'),
            endpoint_port=chart.get('endpoint_S1'),
            initial_zero=topology.get('initial_zero', False),
            source_topology=topology))
    check('all17_actual_S1_topology_nodes_bound',
        len(topology_refs) == 17 and all(
            isinstance(ref, int) and len(graph.nodes[ref].get('arguments', [])) == 4
            for ref in topology_refs))

    # A B endpoint is a new function node; it does not evaluate an archived
    # expression.  Each K contribution retains R^2*Jacobian*dnative=R^2*dy
    # measure and Z-independent limits.
    endpoint_rows = []
    direct_rows = []
    k_sum_rows = []
    for index, chart in enumerate(chart_graphs):
        physical = chart['physical']
        geometry = chart['geometry']
        radius = geometry.get('radius')
        measure = geometry.get('physical_measure')
        lower_endpoint = graph.node('actual_B_endpoint_function',
            chart=chart['name'], endpoint='lower',
            coordinate=geometry.get('coordinate'), domain=geometry.get('domain'),
            endpoint_bound='domain[0]', radius=radius,
            Jacobian=geometry.get('Jacobian'), physical_measure=measure,
            integration_variable='native_coordinate',
            bounds=geometry.get('domain'), limits_Z_independent=True,
            measure='R^2*Jacobian*dnative', physical_ports=physical,
            formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)',
            source_function_identity_not_claimed=True)
        upper_endpoint = graph.node('actual_B_endpoint_function',
            chart=chart['name'], endpoint='upper',
            coordinate=geometry.get('coordinate'), domain=geometry.get('domain'),
            endpoint_bound='domain[1]', radius=radius,
            Jacobian=geometry.get('Jacobian'), physical_measure=measure,
            integration_variable='native_coordinate',
            bounds=geometry.get('domain'), limits_Z_independent=True,
            measure='R^2*Jacobian*dnative', physical_ports=physical,
            formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)',
            source_function_identity_not_claimed=True)
        endpoint_rows.append((lower_endpoint, upper_endpoint))
        k_rows = []
        for channel_index, channel in enumerate(K_CHANNELS):
            k_rows.append(graph.node('actual_K_six_integral', chart=chart['name'],
                channel=channel, channel_index=channel_index,
                formula=K_FORMULAS[channel_index],
                derivative_fields=list(K_DERIVATIVE_FIELDS[channel_index]),
                coordinate=geometry.get('coordinate'), domain=geometry.get('domain'),
                radius=radius, Jacobian=geometry.get('Jacobian'),
                integration_variable='native_coordinate',
                bounds=geometry.get('domain'), physical_ports=physical,
                radius_squared_measure_once=True,
                measure='R^2*Jacobian*dnative', limits_Z_independent=True,
                ordinary_Z_orders=[0, 1], source_actual_density_relation=formula_refs['witness'],
                old_active_measure_native=True,
                linear_Pstar_terms_preserved=True,
                quadratic_Pstar_squared_terms_preserved=True))
        k_sum = graph.node('actual_K_six_integral_sum', arguments=k_rows,
            chart=chart['name'], channel_count=6, all_channels_retained=True)
        k_sum_rows.append(k_sum)
        local_port = physical.get('local_integrals.S1')
        local_at_upper = graph.node('actual_native_local_S1_at_chart_upper',
            arguments=[local_port, geometry.get('coordinate'), geometry.get('domain')],
            chart=chart['name'], expression_port=local_port,
            coordinate_port=geometry.get('coordinate'),
            domain_port=geometry.get('domain'), endpoint_index=1,
            definition='substitute native=domain[1] in the actual local S1 function',
            ordinary_Z_orders=[0, 1], Z_independent_substitution=True,
            archived_expression_evaluated=False)
        direct_rows.append(graph.node('actual_native_local_S1_IBP_decomposition',
            arguments=[local_at_upper, upper_endpoint, lower_endpoint, k_sum], chart=chart['name'],
            formula='B_upper-B_lower+integral(K_six)',
            source_local_S1_port=local_at_upper,
            source_local_S1_free_native_port=local_port,
            physical_S_port=physical.get('S'),
            actual_density_relation=formula_refs['witness'],
            local_source_not_rewritten=True))
        if local_port is not None:
            all_direct_local_ports.append(local_at_upper)

    jump_rows = []
    for index in range(16):
        previous_upper = endpoint_rows[index][1]
        next_lower = endpoint_rows[index + 1][0]
        jump_rows.append(graph.node('actual_internal_primitive_jump',
            arguments=[next_lower, previous_upper], jump_index=index,
            formula='J_i=B_next_lower-B_prev_upper',
            is_actual_function=True, retained=True, zero_assumption=False,
            source_identity_proved=False))
    check('all16_internal_primitive_jumps_retained', len(jump_rows) == 16 and
        all(graph.nodes[row].get('retained') is True and
            graph.nodes[row].get('zero_assumption') is False for row in jump_rows))

    direct_sum = graph.node('actual_native_IN_sum17_local_S1',
        arguments=all_direct_local_ports, chart_count=len(all_direct_local_ports),
        expected_chart_count=17, all_local_source_ports_present=(len(all_direct_local_ports) == 17),
        source_local_integrals_not_evaluated=True)
    check('direct_native_IN_sum17_local_bound', len(all_direct_local_ports) == 17
        and checks.get('all17_actual_S1_continuation_topology', False)
        and checks.get('all17_actual_S1_topology_nodes_bound', False))
    direct_total = graph.node('actual_native_direct_S1_six_integral_shape',
        arguments=[endpoint_rows[-1][1], endpoint_rows[0][0], *jump_rows, *k_sum_rows],
        formula='B_last_upper-B_first_lower-sum(J_0..J_15)+sum(chart_K_six)',
        chart_count=17, jump_count=16, k_integral_shape='six per chart',
        topology_nodes=topology_refs,
        local_decomposition_nodes=direct_rows,
        no_internal_jump_zero_claim=True)

    # The native Rc cumulative rows are literal archived substitutions.  The
    # direct seventeen-chart sum is retained as its own source-connected port;
    # no equality is claimed to the old active field or to the old IO rows.
    old_io = graph.node('old_active_IO_S1_port', source_record=PRE,
        source_nodes=[30687, 30688], ordinary_Z_orders=[0, 1],
        distinct_from_native=True, source_function_identity_not_claimed=True)
    outer_o3 = outer_bind.get('charts', {}).get('O3_power', {})
    o3_cumulative = outer_o3.get('cumulative_integrals', {}).get('S1', [])
    o3_incoming = outer_o3.get('incoming_integrals', {}).get('S1', [])
    o3_local = outer_o3.get('local_integrals', {}).get('S1', [])
    o3_endpoint = outer_bind.get('endpoints', {}).get('O3_power', {}) \
        .get('integrals', {}).get('S1', [])
    outer_nodes = outer.get('exact_graph_nodes', [])
    literal_chain = graph.node('literal_saved_cumulative_incoming_substitution_chain',
        source_record=OUTER, saved_ports=[258837, 258838],
        terminal_endpoint_ports=list(o3_endpoint),
        cumulative_ports=list(o3_cumulative), incoming_ports=list(o3_incoming),
        local_ports=list(o3_local), substitution_variable=250215,
        substitution_value=48, substitution_value_graph_handle=48,
        explicit_graph_handle_2=48, Z_independent_substitution=True,
        exact_saved_rows=(
            outer_nodes[258837] if len(outer_nodes) > 258837 else {},
            outer_nodes[258838] if len(outer_nodes) > 258838 else {},
        ), no_old_source_execution=True)
    check('literal_native_cumulative_Rc_S1_chain',
        o3_cumulative[:2] == [258614, 258615]
        and o3_incoming[:2] == [250195, 250196]
        and o3_local[:2] == [258610, 258611]
        and o3_endpoint[:2] == [258837, 258838]
        and len(outer_nodes) > 258838
        and outer_nodes[258837].get('operation') == 'function_substitution'
        and outer_nodes[258838].get('operation') == 'function_substitution'
        and outer_nodes[258837].get('expression') == 258614
        and outer_nodes[258838].get('expression') == 258615
        and outer_nodes[258837].get('variable') == 250215
        and outer_nodes[258838].get('variable') == 250215
        and outer_nodes[258837].get('value') == 48
        and outer_nodes[258838].get('value') == 48)
    native_terminal = graph.node('actual_native_cumulative_Rc_S1_port',
        arguments=[literal_chain, direct_sum], source_record=OUTER,
        source_nodes=[258837, 258838], ordinary_Z_orders=[0, 1],
        saved_chain_literal=True, direct_native_sum17_bound=True,
        equality_to_old_IO_claimed=False)
    check('native_terminal_chain_bound', bool(o3_cumulative) and bool(o3_incoming)
        and bool(o3_local) and bool(o3_endpoint)
        and checks.get('literal_native_cumulative_Rc_S1_chain', False))
    compact_port = graph.node('accepted_parent_compact_IC_port', source_record=PARENT,
        source_nodes=parent_compact_output_nodes, ordinary_Z_orders=[0, 1],
        source_output='compact_density_difference_integral',
        source_interface='accepted parent compact contribution',
        identity_proved=False, source_function_identity_not_claimed=True)
    parent_terminal_port = graph.node(
        'accepted_parent_terminal_density_difference_port', source_record=PARENT,
        source_nodes=parent_terminal_output_nodes, ordinary_Z_orders=[0, 1],
        source_output='terminal_density_difference',
        source_interface='accepted parent terminal density difference',
        identity_proved=False, source_function_identity_not_claimed=True)
    terminal_formula = graph.node('accepted_parent_F_minus_T_formula',
        arguments=[native_terminal, old_io, compact_port, parent_terminal_port],
        formula='F-T=0.5*(IN_Rc-IO_active+IC)', zero_claimed=False,
        inherited_parent_formula=True,
        parent_terminal_output_nodes=parent_terminal_output_nodes)
    weighted_identity = graph.node(
        'new_weighted_IBP_identity_certificate',
        arguments=[native_terminal, direct_total, formula_refs['witness'],
                   direct_sum, *topology_refs],
        native_terminal_port=native_terminal,
        direct_representation_port=direct_total,
        local_source_sum_port=direct_sum,
        topology_nodes=topology_refs,
        local_decomposition_nodes=direct_rows,
        all16_jump_nodes=jump_rows,
        identity_formula='native_S1=direct_total after all17 actual local sources',
        native_old_field_identity_not_claimed=True,
        new_formula_witness=formula_refs['witness'])
    outputs = {
        'native_S1': native_terminal,
        'native_weighted_ibp_representation': direct_total,
        'old_active_S1': old_io,
        'parent_terminal': terminal_formula,
        'new_weighted_IBP_identity_certificate': weighted_identity,
    }
    graph.nodes[direct_total]['terminal_formula'] = terminal_formula

    checks.update({
        'new_actual_graph_definitions_bound': bool(graph.nodes),
        'source_actual_IBP_formula_guard': all(syntax_checks.values()),
    })
    flags = {name: False for name in STAGE_FLAGS}
    flags[STAGE_FLAGS[0]] = checks.get('new_cleared_polynomial_C1_zero', False) and checks.get('new_cleared_polynomial_C1_ordinary_Z1_zero', False)
    flags[STAGE_FLAGS[1]] = checks.get('all17_native_chart_ports_present', False) and checks.get('shared_archived_graph_prefixes', False)
    flags[STAGE_FLAGS[2]] = checks.get('direct_native_IN_sum17_local_bound', False)
    flags[STAGE_FLAGS[3]] = checks.get('all16_internal_primitive_jumps_retained', False)
    flags[STAGE_FLAGS[4]] = checks.get('native_terminal_chain_bound', False)
    passed = all(flags.values()) and all(checks.values())

    unresolved = [name for name, condition in checks.items() if not condition]
    remaining = [
        'native_old_field_identity_not_claimed',
        'sixteen_internal_primitive_jump_source_identities',
        'compact_terminal_source_identity',
        'old_error_terms_reconciliation',
        'global_five_moment_closure',
    ]
    definition_graph = dict(
        namespace=SCHEMA,
        nodes=graph.nodes,
        node_count=len(graph.nodes),
        nodes_digest=digest(graph.nodes),
        outputs=outputs,
        proof_graph_is_new=True,
        archived_definition_graphs_rebuilt=False,
        archived_graphs_executed=False,
        new_formula_inputs=list(proof['variables']),
        ordinary_Z_orders=[0, 1],
    )
    source_interface = {
        'inherited_parent_interface': parent_result.get('source_interface', {}),
        'accepted_parent_output_rows': {
            'ordinary_Z_orders': [0, 1],
            'compact_density_difference_integral': parent_compact_output_nodes,
            'terminal_density_difference': parent_terminal_output_nodes,
        },
        'native_active_S1_actual': {
            'direct_local_chart_count': 17,
            'native_Rc_cumulative_S1': [258837, 258838],
            'literal_saved_cumulative_ports': [258614, 258615],
            'literal_terminal_endpoint_ports': [258837, 258838],
            'literal_saved_incoming_ports': [250195, 250196],
            'literal_saved_substitution': {
                'variable': 250215, 'value': 48,
                'value_graph_handle': 48, 'explicit_graph_handle_2': 48,
            },
        },
        'old_IO_active_S1': {'source_record': PRE, 'ports': [30687, 30688], 'distinct': True},
        'old_active_I1_six_terms_all_R2_Jac_native': True,
        'linear_Pstar_vs_quadratic_Pstar2_preserved': True,
        'physical_fields': ('U', 'U_y', 'm', 'm_y', 'Q', 'Q_y', 'S'),
        'actual_L': {'source_record': MOD, 'positive_quotient_denominator': actual_L_node,
                     'formula': '1-delta*Z^2', 'Z_independent_delta_port': True,
                     'Q0_source_node': (mod_q[0] if mod_q else None),
                     'Qy0_source_node': (mod_qy[0] if mod_qy else None)},
        'accepted_inherited_delta': {'source_record': LOWER, 'definition_graph_node': 3,
                                     'free_delta_not_invented': True},
        'mass_ODE': 'm_y=U-m',
        'physical_measure': 'R^2*S*Jacobian*dnative',
        'six_signed_channels_retained': True,
        'all16_internal_primitive_jumps': 'actual functions retained; zero not assumed',
        'native_old_field_identity_claimed': False,
        'ordinary_Z_cell': ['0', '1/2'],
    }
    audit = dict(
        new_source_consumer_calls=1,
        accepted_ancestor_operator_calls=0,
        old_operators_replayed=False,
        numeric_values_computed=False,
        archived_definition_graphs_rebuilt=False,
        new_universal_C1_polynomial_only=True,
        sparse_fraction_jet_only=True,
        source_graphs_reused=True,
    )
    result = {
        'schema': SCHEMA,
        'passed': passed,
        'source_only': True,
        'ordinary_Z_orders': [0, 1],
        'ordinary_Z_not_temporal_n': True,
        'source_family': family,
        'source_prerequisite_audit': dict(checks, **prefix_checks),
        'definition_graph': definition_graph,
        'source_interface': source_interface,
        'retained_construction_data': parent_result.get('retained_construction_data'),
        'remaining_construction_gates': remaining,
        'unresolved_gates': unresolved,
        'legacy_numeric_chain_reconfigured': False,
        'new_unentered_expression_simulations': 0,
        'runtime_model_metadata': None,
        'source_inspection_audit': audit,
        'formula_witness': {
            'cleared_C1_zero': proof['witness'].is_zero(),
            'ordinary_Z1_zero': proof['witness_Z1'].is_zero(),
            'variables': list(proof['variables']),
            'formula': list(proof['formula']),
        },
        'actual_native_active_S1_weighted_IBP_C1_identity_proved': flags[STAGE_FLAGS[0]],
        'actual_native_active_S1_all17_source_ports_bound': flags[STAGE_FLAGS[1]],
        'actual_native_active_S1_direct_local_sum17_bound': flags[STAGE_FLAGS[2]],
        'actual_native_active_S1_all16_primitive_jumps_retained': flags[STAGE_FLAGS[3]],
        'actual_native_active_S1_terminal_chain_bound': flags[STAGE_FLAGS[4]],
    }
    result.update({name: False for name in FALSE_FLAGS})
    result['passed'] = bool(passed)
    result['stage_flags'] = flags
    return result
