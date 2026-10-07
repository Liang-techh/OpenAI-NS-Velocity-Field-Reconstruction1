"""Focused independent geometry/C1/transport checks; synthetic values labelled."""
import ast
from fractions import Fraction
import inspect
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_Rc_function_transport as current

HERE, PREFIX, sha = current.HERE, current.PREFIX, current.sha


def geometry_reference():
    """Extract exact dictionaries from the independently accepted binder."""
    spatial = current.current.transfer.spatial
    tree = ast.parse(inspect.getsource(spatial.affine_identity_theorem))
    assignments = []
    for statement in tree.body[0].body:
        if isinstance(statement, ast.For):
            break
        if isinstance(statement, ast.Assign):
            assignments.append(statement)
    namespace = {'sy': sy}
    exec(compile(ast.fix_missing_locations(ast.Module(body=assignments, type_ignores=[])), '<accepted-radius-definitions>', 'exec'), namespace)
    maps, parameters, x = current.exact_radius_maps()
    renames = {namespace[name]: parameters[target] for name, target in (
        ('P', 'logP'), ('C', 'logC'), ('T', 'T'), ('W', 'Tw'), ('B', 'hbB'),
        ('S', 'hbS'), ('sc', 'sc'), ('M', 'Md'))}
    renames.update({namespace[name]: x for name in ('f', 's', 'x')})
    renames.update({namespace['a']: sy.log(4), namespace['b']: sy.log(100), namespace['d']: sy.log(110)})
    for chart in maps:
        original = namespace['rewritten'][chart].subs(renames, simultaneous=True)
        original_jacobian = sy.sympify(namespace['jacobian'][chart]).subs(renames, simultaneous=True)
        if sy.simplify(maps[chart]-original) != 0:
            raise AssertionError('Different exact radius function: '+chart)
        if sy.simplify(sy.diff(maps[chart], x)-original_jacobian) != 0:
            raise AssertionError('Different original dy/dcoordinate: '+chart)
    return dict(passed=True, original_radius_identities=17, original_Jacobian_identities=17,
        independent_source=spatial.__name__+'.affine_identity_theorem', source_hash=sha(Path(spatial.__file__).name))


def symbolic_C1(built):
    """Check the full AST derivatives with each integral a genuine function."""
    Z = sy.Symbol('Z', real=True)
    A = sy.Function('actual_Rc_A')(Z)
    g = built['graph'];atoms = {};integrals = 0
    for index, cell in enumerate(built['cells']):
        for key, row in cell['contributions'].items():
            if row['value'] == g.zero.node:
                if row['Z'] != g.zero.node: raise AssertionError('Flat increment derivative differs')
                continue
            value = sy.Function('I_'+str(index)+'_'+key)(Z)
            atoms[row['value']], atoms[row['Z']] = value, sy.diff(value, Z)
            integrals += 1
    atoms[built['amplitude'].value.node] = A
    atoms[built['amplitude'].Z.node] = sy.diff(A, Z)
    cache = {}
    def walk(index):
        if index in atoms: return atoms[index]
        if index in cache: return cache[index]
        row = g.nodes[index];op = row['operation']
        if op == 'exact_rational': value = sy.Rational(row['numerator'], row['denominator'])
        elif op == 'original_source_parameter': value = sy.Symbol(row['name'], positive=True)
        elif op == 'shared_positive_integer': value = sy.Symbol('N', integer=True, positive=True)
        elif op == 'sum': value = sum(walk(i) for i in row['arguments'])
        elif op == 'product': value = sy.prod(walk(i) for i in row['arguments'])
        elif op == 'negative': value = -walk(row['argument'])
        elif op == 'positive_quotient': value = walk(row['numerator'])/walk(row['denominator'])
        elif op == 'analytic_unary': value = getattr(sy, row['name'])(walk(row['argument']))
        else: raise AssertionError('Unexpected free/bound node in final function: '+op)
        cache[index] = value
        return value
    checked = 0
    for pairs in (built['history'], built['targets'], built['N_scaled_targets']):
        for key, pair in pairs.items():
            if sy.simplify(sy.diff(walk(pair.value.node), Z)-walk(pair.Z.node)) != 0:
                raise AssertionError('Lost first-Z function term: '+key)
            checked += 1
    C = built['numerator']
    if sy.simplify(sy.diff(walk(C.value.node), Z)-walk(C.Z.node)) != 0:
        raise AssertionError('Joint divided numerator derivative differs')
    checked += 1
    m,k,h,e,p = [walk(built['history'][key].value.node) for key in ('m','k','h','e','p')]
    mu = sy.Symbol('mu', positive=True)
    expected = (m/A, (k-A*m)/(mu*A*A), h/A, e/(A*A), p/(A*A))
    for row, value in zip(current.ROWS, expected):
        if sy.simplify(walk(built['targets'][row].value.node)-value) != 0:
            raise AssertionError('Wrong exact target transformation: '+row)
    return dict(passed=True, symbolic_C1_function_derivatives=checked, exact_target_transformations=5,
        genuine_independent_integral_C1_functions=integrals,
        exact_N_dependent_integrals_not_polynomial_cap_coefficients=True)


def structural_sources(owner, built):
    g = built['graph'];source_rows = 0;quiet_rows = 0
    N_nodes = [i for i, row in enumerate(g.nodes) if row['operation'] == 'shared_positive_integer']
    if N_nodes != [built['N'].node]: raise AssertionError('Multiple frequencies')
    previous = {key: dict(value=g.zero.node, Z=g.zero.node) for key in current.RATES}
    for index, cell in enumerate(built['cells']):
        if cell['incoming'] != previous: raise AssertionError('History reset at '+cell['label'])
        previous = cell['outgoing']
        for key, pair in cell['contributions'].items():
            if cell['source_flat_exact_zero']:
                if pair != dict(value=g.zero.node, Z=g.zero.node): raise AssertionError('Nonzero original quiet increment')
                quiet_rows += 1
                continue
            expected = (owner.views[cell['chart']]['five_signed_increment_rate_roots'][key],
                owner.views[cell['chart']]['five_signed_increment_rate_first_derivatives']['Z'][key])
            for root, derivative in zip(expected, ('value', 'Z')):
                integral = g.nodes[pair[derivative]]
                if integral['operation'] != 'definite_integral': raise AssertionError('Source replaced by cover mass')
                ids = g.nodes[integral['integrand']]['arguments']
                rows = [g.nodes[i] for i in ids if g.nodes[i]['operation'] == 'original_function_graph']
                if len(rows) != 1 or rows[0]['source_node'] != root or rows[0]['chart'] != cell['chart']:
                    raise AssertionError('Wrong full signed density root')
                row = rows[0]
                if row['shared_N'] != N_nodes[0] or row['graph_sha256'] != built['source_graph_sha256']:
                    raise AssertionError('Different source graph/frequency')
                phase = g.nodes[row['phase']]
                if phase['operation'] != 'analytic_unary' or phase['name'] != 'fractional_part':
                    raise AssertionError('Free phase used as spatial phase')
                source_rows += 1
    if len(built['cells']) != 24 or len({cell['chart'] for cell in built['cells']}) != 17:
        raise AssertionError('Different original route')
    return dict(passed=True, bound_full_signed_C0_Z_density_roots=source_rows,
        exact_quiet_contribution_rows=quiet_rows, inherited_C0_Z_cell_rows=24*5*2,
        actual_source_graphs_and_single_global_N_bound=True)


def evaluator_reference(built):
    """Manufactured source, not an evaluation of the original physical field."""
    c = mp.mp.clone();c.dps = 40
    g = built['graph']
    # Test a genuine signed cell with variable native Jacobian (dy/dx=1/x)
    # and nonzero original incoming memory through the quiet final power.
    cell = next(row for row in built['cells'] if row['chart'] == 'actual_patch')
    key = 'e';pair = cell['contributions'][key]
    source_nodes = []
    for index in (pair['value'], pair['Z']):
        integrand = g.nodes[g.nodes[index]['integrand']]
        source_nodes.append(next(g.nodes[i] for i in integrand['arguments'] if g.nodes[i]['operation'] == 'original_function_graph'))
    class SyntheticOracle:
        mode = 'synthetic_reference'
        source_family = built['source_family']
        def parameter(self, name):
            return dict(logP=c.exp(c.mpf('0.1'))+11, logC=c.mpf('2'), T=c.mpf('20'), Tw=c.mpf('4'),
                hbB=c.mpf('.01'), hbS=c.mpf('.02'), sc=c.mpf('.2'), Md=c.mpf('.1'), mu=c.mpf('.03'))[name]
        def source(self, row, *, coordinate, phase, Z, N):
            # Derivative of -(1+Z^2), with both signs and nonzero jets.
            if row['source_node'] == source_nodes[0]['source_node']: return -(1+Z*Z)
            if row['source_node'] == source_nodes[1]['source_node']: return -2*Z
            raise AssertionError('Unexpected manufactured function')
        def integrate(self, function, left, right): return c.quad(function, [left, right])
    oracle = SyntheticOracle();Z=c.mpf('.37');N=160
    value = current.evaluate(built, current.FunctionRef(g,pair['value']), oracle=oracle,Z=Z,N=N,ctx=c)
    jet = current.evaluate(built, current.FunctionRef(g,pair['Z']), oracle=oracle,Z=Z,N=N,ctx=c)
    mass = -c.expm1(-1)
    expected, expected_Z = -(1+Z*Z)*mass, -2*Z*mass
    if abs(value-expected)>c.mpf('1e-32') or abs(jet-expected_Z)>c.mpf('1e-32'):
        raise AssertionError('Signed integral or native Jacobian differs')
    no_oracle = False
    try: current.evaluate(built,built['targets']['M'].value,oracle=None,Z=Z,N=N,ctx=c)
    except current.SourceOracleRequired: no_oracle=True
    if not no_oracle: raise AssertionError('Implicit cap fallback admitted')
    cap_rejected = False
    try: g.mul(g.one, {'log_absolute_upper': '3'})
    except TypeError: cap_rejected=True
    if not cap_rejected: raise AssertionError('Cap accepted as function coefficient')
    return dict(passed=True, synthetic_reference_only=True,
        actual_original_field_numerical_evaluation=False, signed_variable_Jacobian_C0_Z_references=2,
        absolute_error=float(max(abs(value-expected),abs(jet-expected_Z))),
        missing_original_oracle_fails_closed=True, cap_as_function_coefficient_rejected=True)


def run(owner, built=None):
    began = time.monotonic()
    built = owner.build() if built is None else built
    result = dict(all_passed=True, **{current.GATE: True}, source_family=owner.family,
        independent_original_geometry=geometry_reference(),
        symbolic_C1_and_targets=symbolic_C1(built),
        source_and_history_structure=structural_sources(owner,built),
        evaluator_reference=evaluator_reference(built),
        **dict.fromkeys(current.current.packets.OPEN,False),
        numerical_original_source_oracle_installed=False, actual_original_numerical_integrals_evaluated=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,
        current_whole_N_selected=False,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Exact original function-integral geometry, signed source binding, first-Z history/target rules and fail-closed evaluator interface. Synthetic quadrature reference is separate from original-field numerical evaluation; controls, terminal closure/global N/recursion are not admitted.')
    (HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Exact original Rc function transport focused check PASS',flush=True)
    return result
