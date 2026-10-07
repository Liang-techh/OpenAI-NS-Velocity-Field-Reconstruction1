"""Exact source-function integrals and Rc targets, separate from cap arithmetic.

This installs the integral representation of the original N-dependent
functions. It does not install a numerical oracle for the original point
phase inverse, choose a globally admissible N, or solve the five controls.
The optional evaluator requires an explicit source oracle; it never reads
coefficient enclosures as values. Synthetic oracle checks are labelled so.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets as current
import lei_ren_part1_paper_compliant_current_generic_loop_function_sources as sources

HERE, PREFIX, sha = current.HERE, current.PREFIX, current.sha
NAME = PREFIX + 'current_native_Rc_function_transport.json'
RECEIPT = PREFIX + 'current_native_Rc_function_transport_check.json'
GATE = 'current_original_exact_N_dependent_Rc_source_integral_and_target_functions_defined'
ROWS = current.repair.ROWS
RATES = {key: Fraction(value) for key, value in current.RATES.items()}


@dataclass(frozen=True)
class FunctionRef:
    """A function expression handle; deliberately not a LogUpper/enclosure."""
    graph: object
    node: int


@dataclass(frozen=True)
class C1Function:
    value: FunctionRef
    Z: FunctionRef


class FunctionTransportGraph:
    def __init__(self):
        self.nodes = []
        self.keys = {}
        self.zero = self.constant(0)
        self.one = self.constant(1)

    def node(self, operation, **attributes):
        row = dict(operation=operation, **attributes)
        key = json.dumps(row, sort_keys=True, separators=(',', ':'))
        if key not in self.keys:
            self.keys[key] = len(self.nodes)
            self.nodes.append(row)
        return FunctionRef(self, self.keys[key])

    def ids(self, arguments):
        if any(type(q) is not FunctionRef or q.graph is not self for q in arguments):
            raise TypeError('Same-graph function handles required; caps are not values')
        return [q.node for q in arguments]

    def constant(self, value):
        q = Fraction(value)
        return self.node('exact_rational', numerator=q.numerator, denominator=q.denominator)

    def add(self, *arguments):
        ids = self.ids(arguments)
        ids = [q for q in ids if q != self.zero.node]
        return self.zero if not ids else FunctionRef(self, ids[0]) if len(ids) == 1 else self.node('sum', arguments=ids)

    def mul(self, *arguments):
        ids = self.ids(arguments)
        if self.zero.node in ids:
            return self.zero
        ids = [q for q in ids if q != self.one.node]
        return self.one if not ids else FunctionRef(self, ids[0]) if len(ids) == 1 else self.node('product', arguments=ids)

    def neg(self, value):
        self.ids([value])
        return self.zero if value == self.zero else self.node('negative', argument=value.node)

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def quotient(self, left, right, certificate):
        self.ids([left, right])
        if right == self.zero:
            raise ValueError('Zero function denominator')
        return self.zero if left == self.zero else self.node('positive_quotient', numerator=left.node,
            denominator=right.node, source_positive_certificate=certificate)

    def unary(self, name, argument):
        self.ids([argument])
        return self.node('analytic_unary', name=name, argument=argument.node)

    def symbol(self, name):
        return self.node('bound_variable', name=name)

    def c1add(self, left, right):
        return C1Function(self.add(left.value, right.value), self.add(left.Z, right.Z))

    def c1mul(self, left, right):
        return C1Function(self.mul(left.value, right.value), self.add(self.mul(left.Z, right.value), self.mul(left.value, right.Z)))

    def c1scale(self, factor, pair):
        return C1Function(self.mul(factor, pair.value), self.mul(factor, pair.Z))


def exact_radius_maps():
    """The accepted radius-minus-inlet formulas, with no cap substitution."""
    P, C, T, W, B, S, sc, Md, x = sy.symbols('logP logC T Tw hbB hbS sc Md x', real=True)
    log4, log100, log110 = sy.log(4), sy.log(100), sy.log(110)
    common = 14*P + 10*C + log110 - log4 + 1000 - B*sc/2
    maps = dict(
        bridge_first=B*(x-sc/2), bridge_second=B*(x-sc/2),
        bridge_macro=4*P*x+(log100-log4+1000)*x+B*(2*(1-x)-sc/2),
        switch_first=4*P+log100-log4+1000+S*x-B*sc/2,
        switch_second=4*P+log100-log4+1000+S*x-B*sc/2,
        switch_power=4*P+log100-log4+1000+(log110-log100)*x+2*S*(1-x)-B*sc/2,
        reshape=4*P+log110-log4+1000+T*x-B*sc/2,
        inner_reference=(4+10*x)*P+10*x*C+(1-x)*T+log110-log4+1000-8*x-B*sc/2,
        axial_restore=common-8+x, restore_buffer=common+x, actual_patch=common-6+sy.log(x),
        Rh_reference=common+x, O2_slope=common+x, O2_axial=common+sy.exp(Md*x),
        O2_buffer=common+sy.exp(Md)+x, O3_slope_mu=common+P+x,
        O3_power=common+P+1+W*x)
    return maps, dict(logP=P, logC=C, T=T, Tw=W, hbB=B, hbS=S, sc=sc, Md=Md), x


def expression(graph, expr, parameters, coordinate=None):
    """Translate only exact symbolic arithmetic, never string eval."""
    if expr.is_Rational:
        return graph.constant(Fraction(int(expr.p), int(expr.q)))
    if expr.is_Symbol:
        if str(expr) == 'x':
            if coordinate is None:
                raise ValueError('Explicit native coordinate binding required')
            return coordinate
        return parameters[str(expr)]
    if expr.func is sy.Add:
        return graph.add(*(expression(graph, q, parameters, coordinate) for q in expr.args))
    if expr.func is sy.Mul:
        return graph.mul(*(expression(graph, q, parameters, coordinate) for q in expr.args))
    if expr.func is sy.Pow and expr.args[1].is_Integer:
        base = expression(graph, expr.args[0], parameters, coordinate)
        power = int(expr.args[1])
        if power < 0:
            return graph.quotient(graph.one, expression(graph, expr.args[0]**(-power), parameters, coordinate), 'original_positive_scale')
        return graph.mul(*([base]*power))
    if expr.func in (sy.exp, sy.log):
        return graph.unary('exp' if expr.func is sy.exp else 'log', expression(graph, expr.args[0], parameters, coordinate))
    if expr == sy.E:
        return graph.unary('exp', graph.one)
    raise ValueError('Unsupported exact radius expression: ' + str(expr))


def exact_route_endpoints(route, parameters):
    result = []
    for label, chart, left, right in route:
        if label == 'initial_flat_collar':
            lo, hi = parameters['sc']/2, 3*parameters['sc']/4
        elif label == 'active_first_bridge':
            lo, hi = 3*parameters['sc']/4, sy.Integer(1)
        elif chart == 'actual_patch':
            lo, hi = sy.Integer(1), sy.E
        elif chart == 'O3_power':
            lo, hi = sy.Rational(str(left))/parameters['Tw'], sy.Rational(str(right))/parameters['Tw']
        else:
            lo, hi = sy.Rational(str(left)), sy.Rational(str(right))
        result.append((label, chart, lo, hi))
    return result


def exact_geometry_theorem():
    maps, parameters, x = exact_radius_maps()
    route = exact_route_endpoints(current.ROUTE, parameters)
    seams = []
    for left, right in zip(route, route[1:]):
        end = maps[left[1]].subs(x, left[3])
        start = maps[right[1]].subs(x, right[2])
        if sy.simplify((end-start).subs(parameters['logP'], sy.exp(parameters['Md'])+11)) != 0:
            raise ArithmeticError('Function-route gap: ' + left[0] + ' -> ' + right[0])
        seams.append(left[0] + ' -> ' + right[0])
    if sy.simplify(maps[route[0][1]].subs(x, route[0][2])) != 0:
        raise ArithmeticError('Function-route inlet differs')
    for chart, value in maps.items():
        if value.has(sy.Symbol('Z')):
            raise ArithmeticError('Z-dependent geometry needs boundary terms')
    return dict(passed=True, true_cells=len(route), charts=len(maps), exact_no_gap_joins=seams,
        original_inlet_offset_exact_zero=True, endpoints_weights_and_global_phase_Z_independent=True,
        positive_Jacobians_and_widths_from_checked_native_geometry=True,
        exact_parameter_relation='logP=exp(Md)+11', velocity_or_higher_derivative_seams_not_admitted=True)


class NativeRcFunctionTransport:
    def __init__(self, owner):
        if type(owner) is not current.NativeRcParameterTargets:
            raise TypeError('Accepted original Rc target owner required')
        self.owner, self.family, self.service = owner, owner.family, owner.service
        self.hashes = {}
        for module in (current, sources):
            receipt = json.loads((HERE/module.RECEIPT).read_bytes())
            manifest = json.loads((HERE/module.NAME).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or manifest['source_family'] != self.family:
                raise ValueError('Same-family checked source/target prerequisites required')
            self.hashes.update(receipt['input_hashes'])
            self.hashes.update({module.NAME: sha(module.NAME), module.RECEIPT: sha(module.RECEIPT)})
        self.views = json.loads(gzip.decompress((HERE/sources.VIEWS).read_bytes()))
        self.hashes.update({sources.VIEWS: sha(sources.VIEWS), Path(__file__).name: sha(Path(__file__).name)})
        self.geometry_theorem = exact_geometry_theorem()
        if set(self.views) != set(exact_radius_maps()[0]):
            raise ValueError('All 17 original function graphs required')
        self.service.bind_hashes(self.hashes)

    def build(self):
        g = FunctionTransportGraph()
        maps, symbols, x = exact_radius_maps()
        bindings = {
            'logP': 'same native seed.logP, exact logP=exp(Md)+11',
            'logC': 'same selected physical_norm_family logCstar',
            'T': 'same actual_long_reshape_mixed_C4 T=400*Abar',
            'Tw': 'same native seed.params.Tw analytic parameter',
            'hbB': 'exp(same native bridge.logh)', 'hbS': 'exp(same native switch.logh)',
            'sc': 'current_inner_exit_strict_collar selected_first_phase_endpoint',
            'Md': 'same native seed.params.Md',
            'mu': 'current_generic_shear_loop_domain right_collar_mu; same checked repair mu'}
        parameters = {name: g.node('original_source_parameter', name=name, exact_source=source)
            for name, source in bindings.items()}
        N = g.node('shared_positive_integer', name='N', lower=current.MIN_N)
        history = {key: C1Function(g.zero, g.zero) for key in RATES}
        cells = []

        def source(chart, node, coordinate, offset):
            phi = g.unary('fractional_part', g.mul(N, offset))
            return g.node('original_function_graph', chart=chart, graph_file=sources.VIEWS,
                graph_sha256=self.hashes[sources.VIEWS], source_node=node, coordinate=coordinate.node,
                phase=phi.node, shared_N=N.node, Z_variable='Z',
                source_coefficients_are_function_recipes_not_cover_values=True)

        for index, (label, chart, left, right) in enumerate(exact_route_endpoints(current.ROUTE, symbols)):
            t = g.symbol('coordinate_' + str(index))
            lo, hi = expression(g, left, parameters), expression(g, right, parameters)
            offset = expression(g, maps[chart], parameters, t)
            yleft = expression(g, sy.simplify(maps[chart].subs(x, left)), parameters)
            yright = expression(g, sy.simplify(maps[chart].subs(x, right)), parameters)
            # Collect width symbolically before conversion; microscopic terms
            # survive even when huge absolute endpoints would round together.
            width = expression(g, sy.simplify(maps[chart].subs(x, right)-maps[chart].subs(x, left)), parameters)
            jacobian = expression(g, sy.diff(maps[chart], x), parameters, t)
            view = self.views[chart]
            incoming = history
            contributions, outgoing = {}, {}
            for key, rate in RATES.items():
                decay = g.unary('exp', g.neg(g.mul(g.constant(rate), width)))
                kernel = g.unary('exp', g.neg(g.mul(g.constant(rate), g.sub(yright, offset))))
                if label == 'initial_flat_collar' or chart == 'O3_power':
                    # Proven original whole flat supports; history is still
                    # propagated. Original pressure/velocity is never zeroed.
                    contribution = C1Function(g.zero, g.zero)
                else:
                    roots = (view['five_signed_increment_rate_roots'][key],
                        view['five_signed_increment_rate_first_derivatives']['Z'][key])
                    pair = []
                    for root in roots:
                        full_density = source(chart, root, t, offset)
                        integrand = g.mul(kernel, full_density, jacobian)
                        pair.append(g.node('definite_integral', integrand=integrand.node,
                            variable='coordinate_' + str(index), lower=lo.node, upper=hi.node,
                            measure='native coordinate; original dy/dcoordinate applied exactly once',
                            exact_function_integral=True, numerical_value_not_installed=True))
                    contribution = C1Function(*pair)
                contributions[key] = contribution
                outgoing[key] = g.c1add(g.c1scale(decay, incoming[key]), contribution)
            history = outgoing
            cells.append(dict(label=label, chart=chart, lower=lo.node, upper=hi.node,
                left_radius_offset=yleft.node, right_radius_offset=yright.node, width=width.node,
                incoming={key: dict(value=q.value.node, Z=q.Z.node) for key, q in incoming.items()},
                contributions={key: dict(value=q.value.node, Z=q.Z.node) for key, q in contributions.items()},
                outgoing={key: dict(value=q.value.node, Z=q.Z.node) for key, q in outgoing.items()},
                source_flat_exact_zero=label == 'initial_flat_collar' or chart == 'O3_power',
                all_incoming_histories_retained=True, original_P0_and_P0_Z_unchanged=True))
        endpoint = expression(g, 2/symbols['Tw'], parameters)
        endpoint_offset = expression(g, maps['O3_power'].subs(x, 2/symbols['Tw']), parameters)
        original = self.views['O3_power']['original_signed_input_graph']['jet_expression_dag']['roots']['E']
        A = C1Function(source('O3_power', original['y0_Z0'], endpoint, endpoint_offset),
            source('O3_power', original['y0_Z1'], endpoint, endpoint_offset))
        ratio = g.quotient(A.Z, A.value, 'same original Rc Ac/S positive theorem')
        den = g.mul(A.value, A.value)
        targets = {}
        for row, key, degree in (('M', 'm', 1), ('I', 'h', 1), ('S', 'e', 2), ('Cp', 'p', 2)):
            divisor = A.value if degree == 1 else den
            value = g.quotient(history[key].value, divisor, 'same original Rc amplitude power positive')
            jet = g.sub(g.quotient(history[key].Z, divisor, 'same original Rc amplitude power positive'),
                g.mul(g.constant(degree), ratio, value))
            targets[row] = C1Function(value, jet)
        numerator = C1Function(g.sub(history['k'].value, g.mul(A.value, history['m'].value)),
            g.sub(history['k'].Z, g.add(g.mul(A.Z, history['m'].value), g.mul(A.value, history['m'].Z))))
        divisor = g.mul(parameters['mu'], den)
        value = g.quotient(numerator.value, divisor, 'same original positive mu and Rc amplitude')
        jet = g.quotient(g.sub(numerator.Z, g.mul(g.constant(2), ratio, numerator.value)), divisor,
            'same original positive mu and Rc amplitude')
        targets[ROWS[1]] = C1Function(value, jet)
        targets = {key: targets[key] for key in ROWS}
        scaled = {key: g.c1scale(N, q) for key, q in targets.items()}
        return dict(graph=g, cells=cells, history=history, amplitude=A, numerator=numerator,
            targets=targets, N_scaled_targets=scaled, N=N, parameters=parameters,
            source_family=self.family, source_graph_sha256=self.hashes[sources.VIEWS])

    def record(self, built):
        encode = lambda pairs: {key: dict(value=q.value.node, Z=q.Z.node) for key, q in pairs.items()}
        return dict(source_family=self.family, **{GATE: True}, function_graph_nodes=built['graph'].nodes,
            exact_original_cells=built['cells'], exact_Rc_history_roots=encode(built['history']),
            exact_Rc_target_roots=encode(built['targets']), exact_N_scaled_target_roots=encode(built['N_scaled_targets']),
            exact_Rc_amplitude_roots=dict(value=built['amplitude'].value.node, Z=built['amplitude'].Z.node),
            exact_joint_divided_numerator_roots=dict(value=built['numerator'].value.node, Z=built['numerator'].Z.node),
            exact_geometry_theorem=self.geometry_theorem, common_N_lower=current.MIN_N,
            original_zero_correction_inlet=True, original_P0_and_P0_Z_unchanged=True,
            exact_N_dependent_integral_functions_defined=True, numerical_original_source_oracle_installed=False,
            actual_original_numerical_integrals_evaluated=False, actual_five_controls_installed=False,
            actual_terminal_Z_function_closure_installed=False, current_whole_N_selected=False,
            **dict.fromkeys(current.packets.OPEN, False), input_hashes=self.hashes,
            scope='Exact integral representations of 24-cell/17-chart original full signed N-dependent source functions and first-Z histories, joint Rc targets and N*r. Function definitions only; no cap coefficients selected as function values, original point oracle/controls/closure/global N/higher jets/recursion remain open.')

    def run(self):
        began = time.monotonic()
        result = self.record(self.build())
        result['execution_seconds'] = time.monotonic()-began
        (HERE/NAME).write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
        print('Exact original Rc source-function integrals and targets defined: 24 cells / 17 charts', flush=True)
        return result


class SourceOracleRequired(RuntimeError):
    pass


def evaluate(built, root, *, oracle, Z, N, ctx):
    """Evaluate with explicitly injected source functions and quadrature.

    No implicit legacy constructor, saved cover midpoint or fallback zero.
    Quadrature is approximate unless the injected oracle certifies it.
    """
    graph = built['graph']
    graph.ids([root])
    if type(N) is not int or N < current.MIN_N:
        raise ValueError('One common integer N>=160 required')
    if oracle is None or not all(callable(getattr(oracle, key, None)) for key in ('parameter', 'source', 'integrate')):
        raise SourceOracleRequired('Original parameter/source/integral oracle required; saved caps cannot supply values')
    if getattr(oracle, 'source_family', None) != built['source_family']:
        raise SourceOracleRequired('Oracle must bind the same original source family')
    mode = getattr(oracle, 'mode', None)
    if mode not in ('original_source', 'synthetic_reference'):
        raise SourceOracleRequired('Oracle must explicitly distinguish original source and synthetic reference')
    if mode == 'original_source' and getattr(oracle, 'source_graph_sha256', None) != built['source_graph_sha256']:
        raise SourceOracleRequired('Original numerical oracle must bind the accepted function graph hash')
    cache = {}

    def walk(index, variables):
        key = (index, tuple(sorted(variables.items())))
        if key in cache:
            return cache[key]
        row = graph.nodes[index];op = row['operation']
        child = lambda i: walk(i, variables)
        if op == 'exact_rational': value = ctx.mpf(row['numerator'])/row['denominator']
        elif op == 'bound_variable': value = variables[row['name']]
        elif op == 'original_source_parameter': value = oracle.parameter(row['name'])
        elif op == 'shared_positive_integer': value = ctx.mpf(N)
        elif op == 'sum': value = sum((child(i) for i in row['arguments']), ctx.mpf(0))
        elif op == 'product':
            value = ctx.mpf(1)
            for i in row['arguments']: value *= child(i)
        elif op == 'negative': value = -child(row['argument'])
        elif op == 'positive_quotient':
            denominator = child(row['denominator'])
            if denominator <= 0: raise ArithmeticError('Oracle violates original positive source denominator')
            value = child(row['numerator'])/denominator
        elif op == 'analytic_unary':
            a = child(row['argument'])
            value = a-ctx.floor(a) if row['name'] == 'fractional_part' else getattr(ctx, row['name'])(a)
        elif op == 'original_function_graph':
            value = oracle.source(row, coordinate=child(row['coordinate']), phase=child(row['phase']), Z=Z, N=N)
        elif op == 'definite_integral':
            def integrand(t):
                return walk(row['integrand'], {**variables, row['variable']: t})
            value = oracle.integrate(integrand, child(row['lower']), child(row['upper']))
        else: raise ValueError('Unknown exact function operation: ' + op)
        if isinstance(value, (dict, current.repair.LogUpper)):
            raise TypeError('Oracle returned a saved cap/record instead of a function value')
        cache[key] = value
        return value
    return walk(root.node, {})


def run(owner):
    return NativeRcFunctionTransport(owner).run()
