"""Bounded, source-preserving partial resolver for Q4 endpoint ports.

Exact scalars are stored as a rational coefficient times 2**exponent.  Huge
exponents are never expanded. Unsupported source operations remain opaque;
an opaque source is not a value, enclosure, zero, or certified error bound.
"""
from fractions import Fraction
import gzip
import hashlib
import json

SCHEMA = 'original_q4_endpoint_resolver_v1'
MAX_NODES = 1500
MAX_DEPTH = 96
MAX_SHIFT = 4096
MAX_BITS = 8192


class Invalid(ValueError):
    pass


def exact(n, d=1, e=0):
    if any(type(v) is not int for v in (n, d, e)) or d <= 0:
        raise Invalid('invalid exact rational')
    if abs(n).bit_length() > MAX_BITS or d.bit_length() > MAX_BITS:
        raise Invalid('coefficient size limit')
    f = Fraction(n, d)
    if f.numerator == 0:
        return {'kind': 'exact_binary_rational', 'numerator': 0, 'denominator': 1, 'exponent': 0}
    return {'kind': 'exact_binary_rational', 'numerator': f.numerator,
            'denominator': f.denominator, 'exponent': e}


def dyadic(value):
    if not isinstance(value, (list, tuple)) or len(value) != 4 or any(type(x) is not int for x in value):
        raise Invalid('invalid dyadic')
    sign, mantissa, exponent, bits = value
    if sign not in (0, 1) or mantissa < 0 or bits != mantissa.bit_length():
        raise Invalid('noncanonical dyadic')
    if mantissa == 0:
        if tuple(value) != (0, 0, 0, 0):
            raise Invalid('noncanonical zero')
        return exact(0)
    if mantissa % 2 == 0:
        raise Invalid('noncanonical mantissa')
    return exact(-mantissa if sign else mantissa, 1, exponent)


def is_exact(obj):
    return isinstance(obj, dict) and obj.get('kind') == 'exact_binary_rational'


def mul(a, b):
    return exact(a['numerator'] * b['numerator'], a['denominator'] * b['denominator'],
                 a['exponent'] + b['exponent'])


def add(a, b):
    e = min(a['exponent'], b['exponent'])
    x, y = a['exponent'] - e, b['exponent'] - e
    if max(x, y) > MAX_SHIFT:
        return None
    return exact((a['numerator'] << x) * b['denominator'] + (b['numerator'] << y) * a['denominator'],
                 a['denominator'] * b['denominator'], e)


def compare(a, b):
    """Exact compare without shifting by an astronomical exponent difference."""
    x = a['numerator'] * b['denominator']
    y = b['numerator'] * a['denominator']
    if x == 0 or y == 0 or (x < 0) != (y < 0):
        return (x > y) - (x < y)
    negative = x < 0
    x, y = abs(x), abs(y)
    lx, ly = x.bit_length() + a['exponent'], y.bit_length() + b['exponent']
    if lx != ly:
        order = (lx > ly) - (lx < ly)
    else:
        delta = a['exponent'] - b['exponent']
        if abs(delta) > MAX_SHIFT:
            raise Invalid('comparison shift limit')
        if delta >= 0:
            x <<= delta
        else:
            y <<= -delta
        order = (x > y) - (x < y)
    return -order if negative else order


def validate_port(port):
    if not isinstance(port, dict) or port.get('schema') != 'q4_physical_Uy_endpoint_source_port_v1':
        raise Invalid('invalid port schema')
    for name in ('expression_node', 'parameter_node', 'coordinate_node', 'ordinary_Z_order'):
        if type(port.get(name)) is not int or port[name] < 0:
            raise Invalid('invalid ' + name)
    if port['ordinary_Z_order'] not in (0, 1, 2):
        raise Invalid('unsupported ordinary Z order')
    h = port.get('archive_sha256')
    if not isinstance(h, str) or len(h) != 64 or any(c not in '0123456789abcdef' for c in h):
        raise Invalid('invalid archive sha')
    sc = port.get('parameter_exact_mpf_tuple')
    scalar = dyadic(sc)
    if scalar['numerator'] <= 0:
        raise Invalid('sc must be positive')
    half = port.get('coordinate_exact_mpf_tuple')
    dyadic(half)
    if half != [sc[0], sc[1], sc[2] - 1, sc[3]]:
        raise Invalid('endpoint is not sc/2')
    if port['parameter_node'] == port['coordinate_node']:
        raise Invalid('parameter and coordinate node collide')


class Resolver:
    def __init__(self, graph, port, env):
        validate_port(port)
        nodes = graph.get('exact_graph_nodes') if isinstance(graph, dict) else None
        if not isinstance(nodes, list):
            raise Invalid('missing graph nodes')
        if not isinstance(env, dict) or any(type(k) is not int for k in env):
            raise Invalid('environment must map integer nodes to exact dyadics')
        self.nodes, self.port = nodes, port
        sc_ref, coord_ref = port['parameter_node'], port['coordinate_node']
        if any(ref >= len(nodes) for ref in (sc_ref, coord_ref, port['expression_node'])):
            raise Invalid('port node out of range')
        sc_row = nodes[sc_ref]
        if not isinstance(sc_row, dict) or sc_row.get('operation') != 'current_original_source_parameter' or sc_row.get('name') != 'sc' or sc_row.get('Z_independent') is not True:
            raise Invalid('port sc is not the source parameter')
        for key in env:
            if key < 0 or key >= len(nodes) or not isinstance(nodes[key], dict):
                raise Invalid('environment node out of range')
            if nodes[key].get('operation') not in ('symbol', 'bound_variable', 'current_original_source_parameter'):
                raise Invalid('environment may bind source variables only')
        self.bindings = {key: dyadic(value) for key, value in env.items()}
        for key, value in ((port['parameter_node'], port['parameter_exact_mpf_tuple']),
                           (port['coordinate_node'], port['coordinate_exact_mpf_tuple'])):
            chosen = dyadic(value)
            if key in self.bindings and compare(self.bindings[key], chosen) != 0:
                raise Invalid('environment conflicts with endpoint port')
            self.bindings[key] = chosen
        self.memo, self.active = {}, set()
        self.visited = 0
        self.opaque = []
        self.flat_branches = []

    def unknown(self, ref, operation, reason):
        entry = {'kind': 'opaque_source_ref', 'node': ref, 'operation': operation, 'reason': reason}
        self.opaque.append(entry)
        return entry

    def get(self, ref, depth=0):
        if type(ref) is not int or ref < 0 or ref >= len(self.nodes):
            raise Invalid('invalid graph reference')
        if ref in self.bindings:
            return self.bindings[ref]
        if ref in self.memo:
            return self.memo[ref]
        if ref in self.active:
            raise Invalid('cycle in source graph')
        row = self.nodes[ref]
        if not isinstance(row, dict):
            raise Invalid('invalid graph node')
        op = row.get('operation')
        if self.visited >= MAX_NODES or depth >= MAX_DEPTH:
            return self.unknown(ref, op, 'traversal budget')
        self.visited += 1
        self.active.add(ref)
        try:
            result = self.interpret(ref, row, depth)
        finally:
            self.active.remove(ref)
        self.memo[ref] = result
        return result

    def interpret(self, ref, row, depth):
        op = row.get('operation')
        child = lambda key: self.get(row[key], depth + 1)
        if op == 'exact_rational':
            return exact(row.get('numerator'), row.get('denominator'))
        if op in ('sum', 'product'):
            args = row.get('arguments')
            if not isinstance(args, list) or not args:
                raise Invalid('missing arithmetic arguments')
            # Preserve source order in the result while visiting the later
            # terms first: the Q4 loop correction branches then have a chance
            # to resolve before a large unresolved background exhausts the
            # budget. This changes no arithmetic ordering or field value.
            values = list(reversed([self.get(i, depth + 1) for i in reversed(args)]))
            if all(is_exact(v) for v in values):
                acc = exact(0 if op == 'sum' else 1)
                for value in values:
                    acc = add(acc, value) if op == 'sum' else mul(acc, value)
                    if acc is None:
                        break
                if acc is not None:
                    return acc
            return {'kind': 'exact_source_expression', 'operation': op, 'node': ref, 'arguments': values}
        if op == 'positive_quotient':
            numerator, denominator = child('numerator'), child('denominator')
            if is_exact(denominator) and denominator['numerator'] <= 0:
                raise Invalid('positive quotient denominator not positive')
            if is_exact(numerator) and is_exact(denominator):
                return exact(numerator['numerator'] * denominator['denominator'],
                             numerator['denominator'] * denominator['numerator'],
                             numerator['exponent'] - denominator['exponent'])
            return {'kind': 'exact_source_expression', 'operation': op, 'node': ref,
                    'numerator': numerator, 'denominator': denominator,
                    'requires_original_positive_denominator': True}
        if op == 'exact_real_comparison':
            left, right = child('left'), child('right')
            operator = row.get('operator')
            if operator not in ('<', '<=', '>', '>=', '=='):
                return self.unknown(ref, op, 'unsupported comparison')
            if not is_exact(left) or not is_exact(right):
                return self.unknown(ref, op, 'comparison operands unresolved')
            c = compare(left, right)
            return {'kind': 'exact_boolean', 'value': {'<': c < 0, '<=': c <= 0,
                    '>': c > 0, '>=': c >= 0, '==': c == 0}[operator]}
        if op == 'logical_or':
            args = row.get('arguments')
            if not isinstance(args, list) or not args:
                raise Invalid('missing logical arguments')
            undecided = False
            # A coordinate comparison can decide the collar without entering
            # the unrelated shear/active expression. Reordering pure OR terms
            # does not evaluate a selected lazy branch's inactive body.
            def coordinate_first(item):
                if type(item) is not int or item < 0 or item >= len(self.nodes):
                    return 1
                candidate = self.nodes[item]
                return 0 if isinstance(candidate, dict) and candidate.get('operation') == 'exact_real_comparison' and candidate.get('left') == self.port['coordinate_node'] else 1
            for item in sorted(args, key=coordinate_first):
                value = self.get(item, depth + 1)
                if value.get('kind') == 'exact_boolean' and value.get('value') is True:
                    return value
                if value.get('kind') != 'exact_boolean':
                    undecided = True
            return self.unknown(ref, op, 'logical operands unresolved') if undecided else {'kind': 'exact_boolean', 'value': False}
        if op == 'original_lazy_flat_branch':
            condition = child('flat_predicate')
            if condition.get('kind') != 'exact_boolean':
                return self.unknown(ref, op, 'flat predicate unresolved; neither body evaluated')
            branch = 'flat_value' if condition['value'] else 'active_body'
            if branch == 'flat_value':
                self.flat_branches.append({'node': ref, 'inactive_body_not_visited': row.get('active_body')})
            return child(branch)
        return self.unknown(ref, op, 'source operation not implemented')


def resolve(graph, port, env):
    """Resolve a caller-supplied graph; archive authentication is external."""
    try:
        resolver = Resolver(graph, port, env)
        expression = resolver.get(port['expression_node'])
        status = 'resolved_exact' if is_exact(expression) else 'unresolved'
        return {'schema': SCHEMA, 'status': status, 'expression': expression,
                'visited_node_count': resolver.visited, 'opaque_observations': resolver.opaque,
                'flat_branches': resolver.flat_branches,
                'scope': {'archive_authenticated': False, 'numeric_field_value_computed': status == 'resolved_exact',
                          'exact_scalar_expression_only': status == 'resolved_exact',
                          'other_source_assumptions_required': True, 'old_errors_removed': False}}
    except (Invalid, KeyError, TypeError, OverflowError) as exc:
        return {'schema': SCHEMA, 'status': 'mismatch', 'error': str(exc),
                'scope': {'archive_authenticated': False, 'numeric_field_value_computed': False,
                          'old_errors_removed': False}}


def resolve_authenticated(raw_archive, ports, env):
    """Authenticate raw gzip and resolve ports authenticated by the caller."""
    if not isinstance(raw_archive, bytes) or not isinstance(ports, list) or not ports:
        raise Invalid('missing raw archive or ports')
    digest = hashlib.sha256(raw_archive).hexdigest()
    for port in ports:
        validate_port(port)
        if port['archive_sha256'] != digest:
            raise Invalid('archive SHA mismatch')
    graph = json.loads(gzip.decompress(raw_archive))
    results = [resolve(graph, port, env) for port in ports]
    for item in results:
        item['scope']['archive_authenticated'] = True
    return {'schema': SCHEMA, 'archive_sha256': digest, 'raw_sha_checked_before_decode': True,
            'port_authentication_required_of_caller': True, 'results': results, 'numeric_field_values_computed': any(
                item['scope'].get('numeric_field_value_computed') is True for item in results),
            'old_errors_removed': False}
