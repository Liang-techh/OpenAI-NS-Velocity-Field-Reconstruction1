"""Retain the original mu in pulse axial curl before interval subtraction."""
import ast
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_differential_field as differential
import lei_ren_part1_paper_compliant_current_original_Rp_centered_scale_arithmetic as centered
import lei_ren_part1_paper_compliant_pulse_mixed_C4 as pulse_mixed
from lei_ren_part1_paper_compliant_pulse_high_jets import CompliantPulseHighJets
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE, PREFIX, sha, ends = differential.HERE, differential.PREFIX, differential.sha, differential.ends
signed, correlated, box = differential.signed, differential.correlated, differential.box
NAME = PREFIX + 'current_original_Rp_axial_vorticity.json.gz'
RECEIPT = PREFIX + 'current_original_Rp_axial_vorticity_check.json'
GATES = ('current_original_Rp_actual_pulse_theta_radial_law_bound',
         'current_original_Rp_nonzero_signed_axial_vorticity_point_installed')
OPEN = differential.OPEN
MU, DELTA = s.symbols('mu delta', positive=True)


def source_law():
    name = PREFIX + 'pulse_mixed_C4.py'
    tree = ast.parse((HERE / name).read_bytes())
    fn = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'transport_mixed')
    assignments = {target.id: node.value for node in fn.body if isinstance(node, ast.Assign)
        for target in node.targets if isinstance(target, ast.Name)}
    expected = ast.parse("c.mpf('.5')+mu", mode='eval').body
    if ast.dump(assignments['bp']) != ast.dump(expected):
        raise ValueError('Original theta radial rate differs')
    physical = assignments['physical']
    theta = next(value for key, value in zip(physical.keys, physical.values)
        if isinstance(key, ast.Constant) and key.value == correlated.physical.mixed.UT)
    expected = ast.parse('[u*(-bp)**k for k in range(5)]', mode='eval').body
    if ast.dump(theta) != ast.dump(expected):
        raise ValueError('Actual original full theta derivative source differs')
    high = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'CompliantPulseMixedC4')
    high = next(node for node in high.body if isinstance(node, ast.FunctionDef) and node.name == '_high_packet')
    calls = [node for node in ast.walk(high) if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name) and node.func.id == 'transport_mixed']
    if len(calls) != 1 or ast.dump(calls[0]) != ast.dump(ast.parse('transport_mixed(self,Z,point,Brows,u)', mode='eval').body):
        raise ValueError('Actual pulse packet does not consume the original theta source')
    radial_name = PREFIX + 'pulse_radial_C4.py'
    radial = ast.parse((HERE / radial_name).read_bytes())
    radial = next(node for node in radial.body if isinstance(node, ast.ClassDef) and node.name == 'CompliantPulseRadialC4')
    data = next(node for node in radial.body if isinstance(node, ast.FunctionDef) and node.name == 'data')
    rows = [node.value for node in ast.walk(data) if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == 'u' for target in node.targets)]
    if len(rows) != 1 or ast.dump(rows[0]) != ast.dump(ast.parse("r*k['U']", mode='eval').body):
        raise ValueError('Original inlet theta coefficient source differs')
    inlet = [node.value for node in ast.walk(high) if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Tuple) and [getattr(v, 'id', None) for v in target.elts] ==
            ['selected', '_', 'u', '_', '_'] for target in node.targets)]
    if len(inlet) != 1 or ast.dump(inlet[0]) != ast.dump(ast.parse('self.data(Z)', mode='eval').body):
        raise ValueError('Actual full theta derivative coefficient must depend only on Z')
    return dict(source_file=name, source_sha256=sha(name), rate_AST=ast.dump(assignments['bp']),
        theta_derivative_AST=ast.dump(theta), source_call_AST=ast.dump(calls[0]),
        inlet_source_file=radial_name, inlet_source_sha256=sha(radial_name),
        inlet_coefficient_AST=ast.dump(rows[0]), inlet_call_AST=ast.dump(inlet[0]),
        identity='Utheta_y=-(1/2+mu)*Utheta', coefficient_independent_of_logR=True)


def operator_expression(graph, node, mu, delta):
    node = node.node if hasattr(node, 'node') else node
    if node == mu.node:
        return MU
    if node == delta.node:
        return DELTA
    data = graph.nodes[node]; op = data['operation']
    rec = lambda value: operator_expression(graph, value, mu, delta)
    if op == 'exact_rational':
        return s.Rational(data['numerator'], data['denominator'])
    if op == 'sum':
        return s.Add(*(rec(value) for value in data['arguments']))
    if op == 'product':
        return s.Mul(*(rec(value) for value in data['arguments']))
    if op == 'negative':
        return -rec(data['argument'])
    if op == 'positive_quotient':
        return rec(data['numerator']) / rec(data['denominator'])
    if op == 'exact_operator_integer_power':
        return rec(data['base']) ** data['exponent']
    if op == 'analytic_unary' and data['name'] in ('sin', 'cos', 'log', 'exp'):
        return getattr(s, data['name'])(rec(data['argument']))
    raise ValueError('Closed original Cartesian operator expression required')


class CurrentOriginalRpAxialVorticity:
    @source_precision
    def __init__(self, before, require_checked=True):
        if type(before) is not differential.CurrentOriginalRpDifferentialField or not before.acceptance_loaded:
            raise ValueError('Accepted original differential field owner required')
        self.before, self.arithmetic, self.product = before, before.before, before.product
        self.graph, self.ctx, self.family_record = before.graph, before.ctx, before.family_record
        self.refined = self.product.before.refined
        self.mu, self.delta = self.product.amplitude.mu, before.physical.delta_function
        self.law = source_law()
        self.hashes = dict(before.hashes)
        for name in (differential.NAME, differential.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        self.acceptance_loaded, self._issued = False, {}
        self.assert_graph()
        if require_checked:
            receipt = json.loads((HERE / RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family'] != self.family_record or \
                    receipt['actual_theta_source_law'] != self.law:
                raise ValueError('Original axial vorticity receipt/source/scope differs')
            for name in (Path(__file__).name, Path(__file__).stem + '_check.py', NAME):
                if receipt['input_hashes'].get(name) != sha(name):
                    raise ValueError('Unbound axial vorticity source ' + name)
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes, receipt['input_hashes'])
            self.acceptance_loaded = True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        native = self.refined.correlated.before
        source_mu = self.product.amplitude.reader().bindings[self.mu.node]
        checks = dict(original_theta_AST_unchanged=source_law() == self.law,
            actual_refined_pulse_class=type(self.refined.pulse) is pulse_mixed.CompliantPulseMixedC4,
            actual_refined_main_method=self.refined.pulse.main.__func__ is CompliantPulseHighJets.main,
            actual_refined_inlet_method=self.refined.pulse.data.__func__ is CompliantPulseRadialC4.data,
            actual_refined_packet_method=self.refined.pulse._high_packet.__func__ is pulse_mixed.CompliantPulseMixedC4._high_packet,
            actual_theta_source_callable=pulse_mixed.CompliantPulseMixedC4._high_packet.__globals__['transport_mixed'] is pulse_mixed.transport_mixed,
            actual_refined_raw_callback=self.refined.raw.pulse_view.__func__ is native._adapter_functions['pulse_view'],
            actual_refined_velocity_callback=self.refined.transport.velocity_rows.__func__ is native._adapter_functions['heat_velocity_rows'],
            actual_refined_dispatch_owner=self.refined.raw.selected is self.refined.dispatch
                and self.refined.dispatch.evaluate.__self__ is self.refined
                and self.refined.dispatch.evaluate.__func__ is type(self.refined)._dispatch,
            actual_source_mu_enclosure=self.refined.pulse.mu._mpi_ == source_mu._mpi_
                and self.refined.pulse.mu.ctx is self.refined.ctx,
            actual_source_delta_enclosure=self.refined.pulse.delta._mpi_ == self.before.physical.delta._mpi_,
            same_original_mu_function=self.mu is self.product.amplitude.mu,
            same_current_source_family=self.family_record == self.refined.family_record)
        if not all(checks.values()):
            raise ValueError('Original axial vorticity source differs: ' + str(checks))
        return checks

    @source_precision
    def evaluate(self, field):
        self.assert_graph()
        record = self.before._require(field)
        if record['chart'] != 'pulse_exit':
            raise ValueError('This admitted axial curl source law requires pulse_exit')
        reader = self.before._reader(record)
        primitives = {}
        for component, index, weight in differential.standard_forms()['vorticity_z']:
            for term in record['entries'][component, index]['raw'].terms:
                log = self.graph.add(*(ref for _, ref in term.log_scale_parts))
                key = (term.source_label, term.source_row.derivative, term.source_row.source_units,
                    term.source_row.powers, tuple(sorted(reader.polynomial(log).items())))
                part = primitives.setdefault(key, dict(row=term.source_row, log=log, operators=[]))
                if part['row'].coefficients[0]._mpi_ != term.source_row.coefficients[0]._mpi_:
                    raise ValueError('Actual common theta primitive required')
                part['operators'].append(self.graph.mul(self.graph.constant(weight), term.operator_function))
        theta, cancelled = {}, []
        for key, part in primitives.items():
            operator = self.graph.add(*part['operators'])
            if not reader.polynomial(operator):
                cancelled.append(key[0]); continue
            if key[0] != correlated.physical.mixed.UT or key[1] not in ((0, 0), (1, 0)) or key[1] in theta:
                raise ValueError('Only the original base/first-radial theta primitives may survive')
            part['operator'] = operator
            theta[key[1]] = (key, part)
        if len(cancelled) != 2 or set(cancelled) != {correlated.physical.mixed.UR} or set(theta) != {(0, 0), (1, 0)}:
            raise ValueError('Original Cartesian radial cancellation and theta pair required')
        key0, base = theta[0, 0]; key1, radial = theta[1, 0]
        if key0[2:] != key1[2:]:
            raise ValueError('Theta radial law must preserve source units, powers and full log scale')
        a = operator_expression(self.graph, base['operator'], self.mu, self.delta)
        b = operator_expression(self.graph, radial['operator'], self.mu, self.delta)
        if s.trigsimp(s.expand(a - (s.Rational(1, 2) + MU) * b + MU * s.sqrt(2))) != 0:
            raise ValueError('Original physical axial curl reduction differs')
        coefficient = -self.ctx.sqrt(2) * self.ctx.mpf(base['row'].coefficients[0])
        if ends(coefficient)[1] >= 0:
            raise ValueError('Original theta coefficient must be strictly positive')
        scale = self.graph.add(base['log'], self.product.logmu)
        result = signed.enclose_factored_sum(self.graph, self.ctx, reader,
            [dict(log_function=scale, log_bound=reader.at(scale), coefficient=coefficient, proofs=[])], record['target'])
        centered_scale = self.product.center(self.product.reader(record['delivery']), scale)
        accuracy = centered.centered_relative_budget(self.ctx, result, record['target'], centered_scale)
        result.update(component='omega_z', physical_derivative=(0, 0, 0), physical_accuracy=accuracy)
        value_record = dict(component='omega_z', unit='velocity/length**1', scale=scale,
            scale_polynomial=reader.polynomial(scale), scale_dag=self.arithmetic._scale_dag(scale),
            coefficient_tuple=coefficient._mpi_, sign=-1, root=object(), multiplier=Fraction(1),
            delivery=record['delivery'], bindings=dict(record['bindings']), binding_snapshot=dict(record['binding_snapshot']),
            aliases=dict(record['aliases']), alias_snapshot=dict(record['alias_snapshot']),
            allowed=record['allowed'], allowed_snapshot=record['allowed_snapshot'],
            absolute_width_satisfied=accuracy['ordinary_numeric_relative_width_satisfied'],
            original_materialized=False, original_row=result)
        value = self.arithmetic._issue(value_record)
        proof = self.graph.node('original_pulse_exit_exact_axial_curl_mu_reduction',
            original_base_theta_source=base['row'].name, base_operator=base['operator'].node,
            radial_operator=radial['operator'].node, exact_original_mu=self.mu.node,
            exact_retained_logmu=self.product.logmu.node, original_theta_source_AST=self.law,
            identity='omega_z=-mu*sqrt(2/R)*lambda**(-2-delta)*Utheta',
            numeric_subtraction_of_near_equal_theta_derivatives_not_used=True)
        self._issued[id(value)] = (value, field, self.before._signature(value_record), dict(
            physical_accuracy=accuracy, exact_source_law_proof=proof.node,
            defining_proof=copy.deepcopy(self.graph.nodes[proof.node]),
            original_base_theta_coefficient=base['row'].coefficients[0],
            original_theta_pair_common_log_scale=base['log'].node,
            exact_radial_primitive_cancellations=2))
        return value

    @source_precision
    def report(self, value):
        self.assert_graph()
        entry = self._issued.get(id(value))
        if entry is None or entry[0] is not value:
            raise ValueError('Axial vorticity value issued by this owner required')
        self.before._require(entry[1])
        issued = self.arithmetic._require(value)
        if self.before._signature(issued) != entry[2] or self.graph.nodes[entry[3]['exact_source_law_proof']] != entry[3]['defining_proof']:
            raise ValueError('Original axial vorticity value/proof changed')
        details = {key: copy.deepcopy(val) for key, val in entry[3].items() if key != 'defining_proof'}
        return dict(value=self.arithmetic.export(value), actual_theta_source_law=self.law, **details,
            original_mu_retained_before_enclosure=True, near_equal_derivatives_not_subtracted=True,
            actual_time_growth_not_measured=True, unrestricted_physical_point_API=False,
            full_certified_physical_accuracy=False, **dict.fromkeys(GATES, self.acceptance_loaded), **dict.fromkeys(OPEN, False))


@source_precision
def run(before, fields):
    began = time.monotonic()
    owner = CurrentOriginalRpAxialVorticity(before, require_checked=False)
    values = {name: owner.evaluate(field) for name, field in fields.items()}
    views = {name: owner.report(value) for name, value in values.items()}
    result = dict(source_family=owner.family_record, actual_axial_vorticity=views,
        source_assertions=owner.assert_graph(), input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began, **dict.fromkeys(GATES + OPEN, False))
    (HERE / NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result), separators=(',', ':'))+'\n').encode(), mtime=0))
    return owner, values, views
