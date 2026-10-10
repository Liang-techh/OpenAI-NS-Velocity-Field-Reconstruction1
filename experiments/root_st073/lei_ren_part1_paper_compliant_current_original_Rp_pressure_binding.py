"""Bind current analytic P0 Taylor rows to both actual pressure datum owners.

The inlet and selected/future datums are distinct live objects. Both enclose
the same fourteen-atom analytic source. Their mass/flatten boxes are bounds,
never definitions of an analytic function. Physical pressure stays factored.
"""
import gzip
import json
import math
from pathlib import Path
import time

import sympy as s
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
from lei_ren_part1_paper_logarithmic_pressure_datum import BETA0, BETA2, LogarithmicPressureDatum
from lei_ren_part1_paper_candidate_pressure_function import q_jets
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
import lei_ren_part1_paper_compliant_current_original_Rp_signed_physical_values as signed

amplitude, correlated = signed.amplitude, signed.correlated
HERE, PREFIX, sha = signed.HERE, signed.PREFIX, signed.sha
NAME = PREFIX+'current_original_Rp_pressure_binding.json.gz'
RECEIPT = PREFIX+'current_original_Rp_pressure_binding_check.json'
GATES = ('current_original_Rp_analytic_P0_six_row_numeric_source_inclusion_installed',
         'current_original_Rp_inlet_and_selected_P0_live_owners_bound',
         'current_original_Rp_P0_lazy_original_pressure_scale_bound')
OPEN = signed.OPEN
ends = correlated.ends


def datum_snapshot(datum):
    """Capture defining recipe and every numerical enclosure, not just labels."""
    return (json.dumps(datum.definition, sort_keys=True), datum.source_sha, datum.datum_sha,
            tuple(sorted(datum.input_hashes.items())), datum.ctx.dps,
            tuple((name, row['beta'], row['mass']._mpi_) for name, row in sorted(datum.stages.items())),
            tuple(getattr(datum, key)._mpi_ for key in
                  ('m0', 'm2', 'tail_upper', 'rho', 'flatten_complex_upper')),
            datum.parameters.Md, datum.parameters.logPstar._mpi_)


def exact_Z_box(ctx, Z):
    """Coordinate interval endpoints are exact rationals, never value oracles."""
    values = Z if isinstance(Z, (tuple, list)) else (Z, Z)
    if len(values) != 2:
        raise ValueError('One exact center or two exact coordinate endpoints required')
    lo, hi = (correlated.inverse.exact_scalar(v) for v in values)
    if not -1 <= lo <= hi <= 1:
        raise ValueError('Real Z coordinate box in [-1,1] required')
    a = ctx.mpf(lo.numerator)/lo.denominator
    b = ctx.mpf(hi.numerator)/hi.denominator
    return ctx.mpf([ends(a)[0], ends(b)[1]])


class CurrentOriginalRpPressureBinding:
    @source_precision
    def __init__(self, before, require_checked=True):
        if type(before) is not signed.CurrentOriginalRpSignedPhysicalValues or not before.acceptance_loaded:
            raise ValueError('Accepted live current signed physical source required')
        self.before = before
        self.amplitude = before.before
        self.inlet, self.selected = self.amplitude.inlet, self.amplitude.selected
        self.datums = dict(inlet=self.inlet.datum, selected=self.selected.datum)
        self.graph, self.ctx, self.family_record = before.graph, before.ctx, before.family_record
        self.frame = self.amplitude.exact.frame
        self.binding = self.frame.reports[amplitude.exact.frame.native.NAME]['current_P0_source_binding']
        self._binding_snapshot = json.dumps(self.binding, sort_keys=True)
        self.hashes = dict(before.hashes)
        for name in (signed.NAME, signed.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        self.acceptance_loaded = False
        self._datum_snapshots = {name: datum_snapshot(datum) for name, datum in self.datums.items()}
        self._datum_objects = dict(self.datums)
        self.P0_nodes = tuple(row.node for row in amplitude.exact.frame.target.rows(
            self.frame.functions['actual_identified_Rp_native_frame_C3']['P0']))
        self._P0_node_snapshots = tuple(json.dumps(self.graph.nodes[node], sort_keys=True) for node in self.P0_nodes)
        z = s.Symbol('Z', real=True)
        m0, m2 = s.symbols('pressure_M0 pressure_M2', real=True)
        flat = s.Function('F_flat')(z)
        analytic = -m0-m2/(1+z*z)**2-flat
        source = s.sympify(self.binding['exact_analytic_function'], locals={
            'Z': z, 'pressure_M0': m0, 'pressure_M2': m2, 'F_flat': s.Function('F_flat')})
        if s.simplify(source-analytic) != 0:
            raise ValueError('Actual current analytic P0 definition required')
        self.analytic_rows = [s.diff(analytic, z, n)/math.factorial(n) for n in range(6)]
        recorded = [s.sympify(expr, locals={'Z': z, 'pressure_M0': m0,
            'pressure_M2': m2, 'F_flat': s.Function('F_flat')})
            for expr in self.binding['first_six_true_Taylor_coefficients']]
        if len(recorded) != 6 or any(s.simplify(a-b) != 0 for a, b in zip(recorded, self.analytic_rows)):
            raise ValueError('Current exact six Taylor rows, with factorial once, required')
        self.q_rows = [s.diff((1+z*z)**-2, z, n)/math.factorial(n) for n in range(6)]
        previous = s.Integer(0)
        for n in range(5):
            following = -(2*z*(n+2)*self.q_rows[n]+(n+3)*previous)/((1+z*z)*(n+1))
            if s.cancel(following-self.q_rows[n+1]) != 0:
                raise ValueError('Numeric q_jets recurrence must give exact ordinary derivatives')
            previous = self.q_rows[n]
        self.scale_log = self.graph.mul(self.graph.constant(2), self.amplitude.logP)
        self.proof = self.graph.node('current_original_Rp_analytic_P0_numeric_source_binding',
            original_native_P0_function_nodes=list(self.P0_nodes),
            source_family=self.family_record, analytic_definition=self.binding['analytic_definition'],
            actual_exact_analytic_function=s.sstr(analytic),
            ordinary_Z_Taylor_rows=[s.sstr(v) for v in self.analytic_rows],
            physical_pressure_scale_log=self.scale_log.node,
            amplitude_parameter_binding=self.amplitude.proof.node,
            two_distinct_live_datum_owners_retained=self.datums['inlet'] is not self.datums['selected'],
            analytic_flatten_function_not_replaced_by_Cauchy_bound=True,
            pressure_Mp_is_separate_from_P0=True)
        self._proof_snapshot = json.dumps(self.graph.nodes[self.proof.node], sort_keys=True)
        self.assert_graph()
        if require_checked:
            receipt = json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family'] != self.family_record:
                raise ValueError('Current P0 inclusion receipt or scope differs')
            correlated.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes, receipt['input_hashes'])
            self.acceptance_loaded = True

    def assert_graph(self):
        self.before.assert_graph()
        checks = dict(same_current_native_frame=self.frame is self.amplitude.exact.frame,
            same_original_P0_definition=json.dumps(self.binding, sort_keys=True) == self._binding_snapshot,
            same_current_inlet_datum=self.inlet.datum is self._datum_objects['inlet'],
            same_current_selected_datum=self.selected.datum is self._datum_objects['selected'],
            actual_selected_future_flatten_P0_owner=self.selected.datum is self.selected.flatten.inlet.datum
                is self.selected.pulse.selection.future.angular.initial.datum,
            actual_typed_numeric_datum_owners=all(type(v) is CompliantPressureDatum and
                v.normalized_jets.__func__ is LogarithmicPressureDatum.normalized_jets for v in self.datums.values()),
            unchanged_live_datum_bounds=all(datum_snapshot(v) == self._datum_snapshots[k] for k, v in self.datums.items()),
            unchanged_original_P0_function_nodes=tuple(json.dumps(self.graph.nodes[n], sort_keys=True)
                for n in self.P0_nodes) == self._P0_node_snapshots,
            same_original_pressure_scale=self.graph.nodes[self.proof.node]['physical_pressure_scale_log'] == self.scale_log.node
                and self.scale_log == self.graph.mul(self.graph.constant(2), self.amplitude.logP),
            unchanged_numeric_binding_proof=json.dumps(self.graph.nodes[self.proof.node], sort_keys=True) == self._proof_snapshot)
        for name, datum in self.datums.items():
            checks[name+'_same_current_source'] = (datum.definition == self.binding['analytic_definition'] and
                datum.source_sha == self.family_record['implicit_source_sha256'] and
                datum.datum_sha == self.family_record['datum_enclosure_sha256'])
            checks[name+'_fourteen_original_atoms'] = set(datum.stages) == set(BETA2+BETA0+('z_flatten',))
            checks[name+'_positive_late_tail'] = ends(datum.tail_upper)[0] > 0
            checks[name+'_original_logP_numeric_inclusion'] = amplitude.contains(
                datum.parameters.logPstar, self.amplitude.logP_box)
        if not all(checks.values()):
            raise ValueError('Current P0 source changed: '+str(checks))
        return checks

    @source_precision
    def _rows(self, Zbox):
        self.assert_graph()
        # The actual inlet callback converts the same raw datum rows to its
        # 240-digit Taylor context. It does not select an enclosure endpoint.
        Z = self.inlet.ctx.mpf(ends(Zbox))
        incoming = self.inlet.incoming(Z)
        owners = {}
        for name, datum in self.datums.items():
            packet = datum.normalized_jets(Z, 5)
            rows = packet['normalized_pressure_coefficients']
            if packet['pressure_units'] != 'P/Pstar^2' or not packet['all_14_pressure_atoms_included'] or \
                    packet['implicit_source_sha256'] != self.family_record['implicit_source_sha256'] or \
                    packet['datum_enclosure_sha256'] != self.family_record['datum_enclosure_sha256']:
                raise ValueError('Actual current complete normalized P0 rows required')
            if any(incoming['P0'][n]._mpi_ != self.inlet.ctx.mpf(ends(row))._mpi_
                   for n, row in enumerate(rows)):
                raise ValueError('Actual inlet P0 rows differ from the live defining datum')
            c = datum.ctx
            q = q_jets(c, c.mpf(Z), 5)
            flatten = []
            for n in range(6):
                if n == 0:
                    term = datum.stages['z_flatten']['mass']
                elif n % 2 and ends(c.mpf(Z)) == (0, 0):
                    term = c.mpf(0)
                else:
                    bound = datum.flatten_complex_upper/datum.rho**n
                    term = c.mpf([ends(-bound)[0], ends(bound)[1]])
                flatten.append(term)
                direct = -(datum.m2*q[n]+(datum.m0 if n == 0 else 0)+term)
                if direct._mpi_ != rows[n]._mpi_:
                    raise ValueError('Fourteen-atom analytic P0 numeric inclusion differs')
            owners[name] = dict(normalized_P0_Taylor_rows=rows, normalized_q_Taylor_rows=q,
                flatten_Taylor_enclosures=flatten, fixed_beta0_mass=datum.m0, fixed_beta2_mass=datum.m2,
                positive_late_tail=datum.tail_upper, flatten_complex_upper=datum.flatten_complex_upper,
                source_sha256=datum.source_sha, datum_sha256=datum.datum_sha,
                actual_numeric_datum_method_called=True, actual_inlet_P0_rows_identical=True,
                all_fourteen_atoms_and_nonzero_late_tail_retained=True)
        return dict(source_family=self.family_record, center_Z=self.inlet.ctx.mpf(Z),
            ordinary_Taylor_order=5, coefficient_convention='dZ^j/j!; factorial applied exactly once',
            exact_analytic_P0_Taylor_rows=[s.sstr(v) for v in self.analytic_rows],
            actual_live_owners=owners, actual_inlet_P0=incoming['P0'], separate_inlet_Mp=incoming['Mp'],
            pressure_is_P0_plus_separate_Mp=True, actual_inlet_pressure=incoming['pressure'],
            exact_pressure_scale_log_function=self.scale_log.node,
            directed_pressure_scale_log=self.amplitude.reader().at(self.scale_log),
            original_P0_numeric_inclusion_proof=self.proof.node,
            mass_enclosures_not_selected_function_values=True, physical_Pstar_squared_not_expanded=True,
            full_physical_accuracy=False, **dict.fromkeys(GATES, self.acceptance_loaded), **dict.fromkeys(OPEN, False))

    @source_precision
    def normalized_jets(self, Z):
        """P0 source bounds for an exact real axial coordinate or coordinate box."""
        return self._rows(exact_Z_box(self.inlet.ctx, Z))

    @source_precision
    def source_binding(self, delivery):
        reader, _, request = self.amplitude.before._validate_delivery(delivery)
        result = self._rows(reader.at(delivery['physical_inverse']['coordinate_functions']['Z']))
        result.update(chart=request.chart, exact_Z=str(request.Z),
            original_physical_inverse_proof=delivery['physical_inverse']['exact_original_inverse_identity'].node,
            same_actual_supported_physical_source=True)
        return result

    @source_precision
    def velocity_pressure(self, delivery, relative_width_target='1/100000000'):
        source = self.source_binding(delivery)
        packet = self.before.velocity_pressure(delivery, relative_width_target)
        packet['current_P0_numeric_inclusion_pending'] = not self.acceptance_loaded
        packet['current_P0_numeric_source_binding'] = source
        packet['values']['p']['current_P0_numeric_inclusion_pending'] = not self.acceptance_loaded
        packet['values']['p']['ordinary_pressure_materialization_not_a_current_P0_binding'] = False
        packet.update(**dict.fromkeys(GATES, self.acceptance_loaded))
        # This gate certifies the pressure source, not absolute pressure error,
        # a general xyz/t API or a completed heat/stress/background field.
        return packet


@source_precision
def run(before, deliveries):
    began = time.monotonic()
    owner = CurrentOriginalRpPressureBinding(before, require_checked=False)
    if not deliveries:
        raise ValueError('Actual live correlated physical source deliveries required')
    views = {name: owner.source_binding(delivery) for name, delivery in deliveries.items()}
    result = dict(candidate_current_P0_numeric_source_binding_constructed=True,
        source_family=owner.family_record, actual_P0_source_bindings=views,
        exact_numeric_binding_proof=owner.graph.nodes[owner.proof.node],
        source_graph_assertions=owner.assert_graph(), **dict.fromkeys(GATES+OPEN, False),
        input_hashes=owner.hashes, execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result), separators=(',', ':'))+'\n').encode(), mtime=0))
    print('Current P0: both actual fourteen-atom datum owners bound to six original Taylor rows', flush=True)
    return owner, views
