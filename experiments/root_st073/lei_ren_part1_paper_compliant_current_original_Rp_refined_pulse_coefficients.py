"""Fresh original pulse coefficients with the proved exact U0 enclosure.

The source functions, future repair, P0, coordinates and positive branch are
unchanged. All U-dependent ratios and selected coefficients are recomputed;
no accepted constant dictionary or selection cache is overwritten.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
from types import SimpleNamespace, MethodType
import time

import lei_ren_part1_paper_compliant_current_original_Rp_physical_accuracy as accuracy
from lei_ren_part1_paper_compliant_axial_high_jets import positive_quadratic_jets
from lei_ren_part1_paper_compliant_pulse_high_jets import _SelectedSource
from lei_ren_part1_paper_compliant_outer_pulse_map import pulse_energy
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import EXPECTED
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

signed, correlated = accuracy.signed, accuracy.correlated
physical, box = signed.physical, correlated.box
HERE, PREFIX, sha, ends = accuracy.HERE, accuracy.PREFIX, accuracy.sha, accuracy.ends
NAME = PREFIX+'current_original_Rp_refined_pulse_coefficients.json.gz'
RECEIPT = PREFIX+'current_original_Rp_refined_pulse_coefficients_check.json'
GATES = ('current_original_Rp_exact_U0_in_actual_pulse_coefficients_installed',
         'current_original_Rp_coherent_refined_positive_selection_installed',
         'current_original_Rp_refined_pulse_physical_operator_rows_installed')
OPEN = accuracy.OPEN


def snapshot_constants(constants):
    return {key: value._mpi_ if hasattr(value, '_mpi_') else value
            for key, value in constants.items()}


def cache_snapshot(selected):
    return tuple((id(owner), name, id(value), tuple(value.keys()))
                 for owner in (selected.axial4, selected.fifth, selected.angular4,
                               selected.energy4, selected.pulse)
                 for name, value in vars(owner).items()
                 if isinstance(value, dict) and 'cache' in name)


class RefinedPositiveSelection:
    """Original C5 implicit equation evaluated on coherent incoming ratios."""
    def __init__(self, owner, future):
        self.owner, self.future, self.ctx, self.cache = owner, future, owner.ctx, {}

    @source_precision
    def select(self, Z):
        c = self.ctx
        Z = c.mpf(Z)
        if ends(Z)[0] < -1 or ends(Z)[1] > 1:
            raise ValueError('Original axial domain [-1,1] required')
        key = Z._mpi_
        if key in self.cache:
            return self.cache[key]
        k, b = self.owner.constants, self.owner.base
        g = IntervalTaylor(c, [Z+Z**3, 1+3*Z**2, 3*Z, 1, 0, 0])
        h = IntervalTaylor(c, [Z**2+2*Z**4+Z**6, 2*Z+8*Z**3+6*Z**5,
            1+12*Z**2+15*Z**4, 8*Z+20*Z**3, 2+15*Z**2, 6*Z])
        moments = [g*k['C1'], g*k['C2']]
        energy = h*k['C_E']+k['C0']
        rows = [jet*factor for jet, factor in zip(moments, b.incoming_factor_caps)]
        affine = b.linear_inverse(rows[0]*(-b.mu), -(rows[1]-rows[0]))
        future = self.future.future(Z)['Section7_34_weighted_future_Taylor']
        target = future-energy*b.mu+b.base
        A2, A1, A0 = b.K, target*0, -target
        for nu, uj, vj in zip(b.nu, affine, b.v):
            A2 += nu*vj*vj
            A1 += uj*(2*nu*vj)
            A0 += uj*uj*nu
        discriminant = A1[0]*A1[0]-4*A2*A0[0]
        if ends(A2)[0] <= 0 or ends(A0[0])[1] >= 0 or ends(discriminant)[0] <= 0:
            raise ArithmeticError('Original positive quadratic branch must be resolved')
        denominator = A1[0]+c.sqrt(discriminant)
        if ends(denominator)[0] <= 0:
            raise ArithmeticError('Stable positive-root denominator must stay positive')
        a0 = -2*A0[0]/denominator
        selected, derivative_denominator = positive_quadratic_jets(c, a0, A2, A1, A0)
        scaled = [uj+selected*vj for uj, vj in zip(affine, b.v)]
        if not ends(scaled[0][0])[1] < 0 < ends(scaled[1][0])[0]:
            raise ArithmeticError('Original selected end-coefficient signs lost')
        result = dict(Z=Z, selected_ap_Taylor=selected,
            selected_scaled_end_coefficient_Taylor=scaled,
            actual_affine_incoming_Taylor=affine,
            actual_row_normalized_incoming_Taylor=rows,
            incoming=dict(moment_Taylor=moments, energy_Taylor=energy,
                Z_independent_constant_definitions=k, ordinary_Taylor_order=5,
                exact_shapes='mi=Ci*(Z+Z^3); e=C0+C_E*(Z^2+2Z^4+Z^6)'),
            energy_target_Taylor=target,
            actual_weighted_incoming_energy_Taylor=energy*b.mu,
            actual_weighted_future_energy_Taylor=future,
            quadratic_coefficients=dict(A2=A2, A1=A1, A0=A0),
            positive_root_derivative_denominator=derivative_denominator,
            common_log_end_coefficient_scale=b.log_end_scale,
            fixed_incoming_row_factor_definitions=b.incoming_factor_log_definitions,
            ordinary_Taylor_order=5, actual_positive_C0_branch_preserved=True,
            selected_root_recomputed_from_original_quadratic=True,
            old_cached_selected_ap_not_used_as_refined_root=True,
            exact_selected_energy_equation='Kpulse*ap^2+sum_j nu_j*(uj+vj*ap)^2=(1-exp(-26))/4-mu*incoming_energy+weighted_full_future_energy',
            actual_selected_ap_c1_c2_C5_available=True,
            admitted_C4_prefix_preserved=False,
            C4_prefix_recomputed_from_same_original_equation=True,
            sixth_derivative_Taylor_remainder_available=False,
            temporal_recursion=False)
        self.cache[key] = result
        return result


class CurrentOriginalRpRefinedPulseCoefficients:
    @source_precision
    def __init__(self, before, require_checked=True):
        if type(before) is not accuracy.CurrentOriginalRpPhysicalAccuracy or not before.acceptance_loaded:
            raise ValueError('Accepted physical accuracy/source owner required')
        self.before, self.amplitude = before, before.before.amplitude
        self.correlated, self.selected = self.amplitude.before, self.amplitude.selected
        self.ctx, self.graph, self.family_record = self.amplitude.ctx, before.graph, before.family_record
        self.hashes = dict(before.hashes)
        for name in (accuracy.NAME, accuracy.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        self.acceptance_loaded, self._sources = False, {}
        self._constants_snapshot = snapshot_constants(self.selected.constants)
        self._cache_snapshot = cache_snapshot(self.selected)
        c = self.ctx
        self.constants = dict(self.selected.constants)
        self.constants['U'] = c.mpf(ends(self.amplitude.U0_box))
        invP = c.exp(-c.mpf(ends(self.amplitude.logP_box)))
        U = self.constants['U']
        self.constants.update(C1=invP*self.constants['M']/U,
            C2=invP*self.constants['K']/(U*U), C0=self.constants['E_Q']/(U*U),
            C_E=self.constants['E_Z']/(U*U))
        # This integral's enclosure, rather than U0, dominates the observed
        # selected-amplitude width. Recompute the identical positive source
        # integral on a finer partition; retain its original defining function.
        self.energy_cells, self.history_cells = 32768, 2048
        self.pulse_map = copy.copy(self.selected.amplitude.pulse)
        self.pulse_map.Kpulse = pulse_energy(c, cells=self.energy_cells)
        self.base = copy.copy(self.selected.amplitude)
        self.base.pulse, self.base.K, self.base.cache = self.pulse_map, self.pulse_map.Kpulse, {}
        # Separate all mutable future/angle caches, retaining their accepted
        # same-source entries. A new Z can use unchanged future algorithms.
        self.angular = copy.copy(self.selected.angular4)
        self.angular.cache = dict(self.selected.angular4.cache)
        self.energy = copy.copy(self.selected.energy4)
        self.energy.cache = dict(self.selected.energy4.cache)
        self.energy.angular = self.angular
        self.fourth = copy.copy(self.selected.axial4)
        self.fourth.constants, self.fourth.energy = self.constants, self.energy
        self.fourth.cache = dict(self.selected.axial4.cache)
        self.future = copy.copy(self.selected.fifth)
        self.future.fourth, self.future.angular4 = self.fourth, self.angular
        for key in ('cache', 'angular_cache', 'energy_cache'):
            setattr(self.future, key, dict(getattr(self.selected.fifth, key)))
        self.selection = RefinedPositiveSelection(self, self.future)
        self.pulse = copy.copy(self.selected.pulse)
        self.pulse.data_cache = {}
        self.pulse.fifth = self.selection
        self.pulse.high = SimpleNamespace(ctx=c, base=self.selected.amplitude, constants=self.constants,
            select=self.selection.select, energy=SimpleNamespace(future=self.future.future))
        self.pulse.high.base, self.pulse.pulse = self.base, self.pulse_map
        self.pulse.selection = _SelectedSource(self.pulse.high)
        self.pulse.Xp = self.pulse.inlet_H/U
        self.constants['Xp'] = self.pulse.Xp
        self._refined_snapshot = snapshot_constants(self.constants)
        self._K_snapshot, self._Xp_snapshot = self.base.K._mpi_, self.pulse.Xp._mpi_
        self.dispatch = SimpleNamespace(evaluate=self._dispatch)
        native = self.correlated.before
        self.raw = box._Proxy(native.pulse.actual)
        self.raw.radius, self.raw.selected = native.box_radius, self.dispatch
        self.raw.raw, self.raw.closed = native.raw, native.heat
        self.raw.assert_graph = self.assert_graph
        self.raw.acceptance_loaded = False
        self.raw.bind('scale', native._adapter_functions['pulse_scale'])
        self.raw.bind('pulse_view', native._adapter_functions['pulse_view'])
        for method in ('factor', 'expression_owner', 'evaluate'):
            self.raw.bind(method, getattr(box.pulse.CurrentOriginalRpRawPulseTransport, method))
        self.transport = box._Proxy(native.actual)
        self.transport.owner, self.transport.assert_graph = self.raw, self.assert_graph
        self.transport.acceptance_loaded = False
        self.transport.bind('velocity_rows', native._adapter_functions['heat_velocity_rows'])
        for method in ('factor', 'evaluate'):
            self.transport.bind(method, getattr(box.mixed.CurrentOriginalRpMixedTransport, method))
        self.proof = self.graph.node('current_original_Rp_coherent_U0_numeric_pulse_refinement',
            original_amplitude_identity=self.amplitude.proof.node,
            original_U0_function=self.amplitude.raw.U0.node,
            original_selected_source_proof=copy.deepcopy(self.selected.current_selection_source_proof),
            original_future_and_analytic_P0_unchanged=True,
            recomputed=['U','C1','C2','C0','C_E','Xp','incoming_C5','positive_ap_C5','end_coefficients','mixed_rows'],
            original_Kpulse_integral_repartitioned=True, energy_cells=self.energy_cells,
            original_positive_history_integrals_repartitioned=True, history_cells=self.history_cells,
            root_formula='-2*A0/(A1+sqrt(A1^2-4*A2*A0)); A2>0,A0<0',
            no_nominal_endpoint_or_cached_selected_root=True)
        self.assert_graph()
        if require_checked:
            receipt = json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family'] != self.family_record:
                raise ValueError('Refined current pulse receipt/scope differs')
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes, receipt['input_hashes'])
            self.acceptance_loaded = True
            self.raw.acceptance_loaded = self.transport.acceptance_loaded = True

    def assert_graph(self):
        self.before.assert_graph()
        checks = dict(accepted_exact_U0=self.amplitude.acceptance_loaded,
            original_constants_unmutated=snapshot_constants(self.selected.constants)==self._constants_snapshot,
            original_cache_keys_unmutated=cache_snapshot(self.selected)==self._cache_snapshot,
            separate_coherent_constants=self.constants is not self.selected.constants and
                self.pulse.high.constants is self.constants is self.fourth.constants,
            current_exact_U0_box=self.constants['U']._mpi_==self.ctx.mpf(ends(self.amplitude.U0_box))._mpi_,
            unchanged_refined_coherent_boxes=snapshot_constants(self.constants)==self._refined_snapshot and
                self.base.K._mpi_==self._K_snapshot and self.pulse.Xp._mpi_==self._Xp_snapshot,
            same_coherent_Xp=self.constants['Xp'] is self.pulse.Xp,
            unchanged_other_source_constants=all(self.constants[key]._mpi_==self.selected.constants[key]._mpi_
                for key in ('M','K','E_Q','E_Z')),
            original_H_and_P=self.pulse.inlet_H is self.selected.pulse.inlet_H and
                self.pulse.inlet_P is self.selected.pulse.inlet_P,
            same_current_P0=self.pulse.selection.future.angular.initial.datum is self.selected.datum,
            same_current_future_repair=self.future.fourth.energy.base is self.selected.future and
                self.angular.repair is self.selected.future.repair,
            separate_future_caches=all(getattr(self.future,key) is not getattr(self.selected.fifth,key)
                for key in ('cache','angular_cache','energy_cache')),
            selection_not_original_cached_root=self.pulse.fifth is self.selection,
            independent_same_source_K_integral=self.base.K is self.pulse_map.Kpulse and
                self.base.pulse is self.pulse_map and self.pulse.high.base is self.base and
                self.base.future is self.selected.future and self.pulse_map.basis is self.selected.amplitude.pulse.basis,
            fixed_numerical_partitions=self.energy_cells==32768 and self.history_cells==2048,
            exact_scales_and_context_unchanged=self.transport.ctx is self.ctx and self.transport.graph is self.graph)
        if not all(checks.values()):
            raise ValueError('Refined current source differs: '+str(checks))
        return checks

    @source_precision
    def _dispatch(self, chart, Z, coordinate):
        if chart not in box.pulse.PULSE:
            raise ValueError('This refinement owns only original active pulse charts')
        method, domain = EXPECTED[chart]
        self.correlated.before._admit(chart, self.ctx.mpf(Z), self.ctx.mpf(coordinate))
        options = {'cells':self.history_cells} if method=='main' else {}
        return dict(source_packet=getattr(self.pulse, method)(Z, coordinate, **options))

    @source_precision
    def source(self, delivery):
        self.assert_graph()
        reader, token, request = self.correlated._validate_delivery(delivery)
        if request.chart not in box.pulse.PULSE:
            raise ValueError('Refined pulse delivery requires an active original pulse chart')
        record = self._sources.get(id(delivery))
        if record is not None:
            if record[0] is not delivery or record[2] != self.correlated._fingerprint(box.report(record[1])):
                raise ValueError('Unchanged derived source view required')
            return record[1]
        Z = reader.at(delivery['physical_inverse']['coordinate_functions']['Z'])
        view = self.transport.evaluate(token.chart, Z, token)
        view.update(exact_axial_coordinate_function=delivery['physical_inverse']['coordinate_functions']['Z'].node,
            directed_axial_coordinate=Z, coherent_exact_U0_refinement_proof=self.proof.node,
            original_delivery_native_identity=delivery['exact_native_inverse_identity'].node,
            whole_native_and_axial_boxes_passed_to_original_algorithms=True,
            source_family=self.family_record, **dict.fromkeys(OPEN, False))
        self._sources[id(delivery)] = (delivery, view, self.correlated._fingerprint(box.report(view)))
        self.assert_graph()
        return view

    @source_precision
    def physical_rows(self, delivery):
        reader, token, request = self.correlated._validate_delivery(delivery)
        view = self.source(delivery)
        Zbox = reader.at(delivery['physical_inverse']['coordinate_functions']['Z'])
        proxy = box._Proxy(self.correlated.physical)
        owner = self
        def source_view(proxy, observed, Z):
            if observed is not owner.source(delivery) or correlated.inverse.exact_scalar(Z)!=request.Z:
                raise ValueError('Live unchanged refined source point required')
            for name in physical.SOURCE.values():
                base = view['original_factorized_values'][name]
                for k in range(5):
                    for j in range(5-k):
                        row = view['log_radius_mixed_rows'][name]['y%d_Z%d'%(k,j)]
                        powers = tuple(map(Fraction,(0,2,2))) if name=='pressure' and k else base.powers
                        expected = owner.transport.factor(name, token.chart, token, view['geometry'], powers, row.coefficients,k,j,False)
                        if type(row) is not physical.mixed.FactorizedMixedSourceRow or \
                                row.derivative!=(k,j) or row.coefficients.ctx is not owner.ctx or row.coefficients.order!=0 or \
                                (row.powers,row.source_units,row.log_scale_parts)!=(expected.powers,expected.source_units,expected.log_scale_parts):
                            raise ValueError('Original units, exact scales and ordinary refined mixed rows required')
            return request.native, request.Z, Zbox
        def coordinates(proxy, geometry, Z, log_tau, theta):
            if correlated.inverse.exact_scalar(Z)!=request.Z or log_tau!=request.forward_coordinates['log_tau'] or \
                    correlated.inverse.exact_scalar(theta)!=request.theta or \
                    owner.correlated._fingerprint(geometry)!=owner.correlated._fingerprint(delivery['actual_source_view']['geometry']):
                raise ValueError('Unchanged original physical coordinates required')
            return request.forward_coordinates
        proxy.source_view = MethodType(source_view, proxy)
        proxy.coordinates = MethodType(coordinates, proxy)
        proxy.bind('map_source_view', physical.CurrentOriginalRpPhysicalSourceMap.map_source_view)
        return proxy.map_source_view(view, request.Z, request.forward_coordinates['log_tau'], request.theta)

    @source_precision
    def evaluate(self, delivery, relative_width_target='1/100000000'):
        self.assert_graph()
        target = Fraction(relative_width_target)
        if not 0 < target < 1:
            raise ValueError('Exact relative width target in (0,1) required')
        mapped, reader = self.physical_rows(delivery), self.amplitude.reader(delivery)
        source_binding = self.before.before.source_binding(delivery)
        values = self.before.before.before
        sections, counts = {}, dict(rows=0, exact_zero_rows=0, nonzero_factored_target_rows=0,
            sign_unresolved_rows=0, ordinary_numeric_target_rows=0)
        identities = dict(original_physical_inverse_identity=delivery['physical_inverse']['exact_original_inverse_identity'].node,
            original_native_inverse_identity=delivery['exact_native_inverse_identity'].node)
        for section in ('Cartesian_spatial_rows','fixed_x_time_rows'):
            sections[section] = {}
            for component, rows in mapped[section].items():
                sections[section][component] = {}
                for label, row in (rows if section=='Cartesian_spatial_rows' else {'dt':rows}).items():
                    groups = values._canonical_groups(row, reader, identities)
                    result = signed.enclose_factored_sum(self.graph,self.ctx,reader,groups,target)
                    result['physical_accuracy'] = self.before.row_accuracy(result,reader,target)
                    result.update(component=component, physical_derivative=row.derivative,
                        coherent_exact_U0_refinement_proof=self.proof.node,
                        original_P0_source_binding_retained=component=='p',
                        signed_canonical_scale_groups=[dict(exact_scale_log_function=g['log_function'].node,
                            signed_coefficient_enclosure=g['coefficient'], exact_scale_identity_proofs=g['proofs']) for g in groups])
                    sections[section][component][label] = result
                    a = result['physical_accuracy']
                    counts['rows'] += 1
                    counts['exact_zero_rows'] += int(result['exact_zero_enclosure'])
                    counts['nonzero_factored_target_rows'] += int(not result['exact_zero_enclosure'] and a['factored_relative_width_satisfied'])
                    counts['sign_unresolved_rows'] += int(a.get('sign') is None)
                    counts['ordinary_numeric_target_rows'] += int(a['ordinary_numeric_delivery_target_satisfied'])
        self.assert_graph()
        return dict(source_family=self.family_record, physical_value_rows=sections, delivery_counts=counts,
            original_P0_source_binding=source_binding, coherent_exact_U0_refinement_proof=self.proof.node,
            full_certified_physical_accuracy=False, unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES,self.acceptance_loaded), **dict.fromkeys(OPEN,False))

    @source_precision
    def velocity_pressure(self, delivery, relative_width_target='1/100000000'):
        _,_,request = self.correlated._validate_delivery(delivery)
        evaluated = self.evaluate(delivery, relative_width_target)
        rows = evaluated['physical_value_rows']['Cartesian_spatial_rows']
        return dict(physical_coordinates={key:ref.node for key,ref in request.forward_coordinates.items()},
            values={name:rows[component]['x0_y0_z0'] for name,component in (('u','ux'),('v','uy'),('w','uz'),('p','p'))},
            source_family=self.family_record, original_P0_source_binding=evaluated['original_P0_source_binding'],
            full_certified_physical_accuracy=False, unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES,self.acceptance_loaded), **dict.fromkeys(OPEN,False))


@source_precision
def run(before, deliveries, relative_width_target='1/100'):
    began = time.monotonic()
    owner = CurrentOriginalRpRefinedPulseCoefficients(before, require_checked=False)
    views = {name:owner.evaluate(delivery,relative_width_target) for name,delivery in deliveries.items()}
    result = dict(source_family=owner.family_record, actual_refined_pulse_physical_values=views,
        source_graph_assertions=owner.assert_graph(), constants=owner.constants, Xp=owner.pulse.Xp,
        accepted_original_constants=owner.selected.constants, accepted_original_Xp=owner.selected.pulse.Xp,
        refined_Kpulse=owner.base.K, accepted_original_Kpulse=owner.selected.amplitude.K,
        numerical_partitions=dict(energy_cells=owner.energy_cells,history_cells=owner.history_cells),
        fresh_selections={name:owner.selection.select(owner.source(delivery)['Z']) for name,delivery in deliveries.items()},
        input_hashes=owner.hashes, execution_seconds=time.monotonic()-began,
        **dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner, views
