"""Saved forward m15 plus one independently sourced reserved response.

No ancestor producer is called. The 212 physical sectors are retained
verbatim, including their numerical and finite-N errors. The complete
pulse is an independent target, never an extra forward-prefix term.
"""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_preRp_signed_inner_m15 import INNER_TERMS
from .original_reserved_m15_response import OriginalReservedM15Response


FORWARD = PREFIX+'compliant_current_original_compact_quiet_m15_numeric.json'
FORWARD_SHA256 = '7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071'


class OriginalRepairedRpM15:
    """Continue only saved axial queries; new queries require a new forward run."""

    def __init__(self, *, response=None, **kwargs):
        if response is not None and kwargs:
            raise ValueError('Supply response owner or constructor options')
        self.response = r = response if response is not None else OriginalReservedM15Response(**kwargs)
        self.ctx, self.iv, self.directory, self.family = r.ctx, r.iv, r.directory, r.family
        self.hashes, self.code_hashes = dict(r.hashes), dict(r.code_hashes)
        raw = (self.directory/FORWARD).read_bytes()
        if hashlib.sha256(raw).hexdigest() != FORWARD_SHA256:
            raise ValueError('Accepted complete 212-term checkpoint required')
        data = json.loads(raw)
        if data['source_family'] != self.family or not data['checks']['passed']:
            raise ValueError('Checked same-source forward checkpoint required')
        self.merge(data['input_hashes'], self.directory, self.hashes)
        self.merge(data['implementation_hashes'], ROOT, self.code_hashes)
        self.hashes[FORWARD] = FORWARD_SHA256
        self.saved = {Fraction(p['Z']): p for p in data['samples']
                      if p['chart'] == 'quiet' and Fraction(p['native_coordinate']) == 1}
        self.pulse = r.load('n1_pulse_numeric_bridge')
        self.pre = r.load('n1_preRp_numeric_bridge')
        self.pre_functions = r.load('n1_preRp_moment_functions')
        self.hashes.update(r.hashes)
        self.full = r.full['samples']['whole_Z']
        identity = r.full['source_identity']
        if not all(identity[k] for k in ('inner_lower_preRp_pulse_pressure_partitions_complete',
                'preRp_pulse_same_original_positive_reference_A_Rv_units',
                'pressure_unit_Z1_exact_zero','m15_minus_already_in_original_signed_components')):
            raise ValueError('Original complete signed m15 partition and fixed units required')
        pulse_identity = self.pulse['source_identity']
        if not all(pulse_identity[k] for k in ('passed',
                'exact_main_gap_end_pressure_C1_graph_and_endpoints_replayed',
                'all_three_charts_and_native_mean_histories_retained',
                'Pstar_applied_once_and_signed_m15_once')):
            raise ValueError('Independent actual complete-pulse source required')
        for name in ('original_repaired_Rp_m15.py','original_preRp_signed_inner_m15.py',
                     'original_n1_inner_m15_numeric.py','original_first_micro_m15.py',
                     'original_restore_m15.py','original_reference_pressure_rows.py'):
            name = 'src/openai_ns_reconstruction/'+name
            digest = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            if self.code_hashes.get(name, digest) != digest:
                raise ValueError('Changed source for endpoint identity: '+name)
            self.code_hashes[name] = digest
        self._bind_partition()
        self._bind_inner()

    def saved_report(self, suffix):
        name = PREFIX+'compliant_current_original_'+suffix+'.json'
        raw = (self.directory/name).read_bytes()
        if self.hashes.get(name) != hashlib.sha256(raw).hexdigest():
            raise ValueError('Required saved source is not a pinned ancestor: '+name)
        data = json.loads(raw)
        if data['source_family'] != self.family or not data['checks']['passed']:
            raise ValueError('Checked saved same-source report required: '+name)
        return data

    def merge(self, group, base, target):
        self.response.merge(group, base, target)

    def _bind_partition(self):
        """Check actual graph handles and normalization, not numerical overlaps."""
        r = self.response
        terminal = r.terminal['terminal_identity_bindings']
        nodes = r.terminal['exact_graph_nodes']
        pulse = self.pulse['samples']['whole_Z']
        if (self.pre_functions['signed_pieces'] != self.pre['true_function_bindings']['pieces']
                or self.pre_functions['signed_pressure_C1_functions'] != self.pre['true_function_bindings']['functions']):
            raise ValueError('Independent legacy preRp piece/integral handles changed')
        if not all(self.pre['source_identity'][key] for key in
                   ('exact_replayed_active_signed_C1_integrals_and_finite_endpoints',
                    'exact_replayed_compact_and_quiet_physical_C1_function_bodies',
                    'all_17_charts_68_views_228_native_cells_retained',
                    'all_228_C3_primitive_handles_identical_to_actual_live_source_functions')):
            raise ValueError('Exact original active/compact/quiet source replay required')
        self.pulse_physical = pulse['true_function_bindings']['physical']['m15']
        full = r.full['true_function_bindings']['physical_five']['m15']
        components = {name: [] for name in ('inner','lower','preRp')}
        inner_lower_pair = []
        for k in (0, 1):
            for handle in (full[k], terminal['terminal_sources']['m15'][k]):
                node = nodes[handle]
                if (node['operation'] != 'sum' or len(node['arguments']) != 2
                        or node['arguments'].count(self.pulse_physical[k]) != 1):
                    raise ValueError('Terminal and complete defect must contain actual pulse once')
            nonpulse = [i for i in nodes[full[k]]['arguments'] if i != self.pulse_physical[k]]
            pre = self.pre['true_function_bindings']['functions']['m15'][k]
            if len(nonpulse) != 1 or nodes[nonpulse[0]]['operation'] != 'sum':
                raise ValueError('Actual full m15 nonpulse sum required')
            children = nodes[nonpulse[0]]['arguments']
            if len(children) != 2 or children[-1] != pre:
                raise ValueError('Exact independent preRp graph component required')
            inner_lower = nodes[children[0]]
            if inner_lower['operation'] != 'sum' or len(inner_lower['arguments']) != 2:
                raise ValueError('Exact inner plus lower m15 graph components required')
            components['inner'].append(inner_lower['arguments'][0])
            components['lower'].append(inner_lower['arguments'][1])
            components['preRp'].append(pre)
            inner_lower_pair.append(children[0])
        a, b = self.full['pressure_normalization'], pulse['normalization']['numerical_source_units']
        for key in ('source_logRc_formula','source_logRp_formula','source_Rv_formula',
                    'source_A','source_logE_formula','normalization','m15'):
            if a[key] != b[key]:
                raise ValueError('Independent complete-pulse m15 units differ: '+key)
        if a['m15'] != dict(A_power=1, Rv_power=2, inverse_mu_numerator=37, inverse_mu_denominator=2):
            raise ValueError('Original p=2 pressure unit required, not m13/P1')
        self.unit = copy.deepcopy(self.full['actual_physical_pressure_unit_logs']['m15'])
        self.bindings = {key: copy.deepcopy(terminal[key]['m15']) for key in
                         ('terminal_sources','reduced_patch','reduced_terminal','patch_endpoint')}
        self.bindings.update(full_unrepaired_m15=copy.deepcopy(full),
            archived_nonpulse_components=components,
            archived_inner_lower_sum=inner_lower_pair,
            complete_pulse_physical=copy.deepcopy(self.pulse_physical),
            complete_pulse_normalized=copy.deepcopy(pulse['true_function_bindings']['normalized']['m15']),
            legacy_preRp_pieces=copy.deepcopy({k: v['m15'] for k, v in self.pre['true_function_bindings']['pieces'].items()}))

    def _bind_inner(self):
        data = self.saved_report('n1_inner_m15_numeric')
        self.inner = {Fraction(p['Z']): p for p in data['samples']
                      if Fraction(p['q']) == 2 and p['cells_per_piece'] == 16}
        self.inner_at_Rin = {Fraction(p['Z']): p for p in data['samples']
                             if Fraction(p['q']) == 1 and p['cells_per_piece'] == 16}
        if not data['checks']['passed']:
            raise ValueError('Checked original six signed inner/lower contributions required')
        first = self.saved_report('first_micro_full_m15_numeric')
        self.first = {Fraction(p['Z']): p for p in first['samples']
                      if p['coordinate'] == dict(kind='exact_decimal_phase',value='1')}
        restore = self.saved_report('restore_m15_numeric')
        self.restore = {Fraction(p['Z']): p for p in restore['samples']
                        if p['chart'] == 'restoration' and p['coordinate']['value'] == '1'}
        self.post_start = {Fraction(p['Z']): p for p in restore['samples']
                           if p['chart'] == 'postrestore' and p['coordinate']['value'] == '0'}

    def crosswalk(self, packet):
        terms = packet['m15_C0_C1_terms']
        if len(terms) != 212 or tuple(t['name'] for t in terms[:6]) != INNER_TERMS:
            raise ValueError('Complete 212-sector forward state required')
        if terms[:6] != self.inner[Fraction(packet['Z'])]['m15_C0_C1_terms']:
            raise ValueError('Original six signed inner/lower sectors changed or duplicated')
        expected = ('compact_exact_original_mean_square_integral','compact_finite_N_endpoint_error',
                    'compact_finite_N_volume_error','compact_finite_N_viscous_error',
                    'quiet_exact_original_mean_square_integral')
        if tuple(t['name'] for t in terms[207:]) != expected:
            raise ValueError('Actual compact four-sector and quiet one-sector continuation required')
        if len({t['name'] for t in terms}) != 212 or any('own_swirl' in t['name'] for t in terms):
            raise ValueError('Forward state must not already contain reserved own_swirl')
        # The old active piece contains every original chart rminus..Rc.
        # The legacy inner/lower cap boxes are never appended to these sectors.
        # The four q=2 collar sectors span BOTH the cutoff to Rin and the
        # lower collar. They cannot be assigned to legacy lower alone.
        groups = (('inner+lower',0,6,'axis..rminus; core plus cutoff/lower q=0..2'),
                  ('active',6,207,'rminus..Rc; all original micro/macro/switch/outer charts'),
                  ('compact',207,211,'Rc..2Rc; dy=dx/x'),
                  ('quiet',211,212,'2Rc..Rp; dy=L dq'))
        return [dict(legacy_component=name, forward_indices=list(range(lo,hi)),
                     forward_names=[t['name'] for t in terms[lo:hi]], physical_domain=scope,
                     legacy_component_appended=False,
                     legacy_graph_pair=copy.deepcopy(self.bindings['archived_inner_lower_sum'] if name=='inner+lower' else
                         self.bindings['legacy_preRp_pieces'].get(name)),
                     exact_original_integrals_and_endpoints=copy.deepcopy(
                         self.pre_functions['defining_bindings'].get(name,{})),
                     radial_weight='R*dR times original signed m15 density; outer pressure part -R*N1p*dR',
                     numeric_approximation_and_all_error_sectors_retained=True)
                for name,lo,hi,scope in groups]

    def inner_lower_split(self, Z):
        """Map the two legacy domains without misassigning collar sectors.

        q=1 ends at Rin, q=2 at rminus. The common core is counted once;
        the lower collar is the signed difference of the same prefix
        functions at2 and1. The saved q=1 enclosures are evidence and do
        not get appended to the forward state.
        """
        key=Fraction(Z);end,rin=self.inner[key],self.inner_at_Rin[key]
        if (rin['chart'] != 'selected_first_cutoff' or end['chart'] != 'selected_lower_collar'
                or rin['exact_sc_tuple'] != end['exact_sc_tuple']
                or rin['exact_phase_tuple'] != end['exact_phase_tuple']
                or rin['core_inlet'] != end['core_inlet']
                or rin['m15_C0_C1_terms'][:2] != end['m15_C0_C1_terms'][:2]):
            raise ValueError('Same core and exact q=1/2 cutoff/lower source split required')
        return dict(common_coordinate='R(q)=4epsilon*exp(hb*(sc/4)*q)',
            core_domain='axis..4epsilon',cutoff_domain='q=0..1; 4epsilon..Rin',
            lower_domain='q=1..2; Rin..rminus',
            legacy_inner_identity='core + collar_prefix(q=1)',
            legacy_lower_identity='collar_prefix(q=2) - collar_prefix(q=1)',
            joint_six_sector_identity='inner+lower = core + collar_prefix(q=2)',
            legacy_inner_graph_pair=copy.deepcopy(self.bindings['archived_nonpulse_components']['inner']),
            legacy_lower_graph_pair=copy.deepcopy(self.bindings['archived_nonpulse_components']['lower']),
            saved_q1_six_sector_enclosures=copy.deepcopy(rin['m15_C0_C1_terms']),
            shared_q1_endpoint_cancels_as_a_source_function=True,
            four_collar_sectors_do_not_belong_to_lower_alone=True,
            saved_q1_terms_appended_to_forward=False,
            enclosure_subtraction_not_selected_as_a_signed_source_value=True)

    def endpoint_reduction(self, terms, Z):
        """Only A0*sigma(sc/2) telescopes; retain both width remainders.

        Both producers use the identical core rho=4 seed and the formula
        A0=32 Z Psi_rho(4)/(1-delta Z^2). Their common sc/2 is fixed in Z.
        Equality of boxes is an additional replay guard, not the proof.
        """
        key = Fraction(Z)
        inner, first = self.inner[key], self.first[key]
        phase = list(inner['exact_sc_tuple']); phase[2] -= 2
        if (inner['exact_phase_tuple'] != phase or first['original_sc_tuple'] != inner['exact_sc_tuple']
                or first['finite_N_full_prefix_cover']['original_sc_tuple'] != inner['exact_sc_tuple']
                or inner['core_inlet']['rho'] != '4'
                or first['incoming_q2_m15_terms'] != inner['m15_C0_C1_terms']):
            raise ValueError('Identical selected sc/phase and shared rho=4 core source required')
        a, b = terms[4], terms[8]
        if (a['name'] != 'first_lower_singular_pressure' or b['name'] != 'first_micro_singular_lower'
                or a['actual_function_binding']['A0'] != '32Z*Psi_core_rho(4)/L'
                or b['actual_function_binding']['primitive'] != 'A0*(sigma(S)-sigma(sc/2))'
                or b['actual_function_binding']['endpoint_sign'] != -1):
            raise ValueError('Source-owned singular endpoint identity required')
        def opposite(raw):
            return [1-raw[0], *raw[1:]] if raw[1] else raw
        if any(x['source_enclosure_exact_mpf_tuples'] != [opposite(v) for v in
                y['source_enclosure_exact_mpf_tuples'][::-1]]
               for x,y in zip(a['coefficient_C0_C1'], b['coefficient_C0_C1'])):
            raise ValueError('Shared singular core coefficient replay changed')
        if a['log_scale_terms'] != b['log_scale_terms']:
            # Serialization metadata can differ; exact directed atoms cannot.
            if [x['source_enclosure_exact_mpf_tuples'] for x in a['log_scale_terms']] != [
                    x['source_enclosure_exact_mpf_tuples'] for x in b['log_scale_terms']]:
                raise ValueError('Same physical epsilon squared and sigma(sc/2) atoms required')
        cancellations = [dict(indices=[4,8], identity='epsilon^2*A0(Z)*sigma(sc/2) - epsilon^2*A0(Z)*sigma(sc/2)=0',
            ordinary_Z1_identity=True, endpoint_Z_independent=True,
            common_core_seed='original core rho=4, OriginalN1CutoffNumeric._inputs',
            exact_sc_tuple=copy.deepcopy(inner['exact_sc_tuple']),
            exact_first_phase_tuple=phase, first_micro_lower_equals_sc_half=True,
            source_functions=['OriginalN1InnerM15Numeric._prepare','OriginalFirstMicroM15.prepare'],
            no_error_caps_cancelled=True, raw_terms_retained=True,
            retained_width_error_names=[terms[i]['name'] for i in (5,9)])]
        restore, post = self.restore[key], self.post_start[key]
        if (post['incoming_m15_terms'] != restore['m15_C0_C1_terms']
                or terms[:140] != restore['m15_C0_C1_terms']):
            raise ValueError('Same full restoration endpoint consumed by postrestore required')
        # The actual row owner chains postrestore.previous_rows=restoration.
        # At this fixed interface all leading U/M/Q/Qy fields are the same:
        # restoration offset -8 plus1 equals postrestore offset -7 plus0.
        if [row['source_enclosure_exact_mpf_tuples'] for row in
                restore['geometry']['atomic_log_offset_terms']] != [row['source_enclosure_exact_mpf_tuples'] for row in
                post['geometry']['atomic_log_offset_terms']]:
            raise ValueError('Actual restoration/postrestore interface radius must be identical')
        pairs = [(124+d,143+d,2) for d in range(3)]+[(133+d,151+d,1) for d in range(2)]
        for ia,ib,power in pairs:
            a,b = terms[ia],terms[ib]
            expected = ('R2F_upper_epsilon_'+str(ia-124),'R2F_lower_epsilon_'+str(ia-124)) if power==2 else (
                'RQy_upper_epsilon_'+str(ia-133),'RQy_lower_epsilon_'+str(ia-133))
            if (a['name'] != 'restoration_'+expected[0] or b['name'] != 'postrestore_'+expected[1]
                    or a['actual_function_binding']['signed_endpoint_orientation'] != -b['actual_function_binding']['signed_endpoint_orientation']
                    or not a['actual_function_binding'].get('actual_changing_source_endpoint',power==2)
                    or not b['actual_function_binding'].get('actual_changing_source_endpoint',power==2)):
                raise ValueError('Original signed leading restoration endpoint pair required')
            if any(x['source_enclosure_exact_mpf_tuples'] != [opposite(v) for v in
                    y['source_enclosure_exact_mpf_tuples'][::-1]] for x,y in zip(a['coefficient_C0_C1'],b['coefficient_C0_C1'])):
                raise ValueError('Leading shared endpoint coefficient replay changed')
            aa,bb = a['log_scale_terms'],b['log_scale_terms']
            if [row['source_enclosure_exact_mpf_tuples'] for row in aa[:3]] != [
                    row['source_enclosure_exact_mpf_tuples'] for row in bb[:3]]:
                raise ValueError('Shared endpoint source scales differ')
            for rows in (aa,bb):
                lo,hi = _bounds(self.ctx, sum((self.response.decode(v) for v in rows[3:]),self.iv.mpf(0)))
                if lo != hi or lo != -7*power:
                    raise ValueError('Exact shared interface offset required')
            cancellations.append(dict(indices=[ia,ib],identity='same leading physical endpoint with opposite orientations',
                source_functions=['OriginalRestoreM15.next_field','OriginalReferencePressureRows.inputs'],
                common_interface='10 logCstar +10 logPstar -8+1 = 10 logCstar +10 logPstar -7+0',
                endpoint_power=power,ordinary_Z1_identity=True,endpoint_Z_independent=True,
                exact_incoming_source_frame_replayed=True,no_error_caps_cancelled=True,raw_terms_retained=True))
        removed = {i for row in cancellations for i in row['indices']}
        return copy.deepcopy([t for i,t in enumerate(terms) if i not in removed]),cancellations

    def target(self, Z):
        """Read the negative complete pulse independently of the forward data."""
        r = self.response
        selected = self.pulse['samples']['whole_Z']
        for sample in self.pulse['samples'].values():
            lo, hi = _bounds(self.ctx, r.decode(sample['Z']))
            zl, zh = _bounds(self.ctx, self.iv.mpf(Z))
            if lo == hi and lo <= zl <= zh <= hi:
                selected = sample
                break
        rows = selected['normalized_true_pulse_C1_enclosures']['m15']
        return dict(name='independent_negative_complete_pulse_m15',
            coefficient_C0_C1=[r.record(-r.decode(row)) for row in rows],
            log_scale_terms=[copy.deepcopy(self.unit)],
            actual_function_binding=dict(physical_handles=self.pulse_physical,
                normalized_handles=selected['true_function_bindings']['normalized']['m15'],
                source=PREFIX+'compliant_current_original_n1_pulse_numeric_bridge.json.gz',
                sign=-1, physical_measure='-R*pressure_density*dR', ordinary_Z1=True),
            selected_source_Z=copy.deepcopy(selected['Z']), independent_of_forward_endpoint=True)

    def compare(self, terms, target):
        """A signed-log comparison ledger; no unmaterializable exp is formed.

        Large positive normalized log caps remain explicit unresolved
        sectors. This routine never substitutes a cap for a function value.
        """
        r, iv = self.response, self.iv
        unit = r.decode(self.unit)
        sectors = []
        for i,term in enumerate(terms):
            log = sum((r.decode(row) for row in term['log_scale_terms']), iv.mpf(0))-unit
            sectors.append(dict(forward_index=i, name=term['name'],
                coefficient_C0_C1=copy.deepcopy(term['coefficient_C0_C1']),
                normalized_log_scale_enclosure=r.record(log),
                error_enclosure_only=term.get('actual_function_binding',{}).get('error_enclosure_only',False)))
        worst = []
        for k in (0,1):
            caps = []
            for sector in sectors:
                coeff = r.decode(sector['coefficient_C0_C1'][k])
                cap = _bounds(self.ctx, abs(coeff))[1]
                if cap:
                    log = r.decode(sector['normalized_log_scale_enclosure'])+iv.log(iv.mpf(cap))
                    caps.append((_bounds(self.ctx,log)[1],sector['forward_index']))
            top,index = max(caps)
            worst.append(dict(ordinary_Z_order=k,largest_individual_log_cap=r.record(iv.mpf(top)),
                              forward_index=index, name=terms[index]['name']))
        # At full support the selected response equals -full_m15. Its
        # independently defined target is -pulse_m15. The checked original
        # graph is full=(inner+lower+preRp)+pulse, so the shared pulse cancels
        # as a source function before enclosure arithmetic. The new forward
        # sectors are still independently approximated: no BG cap is removed.
        components = self.full['pressure_component_enclosures']['m15']
        nonpulse = []
        for name in ('inner','lower','preRp'):
            rows = components[name]['C1_enclosures'] if name != 'preRp' else components[name]
            nonpulse.append(dict(name='minus_archived_'+name,
                physical_handles=copy.deepcopy(self.bindings['archived_nonpulse_components'][name]),
                coefficient_C0_C1=[r.record(-r.decode(row)) for row in rows]))
        return dict(definition='sum(forward212 + reserved_once) - independent_negative_complete_pulse',
            normalized_forward_sectors=sectors,
            target_subtracted_C0_C1=[r.record(-r.decode(row)) for row in target['coefficient_C0_C1']],
            largest_individual_absolute_log_caps=worst,
            source_correlated_difference=dict(
                identity='reserved_once - (-pulse) = -full + pulse = -(inner+lower+preRp)',
                exact_shared_pulse_physical_handles=copy.deepcopy(self.pulse_physical),
                normalized_subtracted_nonpulse_components=nonpulse,
                retained_forward_indices=list(range(212)),
                reserved_and_target_replaced_only_by_exact_source_difference=True,
                independent_forward_212_remain=True,
                no_forward_error_sector_removed=True,independent_numeric_residual_resolved=False),
            caps_are_bounds_not_actual_values=True, numerical_seam_residual_resolved=False,
            interval_overlap_used_as_closure=False, enormous_exp_materialized=False)

    def evaluate(self, Z):
        if type(Z) is not str or not self.ctx.isfinite(self.ctx.mpf(Z)):
            raise ValueError('Saved finite decimal axial query required')
        key = Fraction(Z)
        if key not in self.saved:
            raise ValueError('No saved complete forward state at this Z; ancestors are not replayed')
        packet = self.saved[key]
        crosswalk = self.crosswalk(packet)
        response = self.response.evaluate(Z,'2')
        if response['angular_inverse_response'] != [self.response.record(self.iv.mpf(v)) for v in (0,0,1)]:
            raise ValueError('Exact full-support selected matrix action required')
        own = dict(name='reserved_n1_own_swirl_once',
            coefficient_C0_C1=copy.deepcopy(response['coefficient_C0_C1']),
            log_scale_terms=[copy.deepcopy(self.unit)],
            actual_function_binding=dict(terminal_patch=self.bindings['patch_endpoint'],
                reduced_response=self.bindings['reduced_patch'], physical_measure='-R*theta0*Utheta1 ds',
                angular_matrix_inverse_before_boxing=True, ordinary_Z1_product_rule=True))
        terms = copy.deepcopy(packet['m15_C0_C1_terms'])+[own]
        reduced, cancellations = self.endpoint_reduction(terms,Z)
        target = self.target(Z)
        return dict(Z=Z,source_family=self.family,actual_selected_integer=copy.deepcopy(packet['actual_selected_integer']),
            m15_C0_C1_terms=terms, reduced_m15_C0_C1_terms=reduced,
            original_forward_terms=212, reserved_terms_added=1,
            crosswalk=crosswalk, exact_internal_endpoint_cancellations=cancellations,
            inner_lower_domain_split=self.inner_lower_split(Z),
            reserved_response=response, independent_target=target,
            source_partition_bindings=copy.deepcopy(self.bindings),comparison=self.compare(terms,target),
            contract=dict(saved_forward_212_retained_verbatim=True,six_signed_inner_lower_terms_retained_once=True,
                no_legacy_component_reappended=True,reserved_added_once=True,
                independent_complete_pulse_is_target_only=True,P1_not_used_as_m15=True,
                ordinary_Z1_retained=True,fixed_reference_A_Rv_squared_unit_once=True,
                only_source_proved_internal_endpoint_reduced=True,all_finite_N_error_sectors_retained=True,
                ancestors_not_reexecuted=True),
            limits=dict(matrix_row_reduced_reserved_enclosure=True,complete_213_sector_repaired_Rp_enclosure=True,
                joint_response_input_correlations_resolved=False,
                new_numeric_atoms_registered_as_archived_graph_nodes=False,
                independently_resolved_signed_finite_N_values=False,independent_numeric_Rp_seam=False,
                relative_cancellation_accuracy=False,global_five_moment_numeric_join=False,
                common_global_N=False,exact_heat_exterior=False,stress_cone=False,flat_remainder=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
