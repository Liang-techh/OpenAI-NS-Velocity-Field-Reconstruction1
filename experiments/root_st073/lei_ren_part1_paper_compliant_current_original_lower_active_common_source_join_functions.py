"""Actual lower/first-active traces, common histories and exact source binding.

The two trace graphs remain distinct. The checked same-source integral/IVP
proof identifies their functions; no zero graph, model or interval center
replaces a trace. Source bounds and other global joins remain separate.
"""
import copy
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_n1_intermediate_own_C3_chart_functions as previous
import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions as bridge_source
import lei_ren_part1_paper_compliant_current_core_first_interface as interface
import lei_ren_part1_paper_compliant_current_inner_exit_strict_collar as collar
import lei_ren_part1_paper_compliant_current_generic_left_inlet as left_source

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode=bridge_source.FunctionJet,bridge_source.rows,bridge_source.encode
Algebra=bridge_source.BridgeAlgebra
NAME=PREFIX+'current_original_lower_active_common_source_join_functions.json.gz'
RECEIPT=PREFIX+'current_original_lower_active_common_source_join_functions_check.json'
RECOVERY=PREFIX+'current_original_lower_active_archive_recovery.py'
GATES=('current_original_actual_lower_active_rminus_field_C5_Q4_pressure_C3_join_installed',
       'current_original_actual_lower_active_six_history_inlet_and_pressure_units_bound',
       'current_original_actual_first_micro_flat_collar_C5_source_join_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+(
    'current_original_actual_lower_active_join_quantitative_mixed_bounds_installed',)))


def selected_sc_binding(lower,loop):
    """Resolve the existing semantic leaf from its exact selected source.

    This adds a direct source binding without mutating an accepted ancestor
    graph or changing the original choice of sc.
    """
    raw=parent.linear.read(collar.NAME);rec=parent.linear.read(collar.RECEIPT)
    if not rec['all_passed'] or not rec[collar.GATE]:raise ValueError('Checked original selected collar required')
    record=raw['explicit_current_inner_exit_strict_collar'];saved=record['selected_first_phase_endpoint']
    bits=saved['lower_exact_mpf_tuple']
    if bits!=saved['upper_exact_mpf_tuple'] or bits!=lower.sc_bits:
        raise ValueError('Identical exact original sc singleton required')
    c=MPIntervalContext();c.dps=240
    ledger=parent.linear.read(PREFIX+'K1_ledger.json');norm=parent.linear.read(PREFIX+'physical_norm_family.json')
    gamma=collar.read_interval(c,ledger['gamma']);L=collar.read_interval(c,norm['selected_log_K_upper'])
    selected=c.mpf(collar.endpoints(1/(4*c.sqrt(10*L+c.ln(10/gamma))))[0])
    if list(selected._mpi_[0])!=bits or list(selected._mpi_[1])!=bits:
        raise ValueError('Actual original 240-digit selector differs')
    node=loop.g.nodes[loop.sc.node]
    if (node['operation']!='current_original_source_parameter' or node['name']!='sc'
        or node['definition']!='same original selected_first_phase_endpoint exact singleton'
        or node['source_family']!=lower.family or not node['Z_independent']
        or not node['quantity_is_exact_original_parameter_not_range_endpoint']):
        raise ValueError('The actual original active sc leaf is required')
    source_rec=parent.linear.read(node['current_source_receipt'])
    if sha(node['current_source_receipt'])!=node['current_source_receipt_sha256']:
        raise ValueError('Changed active parameter source receipt')
    for name in (collar.NAME,collar.RECEIPT,Path(collar.__file__).name):
        if source_rec['input_hashes'].get(name)!=sha(name):
            raise ValueError('Active parameter must retain original strict-collar ancestry '+name)
    if not 0<collar.endpoints(selected)[0]<c.mpf('.25'):
        raise ValueError('The actual selected comparison cutoff interval must be below one')
    margin=collar.read_interval(c,record['positive_margins']['subcollar_inside_original_rho4point1_continuation'])
    if collar.endpoints(margin)[0]<=0:raise ValueError('Actual analytic continuation collar required')
    return dict(active_sc_node=loop.sc.node,resolved_exact_sc_node=lower.parameters['sc'].node,
        source_data=collar.NAME,source_receipt=collar.RECEIPT,
        source_data_sha256=sha(collar.NAME),source_receipt_sha256=sha(collar.RECEIPT),
        source_path='explicit_current_inner_exit_strict_collar.selected_first_phase_endpoint',
        definition='lower_endpoint(1/(4*sqrt(10*selected_log_K_upper+log(10/gamma))))',
        original_selector_decimal_precision=240,exact_mpf_tuple=bits,
        exact_singleton_recomputed_not_interval_overlap=True,
        original_active_parameter_projection_binding=parent.active.current.ast_binding(left_source.original_left_inlet),
        comparison_alpha_exactly_one_on_zero_to_sc=True,
        rho4point1_common_analytic_continuation_margin_positive=True)


def physical_active_fields(owner):
    g=owner.g;a=Algebra(g,5);data=owner.modified.charts['first_micro']
    chart=data['source_chart'];P=owner.bridge.inlet.Pstar
    result=dict(data['physical'])
    result['Utheta']=a.scale(chart['modified_fields']['E'],P)
    result['pressure']=Algebra(g,3).scale(chart['modified_fields']['pressure'],g.mul(P,P))
    return result


class CurrentOriginalLowerActiveCommonSourceJoinFunctions:
    def __init__(self,charts_owner,require_checked=True):
        if type(charts_owner) is not previous.CurrentOriginalN1IntermediateOwnC3ChartFunctions or not charts_owner.acceptance_loaded:
            raise TypeError('Checked actual intermediate source owner required')
        charts_owner.assert_graph();self.previous=charts_owner;self.g=g=charts_owner.g
        self.family=copy.deepcopy(charts_owner.family);self.hashes=dict(charts_owner.hashes)
        self.lower=charts_owner.lower;self.core=self.lower.core
        self.modified=charts_owner.modified;self.loop=self.modified.loop
        self.bridge=self.modified.background.bridge;self.inlet=self.bridge.inlet
        if self.inlet.previous is not self.lower or self.bridge.core is not self.core:
            raise ValueError('Actual lower, inlet and bridge must retain the identical core owner')
        for module in (previous,bridge_source,interface,collar):
            rec=parent.linear.read(module.RECEIPT)
            if not rec['all_passed']:raise ValueError('Checked common-source dependency required')
            family=rec.get('source_family') or {key:rec[key] for key in self.family}
            if family!=self.family:raise ValueError('Same actual core/pressure source family required')
            for name,digest in rec['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        for name in (Path(left_source.__file__).name,interface.PREFIX+'current_core_first_interface_check.py',
                     PREFIX+'K1_ledger.json',PREFIX+'physical_norm_family.json',RECOVERY):
            self._pin(name,sha(name))
        theorem=parent.linear.read(interface.RECEIPT)
        self.integrated_theorem=copy.deepcopy(theorem['original_integrated_stress_free_core_equations'])
        self.primitive_theorem=copy.deepcopy(theorem['scaled_atoms_are_the_original_physical_primitives'])
        if not self.integrated_theorem['passed'] or not self.primitive_theorem['passed']:
            raise ValueError('Checked physical integral recovery required')
        self.sc_binding=selected_sc_binding(self.lower.lower,self.loop)
        self.lower_phase=self.lower.lower.fields['query'];self.active_phase=self.loop.charts['first_micro']['coordinate']
        self.left_endpoint=self.lower.lower.parameters['sc_half']
        self.right_endpoint=self.loop.charts['first_micro']['domain'][0]
        a=Algebra(g,5);fields=self.lower.fields
        self.lower_histories=bridge_source.histories(a,self.lower.radial_fields['Phi_rho_rows'][0],fields['U'],
            self.inlet.inlet_moments,self.lower_phase,self.bridge.hb,'lower_active_common')
        self.lower_history_mean=self.lower_histories['M']
        b=Algebra(g,3)
        self.lower_history_pressure=b.add(self.core.functions['P0'],b.scale(
            b.mul(b.mul(self.core.functions['F0'],self.core.functions['F0']),self.lower_histories['C']),
            self.lower.lower.fields['radius']))
        self.left_fields={key:fields[key] for key in ('U','U_y','m','m_y','Utheta','Q','Q_y')}
        self.left_fields['pressure']=self.lower.lower.fields['absolute_swirl_base_pressure']
        self.right_fields=physical_active_fields(self)
        def trace(group,variable,value):
            return {key:Algebra(g,len(rows(jet))-1).at(jet,variable,value) for key,jet in group.items()}
        self.left_trace=trace(self.left_fields,self.lower_phase,self.left_endpoint)
        self.right_trace=trace(self.right_fields,self.active_phase,self.right_endpoint)
        self.left_radius=b.at(b.fixed(self.lower.lower.fields['radius']),self.lower_phase,self.left_endpoint).value
        self.right_radius=b.at(b.fixed(self.loop.charts['first_micro']['radius']),self.active_phase,self.right_endpoint).value
        self.inlet_binding=dict(core_radius=self.core.rho,actual_radius=self.inlet.radius,
            six_core_moment_functions=self.inlet.moment_functions,six_rho4_inlet_functions=self.inlet.inlet_moments,
            comparison_incoming=self.bridge.initial7,lower_history_incoming=self.inlet.inlet_moments,
            normalized_P0=self.bridge.P0,physical_P0=self.core.functions['P0'],
            Pstar=self.inlet.Pstar,physical_axis_scale=g.unary('exp',g.mul(g.constant(2),self.core.parameters['logP'])))
        self.control_binding=dict(comparison=self.bridge.comparison_micro,actual=self.bridge.micro,
            lower_ratio=fields['ratio'],lower_phi=self.lower.radial_fields['Phi_rho_rows'][0],
            lower_psi=self.lower.radial_fields['Psi_rho_rows'][0],
            same_hb=self.bridge.hb,actual_chi=self.lower.lower.fields['chi'],
            active_loop={k:self.loop.charts['first_micro']['loop'][k] for k in
                ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi','flat_predicate','original_collar_upper')},
            active_finite_N=self.loop.charts['first_micro']['finite_N'])
        self.definition=dict(actual_unique_twenty_term_core_and_all_six_inlet_histories=True,
            direct_exact_active_sc_binding_without_changing_accepted_graph=True,
            normalized_P0_times_Pstar_squared_is_physical_core_P0=True,
            comparison_equals_actual_core_on_original_flat_alpha_interval=True,
            nonconstant_actual_chi_and_positive_hb_not_flattened=True,
            lower_and_prescribed_bridge_share_angular_axial_IVP=True,
            actual_finite_N_changes_zero_on_original_two_sided_collar=True,
            leading_U_Uy_mean_mean_y_swirl_C5_Q_Qy_C4_pressure_C3_join=True,
            actual_distinct_trace_functions_not_replaced_by_declared_zero=True,
            quantitative_mixed_bounds_other_joins_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(charts_owner,self.lower,self.core,self.modified,self.loop,self.bridge,self.inlet,g)
        self.acceptance_loaded=False
        if require_checked:
            raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not rec.get('all_passed') or not all(rec.get(k) for k in GATES) or any(rec.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['lower_active_common_source_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual lower/active same-source join required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed lower/active source '+name)
            self.acceptance_loaded=True
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed join source '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(sc_binding=self.sc_binding,lower_histories=self.lower_histories,
            lower_history_mean=self.lower_history_mean,lower_history_pressure=self.lower_history_pressure,
            left_fields=self.left_fields,right_fields=self.right_fields,left_trace=self.left_trace,right_trace=self.right_trace,
            left_radius=self.left_radius,right_radius=self.right_radius,inlet_binding=self.inlet_binding,
            control_binding=self.control_binding,integrated_theorem=self.integrated_theorem,
            primitive_theorem=self.primitive_theorem,definition=self.definition))
    def assert_graph(self):
        self.previous.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.previous,self.previous.lower,self.previous.lower.core,self.previous.modified,
                self.previous.modified.loop,self.previous.modified.background.bridge,self.previous.modified.background.bridge.inlet,self.g)):
            raise ValueError('Actual lower/active common sources changed')
    def resolve_sc(self,ref):
        self.assert_graph()
        if ref!=self.loop.sc:raise ValueError('Only the registered actual sc leaf is resolved here')
        return self.lower.lower.parameters['sc']
    def at_seam(self):
        self.assert_graph()
        return dict(left=self.left_trace,right=self.right_trace,
            actual_C5_Q4_pressure_C3_source_join=self.acceptance_loaded,
            quantitative_mixed_or_global_n1_admission=False)


def run(charts_owner):
    began=time.monotonic();owner=CurrentOriginalLowerActiveCommonSourceJoinFunctions(charts_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,lower_active_common_source_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_LOWER_ACTIVE_COMMON_SOURCE_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(charts_owner) reuses checked actual source owners')
