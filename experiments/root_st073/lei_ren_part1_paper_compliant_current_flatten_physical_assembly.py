"""Add current O5 flatten to the retained twenty physical source owners.

Original Cartesian/time operators and native fixed units are unchanged.
Only the new flatten maps are generated. Output is directed source bounds,
not production point fields or a uniform native pulse interface certificate.
"""
import ast
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_pulse_physical_assembly import (
    CurrentPulsePhysicalAssembly as PULSE, BASE, CHARTS as PRIOR_CHARTS,
    PULSE_CHARTS, UNIFORM, SCOPES, OPEN, sha, HERE, PREFIX,
    pulse_factor_rebase_proof, native_parameter_source_bridge)
from lei_ren_part1_paper_compliant_current_pulse_flatten_source import (
    CurrentPulseFlattenSourceAssembly, CurrentFlattenMixedC4,
    GATE as SOURCE_GATE, RECEIPT as SOURCE_RECEIPT)
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

CHARTS=PRIOR_CHARTS+('flatten',)
PRIOR_REPORT=PREFIX+'current_pulse_physical_assembly.json'
PRIOR_RECEIPT=PREFIX+'current_pulse_physical_assembly_check.json'
RECEIPT=PREFIX+'current_flatten_physical_assembly_check.json'
GATE='current_flatten_cartesian_spatial4_time1_certified'
COMPOSED='current_twenty_one_downstream_physical_source_ownership_certified'


class _CurrentFlattenPhysicalDispatch:
    """One accepted source assembly handles both provider and evaluate calls."""
    def __init__(self,source):self.current=source;self.native=source.dispatch
    def __getattr__(self,name):return getattr(self.native,name)
    def provider(self,chart):return self.current.provider(chart)
    def evaluate(self,chart,Z,coordinate):return self.current.evaluate(chart,Z,coordinate)


def original_flatten_units_binding():
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantGlobalPhysicalAssembly')
    methods={n.name:n for n in cls.body if isinstance(n,ast.FunctionDef)}
    branch=next(n for n in ast.walk(methods['radius']) if isinstance(n,ast.If)
        and ast.unparse(n.test)=="chart == 'flatten'")
    if not any(isinstance(n,ast.Assign) and ast.unparse(n.value)=='v'
            and any(ast.unparse(t)=='offset' for t in n.targets) for n in branch.body):
        raise ValueError('Original flatten offset changed')
    post=next(n for n in ast.walk(methods['radius']) if isinstance(n,ast.If)
        and ast.unparse(n.test)=='chart in POST')
    expected=ast.dump(ast.parse('self.logRp+13/self.params.mu+offset',mode='eval').body)
    if not any(isinstance(n,ast.Return) and isinstance(n.value,ast.Tuple)
            and ast.dump(n.value.elts[0])==expected for n in post.body):
        raise ValueError('Original flatten radius source changed')
    wanted={
        'flatten':"self.dispatch.provider('flatten')",
        'logs[4]':'self.logP+sum(flatten.logEv2_parts.values(),c.mpf(0))/2'}
    rows={}
    for target,expression in wanted.items():
        values=[n.value for n in ast.walk(methods['normalized_sources']) if isinstance(n,ast.Assign)
            and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(v)==ast.dump(ast.parse(expression,mode='eval').body) for v in values)!=1:
            raise ValueError('Original flatten fixed unit changed: '+target)
        rows[target]=expression
    expression='amplitudes.update({label:{4:mp.mpf(1)} for label in (UZ,UT,UR)})'
    if sum(ast.dump(n)==ast.dump(ast.parse(expression,mode='eval').body)
            for n in ast.walk(methods['normalized_sources']) if isinstance(n,ast.Call))!=1:
        raise ValueError('Original all-component Ev0 restoration changed')
    return dict(original_flatten_radius='logRp+13/mu+t',actual_unit_assignments=rows,
        original_velocity_amplitude_update=expression,
        full_derivative_grids_restored_with_fixed_Ev0=True,
        absolute_pressure_Pstar_squared_retained=True,
        terminal_history_metadata_not_used_as_derivative_grid=True,passed=True)


class CurrentFlattenPhysicalAssembly(PULSE):
    @source_precision
    def __init__(self,require_checked=True):
        self.source=CurrentPulseFlattenSourceAssembly()
        self.dispatch=_CurrentFlattenPhysicalDispatch(self.source)
        self.source_owners=dict(self.source.registry);self.hashes=dict(self.source.hashes)
        self.family=self.source.family;self.source_sha=self.source.source;self.datum_sha=self.source.datum_sha
        # The inherited provider graph uses source as the family digest. Keep
        # the actual source assembly separately before setting that field.
        self.source_assembly=self.source;self.source=self.source_sha
        self.core=self.dispatch.anchor.patch.core;self.ctx=c=self.core.ctx
        self.pre=self.dispatch.rh_reference;self.params=self.pre.params;self.pulse=self.dispatch.native_pulse
        self.delta=c.mpf(endpoints(self.params.delta))
        self.logP=c.mpf(endpoints(self.params.logPstar));self.logC=self.core.logC
        self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw
        self.logRv=self.logRp+13/self.params.mu
        self.flatten=self.source_assembly.flatten
        self.native_parameter_bridge=native_parameter_source_bridge(self)
        self.hashes.update(self.native_parameter_bridge['input_hashes'])
        self.rebase_proof=pulse_factor_rebase_proof()
        self.operator_bindings=dict(original_radius=self.radius.__func__ is BASE.radius,
            original_pulse_and_post_normalized_sources=self.normalized_sources.__func__ is PULSE.normalized_sources,
            original_full_evaluate_operator_called=True,
            same_current_provider_and_evaluator=self.dispatch.current is self.source_assembly)
        self.parameter_source_bindings={target:class_assignment('current_flatten_physical_assembly',
            'CurrentFlattenPhysicalAssembly','__init__',target,expression) for target,expression in {
                'self.params':'self.pre.params','self.delta':'c.mpf(endpoints(self.params.delta))',
                'self.logP':'c.mpf(endpoints(self.params.logPstar))','self.logC':'self.core.logC',
                'self.logRref':'c.ln(110)+10*(self.logC+self.logP)',
                'self.logRp':'self.logRref+self.logP+1+self.params.Tw',
                'self.logRv':'self.logRp+13/self.params.mu'}.items()}
        self.flatten_unit_proof=original_flatten_units_binding()
        previous=json.loads((HERE/PRIOR_REPORT).read_bytes())
        receipt=accepted(PRIOR_RECEIPT,self.family,self.source,
            'current_twenty_downstream_physical_source_ownership_certified')
        if (previous['datum_enclosure_sha256']!=self.datum_sha
                or tuple(previous['current_physical_chart_owners'])!=PRIOR_CHARTS
                or previous['current_source_owner_registry']!={k:self.source_owners[k] for k in PRIOR_CHARTS}
                or not receipt['independent_pulse_fixed_unit_and_coordinate_fixture']['passed']
                or not receipt['unchanged_independent_cartesian_coordinate_fixture']['passed']):
            raise ValueError('Retained twenty-owner physical evidence differs')
        for name,digest in receipt['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Retained source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[PRIOR_REPORT]=sha(PRIOR_REPORT);self.hashes[PRIOR_RECEIPT]=sha(PRIOR_RECEIPT)
        self.previous_physical_receipt=receipt;self.previous_physical_manifest=previous
        self.source_identity_proof=previous['exact_source_and_divergence_identity']
        if not all(self.current_provider_graph().values()) or not all(self.operator_bindings.values()):
            raise ValueError('Current flatten physical source graph or operators differ')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.physical_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current flatten physical datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.physical_acceptance_loaded=True

    def current_provider_graph(self):
        return dict(**PULSE.current_provider_graph(self),
            checked_current_flatten_source=self.source_assembly.acceptance_loaded,
            exact_current_flatten_provider=self.dispatch.provider('flatten') is self.flatten,
            current_flatten_class=type(self.flatten) is CurrentFlattenMixedC4,
            same_flatten_native_pulse=self.flatten.pulse is self.pulse,
            same_flatten_native_context=self.flatten.ctx is self.pulse.ctx,
            current_mapper_context_is_core_context=self.ctx is self.core.ctx,
            same_flatten_native_mu=self.flatten.mu is self.pulse.mu,
            same_flatten_native_delta=self.flatten.delta is self.pulse.delta,
            exact_parameter_source_before_copy_checks=self.native_parameter_bridge['passed'],
            mapper_native_mu_copy=endpoints(self.params.mu)==endpoints(self.flatten.mu),
            mapper_native_delta_copy=endpoints(self.delta)==endpoints(self.flatten.delta),
            same_flatten_native_U_object=self.flatten.U is self.pulse.high.constants['U'],
            same_flatten_absolute_datum=self.flatten.inlet.datum is self.pre.datum or (
                self.flatten.inlet.datum.definition==self.pre.datum.definition
                and self.flatten.inlet.datum.datum_sha==self.datum_sha))

    def _annotate(self,packet,owner):
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_source_owner=self.source_owners[owner]['provider'],
            current_source_acceptance_receipt=self.source_owners[owner]['acceptance_receipt'],
            current_flatten_physical_acceptance_receipt=RECEIPT,
            current_flatten_cartesian_spatial4_time1_proved=True,
            **{GATE:self.physical_acceptance_loaded,COMPOSED:self.physical_acceptance_loaded},
            current_twenty_downstream_physical_source_ownership_certified=True,
            current_pulse_terminal_flatten_source_ownership_certified=True,
            current_flatten_physical_owner_installed=self.physical_acceptance_loaded,
            current_source_dispatcher_used=True,current_core_parameter_source_used=True,
            current_Rp_external_pulse_join_certified=True,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))
        if owner in PULSE_CHARTS:
            packet.update(pulse_correlated_radius_velocity_factors_combined_before_bounds=True,
                pulse_source_log_base_meanings=['unused','logPstar','logRp/2','unused',
                    'logPstar-mu*offset','logR','log2'],pulse_original_fixed_units_exactly_rebased=True)
        if owner=='flatten':packet.update(current_flatten_provider_used_for_derivative_grid_and_unit=True,
            original_fixed_Ev0_velocity_and_Pstar_squared_pressure_units_restored=True)
        return packet

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        if chart not in CHARTS or axis:raise ValueError('Twenty-one current downstream physical owners only')
        return self._annotate(BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False),chart)

    def manifest(self):
        old=self.previous_physical_receipt
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,current_provider_graph_identity=self.current_provider_graph(),
            current_physical_parameter_source_bindings=self.parameter_source_bindings,
            native_current_parameter_source_bridge=self.native_parameter_bridge,
            original_physical_operators_retained=self.operator_bindings,
            current_flatten_fixed_unit_source_binding=self.flatten_unit_proof,
            exact_source_and_divergence_identity=self.source_identity_proof,
            retained_twenty_chart_physical_evidence=dict(report=PRIOR_REPORT,receipt=PRIOR_RECEIPT,
                report_sha256=sha(PRIOR_REPORT),owners=list(PRIOR_CHARTS),same_registry_sources_and_datum=True,
                total_regular_source_contributions=old['twenty_owner_regular_source_contributions_checked'],
                supplemental_gap_source_contributions=old['new_current_pulse_spatial_and_time_source_contributions_checked']['pulse_gap_overlap']),
            reused_independent_coordinate_fixture=old['unchanged_independent_cartesian_coordinate_fixture'],
            reused_independent_pulse_unit_fixture=old['independent_pulse_fixed_unit_and_coordinate_fixture'],
            **{GATE:self.physical_acceptance_loaded,COMPOSED:self.physical_acceptance_loaded,UNIFORM:False},
            current_twenty_downstream_physical_source_ownership_certified=True,
            current_pulse_terminal_flatten_source_ownership_certified=True,
            current_flatten_physical_owner_installed=self.physical_acceptance_loaded,
            full_pulse_C4_installed=False,all_profile_source_charts_callable=False,
            output_kind='signed physical Cartesian/time source bounds for twenty-one current downstream owners',
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();packets={}
        for name,t in (('inlet',0),('whole_domain',[0,100]),('exit',100)):
            packets[name]=self.evaluate('flatten',[-1,1],t,theta=None)
            print('Current flatten physically mapped: '+name,flush=True)
        result.update(whole_current_flatten_physical_maps=packets,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentFlattenPhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current flatten physical assembly generated: twenty-one owners, Cartesian spatial4/fixed-x time1',flush=True)
    return result


if __name__=='__main__':run()
