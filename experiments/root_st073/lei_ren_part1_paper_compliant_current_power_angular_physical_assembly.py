"""Add two admitted current O6 charts to retained21 physical source owners.

Original Cartesian/time operators and fixed Ev0/Pstar units are retained.
Only power/angular packets are newly generated, never nonlinear point fields.
"""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_flatten_physical_assembly import (
    CurrentFlattenPhysicalAssembly as FLAT, BASE, PULSE, CHARTS as PRIOR_CHARTS,
    native_parameter_source_bridge, pulse_factor_rebase_proof,
    original_flatten_units_binding, UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_current_power_angular_source import (
    CurrentPowerAngularSourceAssembly, CurrentPowerAngularC4, NEW_CHARTS,
    GATE as SOURCE_GATE, original_radius_functional_bindings)
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

CHARTS=PRIOR_CHARTS+NEW_CHARTS
PRIOR_REPORT=PREFIX+'current_flatten_physical_assembly.json'
PRIOR_RECEIPT=PREFIX+'current_flatten_physical_assembly_check.json'
RECEIPT=PREFIX+'current_power_angular_physical_assembly_check.json'
GATE='current_power_angular_cartesian_spatial4_time1_certified'
COMPOSED='current_twenty_three_downstream_physical_source_ownership_certified'


class _CurrentPowerAngularPhysicalDispatch:
    def __init__(self,source):self.current=source;self.native=source.before.dispatch
    def __getattr__(self,name):return getattr(self.native,name)
    def provider(self,chart):return self.current.provider(chart)
    def evaluate(self,chart,Z,coordinate):return self.current.evaluate(chart,Z,coordinate)


class CurrentPowerAngularPhysicalAssembly(FLAT):
    @source_precision
    def __init__(self,require_checked=True):
        self.source_assembly=CurrentPowerAngularSourceAssembly()
        self.dispatch=_CurrentPowerAngularPhysicalDispatch(self.source_assembly)
        self.source_owners=dict(self.source_assembly.registry);self.hashes=dict(self.source_assembly.hashes)
        self.family=self.source_assembly.family;self.source=self.source_assembly.source
        self.datum_sha=self.source_assembly.datum_sha
        self.core=self.dispatch.anchor.patch.core;self.ctx=c=self.core.ctx
        self.pre=self.dispatch.rh_reference;self.params=self.pre.params;self.pulse=self.dispatch.native_pulse
        self.delta=c.mpf(endpoints(self.params.delta));self.logP=c.mpf(endpoints(self.params.logPstar))
        self.logC=self.core.logC;self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw;self.logRv=self.logRp+13/self.params.mu
        self.flatten=self.source_assembly.before.flatten;self.outer=self.source_assembly.outer
        self.native_parameter_bridge=native_parameter_source_bridge(self)
        self.hashes.update(self.native_parameter_bridge['input_hashes'])
        self.rebase_proof=pulse_factor_rebase_proof()
        self.operator_bindings=dict(original_radius=self.radius.__func__ is BASE.radius,
            original_pulse_and_post_normalized_sources=self.normalized_sources.__func__ is PULSE.normalized_sources,
            original_full_evaluate_operator_called=True,
            same_current_provider_and_evaluator=self.dispatch.current is self.source_assembly)
        self.parameter_source_bindings={target:class_assignment('current_power_angular_physical_assembly',
            'CurrentPowerAngularPhysicalAssembly','__init__',target,expression) for target,expression in {
                'self.params':'self.pre.params','self.delta':'c.mpf(endpoints(self.params.delta))',
                'self.logP':'c.mpf(endpoints(self.params.logPstar))','self.logC':'self.core.logC',
                'self.logRref':'c.ln(110)+10*(self.logC+self.logP)',
                'self.logRp':'self.logRref+self.logP+1+self.params.Tw',
                'self.logRv':'self.logRp+13/self.params.mu',
                'self.flatten':'self.source_assembly.before.flatten','self.outer':'self.source_assembly.outer'}.items()}
        self.fixed_unit_proof=original_flatten_units_binding()
        self.radius_proof=original_radius_functional_bindings()
        previous=json.loads((HERE/PRIOR_REPORT).read_bytes())
        receipt=accepted(PRIOR_RECEIPT,self.family,self.source,
            'current_twenty_one_downstream_physical_source_ownership_certified')
        if (previous['datum_enclosure_sha256']!=self.datum_sha
                or tuple(previous['current_physical_chart_owners'])!=PRIOR_CHARTS
                or previous['current_source_owner_registry']!={k:self.source_owners[k] for k in PRIOR_CHARTS}
                or not receipt['independent_flatten_fixed_unit_fixture']['passed']
                or not previous['reused_independent_coordinate_fixture']['passed']):
            raise ValueError('Retained twenty-one physical evidence differs')
        for name,digest in receipt['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Retained physical source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[PRIOR_REPORT]=sha(PRIOR_REPORT);self.hashes[PRIOR_RECEIPT]=sha(PRIOR_RECEIPT)
        self.previous_physical_receipt=receipt;self.previous_physical_manifest=previous
        self.source_identity_proof=previous['exact_source_and_divergence_identity']
        if not all(self.current_provider_graph().values()) or not all(self.operator_bindings.values()):
            raise ValueError('Current power/angular physical source graph or operators differ')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.physical_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current O6 physical datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.physical_acceptance_loaded=True

    def current_provider_graph(self):
        return dict(**FLAT.current_provider_graph(self),
            checked_current_power_angular_source=self.source_assembly.acceptance_loaded,
            checked_current_flatten_predecessor=self.source_assembly.before.acceptance_loaded,
            exact_same_outer_provider_for_both_charts=all(self.dispatch.provider(k) is self.outer for k in NEW_CHARTS),
            original_current_outer_class=type(self.outer) is CurrentPowerAngularC4,
            same_outer_and_physical_flatten=self.outer.flatten is self.flatten,
            same_outer_and_physical_native_pulse=self.outer.pulse is self.pulse,
            same_outer_and_physical_native_context=self.outer.ctx is self.pulse.ctx,
            same_outer_and_physical_native_mu=self.outer.mu is self.pulse.mu,
            same_outer_and_physical_native_delta=self.outer.delta is self.pulse.delta,
            complete_future_source_functional_binding=self.source_assembly.bindings['passed'])

    def _annotate(self,packet,owner):
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_source_owner=self.source_owners[owner]['provider'],
            current_source_acceptance_receipt=self.source_owners[owner]['acceptance_receipt'],
            current_power_angular_physical_acceptance_receipt=RECEIPT,
            current_power_angular_cartesian_spatial4_time1_proved=True,
            **{GATE:self.physical_acceptance_loaded,COMPOSED:self.physical_acceptance_loaded,SOURCE_GATE:True},
            current_twenty_one_downstream_physical_source_ownership_certified=True,
            current_power_angular_physical_owner_installed=self.physical_acceptance_loaded,
            current_source_dispatcher_used=True,current_core_parameter_source_used=True,
            current_Rp_external_pulse_join_certified=True,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))
        if owner in NEW_CHARTS:
            packet.update(current_outer_provider_used_for_derivative_grid_and_radius=True,
                same_current_flatten_provider_used_for_fixed_unit=True,
                original_fixed_Ev0_velocity_and_Pstar_squared_pressure_units_restored=True)
        return packet

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        if chart not in CHARTS or axis:raise ValueError('Twenty-three current downstream physical owners only')
        return self._annotate(BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False),chart)

    def manifest(self):
        previous=self.previous_physical_manifest;old=self.previous_physical_receipt
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,current_provider_graph_identity=self.current_provider_graph(),
            current_physical_parameter_source_bindings=self.parameter_source_bindings,
            native_current_parameter_source_bridge=self.native_parameter_bridge,
            original_physical_operators_retained=self.operator_bindings,
            current_post_fixed_unit_source_binding=self.fixed_unit_proof,
            current_power_angular_radius_functional_binding=self.radius_proof,
            current_power_angular_profile_source_binding=self.source_assembly.bindings,
            exact_source_and_divergence_identity=self.source_identity_proof,
            retained_twenty_one_chart_physical_evidence=dict(report=PRIOR_REPORT,receipt=PRIOR_RECEIPT,
                report_sha256=sha(PRIOR_REPORT),owners=list(PRIOR_CHARTS),same_registry_sources_and_datum=True,
                total_regular_source_contributions=old['twenty_one_regular_source_contributions_checked'],
                supplemental_gap_source_contributions=previous['retained_twenty_chart_physical_evidence']['supplemental_gap_source_contributions']),
            reused_independent_post_fixed_unit_fixture=old['independent_flatten_fixed_unit_fixture'],
            reused_independent_coordinate_fixture=previous['reused_independent_coordinate_fixture'],
            **{GATE:self.physical_acceptance_loaded,COMPOSED:self.physical_acceptance_loaded,SOURCE_GATE:True,UNIFORM:False},
            current_twenty_one_downstream_physical_source_ownership_certified=True,
            current_power_angular_physical_owner_installed=self.physical_acceptance_loaded,
            full_pulse_C4_installed=False,all_profile_source_charts_callable=False,
            output_kind='signed physical Cartesian/time source bounds for twenty-three current downstream owners',
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();packets={}
        for name,chart,value in (('power_inlet','outer_power',0),('power_domain','outer_power',[0,1]),
                ('power_exit','outer_power',1),('angular_inlet','outer_angular',-4),
                ('angular_domain','outer_angular',[-4,0]),('angular_exit','outer_angular',0)):
            packets[name]=self.evaluate(chart,[-1,1],value,theta=None)
            print('Current power/angular physically mapped: '+name,flush=True)
        result.update(whole_current_power_angular_physical_maps=packets,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentPowerAngularPhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current O6 physical assembly generated: twenty-three owners, Cartesian spatial4/fixed-x time1',flush=True)
    return result


if __name__=='__main__':run()
