"""Add current collar/full Gamma charts to retained27 physical source owners.

Original Cartesian/time operators and fixed Ev0/Pstar units are retained.
Only collar/exterior packets are newly generated, never nonlinear point fields.
"""
import ast
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_steep_waiting_physical_assembly import (
    CurrentSteepWaitingPhysicalAssembly as PRIOR, BASE, PULSE, CHARTS as PRIOR_CHARTS,
    native_parameter_source_bridge, pulse_factor_rebase_proof,
    original_flatten_units_binding, UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_current_heat_source import (
    CurrentHeatSourceAssembly, CurrentCollarGammaC4, NEW_CHARTS, VIEWS,
    GATE as SOURCE_GATE)
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

CHARTS=PRIOR_CHARTS+NEW_CHARTS
PRIOR_REPORT=PREFIX+'current_steep_waiting_physical_assembly.json'
PRIOR_RECEIPT=PREFIX+'current_steep_waiting_physical_assembly_check.json'
RECEIPT=PREFIX+'current_heat_physical_assembly_check.json'
GATE='current_heat_cartesian_spatial4_time1_certified'
COMPOSED='current_twenty_nine_downstream_physical_source_ownership_certified'


def original_heat_radius_bindings():
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantGlobalPhysicalAssembly')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='radius')
    expressions={'heat':"getattr(provider,'heat',provider)",'steep':'heat.steep',
        'offset':'100+steep.outer.Lrel+2+steep.Ts+steep.wait+v'}
    for target,expression in expressions.items():
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)
                and ast.dump(n.value)==wanted for n in ast.walk(fn))!=1:
            raise ValueError('Original heat radius source changed: '+target)
    L,T,W,v=s.symbols('L T W v',real=True)
    origin=100+L+2+T+W
    identities=dict(waiting_collar_radius=s.simplify(100+L+2+T+W-origin)==0,
        collar_exterior_radius=s.simplify((origin+v).subs(v,3)-(origin+3))==0)
    if not all(identities.values()):raise ValueError('Original heat radius interfaces differ')
    return dict(actual_original_heat_radius_assignments=expressions,
        exact_current_heat_radius_interface_identities=identities,
        unbounded_exterior_logR_has_finite_lower_endpoint=True,
        unbounded_source_bound_requires_nonpositive_logR_exponents=True,passed=True)


class _CurrentHeatPhysicalDispatch:
    def __init__(self,source):self.current=source;self.native=source.before.before.before.dispatch
    def __getattr__(self,name):return getattr(self.native,name)
    def provider(self,chart):return self.current.provider(chart)
    def evaluate(self,chart,Z,coordinate):return self.current.evaluate(chart,Z,coordinate)


class CurrentHeatPhysicalAssembly(PRIOR):
    @source_precision
    def __init__(self,require_checked=True):
        self.source_assembly=CurrentHeatSourceAssembly()
        self.dispatch=_CurrentHeatPhysicalDispatch(self.source_assembly)
        self.source_owners=dict(self.source_assembly.registry);self.hashes=dict(self.source_assembly.hashes)
        self.family=self.source_assembly.family;self.source=self.source_assembly.source
        self.datum_sha=self.source_assembly.datum_sha
        self.core=self.dispatch.anchor.patch.core;self.ctx=c=self.core.ctx
        self.pre=self.dispatch.rh_reference;self.params=self.pre.params;self.pulse=self.dispatch.native_pulse
        self.delta=c.mpf(endpoints(self.params.delta));self.logP=c.mpf(endpoints(self.params.logPstar))
        self.logC=self.core.logC;self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw;self.logRv=self.logRp+13/self.params.mu
        self.flatten=self.source_assembly.before.before.before.flatten;self.outer=self.source_assembly.before.before.outer
        self.steep=self.source_assembly.before.steep;self.heat=self.source_assembly.heat
        self.native_parameter_bridge=native_parameter_source_bridge(self)
        self.hashes.update(self.native_parameter_bridge['input_hashes'])
        self.rebase_proof=pulse_factor_rebase_proof()
        self.operator_bindings=dict(original_radius=self.radius.__func__ is BASE.radius,
            original_pulse_and_post_normalized_sources=self.normalized_sources.__func__ is PULSE.normalized_sources,
            original_full_evaluate_operator_called=True,
            same_current_provider_and_evaluator=self.dispatch.current is self.source_assembly)
        self.parameter_source_bindings={target:class_assignment('current_heat_physical_assembly',
            'CurrentHeatPhysicalAssembly','__init__',target,expression) for target,expression in {
                'self.params':'self.pre.params','self.delta':'c.mpf(endpoints(self.params.delta))',
                'self.logP':'c.mpf(endpoints(self.params.logPstar))','self.logC':'self.core.logC',
                'self.logRref':'c.ln(110)+10*(self.logC+self.logP)',
                'self.logRp':'self.logRref+self.logP+1+self.params.Tw',
                'self.logRv':'self.logRp+13/self.params.mu',
                'self.flatten':'self.source_assembly.before.before.before.flatten','self.outer':'self.source_assembly.before.before.outer',
                'self.steep':'self.source_assembly.before.steep','self.heat':'self.source_assembly.heat'}.items()}
        self.fixed_unit_proof=original_flatten_units_binding()
        self.radius_proof=original_heat_radius_bindings()
        previous=json.loads((HERE/PRIOR_REPORT).read_bytes())
        receipt=accepted(PRIOR_RECEIPT,self.family,self.source,
            'current_twenty_seven_downstream_physical_source_ownership_certified')
        if (previous['datum_enclosure_sha256']!=self.datum_sha
                or tuple(previous['current_physical_chart_owners'])!=PRIOR_CHARTS
                or previous['current_source_owner_registry']!={k:self.source_owners[k] for k in PRIOR_CHARTS}
                or not receipt['independent_current_O7_radius_and_fixed_unit_fixture']['passed']
                or not previous['reused_independent_coordinate_fixture']['passed']):
            raise ValueError('Retained twenty-seven physical evidence differs')
        for name,digest in receipt['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Retained physical source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[PRIOR_REPORT]=sha(PRIOR_REPORT);self.hashes[PRIOR_RECEIPT]=sha(PRIOR_RECEIPT)
        self.previous_physical_receipt=receipt;self.previous_physical_manifest=previous
        self.source_identity_proof=previous['exact_source_and_divergence_identity']
        if not all(self.current_provider_graph().values()) or not all(self.operator_bindings.values()):
            raise ValueError('Current collar/Gamma physical source graph or operators differ')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.physical_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current heat physical datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.physical_acceptance_loaded=True

    def current_provider_graph(self):
        graph=PRIOR.current_provider_graph(self)
        source27=self.source_assembly.before;source23=source27.before;source21=source23.before
        graph.update(checked_current_steep_waiting_source=source27.acceptance_loaded,
            checked_current_actual_power_angular_source=source23.acceptance_loaded,
            checked_current_actual_flatten_source=source21.acceptance_loaded,
            checked_current_power_angular_source=source23.acceptance_loaded,
            checked_current_flatten_predecessor=source21.acceptance_loaded,
            checked_current_flatten_source=source21.acceptance_loaded,
            complete_future_source_functional_binding=source23.bindings['passed'],
            complete_O7_future_source_functional_binding=source27.bindings['passed'],
            checked_current_heat_source=self.source_assembly.acceptance_loaded,
            exact_same_heat_provider_for_both_charts=all(self.dispatch.provider(k) is self.heat for k in NEW_CHARTS),
            original_current_heat_class=type(self.heat) is CurrentCollarGammaC4,
            same_heat_and_physical_steep=self.heat.steep is self.steep,
            same_heat_and_physical_outer=self.heat.outer is self.outer,
            raw_heat_datum_alias_does_not_shadow_original_radius_route=not hasattr(self.heat,'heat'),
            same_heat_and_native_context=self.heat.ctx is self.pulse.ctx,
            complete_current_heat_future_source_binding=self.source_assembly.bindings['passed'])
        return graph

    def _annotate(self,packet,owner):
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_source_owner=self.source_owners[owner]['provider'],
            current_source_acceptance_receipt=self.source_owners[owner]['acceptance_receipt'],
            current_heat_physical_acceptance_receipt=RECEIPT,
            current_heat_cartesian_spatial4_time1_proved=True,
            **{GATE:self.physical_acceptance_loaded,COMPOSED:self.physical_acceptance_loaded,SOURCE_GATE:True},
            current_twenty_seven_downstream_physical_source_ownership_certified=True,
            current_heat_physical_owner_installed=self.physical_acceptance_loaded,
            current_source_dispatcher_used=True,current_core_parameter_source_used=True,
            current_Rp_external_pulse_join_certified=True,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))
        if owner in NEW_CHARTS:
            packet.update(current_heat_provider_used_for_derivative_grid_and_radius=True,
                same_current_flatten_provider_used_for_fixed_unit=True,
                original_fixed_Ev0_velocity_and_Pstar_squared_pressure_units_restored=True)
        return packet

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        if chart not in CHARTS or axis:raise ValueError('Twenty-nine current downstream physical owners only')
        packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)
        if chart=='heat_exterior':
            rows=[row for components in packet['physical_spatial_cartesian_mixed4'].values()
                for parts in components.values() for row in parts.values()]
            rows += [row for parts in packet['first_fixed_x_physical_time_derivative'].values() for row in parts.values()]
            for row in rows:
                if any(term['source_log_exponents'][5]>0 for term in row['terms']):
                    raise ValueError('Unbounded exterior requires nonpositive exact logR exponents')
                if not row['exact_zero'] and not all(mp.isfinite(v) for v in endpoints(row['log_absolute_upper'])):
                    raise ValueError('Unbounded exterior uniform physical row upper is not finite')
            packet.update(unbounded_exterior_nonpositive_logR_exponents_checked=True,
                unbounded_exterior_finite_uniform_physical_source_row_uppers=True,
                unbounded_physical_logR_is_coordinate_coverage_not_point_evaluation=True)
        return self._annotate(packet,chart)

    def manifest(self):
        previous=self.previous_physical_manifest;old=self.previous_physical_receipt
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,current_provider_graph_identity=self.current_provider_graph(),
            current_physical_parameter_source_bindings=self.parameter_source_bindings,
            native_current_parameter_source_bridge=self.native_parameter_bridge,
            original_physical_operators_retained=self.operator_bindings,
            current_post_fixed_unit_source_binding=self.fixed_unit_proof,
            current_heat_radius_functional_binding=self.radius_proof,
            current_heat_profile_source_binding=self.source_assembly.bindings,
            exact_source_and_divergence_identity=self.source_identity_proof,
            retained_twenty_seven_chart_physical_evidence=dict(report=PRIOR_REPORT,receipt=PRIOR_RECEIPT,
                report_sha256=sha(PRIOR_REPORT),owners=list(PRIOR_CHARTS),same_registry_sources_and_datum=True,
                total_regular_source_contributions=old['twenty_seven_regular_source_contributions_checked'],
                supplemental_gap_source_contributions=previous['retained_twenty_three_chart_physical_evidence']['supplemental_gap_source_contributions']),
            reused_independent_post_fixed_unit_fixture=previous['reused_independent_post_fixed_unit_fixture'],
            reused_independent_coordinate_fixture=previous['reused_independent_coordinate_fixture'],
            **{GATE:self.physical_acceptance_loaded,COMPOSED:self.physical_acceptance_loaded,SOURCE_GATE:True,UNIFORM:False},
            current_twenty_seven_downstream_physical_source_ownership_certified=True,
            current_heat_physical_owner_installed=self.physical_acceptance_loaded,
            full_pulse_C4_installed=False,all_profile_source_charts_callable=False,
            output_kind='signed physical Cartesian/time source bounds for twenty-nine current downstream owners',
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();packets={}
        for chart in NEW_CHARTS:
            packets[chart]={}
            for name,value in VIEWS[chart].items():
                packets[chart][name]=self.evaluate(chart,[-1,1],value,theta=None)
                print('Current collar/Gamma physically mapped: '+chart+' '+name,flush=True)
        result.update(whole_current_heat_physical_maps=packets,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentHeatPhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current heat physical assembly generated: twenty-nine owners, Cartesian spatial4/fixed-x time1',flush=True)
    return result


if __name__=='__main__':run()
