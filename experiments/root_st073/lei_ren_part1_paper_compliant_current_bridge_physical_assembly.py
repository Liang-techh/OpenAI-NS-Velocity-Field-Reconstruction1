"""Current33 physical source coverage; bridge interfaces admitted separately."""
import ast
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_core_physical_assembly import (
    CurrentCorePhysicalAssembly as PRIOR,CurrentCoreSourceOverlay,
    _CurrentCorePhysicalDispatch,CHARTS as PRIOR_CHARTS,BASE,
    UNIFORM,SCOPES,OPEN,HERE,PREFIX,sha)
from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import (
    CurrentActualBridgeMixedC4,VIEWS as SOURCE_VIEWS)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

BRIDGES={'bridge_first':'first','bridge_second':'second','bridge_macro':'macro'}
CHARTS=PRIOR_CHARTS+tuple(BRIDGES)
RECEIPT=PREFIX+'current_bridge_physical_assembly_check.json'
GATE='current_actual_three_bridge_charts_physical_spatial4_time1_certified'
COMPOSED='current_thirty_three_physical_source_ownership_certified'
COVERAGE='all_33_current_source_charts_physical_spatial4_time1_mapped'
PRIOR_REPORT=PREFIX+'current_core_physical_assembly.json'
PRIOR_RECEIPT=PREFIX+'current_core_physical_assembly_check.json'
SOURCE_RECEIPT=PREFIX+'current_actual_bridge_mixed_C4_check.json'


class CurrentBridgeSourceOverlay(CurrentCoreSourceOverlay):
    def __init__(self,before,core,bridge):
        super().__init__(before,core);self.bridge=bridge
        self.registry.update({chart:dict(provider=PREFIX+'current_actual_bridge_mixed_C4.CurrentActualBridgeMixedC4',
            method='evaluate',coverage_coordinate='s=log(R/Ra)/hb' if phase!='macro' else 'fraction of frozen macro interval',
            domain='[1,2]' if phase=='second' else '[0,1]',acceptance_receipt=SOURCE_RECEIPT)
            for chart,phase in BRIDGES.items()})

    def provider(self,chart):return self.bridge if chart in BRIDGES else super().provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart not in BRIDGES:return super().evaluate(chart,Z,coordinate)
        packet=self.bridge.evaluate(Z,coordinate,BRIDGES[chart])
        return dict(chart=chart,source_packet=packet,
            actual_five_defect_family_sha256=self.bridge.family,implicit_source_sha256=self.bridge.source,
            datum_enclosure_sha256=self.bridge.datum_sha,source_provider=self.registry[chart]['provider'],
            acceptance_receipt=SOURCE_RECEIPT,physical_mixed_grids={k:v for k,v in packet.items()
                if k.startswith('physical_velocity_pressure_') or k.startswith('physical_five_primitive_')})


def original_bridge_mapper_bindings():
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantGlobalPhysicalAssembly')
    methods={n.name:n for n in cls.body if isinstance(n,ast.FunctionDef)}
    expressions={'normalized_sources':("micro_terms(packet,name,phase_to_y=chart in MICRO)",
        "packet.get('actual_parent_axial5_packet')"),
        'evaluate':("self.dispatch.evaluate(chart,Z,coordinate)",
        "self.normalized_sources(chart,Z,coordinate,packet,provider,logR)")}
    for method,calls in expressions.items():
        for expression in calls:
            wanted=ast.dump(ast.parse(expression,mode='eval').body)
            if sum(isinstance(n,ast.Call) and ast.dump(n)==wanted for n in ast.walk(methods[method]))!=1:
                raise ValueError('Original bridge physical call changed: '+expression)
    return dict(current_source_packet_and_parent_consumed_before_mapping=True,
        exact_phase_to_logR_conversion_before_bound=True,
        unchanged_original_radius_normalized_sources_and_full_evaluate=True,passed=True)


class CurrentBridgePhysicalAssembly(PRIOR):
    @source_precision
    def __init__(self,require_checked=True,base=None):
        if base is None:PRIOR.__init__(self,require_checked=True)
        else:
            if type(base) is not PRIOR or not base.core_acceptance_loaded:
                raise ValueError('An already checked exact current30 base required')
            self.__dict__.update(base.__dict__);self.hashes=dict(base.hashes)
            self.source_owners=dict(base.source_owners)
        self.bridge_source=CurrentActualBridgeMixedC4(owner=self)
        self.core_source=CurrentBridgeSourceOverlay(self.source_assembly,self.core,self.bridge_source)
        self.dispatch=_CurrentCorePhysicalDispatch(self.core_source,self.dispatch.native)
        self.source_owners=dict(self.core_source.registry)
        self.hashes.update(self.bridge_source.hashes)
        self.bridge_bindings=original_bridge_mapper_bindings()
        self.overlay_bindings={target:class_assignment('current_bridge_physical_assembly','CurrentBridgePhysicalAssembly',
            '__init__',target,expression) for target,expression in {
                'self.bridge_source':'CurrentActualBridgeMixedC4(owner=self)',
                'self.core_source':'CurrentBridgeSourceOverlay(self.source_assembly,self.core,self.bridge_source)',
                'self.dispatch':'_CurrentCorePhysicalDispatch(self.core_source,self.dispatch.native)'}.items()}
        previous=json.loads((HERE/PRIOR_REPORT).read_bytes())
        old=accepted(PRIOR_RECEIPT,self.family,self.source,'current_thirty_source_physical_ownership_certified')
        if (tuple(previous['current_physical_chart_owners'])!=PRIOR_CHARTS
                or previous['current_source_owner_registry']!={k:self.source_owners[k] for k in PRIOR_CHARTS}
                or previous['datum_enclosure_sha256']!=self.datum_sha):
            raise ValueError('Retained current30 registry or datum differs')
        theorem=accepted(PREFIX+'global_physical_assembly_check.json',self.family,self.source,
            'all_33_original_source_charts_physical_spatial4_time1_mapped')
        for record in (old,theorem):
            for name,digest in record['input_hashes'].items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current physical source conflict: '+name)
                self.hashes[name]=digest
        for name in (PRIOR_REPORT,PRIOR_RECEIPT,PREFIX+'global_physical_assembly_check.json'):
            self.hashes[name]=sha(name)
        self.retained_core_manifest=previous;self.retained_core_receipt=old
        self.original_physical_theorem=theorem
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        failed=[k for k,v in self.current_provider_graph().items() if not v]
        if failed:raise ValueError('Current33 physical object graph differs: '+', '.join(failed))
        self.bridge_physical_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current33 physical datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.bridge_physical_acceptance_loaded=True

    def current_provider_graph(self):
        graph=PRIOR.current_provider_graph(self)
        if not hasattr(self,'bridge_source'):return graph
        graph.update(checked_current_bridge_derivative_source=self.bridge_source.acceptance_loaded,
            exact_same_provider_for_all_three_bridges=all(self.dispatch.provider(k) is self.bridge_source for k in BRIDGES),
            same_physical_bridge_owner=self.bridge_source.owner is self,
            same_nested_core_and_bridge_context=self.bridge_source.core is self.core and self.bridge_source.ctx is self.ctx,
            same_nested_current_history=self.bridge_source.history is self.dispatch.anchor.patch.reference_mixed.long_mixed.history,
            checked_current_core_source=self.core_acceptance_loaded,
            current_bridge_graph_bound=all(self.bridge_source.current_provider_graph().values()),
            original_radius_method=self.radius.__func__ is BASE.radius,
            original_factored_bridge_source_conversion=self.normalized_sources.__func__ is PRIOR.normalized_sources)
        return graph

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        if chart not in CHARTS or (axis and chart!='core'):raise ValueError('Current33 original chart/axis mode required')
        if chart in BRIDGES:
            packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)
            packet.update(datum_enclosure_sha256=self.datum_sha,current_source_owner=self.source_owners[chart]['provider'],
                current_source_acceptance_receipt=SOURCE_RECEIPT,current_source_dispatcher_used=True,
                current_core_parameter_source_used=True,current_actual_bridge_source_used=True,
                **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),**{UNIFORM:False},full_pulse_C4_installed=False,
                core_inner_annulus_interfaces_certified=False)
        else:packet=PRIOR.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)
        packet.update(current_bridge_physical_acceptance_receipt=RECEIPT,
            current_bridge_physical_operators_source_bound=True,
            **{GATE:self.bridge_physical_acceptance_loaded,COMPOSED:self.bridge_physical_acceptance_loaded,
                COVERAGE:self.bridge_physical_acceptance_loaded},current_thirty_source_physical_ownership_certified=True)
        return packet

    def manifest(self):
        old=self.retained_core_receipt
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,current_provider_graph_identity=self.current_provider_graph(),
            current_bridge_overlay_constructor_bindings=self.overlay_bindings,
            original_bridge_mapper_source_bindings=self.bridge_bindings,
            retained_current30_physical_evidence=dict(report=PRIOR_REPORT,receipt=PRIOR_RECEIPT,
                regular_contributions=old['thirty_regular_source_contributions_checked'],
                supplemental_axis_contributions=old['supplemental_axis_contributions_checked'],
                supplemental_gap_contributions=old['retained_twenty_nine_physical_evidence']['supplemental_gap_contributions'],
                same_registry_datum_and_checked_current30=True),
            reused_independent_unchanged_physical_fixtures={k:self.original_physical_theorem[k] for k in (
                'independent_implicit_physical_fixture','independent_micro_scale_fixture','source_and_divergence_identities')},
            **{GATE:self.bridge_physical_acceptance_loaded,COMPOSED:self.bridge_physical_acceptance_loaded,
                COVERAGE:self.bridge_physical_acceptance_loaded},current_thirty_source_physical_ownership_certified=True,
            core_inner_annulus_interfaces_certified=False,full_pulse_C4_installed=False,**{UNIFORM:False},
            **dict.fromkeys(SCOPES,False),**{k:False for k in OPEN if k!=COVERAGE},input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();views={}
        for name,(phase,value) in SOURCE_VIEWS.items():
            chart='bridge_'+phase
            views[name]=self.evaluate(chart,[-1,1],value,theta=None)
            print('Current actual bridge physical map: '+name,flush=True)
        result.update(whole_current_bridge_physical_maps=views,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentBridgePhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current33 physical source coverage generated; full point/global/temporal gates remain open',flush=True)
    return result


if __name__=='__main__':run()
