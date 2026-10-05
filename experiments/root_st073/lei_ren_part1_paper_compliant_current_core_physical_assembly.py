"""Current nested analytic core/axis added to admitted29 source bounds.

Axis is a mode of core, not an additional owner. Original nonlinear core
enclosures and nonsingular Cartesian/time operators are retained. No point
selection, bridge interface, global tensor or temporal recursion is implied.
"""
import ast
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_heat_physical_assembly import (
    CurrentHeatPhysicalAssembly as PRIOR, BASE, CHARTS as PRIOR_CHARTS,
    UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_core_physical_field import CompliantCorePhysicalField as CORE
from lei_ren_part1_paper_compliant_actual_bridge_switch import CoreContext
from lei_ren_part1_paper_compliant_core_physical_field_check import functional_recovery_proof
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

CHARTS=('core',)+PRIOR_CHARTS
RECEIPT=PREFIX+'current_core_physical_assembly_check.json'
GATE='current_core_axis_cartesian_spatial4_time1_certified'
COMPOSED='current_thirty_source_physical_ownership_certified'
SOURCE_GATE='current_core_source_ownership_certified'
PRIOR_REPORT=PREFIX+'current_heat_physical_assembly.json'
PRIOR_RECEIPT=PREFIX+'current_heat_physical_assembly_check.json'
THEOREM=PREFIX+'core_physical_field_check.json'
VIEWS={'whole_core':([0,4],False),'core_axis_inlet':(0,False),'core_exit':(4,False),'whole_axis':(0,True)}


class CurrentCoreSourceOverlay:
    def __init__(self,before,core):
        self.before=before;self.core=core
        self.registry={'core':dict(provider=PREFIX+'core_physical_field.CompliantCorePhysicalField',
            method='profiles',coverage_coordinate='rho=R/epsilon_core',domain='[0,4]',acceptance_receipt=RECEIPT),
            **before.registry}

    def provider(self,chart):return self.core if chart=='core' else self.before.provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart!='core':return self.before.evaluate(chart,Z,coordinate)
        c=self.core.ctx;rho=c.mpf(coordinate);z=c.mpf(Z)
        if endpoints(rho)[0]<0 or endpoints(rho)[1]>4 or endpoints(z)[0]<-1 or endpoints(z)[1]>1:
            raise ValueError('Original core chart rho[0,4], Z[-1,1] required')
        packet=self.core.profiles(rho,z)
        return dict(chart=chart,source_packet=packet,
            actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,
            source_provider=self.registry['core']['provider'],acceptance_receipt=RECEIPT,
            physical_mixed_grids={'ordinary_mixed_profile_grids':packet['ordinary_mixed_profile_grids']})


class _CurrentCorePhysicalDispatch:
    def __init__(self,overlay,native):self.current=overlay;self.native=native
    def __getattr__(self,name):return getattr(self.native,name)
    def provider(self,chart):return self.current.provider(chart)
    def evaluate(self,chart,Z,coordinate):return self.current.evaluate(chart,Z,coordinate)


def original_core_operator_bindings():
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantGlobalPhysicalAssembly')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    wanted=ast.dump(ast.parse('self.core.physical_map(packet,axis=axis)',mode='eval').body)
    if sum(isinstance(n,ast.Call) and ast.dump(n)==wanted for n in ast.walk(fn))!=1:
        raise ValueError('Original nonsingular core map call changed')
    assignments={ast.unparse(t):ast.unparse(n.value) for n in ast.walk(fn)
        if isinstance(n,ast.Assign) for t in n.targets if ast.unparse(t) in ('loglambda_bound','bounds[index][component][label]')}
    target='c.ln(c.mpf(norm))+scale[\'logLambda_term\']+scale[\'amplitude_log_upper\']+scale[\'physical_lambda_exponent\']*logtau/2'
    wanted=ast.dump(ast.parse(target,mode='eval').body)
    if sum(isinstance(n,ast.BinOp) and ast.dump(n)==wanted for n in ast.walk(fn))!=1:
        raise ValueError('Original requested-time core prefactor changed')
    return dict(actual_native_core_physical_map_call_bound=True,
        actual_requested_logtau_scale_expression_bound=True,
        native_original_profiles_and_nonsingular_map_retained=True,
        axis_same_core_coordinate_zero_only=True,passed=True)


class CurrentCorePhysicalAssembly(PRIOR):
    @source_precision
    def __init__(self,require_checked=True):
        # One already checked current29 construction, not another core.
        PRIOR.__init__(self,require_checked=True)
        self.core_source=CurrentCoreSourceOverlay(self.source_assembly,self.core)
        self.dispatch=_CurrentCorePhysicalDispatch(self.core_source,self.dispatch.native)
        self.source_owners=dict(self.core_source.registry)
        self.operator_bindings['same_current_provider_and_evaluator']=self.dispatch.current is self.core_source
        self.core_bindings={target:class_assignment('current_core_physical_assembly','CurrentCoreSourceOverlay',
            '__init__',target,expression) for target,expression in {'self.before':'before','self.core':'core'}.items()}
        self.core_operator_proof=original_core_operator_bindings()
        theorem=accepted(THEOREM,self.family,self.source,'whole_core_and_axis_cartesian_spatial4_enclosures_available')
        self.core_functional_proof=functional_recovery_proof()
        if theorem['functional_recovery_proof']!=self.core_functional_proof:
            raise ValueError('Original core analytic recovery/divergence theorem differs')
        if not theorem['independent_nonsingular_core_fixture']['passed']:
            raise ValueError('Retained independent nonsingular core fixture omitted')
        previous=json.loads((HERE/PRIOR_REPORT).read_bytes())
        old=accepted(PRIOR_RECEIPT,self.family,self.source,'current_twenty_nine_downstream_physical_source_ownership_certified')
        if (tuple(previous['current_physical_chart_owners'])!=PRIOR_CHARTS
                or previous['current_source_owner_registry']!={k:self.source_owners[k] for k in PRIOR_CHARTS}
                or previous['datum_enclosure_sha256']!=self.datum_sha):
            raise ValueError('Retained29 source registry/datum differs')
        for record in (theorem,old):
            for name,digest in record['input_hashes'].items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current core source dependency conflict: '+name)
                self.hashes[name]=digest
        for name in (THEOREM,PRIOR_REPORT,PRIOR_RECEIPT):self.hashes[name]=sha(name)
        self.retained_heat_manifest=previous;self.retained_heat_receipt=old;self.core_theorem=theorem
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        failed=[k for k,v in self.current_provider_graph().items() if not v]
        if failed:raise ValueError('Current nested core/axis object graph differs: '+', '.join(failed))
        self.core_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current core receipt datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.core_acceptance_loaded=True

    def current_provider_graph(self):
        graph=PRIOR.current_provider_graph(self)
        if not hasattr(self,'core_source'):return graph
        graph.update(exact_nested_current_core_provider=self.dispatch.provider('core') is self.dispatch.anchor.patch.core is self.core,
            same_current_core_overlay=self.dispatch.current is self.core_source,
            overlay_keeps_checked29_source_object=self.core_source.before is self.source_assembly,
            original_core_class=type(self.core) is CoreContext and type(self.core.original) is CORE,
            original_core_context_wrapper_is_same_current_source=self.core.original is (
                self.dispatch.anchor.patch.reference_mixed.long_mixed.history.upstream.core),
            original_core_profiles=self.core.profiles.__func__ is CORE.profiles,
            original_core_nonsingular_physical_map=self.core.physical_map.__func__ is CORE.physical_map,
            same_core_family_source_datum=(self.core.family,self.core.source,self.core.datum.datum_sha)==(
                self.family,self.source,self.datum_sha),
            original_core_context_is_mapper_context=self.core.ctx is self.ctx,
            retained_current_native_parameter_source_bridge=self.native_parameter_bridge['passed'])
        return graph

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        if chart not in CHARTS or (axis and (chart!='core' or endpoints(self.ctx.mpf(coordinate))!=(mp.mpf(0),mp.mpf(0)))):
            raise ValueError('Current30 source charts; axis requires core coordinate0')
        if chart!='core':
            packet=PRIOR.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)
        else:
            packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)
            packet.update(datum_enclosure_sha256=self.datum_sha,current_source_owner=self.source_owners['core']['provider'],
                current_source_acceptance_receipt=RECEIPT,current_source_dispatcher_used=True,
                current_core_parameter_source_used=True,**dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),
                **{UNIFORM:False},full_pulse_C4_installed=False,core_inner_annulus_interfaces_certified=False)
        packet.update(current_core_physical_acceptance_receipt=RECEIPT,
            current_core_axis_source_functional_binding_proved=True,
            **{GATE:self.core_acceptance_loaded,COMPOSED:self.core_acceptance_loaded,SOURCE_GATE:self.core_acceptance_loaded},
            current_core_axis_physical_owner_installed=self.core_acceptance_loaded,
            current_twenty_nine_downstream_physical_source_ownership_certified=True)
        return packet

    def manifest(self):
        old=self.retained_heat_receipt
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            current_source_owner_registry=self.source_owners,current_provider_graph_identity=self.current_provider_graph(),
            current_core_constructor_bindings=self.core_bindings,original_core_operator_bindings=self.core_operator_proof,
            retained_canonical_core_functional_proof=self.core_functional_proof,
            reused_independent_nonsingular_core_fixture=self.core_theorem['independent_nonsingular_core_fixture'],
            original_core_axis_is_mode_not_extra_owner=True,
            retained_twenty_nine_physical_evidence=dict(report=PRIOR_REPORT,receipt=PRIOR_RECEIPT,
                regular_contributions=old['twenty_nine_regular_source_contributions_checked'],
                supplemental_gap_contributions=old['retained_twenty_seven_chart_physical_evidence']['supplemental_gap_source_contributions'],
                same_registry_sources_datum_and_checked29=True),
            **{GATE:self.core_acceptance_loaded,COMPOSED:self.core_acceptance_loaded,SOURCE_GATE:self.core_acceptance_loaded},
            current_core_axis_physical_owner_installed=self.core_acceptance_loaded,
            current_twenty_nine_downstream_physical_source_ownership_certified=True,
            core_inner_annulus_interfaces_certified=False,full_pulse_C4_installed=False,**{UNIFORM:False},
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();views={}
        for name,(rho,axis) in VIEWS.items():
            views[name]=dict(source=self.core_source.evaluate('core',[-1,1],rho),
                physical=self.evaluate('core',[-1,1],rho,axis=axis))
            print('Current nested core mapped: '+name,flush=True)
        result.update(whole_current_core_physical_maps=views,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentCorePhysicalAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current30 source owners generated: same analytic core and nonsingular axis',flush=True)
    return result


if __name__=='__main__':run()
