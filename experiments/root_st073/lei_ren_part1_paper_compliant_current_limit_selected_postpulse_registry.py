"""One repaired-limit selected branch through the complete postpulse source chain.

Inject the accepted selected owner into the original postpulse and full exterior
algorithms. No second exact repair is constructed. The new registry records the
actual owner path and its own admission, including the copied pulse routes.
This is a similarity-source installation, not a physical/time assembler.
"""
import copy
import json
import time
from pathlib import Path

import lei_ren_part1_paper_compliant_current_limit_selected_pulse_flatten as selected
import lei_ren_part1_paper_compliant_current_postpulse_energy_history as history
import lei_ren_part1_paper_compliant_current_full_exterior_stress as exterior
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import ROUTES
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE,PREFIX,sha,require=selected.HERE,selected.PREFIX,selected.sha,selected.require
NAME=PREFIX+'current_limit_selected_postpulse_registry.json'
RECEIPT=PREFIX+'current_limit_selected_postpulse_registry_check.json'
GATES=('current_repaired_limit_selected_postpulse_source_registry_installed',
       'current_repaired_limit_selected_full_exterior_five_histories_installed',
       'current_repaired_limit_selected_full_exterior_similarity_stress_mixed4_installed')
# The history-only receipt leaves the exterior theorem open; the checked full
# exterior owner supplies it here. Do not publish that same component gate as
# false alongside its accepted, narrower exterior packet.
OPEN=tuple(key for key in dict.fromkeys(selected.OPEN+history.OPEN+exterior.OPEN)
    if key not in exterior.GATES)
ALIASES={'steep_entry':'steep_in','steep_exit':'steep_out'}
DOMAINS={'flatten':(0,100),'outer_power':(0,1),'outer_angular':(-4,0),
    'steep_entry':(0,1),'steep_power':(0,1),'steep_exit':(0,1),
    'waiting':(0,1),'heat_collar':(0,3),'heat_exterior':(3,None)}


def serialized(value):
    return encode(pack(value))


def transport_proofs(current_history,current_exterior,old_history,old_exterior):
    """Allow only the previously checked native parameter graph substitution."""
    energy=serialized(current_history.proof)
    require(energy==old_history['current_energy_transport_source_proof'],
        'Actual full postpulse energy/meridional source proof changed')
    current=copy.deepcopy(serialized(current_exterior.proof))
    original=copy.deepcopy(old_exterior['current_full_terminal_history_transfer'])
    new_parameter=current['current_original_pressure_terminal_source_proof'].pop(
        'common_exact_parameter_function_bridge')
    old_parameter=original['current_original_pressure_terminal_source_proof'].pop(
        'common_exact_parameter_function_bridge')
    require(current==original,'Full exterior history proof changed beyond exact parameter transport')
    require(new_parameter['passed'] and new_parameter['current_core_compatibility_not_asserted'],
        'Checked native-only exact parameter bridge required')
    require(all(new_parameter['actual_native_only_object_graph'].values()),
        'Actual native parameter object identity lost')
    require(new_parameter['actual_mu_delta_defining_assignments']==
        old_parameter['actual_mu_delta_defining_assignments'],
        'Original exact mu/delta source equations changed')
    pressure=current_history.selected.pressure
    require(pressure.transport['passed'] and
        pressure.transport['only_changed_proof_component']=='common_exact_parameter_function_bridge',
        'Previously admitted pressure parameter transport required')
    return dict(unchanged_full_postpulse_energy_source_proof=True,
        unchanged_full_five_history_and_native_stress_source_proof=True,
        only_changed_proof_path='current_original_pressure_terminal_source_proof.common_exact_parameter_function_bridge',
        actual_exact_mu_delta_defining_assignments=new_parameter['actual_mu_delta_defining_assignments'],
        actual_native_parameter_graph=new_parameter['actual_native_only_object_graph'],
        accepted_current_pressure_parameter_transport=pressure.transport,
        old_common_core_graph_not_imported=True,passed=True)


class CurrentLimitSelectedPostpulseRegistry:
    @source_precision
    def __init__(self,selected_owner=None,require_checked=True):
        self.selected_owner=selected_owner if selected_owner is not None else selected.CurrentLimitSelectedPulseFlatten()
        require(self.selected_owner.acceptance_loaded,'Checked repaired-limit selected pulse required')
        self.family_record=self.selected_owner.family_record
        self.family=self.selected_owner.family;self.source=self.selected_owner.source;self.datum_sha=self.selected_owner.datum_sha
        self.hashes=dict(self.selected_owner.hashes)
        selected.source.inlet.checked('current_limit_selected_pulse_flatten',selected.GATE,self.family_record,self.hashes)
        old_history=selected.source.inlet.checked('current_postpulse_energy_history',history.GATES[0],self.family_record,self.hashes)
        old_exterior=selected.source.inlet.checked('current_full_exterior_stress',exterior.GATES[0],self.family_record,self.hashes)
        require(all(old_history[key] for key in history.GATES) and all(old_exterior[key] for key in exterior.GATES),
            'All original full-history and exterior source gates required')
        self.history=history.CurrentPostpulseEnergyHistory(selected=self.selected_owner.selected)
        self.exterior=exterior.CurrentFullExteriorStress(history=self.history)
        self.ctx=self.history.ctx
        self.hashes.update(self.exterior.hashes)
        self.proof_transport=transport_proofs(self.history,self.exterior,old_history,old_exterior)
        self.registry=self.route_registry()
        self.graph=self.object_graph();require(all(self.graph.values()),'Repaired selected postpulse object graph differs')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=selected.source.inlet.checked('current_limit_selected_postpulse_registry',GATES[0],self.family_record,self.hashes)
            require(all(receipt[key] for key in GATES) and not any(receipt[key] for key in OPEN),
                'Repaired selected postpulse admission scope differs')
            self.acceptance_loaded=True

    def object_graph(self):
        owner=self.selected_owner;sel=owner.selected;h=self.history
        return dict(same_accepted_selected_owner=h.selected is sel,
            same_family_source_and_independent_datum=all((obj.family,obj.source,obj.datum_sha)==
                (self.family,self.source,self.datum_sha) for obj in (owner,sel,h,self.exterior,sel.pressure)),
            same_selected_pulse_in_all_postpulse_sources=h.flatten.pulse is h.outer.pulse is sel.pulse,
            same_selected_flatten_downstream=h.outer.flatten is h.flatten,
            original_selected_flatten_preserved=owner.flatten.pulse is sel.pulse and h.flatten is not owner.flatten,
            same_selected_fifth_callback=h.outer.fifth is sel.fifth,
            same_current_future_in_all_postpulse_sources=h.outer.future is h.steep.future is h.heat.future is sel.future,
            same_unique_exact_repair=h.heat.repair is sel.exact.repair is owner.source_field.exact.repair,
            same_exact_Gamma_source=h.heat.exact_heat is sel.future.heat,
            same_pressure_and_exact_source=self.exterior.pressure is sel.pressure is owner.source_field.pressure and sel.exact is owner.source_field.exact,
            exterior_uses_actual_history=self.exterior.history is h and self.exterior.heat is h.heat,
            same_independent_analytic_P0=h.flatten.inlet.datum is sel.future.angular.initial.datum is sel.exact.repair.angular.initial.datum,
            complete_postpulse_source_graph=all(h.assert_graph().values()),
            actual_full_history_source_transfer=self.exterior.proof['passed'],
            source_proof_transport=self.proof_transport['passed'])

    def route_registry(self):
        registry={}
        owner_path=PREFIX+'current_limit_selected_postpulse_registry.CurrentLimitSelectedPostpulseRegistry.'
        for chart,(method,domain) in selected.source.inlet.EXPECTED.items():
            route=ROUTES[chart]
            require(route[:3]==('pulse_mixed_C4','CompliantPulseMixedC4',method) and route[4]==domain,
                'Original pulse route definition changed: '+chart)
            registry[chart]=dict(provider=owner_path+'selected_owner.pulse',method=method,
                coverage_coordinate=route[3],domain=domain,
                acceptance_receipt=RECEIPT,selected_source_receipt=selected.RECEIPT)
        for chart,(lo,hi) in DOMAINS.items():
            native=ALIASES.get(chart,chart);owner,method=history.METHODS[native]
            if chart=='flatten':coordinate='t=log(R/Rv)'
            elif chart=='outer_power':coordinate='phase=(logR-logRf)/(Lrel-4)'
            elif chart=='outer_angular':coordinate='s=logR-logRrel'
            elif chart in ('steep_entry','steep_exit'):coordinate='t=logR-logRorigin'
            elif chart in ('steep_power','waiting'):coordinate='phase=offset/original_length'
            else:coordinate='t=log(R/Rtail)'
            registry[chart]=dict(provider=owner_path+('exterior' if chart=='heat_exterior' else 'history.'+owner),
                method='exterior' if chart=='heat_exterior' else method,coverage_coordinate=coordinate,
                domain='[3,infinity)' if hi is None else '[%s,%s]'%(lo,hi),
                acceptance_receipt=RECEIPT,selected_source_receipt=selected.RECEIPT)
        for route in registry.values():
            route.update(self.family_record,input_hash_ledger=RECEIPT+'#input_hashes')
        return registry

    def provider(self,chart):
        if chart in selected.source.inlet.PULSE_CHARTS:return self.selected_owner.pulse
        if chart=='heat_exterior':return self.exterior
        if chart not in DOMAINS:raise ValueError('Unknown repaired selected postpulse route')
        return getattr(self.history,history.METHODS[ALIASES.get(chart,chart)][0])

    def underlying_source(self,chart):
        owner=self.provider(chart);method=getattr(owner,self.registry[chart]['method']).__func__
        name=method.__module__+'.py';digest=sha(name)
        require(self.hashes.get(name)==digest,'Actual underlying source method is not in the checked hash ledger')
        return dict(owner_class=type(owner).__module__+'.'+type(owner).__qualname__,
            actual_method=method.__module__+'.'+method.__qualname__,
            actual_method_source_file=name,actual_method_source_sha256=digest,**self.family_record)

    def guard(self,chart,Z,coordinate):
        c=self.ctx;lo,hi=history.endpoints(c.mpf(Z))
        if lo<-1 or hi>1:raise ValueError('Original whole-Z source domain [-1,1] required')
        if chart in selected.source.inlet.PULSE_CHARTS:return
        if chart not in DOMAINS:raise ValueError('Unknown repaired selected postpulse route')
        lower,upper=DOMAINS[chart];lo,hi=history.endpoints(c.mpf(coordinate))
        if lo<lower or (upper is not None and hi>upper):
            raise ValueError('Original route domain required: '+chart+' '+self.registry[chart]['domain'])

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.guard(chart,Z,coordinate);require(all(self.object_graph().values()),'Current source graph changed')
        route=self.registry[chart]
        if chart in selected.source.inlet.PULSE_CHARTS:
            out=dict(self.selected_owner.evaluate(chart,Z,coordinate))
        else:
            point=(self.exterior.exterior(Z,coordinate) if chart=='heat_exterior' else
                self.history.evaluate(ALIASES.get(chart,chart),Z,coordinate)['source_packet'])
            out=dict(chart=chart,source_packet=point)
        out.update(source_provider=route['provider'],acceptance_receipt=RECEIPT,
            underlying_source_provider=self.underlying_source(chart),
            selected_source_receipt=selected.RECEIPT,source_coverage_coordinate=route['coverage_coordinate'],
            source_coordinate_domain=route['domain'],source_family=self.family_record,
            original_analytic_P0_Taylor_retained=self.history.flatten.inlet.incoming(Z)['P0'],
            same_repaired_limit_selected_future_in_all_postpulse_routes=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False),
            output_kind='current repaired selected similarity-source enclosures; no physical point selection')
        return out


@source_precision
def run(field=None):
    began=time.monotonic();field=field if field is not None else CurrentLimitSelectedPostpulseRegistry(require_checked=False)
    result=dict(source_family=field.family_record,**dict.fromkeys(GATES,False),**dict.fromkeys(OPEN,False),
        candidate_selected_postpulse_registry_constructed=True,
        actual_same_object_graph=field.object_graph(),ordered_current_source_registry=field.registry,
        unchanged_energy_source_proof=field.history.proof,
        current_full_terminal_history_transfer=field.exterior.proof,
        strict_original_source_proof_transport=field.proof_transport,
        no_second_exact_repair_constructed=True,
        original_selected_and_native_graphs_preserved=True,
        scope='selected pulse through full Gamma exterior similarity sources on the repaired-limit branch; global physical/absolute post-2Rc closure separate')
    result['actual_source_handshakes']={}
    for chart,Z,coordinate in (('pulse_end','.427',0),('flatten','.427',0),
        ('outer_angular','.427',0),('waiting','.427',1),('heat_collar','.427',0),
        ('heat_exterior','.427','4.23')):
        result['actual_source_handshakes'][chart]=serialized(field.evaluate(chart,Z,coordinate))
        print('Current repaired selected source: '+chart,flush=True)
    result['all_passed']=True;result['input_hashes']=field.hashes
    result['execution_seconds']=time.monotonic()-began
    (HERE/NAME).write_text(json.dumps(serialized(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
