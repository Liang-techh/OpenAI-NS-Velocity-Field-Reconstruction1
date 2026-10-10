"""Current selected Rp owner through every original postpulse source chart.

The existing native algorithms receive the accepted current pulse/flatten.
No repair is re-solved and no pressure constant or moment history is reset.
Outputs remain source enclosures; physical point evaluation is separate.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_Rp_selected_pulse as selected
import lei_ren_part1_paper_compliant_current_limit_heat_pressure_bridge as downstream
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=selected.HERE,selected.PREFIX,selected.sha
NAME=PREFIX+'current_original_Rp_postpulse_source.json'
RECEIPT=PREFIX+'current_original_Rp_postpulse_source_check.json'
GATES=('current_original_Rp_selected_pulse_installed_in_all_postpulse_owners',
       'current_original_Rp_postpulse_source_registry_installed',
       'current_original_Rp_postpulse_five_source_histories_available')
OPEN=tuple(key for key in dict.fromkeys(selected.PENDING+downstream.OPEN+(
    'current_original_Rp_postpulse_absolute_physical_radius_mapper_installed',
    'current_original_Rp_postpulse_pressure_terminal_constants_eliminated',
    'current_original_Rp_uniform_selected_postpulse_seams_certified')) if key not in GATES)
METHODS={**{key:('outer',value) for key,value in
    {'outer_power':'power','outer_angular':'angular'}.items()},
    **{key:('steep',value) for key,value in downstream.steep.METHODS.items()},
    **{key:('heat',value) for key,value in downstream.heat.METHODS.items()}}
VIEWS=(('flatten_exit','flatten',100),('power_inlet','outer_power',0),
    ('power_exit','outer_power',1),('angular_inlet','outer_angular',-4),
    ('angular_active','outer_angular',-3),('angular_exit','outer_angular',0),
    ('steep_entry_inlet','steep_entry',0),('steep_entry_exit','steep_entry',1),
    ('steep_power_inlet','steep_power',0),('steep_power_exit','steep_power',1),
    ('steep_exit_inlet','steep_exit',0),('steep_exit_exit','steep_exit',1),
    ('waiting_inlet','waiting',0),('waiting_exit','waiting',1),
    ('collar_inlet','heat_collar',0),('collar_exit','heat_collar',3),
    ('exterior_inlet','heat_exterior',3),('exterior_fresh','heat_exterior',4))


class CurrentOriginalRpPostpulseSource:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else selected.CurrentOriginalRpSelectedPulse()
        if type(self.before) is not selected.CurrentOriginalRpSelectedPulse or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current original Rp selected owner required')
        self.before.assert_graph()
        self.family_record=self.before.family_record
        self.family,self.source,self.datum_sha=self.before.family,self.before.source,self.before.datum_sha
        self.ctx=self.before.ctx
        # These constructors build only downstream providers, each with fresh
        # mutable caches. Do not invoke the full bridge's repair replay.
        self.power=downstream.CurrentLimitPowerSource(self.before)
        self.steep_source=downstream.CurrentLimitSteepSource(self.power)
        self.heat_source=downstream.CurrentLimitHeatSource(self.steep_source)
        self.outer=self.power.outer;self.steep=self.steep_source.steep;self.heat=self.heat_source.heat
        self.companion=downstream.CurrentLimitHeatCompanion(self.heat_source)
        self.hashes=dict(self.companion.hashes)
        for name in (selected.NAME,selected.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.registry={key:{**value,'acceptance_receipt':RECEIPT}
            for key,value in self.heat_source.registry.items()}
        self.call_trace=[];self.acceptance_loaded=False
        self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[key] for key in GATES) or any(record[key] for key in OPEN):
                raise ValueError('Current postpulse receipt or scope differs')
            if record['source_family']!=self.family_record:raise ValueError('Current postpulse family differs')
            selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        self.before.assert_graph()
        b=self.before;p=self.outer;st=self.steep;h=self.heat
        graph=dict(current_selected_owner_required=b.acceptance_loaded,
            actual_power_consumes_current_owner=self.power.before is b)
        # Spell out every active consumer rather than transporting the old
        # bridge's object graph or an unrelated accepted repair receipt.
        graph['same_actual_pulse']=p.pulse is b.pulse and p.fifth is b.fifth
        graph.update(same_actual_flatten=p.flatten is b.flatten,
            dispatch_stage_owner_aliases=self.power.outer is p and self.steep_source.steep is st
                and self.heat_source.heat is h,
            same_current_future=p.future is st.future is h.future is b.future,
            same_unique_repair=h.repair is p.future.repair is b.angular4.repair is b.seed.exact.repair,
            same_exact_heat=h.exact_heat is b.future.heat is b.future.repair.heat,
            steep_consumes_same_power=st.outer is p and self.steep_source.before is self.power,
            heat_consumes_same_steep=h.steep is st and self.heat_source.before is self.steep_source,
            companion_consumes_same_heat=self.companion.heat is h and self.companion.source_owner is self.heat_source,
            same_current_context=p.ctx is st.ctx is h.ctx is self.companion.ctx is b.ctx,
            independent_current_P0=p.flatten.inlet.datum is b.datum is b.future.angular.initial.datum,
            original_power=p.power.__func__ is downstream.power.BASE.power,
            original_angular=p.angular.__func__ is downstream.power.BASE.angular,
            original_steep=all(getattr(st,method).__func__ is getattr(downstream.steep.BASE,method)
                for method in downstream.steep.METHODS.values()),
            original_heat=all(getattr(h,method).__func__ is getattr(downstream.heat.BASE,method)
                for method in downstream.heat.METHODS.values()),
            same_complete_Gamma_future=self.power.bindings['current_prefix_and_future_decomposition']['passed']
                and self.heat_source.bindings['shared_exact_Gamma_future_binding']['passed'],
            fresh_downstream_cache_owners=len({id(cache) for cache in
                (p.cache,st.cache,st.kernel_cache,h.cache,h.shape_cache,h.tail_cache,h.gamma_cache,
                 self.companion.constants,self.companion.runtime)})==9,
            all_fifteen_routes_present=set(self.registry)==set(b.chain_routes)|{'flatten'}|set(METHODS))
        if not all(graph.values()):raise ValueError('Current postpulse graph differs: '+str(graph))
        return graph

    def provider(self,chart):
        if chart in METHODS:return getattr(self,METHODS[chart][0])
        return self.before.provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.assert_graph()
        if chart not in self.registry:raise ValueError('Unknown current postpulse source chart: '+chart)
        if chart not in METHODS:
            result=dict(self.before.evaluate(chart,Z,coordinate))
        else:
            c=self.ctx;t=c.mpf(coordinate);lo,hi=selected.inlet.endpoints(t)
            if (chart=='outer_angular' and (lo<-4 or hi>0) or
                chart=='heat_collar' and (lo<0 or hi>3) or
                chart=='heat_exterior' and lo<3 or
                chart not in ('outer_angular','heat_collar','heat_exterior') and (lo<0 or hi>1)):
                raise ValueError('Original source chart domain required: '+chart)
            result=dict(chart=chart,source_packet=getattr(self.provider(chart),METHODS[chart][1])(Z,t))
        self.call_trace.append(dict(chart=chart,method=self.registry[chart]['method'],
            actual_current_provider=True))
        result.update(source_family=self.family_record,acceptance_receipt=RECEIPT,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        return result

    @source_precision
    def histories(self,chart,Z,coordinate,packet=None):
        """Five normalized source histories, with explicit raw scaling recipes.

        Mp is defined by original P-P0. This subtraction supplies a directed
        enclosure only; it is never used to assert terminal cancellation.
        """
        if chart not in METHODS and chart!='flatten':raise ValueError('Postpulse history chart required')
        p=packet if packet is not None else self.evaluate(chart,Z,coordinate)['source_packet']
        terminal=self.before.flatten.terminal_histories(Z)
        P0=terminal['P0_over_Pstar_squared']
        if not p['actual_terminal_zero_linear_and_radial_histories_inherited']:
            raise ValueError('Original compact terminal history theorem omitted')
        zero=IntervalTaylor.constant(self.ctx,0,5)
        pressure=p['pressure'] if chart=='flatten' else dict(P0_over_Pstar_squared=P0,
            P_over_Pstar_squared=p['pressure_over_Pstar_squared_Taylor'],
            Mp_over_Pstar_squared=p['pressure_over_Pstar_squared_Taylor']-P0)
        return dict(Mz_over_R_Utheta=zero,
            Mtheta_z_over_sqrt2_R_3half_Utheta_squared=zero,
            Mtheta_over_sqrt2_R_3half_Utheta=p['angular_Taylor'],
            Mztheta_over_R_Utheta_squared=p['energy_Taylor'],
            Mp_over_Pstar_squared=pressure['Mp_over_Pstar_squared'],
            P0_over_Pstar_squared=P0,P_over_Pstar_squared=pressure['P_over_Pstar_squared'],
            raw_moment_scale_recipes=dict(Mz='R*Utheta',Mtheta_z='sqrt(2)*R^(3/2)*Utheta^2',
                Mtheta='sqrt(2)*R^(3/2)*Utheta',Mztheta='R*Utheta^2',Mp='Pstar^2'),
            current_Rv_terminal_source_receipt=selected.RECEIPT,
            source_derived_compact_terminal_linear_zeros=True,
            absolute_pressure_memory_not_reset=True,
            enclosure_subtraction_not_used_as_terminal_zero_proof=True,
            physical_raw_moment_evaluation_installed=False)

    @source_precision
    def heat_constants(self,Z):
        # Keep the actual angular and pressure infinity constants. An
        # enclosure containing zero is not an identity eliminating them.
        return self.companion.terminal_constants(Z)


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpPostpulseSource(before,require_checked=False)
    views={};Z='.413'
    for name,chart,coordinate in VIEWS:
        view=owner.evaluate(chart,Z,coordinate)
        view['five_source_histories']=owner.histories(chart,Z,coordinate,view['source_packet'])
        views[name]=view
        print('Actual current Rp postpulse:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_postpulse_constructed=True,
        actual_same_object_graph=owner.assert_graph(),ordered_current_source_registry=owner.registry,
        actual_source_call_trace=owner.call_trace,fresh_Z=Z,fresh_actual_source_views=views,
        retained_actual_heat_terminal_constants=owner.heat_constants(Z),
        retained_original_power_functional_identities=owner.power.functional_proof,
        retained_original_steep_functional_identities=owner.steep_source.functional_proof,
        retained_original_heat_functional_identities=owner.heat_source.functional_proof,
        current_complete_future_binding=owner.power.bindings['current_prefix_and_future_decomposition'],
        current_heat_source_binding=owner.heat_source.bindings,
        current_absolute_heat_history_binding=owner.companion.history_bindings,
        no_repair_replay_or_reselection=True,physical_radius_handle_only=True,
        current_physical_logRp_node=owner.before.current_Rp_physical_log_radius.node,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(selected.inlet.encode(selected.inlet.pack(result)),indent=2)+'\n',
        encoding='utf8',newline='\n')
    return owner


if __name__=='__main__':run()
