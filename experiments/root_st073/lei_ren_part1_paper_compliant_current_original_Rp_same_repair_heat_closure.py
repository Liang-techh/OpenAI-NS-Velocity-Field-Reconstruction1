"""Absolute heat/pressure closure on the existing current selected repair.

The exact input view installs current pulse/flatten/future/angle providers
without replaying the repair. Original function identities eliminate Dtheta
and Cp; their directed forward subtraction enclosures remain diagnostic.
"""
import json
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_raw_history_transport as raw_source
import lei_ren_part1_paper_compliant_current_limit_heat_pressure_bridge as closure
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=raw_source.HERE,raw_source.PREFIX,raw_source.sha
NAME=PREFIX+'current_original_Rp_same_repair_heat_closure.json'
RECEIPT=PREFIX+'current_original_Rp_same_repair_heat_closure_check.json'
GATES=('current_original_Rp_existing_selected_repair_heat_pressure_closure_runtime_installed',
       'current_original_Rp_heat_Dtheta_Cp_source_function_identities_installed',
       'current_original_Rp_closed_heat_factorized_histories_callable')
OPEN=raw_source.OPEN
VIEWS=(('collar_inlet','heat_collar',0),('collar_active','heat_collar',1),
    ('collar_exit','heat_collar',3),('exterior_inlet','heat_exterior',3),
    ('exterior_fresh','heat_exterior',4))


def absolute_source_binding(before,exact):
    """Bind the current graph to the existing repair's defining functions.

    Only the selected logC is an exact dyadic source datum. Directed radius,
    waiting and amplitude boxes are never substituted for function equality.
    """
    radius=before.radius;native=before.post.before.inlet.exact
    prior=exact.prior;repair=exact.repair;g=radius.graph
    require=closure.require
    source_name=PREFIX+'current_original_O2_source_parameter_frame.json'
    source_check=PREFIX+'current_original_O2_source_parameter_frame_check.json'
    source=json.loads((HERE/source_name).read_bytes());receipt=json.loads((HERE/source_check).read_bytes())
    require(receipt['all_passed'] and receipt['input_hashes'][source_name]==sha(source_name),
        'Accepted exact selected logC source required')
    for name,digest in receipt['input_hashes'].items():
        require(sha(name)==digest,'Changed selected parameter source '+name)
    norm=repair.angular.initial.repair.records['compliant_physical_norm_family']['selected_logCstar']
    chosen=tuple(source['selected_logCstar_exact_mpf_tuple'])
    require(chosen==tuple(norm['lower_exact_mpf_tuple'])==tuple(norm['upper_exact_mpf_tuple'])
        and repair.heat.logC._mpi_==(chosen,chosen),'Repair/current exact selected logC source differs')
    leaf=g.nodes[radius.parameters['logC'].node]
    require(leaf['operation']=='current_original_source_parameter' and leaf['name']=='logC'
        and leaf['definition']=='same selected_logCstar and CurrentLongRadiusPhase exact singleton tuple'
        and leaf['source_family']==before.family_record and source['source_family']==before.family_record
        and leaf['quantity_is_exact_original_parameter_not_range_endpoint']
        and before.hashes[leaf['current_source_receipt']]==leaf['current_source_receipt_sha256'],
        'Absolute logC graph provenance differs')
    ast=closure.angular.class_assignment
    assignments=dict(
        logC=ast('exact_heat_component','SharedExactHeatComponent','__init__','self.logC',
            "read_interval(c,norms['selected_logCstar'])"),
        reference=ast('exact_heat_component','SharedExactHeatComponent','__init__','self.logRref',
            'c.ln(110)+10*(self.logC+self.params.logPstar)'),
        flattenU=ast('current_pulse_flatten_source','CurrentFlattenMixedC4','__init__','self.U',
            "self.inlet.constants['U']"))
    closure.repair.binding('compliant_current_exact_repair_branch','replay_heat','out.tail_finite',
        'p.yd+1+p.Tw+100-30*p.log_mu+2+p.Ts+angular.waiting')
    closure.repair.binding('compliant_current_exact_repair_branch','replay_heat','out.logradius_terms',
        'dict(selected_reference=out.logRref,pulse_term=13/p.mu,finite_offset=out.tail_finite)')
    closure.repair.binding('compliant_current_exact_repair_branch','__init__','angular.loguRp0',
        'c.ln(c.mpf(endpoints(self.flatten.U)))')
    # The prior repair explicitly consumes the old native constructor U.
    # Its admitted current/native function theorem connects that U to U0.
    require(prior.flatten.U is prior.flatten.inlet.constants['U'] is prior.pulse.high.constants['U']
        and prior.graph['current_heat_amplitude_origin_rebound'],'Prior native U source aliases differ')
    require(exact.flatten.U is exact.flatten.inlet.constants['U'] is exact.pulse.high.constants['U'],
        'Current native U source aliases differ')
    native_receipt=json.loads((HERE/raw_source.radius.post.selected.inlet.exact.RECEIPT).read_bytes())
    require(native.acceptance_loaded and native_receipt['all_passed']
        and native_receipt['current_to_original_native_constant_function_identities']==12
        and native_receipt['actual_native_constructor_to_current_constant_function_identities']==12
        and native_receipt['actual_source_functions_identified_not_interval_enclosure_equality'],
        'Accepted current/original native U defining function theorem required')
    q0,constants=native.current_constants();old,kernels=native.original_native_functions(q0)
    field=SimpleNamespace(graph=before.graph,parameters=radius.parameters,
        contracts=radius.frame.bridge.leading.contracts)
    q=raw_source.radius.RadiusInterpreter(field,False,{radius.frame.functions['Tw'].node:s.Symbol('Tw',positive=True)})
    require(s.cancel(q.at(before.U0)-constants['U'])==0 and s.cancel(constants['U']-old['U'])==0,
        'Actual U0 graph/current/native defining functions differ')
    # Expand the actual absolute origin. Other source parameters remain
    # symbols, so cancellation of the microscopic origin is explicit.
    bindings={ref.node:s.Symbol('source_'+name,positive=True) for name,ref in radius.parameters.items()}
    lp=s.exp(40)+11;bindings[radius.parameters['logP'].node]=lp
    mu=bindings[radius.parameters['mu'].node];lc=bindings[radius.parameters['logC'].node]
    W=s.Symbol('same_original_uncorrected_waiting',real=True)
    bindings[radius.functions['waiting'].node]=W
    q=raw_source.radius.RadiusInterpreter(field,False,bindings)
    ref=s.log(110)+10*(lc+lp);Tw=-60*s.log(mu);L=-30*s.log(mu)
    Ts=q.at(radius.functions['Ts']);v=q.at(radius.native)
    desired_Rp=ref+lp+1+Tw
    require(s.simplify(s.expand_log(q.at(radius.logRp)-desired_Rp,force=True))==0,
        'Actual absolute logRp does not equal original selected reference')
    tail=ref+13/mu+lp+1+Tw+100+L+2+Ts+W
    for chart in ('heat_collar','heat_exterior'):
        require(s.simplify(s.expand_log(q.at(radius.maps[chart]['logR'])-tail-v,force=True))==0,
            'Current heat radius differs from existing repair radius '+chart)
    # The already accepted native constant receipt binds yd=logPstar and
    # Tw=-60 log(mu), at the same Md=40 parameter family.
    require(native_receipt['actual_native_constructor_source_binding']['actual_y_d_and_Tw_parameter_source']
        ['same_actual_graph_logP_mu_Tw_definitions'],'Original yd/logP/mu/Tw definitions differ')
    result=dict(passed=True,same_exact_selected_logC_source=True,
        actual_absolute_Rp_matches_existing_repair_reference=True,
        both_heat_radius_functions_match_existing_repair=True,
        same_original_uncorrected_waiting_function_receipt=raw_source.radius.RECEIPT,
        actual_current_U0_equals_original_native_U_function=True,
        current_and_prior_flatten_U_aliases_bound=True,
        repair_loguRp0_source_assignment_is_log_original_native_U=True,
        source_assignments=assignments,
        heat_radius_identity='logR=logRref+13/mu+yd+1+Tw+100-30*logmu+2+Ts+W+t',
        amplitude_identity='U0=original_native_U; loguRp0=log(original_native_U)',
        no_directed_box_equality_used_as_radius_waiting_or_U_function_proof=True,
        input_hashes={source_name:sha(source_name),source_check:sha(source_check)})
    return result


class ExistingSelectedExactRepairView:
    """Checked source-input view, with no repair or future constructor call."""
    def __init__(self,post):
        post.assert_graph();b=post.before;prior=b.seed.exact
        if not prior.acceptance_loaded:raise ValueError('Accepted existing exact repair required')
        self.post=post;self.prior=prior;self.companion=post.companion
        self.heat=post.heat;self.pulse=b.pulse;self.flatten=b.flatten
        self.repair=prior.repair;self.future=b.future;self.angular4=b.angular4
        # Only the active current angular method is exported to closure.
        self.angular5=SimpleNamespace(ctx=b.ctx,angular=b.fifth.angular)
        self.ctx=b.ctx;self.family=b.family;self.source=b.source;self.datum_sha=b.datum_sha
        self.parameter_bridge=post.power.bindings['current_native_parameter_source_bridge']
        self.theorem=prior.theorem;self.uniform=prior.uniform
        self.hashes=dict(post.hashes)
        raw_source.radius.post.selected.inlet.add_hashes(self.hashes,prior.hashes)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.graph=self.assert_graph();self.acceptance_loaded=True

    def assert_graph(self):
        b=self.post.before;r=self.repair
        graph=dict(same_existing_accepted_branch=self.prior.acceptance_loaded and r is self.prior.repair,
            same_current_selected_repair=r is self.future.repair is self.angular4.repair is self.post.heat.repair,
            same_current_future=self.future is b.future is self.post.steep.future is self.post.heat.future,
            same_current_heat_and_companion=self.heat is self.companion.heat is self.post.heat,
            current_pulse_and_flatten=self.pulse is b.pulse and self.flatten is b.flatten and self.flatten.pulse is self.pulse,
            current_angular_C5_callable=self.angular5.angular.__self__ is b.fifth
                and b.fifth.angular4 is self.angular4
                and self.angular5.angular.__func__ is closure.CompliantFifthAxialJets.angular,
            exact_heat_angular_aliases=r.angular is r.heat.angular is self.future.angular
                and r.heat is self.future.heat is self.post.heat.exact_heat,
            actual_parameter_aliases=self.future.params is r.params is r.heat.params is r.angular.params
                is self.post.steep.future.params,
            original_parameter_function_bridge=self.parameter_bridge['passed']
                and self.parameter_bridge['current_core_compatibility_not_asserted'],
            same_independent_P0_object=self.flatten.inlet.datum is b.datum is r.angular.initial.datum,
            same_context=self.ctx is r.ctx is self.future.ctx is self.post.ctx,
            retained_unique_branch_theorem=self.theorem['whole_Z_contraction_and_derivative_inverse_checked'])
        if not all(graph.values()):raise ValueError('Existing selected exact-repair view differs: '+str(graph))
        return graph


class CurrentOriginalRpSameRepairHeatClosure:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.input_raw=before if before is not None else raw_source.CurrentOriginalRpRawHistoryTransport()
        if type(self.input_raw) is not raw_source.CurrentOriginalRpRawHistoryTransport or not self.input_raw.acceptance_loaded:
            raise ValueError('Accepted typed current raw-history transport required')
        # A cheap expression-only copy avoids changing the upstream graph and
        # gives deterministic node IDs irrespective of previous runtime calls.
        self.before=raw_source.CurrentOriginalRpRawHistoryTransport(self.input_raw.before)
        self.before.assert_graph();self.post=self.before.post
        self.family_record=self.post.family_record;self.ctx=self.post.ctx
        self.exact=ExistingSelectedExactRepairView(self.post)
        self.absolute_binding=absolute_source_binding(self.before,self.exact)
        self.collar=closure.collar.CurrentCollarStressMixedC4(companion=self.post.companion)
        self.angular=closure.angular.CurrentAngularTerminalClosure(exact=self.exact,collar=self.collar)
        self.balance=closure.balance.CurrentPressureTerminalBalance(angular=self.angular)
        self.raw=closure.raw.CurrentRawPreheatPressureOperator(angular=self.angular)
        self.pressure=closure.CurrentLimitPressureClosure(balance=self.balance,raw=self.raw,family_record=self.family_record)
        self.hashes=dict(self.before.hashes)
        raw_source.radius.post.selected.inlet.add_hashes(self.hashes,self.pressure.hashes)
        raw_source.radius.post.selected.inlet.add_hashes(self.hashes,self.absolute_binding['input_hashes'])
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.call_trace=[];self.acceptance_loaded=False
        self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[k] for k in GATES) or any(record[k] for k in OPEN):
                raise ValueError('Current same-repair heat closure receipt or scope differs')
            if record['source_family']!=self.family_record:raise ValueError('Same-repair closure family differs')
            if record['actual_absolute_source_binding']!=raw_source.packed(self.absolute_binding):
                raise ValueError('Accepted absolute radius/amplitude source binding differs')
            raw_source.radius.post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        e=self.exact;p=self.post;a=self.angular
        graph=dict(accepted_raw_geometry_source=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            independent_geometry_same_live_providers=self.before is not self.input_raw
                and self.before.graph is not self.input_raw.graph and self.before.post is self.input_raw.post,
            exact_view_bound_to_active_selected_owner=all(e.assert_graph().values()),
            same_collar_companion=self.collar.companion is e.companion is p.companion,
            angular_consumes_same_exact=a.exact is e and a.collar is self.collar,
            balance_raw_pressure_share_one_exact=self.balance.exact is self.raw.exact is self.pressure.exact is e,
            balance_raw_pressure_share_one_angular=self.balance.angular is self.raw.angular is self.pressure.angular is a,
            pressure_consumes_actual_balance_and_raw=self.pressure.balance is self.balance and self.pressure.raw is self.raw,
            unchanged_P0_object=self.raw.flat.inlet.datum is p.before.datum,
            original_function_proofs=self.angular.proof['passed'] and self.balance.proof['passed']
                and self.raw.proof['passed'] and self.pressure.proof['passed'],
            native_only_parameter_transport=self.pressure.transport['passed'],
            independent_closure_caches=len({id(v) for v in
                (self.collar.cache,a.cache,a.terminals,self.balance.cache,self.raw.cache,self.pressure.cache,
                 a.heat.shape_cache,a.heat.tail_cache,a.heat.gamma_cache)})==9,
            no_repair_or_future_replay=e.repair is p.before.seed.exact.repair and e.future is p.before.future)
        graph['same_absolute_radius_and_amplitude_source_functions']=self.absolute_binding['passed']
        if not all(graph.values()):raise ValueError('Same-repair current heat closure graph differs: '+str(graph))
        return graph

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.assert_graph()
        if chart not in ('heat_collar','heat_exterior'):raise ValueError('Actual current heat chart required')
        geometry=self.before.radius.geometry(chart,coordinate)
        packet=self.pressure.evaluate(chart,Z,coordinate)
        c=self.ctx;z=c.mpf(Z);value=raw_source.radius.exact_coordinate(coordinate)
        t=c.mpf(value.numerator)/value.denominator
        forward=self.angular.terminal_constants(z)['current_repaired_forward_terminal']
        if chart=='heat_collar':
            shape=closure.angular.shape_radial5(self.angular.heat,z,t)
            tails=self.angular.heat.collar_tails(z,t)
            remaining=tails['remaining_pressure_in_Rtail_units']
            energy=tails['remaining_energy_in_Rtail_units']*c.exp(self.angular.heat.delta*t)/(2*shape['K_rows'][0]**2)
        else:
            shape=self.angular.heat.local_Gamma(z,t)
            remaining=shape['pressure_numerator']*c.exp(-self.angular.heat.prate*t)
            energy=shape['energy_numerator']/(2*shape['K_rows'][0]**2)
        theta=shape['K_rows'][0]*forward['theta_base']*c.exp(-self.angular.heat.bh*t)
        X=packet['actual_angular_Taylor_after_source_closure'];zero=raw_source.IntervalTaylor.constant(c,0,5)
        P0=packet['original_analytic_P0_Taylor_retained'];P=-remaining*forward['pressure_scale']
        jets=dict(Mz=zero,Mtheta=theta*X*c.sqrt(2),Mtheta_z=zero,
            Mztheta=theta*theta*energy,Mp=P-P0,P0=P0,pressure=P,Utheta=theta)
        histories={name:self.before.factor(name,chart,geometry,raw_source.POWERS[name],value)
            for name,value in jets.items()}
        densities=dict(Mz=zero,Mtheta=theta*c.sqrt(2),Mtheta_z=zero,
            Mztheta=-theta*theta/2,Mp=theta*theta/2,P0=zero,pressure=theta*theta/2)
        powers={**raw_source.POWERS,'Mp':(0,2,0),'pressure':(0,2,0)}
        derivatives={name:self.before.factor('d_native_'+name,chart,geometry,powers[name],value)
            for name,value in densities.items()}
        self.call_trace.append(dict(chart=chart,actual_closed_heat_pressure_called=True,
            actual_current_selected_repair_consumed=True,actual_absolute_radius_mapper_called=True,
            actual_closed_raw_heat_histories_called=True))
        return dict(source_packet=packet,current_absolute_source_geometry=geometry,
            raw_histories=histories,first_native_radial_derivatives=derivatives,
            closed_normalized_histories=dict(theta=theta,angular=X,energy=energy,
                remaining_pressure=remaining,closed_pressure=P,original_P0=P0),
            exact_Ev0_log_parts={key:value.node for key,value in self.before.logEv0_parts.items()},
            source_family=self.family_record,forward_constant_enclosures_retained=True,
            no_second_repair_branch_or_numeric_point_choice=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


def view_report(view):
    return {**{key:value for key,value in view.items() if key not in ('raw_histories','first_native_radial_derivatives')},
        'raw_histories':{key:value.report() for key,value in view['raw_histories'].items()},
        'first_native_radial_derivatives':{key:value.report() for key,value in view['first_native_radial_derivatives'].items()}}


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpSameRepairHeatClosure(before,require_checked=False)
    views={};Z='.457'
    for name,chart,coordinate in VIEWS:
        views[name]=view_report(owner.evaluate(chart,Z,coordinate))
        print('Actual current same-repair heat closure:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_same_repair_closure_constructed=True,
        actual_same_selected_source_graph=owner.assert_graph(),actual_exact_input_view_graph=owner.exact.graph,
        current_parameter_source_bridge=owner.exact.parameter_bridge,
        actual_absolute_radius_and_amplitude_source_binding=owner.absolute_binding,
        reused_actual_angular_function_proof=owner.angular.proof,
        reused_actual_pressure_balance_proof=owner.balance.proof,
        reused_actual_raw_pressure_operator_proof=owner.raw.proof,
        reused_original_pressure_function_identification=owner.pressure.proof,
        current_pressure_parameter_transport=owner.pressure.transport,
        fresh_Z=Z,actual_five_closed_heat_views=views,actual_source_call_trace=owner.call_trace,
        unchanged_selected_repair_branch=True,no_repair_root_or_future_replay=True,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(raw_source.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    return owner


if __name__=='__main__':run()
