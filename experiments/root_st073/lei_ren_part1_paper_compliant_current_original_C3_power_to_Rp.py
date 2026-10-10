"""Portable actual C3 repaired exit -> quiet power -> native pulse input frame.

Consume the genuine compact-repair exit at 2Rc, not a homogeneous shortcut
from Rc. Selected native pulse acceptance remains a separate function bridge.
"""
from fractions import Fraction
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C3_repair_band as band
import lei_ren_part1_paper_compliant_current_original_outer_C3_continuation as outer

target,ranges,current,source=outer.target,outer.ranges,outer.current,outer.source
HERE,PREFIX,sha,ep=outer.HERE,outer.PREFIX,outer.sha,outer.ep
NAME=PREFIX+'current_original_C3_power_to_Rp.json.gz'
RECEIPT=PREFIX+'current_original_C3_power_to_Rp_check.json'
GATES=('current_original_actual_repaired_exit_C3_quiet_power_and_pulse_frame_installed',
    'current_original_quiet_y4_Z3_profile_history_functions_installed',
    'current_original_quiet_C3_directed_magnitude_and_positive_swirl_ranges_installed')
OPEN=band.OPEN


def accepted_inputs():
    reports,hashes={},{}
    for module in (band,band.control,outer):
        report,receipt=outer.rh.read(module.NAME),outer.rh.read(module.RECEIPT)
        if not receipt['all_passed'] or not all(receipt[k] for k in module.GATES):raise ValueError('Accepted actual source required '+module.RECEIPT)
        for name,digest in receipt['input_hashes'].items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Conflicting original source '+name)
            hashes[name]=digest
        hashes[module.NAME]=sha(module.NAME);hashes[module.RECEIPT]=sha(module.RECEIPT)
        reports[module.NAME]=report
    identity=reports[band.NAME]['source_family']
    if any(row['source_family']!=identity for row in reports.values()):raise ValueError('Same original source family required')
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Changed accepted original source '+name)
    return reports,hashes


def positive_archived_amplitude_log(c,record):
    scale=record['formal_positive_scale']
    if any(scale['source_exponents']) or scale['radius_power']!=0:
        raise ValueError('Same normalized endpoint amplitude with factored constant log scale required')
    coefficient=current.packets.interval(c,record['coefficient_interval'])
    if ep(coefficient)[0]<=0 or record['sign']!='positive':raise ValueError('Actual strictly positive amplitude enclosure required')
    logscale=current.packets.interval(c,scale['additional_log_interval'])
    return c.mpf(ep(logscale+c.ln(coefficient))[0])


def constant_interval(c,g,node):
    row=g.nodes[node];op=row['operation'];at=lambda i:constant_interval(c,g,i)
    if op=='exact_rational':return c.mpf(row['numerator'])/row['denominator']
    if op=='sum':return sum((at(i) for i in row['arguments']),c.mpf(0))
    if op=='product':
        result=c.mpf(1)
        for i in row['arguments']:result*=at(i)
        return result
    if op=='negative':return -at(row['argument'])
    if op=='analytic_unary' and row['name'] in ('exp','log'):return getattr(c,row['name'])(at(row['argument']))
    raise ValueError('Exact constant original parameter expression required '+str(row))


def power_flow(g,A,own,mu,s):
    alg=target.C3Algebra(g);alpha=g.add(g.constant('1/2'),mu)
    decay=lambda r:g.unary('exp',g.neg(g.mul(g.constant(r),s)))
    mass=lambda k:g.mul(s,g.unary('exprel',g.neg(g.mul(k,s))))
    E=alg.scale(A,g.unary('exp',g.neg(g.mul(alpha,s))));AA=alg.mul(A,A);zero=alg.fixed(0)
    theta=g.mul(decay('3/2'),s,g.unary('exprel',g.mul(g.sub(g.one,mu),s)))
    hist=dict(m=alg.scale(own['m'],decay(1)),
        h=alg.add(alg.scale(own['h'],decay('3/2')),alg.scale(A,theta)),
        k=alg.scale(own['k'],decay('3/2')),
        e=alg.add(alg.scale(own['e'],decay(1)),alg.scale(AA,g.mul(g.constant('-1/2'),decay(1),mass(g.mul(g.constant(2),mu))))),
        p=alg.add(own['p'],alg.scale(AA,g.mul(g.constant('1/2'),mass(g.add(g.one,g.mul(g.constant(2),mu)))))))
    profiles=[dict(E=alg.scale(E,g.mul(g.constant((-1)**j),*[alpha]*j)),V=zero) for j in range(5)]
    histories=[hist]
    # Analytic density y derivatives: E_y^j=(-alpha)^j E and (E^2)_y^j=(-2alpha)^j E^2.
    EE=alg.mul(E,E)
    for j in range(4):
        densities=dict(m=zero,h=profiles[j]['E'],k=zero,
            e=alg.scale(EE,g.mul(g.constant('-1/2'),g.constant((-2)**j),*[alpha]*j)),
            p=alg.scale(EE,g.mul(g.constant('1/2'),g.constant((-2)**j),*[alpha]*j)))
        histories.append({key:alg.add(densities[key],alg.scale(histories[-1][key],-rate)) for key,rate in current.RATES.items()})
    return dict(profiles_y0_y1_y2_y3_y4_C3=profiles,complete_histories_y0_y1_y2_y3_y4_C3=histories,
        original_alpha=alpha,all_ordinary_y_orders=list(range(5)),ordinary_Z_orders=list(range(4)))


def pulse_input(g,E,own,P0):
    alg=target.C3Algebra(g);EE=alg.mul(E,E)
    return dict(u=E,m1=alg.div(own['m'],E),m2=alg.div(own['k'],EE),X=alg.div(own['h'],E),
        energy=alg.div(own['e'],EE),Mp=own['p'],P0=P0,pressure=alg.add(P0,own['p']),
        exact_common_m_k_already_normalized_by_Pstar=True,extra_Pstar_division_forbidden=True,
        quotient_denominator_is_actual_positive_swirl_not_a_range_value=True,
        native_selected_pulse_consumes_current_frame_not_yet_identified=True)


def build(field):
    g=field.graph;alg=target.C3Algebra(g);ref=lambda i:source.FunctionRef(g,i)
    jet=lambda row:target.C3Function(*map(ref,row));bf=field.band_functions
    x=ref(bf['original_band_variable']);two=g.constant(2);s=g.symbol('current_C3_post_2Rc_power_s')
    at2=lambda q:outer.substitute(g,q,x,two)
    complete2={key:at2(jet(row)) for key,row in bf['complete_histories'].items()}
    leading2={key:at2(jet(row)) for key,row in bf['original_leading_power']['leading_histories'].items()}
    relative2={key:jet(row) for key,row in bf['relative_terminal_actual'].items()}
    profile=bf['profiles_at_band_x']['radial_y_rows'][0]
    actualE2,actualV2=at2(jet(profile['E'])),at2(jet(profile['V']))
    A2=at2(jet(profile['original_E']));P0=jet(bf['original_leading_power']['independent_P0'])
    params={key:ref(node) for key,node in field.control_base['parameters'].items()};mu=params['mu']
    Tw=g.mul(g.constant(-60),g.unary('log',mu));length=g.sub(g.sub(Tw,two),g.unary('log',two))
    actual_flow=power_flow(g,actualE2,complete2,mu,s);leading_flow=power_flow(g,A2,leading2,mu,s)
    endpoint=lambda q:outer.substitute(g,q,s,length)
    terminal={key:endpoint(q) for key,q in actual_flow['complete_histories_y0_y1_y2_y3_y4_C3'][0].items()}
    Ep=endpoint(actual_flow['profiles_y0_y1_y2_y3_y4_C3'][0]['E'])
    Rc_offset=ref(field.control_base['original_Rc_offset']);Rc=g.unary('exp',Rc_offset)
    R2=g.mul(two,Rc);R=g.mul(R2,g.unary('exp',s));Rp=g.mul(R2,g.unary('exp',length))
    return dict(original_band_x=x,quiet_coordinate=s,parameters=params,original_N=ref(field.control_base['N']),
        original_Rc_offset=Rc_offset,Rc=Rc,R2=R2,R=R,Rp=Rp,Tw=Tw,quiet_length=length,
        actual_repaired_exit_complete_C3=complete2,original_leading_exit_C3=leading2,
        actual_repaired_exit_relative_C3=relative2,actual_repaired_exit_E_C3=actualE2,
        actual_repaired_exit_V_C3=actualV2,original_leading_exit_E_C3=A2,original_independent_P0_C3=P0,
        accepted_original_C3_relative_zero_certificates=bf['relative_terminal_zero_certificates'],
        actual_quiet_flow_C3=actual_flow,original_leading_quiet_flow_C3=leading_flow,
        actual_terminal_Rp_histories_C3=terminal,actual_terminal_Rp_swirl_C3=Ep,
        actual_C3_pulse_input_frame=pulse_input(g,Ep,terminal,P0),
        coordinate_identity='s=log(R/(2Rc)); 0<=s<=Tw-2-log2; Rp=Rc*exp(Tw-2)',
        compact_repair_at_Rc_to_2Rc_consumed=True,quiet_relative_zero_is_not_absolute_history_zero=True,
        original_outer_O3_source_ends_at_Rc_not_2Rc=True,
        selected_native_pulse_constructor_consumes_current_frame=False,
        current_outer_leading_endpoint_equals_band_seed_as_functions_not_yet_admitted=True,
        global_time_Cartesian_heat_cone_and_temporal_recursion_not_admitted=True)


class CurrentC3PowerToRp:
    def __init__(self,require_checked=True):
        reports,self.hashes=accepted_inputs();self.raw=reports[band.NAME];self.control_report=reports[band.control.NAME]
        self.outer_report=reports[outer.NAME];self.identity=self.raw['source_family']
        self.graph=outer.rh.restore_graph(self.raw);self.prefix=copy.deepcopy(self.graph.nodes)
        self.band_functions=self.raw['actual_C3_repair_band_functions']
        self.control_base=self.control_report['actual_C3_limit_controls_and_Picard_functions']['exact_C2_limit_functions']['exact_C1_limit_adapter']
        control_nodes=self.control_report['exact_graph_nodes']
        if self.graph.nodes[:len(control_nodes)]!=control_nodes:raise ValueError('Same original control graph prefix required')
        selected=self.raw['actual_selected_repair_integer']
        original_control=outer.rh.read(outer.rh.controls.NAME)['actual_quantitative_C3_bounds']['actual_same_repair_integer']
        if selected!=original_control:raise ValueError('Same original selected repair integer required')
        self.c=outer.rh.MPIntervalContext();self.c.dps=540
        self.rows={tuple(row['exact_Z_cell']):row for row in self.raw['actual_four_Z_C3_band_ranges']}
        range_report=outer.rh.read(ranges.NAME)
        self.source_ranges={tuple(row['exact_Z_cell']):row for row in range_report['actual_four_Z_C3_target_transports']}
        self.range_source_binding=self.bind_range_source()
        self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=outer.rh.read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C3 quiet power required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual quiet power source '+name)
            self.acceptance_loaded=True

    def bind_range_source(self):
        g=self.graph;profile=self.band_functions['profiles_at_band_x'];power=profile['original_power'];amplitudes=[]
        source_report=outer.rh.read(target.NAME)
        source_nodes=source_report['exact_graph_nodes']
        source_amplitudes=source_report['actual_C3_density_transport_targets']['actual_terminal_amplitude']
        identity_keys=('operation','source_family','native_chart','recipe','quantity','Z_order','coordinate','Z_variable')
        for j,node in enumerate(profile['radial_y_rows'][0]['original_E']):
            row=g.nodes[node]
            if row['operation']!='product' or len(row['arguments'])!=2 or power not in row['arguments']:
                raise ValueError('Original E=A*same power function required')
            amplitude=next(i for i in row['arguments'] if i!=power);leaf=g.nodes[amplitude]
            if (leaf['operation']!='current_original_leading_function_recipe' or leaf['source_family']!=self.identity
                or leaf['native_chart']!='O3_power' or leaf['quantity']!='E' or leaf['Z_order']!=j
                or leaf['coordinate']!=g.constant(2).node
                or leaf['recipe']['function']!='WholeZAllNOuterRcFunctions.leading_packet'
                or leaf['recipe']['quantity_paths']['E']!='original_generic_source.common_velocity_E_axial5[Z_order]'
                or leaf.get('Taylor_coefficient_factorial',1)!=(1,1,2,6)[j]):
                raise ValueError('Same original ordinary O3 endpoint amplitude function required')
            canonical=source_nodes[source_amplitudes[j]]
            if (any(leaf.get(k)!=canonical.get(k) for k in identity_keys)
                or leaf.get('Taylor_coefficient_factorial',1)!=canonical.get('Taylor_coefficient_factorial',1)):
                raise ValueError('Band amplitude and original endpoint projection must define the same function row')
            amplitudes.append(amplitude)
        mu=self.control_base['parameters']['mu'];node=g.nodes[mu]
        if node['operation']!='analytic_unary' or node['name']!='exp':raise ValueError('Exact original positive mu expression required')
        logmu=constant_interval(self.c,g,node['argument'])
        for ends,row in self.rows.items():
            p=row['actual_endpoint_source']
            if (p['source_family']!=self.identity or p['original_Rc_source_coordinate']!=2
                or not p['m_k_S_normalization_already_applied'] or not p['third_Taylor_factorial_applied_once']):
                raise ValueError('Same normalized endpoint source projection required')
            saved=current.packets.interval(self.c,self.source_ranges[ends]['actual_positive_mu_log'])
            if ep(saved)[0]>ep(logmu)[0] or ep(logmu)[1]>ep(saved)[1]:raise ValueError('Archived mu must enclose the actual graph expression')
        return dict(actual_amplitude_source_nodes=amplitudes,actual_mu_source_node=mu,
            canonical_original_endpoint_amplitude_nodes=source_amplitudes,
            actual_and_canonical_amplitude_recipe_row_identity=True,
            actual_mu_log_source_node=node['argument'],actual_mu_log=logmu,
            original_endpoint_source_projection=current.ast_binding(band.CurrentC3RepairBand.endpoint_source),
            same_source_recipe_ordinary_endpoint_amplitude_and_normalization=True,
            actual_mu_graph_expression_enclosed_by_every_original_source_cell=True)

    def quantitative_range(self,ends):
        ends=tuple(ends)
        if ends not in self.rows:raise ValueError('Admitted original axial cell required')
        c=self.c;bd=ranges.Bounds(c);row=self.rows[ends];source_row=self.source_ranges[ends]
        read=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(read(rec) for rec in rows[:4]))
        # Accepted endpoint source Taylor rows are already ordinary C3 and normalized.
        A=jet(row['actual_endpoint_source']['ordinary_terminal_amplitude_C3'])
        logA=positive_archived_amplitude_log(c,row['actual_endpoint_source']['ordinary_terminal_amplitude_C3'][0])
        logmu=self.range_source_binding['actual_mu_log'];mu=c.exp(logmu);alpha=c.mpf('.5')+mu
        length=-60*logmu-2-c.ln(2)
        if ep(mu)[0]<=0 or ep(mu)[1]>=ep(c.mpf(1)/6)[0] or ep(length)[0]<=0:raise ValueError('Actual positive quiet width and original small mu required')
        # Relative corrections vanish at 2Rc; absolute histories remain the original leading exit.
        initial={key:jet(q) for key,q in row['actual_complete_history_C3_bounds'].items()}
        AA=bd.product(A,A);half=bd.constant(c.mpf('.5'));zero=bd.fixed(bd.zero)
        profiles=[dict(E=bd.scaled(A,bd.constant(alpha**j)),V=zero) for j in range(5)]
        histories=dict(initial)
        histories['h']=bd.add(initial['h'],bd.scaled(A,bd.constant(c.mpf(2)/3)))
        histories['e']=bd.add(initial['e'],bd.scaled(AA,half))
        # The full exact quiet power square integrates to A2^2/(1+2mu)<=A^2.
        histories['p']=bd.add(initial['p'],bd.scaled(AA,half));history_rows=[histories]
        for j in range(4):
            densities=dict(m=zero,h=profiles[j]['E'],k=zero,
                e=bd.scaled(AA,bd.constant((2*alpha)**j/2)),p=bd.scaled(AA,bd.constant((2*alpha)**j/2)))
            history_rows.append({key:bd.add(densities[key],bd.scaled(history_rows[-1][key],bd.constant(c.mpf(rate.numerator)/rate.denominator)))
                for key,rate in current.RATES.items()})
        P0=jet(row['actual_endpoint_source']['ordinary_independent_P0_C3'])
        # Same actual A>0, E(Rp)=A*exp(-alpha*(Tw-2)). Use outward lower log, never choose a field value.
        logEp=c.mpf(ep(logA-alpha*(-60*logmu-2))[0]);Ecap=profiles[0]['E'];E2=bd.product(Ecap,Ecap)
        frame=dict(u=Ecap,m1=bd.quotient(histories['m'],Ecap,logEp),
            m2=bd.quotient(histories['k'],E2,2*logEp),X=bd.quotient(histories['h'],Ecap,logEp),
            energy=bd.quotient(histories['e'],E2,2*logEp),Mp=histories['p'],P0=P0,pressure=bd.add(P0,histories['p']))
        # A tiny explicit outward slack makes independent operation orderings
        # dominate the same analytic majorants despite directed rounding.
        slack=bd.constant(1+c.mpf('1e-20'))
        profiles=[{key:bd.scaled(q,slack) for key,q in record.items()} for record in profiles]
        history_rows=[{key:bd.scaled(q,slack) for key,q in record.items()} for record in history_rows]
        frame={key:bd.scaled(q,slack) for key,q in frame.items()}
        return dict(source_family=self.identity,exact_Z_cell=ends,actual_positive_mu=mu,actual_alpha=alpha,
            actual_positive_quiet_length=length,actual_terminal_swirl_log_lower=logEp,
            actual_repaired_exit_complete_C3_bounds=ranges.record(initial),
            actual_whole_quiet_profile_y0_y1_y2_y3_y4_C3_bounds=ranges.record(profiles),
            actual_whole_quiet_history_y0_y1_y2_y3_y4_C3_bounds=ranges.record(history_rows),
            actual_independent_P0_C3_bounds=ranges.record(P0),actual_C3_pulse_frame_bounds=ranges.record(frame),
            relative_exit_zero_does_not_zero_complete_histories=True,positive_swirl_denominators_are_actual_source_functions=True,
            whole_quiet_domain_and_terminal_Rp_covered=True,range_caps_not_function_values=True)


def source_bindings():
    return dict(build=current.ast_binding(build),actual_C3_power_semigroup=current.ast_binding(power_flow),
        actual_positive_archived_amplitude=current.ast_binding(positive_archived_amplitude_log),
        exact_constant_parameter_interval=current.ast_binding(constant_interval),
        actual_range_source_node_binding=current.ast_binding(CurrentC3PowerToRp.bind_range_source),
        native_C3_pulse_units=current.ast_binding(pulse_input),quantitative_range=current.ast_binding(CurrentC3PowerToRp.quantitative_range),
        original_C3_repair_band=current.ast_binding(band.build),original_C3_repaired_source_projection=current.ast_binding(band.CurrentC3RepairBand.endpoint_source))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC3PowerToRp(require_checked=False);rows=[field.quantitative_range(ends) for ends in field.rows]
        report=dict(candidate_actual_C3_power_to_Rp_constructed=True,source_family=field.identity,
            accepted_original_C3_band_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.graph.nodes,
            actual_C3_repaired_exit_power_and_pulse_functions=target.encoded(field.functions),
            actual_range_source_node_binding=ranges.record(field.range_source_binding),
            actual_four_Z_C3_power_ranges=ranges.record(rows),source_bindings=source_bindings(),
            **dict.fromkeys(GATES+OPEN,False),selected_native_pulse_constructor_consumes_current_C3_frame=False,
            current_outer_leading_endpoint_band_seed_function_identity_installed=False,
            current_numeric_point_field_oracle_installed=False,global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(target.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual C3 compact-repaired exit, quiet y4/Z3 flow and native Rp frame constructed',flush=True)
    return field


if __name__=='__main__':run()
