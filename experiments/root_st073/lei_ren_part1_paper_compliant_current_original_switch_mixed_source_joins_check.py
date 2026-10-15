"""Actual switch-control jets, original inlet memories and 492 own rows."""
import json
from math import comb
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins as current
import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions_check as prior
import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions as loop_source

rows,encode=current.outer.rows,current.outer.encode
zero,bind,Algebra=prior.zero,prior.bind,prior.Algebra
bridge=prior.prior.bridge


def source_contract(owner):
    owner.assert_graph();g=owner.g;rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and len(rec['checks'])==7
    assert rec['actual_native_to_physical_derivative_requests']==180
    assert rec['checks']['actual_signed_transport_sources']['every_actual_predecessor_pressure_and_five_correction_history_retained']
    assert rec['checks']['actual_recovery_and_pressure_sources']['actual_mean_ODE_and_actual_saved_reduced_S_on_all_ten_charts']
    assert owner.base.acceptance_loaded and tuple(owner.charts)==current.previous.CHARTS
    for key,left,right,endpoint in current.JOINS:
        saved=owner.joins[key]
        assert saved['original_endpoint']==endpoint and saved['left_chart']==left and saved['right_chart']==right
        assert encode(saved['left'])==encode(owner.base.at_own(left,g.one))
        assert encode(saved['right'])==encode(owner.base.at_own(right,g.zero))
        assert owner.loop[right]['incoming']==owner.loop_input['endpoints'][left]
        assert owner.modified[right]['incoming_integrals']==owner.modified_input['endpoints'][left]['integrals']
    for name,refs in owner.controls.items():
        assert refs[0]==current.downstream_source.sigma(g,owner.charts[name]['coordinate'])
        for n in range(1,5):assert owner.factories[name].cache[refs[n-1].node]==refs[n]
    assert owner.charts['first_switch']['jacobian']==owner.charts['second_switch']['jacobian']
    assert owner.charts['post_power']['jacobian']!=owner.charts['second_switch']['jacobian']
    return dict(passed=True,accepted_ten_full_chart_source_receipt=current.previous.RECEIPT,
        all_180_actual_high_derivative_requests_reused_without_rebuilding=True,
        actual_own_endpoint_pairs_and_complete_predecessor_sources_bound=True,
        shared_switch_width_and_distinct_post_power_native_length=True)


def comparison(owner,name):
    g=owner.g;b=bridge.BridgeAlgebra(g,7);bg=owner.backgrounds[name];R=owner.charts[name]['radius']
    old=owner.bridge['comparison_macro'];m=owner.ref(old['phase'])
    phi,U=owner.jet(old['phi']),owner.jet(old['U'])
    incoming={key:b.at(owner.jet(ids),m,g.one) for key,ids in old['moments'].items()}
    distance=g.unary('log',g.quotient(R,g.constant(100),'same positive current radius over100'))
    eq=dict(H=phi,M=U,K=b.mul(phi,U),A=b.mul(U,U),B=b.scale(b.mul(phi,phi),'1/2'),C=b.mul(phi,phi))
    histories={}
    for key,rate in bridge.RATES.items():
        decay=g.unary('exp',g.mul(g.constant(-rate),distance))
        histories[key]=b.add(b.scale(incoming[key],decay),b.scale(eq[key],g.sub(g.one,decay)))
    F0=owner.jet(owner.bridge['core7']['F0']);P0=owner.jet(owner.bridge['P0'])
    known=bridge.comparison_directions(g,phi,U,histories,F0,P0,owner.parameters['delta'])
    assert encode(known)==bg['known_comparison_direction']
    return known,phi,F0,P0


def actual_background_recipes(owner):
    g=owner.g;a=bridge.BridgeAlgebra(g,6);down=current.downstream_source
    h=owner.ref(owner.downstream['hb']);Pstar=owner.ref(owner.downstream['Pstar']);evidence=[]
    for name in ('first_switch','second_switch','post_power'):
        bg=owner.backgrounds[name];x=owner.charts[name]['coordinate'];R=owner.charts[name]['radius'];J=owner.charts[name]['jacobian']
        if name=='first_switch':
            incoming=owner.bridge['R100'];known,barphi,F0,P0=comparison(owner,name)
            actual=down.first_switch_fields(a,known,{key:owner.jet(incoming['actual'][key]) for key in ('phi','U')},
                x,h,R,a.trim(F0),Pstar,barphi)
        elif name=='second_switch':
            incoming=owner.downstream['endpoints']['R1'];known,_,F0,P0=comparison(owner,name)
            actual=down.second_switch_fields(a,known,{key:owner.jet(incoming['actual'][key]) for key in ('phi','U')},x,h,R)
        else:
            incoming=owner.downstream['endpoints']['R2'];y=g.mul(J,x)
            phi=a.scale(owner.jet(incoming['actual']['phi']),g.unary('exp',g.mul(g.constant('-2/5'),y)))
            actual=dict(phi=phi,U=owner.jet(incoming['actual']['U']),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0))
            F0=owner.jet(owner.bridge['core7']['F0']);P0=owner.jet(owner.bridge['P0'])
        assert encode(actual)==bg['actual']
        assert incoming['moments']==bg['incoming_histories']
        histories=bridge.histories(a,actual['phi'],actual['U'],{k:owner.jet(v) for k,v in incoming['moments'].items()},
            x,J,'true_'+('post_switch' if name=='post_power' else name))
        assert encode(histories)==bg['moments']
        fields,radial,generic=bridge.generic_fields(a,actual,histories,a.trim(F0),P0,Pstar,R,owner.parameters['delta'])
        assert encode(fields)==bg['fields'] and encode(radial)==bg['radial'] and encode(generic)==bg['generic']
        evidence.append(name)
    for key,left,right,endpoint in current.JOINS:
        saved=owner.downstream['endpoints'][endpoint];L=owner.backgrounds[left]
        for field in ('phi','U'):
            assert encode(a.at(owner.jet(L['actual'][field]),owner.charts[left]['coordinate'],g.one))==saved['actual'][field]
        for field,ids in L['moments'].items():
            assert encode(a.at(owner.jet(ids),owner.charts[left]['coordinate'],g.one))==saved['moments'][field]
        assert owner.backgrounds[right]['incoming_histories']==saved['moments']
    return dict(passed=True,actual_field_six_history_and_generic_recipes_replayed=evidence,
        actual_first_switch1_and_second_switch1_phi_U_and_six_histories_are_successor_inputs=True,
        same_current_radius_comparison_formula_with_nonzero_inlets_and_B_half_equilibrium=True,
        exact_hydro_pressure_swirl_F0_Pstar_and_physical_measures_preserved=True)


def physical_controls_and_ODE_jets(owner):
    g=owner.g;ev=prior.current.previous.previous.modified_source.lexical.Algebra;flat=[]
    for name,value in (('first_switch',1),('second_switch',0),('second_switch',1)):
        for n,ref in enumerate(owner.controls[name]):
            zero(ev(g,{owner.charts[name]['coordinate'].node:sy.Integer(value)}).at(ref.node)
                -(value if n==0 else 0),('actual sigma flat endpoint',name,value,n))
            flat.append([name,value,n,ref.node])
    # Evaluate the actual stored physical RHS at arbitrary values of its
    # dependent variables. This extracts the source coefficients; it does
    # not assign any derivative or own output a desired seam value.
    y=sy.Symbol('independent_physical_y',real=True);Z=sy.Symbol('Z',real=True)
    H,R0,P,F,BP=sy.symbols('actual_hb actual_seam_R Pstar F0 barPhi',positive=True)
    Phi=sy.Function('unknown_actual_phi')(y);D=sy.Function('current_radius_D')(y)
    drives={k:sy.Function('current_'+k)(y) for k in ('drive_hydro','drive_pressure','drive_swirl')}
    coefficient_rows=[];induction=[]
    for key,left,right,_ in current.JOINS:
        models=[]
        for name,endpoint in ((left,1),(right,0)):
            bg=owner.backgrounds[name];chart=owner.charts[name]
            env={owner.downstream['hb']:H,chart['radius'].node:R0*sy.exp(y),bg['actual']['phi'][0]:Phi,
                owner.downstream['Pstar']:P,owner.bridge['core7']['F0'][0]:F}
            if name in owner.controls:
                # Arbitrary fifth-order tails represent exactly the finite
                # endpoint jets already read above from the actual source.
                env[owner.controls[name][0].node]=endpoint+sy.Symbol('free_flat_tail_'+name)*y**5
                env[bg['known_comparison_direction']['D_over_R'][0]]=D
                for field,body in drives.items():env[bg['known_comparison_direction'][field][0]]=body
                env[owner.bridge['comparison_macro']['phi'][0]]=BP
            alg=prior.geometry_check.Algebra(g,env)
            slope=alg.at(bg['actual']['phi_y'][0])/Phi;axial=alg.at(bg['actual']['U_y'][0])
            models.append(([sy.diff(slope,y,n).subs(y,0) for n in range(5)],
                [sy.diff(axial,y,n).subs(y,0) for n in range(5)]))
        for family in range(2):
            for n,(L,R) in enumerate(zip(models[0][family],models[1][family])):
                zero(L-R,('actual angular or axial physical RHS jets',key,family,n))
                coefficient_rows.append([key,family,n])
        # Finite physical-y Leibniz/ODE induction, arbitrary nonzero inlet
        # values and arbitrary smooth Z functions. Five background orders
        # are required because the original loop roots contain E_y/V_y.
        all_states=[]
        for angular,axial in models:
            phi=[sy.Function('seam_phi')(Z)];U=[sy.Function('seam_U')(Z)]
            for n in range(5):
                phi.append(sum(comb(n,k)*angular[k]*phi[n-k] for k in range(n+1)))
                U.append(axial[n])
            product=lambda A,B,n:sum(comb(n,k)*A[k]*B[n-k] for k in range(n+1))
            histories={field:[sy.Function('nonzero_'+field)(Z)] for field in bridge.RATES}
            for n in range(5):
                density=dict(H=2*phi[n],M=U[n],K=2*product(phi,U,n),A=product(U,U,n),
                    B=product(phi,phi,n),C=product(phi,phi,n))
                for field,rate in bridge.RATES.items():histories[field].append(density[field]-rate*histories[field][-1])
            all_states.append(dict(phi=phi,U=U,**histories))
        for field in all_states[0]:
            for n,(L,R) in enumerate(zip(all_states[0][field],all_states[1][field])):
                zero(L-R,('source physical IVP induction',key,field,n));induction.append([key,field,n])
    return dict(passed=True,actual_native_sigma_flat_rows=flat,
        actual_physical_RHS_coefficient_rows=coefficient_rows,independent_physical_ODE_rows=induction,
        first_axial_driver_turns_off_with_four_flat_derivatives=True,
        second_angular_driver_turns_into_exact_minus_two_fifths_power_slope=True,
        background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets=True,
        ordinary_Z_differentiation_commutes_with_the_finite_smooth_source_induction=True,
        no_dependent_derivative_or_own_pressure_output_bound=True)


def actual_loop_and_phase_transfer(owner):
    g=owner.g;a=current.outer.Algebra(g,5);Pstar=owner.ref(owner.downstream['Pstar']);roots=[]
    for name in ('first_switch','second_switch','post_power'):
        bg=owner.backgrounds[name];chart=owner.loop[name]
        packet={group:{key:owner.jet(ids) for key,ids in bg[group].items()} for group in ('fields','radial','generic')}
        packet['radius']=owner.ref(bg['radius'])
        root=loop_source.generic_roots(g,packet,Pstar)
        assert encode(root)==chart['roots'];inverse=g.nodes[chart['loop']['inverse'][0]]
        assert inverse['operation']=='original_monotone_phase_inverse'
        assert inverse['phase_function']==chart['loop']['Phi'][0] and inverse['target_phase']==chart['actual_phase']
        assert inverse['true_prescribed_source_C5_implicit_inverse']
        offset=g.sub(g.unary('log',g.quotient(owner.ref(chart['radius']),owner.ref(owner.bridge['Ra']),
            'same positive actual R/Ra')),owner.ref(owner.loop_input['phase_origin']))
        assert offset.node==chart['global_offset']
        assert g.unary('fractional_part',g.mul(owner.N,offset)).node==chart['actual_phase']
        D=owner.factories[name];J=owner.charts[name]['jacobian']
        assert D.cache[owner.charts[name]['phase'].node]==g.mul(owner.N,J)
        for factory in D.shared['factories'].values():
            assert factory.cache[owner.charts[name]['phase'].node]==(g.mul(owner.N,J)
                if factory.variable==owner.charts[name]['coordinate'] else g.zero)
        roots.append([name,chart['roots'],chart['actual_phase'],J.node])
    H,Ra=sy.symbols('actual_hb actual_Ra',positive=True);geometry=[]
    for key,left,right,_ in current.JOINS:
        env={owner.downstream['hb']:H,owner.bridge['Ra']:Ra,
            owner.loop_input['phase_origin']:sy.Symbol('same_original_global_phase_origin',real=True),
            owner.charts[left]['coordinate'].node:sy.Integer(1),owner.charts[right]['coordinate'].node:sy.Integer(0)}
        alg=prior.geometry_check.Algebra(g,env)
        zero(alg.at(owner.joins[key]['left_radius'].value.node)-alg.at(owner.joins[key]['right_radius'].value.node),
            ('actual switch radius',key))
        zero(alg.at(owner.loop[left]['global_offset'])-alg.at(owner.loop[right]['global_offset']),
            ('actual same global logarithmic phase',key));geometry.append(key)
    old=current.parent.linear.read(current.previous.previous.RECEIPT)['checks']['actual_smooth_source_seam']
    assert old['original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows']
    rec=current.parent.linear.read(current.previous.RECEIPT)['checks']
    assert rec['actual_recovery_and_pressure_sources']['actual_mean_ODE_and_actual_saved_reduced_S_on_all_ten_charts']
    return dict(passed=True,actual_original_root_inverse_phase_J_bindings=roots,
        actual_radius_and_global_logarithmic_phase_seams=geometry,
        background_y5_supplies_root_y4_and_original_unique_inverse_periodic_lift_U_y4_C5=True,
        actual_nonzero_mean_inlet_is_checked_before_Q_y4_N_y2_pressure_y3_transfer=True,
        actual_Q_N_and_reduced_S_operator_theorem_consumed_from=current.previous.RECEIPT,
        no_raw_fractional_part_differentiation_or_numerical_enclosure_claim=True)


def actual_own_mixed_joins(owner):
    g=owner.g;evidence=[];mean_evidence=[];protected=[]
    for key,left,right,_ in current.JOINS:
        env,Z,_=prior.prior.active.prior.rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
        Lchart,Rchart=owner.charts[left],owner.charts[right]
        x=sy.Symbol(g.nodes[Lchart['coordinate'].node]['name'],real=True)
        t=sy.Symbol(g.nodes[Rchart['coordinate'].node]['name'],real=True)
        H,L=sy.symbols('same_actual_switch_width actual_post_power_length',positive=True)
        env[owner.downstream['hb']]=H
        Jleft=H;Jright=H if right=='second_switch' else L
        if right=='post_power':env[Rchart['jacobian'].node]=L
        Rstar=100*sy.exp(H if left=='first_switch' else 2*H)
        env[Lchart['radius'].node]=Rstar*sy.exp(Jleft*(x-1))
        env[Rchart['radius'].node]=Rstar*sy.exp(Jright*t)
        env[owner.rminus.active_sc.node]=Algebra(owner,env).at(owner.parameters['sc'])
        earlier=current.previous.CHARTS[:current.previous.CHARTS.index(left)+1]
        for name in earlier:
            coordinate=sy.Symbol(g.nodes[owner.charts[name]['coordinate'].node]['name'],real=True)
            for field,ids in owner.loop[name]['contributions'].items():
                bind(env,owner.jet(ids),sy.Function('retained_'+name+'_correction_'+field)(coordinate,Z),Z)
            for field,ids in owner.modified[name]['local_integrals'].items():
                bind(env,owner.jet(ids),sy.Function('retained_'+name+'_signed_'+field)(coordinate,Z),Z)
            if name!=left:
                for field in ('Q','Q_y'):
                    bind(env,owner.jet(owner.modified[name]['physical'][field]),
                        sy.Function('retained_'+name+'_leading_'+field)(coordinate,Z),Z)
        mapped=1+Jright*t/Jleft;mean=sy.Function('common_background_mean_'+key)
        bind(env,owner.jet(owner.backgrounds[left]['fields']['m']),mean(x,Z),Z)
        bind(env,owner.jet(owner.backgrounds[right]['fields']['m']),mean(mapped,Z),Z)
        for name,arg in ((left,x),(right,mapped)):
            for field in ('Q','Q_y'):
                bind(env,owner.jet(owner.modified[name]['physical'][field]),
                    sy.Function('proved_common_leading_'+key+'_'+field)(arg,Z),Z)
            for n,jet in enumerate(owner.charts[name]['N_y']):
                bind(env,jet,sy.Function('proved_common_leading_'+key+'_N'+str(n))(arg,Z),Z)
        # These substitutions are finite seam germs justified by preceding
        # groups, not assertions that the charts agree in a neighborhood.
        # Own pressure and m13/m15 output nodes remain unbound.
        for side in ('left','right'):
            packet=owner.joins[key][side]
            for group,field in (('fields','P1'),('cumulative','m13'),('cumulative','m15')):
                for jet in packet[group][field]:
                    for ref in rows(jet):
                        assert ref.node not in env,('Own pressure output was bound',key,side,field,ref.node)
                        protected.append(ref.node)
        alg=Algebra(owner,env);means=[]
        for name,value in ((left,g.one),(right,g.zero)):
            chart=owner.charts[name]
            jet=current.outer.Algebra(g,5).at(chart['physical']['m'][0],chart['coordinate'],value)
            means.append([alg.at(ref) for ref in rows(jet)])
        for n,(Lvalue,Rvalue) in enumerate(zip(*means)):
            zero(Lvalue-Rvalue,('actual leading mean inlet',key,n));mean_evidence.append([key,n])
        for spec in owner.contract:
            for n in range(spec['physical_y_order']+1):
                left_rows=rows(owner.joins[key]['left'][spec['group']][spec['field']][n])
                right_rows=rows(owner.joins[key]['right'][spec['group']][spec['field']][n])
                for k,(Lref,Rref) in enumerate(zip(left_rows,right_rows)):
                    zero(alg.at(Lref)-alg.at(Rref),('actual switch own mixed join',key,spec['group'],spec['field'],n,k))
                    evidence.append(dict(interface=key,group=spec['group'],field=spec['field'],
                        physical_y_order=n,ordinary_Z_order=k,left_node=Lref.node,right_node=Rref.node))
    assert len(evidence)==492 and len(mean_evidence)==12
    return dict(passed=True,actual_two_switch_own_mixed_rows=evidence,actual_leading_mean_C5_inlet_rows=mean_evidence,
        explicit_unbound_own_pressure_output_nodes=len(protected),
        complete_nonzero_original_correction_S0_S1_pressure_and_five_histories_retained=True,
        finite_common_leading_germs_only_after_source_ODE_inverse_Q_N_proofs=True,
        distinct_post_power_length_and_one_plus_delta_own_recovery_retained=True)


CHECKS=(source_contract,actual_background_recipes,physical_controls_and_ODE_jets,
        actual_loop_and_phase_transfer,actual_own_mixed_joins)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner)
        print('SWITCH_JOIN_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['switch_mixed_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_switch_interfaces=[key for key,*_ in current.JOINS],actual_own_mixed_rows=492,
        exact_graph_nodes=len(raw['exact_graph_nodes']),numerical_oracle_available=False,global_n1=False,
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_TWO_SWITCH_INTERFACES_492_OWN_ROWS_ACCEPTED',flush=True)
    return rec
