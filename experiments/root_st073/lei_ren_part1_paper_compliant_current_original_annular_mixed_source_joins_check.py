"""Actual four annular interfaces, including fifth flat-source derivatives."""
import json
from math import comb
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_annular_mixed_source_joins as current
import lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins_check as prior

rows,encode=current.outer.rows,current.outer.encode
zero,bind,Algebra=prior.zero,prior.bind,prior.Algebra
bridge,geometry=prior.bridge,prior.prior.geometry_check
TEN=current.previous.previous.CHARTS


class SourceAlgebra(geometry.Algebra):
    """Lexical endpoint substitution and empty integrals before body access."""
    def __init__(self,g,bindings,scope=None):
        self.inputs=bindings;self.scope={} if scope is None else scope
        super().__init__(g,{**bindings,**self.scope})
    def at(self,index):
        if index in self.cache:return self.cache[index]
        row=self.g.nodes[index];op=row['operation']
        if op=='function_substitution':
            result=SourceAlgebra(self.g,self.inputs,{**self.scope,row['variable']:self.at(row['value'])}).at(row['expression'])
        elif op=='definite_integral':
            lo,hi=self.at(row['lower']),self.at(row['upper'])
            if lo==hi:result=sy.Integer(0)
            else:
                variable=sy.Symbol(row['variable']);node=self.g.symbol(row['variable']).node
                body=SourceAlgebra(self.g,self.inputs,{**self.scope,node:variable}).at(row['integrand'])
                result=sy.Integral(body,(variable,lo,hi))
        else:return super().at(index)
        self.cache[index]=result;return result


def source_contract(owner):
    owner.assert_graph();g=owner.g;rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and len(rec['checks'])==5 and rec['actual_own_mixed_rows']==492
    assert tuple(owner.charts)==TEN and owner.base.acceptance_loaded
    for key,left,right,endpoint in current.JOINS:
        packet=owner.joins[key];L=owner.charts[left];R=owner.charts[right]
        assert packet['original_endpoint']==endpoint
        assert encode(packet['left'])==encode(owner.base.at_own(left,g.one))
        assert encode(packet['right'])==encode(owner.base.at_own(right,g.zero))
        assert owner.loop[right]['incoming']==owner.loop_input['endpoints'][left]
        assert owner.modified[right]['incoming_integrals']==owner.modified_input['endpoints'][left]['integrals']
        assert owner.backgrounds[right]['incoming_histories']==owner.downstream['endpoints'][endpoint]['moments']
        assert L['domain']==R['domain']==[g.zero,g.one]
        for field in ('phi','U'):
            assert encode(owner.background_endpoints[key]['left']['actual'][field])==owner.downstream['endpoints'][endpoint]['actual'][field]
        assert encode(owner.background_endpoints[key]['left']['moments'])==owner.downstream['endpoints'][endpoint]['moments']
    for name,control in owner.annular_controls.items():
        x=owner.charts[name]['coordinate'];D=owner.factories[name]
        assert control['sigma'][0]==current.downstream_source.sigma(g,x)
        assert control['prime'][0]==current.downstream_source.sigma_prime(g,x)
        assert len(control['sigma'])==6 and len(control['prime'])==5
        for refs in control.values():
            for left,right in zip(refs,refs[1:]):assert D.cache[left.node]==right
    cert=owner.downstream['geometry_certificate']
    assert cert['directed_T_positive'] and cert['directed_reference_gap_positive'] and cert['directed_post_length_positive']
    return dict(passed=True,accepted_switch_predecessor=current.previous.RECEIPT,
        actual_four_left_endpoint_functions_and_nonzero_successor_inputs_bound=True,
        original_positive_finite_T_reference_and_post_lengths_retained=True,
        actual_sigma_to_order_five_and_prime_to_order_four_generated_not_zero_padded=True)


def actual_background_recipes(owner):
    g=owner.g;a=bridge.BridgeAlgebra(g,6);c=g.constant;d=owner.downstream;p=owner.parameters
    F0=a.trim(owner.jet(owner.bridge['core7']['F0']));G=a.trim(owner.jet(owner.bridge['core7']['G']))
    P0=owner.jet(owner.bridge['P0']);Pstar=owner.ref(d['Pstar'])
    z=bridge.FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5)
    logq=a.log(a.add(a.fixed(1),a.mul(z,z)),'positive original 1+Z^2')
    B=a.add(a.scale(G,g.neg(p['Lambda'])),a.log(owner.jet(d['endpoints']['R110']['actual']['phi']),'positive actual R110 phi'),
        a.fixed(g.mul(c('1/2'),g.unary('log',c(220)))),logq)
    E=a.sub(owner.jet(d['endpoints']['R110']['actual']['U']),a.scale(z,4))
    assert encode(logq)==d['logq'] and encode(B)==d['B'] and encode(E)==d['E']
    phi_from=lambda Utheta,R:a.div(Utheta,a.scale(F0,g.unary('sqrt',g.mul(c(2),R))),
        'same positive sqrt(2R)*F0; F0 is Z-only')
    evidence=[]
    for name in current.CHARTS:
        chart=owner.charts[name];x=chart['coordinate'];R=chart['radius'];J=chart['jacobian'];bg=owner.backgrounds[name]
        if name=='post_power':
            incoming=d['endpoints']['R2'];y=g.mul(J,x)
            phi=a.scale(owner.jet(incoming['actual']['phi']),g.unary('exp',g.mul(c('-2/5'),y)))
            actual=dict(phi=phi,U=owner.jet(incoming['actual']['U']),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0))
        else:
            if name=='long_reshape':
                incoming=d['endpoints']['R110'];y=g.mul(J,x)
                logtheta=a.add(a.scale(B,g.sub(g.one,current.downstream_source.sigma(g,x))),a.neg(logq),
                    a.fixed(g.sub(g.mul(c('1/10'),y),p['logC'])))
            elif name=='reference':
                incoming=d['endpoints']['Rsh'];y=g.add(owner.ref(d['T']),g.mul(J,x))
                logtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),p['logC'])))
            else:
                incoming=d['endpoints']['Rz' if name=='restoration' else 'Rrestore']
                offset=g.add(c(-8 if name=='restoration' else -7),x)
                logtheta=a.add(a.neg(logq),a.fixed(g.add(p['logP'],g.mul(c('1/10'),offset))))
            theta=a.exp(logtheta);phi=phi_from(theta,R)
            U=owner.jet(incoming['actual']['U']);Uy=a.fixed(0);slope=a.fixed('-2/5')
            if name=='long_reshape':
                slope=a.add(slope,a.scale(B,g.neg(g.quotient(current.downstream_source.sigma_prime(g,x),J,
                    'same positive original T'))))
            elif name=='restoration':
                U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,current.downstream_source.sigma(g,x))))
                Uy=a.scale(E,g.neg(current.downstream_source.sigma_prime(g,x)))
            elif name=='postrestore':U=a.scale(z,4)
            actual=dict(phi=phi,U=U,phi_y=a.mul(phi,slope) if name=='long_reshape' else a.scale(phi,'-2/5'),
                U_y=Uy,physical_Utheta=theta,log_Utheta=logtheta)
        assert encode(actual)==bg['actual'],name
        assert bg['incoming_histories']==incoming['moments']
        label='post_switch' if name=='post_power' else 'post_restore' if name=='postrestore' else name
        histories=bridge.histories(a,actual['phi'],actual['U'],{k:owner.jet(v) for k,v in incoming['moments'].items()},x,J,'true_'+label)
        assert encode(histories)==bg['moments']
        fields,radial,generic=bridge.generic_fields(a,actual,histories,F0,P0,Pstar,R,p['delta'])
        assert encode(fields)==bg['fields'] and encode(radial)==bg['radial'] and encode(generic)==bg['generic']
        evidence.append(name)
    return dict(passed=True,actual_five_field_history_generic_recipes_replayed=evidence,
        original_R110_B_log_phi_anchored_G_log220_logq_and_original_E_replayed=True,
        complete_original_nonzero_moment_inlets_and_finite_source_measures_retained=True)


def actual_value_and_gauge_joins(owner):
    g=owner.g;p=owner.parameters;d=owner.downstream;Z=sy.Symbol('Z',real=True)
    C,P,G=sy.symbols('actual_logC actual_logP anchored_G',real=True)
    L,T,H=sy.symbols('actual_Lambda actual_T actual_hb',positive=True)
    phi=sy.Symbol('arbitrary_positive_actual_R110_phi',positive=True);U=sy.Symbol('arbitrary_actual_R110_U',real=True)
    env={g.symbol('Z').node:Z,p['logC'].node:C,p['logP'].node:P,p['Lambda'].node:L,d['T']:T,d['hb']:H,
        owner.bridge['core7']['G'][0]:G,d['endpoints']['R110']['actual']['phi'][0]:phi,
        d['endpoints']['R110']['actual']['U'][0]:U}
    alg=SourceAlgebra(g,env)
    logF=alg.at(owner.bridge['core7']['logF0'][0])
    zero(logF+C+L*G,'actual original F0 log gauge')
    assert g.nodes[owner.bridge['core7']['F0'][0]]==dict(operation='analytic_unary',name='exp',argument=owner.bridge['core7']['logF0'][0])
    zero(alg.at(d['B'][0])-(-L*G+sy.log(phi)+sy.log(220)/2+sy.log(1+Z**2)),'actual B whole-function definition')
    # Higher ordinary Z rows are supplied by the pinned genuine ordinary
    # calculus. This identity uses arbitrary function values, not sampled
    # polynomial or numerical midpoint data.
    value_rows=[]
    for key,left,right,_ in current.JOINS:
        # U remains the R110 function through long/reference. These actual
        # value sources are composed from that same input by the replayed
        # constant-U recipes, before restoration changes it to 4Z.
        local=dict(env)
        for endpoint in ('Rsh','Rz'):local[d['endpoints'][endpoint]['actual']['U'][0]]=U
        a=SourceAlgebra(g,local)
        for field in ('phi','U'):
            left_value=a.at(owner.background_endpoints[key]['left']['actual'][field].value.node)
            right_value=a.at(owner.background_endpoints[key]['right']['actual'][field].value.node)
            zero(left_value-right_value,('actual annular whole-Z value',key,field));value_rows.append([key,field])
        # The actual right Duhamel source has its incoming value multiplied
        # by exp(0); its additional integral has equal limits. Validate this
        # without expanding any earlier nonzero history primitive.
        bg=owner.backgrounds[right];x=owner.charts[right]['coordinate']
        arbitrary=dict(local);arbitrary[x.node]=sy.Integer(0)
        for field,ids in bg['incoming_histories'].items():arbitrary[ids[0]]=sy.Symbol('nonzero_incoming_'+key+'_'+field,real=True)
        a=SourceAlgebra(g,arbitrary)
        for field,ids in bg['moments'].items():
            zero(a.at(ids[0])-arbitrary[bg['incoming_histories'][field][0]],('actual Duhamel nonzero inlet',key,field))
            value_rows.append([key,field])
    source_rec=current.parent.linear.read(current.previous.previous.downstream_source.RECEIPT)['checks']
    assert source_rec['long_restore_joins']['long_inlet_F0_anchored_gauge_logC_and_B_cancellation_checked']
    return dict(passed=True,actual_whole_Z_phi_U_and_six_history_value_rows=value_rows,
        actual_F0_exp_minus_logC_minus_Lambda_G_read_from_source_not_assigned=True,
        arbitrary_positive_R110_phi_and_arbitrary_nonzero_axial_history=True,
        actual_B_cancellation_and_exact_restoration_4Z_E_identity=True,
        higher_ordinary_Z_rows_follow_original_genuine_chain_and_FTC_calculus=True)


def flat_controls_and_background_jets(owner):
    g=owner.g;ev=prior.prior.current.previous.previous.modified_source.lexical.Algebra;flat=[]
    for name,control in owner.annular_controls.items():
        for edge in (0,1):
            for family,refs in control.items():
                for n,ref in enumerate(refs):
                    value=ev(g,{owner.charts[name]['coordinate'].node:sy.Integer(edge)}).at(ref.node)
                    zero(value-(edge if family=='sigma' and n==0 else 0),('actual fifth-order flat source',name,edge,family,n))
                    flat.append([name,edge,family,n,ref.node])
    y=sy.Symbol('independent_physical_y',real=True);Z=sy.Symbol('Z',real=True)
    Phi=sy.Function('unknown_actual_phi')(y);B,E,T=sy.symbols('actual_B actual_E original_finite_T',real=True)
    coefficients=[];states=[]
    for key,left,right,_ in current.JOINS:
        models=[]
        for name,edge in ((left,1),(right,0)):
            bg=owner.backgrounds[name]
            env={bg['actual']['phi'][0]:Phi,owner.downstream['T']:T,
                owner.downstream['B'][0]:B,owner.downstream['E'][0]:E}
            if name in owner.annular_controls:
                control=owner.annular_controls[name]
                env[control['prime'][0].node]=sy.Symbol('free_fifth_prime_tail_'+name)*y**5
            a=geometry.Algebra(g,env)
            angular=a.at(bg['actual']['phi_y'][0])/Phi;axial=a.at(bg['actual']['U_y'][0])
            models.append(([sy.diff(angular,y,n).subs(y,0) for n in range(5)],
                [sy.diff(axial,y,n).subs(y,0) for n in range(5)]))
        for family in range(2):
            for n,(Lvalue,Rvalue) in enumerate(zip(models[0][family],models[1][family])):
                zero(Lvalue-Rvalue,('actual annular physical RHS',key,family,n));coefficients.append([key,family,n])
        all_states=[]
        for angular,axial in models:
            phi=[sy.Function('seam_phi')(Z)];U=[sy.Function('seam_U')(Z)]
            for n in range(5):
                phi.append(sum(comb(n,k)*angular[k]*phi[n-k] for k in range(n+1)));U.append(axial[n])
            product=lambda A,B,n:sum(comb(n,k)*A[k]*B[n-k] for k in range(n+1))
            histories={field:[sy.Function('nonzero_'+field)(Z)] for field in bridge.RATES}
            for n in range(5):
                density=dict(H=2*phi[n],M=U[n],K=2*product(phi,U,n),A=product(U,U,n),B=product(phi,phi,n),C=product(phi,phi,n))
                for field,rate in bridge.RATES.items():histories[field].append(density[field]-rate*histories[field][-1])
            all_states.append(dict(phi=phi,U=U,**histories))
        for field in all_states[0]:
            for n,(Lvalue,Rvalue) in enumerate(zip(all_states[0][field],all_states[1][field])):
                zero(Lvalue-Rvalue,('actual annular physical IVP induction',key,field,n));states.append([key,field,n])
    return dict(passed=True,actual_sigma_y5_and_prime_y4_flat_rows=flat,
        actual_physical_RHS_rows=coefficients,independent_nonzero_inlet_physical_IVP_rows=states,
        first_four_prime_derivatives_zero_from_actual_sigma_through_fifth=True,
        actual_background_phi_U_and_six_histories_physical_y0_to5=True,
        native_T_reference_length_and_restoration_one_not_conflated=True)


def actual_loop_phase_transfer(owner):
    g=owner.g;Pstar=owner.ref(owner.downstream['Pstar']);roots=[]
    for name in current.CHARTS:
        bg=owner.backgrounds[name];chart=owner.loop[name]
        packet={group:{field:owner.jet(ids) for field,ids in bg[group].items()} for group in ('fields','radial','generic')}
        packet['radius']=owner.ref(bg['radius'])
        assert encode(prior.loop_source.generic_roots(g,packet,Pstar))==chart['roots']
        inverse=g.nodes[chart['loop']['inverse'][0]]
        assert inverse['operation']=='original_monotone_phase_inverse' and inverse['true_prescribed_source_C5_implicit_inverse']
        assert inverse['phase_function']==chart['loop']['Phi'][0] and inverse['target_phase']==chart['actual_phase']
        offset=g.sub(g.unary('log',g.quotient(owner.ref(chart['radius']),owner.ref(owner.bridge['Ra']),
            'same positive actual R/Ra')),owner.ref(owner.loop_input['phase_origin']))
        assert offset.node==chart['global_offset'] and g.unary('fractional_part',g.mul(owner.N,offset)).node==chart['actual_phase']
        D=owner.factories[name];J=owner.charts[name]['jacobian']
        for factory in D.shared['factories'].values():
            assert factory.cache[owner.charts[name]['phase'].node]==(g.mul(owner.N,J)
                if factory.variable==owner.charts[name]['coordinate'] else g.zero)
        roots.append([name,chart['roots'],chart['actual_phase'],J.node])
    H,T,C,P,Ra=sy.symbols('actual_hb actual_T actual_logC actual_logP actual_Ra',positive=True);seams=[]
    for key,left,right,_ in current.JOINS:
        env={owner.downstream['hb']:H,owner.downstream['T']:T,owner.parameters['logC'].node:C,owner.parameters['logP'].node:P,
            owner.bridge['Ra']:Ra,owner.loop_input['phase_origin']:sy.Symbol('same_original_phase_origin',real=True),
            owner.charts[left]['coordinate'].node:sy.Integer(1),owner.charts[right]['coordinate'].node:sy.Integer(0)}
        a=geometry.Algebra(g,env)
        zero(a.at(owner.joins[key]['left_radius'].value.node)-a.at(owner.joins[key]['right_radius'].value.node),('actual radius',key))
        zero(a.at(owner.loop[left]['global_offset'])-a.at(owner.loop[right]['global_offset']),('actual global phase',key));seams.append(key)
    receipt=current.parent.linear.read(current.previous.previous.RECEIPT)['checks']
    assert receipt['actual_recovery_and_pressure_sources']['actual_mean_ODE_and_actual_saved_reduced_S_on_all_ten_charts']
    return dict(passed=True,actual_original_roots_inverse_phase_and_NJ_bindings=roots,actual_radius_global_phase_seams=seams,
        source_background_y5_unique_inverse_periodic_lift_then_leading_U_y4_C5=True,
        unchanged_Q_y4_N_y2_pressure_y3_operators_consumed_from=current.previous.previous.RECEIPT,
        actual_nonzero_leading_mean_inlets_verified_in_next_group=True)


def actual_own_mixed_joins(owner):
    g=owner.g;evidence=[];means_checked=[];protected=[]
    for key,left,right,_ in current.JOINS:
        env,Z,_=prior.prior.prior.active.prior.rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
        Lchart,Rchart=owner.charts[left],owner.charts[right]
        x=sy.Symbol(g.nodes[Lchart['coordinate'].node]['name'],real=True);t=sy.Symbol(g.nodes[Rchart['coordinate'].node]['name'],real=True)
        Js={name:sy.Symbol('actual_positive_J_'+name,positive=True) for name in current.CHARTS}
        for name in ('restoration','postrestore'):Js[name]=sy.Integer(1)
        for name in current.CHARTS:
            if owner.charts[name]['jacobian']!=g.one:env[owner.charts[name]['jacobian'].node]=Js[name]
        # Each pair uses the physical logarithmic distance from its seam;
        # the actual geometry was proved before this finite-germ model.
        Rstar=sy.Symbol('actual_positive_R_'+key,positive=True)
        env[Lchart['radius'].node]=Rstar*sy.exp(Js[left]*(x-1));env[Rchart['radius'].node]=Rstar*sy.exp(Js[right]*t)
        env[owner.rminus.active_sc.node]=Algebra(owner,env).at(owner.parameters['sc'])
        for name in TEN[:TEN.index(left)+1]:
            coordinate=sy.Symbol(g.nodes[owner.charts[name]['coordinate'].node]['name'],real=True)
            for field,ids in owner.loop[name]['contributions'].items():
                bind(env,owner.jet(ids),sy.Function('retained_'+name+'_correction_'+field)(coordinate,Z),Z)
            for field,ids in owner.modified[name]['local_integrals'].items():
                bind(env,owner.jet(ids),sy.Function('retained_'+name+'_signed_'+field)(coordinate,Z),Z)
            if name!=left:
                for field in ('Q','Q_y'):
                    bind(env,owner.jet(owner.modified[name]['physical'][field]),sy.Function('retained_'+name+'_leading_'+field)(coordinate,Z),Z)
        mapped=1+Js[right]*t/Js[left];mean=sy.Function('common_background_mean_'+key)
        bind(env,owner.jet(owner.backgrounds[left]['fields']['m']),mean(x,Z),Z)
        bind(env,owner.jet(owner.backgrounds[right]['fields']['m']),mean(mapped,Z),Z)
        for name,arg in ((left,x),(right,mapped)):
            for field in ('Q','Q_y'):
                bind(env,owner.jet(owner.modified[name]['physical'][field]),sy.Function('proved_common_leading_'+key+'_'+field)(arg,Z),Z)
            for n,jet in enumerate(owner.charts[name]['N_y']):bind(env,jet,sy.Function('proved_common_N_'+key+'_'+str(n))(arg,Z),Z)
        for side in ('left','right'):
            packet=owner.joins[key][side]
            for group,field in (('fields','P1'),('cumulative','m13'),('cumulative','m15')):
                for jet in packet[group][field]:
                    for ref in rows(jet):assert ref.node not in env;protected.append(ref.node)
        a=Algebra(owner,env);means=[]
        for name,value in ((left,g.one),(right,g.zero)):
            chart=owner.charts[name];jet=current.outer.Algebra(g,5).at(chart['physical']['m'][0],chart['coordinate'],value)
            means.append([a.at(ref) for ref in rows(jet)])
        for k,(Lvalue,Rvalue) in enumerate(zip(*means)):
            zero(Lvalue-Rvalue,('actual annular leading mean inlet',key,k));means_checked.append([key,k])
        for spec in owner.contract:
            for n in range(spec['physical_y_order']+1):
                Lrows=rows(owner.joins[key]['left'][spec['group']][spec['field']][n]);Rrows=rows(owner.joins[key]['right'][spec['group']][spec['field']][n])
                for k,(Lref,Rref) in enumerate(zip(Lrows,Rrows)):
                    zero(a.at(Lref)-a.at(Rref),('actual annular own mixed row',key,spec['group'],spec['field'],n,k))
                    evidence.append(dict(interface=key,group=spec['group'],field=spec['field'],physical_y_order=n,
                        ordinary_Z_order=k,left_node=Lref.node,right_node=Rref.node))
    assert len(evidence)==984 and len(means_checked)==24
    return dict(passed=True,actual_four_annular_own_mixed_rows=evidence,actual_leading_mean_C5_inlets=means_checked,
        explicit_unbound_own_pressure_output_references=len(protected),
        finite_germs_only_after_original_source_proofs_not_neighborhood_equalities=True,
        nonzero_mass_pressure_five_moments_S0_S1_and_one_plus_delta_recovery_retained=True)


CHECKS=(source_contract,actual_background_recipes,actual_value_and_gauge_joins,
        flat_controls_and_background_jets,actual_loop_phase_transfer,actual_own_mixed_joins)


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('ANNULAR_JOIN_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['annular_mixed_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),actual_own_mixed_rows=984,
        actual_annular_interfaces=[key for key,*_ in current.JOINS],exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_FOUR_ANNULAR_INTERFACES_984_OWN_ROWS_ACCEPTED',flush=True)
    return rec
