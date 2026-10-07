"""Same-original active kappa correlations and mixed source-term budgets.

The accepted periodic functions and full inlet-to-Rc geometry are unchanged.
Active derivative formulas replace unrelated a/b extremes only in q support.
All bounds retain original source jets; caps never become field coefficients.
"""
import itertools
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_periodic_mixed_conditioning as preceding
from lei_ren_part1_paper_compliant_current_native_periodic_mixed_conditioning import (
    HERE,PREFIX,sha,previous,packets,ep,native,current,serial,prior,repair,
    LogUpper,SlowJet,ZERO,DY,DZ,DYZ,ORDERS,RATES,MIN_N,
    power,maximum,root_caps,cutoff_caps,active_velocity_correlations,
    velocity_caps,budget_components)

NAME=PREFIX+'current_native_active_kappa_mixed_conditioning.json'
RECEIPT=PREFIX+'current_native_active_kappa_mixed_conditioning_check.json'
GATE='current_original_native_active_kappa_mixed_source_conditioning_and_weighted_attribution_executed'


def exact_kernel_theorem():
    theorem=preceding.exact_kernel_theorem()
    y,Z=sy.symbols('y Z');a=sy.Function('a')(y,Z);b=sy.Function('b')(y,Z)
    t0=-b/a;Delta=a+b*b/a-2
    for i in (y,Z):
        formula=(1-t0*t0)*sy.diff(a,i)-2*t0*sy.diff(b,i)
        if sy.simplify(sy.diff(Delta,i)-formula)!=0:raise ArithmeticError('Active kappa first identity failed')
    formula=(1-t0*t0)*sy.diff(a,y,Z)-2*t0*sy.diff(b,y,Z)+2*(sy.diff(b,y)+t0*sy.diff(a,y))*(sy.diff(b,Z)+t0*sy.diff(a,Z))/a
    if sy.simplify(sy.diff(Delta,y,Z)-formula)!=0:raise ArithmeticError('Active kappa mixed identity failed')
    theorem['original_active_kappa_derivatives']=dict(passed=True,
        Delta='a+b^2/a-2',t0='-b/a',active_constraint='a+b^2/a<2+eta<=5/2',
        active_C0_bounds='a<=3, b^2<=3a, |t0|<=sqrt(3/a_min)',
        first='Delta_i=(1-t0^2)*a_i-2*t0*b_i',
        mixed='Delta_yZ=(1-t0^2)*a_yZ-2*t0*b_yZ+2*(b_y+t0*a_y)*(b_Z+t0*a_Z)/a',
        equivalent_mixed_cross='-2*t0_Z*(b_y+t0*a_y)=-2*t0_y*(b_Z+t0*a_Z)',
        only_modulation_support_bounds_intersected=True,all_original_derivatives_retained=True,
        flat_branch_q_jets_zero_without_zeroing_background=True)
    q,x,s=sy.symbols('q chi s');P=sy.Function('P');H=sy.Function('H');u=q*x
    pu=sy.Subs(sy.diff(P(s),s),s,u);puu=sy.Subs(sy.diff(P(s),s,2),s,u)
    hu=sy.Subs(sy.diff(H(s),s),s,u);huu=sy.Subs(sy.diff(H(s),s,2),s,u)
    tests=((q*P(u),((q,P(u)+u*pu),(x,q*q*pu),(q,q,x*(2*pu+u*puu)),(q,x,q*(2*pu+u*puu)),(x,x,q**3*puu))),
        (q*q*H(u),((q,q*(2*H(u)+u*hu)),(x,q**3*hu),(q,q,2*H(u)+4*u*hu+u*u*huu),(q,x,q*q*(3*hu+u*huu)),(x,x,q**4*huu))))
    for function,rows in tests:
        for row in rows:
            if sy.simplify(sy.diff(function,*row[:-1])-row[-1])!=0:
                raise ArithmeticError('Correlated qP/q2H parameter identity failed')
    theorem['weighted_partial_angle_kernel_bounds']=dict(passed=True,
        angle_first='theta_u=2*sin(theta)/h',
        angle_mixed='theta_uu=4*sin(theta)*cos(theta)/h^2-2*u*sin(theta)/h^3',
        outside_signed_r_abs_ge_half=dict(P='2pi/h',P_u='(2+4pi)/h^2',P_uu='(14+16pi)/h^3',
            H_u='(14+15pi)/h',H_uu='(160+132pi)/h^2'),
        inside_fourier_conversion='h<=2/sqrt(3): previous P_u<=20pi,P_uu<=100pi,H_u<=50pi,H_uu<=600pi',
        uniform_all_finite_signed_u_and_partial_psi=dict(P='4pi/h',P_u='(80pi/3)/h^2',P_uu='160pi/h^3',H='pi',H_u='60pi/h',H_uu='800pi/h^2'),
        no_positive_q_or_p2_lower_required=True,
        exact_functions=dict(F='q*P(chi*q,psi)',G='q^2*H(chi*q,psi)',chi='p2/dstar'),
        F_partials=dict(q='P+u*P_u',chi='q^2*P_u',qq='chi*(2*P_u+u*P_uu)',qchi='q*(2*P_u+u*P_uu)',chichi='q^3*P_uu'),
        G_partials=dict(q='q*(2*H+u*H_u)',chi='q^3*H_u',qq='2*H+4*u*H_u+u^2*H_uu',qchi='q^2*(3*H_u+u*H_uu)',chichi='q^4*H_uu'),
        correlation_bounds='|u|/h^2<=1/2, |u|/h<=1; derivatives of chi remain original source derivatives',
        bounds_not_point_function_values=True)
    return theorem


def active_kappa_caps(c,a,b,t0_original,Delta,logamin):
    C=lambda v:LogUpper.constant(c,v);add=lambda *v:LogUpper.add(c,v)
    inva=LogUpper(c,-logamin)
    t0=previous.minimum(t0_original[ZERO],b[ZERO]*inva,LogUpper(c,(c.ln(3)-logamin)/2))
    oneplus=add(C(1),t0*t0);rows={ZERO:previous.minimum(Delta[ZERO],C(2))};pieces={}
    numerators={i:add(b[i],t0*a[i]) for i in (DY,DZ)}
    for i in (DY,DZ):
        terms=dict(a_source=oneplus*a[i],b_source=C(2)*t0*b[i])
        rows[i]=previous.minimum(Delta[i],add(*terms.values()))
        pieces['y%d_Z%d'%i]={key:value.record() for key,value in terms.items()}
    cross=previous.minimum(C(2)*inva*numerators[DY]*numerators[DZ],
        C(2)*t0_original[DZ]*numerators[DY],C(2)*t0_original[DY]*numerators[DZ])
    terms=dict(a_mixed_source=oneplus*a[DYZ],b_mixed_source=C(2)*t0*b[DYZ],correlated_mixed_cross=cross)
    rows[DYZ]=previous.minimum(Delta[DYZ],add(*terms.values()))
    pieces['y1_Z1']={key:value.record() for key,value in terms.items()}
    selected=SlowJet(c,rows)
    return selected,dict(unrestricted_original_Delta_caps=Delta.record(),active_same_function_Delta_caps=selected.record(),
        active_t0_C0_cap=t0.record(),source_term_caps=pieces,
        exact_identity_theorem='original_active_kappa_derivatives',active_support_only=True,
        derivative_rows_never_clipped_or_zeroed=True)


def primitive_component_caps(c,a,b,E,T1,hat,L,Tslow,curvature,M):
    C=lambda v:LogUpper.constant(c,v);twoPi=C(2*c.pi);fourPi=C(4*c.pi)
    inv2=LogUpper(c,-c.ln(2*c.pi));inv4=LogUpper(c,-c.ln(4*c.pi));half=C('.5')
    cross=dict(slow_y_times_LZ=Tslow[DY]*L[DZ],slow_Z_times_Ly=Tslow[DZ]*L[DY],curvature_times_Ly_LZ=curvature*L[DY]*L[DZ])
    A=dict(a_mixed=half*a[DYZ],a_y_inverse_Z=inv4*a[DY]*L[DZ],a_Z_inverse_y=inv4*a[DZ]*L[DY],a_inverse_mixed=inv4*a[ZERO]*L[DYZ])
    A.update({key:inv4*a[ZERO]*cap for key,cap in cross.items()})
    MM=dict(a_mixed_T1=inv2*a[DYZ]*hat[ZERO],a_y_T1_Z=inv2*a[DY]*hat[DZ],a_Z_T1_y=inv2*a[DZ]*hat[DY],
        a_T1_fixed_mixed=inv2*a[ZERO]*T1[DYZ],a_inverse_fixed_mixed=half*inv2*a[ZERO]*L[DYZ],b_mixed=b[DYZ])
    MM.update({key:inv2*a[ZERO]*cap for key,cap in cross.items()})
    B=dict(E_mixed_M=half*E[DYZ]*M[ZERO],E_y_M_Z=half*E[DY]*M[DZ],E_Z_M_y=half*E[DZ]*M[DY])
    B.update({key:half*E[ZERO]*cap for key,cap in MM.items()})
    return dict(A_yZ={key:cap.record() for key,cap in A.items()},B_over_Pstar_yZ={key:cap.record() for key,cap in B.items()},
        source_product_rules_not_derivatives_of_caps=True,component_sums_precede_old_certificate_intersection=True)


def leading_component_caps(jets):
    c=jets['E'].ctx;C=lambda v:LogUpper.constant(c,v)
    monomials=dict(m=[(1,('B',))],h=[(1,('E','A'))],k=[(1,('E','V','A')),(1,('E','B'))],
        e=[(2,('V','B')),(-1,('E','E','A'))],p=[(1,('E','E','A'))])
    groups={}
    for derivative,order in (('C0',DY),('Z',DYZ)):
        groups[derivative]={}
        for density,terms in monomials.items():
            pieces=[]
            for coefficient,factors in terms:
                for orders in itertools.product(ORDERS,repeat=len(factors)):
                    if tuple(sum(k[i] for k in orders) for i in range(2))!=order:continue
                    cap=C(abs(coefficient));source=[]
                    for factor,k in zip(factors,orders):
                        cap=cap*jets[factor][k];source.append(dict(source=factor,ordinary_y_order=k[0],ordinary_Z_order=k[1]))
                    pieces.append(dict(cap=cap,coefficient=abs(coefficient),sign=1 if coefficient>0 else -1,factors=source))
            groups[derivative][density]=pieces
    return groups


def encode_leading_terms(groups):
    return {derivative:{density:[dict(cap=term['cap'].record(),coefficient=term['coefficient'],sign=term['sign'],factors=term['factors']) for term in terms]
        for density,terms in rows.items()} for derivative,rows in groups.items()}


def joint_weighted_kernel_caps(c,q,chi,P,H):
    C=lambda v:LogUpper.constant(c,v);add=lambda *v:LogUpper.add(c,v)
    CP0,CP,CPP,CH,CHH=(C(4*c.pi),C(80*c.pi/3),C(160*c.pi),C(60*c.pi),C(800*c.pi))
    oldF=q*P;oldG=q*q*H
    F={ZERO:oldF[ZERO]};G={ZERO:oldG[ZERO]}
    for i in (DY,DZ):
        F[i]=add(add(CP0,CP*C('.5'))*q[i],CP*q[ZERO]*q[ZERO]*chi[i])
        G[i]=add(add(C(2*c.pi),CH)*q[ZERO]*q[i],CH*power(q[ZERO],3)*chi[i])
    mixed=add(q[DY]*chi[DZ],q[DZ]*chi[DY])
    F[DYZ]=add(add(CP0,CP*C('.5'))*q[DYZ],CP*power(q[ZERO],2)*chi[DYZ],
        chi[ZERO]*add(C(2)*CP,CPP)*q[DY]*q[DZ],q[ZERO]*add(C(2)*CP,CPP)*mixed,
        CPP*power(q[ZERO],3)*chi[DY]*chi[DZ])
    G[DYZ]=add(add(C(2*c.pi),CH)*q[ZERO]*q[DYZ],CH*power(q[ZERO],3)*chi[DYZ],
        add(C(2*c.pi),C(4)*CH,CHH)*q[DY]*q[DZ],add(C(3)*CH,CHH)*power(q[ZERO],2)*mixed,
        CHH*power(q[ZERO],4)*chi[DY]*chi[DZ])
    f=SlowJet(c,{order:previous.minimum(cap,oldF[order]) for order,cap in F.items()})
    g=SlowJet(c,{order:previous.minimum(cap,oldG[order]) for order,cap in G.items()})
    return f,g,dict(original_chi_slow_caps=chi.record(),qP_slow_caps=f.record(),q2H_slow_caps=g.record(),
        previous_separate_qP_caps=oldF.record(),previous_separate_q2H_caps=oldG.record(),
        same_original_function_cap_intersections=True,partial_angle_and_signed_chi_uniform=True,
        no_q_or_chi_division_used=True)


def primitive_caps(c,roots,positive,eta_log,dstar_log,sigma_caps,primitive):
    C=lambda v:LogUpper.constant(c,v);L=lambda v:LogUpper(c,v);add=lambda *v:LogUpper.add(c,v)
    logamin=packets.interval(c,positive['log_actual_a_positive_lower']);logeta=c.mpf(eta_log);logd=c.mpf(dstar_log)
    if ep(logeta)[1]>ep(c.ln(c.mpf('.5')))[0]:raise ValueError('Original eta<=1/2 required')
    zero=SlowJet.constant(c,0)
    cutoff=serial.p.whole_cutoff_C1(roots,logeta,logamin)
    if cutoff['q'].zero:
        return dict(A=zero,B=zero,record=dict(original_whole_support_flat=True,all_slow_primitive_jets_exact_zero=True))
    a,b,p2,E,Delta=(root_caps(c,roots,k) for k in ('a','b','p2','E','kappa_minus2'))
    velocity_proof=active_velocity_correlations(c,dict(E=E),a,b)
    a.rows[ZERO]=previous.minimum(a[ZERO],C(3));b.rows[ZERO]=previous.minimum(b[ZERO],C('1.5'))
    Delta,active_proof=active_kappa_caps(c,a,b,root_caps(c,roots,'t0'),Delta,logamin)
    q,qproof=cutoff_caps(c,a,Delta,logamin,logeta,sigma_caps,cutoff)
    inva=L(-logamin);t0=SlowJet(c,{k:C(0) for k in ORDERS})
    t0.rows[ZERO]=previous.minimum(b[ZERO]*inva,L((c.ln(3)-logamin)/2),root_caps(c,roots,'t0')[ZERO])
    for i in (DY,DZ):t0.rows[i]=add(b[i],t0[ZERO]*a[i])*inva
    t0.rows[DYZ]=add(b[DYZ],t0[DY]*a[DZ],t0[DZ]*a[DY],t0[ZERO]*a[DYZ])*inva
    original_t0=root_caps(c,roots,'t0')
    for order in ORDERS:t0.rows[order]=previous.minimum(t0[order],original_t0[order])
    u=p2*q;u=SlowJet(c,{k:v.divide_positive(logd) for k,v in u.rows.items()})
    P=SlowJet(c,{ZERO:C(4*c.pi),DY:C(20*c.pi)*u[DY],DZ:C(20*c.pi)*u[DZ],
        DYZ:add(C(20*c.pi)*u[DYZ],C(100*c.pi)*u[DY]*u[DZ])})
    H=SlowJet(c,{ZERO:C(c.pi),DY:C(50*c.pi)*u[DY],DZ:C(50*c.pi)*u[DZ],
        DYZ:add(C(50*c.pi)*u[DYZ],C(600*c.pi)*u[DY]*u[DZ])})
    two,four=SlowJet.constant(c,2),SlowJet.constant(c,4)
    chiJet=SlowJet(c,{k:v.divide_positive(logd) for k,v in p2.rows.items()})
    qP,q2H,weighted_kernel_proof=joint_weighted_kernel_caps(c,q,chiJet,P,H)
    T1=t0*SlowJet.constant(c,2*c.pi)+two*qP
    T2=t0*t0*SlowJet.constant(c,2*c.pi)+four*t0*qP+four*q2H
    nu=SlowJet.constant(c,1)+t0*t0+two*q*q
    numerator={k:add(C(2*c.pi)*nu[k],T2[k]) for k in (DY,DZ,DYZ)}
    hcap=add(C(1),u[ZERO]);Tslow={i:add(t0[i],C(4)*q[i]*hcap,C(12)*q[ZERO]*u[i]) for i in (DY,DZ)}
    chi=p2[ZERO].divide_positive(logd)
    curvature=previous.minimum(C(8)*q[ZERO]*power(hcap,3),
        C(4)*q[ZERO]*power(add(C(1),chi*add(C(1),t0[ZERO])),1.5))
    cross=add(Tslow[DY]*numerator[DZ],Tslow[DZ]*numerator[DY],curvature*numerator[DY]*numerator[DZ])
    psiYZ=add(numerator[DYZ],cross)
    hat={ZERO:T1[ZERO],DY:add(T1[DY],C('.5')*numerator[DY]),DZ:add(T1[DZ],C('.5')*numerator[DZ]),
        DYZ:add(T1[DYZ],cross,C('.5')*numerator[DYZ])}
    A={ZERO:previous.read_cap(c,current.magnitude(primitive['values']['A']).record())}
    for i in (DY,DZ):A[i]=add(C('.5')*a[i],a[ZERO]*numerator[i]*C(1/(4*c.pi)))
    A[DYZ]=add(C('.5')*a[DYZ],C(1/(4*c.pi))*add(a[DY]*numerator[DZ],a[DZ]*numerator[DY],a[ZERO]*psiYZ))
    A[DZ]=previous.minimum(A[DZ],previous.read_cap(c,current.magnitude(primitive['values']['A_Z']).record()))
    M0=C(200);M={ZERO:M0}
    for i in (DY,DZ):
        naive=add(C(1/(2*c.pi))*add(a[i]*hat[ZERO],a[ZERO]*hat[i]),b[i])
        correlated=add(C(2)*t0[i]*A[ZERO],C(2)*t0[ZERO]*A[i],
            C(4)*add(a[i]*q[ZERO],a[ZERO]*q[i]),C(20)*a[ZERO]*q[ZERO]*u[i],
            a[ZERO]*add(C(1),t0[ZERO])*numerator[i]*C(1/(2*c.pi)))
        M[i]=previous.minimum(naive,correlated)
    M[DYZ]=add(C(1/(2*c.pi))*add(a[DYZ]*hat[ZERO],a[DY]*hat[DZ],a[DZ]*hat[DY],a[ZERO]*hat[DYZ]),b[DYZ])
    B=(E*SlowJet(c,M))*SlowJet.constant(c,'.5')
    B.rows[ZERO]=previous.minimum(B[ZERO],previous.read_cap(c,current.magnitude(primitive['values']['B_over_Pstar']).record()))
    B.rows[DZ]=previous.minimum(B[DZ],previous.read_cap(c,current.magnitude(primitive['values']['B_Z_over_Pstar']).record()))
    A=SlowJet(c,A)
    return dict(A=A,B=B,record=dict(original_whole_support_flat=False,native_source_derivative_orders=[list(k) for k in ORDERS],
        normalized_native_a_b_p2_E_Delta_caps={key:root_caps(c,roots,key).record() for key in ('a','b','p2','E','kappa_minus2')},
        active_only_a_and_b_C0_constraints=dict(a='0<a<=3',b='|b|<=3/2',t0='|t0|<=sqrt(3/a_min)'),
        native_active_kappa_correlation=active_proof,original_cutoff_slow_caps=qproof,active_t0_caps=t0.record(),u_caps=u.record(),
        joint_weighted_kernel_proof=weighted_kernel_proof,
        fixed_free_angle_P_caps=P.record(),fixed_free_angle_H_caps=H.record(),
        original_fixed_angle_T1_caps=T1.record(),original_fixed_angle_T2_caps=T2.record(),
        inverse_first_and_mixed_numerator_caps={str(k):v.record() for k,v in numerator.items()},
        normalized_point_curvature_cap=curvature.record(),inverse_yZ_cap=psiYZ.record(),
        correlated_composed_T1_caps={str(k):v.record() for k,v in hat.items()},
        normalized_A_slow_caps=A.record(),normalized_B_over_Pstar_slow_caps=B.record(),active_support_E_correlations=velocity_proof,
        inverse_mixed_components=dict(fixed_angle_mixed=numerator[DYZ].record(),slow_y_times_LZ=(Tslow[DY]*numerator[DZ]).record(),slow_Z_times_Ly=(Tslow[DZ]*numerator[DY]).record(),curvature_times_Ly_LZ=(curvature*numerator[DY]*numerator[DZ]).record()),
        mixed_primitive_components=primitive_component_caps(c,a,b,E,T1,hat,numerator,Tslow,curvature,M),
        all_source_derivatives_retained=True,active_C0_constraints_do_not_change_source_functions=True,
        no_positive_source_scale_or_Poisson_denominator_materialized=True,higher_mixed_jets_not_admitted=True))


class NativeActiveKappaMixedConditioning:
    def __init__(self,owner):
        if type(owner) is not previous.NativeSignedAveraging:raise ValueError('Same original signed averaging owner required')
        checked=json.loads((HERE/preceding.RECEIPT).read_bytes());self.old=json.loads((HERE/preceding.NAME).read_bytes())
        if not checked['all_passed'] or not checked[preceding.GATE] or checked['source_family']!=owner.family:
            raise ValueError('Checked same-original signed averaging route required')
        self.owner=owner;self.target=owner.owner;self.ctx=owner.ctx;self.family=owner.family;self.coordinates=owner.coordinates;self.service=owner.service
        self.service.bind_hashes(checked['input_hashes']);self.service.bind_hashes({preceding.NAME:sha(preceding.NAME),preceding.RECEIPT:sha(preceding.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)})
        manifest=json.loads((HERE/(PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes())
        self.sigma={int(k):previous.read_cap(self.ctx,v) for k,v in manifest['sigma_global_derivative_log_caps'].items()}
        self.theorem=exact_kernel_theorem()

    @native.inlet.source_precision
    def route(self,Z=(-1,1)):
        c=self.ctx;coords=self.coordinates;Z=c.mpf(Z)
        old=next((row for row in self.old['actual_native_periodic_mixed_conditioning_records'].values() if packets.interval(c,row['Z_box'])._mpi_==Z._mpi_),None)
        if old is None:raise ValueError('This replay adapter binds either original accepted Z domain')
        root_owner=self.target.q_owner.owner.owner
        eta=packets.interval(c,root_owner.scales['selected_positive_eta_log']);dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        zero=lambda:{key:coords.scalar(0) for key in RATES}
        running,runningZ=zero(),zero();cells=[]
        for (label,chart,left,right),oldcell in zip(current.ROUTE,old['cells']):
            if (oldcell['label'],oldcell['chart'])!=(label,chart):raise ValueError('Same original ordered cell/function ledger required')
            geom=self.target.geometry(label,chart,left,right);factors={key:current.transfer.true_width_kernel(coords,geom,rate) for key,rate in RATES.items()}
            details={};components={};leading_terms={}
            if label=='initial_flat_collar':values,jets=zero(),zero();source=dict(original_full_initial_flat_receipt=previous.RECEIPT,all_increment_jets_exact_zero=True)
            else:
                source_query=self.target.q_owner.query(chart,Z,geom['coordinate']);roots=source_query['source']['roots'];packet=source_query['source']['packet']
                positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
                primitive=serial.whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
                mixed=primitive_caps(c,roots,positive,eta,dstar,self.sigma,primitive)
                vcaps=velocity_caps(c,packet)
                velocity_proof=active_velocity_correlations(c,vcaps,root_caps(c,roots,'a'),root_caps(c,roots,'b'))
                oldloop=self.owner.loop[chart]
                for key in ('A','B'):
                    for order in ORDERS:mixed[key].rows[order]=previous.minimum(mixed[key][order],previous.read_cap(c,oldloop['slow_phase_held_log_bounds'][key]['y%d_Z%d'%order]))
                loop=dict(slow_phase_held_log_bounds={key:mixed[key].record() for key in ('A','B')})
                bounded=previous.leading_caps(c,{key:value.record() for key,value in vcaps.items()},loop,primitive['values'])
                leading_terms=leading_component_caps(bounded)
                remainder=previous.remainder_caps(bounded);values={};jets={}
                from lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets_check import record_value
                for key in RATES:
                    components[key]=budget_components(bounded['leading'][key],remainder,factors[key],key)
                    cap,proof=previous.ibp_contribution(bounded['leading'][key],remainder['values'][key],factors[key],key)
                    capZ,proofZ=previous.ibp_contribution(bounded['leading'][key],remainder['Z_derivatives'][key],factors[key],key,True)
                    oldcap=current.magnitude(record_value(oldcell['source_C0_Nminus2_covers'][key],coords))
                    oldcapZ=current.magnitude(record_value(oldcell['source_Z_Nminus2_covers'][key],coords))
                    chosen,chosenZ=previous.minimum(cap,oldcap),previous.minimum(capZ,oldcapZ)
                    values[key]=previous.symmetric_cap(chosen,coords);jets[key]=previous.symmetric_cap(chosenZ,coords)
                    details[key]=dict(C0=proof,Z=proofZ,separate_native_budget_components={label:{derivative:v.record() for derivative,v in row.items()} for label,row in components[key].items()},
                        native_C0_coefficient_cap=cap.record(),native_Z_coefficient_cap=capZ.record(),
                        prior_C0_coefficient_cap=oldcap.record(),prior_Z_coefficient_cap=oldcapZ.record(),
                        same_original_uniform_coefficient_cover_intersection=True)
                source=dict(native_provenance=packet.provenance,native_normalized_active_velocity_correlations=velocity_proof,
                    original_periodic_C0_Z_proof=primitive['record'],native_periodic_mixed_proof=mixed['record'],
                    selected_mixed_caps={key:mixed[key].record() for key in ('A','B')},
                    leading_density_caps={key:value.record() for key,value in bounded['leading'].items()},
                    leading_density_slow_term_caps=encode_leading_terms(leading_terms),
                    exact_remainder_contract=remainder['record'],same_original_source_graph_views=previous.VIEWS,
                    formal_mixed_source_definitions_unchanged=True,source_point_values_not_selected_from_caps=True)
            incoming,incomingZ=running,runningZ
            running={key:factors[key]['decay']*incoming[key]+values[key] for key in RATES}
            runningZ={key:factors[key]['decay']*incomingZ[key]+jets[key] for key in RATES}
            record=dict(label=label,chart=chart,original_geometry=geom['record'],original_source=source,per_density_budget=details,
                source_C0_Nminus2_covers={key:v.record() for key,v in values.items()},source_Z_Nminus2_covers={key:v.record() for key,v in jets.items()},
                inherited_C0_Nminus2_covers={key:v.record() for key,v in incoming.items()},inherited_Z_Nminus2_covers={key:v.record() for key,v in incomingZ.items()},
                right_C0_Nminus2_covers={key:v.record() for key,v in running.items()},right_Z_Nminus2_covers={key:v.record() for key,v in runningZ.items()},
                true_mass_once_and_all_endpoints_preserved=True,source_increment_functions_unchanged=True)
            cells.append(dict(record=record,values=values,Z_derivatives=jets,factors=factors,incoming=incoming,incoming_Z=incomingZ,cumulative=running,cumulative_Z=runningZ,budget_components=components,leading_terms=leading_terms))
        endpoint=2/self.target.transfer.geometry.binder.fixed['Tw'];amplitude_query=self.target.q_owner.query('O3_power',Z,endpoint)
        roots=amplitude_query['source']['roots'];logA=packets.interval(c,self.target.owner.reservation['positive_Ac_over_S_log_lower'])
        A=coords.rebase(roots['E'][ZERO],self.family).positive_intersection(logA);AZ=coords.rebase(roots['E'][DZ],self.family)
        mu=c.mpf(packets.interval(c,self.target.owner.domain['right_collar_mu']));logmu=c.mpf(ep(c.ln(mu))[0]);mu_source=coords.scalar(mu)
        if mu._mpi_!=self.target.repair_mu._mpi_:raise ValueError('Same original repair mu required')
        wrap=lambda rows:{key:{-1:coords.scalar(0),-2:value} for key,value in rows.items()}
        targets=current.target_rows(wrap(running),wrap(runningZ),A,AZ,logA,mu_source,logmu);caps=current.target_caps(targets)
        matrix=self.target.repair_manifest['fresh_divided_linear_inverse_and_enclosures'];W=self.target.repair_manifest['fresh_exact_weight_definitions_and_enclosures']
        conditions=repair.contraction_log_conditions(c,caps,matrix,W,logmu,c.ln(MIN_N))
        comparison={}
        for key in repair.ROWS:
            new=previous.read_cap(c,caps['transformed_N_scaled_target_C1_caps'][key]);oldcap=previous.read_cap(c,old['actual_uniform_N_scaled_repair_C1_caps']['transformed_N_scaled_target_C1_caps'][key])
            comparison[key]=dict(new_N_scaled_C1_cap_at_floor=new.record(),prior_averaging_N_scaled_C1_cap_at_floor=oldcap.record(),
                new_certificate_strictly_sharper=(new.log is None or (oldcap.log is not None and ep(new.log)[1]<ep(oldcap.log)[0])),
                logarithmic_upper_cap_reduction=None if new.log is None or oldcap.log is None else oldcap.log-new.log,
                upper_bound_comparison_not_actual_error_measurement=True)
        suffix={key:coords.scalar(1) for key in RATES};weighted=[];dominant={};first_budget={};first_dominant={};factor_budgets=[];factor_dominant={}
        for cell in reversed(cells):
            v=wrap({key:cell['values'][key]*suffix[key] for key in RATES});j=wrap({key:cell['Z_derivatives'][key]*suffix[key] for key in RATES})
            part=current.target_caps(current.target_rows(v,j,A,AZ,logA,mu_source,logmu))
            weighted.append(dict(label=cell['record']['label'],C1_target_caps=part['transformed_N_scaled_target_C1_caps']))
            for key,row in part['transformed_N_scaled_target_C1_caps'].items():
                if row['exact_zero']:continue
                upper=ep(packets.interval(c,row['log_absolute_upper']))[1]
                if key not in dominant or upper>dominant[key][0]:dominant[key]=(upper,cell['record']['label'])
            if cell['record']['label']=='active_first_bridge':
                for component in next(iter(cell['budget_components'].values())):
                    v=wrap({key:previous.symmetric_cap(cell['budget_components'][key][component]['C0'],coords)*suffix[key] for key in RATES})
                    j=wrap({key:previous.symmetric_cap(cell['budget_components'][key][component]['Z'],coords)*suffix[key] for key in RATES})
                    detail=current.target_caps(current.target_rows(v,j,A,AZ,logA,mu_source,logmu))
                    first_budget[component]=detail['transformed_N_scaled_target_C1_caps']
                    for key,row in detail['transformed_N_scaled_target_C1_caps'].items():
                        if row['exact_zero']:continue
                        upper=ep(packets.interval(c,row['log_absolute_upper']))[1]
                        if key not in first_dominant or upper>first_dominant[key][0]:first_dominant[key]=(upper,component)
                for derivative,groups in cell['leading_terms'].items():
                    for density,terms in groups.items():
                        for index,term in enumerate(terms):
                            cap=current.magnitude(cell['factors'][density]['mass'])*LogUpper.constant(c,'.5')*term['cap']
                            v={key:coords.scalar(0) for key in RATES};j=dict(v)
                            (v if derivative=='C0' else j)[density]=previous.symmetric_cap(cap,coords)*suffix[density]
                            detail=current.target_caps(current.target_rows(wrap(v),wrap(j),A,AZ,logA,mu_source,logmu))
                            label=derivative+':'+density+':'+str(index)
                            factor_budgets.append(dict(label=label,slow_derivative=derivative,density=density,
                                signed_source_coefficient=term['sign']*term['coefficient'],factors=term['factors'],
                                weighted_target_caps=detail['transformed_N_scaled_target_C1_caps']))
                            for key,row in detail['transformed_N_scaled_target_C1_caps'].items():
                                if row['exact_zero']:continue
                                upper=ep(packets.interval(c,row['log_absolute_upper']))[1]
                                if key not in factor_dominant or upper>factor_dominant[key][0]:factor_dominant[key]=(upper,label)
            suffix={key:cell['factors'][key]['decay']*suffix[key] for key in RATES}
        record=dict(source_family=self.family,Z_box=Z,all_integer_N_lower=MIN_N,cells=[row['record'] for row in cells],original_true_cells=24,original_charts=17,
            actual_Rc_C0_Nminus2_covers={key:v.record() for key,v in running.items()},actual_Rc_Z_Nminus2_covers={key:v.record() for key,v in runningZ.items()},
            target_C0_orders=current.records(targets['values']),target_Z_orders=current.records(targets['Z_derivatives']),
            actual_uniform_N_scaled_repair_C1_caps=caps,conditional_repair_log_conditions=conditions,
            original_Rc_amplitude=dict(A=A.record(),A_Z=AZ.record(),positive_log_lower=logA,mu=mu,source_provenance=amplitude_query['source']['packet'].provenance),
            same_source_averaging_bound_comparison=comparison,weighted_cell_target_caps=list(reversed(weighted)),
            dominant_target_bound_cells={key:label for key,(_,label) in dominant.items()},
            first_bridge_native_weighted_target_budget_components=first_budget,
            first_bridge_dominant_budget_components={key:label for key,(_,label) in first_dominant.items()},
            first_bridge_weighted_slow_term_budgets=factor_budgets,first_bridge_dominant_slow_terms={key:label for key,(_,label) in factor_dominant.items()},
            slow_term_budgets_are_algebraic_absolute_budget_contributions_not_separate_control_functions=True,
            component_budgets_before_prior_whole_certificate_intersection=True,
            original_all_N_zero_inlet_and_phase_and_P0_retained=True,native_normalization_and_ordinary_width_applied_once=True,
            current_small_actual_target_or_control_or_terminal_closure_certified=False,point_inverse_or_higher_jets_installed=False,
            one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,targets=targets,caps=caps,values=wrap(running),Z_derivatives=wrap(runningZ))


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        target=current.NativeRcParameterTargets(current.preceding.NativeRcC1Histories(current.preceding.preceding.NativeO2C1Histories(current.preceding.preceding.make_middle_owner(bridge))))
        owner=NativeActiveKappaMixedConditioning(previous.NativeSignedAveraging(target));regions={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            regions[name]=owner.route(Z)['record'];print('Native active-kappa mixed conditioning:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_Z_query_count=2,
        actual_native_active_kappa_mixed_conditioning_records=regions,exact_correlated_periodic_kernel_theorem=owner.theorem,
        original_signed_mixed_density_source_graph_views=previous.VIEWS,
        native_C0_y_Z_yZ_caps_and_source_integral_descendants_rebuilt=True,
        actual_controls_or_terminal_closure_installed=False,point_inverse_or_higher_jets_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Same original full24-cell/17-chart route on two Z domains with support-only exact kappa derivative correlations and weighted leading slow-term attribution, with native normalized velocity/input mixed rows, branchwise correlated q jets, uniform signed-r integrated kernel derivatives, normalized inverse curvature and fixed-phi A/B y/yZ caps. Original zero-mean IBP/remainder/endpoints/memory feed N^-2 Rc targets and exact conditional repair recipe. Caps are function covers; no controls, terminal closure, point inverse/higher jets, global N/cone/energy/recursion/full NS.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
