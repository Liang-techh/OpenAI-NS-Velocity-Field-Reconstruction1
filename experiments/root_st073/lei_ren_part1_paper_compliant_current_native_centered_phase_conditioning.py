"""Centered original phase residual and a*nu=v mixed derivative bounds."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_active_kappa_mixed_conditioning as preceding
from lei_ren_part1_paper_compliant_current_native_active_kappa_mixed_conditioning import (
    HERE,PREFIX,sha,previous,packets,ep,native,current,serial,prior,repair,
    LogUpper,SlowJet,ZERO,DY,DZ,DYZ,ORDERS,RATES,MIN_N,
    power,maximum,root_caps,cutoff_caps,active_velocity_correlations,
    velocity_caps,budget_components,active_kappa_caps,primitive_component_caps,
    joint_weighted_kernel_caps,leading_component_caps,encode_leading_terms)

NAME=PREFIX+'current_native_centered_phase_conditioning.json'
RECEIPT=PREFIX+'current_native_centered_phase_conditioning_check.json'
GATE='current_original_native_centered_phase_a_nu_correlation_and_five_target_bounds_executed'


def exact_kernel_theorem():
    theorem=preceding.exact_kernel_theorem()
    y,Z,phi,psi=sy.symbols('y Z phi psi');nu=sy.Function('nu')(y,Z);K=sy.Function('K')(y,Z)
    # Fixed-free-angle residuals: substitute the phase equation only after
    # taking the nu/J derivatives. Both first cross terms remain in L_yZ.
    J=nu*K;phase=(psi+4*K)/(2*sy.pi);T2=(nu-1)*psi+4*J
    for i in (y,Z):
        L=2*sy.pi*phase*sy.diff(nu,i)-sy.diff(T2,i)
        if sy.simplify(L+4*nu*sy.diff(K,i))!=0:raise ArithmeticError('Centered first residual identity failed')
    L=2*sy.pi*phase*sy.diff(nu,y,Z)-sy.diff(T2,y,Z)
    formula=-4*(sy.diff(nu,y)*sy.diff(K,Z)+sy.diff(nu,Z)*sy.diff(K,y)+nu*sy.diff(K,y,Z))
    if sy.simplify(L-formula)!=0:raise ArithmeticError('Centered mixed residual identity failed')
    theorem['centered_original_phase_theorem']=dict(passed=True,
        exact_functions=dict(alpha='sqrt(a)*q=sigma*sqrt((2eta-Delta)/2)',v='a*nu=kappa+2alpha^2',
            Bq='b*q',M2='a*q^2',K='J/nu=(-Bq*P+M2*(H-psi/2))/v'),
        active_v_range='2<=v<=2+2eta<=3',original_smooth_branches_retained=True,
        phase='S=nu*psi+4J=2pi*nu*phi',
        L_i='-4nu*K_i',L_yZ='-4*(nu_y*K_Z+nu_Z*K_y+nu*K_yZ)',
        ordinary_sqrt_a_q_first='sqrt(a)*q_i=alpha_i-alpha*a_i/(2a)',
        ordinary_sqrt_a_q_mixed='sqrt(a)*q_yZ=alpha_yZ-(alpha_y*a_Z+alpha_Z*a_y+alpha*a_yZ)/(2a)+3alpha*a_y*a_Z/(4a^2)',
        A_first='a_i*chi_phase/2+v*K_i/(pi*(1+t^2))',
        A_mixed='a_yZ*chi_phase/2+(v*K_yZ+v_y*K_Z+v_Z*K_y)/(pi*(1+t^2))-2v*t*(t_y*K_Z+t_Z*K_y)/(pi*(1+t^2)^2)+8v*nu*t*t_psi*K_y*K_Z/(pi*(1+t^2)^3)',
        M='-a*hatT1/(2pi)-b*phi; B=E*M/2',
        centered_free_M_first='-b_i*chi_phase-(a*qP)_i/pi',
        centered_free_M_mixed='-b_yZ*chi_phase-(a*qP)_yZ/pi',
        complete_inverse_cross_terms_retained=True,phase_chi_abs_le_one=True,
        C0_bounds='|A|<=a/2<=3/2; |B|/E<=v/2<=3/2 by Cauchy of original T1/T2',
        fixed_free_angle_derivatives_not_total_inverse_derivatives=True,
        bound_products_do_not_define_original_functions=True,no_q_or_chi_or_p2_division=True)
    return theorem


def alpha_caps(c,Delta,logeta,sigma):
    C=lambda v:LogUpper.constant(c,v);L=lambda v:LogUpper(c,v);add=lambda *v:LogUpper.add(c,v)
    inveta=L(-logeta);branches={}
    for name,gamma_lower,g0 in (('body',logeta+c.ln(2),C(c.sqrt(c.mpf('1.5')))),('transition',logeta,L(logeta/2))):
        gg=L(-(c.ln(2)+gamma_lower)/2);gg2=L(-c.ln(2)/2-gamma_lower*c.mpf('1.5'))
        g=SlowJet(c,{ZERO:g0,DY:C('.5')*Delta[DY]*gg,DZ:C('.5')*Delta[DZ]*gg,
            DYZ:add(C('.5')*Delta[DYZ]*gg,C('.25')*Delta[DY]*Delta[DZ]*gg2)})
        if name=='transition':
            s=SlowJet(c,{ZERO:C(1),DY:sigma[1]*Delta[DY]*inveta,DZ:sigma[1]*Delta[DZ]*inveta,
                DYZ:add(sigma[2]*Delta[DY]*Delta[DZ]*inveta*inveta,sigma[1]*Delta[DYZ]*inveta)})
            g=g*s
        branches[name]=g
    result=SlowJet(c,{k:maximum(c,*(b[k] for b in branches.values())) for k in ORDERS})
    return result,dict(original_alpha='sigma*sqrt((2eta-Delta)/2)',branch_caps={k:v.record() for k,v in branches.items()},
        flat_branch_all_alpha_jets_exact_zero=True,no_a_denominator_in_alpha=True)


def centered_phase_caps(c,a,b,E,Delta,q,p2,logamin,logeta,logd,sigma,P,H,qP,Tslow,curvature):
    C=lambda v:LogUpper.constant(c,v);add=lambda *v:LogUpper.add(c,v)
    inva=LogUpper(c,-logamin);half=C('.5');alpha,alpha_proof=alpha_caps(c,Delta,logeta,sigma)
    chi=SlowJet(c,{k:v.divide_positive(logd) for k,v in p2.rows.items()})
    ratio={i:a[i]*inva for i in (DY,DZ,DYZ)}
    qs={ZERO:alpha[ZERO]}
    for i in (DY,DZ):qs[i]=add(alpha[i],half*alpha[ZERO]*ratio[i])
    qs[DYZ]=add(alpha[DYZ],half*alpha[DY]*ratio[DZ],half*alpha[DZ]*ratio[DY],half*alpha[ZERO]*ratio[DYZ],C('.75')*alpha[ZERO]*ratio[DY]*ratio[DZ])
    rho=previous.minimum(C(c.sqrt(3)),b[ZERO]*LogUpper(c,-logamin/2))
    bq={ZERO:previous.minimum(b[ZERO]*q[ZERO],rho*alpha[ZERO])};bqi={i:previous.minimum(b[ZERO]*q[i],rho*qs[i]) for i in (DY,DZ,DYZ)}
    for i in (DY,DZ):bq[i]=add(b[i]*q[ZERO],bqi[i])
    bq[DYZ]=add(b[DYZ]*q[ZERO],b[DY]*q[DZ],b[DZ]*q[DY],bqi[DYZ])
    # v(kappa)=kappa+sigma(1-Delta/eta)^2*(2eta-Delta).
    # Body v is constant; the original transition has gamma<=2eta.
    vk=add(C(1),C(4)*sigma[1]);vkk=C(4)*add(sigma[1],sigma[1]*sigma[1],sigma[2])*LogUpper(c,-logeta)
    v={ZERO:C(3),DY:vk*Delta[DY],DZ:vk*Delta[DZ],DYZ:add(vk*Delta[DYZ],vkk*Delta[DY]*Delta[DZ])}
    m2={ZERO:previous.minimum(alpha[ZERO]*alpha[ZERO],C('1.5'))}
    for i in (DY,DZ):m2[i]=previous.minimum(C(2)*alpha[ZERO]*alpha[i],half*add(v[i],Delta[i]))
    m2[DYZ]=previous.minimum(C(2)*add(alpha[DY]*alpha[DZ],alpha[ZERO]*alpha[DYZ]),half*add(v[DYZ],Delta[DYZ]))
    CP0,CP,CPP,CH,CHH=(C(4*c.pi),C(80*c.pi/3),C(160*c.pi),C(60*c.pi),C(800*c.pi))
    weightedP={};weightedH={}
    for i in (DY,DZ):
        weightedP[i]=CP*add(half*bqi[i],bq[ZERO]*q[ZERO]*chi[i])
        aq_qi=alpha[ZERO]*qs[i]
        weightedH[i]=CH*add(aq_qi,m2[ZERO]*q[ZERO]*chi[i])
    mixed=add(q[DY]*chi[DZ],q[DZ]*chi[DY])
    bq_yqZ=previous.minimum(b[ZERO]*q[DY]*q[DZ],bqi[DY]*q[DZ],bqi[DZ]*q[DY])
    weightedP[DYZ]=add(CP*add(half*bqi[DYZ],bq[ZERO]*mixed,bq[ZERO]*q[ZERO]*chi[DYZ]),
        CPP*add(chi[ZERO]*bq_yqZ,bq[ZERO]*mixed,bq[ZERO]*power(q[ZERO],2)*chi[DY]*chi[DZ]))
    weightedH[DYZ]=add(CH*add(alpha[ZERO]*qs[DYZ],m2[ZERO]*mixed,m2[ZERO]*q[ZERO]*chi[DYZ]),
        CHH*add(qs[DY]*qs[DZ],half*m2[ZERO]*mixed,m2[ZERO]*power(q[ZERO],2)*chi[DY]*chi[DZ]))
    D={ZERO:add(bq[ZERO]*CP0,m2[ZERO]*C(c.pi))}
    for i in (DY,DZ):D[i]=add(bq[i]*CP0,weightedP[i],m2[i]*C(c.pi),weightedH[i])
    D[DYZ]=add(bq[DYZ]*CP0,bq[DY]*P[DZ],bq[DZ]*P[DY],weightedP[DYZ],
        m2[DYZ]*C(c.pi),m2[DY]*H[DZ],m2[DZ]*H[DY],weightedH[DYZ])
    K={ZERO:previous.minimum(D[ZERO].divide_positive(c.ln(2)),C(c.pi/2))}
    for i in (DY,DZ):K[i]=add(D[i],K[ZERO]*v[i]).divide_positive(c.ln(2))
    K[DYZ]=add(D[DYZ],K[DY]*v[DZ],K[DZ]*v[DY],K[ZERO]*v[DYZ]).divide_positive(c.ln(2))
    nu=LogUpper(c,c.ln(3)-logamin);invpi=C(1/c.pi)
    inverse=dict(centered_fixed_mixed=invpi*add(v[ZERO]*K[DYZ],v[DY]*K[DZ],v[DZ]*K[DY]),
        centered_slow_cross=C(2)*invpi*v[ZERO]*add(Tslow[DY]*K[DZ],Tslow[DZ]*K[DY]),
        centered_curvature=C(8)*invpi*v[ZERO]*nu*curvature*K[DY]*K[DZ])
    common=add(*inverse.values())
    aqF={ZERO:previous.minimum(a[ZERO]*qP[ZERO],C(c.sqrt(3))*alpha[ZERO]*CP0)}
    for i in (DY,DZ):aqF[i]=add(add(CP0,half*CP)*C(c.sqrt(3))*qs[i],CP*m2[ZERO]*chi[i])
    aqF[DYZ]=add(add(CP0,half*CP)*C(c.sqrt(3))*qs[DYZ],CP*m2[ZERO]*chi[DYZ],
        add(C(2)*CP,CPP)*chi[ZERO]*qs[DY]*qs[DZ],
        add(C(2)*CP,CPP)*alpha[ZERO]*add(qs[DY]*chi[DZ],qs[DZ]*chi[DY]),CPP*m2[ZERO]*q[ZERO]*chi[DY]*chi[DZ])
    A={ZERO:C('1.5')};M={ZERO:C(3)}
    for i in (DY,DZ):
        A[i]=add(half*a[i],invpi*v[ZERO]*K[i])
        M[i]=add(b[i],invpi*add(a[i]*qP[ZERO],aqF[i],v[ZERO]*K[i]))
    A[DYZ]=add(half*a[DYZ],common)
    M[DYZ]=add(b[DYZ],invpi*add(a[DYZ]*qP[ZERO],a[DY]*qP[DZ],a[DZ]*qP[DY],aqF[DYZ]),common)
    A=SlowJet(c,A);M=SlowJet(c,M);B=(E*M)*SlowJet.constant(c,'.5')
    return dict(A=A,B=B,record=dict(original_centered_K='(-b*q*P+a*q^2*(H-psi/2))/v',
        alpha_slow_caps=alpha.record(),alpha_branch_proof=alpha_proof,
        sqrt_a_q_ordinary_jet_caps={'y%d_Z%d'%k:cap.record() for k,cap in qs.items()},
        original_bq_caps={'y%d_Z%d'%k:cap.record() for k,cap in bq.items()},
        original_aq2_caps={'y%d_Z%d'%k:cap.record() for k,cap in m2.items()},
        original_v_caps={'y%d_Z%d'%k:cap.record() for k,cap in v.items()},v_positive_lower='2',v_upper='3',
        weighted_bq_P_derivative_caps={'y%d_Z%d'%k:cap.record() for k,cap in weightedP.items()},
        weighted_aq2_H_derivative_caps={'y%d_Z%d'%k:cap.record() for k,cap in weightedH.items()},
        centered_numerator_slow_caps={'y%d_Z%d'%k:cap.record() for k,cap in D.items()},
        centered_K_fixed_free_angle_caps={'y%d_Z%d'%k:cap.record() for k,cap in K.items()},
        inverse_mixed_source_components={key:cap.record() for key,cap in inverse.items()},
        normalized_nu_cap_from_original_v_over_a=nu.record(),a_times_qP_slow_caps={'y%d_Z%d'%k:cap.record() for k,cap in aqF.items()},
        normalized_A_slow_caps=A.record(),normalized_M_slow_caps=M.record(),normalized_B_slow_caps=B.record(),
        a_nu_equals_v_used_before_bounding=True,complete_mixed_phase_and_inverse_cross_terms_retained=True,
        original_phase_endpoints_and_all_signed_chi_source_derivatives_retained=True,
        original_functions_not_defined_by_caps=True,no_q_or_chi_division=True))


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
    centered=centered_phase_caps(c,a,b,E,Delta,q,p2,logamin,logeta,logd,sigma_caps,P,H,qP,Tslow,curvature)
    for order in ORDERS:
        A.rows[order]=previous.minimum(A[order],centered['A'][order])
        B.rows[order]=previous.minimum(B[order],centered['B'][order])
    return dict(A=A,B=B,record=dict(original_whole_support_flat=False,native_source_derivative_orders=[list(k) for k in ORDERS],
        normalized_native_a_b_p2_E_Delta_caps={key:root_caps(c,roots,key).record() for key in ('a','b','p2','E','kappa_minus2')},
        active_only_a_and_b_C0_constraints=dict(a='0<a<=3',b='|b|<=3/2',t0='|t0|<=sqrt(3/a_min)'),
        native_active_kappa_correlation=active_proof,original_cutoff_slow_caps=qproof,active_t0_caps=t0.record(),u_caps=u.record(),
        joint_weighted_kernel_proof=weighted_kernel_proof,centered_original_phase_proof=centered['record'],
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


class NativeCenteredPhaseConditioning:
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
        old=next((row for row in self.old['actual_native_active_kappa_mixed_conditioning_records'].values() if packets.interval(c,row['Z_box'])._mpi_==Z._mpi_),None)
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
        owner=NativeCenteredPhaseConditioning(previous.NativeSignedAveraging(target));regions={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            regions[name]=owner.route(Z)['record'];print('Native centered-phase conditioning:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_Z_query_count=2,
        actual_native_centered_phase_conditioning_records=regions,exact_correlated_periodic_kernel_theorem=owner.theorem,
        original_signed_mixed_density_source_graph_views=previous.VIEWS,
        native_C0_y_Z_yZ_caps_and_source_integral_descendants_rebuilt=True,
        actual_controls_or_terminal_closure_installed=False,point_inverse_or_higher_jets_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Same original24-cell/17-chart route with centered free-angle residual K=J/nu, original a*nu=v correlation, full inverse mixed terms and centered B source products; all bounds stay on original functions, with native normalized velocity/input mixed rows, branchwise correlated q jets, uniform signed-r integrated kernel derivatives, normalized inverse curvature and fixed-phi A/B y/yZ caps. Original zero-mean IBP/remainder/endpoints/memory feed N^-2 Rc targets and exact conditional repair recipe. Caps are function covers; no controls, terminal closure, point inverse/higher jets, global N/cone/energy/recursion/full NS.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
