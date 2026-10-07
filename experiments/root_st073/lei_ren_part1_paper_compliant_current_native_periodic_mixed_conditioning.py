"""Native fixed-phi y/yZ caps using correlated periodic kernel identities.

Original fields, loop/inverse and signed density graphs remain unchanged.
Three smooth cutoff branches and native ordinary derivatives replace coarse
global q/Poisson estimates. Old certificates remain available independently.
"""
import ast
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_signed_averaging as previous
import lei_ren_part1_paper_compliant_current_generic_shear_source_bounds as raw_bounds

HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
packets,ep,native,current=previous.packets,previous.ep,previous.native,previous.current
serial,prior,repair=previous.serial,previous.prior,previous.repair
LogUpper,SlowJet=previous.LogUpper,previous.SlowJet
ZERO,DY,DZ,DYZ=previous.ZERO,previous.DY,previous.DZ,previous.DYZ
ORDERS,RATES,MIN_N=previous.ORDERS,previous.RATES,previous.MIN_N
NAME=PREFIX+'current_native_periodic_mixed_conditioning.json'
RECEIPT=PREFIX+'current_native_periodic_mixed_conditioning_check.json'
GATE='current_original_native_parameter_uniform_periodic_y_yZ_conditioning_and_Rc_targets_executed'


def power(cap,exponent):
    return LogUpper(cap.ctx,None if cap.log is None else cap.log*cap.ctx.mpf(str(exponent)))


def maximum(c,*caps):
    nonzero=[v for v in caps if v.log is not None]
    return LogUpper(c,None if not nonzero else c.mpf(max(ep(v.log)[1] for v in nonzero)))


def root_caps(c,roots,name):
    return SlowJet(c,{k:previous.read_cap(c,current.magnitude(roots[name][k]).record()) for k in ORDERS})


def cutoff_caps(c,a,Delta,logamin,logeta,sigma_caps,original_cutoff):
    """Original gamma/a sqrt jets, with correlated g/gamma denominators."""
    C=lambda v:LogUpper.constant(c,v);L=lambda v:LogUpper(c,v);add=lambda *v:LogUpper.add(c,v)
    inva=L(-logamin);inveta=L(-logeta);half=C('.5')
    branches={}
    for name,gamma_lower,g0 in (
        ('body_Delta_le_zero',logeta+c.ln(2),L((c.ln(3)-c.ln(2)-logamin)/2)),
        ('transition_zero_lt_Delta_lt_eta',logeta,L((logeta-logamin)/2))):
        # g=sqrt(gamma/(2a)); these are bounds of correlated functions,
        # not a product of an unrelated g0 cap and 1/gamma cap.
        gg=L(-(c.ln(2)+logamin+gamma_lower)/2)
        gg2=L(-(c.ln(2)+logamin)/2-gamma_lower*c.mpf('1.5'))
        gj={ZERO:g0}
        for i in (DY,DZ):gj[i]=half*add(Delta[i]*gg,a[i]*inva*g0)
        gj[DYZ]=add(half*Delta[DYZ]*gg,C('.25')*Delta[DY]*Delta[DZ]*gg2,
            half*a[DYZ]*inva*g0,C('.75')*a[DY]*a[DZ]*inva*inva*g0,
            C('.25')*add(Delta[DY]*a[DZ],Delta[DZ]*a[DY])*gg*inva)
        if name.startswith('body'):qj=SlowJet(c,gj)
        else:
            sigma=SlowJet(c,{ZERO:C(1),DY:sigma_caps[1]*Delta[DY]*inveta,
                DZ:sigma_caps[1]*Delta[DZ]*inveta,
                DYZ:add(sigma_caps[2]*Delta[DY]*Delta[DZ]*inveta*inveta,sigma_caps[1]*Delta[DYZ]*inveta)})
            qj=SlowJet(c,gj)*sigma
        branches[name]=qj
    q=SlowJet(c,{k:maximum(c,*(v[k] for v in branches.values())) for k in ORDERS})
    # Retain the tighter existing native three-branch C0/Z function caps.
    for order,name in ((ZERO,'q'),(DZ,'q_Z')):
        q.rows[order]=previous.minimum(q[order],previous.read_cap(c,current.magnitude(original_cutoff[name]).record()))
    return q,dict(original_gamma='2*eta-Delta',original_q='sigma(1-Delta/eta)*sqrt(gamma/(2a))',
        body_sigma_exact_one=True,flat_Delta_ge_eta_all_q_jets_exact_zero=True,
        gamma_body_positive_lower='2*eta',gamma_transition_range='[eta,2eta]',
        correlated_sqrt_ratio_bounds=dict(g_over_gamma='1/sqrt(2*a*gamma)',g_over_gamma_squared='1/(sqrt(2*a)*gamma^(3/2))'),
        exact_mixed_sqrt_rule='g_yZ=-g*Delta_yZ/(2gamma)-g*Delta_y*Delta_Z/(4gamma^2)-g*a_yZ/(2a)+3g*a_y*a_Z/(4a^2)+g*(Delta_y*a_Z+Delta_Z*a_y)/(4gamma*a)',
        branch_caps={key:value.record() for key,value in branches.items()},union_caps=q.record(),
        source_derivatives_not_clipped_or_erased=True,only_function_bounds_intersected=True)


def exact_kernel_theorem():
    r,z=sy.symbols('r cospsi',real=True);D=1-2*r*z+r*r;s=1-r*r
    cos_theta=((1+r*r)*z-2*r)/D
    if sy.cancel((cos_theta+r)/s-(z-r)/D)!=0:raise ArithmeticError('Correlated original point kernel identity failed')
    # d=beta*w=h*cos(theta)+u; theta_psi=1+2u*d and
    # h^2 sin(theta)^2=1+2u*d-d^2.
    h,u,d=sy.symbols('h u d',real=True)
    if sy.expand((h*h-(d-u)**2).subs(h*h,1+u*u)-(1+2*u*d-d*d))!=0:
        raise ArithmeticError('Original normalized kernel curvature identity failed')
    y,Z=sy.symbols('y Z');gam=sy.Function('gamma')(y,Z);a=sy.Function('a')(y,Z)
    g=sy.sqrt(gam/(2*a));gy,gZ=sy.diff(gam,y),sy.diff(gam,Z)
    formula=g*(sy.diff(gam,y,Z)/(2*gam)-gy*gZ/(4*gam**2)-sy.diff(a,y,Z)/(2*a)
        +3*sy.diff(a,y)*sy.diff(a,Z)/(4*a*a)-(gy*sy.diff(a,Z)+gZ*sy.diff(a,y))/(4*gam*a))
    if sy.simplify(sy.diff(g,y,Z)-formula)!=0:raise ArithmeticError('Original correlated sqrt mixed rule failed')
    signed=previous.functions.signed
    tree=ast.parse((HERE/Path(signed.__file__).name).read_text(encoding='utf8'))
    source=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='source_leaves')
    rows=next(n.value for n in ast.walk(source) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='rows' for t in n.targets))
    words={kw.arg:kw.value for kw in rows.keywords}
    for key,expression in (('C','tuple(E[j]-2*E[j+1] for j in range(3))'),('B','tuple(2*V[j+1] for j in range(3))')):
        if ast.dump(words[key])!=ast.dump(ast.parse(expression,mode='eval').body):raise ValueError('Original shear/velocity source identity changed')
    source=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='from_packet')
    roots=next(n.value for n in ast.walk(source) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='roots' for t in n.targets))
    words={kw.arg:kw.value for kw in roots.keywords}
    for key,expression in (('a',"dag.quotient_jet(C,E,'E')"),('b',"dag.quotient_jet(B,E,'E')")):
        if ast.dump(words[key])!=ast.dump(ast.parse(expression,mode='eval').body):raise ValueError('Original normalized shear quotient changed')
    return dict(passed=True,original_signed_r_kernel_identity=True,original_correlated_sqrt_mixed_rule=True,
        normalized_integrals=dict(P='h^-1*W1',H='h^-2*W2',d='h^-1*w'),
        exact_kernel_relations=dict(d='h*cos(theta)+u',theta_psi='1+2u*d',d_psi_squared='(1+2u*d)^2*(1+2u*d-d^2)'),
        signed_parameter_u_and_both_r_sectors=True,
        whole_parameter_caps=dict(W1= '4pi',W1_u='16pi',W1_uu='40pi',P='4pi',P_u='20pi',P_uu='100pi',
            H='pi',H_u='50pi',H_uu='600pi',d='2h',d_u=6,d_uu=50,theta_u=2,theta_uu=14),
        inside_r_abs_le_half='Fourier sums bound |w|<=2, |w_r|<=4, |w_rr|<=16 and |W1_r|<=4, |W1_rr|<=16. H_r<=40pi, H_rr<=272pi, hence H_uu<=392pi.',
        outside_r_abs_ge_half='theta_u<=2,theta_uu<=14,r_u<=1,|r_uu|<=3,s_u<=2,|s_uu|<=10. W1=(theta-psi)/(2r) gives W1_uu<=22+28pi<40pi. H=K/(4r^2), |K|<=6pi+2,|K_u|<=16pi+12,|K_uu|<=80pi+104; reciprocal derivatives<=1,4,36 give H_u<=40pi+20<50pi and H_uu<=424pi+272<600pi.',
        H_C0='H is nonnegative cumulative integral and H(2pi)=pi by original Fourier orthogonality',
        point_d_derivatives='|w|<=2h^2,|w_r|<=4h^4,|w_rr|<=16h^6; beta=h^-1, r_u=h^-3 imply |d_u|<=6,|d_uu|<=44<50',
        curvature_conditioning='u=(p2/dstar)*q and t=t0+2q*d imply |t_psi|<=2q*(1+abs(p2/dstar)*abs(t-t0))^(3/2). Dividing by (1+t^2)^2 or multiplying by 2|t|/(1+t^2)^3 gives <=4q*(1+abs(p2/dstar)*(1+|t0|))^(3/2). No division by q is used.',
        slow_derivatives_hold_phi_fixed=True,Poisson_denominator_inverse_powers_not_used=True,
        native_point_inverse_or_higher_mixed_jets_installed=False)


def active_velocity_correlations(c,velocity,a,b):
    """Bounds restricted to modulation support, never global field caps."""
    E=velocity['E'];V=velocity.get('V');C=lambda v:LogUpper.constant(c,v);add=lambda *v:LogUpper.add(c,v)
    before={key:value.record() for key,value in velocity.items()}
    E.rows[DY]=previous.minimum(E[DY],E[ZERO])
    E.rows[DYZ]=previous.minimum(E[DYZ],add(E[DZ],C('.5')*a[DZ]*E[ZERO]))
    b0=previous.minimum(b[ZERO],C('1.5'))
    if V is not None:
        V.rows[DY]=previous.minimum(V[DY],C('.5')*b0*E[ZERO])
        V.rows[DYZ]=previous.minimum(V[DYZ],C('.5')*add(b[DZ]*E[ZERO],b0*E[DZ]))
    return dict(unrestricted_native_velocity_caps=before,active_support_only_velocity_caps={key:value.record() for key,value in velocity.items()},
        exact_same_original_source_identities=dict(E_y='(1-a)*E/2',E_yZ='(1-a)*E_Z/2-a_Z*E/2',
            V_y='b*E/2',V_yZ='(b_Z*E+b*E_Z)/2'),
        active_only_a_le3_and_b_abs_le1point5=True,
        restrictions_used_only_in_modulation_densities_and_primitive_jets=True,
        original_velocity_outside_q_support_not_clipped_or_zeroed=True)


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
    T1=t0*SlowJet.constant(c,2*c.pi)+two*q*P
    T2=t0*t0*SlowJet.constant(c,2*c.pi)+four*t0*q*P+four*q*q*H
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
        original_cutoff_slow_caps=qproof,active_t0_caps=t0.record(),u_caps=u.record(),
        fixed_free_angle_P_caps=P.record(),fixed_free_angle_H_caps=H.record(),
        original_fixed_angle_T1_caps=T1.record(),original_fixed_angle_T2_caps=T2.record(),
        inverse_first_and_mixed_numerator_caps={str(k):v.record() for k,v in numerator.items()},
        normalized_point_curvature_cap=curvature.record(),inverse_yZ_cap=psiYZ.record(),
        correlated_composed_T1_caps={str(k):v.record() for k,v in hat.items()},
        normalized_A_slow_caps=A.record(),normalized_B_over_Pstar_slow_caps=B.record(),active_support_E_correlations=velocity_proof,
        all_source_derivatives_retained=True,active_C0_constraints_do_not_change_source_functions=True,
        no_positive_source_scale_or_Poisson_denominator_materialized=True,higher_mixed_jets_not_admitted=True))


def velocity_caps(c,packet):
    return {name:SlowJet(c,{k:previous.read_cap(c,raw_bounds.modal_partial_bound(packet,packet.velocity[component][k[0]],k[1]).record())
        for k in ORDERS}) for name,component in (('E','theta'),('V','axial'))}


def budget_components(leading,remainder,factors,key):
    c=leading.ctx;C=lambda v:LogUpper.constant(c,v)
    decay=current.magnitude(factors['decay']);mass=current.magnitude(factors['mass'])
    rate=RATES[key];rate=C(c.mpf(rate.numerator)/rate.denominator);components={}
    for derivative,order,slow in (('C0',ZERO,DY),('Z',DZ,DYZ)):
        G=leading[order]*C('.5');Gy=leading[slow]*C('.5')
        for label,cap in (('right_endpoint',G),('decayed_left_endpoint',decay*G),
            ('slow_variable_integral',mass*Gy),('kernel_derivative_integral',mass*rate*G),
            ('quadratic_remainder_integral',mass*remainder['Z_derivatives' if derivative=='Z' else 'values'][key])):
            components.setdefault(label,{})[derivative]=cap
    return components


class NativePeriodicMixedConditioning:
    def __init__(self,owner):
        if type(owner) is not previous.NativeSignedAveraging:raise ValueError('Same original signed averaging owner required')
        checked=json.loads((HERE/previous.RECEIPT).read_bytes());self.old=json.loads((HERE/previous.NAME).read_bytes())
        if not checked['all_passed'] or not checked[previous.GATE] or checked['source_family']!=owner.family:
            raise ValueError('Checked same-original signed averaging route required')
        self.owner=owner;self.target=owner.owner;self.ctx=owner.ctx;self.family=owner.family;self.coordinates=owner.coordinates;self.service=owner.service
        self.service.bind_hashes(checked['input_hashes']);self.service.bind_hashes({previous.NAME:sha(previous.NAME),previous.RECEIPT:sha(previous.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)})
        manifest=json.loads((HERE/(PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes())
        self.sigma={int(k):previous.read_cap(self.ctx,v) for k,v in manifest['sigma_global_derivative_log_caps'].items()}
        self.theorem=exact_kernel_theorem()

    @native.inlet.source_precision
    def route(self,Z=(-1,1)):
        c=self.ctx;coords=self.coordinates;Z=c.mpf(Z)
        old=next((row for row in self.old['actual_original_native_signed_averaging_records'].values() if packets.interval(c,row['Z_box'])._mpi_==Z._mpi_),None)
        if old is None:raise ValueError('This replay adapter binds either original accepted Z domain')
        root_owner=self.target.q_owner.owner.owner
        eta=packets.interval(c,root_owner.scales['selected_positive_eta_log']);dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        zero=lambda:{key:coords.scalar(0) for key in RATES}
        running,runningZ=zero(),zero();cells=[]
        for (label,chart,left,right),oldcell in zip(current.ROUTE,old['cells']):
            if (oldcell['label'],oldcell['chart'])!=(label,chart):raise ValueError('Same original ordered cell/function ledger required')
            geom=self.target.geometry(label,chart,left,right);factors={key:current.transfer.true_width_kernel(coords,geom,rate) for key,rate in RATES.items()}
            details={};components={}
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
                remainder=previous.remainder_caps(bounded);values={};jets={}
                from lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets_check import record_value
                for key in RATES:
                    components[key]=budget_components(bounded['leading'][key],remainder,factors[key],key)
                    cap,proof=previous.ibp_contribution(bounded['leading'][key],remainder['values'][key],factors[key],key)
                    capZ,proofZ=previous.ibp_contribution(bounded['leading'][key],remainder['Z_derivatives'][key],factors[key],key,True)
                    oldcap=current.magnitude(record_value(oldcell['true_width_source_Nminus2_C0_covers'][key],coords))
                    oldcapZ=current.magnitude(record_value(oldcell['true_width_source_Nminus2_Z_covers'][key],coords))
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
            cells.append(dict(record=record,values=values,Z_derivatives=jets,factors=factors,incoming=incoming,incoming_Z=incomingZ,cumulative=running,cumulative_Z=runningZ,budget_components=components))
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
        suffix={key:coords.scalar(1) for key in RATES};weighted=[];dominant={};first_budget={};first_dominant={}
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
        owner=NativePeriodicMixedConditioning(previous.NativeSignedAveraging(target));regions={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            regions[name]=owner.route(Z)['record'];print('Native correlated periodic mixed conditioning:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_Z_query_count=2,
        actual_native_periodic_mixed_conditioning_records=regions,exact_correlated_periodic_kernel_theorem=owner.theorem,
        original_signed_mixed_density_source_graph_views=previous.VIEWS,
        native_C0_y_Z_yZ_caps_and_source_integral_descendants_rebuilt=True,
        actual_controls_or_terminal_closure_installed=False,point_inverse_or_higher_jets_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Same original full24-cell/17-chart route on two Z domains, with native normalized velocity/input mixed rows, branchwise correlated q jets, uniform signed-r integrated kernel derivatives, normalized inverse curvature and fixed-phi A/B y/yZ caps. Original zero-mean IBP/remainder/endpoints/memory feed N^-2 Rc targets and exact conditional repair recipe. Caps are function covers; no controls, terminal closure, point inverse/higher jets, global N/cone/energy/recursion/full NS.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
