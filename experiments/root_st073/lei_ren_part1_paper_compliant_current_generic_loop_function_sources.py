"""Section 11 analytic loop, changed velocities and five signed source rates.

Functional graph nodes bind the original signed source functions, not
selected values of saved covers. eta,d_star,R,N and the inverse remain
formal. Integral nodes have explicit bound variables and monotonicity.
This stage defines candidates and rates; whole changed transport, repair,
N, high physical derivative admission and a point field remain separate.
"""
import ast
import gzip
import json
import math
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_signed_jets as signed

packets,current=signed.packets,signed.current
HERE,PREFIX,sha=signed.HERE,signed.PREFIX,signed.sha
NAME=PREFIX+'current_generic_loop_function_sources.json'
RECEIPT=PREFIX+'current_generic_loop_function_sources_check.json'
VIEWS=PREFIX+'current_generic_loop_function_sources_views.json.gz'
GATE='current_original_generic_loop_velocity_and_five_increment_function_graphs_defined'
OPEN=signed.OPEN
FIRST=((1,0),(0,1));ZERO=(0,0)


class FrequencyLogBound:
    """Positive majorant sum_p exp(coefficient_log[p])*N^p, N unchosen."""
    def __init__(self,c,terms):
        self.ctx=c;self.terms={p:v for p,v in terms.items() if v.log is not None}
        if any(type(p) is not int or v.ctx is not c for p,v in self.terms.items()):raise ValueError('Integer N powers in the same source context required')
    def __add__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Same frequency source context required')
        terms={}
        for p in set(self.terms)|set(other.terms):
            terms[p]=current.bounds.LogUpper.add(self.ctx,[q for q in (self.terms.get(p),other.terms.get(p)) if q is not None])
        return FrequencyLogBound(self.ctx,terms)
    def __mul__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Same frequency source context required')
        terms={}
        for p,left in self.terms.items():
            for q,right in other.terms.items():
                key=p+q;value=left*right
                terms[key]=current.bounds.LogUpper.add(self.ctx,[terms[key],value]) if key in terms else value
        return FrequencyLogBound(self.ctx,terms)
    def scale(self,cap):return FrequencyLogBound(self.ctx,{p:v*cap for p,v in self.terms.items()})
    def evaluate(self,logN):
        return current.bounds.LogUpper.add(self.ctx,[current.bounds.LogUpper(self.ctx,cap.log+p*logN) for p,cap in self.terms.items()])
    def record(self):return dict(exact_zero=not self.terms,terms=[dict(N_power=p,coefficient=v.record()) for p,v in sorted(self.terms.items())],
        actual_N_or_coefficient_exponential_not_materialized=True)


def frequency_source_bounds(c,raw,loop,*,whole_q_flat=False):
    """For logN>=max(0,log|A|), bound expm1 and retain O(1) fast y terms."""
    read=lambda row:current.bounds.LogUpper(c,None if row['exact_zero'] else packets.interval(c,row['log_absolute_upper']))
    C=lambda value:current.bounds.LogUpper.constant(c,value)
    poly=lambda value,p=0:FrequencyLogBound(c,{p:value})
    W=lambda name,j=0,k=0:read(raw[name]['y%d_Z%d'%(j,k)])
    A=lambda j=0,k=0:read(loop['slow_phase_held_log_bounds']['A']['y%d_Z%d'%(j,k)])
    B=lambda j=0,k=0:read(loop['slow_phase_held_log_bounds']['B']['y%d_Z%d'%(j,k)])
    expone=current.bounds.LogUpper(c,c.mpf(1));two=C(2);half=C('.5')
    ep=packets.recovery.endpoints
    minimum=c.mpf(max(ep(c.mpf(0))[1],ep(A().log)[1] if A().log is not None else 0))
    E,V=poly(W('E')),poly(W('V'))
    if whole_q_flat:
        z=FrequencyLogBound(c,{})
        return dict(required_positive_log_N_lower=minimum,increment_velocity_log_polynomials=dict(delta_E=z.record(),delta_V=z.record()),
            five_increment_rate_log_polynomials={k:z.record() for k in ('m','h','k','e','p')},
            five_increment_rate_first_log_polynomials={i:{k:z.record() for k in ('m','h','k','e','p')} for i in ('y','Z')},
            whole_source_q_flat_by_original_power_and_eta_boundary_proof=True,
            incoming_five_history_defects_not_zeroed=True)
    rr=poly(W('E')*A()*expone,-1);ss=poly(B(),-1)
    rr_i={};ss_i={}
    for key,j,k in (('y',1,0),('Z',0,1)):
        rr_i[key]=poly(W('E',j,k)*A()*expone,-1)+poly(W('E')*A(j,k)*expone,-1)
        ss_i[key]=poly(B(j,k),-1)
        if key=='y':
            rr_i[key]=rr_i[key]+poly(W('E')*read(loop['phase_and_first_slow_phase_log_bounds']['A_phi'])*expone)
            ss_i[key]=ss_i[key]+poly(read(loop['phase_and_first_slow_phase_log_bounds']['B_phi']))
    inc=dict(m=ss,h=rr,k=V*rr+E*ss+rr*ss,
        e=(V*ss).scale(two)+ss*ss+E*rr+(rr*rr).scale(half),p=E*rr+(rr*rr).scale(half))
    derivatives={}
    for key,j,k in (('y',1,0),('Z',0,1)):
        Ei,Vi=poly(W('E',j,k)),poly(W('V',j,k));ri,si=rr_i[key],ss_i[key]
        derivatives[key]=dict(m=si,h=ri,k=Vi*rr+V*ri+Ei*ss+E*si+ri*ss+rr*si,
            e=(Vi*ss+V*si+ss*si).scale(two)+Ei*rr+E*ri+rr*ri,
            p=Ei*rr+E*ri+rr*ri)
    return dict(required_positive_log_N_lower=minimum,
        increment_velocity_log_polynomials=dict(delta_E=rr.record(),delta_V=ss.record()),
        five_increment_rate_log_polynomials={k:v.record() for k,v in inc.items()},
        five_increment_rate_first_log_polynomials={key:{k:v.record() for k,v in row.items()} for key,row in derivatives.items()},
        expm1_bound='|A/N|<=1 implies |expm1(A/N)|<=e*|A|/N and exp(A/N)<=e',
        fast_y_derivatives_have_N_power0_terms=True,Z_and_value_sources_have_only_negative_N_powers=True,
        nonzero_original_V_and_all_signed_density_cross_terms_bounded=True,
        whole_source_q_flat_by_original_power_and_eta_boundary_proof=False,
        actual_N_not_selected=True)


class FunctionGraph:
    """Acyclic functional instructions, including lazy branch and binders."""
    def __init__(self,inputs):
        self.nodes=list(inputs['jet_expression_dag']['nodes'])
        self.keys={json.dumps(value,sort_keys=True):i for i,value in enumerate(self.nodes)}
        self.source=inputs;self.zero=self.integer(0);self.one=self.integer(1)
    def node(self,operation,**attributes):
        value=dict(operation=operation,**attributes);key=json.dumps(packets.encode(value),sort_keys=True)
        if key not in self.keys:self.keys[key]=len(self.nodes);self.nodes.append(value)
        return self.keys[key]
    def integer(self,n):
        if type(n) is not int:raise ValueError('Exact integer graph weight required')
        return self.node('constant',value=n)
    def add(self,*args):
        args=[q for q in args if q!=self.zero]
        return self.zero if not args else args[0] if len(args)==1 else self.node('sum',arguments=args)
    def neg(self,arg):return self.zero if arg==self.zero else self.node('negative',argument=arg)
    def sub(self,left,right):return self.add(left,self.neg(right))
    def mul(self,*args):
        if self.zero in args:return self.zero
        args=[q for q in args if q!=self.one]
        return self.one if not args else args[0] if len(args)==1 else self.node('product',arguments=args)
    def div(self,left,right,condition):
        return self.zero if left==self.zero else self.node('positive_function_quotient',
            numerator=left,denominator=right,positive_certificate=condition)
    def function(self,name,arg):return self.node('analytic_unary',name=name,argument=arg)
    def root(self,name,k=ZERO):return self.source['jet_expression_dag']['roots'][name]['y%d_Z%d'%k]
    def integral(self,body,angle):
        return self.node('angle_integral',integrand=body,angle_variable='psi',upper_angle=angle,
            lower_angle_exact_zero=True,ordinary_slow_parameters_fixed=True)
    def at_inverse(self,body,inverse):
        return self.node('substitute_inverse_angle',body=body,inverse_angle=inverse,angle_variable='psi')
    def flat(self,body,Delta,eta):
        return self.node('flat_zero_branch',active_body=body,Delta=Delta,eta=eta,
            flat_when_Delta_ge_eta=True,flat_value=self.zero,
            active_body_not_evaluated_on_flat_branch=True)


def build(inputs,V_rows,log_scales,loop_bounds):
    """First slow derivatives plus actual fast phase for E_N,V_N and rates."""
    graph=FunctionGraph(inputs);g=graph
    one,two,four=g.one,g.integer(2),g.integer(4)
    pi=g.node('mathematical_pi');twopi=g.mul(two,pi)
    inv=lambda node,certificate:g.div(one,node,certificate)
    half=inv(two,'exact_positive_integer2');inv2pi=inv(twopi,'exact_positive_2pi')
    psi=g.node('function_variable',name='psi',domain='closed[0,2pi]')
    phi=g.node('function_variable',name='phi',domain='periodic_modulo1')
    N=g.node('shared_positive_integer_parameter',name='common_N',domain='N>=1',same_for_all_charts=True)
    invN=inv(N,'common_N_positive_integer')
    scale=lambda name:g.node('original_positive_log_scale',name=name,
        selected_source_log=log_scales[name],source_scale_never_materialized=True)
    eta=scale('eta');dstar=scale('d_star')
    a,b,E,p1,p2,t0=(g.root(key) for key in ('a','b','E','p1','p2','t0'))
    # Delta is a signed expression, with correlated original mu cutoffs on O3.
    Delta=g.root('kappa_minus2')
    if inputs['original_O3_small_excess_definition']:
        Delta=g.node('original_correlated_O3_excess',definition=inputs['original_O3_small_excess_definition'],
            original_positive_mu=inputs['original_O3_positive_mu'],
            coordinate_source=inputs['source_provenance']['coordinate_contract'],
            original_source_provenance=inputs['source_provenance'],generic_equal_expression=Delta,
            derivative_y_order=0,derivative_Z_order=0,rounded_kappa_minus2_forbidden=True)
    delta_i={}
    for i in FIRST:
        value=g.root('kappa_minus2',i)
        if inputs['original_O3_small_excess_definition']:
            value=g.node('original_correlated_O3_excess',definition=inputs['original_O3_small_excess_definition'],
                original_positive_mu=inputs['original_O3_positive_mu'],
                coordinate_source=inputs['source_provenance']['coordinate_contract'],
                original_source_provenance=inputs['source_provenance'],generic_equal_expression=value,
                derivative_y_order=i[0],derivative_Z_order=i[1],rounded_kappa_minus2_forbidden=True)
        delta_i[i]=value
    x=g.sub(one,g.div(Delta,eta,'actual_eta_positive'))
    gamma=g.sub(g.mul(two,eta),Delta)
    sigma=g.function('original_flat_sigma',x)
    root=g.function('positive_sqrt',g.div(gamma,g.mul(two,a),'active_gamma_and_a_positive'))
    q=g.flat(g.mul(sigma,root),Delta,eta)
    q_i={}
    for i in FIRST:
        sig_i=g.neg(g.mul(g.function('original_flat_sigma_prime',x),g.div(delta_i[i],eta,'actual_eta_positive')))
        Li=g.neg(g.mul(half,g.add(g.div(delta_i[i],gamma,'active_gamma_ge_eta_positive'),
            g.div(g.root('a',i),a,'a'))))
        q_i[i]=g.flat(g.mul(root,g.add(sig_i,g.mul(sigma,Li))),Delta,eta)
    u=g.div(g.mul(p2,q),dstar,'actual_d_star_positive')
    h=g.function('positive_sqrt',g.add(one,g.mul(u,u)))
    ih=inv(h,'h_ge1');r=g.mul(u,ih);alpha=g.mul(two,q,ih)
    ri={};alphai={}
    for i in FIRST:
        ui=g.div(g.add(g.mul(g.root('p2',i),q),g.mul(p2,q_i[i])),dstar,'actual_d_star_positive')
        ri[i]=g.mul(ui,ih,ih,ih)
        hi_inv_i=g.neg(g.mul(u,ui,ih,ih,ih))
        alphai[i]=g.mul(two,g.add(g.mul(q_i[i],ih),g.mul(q,hi_inv_i)))
    # Correlated positive denominator: 1-|r| never subtracts rounded units.
    rho=g.div(one,g.mul(h,g.add(h,g.function('absolute',u))),'h_h_plus_abs_u_positive')
    trig=g.node('Poisson_half_angle_squared',angle=psi,signed_r=r,
        expression='sin(psi/2)^2 if r>=0 else cos(psi/2)^2')
    D=g.add(g.mul(rho,rho),g.mul(four,g.function('absolute',r),trig))
    nn=g.node('correlated_cos_minus_r',angle=psi,signed_r=r,one_minus_abs_r=rho)
    w=g.div(nn,D,'correlated_Poisson_D_positive')
    wr=g.div(g.add(g.neg(D),g.mul(two,nn,nn)),g.mul(D,D),'correlated_Poisson_D_squared_positive')
    wpsi=g.neg(g.div(g.mul(ih,ih,g.function('sin',psi)),g.mul(D,D),'correlated_Poisson_D_squared_positive'))
    t=g.add(t0,g.mul(alpha,w))
    ti={i:g.add(g.root('t0',i),g.mul(alphai[i],w),g.mul(alpha,wr,ri[i])) for i in FIRST}
    tpsi=g.mul(alpha,wpsi)
    ratio=g.add(one,g.mul(t0,t0),g.mul(two,q,q))
    K=g.div(one,g.mul(twopi,ratio),'correlated_ratio_ge1')
    Ki={i:g.neg(g.div(g.mul(K,g.add(g.mul(two,t0,g.root('t0',i)),g.mul(four,q,q_i[i]))),ratio,'correlated_ratio_ge1')) for i in FIRST}
    T1=g.integral(t,psi);T2=g.integral(g.mul(t,t),psi)
    P=g.add(psi,T2);Phi=g.mul(K,P)
    lambda0=g.mul(K,g.add(one,g.mul(t,t)))
    inverse=g.node('monotone_phase_inverse',phase_function=Phi,angle_variable='psi',target_phase=phi,
        angle_endpoints_exact=['0','2pi'],phase_endpoints_exact=['0','1'],
        positive_derivative_log_lower=loop_bounds['positive_conditioning']['log_lambda_positive_lower'],
        flat_Delta=Delta,flat_eta=eta,
        exact_flat_inverse='psi=2pi*phi')
    at=lambda body:g.at_inverse(body,inverse)
    ll=at(lambda0);ilambda=inv(ll,'actual_monotone_phase_lambda_positive')
    psi_i={}
    for i in FIRST:
        T2i=g.integral(g.mul(two,t,ti[i]),psi)
        Phii=g.add(g.mul(Ki[i],P),g.mul(K,T2i))
        psi_i[i]=g.neg(g.mul(at(Phii),ilambda))
    chi=g.sub(phi,g.mul(inverse,inv2pi))
    T1hat=at(T1);tt=at(t)
    M=g.sub(g.neg(g.mul(a,T1hat,inv2pi)),g.mul(b,phi))
    A=g.flat(g.mul(half,a,chi),Delta,eta)
    BB=g.flat(g.mul(half,E,M),Delta,eta)
    Ai={};Bi={}
    for i in FIRST:
        T1hati=g.add(at(g.integral(ti[i],psi)),g.mul(tt,psi_i[i]))
        Mi=g.sub(g.neg(g.mul(g.add(g.mul(g.root('a',i),T1hat),g.mul(a,T1hati)),inv2pi)),g.mul(g.root('b',i),phi))
        Ai[i]=g.flat(g.mul(half,g.sub(g.mul(g.root('a',i),chi),g.mul(a,psi_i[i],inv2pi))),Delta,eta)
        Bi[i]=g.flat(g.mul(half,g.add(g.mul(g.root('E',i),M),g.mul(E,Mi))),Delta,eta)
    Aphi=g.flat(g.mul(half,g.sub(a,g.mul(a,ilambda,inv2pi))),Delta,eta)
    Bphi=g.flat(g.mul(half,E,g.sub(g.neg(g.mul(a,tt,ilambda,inv2pi)),b)),Delta,eta)
    V={}
    for i in (ZERO,*FIRST):
        label='V_y%d_Z%d'%i
        leaf=packets.encode(signed.expressions.RadiusPolynomial(V_rows[i].algebra,{0:V_rows[i]}).record())
        g.source['jet_expression_dag']['source_derivative_leaves'][label]=leaf
        V[i]=g.node('source_derivative',name=label)
    exparg=g.mul(A,invN);expfactor=g.function('exp',exparg)
    deltaE=g.mul(E,g.function('expm1',exparg));deltaV=g.mul(BB,invN)
    EN=g.mul(E,expfactor);VN=g.add(V[ZERO],deltaV)
    xpartials={(1,0):g.add(g.mul(Ai[(1,0)],invN),Aphi),(0,1):g.mul(Ai[(0,1)],invN)}
    Bpartials={(1,0):g.add(g.mul(Bi[(1,0)],invN),Bphi),(0,1):g.mul(Bi[(0,1)],invN)}
    ENi={i:g.mul(expfactor,g.add(g.root('E',i),g.mul(E,xpartials[i]))) for i in FIRST}
    VNi={i:g.add(V[i],Bpartials[i]) for i in FIRST}
    riE={i:g.add(g.mul(g.root('E',i),g.function('expm1',exparg)),g.mul(E,expfactor,xpartials[i])) for i in FIRST}
    riV=Bpartials
    def rates(E0,V0,r0,s0):
        return dict(m=s0,h=r0,k=g.add(g.mul(V0,r0),g.mul(E0,s0),g.mul(r0,s0)),
            e=g.sub(g.add(g.mul(two,V0,s0),g.mul(s0,s0)),g.add(g.mul(E0,r0),g.mul(half,r0,r0))),
            p=g.add(g.mul(E0,r0),g.mul(half,r0,r0)))
    flux=rates(E,V[ZERO],deltaE,deltaV);flux_i={}
    for i in FIRST:
        Ei,Vi,rr,ss=g.root('E',i),V[i],riE[i],riV[i]
        flux_i[i]=dict(m=ss,h=rr,
            k=g.add(g.mul(Vi,deltaE),g.mul(V[ZERO],rr),g.mul(Ei,deltaV),g.mul(E,ss),g.mul(rr,deltaV),g.mul(deltaE,ss)),
            e=g.sub(g.add(g.mul(two,Vi,deltaV),g.mul(two,V[ZERO],ss),g.mul(two,deltaV,ss)),
                g.add(g.mul(Ei,deltaE),g.mul(E,rr),g.mul(deltaE,rr))),
            p=g.add(g.mul(Ei,deltaE),g.mul(E,rr),g.mul(deltaE,rr)))
    v=g.mul(a,ratio);aL=g.div(v,g.add(one,g.mul(tt,tt)),'one_plus_t_squared_positive')
    bL=g.neg(g.mul(aL,tt));vminus2=g.add(Delta,g.mul(two,a,q,q))
    H=g.add(p1,g.mul(p2,tt));J=g.sub(p2,g.mul(p1,tt));frozenD=g.sub(H,v)
    frozenQ=g.sub(g.mul(two,frozenD,frozenD),g.mul(vminus2,J,J))
    aN=g.sub(aL,g.mul(two,Ai[(1,0)],invN))
    bN=g.mul(g.function('exp',g.neg(exparg)),g.add(bL,g.div(g.mul(two,Bi[(1,0)],invN),E,'E')))
    return dict(chart=inputs['chart'],source_family=inputs['source_family'],
        function_graph_nodes=g.nodes,original_signed_input_graph=inputs,
        roots=dict(q=q,Delta=Delta,eta=eta,d_star=dstar,N=N,Phi=Phi,inverse_angle=inverse,
            A=A,B_over_Pstar=BB,A_phi=Aphi,B_over_Pstar_phi=Bphi,
            A_y_slow=Ai[(1,0)],A_Z_slow=Ai[(0,1)],B_y_slow=Bi[(1,0)],B_Z_slow=Bi[(0,1)],
            E_N=EN,V_N=VN,delta_E=deltaE,delta_V=deltaV,a_N=aN,b_N=bN,
            direction_at_free_angle=t,T1_at_free_angle=T1,T2_at_free_angle=T2,lambda_at_free_angle=lambda0,
            a_L=aL,b_L=bL,v_minus2=vminus2,frozen_D=frozenD,frozen_Q=frozenQ),
        changed_velocity_first_ordinary_derivatives={
            'theta':{'y':ENi[(1,0)],'Z':ENi[(0,1)]},'axial':{'y':VNi[(1,0)],'Z':VNi[(0,1)]}},
        five_signed_increment_rate_roots=flux,
        five_signed_increment_rate_first_derivatives={('y' if i==(1,0) else 'Z'):row for i,row in flux_i.items()},
        common_phase_binding=dict(phi='fractional_part(N*log(R/r_minus))',
            N_parameter='common_N, same positive integer on all charts',
            radius_source=inputs['original_radius_source'],
            microscopic_left_offset_and_actual_log_r_minus_source=PREFIX+'current_generic_shear_loop_domain.json/current_original_generic_loop_domain',
            slow_derivatives_hold_phi_fixed=True,total_y_adds_N_times_phase_derivative=True,total_Z_phase_derivative_exact_zero=True,
            local_chart_selector_not_used_as_fast_phase=True),
        own_defect_transport_contract=dict(rates=dict(m=1,h='1.5',k='1.5',e=1,p=0),
            formula='D_j(y,Z)=exp(-rate_j*(y-y_in))*D_j(y_in,Z)+integral_y_in^y exp(-rate_j*(y-s))*delta_density_j(s,Z)ds',
            unchanged_original_P0=True,actual_original_inlet_histories_required=True,
            quiet_gap_rule='delta_density=0 but D_out=exp(-rate*width)*D_in; pressure rate0 preserves memory',
            same_global_source_family_and_N_across_interfaces=True),
        scope='Actual analytic source-bound Section11 function graph and changed five signed source rates/first y,Z derivatives, before own transport integration/new repair/N and high-order physical admission.',
        actual_scale_radius_or_saved_coefficient_value_materialized=False,
        signed_function_point_values_not_taken_from_saved_covers=True)


def exact_theorem():
    E,V,r,u,a,b,Ay,By,Aphi,Bphi,N=s.symbols('E V r u a b Ay By Aphi Bphi N',real=True)
    checks={}
    densities=lambda ee,vv:dict(m=vv,h=ee,k=ee*vv,e=vv*vv-ee*ee/2,p=ee*ee/2)
    old=densities(E,V);new=densities(E+r,V+u)
    inc=dict(m=u,h=r,k=V*r+E*u+r*u,e=2*V*u+u*u-E*r-r*r/2,p=E*r+r*r/2)
    for key in inc:
        if s.expand(new[key]-old[key]-inc[key])!=0:raise ArithmeticError('Changed signed density identity failed')
        checks['full_signed_increment_'+key]=True
    y,z,phi=s.symbols('y Z phi');AA=s.Function('A')(y,z,phi);BB=s.Function('B')(y,z,phi)
    ef=s.Function('E')(y,z);vf=s.Function('V')(y,z)
    total=lambda f:s.diff(f,y)+N*s.diff(f,phi)
    EN=ef*s.exp(AA/N);VN=vf+BB/N
    expectedE=s.exp(AA/N)*(s.diff(ef,y)+ef*(s.diff(AA,y)/N+s.diff(AA,phi)))
    expectedV=s.diff(vf,y)+s.diff(BB,y)/N+s.diff(BB,phi)
    for label,left,right in (('changed_E_total_y',total(EN),expectedE),('changed_V_total_y',total(VN),expectedV),
        ('changed_E_Z',s.diff(EN,z),s.exp(AA/N)*(s.diff(ef,z)+ef*s.diff(AA,z)/N)),
        ('changed_V_Z',s.diff(VN,z),s.diff(vf,z)+s.diff(BB,z)/N)):
        if s.simplify(left-right)!=0:raise ArithmeticError('Common fast phase derivative failed')
        checks[label]=True
    sourcefile=PREFIX+'current_generic_shear_loop.py';tree=ast.parse((HERE/sourcefile).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='GenericShearLoop')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='modulate')
    wanted={'x':"point['A']/N",'delta_theta':'self.Utheta*c.expm1(x)','delta_z':"point['B']/N"}
    for target,expression in wanted.items():
        found=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        if len(found)!=1 or ast.dump(found[0])!=ast.dump(ast.parse(expression,mode='eval').body):raise ValueError('Original velocity modulation changed')
        checks['original_modulate_'+target]=True
    result=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Return) and isinstance(n.value,ast.Call))
    keywords={kw.arg:kw.value for kw in result.keywords}
    for target,expression in {'Utheta':'self.Utheta*c.exp(x)','Uz':'Uz+delta_z',
        'a_N':"point['aL']-2*Ay/N",'b_N':"c.exp(-x)*(point['bL']+2*By/(N*self.Utheta))"}.items():
        if ast.dump(keywords[target])!=ast.dump(ast.parse(expression,mode='eval').body):raise ValueError('Original changed velocity/shear keyword changed')
        checks['original_modulate_'+target]=True
    densityfile=PREFIX+'current_generic_shear_moment_recovery.py'
    densitytree=ast.parse((HERE/densityfile).read_text(encoding='utf8'))
    for name,expected in (
        ('history_densities',{'m':'V','h':'E','k':'E*V','e':'V*V-E*E/2','p':'E*E/2'}),
        ('increment_densities',{'m':'u','h':'e','k':'V*e+E*u+e*u','e':'2*V*u+u*u-E*e-e*e/2','p':'E*e+e*e/2'})):
        fn=next(n for n in densitytree.body if isinstance(n,ast.FunctionDef) and n.name==name)
        result=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Return));keywords={kw.arg:kw.value for kw in result.keywords}
        for key,expression in expected.items():
            if ast.dump(keywords[key])!=ast.dump(ast.parse(expression,mode='eval').body):raise ValueError('Original five density source keyword changed')
            checks['original_'+name+'_'+key]=True
    rho,angle=s.symbols('abs_r angle',real=True)
    for sign in (1,-1):
        rr=sign*rho;gap=1-rho
        trig=s.sin(angle/2)**2 if sign==1 else s.cos(angle/2)**2
        den=gap**2+4*rho*trig
        numerator=gap-2*s.sin(angle/2)**2 if sign==1 else 2*s.cos(angle/2)**2-gap
        if s.trigsimp(s.expand(den-(1-2*rr*s.cos(angle)+rr*rr)))!=0 or s.trigsimp(numerator-(s.cos(angle)-rr))!=0:
            raise ArithmeticError('Correlated Poisson half-angle semantic identity failed')
        checks['correlated_Poisson_half_angle_sign_'+str(sign)]=True
    return dict(passed=True,exact_identities=checks,
        original_full_signed_five_density_cross_terms_retained=True,
        original_global_N_and_phase_held_derivative_convention_bound=True,
        actual_N_or_whole_transport_not_admitted_by_identities=True)


class CurrentLoopFunctionSources:
    def __init__(self):
        self.signed=signed.CurrentSignedInputJets();self.service=self.signed.service;self.family=self.signed.family
        checked=json.loads((HERE/signed.RECEIPT).read_bytes())
        if not checked['all_passed'] or not checked[signed.GATE] or checked['source_family']!=self.family or any(checked.get(k) for k in OPEN):
            raise ValueError('Checked original signed input/source scope required')
        self.service.bind_hashes(checked['input_hashes']);self.service.bind_hashes({signed.RECEIPT:sha(signed.RECEIPT)})
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name),PREFIX+'current_generic_shear_loop.py':sha(PREFIX+'current_generic_shear_loop.py'),
            PREFIX+'current_generic_shear_moment_recovery.py':sha(PREFIX+'current_generic_shear_moment_recovery.py')})
        manifest=json.loads((HERE/current.NAME).read_bytes());self.bounds=manifest['current_actual_loop_jet_log_bounds_by_chart']
        scales=manifest['current_actual_logarithmic_scales']
        self.scales=dict(eta=scales['selected_positive_eta_log'],d_star=scales['logarithmic_selected_positive_lower_constants']['d_star'])
        self.raw_prior=json.loads((HERE/(PREFIX+'current_generic_shear_source_bounds.json')).read_bytes())['current_original_source_log_bound_charts']
        self.raw_O3=json.loads((HERE/(PREFIX+'current_generic_shear_O3_sources.json')).read_bytes())['original_O3_quotient_log_norms']
        self.theorem=exact_theorem()
    def chart(self,chart):
        inputs=self.signed.saved(chart)
        packet=self.service.saved(chart) if chart in packets.CHARTS else self.signed.provider.saved(chart)
        V={(j,k):signed.ordinary_axial_coefficient(packet.velocity['axial'][j],k) for j,k in (ZERO,*FIRST)}
        result=build(inputs,V,self.scales,self.bounds[chart])
        raw=(self.raw_prior[chart] if chart in self.raw_prior else self.raw_O3[chart])['ordinary_mixed_source_log_norms']
        if chart=='O3_power' and inputs['original_O3_small_excess_definition']!='2*mu':raise ValueError('Canonical original constant power excess required for zero source branch')
        result['actual_source_frequency_majorants']=frequency_source_bounds(self.service.ctx,raw,self.bounds[chart],whole_q_flat=chart=='O3_power')
        return result
    def run(self):
        records={chart:self.chart(chart) for chart in self.bounds}
        (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(packets.encode(records),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
        self.service.bind_hashes({VIEWS:sha(VIEWS)})
        result=dict(source_family=self.family,function_source_cover_charts=len(records),
            actual_generic_loop_function_graph_views=VIEWS,
            source_graph_inventory={chart:dict(function_nodes=len(v['function_graph_nodes']),
                signed_five_increment_roots=list(v['five_signed_increment_rate_roots']),
                first_rate_derivatives=['y','Z'],same_global_N_and_phase=True) for chart,v in records.items()},
            exact_loop_modulation_and_signed_increment_theorem=self.theorem,
            **{GATE:True},**dict.fromkeys(OPEN,False),
            signed_current_point_loop_or_inverse_jets_installed=False,
            actual_changed_five_moment_transport_integrated=False,
            current_whole_N_selected=False,full_velocity_mixed4_or_stress_mixed3_admitted=False,
            source_graph_ancestor_constructors_called=False,
            scope='Actual source-function graph attachment for Section11 q/Poisson phase/inverse/A/B and E_N,V_N, five full signed increment rates and total first y/Z derivatives. No numerical point values selected from covers; changed own histories/repair/commonN/global cone/recursion remain open.',
            input_hashes=self.service.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        print('Actual generic loop and changed five signed source functions defined:17 charts',flush=True)
        return result


def run():return CurrentLoopFunctionSources().run()


if __name__=='__main__':run()
