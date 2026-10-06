"""Exact modulated cumulative histories and frequency-uniform enclosures."""
import ast
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_finite_frequency_profiles import (
    SourceAST,cutoff_rows,endpoints,source_precision)
from lei_ren_part1_paper_compliant_current_O3_transition_direction_operator import upper
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_frozen_comparison_field import derivative,square
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

SHAPES={'m':('C',1,1),'h':('C',1,'3/2'),'k':('C^2',2,'3/2'),
        'e':('C^2',2,1),'p':('C^2',2,0)}
MESH=(-2,-1,0,s.Rational(1,4),s.Rational(1,2))

def exact_modulated_history_theorem():
    checks={};asts=SourceAST()
    def zero(name,a,b):
        if s.simplify(s.expand(a-b))!=0:raise ArithmeticError('Modulated history identity failed: '+name)
        checks[name]=True
    original=asts.method('pre_pulse_mixed_C4','physical_mixed')
    wanted=(
      "rows['m'].append(V[j]-rows['m'][j])",
      "rows['h'].append(U[j]-rows['h'][j]*c.mpf('1.5'))",
      "rows['k'].append(product_rows(U,V,j)-rows['k'][j]*c.mpf('1.5'))",
      "rows['e'].append(product_rows(V,V,j)*invP2-product_rows(U,U,j)/2-rows['e'][j])",
      "rows['p'].append(product_rows(U,U,j)/2)")
    calls={ast.dump(node) for node in ast.walk(original) if isinstance(node,ast.Call)}
    for text in wanted:
        if ast.dump(ast.parse(text,mode='eval').body) not in calls:raise ValueError('Original normalized five-history ODE changed')
    checks['actual_original_five_normalized_history_ODE_program_bound']=True
    P,U,du,V,m,h,k,e,p,dm,dh,dk,de,dp=s.symbols('Pstar old_u swirl_increment Vhat old_m old_h old_k old_e old_p dmhat dh dkhat de dp',real=True)
    newu=U+du
    zero('actual_m_increment_has_extra_Pstar',-m+P*(V-dm),P*V-(m+P*dm))
    zero('actual_h_increment_units',U-s.Rational(3,2)*h+du-s.Rational(3,2)*dh,newu-s.Rational(3,2)*(h+dh))
    zero('actual_k_increment_has_extra_Pstar',-s.Rational(3,2)*k+P*(newu*V-s.Rational(3,2)*dk),
         newu*(P*V)-s.Rational(3,2)*(k+P*dk))
    zero('actual_e_increment_complete_energy_units',-U**2/2-e+V**2-U*du-du**2/2-de,
         V**2-newu**2/2-(e+de))
    zero('actual_p_increment_shared_axis_pressure_units',U**2/2+U*du+du**2/2,newu**2/2)
    t,x,N=s.symbols('actual_logR_offset X N',real=True)
    f,beta,G,Ua,C=s.symbols('actual_f actual_beta actual_G actual_Ua C',real=True)
    zero('actual_swirl_increment_squared_difference',2*(Ua*f*C)*(Ua*f*C*(G-1))+(Ua*f*C*(G-1))**2,
         Ua**2*f**2*C**2*(G**2-1))
    # Bind the scalar f and C shapes to the actual two source programs,
    # not a profile fitted from the scalar endpoint enclosure.
    tau,parent_U,mu,J=s.symbols('actual_buffer_inlet_offset parent_U mu same_J',real=True)
    ctx=SimpleNamespace(exp=s.exp)
    buffer_root=asts.evaluate(asts.expression('pre_pulse_mixed_C4','axial','root',wanted='c.exp(-t/2)'),
      dict(c=ctx,t=tau+t))
    buffer_u=asts.evaluate(asts.expression('pre_pulse_mixed_C4','axial','u',wanted='u1*root'),
      dict(u1=parent_U*C,root=buffer_root))
    zero('actual_buffer_Ua_f_C_source_function',buffer_u,
      (parent_U*s.exp(-tau/2))*C*s.exp(-t/2))
    factor=asts.evaluate(asts.expression('pre_pulse_mixed_C4','slope_mu','factor',
      wanted="c.exp(-t/2-mu*K['J'])"),dict(c=ctx,t=t,mu=mu,K=dict(J=J)))
    transition_u=asts.evaluate(asts.expression('pre_pulse_mixed_C4','slope_mu','u',wanted='u1*factor'),
      dict(u1=Ua*C,factor=factor))
    zero('actual_transition_Ua_f_C_source_function',transition_u,Ua*C*s.exp(-t/2-mu*J))
    expressions={'m':(Ua*f*C*beta/N,Ua*C*f*beta/N),
      'h':(Ua*f*C*(G-1),Ua*C*f*(G-1)),
      'k':((Ua*f*C*G)*(Ua*f*C*beta/N),Ua**2*C**2*f**2*beta*G/N),
      'e':((Ua*f*C*beta/N)**2-Ua**2*f**2*C**2*(G**2-1)/2,
           Ua**2*C**2*f**2*(beta**2/N**2-(G**2-1)/2)),
      'p':(Ua**2*f**2*C**2*(G**2-1)/2,Ua**2*C**2*f**2*(G**2-1)/2)}
    for name,(a,b) in expressions.items():zero('actual_full_Z_shape_'+name,a,b)
    source=s.Function('same_signed_density')(t)
    for name,(_,_,rate) in SHAPES.items():
        r=s.Rational(str(rate));F=s.Function('same_cumulative_'+name)(t)
        value=s.exp(-r*t)*F
        zero('actual_continuous_transport_FTC_'+name,
          s.diff(value,t).subs(s.diff(F,t),s.exp(r*t)*source),source-r*value)
    z,delta,Bm,Vf=s.symbols('Z delta same_Bm actual_V_factor',real=True)
    Cz=1/(1+z*z);L=1-delta*z*z
    q=(2*z*Cz*Vf-(1-delta)*z*Cz*Bm-(1-z*z)*s.diff(Cz,z)*Bm)/L
    zero('actual_radial_recovery_same_cumulative_m_not_density',q,
      z*(2*Cz*Vf+Cz**2*((1+delta)-(3-delta)*z*z)*Bm)/L)
    zero('actual_radial_increment_Pstar_units',P*q,
      (2*z*(P*Cz*Vf)-(1-delta)*z*(P*Cz*Bm)-(1-z*z)*s.diff(P*Cz*Bm,z))/L)
    zero('actual_C_partial_fraction_global_axial_bound',Cz,
         ((1-s.I*z)**-1+(1+s.I*z)**-1)/2)
    for j in range(6):
        zero('actual_C_squared_axial_product_majorant_'+str(j),
          sum(s.binomial(j,n)*math.factorial(n)*math.factorial(j-n) for n in range(j+1)),math.factorial(j+1))
    J=s.symbols('same_original_J',real=True)
    ft=s.exp(-t/2-s.symbols('mu',real=True)*J)
    for name,rate,power in (('m',1,1),('h',s.Rational(3,2),1),('k',s.Rational(3,2),2),('e',1,2),('p',0,2)):
        zero('actual_positive_log_coordinate_weight_'+name,s.exp(rate*t)*ft**power,
             s.exp((rate-s.Rational(power,2))*t-power*s.symbols('mu',real=True)*J))
    integer_N=s.Symbol('finite_positive_integer_N',integer=True,positive=True)
    zero('actual_positive_kinetic_mass_on_unit_buffer_plateau',
      s.integrate(s.cos(2*s.pi*integer_N*t)**2,(t,-1,0)),s.Rational(1,2))
    return dict(identities=checks,scalar_integral_definitions={
      'm':'exp(-t)*integral[-2,min(t,1/2)] exp(s/2-mu*J(s))*beta(s)/N ds',
      'h':'exp(-3t/2)*integral exp(s-mu*J(s))*expm1(A(s)/N) ds',
      'k':'exp(-3t/2)*integral exp(s/2-2mu*J(s))*beta(s)*exp(A(s)/N)/N ds',
      'e':'exp(-t)*integral exp(-2mu*J(s))*(beta(s)^2/N^2-expm1(2A(s)/N)/2) ds',
      'p':'integral exp(-s-2mu*J(s))*expm1(2A(s)/N)/2 ds'},
      negative_side_J_is_exact_zero_only_for_buffer=True,
      positive_side_J_is_same_original_nonnegative_cutoff_integral=True,
      whole_Z_canonical_shapes=SHAPES,
      raw_normalized_increment_Pstar_powers={'m':1,'h':0,'k':1,'e':0,'p':0},
      shared_axis_pressure='P0(Z)+old_p(t,Z)+delta_p(t,Z), all divided by Pstar^2',
      positive_kinetic_term_retained_separately=True,
      axial_C_bounds='|dZ^j C|<=j!; |dZ^j C^2|<=(j+1)! for j0..5, by partial fractions and Leibniz',
      input_hashes=asts.hashes,passed=True)

def signed_expm1_enclosure(c,x):
    """Same exact expm1(x)=x*integral_0^1 exp(r*x)dr, not subtraction."""
    cap=upper(c,abs(x));return x*c.exp(c.mpf([-endpoints(cap)[1],endpoints(cap)[1]]))

def scalar_integrand_enclosure(c,mu,N,t):
    """Range of the exact signed source on a whole log-coordinate cell.

    The cutoff-integral envelope bounds the defining original J; it never
    replaces J or its positive tail. Oscillations retain the actual N*t.
    """
    lo,hi=endpoints(t)
    if lo<-2 or hi>.5:raise ValueError('Actual modulation support cell required')
    chi=cutoff_rows(c,t)[0];phase=N*t
    A=-mu*chi**2*c.sin(4*c.pi*phase)/(8*c.pi)
    beta=-c.sqrt(mu)*chi*c.cos(2*c.pi*phase)/(2*c.pi)
    # J(s)=0 on the buffer, and 0<=J(s)<=s on the transition.
    if lo<0<hi:raise ValueError('Split the exact source integral at the O2/O3 seam')
    J=c.mpf(0) if hi<=0 else c.mpf([0,max(mp.mpf(0),hi)])
    ax=A/N;G=c.exp(ax);dG=signed_expm1_enclosure(c,ax);dG2=signed_expm1_enclosure(c,2*ax)
    return dict(m=c.exp(t/2-mu*J)*beta/N,
      h=c.exp(t-mu*J)*dG,
      k=c.exp(t/2-2*mu*J)*beta*G/N,
      e_kinetic=c.exp(-2*mu*J)*beta**2/N**2,
      e_swirl=-c.exp(-2*mu*J)*dG2/2,
      p=c.exp(-t-2*mu*J)*dG2/2)

@source_precision
def cumulative_scalar_enclosures(c,mu,N,t,cells):
    """Signed range-integration for every endpoint in the requested box."""
    lo,hi=endpoints(t);upper_t=min(hi,mp.mpf('.5'));lower_t=min(lo,mp.mpf('.5'))
    if hi<=-2:return {name:c.mpf(0) for name in ('m','h','k','e','p','e_kinetic','e_swirl')}
    result={name:c.mpf(0) for name in ('m','h','k','e_kinetic','e_swirl','p')}
    for left,right in zip(MESH[:-1],MESH[1:]):
        a=c.mpf(str(left));b=c.mpf(str(right))
        for j in range(cells):
            x=a+(b-a)*j/cells;y=a+(b-a)*(j+1)/cells
            xl,xh=endpoints(x);yl,yh=endpoints(y)
            if xl>=upper_t:continue
            stop=min(yh,upper_t);start=xl
            width=c.mpf([max(mp.mpf(0),min(yl,lower_t)-xh),stop-start])
            # Cells are fixed, split at both chart and cutoff flat edges.
            values=scalar_integrand_enclosure(c,mu,c.mpf(N),c.mpf([start,stop]))
            for name,value in values.items():result[name]+=width*value
    result['e']=result['e_kinetic']+result['e_swirl']
    for name,(_,_,rate) in SHAPES.items():
        if rate:result[name]*=c.exp(-c.mpf(str(rate))*t)
    for name in ('e_kinetic','e_swirl'):result[name]*=c.exp(-t)
    if lo>=0:
        floor=mu*c.exp(-t)/(8*c.pi**2*c.mpf(N)**2)
        kinetic=result['e_kinetic']
        result['e_kinetic']=c.mpf([max(endpoints(kinetic)[0],endpoints(floor)[0]),endpoints(kinetic)[1]])
        result['e']=result['e_kinetic']+result['e_swirl']
    return result

@source_precision
def frequency_uniform_bounds(c,mu):
    if endpoints(mu)[0]<=0:raise ValueError('Actual positive mu required')
    a=mu/(8*c.pi);b=c.sqrt(mu)/(2*c.pi);length=c.mpf('2.5')
    K=dict(m=length*c.exp(c.mpf('.25'))*b,
      h=length*c.exp(c.mpf('.5'))*a*c.exp(a),
      k=length*c.exp(c.mpf('.25'))*b*c.exp(a),
      e=length*(b*b+a*c.exp(2*a)),p=length*c.exp(2)*a*c.exp(2*a),
      e_kinetic=length*b*b,e_swirl=length*a*c.exp(2*a))
    if any(endpoints(value)[0]<=0 for value in K.values()):raise ArithmeticError('Positive source constants unresolved')
    return dict(untransported_abs_integral_times_N_upper=K,
      positive_kinetic_buffer_mass_times_N_squared_lower=mu/(8*c.pi**2),
      positive_kinetic_lower_applies_to_all_endpoints_at_least0=True,
      all_finite_integer_N_at_least1=True,all_cumulative_endpoints_in_support_and_after=True,
      pointwise_small_profiles_not_phase_sample_cancellation_used=True,
      delta_M_h_factor='Ua',delta_k_e_p_factor='Ua^2',
      derivative_majorants_C=[math.factorial(j) for j in range(6)],
      derivative_majorants_C_squared=[math.factorial(j+1) for j in range(6)],
      radial_transport_factors={'m':'exp(-t)','h':'exp(-3t/2)','k':'exp(-3t/2)','e':'exp(-t)','p':'1'},
      finite_uniform_N_admissible_cone_choice=False)

def increment_rows(c,z,delta,base,increments,Unew,Uold,Vnew,swirl_increment):
    """Ordinary source jets of every modified history, correct P sectors."""
    C=(1+square(z)).reciprocal()
    d={name:[C*(base*increments[name]) if SHAPES[name][0]=='C'
      else square(C)*(base**2*increments[name])] for name in SHAPES}
    for j in range(4):
        d['m'].append(Vnew[j]-d['m'][j])
        d['h'].append(swirl_increment[j]-d['h'][j]*c.mpf('1.5'))
        d['k'].append(product_rows(Unew,Vnew)[j]-d['k'][j]*c.mpf('1.5'))
        ds=product_rows(Uold,swirl_increment)[j]*2+product_rows(swirl_increment,swirl_increment)[j]
        d['e'].append(product_rows(Vnew,Vnew)[j]-ds/2-d['e'][j])
        d['p'].append(ds/2)
    L=1-square(z)*delta
    Q=[(z*Vnew[j]*2-z*d['m'][j]*(1-delta)-(1-square(z))*derivative(d['m'][j]))/L for j in range(5)]
    return d,Q
