"""Whole closed O3 cone for the actual finite-N modified source.

Weighted axial stress preserves microscopic shear factors and full native
stress sectors. The proof is uniform in all phases, Z and offset, not a
phase grid. Left O2 taper, quiet signed cones and the common global N remain
open. No original graph or fixed-frequency control vector is reconstructed.
"""
import functools
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_completed_tensor_error_majorants as errors
import lei_ren_part1_paper_compliant_current_O3_correlated_finite_N_shear as primitive
from lei_ren_part1_paper_compliant_current_original_cone_operator import original_cone_theorem

HERE,PREFIX,sha=errors.HERE,errors.PREFIX,errors.sha
NAME=PREFIX+'current_O3_modified_transition_cone.json'
RECEIPT=PREFIX+'current_O3_modified_transition_cone_check.json'
BASE_NAME=PREFIX+'current_O3_transition_direction.json'
BASE_RECEIPT=PREFIX+'current_O3_transition_direction_check.json'
MINIMUM_N=22
parameters,endpoints,upper=errors.parameters,errors.endpoints,errors.upper
OPEN=('current_modified_O2_whole_cone_certified','current_modified_quiet_signed_cone_certified',
    'common_N_modified_cones_certified','global_admissible_tensor_certified',
    'global_physical_velocity_interface_composition_certified','actual_coefficient_recursion_certified',
    'finite_energy_certified','global_flat_remainder_certified','full_corrected_NS_certified')


@functools.lru_cache(maxsize=1)
def exact_theorem():
    """Source-bound weighted cone and all-N monotone error transfer."""
    asts=errors.recovered.histories.SourceAST();checks={}
    def zero(name,a,b):
        diff=s.expand(a-b)
        if diff!=0 and s.cancel(diff)!=0:raise ArithmeticError('Modified cone identity: '+name)
        checks[name]=True
    mu=s.Symbol('mu',positive=True);root=s.Symbol('sqrtmu',positive=True)
    g,h,theta,zeta=s.symbols('g h theta weighted_axial',real=True)
    a=2+root**2*g;b=root*h;tz=zeta/root
    D=theta-b*tz/a;J=tz+b*theta/a;km=a-2+b*b/a
    H=a*theta-h*zeta;K=a*zeta+root**2*h*theta
    mhat=2*g+h*h+root**2*g*g
    W=(2*a-root**2*h*h)*theta*theta-2*a*h*theta*zeta-a*g*zeta*zeta
    zero('same_signed_direction_without_cutoff_division',H,a*D)
    zero('same_signed_weighted_cross',K,root*a*J)
    zero('same_factored_actual_primitive_excess',mhat,a*km/root**2)
    zero('same_full_signed_two_vector_quadratic',2*a*H*H-mhat*K*K,a**3*(2*D*D-km*J*J))
    zero('reduced_weighted_direction_polynomial',2*a*H*H-mhat*K*K,(a*a+root**2*h*h)*W)
    # Current profile derivatives retain physical Pstar in both components;
    # that common factor cancels in bs, with its original positive sign.
    P,U,Uy,Vy=s.symbols('Pstar Utheta_over_Pstar Utheta_y_over_Pstar Vhat_y',positive=True)
    zero('actual_physical_Pstar_cancels_in_signed_bs',2*(P*Vy)/(P*U),2*Vy/U)
    source=original_cone_theorem()
    if not source['passed']:raise ValueError('Original signed leading cone theorem required')
    # Source-bound finite-N shear, in dimensionless g/h units. No 2+tiny
    # subtraction is used to recover g. sigma remains the actual cutoff.
    chi,chiy,N,L,sig,S,C=s.symbols('chi chi_y N log_slope sigma sinphase cosphase',real=True)
    A=-mu*chi**2*S*C/(4*s.pi)
    gn=2*sig+chi**2*(C*C-S*S)+chi*chiy*S*C/(s.pi*N)
    hn=s.exp(-A/N)*(2*chi*S-(chi*L+chiy)*C/(s.pi*N))
    env=dict(c=SimpleNamespace(pi=s.pi,sqrt=s.sqrt,exp=s.exp,
        sin=lambda x:2*S*C if x==4*s.pi*s.Symbol('phase',real=True) else S,
        cos=lambda x:C*C-S*S if x==4*s.pi*s.Symbol('phase',real=True) else C),
        mu=mu,chi=[chi,chiy],phase=s.Symbol('phase',real=True),n=N,sig=sig,logEy=L,A=[A])
    for name in ('da','slowA','beta','slowB_over_E','aexcess','b'):
        env[name]=asts.evaluate(asts.expression('current_O3_finite_frequency_profiles','profile',name),env)
    zero('same_actual_finite_N_g_without_dividing_tiny_numbers',env['aexcess'],mu*gn)
    zero('same_actual_finite_N_h_with_original_bs_sign',env['b'],s.sqrt(mu)*hn)
    # Native stress order0 uses local velocity through y1 and history values
    # only. Their positive cap formulas increase with q=1/N.
    q=s.Symbol('inverse_N',positive=True)
    c0,c1,F0,F1=s.symbols('C0 C1 F0 F1',positive=True)
    cc=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),sqrt=s.sqrt,exp=s.exp,pi=s.pi)
    cutoff=(c0,c1,s.Integer(1),s.Integer(1),s.Integer(1))
    fs=[[F0]+[s.Integer(1)]*5,[F1]+[s.Integer(1)]*5]+[[s.Integer(1)]*6 for _ in range(3)]
    profiles=parameters.normalized_envelopes(cc,mu,1/q,cutoff,fs,s.Integer(1))
    exponent=mu*c0*c0*q/(8*s.pi);h1=mu*c0*c0/2+mu*c0*c1*q/(4*s.pi)
    targets=dict(theta0=F0*s.exp(exponent)*exponent,
        theta1=s.exp(exponent)*(F1*exponent+F0*h1),
        axial0=s.sqrt(mu)*F0*c0*q/(2*s.pi),
        axial1=s.sqrt(mu)*(2*s.pi*F0*c0+(F1*c0+F0*c1)*q)/(2*s.pi))
    for label,key in (('theta','theta_increment_over_Pstar_over_Ad_majorants'),
                      ('axial','modified_axial_over_Pstar_over_Ad_majorants')):
        for j in (0,1):
            name=label+str(j);zero('same_positive_local_cap_'+name,profiles[key][j][0],targets[name])
            derivative=s.diff(targets[name],q)
            if label=='theta':derivative=s.expand(derivative/s.exp(exponent))
            poly=s.Poly(s.expand(derivative),q)
            if not all(v.is_nonnegative for v in poly.all_coeffs()):raise ValueError('Nonmonotone native0 profile cap')
            checks['positive_inverse_N_derivative_'+name]=True
    model,_=errors.source_model()
    for label,parts in model['stress'].items():
        for name,p in parts.items():
            for _,monomial in errors.polynomial_terms(p['expression'],0,0,p['mode'][0]):
                for (fn,j,k),power in monomial:
                    if fn in ('du','V') and j>1 or fn in ('dm','dh','dk','de','dp') and j!=0 or fn=='dr':
                        raise ValueError('Higher growing N derivatives used in leading cone')
            checks['same_native0_monotone_source_dependencies_'+label+'/'+name]=True
    # All cumulative values are nonnegative multiples of q or min of two
    # increasing bounds K*q and Kswirl*q+Kkinetic*q².
    K1,K2=s.symbols('positive_K1 positive_K2',positive=True)
    for name,expr in (('history',K1*q),('kinetic_plus_swirl',K1*q+K2*q*q)):
        if s.diff(expr,q).is_positive is not True:raise ValueError('Nonmonotone history value')
        checks['positive_inverse_N_derivative_'+name]=True
    for target in ('initial','kinetic','swirl'):
        asts.expression('current_O3_recovered_error_majorants','modulation',target)
    # exp(x)<=1/(1-x), for 0<=x<1: derivative of (1-x)exp(x) is -x exp(x).
    x=s.Symbol('nonnegative_x',nonnegative=True)
    zero('positive_exponential_rational_upper_initial',((1-x)*s.exp(x)).subs(x,0),1)
    zero('positive_exponential_rational_upper_derivative',s.diff((1-x)*s.exp(x),x),-x*s.exp(x))
    mumax=s.Rational(1,1000);chi1=128;n0=MINIMUM_N
    G=2+s.Rational(chi1,6*n0)
    Hcap=(2+(s.Rational(1,2)+mumax+chi1)/(3*n0))/(1-mumax/(24*n0))
    if G>=3 or Hcap>=4:raise ValueError('Uniform supported shear bounds too large')
    checks['all_N_ge22_actual_absolute_g_less3']=True
    checks['all_N_ge22_actual_absolute_h_less4']=True
    aa,bb=s.symbols('actual_positive_a actual_signed_b',real=True)
    zero('actual_vs_excess_above_correlated_primitive_when_a_bounded',
        aa-2+bb*bb/aa-(aa-2+bb*bb/(2+3*mu)),bb*bb*(2+3*mu-aa)/(aa*(2+3*mu)))
    asts.expression('pre_pulse_mixed_C4','slope_mu','logU')
    hashes={**source['input_hashes'],**asts.hashes,Path(__file__).name:sha(Path(__file__).name)}
    name=PREFIX+'current_original_cone_operator.py';hashes[name]=sha(name)
    return dict(passed=True,identities=checks,input_hashes=hashes,
        stress_normalization='theta=Ttheta/[Pstar*sqrt(R/2)]; zeta=sqrt(mu)*Tz/[Pstar*sqrt(R/2)]',
        shear_normalization='a=2+mu*g; bs=sqrt(mu)*h; signed bs=+2*Uz_y/Utheta',
        leading_weighted_cone='a>0; mhat=2g+h²+mu*g²>0; H=a*theta-h*zeta>0; W=(2a-mu*h²)theta²-2a*h*theta*zeta-a*g*zeta²>0',
        all_phase_source_bounds=dict(minimum_integer_N=n0,absolute_g_upper=3,absolute_h_upper=4,
            proof='pi>=3; mu<=.001; chi<=1; chi_y<=128; sigma<=1/2 while chi is active'),
        positive_native0_caps_nonincreasing_for_all_integer_N_ge_minimum=True,
        original_quiet_negative_b_label_not_used_in_signed_cone=True)


class CurrentModifiedO3TransitionCone:
    def __init__(self,require_checked=True):
        self.errors=errors.CurrentCompletedTensorErrorMajorants();self.ctx=self.errors.ctx
        self.data=self.errors.data;self.mu=self.errors.mu;self.theorem=exact_theorem()
        self.primitive=primitive.CurrentCorrelatedFiniteNShear();self.hashes=dict(self.errors.hashes)
        read=parameters.numeric.transport.read_interval
        base=json.loads((HERE/BASE_NAME).read_bytes());receipt=json.loads((HERE/BASE_RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt['current_whole_O3_variable_transition_direction_bound_certified']:
            raise ValueError('Checked whole actual baseline direction required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed whole baseline cone source: '+name)
        for name,value in self.data['source']['accepted']['source_family'].items():
            if base.get(name)!=value:raise ValueError('Different baseline current family/pressure datum')
        self.base=base['current_whole_variable_transition_bounds']
        if not self.base['whole_closed_direction_and_positive_stress_certified'] or not self.base['complete_original_signed_stress_pressure_velocity_and_remainder_preserved']:
            raise ValueError('Incomplete original signed baseline vector')
        self.read=read
        c=self.ctx;floor=read(c,self.base['full_theta_uniform_lower'])
        self.logfloor=self.data['logAd']-c.mpf('.5')-self.mu/2+c.ln(floor)
        # This is the accepted bound for 2*mu*(Tz0/Ttheta0)^2,
        # valid on the closed transition including its old zero-shear seam.
        self.rho0=upper(c,c.sqrt(read(c,self.base['full_directional_expression_upper'])/2))
        self.error_rows=self.leading_error_rows(MINIMUM_N)
        self.proof=self.prove()
        self.hashes.update({**receipt['input_hashes'],BASE_RECEIPT:sha(BASE_RECEIPT),
            **self.primitive.hashes,primitive.RECEIPT:sha(primitive.RECEIPT),
            errors.RECEIPT:sha(errors.RECEIPT),**self.theorem['input_hashes']})
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed']:raise ValueError('Checked actual modified O3 cone required')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed modified closed O3 cone input: '+name)

    def leading_error_rows(self,N):
        c=self.ctx;q,caps=self.errors.inputs('O3',(0,1),N);result={}
        for label,parts in self.errors.model['stress'].items():
            result[label]={}
            for name,p in parts.items():
                cap=errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),self.errors.delta,0,p['mode'][0])[0,0]
                rp,pp,_,_=p['mode'];extra=self.data['log_mu']/2 if label=='axial' else c.mpf(0)
                partslog=dict(relative_radius=(c.mpf(str(rp))-c.mpf('.5'))*q['native_logR'],
                    relative_Pstar=(pp-1)*self.data['logP'],external_Ad=p['error_Ad_power']*self.data['logAd'],
                    weighted_axial=extra,coefficient=c.ln(cap),baseline_floor=-self.logfloor)
                result[label][name]=dict(source_mode=p['mode'],error_Ad_power=p['error_Ad_power'],
                    native_coefficient_absolute_cap=cap,relative_error_log_parts=partslog,
                    log_absolute_error_over_baseline_theta_lower=sum(partslog.values(),c.mpf(0)))
        return result

    def prove(self):
        c=self.ctx;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Whole modified O3 cone unresolved: '+name)
            positive[name]=value
        for label,parts in self.error_rows.items():
            for name,p in parts.items():require('native_error_log_gap/'+label+'/'+name,
                -1010-p['log_absolute_error_over_baseline_theta_lower'])
            require('sum_of_native_errors_below_exp_minus1000/'+label,10-c.ln(len(parts)))
        eps=c.exp(-1000);rho=(self.rho0+eps)/(1-eps)
        require('positive_modified_theta_reserve',1-eps)
        require('weighted_axial_ratio_below_exp_minus390',c.exp(-390)-rho)
        # Fixed rational caps avoid assembling 2+microscopic mu in an interval
        # and then subtracting 2. Actual a is kept in its source factorization.
        alo=c.mpf('1.997');ahi=c.mpf('2.003');G=c.mpf(3);H=c.mpf(4);mumax=c.mpf('.001')
        direction=alo-H*rho
        W=2*alo-mumax*H*H-2*ahi*H*rho-ahi*G*rho*rho
        require('weighted_H_over_modified_theta',direction)
        require('weighted_direction_polynomial_over_modified_theta_squared',W)
        Dlower=direction/ahi;Qlower=alo*alo*W/ahi**3
        require('original_signed_D_over_modified_theta',Dlower)
        require('original_signed_Q_over_modified_theta_squared',Qlower)
        sig=primitive.sigma_jets(c,c.mpf('.25'))[0]
        f0=c.mpf(min(endpoints(c.mpf('.25'))[0],endpoints(2*sig)[0]))
        require('closed_O3_correlated_shear_floor_over_mu',f0)
        require('whole_actual_vs_minus2_lower',self.mu*f0)
        require('whole_a_positive_lower',alo)
        return dict(minimum_integer_N=MINIMUM_N,all_finite_integer_N_at_least_minimum=True,
            current_O3_offset_domain=(0,1),whole_Z_domain=(-1,1),all_actual_phases=True,
            exact_actual_phase='N*offset; same translated N*logR source',positive_margins=positive,
            baseline_theta_log_lower=self.logfloor,baseline_weighted_axial_ratio_upper=self.rho0,
            relative_theta_error_upper=eps,relative_weighted_axial_error_upper=eps,
            modified_weighted_axial_ratio_upper=rho,actual_a_lower=alo,actual_a_upper=ahi,
            actual_g_absolute_upper=G,actual_h_absolute_upper=H,
            actual_vs_minus2_lower=self.mu*f0,weighted_H_over_theta_lower=direction,
            original_signed_D_over_theta_lower=Dlower,original_signed_Q_over_theta_squared_lower=Qlower,
            reduced_weighted_quadratic_over_theta_squared_lower=W,
            modified_theta_log_lower=self.logfloor+c.ln(1-eps),
            original_pressure_energy_moments_and_all_cross_terms_retained=True,
            no_cutoff_or_axial_coordinate_division_used=True,
            resolved_signed_coefficients_or_control_midpoints_used=False,
            current_modified_closed_O3_signed_two_vector_cone_certified=True,
            **{key:False for key in OPEN})

    def query(self,offset,N):
        parameters.positive_integer_N(N);c=self.ctx;t=c.mpf(offset);lo,hi=endpoints(t)
        if N<MINIMUM_N or lo<0 or hi>1 or not all(mp.isfinite(v) for v in (lo,hi)):
            raise ValueError('Current closed O3 offset subset[0,1], finite integer N>=22 required')
        current=self.primitive.query(t,N)
        # g/h are ratios of factored source quantities, not differences of
        # rounded 2+mu terms. Their point/box enclosures are diagnostic only;
        # uniform cone proof above is independent of phase-box width.
        g=current['finite_N_theta_shear_excess']/self.mu
        h=current['finite_N_axial_shear']/c.sqrt(self.mu)
        return dict(offset=t,finite_integer_N=N,actual_phase=c.mpf(N)*t,
            signed_current_g_enclosure=g,signed_current_h_enclosure=h,
            actual_bs_convention='+2*Uz_y/Utheta=sqrt(mu)*h',
            source_correlated_primitive_lower=current['source_certified_joint_margin_lower_envelope'],
            uniform_whole_source_certificate=self.proof,source_family=current['source_family'],
            current_modified_closed_O3_signed_two_vector_cone_certified=True,
            **{key:False for key in OPEN})


def run():
    field=CurrentModifiedO3TransitionCone(require_checked=False)
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_weighted_source_cone_theorem=field.theorem,whole_modified_closed_O3_cone=field.proof,
        complete_leading_stress_error_log_ledger=field.error_rows,
        examples={name:field.query(t,N) for name,t,N in (
            ('closed_whole',(0,1),22),('old_zero_shear_seam','0',22),
            ('right_flat_taper',('.25','.5'),22),('right_flat_edge','.5',22),
            ('terminal','1',22),('current_frequency',('0','.5'),10**12))},
        current_modified_closed_O3_signed_two_vector_cone_certified=True,
        scoped_source_cone_does_not_increment_unmodified_registry_counts=True,
        **{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Modified whole closed O3 cone generated; all phases/Z/offset, all integer N>=22',flush=True)
    return result


if __name__=='__main__':run()
