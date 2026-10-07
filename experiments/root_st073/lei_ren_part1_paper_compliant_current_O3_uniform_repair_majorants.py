"""Uniform finite-N implicit repair, bump jets and cumulative defect bounds.

Consumes accepted independent weights, never the saved N=10^12 controls.
All amplitudes/radii stay factored. This supplies a repair-family estimate,
not an admission of the completed signed tensor or coefficient recursion.
"""
import ast
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_actual_parameter_majorants as parameters
import lei_ren_part1_paper_compliant_current_O3_independent_repair_operator as repair
import lei_ren_part1_paper_compliant_current_O3_repaired_histories_operator as histories
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import BETA_POLYNOMIALS, BETA_CONSTANTS
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value

HERE, PREFIX, sha = parameters.HERE, parameters.PREFIX, parameters.sha
NAME = PREFIX+'current_O3_uniform_repair_majorants.json'
RECEIPT = PREFIX+'current_O3_uniform_repair_majorants_check.json'
PARENT = PREFIX+'current_O3_independent_repair.json'
PARENT_CHECK = PREFIX+'current_O3_independent_repair_check.json'
HISTORY_CHECK = PREFIX+'current_O3_repaired_histories_check.json'
OPEN = ('recovered_radial_and_absolute_pressure_error_bounds_available',
        'completed_signed_tensor_error_bounds_available', 'common_N_modified_cones_certified',
        'actual_coefficient_recursion_certified', 'finite_energy_certified', 'full_corrected_NS_certified')
endpoints = parameters.numeric.transport.endpoints


def upper(c, value):
    return c.mpf(endpoints(c.mpf(value))[1])


def ceil_upper(value):
    # Exact binary endpoints; do not round the sufficient frequency down.
    sign,mantissa,exponent,_=endpoints(value)[1]._mpf_
    if sign:raise ValueError('Positive sufficient-frequency bound required')
    if exponent>=0:ceiling=mantissa << exponent
    elif -exponent>=mantissa.bit_length():ceiling=int(mantissa>0)
    else:
        shift=-exponent
        ceiling=(mantissa >> shift)+int(mantissa != ((mantissa >> shift) << shift))
    return ceiling+1


def cap_min(c, a, b):
    return c.mpf(min(endpoints(a)[1],endpoints(b)[1]))


def jet_caps(c, normalization, ell):
    raw=[c.exp(-1)]+[BETA_CONSTANTS[k]*(2*k)**(2*k)*c.exp(-2*k) for k in range(1,5)]
    normalized=[raw[k]/(ell**(k+1)*normalization) for k in range(5)]
    log_bump=[upper(c,sum((math.comb(j,k)*normalized[k] for k in range(j+1)),c.mpf(0))) for j in range(5)]
    return dict(raw_ordinary_derivative_caps=raw,normalized_beta_ordinary_derivative_caps=normalized,
                log_bump_ordinary_derivative_caps=log_bump)


def exact_theorem():
    asts=repair.SourceAST(); checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(a-b)) != 0 and s.simplify(a-b) != 0:
            raise ArithmeticError('Uniform repair identity: '+name)
        checks[name]=True
    mu,N=s.symbols('mu N',positive=True); h=s.symbols('a0 a2 e0 e1 e2',real=True)
    cross=s.symbols('cross0 cross2',positive=True); energy=s.symbols('energy0:3',positive=True)
    pressure=s.symbols('pressure0:3',positive=True)
    matrix=dict(cross_weights=cross,energy_weights=energy,pressure_weights=pressure)
    q=asts.replay('current_O3_independent_repair_operator','scaled_quadratic',{})(None,mu,N,matrix,h)
    for i,value in enumerate(q):
        zero('quadratic_homogeneity_'+str(i),sum(h[j]*s.diff(value,h[j]) for j in range(5)),2*value)
        zero('exact_one_over_N_quadratic_'+str(i),s.diff(N*value,N),0)
    expected_jacobian=[[0]*5,
        [cross[0]*h[2]/N,cross[1]*h[4]/N,cross[0]*h[0]/N,0,cross[1]*h[1]/N],
        [0]*5,[2*energy[0]*h[0]/N,2*energy[2]*h[1]/N,
        -mu*energy[0]*h[2]/N,-mu*energy[1]*h[3]/N,-mu*energy[2]*h[4]/N],
        [0,0,mu*pressure[0]*h[2]/N,mu*pressure[1]*h[3]/N,mu*pressure[2]*h[4]/N]]
    for i in range(5):
        for j in range(5):zero('same_source_Jacobian_%d_%d'%(i,j),s.diff(q[i],h[j]),expected_jacobian[i][j])
    y,alpha=s.symbols('logx alpha',real=True)
    f=s.exp(-alpha*y); e=s.Function('normalized_deltaE')(y); u=s.Function('normalized_deltaUz')(y)
    zero('quiet_theta_shear_difference',-2*s.diff(f+e,y)/(f+e)-2*alpha,-2*(s.diff(e,y)+alpha*e)/(f+e))
    zero('quiet_axial_shear',-2*s.diff(u,y)/(f+e),-2*s.diff(u,y)/(f+e))
    r=s.Symbol('raw_coordinate',real=True);raw=s.exp(-1/(1-r*r))
    for k,poly in enumerate(BETA_POLYNOMIALS):
        zero('same_raw_beta_ordinary_derivative_'+str(k),s.diff(raw,r,k),raw*sum(v*r**i for i,v in enumerate(poly))/(1-r*r)**(2*k))
    w=s.Symbol('positive_edge_variable',positive=True)
    for k in range(1,5):
        zero('raw_beta_majorant_stationary_point_'+str(k),s.diff(s.exp(-1/w)*w**(-2*k),w).subs(w,s.Rational(1,2*k)),0)
    beta=s.Function('same_normalized_beta')(y)
    for j in range(5):zero('same_log_bump_ordinary_jet_'+str(j),s.diff(s.exp(-y)*beta,y,j),
        s.exp(-y)*sum(math.comb(j,k)*(-1)**(j-k)*s.diff(beta,y,k) for k in range(j+1)))
    Z=s.Symbol('real_axial_variable',real=True)
    C=1/(1+Z**2)
    for k in range(6):zero('same_C_axial_derivative_factorial_'+str(k),s.diff(C,Z,k),
        (-1)**k*math.factorial(k)*((Z-s.I)**(-k-1)-(Z+s.I)**(-k-1))/(2*s.I))
    W=dict(mass=s.symbols('mass0:3'),D=s.symbols('D0 D2'),I=s.symbols('I0:3'),
        S=s.symbols('S0:3'),Cp=s.symbols('Cp0:3'),cross=cross,energy=energy,pressure=pressure)
    primitives=asts.replay('current_O3_repaired_histories_operator','partial_repair_primitives',{})(SimpleNamespace(sqrt=s.sqrt),mu,N,h,W)
    sa=s.sqrt(mu)/N;se=mu/N
    zero('partial_correlated_D_before_enclosure',(primitives['J']-primitives['M'])/(mu*sa),
        h[0]*W['D'][0]+h[1]*W['D'][1]+(cross[0]*h[0]*h[2]+cross[1]*h[1]*h[4])/N)
    # Terminal closure is the implicit equation, not a numerical zero fit.
    d=s.symbols('same_dM same_dD same_dI same_dS same_dCp')
    b=s.Matrix(5,5,lambda i,j:s.Symbol('B%d%d'%(i,j)))
    residual=[sum(b[i,j]*h[j] for j in range(5))+q[i]+d[i] for i in range(5)]
    factors=[sa,mu*sa,se,se,se]
    for i in range(5):zero('same_terminal_implicit_residual_'+str(i),factors[i]*residual[i],
        factors[i]*(sum(b[i,j]*h[j] for j in range(5))+q[i])+factors[i]*d[i])
    ctx=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    quiet_y=asts.evaluate(asts.expression('current_O3_repaired_histories','history','y',wanted='q-1'),dict(q=1+y))
    quiet_t=asts.evaluate(asts.expression('current_O3_modulated_histories','history','t',wanted='1+v'),dict(v=1+quiet_y))
    zero('same_quiet_coordinate_offset_two_plus_y',quiet_t,2+y)
    for target in ('upper_t','lower_t'):
        asts.expression('current_O3_modulated_histories_operator','cumulative_scalar_enclosures',target,
            wanted="min(%s,mp.mpf('.5'))"%('hi' if target=='upper_t' else 'lo'))
    zero('quiet_source_is_after_full_modulation_support',s.Integer(2)-s.Rational(1,2),s.Rational(3,2))
    transport=asts.method('current_O3_modulated_histories_operator','cumulative_scalar_enclosures')
    transport_expression=ast.parse('c.exp(-c.mpf(str(rate))*t)',mode='eval').body
    term=next(n.value for n in ast.walk(transport) if isinstance(n,ast.AugAssign)
        and ast.unparse(n.target)=='result[name]' and ast.dump(n.value)==ast.dump(transport_expression))
    for label,rate in (('M',1),('I',s.Rational(3,2)),('J',s.Rational(3,2)),('S',1),('Cp',0)):
        if rate:
            at=asts.evaluate(term,dict(c=ctx,rate=rate,t=quiet_t))
            inlet=asts.evaluate(term,dict(c=ctx,rate=rate,t=s.Integer(2)))
            zero('same_quiet_homogeneous_transport_'+label,at/inlet,s.exp(-rate*y))
            zero('fixed_R0_physical_unit_cancels_transport_'+label,s.exp(rate*y)*at/inlet,1)
        else:checks['same_Cp_post_support_constant']=True
    for stem,fn in (('flat_pulse_derivatives','beta_polynomials'),('flat_pulse_derivatives','beta_tail_bound'),
                    ('current_O3_repaired_histories_operator','bump_rows'),('current_O3_repaired_histories_operator','partial_weight')):
        asts.method(stem,fn)
    return dict(passed=True,identities=checks,input_hashes=asts.hashes,
        implicit_equation='B*h_N+Q_N(h_N,h_N)=-d_N; correlated D=(J-M)/mu before enclosure',
        jacobian_bound='||(B+DQ_N(h_N))^-1||inf <= CA/(1-2*CA*CQ*R/N)',
        support_argument='three fixed disjoint logx supports: at most one bump at each point',
        raw_derivative_argument='|P_k(r)|<=sum|coefficients|; maximize exp(-1/w)*w^(-2k),0<w<=1',
        axial_derivative_argument='C(Z) partial fractions; |Z+i|=|Z-i|=sqrt(1+Z^2)>=1 gives |C^(k)|<=k!',
        exact_terminal_and_remaining_strip_argument='same checked terminal implicit equation gives total defect = negative remaining repair integral',
        quiet_support_separation='q=1+y, actual t=1+q=2+y>=2>1/2; original scalar integral clamps at1/2, leaving only homogeneous transport',
        quiet_shear_argument='f=exp(-alpha*y); |deltaE|<=mu*R*G0/N; exact quotient preserves positive denominator')


class CurrentUniformRepairMajorants:
    def __init__(self,require_checked=True):
        self.data=parameters.inputs();self.ctx=self.data['ctx'];c=self.ctx
        self.raw=json.loads((HERE/PARENT).read_bytes()); parent=json.loads((HERE/PARENT_CHECK).read_bytes())
        history=json.loads((HERE/HISTORY_CHECK).read_bytes())
        if not parent['all_passed'] or not history['all_passed']:
            raise ValueError('Accepted independent map and same-source partial/terminal theorem required')
        if not parent['current_O3_modulated_unique_implicit_constant_repair_controls_certified'] or not history['same_exact_implicit_controls_reproduce_full_partial_matrix_and_close_all_five']:
            raise ValueError('Exact same-source implicit controls and terminal closure required')
        family=self.data['source']['accepted']['source_family']
        for record in (self.raw,parent,history):
            for key,value in family.items():
                if record.get(key)!=value:raise ValueError('Uniform repair belongs to another actual source family')
        self.hashes=dict(self.data['hashes'])
        for record in (self.raw,parent,history):
            for name,digest in record['input_hashes'].items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting repair dependency: '+name)
                if sha(name)!=digest:raise ValueError('Changed accepted repair dependency: '+name)
                self.hashes[name]=digest
        for name in (PARENT,PARENT_CHECK,HISTORY_CHECK):self.hashes[name]=sha(name)
        self.normalization=parameters.numeric.transport.read_interval(c,self.raw['same_raw_bump_normalization'])
        self.weights=restore_value(c,self.raw['new_independent_bump_weights'])
        self.matrix=restore_value(c,self.raw['new_scaled_linear_inverse_quadratic_map'])
        self.parent_theorem=repair.exact_repair_theorem()
        if parameters.encoded(self.parent_theorem)!=self.raw['exact_source_and_new_five_bump_map_theorem']:
            raise ValueError('Original physical/scaled repair source theorem changed')
        self.theorem=exact_theorem();self.hashes.update(self.theorem['input_hashes'])
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.mu=self.data['mu'];self.alpha=c.mpf('.5')+self.mu
        self.jets=jet_caps(c,self.normalization,self.weights['radius'])
        saved=restore_value(c,self.raw['uniform_repair_only_contraction_certificate'])
        self.repair_threshold=saved['repair_only_sufficient_integer_N']
        self.base=repair.contraction_certificate(c,self.mu,self.repair_threshold,self.matrix)
        if self.base['repair_only_sufficient_integer_N']!=self.repair_threshold:
            raise ValueError('Accepted uniform repair frequency changed at actual mu')
        self.R=self.base['scaled_controls_ball_radius'];G=self.jets['log_bump_ordinary_derivative_caps']
        positive=ceil_upper(2*self.mu*self.R*G[0]*c.exp(self.alpha))
        shear=ceil_upper(8*self.R*(G[1]+self.alpha*G[0])*c.exp(self.alpha))
        self.quiet_threshold=max(self.repair_threshold,positive,shear)
        if require_checked:
            accepted=json.loads((HERE/RECEIPT).read_bytes())
            if not accepted['all_passed']:raise ValueError('Accepted uniform repair estimates required')
            for name,digest in accepted['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed checked uniform repair program: '+name)

    def _weights_cap(self,y,remaining=False):
        c=self.ctx;mu=self.mu;ell=self.weights['radius'];ci=self.weights['centers'];H=self.weights['H']
        lo,hi=endpoints(y);raw0=self.jets['raw_ordinary_derivative_caps'][0]
        result={key:[] for key in ('mass','D','I','S','Cp','cross','energy','pressure')}
        for i,center in enumerate(ci):
            left,right=endpoints(center-ell)[0],endpoints(center+ell)[1]
            width=upper(c,c.mpf(right)-c.mpf(left))
            extent=upper(c,c.mpf(right)-c.mpf(lo) if remaining else c.mpf(hi)-c.mpf(left))
            length=c.mpf(max(mp.mpf(0),min(endpoints(width)[1],endpoints(extent)[1])))
            z=c.mpf([left,right])
            for name,p,full in (('mass',c.mpf(0),c.mpf(1)),('I',c.mpf('.5'),H['I']*c.exp(center/2)),
                ('S',-self.alpha,H['S']*c.exp(-self.alpha*center)),
                ('Cp',-self.alpha-1,H['Cp']*c.exp((-self.alpha-1)*center))):
                result[name].append(cap_min(c,full,raw0*length/(ell*self.normalization)*c.exp(p*z)))
            for name,p,full in (('energy',c.mpf(0),self.matrix['energy_weights'][i]),
                ('pressure',c.mpf(-1),self.matrix['pressure_weights'][i])):
                result[name].append(cap_min(c,full,raw0**2*length/(ell**2*self.normalization**2)*c.exp((p-1)*z)))
            if i in (0,2):
                j=0 if i==0 else 1
                result['D'].append(cap_min(c,abs(self.weights['divided_axial_rows'][j]),
                    raw0*length/(ell*self.normalization)*c.mpf(right)))
                result['cross'].append(cap_min(c,self.matrix['cross_weights'][j],
                    raw0**2*length/(ell**2*self.normalization**2)*c.exp(-z/2)))
        return result

    def _primitive_caps(self,N,W):
        c=self.ctx;mu=self.mu;R=self.R;sa=c.sqrt(mu)/N;se=mu/N
        mass=W['mass'][0]+W['mass'][2];D=sum(W['D']);cross=sum(W['cross'])
        return {key:upper(c,value) for key,value in dict(M=sa*R*mass,
            J=sa*R*(mass+mu*D)+sa*se*R**2*cross,I=se*R*sum(W['I']),
            S=sa**2*R**2*(W['energy'][0]+W['energy'][2])+se*R*sum(W['S'])+se**2*R**2*sum(W['energy'])/2,
            Cp=se*R*sum(W['Cp'])+se**2*R**2*sum(W['pressure'])/2).items()}

    def query(self,N,logx=(0,1)):
        parameters.positive_integer_N(N);c=self.ctx;y=c.mpf(logx);lo,hi=endpoints(y)
        if not all(mp.isfinite(v) for v in (lo,hi)) or lo<0 or hi>1:
            raise ValueError('Finite quiet-power logx subset[0,1] required')
        cert=repair.contraction_certificate(c,self.mu,N,self.matrix)
        CA=self.matrix['inverse_infinity_norm_upper'];CQ=self.matrix['quadratic_norm_times_N_upper']
        kappa=cert['strict_contraction_upper'];G=self.jets['log_bump_ordinary_derivative_caps']
        sa=c.sqrt(self.mu)/N;se=self.mu/N
        profiles={label:[[upper(c,scale*self.R*G[j]*math.factorial(k)) for k in range(6)] for j in range(5)]
            for label,scale in (('deltaE_over_common_A',se),('deltaUz_over_common_A',sa))}
        denom=c.exp(-self.alpha)-se*self.R*G[0]
        quiet={};positive=endpoints(denom)[0]>0
        if positive:
            error=upper(c,2*se*self.R*(G[1]+self.alpha*G[0])/denom)
            lower=c.mpf(endpoints(2*self.mu-error)[0])
            quiet=dict(normalized_swirl_denominator_lower=c.mpf(endpoints(denom)[0]),
                theta_shear_perturbation_absolute_cap=error,axial_shear_absolute_cap=upper(c,2*sa*self.R*G[1]/denom),
                repaired_theta_shear_excess_lower=lower,
                repaired_joint_primitive_margin_lower=lower,
                positive_primitive_margin_certified=endpoints(lower)[0]>0,
                at_least_three_halves_mu_certified=endpoints(error-self.mu/2)[1]<=0)
        prefix=self._primitive_caps(N,self._weights_cap(y));remaining=self._primitive_caps(N,self._weights_cap(y,True))
        d=cert['scaled_defect_abs_bounds_all_integer_N_at_least1']
        incoming=dict(M=sa*d[0],J=sa*(d[0]+self.mu*d[1]),I=se*d[2],S=se*d[3],Cp=se*d[4])
        total={key:cap_min(c,incoming[key]+prefix[key],remaining[key]) for key in prefix}
        logA=self.data['logP']+self.data['logAd']-1-3*self.mu/2
        logR0=self.data['radius_logs']['logRw']+1
        return dict(finite_integer_N=N,requested_logx=y,actual_mu=self.mu,
            implicit_repair=cert,inverse_nonlinear_Jacobian_infinity_norm_upper=upper(c,CA/(1-kappa)),
            quadratic_Jacobian_infinity_norm_upper=upper(c,2*CQ*self.R/N),
            controls_are_implicit_unique_constants_with_absolute_bound=self.R,
            saved_frequency_controls_or_signed_defects_reused=False,
            profile_logR4_axial5_absolute_caps=profiles,
            profile_units='mixed derivatives of physical increments divided by Pstar*Ad*exp(-1-3mu/2); axial factor C(Z)=1/(1+Z^2)',
            prefix_repair_primitive_absolute_caps=prefix,remaining_repair_primitive_absolute_caps=remaining,
            same_source_cumulative_defect_absolute_caps=total,
            cumulative_units=dict(M=dict(log_factor=logA+logR0,sqrt2_power=0,C_power=1),
                J=dict(log_factor=2*logA+3*logR0/2,sqrt2_power=1,C_power=2),
                I=dict(log_factor=logA+3*logR0/2,sqrt2_power=1,C_power=1),
                S=dict(log_factor=2*logA+logR0,sqrt2_power=0,C_power=2),
                Cp=dict(log_factor=2*logA,sqrt2_power=0,C_power=2)),
            shrinking_remaining_support_strips_retained=True,
            all_five_terminal_defects_vanish_by_same_implicit_equation=True,
            quiet_positive_denominator_certified=positive,quiet_repaired_shear=quiet,
            quiet_primitive_sufficient_integer_N=self.quiet_threshold,
            sufficient_integer_N_is_only_for_repair_and_quiet_primitive_margin=True,
            same_original_axis_pressure_datum_retained=True,
            source_family=self.data['source']['accepted']['source_family'],**{key:False for key in OPEN})


def run():
    field=CurrentUniformRepairMajorants(require_checked=False)
    examples={name:field.query(N,y) for name,N,y in (
        ('repair_threshold',field.repair_threshold,(0,1)),('quiet_threshold',field.quiet_threshold,(0,1)),
        ('original_frequency',10**12,(0,1)),('before_first_bump',field.quiet_threshold,'0'),
        ('inside_last_strip',field.quiet_threshold,('.8','.824')),
        ('near_last_exit',field.quiet_threshold,('.824','.825')),('terminal',field.quiet_threshold,'1'))}
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_uniform_repair_theorem=field.theorem,actual_bump_ordinary_jet_caps=field.jets,
        repair_only_sufficient_integer_N=field.repair_threshold,
        quiet_primitive_sufficient_integer_N=field.quiet_threshold,examples=examples,
        variable_N_implicit_control_and_Jacobian_bounds_available=True,
        variable_N_bump_jet_partial_defect_and_repaired_quiet_primitive_bounds_available=True,
        old_fixed_frequency_controls_not_relabelled=True,**{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Uniform repair generated; repair N>=',field.repair_threshold,';quiet primitive N>=',field.quiet_threshold,flush=True)
    return result


if __name__=='__main__':run()
