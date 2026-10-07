"""Own-moment radial and cumulative absolute-pressure error mixed jets.

Bounds the changed O2/O3 modulation and quiet repair source functions for
variable N, with original datum and Pstar sectors. No ancestor constructors
or fixed-frequency signed controls are used. Completed tensor cones remain
open. Returned native radial rows include their sqrt(R/2) derivative shift.
"""
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_uniform_repair_majorants as uniform
import lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator as histories
from lei_ren_part1_paper_compliant_current_O3_finite_frequency_profiles import cutoff_rows

parameters=uniform.parameters
HERE,PREFIX,sha=uniform.HERE,uniform.PREFIX,uniform.sha
NAME=PREFIX+'current_O3_recovered_error_majorants.json'
RECEIPT=PREFIX+'current_O3_recovered_error_majorants_check.json'
endpoints,upper=uniform.endpoints,uniform.upper
OPEN=('completed_signed_tensor_error_bounds_available','common_N_modified_cones_certified',
      'global_physical_velocity_interface_composition_certified','actual_coefficient_recursion_certified',
      'finite_energy_certified','full_corrected_NS_certified')


def products(a,b,order=4):
    return [sum(math.comb(j,i)*a[i]*b[j-i] for i in range(j+1)) for j in range(order+1)]


def scalar_history_caps(initial,old,increment,axial):
    """Absolute ordinary logR bounds from the defining history ODEs."""
    d={label:[initial[label]] for label in ('m','h','k','e','p')}
    modified=[old[j]+increment[j] for j in range(5)]
    ksource=products(modified,axial);kinetic=products(axial,axial)
    cross=products(old,increment);square=products(increment,increment)
    for j in range(4):
        d['m'].append(axial[j]+d['m'][j])
        d['h'].append(increment[j]+3*d['h'][j]/2)
        d['k'].append(ksource[j]+3*d['k'][j]/2)
        d['e'].append(kinetic[j]+cross[j]+square[j]/2+d['e'][j])
        d['p'].append(cross[j]+square[j]/2)
    return d


def radial_factor_caps(c,delta):
    """Full axial derivative caps for the actual rational recovery factors."""
    lo,hi=endpoints(delta)
    if lo<0 or hi>=1:raise ValueError('Actual 0<=delta<1 required for radial recovery')
    H=[upper(c,1/(1-delta))]
    for k in range(1,6):
        term=2*delta*k*H[k-1]
        if k>1:term+=delta*k*(k-1)*H[k-2]
        H.append(upper(c,term/(1-delta)))
    Z=[c.mpf(1),c.mpf(1)]+[c.mpf(0)]*4
    one_minus_Z2=[c.mpf(1),c.mpf(2),c.mpf(2)]+[c.mpf(0)]*3
    C=[c.mpf(math.factorial(k)) for k in range(7)]
    zc=products(Z,C,5);dc=products(one_minus_Z2,C[1:],5)
    V=[upper(c,2*v) for v in products(zc,H,5)]
    M=[upper(c,v) for v in products([abs(1-delta)*zc[k]+dc[k] for k in range(6)],H,5)]
    return dict(inverse_L_axial_derivative_caps=H,axial_velocity_factor_caps=V,
                own_M_factor_caps=M,whole_axial_domain=(-1,1),minimum_L=1-delta)


def recover(c,d,axial,factors):
    # m and axial are Ad*C times scalar rows. Pressure is Ad²*C².
    Q=[[factors['axial_velocity_factor_caps'][k]*axial[j]+factors['own_M_factor_caps'][k]*d['m'][j]
        for k in range(6)] for j in range(5)]
    radial=[[upper(c,sum((math.comb(j,i)*c.mpf('.5')**(j-i)*Q[i][k] for i in range(j+1)),c.mpf(0)))
        for k in range(6)] for j in range(5)]
    pressure=[[upper(c,d['p'][j]*math.factorial(k+1)) for k in range(6)] for j in range(5)]
    all_histories={label:[[upper(c,value*math.factorial(k+(label in ('k','e','p')))) for k in range(6)]
        for value in rows] for label,rows in d.items()}
    return dict(radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5=radial,
        absolute_pressure_error_over_Pstar2_Ad2_ordinary_logR4_axial5=pressure,
        five_normalized_history_error_ordinary_logR4_axial5=all_histories)


class _Formula:
    """Small symbolic field for replay of the original pure recovery program."""
    def __init__(self,value):self.value=s.sympify(value)
    @staticmethod
    def unwrap(value):return value.value if isinstance(value,_Formula) else s.sympify(value)
    def __add__(self,other):return _Formula(self.value+self.unwrap(other))
    __radd__=__add__
    def __sub__(self,other):return _Formula(self.value-self.unwrap(other))
    def __rsub__(self,other):return _Formula(self.unwrap(other)-self.value)
    def __mul__(self,other):return _Formula(self.value*self.unwrap(other))
    __rmul__=__mul__
    def __truediv__(self,other):return _Formula(self.value/self.unwrap(other))
    def __pow__(self,other):return _Formula(self.value**other)
    def reciprocal(self):return _Formula(1/self.value)


def exact_theorem():
    asts=histories.SourceAST();checks={}
    def zero(name,a,b):
        difference=s.cancel(s.together(_Formula.unwrap(a)-_Formula.unwrap(b)))
        if difference!=0:raise ArithmeticError('Recovered error source identity: '+name)
        checks[name]=True
    Z,delta,Ad=s.symbols('Z delta Ad',real=True)
    C=1/(1+Z**2);L=1-delta*Z**2
    old=s.symbols('old0:5');inc=s.symbols('increment0:5');axial=s.symbols('axial0:5')
    initial=dict(zip(('m','h','k','e','p'),s.symbols('m0 h0 k0 e0 p0')))
    program=asts.replay('current_O3_modulated_histories_operator','increment_rows',dict(
        square=lambda x:x*x,derivative=lambda x:_Formula(s.diff(x.value,Z)),
        product_rows=lambda a,b:products(a,b),SHAPES=histories.SHAPES))
    d,Q=program(SimpleNamespace(mpf=lambda v:s.Rational(str(v))),_Formula(Z),delta,Ad,initial,
        [_Formula(Ad*C*(old[j]+inc[j])) for j in range(5)],[_Formula(Ad*C*v) for v in old],
        [_Formula(Ad*C*v) for v in axial],[_Formula(Ad*C*v) for v in inc])
    signed={key:[value] for key,value in initial.items()}
    ksource=products([old[j]+inc[j] for j in range(5)],axial)
    kinetic=products(axial,axial);cross=products(old,inc);square=products(inc,inc)
    for j in range(4):
        signed['m'].append(axial[j]-signed['m'][j])
        signed['h'].append(inc[j]-s.Rational(3,2)*signed['h'][j])
        signed['k'].append(ksource[j]-s.Rational(3,2)*signed['k'][j])
        signed['e'].append(kinetic[j]-cross[j]-square[j]/2-signed['e'][j])
        signed['p'].append(cross[j]+square[j]/2)
    for key in d:
        power=1 if key in ('m','h') else 2
        for j in range(5):zero('same_actual_history_source_%s_%d'%(key,j),d[key][j],Ad**power*C**power*signed[key][j])
    rV=2*Z*C/L;rM=(-(1-delta)*Z*C-(1-Z**2)*s.diff(C,Z))/L
    for j in range(5):zero('same_actual_radial_own_M_factorization_'+str(j),Q[j],
        Ad*(rV*axial[j]+rM*signed['m'][j]))
    zero('pressure_value_is_cumulative_history',d['p'][0],Ad**2*C**2*initial['p'])
    zero('pressure_first_logR_derivative_has_cross_term',d['p'][1],Ad**2*C**2*(old[0]*inc[0]+inc[0]**2/2))
    shifted=asts.replay('collar_stress_C3','shifted_rows',dict(math=math))(Q,s.Rational(1,2),4)
    for j in range(5):zero('same_single_sqrtR_radial_derivative_shift_'+str(j),shifted[j],
        sum(math.comb(j,i)*s.Rational(1,2)**(j-i)*Q[i].value for i in range(j+1)))
    for k in range(7):zero('C_global_partial_fraction_derivative_'+str(k),s.diff(C,Z,k),
        (-1)**k*math.factorial(k)*((Z-s.I)**(-k-1)-(Z+s.I)**(-k-1))/(2*s.I))
    for k in range(6):zero('C_squared_factorial_derivative_cap_'+str(k),
        sum(math.comb(k,i)*math.factorial(i)*math.factorial(k-i) for i in range(k+1)),math.factorial(k+1))
    # This exact differentiated identity proves the inverse-L cap recurrence.
    H=s.Function('same_inverse_L')(Z)
    for k in range(1,6):zero('inverse_L_Leibniz_recurrence_'+str(k),s.diff(L*H,Z,k),
        L*s.diff(H,Z,k)-2*delta*k*Z*s.diff(H,Z,k-1)-delta*k*(k-1)*s.diff(H,Z,max(0,k-2)))
    y,f2=s.symbols('quiet_logx same_f2',real=True)
    quiet_y=asts.evaluate(asts.expression('current_O3_repaired_histories','history','y',wanted='q-1'),dict(q=1+y))
    quiet_t=asts.evaluate(asts.expression('current_O3_modulated_histories','history','t',wanted='1+v'),dict(v=1+quiet_y))
    zero('same_quiet_coordinate_after_modulation_support',quiet_t,2+y)
    for target,endpoint in (('upper_t','hi'),('lower_t','lo')):
        asts.expression('current_O3_modulated_histories_operator','cumulative_scalar_enclosures',target,
            wanted="min(%s,mp.mpf('.5'))"%endpoint)
    incoming=dict(zip(('M','I','J','S','Cp'),s.symbols('inM inI inJ inS inCp')))
    partial=dict(zip(('M','I','J','S','Cp'),s.symbols('partialM partialI partialJ partialS partialCp')))
    recipes=(('m','M',1,1),('h','I',1,s.Rational(3,2)),('k','J',2,s.Rational(3,2)),('e','S',2,1),('p','Cp',2,0))
    previous={key:f2**power*s.exp(-rate*y)*incoming[label] for key,label,power,rate in recipes}
    scalar=asts.evaluate(asts.expression('current_O3_repaired_histories','history','scalar',wanted=
        "dict(m=previous['m']+self.f2*primitives['M']/x,h=previous['h']+self.f2*primitives['I']/x**c.mpf('1.5'),k=previous['k']+self.f2**2*primitives['J']/x**c.mpf('1.5'),e=previous['e']+self.f2**2*primitives['S']/x,p=previous['p']+self.f2**2*primitives['Cp'])"),
        dict(c=SimpleNamespace(mpf=lambda v:s.Rational(str(v))),self=SimpleNamespace(f2=f2),previous=previous,primitives=partial,x=s.exp(y)))
    for key,label,power,rate in recipes:zero('same_quiet_fixed_units_to_scalar_history_'+key,scalar[key],
        f2**power*s.exp(-rate*y)*(incoming[label]+partial[label]))
    asts.expression('current_O3_modulated_histories','history','pressure',wanted='[a+b for a,b in zip(original_P,delta_rows[\'p\'])]')
    asts.expression('current_O3_repaired_histories','history','original_pressure',wanted="parent['original_absolute_pressure_over_Pstar2_ordinary_logR_rows']")
    asts.expression('current_O3_repaired_histories','history','Ahat',wanted="(1+square(z)).reciprocal()*(self.histories.Ua*self.f2)")
    asts.expression('current_O3_modulated_histories','history','radial_increment',wanted="shifted_rows(Q,c.mpf('.5'),4)")
    for fn in ('frequency_uniform_bounds','exact_modulated_history_theorem'):
        asts.method('current_O3_modulated_histories_operator',fn)
    return dict(passed=True,identities=checks,input_hashes=asts.hashes,
        pressure_identity='delta_p/Pstar^2 is the own cumulative Cp history; its first logR derivative is Uold*deltaUtheta+(deltaUtheta)^2/2',
        radial_identity='delta_ur=Pstar*Ad*sqrt(R/2)*(rV*axial_scalar+rM*own_M_scalar)',
        original_axis_pressure_datum_is_unchanged=True,ordinary_derivatives_not_Taylor_coefficients=True,
        quiet_support_and_homogeneous_transport_consumed_from_checked_uniform_repair=True,
        whole_Z_factor_argument='C partial fractions through6, C² Leibniz through5, inverse(1-delta Z²) recurrence through5',
        completed_tensor_cones_not_implied=True)


class CurrentRecoveredErrorMajorants:
    def __init__(self,require_checked=True):
        self.repair=uniform.CurrentUniformRepairMajorants();self.data=self.repair.data;self.ctx=self.repair.ctx
        self.mu=self.data['mu'];self.delta=self.data['source']['delta']
        self.factors=radial_factor_caps(self.ctx,self.delta)
        self.norms=parameters.normalized_source_norms(self.ctx,self.mu)
        self.frequency=histories.frequency_uniform_bounds(self.ctx,self.mu)
        self.theorem=exact_theorem()
        self.hashes={**self.repair.hashes,uniform.RECEIPT:sha(uniform.RECEIPT),**self.theorem['input_hashes'],Path(__file__).name:sha(Path(__file__).name)}
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed']:raise ValueError('Accepted recovered error jets required')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed recovered error dependency: '+name)

    def _finish(self,N,coordinate,d,axial,logR,details):
        c=self.ctx;bounds=recover(c,d,axial,self.factors)
        signed={key:[[c.mpf([-1,1])*upper(c,v) for v in row] for row in rows]
            for key,rows in bounds['five_normalized_history_error_ordinary_logR4_axial5'].items()}
        return dict(finite_integer_N=N,coordinate=coordinate,actual_mu=self.mu,actual_delta=self.delta,
            normalized_scalar_history_error_ordinary_logR_caps=d,
            normalized_signed_history_error_enclosures=signed,**bounds,
            radial_factor_axial_caps=self.factors,
            native_cylindrical_radial_error_log_prefactor=self.data['logP']+self.data['logAd']+logR/2-c.ln(2)/2,
            native_absolute_pressure_error_log_prefactor=2*(self.data['logP']+self.data['logAd']),
            native_logR=logR,radial_sqrtR_derivative_shift_already_applied=True,
            pressure_is_cumulative_value_not_local_cross_plus_square=True,
            same_original_axis_pressure_datum_retained=True,
            history_units=dict(m='delta(Mz/R)/(Pstar*Ad)',h='delta(Mtheta/(sqrt2*R^1.5*Pstar))/Ad',
                k='delta(Mtheta_z/(sqrt2*R^1.5*Pstar))/(Pstar*Ad^2)',
                e='delta(Mztheta/(R*Pstar^2))/Ad^2',p='delta(Mp/Pstar^2)/Ad^2'),
            actual_current_changed_O2_O3_native_radial_pressure_error_jets_available=True,
            source_family=self.data['source']['accepted']['source_family'],**details,**{key:False for key in OPEN})

    def modulation(self,offset,N):
        parameters.positive_integer_N(N);c=self.ctx;t=c.mpf(offset);lo,hi=endpoints(t)
        if not all(mp.isfinite(v) for v in (lo,hi)) or lo< -11 or hi>1:
            raise ValueError('Finite current buffer/transition logR offset subset[-11,1] required')
        K=self.frequency['untransported_abs_integral_times_N_upper']
        zero=hi<=-2
        initial={key:upper(c,K[key]*c.exp(-c.mpf(str(rate))*t)/N) if not zero else c.mpf(0)
            for key,rate in (('m',1),('h','1.5'),('k','1.5'),('e',1),('p',0))}
        kinetic=upper(c,K['e_kinetic']*c.exp(-t)/N**2) if not zero else c.mpf(0)
        swirl=upper(c,K['e_swirl']*c.exp(-t)/N) if not zero else c.mpf(0)
        initial['e']=uniform.cap_min(c,initial['e'],kinetic+swirl)
        old=[upper(c,parameters.maximum(c,self.norms['T'][j],c.exp(-t/2)/2**j)) for j in range(5)]
        C=[uniform.cap_min(c,self.norms['C'][j],abs(v)) for j,v in enumerate(cutoff_rows(c,t))]
        F=[[old[j]*math.factorial(k) for k in range(6)] for j in range(5)]
        profile=parameters.normalized_envelopes(c,self.mu,N,C,F,self.norms['L1'])
        inc=[upper(c,row[0]) for row in profile['theta_increment_over_Pstar_over_Ad_majorants']]
        axial=[upper(c,row[0]) for row in profile['modified_axial_over_Pstar_over_Ad_majorants']]
        d=scalar_history_caps(initial,old,inc,axial)
        lower=c.mpf(0)
        if lo>=0:lower=c.mpf(endpoints(self.frequency['positive_kinetic_buffer_mass_times_N_squared_lower']*c.exp(-t)/N**2)[0])
        return self._finish(N,t,d,axial,self.data['radius_logs']['logRd']+t,dict(
            chart='O2_O3_modulation',positive_kinetic_history_scalar_enclosure=c.mpf([endpoints(lower)[0],endpoints(kinetic)[1]]),
            absolute_swirl_history_scalar_cap=swirl,local_cutoff_derivative_caps=C,
            original_source_errors_vanish_before_flat_left_support=zero,
            independent_repair_not_yet_active_on_this_chart=True))

    def quiet_repair(self,N,logx=(0,1)):
        inherited=self.repair.query(N,logx);c=self.ctx;y=inherited['requested_logx'];mu=self.mu
        f2=c.exp(-1-3*mu/2);alpha=c.mpf('.5')+mu
        local=[]
        for center in self.repair.weights['centers']:
            local.append(uniform.histories.bump_rows(c,y,center,self.repair.weights['radius'],self.repair.normalization))
        G=[uniform.cap_min(c,self.repair.jets['log_bump_ordinary_derivative_caps'][j],
            c.mpf(max(endpoints(abs(rows[j]))[1] for rows in local))) for j in range(5)]
        inc=[upper(c,f2*mu*self.repair.R*g/N) for g in G]
        axial=[upper(c,f2*c.sqrt(mu)*self.repair.R*g/N) for g in G]
        old=[upper(c,f2*alpha**j*c.exp(-alpha*y)) for j in range(5)]
        total=inherited['same_source_cumulative_defect_absolute_caps']
        initial={key:upper(c,f2**power*c.exp(-c.mpf(str(rate))*y)*total[label])
            for key,label,power,rate in (('m','M',1,1),('h','I',1,'1.5'),('k','J',2,'1.5'),('e','S',2,1),('p','Cp',2,0))}
        d=scalar_history_caps(initial,old,inc,axial)
        zero=all(v._mpi_==c.mpf(0)._mpi_ for rows in d.values() for v in rows)
        return self._finish(N,y,d,axial,self.data['radius_logs']['logRw']+1+y,dict(
            chart='quiet_O3_independent_repair',local_disjoint_bump_logR_derivative_caps=G,
            same_source_shrinking_terminal_defect_caps=total,
            all_error_logR4_axial5_rows_vanish_after_last_bump=zero,
            finite_N_controls_remain_same_unique_implicit_vector=True,
            quiet_primitive_sufficient_integer_N=self.repair.quiet_threshold))


def run():
    field=CurrentRecoveredErrorMajorants(require_checked=False)
    examples={name:field.modulation(t,N) for name,t,N in (
        ('before_modulation','-3',1),('left_taper','-1.99',37),('old_seam','0',37),
        ('after_cutoff','.75',37),('crossing_seam',('-.1','.1'),10**12))}
    examples.update({name:field.quiet_repair(N,y) for name,N,y in (
        ('quiet_whole',field.repair.quiet_threshold,(0,1)),('first_bump',10**12,('.19','.21')),
        ('near_last_exit',field.repair.quiet_threshold,('.824','.825')),('terminal',field.repair.quiet_threshold,'1'))})
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_recovered_error_source_theorem=field.theorem,examples=examples,
        current_O2_O3_and_quiet_radial_pressure_mixed4_axial5_error_bounds_available=True,
        original_pressure_datum_and_single_radial_shift_retained=True,
        old_fixed_frequency_controls_not_reused=True,**{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Recovered own-M radial and cumulative absolute-pressure error jets generated; native logR4/Z5',flush=True)
    return result


if __name__=='__main__':run()
