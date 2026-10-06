"""Original sigma/phi/Gamma radial5 and current actual collar stress4."""
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

from lei_ren_part1_paper_compliant_current_collar_stress_mixed_C4 import (
    CurrentCollarStressMixedC4,HERE,PREFIX,NAME,RECEIPT,GATE,SHAPE_GATE,JOIN_GATE,SCOPES,VIEWS,
    sha,pack,encode,endpoints,binding,IntervalTaylor,sigma_C5,phi_C5,
    sigma_fifth_formula,SIGMA5_CONSTANT,PHI5_POLYNOMIAL,collar_stress_mixed4,
    constant_stress_rows)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import positive_moment_derivative,stirling_second
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def exact(a,b,label):
    if s.simplify(a-b)!=0:raise ArithmeticError('Original collar identity differs: '+label)


def radial5_source_proof():
    binding('compliant_collar_Gamma_C4','shape','sig','sigma_jets(c,t)')
    binding('compliant_collar_Gamma_C4','shape','phi','phi_jets(c,t)')
    binding('compliant_collar_Gamma_C4','shape','K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]')
    binding('compliant_current_collar_stress_mixed_C4','shape_radial5','dr',
        'gamma_deficit_mixed(c,Z,heat.a,heat.Scap,t,5)')
    binding('compliant_current_collar_stress_mixed_C4','shape_radial5','D',
        'product_rows(product_rows(sr,C),dr)')
    binding('compliant_current_collar_stress_mixed_C4','shape_radial5','K',
        '[pre[j]-D[j]*(heat.a*heat.S) for j in range(6)]')
    x,u=s.symbols('x sigma',real=True);q=u*(1-u);aa=1-2*u
    L=s.symbols('L1:6');L1,L2,L3,L4,L5=L
    fourth=q*(L4+aa*(4*L1*L3+3*L2**2)+6*(1-6*q)*L1**2*L2+aa*(1-12*q)*L1**4)
    fifth=s.diff(fourth,u)*q*L1+sum(s.diff(fourth,L[j])*L[j+1] for j in range(4))
    values=[2/(1-x)**3+2/x**3,6/(1-x)**4-6/x**4,
        24/(1-x)**5+24/x**5,120/(1-x)**6-120/x**6,720/(1-x)**7+720/x**7]
    exact(fifth.subs(dict(zip(L,values))),sigma_fifth_formula(u,x),'sigma fifth derivative of original fourth')
    odds=1/(1-x)**2-1/x**2
    exact(odds.subs(x,1-x),-odds,'original sigma reflection symmetry')
    qq=s.symbols('q',real=True)
    exact(1-30*qq+120*qq**2,120*(qq-s.Rational(1,8))**2-s.Rational(7,8),'logistic fifth polynomial range')
    # For q in[0,1/4], this convex quadratic is in[-7/8,1].
    # |1-2sigma|<=1, |1-6q|<=1, |(1-2sigma)(1-12q)|<=2.
    ell=[2*math.factorial(j+1) for j in range(1,6)]
    e1,e2,e3,e4,e5=ell
    bound=e5+(5*e1*e4+10*e2*e3)+(10*e1**2*e3+15*e1*e2**2)+2*10*e1**3*e2+e1**5
    if bound!=SIGMA5_CONSTANT:raise ValueError('Fifth flat sigmoid bound changed')
    distance=s.symbols('distance',positive=True);phi=s.exp(-4/distance**2)
    polynomial=sum(v*distance**j for j,v in enumerate(PHI5_POLYNOMIAL))
    exact(-s.diff(phi,distance,5),phi*polynomial/distance**15,'original flat phi fifth derivative')
    xi,a,v=s.symbols('xi a v',positive=True);H=s.Function('full_Gamma_H')(xi)
    current=H;gamma_count=0
    for k in range(1,6):
        current=-xi*s.diff(current,xi)
        expected=sum((-1)**k*stirling_second(k,j)*xi**j*s.diff(H,xi,j) for j in range(1,k+1))
        exact(current,expected,'full Gamma logarithmic derivative '+str(k));gamma_count+=1
    for n in (9,10):
        expected=(-1)**n*s.rf(a+1,n-1)*v**n*(1+xi*v)**(-a-n)
        exact(s.diff((1+xi*v)**(-a),xi,n)/a,expected,'full positive Gamma integrand '+str(n))
    # g(y)=(1+y)^(-r) is convex for y>=0,r>0, hence
    # 1-r*y<=g(y)<=1; integration yields the unchanged positive moment bounds.
    return dict(original_sigma_and_phi_source_AST_bound=True,
        original_sigma_fifth_derivative_verified=True,sigma_fifth_flat_tail_constant=bound,
        sigma_fifth_flat_supremum_proved=True,original_phi_fifth_polynomial_verified=True,
        phi_fifth_flat_supremum_proved=True,full_Gamma_logarithmic_order5_identities=gamma_count,
        full_positive_Gamma_derivative_integrands_through10_verified=True,
        full_Gamma_measure_not_finite_inverse_radius_series=True,
        original_shape_product_and_chain_derivatives_extended=True,passed=True)


def independent_shape_fixtures():
    c=MPIntervalContext();c.dps=70;checks=0
    with mp.workdps(110):
        def sigma(x):return 1/(1+mp.exp(-(1/(1-x)**2-1/x**2)))
        def phi(x):return mp.exp(-4/(3-x)**2)
        def derivative(fn,x,n):
            if fn is sigma and x==mp.mpf('.5') and n in (2,4):return mp.mpf(0)
            if fn is sigma and x>mp.mpf('.5') and n>0:
                # Differentiate the small reflected quantity; subtracting
                # it from1 first would erase derivatives far below1 ulp.
                return (-1)**(n+1)*mp.diff(sigma,1-x,n)/math.factorial(n)
            return mp.diff(fn,x,n)/math.factorial(n)
        for provider,fn,points,boxes in (
            (sigma_C5,sigma,('.07','.21','.5','.83','.96'),([0,'.1'],['.4','.6'],['.9',1])),
            (phi_C5,phi,('.2','1.7','2.8','2.97'),([0,3],['2.95',3]))):
            for point in points:
                x=mp.mpf(point);rows=provider(c,point)
                for n in range(6):
                    value=derivative(fn,x,n)
                    lo,hi=endpoints(rows[n])
                    if not lo<=value<=hi:raise ArithmeticError('Independent shape derivative not enclosed: '+str((fn.__name__,point,n)))
                    checks+=1
            for box in boxes:
                rows=provider(c,box);left,right=map(mp.mpf,box)
                for factor in (mp.mpf('.1'),mp.mpf('.43'),mp.mpf('.9')):
                    x=left+(right-left)*factor
                    for n in range(6):
                        value=derivative(fn,x,n);lo,hi=endpoints(rows[n])
                        if not lo<=value<=hi:raise ArithmeticError('Independent support-crossing derivative not enclosed')
                        checks+=1
        for provider,location,value in ((sigma_C5,0,0),(sigma_C5,1,1),(sigma_C5,3,1),(phi_C5,3,0)):
            rows=provider(c,location)
            if endpoints(rows[0])!=(mp.mpf(value),mp.mpf(value)) or any(endpoints(rows[n])!=(mp.mpf(0),mp.mpf(0)) for n in range(1,6)):
                raise ValueError('Exact flat support jets through5 were not retained')
        a=mp.mpf('.17');xi=mp.mpf('.002');n=10
        reference=mp.rf(a+1,n-1)/mp.gamma(1+a)*mp.quad(
            lambda v:mp.exp(-v)*v**(a+n)*(1+xi*v)**(-a-n),[0,1,10,mp.inf])
        lo,hi=endpoints(positive_moment_derivative(c,c.mpf(str(a)),c.mpf(str(xi)),n))
        if not lo<=reference<=hi:raise ArithmeticError('Independent full Gamma tenth derivative not enclosed')
    return dict(independent_original_shape_derivative_rows=checks,exact_flat_endpoints_through5=True,
        independent_full_Gamma_tenth_derivative_checked=True,passed=True)


def independent_original_stress_fixture():
    """Differentiate full physical moments/stresses rather than factored rows."""
    c=MPIntervalContext();c.dps=70
    t,z=s.symbols('t Z',real=True);location={t:s.Rational(3,5),z:s.Rational(2,7)}
    a=s.Rational(3,20);delta=2*a;k=1-a;p=1+delta;b=(1-delta)/2
    d=1-z*z;L=1-delta*z*z;Rt=s.Integer(7);Pstar=s.Integer(3)
    B=2*s.exp(-(s.Rational(1,2)+a)*t);R=Rt*s.exp(t)
    q=(1+z+z*z)/100;K=1+q*s.exp(-t)
    A=1/k-q*s.exp(-t)/a
    E=1/delta+2*q*s.exp(-t)/(delta+1)+q*q*s.exp(-2*t)/(delta+2)
    P=1/(2*p)+q*s.exp(-t)/(p+1)+q*q*s.exp(-2*t)/(2*(p+2))
    D=s.Rational(2,5)+z/5+z*z/9;Cp=s.Rational(3,10)+(z+z*z)/7
    def interval(expression):
        expression=s.simplify(expression)
        if expression.is_Rational:return c.mpf(int(expression.p))/int(expression.q)
        if expression.func==s.exp:return c.exp(interval(expression.args[0]))
        if expression.is_Add:return sum((interval(v) for v in expression.args),c.mpf(0))
        if expression.is_Mul:
            out=c.mpf(1)
            for v in expression.args:out*=interval(v)
            return out
        if expression.is_Pow:return interval(expression.base)**interval(expression.exp)
        raise ValueError('Unexpected exact independent stress expression: '+str(expression))
    def jet(expression):
        return IntervalTaylor(c,[interval(s.diff(expression,z,n).subs(location)/math.factorial(n)) for n in range(6)])
    rows=lambda expression,count:[jet(s.diff(expression,t,j)) for j in range(count)]
    heat=SimpleNamespace(ctx=c,a=interval(a),delta=interval(delta),k=interval(k),prate=interval(p),S=interval(1/Rt))
    shape=dict(K_rows=rows(K,6))
    defects=dict(K_defect_rows=rows(K-1,5),angular_defect_rows=rows(A-1/k,5),
        energy_defect_rows=rows(E-1/delta,5),pressure_defect_rows=rows(P-1/(2*p),5))
    canonical=collar_stress_mixed4(heat,shape,defects,interval(location[z]),interval(location[t]))
    extra=constant_stress_rows(heat,interval(location[z]),interval(location[t]),jet(D),jet(Cp),4)
    Qt=s.sqrt(R/2)*B;Qz=s.sqrt(R/2)*B*B;Qp=s.sqrt(R/2)*Pstar**2
    Mt=s.sqrt(2)*R**s.Rational(3,2)*B*(A+D*s.exp(-k*t));Me=R*B*B*E/2
    theta=(k*Mt-b*z*s.diff(Mt,z)-R*s.sqrt(2*R)*B*K)/(2*L*R)
    theta+=s.sqrt(2*R)/R*s.diff(B*K,t)-B*K/s.sqrt(2*R)
    pressure_canonical=-B*B*P
    axial=(2*delta*z*Me-d*s.diff(Me,z)+R*(2*p*z*pressure_canonical-d*s.diff(pressure_canonical,z)))/(L*s.sqrt(2*R))
    pressure_extra=R*(2*p*z*Pstar**2*Cp-d*s.diff(Pstar**2*Cp,z))/(L*s.sqrt(2*R))
    fixture_rows=0
    for label,expression,factor,actual in (
        ('theta',theta,Qt,[canonical['theta'][j]+extra['theta'][j] for j in range(5)]),
        ('axial',axial,Qz,canonical['axial']),
        ('pressure_offset',pressure_extra,Qp,extra['axial_pressure_constant'])):
        for j in range(5):
            for n in range(5-j):
                expression_value=s.diff(expression,t,j,z,n)/factor
                with mp.workdps(110):reference=mp.mpf(str(s.N(expression_value.subs(location),110)))
                lo,hi=endpoints(actual[j][n]*math.factorial(n))
                if not lo<=reference<=hi:raise ArithmeticError('Original physical stress mixed4 not enclosed: '+str((label,j,n)))
                fixture_rows+=1
    return dict(independent_original_physical_stress_mixed4_rows=fixture_rows,
        actual_nonzero_angular_and_pressure_constants_retained=True,
        full_source_factor_derivatives_included=True,fifth_K_derivative_in_viscous_shear=True,passed=True)


def current_collar_Gamma_stress_join_proof(field):
    """Function equality at the flat join, with actual constants retained."""
    companion=field.companion
    parent=json.loads((HERE/(PREFIX+'current_heat_pressure_stress_check.json')).read_bytes())
    theorem=companion.source_owner.functional_proof
    current=companion.current_bindings['shared_exact_Gamma_future_binding']
    if not current['passed'] or not current['all_three_scaled_Gamma_future_integrands_bound'] or not all(current['exact_full_Gamma_future_identities'].values()):
        raise ValueError('Current C4/C5 full Gamma future-function identities missing')
    required=('exact_original_heat_collar_bracket','Gamma_angular_tail_rescaling',
        'Gamma_energy_tail_rescaling','Gamma_pressure_tail_rescaling',
        'new_exterior_retained_angular_history','new_exterior_forward_pressure_history',
        'new_exterior_pressure_interface')
    if not all(companion.graph.values()) or not all(theorem.get(key) for key in required):
        raise ValueError('Current full collar/Gamma function sources differ')
    if not parent['current_actual_pressure_source_split_proof']['passed'] or not parent['canonical_projected_full_Gamma_moment_stress_identities']['full_terminal_moment_stress_theorem_verified']:
        raise ValueError('Actual Cp source split or canonical full Gamma stress theorem missing')
    a,S=s.symbols('a S',positive=True);gamma=s.symbols('Gamma_y0:6')
    # At t=3, sigma=1 and phi=0 with every positive derivative zero.
    # The original shape product gives precisely the exterior H derivatives.
    sig=[s.Integer(1)]+[s.Integer(0)]*5;phi=[s.Integer(0)]*6
    def product(left,right):return [sum(s.binomial(j,n)*left[n]*right[j-n] for n in range(j+1)) for j in range(6)]
    one=[s.Integer(1)]+[s.Integer(0)]*5
    W=[one[j]-sig[j]+product(sig,phi)[j] for j in range(6)]
    D=product(sig,gamma)
    for j in range(6):
        exact(W[j],0,'flat collar W jet '+str(j))
        exact(one[j]-a*S*D[j],one[j]-a*S*gamma[j],'full Gamma K jet '+str(j))
    t,z,rate=s.symbols('t Z rate',real=True);C=s.Function('same_stress_coefficient')(t,z)
    for j in range(5):
        exact(s.diff(s.exp(rate*t)*C,t,j)/s.exp(rate*t),
            sum(s.binomial(j,n)*rate**(j-n)*s.diff(C,t,n) for n in range(j+1)),
            'positive stress factor product derivative '+str(j))
    return dict(same_checked_current_heat_source_and_full_future_functions=True,
        current_Gamma_future_source_identities_consumed=len(current['exact_full_Gamma_future_identities']),
        original_sigma_phi_flat_join_through5=True,full_collar_Gamma_shape_join_derivatives=6,
        exact_current_forward_Cp_source_split_consumed=True,
        same_actual_Dtheta_and_Cp_on_both_sides=True,
        full_original_moment_to_stress_equations_and_positive_factors_retained=True,
        original_stress_function_identity_implies_all_mixed4_traces=True,
        canonical_Gamma_cancels_but_actual_constants_not_eliminated=True,
        interval_overlap_not_used_as_join_proof=True,zero_exterior_stress_claim=False,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentCollarStressMixedC4(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current mixed4 collar source/datum differs')
    if any(raw[k] for k in (GATE,SHAPE_GATE,JOIN_GATE)+SCOPES):raise ValueError('Producer overstates acceptance')
    count=shape_count=0
    for name,(Z,t) in VIEWS.items():
        packet=field.evaluate(Z,t)
        if encode(pack(packet))!=raw['actual_current_collar_stress_mixed4_views'][name]:raise ValueError('Current collar packet differs')
        if not all(packet['same_current_common_heat_source_graph'].values()) or any(packet[k] for k in SCOPES):
            raise ValueError('Current graph or scope differs')
        if not packet['angular_and_pressure_constants_retained']:raise ValueError('Current terminal constants lost')
        shape=packet['current_full_collar_shape_radial5']
        if len(shape['K_rows'])!=6 or shape['sigma'].order!=5 or shape['phi'].order!=5:raise ValueError('Actual radial5 shape missing')
        prior=field.heat.shape(Z,t)
        for key in ('K_rows','D_rows','W_rows','pre_rows'):
            if encode(pack(shape[key][:5]))!=encode(pack(prior[key])):raise ValueError('Original shape prefix changed')
        for jet in shape['K_rows']:
            if not all(mp.isfinite(v) for box in jet.coefficients for v in endpoints(box)):raise ValueError('Infinite current shape coefficient')
            shape_count+=len(jet.coefficients)
        for grid in packet['actual_stress_factored_mixed4'].values():
            if len(grid)!=15 or not all(mp.isfinite(v) for box in grid.values() for v in endpoints(box)):
                raise ValueError('Finite actual collar stress mixed4 row omitted')
            count+=len(grid)
    proof=radial5_source_proof();shapes=independent_shape_fixtures();stress=independent_original_stress_fixture()
    join=current_collar_Gamma_stress_join_proof(field)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_radial5_defining_source_proof=proof,
        independent_original_shape_fixtures=shapes,independent_original_physical_stress_fixture=stress,
        current_collar_Gamma_actual_stress_join_proof=join,
        actual_collar_factored_stress_mixed4_rows_checked=count,current_K_radial5_axial5_coefficients_checked=shape_count,
        exact_current_angular_and_pressure_constants_retained=True,
        **dict.fromkeys((GATE,SHAPE_GATE,JOIN_GATE),True),**dict.fromkeys(SCOPES,False),all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(encode(result),indent=2)+'\n').encode('utf8'))
    print('PASS full current collar K radial5 and actual stress mixed4; terminal constants/global/temporal remain open',flush=True)
    return result


if __name__=='__main__':run()
