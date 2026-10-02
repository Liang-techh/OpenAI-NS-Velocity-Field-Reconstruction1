"""Independent full-kernel quadrature and source moment normalization checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_long_reshape_profiles import (
    CompliantLongReshapeProfiles,IntervalTaylor,backward_kernel,
    normalized_amplitude_decay,relative_amplitude_jet)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_long_reshape_profiles.json'


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    a=mp.exp(-1/t**2);b=mp.exp(-1/(1-t)**2)
    return a/(a+b)


def finite_kernel_fixture():
    """Full independent quadrature, both signs of B and all axial orders5."""
    with mp.workdps(80):
        c=MPIntervalContext();c.dps=110;T=mp.mpf(4000);tol=mp.mpf('1e-65')
        Z,delta=s.symbols('Z delta',real=True)
        counts=dict(kernels=0,inherited=0,amplitudes=0)
        for sign in (-1,1):
            polynomial=sign*3+s.Rational(7,10)*Z-Z*Z/5+Z**5/200
            fn=s.lambdify(Z,polynomial,'mpmath')
            initial=lambda z:mp.mpf('.3')+z/20+z*z/100+z**4/200
            z=mp.mpf('.3')
            def jet(f):
                return IntervalTaylor(c,[c.mpf([q-tol,q+tol])/math.factorial(k)
                    for k in range(6) for q in (mp.diff(f,z,k),)])
            B=jet(fn);inlet=jet(initial)
            for y in (mp.mpf(0),mp.mpf(1),T/2,T):
                sig=sigma(y/T)
                # Direct symbolic Z differentiation of the original
                # integrand, not the bound recurrence used by the provider.
                for kind,kappa,m,low,high in (('theta','1.6',1,'1.55','1.65'),
                        ('pressure','.2',2,'.1','.3'),('swirl','1.2',2,'1.1','1.3')):
                    bounded=backward_kernel(c,B,c.mpf(T),c.mpf(y),m,low,high,[])
                    expression=s.exp(m*polynomial*delta)
                    derivatives=[s.lambdify((Z,delta),s.simplify(s.diff(expression,Z,n)/expression),'mpmath')
                                 for n in range(6)]
                    cuts=[mp.mpf(0)]+[mp.mpf(q) for q in (1,5,20,80,400) if q<y]+([y] if y else [])
                    for n in range(6):
                        def integrand(t):
                            change=sigma((y-t)/T)-sig
                            return mp.exp(-mp.mpf(kappa)*t+m*fn(z)*change)*derivatives[n](z,change)
                        expected=mp.quad(integrand,cuts) if y else mp.mpf(0)
                        a,b=endpoints(bounded[n]*math.factorial(n))
                        if not a-tol*100<=expected<=b+tol*100:
                            raise ArithmeticError('Independent full '+kind+' kernel derivative failed')
                        counts['kernels']+=1
                    logjet=B*(m*c.mpf(sig))-c.mpf(y)*c.mpf(kappa)
                    inherited=normalized_amplitude_decay(c,logjet,low,high,c.mpf(y),[],initial=inlet)
                    exact=lambda zz:initial(zz)*mp.exp(-mp.mpf(kappa)*y+m*fn(zz)*sig)
                    for n in range(6):
                        expected=mp.diff(exact,z,n);a,b=endpoints(inherited[n]*math.factorial(n))
                        if not a-tol*100<=expected<=b+tol*100:
                            raise ArithmeticError('Independent nonzero inlet decay derivative failed')
                        counts['inherited']+=1
                # Actual amplitude derivatives are distinct from shape
                # derivatives; verify normalized Utheta and Utheta^2 jets.
                logu=B*(1-c.mpf(sig))+c.mpf(y)/10
                for m in (1,2):
                    jetu=relative_amplitude_jet(logu,m)
                    actual=lambda zz:mp.exp(m*(y/10+(1-sig)*fn(zz)))
                    base=actual(z)
                    for n in range(6):
                        expected=mp.diff(actual,z,n)/base;a,b=endpoints(jetu[n]*math.factorial(n))
                        if not a-tol*100<=expected<=b+tol*100:
                            raise ArithmeticError('Current amplitude Bell derivative failed')
                        counts['amplitudes']+=1
        return dict(independent_full_kernel_axial_derivatives=counts['kernels'],
            independent_nonzero_history_decay_axial_derivatives=counts['inherited'],
            independent_current_amplitude_derivatives=counts['amplitudes'],
            both_signs_of_B_tested=True,finite_parameter_fixture_only=True,
            actual_source_admission=False,passed=True)


def structural_identities():
    R=s.symbols('R',positive=True);u=s.Function('u')(R);V=s.symbols('V',real=True);q=s.symbols('q',real=True)
    H,K,m,A,b,p=[s.Function(name)(R) for name in ('H','K','m','A','b','p')]
    # Differentiate the PHYSICAL primitives directly. Normalize only after
    # differentiation, so a misplaced root2 or R power cannot pass.
    moments=(s.sqrt(2)*R**s.Rational(3,2)*u*H,R*m,
             s.sqrt(2)*R**s.Rational(3,2)*u*K,R*A-R*u*u*b/2,u*u*p/2)
    updates={s.diff(u,R):q*u/R,
        s.diff(H,R):(1-(s.Rational(3,2)+q)*H)/R,
        s.diff(K,R):(V-(s.Rational(3,2)+q)*K)/R,
        s.diff(m,R):(V-m)/R,s.diff(A,R):(V*V-A)/R,
        s.diff(b,R):(1-(1+2*q)*b)/R,s.diff(p,R):(1-2*q*p)/R}
    rhs=(s.sqrt(2*R)*u,V,s.sqrt(2*R)*u*V,V*V-u*u/2,u*u/(2*R))
    for moment,expected in zip(moments,rhs):
        if s.simplify(s.diff(moment,R).xreplace(updates)-expected)!=0:
            raise ArithmeticError('Five physical primitive normalization identity failed')
    z,logC,logP,y,T,S=s.symbols('z logC logP y T S',real=True)
    B=s.symbols('B',real=True);Q=1+z*z
    logfield=y/10-logC-logP-s.log(Q)+(1-S)*B
    if s.simplify(logfield.subs(S,1)-(y/10-logC-logP-s.log(Q)))!=0:
        raise ArithmeticError('Reshape reference endpoint log identity failed')
    offset=y-10*(logC+logP)
    if s.simplify(logfield.subs(S,1)-(offset/10-s.log(Q)))!=0:
        raise ArithmeticError('Reference radius compatibility identity failed')
    return dict(exact_five_physical_primitive_RHS=5,exact_reference_log_endpoint=1,
                exact_reference_radius_compatibility=1,passed=True)


def run():
    with mp.workdps(280):
        receipt=json.loads((HERE/NAME).read_bytes())
        for name,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Reshape source changed: '+name)
        c=MPIntervalContext();c.dps=240;capcount=0;bounds=0;moments=0
        for proof in receipt['factored_exponential_cap_proofs']:
            logvalue=read_interval(c,proof['log_magnitude_upper']);cap=read_interval(c,proof['log_cap'])
            if endpoints(logvalue)[1]>endpoints(cap)[0] or not proof['exact_source_not_replaced']:
                raise ArithmeticError('Factored actual-source exponential cap was not proved')
            capcount+=1
        join=json.loads((HERE/'lei_ren_part1_paper_compliant_reference_join_bounds.json').read_bytes())
        if not (join['reference_join_family_sha256']==receipt['reference_join_family_sha256']
                and join['admitted_inner_parameter_family_sha256']==receipt['admitted_inner_parameter_family_sha256']
                and join['same_family_R110_Rh_relaxed_cone_analytically_certified']
                and join['actual_moments_continuously_inherited']):
            raise ValueError('Actual B C2/shear source theorem is not bound')
        A=read_interval(c,receipt['Abar']);T=read_interval(c,receipt['T'])
        selected=json.loads((HERE/'lei_ren_part1_paper_compliant_physical_norm_family.json').read_bytes())
        if endpoints(A)!=endpoints(read_interval(c,selected['A_upper'])) or endpoints(T)!=endpoints(400*A):
            raise ArithmeticError('Original selected A/T changed')
        if endpoints(read_interval(c,join['B_C2_upper']))!=endpoints(2*A):raise ArithmeticError('Actual B C2 theorem changed')
        for stored,expected in zip(join['shape_shear_interval'],('.7','.9')):
            lo,hi=endpoints(read_interval(c,stored));a,b=endpoints(c.mpf(expected))
            if not lo<=a<=b<=hi:raise ArithmeticError('Original shear certificate changed')
        step=json.loads((HERE/'lei_ren_part1_paper_shared_fixed_step_bound.json').read_bytes())
        if not step['global_derivative_bound_certified'] or step['global_derivative_upper']!=8:
            raise ValueError('Sigma derivative source bound missing')
        radius=read_interval(c,join['B_C2_upper'])*step['global_derivative_upper']/T
        q=c.mpf([endpoints(c.mpf('.1')-radius)[0],endpoints(c.mpf('.1')+radius)[1]])
        if endpoints(q)!=endpoints(read_interval(c,receipt['direct_same_source_log_u_slope_interval'])):
            raise ArithmeticError('Source log-u rate was not recomputed correctly')
        for name,shift,multiplier in (('theta','1.5',1),('pressure','0',2),('swirl','1',2)):
            rates=q*multiplier+c.mpf(shift);low,high=map(c.mpf,receipt['kernel_rate_bounds'][name])
            if endpoints(rates)[0]<endpoints(low)[1] or endpoints(rates)[1]>endpoints(high)[0]:
                raise ArithmeticError('Full kernel rate bounds do not follow from actual source '+name)
        packets=[receipt[n] for n in ('whole_reshape','actual_R110_inlet','actual_Rsh_exit')]+receipt['samples']
        for packet in packets:
            for name in ('log_Utheta_over_Pstar_axial5_coefficients','Utheta_true_axial5_divided_by_current_Utheta',
                         'Uz_actual_axial5_coefficients','actual_Q_axial4_coefficients',
                         'pressure_axis_axial5_coefficients','angular_shear_axial5_coefficients'):
                row=packet[name]
                if len(row)!=(5 if name=='actual_Q_axial4_coefficients' else 6):raise ValueError('Reshape axial order lost')
                for value in row:
                    lo,hi=endpoints(read_interval(c,value))
                    if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite actual reshape profile')
                    bounds+=1
            for group in ('actual_normalized_moment_shape_axial5_coefficients',
                          'actual_true_moment_derivatives_divided_by_current_prefactor',
                          'full_backward_kernel_axial5_coefficients',
                          'actual_inherited_R110_moment_contribution_axial5_coefficients'):
                for row in packet[group].values():
                    if len(row)!=6:raise ValueError('Moment/kernel axial5 order lost')
                    for value in row:
                        lo,hi=endpoints(read_interval(c,value))
                        if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite normalized moment/kernel')
                        moments+=1
            a,b=endpoints(read_interval(c,packet['angular_shear_axial5_coefficients'][0]))
            if a<endpoints(c.mpf('.7'))[0] or b>endpoints(c.mpf('.9'))[1]:raise ArithmeticError('Original shear range lost')
            if packet['moment_reset'] or packet['numerical_kernel_truncation'] or packet['reshape_radial_mixed4_certified']:
                raise ValueError('Source/history or scope changed')
        provider=CompliantLongReshapeProfiles();inputs=provider.inputs([-1,1])
        start=receipt['actual_R110_inlet'];end=receipt['actual_Rsh_exit']
        # Compare all exact y=0 normalized histories to the parent, rather
        # than accepting a broad interval overlap as an interface proof.
        for name,jet in inputs['moments'].items():
            stored=start['actual_normalized_moment_shape_axial5_coefficients'][name]
            if any(endpoints(read_interval(c,q))!=endpoints(value) for q,value in zip(stored,jet.coefficients)):
                raise ArithmeticError('Actual R110 normalized histories not preserved exactly')
        parent=json.loads((HERE/'lei_ren_part1_paper_compliant_inner_switch_profiles.json').read_bytes())['actual_R110_inlet']
        if start['Uz_actual_axial5_coefficients']!=parent['Uz_actual_axial5_coefficients'] or start['pressure_axis_axial5_coefficients']!=parent['pressure_axis_axial5_coefficients']:
            raise ArithmeticError('Actual axial velocity/pressure inlet changed')
        if any(p['Uz_actual_axial5_coefficients']!=start['Uz_actual_axial5_coefficients'] for p in (end,receipt['whole_reshape'])):
            raise ArithmeticError('Long reshape introduced axial velocity change')
        if any(endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0))
               for row in start['full_backward_kernel_axial5_coefficients'].values() for value in row):
            raise ArithmeticError('Reshape kernel did not vanish exactly at inlet')
        result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
            implicit_source_sha256=receipt['implicit_source_sha256'],reference_join_family_sha256=receipt['reference_join_family_sha256'],
            kernel_fixture=finite_kernel_fixture(),structural_identities=structural_identities(),
            actual_factored_source_exponential_caps_checked=capcount,actual_log_velocity_and_pressure_bounds_checked=bounds,
            actual_normalized_moment_and_kernel_bounds_checked=moments,
            original_selected_A_T_retained=True,exact_R110_histories_retained=True,
            actual_B_C2_shear_sigma_and_kernel_rate_theorems_directly_bound=True,
            actual_long_reshape_log_velocity_axial5_available=True,
            actual_long_reshape_moment_axial5_enclosures_available=True,
            actual_long_reshape_radial_recovery_axial4_available=True,actual_Rsh_exit_axial5_available=True,
            reshape_radial_mixed4_certified=False,reference_continuation_installed=False,
            actual_five_moment_patch_connected=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,
            input_hashes={**receipt['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Full reshape kernel axial5, physical primitive normalization and actual R110 history checks PASS',flush=True)
    return result


if __name__=='__main__':run()
