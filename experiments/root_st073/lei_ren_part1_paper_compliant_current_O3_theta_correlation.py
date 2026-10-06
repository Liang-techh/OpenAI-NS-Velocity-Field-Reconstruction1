"""Actual O3 power theta source correlation; no O3 cone admission yet."""
import importlib
import json
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_pulse_entrance_cone import (
    CurrentPulseEntranceCone,HERE,PREFIX,sha,pack,encode,source_precision,OPEN)
from lei_ren_part1_paper_compliant_current_original_cone import SourceAST

NAME=PREFIX+'current_O3_theta_correlation.json'


def original_full_theta_identity():
    """Replay the actual coefficient operator, keeping incoming M and K."""
    asts=SourceAST();z,mu,delta,t=s.symbols('Z mu delta t',real=True)
    U,Ps,R=s.symbols('actual_U0 Pstar R',positive=True)
    M,K,X0=s.symbols('actual_M0 actual_K0 actual_Xpre',real=True)
    r=1-mu;C=1/(1+z*z);L=1-delta*z*z;b=(1-delta)/2;k=1-delta/2
    X=1/r+(X0-1/r)*s.exp(-r*t);u=U*s.exp(-(s.Rational(1,2)+mu)*t)
    B=Ps*u;zero=s.Integer(0)
    m=M*z*s.exp(-t)/(Ps*u*C);n=K*z*C*s.exp(-3*t/2)/(Ps*(u*C)**2)
    shifted=importlib.import_module(PREFIX+'collar_stress_C3').shifted_rows
    product=importlib.import_module(PREFIX+'collar_Gamma_C4').product_rows
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    mrows=[m*(-(s.Rational(1,2)-mu))**j for j in range(5)]
    nrows=[n*(-(s.Rational(1,2)-2*mu))**j for j in range(5)]
    coefficient=asts.replay('pulse_end_stress_C3','pulse_coefficients',dict(
        axial_derivative=lambda v:s.diff(v,z),product_rows=product,shifted_rows=shifted,
        mp=c))(delta,mu,z,C,X,[zero]*5,mrows,nrows,[zero]*5,[zero]*5,[zero]*5)
    # The replay concerns theta only; energy, loss and pressure slots are
    # absent from that operator. No axial/remainder row is inferred from
    # their zero placeholders here.
    checks={}
    def zero_identity(name,a,b):
        if s.cancel(s.expand_power_exp(a-b))!=0:raise ArithmeticError('O3 full theta identity failed: '+name)
        checks[name]=True
    theta=coefficient['theta']
    A=k*C+2*b*z*z*C*C
    N=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
    original=theta['equilibrium']['shape'][0]+theta['signed_original_memory']['shape'][0] \
        +B*theta['meridional_transport']['shape'][0]
    reduced=(A*X-C)/L+C*M*s.exp(-t)+N*(K/U)*s.exp(-r*t)
    zero_identity('actual_three_inertial_sectors_reduce_together',original,reduced)
    equilibrium=(A/r-C)/L
    zero_identity('actual_positive_equilibrium_numerator',equilibrium,
        ((mu-delta/2)*C+(1-delta)*z*z*C*C)/(r*L))
    initial=(A*X0-C)/L+C*M+N*K/U
    memory_transport=A*(X0-1/r)/L+N*K/U
    drift=equilibrium+(initial-equilibrium)*s.exp(-t) \
        +s.exp(-t)*(s.exp(mu*t)-1)*memory_transport
    zero_identity('actual_initial_theta_and_expm1_drift_source_identity',reduced,drift)
    zero_identity('actual_axis_full_inertial_theta',reduced.subs(z,0),
        k*X-1+M*s.exp(-t)-(K/U)*s.exp(-r*t))
    radial=theta['radial_shear']['shape'][0]/R
    zero_identity('actual_radial_shear_remains_separate_inverse_R_mode',radial,-2*(1+mu)*C/R)
    zero_identity('actual_full_theta_including_radial_shear',original+radial,reduced-2*(1+mu)*C/R)
    # Original O2 axial history: K/U=M+4*(X-1). Its zero-mu
    # axis reduction shows why throwing away incoming transport fails.
    zero_identity('O2_axial_axis_correlation_at_zero_mu_delta',
        initial.subs({z:0,mu:0,delta:0,K:U*(M+4*(X0-1))}),3*(1-X0))
    zero_identity('actual_expm1_drift_zero_at_power_inlet',drift.subs(t,0),initial)
    for name,part in theta.items():
        if name=='radial_shear':expected=(-.5,1,0,0)
        elif name=='meridional_transport':expected=(.5,2,1,0)
        else:expected=(.5,1,0,1 if name=='signed_original_memory' else 0)
        if tuple(part['mode'])!=expected:raise ValueError('Actual original O3 theta mode changed')
        checks['actual_'+name+'_mode']=True
    return dict(identities=checks,
        coefficient_symbols=dict(U='Actual scalar swirl at O3 power phase0',M='Actual raw m=M*Z at phase0',
            K='Actual raw k=K*Z/(1+Z^2) at phase0',Xpre='Actual scalar h/u at phase0',t='Tw*phase'),
        formulas=dict(C=str(C),L=str(L),A=str(A),N=str(N),
            equilibrium=str(equilibrium),initial_full_inertial_theta=str(initial),
            correlated_memory_transport=str(memory_transport),
            full_inertial_theta=str(reduced),correlated_expm1_form=str(drift),radial_shear=str(radial)),
        zero_axial_placeholders_not_claimed_as_current_energy_pressure_or_remainder=True,
        order_zero_full_theta_source_only=True,input_hashes=asts.hashes,passed=True)


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseEntranceCone()
    if type(field) is not CurrentPulseEntranceCone or not field.acceptance_loaded:
        raise ValueError('Checked current entrance source graph required')
    field.assert_graph();owner=field.registry.owners['incoming'];owner.assert_graph()
    original=importlib.import_module(PREFIX+'current_pulse_entrance_incoming_background_tensor')
    source=field.bindings['current_actual_incoming_function_proof']
    if owner.actual_power.__func__ is not original.CurrentPulseEntranceIncomingBackgroundTensor.actual_power \
            or owner.actual_power.__func__.__globals__['pulse_coefficients'] is not original.pulse_coefficients:
        raise ValueError('Actual current O3 power coefficient operator differs')
    if not source['passed'] or not all(source['identities'].values()):raise ValueError('Current complete production function proof required')
    asts=SourceAST()
    asts.expression('current_pulse_entrance_incoming_background_tensor','actual_power','pre',wanted='self.physical.pre.power(Z,phase)')
    asts.expression('current_pulse_entrance_incoming_background_tensor','actual_power','rows',
        wanted='pulse_coefficients(delta,mu,z,C,X,zeros,m,n,e,zeros,P)')
    asts.expression('current_pulse_entrance_incoming_background_tensor','actual_power','m1',wanted="copy_jet(c,raw['m'][0])*invP/u")
    asts.expression('current_pulse_entrance_incoming_background_tensor','actual_power','m2',wanted="copy_jet(c,raw['k'][0])*invP/(u*u)")
    theorem=original_full_theta_identity()
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_actual_full_O3_theta_source_identity=theorem,
        current_pre_power_function_and_native_coefficient_binding=source,
        actual_current_O3_power_input_and_normalization_AST_bound=True,
        actual_full_incoming_M_K_transport_not_dropped=True,
        original_P0_energy_axial_stress_and_remainder_unchanged=True,
        current_O3_power_whole_cone_certified=False,current_O3_transition_whole_cone_certified=False,
        uniform_positive_initial_theta_and_drift_bounds_still_required=True,
        initial_Xpre_not_confused_with_final_pulse_Xp=True,
        actual_current_tensor_region_count=33,current_strict_nonzero_regions=14,
        scope='Current source-functional order-zero O3 power theta reduction only. Complete incoming M/K, actual memory and separate inverse-R shear retained. No sign bound, whole O3 cone, higher derivative estimate, wave field or coefficient-recursion admission.',
        input_hashes={**field.hashes,**asts.hashes,**theorem['input_hashes'],
            PREFIX+'current_O3_theta_correlation.py':sha(PREFIX+'current_O3_theta_correlation.py')},
        **dict.fromkeys(OPEN,False),all_identities_passed=True)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current full O3 power theta correlation constructed; whole O3 cone still open',flush=True)
    return result


if __name__=='__main__':run()
