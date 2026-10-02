"""Independent original-switch integration and actual-history transport checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_inner_switch_profiles import (
    CompliantInnerSwitchProfiles, power_transport, IntervalTaylor, MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_inner_switch_profiles.json'


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    a=mp.exp(-1/t**2);b=mp.exp(-1/(1-t)**2)
    return a/(a+b)


def switch_integration_fixture():
    """Integrate source a,b directly for nonconstant current-R Dbar/Ebar.

    Finite fixture hb is NOT admitted for the actual source. In particular
    test the O(hb) term that is lost by omitting sigma in the second a.
    """
    with mp.workdps(75):
        h=mp.mpf('.025');tol=mp.mpf('1e-60');count=0;axial=0
        mass=mp.quad(sigma,[0,mp.mpf('.5'),1])
        if abs(mass-mp.mpf('.5'))>tol:raise ArithmeticError('Original sigma half integral failed')
        for t in (mp.mpf('.01'),mp.mpf('.2'),mp.mpf('.5'),mp.mpf('.85')):
            if abs(sigma(t)+sigma(1-t)-1)>tol:raise ArithmeticError('Original sigma symmetry failed')
        for z in (mp.mpf(0),mp.mpf('.4')):
            dz=mp.mpf('.03')+mp.mpf('.004')*z+mp.mpf('.002')*z*z
            D=lambda R:R*dz
            def a(q):
                radius=100*mp.exp(h*q)
                if q<=1:return h*D(radius)
                if q<=2:return h*D(radius)*(1-sigma(q-1))+mp.mpf('.8')*sigma(q-1)
                return mp.mpf('.8')
            JD=mp.quad(lambda q:D(100*mp.exp(h*q)),[0,1])
            JD+=mp.quad(lambda q:(1-sigma(q-1))*D(100*mp.exp(h*q)),[1,mp.mpf('1.5'),2])
            for phase in (mp.mpf('.5'),mp.mpf(1),mp.mpf('1.5'),mp.mpf(2),mp.mpf(3)):
                cuts=[mp.mpf(0)]+[mp.mpf(q) for q in (1,2) if q<phase]+[phase]
                direct=-h*mp.quad(a,cuts)/2
                stop=min(phase,1)
                dmass=mp.quad(lambda q:D(100*mp.exp(h*q)),[0,stop])
                smass=mp.mpf(0)
                if phase>1:
                    upper=min(phase,2)
                    dmass+=mp.quad(lambda q:(1-sigma(q-1))*D(100*mp.exp(h*q)),[1,upper])
                    smass=mp.quad(sigma,[0,upper-1])
                recovered=-h*h*dmass/2-mp.mpf('.4')*h*smass-mp.mpf('.4')*h*max(0,phase-2)
                if abs(direct-recovered)>tol:raise ArithmeticError('Original current-R switch integration failed')
                count+=1
                if phase>=2:
                    corrected=-mp.mpf('.4')*h*phase+mp.mpf('.6')*h-h*h*JD/2
                    if abs(direct-corrected)>tol:raise ArithmeticError('Post-switch .6hb term failed')
                    wrong=-mp.mpf('.4')*h*phase-h*h*JD/2
                    if abs(direct-wrong)<h/2:raise ArithmeticError('Fixture did not distinguish dropped sigma term')
                    count+=1
            F100=lambda z0:mp.exp(mp.mpf('.13')*z0+mp.mpf('.02')*z0*z0)
            barphi=lambda z0:1+mp.mpf('.03')*z0
            F0=lambda z0:mp.exp(mp.mpf('.07')*z0)
            E=lambda R,z0:(mp.mpf('.02')+mp.mpf('.003')*z0)*mp.sqrt(R)
            # On the first switch D=R*d(z) integrates exactly in logR.
            def F(q,z0):
                radius=100*mp.exp(h*q);d0=mp.mpf('.03')+mp.mpf('.004')*z0+mp.mpf('.002')*z0*z0
                return F100(z0)*mp.exp(-h*d0*(radius-100)/2)
            for phase in (mp.mpf('.5'),mp.mpf(1),mp.mpf('1.5'),mp.mpf(2),mp.mpf(3)):
                stop=min(phase,1)
                physical=-h*h*mp.quad(lambda q:(1-sigma(q))*mp.sqrt(100*mp.exp(h*q)/2)*F(q,z)*E(100*mp.exp(h*q),z),[0,stop])
                factored=-h*h*mp.quad(lambda q:(1-sigma(q))*(F(q,z)/(F0(z)*barphi(z)))*
                    (mp.sqrt(100*mp.exp(h*q)/2)*F0(z)*barphi(z)*E(100*mp.exp(h*q),z)),[0,stop])
                if abs(physical-factored)>tol:raise ArithmeticError('First-switch physical/factored axial integral disagrees')
                if phase==1:V1=physical
                if phase>1 and abs(physical-V1)>tol:raise ArithmeticError('b=0 did not retain exact first-switch velocity')
                axial+=1
            frozenJD=mp.mpf('1.5')*D(100)
            if abs(JD-frozenJD)<mp.mpf('.001'):
                raise ArithmeticError('Fixture does not distinguish current-R from frozen incoming Dbar')
        return dict(independent_angular_switch_integrations=count,
            independent_physical_factored_axial_integrations=axial,
            original_sigma_symmetries=4,original_sigma_half_integrals=1,
            dropped_sigma_and_frozen_Dbar_variants_rejected=True,
            finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def power_moment_fixture():
    with mp.workdps(90):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-75');count=0
        phi=lambda z:mp.exp(mp.mpf('.1')*z+mp.mpf('.02')*z*z)
        velocity=lambda z:mp.mpf('.3')+2*z+mp.mpf('.07')*z*z
        names=('H','mean','K','A','B','C')
        initials={n:(lambda z,j=j:mp.mpf(j+2)/7+mp.mpf(j+1)*z/50+z*z/100)
                  for j,n in enumerate(names)}
        def jet(fn,z):
            return IntervalTaylor(c,[c.mpf([d-tol,d+tol])/math.factorial(k)
                for k in range(6) for d in (mp.diff(fn,z,k),)])
        for z in (mp.mpf(0),mp.mpf('.4')):
            initial={n:jet(fn,z) for n,fn in initials.items()}
            for x in (mp.mpf(1),mp.mpf('1.005'),mp.mpf('1.1'),mp.mpf(2)):
                transported=power_transport(c,c.mpf(1/x),initial,jet(phi,z),jet(velocity,z))
                actual=dict(H=transported[MTH],mean=transported[MZ],K=transported[MTHZ],
                    A=transported[MZT]['axial'],B=transported[MZT]['swirl'],C=transported[MP])
                # Direct positive radial quadrature, without using the
                # closed weights in power_transport, and nonzero histories.
                angular=mp.quad(lambda q:2*q*q**(-mp.mpf(2)/5),[1,x])/x**2
                swirl=mp.quad(lambda q:q*q**(-mp.mpf(4)/5),[1,x])/x**2
                pressure=mp.quad(lambda q:q**(-mp.mpf(4)/5),[1,x])/x
                linear=mp.quad(lambda q:mp.mpf(1),[1,x])/x
                target=dict(H=lambda zz:initials['H'](zz)/x**2+phi(zz)*angular,
                    mean=lambda zz:initials['mean'](zz)/x+velocity(zz)*linear,
                    K=lambda zz:initials['K'](zz)/x**2+phi(zz)*velocity(zz)*angular,
                    A=lambda zz:initials['A'](zz)/x+velocity(zz)**2*linear,
                    B=lambda zz:initials['B'](zz)/x**2+phi(zz)**2*swirl,
                    C=lambda zz:initials['C'](zz)/x+phi(zz)**2*pressure)
                for name,row in actual.items():
                    for k in range(6):
                        expected=mp.diff(target[name],z,k);lo,hi=endpoints(row[k]*math.factorial(k))
                        if not lo-tol*1000<=expected<=hi+tol*1000:
                            raise ArithmeticError('Independent cumulative '+name+' derivative failed')
                        count+=1
        return dict(independently_integrated_nonzero_history_axial_derivatives=count,
            finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def structural_identities():
    x=s.symbols('x',positive=True);f,v=s.symbols('f v',real=True)
    incoming=s.symbols('H m K A B C');theta=1/x
    shapes=[theta**2*incoming[0]+s.Rational(5,4)*f*(theta**s.Rational(2,5)-theta**2),
            theta*incoming[1]+v*(1-theta),
            theta**2*incoming[2]+s.Rational(5,4)*f*v*(theta**s.Rational(2,5)-theta**2),
            theta*incoming[3]+v*v*(1-theta),
            theta**2*incoming[4]+s.Rational(5,6)*f*f*(theta**s.Rational(4,5)-theta**2),
            theta*incoming[5]+5*f*f*(theta**s.Rational(4,5)-theta)]
    powers=(2,1,2,1,2,1)
    rhs=(2*x*f*x**(-s.Rational(2,5)),v,2*x*f*v*x**(-s.Rational(2,5)),v*v,
         x*f*f*x**(-s.Rational(4,5)),f*f*x**(-s.Rational(4,5)))
    for i,(shape,power,source) in enumerate(zip(shapes,powers,rhs)):
        if s.simplify(s.diff(x**power*shape,x)-source)!=0 or s.simplify(shape.subs(x,1)-incoming[i])!=0:
            raise ArithmeticError('Power moment source or nonzero history identity failed')
    Cstar,Lambda,G,phi,z=s.symbols('Cstar Lambda G phi z',positive=True)
    F0=s.exp(-s.log(Cstar)-Lambda*G);u=s.sqrt(220)*F0*phi
    B=-Lambda*G+s.log(phi)+s.log(220)/2+s.log(1+z*z)
    if s.simplify(s.expand_log(s.log(Cstar*u*(1+z*z)),force=True)-B)!=0:
        raise ArithmeticError('Actual R110 logCstar cancellation failed')
    return dict(exact_power_primitive_RHS=6,exact_inherited_R2_histories=6,
                exact_R110_log_amplitude_cancellation=1,passed=True)


def run():
    with mp.workdps(280):
        receipt=json.loads((HERE/NAME).read_bytes())
        for name,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Switch source changed: '+name)
        c=MPIntervalContext();c.dps=240;proof_count=0;bound_count=0;moment_count=0
        sign=receipt['angular_sign_and_shear_source']
        gate=json.loads((HERE/sign['certificate']).read_bytes());ledger=json.loads((HERE/sign['ledger']).read_bytes())
        if not (gate['actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified']
                and gate['exact_implicit_exit_field_specified']
                and gate['admitted_inner_parameter_family_sha256']==receipt['admitted_inner_parameter_family_sha256']
                and ledger['admitted_inner_parameter_family_sha256']==receipt['admitted_inner_parameter_family_sha256']
                and endpoints(read_interval(c,ledger['decreasing_K_smallness_bounds']['initial_switch_a']))[1]<endpoints(c.mpf(4)/5)[0]):
            raise ValueError('Same-source positive Dbar and short-switch bounds missing')
        logh=read_interval(c,receipt['source_log_hb_enclosure'])
        if endpoints(logh)[1]>=endpoints(c.ln(c.ln(c.mpf(110)/100)/2))[0]:raise ArithmeticError('R2 is not below110')
        for proof in receipt['source_log_product_cap_proofs']:
            upper=read_interval(c,proof['input_absolute_upper']);extra=read_interval(c,proof['additional_source_log'])
            cap=read_interval(c,proof['cap']);computed=logh+extra+c.ln(upper)
            if endpoints(cap)[0]<=0 or endpoints(computed)[1]>endpoints(c.ln(cap))[0]:
                raise ArithmeticError('Actual source h/product logarithmic enclosure failed')
            proof_count+=1
        packets=[receipt[n] for n in ('actual_R100_inlet','first_switch','second_switch','whole_short_switches','actual_R110_inlet')]+receipt['samples']
        for packet in packets:
            arrays=[(k,packet[k]) for k in ('F_actual_over_F0_axial5_coefficients','F_actual_true_axial5_divided_by_F0',
                'Uz_actual_axial5_coefficients','actual_Q_axial4_coefficients','pressure_axis_axial5_coefficients',
                'pressure_increment_true_axial5_divided_by_R_F0_squared')]
            if 'actual_R110_log_shape_axial5_coefficients' in packet:
                arrays.append(('actual_R110_log_shape_axial5_coefficients',packet['actual_R110_log_shape_axial5_coefficients']))
            for name,row in arrays:
                if len(row)!=(5 if name=='actual_Q_axial4_coefficients' else 6):raise ValueError('Switch derivative order lost')
                for value in row:
                    lo,hi=endpoints(read_interval(c,value))
                    if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite switch profile')
                    bound_count+=1
            if endpoints(read_interval(c,packet['F_actual_over_F0_axial5_coefficients'][0]))[0]<=0:
                raise ArithmeticError('Actual switch positive swirl lost')
            for name,row in packet['actual_moment_shape_axial5_coefficients'].items():
                for values in (row.values() if name==MZT else (row,)):
                    if len(values)!=6:raise ValueError('Actual inherited moment order lost')
                    for value in values:
                        lo,hi=endpoints(read_interval(c,value))
                        if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite actual switch moment')
                        moment_count+=1
            if packet['switch_radial_mixed4_certified'] or packet['comparison_moments_substituted']:
                raise ValueError('Unbuilt scope or comparison history promoted')
        # Compare the exact inlet to the accepted parent receipt, retaining
        # its nonzero moments, pressure, velocity and derivative conventions.
        parent=json.loads((HERE/'lei_ren_part1_paper_compliant_inner_bridge_profiles.json').read_bytes())['terminal']
        start=receipt['actual_R100_inlet']
        for key in ('F_actual_over_F0_axial5_coefficients','Uz_actual_axial5_coefficients',
                    'actual_moment_shape_axial5_coefficients','pressure_axis_axial5_coefficients'):
            if start[key]!=parent[key]:raise ArithmeticError('R100 actual bridge history was changed')
        for value in start['log_F_actual_over_F100_coefficients']:
            if endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('R100 switch log not exactly0')
        if receipt['first_switch']['Uz_actual_axial5_coefficients']!=receipt['second_switch']['Uz_actual_axial5_coefficients']:
            raise ArithmeticError('Second switch introduced axial change')
        if receipt['second_switch']['Uz_actual_axial5_coefficients']!=receipt['actual_R110_inlet']['Uz_actual_axial5_coefficients']:
            raise ArithmeticError('Power interval introduced axial change')
        # Direct source h^2 check of the angular operator before capping.
        provider=CompliantInnerSwitchProfiles();inp=provider.inputs([-1,1]);order_count=0
        for k in range(6):
            operand=inp['comparison']['direction']['D_over_R'][k]*c.mpf('82.5')
            absolute=c.mpf(max(abs(v) for v in endpoints(operand)))
            if endpoints(2*logh+c.ln(absolute))[1]>endpoints(c.ln(inp['h2_bounds'][k]))[1]:
                raise ArithmeticError('Actual angular h^2 factor not enclosed')
            order_count+=1
        result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
            implicit_source_sha256=receipt['implicit_source_sha256'],
            switch_fixture=switch_integration_fixture(),power_moment_fixture=power_moment_fixture(),
            structural_identities=structural_identities(),
            actual_source_width_product_log_proofs_checked=proof_count,
            actual_switch_profile_bounds_checked=bound_count,actual_switch_moment_bounds_checked=moment_count,
            actual_angular_h_squared_orders_checked=order_count,exact_R100_bridge_history_retained=True,
            exact_zero_axial_change_after_first_switch=True,
            actual_100_110_switch_axial5_enclosures_available=True,
            actual_100_110_radial_recovery_axial4_available=True,actual_R110_inlet_axial5_available=True,
            switch_radial_mixed4_certified=False,long_reshape_with_this_inlet_installed=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,all_passed=True,
            input_hashes={**receipt['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original switch integration, inherited moment transport and actual R110 inlet checks PASS',flush=True)
    return result


if __name__=='__main__':run()
