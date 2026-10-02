"""Independent centered physical transport and exact five-row normalization."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_reference_restore_profiles import (
    CompliantReferenceRestoreProfiles,IntervalTaylor,restoration_kernels,
    restore_centered,scale_history)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_frozen_comparison_field import axial_jet
from lei_ren_part1_paper_compliant_inner_bridge_profiles import symmetric
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_reference_restore_profiles.json'


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    a=mp.exp(-1/t**2);b=mp.exp(-1/(1-t)**2)
    return a/(a+b)


def restoration_fixture():
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-68');kernel_count=0;derivative_count=0
        rates={'mean':(mp.mpf(1),1),'mixed':(mp.mpf('1.6'),1),'square':(mp.mpf(1),2)}
        for box in (0,mp.mpf('.5'),1,[mp.mpf('.3'),mp.mpf('.8')]):
            bounded=restoration_kernels(c,c.mpf(box))
            points=(mp.mpf('.3'),mp.mpf('.55'),mp.mpf('.8')) if isinstance(box,list) else (mp.mpf(box),)
            for t in points:
                for name,(rate,power) in rates.items():
                    direct=mp.quad(lambda q:mp.exp(-rate*(t-q))*(1-sigma(q))**power,[0,t]) if t else mp.mpf(0)
                    lo,hi=endpoints(bounded[name])
                    if not lo-tol<=direct<=hi+tol:raise ArithmeticError('Full restoration kernel enclosure failed')
                    kernel_count+=1
        z=mp.mpf('.3');Efn=lambda zz:mp.mpf('.02')+zz/300+zz**3/2000
        names=('mean_error','angular_error','mixed_error','axial_square','swirl_error','pressure_error')
        functions={name:(lambda zz,i=i:mp.mpf(i+1)/100+zz/1000+zz**2/2000+zz**5/10000)
                   for i,name in enumerate(names)}
        def jet(fn):
            return IntervalTaylor(c,[c.mpf([q-tol,q+tol])/math.factorial(k)
                for k in range(6) for q in (mp.diff(fn,z,k),)])
        initial={n:jet(fn) for n,fn in functions.items()};E=jet(Efn)
        for t in (mp.mpf(0),mp.mpf('.5'),mp.mpf(1)):
            kernels=restoration_kernels(c,c.mpf(t));actual=restore_centered(c,initial,E,c.mpf(t),kernels)
            # Direct physical radial increments, Rz=1, u_ref=A(z)R^.1,
            # V=4z+E(z)*alpha(logR). The cancellation is done only AFTER
            # integrating these physical increments, independently of the
            # normalized recurrence used by restore_centered.
            R=mp.exp(t);amplitude=lambda zz:mp.mpf('1.7')*mp.exp(mp.mpf('-.8'))/(1+zz*zz)
            mass=mp.quad(lambda q:mp.exp(q)*(1-sigma(q)),[0,t]) if t else mp.mpf(0)
            mixed=mp.quad(lambda q:mp.exp(mp.mpf('1.6')*q)*(1-sigma(q)),[0,t]) if t else mp.mpf(0)
            quadratic=mp.quad(lambda q:mp.exp(q)*(1-sigma(q))**2,[0,t]) if t else mp.mpf(0)
            def centered_targets(zz):
                u0=amplitude(zz);uR=u0*R**mp.mpf('.1');v0=4*zz;e=Efn(zz)
                mz0=v0+functions['mean_error'](zz)
                h0=mp.sqrt(2)*u0*(mp.mpf('0.625')+functions['angular_error'](zz))
                k0=4*zz*h0+mp.sqrt(2)*u0*functions['mixed_error'](zz)
                axial0=v0*v0+8*zz*functions['mean_error'](zz)+functions['axial_square'](zz)
                mz=mz0+v0*(R-1)+e*mass
                h=h0+mp.quad(lambda q:mp.sqrt(2)*u0*mp.exp(mp.mpf('1.6')*q),[0,t]) if t else h0
                k=k0+4*zz*(h-h0)+mp.sqrt(2)*u0*e*mixed
                axial=axial0+v0*v0*(R-1)+8*zz*e*mass+e*e*quadratic
                return dict(mean_error=mz/R-4*zz,angular_error=h/(mp.sqrt(2)*R**mp.mpf('1.5')*uR)-mp.mpf('.625'),
                    mixed_error=(k-4*zz*h)/(mp.sqrt(2)*R**mp.mpf('1.5')*uR),
                    axial_square=axial/R-8*zz*mz/R+16*zz*zz,
                    swirl_error=functions['swirl_error'](zz)*R**mp.mpf('-1.2'),
                    pressure_error=functions['pressure_error'](zz)*R**mp.mpf('-.2'))
            for name,row in actual.items():
                for n in range(6):
                    expected=mp.diff(lambda zz:centered_targets(zz)[name],z,n)
                    lo,hi=endpoints(row[n]*math.factorial(n))
                    if not lo-tol*1000<=expected<=hi+tol*1000:
                        raise ArithmeticError('Independent physical centered '+name+' derivative failed')
                    derivative_count+=1
        return dict(independent_full_restoration_kernel_integrals=kernel_count,
            independently_integrated_physical_centered_axial_derivatives=derivative_count,
            finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def structural_identities():
    R,Rm,Am,z=s.symbols('R Rm Am z',positive=True);u=s.Function('u')(R)
    V=s.symbols('V',real=True);em,eh,ek,aa,eb,ep=[s.Function(n)(R) for n in ('em','eh','ek','aa','eb','ep')]
    m=4*z+em;H=s.Rational(5,8)+eh;K=4*z*H+ek
    A=16*z*z+8*z*em+aa;b=s.Rational(5,6)+eb;p=5+ep
    moments=[R*m,s.sqrt(2)*R**s.Rational(3,2)*u*H,
             s.sqrt(2)*R**s.Rational(3,2)*u*K,R*A-R*u*u*b/2,u*u*p/2]
    updates={s.diff(u,R):u/(10*R),s.diff(em,R):(V-4*z-em)/R,
        s.diff(eh,R):-s.Rational(8,5)*eh/R,s.diff(ek,R):(V-4*z-s.Rational(8,5)*ek)/R,
        s.diff(aa,R):((V-4*z)**2-aa)/R,s.diff(eb,R):-s.Rational(6,5)*eb/R,
        s.diff(ep,R):-ep/(5*R)}
    expected=[V,s.sqrt(2*R)*u,s.sqrt(2*R)*u*V,V*V-u*u/2,u*u/(2*R)]
    for row,rhs in zip(moments,expected):
        if s.simplify(s.diff(row,R).xreplace(updates)-rhs)!=0:
            raise ArithmeticError('Centered physical moment source identity failed')
    physical_centered=moments[3]-(16*z*z*R-s.Rational(5,12)*R*u*u)-8*z*(moments[0]-4*z*R)
    if s.simplify(physical_centered-(R*aa-R*u*u*eb/2))!=0:
        raise ArithmeticError('Centered d4 omitted reference/swirl term')
    x=s.symbols('x',real=True);ratio=s.exp(x);localu=Am*s.exp(x/10)
    # All normalization factors come from the physical (4.42) definitions.
    definitions=[(R*em)/Rm,(s.sqrt(2)*R**s.Rational(3,2)*u*ek)/(s.sqrt(2)*Rm**s.Rational(3,2)*Am),
        (s.sqrt(2)*R**s.Rational(3,2)*u*eh)/(s.sqrt(2)*Rm**s.Rational(3,2)*Am),
        (R*aa-R*u*u*eb/2)/(Rm*Am*Am),(u*u*ep/2)/(Am*Am)]
    targets=[s.exp(x)*em,s.exp(s.Rational(8,5)*x)*ek,s.exp(s.Rational(8,5)*x)*eh,
             s.exp(x)*aa/Am**2-s.exp(s.Rational(6,5)*x)*eb/2,s.exp(x/5)*ep/2]
    # The center functions refer to the same current radius on both sides.
    # Replace them by scalars before changing only the radial prefactors.
    center_symbols=dict(zip((em,eh,ek,aa,eb,ep),s.symbols('em0 eh0 ek0 aa0 eb0 ep0')))
    for row,target in zip(definitions,targets):
        scaled=row.xreplace(center_symbols).subs({R:Rm*ratio,u:localu},simultaneous=True)
        if s.simplify(scaled-target.xreplace(center_symbols))!=0:
            raise ArithmeticError('Original normalized defect row scale failed')
    # After restoration V=4z, physical defect functions are constant in R.
    for row in (R*em,s.sqrt(2)*R**s.Rational(3,2)*u*eh,
                s.sqrt(2)*R**s.Rational(3,2)*u*ek,physical_centered,u*u*ep/2):
        if s.simplify(s.diff(row,R).xreplace(updates).subs(V,4*z))!=0:
            raise ArithmeticError('Unpatched defect was not constant after restoration')
    return dict(exact_centered_physical_primitive_RHS=5,exact_reference_and_swirl_d4_cancellation=1,
        exact_original_defect_row_normalizations=5,exact_unpatched_defect_radial_constants=5,passed=True)


def source_binding(provider,inp):
    c=provider.ctx;Z=c.mpf([-1,1]);switch=provider.reshape.switch.inputs(Z)
    core=provider.core.normalized_jets(4,Z);z=inp['z'];base=z*4+provider.core.j
    if any(endpoints(a)!=endpoints(b) for a,b in zip(base.coefficients,core['source']['u'].coefficients)):
        raise ArithmeticError('Actual core axial baseline changed')
    core_v=axial_jet(c,core['Uz']);psi=axial_jet(c,core['Psi'])
    # Bind the ordinary derivatives in the parent's original arithmetic
    # order, before conversion to Taylor coefficients introduces rounding.
    for k in range(6):
        key='rho0_Z'+str(k)
        expected=base[k]*math.factorial(k)+provider.core.epsilon*core['Psi'][key]
        if endpoints(core['Uz'][key])!=endpoints(expected):
            raise ArithmeticError('Actual core axial velocity is not 4Z+j+epsilon*Psi')
    cap=provider.reshape.switch.bridge.cap
    three=symmetric(c,sum((c.mpf([0,endpoints(cap)[1]]) for _ in range(3)),c.mpf(0)))
    for k in range(6):
        if endpoints(switch['v'][k])!=endpoints(core_v[k]+three):
            raise ArithmeticError('Actual bridge three-component drive cover changed')
        lo,hi=endpoints(switch['velocity_increment'][k]);boundlo,boundhi=endpoints(three)
        if lo<boundlo or hi>boundhi:raise ArithmeticError('First switch exceeds three drive caps')
        if endpoints(switch['v_cover'][k])!=endpoints(switch['v'][k]+switch['velocity_increment'][k]):
            raise ArithmeticError('Actual first-switch axial history changed')
        if endpoints(inp['parent']['Uz_actual_axial5_coefficients'][k])!=endpoints(switch['v_cover'][k]):
            raise ArithmeticError('Actual Rsh velocity did not retain first-switch history')
    raw=psi*provider.core.epsilon+provider.core.j
    raw+=IntervalTaylor(c,[symmetric(c,cap*6)]*6)
    for k in range(6):
        lo,hi=endpoints(inp['E'][k]);boundlo,boundhi=endpoints(raw[k])
        if lo<boundlo or hi>boundhi:raise ArithmeticError('Centered E is outside actual six-cap cover')
    # Original alpha_box uses b/(a+b), sigma_jets uses a/(a+b),
    # with a=exp(-1/t^2), b=exp(-1/(1-t)^2). This exact identity
    # holds throughout the open interval; both use exact flat endpoints.
    a,b=s.symbols('a b',positive=True)
    if s.simplify(b/(a+b)-(1-a/(a+b)))!=0:raise ArithmeticError('Original cutoff complement identity failed')
    cutoff_count=0
    for t in ('0','.125','.25','.5','.75','.875','1'):
        point=c.mpf(t);old=alpha_box(c,point+1,c.mpf(1));new=1-sigma_jets(c,point)[0]
        reference=1-sigma(mp.mpf(t));ol,oh=endpoints(old);nl,nh=endpoints(new)
        if not (ol-mp.mpf('1e-230')<=reference<=oh+mp.mpf('1e-230') and
                nl-mp.mpf('1e-230')<=reference<=nh+mp.mpf('1e-230')):
            raise ArithmeticError('Original restoration cutoff implementations differ')
        cutoff_count+=1
    return dict(exact_actual_core_bridge_switch_Rsh_axial_source_chain=True,
        actual_six_drive_caps_checked_through_axial_order5=True,
        same_original_cutoff_complement_identity=True,original_cutoff_points_checked=cutoff_count,passed=True)


def run():
    with mp.workdps(280):
        receipt=json.loads((HERE/NAME).read_bytes());c=MPIntervalContext();c.dps=240
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Reference restore source changed: '+path)
        proof_count=0;bounds=0;moments=0
        core=json.loads((HERE/'lei_ren_part1_paper_compliant_core_transfer.json').read_bytes())
        logP=read_interval(c,core['logPstar']);logcap=-1000*c.ln(10)-2*logP
        if endpoints(logcap)!=endpoints(read_interval(c,receipt['positive_log_cap'])):
            raise ArithmeticError('Source-sensitive transport cap changed')
        for proof in receipt['retained_source_transport_cap_proofs']:
            rate=read_interval(c,proof['rate']);gap=read_interval(c,proof['source_gap_enclosure'])
            upper=read_interval(c,proof['input_absolute_upper']);computed=-rate*gap+c.ln(upper)
            if endpoints(rate)[0]<=0 or endpoints(gap)[0]<0 or endpoints(computed)[1]>endpoints(logcap)[0]:
                raise ArithmeticError('Actual factored transport cap inequality failed')
            proof_count+=1
        packets=[receipt[n] for n in ('whole_reference','actual_Rsh_inlet','actual_Rz_inlet',
            'whole_axial_restore','restore_start','restore_end','unpatched_Rm','unpatched_Rh')]+receipt['samples']
        for packet in packets:
            for name in ('log_Utheta_over_Pstar_axial5_coefficients','Utheta_true_axial5_divided_by_current_Utheta',
                         'Uz_actual_axial5_coefficients','actual_Q_axial4_coefficients','pressure_axis_axial5_coefficients',
                         'actual_E_V110_minus_4Z_axial5_coefficients'):
                row=packet[name]
                if len(row)!=(5 if name=='actual_Q_axial4_coefficients' else 6):raise ValueError('Reference restore derivative order lost')
                for value in row:
                    lo,hi=endpoints(read_interval(c,value))
                    if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite reference restore profile')
                    bounds+=1
            for group in ('actual_centered_moment_axial5_coefficients','actual_normalized_moment_shape_axial5_coefficients',
                          'actual_true_moment_derivatives_divided_by_current_prefactor'):
                for row in packet[group].values():
                    if len(row)!=6:raise ValueError('Actual centered moment derivative order lost')
                    for value in row:
                        lo,hi=endpoints(read_interval(c,value))
                        if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite centered moment')
                        moments+=1
            if packet['moment_reset'] or packet['radial_mixed4_certified'] or packet['actual_five_moment_patch_connected']:
                raise ValueError('Unbuilt scope or reference moment reset promoted')
        provider=CompliantReferenceRestoreProfiles();inp=provider.inputs([-1,1])
        start=receipt['actual_Rsh_inlet'];parent=json.loads((HERE/'lei_ren_part1_paper_compliant_long_reshape_profiles.json').read_bytes())['actual_Rsh_exit']
        if start['pressure_axis_axial5_coefficients']!=parent['pressure_axis_axial5_coefficients']:
            raise ArithmeticError('Original analytic pressure datum changed')
        for name,jet in inp['source_centered_Rsh'].items():
            if any(endpoints(read_interval(c,row))!=endpoints(value) for row,value in
                   zip(start['actual_centered_moment_axial5_coefficients'][name],jet.coefficients)):
                raise ArithmeticError('Actual Rsh centered histories were not copied exactly')
        for name in ('actual_centered_moment_axial5_coefficients','Uz_actual_axial5_coefficients'):
            # Compare exact interval endpoints, not overlap, at restore inlet.
            left=receipt['actual_Rz_inlet'][name];right=receipt['restore_start'][name]
            groups=left.keys() if isinstance(left,dict) else (None,)
            for group in groups:
                a=left[group] if group is not None else left;b=right[group] if group is not None else right
                if any(endpoints(read_interval(c,x))!=endpoints(read_interval(c,y)) for x,y in zip(a,b)):
                    raise ArithmeticError('Restoration changed exact Rz inlet history')
        for key in ('restore_end','unpatched_Rm','unpatched_Rh'):
            row=receipt[key]['Uz_actual_axial5_coefficients']
            exact=(c.mpf([-4,4]),c.mpf(4),c.mpf(0),c.mpf(0),c.mpf(0),c.mpf(0))
            if any(endpoints(read_interval(c,a))!=endpoints(b) for a,b in zip(row,exact)):
                raise ArithmeticError('Original restoration did not end at exact4Z')
        defects=receipt['actual_five_defects_at_Rh']['actual_normalized_five_defect_axial5_coefficients'];defect_count=0
        if len(defects)!=5 or any(len(row)!=6 for row in defects):raise ValueError('Actual defect five-row axial5 shape lost')
        for i,row in enumerate(defects):
            for value in row:
                lo,hi=endpoints(read_interval(c,value))
                if not (lo<=hi and mp.isfinite(lo) and mp.isfinite(hi)):raise ArithmeticError('Nonfinite actual normalized defect')
                defect_count+=1
            if i in (0,1,3) and endpoints(read_interval(c,row[0]))[0]<=0:
                raise ArithmeticError('Actual signed d1/d2/d4 source positivity lost')
        E=inp['E'];rho=inp['rho']
        ledger=provider.reshape.switch.bridge.records['K1_ledger']
        gate=provider.reshape.switch.bridge.records['global_exit_certificate']
        source_rho=read_interval(c,ledger['rho_core_C2_bound'])+read_interval(c,gate['rho_bridge_C2_upper'])
        if endpoints(rho)!=endpoints(source_rho):raise ArithmeticError('Centered E C2 operands changed')
        # Compare the actual cover with the source theorem cover directly.
        # Subtracting another uncertain copy of j loses its correlation.
        for k in range(3):
            source_cover=c.mpf([-endpoints(rho)[1],endpoints(rho)[1]])/math.factorial(k)
            if k==0:source_cover+=provider.core.j
            lo,hi=endpoints(E[k]);boundlo,boundhi=endpoints(source_cover)
            if lo<boundlo or hi>boundhi:
                raise ArithmeticError('Centered E C2 source theorem lost')
        result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
            implicit_source_sha256=receipt['implicit_source_sha256'],reference_join_family_sha256=receipt['reference_join_family_sha256'],
            restoration_fixture=restoration_fixture(),structural_identities=structural_identities(),
            actual_source_binding=source_binding(provider,inp),
            actual_source_factored_transport_caps_checked=proof_count,
            actual_velocity_pressure_profile_bounds_checked=bounds,actual_centered_moment_bounds_checked=moments,
            actual_normalized_defect_bounds_checked=defect_count,exact_Rsh_Rz_histories_and_P0_retained=True,
            original_axial_restoration_reaches_exact_4Z=True,negative_inherited_defect_terms_retained=True,
            actual_reference_continuation_axial5_enclosures_available=True,actual_axial_restoration_axial5_enclosures_available=True,
            actual_unpatched_Rm_Rh_moment_axial5_enclosures_available=True,
            actual_five_defect_profile_axial5_enclosures_available=True,actual_radial_recovery_axial4_available=True,
            radial_mixed4_certified=False,all_radial_endpoint_joins_certified=False,actual_five_moment_patch_connected=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,all_passed=True,input_hashes={**receipt['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual centered reference/restore history, normalized five defects and physical primitive checks PASS',flush=True)
    return result


if __name__=='__main__':run()
