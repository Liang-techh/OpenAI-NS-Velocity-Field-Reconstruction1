"""Independent implicit-coordinate fixture and source-scale/solenoidal checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    cartesian_source_row,time_source_row,ordinary_terms,micro_terms,log_row,
    zero_powers,shift,INDICES,UZ,UT,UR,P,BASES)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import positive_exp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'global_physical_assembly.json'


def source_identities():
    z,delta=s.symbols('Z delta',real=True)
    m,V,Vy,mz,Vz=s.symbols('m V Vy mz Vz',real=True)
    d=1-z*z;L=1-delta*z*z
    Q=(2*z*V-(1-delta)*z*m-d*mz)/L
    Qy=(2*z*Vy-(1-delta)*z*(V-m)-d*(Vz-mz))/L
    divergence=s.simplify(Qy+Q+((-1-delta)*z*V+d*Vz-2*z*Vy)/L)
    if divergence!=0:raise ArithmeticError('Exact original Mz recovery is not solenoidal')
    h,p,f,u,k,R,phi,lp,lf,lh=s.symbols('h p f u k R phi lp lf lh',real=True)
    old=h*lh+2*p*lp+2*f*lf+2*u*(s.log(2)/2+R/2+lf+phi-lp)-k*lh
    new=(h-k)*lh+(2*p-2*u)*lp+(2*f+2*u)*lf+2*u*phi+u*R+u*s.log(2)
    if s.expand(old-new)!=0:raise ArithmeticError('Microscopic linear/quadratic source-factor conversion failed')
    # This is FIXED-basepoint normalization. All full derivatives of the
    # source are in the input rows. Re-differentiating the unit is wrong.
    x=s.symbols('x',positive=True);G=s.Function('G')(x);A0=s.symbols('A0',nonzero=True)
    for n in range(5):
        if s.simplify(A0*s.diff(G/A0,x,n)-s.diff(G,x,n))!=0:raise ArithmeticError('Fixed source-unit derivative identity failed')
    return dict(exact_original_Mz_recovery_divergence_identity=True,
        micro_logu_Pstar_F0_phi_radial_factor_identity=True,fixed_basepoint_unit_derivative_identities=5,
        no_coordinate_rescaling_claimed_as_temporal_recursion=True,passed=True)


def materialize_finite(c,rows,logs,loglambda):
    total=c.mpf(0)
    for powers,coefficient in rows[0]:
        exponent=sum((logs[index]*power for index,power in enumerate(powers) if power),c.mpf(0))+rows[1]*loglambda
        total+=coefficient*c.exp(exponent)
    return total


def implicit_physical_fixture():
    """Finite streamfunction-derived field, independent lambda root, all C4."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-60')
        delta=mp.mpf('.03');tau=mp.mpf('.7');Z0=mp.mpf('.4');Y0=mp.mpf('.2');theta=mp.mpf('.7')
        lam=mp.sqrt(tau/(1-Z0*Z0));R=mp.exp(Y0);r=lam*mp.sqrt(2*R)
        x0=r*mp.cos(theta);y0=r*mp.sin(theta);z0=lam**(1-delta)*Z0
        def mass(y,z):return 4*z+mp.mpf('.3')*z*mp.exp(mp.mpf('.17')*y)+mp.mpf('.2')*z*z*mp.exp(-mp.mpf('.11')*y)
        def axial(y,z):return 4*z+mp.mpf('.351')*z*mp.exp(mp.mpf('.17')*y)+mp.mpf('.178')*z*z*mp.exp(-mp.mpf('.11')*y)
        def massz(y,z):return 4+mp.mpf('.3')*mp.exp(mp.mpf('.17')*y)+mp.mpf('.4')*z*mp.exp(-mp.mpf('.11')*y)
        def radial(y,z):
            Q=(2*z*axial(y,z)-(1-delta)*z*mass(y,z)-(1-z*z)*massz(y,z))/(1-delta*z*z)
            return mp.sqrt(mp.exp(y)/2)*Q
        functions={UR:radial,UZ:axial,
            UT:lambda y,z:mp.exp(-mp.mpf('.21')*y)*(1-mp.mpf('.2')*z+mp.mpf('.15')*z**3)+mp.exp(mp.mpf('.13')*y)*mp.cos(z),
            P:lambda y,z:mp.exp(-mp.mpf('.23')*y)*(1+z**4)+mp.exp(mp.mpf('.04')*y)*mp.sin(z)}
        beta={UR:mp.mpf(-1),UT:-1-delta,UZ:-1-delta,P:-2-2*delta}
        # Distinct fixed source units exercise all unit placement, not only
        # an implementation with every normalization set to one.
        amp_powers={UR:2,UT:1,UZ:-1,P:3};logs=(c.mpf(0),)*4+(c.ln(2),c.mpf(Y0),c.ln(2))
        amplitudes={label:{4:mp.mpf(power)} for label,power in amp_powers.items()};grids={}
        for label,fn in functions.items():
            unit=mp.mpf(2)**amp_powers[label];grid={}
            for k in range(5):
                for n in range(5-k):
                    v=mp.diff(fn,(Y0,Z0),(k,n))/unit
                    grid['y'+str(k)+'_Z'+str(n)]=c.mpf([v-tol,v+tol])
            grids[label]=ordinary_terms(grid)
        roots={}
        def physical(component,x,y,z,remaining=tau):
            key=(z,remaining,mp.mp.prec)
            if key not in roots:roots[key]=mp.findroot(lambda ll:ll*ll-ll**(2*delta)*z*z-remaining,lam,tol=mp.eps*16,verify=True)
            ll=roots[key];rr=mp.sqrt(x*x+y*y);Z=z/ll**(1-delta);Y=mp.log(rr*rr/(2*ll*ll))
            ur=ll**beta[UR]*functions[UR](Y,Z);ut=ll**beta[UT]*functions[UT](Y,Z)
            return {'ux':lambda:ur*x/rr-ut*y/rr,'uy':lambda:ur*y/rr+ut*x/rr,
                'uz':lambda:ll**beta[UZ]*functions[UZ](Y,Z),'p':lambda:ll**beta[P]*functions[P](Y,Z)}[component]()
        mapped={};count=0
        for i,j,b in INDICES:
            for component in ('ux','uy','uz','p'):
                parts=cartesian_source_row(c,grids,component,i,j,b,c.mpf(Z0),c.mpf(delta),c.mpf(mp.cos(theta)),c.mpf(mp.sin(theta)),amplitudes)
                value=sum((materialize_finite(c,rows,logs,c.ln(c.mpf(lam))) for rows in parts.values()),c.mpf(0))
                actual=mp.diff(lambda x,y,z:physical(component,x,y,z),(x0,y0,z0),(i,j,b))
                if not endpoints(value)[0]<=actual<=endpoints(value)[1]:raise ArithmeticError('Independent new Cartesian source map failed: '+str((component,i,j,b)))
                mapped[component,i,j,b]=value;count+=1
        time={label:materialize_finite(c,time_source_row(c,grids,label,c.mpf(Z0),c.mpf(delta),amplitudes),logs,c.ln(c.mpf(lam))) for label in functions}
        times={'ux':time[UR]*c.mpf(mp.cos(theta))-time[UT]*c.mpf(mp.sin(theta)),
            'uy':time[UR]*c.mpf(mp.sin(theta))+time[UT]*c.mpf(mp.cos(theta)),'uz':time[UZ],'p':time[P]}
        for component,value in times.items():
            actual=-mp.diff(lambda remaining:physical(component,x0,y0,z0,remaining),tau)
            if not endpoints(value)[0]<=actual<=endpoints(value)[1]:raise ArithmeticError('Independent fixed-x time source map failed')
        divergence=mapped['ux',1,0,0]+mapped['uy',0,1,0]+mapped['uz',0,0,1]
        if not endpoints(divergence)[0]<=0<=endpoints(divergence)[1]:raise ArithmeticError('Mapped divergence lost zero')
        actual_div=sum(mp.diff(lambda q:physical(component,*(q if i==axis else (x0,y0,z0)[i] for i in range(3))),
            (x0,y0,z0)[axis]) for axis,component in enumerate(('ux','uy','uz')))
        if abs(actual_div)>mp.mpf('1e-70'):raise ArithmeticError('Independent physical streamfunction fixture is not divergence-free')
        return dict(independent_implicit_root_cartesian_derivatives=count,independent_fixed_x_time_derivatives=4,
            independent_physical_divergence_checked=True,distinct_fixed_source_units_checked=True,
            finite_fixture_only=True,actual_core_point_coefficients_reconstructed=False,passed=True)


def microscopic_scale_fixture():
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110
        hb=c.mpf(-3000);lp=c.mpf(4000);lf=c.mpf(-5000);lr=c.mpf('.7');phi=c.mpf('.03')
        lu=c.ln(2)/2+lr/2+lf+phi-lp;oldlogs=(hb,2*lp,2*lf,2*lu)
        logs=(hb,lp,lf,phi,c.mpf(0),lr,c.ln(2));count=0
        for key in ((2,1,1,0),(3,0,0,1),(4,-1,1,1),(0,1,0,0)):
            ledger=[]
            for k in range(5):
                for n in range(5-k):ledger.append(dict(physical_row='Uz/s'+str(k)+'_Z'+str(n),terms=[dict(source_exponents=key,final_ordinary_coefficient=c.mpf('1.25'))]))
            for phase_to_y in (True,False):
                converted=micro_terms(dict(final_factored_physical_row_ledgers=ledger),'Uz',phase_to_y=phase_to_y)
                for (k,n),terms in converted.items():
                    true_log=sum((oldlogs[j]*power for j,power in enumerate(key)),c.mpf(0))-(k*hb if phase_to_y else 0)+c.ln(c.mpf('1.25'))
                    row=log_row(c,terms,logs,c.mpf(-1),c.mpf('-.5'));expected=true_log+c.mpf('.5')
                    actual=endpoints(expected);bound=endpoints(row['log_absolute_upper'])
                    # A directed upper bound can be strictly above the source
                    # log interval; interval overlap is not its acceptance rule.
                    if bound[1]<actual[0] or abs(bound[1]-actual[0])>mp.mpf('1e-65'):
                        raise ArithmeticError('Extreme phase/macro factor log conversion failed')
                    count+=1
        # The original phase source here is exp(-8000), but its fourth logR
        # derivative is exp(4000). Preserve that source through conversion.
        source_phase=-8000;resolved_cap=positive_exp(c,c.mpf(source_phase))
        if endpoints(resolved_cap)[0]!=0 or endpoints(c.mpf(source_phase)-4*hb)[0]!=4000:
            raise ArithmeticError('Extreme source fixture does not exercise pre-cap conversion')
        return dict(independent_signed_micro_and_macro_source_scale_rows=count,
            micro_logR_conversion_rows=60,macro_without_extra_width_conversion_rows=60,phase_cap_not_used_as_source=True,
            exact_logR_conversion_before_bound=True,passed=True)


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Global physical receipt source changed: '+name)
    identities=source_identities();fixture=implicit_physical_fixture();micro=microscopic_scale_fixture()
    c=MPIntervalContext();c.dps=240;count=0;terms=0
    charts=receipt['whole_source_chart_physical_maps']
    if len(charts)!=33:raise ValueError('Incomplete leading source chart coverage')
    for chart,packet in charts.items():
        if chart=='core':
            native=packet['core_axis_nonsingular_native_map']
            if not native['includes_axis_without_inverse_radius']:raise ValueError('Missing nonsingular core map')
            for components in packet['requested_core_spatial_log_bounds'].values():
                for parts in components.values():
                    for value in parts.values():
                        if value is not None and not mp.isfinite(endpoints(read_interval(c,value))[1]):raise ArithmeticError('Nonfinite core physical log bound')
                        count+=1
            for parts in packet['requested_core_time_log_bounds'].values():
                for value in parts.values():
                    if value is not None and not mp.isfinite(endpoints(read_interval(c,value))[1]):raise ArithmeticError('Nonfinite requested core time log bound')
                    count+=1
            continue
        if len(packet['physical_spatial_cartesian_mixed4'])!=35:raise ValueError('Cartesian multiindex coverage incomplete')
        if not packet['normalization_not_differentiated_twice'] or not packet['moving_cylindrical_basis_differentiated']:
            raise ValueError('Physical source units/basis were lost')
        if chart in ('bridge_first','bridge_second','switch_first','switch_second') and not packet['microscope_phase_to_logR_applied_before_source_bound']:
            raise ValueError('Micro physical derivatives used capped phase values')
        factored=chart in ('bridge_first','bridge_second','bridge_macro','switch_first','switch_second')
        mode='uncapped_factored_rows' if factored else 'provider_prebounded_mixed_rows'
        if packet['global_source_row_mode']!=mode:raise ValueError('Source row precision mode was overstated: '+chart)
        if chart.startswith(('bridge_','switch_')):
            source=packet['original_radius_source']
            if not (source['microscopic_radius_variation_retained_formally'] and source['absolute_radius_not_rounded_as_source'] and source['exact_log_source']):
                raise ValueError('Original positive microscopic radius source lost: '+chart)
            if not source['numeric_logR_bound_uses_positive_hb_enclosure'] or source['formal_hb_radius_correlation_evaluated']:
                raise ValueError('Conservative radius enclosure was promoted to exact correlated evaluation: '+chart)
        for components in list(packet['physical_spatial_cartesian_mixed4'].values())+[packet['first_fixed_x_physical_time_derivative']]:
            for parts in components.values():
                for row in parts.values():
                    if row['source_row_mode']!=mode or row['original_source_factors_combined_before_enclosure']!=factored:
                        raise ValueError('Provider caps were hidden in physical source metadata: '+chart)
                    if not row['exact_zero']:
                        if row['log_absolute_upper'] is None or not mp.isfinite(endpoints(read_interval(c,row['log_absolute_upper']))[1]):raise ArithmeticError('Nonfinite physical source log upper bound: '+chart)
                        if endpoints(read_interval(c,row['physical_lambda_exponent']))[1]>=0:raise ArithmeticError('Uniform lambda bound direction is invalid')
                    if row['exact_zero']!= (not row['terms']):raise ValueError('Exact-zero source semantics changed')
                    count+=1;terms+=len(row['terms'])
    axis=receipt['nonsingular_axis_map']['core_axis_nonsingular_native_map']
    if not axis['axis_only'] or not axis['includes_axis_without_inverse_radius']:raise ValueError('Axis was mapped by singular cylindrical templates')
    for flag in ('full_point_physical_field_evaluation','full_background_NS_validation','physical_energy_integral_certified',
        'admissible_stress_lift_constructed','independently_bounded_flat_remainder','temporal_recursion'):
        if receipt[flag]:raise ValueError('Physical assembly scope promoted: '+flag)
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],implicit_source_sha256=receipt['implicit_source_sha256'],
        source_and_divergence_identities=identities,independent_implicit_physical_fixture=fixture,independent_micro_scale_fixture=micro,
        whole_original_charts_checked=33,finite_upper_physical_source_rows_checked=count,retained_signed_source_terms_checked=terms,
        axis_nonsingular_source_map_checked=True,all_33_original_source_charts_physical_spatial4_time1_mapped=True,
        exact_source_divergence_identity_checked=True,full_point_physical_field_evaluation=False,
        full_background_NS_validation=False,physical_energy_integral_certified=False,admissible_stress_lift_constructed=False,
        independently_bounded_flat_remainder=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Common physical assembly PASS:140 independent Cartesian derivatives,4 time derivatives,120 micro/macro factor rows,source divergence identity',flush=True)
    return result


if __name__=='__main__':run()
