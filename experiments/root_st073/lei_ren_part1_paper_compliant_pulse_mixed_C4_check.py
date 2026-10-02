"""Independent closed-form mixed physical derivatives and actual O.4 gates."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import transport_mixed
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_pulse_mixed_C4.json'
UZ='Uz_over_Pstar_without_common_theta_radial_factor'
UT='Utheta_over_Pstar_without_common_theta_radial_factor'
UR='Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor'
P='P_over_Pstar_squared'


def independent_closed_form_fixture():
    """Differentiate actual physical components, not the transport algorithm."""
    with mp.workdps(90):
        c=MPIntervalContext(); c.dps=100; radius=mp.mpf('1e-70')
        mu,delta,a=mp.mpf('.04'),mp.mpf('.03'),mp.mpf('.17')
        bp=mp.mpf('.5')+mu; rate=1-mu; prate=1+2*mu
        ap=lambda z:mp.sqrt(1+z*z)+mp.mpf('.02')*z
        apz=lambda z:z/mp.sqrt(1+z*z)+mp.mpf('.02')
        theta=lambda z:1/(1+z*z)
        b=lambda y,z:ap(z)*mp.exp(a*y)
        def moment(y,z,row):
            lam=mp.mpf('.5')-row*mu
            initial=mp.mpf('.03')*(z+z**3) if row==1 else mp.mpf('.02')*(1+z*z)
            return initial*mp.exp(-lam*y)+ap(z)*(mp.exp(a*y)-mp.exp(-lam*y))/(a+lam)
        def momentz(y,z):
            lam=mp.mpf('.5')-mu
            return mp.mpf('.03')*(1+3*z*z)*mp.exp(-lam*y)+apz(z)*(mp.exp(a*y)-mp.exp(-lam*y))/(a+lam)
        X=lambda y,z:1/rate+(mp.mpf('.7')-1/rate)*mp.exp(-rate*y)
        energy=lambda y,z:(1+z*z)*mp.exp(2*mu*y)+ap(z)**2*(mp.exp(2*a*y)-mp.exp(2*mu*y))/(2*a-2*mu)+(1-mp.exp(2*mu*y))/(4*mu)
        pressure=lambda y,z:mp.mpf('.1')/(1+z*z)**2+theta(z)**2*(1-mp.exp(-prate*y))/(2*prate)
        A=lambda y,z:(2*z*b(y,z)-(1-delta)*z*moment(y,z,1)
            -(1-z*z)*(momentz(y,z)-2*z*moment(y,z,1)/(1+z*z)))/(1-delta*z*z)
        physical={UZ:lambda y,z:theta(z)*mp.exp(-bp*y)*b(y,z),
            UT:lambda y,z:theta(z)*mp.exp(-bp*y),
            UR:lambda y,z:theta(z)*mp.exp(-mu*y)*A(y,z),P:pressure}
        field=CompliantPulseRadialC4.__new__(CompliantPulseRadialC4)
        field.ctx=c; field.mu=c.mpf(mu); field.delta=c.mpf(delta); field.prate=c.mpf(prate)
        counts={}; y=mp.mpf('.3')
        for Z in (mp.mpf(-1),mp.mpf(0),mp.mpf('.3'),mp.mpf(1)):
            def jet(f):
                return IntervalTaylor(c,[c.mpf([v-radius,v+radius]) for v in mp.taylor(f,Z,5)])
            u=jet(theta); B=jet(lambda z:b(y,z)); Brows=[B*a**k for k in range(5)]
            point=dict(Mz_over_R_Utheta=jet(lambda z:moment(y,z,1)),
                Mtheta_z_over_sqrt2_R_3half_Utheta_squared=jet(lambda z:moment(y,z,2)),
                Mtheta_over_sqrt2_R_3half_Utheta=jet(lambda z:X(y,z)),
                Mztheta_over_R_Utheta_squared=jet(lambda z:energy(y,z)),
                pressure=dict(P_over_Pstar_squared=jet(lambda z:pressure(y,z)),
                    P_y_over_Pstar_squared=u*u*(mp.exp(-prate*y)/2)))
            actual=transport_mixed(field,c.mpf(Z),point,Brows,u)
            for label,function in physical.items():
                factor=mp.mpf(1) if label==P else mp.exp((-mu if label==UR else -bp)*y)
                for k in range(5):
                    for n in range(5-k):
                        target=mp.diff(function,(y,Z),(k,n))/factor
                        value=actual['physical_mixed_derivatives_total_order_le4'][label]['y'+str(k)+'_Z'+str(n)]
                        lo,hi=endpoints(value)
                        if not lo<=target<=hi:raise ArithmeticError('Independent physical mixed derivative failed: '+label+' '+str((Z,k,n)))
                        counts[label]=counts.get(label,0)+1
            primitives=actual['primitive_y_derivative_Taylor']
            for label,function in (
                ('Mz_over_R_Utheta',lambda yy,z:moment(yy,z,1)),
                ('Mtheta_z_over_sqrt2_R_3half_Utheta_squared',lambda yy,z:moment(yy,z,2)),
                ('Mtheta_over_sqrt2_R_3half_Utheta',X),('Mztheta_over_R_Utheta_squared',energy)):
                for k in range(5):
                    for n in range(5-k):
                        target=mp.diff(function,(y,Z),(k,n))/math.factorial(n)
                        lo,hi=endpoints(primitives[label][k][n])
                        if not lo<=target<=hi:raise ArithmeticError('Independent mixed primitive derivative failed: '+label)
                        counts[label]=counts.get(label,0)+1
        return dict(independent_closed_form_mixed_derivative_counts=counts,
            actual_Md40_source_admission=False,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Mixed C4 source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    old=json.loads((HERE/'lei_ren_part1_paper_compliant_pulse_radial_C4.json').read_bytes())
    if endpoints(read(r['whole_gap_box_overlap_xi']))[1]>mp.mpf('12.0001'):
        raise ArithmeticError('Uniform gap boxes do not cover the coordinate overlap')
    groups=('main_samples','entrance_samples','gap_samples','gap_end_samples','end_samples')
    samples=[p for group in groups for p in r[group]]+[r[key] for key in
        ('whole_Z_terminal','whole_Z_main','whole_Z_end_support','whole_Z_support_crossing',
         'whole_Z_entrance_box','whole_Z_main_box','whole_Z_exit_box',
         'whole_Z_gap_main_box','whole_Z_gap_end_box','whole_Z_end_box')]
    checks=0
    def overlaps(a,b,label):
        a,b=read(a),read(b)
        if max(endpoints(a)[0],endpoints(b)[0])>min(endpoints(a)[1],endpoints(b)[1]):
            raise ArithmeticError('Same-source mixed derivatives disjoint: '+label)
    for point in samples:
        grid=point['physical_mixed_derivatives_total_order_le4']
        for label in (UZ,UT,UR,P):
            if len(grid[label])!=15:raise ValueError('Mixed multiindex coverage incomplete')
            for value in grid[label].values():
                if any(not mp.isfinite(v) for v in endpoints(read(value))):raise ArithmeticError('Mixed derivative nonfinite')
                checks+=1
        if not point['pulse_all_mixed_derivatives_total_order_le4_available'] or not point['divergence_preserved_by_Mz_recovery']:
            raise ValueError('Original source mixed recovery missing')
        for key in ('full_pulse_C4_installed','full_outer_C4_certified','whole_outer_cone_certified','temporal_recursion',
                    'uniform_pulse_C4_chart_interface_certificate_available','cartesian_spatial_derivatives_certified'):
            if point[key]:raise ValueError('Mixed coordinate scope promoted to incomplete claim: '+key)
    comparisons=0
    for group in groups:
        for point,prior in zip(r[group],old[group]):
            for key in ('Uz_over_Utheta','Uz_y_over_Utheta','Ur_over_sqrt_R_over_2_Utheta',
                        'Ur_y_over_sqrt_R_over_2_Utheta','Ur_Z_over_sqrt_R_over_2_Utheta','Ur_yZ_over_sqrt_R_over_2_Utheta'):
                for a,b in zip(point[key]['coefficients'],prior[key]['coefficients']):
                    overlaps(a,b,key); comparisons+=1
    interfaces={}
    for label,left,right in (
        ('main_gap_xi11',r['main_samples'][-1],r['gap_samples'][0]),
        ('gap_coordinate_overlap',r['gap_samples'][2],r['gap_end_samples'][0]),
        ('gap_end_s_minus4',r['gap_end_samples'][-1],r['end_samples'][0])):
        for component in (UZ,UT,UR,P):
            for index,a in left['physical_mixed_derivatives_total_order_le4'][component].items():
                overlaps(a,right['physical_mixed_derivatives_total_order_le4'][component][index],label+' '+component+' '+index)
        interfaces[label]=True
    whole=r['whole_Z_terminal']['physical_mixed_derivatives_total_order_le4']
    for component in (UZ,UR):
        for value in whole[component].values():
            if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Terminal mixed zero not exact')
    # Terminal energy remains the selected positive complete future energy;
    # taking more radial derivatives must not zero that accumulated history.
    point=r['whole_Z_terminal']
    if endpoints(read(point['Mztheta_over_R_Utheta_squared']['coefficients'][0]))[0]<=0:
        raise ArithmeticError('Positive terminal future energy lost')
    for a,b in zip(point['Mztheta_over_R_Utheta_squared']['coefficients'],point['positive_terminal_future_energy_Taylor']['coefficients']):
        overlaps(a,b,'terminal selected energy')
    fixture=independent_closed_form_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        independent_closed_form_mixed_fixture=fixture,finite_actual_mixed_derivative_enclosures_checked=checks,
        earlier_actual_source_derivative_comparisons=comparisons,directed_mixed_chart_overlaps=interfaces,
        whole_Z_terminal_mixed_zero_and_positive_energy_checked=True,
        pulse_all_mixed_derivatives_total_order_le4_available=True,original_flat_pulse_support_derivatives_installed=True,
        factored_all_chart_mixed_derivative_boxes_available=True,
        mixed_coordinate_scope='y=log similarity R and paper axial Z; not Cartesian spatial derivatives',
        all_passed=True,full_pulse_C4_installed=False,full_outer_C4_certified=False,whole_outer_cone_certified=False,
        uniform_pulse_C4_chart_interface_certificate_available=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Actual mixed order<=4: independent physical/primitive derivatives, original-source agreement, chart overlaps and terminal closure PASS',flush=True)
    return out


if __name__=='__main__':run()
