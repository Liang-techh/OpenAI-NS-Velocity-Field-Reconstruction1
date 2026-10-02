"""Independent physical radial C4 and actual fifth-moment chart checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_pulse_radial_C4.json'


def independent_physical_fixture():
    with mp.workdps(85):
        c=MPIntervalContext(); c.dps=75; radius=mp.mpf('1e-65')
        delta,mu=mp.mpf('.03'),mp.mpf('.04')
        ap=lambda z:mp.sqrt(1+z*z)+mp.mpf('.02')*z
        b=lambda z:mp.mpf('1.1')*ap(z); by=lambda z:mp.mpf('.6')*ap(z)
        m=lambda z:mp.mpf('.03')*(z+z**3)+mp.mpf('.3')*ap(z)
        theta=lambda z:1/(1+z*z)
        ur=lambda z:(2*z*theta(z)*b(z)-(1-delta)*z*theta(z)*m(z)
            -(1-z*z)*mp.diff(lambda x:theta(x)*m(x),z))/(1-delta*z*z)
        ury=lambda z:((1+delta)*z*theta(z)*b(z)-(1-z*z)*mp.diff(lambda x:theta(x)*b(x),z)
            +2*z*theta(z)*(by(z)-(mp.mpf('.5')+mu)*b(z)))/(1-delta*z*z)-ur(z)/2
        field=CompliantPulseRadialC4.__new__(CompliantPulseRadialC4)
        field.ctx=c; field.delta=c.mpf(delta); field.mu=c.mpf(mu); checks={}
        for Z in (mp.mpf(-1),mp.mpf(0),mp.mpf('.3'),mp.mpf(1)):
            jet=lambda f:IntervalTaylor(c,[c.mpf([v-radius,v+radius]) for v in mp.taylor(f,Z,5)])
            u=jet(theta); B=jet(b); M=jet(m); field.data=lambda z:(None,None,u,None,None)
            point=dict(Z=c.mpf(Z),Uz_over_Utheta=B,Uz_y_over_Utheta=jet(by)-B*(c.mpf('.5')+field.mu),
                Mz_over_R_Utheta=M,Mtheta_over_sqrt2_R_3half_Utheta=c.mpf(1),
                Ur_over_sqrt_R_over_2_Utheta=field.radial(c.mpf(Z),B,M))
            point=field._high_packet(point)
            cases=[('Ur',point['factored_Ur_over_sqrt_R_over_2_Pstar_Taylor'],ur),
                ('Ur_y',point['Ur_y_over_sqrt_R_over_2_Utheta']*u.truncate(4),ury),
                ('Ur_Z',point['Ur_Z_over_sqrt_R_over_2_Utheta']*u.truncate(3),lambda z:mp.diff(ur,z)),
                ('Ur_yZ',point['Ur_yZ_over_sqrt_R_over_2_Utheta']*u.truncate(3),lambda z:mp.diff(ury,z))]
            for label,actual,function in cases:
                for n,value in enumerate(mp.taylor(function,Z,actual.order)):
                    lo,hi=endpoints(actual[n])
                    if not lo<=value<=hi:raise ArithmeticError('Independent physical radial C4 failed: '+label+' '+str(n))
                    checks[label+'_Z'+str(Z)+'_order'+str(n)]=True
        return dict(checks=checks,actual_Md40_source_admission=False,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Radial C4 source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    old=json.loads((HERE/'lei_ren_part1_paper_compliant_pulse_high_jets.json').read_bytes())
    keys=('Uz_over_Utheta','Uz_y_over_Utheta','Mz_over_R_Utheta',
        'Mtheta_z_over_sqrt2_R_3half_Utheta_squared','Mtheta_over_sqrt2_R_3half_Utheta','Mztheta_over_R_Utheta_squared')
    radial_keys=('Ur_over_sqrt_R_over_2_Utheta','Ur_y_over_sqrt_R_over_2_Utheta')
    mixed_keys=('Ur_Z_over_sqrt_R_over_2_Utheta','Ur_yZ_over_sqrt_R_over_2_Utheta')
    groups=('main_samples','entrance_samples','gap_samples','gap_end_samples','end_samples')
    samples=[p for name in groups for p in r[name]]+[r[name] for name in
        ('whole_Z_terminal','whole_Z_main','whole_Z_end_support','whole_Z_support_crossing')]
    for point in samples:
        for key in keys:
            if len(point[key]['coefficients'])!=6:raise ValueError('Fifth pulse primitive lost: '+key)
        for key in radial_keys:
            if len(point[key]['coefficients'])!=5:raise ValueError('Physical radial axial C4 lost')
        for key in mixed_keys:
            if len(point[key]['coefficients'])!=4:raise ValueError('Mixed physical derivative C3 lost')
        p=point['pressure']
        for key in ('Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared','P_y_over_Pstar_squared'):
            if len(p[key]['coefficients'])!=6:raise ValueError('True pressure fifth-order data lost')
        for mv,p0v,pv in zip(p['Mp_over_Pstar_squared']['coefficients'],p['P0_over_Pstar_squared']['coefficients'],p['P_over_Pstar_squared']['coefficients']):
            expected=read(mv)+read(p0v); actual=read(pv)
            if not endpoints(actual)[0]<=endpoints(expected)[0]<=endpoints(expected)[1]<=endpoints(actual)[1]:
                raise ArithmeticError('Original analytic pressure not retained')
        if not point['radial_velocity_axial_C4_available'] or not point['divergence_preserved_by_Mz_recovery']:
            raise ValueError('Actual radial recovery omitted')
        for key in ('full_pulse_C4_installed','full_outer_C4_certified','whole_outer_cone_certified','temporal_recursion','fifth_derivative_Taylor_remainder_available'):
            if point[key]:raise ValueError('Axial C4 promoted to full field C4: '+key)
    comparisons=0
    for name in groups:
        for actual,prior in zip(r[name],old[name]):
            for key in keys+radial_keys+mixed_keys:
                for av,bv in zip(actual[key]['coefficients'],prior[key]['coefficients']):
                    a,b=endpoints(read(av)); d,e=endpoints(read(bv))
                    if max(a,d)>min(b,e):raise ArithmeticError('Earlier pulse source derivatives changed: '+key)
                    comparisons+=1
    whole=r['whole_Z_terminal']
    for key in ('Uz_over_Utheta','Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared')+radial_keys+mixed_keys:
        for v in whole[key]['coefficients']:
            if endpoints(read(v))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Terminal higher derivative closure lost: '+key)
    if endpoints(read(whole['Mztheta_over_R_Utheta_squared']['coefficients'][0]))[0]<=0:
        raise ArithmeticError('Terminal future energy not positive')
    for actual,expected in zip(whole['Mztheta_over_R_Utheta_squared']['coefficients'],whole['positive_terminal_future_energy_Taylor']['coefficients']):
        a,b=endpoints(read(actual)); d,e=endpoints(read(expected))
        if not a<=d<=e<=b:raise ArithmeticError('Fifth terminal selected energy identity lost')
    interfaces={}
    for label,left,right in (('main_gap_xi11',r['main_samples'][-1],r['gap_samples'][0]),
        ('gap_coordinate_overlap',r['gap_samples'][2],r['gap_end_samples'][0]),
        ('gap_end_s_minus4',r['gap_end_samples'][-1],r['end_samples'][0])):
        for key in keys+radial_keys+mixed_keys:
            for av,bv in zip(left[key]['coefficients'],right[key]['coefficients']):
                a,b=endpoints(read(av)); d,e=endpoints(read(bv))
                if max(a,d)>min(b,e):raise ArithmeticError('Radial C4 chart interface disjoint: '+label+' '+key)
        for key in ('Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared','P_y_over_Pstar_squared'):
            for av,bv in zip(left['pressure'][key]['coefficients'],right['pressure'][key]['coefficients']):
                a,b=endpoints(read(av)); d,e=endpoints(read(bv))
                if max(a,d)>min(b,e):raise ArithmeticError('Fifth pressure chart interface disjoint')
        interfaces[label]=True
    fixture=independent_physical_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        independent_physical_radial_C4_fixture=fixture,unchanged_lower_order_source_comparisons=comparisons,
        directed_fifth_primitive_radial_C4_interfaces=interfaces,whole_Z_terminal_fifth_derivatives_checked=True,
        pulse_primitives_axial_C5_installed=True,pulse_radial_velocity_axial_C4_installed=True,
        radial_velocity_Z_and_yZ_axial_C3_available=True,all_passed=True,full_pulse_C4_installed=False,
        full_outer_C4_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Actual radial axial C4: independent physical derivatives, all-chart lower-order agreement, fifth terminal and chart interfaces PASS',flush=True)
    return out


if __name__=='__main__':run()
