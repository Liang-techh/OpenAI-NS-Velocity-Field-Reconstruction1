"""Inner exit stress coefficients through width one, with shared pressure.

Differentiating width-two velocity supplies width-one shear, not width-two
stress. No omitted-order or admissibility certificate is implied.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_pressure_width_jet import IntervalPressureWidthJet
from lei_ren_part1_paper_collar_interval_receipt_reader import read_inlet
from lei_ren_part1_paper_collar_first_width_interval_endpoint import first_width_coefficients
from lei_ren_part1_paper_collar_second_width_moments import second_width_spatial_coefficients
from lei_ren_part1_paper_collar_width_physical_endpoint import physical_endpoint
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_collar_low_width_spatial import primitive_nodes


def stress_coefficients(inlet, *, s, chi, J, K):
    first=first_width_coefficients(inlet,s=s,switch_integral=J)
    second=second_width_spatial_coefficients(inlet,s=s,switch_primitive=J,integrated_primitive=K)
    field=physical_endpoint(inlet,first,second,s=s)
    ctx=inlet['ctx'];template=inlet['F'].value
    W=IntervalPressureWidthJet(ctx,{(0,1):1},pressure_order=template.pressure_order,width_order=2)
    zero=template*0;s=zero+s;chi=zero+chi;J=zero+J
    r=field['R'].value;F=field['F'].value;u=field['Uz'].value
    z=inlet['z'].value;dt=zero+inlet['delta'];d=1-z*z;L=1-dt*z*z
    m={k:v.value for k,v in field['moments'].items()}
    mz={k:v.tangent for k,v in field['moments'].items()}
    root=(2*r).sqrt();V=-r+(1-dt)*z*m['z']+d*mz['z']
    It=F*V/L+((1-dt/2)*m['theta']-(1-dt)*z*mz['theta']/2
        -d*mz['theta_z']+(2*dt-1)*z*m['theta_z'])/(2*L*r)
    Iz=(V*u+(1-dt)*(m['z']-z*mz['z'])/2+2*dt*z*m['z_theta']
        -d*mz['z_theta']+r*(2*(1+dt)*z*field['P'].value-d*field['P'].tangent))/(L*root)
    Ra=inlet['R'].value;Fa=inlet['F'].value;D=inlet['D'].value;I=inlet['I_z'].value
    a=Ra*inlet['F_R'].value/Fa;b=(Ra/2).sqrt()
    gp1=-D*chi/2
    gp2=-(Ra*inlet['D_R'].value*s*chi+D*(1-chi))/2
    up1=-b*I*chi
    up2=-b*(Ra*inlet['I_z_R'].value*s*chi+I*((1/2-a)*s*chi-D*J*chi/2+1-chi))
    St=2*Fa*(gp1+W*(gp2+first['g'].value*gp1))
    Sz=(2/Ra).sqrt()*(up1+W*(up2-s*up1/2))
    values=dict(I_theta=It,I_z=Iz,S_theta=St,S_z=Sz,T_theta=It+St,T_z=Iz+Sz)
    def coefficient(ring,w):
        result=ring.component(0,w)*0
        for p in range(ring.pressure_order+1):result+=ring.component(p,w)
        return dict(nominal=result.nominal,pressure_error=result.difference)
    result={str(w):{k:coefficient(v,w) for k,v in values.items()} for w in (0,1)}
    # Separate algebraic identities catch radial shear factors and switch
    # derivatives; they do not test an omitted-order nonlinear remainder.
    expected={'T_theta':(1-chi)*inlet['I_theta'].value,
              'T_z':(1-chi)*I,'S_theta':-chi*inlet['I_theta'].value,'S_z':-chi*I}
    for name,target in expected.items():
        discrepancy=coefficient(values[name]-target,0)['nominal']
        lo,hi=endpoints(discrepancy)
        if not lo<=0<=hi:raise AssertionError('Leading stress identity: '+name)
    if endpoints(chi.evaluate().nominal)==(0,0):
        for name,target in (('S_theta',inlet['I_theta'].value),('S_z',I)):
            got=coefficient(values[name],1)['nominal']
            discrepancy=got+target.evaluate(width=0).nominal
            lo,hi=endpoints(discrepancy)
            if not lo<=0<=hi:raise AssertionError('Exit first-width shear identity: '+name)
    return result


def run():
    inlet,source=read_inlet();profile,_=accepted_profile();e=ScheduleEndpointEnclosures(profile.schedule)
    with mp.workdps(e.precision+40):
        nodes=primitive_nodes(e);rows=[]
        for text in ('0','.25','.5','.75','1','2'):
            s=e.iv.mpf(text);_,J,K=nodes[int(mp.mpf(text)*32)]
            chi=e.sigma_interval(1-s) if mp.mpf(text)<1 else e.iv.mpf(0)
            values=stress_coefficients(inlet,s=s,chi=chi,J=J,K=K)
            rows.append(dict(s=text,chi=chi,coefficients=values))
        report=dict(source=source,Z='.3',width_orders=[0,1],samples=rows,
            same_pressure_and_five_moments=True,second_width_stress_available=False,
            missing_for_second_width_stress='third-width velocity derivatives',
            leading_stress_identity_checks_passed=True,exit_first_width_shear_checks_passed=True,
            cone_certified=False,omitted_width_orders_enclosed=False,
            source_errors_enclosed=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        for row in rows:
            print('s',row['s'],'Tz0',mp.nstr(endpoints(row['coefficients']['0']['T_z']['nominal'])[0],12),flush=True)


if __name__=='__main__':run()
