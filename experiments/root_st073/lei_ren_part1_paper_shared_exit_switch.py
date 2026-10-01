"""Same-source actual Section 9.4 switch from physical R100 to R110.

Full scaled pressure/swirl inputs, directed switch integrals and exact
terminal constant-power moments. Stress-cone certification is separate.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_shared_comparison import SharedComparison
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_state,_pack
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_exit_switch_enclosure import constant_power
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_shared_inlet import diff
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def run(cells_per_half=16):
    if not isinstance(cells_per_half,int) or cells_per_half<1:raise ValueError('positive cells required')
    calc=SharedComparison();c=calc.ctx
    names=('lei_ren_part1_paper_shared_comparison.json',
        'lei_ren_part1_paper_shared_exit_continuation.json')
    comp,continuation=[json.loads((HERE/n).read_bytes()) for n in names]
    for record in (comp,continuation):
        if record['identity']!=calc.identity or record['parameter_family_sha256']!=calc.parameters.sha:raise ValueError('switch source/parameter changed')
        for n,d in record['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('dependency changed: '+n)
    if continuation['target_physical_R']!='100':raise ValueError('switch requires physical R100')
    with mp.workdps(c.dps+40):
        state0=restore_state(c,continuation['target']['normalized_state'],1)
        phi0,u0=state0['phi'],state0['U'];g=phi0*0;J=g
        moments={n:v for n,v in state0.items() if n not in ('phi','U')}
        frozen=restore_state(c,comp['endpoint_state'],2);f,u=frozen['phi'],frozen['U']
        se=read_interval(c,comp['endpoint_scaled_radius'])
        gamma_upper=read_interval(c,continuation['shared_epsilon_upper'])
        gamma=c.mpf([0,endpoints(gamma_upper)[1]])
        hb=calc.h;step=hb/cells_per_half;cells=[]
        L=1-calc.z*calc.z*calc.delta
        for n in range(2*cells_per_half):
            left=c.mpf(n)/cells_per_half;right=c.mpf(n+1)/cells_per_half
            xi=c.mpf([endpoints(left)[0],endpoints(right)[1]])
            x=hb*xi
            s=c.exp(x)*100/calc.eps;ds=s-se;ds2=s*s-se*se
            increments=dict(theta=f*ds2,z=u*ds,theta_z=f*u*ds2,p=f*f*ds,
                u_squared=u*u*ds,weighted_phi_squared=f*f*ds2/2)
            bar=dict(phi=f,U=u,**{key:frozen[key]+v for key,v in increments.items()})
            base=calc.transfer(bar,s)
            if n<cells_per_half:
                blend=alpha_box(c,xi+1,c.mpf(1))
                A=base['unmodulated_driver_A']*gamma
            else:
                blend=c.mpf(0);sigma=1-alpha_box(c,xi,c.mpf(1))
                A=base['unmodulated_driver_A']*(gamma*(1-sigma))-sigma*c.mpf('.4')
            # Actual switch uses log(F/F100), so its axial driver uses Phi100.
            B=-(phi0/f)*base['scaled_axial_numerator']*gamma*blend/(2*L)
            partial=c.mpf([0,endpoints(step)[1]])
            gr=g+A*partial;eg=gr.exp();phi=phi0*eg
            U_rhs=eg*B;uv=u0+J+U_rhs*partial
            rhs=dict(theta=phi*(2*s*s),z=uv*s,theta_z=phi*uv*(2*s*s),
                p=phi*phi*s,u_squared=uv*uv*s,weighted_phi_squared=phi*phi*s*s)
            box=dict(phi=phi,U=uv,**{key:v+rhs[key]*partial for key,v in moments.items()})
            cells.append(dict(phase_xi=xi,x_interval=x,stage=1 if n<cells_per_half else 2,
                actual_A=A,actual_B=B,normalized_state=box))
            moments={key:v+rhs[key]*step for key,v in moments.items()}
            g,J=g+A*step,J+U_rhs*step
        s2=100*c.exp(2*hb)/calc.eps
        state2=dict(phi=phi0*g.exp(),U=u0+J,**moments)
        length=c.ln(c.mpf('1.1'))-2*hb
        if endpoints(length)[0]<=0:raise ValueError('constant-power interval absent')
        target=constant_power(c,state2,s2,length);st=c.mpf(110)/calc.eps
        if endpoints(target['phi'][0])[0]<=0:raise ValueError('R110 normalized swirl not positive')
        pressure=calc.p0_scaled.truncate(1)+calc.S_scaled.truncate(1)*target['p']
        z=calc.z.truncate(1);d=1-z*z;L=1-z*z*calc.delta;mz=target['z']
        Ur=(z*target['U']*(2*st)-z*mz*(1-calc.delta)-d*diff(mz))*c.sqrt(calc.eps)/(L*c.sqrt(2*st))
        hashes={**calc.input_hashes,**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
            Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'lei_ren_part1_paper_interval_exit_switch_enclosure.py':hashlib.sha256((HERE/'lei_ren_part1_paper_interval_exit_switch_enclosure.py').read_bytes()).hexdigest()}
        result=dict(identity=calc.identity,input_hashes=hashes,
            implicit_source_sha256=calc.identity['implicit_source_sha256'],
            datum_enclosure_sha256=calc.identity['datum_enclosure_sha256'],
            physical_R=110,scaled_radius=st,cells_per_half=cells_per_half,
            switch_cells=cells,switch_endpoint_state=state2,constant_power_log_length=length,
            normalized_actual_R110_state=target,restored_pressure_scaled=pressure,physical_Ur=Ur,
            terminal_log_swirl_slope=target['phi']*0-c.mpf('.4'),terminal_axial_slope=target['U']*0,
            controlled_actual_switch_integrals=True,constant_power_moments_exact=True,
            same_source_core_comparison_actual_exit_chain=True,
            retained_axial_order=1,physical_F0_not_materialized=True,
            reference_inlet_five_moment_identities_repaired=False,
            switch_stress_cone_certified=False,admissible_stress_lift_constructed=False,
            parameter_family_sha256=calc.parameters.sha,parameters=calc.parameters.report(),
            h_b_equals_exit_epsilon_by_definition=True,full_Section9_parameter_admission=False,
            whole_axis=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Shared-parameter controlled actual switch to R110: 32 cells plus exact constant-power moments')
    return result

if __name__=='__main__':run()
