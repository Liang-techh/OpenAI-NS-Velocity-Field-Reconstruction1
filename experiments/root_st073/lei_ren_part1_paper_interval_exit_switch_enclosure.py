"""Controlled actual Section 9.4 source switches R=100..110.

Two cell-range integrations followed by exact constant-power moments.
No terminal moment repair, whole-axis cone, or temporal recursion claim.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box, pack
from lei_ren_part1_paper_interval_exit_continuation_enclosure import (
    load_receipts, continue_exit, _frozen_state, _driver, _physical_packet)
from lei_ren_part1_paper_interval_exit_stress_enclosure import evaluate
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def zero_axial_source_cone(calc, stress, a):
    # With Sz=0 and st=-a, kappa=a and -dot/(-st)=tt exactly.
    # Cancel the common nonzero st before interval division.
    c=calc.ctx
    tt=stress['normalized_stress']['T_theta_over_F']
    if endpoints(stress['normalized_stress']['S_z_over_F']) != endpoints(c.mpf(0)):
        raise ValueError('zero axial shear required for cancellation')
    av=a[0] if hasattr(a,'coefficients') else a
    lo,hi=endpoints(av)
    if not 0<lo or hi>2:
        return stress['cone']
    margin=tt-(2-av)
    passed=endpoints(tt)[0]>0 and endpoints(margin)[0]>0
    return dict(prerequisites_certified=True,kappa=av,branch='kappa<=2',
        direction_certified=endpoints(tt)[0]>0,margin=margin,
        relaxed_cone_certified=passed,admissible_cone_certified=False,
        family_cone_ruled_out=endpoints(tt)[1]<=0 or endpoints(margin)[1]<=0,
        zero_axial_source_cancelled_algebraically=True,
        status='certified' if passed else 'not certified: zero-axial-source margin unresolved')


def constant_power(c, state, s2, length):
    """phi falls as exp(-2x/5), U stays fixed, exact six moments."""
    f,u = state['phi'],state['U']
    k = c.mpf('0.4')
    phi = f*c.exp(-k*length)
    inc1 = s2*(c.exp(length)-1)
    angular = s2*s2*(c.exp((2-k)*length)-1)*2/(2-k)
    quadratic = s2*(c.exp((1-2*k)*length)-1)/(1-2*k)
    weighted = s2*s2*(c.exp((2-2*k)*length)-1)/(2-2*k)
    increments = dict(theta=f*angular,z=u*inc1,theta_z=f*u*angular,
        p=f*f*quadratic,u_squared=u*u*inc1,weighted_phi_squared=f*f*weighted)
    return dict(phi=phi,U=u,**{key:state[key]+value for key,value in increments.items()})


def integrate_switch(calc, continuation, cells_per_half=16):
    if not isinstance(cells_per_half,int) or cells_per_half<1:
        raise ValueError('positive integer cell count required')
    c = calc.ctx
    with mp.workdps(calc.precision+40):
        if endpoints(continuation['target_physical_R']) != endpoints(c.mpf(100)):
            raise ValueError('switch requires R100 actual endpoint')
        state0 = continuation['normalized_exit_state']
        phi0,u0 = state0['phi'],state0['U']
        bar = continuation['endpoint_comparison_state']
        se = continuation['endpoint_scaled_radius']
        eps_exit = continuation['exit_shear_epsilon']
        hb = c.mpf('.005')
        zero = phi0*0
        g,j = zero,zero
        moments = {key:value for key,value in state0.items() if key not in ('phi','U')}
        stress_cells=[]
        def comparison(s):
            ds=c.mpf([max(mp.mpf(0),endpoints(s-se)[0]),endpoints(s-se)[1]])
            ds2=c.mpf([max(mp.mpf(0),endpoints(s*s-se*se)[0]),endpoints(s*s-se*se)[1]])
            return _frozen_state(bar,bar['phi'],bar['U'],ds,ds2)
        step=hb/cells_per_half
        for index in range(2*cells_per_half):
            left=step*index
            x=c.mpf([endpoints(left)[0],endpoints(left+step)[1]])
            s=c.exp(x)*100/calc.eps
            driver=_driver(calc,bar['phi'],bar['U'],comparison(s),s)
            # alpha_box is decreasing on its second half. Shifting x by hb
            # gives sigma(x/hb)=1-alpha_box(x+hb,hb).
            if index<cells_per_half:
                blend=alpha_box(c,x+hb,hb)
                a=driver['D']*eps_exit
            else:
                blend=c.mpf(0)
                sigma=1-alpha_box(c,x,hb)
                a=driver['D']*(eps_exit*(1-sigma))+(phi0*0+c.mpf('.8'))*sigma
            A=-a/2
            B=-(phi0/bar['phi'].truncate(1))*driver['I_z']*eps_exit*c.sqrt(s*calc.eps/2)*blend
            partial=c.mpf([0,endpoints(step)[1]])
            gr=g+A*partial
            phi=phi0*gr.exp()
            urhs=gr.exp()*B
            u=u0+j+urhs*partial
            rhs=dict(theta=phi*(2*s*s),z=u*s,theta_z=phi*u*(2*s*s),
                p=phi*phi*s,u_squared=u*u*s,weighted_phi_squared=phi*phi*s*s)
            range_m={key:value+rhs[key]*partial for key,value in moments.items()}
            packet=dict(normalized_exit_state=dict(phi=phi,U=u,**range_m),
                physical_R=s*calc.eps,actual_terminal_g_y=A,actual_terminal_U_y=urhs)
            stress=evaluate(calc,packet)
            if index>=cells_per_half:
                stress['cone']=zero_axial_source_cone(calc,stress,a)
            stress_cells.append(dict(index=index,stage=1 if index<cells_per_half else 2,
                x=x,cone=stress['cone']))
            moments={key:value+rhs[key]*step for key,value in moments.items()}
            g,j=g+A*step,j+urhs*step
        s2=c.exp(2*hb)*100/calc.eps
        state2=dict(phi=phi0*g.exp(),U=u0+j,**moments)
        length=c.ln(c.mpf('1.1'))-2*hb
        if endpoints(length)[0]<=0:
            raise ValueError('R110 must follow both switches')
        target=constant_power(c,state2,s2,length)
        sT=c.mpf(110)/calc.eps
        physical=_physical_packet(calc,target,sT)
        gy=target['phi']*0-c.mpf('.4')
        uy=target['U']*0
        target_packet=dict(normalized_exit_state=target,physical_R=c.mpf(110),
            actual_terminal_g_y=gy,actual_terminal_U_y=uy)
        terminal_stress=evaluate(calc,target_packet)
        terminal_stress['cone']=zero_axial_source_cone(calc,terminal_stress,c.mpf('.8'))
        range_length=c.mpf([0,endpoints(length)[1]])
        constant_state=constant_power(c,state2,s2,range_length)
        constant_packet=dict(normalized_exit_state=constant_state,
            physical_R=s2*calc.eps*c.exp(range_length),actual_terminal_g_y=gy,actual_terminal_U_y=uy)
        constant_stress=evaluate(calc,constant_packet)
        constant_stress['cone']=zero_axial_source_cone(calc,constant_stress,c.mpf('.8'))
        return dict(center_family=calc.center_family,state_sha256=calc.state_hash,
            cells_per_half=cells_per_half,hb=hb,switch_endpoint_R=s2*calc.eps,
            switch_endpoint_state=state2,normalized_exit_state=target,
            actual_terminal_g_y=gy,actual_terminal_U_y=uy,**physical,
            switch_cell_cones=stress_cells,R110_stress=terminal_stress,
            constant_power_range_stress=constant_stress,
            retained_axial_order=1,controlled_switch_cell_integrals=True,
            final_constant_power_moments_exact=True,midpoint_projection_used=False,
            original_construction_parameter_errors_enclosed=False,
            whole_axis_transition=False,whole_switch_cone_certified=False,
            terminal_five_moment_closure=False,Ur_Z_available=False,temporal_recursion=False)


def run(cells_per_half=16):
    calc=IntervalComparisonJets(4)
    bridge,cells=load_receipts(calc)
    continuation=continue_exit(calc,bridge,cells['endpoint'],cells['initial_phi'])
    result=integrate_switch(calc,continuation,cells_per_half)
    here=Path(__file__).parent
    names=('lei_ren_part1_paper_interval_exit_switch_enclosure.py',
        'lei_ren_part1_paper_interval_exit_continuation_enclosure.py',
        'lei_ren_part1_paper_interval_exit_continuation_enclosure.json',
        'lei_ren_part1_paper_interval_exit_stress_enclosure.py',
        'lei_ren_part1_paper_interval_comparison_enclosure.py')
    result['input_hashes']={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    counts={}
    for row in result['switch_cell_cones']:
        status=row['cone']['status'];counts[status]=counts.get(status,0)+1
    print('Actual controlled switch to R110 complete;',counts,'R110:',result['R110_stress']['cone']['status'],flush=True)
    return result

if __name__=='__main__':
    run()
