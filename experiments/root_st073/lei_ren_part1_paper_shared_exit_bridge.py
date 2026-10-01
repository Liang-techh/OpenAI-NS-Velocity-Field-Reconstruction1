"""Actual exit bridge with exactly shared h_b=epsilon_b.

Phase xi=y/h_b preserves the paper cutoffs even when the positive h value
is implicit. This is conditional on finite full paper K and input gates.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_shared_comparison import SharedComparison
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_candidate_shared_inlet import diff
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
HERE=Path(__file__).parent

def integrate(calc,comparison,core_cells=16):
    if comparison['parameter_family_sha256']!=calc.parameters.sha or comparison['identity']!=calc.identity:
        raise ValueError('shared parameter/source mismatch')
    if comparison['phase_domain']!=['1','2'] or not comparison['comparison_discretization_enclosed']:
        raise ValueError('phase coverage missing')
    if core_cells<1 or not isinstance(core_cells,int):raise ValueError('positive core cell count')
    c=calc.ctx
    with mp.workdps(c.dps+40):
        h=calc.h;step=h/core_cells;packets=[]
        for n in range(core_cells):
            left=c.mpf(n)/core_cells;right=c.mpf(n+1)/core_cells
            xi=c.mpf([endpoints(left)[0],endpoints(right)[1]])
            s=4*c.exp(h*xi);state=calc.core_state(s)
            packets.append(dict(phase_xi=xi,step=step,scaled_radius=s,state=state,
                transfer=calc.transfer(state,s),stage='core_comparison'))
        packets += [dict(row,stage='switched_comparison') for row in comparison['cell_enclosures']]
        start={n:v.truncate(1) for n,v in calc.core_state(c.mpf(4)).items()}
        phi_a,u_a=start['phi'],start['U'];g=phi_a*0;J=g
        moments={n:v for n,v in start.items() if n not in ('phi','U')};saved=[]
        for row in packets:
            xi,s,dy=row['phase_xi'],row['scaled_radius'],row['step']
            chi=h+(1-h)*alpha_box(c,xi+1,c.mpf(1))
            lo,hi=endpoints(chi);chi=c.mpf([max(mp.mpf(0),lo),min(mp.mpf(1),hi)])
            base=row['transfer'];A=base['unmodulated_driver_A'].truncate(1)*chi
            B=base['unmodulated_driver_B'].truncate(1)*chi
            partial=c.mpf([0,endpoints(dy)[1]])
            gr=g+A*partial;eg=gr.exp();phi=phi_a*eg;urhs=eg*B
            u=u_a+J+urhs*partial
            rhs=dict(theta=phi*(2*s*s),z=u*s,theta_z=phi*u*(2*s*s),
                p=phi*phi*s,u_squared=u*u*s,weighted_phi_squared=phi*phi*s*s)
            state=dict(phi=phi,U=u,**{n:v+rhs[n]*partial for n,v in moments.items()})
            saved.append(dict(phase_xi=xi,stage=row['stage'],actual_chi=chi,
                actual_driver_A=A,actual_driver_B=B,state=state))
            moments={n:v+rhs[n]*dy for n,v in moments.items()};g,J=g+A*dy,J+urhs*dy
        state=dict(phi=phi_a*g.exp(),U=u_a+J,**moments)
        if endpoints(state['phi'][0])[0]<=0:raise ValueError('shared actual swirl not positive')
        s=comparison['endpoint_scaled_radius'];z=calc.z.truncate(1);d=1-z*z;L=1-z*z*calc.delta
        pressure=calc.p0_scaled.truncate(1)+calc.S_scaled.truncate(1)*state['p']
        Ur=(z*state['U']*(2*s)-z*state['z']*(1-calc.delta)-d*diff(state['z']))*c.sqrt(calc.eps)/(L*c.sqrt(2*s))
        compname='lei_ren_part1_paper_shared_comparison.json'
        hashes={**calc.input_hashes,compname:hashlib.sha256((HERE/compname).read_bytes()).hexdigest(),
            Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        return dict(identity=calc.identity,input_hashes=hashes,parameter_family_sha256=calc.parameters.sha,
            parameters=calc.parameters.report(),phase_domain=['0','2'],closed_cell_count=len(packets),
            endpoint_scaled_radius=s,physical_R=s*calc.eps,
            normalized_actual_exit_state=state,log_F_over_Fa=g,Uz_increment=J,
            restored_pressure_scaled=pressure,physical_Ur=Ur,cell_enclosures=saved,
            actual_endpoint_driver_A=comparison['endpoint_transfer']['unmodulated_driver_A']*h,
            actual_endpoint_driver_B=comparison['endpoint_transfer']['unmodulated_driver_B']*h,
            shared_epsilon_box=h,shared_epsilon_upper=calc.parameters.h_upper,
            h_b_equals_exit_epsilon_by_definition=True,phase_cutoff_exact=True,
            actual_exit_ODE_integral_enclosed=True,actual_exit_not_frozen=True,
            positive_h_family_enclosed_without_zero_selection=True,retained_axial_order=1,
            full_Section9_parameter_admission=False,full_K_norm_certificate_missing=True,
            source_j_eta_tol_relation_verified=False,whole_axis=False,
            terminal_five_moment_repair_completed=False,temporal_recursion=False)
def run():
    calc=SharedComparison();comp=calc.integrate();r=integrate(calc,comp)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(r)),indent=2)+'\n',encoding='utf-8')
    print('Shared-parameter actual bridge: 48 phase cells; width and plateau identical')
    return r
if __name__=='__main__':run()
