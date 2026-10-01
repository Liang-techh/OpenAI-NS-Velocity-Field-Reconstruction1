"""Actual initial exit ODE, with implicit positive microscopic shear.

The comparison drivers are cutoff-modulated; this evolves the actual
exit field rather than freezing the comparison. C1 axial family only.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_comparison import LogarithmicComparison
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_candidate_shared_inlet import diff

HERE=Path(__file__).parent

def integrate(calc,comparison,core_cells=16):
    if not isinstance(core_cells,int) or core_cells<1:raise ValueError('positive core cell count')
    if comparison['identity']!=calc.identity or not comparison['cell_integral_discretization_enclosed'] or not comparison['core_analytic_tail_propagated']:
        raise ValueError('comparison admission missing')
    if comparison['y_interval']!=['.005','.01'] or comparison['core_branch_y_interval']!=['0','.005']:
        raise ValueError('comparison coverage changed')
    c=calc.ctx
    with mp.workdps(c.dps+40):
        major=json.loads((HERE/'lei_ren_part1_paper_logarithmic_core_majorant.json').read_bytes())
        loglam=read_interval(c,major['logLambda']);lam=read_interval(c,major['Lambda'])
        if endpoints(lam)[0]<endpoints(4*loglam+1000)[1]:
            raise ValueError('microscopic shear upper envelope unproved')
        j=c.mpf('1e-14');sigma=j/500
        realG=(1+j)*(1+c.mpf('1e-200'))/(2*sigma)
        Grealbar=c.mpf(endpoints(realG)[1])
        shear_upper=c.exp(-4*loglam-1000)
        shear=c.mpf([0,endpoints(shear_upper)[1]])
        # Actual gamma=exp(-Lambda*(2Grealbar+1))>0. gamma<=exp(-Lambda)
        # <=exp(-4logLambda-1000), with the above explicit guard.
        hb=c.mpf('.005');step=hb/core_cells;packets=[]
        for n in range(core_cells):
            lo=step*n;hi=step*(n+1)
            y=c.mpf([endpoints(lo)[0],endpoints(hi)[1]]);s=4*c.exp(y)
            state=calc.core_state(s)
            packets.append(dict(y_interval=y,step=step,scaled_radius=s,state=state,
                transfer=calc.transfer(state,s),phase='core_comparison'))
        for row in comparison['cell_enclosures']:
            packets.append(dict(row,step=hb/comparison['cells'],phase='switched_comparison'))
        start={n:v.truncate(1) for n,v in calc.core_state(c.mpf(4)).items()}
        phi_a,u_a=start['phi'],start['U'];g=phi_a*0;J=g
        moments={n:v for n,v in start.items() if n not in ('phi','U')}
        saved=[]
        for n,packet in enumerate(packets):
            y,s,dy=packet['y_interval'],packet['scaled_radius'],packet['step']
            base=packet['transfer']
            chi=shear+(1-shear)*alpha_box(c,y+hb,hb)
            # Convex combination: exact mathematical chi lies in [0,1].
            a,b=endpoints(chi);chi=c.mpf([max(mp.mpf(0),a),min(mp.mpf(1),b)])
            A=base['unmodulated_driver_A'].truncate(1)*chi
            B=base['unmodulated_driver_B'].truncate(1)*chi
            partial=c.mpf([0,endpoints(dy)[1]])
            grange=g+A*partial;eg=grange.exp();urhs=eg*B
            u=u_a+J+urhs*partial;phi=phi_a*eg
            rhs=dict(theta=phi*(2*s*s),z=u*s,theta_z=phi*u*(2*s*s),
                p=phi*phi*s,u_squared=u*u*s,weighted_phi_squared=phi*phi*s*s)
            box=dict(phi=phi,U=u,**{key:value+rhs[key]*partial for key,value in moments.items()})
            saved.append(dict(y_interval=y,phase=packet['phase'],cutoff_chi=chi,
                actual_driver_A=A,actual_driver_B=B,actual_state=box))
            moments={key:value+rhs[key]*dy for key,value in moments.items()}
            g,J=g+A*dy,J+urhs*dy
        state=dict(phi=phi_a*g.exp(),U=u_a+J,**moments)
        if endpoints(state['phi'][0])[0]<=0:raise ValueError('normalized actual swirl positivity lost')
        s=comparison['endpoint_scaled_radius'];z=calc.z.truncate(1)
        pressure=calc.p0_scaled.truncate(1)+calc.S_scaled.truncate(1)*state['p']
        mz=state['z'];L=1-z*z*calc.delta;d=1-z*z
        Ur=(z*state['U']*(2*s)-z*mz*(1-calc.delta)-d*diff(mz))*c.sqrt(calc.eps)/(L*c.sqrt(2*s))
        Aend=comparison['endpoint_transfer']['unmodulated_driver_A']*shear
        Bend=comparison['endpoint_transfer']['unmodulated_driver_B']*shear
        return dict(identity=calc.identity,implicit_source_sha256=calc.identity['implicit_source_sha256'],
            datum_enclosure_sha256=calc.identity['datum_enclosure_sha256'],
            input_hashes={**calc.input_hashes,'lei_ren_part1_paper_logarithmic_comparison.json':hashlib.sha256((HERE/'lei_ren_part1_paper_logarithmic_comparison.json').read_bytes()).hexdigest(),Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            core_cells=core_cells,switch_cells=comparison['cells'],closed_cell_count=len(packets),
            y_interval=['0','.01'],endpoint_scaled_radius=s,physical_R=s*calc.eps,
            shear_definition='gamma=exp(-Lambda*(2Grealbar+1))',Grealbar=Grealbar,
            actual_shear_strictly_positive_by_definition=True,
            shear_parameter_policy='candidate conservative uniform real-G bound; not paper-derived cstar K^-100',
            comparison_hb_candidate='.005',paper_hb_cstar_K_inverse100_relation_verified=False,
            full_Section9_parameter_admission=False,
            shear_upper_envelope=shear_upper,shear_zero_lower_is_bound_not_selected_value=True,
            log_F_over_Fa=g,Uz_increment=J,normalized_actual_exit_state=state,
            restored_pressure_scaled=pressure,physical_Ur=Ur,
            actual_endpoint_driver_A=Aend,actual_endpoint_driver_B=Bend,
            actual_endpoint_phi_y=state['phi']*Aend,actual_endpoint_Uz_y=g.exp()*Bend,
            cell_enclosures=saved,retained_axial_order=1,
            actual_exit_ODE_integral_enclosed=True,frozen_comparison_not_used_as_actual_field=True,
            physical_F0_not_materialized=True,actual_F0_times_phi_strictly_positive_by_definition=True,
            terminal_reference_inlet_constructed=False,terminal_five_moment_repair_completed=False,
            whole_axis=False,admissible_stress_lift_constructed=False,temporal_recursion=False)

def run():
    calc=LogarithmicComparison();comparison=calc.integrate();result=integrate(calc,comparison)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Md40 actual exit bridge:',result['closed_cell_count'],'closed cells; normalized positive swirl; C1 axial family')
    return result

if __name__=='__main__':run()
