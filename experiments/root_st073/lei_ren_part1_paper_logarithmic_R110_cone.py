"""Local-axis relaxed cone at the actual R110 zero-axial-source exit.

This checks one required interface, not the full transition/stress lift.
Positive physical F0 is implicit, allowing exact normalized cancellation.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_comparison import LogarithmicComparison,scaled_transfer
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_state,_pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def run():
    calc=LogarithmicComparison();c=calc.ctx;name='lei_ren_part1_paper_logarithmic_exit_switch.json'
    raw=json.loads((HERE/name).read_bytes())
    if raw['identity']!=calc.identity:raise ValueError('cone source changed')
    for n,d in raw['input_hashes'].items():
        if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('dependency changed: '+n)
    with mp.workdps(c.dps+40):
        m=restore_state(c,raw['normalized_actual_R110_state'],1)
        s=c.mpf(110)/calc.eps;z=calc.z.truncate(1)
        gy=-c.mpf('.4');Sz=c.mpf(0);St=2*gy
        transfer=scaled_transfer(c,m['phi'],m['U'],m,calc.ell_scaled.truncate(1),
            calc.S_scaled.truncate(1),calc.p0_scaled.truncate(1),calc.eps,s,z,calc.delta,calc.initial_phi.truncate(1))
        D=transfer['D'][0];Ttheta=D+St
        # kappa=.8<=2, Sz=0. Exact nonzero St cancellation reduces (3.23)
        # to Ttheta/F>0 and Ttheta/F > 2-kappa, i.e. D>2.
        kappa=-St;margin=Ttheta-(2-kappa)
        passed=endpoints(Ttheta)[0]>0 and endpoints(margin)[0]>0
        result=dict(identity=calc.identity,physical_R=110,axial_family=['.49','.51'],
            input_hashes={**calc.input_hashes,name:hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            S_theta_over_F=St,S_z_over_F=Sz,kappa=kappa,I_theta_over_F=D,T_theta_over_F=Ttheta,
            weak_branch_margin=margin,zero_axial_source_cancelled_exactly=True,
            positive_F0_and_Phi_cancelled_by_analytic_definition=True,
            relaxed_3_23_cone_certified_at_R110_local_axis=passed,
            admissible_stress_lift_constructed=False,whole_transition_cone_certified=False,
            paper_hb_cstar_K_inverse100_relation_verified=False,full_Section9_parameter_admission=False,
            whole_axis=False,terminal_five_moment_repaired=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Md40 actual R110 local-axis relaxed cone:',passed)
    return result

if __name__=='__main__':run()
