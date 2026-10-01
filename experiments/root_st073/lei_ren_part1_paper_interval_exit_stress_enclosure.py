"""Directed stress and cone tests from actual exit fields and moments.

The angular amplitude is cancelled algebraically before normalizing angular
stress. The prescribed tiny shear is evaluated directly, never by subtraction
of two larger velocity derivatives. Inconclusive intervals are not failures.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet
from lei_ren_part1_paper_interval_comparison_enclosure import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def square(c,value):
    lo,hi = endpoints(value)
    left,right = c.mpf(lo)**2,c.mpf(hi)**2
    lower = c.mpf(0) if lo <= 0 <= hi else c.mpf(min(endpoints(left)[0],endpoints(right)[0]))
    upper = c.mpf(max(endpoints(left)[1],endpoints(right)[1]))
    return c.mpf([endpoints(lower)[0],endpoints(upper)[1]])


def cone(c,st,sz,tt,tz):
    prerequisites = endpoints(st)[1] < 0
    if not prerequisites:
        return dict(prerequisites_certified=False,relaxed_cone_certified=False,
            admissible_cone_certified=False,status='not certified: angular shear interval reaches zero')
    kappa = -(square(c,st)+square(c,sz))/st
    dot,cross = tt*st+tz*sz,-tt*sz+tz*st
    direction = endpoints(dot)[1] < 0
    if endpoints(kappa)[0] > 2:
        branch = 'kappa>2'
        margin = square(c,dot)*2-(kappa-2)*square(c,cross)
    elif endpoints(kappa)[1] <= 2:
        branch = 'kappa<=2'
        margin = -dot/(-st)-(2-kappa)
    else:
        return dict(prerequisites_certified=True,kappa=kappa,dot=dot,cross=cross,
            direction_certified=direction,branch='not isolated',
            relaxed_cone_certified=False,admissible_cone_certified=False,
            status='not certified: kappa interval crosses 2')
    passed = direction and endpoints(margin)[0] > 0
    ruled_out = endpoints(dot)[0] >= 0 or endpoints(margin)[1] <= 0
    return dict(prerequisites_certified=True,kappa=kappa,dot=dot,cross=cross,
        direction_certified=direction,branch=branch,margin=margin,
        relaxed_cone_certified=passed,admissible_cone_certified=passed and branch=='kappa>2',
        family_cone_ruled_out=ruled_out,
        status='certified' if passed else ('ruled out over enclosed family' if ruled_out else 'not certified: interval margin/direction unresolved'))


def evaluate(calc,packet):
    c = calc.ctx
    with mp.workdps(calc.precision+40):
        state = packet['normalized_exit_state']
        phi,u = state['phi'],state['U']
        R = packet['physical_R'] if 'physical_R' in packet else packet['target_physical_R']
        s = R*calc.lam
        from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_jet
        gy = restore_jet(c,packet['actual_terminal_g_y'],order=1)
        uy = restore_jet(c,packet['actual_terminal_U_y'],order=1)
        if endpoints(phi[0])[0] <= 0 or endpoints(calc.F0)[0] <= 0:
            raise ValueError('Positive angular field needed for normalized stress')
        inlet = normalized_inlet([phi],[u],state,calc.S.truncate(1),calc.ell.truncate(1),
            calc.p0.truncate(1),calc.lam,s,calc.z.truncate(1),calc.delta,
            endpoints_override=dict(phi_exit=phi,u_exit=u,phi_s=phi*gy/s,u_s=uy/s))
        F = calc.F0*phi[0]
        root = c.sqrt(2*R)
        Stheta = F*gy[0]*2
        Sz = uy[0]*root/R
        Itheta = calc.F0*inlet['itheta'][0]
        Iz = inlet['iz'][0]*c.sqrt(calc.eps)
        Ttheta = calc.F0*inlet['ttheta'][0]
        Tz = inlet['tz'][0]*c.sqrt(calc.eps)
        # S_theta/F is known exactly from the ODE normalization, without
        # interval division of two correlated copies of the tiny F0.
        st = gy[0]*2
        sz = Sz/F
        tt = inlet['ttheta'][0]/phi[0]
        tz = Tz/F
        return dict(center_family=calc.center_family,physical_R=R,
            actual_exit_fields_and_moments_used=True,
            angular_amplitude_cancelled_algebraically=True,
            tiny_shear_subtraction_used=False,
            stress=dict(I_theta=Itheta,I_z=Iz,S_theta=Stheta,S_z=Sz,T_theta=Ttheta,T_z=Tz),
            normalized_stress=dict(S_theta_over_F=st,S_z_over_F=sz,T_theta_over_F=tt,T_z_over_F=tz),
            cone=cone(c,st,sz,tt,tz),
            original_construction_parameter_errors_enclosed=False,
            whole_axis_or_whole_annulus_cone_certified=False,
            terminal_five_moment_closure=False,temporal_recursion=False)


def run():
    from lei_ren_part1_paper_interval_exit_continuation_enclosure import load_receipts
    from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
    calc = IntervalComparisonJets(4)
    bridge,cells = load_receipts(calc)
    comparison,initial_phi = cells['endpoint'],cells['initial_phi']
    from lei_ren_part1_paper_interval_exit_continuation_enclosure import continue_exit
    endpoint = continue_exit(calc,bridge,comparison,initial_phi,target_R='100')
    range_packet = dict(normalized_exit_state=endpoint['normalized_range_state'],
        physical_R=endpoint['physical_range']['physical_R'],
        actual_terminal_g_y=-(endpoint['comparison_range_driver']['D']*endpoint['exit_shear_epsilon']/2))
    range_packet['actual_terminal_U_y'] = -(endpoint['normalized_range_state']['phi']
        /endpoint['endpoint_comparison_state']['phi'].truncate(1)
        *endpoint['exit_shear_epsilon']*endpoint['comparison_range_driver']['I_z']
        *calc.ctx.sqrt(range_packet['physical_R']/2))
    result = dict(initial_exit=evaluate(calc,bridge),R100=evaluate(calc,endpoint),
                  full_exit_range=evaluate(calc,range_packet),state_sha256=calc.state_hash)
    result['input_hashes'] = {'lei_ren_part1_paper_interval_exit_continuation_enclosure.json':
        hashlib.sha256(Path(__file__).with_name('lei_ren_part1_paper_interval_exit_continuation_enclosure.json').read_bytes()).hexdigest()}
    here = Path(__file__).parent
    result['input_hashes'].update({name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
        ('lei_ren_part1_paper_interval_exit_stress_enclosure.py',
         'lei_ren_part1_paper_interval_exit_continuation_enclosure.py',
         'lei_ren_part1_paper_candidate_shared_inlet.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual exit cone enclosure:',result['initial_exit']['cone']['status'],
          '; R100:',result['R100']['cone']['status'],
          '; full range:',result['full_exit_range']['cone']['status'],flush=True)
    return result


if __name__ == '__main__':
    run()
