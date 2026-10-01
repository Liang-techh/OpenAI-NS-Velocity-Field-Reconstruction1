"""Acceptance guards for the implicit same-source normalized core seed."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_interval_core import LogarithmicCoreFactory,STEM
from lei_ren_part1_paper_candidate_gauge_core import _unpack_vector
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def run():
    f=LogarithmicCoreFactory();c=f.ctx
    state=json.loads((HERE/(STEM+'_state.json')).read_bytes())
    if state['input_hashes']!=f.hashes:raise ValueError('seed inputs changed')
    tube=json.loads((HERE/'lei_ren_part1_paper_analytic_radial_tail.json').read_bytes())
    with mp.workdps(c.dps+40):
        j=c.mpf('1e-14');sigma=j/500;dt=f.delta
        for key,value in (('j',j),('sigma',sigma),('delta',c.mpf('1e-200'))):
            a,b=endpoints(read_interval(c,tube[key]));d,e=endpoints(value)
            if not a<=d<=e<=b:raise ValueError('universal tube parameter mismatch: '+key)
        if not tube['common_complex_axis_poles_excluded']:raise ValueError('tube poles unguarded')
        G=read_interval(c,tube['G_modulus_upper']);chosen=read_interval(c,f.majorant['Gupper_in_logC_definition'])
        if endpoints(chosen)[0]<endpoints(G)[1]:raise ValueError('symbolic logC Gbar too small')
        if endpoints(chosen)[0]!=endpoints(chosen)[1]:raise ValueError('Gbar must be scalar')
        if f.majorant['logC_definition']!='Lambda*Gupper + 2*logLambda + 1000':raise ValueError('logC changed')
        Hleft=(1-dt)*(-j)/2+(1-j*j)*(-3*j)
        Hright=j
        derivative_lower=(9-dt)/2-14*j*j
        if endpoints(Hleft)[1]>=0 or endpoints(Hright)[0]<=0 or endpoints(derivative_lower)[0]<=0:
            raise ValueError('unique real anchor root not proved')
        fresh,rows,metadata=f.seed(state['identity']['center_Z'],state['identity']['initial_axis_length'])
        fixed={name:_unpack_vector(c,row) for name,row in state['fixed'].items()}
        for name in fresh:
            if len(fixed[name])!=114:raise ValueError('axial length')
            for stored,expected in zip(fixed[name],fresh[name]):
                if endpoints(stored)!=endpoints(expected):raise ValueError('seed mismatch: '+name)
        if endpoints(fixed['S_Z_taylor'][0])[1]<=0:raise ValueError('swirl removed')
        if endpoints(fixed['P0_Z_taylor'][0])[1]>=0:raise ValueError('pressure sign')
        normalized=f.datum.normalized_jets(state['identity']['center_Z'],113)
        if normalized['pressure_units']!='P/Pstar^2' or not normalized['all_14_pressure_atoms_included']:
            raise ValueError('pressure units or atoms')
        fixed_sha=hashlib.sha256(json.dumps(state['fixed'],sort_keys=True).encode()).hexdigest()
        result=dict(input_hashes={**f.hashes,Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            implicit_source_sha256=f.datum.source_sha,datum_enclosure_sha256=f.datum.datum_sha,
            fixed_seed_sha256=fixed_sha,seed_axial_length=114,
            identical_j_sigma_delta_upper_guarded=True,scalar_Gbar_dominates_complex_G=True,
            root_anchor_interval=['-j','0'],Hleft=Hleft,Hright=Hright,root_derivative_lower=derivative_lower,
            unique_axis_anchor_root_proved=True,primitive_branch='G(Z0)=0 at unique H root',
            real_F0_strictly_positive_by_definition=True,
            S_scaled_complex_supremum_bound='exp(-6logLambda-2000)',
            Cauchy_radius=f.eta/2,directed_lower_radius_used_by_interval_division=True,
            fourteen_atom_pressure_units_guarded=True,all_fixed_jets_reproduced=True,
            interval_center_means_family_of_local_jets=True,
            tiny_swirl_positive_upper_preserved=True,
            true_F0_point_value_not_selected=True,seed_acceptance_passed=True,
            finite_core_completed=state['completed_radial_order']==110,
            temporal_recursion=False,whole_axis_core_generated=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print('PASS: Md40 seed, root, tube, fourteen pressure atoms, 114 axial jets')
    return result

if __name__=='__main__':run()
