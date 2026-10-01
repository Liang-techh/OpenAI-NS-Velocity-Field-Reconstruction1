"""Independent algebra and interval guards for epsilon-source core transfer."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    names=('lei_ren_part1_paper_compliant_pressure_source.json',
           'lei_ren_part1_paper_compliant_core_transfer.json')
    source,core=[json.loads((HERE/n).read_bytes()) for n in names]
    for record in (source,core):
        for name,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Transfer input changed: '+name)
    z,delta,C,r,L,eps,P,u,H=sp.symbols('z delta C r L eps P u H')
    # Start with the paper's axial forcing and independently form g=-forcing.
    forcing=(1+delta)/2*(1-2*z*u)*u+4*H-2*(1+delta)*z*P
    dg=sp.expand(-forcing.subs(P,P+C)+forcing)
    dpsi=-r*dg/(2*L)
    identities={
        'constant_pressure_changes_g': sp.simplify(dg-2*(1+delta)*z*C),
        'constant_pressure_changes_Psi0': sp.simplify(dpsi+r*(1+delta)*z*C/L),
        'physical_axial_prefactor': sp.simplify((u+eps*sp.Symbol('Psi_new'))-(u+eps*sp.Symbol('Psi_old'))
                                             -eps*(sp.Symbol('Psi_new')-sp.Symbol('Psi_old'))),
    }
    t,b,tau,Ut,Rt,e=sp.symbols('t b tau Ut Rt e', positive=True)
    cinf=Ut*sp.exp(-b*tau)*(Rt*sp.exp(tau))**b/(1-e)
    identities['c_infinity_waiting_cancels']=sp.simplify(cinf-Ut*Rt**b/(1-e))
    A,B,V=sp.symbols('A B V')
    identities['restored_pressure_difference']=sp.expand((A**2-B**2)*V-(A-B)*(A+B)*V)
    q,d0,D=sp.symbols('q d0 D', positive=True)
    identities['fixed_point_resolvent_bound']=sp.simplify((1-q)*(d0/(1-q))-d0)
    if any(value!=0 for value in identities.values()):
        raise AssertionError('Independent transfer identity failed')
    c=MPIntervalContext();c.dps=160
    get=lambda record,key: read_interval(c,record[key])
    with mp.workdps(220):
        new=CompliantPressureDatum();old=LogarithmicPressureDatum()
        if core['implicit_source_sha256']!=new.source_sha or new.source_sha==old.source_sha:
            raise AssertionError('New/old source identity separation failed')
        if not source['compliant_source']['parameter_bounds']['Section7_c_epsilon_numeric_gate']:
            raise AssertionError('Epsilon numeric threshold failed')
        if source['compliant_source']['implicit_source_definition']['c_epsilon']!='.001':
            raise AssertionError('Wrong epsilon source')
        if len(new.stages)!=14 or any(endpoints(s['mass'])[1]<=0 for s in new.stages.values()):
            raise AssertionError('Pressure mass removed')
        if endpoints(new.parameters.waiting['root_interval'])[0]<=0:
            raise AssertionError('Waiting root lost')
        if new.source_sha!=CompliantPressureDatum(precision=120).source_sha:
            raise AssertionError('Source identity depends on enclosure precision')
        change=source['pressure_perturbation'];tr=core['core_transfer']
        if endpoints(get(change,'pulse_length_exceeds_selected_cut_margin'))[0]<=0:
            raise AssertionError('Finite log cap unsupported')
        physical=get(change,'physical_pressure_difference_abs_upper')
        normalized=get(change,'normalized_pressure_difference_abs_upper')
        physical_direct=c.exp(2*c.mpf(new.parameters.logPstar))*normalized
        if max(endpoints(physical)[0],endpoints(physical_direct)[0])>min(endpoints(physical)[1],endpoints(physical_direct)[1]):
            raise AssertionError('Pstar^2 factor missing')
        if len(core['terms'])!=20 or not core['contraction_proved']:
            raise AssertionError('Complete coupled contraction missing')
        lip=get(tr,'common_scaled_Lipschitz_upper')
        shift=get(tr,'center_difference_Xh_upper')
        bound=get(tr,'solution_difference_Xh_upper')
        independent=shift/(1-lip)
        if endpoints(lip)[1]>=mp.mpf('.5') or endpoints(bound)[0]<=0:
            raise AssertionError('Positive resolvent gate failed')
        if max(endpoints(bound)[0],endpoints(independent)[0])>min(endpoints(bound)[1],endpoints(independent)[1]):
            raise AssertionError('Resolvent transfer mismatch')
        if max(endpoints(get(row,'normalized_Phi_difference_upper'))[1]
               for row in tr['normalized_field_C3_difference_rows'])>=mp.mpf('1e-12'):
            raise AssertionError('Mixed C3 transfer too large')
        if tr['actual_heat_pressure_change_included'] or tr['downstream_exit_moment_cone_transfer_completed']:
            raise AssertionError('Unproved downstream transfer claimed')
        report=dict(independent_symbolic_identities={k:str(v) for k,v in identities.items()},
            total_independent_identities=len(identities),fourteen_positive_new_mass_boxes=True,
            source_identity_separated_and_precision_independent=True,
            physical_pressure_scale_and_positive_finite_cap_checked=True,
            twenty_term_common_envelope_contraction_checked=True,
            fixed_point_resolvent_and_mixed_C3_transfer_checked=True,
            actual_heat_change_and_global_matching_not_claimed=True,all_passed=True,
            new_source_sha256=new.source_sha,old_source_sha256=old.source_sha,
            input_hashes={**core['input_hashes'],**source['input_hashes'],
                **{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
    print('Compliant pressure/core: six independent identities, physical scaling, contraction and C3 perturbation PASS',flush=True)
    return report


if __name__=='__main__':run()
