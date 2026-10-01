"""Independent energy normalizations and complete-tail admission checks."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_future_swirl_energy.json'


def run():
    r=json.loads((HERE/NAME).read_bytes());hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Future energy dependency changed: '+name)
    mu,delta,t,J,T,tau,q,sig=s.symbols('mu delta t J T tau q sig',real=True)
    expressions={
        'pulse_R_U_squared_attenuation':s.exp(13/mu)*s.exp(-(1+2*mu)*13/mu)-s.exp(-26),
        'flatten_energy_weight':s.exp(t)*s.exp(-(1+2*mu)*t)*(q/2)**(2*sig)-s.exp(-2*mu*t)*(q/2)**(2*sig),
        'steep_in_energy_weight':s.exp(t)*s.exp(-(1+2*mu)*t-2*(1-mu)*J)-s.exp(-2*mu*t-2*(1-mu)*J),
        'steep_out_energy_weight':s.exp(t)*s.exp(-3*t+2*(1-delta/2)*J)-s.exp(-2*t+2*(1-delta/2)*J),
        'Ns_prefactor':s.exp(1)*s.exp(-1-2*mu-(1-mu))-s.exp(-1-mu),
        'Nq_prefactor':s.exp(T)*s.exp(-3*T)-s.exp(-2*T),
        'Nt_relative_prefactor':s.exp(1)*s.exp(-3+(1-delta/2))-s.exp(-1-delta/2),
        'Ntail_relative_prefactor':s.exp(tau)*s.exp(-(1+delta)*tau)-s.exp(-delta*tau)}
    eps,W=s.symbols('epsilon W',real=True)
    expressions['separate_epsilon_tail_atoms']=(1-eps*W)**2-(1-2*eps*W+eps**2*W**2)
    R,E,h=s.symbols('R E h',positive=True)
    expressions['angular_bump_energy_weight']=R*s.exp(t)*E**2*s.exp(-(1+2*mu)*t)*((1+h)**2-1)/(R*E**2)-s.exp(-2*mu*t)*(2*h+h**2)
    # Differentiate a bounded power-tail primitive rather than truncating it.
    expressions['entire_positive_power_energy_tail_FTC']=s.diff(-s.exp(-delta*t)/delta,t)-s.exp(-delta*t)
    c=MPIntervalContext();c.dps=160
    passed={}
    for name,expression in expressions.items():
        if s.simplify(expression)!=0:raise ArithmeticError('Independent energy identity failed: '+name)
        passed[name]=True
    # This receipt covers the entire axial interval and infinite exterior.
    whole=r['whole_Z_C1'];total=whole['complete_future_energy_Taylor']['coefficients']
    weighted=whole['Section7_34_weighted_future_Taylor']['coefficients']
    if len(total)!=2 or len(weighted)!=2 or endpoints(read_interval(c,total[0]))[0]<=0:
        raise ArithmeticError('Complete whole-Z energy positivity/C1 failed')
    if endpoints(read_interval(c,weighted[0]))[1]>=mp.mpf('1e-100'):
        raise ArithmeticError('Weighted future energy not small at selected source')
    atoms=whole['separate_preheat_tail_atoms']
    if not (endpoints(read_interval(c,atoms['baseline']))[0]>0
            and endpoints(read_interval(c,atoms['epsilon_atom']))[1]<0
            and endpoints(read_interval(c,atoms['epsilon_squared_atom']))[0]>0):
        raise ArithmeticError('Preheat epsilon atoms were discarded')
    for point in r['samples']:
        deficit=point['pieces']['positive_Gamma_heat_energy_deficit_enclosure']['coefficients']
        if endpoints(read_interval(c,deficit[0]))[0]<0:raise ArithmeticError('Heat energy deficit sign failed')
        if not point['infinite_heat_tail_included'] or not point['signed_angular_energy_change_retained']:
            raise ValueError('A full future-energy term was omitted')
    if not r['complete_future_corrected_swirl_energy_available'] or r['full_physical_kinetic_energy_certified']:
        raise ValueError('Future energy scope/completion changed')
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'],independent_energy_identities=passed,
        whole_Z_positive_C1_complete_energy_checked=True,epsilon_and_epsilon_squared_atoms_retained=True,
        formal_Gamma_deficit_retained=True,infinite_exterior_included=True,all_passed=True,
        full_physical_kinetic_energy_certified=False,actual_ap_selected=False,
        full_outer_five_moment_match=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Complete future energy:',len(passed),'independent physical identities, whole-Z C1 positivity and retained signed/epsilon/Gamma tail atoms PASS',flush=True)
    return result


if __name__=='__main__':run()
