"""Focused exact radius ordering and original two-sided collar admission."""
import json
from pathlib import Path
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_loop_domain as source
from lei_ren_part1_paper_compliant_current_correlated_radius_operator import original_equalities,SYMBOLS


def radius_identities():
    lc,lp,lh,sc=s.symbols('logCstar logPstar loghb s_c',real=True)
    ra=s.log(4)-4*lp-1000
    rw=s.log(110)+10*lc+11*lp+1
    logminus=ra+s.exp(lh)*sc/2;logplus=rw+1
    a=SYMBOLS;equalities=original_equalities()
    replacements={a['log_epsilon']:-4*lp-1000,a['log_C']:lc,a['log_P']:lp}
    checks=dict(
        original_Ra_exact_source=s.simplify(equalities[a['log_Ra']].subs(replacements)-ra)==0,
        original_Rw_exact_source=s.simplify((equalities[a['log_Rref']]+a['log_P']+1).subs(replacements)-rw)==0,
        chosen_positive_log_radius_gap=s.simplify(logplus-logminus-(s.log(s.Rational(110,4))+10*lc+15*lp+1002-s.exp(lh)*sc/2))==0,
        chosen_left_to_R110_gap=s.simplify(s.log(110)-logminus-(s.log(s.Rational(110,4))+4*lp+1000-s.exp(lh)*sc/2))==0,
        right_to_Rc_gap=s.simplify((rw+2)-(rw+1)-1)==0,
        Rc_to_twice_Rc_gap=s.simplify((rw+2+s.log(2))-(rw+2)-s.log(2))==0)
    if not all(checks.values()):raise ArithmeticError('Original relative radius identities failed')
    return checks


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed original loop-domain dependency: '+name)
    field=source.CurrentLoopDomain();domain=field.domain();c=field.ctx;endpoints=source.packets.recovery.endpoints
    if source.packets.encode(domain)!=manifest['current_original_generic_loop_domain']:
        raise ValueError('Original collar/reservation domain differs')
    fraction=domain['left_collar_fraction_box'];center=domain['left_center_fraction']
    if not endpoints(fraction)[0]<endpoints(center)[0]<=endpoints(center)[1]<endpoints(fraction)[1]:
        raise ArithmeticError('Actual left collar must be two-sided at r_minus')
    right=domain['right_collar_original_power_offset_box'];middle=domain['right_center_original_power_offset']
    if not endpoints(right)[0]<endpoints(middle)[0]<=endpoints(middle)[1]<endpoints(right)[1]:
        raise ArithmeticError('Actual right collar must be two-sided at r_plus')
    for key in ('positive_log_r_plus_over_r_minus_lower','positive_log_R110_over_r_minus_lower',
        'left_collar_kappa_minus2_lower','right_collar_kappa_minus2'):
        if endpoints(domain[key])[0]<=0:raise ArithmeticError('Strict source geometry/collar margin lost')
    boundary=domain['both_boundary_kappa_excess_log_lower']
    for key in ('left_collar_kappa_minus2_lower','right_collar_kappa_minus2'):
        if endpoints(boundary-c.ln(domain[key]))[1]>0:
            raise ArithmeticError('Chosen boundary log lower exceeds an original source excess')
    if any(domain[k] for k in ('complete_modification_box_H0_minus2_certified',
        'original_outer_admission_twice_Rc_through_Rb_certified','complete_Section11_loop_domain_certified',
        'whole_generic_scales_instantiated')) or any(domain[k] for k in source.OPEN):
        raise ArithmeticError('Endpoint/reservation certificate exceeds scope')
    if not domain['formal_left_log_radius_offset']['not_added_to_huge_log_Ra']:
        raise ArithmeticError('Microscopic source offset lost by radius arithmetic')
    exact=radius_identities()
    query=domain['left_two_sided_original_source_query']
    if not query['strict_nonzero_stress_cone_certified_for_entire_query_box'] or query['exact_zero_core_inlet_case_only']:
        raise ArithmeticError('Strict collar confused with zero Ra endpoint')
    result=dict(all_passed=True,source_family=field.family,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        exact_relative_original_radius_identities=exact,two_sided_current_original_strict_collars=2,
        formal_positive_left_offset_retained=True,actual_positive_boundary_excesses_preserved_separately=True,
        whole_modification_margin_or_full_scales_promoted=False,ancestor_constructors_called=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),
            Path(__file__).name:source.sha(Path(__file__).name),
            source.PREFIX+'current_correlated_radius_operator.py':source.sha(source.PREFIX+'current_correlated_radius_operator.py')})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original loop domain PASS: six exact radius identities, two strict two-sided collars; whole margins/scales open',flush=True)
    return result


if __name__=='__main__':run()
