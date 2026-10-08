"""Check the original parameter correlation and safe compact source constants."""
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame as source


def run():
    began=time.monotonic();frame=source.OriginalO2SourceParameterFrame();d=frame.definitions
    report=json.loads((source.HERE/source.NAME).read_bytes())
    assert report[source.GATE] and report['source_family']==frame.family
    for name,digest in report['input_hashes'].items():assert source.sha(name)==digest,name
    # Independent parameter identities, with exp represented by its exact
    # protected source operation instead of a hardware-sized scalar.
    assert d['Md']==40 and d['logPstar']==s.exp(40)+11 and d['yd']==d['logPstar']
    assert d['Td']+1==d['yd']
    assert d['log_mu']==-s.log(1000)-4*d['logPstar']
    assert d['log_delta']==-4*d['logPstar']-30
    assert d['log_epsilon']==d['log_delta']-s.log(1000)
    assert d['Tw']==-60*d['log_mu'] and d['L']==-30*d['log_mu']
    assert d['Ts']==4*(s.log(2)-d['log_delta']) and d['Tf']==100
    for name,log_name in (('mu','log_mu'),('delta','log_delta'),('epsilon','log_epsilon'),
        ('Pstar','logPstar'),('Rref','logRref')):
        assert d[name]==source.ExactSourceExponential(d[log_name])
    assert d['logRref']==s.log(110)+10*(d['logCstar']+d['logPstar'])
    assert frame.radius_log('.7')-frame.radius_log('.2')==s.Rational(1,2)
    sign,mantissa,exponent,bits=frame.selected_logCstar_mpf_tuple
    assert sign==0 and mantissa.bit_length()==bits
    assert d['logCstar'].args==(s.Integer(mantissa),s.Integer(exponent))
    assert d['logCstar'].evalf(20)==d['logCstar']
    assert source.ExactSourceExponential(s.Integer(2)).evalf(20)==source.ExactSourceExponential(s.Integer(2))
    assert frame.delta_branch_proof['four_logP_plus30_lower']>frame.delta_branch_proof['two_hundred_log10_strict_upper']
    assert s.Rational(1)+3+s.Rational(9,2)+s.Rational(27,6)>10
    assert 1+40+40**2//2==841
    p=frame.owner.pressure.partition;z=frame.owner.template['z']
    assert not d['W'].free_symbols and not d['W'].has(p['W'])
    assert all(s.diff(value,z)==0 for value in frame.substitutions.values())
    for order in (0,1,2):
        jet=frame.pressure_at(s.Rational(1,3),order)
        assert not jet.free_symbols
        assert not jet.has(*(p[key] for key in ('mu','delta','epsilon','yd','Tw','L','Ts','W')))
    got=frame.functions('.53')
    assert got['original_y_exact']==s.Rational(53,100)
    assert got['original_R']==source.ExactSourceExponential(d['logRref']+s.Rational(53,100))
    assert got['original_Pstar']==d['Pstar'] and got['original_delta']==d['delta']
    for key in ('p1','p2'):
        for row in got[key]:assert row.free_symbols<={z}
    assert not got['original_pressure_integrals_evaluated_numerically']
    assert not got['conditioned_native_phase_or_scalar_oracle_installed']
    rejected=0
    for function in (lambda:source.exact_coordinate(-1),lambda:source.exact_coordinate(2),
        lambda:source.ExactPositiveDyadic(0,5),lambda:source.ExactPositiveDyadic(1,s.Rational(1,2)),
        lambda:frame.pressure_at(2),lambda:source.exact_coordinate(frame.owner.profiles.ctx.nan),
        lambda:source.exact_coordinate(frame.owner.profiles.ctx.inf)):
        try:function()
        except ValueError:rejected+=1
    assert rejected==7
    flags=('original_p1_p2_scalar_point_values_installed','numerical_original_source_point_or_integral_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*source.inertial.profiles.loop.OPEN)
    assert all(report[key] is False for key in flags)
    result=dict(all_passed=True,**{source.GATE:True},source_family=frame.family,
        selected_physical_family_sha256=frame.selected_physical_family_sha256,
        source_definition_binding=frame.source_bindings,delta_branch_proof=frame.delta_branch_proof,
        independent_parameter_radius_and_length_identities=True,
        selected_dyadic_builder_parameter_used_not_pressure_or_radius_enclosure_endpoints=True,
        astronomical_constants_remain_protected_under_evalf=True,
        all_pressure_parameter_symbols_bound_at_three_Z_orders=True,
        full_O2_inertial_functions_have_only_Z_free=True,
        Z_independent_parameter_substitution_preserves_previous_checked_C1_identities=True,
        shared_exact_requested_radial_coordinate=True,invalid_source_domain_rejections=rejected,
        **dict.fromkeys(flags,False),input_hashes={**frame.hashes,
            source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Exact source parameter correlations and original selected radius definition only. Compact protected exponentials prevent hidden astronomical expansion. Pressure and native phase are not numerically evaluated; approximate radial coefficients and all global field/recursion limitations remain.')
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original O2 parameter frame: independent correlations, compact constants and full binding PASS',flush=True)
    return result


if __name__=='__main__':run()
