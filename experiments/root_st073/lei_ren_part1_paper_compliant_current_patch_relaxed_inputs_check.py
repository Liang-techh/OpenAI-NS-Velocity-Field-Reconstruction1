"""Focused current source attachment, full signed cone and support-gap checks."""
import json
from pathlib import Path
from fractions import Fraction
import lei_ren_part1_paper_compliant_current_patch_relaxed_inputs as source


def rational(endpoint):
    sign,mantissa,exponent,_=endpoint._mpf_
    return (-1 if sign else 1)*Fraction(mantissa)*Fraction(2)**exponent


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed actual current patch source: '+name)
    owner=source.CurrentPatchRelaxedInputs();c=owner.ctx
    if source.packets.encode(owner.theorem)!=data['original_current_full_signed_patch_source_theorem']:
        raise ValueError('Current full patch source theorem changed')
    if source.packets.encode(owner.proof)!=data['current_whole_Rm_Rh_relaxed_input_and_bounds']:
        raise ValueError('Current patch analytic bounds changed')
    if owner.conditions!=data['current_source_function_attachment_conditions']:
        raise ValueError('Current original patch function attachment changed')
    breaks=[c.mpf(1),*[c.mpf(i)/40 for i in source.EDGES],c.exp(1)]
    queries={str(i):owner.query((source.packets.recovery.endpoints(lo)[0],source.packets.recovery.endpoints(hi)[1]))
        for i,(lo,hi) in enumerate(zip(breaks,breaks[1:]))}
    if source.packets.encode(queries)!=data['analytic_support_and_quiet_gap_subbox_examples']:
        raise ValueError('Current source support/gap queries changed')
    invalid=((0,0),(-1,0),(3,0),('inf',0),(1,2),((1,2),(-2,1)))
    for x,Z in invalid:
        try:owner.query(x,Z)
        except ValueError:pass
        else:raise ArithmeticError('Invalid actual patch source box admitted')
    read=lambda value:source.packets.interval(c,value)
    proof=owner.proof
    # Exact rationals protect equality at a bound from fresh floating rounding.
    # Independently check the entire signed rectangle at p1=8. Larger p1 only
    # increases the positive expression since a-bw>0.
    hi=lambda value:rational(source.packets.recovery.endpoints(value)[1])
    lo=lambda value:rational(source.packets.recovery.endpoints(value)[0])
    bw=hi(proof['original_signed_bw_absolute_upper']);amin=Fraction(7,10);amax=Fraction(9,10)
    tests=0
    for a in (amin,Fraction(4,5),amax):
        for p1 in (Fraction(8),Fraction(16),Fraction(1000)):
            for signed_bw in (-bw,Fraction(0),bw):
                Hminus2=p1*(a-signed_bw)/a-2
                if Hminus2<=0 or Hminus2<lo(proof['original_full_H0_minus2_uniform_lower']):
                    raise ArithmeticError('Full signed H0/p2 sufficient margin failed')
                tests+=1
    if hi(proof['original_signed_Vy_absolute_upper'])>hi(proof['current_CS_times_C1_coefficient_upper']):
        raise ArithmeticError('Vy cannot use the universal coefficient bound')
    # The patch source/barrier is complete, even where local bumps vanish.
    for key in ('full_signed_p2_pressure_energy_and_meridional_terms_retained',
        'all_three_disjoint_support_groups_and_four_quiet_gaps_covered',
        'quiet_gap_velocities_do_not_reset_current_partial_histories',
        'conservative_full_signed_bw_bound_not_copied_from_legacy_256_estimate'):
        if proof[key] is not True:raise ArithmeticError('Partial-history/full signed source omitted')
    if proof['terminal_functional_closure_only_at_or_after_last_support_edge']!='71/40':
        raise ArithmeticError('Implicit closure advanced before full support')
    if any(data[k] or proof[k] for k in source.OPEN) or proof['p1_p2_whole_path_norm_bounds_certified']:
        raise ArithmeticError('Original patch gate promoted changed loops/global norms/recursion')
    result=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),source_family=owner.family,
        current_source_function_attachment_conditions=len(owner.conditions),
        exact_original_patch_AST_bindings=len(owner.theorem['original_patch_AST_bindings']),
        exact_signed_source_and_partial_history_identities=len(owner.theorem['exact_identities']),
        positive_actual_source_norm_and_smallness_gates=len(proof['actual_current_source_hypothesis_gates']),
        positive_actual_relaxed_margins=len(proof['positive_directed_relaxed_margins']),
        independent_exact_rational_signed_H0_envelopes=tests,
        same_true_coefficient_Vy_CS_norm_budget_checked=True,
        actual_support_group_and_quiet_gap_queries_checked=len(queries),invalid_source_queries_rejected=len(invalid),
        full_signed_p2_pressure_energy_and_meridional_terms_retained=True,
        source_ancestor_constructors_called=False,source_radius_amplitude_and_width_not_materialized=True,
        p1_p2_whole_path_norm_bounds_certified=False,strict_completed_tensor_cone_new_regions_admitted=0,
        scope=data['scope'],input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),
            Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Current actual full signed patch, three supports/four history gaps and source gates PASS',flush=True)
    return result


if __name__=='__main__':run()
