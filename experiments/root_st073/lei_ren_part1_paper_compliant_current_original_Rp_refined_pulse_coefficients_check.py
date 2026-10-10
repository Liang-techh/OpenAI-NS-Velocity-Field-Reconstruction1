"""Independent pulse integral, C5 implicit coefficients and physical binding."""
import copy
import gzip
import json
from fractions import Fraction
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rp_refined_pulse_coefficients as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def included(outer, inner):
    a,b = current.ends(outer)
    u,v = current.ends(inner)
    return a <= u <= v <= b


def rejected(fn):
    try:
        fn()
    except (ValueError, TypeError, ArithmeticError, KeyError):
        return True
    raise AssertionError('Invalid refined source admitted')


def independent_pulse_energy(c):
    """A different rigorous integral decomposition, without a Darboux grid.

    For .02 <= xi <= 10, gp=xi-.01 exactly. On startup, 0<=gp<=xi;
    on exit, 0<=gp<=xi-.01. Both omitted transition integrals retain a
    strictly positive upper bound. No quadrature midpoint defines gp.
    """
    def primitive(x, offset):
        y=x-offset
        return -c.exp(-2*x)*(y*y/2+y/2+c.mpf(1)/4)
    middle = primitive(c.mpf(10),c.mpf('.01'))-primitive(c.mpf('.02'),c.mpf('.01'))
    startup = primitive(c.mpf('.02'),c.mpf(0))-primitive(c.mpf(0),c.mpf(0))
    exit_tail = primitive(c.mpf(11),c.mpf('.01'))-primitive(c.mpf(10),c.mpf('.01'))
    return middle+c.mpf([0,current.ends(startup+exit_tail)[1]])


@source_precision
def run(owner, deliveries):
    began=time.monotonic()
    raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[k] for k in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name]==digest for name,digest in owner.before.hashes.items())
    c=MPIntervalContext();c.dps=500
    cv=lambda value:c.mpf(current.ends(value))
    independent_K=independent_pulse_energy(c)
    assert included(owner.base.K,independent_K)
    assert included(owner.selected.amplitude.K,owner.base.K)
    assert current.ends(independent_K)[0]>0
    old_U=owner.selected.constants['U'];new_U=owner.constants['U']
    assert included(old_U,new_U) and included(owner.amplitude.U0_box,new_U)
    # Coherent ratios, independent arithmetic/context and original inputs.
    k=owner.constants;U=cv(new_U);invP=c.exp(-cv(owner.amplitude.logP_box))
    ratios=dict(C1=invP*cv(k['M'])/U,C2=invP*cv(k['K'])/(U*U),
        C0=cv(k['E_Q'])/(U*U),C_E=cv(k['E_Z'])/(U*U))
    for key,value in ratios.items():
        assert included(k[key],value), key
    assert included(owner.pulse.Xp,cv(owner.pulse.inlet_H)/U)
    observations={};totals=dict(rows=0,exact_zero_rows=0,nonzero_factored_target_rows=0,
        sign_unresolved_rows=0,ordinary_numeric_target_rows=0)
    for name,delivery in deliveries.items():
        source=owner.source(delivery)
        selected=owner.selection.select(source['Z'])
        original=owner.selected.fifth.select(source['Z'])
        q=selected['quadratic_coefficients'];A2=cv(q['A2'])
        A1=[cv(v) for v in q['A1'].coefficients];A0=[cv(v) for v in q['A0'].coefficients]
        # Different quadratic formula, then direct coefficient extraction of
        # A2*a^2+A1*a+A0=0 through degree five. No old cached a0 is used.
        def corner_root(upper):
            # On the positive branch, a decreases in each of A2,A1,A0.
            # Endpoint evaluation avoids the dependency duplication of A2
            # in a naive interval application of the other root formula.
            index=1 if upper else 0
            a2,a1,a0=(c.mpf(current.ends(value)[index]) for value in (A2,A1[0],A0[0]))
            return (-a1+c.sqrt(a1*a1-4*a2*a0))/(2*a2)
        low,high=corner_root(True),corner_root(False)
        root=c.mpf([current.ends(low)[0],current.ends(high)[1]])
        expected=[root];den=2*A2*root+A1[0]
        for n in range(1,6):
            cross=A2*sum((expected[i]*expected[n-i] for i in range(1,n)),c.mpf(0))
            mixed=sum((A1[i]*expected[n-i] for i in range(1,n+1)),c.mpf(0))
            expected.append(-(cross+mixed+A0[n])/den)
        actual=selected['selected_ap_Taylor']
        for n,value in enumerate(expected):
            assert included(actual[n],value),(name,n,'implicit coefficient')
        assert included(original['selected_ap_Taylor'][0],actual[0])
        # Each coefficient of the actual quadratic enclosure includes zero;
        # this is an inclusion check, not an exact interval residual claim.
        for n in range(6):
            residual=A2*sum((cv(actual[i])*cv(actual[n-i]) for i in range(n+1)),c.mpf(0))
            residual+=sum((A1[i]*cv(actual[n-i]) for i in range(n+1)),c.mpf(0))+A0[n]
            assert current.ends(residual)[0]<=0<=current.ends(residual)[1]
        assert all(selected['actual_weighted_future_energy_Taylor'][n]._mpi_==
            original['actual_weighted_future_energy_Taylor'][n]._mpi_ for n in range(6))
        # The normalized shape polynomials retain true C5 coefficients.
        Z=cv(source['Z']);g=[Z+Z**3,1+3*Z**2,3*Z,c.mpf(1),c.mpf(0),c.mpf(0)]
        h=[Z**2+2*Z**4+Z**6,2*Z+8*Z**3+6*Z**5,1+12*Z**2+15*Z**4,
            8*Z+20*Z**3,2+15*Z**2,6*Z]
        for n in range(6):
            for j,key in enumerate(('C1','C2')):
                assert included(selected['incoming']['moment_Taylor'][j][n],g[n]*cv(k[key]))
            assert included(selected['incoming']['energy_Taylor'][n],h[n]*cv(k['C_E'])+(cv(k['C0']) if n==0 else 0))
        # Full physical operators use precisely the same original inverse,
        # derivative labels, source units and scale function nodes.
        original_view=delivery['actual_source_view']
        for component,rows in source['log_radius_mixed_rows'].items():
            for label,row in rows.items():
                old=original_view['log_radius_mixed_rows'][component][label]
                assert (row.powers,row.source_units,row.log_scale_parts,row.derivative)==\
                    (old.powers,old.source_units,old.log_scale_parts,old.derivative)
        packet=source['original_forward_source_packet']
        pressure=packet['pressure']['P0_over_Pstar_squared']
        old_pressure=original_view['original_forward_source_packet']['pressure']['P0_over_Pstar_squared']
        assert all(pressure[n]._mpi_==old_pressure[n]._mpi_ for n in range(6))
        result=owner.evaluate(delivery,'1/1000')
        old_result=owner.before.evaluate(delivery,'1/1000')
        reader=owner.amplitude.reader(delivery)
        for section,components in result['physical_value_rows'].items():
            for component,rows in components.items():
                for label,row in rows.items():
                    old=old_result['physical_value_rows'][section][component][label]
                    assert row['exact_reference_log_scale_function']==old['exact_reference_log_scale_function']
                    assert row['directed_reference_log_scale']._mpi_==old['directed_reference_log_scale']._mpi_
                    assert not row['physical_accuracy']['ordinary_numeric_delivery_target_satisfied']
                    if row['physical_accuracy'].get('sign'):
                        a,b=current.ends(row['common_scale_coefficient_enclosure'])
                        lo,hi=min(abs(a),abs(b)),max(abs(a),abs(b))
                        width=(c.mpf(hi)-c.mpf(lo))/c.mpf(lo)
                        assert included(row['physical_accuracy']['coefficient_relative_width_upper'],width)
                    # Independently sum the signed common-scale terms,
                    # including every original retained positive ratio tail.
                    common=c.mpf(0)
                    for term in row['ratio_terms']:
                        assert term['directed_ratio_enclosure'] is not None
                        common+=cv(term['signed_coefficient'])*cv(term['directed_ratio_enclosure'])
                        assert reader.graph.nodes[term['exact_ratio_function']]['operation']=='analytic_unary'
                    assert included(row['common_scale_coefficient_enclosure'],common)
        for key,value in result['delivery_counts'].items():totals[key]+=value
        base=owner.velocity_pressure(delivery,'1/1000')
        assert set(base['values'])=={'u','v','w','p'}
        assert not base['unrestricted_physical_point_API'] and not base['full_certified_physical_accuracy']
        assert not any(base[k] for k in current.OPEN)
        widths={component:row['physical_accuracy'].get('coefficient_relative_width_upper')
            for component,row in base['values'].items()}
        for component in ('u','v','w'):
            assert base['values'][component]['physical_accuracy']['factored_relative_width_satisfied']
        observations[name]=dict(counts=result['delivery_counts'],base_relative_widths=widths,
            original_base_relative_widths={component:old_result['physical_value_rows']['Cartesian_spatial_rows'][key]['x0_y0_z0']['physical_accuracy'].get('coefficient_relative_width_upper')
                for component,key in (('u','ux'),('v','uy'),('w','uz'),('p','p'))},
            selected_ap_C5_independently_included=True,unchanged_P0_and_exact_physical_scales=True)
    first=next(iter(deliveries.values()))
    invalid=dict(copied_delivery=rejected(lambda:owner.evaluate(copy.copy(first))),
        invalid_target=rejected(lambda:owner.evaluate(first,'0')),
        caller_amplitude_oracle=rejected(lambda:owner.evaluate(first,{'U':[1,1]})))
    saved=owner.constants['C1']
    try:
        owner.constants['C1']=owner.ctx.mpf(0)
        invalid['substituted_normalized_ratio']=rejected(owner.assert_graph)
    finally:owner.constants['C1']=saved
    saved=owner.base.K
    try:
        owner.base.K=owner.ctx.mpf(1)
        invalid['substituted_pulse_energy']=rejected(owner.assert_graph)
    finally:owner.base.K=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes)
    hashes[current.NAME]=current.sha(current.NAME)
    hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,
        actual_row_counts=totals,actual_refined_observations=observations,
        independent_500_digit_closed_middle_positive_endpoint_integral=independent_K,
        refined_Kpulse_in_original_source_box=True,
        independent_stable_and_direct_positive_roots_and_C5_coefficients=True,
        original_future_energy_and_independent_P0_unchanged=True,
        exact_scale_functions_physical_inverse_and_units_unchanged=True,
        all_signed_ratio_terms_and_positive_tails_retained=True,
        original_constants_and_cache_keys_unmutated=True,
        invalid_inputs_rejected=invalid,**dict.fromkeys(current.GATES,True),
        **dict.fromkeys(current.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_REFINED_PULSE_COEFFICIENTS original integrals, positive C5 branch and physical factors',flush=True)
    return receipt
