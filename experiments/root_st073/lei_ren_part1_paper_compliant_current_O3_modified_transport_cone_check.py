"""Focused transport history retention, signed cone and common-N checks."""
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_O3_modified_transport_cone as source


def independent_transport_fixtures(field):
    """Moderate signed histories test derivative and radius units separately."""
    count=0
    with mp.workdps(100):
        for mu in (mp.mpf('.001'),mp.mpf('1e-20')):
            alpha=mp.mpf('.5')+mu
            for v in (mp.mpf(0),mp.mpf('.537'),mp.mpf(1)):
                def u(x):return mp.mpf('1.7')*mp.exp(-alpha*x)
                a=1-2*mp.diff(u,v)/u(v)
                if abs(a-(2+2*mu))>mp.mpf('1e-90'):
                    raise ArithmeticError('Independent actual transport shear differs')
                count+=1
                for initial,rate in ((mp.mpf('-.8'),1),(mp.mpf('.4'),mp.mpf('1.5')),
                    (mp.mpf('-.7'),mp.mpf('1.5')),(mp.mpf('.9'),1),(mp.mpf('-.3'),0)):
                    def history(x):return initial*mp.exp(-rate*(1+x))
                    for j in range(5):
                        value=mp.diff(history,v,j);expected=(-rate)**j*history(v)
                        if abs(value-expected)>mp.mpf('1e-90'):
                            raise ArithmeticError('Independent ordinary history row differs')
                        count+=1
                # The pure-m radial source Q=-M acquires sqrt(R) once.
                def radial(x):return -mp.mpf('.8')*mp.exp(-(1+x))
                def physical_radial(x):return mp.exp(x/2)*radial(x)
                for j in range(5):
                    shifted=sum(mp.binomial(j,i)*mp.mpf('.5')**(j-i)*mp.diff(radial,v,i)
                        for i in range(j+1))
                    if abs(mp.diff(physical_radial,v,j)/mp.exp(v/2)-shifted)>mp.mpf('1e-90'):
                        raise ArithmeticError('Independent radial radius shift differs')
                    count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed transport/common-N source: '+name)
    field=source.CurrentModifiedTransportCone(require_checked=False);c=field.ctx
    for key,value in (('exact_transport_source_and_native_units_theorem',field.theorem),
        ('original_local_power_baseline',field.baseline),('whole_modified_pre_repair_transport_cone',field.proof),
        ('complete_leading_stress_error_log_ledger',field.error_rows),('scoped_common_frequency_ledger',field.common)):
        if source.parameters.encoded(value)!=data[key]:raise ValueError('Transport source proof/ledger changed: '+key)
    for name,value in field.proof['positive_margins'].items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('Transport cone reserve unresolved: '+name)
    if field.proof['original_signed_D_over_theta_exact']!=1:
        raise ArithmeticError('Actual bs=0 signed D quotient differs')
    if source.endpoints(field.proof['original_signed_Q_over_theta_squared_lower'])[0]<=mp.mpf('1.9'):
        raise ArithmeticError('Transport signed quadratic reserve missing')
    specs=(('whole_bridge',(0,1),field.minimum_N),('O3_to_power_seam','0',field.minimum_N),
        ('interior','.537',field.minimum_N),('repair_entry','1',field.minimum_N),('current_frequency',(0,1),10**12))
    histories=0
    for name,v,N in specs:
        query=field.query(v,N)
        if source.parameters.encoded(query)!=data['examples'][name]:raise ValueError('Transport seam/query changed: '+name)
        q,caps=field.inputs(v,N)
        for key in ('m','h','k','e','p'):
            if source.endpoints(q['normalized_scalar_history_caps'][key][0])[1]<=0:
                raise ArithmeticError('Surviving transport history erased: '+key)
            histories+=1
        for key in ('du','V'):
            if any(source.endpoints(value)!=(0,0) for value in caps[key].values()):
                raise ArithmeticError('Post-cutoff local modulation is not exact zero')
        if source.endpoints(caps['dr'][0,0])[1]<=0:
            raise ArithmeticError('Own-moment radial error erased with local velocity')
        if any(source.endpoints(q['normalized_scalar_history_caps']['p'][j])!=(0,0) for j in range(1,5)):
            raise ArithmeticError('Constant pressure history acquired logR derivatives')
    zeros=nonzeros=comparisons=0
    for label,parts in field.error_rows.items():
        for name,row in parts.items():
            if row['source_error_exactly_zero']:
                if source.endpoints(row['native_coefficient_absolute_cap'])!=(0,0) or row['relative_error_log_parts'] is not None or row['log_absolute_error_over_baseline_theta_lower'] is not None:
                    raise ArithmeticError('Exact zero stress error uses a finite logarithm')
                zeros+=1
            else:
                if source.endpoints(-1010-row['log_absolute_error_over_baseline_theta_lower'])[0]<=0:
                    raise ArithmeticError('Transport full error lacks local baseline margin')
                nonzeros+=1
    if not zeros or not nonzeros:raise ArithmeticError('Transport must retain local zeros and cumulative stress changes')
    for N in (2*field.minimum_N,10**12):
        for label,parts in field.leading_error_rows(N).items():
            for name,row in parts.items():
                if source.endpoints(row['native_coefficient_absolute_cap'])[1]>source.endpoints(field.error_rows[label][name]['native_coefficient_absolute_cap'])[1]:
                    raise ArithmeticError('Larger N enlarged transported native error cap')
                comparisons+=1
    common=field.common
    if common['sufficient_common_integer_N']!=max(common['regional_threshold_ledger'].values()):
        raise ArithmeticError('Common N did not include every changed-region threshold')
    if common['sufficient_common_integer_N']!=field.quiet.minimum_N or not common['earlier_O2_buffer_and_exact_degenerate_left_edge_excluded_from_strict_certificate']:
        raise ArithmeticError('Changed-region common N incorrectly covered degenerate O2')
    if len(common['accepted_regional_certificates'])!=2:
        raise ArithmeticError('Checked modulation certificates missing')
    invalid=((-1,field.minimum_N),(2,field.minimum_N),(('.99','1.01'),field.minimum_N),
        (0,field.minimum_N-1),(0,True),(0,1.5),(mp.inf,field.minimum_N))
    for v,N in invalid:
        try:field.query(v,N)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid transport/common-N query admitted')
    for key in source.OPEN:
        if data.get(key) is not False or common.get(key) is not False:
            raise ValueError('Scoped transport/common N promoted global gate')
    fixtures=independent_transport_fixtures(field)
    result=dict(all_passed=True,source_family=data['source_family'],minimum_integer_N=field.minimum_N,
        exact_source_transport_native_unit_and_all_N_identities=len(field.theorem['identities']),
        strict_transport_cone_inequalities=len(field.proof['positive_margins']),
        independent_signed_transport_shear_history_and_radial_unit_fixtures=fixtures,
        surviving_cumulative_history_queries=histories,exact_local_zero_stress_pieces=zeros,
        surviving_complete_stress_error_pieces=nonzeros,larger_N_native_error_comparisons=comparisons,
        source_seam_and_bridge_queries_recomputed=len(specs),invalid_queries_rejected=len(invalid),
        scoped_common_frequency_includes_all_changed_region_constraints=True,
        current_modified_pre_repair_power_transport_cone_certified=True,
        current_modified_modulation_transport_and_quiet_common_N_certified=True,
        **{key:False for key in source.OPEN},input_hashes={**data['input_hashes'],
            source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.parameters.encoded(result),indent=2)+'\n').encode())
    print('Transport cone and scoped common N PASS: N>=',field.minimum_N,';',nonzeros,'retained stress pieces;',fixtures,'independent fixtures',flush=True)
    return result


if __name__=='__main__':run()
