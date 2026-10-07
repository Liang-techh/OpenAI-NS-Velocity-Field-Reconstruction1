"""Focused common-unit, source-row, factor preservation and recovery checks."""
import gzip
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_source_packets as source


def overlap(a,b,label):
    al,ah=source.recovery.endpoints(a);bl,bh=source.recovery.endpoints(b)
    if ah<bl or bh<al:raise ArithmeticError('Independent saved source covers disagree: '+label)


def modal_overlap(a,b,label):
    if a.algebra is not b.algebra:raise ValueError('Comparison must retain identical source bases')
    zero=a.ctx.mpf(0);count=0
    for key in set(a.terms)|set(b.terms):
        left=a.terms.get(key);right=b.terms.get(key)
        for j in range(min(a.order,b.order)+1):
            overlap(left[j] if left is not None else zero,right[j] if right is not None else zero,label)
            count+=1
    return count


def plain_cover_overlap(a,b,label):
    """Compare total intervals where older plain caches already combined modes.

    Plain energy covers contain the previously enclosed V^2/S^2 term in
    their coefficient interval. Its mode cannot be compared to a distinct
    formal S^-2 term coefficient by coefficient. Only reciprocal Pstar
    normalization (representable arbitrary-precision exponent) is used in
    this range check; microscopic width/F0/swirl factors are never resolved.
    The normalization does not select any point/source value from a cover.
    """
    def combined(row):
        out=source.IntervalTaylor.constant(row.ctx,0,row.order)
        for powers,coefficient in row.terms.items():
            if any(powers[j] for j in (0,2,3)) or powers[1]>0:
                raise ValueError('Plain-cover comparison permits reciprocal Pstar modes only')
            out=out+coefficient*row.ctx.exp(row.algebra.logs[1]*powers[1])
        return out
    left,right=combined(a),combined(b)
    for j in range(min(left.order,right.order)+1):overlap(left[j],right[j],label)
    return min(left.order,right.order)+1


def source_checks(packet,field):
    count=0;algebra=packet.algebra
    compare=plain_cover_overlap if packet.provenance['raw_record_form']=='converted_plain_raw_cover' else modal_overlap
    for name,rows in packet.velocity.items():
        for j,row in enumerate(rows):
            restored=algebra.shift(row,(0,.5,0,0)) if name in ('axial','radial') else row
            if source.encode(restored)!=source.encode(packet.native_velocity[name][j]):
                raise ArithmeticError('Native velocity units or source factors lost')
            count+=compare(row,field['physical_velocity_pressure_ordinary_y_rows'][name][j],packet.chart+' velocity '+name)
    for name,rows in packet.histories.items():
        for j,row in enumerate(rows):
            restored=algebra.shift(row,(0,.5,0,0)) if name in ('m','k') else row
            if source.encode(restored)!=source.encode(packet.native_histories[name][j]):
                raise ArithmeticError('Native history units or source factors lost')
            count+=compare(row,field['own_normalized_history_ordinary_y_rows'][name][j],packet.chart+' history '+name)
    for j,row in enumerate(packet.absolute_pressure):
        count+=compare(row,field['physical_velocity_pressure_ordinary_y_rows']['pressure'][j],packet.chart+' pressure')
    if field['original_axis_pressure_over_S_squared'] is not packet.P0:
        raise ArithmeticError('Recovery replaced the original axis datum')
    if any(field[k] is not False for k in source.OPEN):raise ArithmeticError('Original covers promoted missing stages')
    if packet.algebra.proofs or packet.algebra.final_rows:raise ArithmeticError('Source factors were resolved')
    return count


def changed_factored_fixture():
    """Nonzero own histories, mixed source factors and unchanged P0; no actual-factor evaluation."""
    c=MPIntervalContext();c.dps=100
    algebra=source.FactoredAlgebra(c,(c.mpf('-2'),2*c.ln(3),2*c.ln(2),c.mpf('-1')),[])
    Z=source.IntervalTaylor.variable(c,'.2',5)
    family=dict(zip(source.FAMILY_KEYS,('fixture','fixture_source','fixture_datum')))
    def mode(value,powers):return algebra.shift(value,powers)
    original={k:mode(1+Z/10, (j%2,-.5,j%2,0)) for j,k in enumerate(source.recovery.RATES)}
    defect={k:mode((j+1)/100+Z*Z/50,(0,-1,0,.5)) for j,k in enumerate(source.recovery.RATES)}
    P0=algebra.lift(-2+Z/3)
    E=[mode(1+Z*Z/10,(0,0,0,.5))]+[mode((j+1)/50+Z/100,(0,0,0,.5)) for j in range(4)]
    V=[mode('.7' if j==0 else '.03',(0,-.5,0,0)) for j in range(5)]
    packet=source.CurrentSourcePacket('fixture',family,algebra,Z,{}, {},(),P0,{}, {},{})
    state=source.FactoredRecoveryState(packet,original=original,defect=defect)
    factored=state.field_rows(Z=Z,delta='.001',E_rows=E,V_rows=V)
    def fixture_value(row):
        # Only modest artificial fixture logs are exponentiated. Production
        # packets retain all source modes and never call this evaluator.
        value=source.IntervalTaylor.constant(c,0,row.order)
        for key,coefficient in row.terms.items():
            factor=c.exp(sum((log*p for log,p in zip(algebra.logs,key)),c.mpf(0)))
            value=value+coefficient*factor
        return value
    plain=source.recovery.GenericMomentRecovery(c,source_family=family,P0=fixture_value(P0),
        original={k:fixture_value(v) for k,v in original.items()},defect={k:fixture_value(v) for k,v in defect.items()})
    expected=plain.field_rows(Z=Z,delta='.001',E_rows=[fixture_value(v) for v in E],V_rows=[fixture_value(v) for v in V])
    count=0
    for group in ('own_normalized_history_ordinary_y_rows','physical_velocity_pressure_ordinary_y_rows',
                  'physical_primitive_ordinary_y_rows','full_signed_stress_ordinary_y_rows'):
        for name,rows in factored[group].items():
            for j,row in enumerate(rows):
                actual=fixture_value(row);reference=expected[group][name][j]
                for k in range(min(actual.order,reference.order)+1):
                    overlap(actual[k],reference[k],'changed factored '+group+' '+name);count+=1
    if factored['original_axis_pressure_over_S_squared'] is not P0 or algebra.final_rows or algebra.proofs:
        raise ArithmeticError('Factored own-field recovery changed P0 or resolved factors')
    return count


def guards(service):
    view,origin=service.saved_view('bridge_first');algebra=service.saved('bridge_first').algebra
    calls=[lambda:service.saved('unknown'),lambda:service.saved('bridge_first','macro_whole'),
           lambda:service.query('bridge_first',0,'.5'),
           lambda:source.CurrentSourcePackets(owners={'unknown':object()})]
    wrong=source.FactoredAlgebra(algebra.ctx,algebra.logs,[])
    calls.append(lambda:source.factored_check_jets(algebra.ctx,[algebra.lift(1),wrong.lift(1)]))
    calls.append(lambda:source.decode_row(algebra,dict(axial_order=5,fixed_basepoint_log_factors=True,
        terms=[dict(source_exponents=[0,.3,0,0],axial_Taylor_coefficients=[0]*6)])))
    bad={**view,'current_unresolved_raw_source_rows':{**view['current_unresolved_raw_source_rows'],
         'inverse_width_source_shift_applied_before_any_resolution':False}}
    calls.append(lambda:service.adapt('bridge_first',bad,service.ctx,origin))
    bad_domain={**view,'current_actual_source_stress_packet':{**view['current_actual_source_stress_packet'],'Z':2}}
    calls.append(lambda:service.adapt('bridge_first',bad_domain,service.ctx,origin))
    for call in calls:
        try:call()
        except (ValueError,TypeError,KeyError):pass
        else:raise ArithmeticError('Invalid source/basis/coordinate input admitted')
    return len(calls)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed typed source dependency: '+name)
    covers=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    service=source.CurrentSourcePackets();count=0
    if list(covers)!=list(source.CHARTS) or manifest['current_original_chart_count']!=len(covers):
        raise ValueError('Original source chart inventory differs')
    for chart in source.CHARTS:
        packet=service.saved(chart);field=packet.recover_original(service.data['delta'])
        if source.encode(dict(packet=packet.record(),original_recovered_field=field))!=covers[chart]:
            raise ValueError('Saved factored source/recovery differs: '+chart)
        count+=source_checks(packet,field)
    if source.encode(source.unit_theorem())!=manifest['common_unit_theorem']:
        raise ValueError('Source unit theorem differs')
    result=dict(all_passed=True,**{source.GATE:True,source.RECOVERY_GATE:True},**dict.fromkeys(source.OPEN,False),
        source_family=service.family,current_original_chart_count=len(covers),
        source_unit_identities=len(manifest['common_unit_theorem']['identities']),
        actual_source_recovery_coefficients_checked=count,
        factored_cache_modes_compared_before_resolution=True,
        older_plain_energy_covers_compared_after_reciprocal_unit_normalization=True,
        nonzero_changed_factored_fixture_coefficients=changed_factored_fixture(),
        invalid_packet_basis_coordinate_guards=guards(service),
        source_factors_resolved_in_production=False,ancestor_constructors_called=False,
        successful_arbitrary_coordinate_live_owner_query_tested=False,
        scope='Receipt-bound original saved box covers and factored recovery algebra; cache not a point function or a new generic-loop whole-current installation.',
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Typed current source covers and factored recovery PASS:',len(covers),'charts,',count,'source coefficients',flush=True)
    return result


if __name__=='__main__':run()
