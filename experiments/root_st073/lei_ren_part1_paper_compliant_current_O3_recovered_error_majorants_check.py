"""Independent original recovery-row fixtures and current variable-N queries."""
import json
import math
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import mpmath as mp
import lei_ren_part1_paper_compliant_current_O3_recovered_error_majorants as source
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows


def sup(value):
    # mpmath scalar abs() can round to the unrelated global precision.
    return max(mp.make_mpf((0,v._mpf_[1],v._mpf_[2],v._mpf_[3])) for v in source.endpoints(value))


def original_interval_row_fixtures():
    c=MPIntervalContext();c.dps=90;count=0;comparisons=0
    for delta in ('0','.03','.5','.9'):
        factors=source.radial_factor_caps(c,c.mpf(delta))
        for zv in ('-1','-.7','0','.3','1'):
            z=IntervalTaylor.variable(c,c.mpf(zv),6)
            C=(1+source.histories.square(z)).reciprocal()
            for sign in (1,-1):
                old=[c.mpf(v) for v in ('2','-1','.5','-.25','.125')]
                inc=[sign*c.mpf(v) for v in ('.13','-.22','.31','-.4','.49')]
                axial=[-sign*c.mpf(v) for v in ('.6','-.5','.4','-.3','.2')]
                initial={key:sign*c.mpf(v) for key,v in zip(('m','h','k','e','p'),('.7','-.6','.5','-.4','.3'))}
                # Execute the original directed IntervalTaylor program, not
                # the new positive-cap recurrence or symbolic replay field.
                d,Q=source.histories.increment_rows(c,z,c.mpf(delta),c.mpf(1),initial,
                    [C*(old[j]+inc[j]) for j in range(5)],[C*v for v in old],
                    [C*v for v in axial],[C*v for v in inc])
                native_radial=shifted_rows(Q,c.mpf('.5'),4)
                caps=source.scalar_history_caps({key:abs(v) for key,v in initial.items()},
                    [abs(v) for v in old],[abs(v) for v in inc],[abs(v) for v in axial])
                bounded=source.recover(c,caps,[abs(v) for v in axial],factors)
                for j in range(5):
                    for k in range(6):
                        for label,jet,cap in (
                            ('radial',native_radial[j],bounded['radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5'][j][k]),
                            ('pressure',d['p'][j],bounded['absolute_pressure_error_over_Pstar2_Ad2_ordinary_logR4_axial5'][j][k])):
                            if sup(jet[k]*math.factorial(k))>source.endpoints(cap)[1]:
                                raise ArithmeticError('Original IntervalTaylor '+label+' row exceeds mixed cap')
                            comparisons+=1
                        for key in d:
                            if sup(d[key][j][k]*math.factorial(k))>source.endpoints(bounded['five_normalized_history_error_ordinary_logR4_axial5'][key][j][k])[1]:
                                raise ArithmeticError('Original signed history %s exceeds mixed cap:delta=%s,Z=%s,sign=%s,j=%s,k=%s;source=%s,cap=%s'%
                                    (key,delta,zv,sign,j,k,d[key][j][k]*math.factorial(k),bounded['five_normalized_history_error_ordinary_logR4_axial5'][key][j][k]))
                            comparisons+=1
                count+=1
    return count,comparisons


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed recovered error input: '+name)
    field=source.CurrentRecoveredErrorMajorants(require_checked=False)
    if source.parameters.encoded(field.theorem)!=data['exact_recovered_error_source_theorem']:
        raise ValueError('Original recovery source theorem changed')
    specs=(('before_modulation','modulation','-3',1),('left_taper','modulation','-1.99',37),
        ('old_seam','modulation','0',37),('after_cutoff','modulation','.75',37),
        ('crossing_seam','modulation',('-.1','.1'),10**12),
        ('quiet_whole','quiet_repair',(0,1),field.repair.quiet_threshold),
        ('first_bump','quiet_repair',('.19','.21'),10**12),
        ('near_last_exit','quiet_repair',('.824','.825'),field.repair.quiet_threshold),
        ('terminal','quiet_repair','1',field.repair.quiet_threshold))
    for name,method,coordinate,N in specs:
        fresh=getattr(field,method)(coordinate,N) if method=='modulation' else field.quiet_repair(N,coordinate)
        if source.parameters.encoded(fresh)!=data['examples'][name]:raise ValueError('Changed recovered query: '+name)
        for label in ('radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5',
                      'absolute_pressure_error_over_Pstar2_Ad2_ordinary_logR4_axial5'):
            rows=fresh[label]
            if len(rows)!=5 or any(len(row)!=6 for row in rows):raise ValueError('Incomplete mixed recovery rows')
        for key,rows in fresh['normalized_signed_history_error_enclosures'].items():
            for j,row in enumerate(rows):
                for k,v in enumerate(row):
                    lo,hi=source.endpoints(v)
                    cap=fresh['five_normalized_history_error_ordinary_logR4_axial5'][key][j][k]
                    if lo>0 or hi<source.endpoints(cap)[1] or sup(field.ctx.mpf(lo))<source.endpoints(cap)[1]:
                        raise ArithmeticError('Signed enclosure lost one endpoint at global scalar precision')
    def is_zero(query):
        return all(v._mpi_==field.ctx.mpf(0)._mpi_ for label in (
            'radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5',
            'absolute_pressure_error_over_Pstar2_Ad2_ordinary_logR4_axial5') for row in query[label] for v in row)
    before=field.modulation('-11',1);terminal=field.quiet_repair(field.repair.quiet_threshold,'1')
    after_support=field.quiet_repair(field.repair.quiet_threshold,('.9','1'))
    if not all(is_zero(q) for q in (before,terminal,after_support)):
        raise ArithmeticError('Original left/terminal source error jets do not vanish')
    after=field.modulation('.75',37)
    if any(v._mpi_!=field.ctx.mpf(0)._mpi_ for v in after['local_cutoff_derivative_caps']):
        raise ArithmeticError('Post-cutoff local profile is not exactly original')
    if source.endpoints(after['normalized_scalar_history_error_ordinary_logR_caps']['p'][0])[0]<=0:
        raise ArithmeticError('Pressure history was replaced by zero local pressure density')
    if any(v._mpi_!=field.ctx.mpf(0)._mpi_ for v in after['normalized_scalar_history_error_ordinary_logR_caps']['p'][1:]):
        raise ArithmeticError('Post-cutoff cumulative pressure should be radially constant')
    if source.endpoints(field.modulation('0',37)['positive_kinetic_history_scalar_enclosure'])[0]<=0:
        raise ArithmeticError('Separate nonzero kinetic buffer mass was dropped')
    whole=field.quiet_repair(field.repair.quiet_threshold,(0,1))
    near=field.quiet_repair(field.repair.quiet_threshold,('.824','.825'))
    for label in ('radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5',
                  'absolute_pressure_error_over_Pstar2_Ad2_ordinary_logR4_axial5'):
        if not source.endpoints(near[label][0][0])[1]<source.endpoints(whole[label][0][0])[1]:
            raise ArithmeticError('Near-exit recovered value bound fails to shrink')
    invalid=(('modulation',-12,1),('modulation',2,1),('modulation',mp.inf,1),('modulation',0,True),
             ('quiet_repair',0,field.repair.repair_threshold-1),('quiet_repair',(-1,0),10**12))
    for method,coordinate,N in invalid:
        try:
            if method=='modulation':field.modulation(coordinate,N)
            else:field.quiet_repair(N,coordinate)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid recovery query accepted')
    for key in source.OPEN:
        if data.get(key) is not False:raise ValueError('Recovery bounds promoted a whole-cone/global gate')
    fixtures,comparisons=original_interval_row_fixtures()
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_original_recovery_radial_shift_and_inverse_L_identities=len(field.theorem['identities']),
        actual_variable_N_modulation_and_quiet_queries_recomputed=len(specs),
        independent_original_IntervalTaylor_recovery_fixtures=fixtures,
        independent_signed_history_radial_pressure_mixed_jet_comparisons=comparisons,
        invalid_queries_rejected=len(invalid),
        cumulative_pressure_value_not_replaced_by_local_derivative=True,
        separate_positive_kinetic_mass_and_same_axis_datum_retained=True,
        left_and_post_repair_native_error_mixed_jets_exactly_zero=True,
        local_O2_O3_quiet_radial_pressure_error_jets_certified=True,
        **{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Recovered native errors PASS:',len(field.theorem['identities']),'source identities;',fixtures,
        'original interval fixtures;',comparisons,'mixed row comparisons',flush=True)
    return result


if __name__=='__main__':run()
