"""Independent signed original-program fixtures and physical error checks."""
import gzip
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_O3_completed_tensor_error_majorants as source
import lei_ren_part1_paper_compliant_current_modified_pre_stress_operator as original
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket


def sup(v):
    return max(mp.make_mpf((0,x._mpf_[1],x._mpf_[2],x._mpf_[3])) for x in source.endpoints(v))


def original_signed_fixtures():
    """Execute the original IntervalTaylor programs, not symbolic caps."""
    c=MPIntervalContext();c.dps=90;model,_=source.source_model();count=comparisons=physical_comparisons=0
    for delta,zv,sign in (('0','-1',1),('.03','-.7',-1),('.5','0',1),('.9','.3',-1),('.9','1',1)):
        z=IntervalTaylor.variable(c,c.mpf(zv),5);C=(1+z*z).reciprocal();y=c.mpf('.17')
        row={};caps={}
        for n,name in enumerate(source.FUNCTIONS):
            # Independent ordinary logR derivative rows of a polynomial,
            # with signed rational/polynomial axial factors.
            coefficients=[c.mpf(sign*(-1)**(n+j)*(n+1))/(10*(j+1)) for j in range(5)]
            shape=C if n%3==0 else C*C if n%3==1 else 1+z/3-z*z/7
            raw=[shape*sum(coefficients[i]*math.factorial(i)/math.factorial(i-j)*y**(i-j)
                for i in range(j,5)) for j in range(5)]
            caps[name]={(j,k):source.upper(c,abs(raw[j][k]*math.factorial(k))) for j,k in source.indices(4)}
            row[name]=shifted_rows(raw,c.mpf('.5'),4) if name in ('R0','dr') else raw
        for Ad in ('.03','1'):
            modified=source.source_view(row,c.mpf(Ad));base=source.source_view(row,c.mpf(0))
            for kind,fn,total in (('stress',original.modified_pre_stress_rows,3),
                                  ('remainder',original.modified_pre_remainder_sectors,2)):
                changed=fn(c,c.mpf(delta),z,modified);unchanged=fn(c,c.mpf(delta),z,base)
                bounded={label:{} for label in changed}
                for label,parts in model[kind].items():
                    for name,part in parts.items():
                        grid=source.bound_expression(c,part['expression'],caps,c.mpf(zv),c.mpf(delta),total,part['mode'][0])
                        oldname=part['original_sector']
                        for jk,cap in grid.items():
                            bounded[label].setdefault(oldname,{})[jk]=bounded[label].get(oldname,{}).get(jk,c.mpf(0))+c.mpf(Ad)**part['error_Ad_power']*cap
                for label,parts in changed.items():
                    for name,part in parts.items():
                        if name not in bounded[label]:continue # exact AST-proved baseline cancellation
                        key='full_derivative_rows' if kind=='stress' else 'rows'
                        for j,k in source.indices(total):
                            value=(part[key][j]-unchanged[label][name][key][j])[k]*math.factorial(k)
                            # Directed subtraction of two baseline enclosures
                            # adds numerical width, even for a tight triangle
                            # equality. Allow only the working-precision floor;
                            # the exact source proof establishes the real bound.
                            slack=c.mpf(2)**(-c.prec+20)*(1+sup(part[key][j][k]*math.factorial(k))+
                                sup(unchanged[label][name][key][j][k]*math.factorial(k)))
                            if sup(value)>source.endpoints(bounded[label][name][j,k]+slack)[1]:
                                raise ArithmeticError('Signed original source exceeds error cap: %s/%s/%s/y%d/Z%d'%(kind,label,name,j,k))
                            comparisons+=1
                        actualgrid={'y%d_Z%d'%(j,k):(part[key][j]-unchanged[label][name][key][j])[k]*math.factorial(k)
                            for j,k in source.indices(total)}
                        capgrid={'y%d_Z%d'%jk:c.mpf([-1,1])*source.upper(c,v) for jk,v in bounded[label][name].items()}
                        beta=-2-c.mpf(delta) if kind=='stress' else part['beta']
                        for a,b in source.indices(total):
                            actual=physical_bracket(c,actualgrid,a,b,c.mpf(zv),c.mpf(delta),beta)
                            cap=physical_bracket(c,capgrid,a,b,c.mpf(zv),c.mpf(delta),beta)
                            slack=c.mpf(2)**(-c.prec+20)*(1+sup(actual)+sup(cap))
                            if sup(actual)>source.endpoints(c.mpf(sup(cap))+slack)[1]:
                                raise ArithmeticError('Original physical derivative exceeds mapped error bound')
                            physical_comparisons+=1
            count+=1
    return count,comparisons,physical_comparisons


def original_shift_checks(field):
    """Compare recovered normalized M/R and Q jets with independent raw rows."""
    c=field.ctx;read=source.parameters.numeric.transport.read_interval
    views=json.loads(gzip.decompress((source.HERE/source.BASE_VIEWS).read_bytes()));count=0
    for chart,key in (('O2','O2_full_buffer'),('O3','O3_full_transition'),('quiet','quiet_full_repair')):
        pre=views[key]['actual_source']['current_original_pre_source']
        for n,label,rate,rawkey in (
            ('M0','Mz_over_current_R',c.mpf(1),'actual_normalized_primitive_y_derivative_axial5'),
            ('R0','Ur_over_current_sqrt_R_over_2',c.mpf('.5'),'actual_Q_y_derivative_axial4')):
            grid=(pre['physical_five_primitive_y_Z_mixed4'] if n=='M0' else pre['physical_velocity_pressure_y_Z_mixed4'])[label]
            for j,k in source.indices(4):
                # Shift and inverse shift compose as binomial polynomials.
                terms=[sum(math.comb(i,a)*(-rate)**(i-a)*read(c,grid['y%d_Z%d'%(a,k)])
                    for a in range(i+1)) for i in range(j+1)]
                restored=sum(math.comb(j,i)*rate**(j-i)*terms[i] for i in range(j+1))
                original_value=read(c,grid['y%d_Z%d'%(j,k)])
                lo,hi=source.endpoints(restored);a,b=source.endpoints(original_value)
                if lo>a or hi<b:raise ArithmeticError('Physical radius shift failed to round-trip')
                count+=1
    return count


def check_physical_rows(query):
    p=query['physical_error_bounds'];count=0
    for key in ('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_divergence_mixed2',
                'physical_three_component_remainder_mixed2'):
        for parts in p[key].values():
            for grid in parts.values():
                for row in grid.values():
                    if row['exact_zero']:
                        if row['log_absolute_upper'] is not None or sup(row['signed_coefficient'])!=0:
                            raise ArithmeticError('Zero physical error sector has a nonzero bound')
                    else:
                        if row['log_absolute_upper'] is None:raise ArithmeticError('Physical error log bound missing')
                        if not row['positive_source_factors_not_materialized']:raise ArithmeticError('Huge factors materialized')
                        parts=row['actual_source_log_parts']
                        if 'original_Rd_amplitude' not in parts:raise ArithmeticError('Additional Ad factor missing')
                    count+=1
    if p['exact_completed_tensor_radial_divergence']._mpi_!=mp.iv.mpf(0)._mpi_:
        raise ArithmeticError('Completed radial divergence must be exactly zero')
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    examples=json.loads(gzip.decompress((source.HERE/source.VIEWS_NAME).read_bytes()))
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed completed error input: '+name)
    field=source.CurrentCompletedTensorErrorMajorants(require_checked=False);c=field.ctx
    if source.parameters.encoded(field.theorem)!=data['exact_original_source_difference_theorem']:
        raise ValueError('Original source difference theorem changed')
    fixtures,comparisons,physical_comparisons=original_signed_fixtures();shifts=original_shift_checks(field)
    specs=(('before_support','O2','-3',1),('left_taper','O2','-1.99',37),
        ('transition','O3',('0','.5'),10**12),('after_cutoff','O3','.75',37),
        ('quiet_whole','quiet',(0,1),field.recovered.repair.quiet_threshold),
        ('terminal','quiet',('.9','1'),field.recovered.repair.quiet_threshold))
    checks=0;fresh={}
    for name,chart,coordinate,N in specs:
        q=field.query(chart,coordinate,N);fresh[name]=q
        if source.parameters.encoded(q)!=examples[name]:raise ValueError('Changed completed error query: '+name)
        checks+=check_physical_rows(q)
    timed=field.at_physical_time('O3','.2',10**12,Z=('.3','.6'),log_tau=('-3','-1'))
    if source.parameters.encoded(timed)!=examples['actual_time']:raise ValueError('Actual lambda query changed')
    p=timed['physical_error_bounds'];expected=(c.mpf([-3,-1])-c.ln(1-c.mpf(['.3','.6'])**2))/2
    if p['actual_log_lambda']._mpi_!=expected._mpi_:raise ArithmeticError('Physical time used conservative lambda')
    if p['requested_log_tau']._mpi_!=c.mpf([-3,-1])._mpi_:raise ArithmeticError('Mapper time leaked into physical metadata')
    checks+=check_physical_rows(timed)
    for name in ('before_support','terminal'):
        if any(v._mpi_!=c.mpf(0)._mpi_ for labels in fresh[name]['native_error_absolute_mixed_caps'].values()
            for parts in labels.values() for grid in parts.values() for v in grid.values()):
            raise ArithmeticError('Flat support/implicit repair exit error did not vanish')
    after=fresh['after_cutoff']['native_error_absolute_mixed_caps']
    if not any(source.endpoints(v)[1]>0 for v in after['stress']['axial']['actual_absolute_pressure_P2_Ad2'].values()):
        raise ArithmeticError('Post-cutoff cumulative absolute pressure was dropped')
    invalid=(lambda:field.query('O3',('-.1','.1'),37),lambda:field.query('O2','.1',37),
        lambda:field.query('quiet',0,field.recovered.repair.repair_threshold-1),
        lambda:field.query('O3',0,True),lambda:field.query('O3',0,37,Z=(-2,1)),
        lambda:field.at_physical_time('O3',0,37,Z=(-1,1)),
        lambda:field.query('O3',0,37,log_lambda=mp.inf),lambda:field.query('O3',0,37,viscosity=0))
    for fn in invalid:
        try:fn()
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid tensor error query accepted')
    for flag in source.OPEN:
        if data.get(flag) is not False:raise ValueError('Error bounds promoted signed/global gate')
    receipt=dict(all_passed=True,source_family=data['source_family'],
        exact_actual_source_difference_and_shift_identities=len(field.theorem['identities']),
        original_signed_IntervalTaylor_fixtures=fixtures,original_signed_source_mixed_comparisons=comparisons,
        independent_original_physical_derivative_comparisons=physical_comparisons,
        signed_fixture_roundoff_floor='2^(-ctx.prec+20) times ordinary-row magnitude; production caps are not enlarged',
        baseline_physical_radius_shift_roundtrip_checks=shifts,physical_error_log_rows_checked=checks,
        variable_N_current_queries_recomputed=len(specs)+1,invalid_queries_rejected=len(invalid),
        actual_lambda_time_correlation_and_metadata_retained=True,
        cached_modified_fixed_N_values_used=False,original_nonzero_cross_terms_and_cumulative_pressure_retained=True,
        completed_signed_tensor_error_bounds_available=True,**{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Completed tensor errors PASS:',fixtures,'signed fixtures;',comparisons,'mixed comparisons;',checks,'physical rows',flush=True)
    return receipt


if __name__=='__main__':run()
