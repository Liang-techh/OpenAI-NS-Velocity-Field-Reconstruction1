"""Integrate the original inner-reference E squared weight before enclosing p.

This reads accepted same-N source archives. It changes one pressure integral
cover, not its defining density, then retains all actual upstream memory.
"""
import ast
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_upstream_active_kappa_C1 as upstream
import lei_ren_part1_paper_compliant_current_native_upstream_paired_branch_C1_check as saved

common=upstream.common;accepted=saved.accepted;replay=accepted.replay
HERE,PREFIX,sha=upstream.HERE,upstream.PREFIX,upstream.sha
NAME=PREFIX+'current_inner_reference_weighted_pressure.json'
RECEIPT=PREFIX+'current_inner_reference_weighted_pressure_check.json'
GATE='original_inner_reference_E_squared_weight_integrated_and_actual_pressure_memory_replayed'
KEYS=upstream.KEYS;N=upstream.N;ep=upstream.ep;iv=common.interval
encode=upstream.encode;same=accepted.same


def source_identity():
    """Bind the radial weight and normalization to original source ASTs."""
    guards=[]
    def tree(stem):return ast.parse((HERE/(PREFIX+stem+'.py')).read_text(encoding='utf8'))
    def fn(stem,name):
        matches=[n for n in ast.walk(tree(stem)) if isinstance(n,ast.FunctionDef) and n.name==name]
        if len(matches)!=1:raise ValueError('Unique original callable required: '+stem+'.'+name)
        return matches[0]
    def expr(stem,name,target,expected):
        node=fn(stem,name);wanted=ast.dump(ast.parse(expected,mode='eval').body,include_attributes=False)
        found=[n for n in ast.walk(node) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets)
            and ast.dump(n.value,include_attributes=False)==wanted]
        if len(found)!=1:raise ValueError('Original source expression changed: '+stem+'.'+name+'.'+target)
        guards.append(dict(module=PREFIX+stem+'.py',callable=name,target=target,expression=expected))
    expr('reference_restore_profiles','reference','gap','(self.loggap-8)*phase')
    expr('reference_restore_profiles','reference','logu','-logq+(self.reshape.T/10-self.core.logC-self.core.logP+gap/10)')
    # __init__ is class scoped because a module may contain helper classes.
    for stem,klass,target,expected in (
        ('reference_restore_profiles','CompliantReferenceRestoreProfiles','self.loggap','self.reshape.logref-self.reshape.T'),
        ('long_reshape_profiles','CompliantLongReshapeProfiles','self.logref','10*(self.core.logC+self.core.logP)')):
        cls=next(n for n in tree(stem).body if isinstance(n,ast.ClassDef) and n.name==klass)
        node=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
        wanted=ast.dump(ast.parse(expected,mode='eval').body,include_attributes=False)
        if not any(isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)
            and ast.dump(n.value,include_attributes=False)==wanted for n in ast.walk(node)):
            raise ValueError('Original logarithmic-radius length changed')
        guards.append(dict(module=PREFIX+stem+'.py',callable=klass+'.__init__',target=target,expression=expected))
    expr('current_restore_stress_operator','raw_restore_rows','logu',
        "IntervalTaylor(c,[c.mpf(endpoints(value)) for value in parent['log_Utheta_over_Pstar_axial5_coefficients']])")
    expr('current_restore_stress_operator','raw_restore_rows','u0','positive_exp(c,logu[0])')
    expr('current_restore_stress_operator','raw_restore_rows','u','amp*u0')
    raw=fn('current_restore_stress_operator','raw_restore_rows')
    theta='[u*c.mpf(".1")**j for j in range(5)]'
    wanted=ast.dump(ast.parse(theta,mode='eval').body,include_attributes=False)
    if not any(isinstance(n,ast.keyword) and n.arg=='theta'
        and ast.dump(n.value,include_attributes=False)==wanted for n in ast.walk(raw)):
        raise ValueError('Original theta ordinary-y derivative compiler changed')
    guards.append(dict(module=PREFIX+'current_restore_stress_operator.py',callable='raw_restore_rows',target='theta',expression=theta))
    expr('current_native_signed_input_enclosures','packet','algebra',
        'packets.FactoredAlgebra(self.ctx,(self.ctx.mpf(0),2*self.native.seed.logP,self.ctx.mpf(0),2*logu),[])')
    expr('current_native_signed_input_enclosures','packet','vel',
        "{key:tuple(algebra.shift(v,packets.INVERSE_S) if key in ('axial','radial') else v for v in values) for key,values in nv.items()}")
    expr('current_generic_loop_function_sources','build','deltaE','g.mul(E,g.function("expm1",exparg))')
    # The original flux AST itself is checked; no replacement density is installed.
    source_tree=tree('current_generic_loop_function_sources')
    wanted=ast.dump(ast.parse('g.add(g.mul(E0,r0),g.mul(half,r0,r0))',mode='eval').body,include_attributes=False)
    if not any(isinstance(n,ast.keyword) and n.arg=='p' and ast.dump(n.value,include_attributes=False)==wanted
        for n in ast.walk(source_tree)):
        raise ValueError('Original pressure-density expression changed')
    actual=json.loads((HERE/(PREFIX+'actual_reference_restore_mixed_C4_check.json')).read_bytes())
    fixture=json.loads((HERE/(PREFIX+'reference_restore_mixed_C4_check.json')).read_bytes())
    if not (actual['all_passed'] and actual['current_correlated_E_installed']
        and actual['current_exact_Rz_and_restore_exit_two_sided_rows_checked']==270
        and actual['current_centered_E_source_bindings']['unchanged_original_callables']['reference_branch']
        and fixture['all_passed'] and fixture['exact_positive_source_amplitudes_retained']
        and fixture['independent_physical_fixture']['independent_physical_closed_primitive_mixed_derivatives']==135
        and fixture['independent_physical_fixture']['actual_source_admission'] is False):
        raise ValueError('Original mixed-row receipts required without global admission')
    E,x=sy.symbols('E x',real=True);d=E*(sy.exp(x)-1)
    if sy.expand(E*d+d*d/2-E*E*(sy.exp(2*x)-1)/2)!=0:
        raise ArithmeticError('Original pressure factorization failed')
    return dict(passed=True,original_source_AST_guards=guards,
        ordinary_y_coordinate='y=log R; gap=log(R/Rsh), dy=dgap',
        original_theta_unit='E=Utheta/Pstar; theta not shifted by Pstar again',
        exact_original_identity='E_y=E/10; (E²)_y=E²/5',
        exact_original_pressure_density='p=E²*expm1(2*A/N)/2',
        exact_weight_integral='integral E² dy=5*E(right)²*(1-exp(-true_width/5))',
        weight_upper='integral E² dy<=5*whole_cell_E_absolute_upper²',
        source_function_A_can_vary_in_y_and_phase=True,
        magnitude_majorant='|integral p dy|<=5/2*Emax²*(exp(2*Amax/N)-1)',
        true_width_already_integrated_not_multiplied_again=True,
        fixed_Z_geometry_and_rate_zero_p_kernel_required=True,
        no_oscillatory_cancellation_or_selector_derivative_assumed=True)


def pressure_bound(c,payload):
    coordinates,*_=replay.inputs(c,payload)
    route=payload['complete_actual_pre_O2_tile_route']
    row=route['actual_serial_chart_C1_history_records']['inner_reference']
    source=row['original_full_box_signed_C1_source'];prov=source['source_provenance']
    if prov['chart']!='inner_reference' or prov['native_provider_group']!='restore':
        raise ValueError('Original inner-reference source required')
    if not prov['original_capped_amplitude_covers_replaced_by_exact_formal_unit']:
        raise ValueError('Actual signed source amplitude basis required')
    adapter=prov['amplitude_adapter']
    if adapter['source_module']!='current_restore_stress_operator' or adapter['source_function']!='raw_restore_rows':
        raise ValueError('Original restore row compiler required')
    logu=iv(c,adapter['frozen_positive_source_log'])
    # The signed source compiler replaces the generic plain basis with 2 logu.
    bases=(c.mpf(0),coordinates.logP_squared,c.mpf(0),2*logu,iv(c,prov['logR_cover']))
    roots=source['original_correlated_shear_and_q']['correlated_shear_and_signed_root_enclosures']
    E=saved.native_restore(roots['E']['y0_Z0'],bases,coordinates.ledger)
    A=saved.native_restore(row['original_whole_period_C1_cover']['selected_primitive_C0_Z_covers']['A'],bases,coordinates.ledger)
    if E.scale.powers!=(0,0,0,.5,0) or ep(E.coefficient)!=(1,1):
        raise ValueError('Same exact formal E amplitude required')
    if ep(common.absolute_log_upper(E))!=ep(iv(c,roots['E']['y0_Z0']['log_absolute_upper'])):
        raise ValueError('Signed source basis does not reproduce original E log cover')
    if any(A.scale.powers):raise ValueError('Original whole-period A cover must remain dimensionless')
    Amax=upstream.baseline.serial.absolute_upper(A)
    logA=upstream.baseline.magnitude_log(Amax)
    if logA is None:amax=c.mpf(0)
    else:
        if logA>ep(c.ln(159))[1]:raise ValueError('Accepted bounded original primitive required')
        amax=c.exp(c.mpf(logA))
    tmax=amax*2/N
    if ep(tmax)[1]>1:raise ValueError('Only finite small original pressure exponent materialized')
    exponential=c.exp(tmax)-1
    if ep(exponential)[0]<0:raise ArithmeticError('Positive directed exponential increment required')
    EU=upstream.baseline.serial.absolute_upper(E)
    candidate=upstream.baseline.serial.symmetric_bound(coordinates.rebase(EU*EU,payload['source_family'])*(c.mpf(5)/2)*coordinates.scalar(exponential))
    old=common.restore_common_source(row['true_log_radius_signed_C0_contributions']['p'],coordinates)
    before=upstream.baseline.magnitude_log(old);after=upstream.baseline.magnitude_log(candidate)
    tighter=after is None or before is not None and after<before
    value=candidate if tighter else old
    proof=dict(original_source_provenance=prov,native_signed_source_log_bases=bases,
        original_E_C0=E.record(),original_whole_period_A_C0=A.record(),original_E_log_cover=roots['E']['y0_Z0']['log_absolute_upper'],
        whole_cell_E_absolute_upper=EU.record(),whole_period_A_absolute_upper=Amax.record(),
        finite_Amax=amax,twice_Amax_over_N=tmax,directed_positive_exponential_increment=exponential,
        candidate_pressure_integral=candidate.record(),accepted_pressure_integral=old.record(),
        selected_pressure_integral=value.record(),before_absolute_upper_log=before,after_absolute_upper_log=after,
        strict_upper_reduction=tighter and before!=after,
        original_true_width_geometry=row['actual_whole_native_chart_geometry'],
        only_integral_enclosure_changed_original_density_unchanged=True,
        Emax_and_Amax_are_function_upper_bounds_not_chosen_field_values=True,
        exact_signed_amplitude_basis_reproduces_archived_original_log_cover=True)
    return coordinates,value,proof


def pressure_path(c,payload,refined):
    coordinates,original,originalZ,P0,P0Z=replay.inputs(c,payload)
    refined=common.restore_common_source(encode(refined.record()),coordinates)
    route=payload['complete_actual_pre_O2_tile_route'];flat,rows=accepted.route_rows(route)
    prior=common.restore_common_source(flat['actual_correction_C0_enclosures']['p'],coordinates)
    oldprior=prior;path=[]
    if not prior.zero or common.previous.five.RATES['p']!=0:
        raise ValueError('Original sc/2 correction and pressure rate required')
    for label,row in rows:
        geometry_key,source_key,out,background_key,own,inherited=accepted.row_fields(label)
        g=row[geometry_key]
        if not g['width_and_endpoints_independent_of_Z']:raise ValueError('Original fixed-Z interval required')
        width=common.restore_common_source(g['positive_true_log_radius_width'],coordinates)
        if ep(width.coefficient)[0]<=0:raise ValueError('True positive original width required')
        old=common.restore_common_source(row['true_log_radius_signed_C0_contributions']['p'],coordinates)
        addition=refined if label=='inner_reference' else old
        incoming=prior;oldprior=oldprior+old;prior=prior+addition
        same(c,row[out+'C0']['p'],oldprior)
        background=common.restore_common_source(row[background_key]['original_normalized_history_C0_enclosures']['p'],coordinates)
        path.append(dict(label=label,actual_inherited_pressure_correction=incoming.record(),
            accepted_original_pressure_integral=old.record(),selected_original_pressure_integral=addition.record(),
            actual_right_pressure_correction=prior.record(),original_right_pressure_background=background.record(),
            actual_right_pressure_own_history=(background+prior).record(),
            exact_original_pressure_decay=1,original_true_width_retained=width.record(),
            only_inner_reference_integral_cover_changed=label=='inner_reference'))
    same(c,route['actual_original_inlet_to_O2_inlet_correction_C0']['p'],oldprior)
    values=dict(original,p=prior)
    return values,path


def pressure_replay(c,payload,values,record,archive):
    result=replay.replay(c,payload,record,archive)
    coordinates,old,jets,P0,P0Z=replay.inputs(c,payload)
    values={k:common.restore_common_source(encode(v.record()),coordinates) for k,v in values.items()}
    increment=coordinates.scalar(iv(c,record['five_original_C0_integral_contributions']['p']))
    original=coordinates.scalar(iv(c,record['source_defined_original_histories_at_y1']['p']))
    pressure=values['p']+increment;own=original+pressure
    result['actual_refined_upstream_C0']['p']=values['p'].record()
    result['propagated_actual_correction_C0']['p']=pressure.record()
    result['actual_original_plus_propagated_correction_C0']['p']=own.record()
    result['actual_absolute_pressure']=(P0+own).record()
    result['actual_upstream_binding'].update(manifest=upstream.NAME,manifest_sha256=sha(upstream.NAME),
        source_route_retained_with_pressure_integral_cover_refinement=True)
    result['weighted_pressure_cover_binding']=dict(manifest=NAME,source=Path(__file__).name,source_sha256=sha(Path(__file__).name),
        original_source_archive=archive['filename'],original_source_archive_sha256=sha(archive['filename']),
        report_hash_bound_by_separate_receipt_no_circular_self_hash=True)
    result['only_pressure_C0_integral_and_its_actual_serial_outputs_refined']=True
    result['pressure_C0_log_absolute_upper']=upstream.baseline.magnitude_log(pressure)
    return result


def same_Z(left,right):
    return tuple(replay.Fraction(v) for v in left)==tuple(replay.Fraction(v) for v in right)


def run():
    begin=time.monotonic();manifest=json.loads((HERE/upstream.NAME).read_bytes())
    hashes={};common.attach_receipt(hashes,upstream,manifest['source_family'])
    identity=source_identity();rows=[];outputs=[]
    originals=json.loads((HERE/common.NAME).read_bytes())
    c=MPIntervalContext();c.dps=240
    with mp.workdps(300):
        for archive in manifest['genuine_active_kappa_upstream_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            if hashlib.sha256(raw).hexdigest()!=archive['lossless_original_json_sha256']:
                raise ValueError('Accepted full original source archive changed')
            payload=json.loads(raw);coordinates,refined,proof=pressure_bound(c,payload)
            values,path=pressure_path(c,payload,refined)
            if not proof['strict_upper_reduction']:raise ArithmeticError('No pressure integral bound improvement')
            rows.append(dict(exact_Z_range=payload['exact_Z_range'],original_source_archive=archive,
                pressure_integral_proof=proof,actual_twelve_chart_pressure_path=path,
                actual_refined_upstream_C0={k:v.record() for k,v in values.items()},
                original_Z_covers_backgrounds_and_P0_retained=True))
            for record in originals['original_complete_same_N_source_incoming_O2_integral_refinements']:
                if same_Z(record['exact_Z_range'],payload['exact_Z_range']):
                    outputs.append(pressure_replay(c,payload,values,record,archive))
            print('Original E-squared weighted pressure and actual memory:',payload['exact_Z_range'],flush=True)
    if len(rows)!=2 or len(outputs)!=4:raise ValueError('Both exact axial tiles and all four accepted O2 records required')
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=manifest['source_family'],candidate_N=N,
        original_weighted_pressure_theorem=identity,actual_weighted_pressure_tiles=rows,
        genuine_original_O2_weighted_pressure_replays=outputs,input_hashes=hashes,execution_seconds=time.monotonic()-begin,
        original_other_four_C0_and_all_five_Z_covers_backgrounds_P0_unchanged=True,
        accepted_O2_integrals_reused_without_recomputation=True,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False),
        scope='Same original source/N1024 and strict-sign axial tiles. Exact inner_reference E_y=E/10 removes width times max pressure overestimate; complete rate-zero memory retained. Not moment closure, repaired field, stress, scale recursion, or NS residual.')
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
