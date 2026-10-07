"""Actual whole input coverage and independently evaluated finite scale fixtures."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_uniform_inputs as source
from lei_ren_part1_paper_compliant_current_generic_shear_loop import GenericLoopScales

packets=source.packets


def scale_fixtures():
    """Only modest artificial constants are exponentiated for independent comparison."""
    count=0
    for amin,m,boundary,t0,p1,p2 in (
        ('.8','1','.019','1','8','2'),('.1','.01','.004','3','12','5'),
        ('2','1','.1','0','5','0'),('.005','2','.0002','7','100','20')):
        expected=GenericLoopScales(a_min=amin,margin_min=m,boundary_kappa_excess_min=boundary,
            t0_abs_max=t0,p1_abs_max=p1,p2_abs_max=p2,dps=100)
        c=MPIntervalContext();c.dps=110;ep=packets.recovery.endpoints
        def cap(value):
            exact=c.mpf(value)
            return dict(exact_zero=value=='0',log_absolute_upper=None if value=='0' else c.ln(exact))
        previous=dict(positive_noncore_source_denominator_theorem=dict(source_charts={
            'fixture':dict(log_actual_a_positive_lower=c.ln(c.mpf(amin)))}),
            current_original_source_log_bound_charts={'fixture':dict(admitted_original_quotient_log_norms={
                key:{'y0_Z0':cap(value)} for key,value in zip(('t0','p1','p2'),(t0,p1,p2))})})
        fake=SimpleNamespace(ctx=c,previous=previous,O3=dict(original_O3_quotient_log_norms={}),
            geometry=dict(both_boundary_kappa_excess_log_lower=c.ln(c.mpf(boundary)),
                allowed_eta_log_upper_for_automatic_constant_edges=c.ln(c.mpf(boundary)/2)),scale_theorem={})
        result=source.CurrentUniformInputs.log_scales(fake,dict(uniform_log_H0_minus2_positive_lower=c.ln(c.mpf(m))))
        for key in ('q_star','B_star','J_star'):
            value=c.exp(result['logarithmic_conservative_upper_constants'][key])
            # Independently evaluated original formulas are slightly rounded.
            if ep(value)[1]<getattr(expected,key)*(1-expected.ctx.mpf('1e-90')):
                raise ArithmeticError('Conservative logarithmic scale below independent formula: '+key)
            count+=1
        eta=c.exp(result['selected_positive_eta_log'])
        if ep(eta)[0]>expected.eta*(1+expected.ctx.mpf('1e-90')):
            raise ArithmeticError('Logarithmic eta exceeds original independent recipe')
        count+=1
        actualB=c.exp(result['logarithmic_conservative_upper_constants']['B_star'])
        actualJ=c.exp(result['logarithmic_conservative_upper_constants']['J_star'])
        if ep(actualJ-c.mpf(p1)*actualB-c.mpf(p2))[1]<0:
            raise ArithmeticError('Actual selected logarithmic B/J relation failed')
        count+=1
    return count


def source_guards(field):
    calls=[];original=field.rows
    def test(stem,path,value):
        bad=copy.deepcopy(original);row=bad[stem]
        for key in path[:-1]:row=row[key]
        row[path[-1]]=value;field.rows=bad
        try:field.margins()
        finally:field.rows=original
    calls.append(lambda:test('current_reshape_relaxed_inputs',
        ('current_whole_R110_Rz_relaxed_input_and_bounds','whole_path_H0_minus2_strictly_above1'),False))
    calls.append(lambda:test('current_O2_modified_taper_cone',
        ('original_whole_closed_O2_taper_baseline_bounds','same_original_O2_pressure_energy_moments_and_all_stress_sectors'),False))
    calls.append(lambda:test('current_O3_transition_direction',
        ('current_whole_variable_transition_bounds','full_directional_expression_upper'),2))
    for call in calls:
        try:call()
        except (ValueError,ArithmeticError):pass
        else:raise ArithmeticError('Invalid original source/margin admitted')
    return len(calls)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed whole original input prerequisite: '+name)
    field=source.CurrentUniformInputs();c=field.ctx;ep=packets.recovery.endpoints
    margins=field.margins();scales=field.log_scales(margins)
    for key,actual in (('whole_actual_original_generic_input_margin',margins),
        ('current_actual_logarithmic_loop_scales',scales),
        ('exact_full_signed_source_normalization_theorem',field.theorem)):
        if packets.encode(actual)!=manifest[key]:raise ValueError('Current whole original input record differs: '+key)
    charts=set(scales['whole_source_chart_inventory']);groups=margins['chart_groups'];covered=set()
    for record in groups.values():
        covered.update(record['source_charts'])
        if ep(margins['uniform_log_H0_minus2_positive_lower']-record['log_H0_minus2_positive_lower'])[1]>0:
            raise ArithmeticError('Uniform log margin exceeds a current source lower')
    if charts!=covered or len(charts)!=17:raise ArithmeticError('Whole chosen noncore loop source coverage incomplete')
    first=groups['O2_buffer_first']['original_source_proof'];last=groups['O2_buffer_last']['original_source_proof']
    if tuple(first['original_selector_domain'])!=(0,9) or tuple(last['original_selector_domain'])!=(9,11):
        raise ArithmeticError('Full original buffer must cover both actual source pieces')
    if tuple(first['shared_offset_domain'])!=(-11,-2) or tuple(last['shared_offset_domain'])!=(-2,0):
        raise ArithmeticError('Original pure buffer source selector rebasing failed')
    for value in margins['source_margin_attachment_conditions'].values():
        if ep(value)[0]<=0:raise ArithmeticError('Actual whole source attachment gap lost')
    for limit in scales['eta_required_log_upper_constraints'].values():
        if ep(scales['selected_positive_eta_log']-limit)[1]>0:
            raise ArithmeticError('Selected formal eta exceeds a required actual source constraint')
    if not scales['eta_and_two_plus_eta_stay_formal_separate_sources'] or scales['old_scalar_fixture_constants_not_consumed'] is not True:
        raise ArithmeticError('Actual eta/source bounds replaced by scalar fixture')
    if any(manifest.get(k) for k in source.OPEN) or manifest['phase_held_loop_primitive_derivative_bounds_certified']:
        raise ArithmeticError('Whole source constants promoted missing physical loop stages')
    result=dict(all_passed=True,source_family=field.family,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        whole_modification_original_relaxed_input_certified=True,
        current_whole_source_loop_scales_in_log_form_certified=True,
        whole_original_source_charts=17,full_O2_buffer_selector_intervals=[[0,9],[9,11]],
        full_signed_H0_D_Q_normalization_identities=len(field.theorem['exact_identities']),
        original_Section11_scale_recipe_AST_bindings=len(field.scale_theorem['original_scale_recipe_AST_bindings']),
        actual_source_margin_attachment_conditions=len(margins['source_margin_attachment_conditions']),
        selected_positive_eta_log_constraints=len(scales['eta_required_log_upper_constraints']),
        independent_modest_original_scale_formula_comparisons=scale_fixtures(),
        invalid_original_source_margin_guards=source_guards(field),
        full_energy_absolute_pressure_and_post_radial_theta_denominator_retained=True,
        actual_eta_amplitude_radius_width_or_inverse_not_materialized=True,
        no_new_global_or_modified_strict_region_claimed=True,ancestor_constructors_called=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),
            Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole actual original generic input PASS:17 charts, full11-unit buffer, current H0 margin/log scales; independent fixtures PASS',flush=True)
    return result


if __name__=='__main__':run()
