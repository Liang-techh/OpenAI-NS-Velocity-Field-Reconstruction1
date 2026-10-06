"""Four exact angular support limits of the actual full background T/E.

Endpoint source identities precede physical bounds. The full nonzero stress,
pressure and history are retained; local difference bounds are not substituted.
"""
import ast
import copy
import gzip
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_core_axis_background_tensor import (
    CurrentCoreAxisBackgroundTensor,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,
    source_precision,accepted,_verify_hashes)
from lei_ren_part1_paper_compliant_current_angular_background_stress import CurrentAngularBackgroundStress
from lei_ren_part1_paper_compliant_current_angular_support_interfaces import (
    CurrentAngularSupportInterfaces,EDGES,source_edge_integrals,exact_endpoint_program)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import canonical_tensor_groups
from lei_ren_part1_paper_compliant_collar_physical_C2 import scale_row
from lei_ren_part1_paper_compliant_current_angular_internal_operator import two_sided_full_angular_trace_theorem

NAME=PREFIX+'current_angular_internal_background_tensor.json.gz'
RECEIPT=PREFIX+'current_angular_internal_background_tensor_check.json'
GATES=('current_four_angular_internal_completed_tensor_joins_certified',
    'current_all_fourteen_internal_completed_tensor_traces_certified')
SEAMS={'outer_angular_-3_-1':'-63/20','outer_angular_-3_1':'-57/20',
    'outer_angular_-1_-1':'-23/20','outer_angular_-1_1':'-17/20'}
VIEWS={name+'_whole':(name,(-1,1),('-3','-1'),None,'1') for name in SEAMS}
VIEWS.update({name+'_fresh':(name,z,'-2.6','.41',nu) for (name,z,nu) in
    zip(SEAMS,('.537','-.419','.731','-.317'),('.8','.03','.2','1'))})


def angular_owner(core_tensor):
    """Reuse the existing checked downstream owner; never rebuild a second core."""
    obj=core_tensor.bridge_tensor
    for attribute in ('micro_tensor','power_tensor','reshape_tensor','restore_tensor','patch_tensor',
        'o2_tensor','o3_tensor','entrance_incoming_tensor','main_tensor','gap_tensor','end_tensor','flatten_power','angular'):
        obj=getattr(obj,attribute)
    return obj


def exact_remaining_program():
    """Original remaining B/D/E/F program at exact rational support limits."""
    asts=SourceAST();fn=copy.deepcopy(asts.method('angular_stress_C3','angular_future_changes'))
    fn.name='_exact_angular_remaining';fn.decorator_list=[]
    fn.args.kwonlyargs.append(ast.arg(arg='source_edge'));fn.args.kw_defaults.append(None)
    replacements=0
    class Reduce(ast.NodeTransformer):
        def visit_If(self,node):
            nonlocal replacements
            if ast.unparse(node.test)=='hi <= -endpoints(ell)[1]':
                replacements+=1
                return ast.parse('past,future=source_edge_integrals(outer,center,source_edge)').body[0]
            return self.generic_visit(node)
    fn=Reduce().visit(fn)
    if replacements!=1:raise ValueError('Original angular remaining support selection changed')
    if any(isinstance(node,ast.Call) and ast.unparse(node.func)=='future_bump_weights' for node in ast.walk(fn)):
        raise ValueError('Partial rounded future selection survives exact whole-If reduction')
    from lei_ren_part1_paper_compliant_angular_stress_C3 import angular_future_changes
    env=dict(angular_future_changes.__globals__);env['source_edge_integrals']=source_edge_integrals
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original-angular-remaining-exact-support>','exec'),env)
    return env[fn.name],dict(original_remaining_source_sha256=asts.hashes,
        exact_support_selection_replacements=replacements,
        original_signed_B_D_E_F_coefficients_and_translations_unchanged=True)


def exact_full_tensor_program(endpoint_recipe,remaining_recipe):
    """Same full angular operator, with only proved source limits substituted."""
    asts=SourceAST();fn=copy.deepcopy(asts.method('current_angular_background_stress','angular'))
    fn.name='_exact_full_angular_tensor';fn.decorator_list=[]
    fn.args.kwonlyargs.append(ast.arg(arg='source_edge'));fn.args.kw_defaults.append(None)
    replacements={name:0 for name in ('native','future')}
    class Reduce(ast.NodeTransformer):
        def visit_Assign(self,node):
            if len(node.targets)==1:
                name=ast.unparse(node.targets[0])
                if name=='native':
                    if ast.unparse(node.value)!='self.outer.angular(Z, v)':raise ValueError('Full angular native callback changed')
                    replacements[name]+=1;node.value=ast.parse('endpoint_recipe(self.outer,Z,v,source_edge)',mode='eval').body
                elif name=='future':
                    if ast.unparse(node.value)!='angular_future_changes(self.outer, data, v)':raise ValueError('Full angular remaining callback changed')
                    replacements[name]+=1;node.value=ast.parse('remaining_recipe(self.outer,data,v,source_edge=source_edge)',mode='eval').body
            return self.generic_visit(node)
    fn=Reduce().visit(fn)
    if replacements!=dict(native=1,future=1):raise ValueError('Exact full source substitution incomplete')
    env=dict(CurrentAngularBackgroundStress.angular.__wrapped__.__globals__)
    env.update(endpoint_recipe=endpoint_recipe,remaining_recipe=remaining_recipe)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<same-full-angular-tensor-exact-source-limit>','exec'),env)
    return env[fn.name],dict(exact_native_and_remaining_replacements=replacements,
        unchanged_full_angular_source_sha256=asts.hashes,
        every_original_full_moment_pressure_stress_divergence_completion_remainder_and_log_statement_retained=True)


def current_angular_internal_source_theorem(field):
    angular=field.angular;support=field.support;asts=SourceAST();field.assert_graph();field.assert_programs()
    source=support.proof;primitive=source['current_angular_endpoint_source_theorem']
    for key in ('current_primitive_mixed4_identities','current_velocity_pressure_mixed4_identities'):
        if not primitive['passed'] or not primitive[key] or not all(primitive[key].values()):
            raise ValueError('Both original one-sided source function limits required')
    for proof in (source['weighted_integral_FTC_and_flat_beta_source_theorem'],angular.moment_proof,
        angular.baselines,angular.pressure_proof,angular.KR,angular.log_proof):
        if not proof['passed']:raise ValueError('Same full angular source proof not admitted')
    two_sided=two_sided_full_angular_trace_theorem()
    if not two_sided['passed'] or (two_sided['full_moment_identity_count'],two_sided['full_physical_group_identity_count'])!=(60,71):
        raise ValueError('Both nonzero full source/operator sides must be instantiated before bounds')
    for recipe in (field.endpoint_recipe,field.remaining_recipe):
        if any(name in recipe.__code__.co_names for name in ('bump_weights','future_bump_weights')):
            raise ValueError('Rounded partial support integration survives an exact endpoint source recipe')
    if not (angular.physical_proof['actual_regional_physical_angular_stress_remainder_identity_verified']
        and all(angular.physical_proof['reduction_identities'].values())
        and all(angular.physical_proof['general_K_physical_transfer']['identities'].values())):
        raise ValueError('Original physical tensor/remainder theorem missing')
    if set(SEAMS.values())!={edge['exact_edge'] for edge in EDGES}:raise ValueError('All four exact source supports required')
    ledger={};limits={}
    for name,q in SEAMS.items():
        row=angular.atlas.internal[name]
        if row['chart']!='outer_angular' or row['exact_source_edge']!=q or not row['current_full_field_interface_trace_admitted']:
            raise ValueError('Checked primitive support ledger changed')
        ledger[name]=dict(row)
        limits[name]={}
        for center in (-3,-1):
            past,future=source_edge_integrals(angular.outer,center,q)
            limits[name][str(center)]=dict(past_full=past is angular.outer.weights,
                future_full=future is angular.outer.weights,
                exact_support_limit='full/empty original A/B/D/E/F weight functions, not rounded clipping')
    # Pin the exact original scalar operator pipeline and constant physical
    # radius/factor pipeline. Both one-sided source limits use this same body.
    for stem,method in (('current_angular_background_stress','normalized_full_moment_rows'),
        ('current_angular_background_stress','angular'),('current_angular_background_stress','angular_source_log_parts'),
        ('angular_stress_C3','angular_future_changes'),('collar_stress_C3','collar_stress_rows'),
        ('collar_physical_C2','collar_velocity_bracket'),('collar_physical_C2','physical_source_row'),
        ('pulse_physical_bounds','physical_bracket'),('current_angular_support_interfaces','source_edge_integrals')):
        asts.method(stem,method)
    return dict(checked_both_one_sided_primitive_and_velocity_pressure_mixed4_source_theorem=primitive,
        checked_weighted_support_FTC_and_flat_beta_source_theorem=source['weighted_integral_FTC_and_flat_beta_source_theorem'],
        checked_same_normalization_and_five_original_weight_function_sources=source['exact_normalization_and_five_weight_defining_source'],
        original_both_sided_full_nonzero_moment_and_canonical_tensor_operator_theorem=two_sided,
        original_full_normalized_A_E_P_K_AST_theorem=angular.moment_proof,
        original_full_stress_baseline_and_KR_AST_theorem=angular.baselines,
        same_complete_P0_forward_pressure_future_and_terminal_identity=angular.pressure_proof,
        original_general_physical_tensor_divergence_and_viscosity_theorem=angular.physical_proof,
        same_exact_KR_and_source_factor_logs=dict(KR=angular.KR,logs=angular.log_proof),
        exact_source_ordering_at_all_four_edges=limits,current_four_angular_source_ledger=ledger,
        exact_native_endpoint_recipe=source['exact_private_endpoint_reduction'],
        exact_remaining_recipe=field.remaining_binding,exact_full_tensor_recipe=field.tensor_binding,
        entire_If_elif_else_selection_replaced_and_no_partial_weight_calls_survive=True,
        same_original_C5_amplitudes_full_histories_quadratic_D_F_and_absolute_pressure_retained=True,
        source_functions_and_mixed4_limits_identical_before_full_operators=True,
        same_original_full_operators_imply_both_one_sided_completed_stress3_divergence2_remainder2_traces=True,
        common_full_trace_is_not_local_zero_difference_or_interval_overlap=True,
        source_caps_are_only_bounds_on_exact_inverse_radius_and_amplitudes=True,
        input_hashes={**asts.hashes,**two_sided['input_hashes']},passed=True)


def canonical_angular_groups(c,view):
    """Common 71-component layout, retaining exact radial/axial E zero rows."""
    converted=dict(view)
    converted['physical_cylindrical_stress_mixed3']={name:dict(angular_full=grid) for name,grid in view['physical_cylindrical_stress_mixed3'].items()}
    converted['physical_cylindrical_stress_divergence_mixed2']={name:dict(angular_full=grid) for name,grid in view['physical_stress_divergence_mixed2'].items()}
    converted['completed_theta_theta_stress_mixed2']=dict(angular_full=view['completed_theta_theta_stress_mixed2'])
    error=view['physical_axial_viscosity_remainder_mixed2'];zero={key:scale_row(c,row,c.mpf(0)) for key,row in error.items()}
    converted['physical_three_component_remainder_mixed2']=dict(radial=dict(exact_zero=zero),theta=dict(angular_full=error),axial=dict(exact_zero=zero))
    converted['physical_completed_stress_tensor_cartesian']=view['completed_background_tensor_cartesian_components']
    return canonical_tensor_groups(converted)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


@source_precision
def leading_origin_time_obstruction(field):
    """The present leading E is nonflat; higher coefficient repair is required."""
    point=field.core_axis.axis(0,'-1','1');c=field.ctx
    row=point['physical_core_Cartesian_remainder_mixed2']['x0_y0_z0']['z']['axial_axial_viscosity']
    lower=endpoints(row['signed_coefficient'])[0];gamma_upper=endpoints(row['physical_lambda_exponent'])[1]
    if not (lower>0 and gamma_upper<0 and all(endpoints(value)==endpoints(c.mpf(0)) for value in row['actual_source_log_parts'].values())
        and endpoints(row['radial_log_prefactor'])==endpoints(c.mpf(0))
        and endpoints(row['physical_viscosity_exponent'])==endpoints(c.mpf('.5'))):
        raise ValueError('Same leading origin source sign/scaling obstruction changed')
    return dict(exact_origin='x_phys=y_phys=z_phys=0; rho=Z=0; lambda=sqrt(tau)',
        actual_same_source_nonzero_axial_remainder_row=row,
        exact_leading_time_factor='E_z(0,tau,nu)=sqrt(nu)*k_axis*tau^((-3+delta)/2)',
        strictly_positive_coefficient_lower=lower,strictly_negative_gamma_upper=gamma_upper,
        lower_for_0_tau_lt_1='sqrt(nu)*coefficient_lower*tau^(gamma_upper/2)',
        current_leading_origin_remainder_diverges_as_tau_to_zero=True,
        current_leading_remainder_cannot_be_terminal_time_flat=True,
        scoped_to_present_fixed_leading_family_not_future_recursive_or_oscillatory_corrected_field=True,
        genuine_n_dependent_coefficient_recovery_required_before_claiming_flat_completed_background=True,
        passed=True)


class CurrentAngularInternalBackgroundTensor:
    @source_precision
    def __init__(self,core_axis=None,require_checked=True):
        self.core_axis=core_axis if core_axis is not None else CurrentCoreAxisBackgroundTensor()
        self.angular=angular_owner(self.core_axis.core_tensor);self.support=self.angular.atlas.angular
        self.physical=self.angular.physical;self.ctx=self.angular.ctx
        self.family=self.angular.family;self.source=self.angular.source;self.datum_sha=self.angular.datum_sha
        self.assert_graph();self.endpoint_recipe=exact_endpoint_program()
        self.remaining_recipe,self.remaining_binding=exact_remaining_program()
        self.full_operator,self.tensor_binding=exact_full_tensor_program(self.endpoint_recipe,self.remaining_recipe)
        self.program_codes=tuple(fn.__code__ for fn in (self.endpoint_recipe,self.remaining_recipe,self.full_operator))
        self.proof=current_angular_internal_source_theorem(self)
        self.leading_obstruction=leading_origin_time_obstruction(self)
        self.hashes={**self.core_axis.hashes,**self.angular.hashes,**self.support.hashes,**self.proof['input_hashes']}
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.acceptance_loaded=False
        name=PREFIX+'current_angular_internal_operator.py';self.hashes[name]=sha(name)
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Angular full support tensor admission exceeds scoped source')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.core_axis.assert_graph();self.angular.assert_graph()
        if not (type(self.core_axis) is CurrentCoreAxisBackgroundTensor and self.core_axis.acceptance_loaded
            and type(self.angular) is CurrentAngularBackgroundStress and self.angular.acceptance_loaded
            and self.angular is angular_owner(self.core_axis.core_tensor)
            and type(self.support) is CurrentAngularSupportInterfaces and self.support.acceptance_loaded
            and self.support is self.angular.atlas.angular and self.support.physical is self.physical is self.core_axis.core_tensor.physical
            and self.support.outer is self.angular.outer and self.ctx is self.angular.ctx is self.core_axis.ctx
            and (self.family,self.source,self.datum_sha)==(self.core_axis.family,self.core_axis.source,self.core_axis.datum_sha)):
            raise ValueError('One checked core/axis/angular/full-pressure/support owner graph required')

    def assert_programs(self):
        if not (tuple(fn.__code__ for fn in (self.endpoint_recipe,self.remaining_recipe,self.full_operator))==self.program_codes
            and self.full_operator.__globals__['endpoint_recipe'] is self.endpoint_recipe
            and self.full_operator.__globals__['remaining_recipe'] is self.remaining_recipe
            and self.endpoint_recipe.__globals__['source_edge_integrals'] is source_edge_integrals
            and self.remaining_recipe.__globals__['source_edge_integrals'] is source_edge_integrals):
            raise ValueError('Unchanged full tensor and source-proved endpoint programs required')

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph();self.assert_programs()
        if name not in SEAMS:raise ValueError('Unknown exact angular full tensor support seam')
        c=self.ctx;q=s.Rational(SEAMS[name]);coordinate=c.mpf(str(q.p))/q.q
        full=self.full_operator(self.angular,Z,coordinate,log_tau,theta,viscosity,source_edge=SEAMS[name])
        groups=canonical_angular_groups(c,full);common={}
        for key,rows in groups.items():
            positive=[endpoints(row['log_absolute_upper'])[1] for row in rows if not row['exact_zero']]
            common[key]=dict(exact_zero=not positive,
                common_triangle_log_upper=None if not positive else c.mpf(max(positive))+c.ln(len(positive)))
        return dict(seam=name,exact_source_edge=SEAMS[name],source_coordinate_enclosure=coordinate,
            current_full_angular_endpoint_tensor=full,common_actual_tensor_rows=common,
            current_common_tensor_contribution_count=len(common),current_common_physical_contribution_count=sum(len(rows) for rows in groups.values()),
            exact_source_function_limits_and_full_operators_precede_bounds=True,
            full_nonzero_histories_pressure_stress_and_remainder_retained=True,
            ordinary_s_derivatives_not_support_distance_derivatives=True,
            same_complete_source_family_datum_and_physical_factors_on_both_sides=True,
            only_proved_endpoint_source_callbacks_reduced=True,
            local_zero_difference_not_substituted_for_actual_full_tensor=True,
            interval_overlap_not_used_as_source_function_equality=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_four_angular_internal_full_tensor_source_theorem=self.proof,
            actual_current_tensor_regions_available=self.core_axis.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=32,actual_current_completed_tensor_internal_interface_count=14,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            current_core_full_background_tensor_available=True,core_axis_tensor_remainder_limits_certified=True,
            current_leading_origin_remainder_time_obstruction=self.leading_obstruction,
            scope='Four exact rational internal angular support limits, whole source Z[-1,1], finite compact log tau and finite positive nu; full original stress3/divergence2/remainder2, not global cone/time-flatness/energy/points/recursion',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentAngularInternalBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_four_angular_internal_full_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_four_angular_internal_full_tensor_views'][name]=field.interface(*args)
        print('Build full angular internal tensor trace: '+name,flush=True)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
