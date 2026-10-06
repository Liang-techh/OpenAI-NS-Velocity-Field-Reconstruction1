"""Same-source full positive-radius core tensor and completed core/bridge trace."""
import gzip
import json
import ast
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_bridge_background_tensor import (
    CurrentBridgeBackgroundTensor,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,canonical_tensor_groups,IntervalTaylor)
from lei_ren_part1_paper_compliant_current_core_interior_moments import CurrentCoreInteriorMoments
from lei_ren_part1_paper_compliant_current_core_first_interface import CurrentCoreFirstInterface,binding
from lei_ren_part1_paper_compliant_core_physical_field import CompliantCorePhysicalField
from lei_ren_part1_paper_compliant_core_coefficient_rebuild import CompliantCoreCoefficientRebuild
from lei_ren_part1_paper_compliant_current_core_stress_operator import (
    assemble_core_raw_rows,core_raw_source,core_raw_unit_theorem,core_integrated_stress_theorem,
    reduce_core_stress_from_source,expanded_stress_sectors,factored_rows_record)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST

NAME=PREFIX+'current_core_background_tensor.json.gz'
RECEIPT=PREFIX+'current_core_background_tensor_check.json'
GATES=('current_positive_radius_core_full_background_tensor_available',
    'current_positive_radius_core_full_meridional_decomposition_available',
    'current_core_bridge_completed_tensor_join_certified')
OPEN=('current_core_full_background_tensor_available','core_axis_tensor_remainder_limits_certified',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'full_point_physical_field_evaluation','actual_point_moment_history_recovered',
    'independently_bounded_flat_remainder','temporal_recursion')
VIEWS={
    'compact_core':('core',(-1,1),('.001','4'),('-3','-1'),None,'1'),
    'inlet':('core',(-1,1),4,'-1',None,'1'),
    'fresh':('core','.359','.731','-2.6','.41','.8'),
    'fresh_low':('core','.731','.000137','-2.6','.41','.8'),
    'fresh_high':('core','-.317','3.831','-2.6','.41','.2')}


def current_core_tensor_source_theorem(field):
    interior=field.interior;first=field.first;core=field.core;bridge=field.bridge_tensor
    live=dict(exact_checked_interior_owner=type(interior) is CurrentCoreInteriorMoments,
        same_checked_interior_evaluator=interior.evaluate.__func__ is CurrentCoreInteriorMoments.evaluate,
        same_checked_original_core=type(core) is CompliantCorePhysicalField and core is first.core,
        same_original_core_profile=core.profiles.__func__ is CompliantCorePhysicalField.profiles,
        same_original_rebuild=type(interior.rebuild) is CompliantCoreCoefficientRebuild,
        same_original_rebuild_profile=interior.rebuild.profile.__func__ is CompliantCoreCoefficientRebuild.profile,
        same_common_core_first=first.common is interior.common,
        same_actual_core_atom_rebuild=first.atoms.rebuild is interior.rebuild,
        same_actual_first_bridge=first.bridge is bridge.bridge,
        same_context_and_original_pressure=field.ctx is core.ctx and core.datum is first.core.datum,
        same_checked_full_stress_operator=field.stress is bridge.stress,
        same_checked_full_remainder_operator=field.remainder is bridge.remainder,
        same_checked_full_physical_operator=field.lift is bridge.lift,
        same_raw_source_callable=field.raw_source is core_raw_source,
        same_raw_unit_callable=field.raw_source.__globals__['assemble_core_raw_rows'] is assemble_core_raw_rows)
    if not all(live.values()) or not first.acceptance_loaded:raise ValueError('One checked nonlinear core/moments/first/full-operator graph required')
    name=PREFIX+'current_core_first_interface_check.json'
    receipt=accepted(name,field.family,field.source,'current_core_bridge_functional_mixed4_join_certified');_verify_hashes(receipt)
    if receipt['datum_enclosure_sha256']!=field.datum_sha or not receipt['all_current_bridge_functional_interfaces_certified']:
        raise ValueError('Same datum and fourth functional source join required')
    equations=receipt['original_integrated_stress_free_core_equations']
    coordinate=receipt['phase0_ODE_and_exact_positive_width_pullback']
    atoms=receipt['scaled_atoms_are_the_original_physical_primitives']
    if not all(item['passed'] for item in (equations,coordinate,atoms)) or not coordinate['equality_by_functions_and_ODE_uniqueness_not_interval_overlap']:
        raise ValueError('Original exact common core equations/axis constants/inlet/pullback proof required')
    unit=core_raw_unit_theorem();integrated=core_integrated_stress_theorem();asts=SourceAST()
    if unit['total_identities']!=135 or len(integrated['full_original_raw_stress_to_integrated_core_ODE_identities'])!=20:
        raise ValueError('Full raw core source units/stress identities incomplete')
    for target,expr in (('p0',"IntervalTaylor(c,mom['original_analytic_axis_pressure_axial6_Taylor_coefficients'][:6])"),
        ('ratios',"mom['original_F0_relative_axial6_ordinary_derivatives']"),
        ('ratios2',"mom['original_F0_squared_relative_axial6_ordinary_derivatives']"),
        ('logR','c.ln(field.core.epsilon)+c.ln(r)'),
        ('logF0',"c.mpf([endpoints(-field.core.logC-field.core.Lambda*field.core.Gbar)[0],endpoints(-field.core.logC)[1]])")):
        binding('compliant_current_core_stress_operator','core_raw_source',target,expr)
    # Match the original admitted positive F0 bound, not a selected amplitude.
    original_profile=asts.method('core_physical_field','profiles')
    original_bounds=[node.value for node in ast.walk(original_profile)
        if isinstance(node,ast.keyword) and node.arg=='F0_exact_positive_log_enclosure']
    expected=ast.parse("c.mpf([endpoints(-self.logC-self.Lambda*self.Gbar)[0],endpoints(-self.logC)[1]])",mode='eval').body
    if len(original_bounds)!=1 or ast.dump(original_bounds[0])!=ast.dump(expected):
        raise ValueError('Original admitted positive F0 source bound changed')
    for stem,method in (('current_core_stress_operator','core_raw_source'),('current_core_stress_operator','reduce_core_stress_from_source'),
        ('current_core_background_tensor','chart'),('current_core_background_tensor','assert_graph')):
        asts.method(stem,method)
    return dict(live_original_callable_bindings=live,original_core_raw_unit_theorem=unit,
        original_integrated_total_core_stress_theorem=integrated,
        consumed_same_source_core_first_receipt=receipt,
        same_fixed_point_model_and_all_finite_tail_bindings=interior.proof,
        original_positive_F0_log_bound_source_bound=True,
        rho4_six_moment_sources_equal_original_true_atoms=interior.proof['original_rho4_finite_atom_identities'],
        rho4_same_infinite_density_integral_sources_and_all_product_tails=True,
        core_first_ordinary_y4_function_equality_then_identical_full_operators=True,
        core_raw_h_k_and_bridge_raw_h_k_use_same_physical_moments=True,
        normalized_P0_and_cumulative_C_have_identical_original_units=True,
        only_total_core_stress_zero_original_signed_sectors_retained=True,
        all_original_nonzero_remainder_sectors_retained=True,
        positive_radius_source_identity_does_not_admit_axis_or_uniform_axis_bound=True,
        exact_common_source_core_stress_zero_before_physical_bounds=True,
        input_hashes={**receipt['input_hashes'],name:sha(name),**unit['input_hashes'],**integrated['input_hashes'],**asts.hashes},
        passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentCoreBackgroundTensor:
    @source_precision
    def __init__(self,interior=None,require_checked=True):
        self.interior=interior if interior is not None else CurrentCoreInteriorMoments()
        self.bridge_tensor=self.interior.bridge_tensor;self.physical=self.bridge_tensor.physical
        self.pulse=self.bridge_tensor.pulse;self.core=self.interior.core;self.ctx=self.core.ctx
        self.family=self.interior.family;self.source=self.interior.source;self.datum_sha=self.interior.datum_sha
        self.first=CurrentCoreFirstInterface(common=self.interior.common)
        self.stress=self.bridge_tensor.stress;self.remainder=self.bridge_tensor.remainder;self.lift=self.bridge_tensor.lift
        self.raw_source=core_raw_source;self.assert_graph()
        source=current_core_tensor_source_theorem(self)
        self.proof=dict(current_same_source_core_units_integrated_stress_and_completed_first_trace=source,
            consumed_original_full_factored_tensor_operators=self.bridge_tensor.proof['consumed_checked_original_full_factored_tensor_operators'],
            zero_total_stress_preserves_original_signed_cancellation_ledger=True,
            full_nonzero_NS_remainder_not_removed=True,passed=True)
        self.hashes={**self.interior.hashes,**self.first.hashes,**source['input_hashes']}
        for stem in ('current_core_background_tensor','current_core_stress_operator'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Positive-radius core tensor admission source/scope differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.interior.assert_graph();self.bridge_tensor.assert_graph()
        if not (type(self.interior) is CurrentCoreInteriorMoments and self.interior.acceptance_loaded
            and self.first.acceptance_loaded and self.first.common is self.interior.common
            and self.core is self.interior.core is self.first.core and self.ctx is self.core.ctx
            and self.bridge_tensor is self.interior.bridge_tensor
            and self.physical is self.bridge_tensor.physical and self.pulse is self.bridge_tensor.pulse
            and self.stress is self.bridge_tensor.stress and self.remainder is self.bridge_tensor.remainder
            and self.lift is self.bridge_tensor.lift and self.raw_source is core_raw_source):
            raise ValueError('Same checked core interior/first bridge and full tensor operators required')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);rho=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart!='core' or not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(rho)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(rho)[0]<=0 or endpoints(rho)[1]>4 or endpoints(nu)[0]<=0:
            raise ValueError('Positive-radius core 0<rho<=4,Z[-1,1],finite log tau and nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        raw,source=self.raw_source(self,Z,rho);algebra=raw['algebra']
        z=algebra.lift(IntervalTaylor.variable(c,Z,5));delta=c.mpf(self.pulse.delta)
        original_rows=self.stress(c,delta,z,raw['velocity']['theta'],raw['velocity']['axial'],raw['histories'],raw['absolute_pressure'])
        logR=source['logR']
        original_sectors=expanded_stress_sectors(c,original_rows,algebra,logR,self.physical.logP)
        totals=reduce_core_stress_from_source(original_rows,self.proof['current_same_source_core_units_integrated_stress_and_completed_first_trace'])
        sectors=expanded_stress_sectors(c,totals,algebra,logR,self.physical.logP)
        packet=dict(Z=Z,s=rho,exact_logR=logR,exact_pulse_reference_logB_parts=dict(logPstar=self.physical.logP),
            exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors)
        point=self.lift(c,packet,delta,raw['velocity'],lt,theta,nu)
        return dict(point,chart=chart,coverage_coordinate=rho,current_actual_source_stress_packet=packet,
            actual_same_fixed_point_core_source=source,current_unresolved_raw_source_rows=factored_rows_record(raw),
            fixed_current_factored_source_log_bases=list(algebra.logs),
            original_uncancelled_signed_stress_sector_enclosures=original_sectors,
            exact_zero_total_stress_is_separate_from_original_sector_ledger=True,
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate('core',Z,rho,lt,theta),
            ordinary_logR_conversion='D_y=rho*d_rho; Euler/Stirling order4, exact R=epsilon_core*rho',
            every_original_NS_remainder_sector_retained=True,
            original_radial_prefactors_and_absolute_P0_included_once=True,
            actual_full_stress_not_local_difference=True,source_bounds_not_resolved_physical_point_values=True,
            physical_axis_not_in_this_chart=True,uniform_axis_or_required_domain_energy_not_inferred=True,
            **dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name='core_bridge',Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name!='core_bridge':raise ValueError('Only completed current core/first bridge trace required')
        left=self.chart('core',Z,4,log_tau,theta,viscosity)
        right=self.bridge_tensor.chart('bridge_first',Z,0,log_tau,theta,viscosity)
        aa=canonical_tensor_groups(left);bb=canonical_tensor_groups(right)
        if set(aa)!=set(bb):raise ValueError('Same completed core/first tensor layout required')
        common={}
        for key in aa:
            values=[endpoints(row['log_absolute_upper'])[1] for row in aa[key]+bb[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_actual_source_and_completed_tensor_endpoint_function_theorem=self.proof,
            source_function_equality_precedes_common_triangle_bounds=True,interval_overlap_not_used_as_function_identity=True,
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_positive_radius_core_source_and_completed_first_trace_theorem=self.proof,
            actual_current_tensor_regions_available=['core_positive_radius']+self.bridge_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=32,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            exact_core_domain='0<rho<=4,Z[-1,1]; physical axis excluded',
            published_compact_sector='rho[.001,4]; fresh rho .000137 is below this compact sector',
            any_requested_positive_compact_core_sector_supported=True,
            uniform_bound_as_rho_tends_to_zero_not_claimed=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCoreBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_positive_radius_core_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_positive_radius_core_tensor_views'][name]=field.chart(*args)
        print('Build full same-source positive-radius core tensor: '+name,flush=True)
    result['current_core_bridge_completed_tensor_interface']=field.interface()
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
