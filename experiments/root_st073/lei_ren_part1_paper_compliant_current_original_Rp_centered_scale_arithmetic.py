"""Exact fixed binary-rational source origins and centered scale error.

L=q*logC+R is a proved identity of the unchanged original scale functions.
The exact selected logC singleton is not rounded into the error width of R.
No original physical exponential is expanded or replaced by a practical scale.
"""
from fractions import Fraction
import copy
import gzip
import json
import math
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rp_refined_pulse_coefficients as refined
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

accuracy, signed, correlated = refined.accuracy, refined.signed, refined.correlated
HERE, PREFIX, sha, ends = refined.HERE, refined.PREFIX, refined.sha, refined.ends
NAME = PREFIX+'current_original_Rp_centered_scale_arithmetic.json.gz'
RECEIPT = PREFIX+'current_original_Rp_centered_scale_arithmetic_check.json'
GATES = ('current_original_Rp_exact_binary_rational_logC_origin_installed',
         'current_original_Rp_same_source_finite_parameter_log_refinement_installed',
         'current_original_Rp_centered_physical_scale_width_budget_installed')
OPEN = refined.OPEN


def exact_scaled_binary_rational(singleton, coefficient):
    """Encode q*m*2**e exactly, without an exponent-sized integer shift."""
    if singleton[0] != singleton[1]:
        raise ValueError('An exact original source singleton is required')
    sign, mantissa, exponent, bits = singleton[0]
    coefficient = Fraction(coefficient)
    numerator = (-1 if sign else 1)*mantissa*coefficient.numerator
    denominator = coefficient.denominator
    if not numerator:
        return dict(numerator=0, denominator=1, binary_exponent=0)
    divisor = math.gcd(abs(numerator), denominator)
    numerator //= divisor
    denominator //= divisor
    while numerator % 2 == 0:
        numerator //= 2
        exponent += 1
    while denominator % 2 == 0:
        denominator //= 2
        exponent -= 1
    return dict(numerator=numerator, denominator=denominator, binary_exponent=exponent)


@source_precision
def centered_relative_budget(ctx, row, target, center):
    """Only the true variation of R enters width(exp(q*logC+R))."""
    result = accuracy.relative_error_budget(ctx,row,target)
    result.update(original_absolute_scale_log_function_retained=True,
        exact_fixed_origin_is_not_a_numeric_scale_approximation=True,
        centered_scale_error_identity='width(L)=width(R), L=q*exact_logC+R',
        exact_fixed_origin=center['exact_fixed_origin'],
        directed_residual_log_scale=center['directed_residual_log_scale'],
        centered_reference_log_width_bound=center['residual_log_width_bound'],
        direct_rounded_reference_log_width_bound=result.get('reference_log_width_bound'),
        centered_scale_proof=center['exact_centered_scale_identity_proof'])
    if row['exact_zero_enclosure'] or not result.get('sign'):
        return result
    common=row['common_scale_coefficient_enclosure']
    a,b=ends(common)
    # mpf abs rounds to the ambient scalar precision; retain the exact
    # directed interval endpoints before computing the ratio in ctx.
    magnitudes=[mp.make_mpf((0,*endpoint._mpf_[1:])) for endpoint in (a,b)]
    minimum,maximum=min(magnitudes),max(magnitudes)
    coefficient_log_ratio=ctx.ln(ctx.mpf(maximum)/ctx.mpf(minimum))
    ordinary_log_ratio=center['residual_log_width_bound']+coefficient_log_ratio
    eps=Fraction(target)
    allowed=ctx.ln(1+ctx.mpf(eps.numerator)/eps.denominator)
    result.update(reference_log_width_bound=center['residual_log_width_bound'],
        log_ordinary_magnitude_ratio_bound=ordinary_log_ratio,
        log_allowed_magnitude_ratio_bound=allowed,
        ordinary_numeric_relative_width_satisfied=ends(ordinary_log_ratio)[1]<=ends(allowed)[0],
        remaining_scale_error_comes_from_directed_residual_not_exact_origin=True)
    result.pop('ordinary_relative_width_bound',None)
    result.pop('ordinary_relative_width_upper_representation',None)
    if ends(ordinary_log_ratio)[1]<=signed.EXP_LOG_LIMIT:
        result['ordinary_relative_width_bound']=ctx.exp(ordinary_log_ratio)-1
    else:
        result['ordinary_relative_width_upper_representation']=dict(
            kind='expm1_of_directed_centered_log_ratio_bound',argument=ordinary_log_ratio,
            physical_exponential_not_materialized=True)
    result['ordinary_numeric_delivery_target_satisfied']=(result['ordinary_numeric_relative_width_satisfied']
        and row['ordinary_numeric_materialized'] and row.get('requested_relative_width_satisfied',False))
    if result['ordinary_numeric_relative_width_satisfied']:
        result.pop('ordinary_unresolved_reason',None)
    else:
        result['ordinary_unresolved_reason']='Residual source/arithmetic enclosure still exceeds ordinary width target; exact origin was removed from error, not from the field'
    return result


class CurrentOriginalRpCenteredScaleArithmetic:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not refined.CurrentOriginalRpRefinedPulseCoefficients or not before.acceptance_loaded:
            raise ValueError('Accepted refined current pulse/source owner required')
        self.before,self.amplitude=before,before.amplitude
        self.graph,self.family_record=before.graph,before.family_record
        self.ctx=MPIntervalContext();self.ctx.dps=400
        self.hashes=dict(before.hashes)
        for name in (refined.NAME,refined.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        self._readers={}
        self.pivot=self.amplitude.radius.parameters['logC']
        base=self.amplitude.reader()
        self.exact_singleton=base.bindings[self.pivot.node]._mpi_
        if self.exact_singleton[0]!=self.exact_singleton[1]:
            raise ValueError('Actual selected logC must stay an exact singleton')
        # Same elementary source graph, independent directed arithmetic.
        # Finite parameter exponentials are the existing whitelist only.
        pure=correlated.locator.CancelledGraphBounds(self.graph,self.ctx,{},base.allowed_exponentials)
        self.finite_bindings={}
        self.inclusion={}
        for name,ref in (('mu',self.amplitude.mu),('logP',self.amplitude.logP),
                         ('Tw',self.amplitude.Tw),('logU0',self.amplitude.logU0)):
            value=pure.at(ref)
            previous=base.bindings[ref.node]
            if not accuracy.pressure.amplitude.contains(previous,value):
                raise ArithmeticError('Fresh elementary source must be included in old '+name+' binding')
            self.finite_bindings[ref.node]=value
            self.inclusion[name]=dict(exact_source_function=ref.node,original_box=previous,
                fresh_box=value,original_graph_recomputed=True,included_in_original_box=True)
        U=self.ctx.exp(self.finite_bindings[self.amplitude.logU0.node])
        if not accuracy.pressure.amplitude.contains(base.bindings[self.amplitude.raw.U0.node],U):
            raise ArithmeticError('Fresh exact U0 must stay inside accepted U0 source binding')
        self.finite_bindings[self.amplitude.raw.U0.node]=U
        self.inclusion['U0']=dict(exact_source_function=self.amplitude.raw.U0.node,
            original_box=base.bindings[self.amplitude.raw.U0.node],fresh_box=U,
            original_graph_and_accepted_compact_identity_recomputed=True,included_in_original_box=True)
        self._snapshots={node:value._mpi_ for node,value in self.finite_bindings.items()}
        self._pivot_node=copy.deepcopy(self.graph.nodes[self.pivot.node])
        self.parameter_proof=self.graph.node('current_original_Rp_independent_finite_parameter_log_arithmetic',
            original_amplitude_identity=self.amplitude.proof.node,
            refined_source_nodes=sorted(self.finite_bindings),directed_decimal_precision=400,
            source_parameter_boxes_recomputed_from_original_elementary_definitions=True,
            no_original_parameter_or_radius_changed=True)
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family']!=self.family_record:
                raise ValueError('Centered original-scale receipt/scope differs')
            refined.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        base=self.amplitude.reader()
        checks=dict(accepted_current_refined_coefficient_source=self.before.acceptance_loaded,
            same_exact_original_pivot=self.graph.nodes[self.pivot.node]==self._pivot_node and
                base.bindings[self.pivot.node]._mpi_==self.exact_singleton,
            original_logC_is_exact_singleton=self.exact_singleton[0]==self.exact_singleton[1],
            unchanged_fresh_source_bindings=all(value._mpi_==self._snapshots[node] for node,value in self.finite_bindings.items()),
            same_source_alias_and_whitelist=set(self.finite_bindings)<=set(base.bindings),
            independent_bounded_arithmetic_precision=self.ctx.dps==400,
            original_physical_exponential_guard=signed.EXP_LOG_LIMIT==1000)
        if not all(checks.values()):raise ValueError('Centered original source differs: '+str(checks))
        return checks

    @source_precision
    def reader(self,delivery):
        self.assert_graph()
        original=self.amplitude.reader(delivery)
        base=correlated.locator.CancelledGraphBounds(self.graph,self.ctx,
            {**original.bindings,**self.finite_bindings},original.allowed_exponentials)
        reader=correlated.CorrelatedGraphBounds(base,original.aliases)
        if reader.bindings[self.pivot.node]._mpi_!=self.exact_singleton:
            raise ValueError('No rounding/substitution of the exact selected source origin')
        self._readers[id(reader)]=(reader,{node:value._mpi_ for node,value in reader.bindings.items()},
            dict(reader.aliases),reader.allowed_exponentials)
        return reader

    @source_precision
    def center(self,reader,log_function):
        if not hasattr(log_function,'node'):
            log_function=refined.box.pulse.radius.FunctionRef(self.graph,log_function)
        if getattr(log_function,'graph',None) is not self.graph:
            raise ValueError('Same original source graph function required')
        record=self._readers.get(id(reader))
        if record is None or record[0] is not reader or \
                {node:value._mpi_ for node,value in reader.bindings.items()}!=record[1] or \
                reader.aliases!=record[2] or reader.allowed_exponentials!=record[3]:
            raise ValueError('Unchanged source-centered reader issued by this owner required')
        if reader.graph is not self.graph or reader.ctx is not self.ctx or \
                reader.bindings[self.pivot.node]._mpi_!=self.exact_singleton:
            raise ValueError('Current exact source-centered reader required')
        # Mutable caller memo/polynomial caches never supply an enclosure.
        reader=correlated.CorrelatedGraphBounds(correlated.locator.CancelledGraphBounds(
            self.graph,self.ctx,dict(reader.bindings),reader.allowed_exponentials),dict(reader.aliases))
        polynomial=reader.polynomial(log_function)
        coefficient=polynomial.get((self.pivot.node,),Fraction(0))
        residual=dict(polynomial)
        residual.pop((self.pivot.node,),None)
        origin=self.graph.mul(self.graph.constant(coefficient),self.pivot)
        rest=signed.polynomial_function(self.graph,residual)
        if reader.polynomial(self.graph.sub(log_function,self.graph.add(origin,rest))):
            raise ValueError('Exact original fixed-origin identity required')
        encoding=exact_scaled_binary_rational(self.exact_singleton,coefficient)
        value=reader.at(rest)
        a,b=ends(value)
        width=self.ctx.mpf(b)-self.ctx.mpf(a)
        proof=self.graph.node('exact_original_Rp_fixed_binary_rational_scale_origin',
            original_log_scale_function=log_function.node if hasattr(log_function,'node') else log_function,
            exact_original_logC_source=self.pivot.node,exact_coefficient=signed.rational_record(coefficient),
            exact_center_function=origin.node,exact_center_encoding=encoding,
            exact_residual_function=rest.node,
            identity='L=q*exact_original_logC+R; exact fixed center has zero enclosure variation',
            proof_method='exact rational graph polynomial and normalized binary-rational source datum',
            remaining_logC_products_in_residual=any(self.pivot.node in atoms for atoms in residual),
            original_scale_function_not_replaced=True)
        return dict(exact_fixed_origin=dict(kind='exact_scaled_binary_rational',**encoding,
            exact_source_logC_node=self.pivot.node,coefficient=signed.rational_record(coefficient),
            defining_source_tuple=list(self.exact_singleton[0]),source_function=origin.node,
            zero_error_variation=True,no_exponent_sized_integer_shift=True),
            exact_residual_log_scale_function=rest.node,directed_residual_log_scale=value,
            residual_log_width_bound=width,exact_centered_scale_identity_proof=proof.node,
            residual_still_contains_nonconstant_logC_products=any(self.pivot.node in atoms for atoms in residual),
            original_full_log_scale_function=log_function.node if hasattr(log_function,'node') else log_function,
            ordinary_numeric_scale_not_materialized=True)

    @source_precision
    def evaluate(self,delivery,relative_width_target='1/100000000'):
        self.assert_graph()
        target=Fraction(relative_width_target)
        if not 0<target<1:raise ValueError('Exact relative width target in (0,1) required')
        _,_,request=self.amplitude.before._validate_delivery(delivery)
        if request.chart in refined.box.pulse.PULSE:
            mapped=self.before.physical_rows(delivery);provider='accepted_refined_current_pulse'
        else:
            mapped=self.amplitude.before.physical_rows(delivery);provider='accepted_original_postpulse_or_heat'
        reader=self.reader(delivery)
        values=self.before.before.before.before
        error_owner=self.before.before
        identities=dict(original_physical_inverse_identity=delivery['physical_inverse']['exact_original_inverse_identity'].node,
            original_native_inverse_identity=delivery['exact_native_inverse_identity'].node)
        sections={};counts=dict(rows=0,exact_zero_rows=0,nonzero_factored_target_rows=0,
            centered_nonzero_rows=0,remaining_nonconstant_origin_rows=0,ordinary_width_target_rows=0,
            ordinary_numeric_target_rows=0,sign_unresolved_rows=0)
        for section in ('Cartesian_spatial_rows','fixed_x_time_rows'):
            sections[section]={}
            for component,rows in mapped[section].items():
                sections[section][component]={}
                for label,row in (rows if section=='Cartesian_spatial_rows' else {'dt':rows}).items():
                    groups=values._canonical_groups(row,reader,identities)
                    result=signed.enclose_factored_sum(self.graph,values.ctx,reader,groups,target)
                    if result['exact_zero_enclosure']:
                        budget=accuracy.relative_error_budget(self.ctx,result,target)
                    else:
                        center=self.center(reader,result['exact_reference_log_scale_function'])
                        budget=centered_relative_budget(self.ctx,result,target,center)
                        result['centered_original_log_scale']=center
                        counts['centered_nonzero_rows']+=1
                        counts['remaining_nonconstant_origin_rows']+=int(center['residual_still_contains_nonconstant_logC_products'])
                    # Keep all signed coefficient/ratio/positive-tail ledgers.
                    legacy=error_owner.row_accuracy(result,reader,target)
                    for key in ('signed_coefficient_and_ratio_variation_ledger',
                                'coefficient_ledger_keeps_shared_source_dependencies',
                                'normalized_error_does_not_replace_absolute_physical_error',
                                'actual_common_scale_source_ledger'):
                        if key in legacy:budget[key]=legacy[key]
                    result.update(physical_accuracy=budget,component=component,physical_derivative=row.derivative,
                        original_parameter_arithmetic_proof=self.parameter_proof.node,
                        signed_canonical_scale_groups=[dict(exact_scale_log_function=g['log_function'].node,
                            signed_coefficient_enclosure=g['coefficient'],exact_scale_identity_proofs=g['proofs']) for g in groups])
                    sections[section][component][label]=result
                    counts['rows']+=1;counts['exact_zero_rows']+=int(result['exact_zero_enclosure'])
                    counts['nonzero_factored_target_rows']+=int(not result['exact_zero_enclosure'] and budget['factored_relative_width_satisfied'])
                    counts['ordinary_width_target_rows']+=int(budget['ordinary_numeric_relative_width_satisfied'])
                    counts['ordinary_numeric_target_rows']+=int(budget['ordinary_numeric_delivery_target_satisfied'])
                    counts['sign_unresolved_rows']+=int(budget.get('sign') is None)
        return dict(source_family=self.family_record,chart=request.chart,physical_value_rows=sections,
            delivery_counts=counts,source_coefficient_provider=provider,
            coefficient_history_partition=dict(main_pulse_charts_2048=request.chart in ('pulse_main','pulse_exit'),
                other_pulse_chart_history_partitions_unchanged=request.chart not in ('pulse_main','pulse_exit'),
                original_full_energy_integral_cells=32768 if request.chart in refined.box.pulse.PULSE else None,
                whole_chart_uniform_accuracy_not_certified=True),
            original_P0_source_binding=error_owner.before.source_binding(delivery),
            exact_source_origin_and_directed_residual_are_separate=True,
            original_physical_source_and_parameters_unchanged=True,
            full_certified_physical_accuracy=False,unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def velocity_pressure(self,delivery,relative_width_target='1/100000000'):
        _,_,request=self.amplitude.before._validate_delivery(delivery)
        view=self.evaluate(delivery,relative_width_target)
        rows=view['physical_value_rows']['Cartesian_spatial_rows']
        return dict(physical_coordinates={key:ref.node for key,ref in request.forward_coordinates.items()},
            values={name:rows[component]['x0_y0_z0'] for name,component in (('u','ux'),('v','uy'),('w','uz'),('p','p'))},
            original_P0_source_binding=view['original_P0_source_binding'],source_family=self.family_record,
            full_certified_physical_accuracy=False,unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


@source_precision
def run(before,deliveries,relative_width_target='1/1000'):
    began=time.monotonic()
    owner=CurrentOriginalRpCenteredScaleArithmetic(before,require_checked=False)
    views={name:owner.evaluate(delivery,relative_width_target) for name,delivery in deliveries.items()}
    result=dict(source_family=owner.family_record,actual_centered_physical_scale_views=views,
        source_graph_assertions=owner.assert_graph(),finite_parameter_source_inclusion=owner.inclusion,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,views
