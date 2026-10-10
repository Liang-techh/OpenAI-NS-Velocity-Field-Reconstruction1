"""Actual original physical rows as signed common-scale value enclosures.

Exact scale equality is proved before signed coefficients are combined.
Distinct scales retain their exact ratio functions and positive tail bounds.
Ordinary numbers are optional representations, never substitute parameters.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_amplitude_binding as amplitude
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

correlated=amplitude.correlated
physical=correlated.physical
HERE,PREFIX,sha=amplitude.HERE,amplitude.PREFIX,amplitude.sha
NAME=PREFIX+'current_original_Rp_signed_physical_values.json.gz'
RECEIPT=PREFIX+'current_original_Rp_signed_physical_values_check.json'
GATES=('current_original_Rp_proved_canonical_physical_scale_groups_installed',
       'current_original_Rp_signed_common_scale_value_enclosures_installed',
       'current_original_Rp_positive_ratio_tail_and_cancellation_bounds_installed',
       'current_original_Rp_supported_family_uvw_pressure_value_interface_installed')
OPEN=amplitude.OPEN
ends=amplitude.ends
# This is an evaluation resource guard, not a radius/amplitude parameter.
EXP_LOG_LIMIT=1000


def rational_record(value):
    value=Fraction(value)
    return dict(numerator=value.numerator,denominator=value.denominator)


def polynomial_function(graph,polynomial):
    terms=[]
    for atoms,coefficient in sorted(polynomial.items()):
        refs=[correlated.box.pulse.radius.FunctionRef(graph,node) for node in atoms]
        terms.append(graph.mul(graph.constant(coefficient),*refs))
    return graph.add(*terms) if terms else graph.zero


def signed_log_enclosure(ctx,coefficient,scale_log):
    """A coefficient interval crossing zero has no finite log lower bound."""
    lo,hi=ends(coefficient)
    if lo==hi==0:return dict(sign=0,exact_zero_enclosure=True,log_magnitude=None)
    if lo>0 or hi<0:
        sign=1 if lo>0 else -1
        magnitude=coefficient if sign>0 else -coefficient
        return dict(sign=sign,exact_zero_enclosure=False,
            log_magnitude=scale_log+ctx.ln(magnitude),coefficient_sign_certified=True)
    maximum=max(abs(lo),abs(hi))
    upper=scale_log+ctx.ln(ctx.mpf(maximum))
    return dict(sign=None,exact_zero_enclosure=False,coefficient_sign_certified=False,
        log_magnitude_lower='-infinity',log_magnitude_upper=ctx.mpf(ends(upper)[1]),
        reason='The signed coefficient enclosure includes zero; no nonzero magnitude is inferred')


def ratio_enclosure(ctx,log_ratio):
    """Monotonic exponential bounds; retain positive tails without exp(huge)."""
    lo,hi=ends(log_ratio)
    if hi>EXP_LOG_LIMIT:
        return None,dict(kind='unresolved_ratio',reason='Positive log-ratio upper bound exceeds the evaluation guard')
    if hi<=-EXP_LOG_LIMIT:
        upper=ends(ctx.exp(ctx.mpf(-EXP_LOG_LIMIT)))[1]
        return ctx.mpf([0,upper]),dict(kind='retained_positive_tail',
            exact_ratio_is_strictly_positive=True,upper_log_cutoff=-EXP_LOG_LIMIT,
            zero_lower_is_only_an_enclosure=True,tail_not_replaced_by_zero=True)
    if lo < -EXP_LOG_LIMIT:
        upper=ends(ctx.exp(ctx.mpf(hi)))[1]
        return ctx.mpf([0,upper]),dict(kind='monotone_upper_with_unresolved_positive_lower',
            exact_ratio_is_strictly_positive=True,zero_lower_is_only_an_enclosure=True)
    return ctx.exp(log_ratio),dict(kind='directed_finite_exponential',exact_ratio_is_strictly_positive=True)


@source_precision
def enclose_factored_sum(graph,ctx,reader,groups,relative_width_target):
    """Factor around an existing exact scale; reference choice changes no field."""
    active=[group for group in groups if ends(group['coefficient'])!=(0,0)]
    if not active:
        return dict(exact_zero_enclosure=True,ordinary_numeric_enclosure=ctx.mpf(0),
            ordinary_numeric_materialized=True,requested_relative_width_satisfied=True,
            signed_log_value=dict(sign=0,exact_zero_enclosure=True,log_magnitude=None),ratio_terms=[])
    # Choose one real source scale as the reference. Its bound endpoint only
    # chooses which exact function to factor out; it never defines a value.
    reference=max(active,key=lambda group:(ends(group['log_bound'])[0],-group['log_function'].node))
    common=ctx.mpf(0);ratios=[];unresolved=[]
    for group in active:
        difference=graph.sub(group['log_function'],reference['log_function'])
        bound=reader.at(difference)
        equal=not reader.polynomial(difference)
        value,method=(ctx.mpf(1),dict(kind='exact_common_scale',exact_ratio_is_one=True)) if equal else ratio_enclosure(ctx,bound)
        ratio=graph.unary('exp',difference)
        record=dict(source_log_scale=group['log_function'].node,
            exact_ratio_function=ratio.node,exact_log_ratio_function=difference.node,
            directed_log_ratio=bound,directed_ratio_enclosure=value,
            signed_coefficient=group['coefficient'],method=method)
        if value is None:unresolved.append(record)
        else:common+=group['coefficient']*value
        ratios.append(record)
    result=dict(exact_reference_log_scale_function=reference['log_function'].node,
        directed_reference_log_scale=reference['log_bound'],ratio_terms=ratios,
        reference_is_an_existing_exact_source_scale=True,coefficient_or_scale_midpoint_not_selected=True,
        identity='value=exp(Lref)*sum(c_i*exp(L_i-Lref))',
        exact_zero_enclosure=False,ordinary_numeric_materialized=False,
        requested_relative_width_satisfied=False)
    if unresolved:
        result.update(common_scale_coefficient_enclosure=None,signed_log_value=None,
            unresolved_accuracy_reason='Some distinct-scale ratios lack a feasible numeric enclosure; original signed factored groups remain the value representation')
        return result
    result['common_scale_coefficient_enclosure']=common
    result['signed_log_value']=signed_log_enclosure(ctx,common,reference['log_bound'])
    if ends(common)==(0,0):
        result.update(exact_zero_enclosure=True,ordinary_numeric_enclosure=ctx.mpf(0),
            ordinary_numeric_materialized=True,requested_relative_width_satisfied=True,
            exact_cancellation_before_scale_expansion=True)
        return result
    lo,hi=ends(reference['log_bound'])
    if min(lo,hi)<-EXP_LOG_LIMIT or max(lo,hi)>EXP_LOG_LIMIT:
        result['unresolved_accuracy_reason']='Original common scale cannot be expanded within the evaluation guard; keep its exact function, signed coefficient enclosure and directed log magnitude'
        return result
    value=common*ctx.exp(reference['log_bound']);a,b=ends(value)
    result.update(ordinary_numeric_enclosure=value,ordinary_numeric_materialized=True)
    if a==b==0:
        result['requested_relative_width_satisfied']=True
    elif a>0 or b<0:
        width=(ctx.mpf(b)-ctx.mpf(a))/ctx.mpf(min(abs(a),abs(b)))
        result['relative_width_bound']=width
        target=Fraction(relative_width_target)
        target_bound=ctx.mpf(target.numerator)/target.denominator
        result['requested_relative_width_satisfied']=ends(width)[1]<=ends(target_bound)[0]
    else:
        result['unresolved_accuracy_reason']='Signed numeric enclosure includes zero, so relative error is unresolved'
    if not result['requested_relative_width_satisfied'] and 'unresolved_accuracy_reason' not in result:
        result['unresolved_accuracy_reason']='Full source, coordinate and scale enclosure width exceeds the requested relative width'
    return result


class CurrentOriginalRpSignedPhysicalValues:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else amplitude.CurrentOriginalRpAmplitudeBinding()
        if type(self.before) is not amplitude.CurrentOriginalRpAmplitudeBinding or not self.before.acceptance_loaded:
            raise ValueError('Accepted current exact amplitude/source binding required')
        self.before.assert_graph();self.graph=self.before.graph;self.ctx=self.before.ctx
        self.family_record=self.before.family_record;self.hashes=dict(self.before.hashes)
        self.acceptance_loaded=False
        for name in (amplitude.NAME,amplitude.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN) or \
                    receipt['source_family']!=self.family_record:
                raise ValueError('Current signed physical value receipt or scope differs')
            correlated.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        self.before.assert_graph()
        checks=dict(accepted_exact_current_amplitude=self.before.acceptance_loaded,
            same_current_graph_context=self.graph is self.before.graph and self.ctx is self.before.ctx,
            same_actual_correlated_source=self.before.before.acceptance_loaded,
            unchanged_current_source_family=self.family_record==self.before.family_record)
        if not all(checks.values()):raise ValueError('Current signed physical source differs')
        return checks

    def _canonical_groups(self,row,reader,identity_nodes):
        if type(row) is not physical.PhysicalSourceRow:
            raise TypeError('Actual original physical source row required')
        groups={};g=self.graph
        for source in row.groups():
            logref=g.add(*[ref for _,ref in source['log_scale_parts']]);poly=reader.polynomial(logref)
            key=(source['scale_units'],source['exact_source_scale_powers'],tuple(sorted(poly.items())))
            if key not in groups:
                canonical=polynomial_function(g,poly)
                if reader.polynomial(g.sub(logref,canonical)):
                    raise ValueError('Exact current scale canonicalization identity required')
                groups[key]=dict(log_function=canonical,log_bound=reader.at(canonical),
                    coefficient=self.ctx.mpf(0),sources=[],proofs=[],polynomial=poly,
                    source_units=source['scale_units'],source_powers=source['exact_source_scale_powers'])
            group=groups[key]
            if reader.polynomial(g.sub(logref,group['log_function'])):
                raise ValueError('Only exactly equal current positive scales may be merged')
            coefficient=source['signed_coefficient']
            if coefficient.ctx is not self.ctx or coefficient.order!=0:
                raise ValueError('Same current ordinary signed coefficient enclosure required')
            group['coefficient']+=coefficient[0]
            proof=g.node('exact_original_Rp_canonical_physical_positive_scale_identity',
                original_scale_log_function=logref.node,canonical_scale_log_function=group['log_function'].node,
                amplitude_identity=self.before.proof.node,**identity_nodes,
                identical_source_units=source['scale_units'],identical_source_scale_powers=[
                    rational_record(v) for v in source['exact_source_scale_powers']],
                proof_method='exact rational polynomial under validated original source/root/native/amplitude identities',
                same_physical_component=row.component,ordinary_physical_derivative=row.derivative)
            group['proofs'].append(proof.node)
            group['sources'].append(dict(original_scale_log_function=logref.node,
                signed_original_group_coefficient=coefficient,contributors=len(source['contributors']),
                source_units=source['scale_units'],source_powers=[rational_record(v) for v in source['exact_source_scale_powers']]))
        return list(groups.values())

    @source_precision
    def evaluate(self,delivery,relative_width_target='1/100000000'):
        self.assert_graph()
        target=correlated.inverse.exact_scalar(relative_width_target)
        if not 0<target<1:raise ValueError('Exact relative width target in (0,1) required')
        reader=self.before.reader(delivery);mapped=self.before.before.physical_rows(delivery)
        identities=dict(original_physical_inverse_identity=mapped['correlated_original_physical_input_identity'].node,
            original_native_inverse_identity=mapped['correlated_original_native_inverse_identity'].node)
        sections={};counts=dict(rows=0,literal_groups=0,canonical_groups=0,numeric_rows=0,
            relative_width_satisfied_rows=0,retained_positive_ratio_tails=0)
        for section in ('Cartesian_spatial_rows','fixed_x_time_rows'):
            sections[section]={}
            for component,rows in mapped[section].items():
                values=rows if section=='Cartesian_spatial_rows' else {'dt':rows}
                sections[section][component]={}
                for label,row in values.items():
                    groups=self._canonical_groups(row,reader,identities)
                    result=enclose_factored_sum(self.graph,self.ctx,reader,groups,target)
                    reported=[]
                    for group in groups:
                        reported.append(dict(exact_scale_log_function=group['log_function'].node,
                            directed_scale_log=group['log_bound'],signed_coefficient_enclosure=group['coefficient'],
                            signed_log_group=signed_log_enclosure(self.ctx,group['coefficient'],group['log_bound']),
                            exact_scale_identity_proofs=group['proofs'],original_groups=group['sources'],
                            same_source_units=group['source_units'],same_source_scale_powers=[
                                rational_record(v) for v in group['source_powers']],
                            exact_polynomial=[dict(atoms=list(atoms),coefficient=rational_record(value))
                                for atoms,value in sorted(group['polynomial'].items())]))
                    result.update(component=component,physical_derivative=row.derivative,
                        signed_canonical_scale_groups=reported,actual_source_and_original_physical_operators_retained=True,
                        current_P0_numeric_inclusion_pending=component=='p',
                        ordinary_pressure_materialization_not_a_current_P0_binding=component=='p')
                    sections[section][component][label]=result
                    counts['rows']+=1;counts['literal_groups']+=len(row.groups());counts['canonical_groups']+=len(groups)
                    counts['numeric_rows']+=int(result['ordinary_numeric_materialized'])
                    counts['relative_width_satisfied_rows']+=int(result['requested_relative_width_satisfied'])
                    counts['retained_positive_ratio_tails']+=sum(v['method']['kind']=='retained_positive_tail' for v in result['ratio_terms'])
        return dict(source_family=self.family_record,chart=mapped['chart'],exact_Z=mapped['exact_Z'],
            exact_native_coordinate=mapped['exact_native_coordinate'],relative_width_target=rational_record(target),
            actual_current_amplitude_identity=self.before.proof.node,**identities,
            physical_value_rows=sections,delivery_counts=counts,exponential_log_evaluation_guard=EXP_LOG_LIMIT,
            exact_zero_and_sign_uncertainty_preserved=True,all_original_source_owners_unmutated=True,
            original_astronomical_scale_not_replaced=True,full_physical_accuracy_and_P0_source_binding_still_open=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def velocity_pressure(self,delivery,relative_width_target='1/100000000'):
        """u/v/w/p at the validated exact physical point, with honest accuracy.

        This supported-family interface is not an unrestricted x/y/z/t API.
        Components retain factored/log values when ordinary numbers fail.
        """
        _,_,request=self.before.before._validate_delivery(delivery)
        evaluated=self.evaluate(delivery,relative_width_target)
        rows=evaluated['physical_value_rows']['Cartesian_spatial_rows']
        return dict(physical_coordinates={key:ref.node for key,ref in request.forward_coordinates.items()},
            values={name:rows[component]['x0_y0_z0'] for name,component in
                (('u','ux'),('v','uy'),('w','uz'),('p','p'))},
            convention='u=ux, v=uy, w=uz; original physical Cartesian components',
            supported_domain='current exact source-parameterized correlated physical input family',
            chart=request.chart,exact_native_coordinate=str(request.native),exact_Z=str(request.Z),
            source_family=self.family_record,actual_current_amplitude_identity=self.before.proof.node,
            source_and_coordinate_enclosures_retained=True,
            unrestricted_physical_point_API=False,current_P0_numeric_inclusion_pending=True,
            full_certified_physical_accuracy=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


@source_precision
def run(before=None,deliveries=None):
    began=time.monotonic();owner=CurrentOriginalRpSignedPhysicalValues(before,require_checked=False)
    if not deliveries:raise ValueError('Actual live current correlated physical source deliveries required')
    views={name:owner.evaluate(delivery) for name,delivery in deliveries.items()}
    result=dict(candidate_current_signed_physical_value_delivery_constructed=True,source_family=owner.family_record,
        actual_current_signed_physical_values=views,source_graph_assertions=owner.assert_graph(),
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(correlated.report(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('Current original signed canonical physical values generated',flush=True)
    return owner,views


if __name__=='__main__':run()
