"""Current exact scaled swirl source inside the original Cauchy seed boxes.

Uses the existing current bridge's actual atom/amplitude/rebuild objects.
The selected Cstar surplus is kept in logarithmic inequalities and connects
the transfer's weighted norm to the old unweighted S consumer. No source
exponential, parameter representative or nonlinear point solution is chosen.
"""
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import CurrentActualBridgeMixedC4
from lei_ren_part1_paper_compliant_current_core_recurrence_source import binding,function,sha
from lei_ren_part1_paper_compliant_core_physical_field import relative_amplitude_jet
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'current_core_scaled_swirl_source.json'
RECEIPT=PREFIX+'current_core_scaled_swirl_source_check.json'
GATE='current_core_exact_scaled_swirl_source_jets_admitted'
PARENTS=(PREFIX+'current_core_recurrence_source_check.json',
    PREFIX+'anchored_axis_amplitude_check.json',PREFIX+'current_actual_bridge_mixed_C4_check.json')
OPEN=('finite_rows_and_tails_bound_to_same_nonlinear_fixed_point',
    'current_core_bridge_functional_mixed4_join_certified','all_current_bridge_functional_interfaces_certified',
    'actual_point_moment_history_recovered','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','global_completed_tensor_admissibility',
    'admissible_stress_lift_constructed','temporal_recursion')


def defining_source_bindings():
    binding('compliant_anchored_axis_amplitude','evaluate','logF0','-self.core.logC-self.core.Lambda*actual_G')
    binding('compliant_core_coefficient_rebuild','__init__','required','self.core.Lambda*self.core.Gbar+2*self.core.logLambda+1000')
    binding('compliant_core_coefficient_rebuild','__init__','self.S_bound','c.exp(-6*self.core.logLambda-2000)')
    binding('compliant_core_coefficient_rebuild','seed','bound','self.S_bound/(self.eta/2)**k')
    binding('compliant_rooted_core_field','build_root_rows','bound','self.rebuild.S_bound/(self.rebuild.eta/2)**k')
    binding('compliant_core_transfer','run','weight',"get(linear,'Cauchy_weight_supremum_upper')")
    binding('compliant_core_transfer','run','Fsquare','ctx.exp(-4*logLambda-2000)*weight')
    binding('compliant_core_transfer','run','Gupper',"ctx.mpf(endpoints(get(tube,'G_modulus_upper'))[1])")
    binding('compliant_physical_norm_family','run','min_logC','b.hi(lam*G+2*logLam+1000)')
    binding('compliant_physical_norm_family','run','selected_logC','b.hi(2*max(endpoints(v)[1] for v in restrictions.values()))')
    binding('shared_linear_resolvent','run','weight','ctx.mpf([1,max(mp.mpf(1),endpoints(4*ratio)[1])])')
    binding('shared_linear_resolvent','run','ratio','h/(eta/2)')
    binding('shared_analytic_tube','run','h','eta/8')
    import ast
    seed=function('compliant_core_coefficient_rebuild','seed')
    expected=ast.dump(ast.parse('c.mpf([0,endpoints(bound)[1]]) if k==0 else symmetric(c,bound)',mode='eval').body)
    calls=[n for n in ast.walk(seed) if isinstance(n,ast.Call) and ast.unparse(n.func)=='S.append']
    if len(calls)!=1 or ast.dump(calls[0].args[0])!=expected:
        raise ValueError('Original positive/symmetric Cauchy source seed changed')
    z,j,delta,logLam,logC,G=s.symbols('Z j delta logLambda logC G',real=True)
    H=(1-delta)*z/2+(1-z*z)*(4*z+j)
    if s.expand(H-(-4*z**3-j*z*z+(9-delta)*z/2+j))!=0:
        raise ArithmeticError('Anchored and seed H source polynomials differ')
    if s.expand(-2*logLam+2*(-logC-s.exp(logLam)*G)
            -(-2*logLam-2*logC-2*s.exp(logLam)*G))!=0:
        raise ArithmeticError('Exact scaled swirl logarithm identity failed')
    return dict(original_anchored_amplitude_log_definition_bound=True,
        transfer_actual_Cauchy_weight_consumer_and_seed_bound_AST_bound=True,
        same_cubic_H_in_amplitude_and_core_seed=True,selected_Cstar_surplus_definition_bound=True,
        original_positive_value_and_symmetric_derivative_boxes_retained=True,
        exact_log_S='-2logLambda-2selected_logCstar-2Lambda*G(Z)',
        all_Cauchy_orders_rule='|S^(k)(x)/k!| <= selected_transfer_scaled_upper/(eta/2)^k',
        passed=True)


def selected_transfer_upper_log(c,logLambda,weighted_Fsquare_upper,selected_surplus):
    """Consume the weighted transfer bound; never exponentiate the surplus."""
    return -2*logLambda+c.ln(weighted_Fsquare_upper)-2*selected_surplus


class CurrentCoreScaledSwirlSource:
    @source_precision
    def __init__(self,bridge=None,require_checked=True):
        self.bridge=bridge if bridge is not None else CurrentActualBridgeMixedC4()
        if not self.bridge.acceptance_loaded:raise ValueError('Checked current actual bridge required')
        self.atoms=self.bridge.upstream.comparison.atoms
        self.rebuild=self.atoms.rebuild;self.amplitude=self.atoms.field.peak.amplitude
        self.core=self.rebuild.core;self.ctx=c=self.core.ctx
        self.family=self.core.family;self.source=self.core.source;self.datum_sha=self.core.datum.datum_sha
        self.hashes=dict(self.bridge.hashes);self.parents={}
        gates=('current_core_radial_production_scaling_certified',
            'actual_anchored_G_resolved_with_directed_error','current_actual_three_bridge_charts_mixed4_available')
        for name,gate in zip(PARENTS,gates):
            row=accepted(name,self.family,self.source,gate);_verify_hashes(row)
            if not row['all_passed'] or row['datum_enclosure_sha256']!=self.datum_sha:
                raise ValueError('Same current core source admission required')
            for key,digest in row['input_hashes'].items():
                if key in self.hashes and self.hashes[key]!=digest:raise ValueError('Current source conflict: '+key)
                self.hashes[key]=digest
            self.hashes[name]=sha(name);self.parents[name]=row
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.graph=self.current_source_graph()
        if not all(self.graph.values()):raise ValueError('Existing current amplitude/rebuild graph required')
        self.bindings=defining_source_bindings()
        major=self.core.records['core_transfer'];norm=self.core.records['physical_norm_family']
        tube=self.core.records['shared_analytic_tube'];linear=self.core.records['shared_linear_resolvent']
        base_family=major['analytic_core_family_sha256']
        if (norm['base_analytic_core_family_sha256']!=base_family or tube['analytic_core_family_sha256']!=base_family
                or linear['input_sha256']!=sha('lei_ren_part1_paper_shared_analytic_tube.json')
                or norm['definition']['pressure_source_sha256']!=self.source
                or norm['definition']['pressure_datum_sha256']!=self.datum_sha
                or norm['definition']['amplitude']!='F0=exp(-logCstar-Lambda*G)'
                or major['logC_definition']!='Lambda*Gupper + 2*logLambda + 1000'
                or not tube['common_complex_axis_poles_excluded']):
            raise ValueError('Original amplitude/tube/pressure definitions differ')
        read=lambda row,key:read_interval(c,row[key])
        self.eta=read(tube,'complex_tube_radius');self.radius=self.eta/2
        self.weight=read(linear,'Cauchy_weight_supremum_upper')
        self.weighted_Fsquare=read(major,'F0_squared_Xh_norm_upper')
        if (self.radius._mpi_!=read(linear,'Cauchy_disk_radius')._mpi_
                or self.eta._mpi_!=self.rebuild.eta._mpi_
                or self.core.Gbar._mpi_!=read(tube,'Gbar_for_symbolic_logC')._mpi_
                or self.core.logC._mpi_!=read(norm,'selected_logCstar')._mpi_
                or endpoints(self.weight)[0]<1 or endpoints(self.radius)[0]<=0):
            raise ValueError('Same source bounds and half-tube radius required')
        # Fsquare is the baseline weighted bound, before the additionally
        # selected Cstar surplus. Its actual weight is consumed, including
        # outward rounding; it is not silently replaced by the number one.
        # Replay the producer's declared precision, not the downstream
        # context's different rounding precision. No interval overlap test.
        producer_ctx=MPIntervalContext();producer_ctx.dps=major['precision']
        producer_logLambda=read_interval(producer_ctx,major['logLambda'])
        producer_weight=read_interval(producer_ctx,linear['Cauchy_weight_supremum_upper'])
        expected=producer_ctx.exp(-4*producer_logLambda-2000)*producer_weight
        if expected._mpi_!=read_interval(producer_ctx,major['F0_squared_Xh_norm_upper'])._mpi_:
            raise ValueError('Transfer weighted source formula changed')
        guard=self.core.Lambda*self.core.Gbar+2*self.core.logLambda+1000
        self.surplus=self.core.logC-guard
        if endpoints(self.surplus)[0]<1000:raise ValueError('Selected complex-amplitude surplus unavailable')
        self.required_log=selected_transfer_upper_log(c,self.core.logLambda,self.weighted_Fsquare,self.surplus)
        self.consumer_log=c.ln(self.rebuild.S_bound)
        self.margin=self.consumer_log-self.required_log
        if endpoints(self.margin)[0]<=0:raise ValueError('Existing S consumer fails weighted selected-source dominance')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Scaled swirl datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def current_source_graph(self):
        return dict(same_actual_atom_rebuild=self.rebuild is self.atoms.rebuild is self.amplitude.rebuild,
            same_current_original_core=self.core is self.amplitude.core is self.bridge.upstream.core is self.bridge.core.original,
            same_actual_anchored_amplitude=self.amplitude is self.atoms.field.peak.amplitude,
            same_current_records=self.core.records is self.bridge.core.records,
            same_interval_context=self.ctx is self.bridge.ctx is self.amplitude.ctx,
            same_current_family_source_datum=(self.family,self.source,self.datum_sha)==(
                self.bridge.family,self.bridge.source,self.bridge.datum_sha))

    @source_precision
    def cauchy_rows(self,fixed):
        c=self.ctx;rows=[]
        for k,box in enumerate(fixed['S_Z_taylor']):
            lo,hi=endpoints(box);upper=max(abs(lo),abs(hi))
            if upper<=0 or (k==0 and lo!=0) or (k>0 and lo!=-hi):
                raise ValueError('Original positive/symmetric Cauchy seed semantics changed')
            if upper!=endpoints(self.rebuild.S_bound/self.radius**k)[1]:
                raise ValueError('Original Cauchy seed does not use the same source bound/radius')
            required=self.required_log-k*c.ln(self.radius)
            margin=c.ln(c.mpf(upper))-required
            if endpoints(margin)[0]<=0:raise ValueError('Exact source Taylor coefficient is not admitted')
            rows.append(dict(axial_Taylor_order=k,ordinary_derivative_divisor_factorial=math.factorial(k),
                original_seed_box=box,selected_source_Taylor_log_upper=required,
                consumer_log_dominance_margin=margin,exact_positive_source_not_zero_or_cap=True))
        return rows

    @source_precision
    def certified_seed(self,Z,degree=24,depth=6):
        fixed,rows=self.rebuild.seed(Z,degree,depth)
        proof=self.cauchy_rows(fixed)
        return dict(Z=self.ctx.mpf(Z),radial_degree=degree,axial_depth=depth,
            fixed=fixed,rows=rows,exact_S_source='epsilon_core^2 exp(-2selected_logCstar-2Lambda*G(Z))',
            admitted_S_axial_Taylor_coefficients=proof,seed_rebuilt_from_current_original_callable=True,
            exact_positive_S_retained_formally=True,source_value_or_derivative_not_selected=True,
            **{GATE:self.acceptance_loaded,'exact_scaled_S_source_jet_admitted':self.acceptance_loaded},
            **dict.fromkeys(OPEN,False))

    @source_precision
    def logarithmic_source(self,Z):
        packet=self.amplitude.evaluate(Z);c=self.ctx
        logS=-2*self.core.logLambda+2*packet['logF0']
        relative=self.core.axis_inputs(packet['Z'])['relative_F0_squared_derivatives']
        if endpoints(logS-self.consumer_log)[1]>=0:
            raise ValueError('Actual anchored scaled source exceeds admitted Cauchy value bound')
        return dict(Z=packet['Z'],log_S_source_enclosure=logS,
            relative_ordinary_S_derivatives_through5=relative,
            exact_source_definition='S^(k)=exp(log_S)*relative_ordinary_S_derivative[k]',
            original_anchored_amplitude_used=True,exp_log_S_not_materialized=True,
            source_not_replaced_by_Cauchy_box=True,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))

    @source_precision
    def report(self):
        seeds={'whole_axis':self.certified_seed([-1,1]),'fresh_Z':self.certified_seed('.271'),
            'source_anchor':self.certified_seed(self.amplitude.roots['anchor_interval'])}
        root=self.atoms.field.root_rows
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_source_graph=self.current_source_graph(),
            defining_source_bindings=self.bindings,complex_tube_receipt_sha256=sha('lei_ren_part1_paper_shared_analytic_tube.json'),
            Cauchy_radius=self.radius,actual_transfer_Cauchy_weight=self.weight,
            actual_transfer_weighted_F0_squared_norm_upper=self.weighted_Fsquare,
            selected_logCstar_surplus=self.surplus,selected_scaled_source_log_upper=self.required_log,
            original_S_consumer_log_upper=self.consumer_log,selected_source_dominance_margin=self.margin,
            current_original_seed_packets=seeds,existing_actual_root_seed_S_rows=self.cauchy_rows(root['fixed']),
            actual_logarithmic_source_packets={Z:self.logarithmic_source(Z) for Z in ('.271','anchor')},
            baseline_weighted_source_bound_adjusted_by_exact_selected_Cstar_surplus=True,
            same_original_Cauchy_seed_consumers_proved_to_enclose_exact_source=True,
            all_nonnegative_Taylor_orders_by_Cauchy_theorem=True,
            no_source_exponential_materialized=True,no_old_seed_rows_changed=True,
            **{GATE:False,'exact_scaled_S_source_jet_admitted':False},**dict.fromkeys(OPEN,False),input_hashes=self.hashes)


def run(bridge=None):
    result=CurrentCoreScaledSwirlSource(bridge=bridge,require_checked=False).report()
    (HERE/NAME).write_bytes((json.dumps(encode(result),indent=2)+'\n').encode('utf8'))
    print('Current exact scaled swirl source and original Cauchy consumer packets generated',flush=True)
    return result


if __name__=='__main__':run()
