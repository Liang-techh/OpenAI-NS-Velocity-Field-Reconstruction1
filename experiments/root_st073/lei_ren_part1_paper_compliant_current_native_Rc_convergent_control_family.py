"""Exact source-bound C1 fixed-point family and cached finite evaluation.

For each source/local-admitted finite integer N the realized original C1
target feeds one Banach map. Its exact C1 limit and terminal identities are
defined here. Numerical evaluation remains explicit and approximate: the
mathematical Picard tail does not cover oracle/quadrature/roundoff errors.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_C1_integral_realization as integrals
import lei_ren_part1_paper_compliant_current_native_band_function_evaluator as evaluator

allN,band,controls,source,packets=integrals.allN,integrals.band,integrals.controls,integrals.source,integrals.packets
HERE,PREFIX,sha,ep,require=integrals.HERE,integrals.PREFIX,integrals.sha,integrals.ep,integrals.require
NAME=PREFIX+'current_native_Rc_convergent_control_family.json'
RECEIPT=PREFIX+'current_native_Rc_convergent_control_family_check.json'
GATE='current_original_realized_C1_target_convergent_control_family_and_terminal_limit_defined'


def generic_control_operator(built):
    """One exact map template, not a fixed-depth or cap-valued solution."""
    g=built['graph'];W=controls.exact_weights(g,built['parameters']['mu']);N=built['N']
    h=[source.C1Function(g.symbol('repair_control_'+str(i)),g.symbol('repair_control_Z_'+str(i))) for i in range(5)]
    d=[built['N_scaled_targets'][key] for key in controls.ROWS]
    overN=lambda v:g.quotient(v,N,'one common original positive integer N')
    Q=controls.quadratic(g,W,[q.value for q in h])
    QZ=controls.quadratic_action(g,W,[q.value for q in h],[q.Z for q in h])
    value=controls.inverse_action(g,W,[g.add(q.value,overN(r)) for q,r in zip(d,Q)])
    jet=controls.inverse_action(g,W,[g.add(q.Z,overN(r)) for q,r in zip(d,QZ)])
    next_h=controls.c1_vectors(g,[g.neg(v) for v in value],[g.neg(v) for v in jet])
    residual=controls.c1_vectors(g,
        [g.add(q.value,g.add(*(g.mul(b,v.value) for b,v in zip(row,h))),overN(r)) for q,row,r in zip(d,W['B'],Q)],
        [g.add(q.Z,g.add(*(g.mul(b,v.Z) for b,v in zip(row,h))),overN(r)) for q,row,r in zip(d,W['B'],QZ)])
    column=[controls.quadratic_action(g,W,[q.value for q in h],[g.one if j==i else g.zero for j in range(5)]) for i in range(5)]
    jacobian=[[g.add(W['B'][i][j],overN(column[j][i])) for j in range(5)] for i in range(5)]
    invN=g.quotient(g.one,N,'one common original positive integer N')
    history={key:g.c1add(g.c1scale(invN,built['coefficient_history'][-1][key]),
        g.c1scale(g.mul(invN,invN),built['coefficient_history'][-2][key])) for key in allN.RATES}
    family=dict(**built,history=history,repair_weights=W,finite_picard_sequence=[h],control_residual=residual,
        implicit_Z_matrix=jacobian,implicit_Z_rhs=[g.neg(q.Z) for q in d],
        band_variable='repair_x',generic_controls=h,generic_Picard_map=next_h)
    return band.exact_partial_band_graph(family)


@dataclass(frozen=True)
class FiniteControlEvaluation:
    depth: int
    values: tuple
    Z_derivatives: tuple
    residual_values: tuple
    residual_Z_derivatives: tuple
    approximate_only: bool=True
    original_numeric_fixed_point_certified: bool=False


class CachedC1ControlEvaluator:
    """Evaluate the exact same map at arbitrary finite depth with an oracle.

    The cached evaluator retains actual original target/source integrals.
    It never obtains a coefficient from a range or log-norm endpoint.
    """
    def __init__(self,built,*,oracle,Z,N,ctx,original_log_N_required=None):
        if getattr(oracle,'mode',None)=='original_source':
            if original_log_N_required is None:raise ValueError('Original source/local frequency contract required')
            actual=ctx.log(N);required=ep(original_log_N_required)[1]
            require(actual>=required,'Original candidate N does not satisfy the current source/local C1 contract')
        self.built,self.ctx=built,ctx
        self.value=evaluator.BandFunctionEvaluator(built,oracle=oracle,Z=Z,N=N,ctx=ctx)

    def variables(self,h,hZ):
        return {'Z':self.value.Z,**{'repair_control_'+str(i):v for i,v in enumerate(h)},
                **{'repair_control_Z_'+str(i):v for i,v in enumerate(hZ)}}

    def evaluate(self,depth):
        if type(depth) is not int or not 1<=depth<=4096:raise ValueError('Finite Picard depth must lie in[1,4096]')
        h=hZ=(self.ctx.mpf(0),)*5
        for unused in range(depth):
            variables=self.variables(h,hZ)
            next_h=tuple(self.value(q.value,variables) for q in self.built['generic_Picard_map'])
            next_Z=tuple(self.value(q.Z,variables) for q in self.built['generic_Picard_map'])
            h,hZ=next_h,next_Z
        variables=self.variables(h,hZ)
        residual=tuple(self.value(q.value,variables) for q in self.built['control_residual'])
        residual_Z=tuple(self.value(q.Z,variables) for q in self.built['control_residual'])
        return FiniteControlEvaluation(depth,h,hZ,residual,residual_Z)


class NativeRcConvergentControlFamily:
    def __init__(self):
        self.integrals=integrals.NativeRcC1IntegralRealization()
        self.ctx,self.family=self.integrals.ctx,self.integrals.family;self.hashes=dict(self.integrals.hashes)
        checked=json.loads((HERE/integrals.RECEIPT).read_bytes())
        require(checked['all_passed'] and checked[integrals.GATE] and checked['source_family']==self.family,
                'Realized original C1 target functions required before the fixed-point family')
        self.hashes.update(checked['input_hashes']);self.hashes[integrals.RECEIPT]=sha(integrals.RECEIPT)
        self.conditions=self.integrals.data['current_native_Rc_all_N_function_controls']['actual_uniform_repair_C1_log_N_conditions']
        prior=integrals.previous
        row=json.loads((HERE/prior.NAME).read_bytes())
        require(self.hashes[prior.NAME]==sha(prior.NAME),'Original source/local frequency budget changed')
        self.logN=packets.interval(self.ctx,row['source_repair_active_flat_quiet_band_C0_sufficient_log_N_lower'])
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.hashes[Path(evaluator.__file__).name]=sha(Path(evaluator.__file__).name)

    def contract(self):
        c=self.ctx;read=lambda key:packets.interval(c,self.conditions[key]['log_absolute_upper'])
        CA,CQ,rho,D=[read(key) for key in ('exact_integral_matrix_inverse_log_cap','general_transformed_quadratic_C1_log_cap',
            'formal_control_C1_ball_radius_log','actual_target_C1_log_cap')]
        require(self.conditions['B_mu_Z_independent_and_C1_product_norm_submultiplicative'] and
                self.conditions['unique_C1_control_functions_in_certified_ball_for_any_compatible_actual_C1_target'],
                'Same original Z-independent inverse and full C1 quadratic map required')
        require(ep(self.logN-packets.interval(c,self.conditions['repair_sufficient_common_log_N_lower']))[0]>=0,
                'Current source/local frequency does not contain repair conditions')
        logL=c.mpf(ep(c.ln(2)+CA+CQ+rho-self.logN)[1])
        require(ep(logL+c.ln(2))[1]<=0,'Current original C1 operator does not contract by at most half')
        require(ep(rho-CA-D-c.ln(2))[1]>=0,'Original C1 ball does not contain twice the first Picard cap')
        return dict(source_family=self.family,whole_Z_domain=[-1,1],
            actual_source_local_sufficient_log_N=self.logN,
            actual_uniform_C1_contraction_log_upper=logL,uniform_contraction_at_most='1/2',
            original_common_control_C1_ball_radius=self.conditions['formal_control_C1_ball_radius_log'],
            exact_target_C1_functions_receipt=integrals.RECEIPT,
            exact_limit='h_star(N,Z)=lim_C1 h_n; h0=0; h_next=-B_exact(mu)^-1*(d_exact(N,Z)+Q_exact(mu,h)/N)',
            target_is_the_realized_original_integral_function_not_an_arbitrary_C1_input=True,
            actual_global_integer_N_not_selected=True,integer_or_exp_logN_not_materialized=True,
            tail='norm_C1(h_star-h_n)<=rho*2^-n',
            preconditioned_residual='norm_C1(B^-1*(B*h_n+d+Q(h_n)/N))<=3*rho*2^(-n-1)',
            implicit_Z_identity='(B+DQ(h_star)/N)*h_star_Z=-d_Z',
            implicit_Z_inverse='(I+B^-1*DQ(h_star)/N)^-1*B^-1, Neumann norm<=2*CA',
            exact_limit_exists_unique_in_original_C1_ball=True,
            uniform_C1_convergence_not_numerical_iteration_error_bound=True)

    def factored_tail(self,depth):
        if type(depth) is not int or depth<0:raise ValueError('Nonnegative exact Picard depth required')
        # Retain the small dyadic factor separately: adding -n*log2 to
        # the astronomical log-radius at fixed precision would erase it.
        return dict(depth=depth,original_radius=self.conditions['formal_control_C1_ball_radius_log'],
            relative_control_C1_tail=dict(numerator=1,denominator=2**depth),
            relative_preconditioned_residual_C1_tail=dict(numerator=3,denominator=2**(depth+1)),
            exact_mathematical_iteration_tail=True,
            numerical_oracle_quadrature_and_roundoff_errors_included=False)

    def build(self):return generic_control_operator(self.integrals.build())


def run():
    began=time.monotonic();owner=NativeRcConvergentControlFamily();built=owner.build();contract=owner.contract()
    rows=controls.pair_roots
    record=dict(source_family=owner.family,**{GATE:True},
        exact_function_graph_nodes=built['graph'].nodes,exact_generic_control_roots=rows(built['generic_controls'],controls.CONTROLS),
        exact_generic_Picard_C1_map_roots=rows(built['generic_Picard_map'],controls.CONTROLS),
        exact_control_residual_roots=rows(built['control_residual'],controls.ROWS),
        exact_implicit_Z_matrix_roots=[[v.node for v in row] for row in built['implicit_Z_matrix']],
        exact_implicit_Z_rhs_roots=[v.node for v in built['implicit_Z_rhs']],
        exact_limit_family_partial_history_roots=rows(built['partial_band_histories'].values(),built['partial_band_histories'].keys()),
        exact_limit_family_terminal_residual_identity_roots=rows(built['finite_terminal_residual_identities'].values(),built['finite_terminal_residual_identities'].keys()),
        exact_original_C1_limit_contract=contract,factored_iteration_tail_examples=[owner.factored_tail(n) for n in (0,3,16,64)],
        exact_same_source_C1_fixed_point_family_defined=True,
        terminal_C1_limit_identities='Dm=A*resM/(2N),Dh=A*resI/(2^(3/2)N),De=A²*resS/(2N),Dp=A²*resCp/N,Dk=A²*(resM+mu*resD)/(2^(3/2)N); all five residual functions and their Z derivatives tend uniformly to zero',
        exact_family_terminal_functions_zero_in_C1_limit_for_each_source_local_admitted_N=True,
        finite_Picard_terminal_functions_declared_zero=False,
        cached_arbitrary_depth_C1_evaluator_implemented=True,
        evaluator_is_approximate_unless_oracle_and_roundoff_are_separately_certified=True,
        mathematical_tail_does_not_certify_numeric_quadrature_errors=True,
        actual_source_ancestor_constructors_called=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Genuine exact C1 fixed-point family driven by realized original all24 integrals for each source/local-admitted finite N, arbitrary-depth cached value/Z iteration, factored rigorous mathematical tails and terminal C1 limit identities. No selected global N, evaluated/certified original numeric oracle or installed numeric controls/tail/terminal field, spatial higher seams, global outer cone or recursion.')
    (HERE/NAME).write_text(json.dumps(packets.encode(record),indent=2)+'\n',encoding='utf8')
    print('Realized original C1 target feeds exact convergent control family and terminal C1 limit',flush=True)
    return record


if __name__=='__main__':run()
