"""Current exact inlet amplitude, with directed source-bound log enclosures.

Both original sigma primitives have endpoint value 1/2 by reflection.
This permits an exact compact log amplitude before interval arithmetic.
The copied production U box is admitted only after an explicit inclusion.
No source owner is mutated and no enormous physical exponential is taken.
"""
import json
from pathlib import Path
import time
import sympy as s
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_compliant_current_original_Rp_correlated_physical_source as correlated
import lei_ren_part1_paper_compliant_current_original_Rp_native_constants as exact
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=correlated.HERE,correlated.PREFIX,correlated.sha
NAME=PREFIX+'current_original_Rp_amplitude_binding.json'
RECEIPT=PREFIX+'current_original_Rp_amplitude_binding_check.json'
GATES=('current_original_Rp_exact_U0_numeric_source_inclusion_installed',
       'current_original_Rp_exact_logU0_directed_source_binding_installed',
       'current_original_Rp_proved_amplitude_log_scale_reader_installed')
OPEN=correlated.OPEN
ends=correlated.ends
zero=exact.frame.bridge.leading_check.zero


def contains(outer,inner):
    a,b=ends(outer);x,y=ends(inner)
    return a<=x<=y<=b


def sigma_endpoint_proofs(constants,q,kernels):
    """Check the actual primitive definitions before applying reflection."""
    graph=constants.frame.graph;proofs={}
    contract=constants.frame.bridge.leading.contracts['original_sigma']
    for key in ('J1','J_transition'):
        expr=kernels[key];name=str(expr.func)
        if not name.startswith('original_kernel_') or expr.args!=(s.Integer(1),):
            raise ValueError('Actual current sigma primitive at its endpoint required')
        node=int(name[len('original_kernel_'):]);n=graph.nodes[node]
        if n['operation']!='definite_integral' or not all(n[k] for k in
                ('exact_function_integral','original_kernel_function','numerical_value_not_installed',
                 'coordinate_and_kernel_Z_independent')) or q.at(n['lower'])!=0:
            raise ValueError('Exact original integral from zero required')
        body=graph.nodes[n['integrand']];variable=graph.nodes[body['argument']]
        if body['operation']!='exact_original_source_flat_sigma' or body['source_binding']!=contract or \
                body['definition']!='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside' or \
                not body['smooth_flat_endpoints'] or not body['Z_independent_argument'] or \
                variable!={'operation':'bound_variable','name':n['variable']}:
            raise ValueError('Same original reflected flat sigma integrand required')
        proofs[key]=dict(integral_node=node,integrand_node=n['integrand'],
            exact_endpoint='1',exact_value='1/2',same_original_sigma_contract=contract,
            proof='sigma(1-s)=1-sigma(s); reflect integral_0^1; 2*J(1)=1')
    return proofs


class CurrentOriginalRpAmplitudeBinding:
    @source_precision
    def __init__(self,before=None,constants=None,require_checked=True):
        self.before=before if before is not None else correlated.CurrentOriginalRpCorrelatedPhysicalSource()
        if type(self.before) is not correlated.CurrentOriginalRpCorrelatedPhysicalSource or not self.before.acceptance_loaded:
            raise ValueError('Accepted current correlated physical source required')
        self.before.assert_graph();self.graph=self.before.graph;self.ctx=self.before.ctx
        self.family_record=self.before.family_record;self.radius=self.before.radius
        self.raw=self.before.physical.owner.raw;self.selected=self.before.physical.owner.selected
        self.inlet=self.selected.inlet
        self.exact=constants if constants is not None else exact.CurrentOriginalRpNativeConstants(self.radius.frame)
        if type(self.exact) is not exact.CurrentOriginalRpNativeConstants or not self.exact.acceptance_loaded or \
                self.exact.frame is not self.radius.frame or self.exact.identity!=self.family_record or \
                self.inlet.exact.frame is not self.exact.frame or self.raw.graph is not self.graph:
            raise ValueError('Same current exact amplitude frame and production source required')
        self.hashes=dict(self.before.hashes)
        correlated.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,self.exact.hashes)
        for name in (correlated.NAME,correlated.RECEIPT,exact.NAME,exact.RECEIPT,Path(__file__).name):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        q,constants=self.exact.current_constants();_,kernels=self.exact.original_native_functions(q)
        self.sigma_proofs=sigma_endpoint_proofs(self.exact,q,kernels)
        Tw=q.at(self.exact.frame.functions['Tw'])
        compact=-s.Rational(57,10)-s.exp(40)/2-Tw/2-q.mu*Tw-q.mu/2
        substituted=constants['U'].subs({kernels[key]:s.Rational(1,2) for key in self.sigma_proofs})
        if not zero(substituted-s.exp(compact)):
            raise ValueError('Current exact U0 must equal the reflected-source compact exponential')
        self.symbolic_logU0=s.sstr(compact)
        self.source_U0_expression=s.sstr(constants['U'])
        g=self.graph;f=self.radius.functions
        self.Tw=correlated.box.pulse.radius.FunctionRef(g,self.exact.frame.functions['Tw'].node)
        self.exp40=g.unary('exp',g.constant(40));self.logP=f['logP'];self.mu=f['mu']
        self.logU0=g.add(g.constant('-57/10'),g.mul(g.constant('-1/2'),self.exp40),
            g.mul(g.constant('-1/2'),self.Tw),g.neg(g.mul(self.mu,self.Tw)),g.mul(g.constant('-1/2'),self.mu))
        self.proof=g.node('exact_original_Rp_U0_reflected_sigma_log_identity',
            original_U0_function=self.raw.U0.node,original_logU0_function=self.raw.logU0.node,
            exact_compact_logU0_function=self.logU0.node,
            current_native_frame_u_row=exact.frame.target.rows(
                self.exact.frame.functions['actual_identified_Rp_native_frame_C3']['u'])[0].node,
            exact_Z_substitution=g.zero.node,sigma_endpoint_proofs=self.sigma_proofs,
            identity='U0=exp(-57/10-exp(40)/2-(1/2+mu)*Tw-mu/2)>0',
            Tw_definition='-60*log(mu); original current parameter, not an endpoint')
        # Evaluate the original elementary parameter graph afresh at higher
        # precision. Copied coarse mu/Tw boxes cannot prove containment of a
        # previously rounded amplitude box merely by interval overlap.
        self.amplitude_ctx=MPIntervalContext();self.amplitude_ctx.dps=300
        reader=correlated.locator.CancelledGraphBounds(g,self.amplitude_ctx,{},
            self.before.locator.source_exponentials)
        self.logP_box=reader.at(self.logP);self.Tw_box=reader.at(self.Tw)
        self.logU0_box=reader.at(self.logU0);self.U0_box=self.amplitude_ctx.exp(self.logU0_box)
        if ends(self.U0_box)[0]<=0 or not contains(self.inlet.constants['U'],self.U0_box) or \
                not contains(self.selected.constants['U'],self.U0_box):
            raise ArithmeticError('Copied production U boxes must enclose the actual exact current U0')
        self._snapshots=dict(inlet_U=self.inlet.constants['U']._mpi_,
            selected_U=self.selected.constants['U']._mpi_,U0=self.U0_box._mpi_,
            logU0=self.logU0_box._mpi_,logP=self.logP_box._mpi_,Tw=self.Tw_box._mpi_)
        self._proof_record=json.dumps(self.graph.nodes[self.proof.node],sort_keys=True)
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN) or \
                    receipt['source_family']!=self.family_record:
                raise ValueError('Current amplitude receipt or scope differs')
            correlated.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        self.before.assert_graph();g=self.graph
        source=exact.frame.target.rows(self.exact.frame.functions['actual_identified_Rp_native_frame_C3']['u'])[0]
        actual=g.nodes[self.raw.U0.node]
        checks=dict(same_current_native_source=self.inlet.exact.frame is self.exact.frame is self.radius.frame,
            actual_current_U0_substitution=actual==dict(operation='function_substitution',expression=source.node,
                variable=g.symbol('Z').node,value=g.zero.node,Z_independent_substitution=True),
            same_raw_logU0=g.nodes[self.raw.logU0.node]==dict(operation='analytic_unary',name='log',argument=self.raw.U0.node),
            unchanged_source_identity_proof=json.dumps(g.nodes[self.proof.node],sort_keys=True)==self._proof_record,
            same_compact_log_and_parameter_nodes=g.nodes[self.proof.node]['exact_compact_logU0_function']==self.logU0.node
                and self.logP==self.radius.functions['logP'] and self.mu==self.radius.functions['mu']
                and self.Tw.node==self.exact.frame.functions['Tw'].node,
            unchanged_full_production_U_boxes=self.inlet.constants['U']._mpi_==self._snapshots['inlet_U'] and
                self.selected.constants['U']._mpi_==self._snapshots['selected_U'],
            original_current_U0_numeric_inclusion=contains(self.inlet.constants['U'],self.U0_box) and
                contains(self.selected.constants['U'],self.U0_box),
            unchanged_directed_source_bindings=all(value._mpi_==self._snapshots[key] for key,value in
                (('U0',self.U0_box),('logU0',self.logU0_box),('logP',self.logP_box),('Tw',self.Tw_box))),
            same_current_source_graph=self.raw.graph is g is self.before.graph,
            unchanged_independent_precision=self.amplitude_ctx.dps==300,
            same_shared_production_amplitude=self.selected.constants is self.selected.axial4.constants is
                self.selected.pulse.high.constants is self.selected.flatten.inlet.constants)
        if not all(checks.values()):raise ValueError('Current amplitude source changed: '+str(checks))
        return checks

    @source_precision
    def reader(self,delivery=None):
        """No caller-supplied bound or equality alias is accepted."""
        self.assert_graph()
        base=(self.before.locator.reader() if delivery is None else self.before._validate_delivery(delivery)[0])
        bindings={**base.bindings,self.raw.U0.node:self.U0_box,self.logP.node:self.logP_box,
            self.Tw.node:self.Tw_box,self.logU0.node:self.logU0_box}
        seeded=correlated.locator.CancelledGraphBounds(self.graph,self.ctx,bindings,base.allowed_exponentials)
        aliases={**getattr(base,'aliases',{}),self.raw.logU0.node:self.logU0.node}
        return correlated.CorrelatedGraphBounds(seeded,aliases)

    @source_precision
    def physical_log_scale_rows(self,delivery):
        """Actual original physical rows with evaluated positive-scale logs.

        Signed coefficient groups are retained individually. A bounded scale
        logarithm does not certify the accuracy of a sum of different scales.
        """
        reader=self.reader(delivery);mapped=self.before.physical_rows(delivery);out={}
        for section in ('Cartesian_spatial_rows','fixed_x_time_rows'):
            out[section]={}
            for component,rows in mapped[section].items():
                values=rows if section=='Cartesian_spatial_rows' else {'dt':rows}
                out[section][component]={}
                for label,row in values.items():
                    groups=[]
                    for group in row.groups():
                        logscale=self.graph.add(*[ref for _,ref in group['log_scale_parts']])
                        bound=reader.at(logscale)
                        groups.append(dict(exact_positive_scale_log_function=logscale.node,
                            directed_positive_scale_log=bound,signed_source_coefficient=group['signed_coefficient'],
                            source_scale_units=group['scale_units'],source_scale_powers=[dict(
                                numerator=v.numerator,denominator=v.denominator) for v in group['exact_source_scale_powers']],
                            original_log_parts={name:ref.node for name,ref in group['log_scale_parts']},
                            source_contributors=len(group['contributors'])))
                    out[section][component][label]=dict(physical_derivative=row.derivative,groups=groups)
        return dict(source_family=self.family_record,chart=mapped['chart'],exact_Z=mapped['exact_Z'],
            exact_native_coordinate=mapped['exact_native_coordinate'],actual_current_U0_identity=self.proof.node,
            actual_original_correlated_physical_inverse_identity=mapped['correlated_original_physical_input_identity'].node,
            actual_original_native_inverse_identity=mapped['correlated_original_native_inverse_identity'].node,
            evaluated_physical_log_scale_rows=out,actual_current_source_and_original_operators_retained=True,
            source_and_coordinate_enclosures_retained=True,distinct_signed_scales_not_summed=True,
            absolute_physical_values_and_accuracy_not_materialized=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    def report(self):
        return dict(source_family=self.family_record,actual_current_source_U0_expression=self.source_U0_expression,
            exact_compact_logU0_expression=self.symbolic_logU0,
            actual_current_amplitude_identity_proof=self.proof.node,
            original_U0_node=self.raw.U0.node,original_logU0_node=self.raw.logU0.node,
            compact_logU0_node=self.logU0.node,sigma_endpoint_identity_proofs=self.sigma_proofs,
            directed_U0=self.U0_box,directed_logU0=self.logU0_box,directed_logPstar=self.logP_box,
            actual_current_Tw=self.Tw_box,original_inlet_U=self.inlet.constants['U'],
            actual_selected_production_U=self.selected.constants['U'],source_graph_assertions=self.assert_graph(),
            independent_elementary_source_evaluation_precision=self.amplitude_ctx.dps,
            full_original_production_U_boxes_retained=True,all_source_owners_unmutated=True,
            unrelated_moment_integrals_not_needed_for_U0=True,
            P0_numeric_source_inclusion_and_physical_materialization_still_open=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


@source_precision
def run(before=None,constants=None):
    began=time.monotonic();owner=CurrentOriginalRpAmplitudeBinding(before,constants,require_checked=False)
    report=owner.report();report.update(candidate_current_U0_source_binding_constructed=True,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(correlated.report(report),indent=2)+'\n',encoding='utf8',newline='\n')
    print('Current original U0 compact source logarithm and production inclusion generated',flush=True)
    return owner


if __name__=='__main__':run()
