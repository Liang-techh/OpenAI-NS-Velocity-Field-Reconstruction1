"""Directed original delta/logC products with unchanged physical functions.

Source delta is exp(-4*logP-30), not a pulse or transport decay atom.
Its min branch is proved before recomputation. Mixed origins retain their
nonzero uncertainty; the exact binary-rational fixed origin is unchanged.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rp_remaining_pressure_tail as tail
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

centered,correlated,signed,box=tail.centered,tail.correlated,tail.signed,tail.box
HERE,PREFIX,sha,ends=tail.HERE,tail.PREFIX,tail.sha,tail.ends
NAME=PREFIX+'current_original_Rp_source_product_arithmetic.json.gz'
RECEIPT=PREFIX+'current_original_Rp_source_product_arithmetic_check.json'
GATES=('current_original_Rp_actual_delta_source_min_branch_bound',
       'current_original_Rp_same_source_mixed_parameter_products_refined',
       'current_original_Rp_refined_physical_relative_scale_width_installed')
OPEN=tail.OPEN


def contains(outer,inner):
    a,b=ends(outer);u,v=ends(inner)
    return a<=u<=v<=b


class CurrentOriginalRpSourceProductArithmetic:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not tail.CurrentOriginalRpRemainingPressureTail or not before.acceptance_loaded:
            raise ValueError('Accepted original remaining-pressure source required')
        self.before,self.amplitude=before,before.amplitude
        self.graph,self.family_record=before.graph,before.family_record
        self.ctx=MPIntervalContext();self.ctx.dps=700
        self.pivot=before.before.pivot
        self.exact_singleton=before.before.exact_singleton
        self.delta=self.amplitude.radius.functions['delta']
        self.logdelta=self.amplitude.radius.functions['logdelta']
        self.logmu=self.amplitude.radius.functions['logmu']
        self.hashes=dict(before.hashes)
        for filename in (tail.NAME,tail.RECEIPT,Path(__file__).name):self.hashes[filename]=sha(filename)
        self.acceptance_loaded=False;self._readers={}
        g=self.graph
        bare=correlated.locator.CancelledGraphBounds(g,self.ctx,{},())
        exp40=g.unary('exp',g.constant(40))
        expected_logP=g.add(exp40,g.constant(11))
        expected_logdelta=g.add(g.mul(g.constant(-4),self.amplitude.logP),g.constant(-30))
        expected_logmu=g.add(g.unary('log',g.constant('1/1000')),
            g.mul(g.constant(-4),self.amplitude.logP))
        mu_node=g.nodes[self.amplitude.mu.node];delta_node=g.nodes[self.delta.node]
        if (mu_node.get('operation'),mu_node.get('name'))!=('analytic_unary','exp') or \
                (delta_node.get('operation'),delta_node.get('name'))!=('analytic_unary','exp'):
            raise ValueError('Actual original mu and delta exponential definitions required')
        for actual,expected in ((self.amplitude.logP,expected_logP),(self.logdelta,expected_logdelta),
                                (mu_node['argument'],expected_logmu)):
            if not hasattr(actual,'node'):actual=box.pulse.radius.FunctionRef(g,actual)
            if bare.polynomial(g.sub(actual,expected)):
                raise ValueError('Original exp(40)+11 / mu / delta source hierarchy differs')
        if delta_node['argument']!=self.logdelta.node:
            raise ValueError('Delta atom must use its actual original logdelta child DAG')
        self._source_nodes={node:copy.deepcopy(g.nodes[node]) for node in
            (self.delta.node,self.logdelta.node,self.amplitude.mu.node,self.amplitude.logP.node)}
        self.native_definition=copy.deepcopy(before.terminal.raw.flat.inlet.datum.definition)
        if self.native_definition['delta']!='min(1e-200,exp(-4logPstar-30))' or \
                self.native_definition['logPstar']!='exp(Md)+11' or self.native_definition['Md']!='40':
            raise ValueError('Accepted actual Md=40 native delta min definition required')
        original=self.amplitude.reader()
        self._allowed=original.allowed_exponentials
        # Only the defining original delta is additionally admitted, locally
        # for this parameter computation; physical exp guards stay unchanged.
        pure=correlated.locator.CancelledGraphBounds(g,self.ctx,{},self._allowed|{self.delta.node})
        self.refs=dict(mu=self.amplitude.mu,logP=self.amplitude.logP,Tw=self.amplitude.Tw,
            logU0=self.amplitude.logU0,delta=self.delta,logmu=self.logmu,logdelta=self.logdelta)
        native_delta=before.terminal.raw.flat.delta
        self.finite_bindings={};self.inclusion={}
        for name,ref in self.refs.items():
            value=pure.at(ref)
            old=native_delta if name=='delta' else original.at(ref)
            if not contains(old,value):raise ArithmeticError('Fresh original source must remain inside '+name)
            self.finite_bindings[ref.node]=value
            self.inclusion[name]=dict(source_function=ref.node,original_box=old,fresh_box=value,
                original_defining_graph_recomputed=True,included_in_original_box=True)
        U=self.ctx.exp(self.finite_bindings[self.amplitude.logU0.node])
        if not contains(original.at(self.amplitude.raw.U0),U):
            raise ArithmeticError('Fresh exact U0 must stay inside its accepted source box')
        self.finite_bindings[self.amplitude.raw.U0.node]=U
        self.inclusion['U0']=dict(source_function=self.amplitude.raw.U0.node,
            original_box=original.at(self.amplitude.raw.U0),fresh_box=U,
            original_compact_identity_recomputed=True,included_in_original_box=True)
        if not contains(native_delta,self.finite_bindings[self.delta.node]):
            raise ArithmeticError('Fresh delta must stay inside the actual native delta source')
        threshold=-200*self.ctx.ln(10)
        if ends(self.finite_bindings[self.logdelta.node])[1]>=ends(threshold)[0]:
            raise ValueError('Original exponential delta min branch must be strict')
        self.native_delta_snapshot=native_delta._mpi_
        self._snapshots={node:value._mpi_ for node,value in self.finite_bindings.items()}
        self.mixed_product=g.mul(self.pivot,self.delta)
        self.ratio=g.quotient(self.delta,self.amplitude.mu,'original positive source mu')
        self.ratio_constant=g.mul(g.constant(1000),g.unary('exp',g.constant(-30)))
        self.proof=g.node('current_original_Rp_source_delta_product_arithmetic',
            actual_delta_function=self.delta.node,actual_logdelta_function=self.logdelta.node,
            actual_mu_function=self.amplitude.mu.node,actual_logP_function=self.amplitude.logP.node,
            source_parameter_identities='logP=exp(40)+11; delta=exp(-4logP-30); mu=exp(-4logP)/1000',
            strict_actual_native_min_branch='-4logP-30 < -200log(10)',
            same_exact_logC=self.pivot.node,exact_mixed_product=self.mixed_product.node,
            mixed_product_is_logC_times_delta_not_C_times_delta=True,
            exact_delta_mu_ratio=self.ratio.node,exact_ratio_constant=self.ratio_constant.node,
            ratio_identity='delta/mu=1000*exp(-30), real positive source exponentials',
            directed_decimal_precision=700,refined_source_nodes=sorted(self.finite_bindings),
            original_physical_exponential_guard=1000,
            nonconstant_mixed_origin_has_directed_uncertainty=True,
            original_parameters_coordinates_P0_velocity_pressure_functions_unchanged=True)
        self._proof_snapshot=copy.deepcopy(g.nodes[self.proof.node])
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family']!=self.family_record:
                raise ValueError('Original source-product receipt/scope differs')
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        original=self.amplitude.reader()
        checks=dict(accepted_original_remaining_pressure_source=self.before.acceptance_loaded,
            same_original_source_defining_nodes=all(self.graph.nodes[node]==value for node,value in self._source_nodes.items()),
            same_native_delta_definition=self.before.terminal.raw.flat.inlet.datum.definition==self.native_definition,
            same_actual_native_delta_box=self.before.terminal.raw.flat.delta._mpi_==self.native_delta_snapshot,
            unchanged_fresh_parameter_boxes=set(self.finite_bindings)==set(self._snapshots)
                and all(value._mpi_==self._snapshots[node] for node,value in self.finite_bindings.items()),
            same_exact_singleton_origin=original.bindings[self.pivot.node]._mpi_==self.exact_singleton,
            unchanged_parameter_product_proof=self.graph.nodes[self.proof.node]==self._proof_snapshot,
            bounded_source_precision=self.ctx.dps==700,
            unchanged_original_source_exponential_whitelist=original.allowed_exponentials==self._allowed,
            unchanged_physical_exponential_guard=signed.EXP_LOG_LIMIT==1000)
        if not all(checks.values()):raise ValueError('Original source-product arithmetic differs: '+str(checks))
        return checks

    @source_precision
    def reader(self,delivery):
        self.assert_graph()
        original=self.before.before.reader(delivery)
        reader=correlated.CorrelatedGraphBounds(correlated.locator.CancelledGraphBounds(
            self.graph,self.ctx,{**original.bindings,**self.finite_bindings},original.allowed_exponentials),original.aliases)
        self._readers[id(reader)]=(reader,{node:value._mpi_ for node,value in reader.bindings.items()},
            dict(reader.aliases),reader.allowed_exponentials)
        return reader

    @source_precision
    def center(self,reader,log_function):
        return centered.CurrentOriginalRpCenteredScaleArithmetic.center(self,reader,log_function)

    @source_precision
    def evaluate(self,delivery,relative_width_target='1/100000000'):
        self.before.source(delivery)
        proxy=box._Proxy(self.before.before)
        sourceproxy=box._Proxy(self.before.refined)
        sourceproxy.physical_rows=self.before.physical_rows
        proxy.before=sourceproxy;proxy.ctx=self.ctx
        proxy.reader=self.reader;proxy.center=self.center;proxy.assert_graph=self.assert_graph
        proxy.parameter_proof=self.proof;proxy.acceptance_loaded=self.acceptance_loaded
        result=centered.CurrentOriginalRpCenteredScaleArithmetic.evaluate(proxy,delivery,relative_width_target)
        result.update(original_source_product_arithmetic_proof=self.proof.node,
            exact_original_pressure_tail_identity=self.before.proof.node,
            source_coefficient_provider='accepted_refined_velocity_and_remaining_pressure_with_same_source_parameter_arithmetic',
            mixed_origin_uncertainty_not_discarded=True,
            parameter_product_directed_precision=700,
            **dict.fromkeys(tail.GATES,self.before.acceptance_loaded),
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        return result

    @source_precision
    def velocity_pressure(self,delivery,relative_width_target='1/100000000'):
        _,_,request=self.amplitude.before._validate_delivery(delivery)
        result=self.evaluate(delivery,relative_width_target)
        rows=result['physical_value_rows']['Cartesian_spatial_rows']
        return dict(values={name:rows[component]['x0_y0_z0'] for name,component in
            (('u','ux'),('v','uy'),('w','uz'),('p','p'))},
            physical_coordinates={name:ref.node for name,ref in request.forward_coordinates.items()},
            source_family=self.family_record,original_P0_source_binding=result['original_P0_source_binding'],
            original_source_product_arithmetic_proof=self.proof.node,
            exact_original_pressure_tail_identity=self.before.proof.node,
            full_certified_physical_accuracy=False,unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


@source_precision
def run(before,deliveries,relative_width_target='1/1000'):
    began=time.monotonic();owner=CurrentOriginalRpSourceProductArithmetic(before,require_checked=False)
    views={name:owner.evaluate(delivery,relative_width_target) for name,delivery in deliveries.items()}
    result=dict(source_family=owner.family_record,actual_source_product_views=views,
        source_graph_assertions=owner.assert_graph(),fresh_source_parameter_inclusion=owner.inclusion,
        source_product_arithmetic_proof=owner.graph.nodes[owner.proof.node],input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,views
