"""Original cutoff-local density evaluation and supported exponential factors.

Original A/N and A_Z remain formal. Only exponential factors consume the
exact support bound. Cutoff and signed-u density tuples are hulled after
their nonlinear calculation; none of the overlapping branches is added.
"""
import ast
import inspect
import json
from pathlib import Path
import textwrap
import time
import lei_ren_part1_paper_compliant_current_native_cutoff_q_cover as cutoff
import lei_ren_part1_paper_compliant_current_native_signed_u_factored_oracle as original

cover=cutoff.cover;density=cover.density;native=cutoff.native;packets=cutoff.packets
HERE,PREFIX,sha=cutoff.HERE,cutoff.PREFIX,cutoff.sha
NAME=PREFIX+'current_native_cutoff_density_oracle.json'
RECEIPT=PREFIX+'current_native_cutoff_density_oracle_check.json'
GATE='current_original_cutoff_local_supported_density_all17_declared_queries_executed'


def replace_expressions(tree,changes):
    patterns=[ast.dump(ast.parse(old,mode='eval').body) for old,new in changes];counts=[0]*len(changes)
    class Adapter(ast.NodeTransformer):
        def visit(self,node):
            for i,pattern in enumerate(patterns):
                if ast.dump(node)==pattern:
                    counts[i]+=1;return ast.copy_location(ast.parse(changes[i][1],mode='eval').body,node)
            return super().visit(node)
    tree=Adapter().visit(tree);ast.fix_missing_locations(tree)
    if counts!=[1]*len(changes):raise ValueError('Accepted original body changed; explicit AST adapter review required')
    return tree,counts


def compile_supported_density(support):
    tree=ast.parse(inspect.getsource(density.density_Z_kernels))
    tree,counts=replace_expressions(tree,[('density.factored_expm1(A)','supported_expm1(A,N)'),
        ('c.exp(phase.bounded_value(A))','supported_exponential(A,N)')])
    scope=dict(vars(density));scope.update(supported_expm1=support.supported_expm1,
        supported_exponential=support.supported_exponential)
    exec(compile(tree,'<original-density-with-source-supported-exponential-factors>','exec'),scope)
    return scope['density_Z_kernels'],counts


def compile_conditional_query(source_query,kernel):
    tree=ast.parse(textwrap.dedent(inspect.getsource(cover.NativeSignedUDensityCover.spatial_query)))
    tree.body[0].name='conditional_original_spatial_query'
    tree,counts=replace_expressions(tree,[
        ("self.first_owner.owner.query(chart,Z,geometry['raw']['coordinate'])","conditional_source_query(chart,Z,geometry['raw']['coordinate'])"),
        ("density.density_Z_kernels(E,E_Z,V,V_Z,got['values'],N)","supported_density_kernels(E,E_Z,V,V_Z,got['values'],N)")])
    scope=dict(vars(cover));scope.update(conditional_source_query=source_query,supported_density_kernels=kernel)
    exec(compile(tree,'<original-spatial-density-with-conditional-cutoff-source>','exec'),scope)
    return scope['conditional_original_spatial_query'],counts


class NativeCutoffDensityCover(cover.NativeSignedUDensityCover):
    def __init__(self,owner):
        super().__init__(owner)
        receipt=json.loads((HERE/cutoff.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(cutoff.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Checked same-original-family cutoff q jets and support identity required')
        self.qcover=cutoff.NativeCutoffQCover(self.first_owner.owner)
        self.kernel,self.kernel_adapter_counts=compile_supported_density(self.qcover)
        self.hashes={**self.hashes,**receipt['input_hashes'],cutoff.RECEIPT:sha(cutoff.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}
        self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        N=density.density.candidate_integer(N)
        if N<160:raise ValueError('Original candidate N>=160 required')
        geometry=self.first_owner.binder.query(chart,Z,coordinate,N)
        query=self.qcover.query(chart,Z,geometry['raw']['coordinate']);cells=[];records=[]
        for branch in query['branches']:
            def source_query(source_chart,source_Z,source_coordinate):
                if source_chart!=chart or source_coordinate._mpi_!=geometry['raw']['coordinate']._mpi_ or self.ctx.mpf(source_Z)._mpi_!=self.ctx.mpf(Z)._mpi_:
                    raise ValueError('Conditional q source must use the exact original chart/Z/coordinate')
                return branch['query']
            evaluator,adapter_counts=compile_conditional_query(source_query,self.kernel)
            got=evaluator(self,chart,Z,coordinate,N)
            for cell in got['cells']:
                cell['record']['original_cutoff_condition']=branch['record']['condition']
                cell['record']['original_cutoff_branch']=branch['record']['conditional_branch']
                cells.append(cell)
            records.append(dict(condition=branch['record']['condition'],branch=branch['record']['conditional_branch'],
                spatial_query_status=got['record']['status'],original_conditional_spatial_density=got['record'],
                original_source_query_AST_replacements=adapter_counts))
        enclosed=bool(cells) and all(cell['values'] is not None for cell in cells)
        record=dict(status='enclosed' if enclosed else 'requires_original_source_or_phase_refinement',chart=chart,
            source_family=self.family,candidate_N=N,original_q_slow_jet_source=query['record'],
            actual_original_radius_phase=geometry['record'],original_cutoff_branch_density_queries=records,
            original_supported_density_AST_replacements=self.kernel_adapter_counts,
            exact_original_primitive_support_argument=cutoff.primitive_identity_checks(),
            exponent_absolute_upper=self.ctx.mpf(5)/(4*N),
            original_A_over_N_formal_factor_and_A_Z_preserved=True,
            only_exponential_factor_ranges_intersected_with_original_support=True,
            nonlinear_density_precedes_cutoff_and_signed_u_branch_unions=True,
            overlapping_cutoff_or_signed_u_branches_hulled_not_added=True,
            selected_provider_query_not_whole_chart_admission=True,
            global_histories_or_actual_targets_controls_or_recursion_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,source=query['source'],geometry=geometry,qsource=query)


class NativeCutoffFactoredOracle(original.SignedUCoverFactoredOracle):
    def __init__(self,role_owner,built=None):
        super().__init__(role_owner,built)
        self.owner=NativeCutoffDensityCover(self.original_density_owner)
        self.hashes={**self.hashes,**self.owner.hashes}
        self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        frame=super().density_frame(chart=chart,Z=Z,coordinate=coordinate,N=N)
        frame.record['phase_solver_backend']='original cutoff-local q jets and signed-u first jets with source-supported exponential factors'
        frame.record['original_cutoff_q_receipt']=cutoff.RECEIPT
        frame.record['original_A_A_Z_and_nonzero_E_V_cross_terms_retained']=True
        return frame


@native.inlet.source_precision
def run(role_owner):
    began=time.monotonic();built=role_owner.build();oracle=NativeCutoffFactoredOracle(role_owner,built)
    rows={(row.get('chart'),row.get('function_role')):row for row in built['graph'].nodes if row['operation']=='original_function_graph'}
    records={};count=0;dispatches=0
    for chart,coordinate in original.POINTS.items():
        frame=oracle.density_frame(chart=chart,Z=(-1,1),coordinate=coordinate,N=2048)
        records[chart]=frame.record
        if frame.values is not None:
            count+=1
            for key in density.RATES:
                for order in ('C0','Z'):
                    row=rows.get((chart,'density_'+key+'_'+order))
                    if row is not None:oracle.dispatch_function_range(row,frame);dispatches+=1
        print('Original cutoff-local density provider:',chart,frame.record['status'],flush=True)
    if count!=17 or dispatches!=160:
        raise ArithmeticError('All17 declared original provider queries and160 density dispatches required; no completed receipt written')
    result=dict(source_family=oracle.source_family,**{GATE:True},actual_original_declared_provider_queries=records,
        original_declared_query_count=17,enclosed_original_declared_queries=count,unresolved_original_declared_queries=17-count,
        actual_original_role_dispatches=dispatches,candidate_N=2048,full_Z_interval=[-1,1],
        exact_original_primitive_support_bound=cutoff.primitive_identity_checks(),
        original_A_over_N_and_A_Z_not_replaced_by_bound=True,
        all17_whole_chart_or_continuous_route_ranges_admitted=False,full_factored_function_graph_evaluator_installed=False,
        actual_original_full_route_numerical_integrals_evaluated=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=oracle.hashes,
        scope='Cutoff-local original q jets, signed-u phase/density ranges and exact source-support exponential factors at all seventeen declared full-Z providers. Formal A/A_Z and all source cross terms remain intact; whole-chart/full-route integral targets, controls/global N/closure/recursion remain separate.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result
