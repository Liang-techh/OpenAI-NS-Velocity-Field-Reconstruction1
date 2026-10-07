"""Original branch-local q squared jets with eta collected before evaluation.

The linear q and q_Z functions remain unchanged. Five quadratic sites in
the original first-jet body consume the equivalent direct q squared jet.
This improves arithmetic coordinates, not the defining field or cutoff.
"""
import ast
import hashlib
import inspect
import json
import math
from pathlib import Path
import textwrap
import time
import lei_ren_part1_paper_compliant_current_native_bounded_primitive_transport as accepted

density=accepted.density;cutoff=accepted.cutoff;cover=cutoff.cover;first=cover.first
current=cutoff.current;prior=cutoff.prior;native=cutoff.native;packets=cutoff.packets
HERE,PREFIX,sha=cutoff.HERE,cutoff.PREFIX,cutoff.sha;ep=cutoff.ep;ZERO=cutoff.ZERO;ORDERS=cutoff.ORDERS
NAME=PREFIX+'current_native_collected_q2_transport.json'
RECEIPT=PREFIX+'current_native_collected_q2_transport_check.json'
GATE='current_original_collected_q_squared_ordinary_jets_and_fixed_N_transport_executed'


def conditional_q2_jet(roots,eta_log,log_a_lower,branch):
    """Ordinary y<=2/Z<=1 jet of the original q^2, without squaring q jets."""
    a=roots['a'];c=a[ZERO].ctx;scalar=a[ZERO].scalar;zero=scalar(0)
    if branch['name']=='flat':
        return {order:scalar(0) for order in ORDERS}
    eta=prior.ScaledEnclosure(prior.FormalScale(a[ZERO].scale.bases,offset=eta_log),1,a[ZERO].ledger)
    Delta=dict(roots['kappa_minus2']);Delta[ZERO]=branch['Delta']
    if branch['name']=='negative':
        numerator={order:(eta*2 if order==ZERO else zero)-value for order,value in Delta.items()}
        numerator[ZERO]=numerator[ZERO].positive_intersection(eta_log+c.ln(2))
    elif branch['name']=='transition':
        theta=branch['theta'];s=[value*math.factorial(n) for n,value in enumerate(prior.sigma_jets(c,1-theta))]
        # H(theta)=sigma(1-theta)^2*(2-theta). These are ordinary
        # derivatives with respect to theta, including its minus signs.
        G=[s[0]*s[0],-2*s[0]*s[1],2*s[1]*s[1]+2*s[0]*s[2],
            -6*s[1]*s[2]-2*s[0]*s[3]]
        H=[scalar((2-theta)*G[n]-(n*G[n-1] if n else 0)) for n in range(4)]
        dy,dyy,dz,dyz,dyyz=[Delta[order] for order in ((1,0),(2,0),(0,1),(1,1),(2,1))]
        inverse_eta=lambda value:value.positive_divide(eta,eta_log)
        eta2=prior.ScaledEnclosure(prior.FormalScale(a[ZERO].scale.bases,offset=eta_log*2),1,a[ZERO].ledger)
        numerator={ZERO:eta*H[0],(1,0):H[1]*dy,(0,1):H[1]*dz,
            (2,0):H[1]*dyy+inverse_eta(H[2]*current.square(dy)),
            (1,1):H[1]*dyz+inverse_eta(H[2]*dy*dz),
            (2,1):H[1]*dyyz+inverse_eta(H[2]*(dyy*dz+dy*dyz*2))
                +(H[3]*current.square(dy)*dz).positive_divide(eta2,eta_log*2)}
    else:raise ValueError('Original negative/transition/flat cutoff branch required')
    return current.quotient_jet(numerator,{order:value*2 for order,value in a.items()},log_a_lower+c.ln(2))


def q2_record(rows,branch):
    return dict(original_cutoff_branch=branch['name'],original_cutoff_condition=branch['condition'],
        ordinary_q_squared_rows={'y%d_Z%d'%order:value.record() for order,value in rows.items()},
        transition_first_eta_over_eta_cancelled_before_interval_evaluation=branch['name']=='transition',
        second_and_third_transition_inverse_eta_factors_retained=True,
        original_a_positive_denominator_only=True,q_squared_positive_lower_not_required=True,
        sigma_Taylor_coefficients_converted_to_ordinary_derivatives=True,
        original_linear_q_and_q_derivative_objects_retained=True,
        branch_range_not_field_selection=True)


class NativeCollectedQ2Cover(cutoff.NativeCutoffQCover):
    @native.inlet.source_precision
    def query(self,chart,Z,coordinate):
        query=super().query(chart,Z,coordinate)
        positive=self.owner.owner.owner.decode(self.owner.owner.owner.inventory[chart]['actual_positive_denominator_theorem'])
        candidates,_=cutoff.conditional_cutoff_branches(query['source']['roots']['kappa_minus2'][ZERO],self.eta_log)
        conditional={row['name']:row for row in candidates}
        for item in query['branches']:
            original=item['query'];name=item['record']['conditional_branch']
            branch=dict(name=name,condition=item['record']['condition'],Delta=original['source']['roots']['kappa_minus2'][ZERO])
            if name=='transition':
                # Repeat the original range restriction from the unconditioned
                # parent. Do not divide two copies of the huge eta log again.
                branch['theta']=conditional[name]['theta']
            rows=conditional_q2_jet(original['source']['roots'],self.eta_log,positive['log_actual_a_positive_lower'],branch)
            item['query']=dict(original,q2_rows=rows)
            item['record']['original_direct_q_squared_jet']=q2_record(rows,branch)
        query['record']['original_branch_local_q_squared_jets']=[item['record']['original_direct_q_squared_jet'] for item in query['branches']]
        query['record']['original_q_and_six_ordinary_q_rows_unchanged']=True
        return query


def compile_collected_first_jets(source,qrows,q2rows,dstar_log,phi,branch):
    if source['q'] is not qrows[ZERO]:raise ValueError('Original branch q/jet C0 identity required')
    for value in q2rows.values():source['q'].coerce(value)
    if source['q'].zero and any(not value.zero for value in q2rows.values()):
        raise ValueError('Flat original q requires exact zero q squared jets')
    tree=ast.parse(inspect.getsource(first.conditioned_first_jets));original=ast.dump(tree)
    changes=[('phase.ConditionedPhase(source, dstar_log)','branch_factory(source,dstar_log)'),
        ('-loop.u * loop.hinv * loop.hinv * loop.hinv * u.d[key]','-loop.uhinvcube*u.d[key]'),
        ('loop.hinv * loop.hinv * loop.hinv * u.d[key]','loop.hinvcube*u.d[key]'),
        ('q.square()','original_q2_dual')]
    patterns=[ast.dump(ast.parse(old,mode='eval').body) for old,new in changes];counts=[0]*len(changes)
    class Adapter(ast.NodeTransformer):
        def visit(self,node):
            for i,pattern in enumerate(patterns):
                if ast.dump(node)==pattern:
                    counts[i]+=1;return ast.copy_location(ast.parse(changes[i][1],mode='eval').body,node)
            return super().visit(node)
    tree=Adapter().visit(tree);ast.fix_missing_locations(tree)
    if counts!=[1,1,1,5]:raise ValueError('Original first-jet body changed; reviewed five-site q squared adapter required')
    scope=dict(vars(first));scope.update(branch_factory=lambda src,log:cover.BranchConditionedPhase(src,log,branch),
        original_q2_dual=first.source_dual(q2rows[ZERO],q2rows))
    exec(compile(tree,'<original-first-jets-with-branch-local-direct-q-squared>','exec'),scope)
    got=scope['conditioned_first_jets'](source,qrows,dstar_log,phi)
    got['record'].update(signed_u_branch=branch['name'],original_direct_q_squared_first_jet_adapter=dict(
        original_body_sha256=hashlib.sha256(original.encode()).hexdigest(),AST_replacement_counts=counts,
        original_linear_q_q_Z_u_h_and_implicit_phase_rules_unchanged=True,
        direct_q_squared_rows_from_same_cutoff_branch=True,field_or_derivative_cap_not_substituted=True))
    return got


def compile_collected_query(source_query,kernel):
    tree=ast.parse(textwrap.dedent(inspect.getsource(cover.NativeSignedUDensityCover.spatial_query)))
    tree.body[0].name='collected_original_spatial_query'
    tree,counts=density.replace_expressions(tree,[
        ("self.first_owner.owner.query(chart,Z,geometry['raw']['coordinate'])","conditional_source_query(chart,Z,geometry['raw']['coordinate'])"),
        ("branch_first_jets(source,qsource['rows'],self.first_owner.dstar,phi,branch)","compile_collected_first_jets(source,qsource['rows'],qsource['q2_rows'],self.first_owner.dstar,phi,branch)"),
        ("density.density_Z_kernels(E,E_Z,V,V_Z,got['values'],N)","supported_density_kernels(E,E_Z,V,V_Z,got['values'],N)")])
    scope=dict(vars(cover));scope.update(conditional_source_query=source_query,supported_density_kernels=kernel,
        compile_collected_first_jets=compile_collected_first_jets)
    exec(compile(tree,'<original-spatial-density-with-collected-q-squared-first-jets>','exec'),scope)
    return scope['collected_original_spatial_query'],counts


class NativeCollectedQ2Density(accepted.NativeBoundedPrimitiveDensity):
    def __init__(self,owner):
        super().__init__(owner);self.qcover=NativeCollectedQ2Cover(self.first_owner.owner)
        self.hashes={**self.hashes,Path(__file__).name:sha(Path(__file__).name)};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        tree=ast.parse(textwrap.dedent(inspect.getsource(density.NativeCutoffDensityCover.spatial_query)))
        tree.body[0].name='collected_cutoff_spatial_query'
        tree,counts=density.replace_expressions(tree,[
            ('compile_conditional_query(source_query,self.kernel)','compile_collected_query(source_query,self.kernel)')])
        scope=dict(vars(density));scope['compile_collected_query']=compile_collected_query
        exec(compile(tree,'<original-cutoff-spatial-query-with-collected-q-squared>','exec'),scope)
        start=len(self.support_rows);got=scope['collected_cutoff_spatial_query'](self,chart,Z,coordinate,N)
        rows=self.support_rows[start:];got['record'].update(
            original_C0_A_B_source_support_ranges_before_nonlinear_density=rows,
            original_A_Z_B_Z_are_derived_from_same_original_field_with_direct_q_squared=True,
            original_linear_q_derivatives_retained=True,original_spatial_query_AST_replacements=counts,
            exact_original_C0_primitive_support_theorem=accepted.primitive_support_theorem())
        got['record']['original_A_over_N_formal_factor_and_A_Z_preserved']=all(row['A']['original_tighter_formal_range_retained'] for row in rows)
        return got


class NativeCollectedQ2Oracle(accepted.NativeBoundedPrimitiveOracle):
    def __init__(self,role_owner,built=None):
        super().__init__(role_owner,built);self.owner=NativeCollectedQ2Density(self.original_density_owner)
        self.hashes={**self.hashes,**self.owner.hashes};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        frame=super().density_frame(chart=chart,Z=Z,coordinate=coordinate,N=N)
        frame.record['phase_solver_backend']='original cutoff/signed-u phase body with branch-local direct q squared jets and C0 primitive supports'
        return frame


class NativeCollectedQ2Transport(accepted.NativeBoundedPrimitiveTransport):
    def __init__(self,role_owner):
        super().__init__(role_owner)
        checked=json.loads((HERE/accepted.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(accepted.GATE) or checked['source_family']!=self.family:
            raise ValueError('Checked original C0 primitive support baseline required')
        self.oracle=NativeCollectedQ2Oracle(role_owner,self.built)
        self.hashes={**self.hashes,**self.oracle.hashes,**checked['input_hashes'],accepted.RECEIPT:sha(accepted.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)};self.service.bind_hashes(self.hashes)


@native.inlet.source_precision
def run(role_owner,*,return_live=False):
    began=time.monotonic();owner=NativeCollectedQ2Transport(role_owner);got=owner.route()
    if got['history'] is None or len(got['cells'])!=24:raise ArithmeticError('All24 original continuous cells required')
    baseline=json.loads((HERE/accepted.NAME).read_bytes());comparison=accepted.compare_targets(owner.ctx,got['record'],baseline)
    result=dict(**got['record'],**{GATE:True},ordinary_q_squared_orders=[list(order) for order in ORDERS],
        comparison_with_checked_original_C0_support_target_ranges=comparison,
        strict_target_absolute_upper_reductions=sum(row['strict_absolute_upper_reduction'] for row in comparison.values()),
        defining_original_q_and_linear_q_Z_unchanged=True,only_five_original_quadratic_first_jet_sites_adapted=True,
        useful_repair_contraction_or_actual_controls_established=False,
        execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Original branch-local q squared ordinary y2/Z1 jet with exact first-order eta cancellation; five equivalent quadratic first-jet sites; same full24 fixed-N original C1 transport. Genuine linear q_Z and higher inverse eta sensitivity remain. No actual controls/global N/closure/recursion admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original collected q squared transported through24 cells;strict target upper reductions:',result['strict_target_absolute_upper_reductions'],flush=True)
    return (result,owner,got) if return_live else result
