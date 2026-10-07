"""Overlapping signed-u covers for the original spatial C0/Z density source.

Only the geometric C0 range of u is conditional. All original root/q jets,
radius phase, velocity terms, units and directed ledger remain unchanged.
Overlapping branches are hulled after nonlinear density evaluation, not added.
"""
import ast
import copy
import inspect
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals as density

first=density.first;phase=first.phase;current=first.current;prior=first.prior
native=first.native;packets=first.packets;ep=first.ep;ZERO=first.ZERO
HERE,PREFIX,sha=first.HERE,first.PREFIX,first.sha
NAME=PREFIX+'current_native_signed_u_phase_cover.json'
RECEIPT=PREFIX+'current_native_signed_u_phase_cover_check.json'
GATE='current_original_overlapping_signed_u_spatial_density_C0_Z_cover_executed'


def unit_hull(value):
    """Directed range intersection, never a source point or a zero tail."""
    return phase.clipped(value.ctx,value.coefficient*value.bounded_exp(value.scale.evaluate()),0,1)


def bind_source_context(query,c):
    """Outward-copy original ranges into one context without choosing values."""
    source=query['source'];q=source['q'];old_bases=q.scale.bases
    original_values=[q]+[value for rows in source['roots'].values() for value in rows.values()]
    if query['rows'] is not None:original_values+=list(query['rows'].values())
    if any(value.scale.bases is not old_bases or value.ledger is not q.ledger or value.ctx is not q.ctx for value in original_values):
        raise ValueError('Original root/q ranges must share one basis, context and ledger')
    if q.ctx is c:return query
    copied=lambda value:c.mpf(ep(value))
    bases=tuple(copied(value) for value in old_bases)
    def rebind(value):
        if value.scale.bases is not old_bases or value.ledger is not q.ledger:
            raise ValueError('Original root/q ranges must share a source basis and ledger')
        scale=prior.FormalScale(bases,value.scale.powers,copied(value.scale.offset))
        result=prior.ScaledEnclosure(scale,copied(value.coefficient),value.ledger)
        if result.ctx is not c:raise ValueError('Directed source context copy failed')
        return result
    before=q.ledger.get('coefficient_log_coordinate_rescalings',0)
    roots={name:{order:rebind(value) for order,value in rows.items()} for name,rows in source['roots'].items()}
    rows=None if query['rows'] is None else {order:rebind(value) for order,value in query['rows'].items()}
    rebound_q=rebind(q)
    record=dict(query['record'],original_source_ranges_directed_context_bridge=True,
        source_log_values_and_exact_factor_powers_unchanged=True,original_source_ledger_unchanged=True,
        coefficient_and_log_endpoints_outward_copied=True,source_function_values_not_selected=True,
        context_bridge_coefficient_log_coordinate_rescalings=q.ledger.get('coefficient_log_coordinate_rescalings',0)-before,
        conservative_range_conversion_not_added_source_correlations=True,
        packet_metadata_kept_separate_from_rebound_density_arithmetic=True)
    return dict(query,source=dict(source,q=rebound_q,roots=roots),rows=rows,record=record)


def signed_u_branches(source,dstar_log):
    q=source['q'];c=q.ctx
    dstar=prior.ScaledEnclosure(prior.FormalScale(q.scale.bases,offset=dstar_log),1,q.ledger)
    u=(source['roots']['p2'][ZERO]*q).positive_divide(dstar,dstar_log)
    if q.zero:return u,[dict(name='flat_q',u=u,condition='q identically zero')],[]
    lo,hi=ep(u.coefficient);logs=ep(u.scale.evaluate())
    threshold=ep(c.ln(c.mpf('.125')))[0]
    central_upper=ep(c.ln(c.mpf('.25')))[1]
    branches=[];empty=[]
    for name,sign,magnitude,lower_magnitude in (
        ('negative_tail',-1,max(-lo,0),-hi if hi<0 else 0),
        ('positive_tail',1,max(hi,0),lo if lo>0 else 0)):
        upper=logs[1]+ep(c.ln(c.mpf(magnitude)))[1] if magnitude else -mp.inf
        if upper<threshold:
            empty.append(dict(name=name,proof='original signed range has no magnitude >=1/8'));continue
        lower=logs[0]+ep(c.ln(c.mpf(lower_magnitude)))[0] if lower_magnitude else -mp.inf
        L=c.mpf((max(lower,threshold),upper))
        conditional=prior.ScaledEnclosure(prior.FormalScale(u.scale.bases,offset=L),sign,u.ledger)
        branches.append(dict(name=name,sign=sign,u=conditional,log_magnitude=L,
            condition='u<=-1/8' if sign<0 else 'u>=1/8'))
    fixed_sign_lower=0 if lo<=0<=hi else min(abs(lo),abs(hi))
    lower=logs[0]+ep(c.ln(c.mpf(fixed_sign_lower)))[0] if fixed_sign_lower else -mp.inf
    if lower>central_upper:
        empty.append(dict(name='central',proof='original signed magnitude strictly greater than1/4'))
    else:
        # An outer C0 range of the restricted original function, not a new u.
        interval=c.mpf((0 if lo>=0 else '-.25',0 if hi<=0 else '.25'))
        branches.append(dict(name='central',u=u.scalar(interval),condition='|u|<=1/4'))
    if not branches:raise ArithmeticError('Overlapping u cover unexpectedly empty')
    return u,branches,empty


class BranchConditionedPhase(phase.ConditionedPhase):
    def __init__(self,source,dstar_log,branch):
        self.source=source;self.q=source['q'];self.roots=source['roots'];self.c=c=self.q.ctx
        self.point=mp.mp.clone();self.point.dps=c.dps+20
        self.a=self.roots['a'][ZERO];self.t0=self.roots['t0'][ZERO];self.E=self.roots['E'][ZERO]
        self.scalar=self.q.scalar
        self.dstar=prior.ScaledEnclosure(prior.FormalScale(self.q.scale.bases,offset=dstar_log),1,self.q.ledger)
        self.u=branch['u'];self.branch=branch
        self.nu=self.scalar(1)+current.square(self.t0)+current.square(self.q)*2
        self.flat=self.q.zero
        if self.flat:self.geometry='flat';return
        if branch['name']=='central':
            self.geometry='small_r_series';uf=phase.bounded_value(self.u);h=c.sqrt(1+uf**2)
            self.r=uf/h;self.s=1/(1+uf**2);self.hinv=self.scalar(1/h)
            self.rho=self.scalar(1/(h*(h+abs(uf))))
            self.uhinvcube=self.u*self.hinv*self.hinv*self.hinv
            self.hinvcube=self.hinv*self.hinv*self.hinv
        else:
            self.geometry='signed_Mobius';self.sign=branch['sign']
            L=self.u.scale;z=self.q.bounded_exp((-L-L).evaluate());H=c.sqrt(1+z)
            self.hinv=prior.ScaledEnclosure(-L,1/H,self.q.ledger)
            self.rho=prior.ScaledEnclosure(-L-L,1/(H*(H+1)),self.q.ledger)
            self.s_source=prior.ScaledEnclosure(-L-L,1/(1+z),self.q.ledger)
            self.r=phase.clipped(c,1-unit_hull(self.rho),c.mpf('.125')/c.sqrt(1+c.mpf('.125')**2),1)*self.sign
            self.s=unit_hull(self.s_source)
            # Collect dependent L factors before interval evaluation.
            # u*h^-3=sign*exp(-2L)/(1+exp(-2L))^(3/2).
            self.uhinvcube=prior.ScaledEnclosure(-L-L,self.sign/(H**3),self.q.ledger)
            self.hinvcube=prior.ScaledEnclosure(-L-L-L,1/(H**3),self.q.ledger)
        self.hinv_finite=unit_hull(self.hinv)
        self.n0=self.normalized(self.scalar(1),1,True)
        self.nt=self.normalized(current.square(self.t0),1,True)
        self.nq=self.normalized(current.square(self.q),c.mpf('.5'),True)
        self.ntq=self.normalized(self.t0*self.q,1/(2*c.sqrt(2)),False)

    def normalized(self,value,cap,positive):
        # The inherited theorem fallback reads directed endpoints of cap.
        # Convert literal caps before it is needed on wider chart sources.
        return super().normalized(value,self.c.mpf(cap),positive)

    def angle_endpoint(self,x,inverse=False):
        c=self.c
        if x in (0,c.mpf('.5'),1):return c.mpf(x)
        reflect=x>c.mpf('.5');argument=c.mpf(x)
        if reflect:argument=1-argument
        rho=unit_hull(self.rho);signed=self.sign*(-1 if inverse else 1)
        plus,minus=(2-rho,rho) if signed>0 else (rho,2-rho)
        angle=phase.clipped(c,c.atan2(plus*c.sin(c.pi*argument),minus*c.cos(c.pi*argument))/c.pi,0,1)
        return 1-angle if reflect else angle

    def geometry_record(self):
        record=super().geometry_record()
        record.update(conditional_branch=self.branch['name'],conditional_original_u_condition=self.branch['condition'],
            branch_is_function_range_restriction_not_source_selection=True,
            original_root_and_q_derivative_rows_unchanged=True,branch_boundary_derivatives_not_introduced=True)
        if not self.flat:
            record.update(original_signed_u_times_h_inverse_cubed=self.uhinvcube.record(),
                h_inverse_cubed=self.hinvcube.record(),dependent_log_factors_collected_before_evaluation=True)
        return record


def branch_first_jets(source,qrows,dstar_log,phi,branch):
    """Explicit AST adapter: constructor and two equivalent collected factors."""
    tree=ast.parse(inspect.getsource(first.conditioned_first_jets));original=copy.deepcopy(tree)
    expected=[
        'phase.ConditionedPhase(source, dstar_log)',
        '-loop.u * loop.hinv * loop.hinv * loop.hinv * u.d[key]',
        'loop.hinv * loop.hinv * loop.hinv * u.d[key]']
    replacement=['branch_factory(source, dstar_log)','-loop.uhinvcube * u.d[key]','loop.hinvcube * u.d[key]']
    counts=[0,0,0];patterns=[ast.dump(ast.parse(x,mode='eval').body) for x in expected]
    class Adapter(ast.NodeTransformer):
        def visit(self,node):
            for i,pattern in enumerate(patterns):
                if ast.dump(node)==pattern:
                    counts[i]+=1;return ast.copy_location(ast.parse(replacement[i],mode='eval').body,node)
            return super().visit(node)
    tree=Adapter().visit(tree);ast.fix_missing_locations(tree)
    if counts!=[1,1,1]:raise ValueError('Accepted original first-jet body changed; explicit adapter review required')
    scope=dict(vars(first));scope['branch_factory']=lambda src,log:BranchConditionedPhase(src,log,branch)
    exec(compile(tree,'<original-first-jets-with-signed-u-cover>','exec'),scope)
    got=scope['conditioned_first_jets'](source,qrows,dstar_log,phi)
    got['record']['signed_u_branch']=branch['name']
    got['record']['original_first_jet_AST_adapter']=dict(constructor_changes=counts[0],equivalent_collected_geometry_factors=counts[1:],
        original_body_sha256=__import__('hashlib').sha256(ast.dump(original).encode()).hexdigest(),
        original_source_roots_derivatives_and_implicit_rules_unchanged=True)
    return got


class NativeSignedUDensityCover:
    def __init__(self,owner):
        if type(owner) is not density.NativeDensityC1LocalIntegrals:
            raise TypeError('Accepted original density owner required')
        receipt=json.loads((HERE/density.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(density.GATE) or receipt['source_family']!=owner.family:
            raise ValueError('Accepted same-family original density source required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.first_owner=owner.owner
        self.hashes={**receipt['input_hashes'],density.RECEIPT:sha(density.RECEIPT),Path(__file__).name:sha(Path(__file__).name)}
        self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        N=density.density.candidate_integer(N)
        if N<160:raise ValueError('Original integer N>=160 required')
        geometry=self.first_owner.binder.query(chart,Z,coordinate,N)
        qsource=self.first_owner.owner.query(chart,Z,geometry['raw']['coordinate'])
        qsource=bind_source_context(qsource,self.ctx)
        source=qsource['source'];packet=source['packet']
        original_u,branches,empty=signed_u_branches(source,self.first_owner.dstar)
        E=source['roots']['E'][ZERO];E_Z=source['roots']['E'][(0,1)]
        signed_owner=self.first_owner.owner.owner.owner
        def leaf(k):
            row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
            return prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row})
        V,V_Z=[signed_owner.leaf(leaf(k),E.scale.bases,E.ledger) for k in (0,1)]
        for value in (E_Z,V,V_Z):
            E.coerce(value)
            if value.ctx is not self.ctx:raise ValueError('Original E/V density context bridge incomplete')
        cells=[]
        for phi in geometry['phase_boxes']:
            local=[]
            for branch in branches:
                try:
                    got=branch_first_jets(source,qsource['rows'],self.first_owner.dstar,phi,branch)
                except ArithmeticError as error:
                    local.append(dict(record=dict(status='requires_original_phase_range_refinement',
                        signed_u_branch=branch['name'],arithmetic_obstruction=str(error),zero_or_midpoint_fallback_used=False),values=None))
                    continue
                values=None
                if got['values'] is not None:
                    try:values=density.density_Z_kernels(E,E_Z,V,V_Z,got['values'],N)
                    except ArithmeticError as error:
                        got['record']=dict(status='requires_original_density_exponent_refinement',
                            signed_u_branch=branch['name'],arithmetic_obstruction=str(error),
                            original_phase_first_jet_source=got['record'],zero_or_midpoint_fallback_used=False)
                local.append(dict(record=got['record'],values=values))
            enclosed=all(cell['values'] is not None for cell in local)
            values=None
            if enclosed:
                union=density.local.same_source_union
                values={group:{key:union([cell['values'][group][key] for cell in local]) for key in density.RATES}
                    for group in ('kernels','Z_derivatives')}
            cells.append(dict(record=dict(status='enclosed' if enclosed else 'requires_original_q_or_phase_refinement',
                actual_global_radius_phase_cell=phi,original_conditional_branch_first_jets=[cell['record'] for cell in local],
                original_nonlinear_density_computed_before_branch_union=True,
                overlapping_branch_ranges_hulled_not_added=True),values=values,branches=local))
        enclosed=bool(cells) and all(cell['values'] is not None for cell in cells)
        record=dict(status='enclosed' if enclosed else 'requires_original_q_or_phase_refinement',source_family=self.family,
            chart=chart,candidate_N=N,original_q_slow_jet_source=qsource['record'],actual_original_radius_phase=geometry['record'],
            original_unconditioned_u=original_u.record(),conditional_branches=[dict(name=b['name'],condition=b['condition'],u=b['u'].record()) for b in branches],
            branches_proved_empty=empty,spatial_signed_density_branch_cells=[cell['record'] for cell in cells],
            cover_theorem='(-infinity,-1/8] union [-1/4,1/4] union [1/8,infinity)=R',
            all_nonempty_conditional_branch_ranges_required=True,original_nonzero_E_V_and_all_cross_terms_retained=True,
            branch_is_range_restriction_not_field_definition=True,original_root_q_jet_rows_unchanged=True,
            q_zero_only_flat_shortcut=True,no_branch_boundary_derivatives=True,
            actual_common_N_radius_phase_preserved=True,source_midpoint_or_zero_fallback_used=False,
            whole_chart_or_full_numerical_Rc_route_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,source=source,geometry=geometry,qsource=qsource)

    @native.inlet.source_precision
    def contribution(self,*,Z,left,right,N):
        lo=density.spatial.exact_coordinate(left);hi=density.spatial.exact_coordinate(right)
        if lo is None or hi is None or not 0<=lo<hi<=1:raise ValueError('Exact ordered O2_slope endpoints required')
        c=self.ctx;cv=lambda f:c.mpf(f.numerator)/f.denominator
        box=c.mpf((ep(cv(lo))[0],ep(cv(hi))[1]));query=self.spatial_query('O2_slope',Z,box,N)
        if any(cell['values'] is None for cell in query['cells']):raise ArithmeticError('Original covered cell unresolved; no integral fallback')
        width=cv(hi-lo);union=density.local.same_source_union
        ranges={group:{key:union([cell['values'][group][key] for cell in query['cells']]) for key in density.RATES}
            for group in ('kernels','Z_derivatives')}
        masses={key:density.local.positive_kernel_mass(c,width,rate) for key,rate in density.RATES.items()}
        contributions={key:ranges['kernels'][key]*masses[key] for key in density.RATES}
        derivatives={key:ranges['Z_derivatives'][key]*masses[key] for key in density.RATES}
        record=dict(original_whole_cell_signed_u_density_cover=query['record'],
            exact_O2_log_radius_width_fraction=density.spatial.fractional_record(hi-lo),
            positive_Duhamel_masses=masses,C0_contributions={k:v.record() for k,v in contributions.items()},
            Z_derivative_contributions={k:v.record() for k,v in derivatives.items()},
            phase_and_overlapping_u_unions_before_single_integral_mass=True,
            original_C1_derivative_under_fixed_Z_independent_integral=True,pressure_zero_rate_memory_not_reset=True,
            original_incoming_histories_not_replaced_by_local_contributions=True,
            full_numerical_Rc_targets_or_controls_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,contributions=contributions,Z_derivatives=derivatives,source=query['source'],geometry=query['geometry'])


@native.inlet.source_precision
def run(original_owner):
    began=time.monotonic();owner=NativeSignedUDensityCover(original_owner);c=owner.ctx
    box=c.mpf((ep(c.mpf('.12'))[0],ep(c.mpf('.15'))[1]))
    query=owner.spatial_query('O2_slope',(-1,1),box,2048)
    if query['record']['status']!='enclosed':raise ArithmeticError('Declared signed-crossing original O2 cover unresolved')
    integral=owner.contribution(Z=(-1,1),left='.12',right='.15',N=2048)
    result=dict(source_family=owner.family,**{GATE:True},actual_broad_original_O2_query=query['record'],
        actual_broad_original_O2_local_integral=integral['record'],candidate_N=2048,
        full_Z_interval=[-1,1],O2_coordinate_interval=['.12','.15'],
        original_broad_crossing_source_resolved=True,actual_local_C0_Z_integral_rows=10,
        scalar_values_or_actual_five_controls_or_terminal_closure_admitted=False,
        all_17_chart_original_range_oracle_admitted=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Original broad signed-crossing O2 source/radius-phase density C0/Z cover and ten local Duhamel contribution ranges. No selected field values, whole-route targets, controls, global closure or recursion admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result
