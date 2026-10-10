"""Actual physical-coordinate inverse on the accepted original Rp graph.

Reuse the original directed implicit solver, not its old registry/owners.
Both exact Cartesian inputs and logarithmic cylindrical inputs produce true
inverse function nodes plus directed coordinate enclosures. A log-z adapter
copies the same root iteration without ever materializing exp(log_abs_z).
Chart admission and sufficiently accurate velocity point evaluation remain
separate tasks; no uncertain root or source coefficient is selected.
"""
import ast
import copy
from fractions import Fraction
import functools
import gzip
import hashlib
import inspect
import json
import math
from pathlib import Path
import time
import textwrap
import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_Rp_physical_source_map as physical
import lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator as original
from lei_ren_part1_paper_compliant_current_global_tensor_cover_operator import independent_implicit_cover_theorem
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=physical.HERE,physical.PREFIX,physical.sha
NAME=PREFIX+'current_original_Rp_physical_inverse.json.gz'
RECEIPT=PREFIX+'current_original_Rp_physical_inverse_check.json'
GATES=('current_original_Rp_exact_physical_inverse_function_nodes_installed',
       'current_original_Rp_directed_Cartesian_coordinate_inverse_installed',
       'current_original_Rp_directed_log_cylindrical_inverse_installed')
OPEN=tuple(dict.fromkeys(physical.OPEN+('current_original_Rp_native_chart_location_installed',
    'current_original_Rp_inverse_intervals_consumed_by_point_source',)))
ends=original.endpoints


def exact_scalar(value):
    """A binary float means that exact binary number, not a guessed decimal."""
    if type(value) is float:
        if not math.isfinite(value):raise ValueError('Finite physical scalar required')
        return Fraction.from_float(value)
    if type(value) in (int,str) or type(value) is Fraction:
        try:return Fraction(value)
        except (ValueError,OverflowError):raise ValueError('Finite exact physical scalar required') from None
    raise TypeError('Exact rational/string or finite binary float required; intervals are not input values')


def box(c,q):return c.mpf(q.numerator)/q.denominator


def current_forward_unit_binding():
    tree=ast.parse(textwrap.dedent(inspect.getsource(inspect.unwrap(physical.CurrentOriginalRpPhysicalSourceMap.coordinates))))
    statements={
        'logr':"g.add(loglambda, g.mul(half, g.add(g.unary('log', g.constant(2)), lr)))",
        'axial':"g.mul(z, g.unary('exp', g.mul(g.sub(g.one, delta), loglambda)))"}
    bound={}
    for name,expression in statements.items():
        wanted=ast.dump(ast.parse(expression,mode='eval').body,include_attributes=False)
        found=[node for node in ast.walk(tree) if isinstance(node,ast.Assign)
            and any(isinstance(target,ast.Name) and target.id==name for target in node.targets)]
        if len(found)!=1 or ast.dump(found[0].value,include_attributes=False)!=wanted:
            raise ValueError('Current original unit-viscosity coordinate source changed: '+name)
        bound[name]=expression
    return dict(passed=True,current_original_nu1_forward_assignments=bound,
        current_forward_coordinates_AST_sha256=hashlib.sha256(ast.dump(tree,include_attributes=False).encode()).hexdigest())


@functools.lru_cache(maxsize=1)
def log_axial_program():
    """Copy the admitted root loop; change only how log(abs(z)) enters it."""
    source=ast.parse(inspect.getsource(inspect.unwrap(original.implicit_log_coordinate_map)))
    function=source.body[0]
    branches=[node for node in function.body if isinstance(node,ast.If) and ast.unparse(node.test)=='zl == zh == 0']
    if len(branches)!=1:raise ValueError('Original zero/nonzero inverse split changed')
    nonzero=branches[0].orelse
    start=next(i for i,node in enumerate(nonzero) if isinstance(node,ast.Assign)
        and any(isinstance(target,ast.Name) and target.id=='axial_log' for target in node.targets))
    if start!=2 or ast.unparse(nonzero[start])!='axial_log = 2 * logabs - ctx.ln(nu)':
        raise ValueError('Original axial logarithm source binding changed')
    body=copy.deepcopy(nonzero[start:])
    required={
        'lo':'endpoints(lt / 2)[0]',
        'hi':'max(endpoints((lt + ctx.ln(2)) / 2)[1], endpoints((axial_log + ctx.ln(2)) / (2 * (1 - d)))[1])',
        'derivative':'ctx.mpf([endpoints(2 * (1 - d))[0], 2])',
        'f':'2 * ctx.mpf(mid) - logaddexp(ctx, lt, axial_log + 2 * d * ctx.mpf(mid))',
        'contracted':'intersect(ctx, ctx.mpf([lo, hi]), ctx.mpf(mid) - f / derivative)'}
    for name,wanted in required.items():
        assignments=[node for node in ast.walk(ast.Module(body=body,type_ignores=[]))
            if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in node.targets)]
        if sum(ast.unparse(node.value)==wanted for node in assignments)!=1:
            raise ValueError('Original directed log root statement changed: '+name)
    prefix=ast.parse('''lt=ctx.mpf(log_tau)
logabs=ctx.mpf(log_abs_z)
d=ctx.mpf(delta)
nu=ctx.mpf(1)
zl=zh=sign_z
tolerance=mp.mpf(relative_tolerance)
''').body
    suffix=ast.parse('''return dict(actual_log_lambda=q,Z=Z,exact_log_one_minus_Z_squared_source=log_complement,
    log_root_enclosure_width=ctx.mpf(endpoints(q)[1]-endpoints(q)[0]),
    solver_status=status,iterations=iterations,relative_log_root_tolerance=str(relative_tolerance),
    physical_viscosity=nu,source_delta=d,requested_log_tau=lt,
    requested_log_abs_z=logabs,physical_z_sign=sign_z,
    true_finite_point_has_abs_Z_strictly_below_one=True,
    Z_enclosure_can_touch_infinity_limit_due_to_precision=True,
    midpoint_used_only_as_root_iteration_trial=True,actual_log_tau_retained=True,
    exp_truncation_is_enclosure_only_not_source=True)
''').body
    generated=ast.FunctionDef(name='_original_log_axial_inverse',args=ast.arguments(
        posonlyargs=[],args=[ast.arg(arg=name) for name in
            ('ctx','log_abs_z','log_tau','delta','sign_z','relative_tolerance','max_steps')],
        kwonlyargs=[],kw_defaults=[],defaults=[]),body=prefix+body+suffix,decorator_list=[])
    env=dict(inspect.unwrap(original.implicit_log_coordinate_map).__globals__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[generated],type_ignores=[])),
        '<original-directed-log-root-with-explicit-log-z>','exec'),env)
    return env[generated.name],dict(passed=True,
        original_directed_implicit_solver_AST_sha256=hashlib.sha256(ast.dump(function,include_attributes=False).encode()).hexdigest(),
        copied_original_nonzero_root_and_Z_program_AST_sha256=hashlib.sha256(ast.dump(ast.Module(body=body,type_ignores=[]),include_attributes=False).encode()).hexdigest(),
        exact_adapted_program_AST_sha256=hashlib.sha256(ast.dump(generated,include_attributes=False).encode()).hexdigest(),
        original_root_loop_derivative_bracket_and_exp_helpers_unchanged=True,
        only_explicit_log_abs_z_input_and_exact_nu1_prefix_added=True,
        no_physical_z_or_lambda_exponential_materialized=True,
        original_defining_statements=required)


@source_precision
def log_coordinate_map(c,log_abs_z,log_tau,delta,sign_z,relative_tolerance='1e-60',max_steps=512):
    if type(sign_z) is not int or sign_z not in (-1,0,1):raise ValueError('Explicit exact physical z sign required')
    lt=c.mpf(log_tau);d=c.mpf(delta)
    if not all(mp.isfinite(v) for v in ends(lt)+ends(d)) or not 0<=ends(d)[0]<=ends(d)[1]<1:
        raise ValueError('Finite log time and admitted delta interval required')
    tolerance=mp.mpf(relative_tolerance)
    if tolerance<=0 or not mp.isfinite(tolerance) or type(max_steps) is not int or max_steps<1:
        raise ValueError('Positive finite tolerance and integer step count required')
    if sign_z==0:
        if log_abs_z is not None:raise ValueError('Zero axial coordinate has no finite log(abs(z))')
        return original.implicit_log_coordinate_map(c,0,lt,1,d,relative_tolerance,max_steps)
    if log_abs_z is None:raise ValueError('Nonzero axial coordinate requires finite log(abs(z))')
    lz=c.mpf(log_abs_z)
    if not all(mp.isfinite(v) for v in ends(lz)):raise ValueError('Finite axial logarithm required')
    fn,_=log_axial_program()
    return fn(c,lz,lt,d,sign_z,relative_tolerance,max_steps)


class CurrentOriginalRpPhysicalInverse:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else physical.CurrentOriginalRpPhysicalSourceMap()
        if type(self.before) is not physical.CurrentOriginalRpPhysicalSourceMap or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current physical source map required')
        self.graph=self.before.graph;self.ctx=self.before.ctx;self.radius=self.before.radius
        self.delta=self.before.delta;self.delta_function=self.before.delta_function
        self.family_record=self.before.family_record;self.original_N_definition=self.before.original_N_definition
        self.theorem=original.implicit_physical_source_theorem();self.independent_theorem=independent_implicit_cover_theorem()
        _,self.program_binding=log_axial_program()
        self.forward_unit_binding=current_forward_unit_binding()
        self.hashes=dict(self.before.hashes)
        physical.mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,self.theorem['input_hashes'])
        for name in (physical.NAME,physical.RECEIPT,Path(__file__).name,
                PREFIX+'current_physical_tensor_locator_operator.py',PREFIX+'current_global_tensor_cover_operator.py'):
            self.hashes[name]=sha(name)
        # The old receipt supplies only the independently checked pure root
        # program, never its old 33-region routing or archived scalar values.
        name=PREFIX+'current_physical_tensor_locator_check.json'
        old=json.loads((HERE/name).read_bytes())
        if not old['all_passed'] or not old['current_actual_implicit_physical_coordinates_directed']:
            raise ValueError('Accepted original directed implicit root program required')
        if (old['actual_five_defect_family_sha256'],old['implicit_source_sha256'])!=(self.before.before.post.family,self.before.before.post.source):
            raise ValueError('Original pure inverse theorem family differs')
        self.original_root_fixture_evidence=old['independent_implicit_physical_root_fixtures']
        physical.mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,old['input_hashes'])
        self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self.call_trace=[];self._solve_witnesses={};self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current original inverse receipt/scope differs')
            if receipt['source_family']!=self.family_record or receipt['original_log_axial_program_binding']!=self.program_binding:
                raise ValueError('Current inverse family or directed source program differs')
            if receipt['current_forward_unit_viscosity_binding']!=self.forward_unit_binding:
                raise ValueError('Current original forward unit binding differs')
            physical.mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        result=dict(accepted_same_current_physical_map=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            same_original_inverse_graph_context_and_delta=self.graph is self.before.graph is self.radius.graph
                and self.ctx is self.before.ctx and self.delta_function==self.radius.functions['delta']
                and self.delta is self.before.delta,
            exact_unit_viscosity_matches_current_forward_map=self.forward_unit_binding==current_forward_unit_binding()
                and self.forward_unit_binding['passed'],
            original_directed_solver_theorem=self.theorem['passed'] and self.independent_theorem['passed']
                and self.independent_theorem['unique_global_axial_inverse'],
            original_log_input_program_recomputed=self.program_binding==log_axial_program()[1],
            original_finite_N_and_analytic_pressure_history_retained=self.original_N_definition==self.before.original_N_definition)
        if not all(result.values()):raise ValueError('Current original inverse source graph differs: '+str(result))
        return result

    def inverse_functions(self,logr,logz,lt,theta,sign_z):
        self.assert_graph();g=self.graph
        for ref in (logr,lt,theta)+(() if logz is None else (logz,)):
            if type(ref) is not physical.mixed.pulse.radius.FunctionRef or ref.graph is not g:
                raise TypeError('Same current graph exact physical input functions required')
        if type(sign_z) is not int or sign_z not in (-1,0,1) or (sign_z==0)!=(logz is None):
            raise ValueError('Consistent exact axial sign/log input required')
        q=(g.mul(g.constant('1/2'),lt) if sign_z==0 else g.node(
            'exact_original_physical_log_lambda_inverse',log_tau=lt.node,log_abs_z=logz.node,
            physical_z_sign=sign_z,delta=self.delta_function.node,viscosity=g.one.node,
            definition='unique finite q: 2q=logaddexp(log_tau,2log_abs_z+2delta*q)',
            positive_root_certificate='tau=exp(log_tau)>0; 0<delta<1; F_prime>=2(1-delta)>0',
            original_source_module=PREFIX+'current_physical_tensor_locator_operator.py',
            original_source_function='implicit_log_coordinate_map'))
        lr=g.add(g.mul(g.constant(2),logr),g.neg(g.unary('log',g.constant(2))),g.mul(g.constant(-2),q))
        lz=None if sign_z==0 else g.sub(logz,g.mul(g.sub(g.one,self.delta_function),q))
        Z=g.zero if lz is None else g.mul(g.constant(sign_z),g.unary('exp',lz))
        complement=g.sub(lt,g.mul(g.constant(2),q))
        return dict(log_r=logr,log_tau=lt,theta=theta,log_lambda=q,logR=lr,R=g.unary('exp',lr),
            lambda_value=g.unary('exp',q),Z=Z,log_abs_Z=lz,
            log_one_minus_Z_squared=complement,delta=self.delta_function,
            tau=g.unary('exp',lt),t=g.sub(g.one,g.unary('exp',lt)))

    def _solve_fingerprint(self,kind,exact_inputs,refs,mapping,lr,theta_box):
        if any(ref is not None and ref.graph is not self.graph for ref in refs.values()):
            raise ValueError('Current inverse result functions required')
        data=dict(kind=kind,exact_inputs=exact_inputs,
            refs={key:None if ref is None else ref.node for key,ref in refs.items()},
            mapping=mapping,log_radius=lr,theta=theta_box)
        return hashlib.sha256(json.dumps(physical.mixed.pulse.raw.packed(data),sort_keys=True,separators=(',',':')).encode()).hexdigest()

    def _record_solve(self,kind,exact_inputs,refs,mapping,lr,theta_box):
        witness=object()
        self._solve_witnesses[id(witness)]=(witness,mapping,refs,
            self._solve_fingerprint(kind,exact_inputs,refs,mapping,lr,theta_box))
        return witness

    @source_precision
    def assemble(self,kind,exact_inputs,refs,mapping,lr,theta_box,*,witness=None):
        self.assert_graph();c=self.ctx;q=mapping['actual_log_lambda']
        proof=self._solve_witnesses.get(id(witness))
        if proof is None or proof[0] is not witness or proof[1] is not mapping or proof[2] is not refs \
                or proof[3]!=self._solve_fingerprint(kind,exact_inputs,refs,mapping,lr,theta_box):
            raise ValueError('Unchanged live solve from this current inverse owner required')
        logR=2*lr-c.ln(2)-2*q
        lo,hi=ends(mapping['Z']);qlo,qhi=ends(q)
        result=dict(input_kind=kind,exact_input_record=exact_inputs,source_family=self.family_record,
            coordinate_functions=refs,directed_inverse_mapping=mapping,
            physical_log_radius_enclosure=lr,theta_enclosure=theta_box,source_logR_enclosure=logR,
            inverse_log_lambda_absolute_width=c.mpf(qhi)-c.mpf(qlo),
            strict_numerical_Z_interior_resolved=(-1<lo<=hi<1),
            unique_exact_finite_point_has_abs_Z_strictly_below_one=True,
            numerical_absolute_radius_or_lambda_materialized=False,
            numerical_bounds_are_not_defining_coordinate_functions=True,
            arbitrary_root_or_source_coefficient_not_selected=True,
            original_N_definition=self.original_N_definition,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self.call_trace.append(dict(kind=kind,actual_current_delta_enclosure_consumed=True,
            original_exact_delta_function=refs['delta'].node,
            source_intervals_not_midpoint_values=True,
            solver_status=mapping['solver_status'],root_iterations=mapping['iterations']))
        return result

    @source_precision
    def cartesian(self,x,y,z,t,relative_tolerance='1e-60',max_steps=512):
        self.assert_graph();g=self.graph;c=self.ctx
        xx,yy,zz,tt=map(exact_scalar,(x,y,z,t));tau=1-tt
        if tau<=0:raise ValueError('Physical t<1 and strictly positive tau required')
        if xx==yy==0:raise ValueError('Current Rp inverse is off axis; actual core limit remains open')
        xref,yref,zref,tref=(g.constant(v) for v in (xx,yy,zz,tt))
        norm2=g.add(g.mul(xref,xref),g.mul(yref,yref))
        lrref=g.mul(g.constant('1/2'),g.unary('log',norm2))
        ltref=g.unary('log',g.constant(tau));sign=1 if zz>0 else -1 if zz<0 else 0
        lzref=None if not sign else g.unary('log',g.constant(abs(zz)))
        angle=g.node('exact_physical_atan2',x=xref.node,y=yref.node,
            branch='(-pi,pi]; negative x with exact y=0 uses +pi',positive_radius_certificate='x^2+y^2>0')
        xb,yb,zb=map(lambda value:box(c,value),(xx,yy,zz))
        lr=c.ln(xb**2+yb**2)/2;lt=c.ln(box(c,tau));theta=c.atan2(yb,xb)
        mapping=original.implicit_log_coordinate_map(c,zb,lt,1,self.delta,relative_tolerance,max_steps)
        refs=self.inverse_functions(lrref,lzref,ltref,angle,sign)
        refs.update(x=xref,y=yref,z=zref,requested_t=tref)
        inputs=dict(x=str(xx),y=str(yy),z=str(zz),t=str(tt))
        witness=self._record_solve('exact_Cartesian',inputs,refs,mapping,lr,theta)
        return self.assemble('exact_Cartesian',inputs,refs,mapping,lr,theta,witness=witness)

    @source_precision
    def log_cylindrical(self,log_r,log_abs_z,sign_z,log_tau,theta='0',relative_tolerance='1e-60',max_steps=512):
        self.assert_graph();c=self.ctx;g=self.graph
        lr,lt,angle=map(exact_scalar,(log_r,log_tau,theta))
        lz=None if log_abs_z is None else exact_scalar(log_abs_z)
        mapping=log_coordinate_map(c,None if lz is None else box(c,lz),box(c,lt),self.delta,
            sign_z,relative_tolerance,max_steps)
        refs=self.inverse_functions(g.constant(lr),None if lz is None else g.constant(lz),g.constant(lt),g.constant(angle),sign_z)
        inputs=dict(log_r=str(lr),log_abs_z=None if lz is None else str(lz),sign_z=sign_z,log_tau=str(lt),theta=str(angle))
        lrb,ab=box(c,lr),box(c,angle)
        witness=self._record_solve('exact_log_cylindrical',inputs,refs,mapping,lrb,ab)
        return self.assemble('exact_log_cylindrical',inputs,refs,mapping,lrb,ab,witness=witness)


def report(view):
    return {**{key:value for key,value in view.items() if key!='coordinate_functions'},
        'exact_coordinate_function_nodes':{key:None if value is None else value.node
            for key,value in view['coordinate_functions'].items()}}


CASES={'ordinary':('0.7','-0.2','0.3','0.9'),
    'negative_axial':('-0.7','0.2','-0.3','0.9'),
    'zero_axial':('0.7','0.2','0','0.9'),
    'positive_y_axis':('0','1','0.013','0.999'),
    'negative_x_angle_cut':('-1','0','0.013','0.999')}
LOG_CASES={'near_time_zero_z':('0',None,0,'-1e1000','7/10'),
    'near_time_nonzero_z':('0','-2',1,'-1e1000','7/10'),
    'tiny_axial_log':('0','-1e1000',-1,'-10','7/10'),
    'large_axial_log':('0','1e1000',1,'-10','7/10')}


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpPhysicalInverse(before,require_checked=False)
    views={name:owner.cartesian(*args) for name,args in CASES.items()}
    views.update({name:owner.log_cylindrical(*args) for name,args in LOG_CASES.items()})
    result=dict(source_family=owner.family_record,candidate_current_physical_inverse_constructed=True,
        actual_source_graph=owner.assert_graph(),original_log_axial_program_binding=owner.program_binding,
        current_forward_unit_viscosity_binding=owner.forward_unit_binding,
        original_implicit_inverse_theorem=owner.theorem,independent_unique_axial_inverse_theorem=owner.independent_theorem,
        accepted_original_root_fixture_evidence=owner.original_root_fixture_evidence,
        actual_current_inverse_views={name:report(view) for name,view in views.items()},
        actual_source_call_trace=owner.call_trace,original_N_definition=owner.original_N_definition,
        exact_current_inverse_expression_graph=owner.graph.nodes,
        scope='Actual original graph delta; exact Cartesian and log cylindrical input functions; directed implicit lambda/Z/logR/angle enclosures. Off-axis, tau>0, nu=1. No old registry/locator owner. Native chart location, source point accuracy, global/axis field, stress, energy and recursion remain open',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(physical.mixed.pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('CURRENT_ORIGINAL_RP_PHYSICAL_INVERSE Cartesian and log-scale queries constructed',flush=True)
    return owner,views


if __name__=='__main__':run()
