"""True original actual-patch terminal source points, including the Rh seam.

After the last compact support, the accepted unique implicit map supplies
exact background moment identities. The closed source is evaluated in the
original patch coordinate; no active-patch coefficient box defines a value.
This point adapter does not supply incoming correction histories or integrals.
"""
import ast
import copy
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_point_source_leaves as leaves
import lei_ren_part1_paper_compliant_actual_Rh_source_join as join

reference=leaves.reference; base=leaves.base; point=base.point; slow=leaves.slow
precise=leaves.precise; density=leaves.density; ep=leaves.ep
HERE,PREFIX,sha=leaves.HERE,leaves.PREFIX,leaves.sha
NAME=PREFIX+'current_original_patch_terminal_points.json'
RECEIPT=PREFIX+'current_original_patch_terminal_points_check.json'
GATE='actual_original_patch_terminal_true_source_C0_Z_point_callbacks_installed'
EDGE=s.Rational(71,40)


def log_coordinate(x):
    if x==s.E:return s.Integer(1)
    if x.func==s.exp:return x.args[0]
    return s.log(x)


def terminal_coordinate(value):
    if isinstance(value,dict):
        if set(value)!={'original_patch_log_offset'}:raise ValueError('Original patch coordinate recipe required')
        q=precise.rational(value['original_patch_log_offset'])
        if not 0<=q<=1:raise ValueError('Original patch log offset in[0,1] required')
        x=s.exp(q)
    elif isinstance(value,s.Basic) and (value==s.E or value.func==s.exp):
        q=log_coordinate(value)
        if not q.is_Rational or not 0<=q<=1:raise ValueError('Exact original patch exponential coordinate required')
        precise.rational(q);x=value
    else:x=precise.rational(value)
    if x.is_Rational:
        if x<EDGE:raise ValueError('Unresolved active patch excluded; x>=71/40 required')
    c=base.MPIntervalContext();c.dps=100
    with mp.workdps(140):
        v=numeric_coordinate(c,x)
        if x.is_Rational and ep(v)[1]>ep(c.exp(1))[0]:raise ValueError('Original patch x<=e required')
        if not x.is_Rational and ep(v)[0]<ep(c.mpf(71)/40)[1]:raise ValueError('Unresolved active patch excluded')
    return x


def numeric_coordinate(c,x):
    if x.is_Rational:return c.mpf(int(x.p))/int(x.q)
    q=log_coordinate(x);return c.exp(c.mpf(int(q.p))/int(q.q))


def patch_offset(c,x):
    if x.is_Rational:return -6+c.ln(numeric_coordinate(c,x))
    q=log_coordinate(x);return -6+c.mpf(int(q.p))/int(q.q)


def phase_coordinate(x):
    return str(x) if x.is_Rational else {'original_patch_log_offset':str(log_coordinate(x))}


def terminal_recipe_binding():
    # This theorem binds the actual terminal expressions and their physical
    # unit conversions, not just the values of saved chart ranges.
    theorem=join.Rh_functional_identity()
    if not theorem['passed'] or not theorem['current_unique_implicit_map_full_weight_closure_used']:
        raise ValueError('Original full-support unique implicit-map identity required')
    tree=ast.parse((HERE/(PREFIX+'actual_moment_patch.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantActualMomentPatch')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    terminal=next(n.value for n in fn.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='terminal' for t in n.targets))
    if ast.dump(terminal)!=ast.dump(ast.parse('endpoints(x)[0]>=mp.mpf(71)/40',mode='eval').body):
        raise ValueError('Original terminal support domain changed')
    branch=next(n for n in fn.body if isinstance(n,ast.If) and isinstance(n.test,ast.Name) and n.test.id=='terminal')
    reset=next(n.value for n in branch.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='defects' for t in n.targets))
    if ast.dump(reset)!=ast.dump(ast.parse('[zero]*5',mode='eval').body):raise ValueError('Original terminal functional identity changed')
    x,z,ps=s.symbols('x Z Pstar',positive=True)
    am=s.exp(-s.Rational(3,5))/(1+z*z);u=am*x**s.Rational(1,10);V=4*z
    # Convert original Rm units to the current-R units of the inertial compiler.
    patch=dict(m=V,h=am*s.Rational(5,8)*x**s.Rational(8,5)/x**s.Rational(3,2),
        k=V*am*s.Rational(5,8)*x**s.Rational(8,5)/x**s.Rational(3,2),
        e=V*V/ps**2-am**2*s.Rational(5,12)*x**s.Rational(6,5)/x,
        p=am**2*s.Rational(5,2)*x**s.Rational(1,5))
    closed=dict(m=V,h=s.Rational(5,8)*u,k=V*s.Rational(5,8)*u,
        e=V*V/ps**2-s.Rational(5,12)*u*u,p=s.Rational(5,2)*u*u)
    identities={name:s.simplify(patch[name]-closed[name])==0 for name in patch}
    if not all(identities.values()):raise ArithmeticError('Patch current-R normalization failed')
    return dict(passed=True,accepted_terminal_physical_function_theorem=theorem,
        exact_last_support_edge='71/40',terminal_branch_bound_to_original_source=True,
        current_R_five_history_normalization_identities=identities,
        true_offset_recipe='log(x)-6',exact_velocity_recipes={'E':'exp(-3/5)*x^(1/10)/(1+Z^2)','V':'4Z','a':'4/5','b':'0'},
        inertial_compiler_uses_current_R_normalized_histories=True,
        P0_separate_same_analytic_function_added_once=True,
        zero_terminal_background_defect_is_implicit_function_identity_not_incoming_correction_zero=True)


def project_method(module,clsname,method):
    """Change coordinate validation/physical radius only; reverse AST check."""
    tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==clsname)
    original=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
    fn=copy.deepcopy(original);changed=[]
    targets={
        'source.exact_coordinate(y)':'terminal_coordinate(y)',
        'point.source.exact_coordinate(y)':'terminal_coordinate(y)',
        'self.frame.radius_log(y)':"self.frame.definitions['logRref']-6+s.log(y)",
        'c.mpf(int(y.p)) / int(y.q)':'patch_offset(c,y)'}
    expected={'evaluate':2,'query':2}[method]
    maps={ast.dump(ast.parse(k,mode='eval').body):v for k,v in targets.items()}
    class Project(ast.NodeTransformer):
        def visit(self,node):
            key=ast.dump(node)
            if key in maps:
                new=ast.parse(maps[key],mode='eval').body;changed.append((key,ast.dump(new)))
                return new
            return super().visit(node)
    Project().visit(fn)
    if len(changed)!=expected:raise ValueError('Original point coordinate sites changed')
    restored=copy.deepcopy(fn)
    reverse={b:next(ast.parse(k,mode='eval').body for k in targets if ast.dump(ast.parse(k,mode='eval').body)==a) for a,b in changed}
    class Restore(ast.NodeTransformer):
        def visit(self,node):
            key=ast.dump(node)
            if key in reverse:return copy.deepcopy(reverse[key])
            return super().visit(node)
    Restore().visit(restored)
    if ast.dump(restored)!=ast.dump(original):raise ValueError('Original inertial/pressure/phase math changed')
    env=dict(vars(module),terminal_coordinate=terminal_coordinate,patch_offset=patch_offset)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original point method: actual patch terminal coordinate/radius>','exec'),env)
    return env[method]


class PatchTerminalFactoredInputs(point.OriginalO2FactoredPointInputs):
    mode='original_actual_patch_terminal_closed_point_coefficients_with_directed_errors'
    def radial(self,y):
        y=terminal_coordinate(y)
        if y not in self.radial_cache:
            c=self.ctx;iv=self.interval
            f=c.exp(patch_offset(c,y)/10);F=iv.exp(patch_offset(iv,y)/10)
            values=dict(f=f,H=c.mpf(5)/8*f,D=c.mpf(5)/12*f*f,P=c.mpf(5)/2*f*f,a=c.mpf(4)/5)
            boxes=dict(f=F,H=iv.mpf(5)/8*F,D=iv.mpf(5)/12*F**2,P=iv.mpf(5)/2*F**2,a=iv.mpf(4)/5)
            self.radial_cache[y]=(values,boxes)
        return self.radial_cache[y]


PatchTerminalFactoredInputs.evaluate=project_method(point,'OriginalO2FactoredPointInputs','evaluate')
terminal_query=project_method(base,'OriginalO2ConditionedPrimitives','query')


class OriginalPatchTerminalPointLeaves:
    mode='original_patch_terminal_point_partial'
    def __init__(self,provider=None):
        self.provider=provider if provider is not None else leaves.OriginalPointSourceLeaves()
        if type(self.provider) is not leaves.OriginalPointSourceLeaves:raise TypeError('Genuine accepted source provider required')
        p=self.provider
        self.family=self.source_family=p.family;self.source_graph_sha256=p.source_graph_sha256
        self.built=p.built;self.graph_digest=p.graph_digest;self.phase=p.phase;self.views=p.views
        self.rows=p.rows;self.hashes=dict(p.hashes);self.recipe=terminal_recipe_binding()
        self.owner=copy.copy(p.reference.owner);old=self.owner.inputs
        inputs=PatchTerminalFactoredInputs.__new__(PatchTerminalFactoredInputs)
        inputs.__dict__=dict(old.__dict__,radial_cache={},point_cache={});self.owner.inputs=inputs
        self._frames={};self._issued={};self.materialization_records=[]
        for stem,flags in (('actual_feedback_patch_mixed_C4_check',('current_Rh_same_unique_implicit_function_and_full_support_retained','current_five_functional_terminal_identities_connected')),
            ('actual_Rh_source_join_check',('current_Rh_external_neighbor_join_certified',))):
            name=PREFIX+stem+'.json';record=json.loads((HERE/name).read_bytes())
            if not record['all_passed'] or not all(record[key] for key in flags):raise ValueError('Accepted current terminal identity required')
            for key in ('implicit_source_sha256','datum_enclosure_sha256'):
                if record[key]!=self.family[key]:raise ValueError('Current terminal source family differs')
            p.bind_hashes({**record['input_hashes'],name:sha(name)})
        self.hashes=dict(p.hashes)
        for name in (Path(__file__).name,Path(join.__file__).name):self.hashes[name]=sha(name)

    def require_row(self,row):
        g=self.built['graph']
        if leaves.digest_rows(g.nodes)!=self.graph_digest or not any(row is item for item in g.nodes):raise ValueError('Unchanged issued original row required')
        if row.get('operation')!='original_function_graph' or row.get('graph_sha256')!=self.source_graph_sha256:raise ValueError('Same original source graph required')
        if row.get('chart')!='actual_patch':raise NotImplementedError('This adapter supports original actual_patch terminal only')
        role=row.get('function_role','');parts=role.split('_')
        if row.get('source_graph_namespace')!='function_graph_nodes' or len(parts)!=3 or parts[0]!='density':raise ValueError('Original density role required')
        key,order=parts[1:]
        if key not in leaves.transport.RATES or order not in ('C0','Z'):raise ValueError('Original C0/Z density role required')
        view=self.views['actual_patch'];expected=view['five_signed_increment_rate_roots'][key] if order=='C0' else view['five_signed_increment_rate_first_derivatives']['Z'][key]
        if row['source_node']!=expected or row['shared_N']!=self.built['N'].node:raise ValueError('Original root/N role differs')
        return key,order

    def frame(self,*,chart='actual_patch',coordinate,Z,N,bits=80):
        if chart!='actual_patch':raise NotImplementedError('Only original patch terminal installed')
        x=terminal_coordinate(coordinate);z=point.pressure.exact_Z(Z)
        if type(N) is not int or N<160 or N.bit_length()>4096:raise ValueError('Candidate integer N>=160 required')
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        cache=(x,z,N,bits)
        if cache in self._frames:return self._frames[cache]
        row=self.rows['actual_patch','density_m_C0'];self.require_row(row)
        geometry=self.phase.phase_for_source(row,self.built,coordinate=phase_coordinate(x),N=N)
        c=self.owner.ctx;pieces=[]
        with mp.workdps(c.dps+40):
            query=terminal_query(self.owner,y=x,Z=z);kernel=query['kernel']
            query['basis_contract'].update(source_R='Rref*exp(-6)*x',exact_a1_equals2_intersection_used=False,
                coordinate_kind='original actual_patch x=R/Rm',patch_current_R_normalization_bound=True)
            a=base.conditioned.bounded_value(query['roots']['a'][0,0])
            if not 0<ep(a)[0]<=ep(a)[1]<=2:raise ArithmeticError('Actual terminal modulation bound failed')
            for box in geometry['actual_phase_directed_boxes']:
                inverse=kernel.evaluate(c.mpf(box),bits=bits)
                if inverse['status']!='enclosed':raise ArithmeticError('Original terminal source inverse needs refinement')
                selected=inverse['selected_inverse'];primitives=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                jets,derivative=slow.slow_values(kernel,query['roots'],selected['coordinate_interval'],selected['chart'])
                interpreter=density.BoundDensityGraph(self.views['actual_patch'],query,primitives,jets,N);values=interpreter.outputs()
                contract=interpreter.record();contract.update(point_chart='actual_patch',actual_point_a_enclosure=a,
                    bound_proof='original terminal a=4/5; same A=a*(phi-psi/(2*pi))/2; N>=160',
                    integration_coordinate='log(R/Rref)=log(x)-6; density per log radius, not per x',
                    B_over_Pstar_already_normalized_not_multiplied_again=True)
                inverse.update(free_phase_parameter_not_spatial_phase=False,original_common_N_and_radius_phase_bound=True)
                pieces.append(dict(values=values,primitives=primitives,jets=jets,inverse=inverse,derivative=derivative,interpreter=interpreter,execution_contract=contract))
            union=lambda rows:leaves.unions.same_source_union(rows)
            values={order:{key:union([piece['values'][label][key] for piece in pieces]) for key in leaves.transport.RATES}
                for order,label in (('C0','densities'),('Z','density_Z'))}
            if any(v.ctx is not c or v.ledger is not query['ledger'] or v.scale.bases is not kernel.q.scale.bases for row in values.values() for v in row.values()):raise ValueError('One exact source point basis/context/ledger required')
            record=dict(source_family=self.family,chart=chart,original_coordinate_exact=str(x),original_Z_exact=str(z),explicit_candidate_N=N,
                source_graph_sha256=self.source_graph_sha256,actual_original_point_phase=geometry,
                actual_defining_point_inputs_and_errors=point.record_query(query['point']),source_factor_basis=query['basis_contract'],
                native_point_inverse_and_derivative_records=[dict(inverse=p['inverse'],derivative=p['derivative'],graph=p['execution_contract']) for p in pieces],
                actual_original_five_C0_Z_point_source_enclosures={order:{key:v.record() for key,v in row.items()} for order,row in values.items()},
                numerical_arithmetic_ledger=query['ledger'],original_P0_P0_Z_not_modified=True,
                ordinary_Z_at_fixed_true_radius_phase=True,actual_graph_nodes_executed_not_cover_values_selected=True,
                original_pressure_coefficient_and_late_flatten_errors_retained=True,
                terminal_background_identity_not_unknown_incoming_correction_zero=True,
                point_enclosure_not_a_continuous_integral=True)
            frame=leaves.OriginalPointSourceFrame(chart,x,z,N,c,values,record,query,tuple(pieces))
        self._frames[cache]=frame;self._issued[id(frame)]=frame;return frame

    def dispatch(self,row,frame):
        key,order=self.require_row(row)
        if type(frame) is not leaves.OriginalPointSourceFrame or self._issued.get(id(frame)) is not frame or frame.chart!='actual_patch':raise ValueError('Issued matching terminal frame required')
        self.phase.phase_for_source(row,self.built,coordinate=phase_coordinate(frame.coordinate),N=frame.N)
        return frame.values[order][key]

    def source_factored(self,row,*,coordinate,Z,N,bits=80):
        self.require_row(row);return self.dispatch(row,self.frame(coordinate=coordinate,Z=Z,N=N,bits=bits))

    source=leaves.OriginalPointSourceLeaves.source
    parameter=leaves.OriginalPointSourceLeaves.parameter
    integrate=leaves.OriginalPointSourceLeaves.integrate


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();owner=OriginalPatchTerminalPointLeaves();frames=[];ordinary=0;unmaterialized=[]
    samples=[('71/40','.37',160),('1.831','-.37',257),('2.337','0',257),('2.337','.37',257),
        ({'original_patch_log_offset':'1'},'.37',257)]
    for x,z,N in samples:
        frame=owner.frame(coordinate=x,Z=z,N=N);frames.append(frame)
        for key in leaves.transport.RATES:
            for order in ('C0','Z'):
                row=owner.rows['actual_patch','density_'+key+'_'+order]
                assert owner.dispatch(row,frame) is frame.values[order][key]
                try:owner.source(row,coordinate=phase_coordinate(frame.coordinate),Z=z,N=N);ordinary+=1
                except ArithmeticError:unmaterialized.append(dict(coordinate=str(frame.coordinate),Z=z,N=N,role=row['function_role'],source_retained_factored=True))
        print('Original actual-patch terminal source:',str(frame.coordinate),z,N,flush=True)
    result=dict(**{GATE:True},mode=owner.mode,source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        exact_original_transport_graph_sha256=owner.graph_digest,terminal_source_recipe=owner.recipe,
        actual_original_point_source_frames=[f.record for f in frames],original_role_bound_point_dispatches=10*len(frames),
        ordinary_interval_callbacks=ordinary,actual_ordinary_interval_point_callbacks=owner.materialization_records,
        unmaterializable_source_values_retained_factored=unmaterialized,
        actual_patch_active_support_source_installed=False,actual_patch_terminal_integral_installed=False,
        full_17_chart_source_or_integral_oracle_installed=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(precise.phase.packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Genuine actual_patch terminal C0/Z point source for71/40<=x<=e, same original graph/radius phase/pressure errors/current-R inertial compiler. Active patch, full incoming correction history, integrals, controls, global N and recursion remain open.')
    (HERE/NAME).write_bytes(json.dumps(precise.encode(result),indent=2).encode()+b'\n')
    return (result,owner,frames) if return_live else result


if __name__=='__main__':run()
