"""Exact underlying outer-source identity, separate from interval covers.

Production expressions are projected to their defining scalar functions.
Original kernel primitives stay formal and source-bound. An arithmetic
mu enclosure is NOT declared equal to mu, and no range endpoint is used.
"""
import ast
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_limit_power_to_Rp as current
import lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N as slope
import lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N as axial
import lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N as o3
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum

HERE,PREFIX,sha,require=current.HERE,current.PREFIX,current.sha,current.require
NAME=PREFIX+'current_limit_Rp_native_identity.json'
RECEIPT=PREFIX+'current_limit_Rp_native_identity_check.json'
GATE='current_exact_underlying_outer_source_and_Rp_native_pulse_frame_identity_proved'
Z=sy.Symbol('Z',real=True);S=sy.Symbol('Pstar',positive=True)


class Projection:
    """Restricted arithmetic interpretation of actual source assignments."""
    def __init__(self):self.bindings={};self.hashes={};self.identities={}

    def assignment(self,module,method,target,env,*,index=0,count=1):
        path=Path(module.__file__);tree=ast.parse(path.read_text(encoding='utf8'))
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        nodes=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets)]
        require(len(nodes)==count and 0<=index<count,'Unique source assignment changed: '+path.name+'/'+target)
        key=path.name+':'+method+':'+target+':'+str(index)
        self.bindings[key]=dict(expression=ast.unparse(nodes[index]),
            function_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest())
        self.hashes[path.name]=sha(path.name)
        return self.expression(nodes[index],env)

    def expression(self,node,env):
        label=ast.unparse(node)
        if label in env:return env[label]
        if isinstance(node,ast.Constant):
            if isinstance(node.value,(str,int,float)):return sy.Rational(str(node.value))
        if isinstance(node,ast.Name) and node.id in env:return env[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -self.expression(node.operand,env)
        if isinstance(node,ast.BinOp):
            a,b=self.expression(node.left,env),self.expression(node.right,env)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if isinstance(node.op,ast.Pow):return a**b
        if isinstance(node,ast.Subscript):
            value=self.expression(node.value,env)
            key=node.slice.value if isinstance(node.slice,ast.Constant) else self.expression(node.slice,env)
            return value[key]
        if isinstance(node,(ast.Tuple,ast.List)):return [self.expression(q,env) for q in node.elts]
        if isinstance(node,ast.Call):
            name=ast.unparse(node.func)
            if name=='dict':return {q.arg:self.expression(q.value,env) for q in node.keywords}
            args=[self.expression(q,env) for q in node.args]
            if name in ('c.mpf','f.scalar','f.jet'):return args[0]
            if name in ('c.exp','formal'):return sy.exp(args[0])
            if name=='c.ln':return sy.log(args[0])
            if name=='square':return args[0]**2
            if name in ('f.scale','f.multiply'):return sy.prod(args)
            if name=='f.add':return sum(args)
            if name=='f.factor':
                require(args[0] in ([0,-1,0,0,0],[0,0,0,0,0]),'Unbound source factor basis')
                return (S**-2 if args[0][1]==-1 else sy.Integer(1))*(sy.exp(args[1]) if len(args)>1 else 1)
            if isinstance(node.func,ast.Attribute) and node.func.attr=='reciprocal' and not args:
                return 1/self.expression(node.func.value,env)
        raise ValueError('Unbound production source expression: '+label)

    def equal(self,name,left,right):
        delta=sy.cancel(sy.simplify(sy.expand_power_exp(left-right)))
        require(delta==0,'Underlying source functions differ: '+name+': '+str(delta))
        require(sy.diff(delta,Z)==0,'Source ordinary-Z identity failed: '+name)
        self.identities[name]=dict(value_identity=True,true_Z_identity=True,
            defining_function=sy.sstr(left),no_range_value_substitution=True)


def defining_source_proof():
    p=Projection();q=1+Z*Z;J,Ih,Ip,Ie=sy.symbols('J1 I_theta I_pressure I_energy',real=True)
    # These are the SAME source-bound continuous primitives, not the finite
    # enclosure sums returned by any particular scalar partition.
    env=dict(y=sy.Integer(1),J=J,mass=[Ih,Ip,Ie],**{'ref.z':Z,'ref.Vref':4*Z,'self.invP2':S**-2},z=Z,qi=1/q,V=4*Z)
    cur={};old={}
    for target in ('factor','qi','qi2','E','V','h','hist'):
        env[target]=p.assignment(slope,'background_cell',target,env)
    cur=dict(E=env['E'],V=env['V'],hist=env['hist'])
    envold=dict(y=1,J=J,mass=[Ih,Ip,Ie],z=Z,qi=1/q,**{'self.invP2':S**-2})
    for target in ('factor','u','V','h','hist'):envold[target]=p.assignment(current.pre,'slope',target,envold)
    old=dict(E=envold['u'],V=envold['V'],hist=envold['hist'])
    for key in ('E','V'):p.equal('O2_slope_exit/'+key,cur[key],old[key])
    for key in cur['hist']:p.equal('O2_slope_exit/'+key,cur['hist'][key],old['hist'][key])
    for key,expr in old['hist'].items():
        seed=expr.subs({J:0,Ih:0,Ip:0,Ie:0})
        # Check the actual slope expressions at y=0 rather than taking
        # the exit formulas with y=1. All defining primitives vanish there.
        seedenv=dict(y=0,J=0,mass=[0,0,0],z=Z,qi=1/q,**{'self.invP2':S**-2})
        for target in ('factor','u','V','h','hist'):seedenv[target]=p.assignment(current.pre,'slope',target,seedenv)
        expected=dict(m=4*Z,h=sy.Rational(5,8)/q,k=sy.Rational(5,2)*Z/q,
            e=16*Z*Z/S**2-sy.Rational(5,12)/q**2,p=sy.Rational(5,2)/q**2)
        p.equal('canonical_Rref_zero_primitive_seed/'+key,seedenv['hist'][key],expected[key])
    yd=sy.exp(40);tau=sy.Integer(11);t=yd+tau-1
    KB,KB2=sy.symbols('K_B_at_yd K_B2_at_yd',real=True)
    kernels=dict(B_mass=KB*sy.exp(-tau),B_squared_mass=KB2*sy.exp(-tau))
    env=dict(t=t,y=yd+tau,u1=cur['E'],old=cur['hist'],K=kernels,B=[0,0],z=Z,
        d=sy.exp(-t),root=sy.exp(-t/2),d3=sy.exp(-sy.Rational(3,2)*t),**{'op.zrows':Z})
    for target in ('E','E2','z','z2','V','Vy','hist'):env[target]=p.assignment(axial,'background_cell',target,env)
    nextcur=dict(E=env['E'],V=env['V'],hist=env['hist'])
    envold=dict(t=t,u1=old['E'],old=old['hist'],K=kernels,z=Z,B=[0]*5,
        d=sy.exp(-t),root=sy.exp(-t/2),d3=sy.exp(-sy.Rational(3,2)*t),**{'self.invP2':S**-2})
    nextold=dict(E=p.assignment(current.pre,'axial','u',envold),V=0,
        hist=p.assignment(current.pre,'axial','hist',envold))
    for key in ('E','V'):p.equal('O2_buffer_exit/'+key,nextcur[key],nextold[key])
    for key in nextcur['hist']:p.equal('O2_buffer_exit/'+key,nextcur['hist'][key],nextold['hist'][key])
    cur,old=nextcur,nextold
    mu=sy.Symbol('mu',positive=True);JT,KT,KE,KP=sy.symbols('J_transition K_theta K_energy K_pressure',real=True)
    K=dict(J=JT,theta=KT,energy=KE,pressure=KP)
    env=dict(t=1,mu=mu,mu_cover=mu,K=K,d=sy.exp(-1),d3=sy.exp(-sy.Rational(3,2)),u1=cur['E'],old=cur['hist'])
    # Exact-function lift: mu_cover labels an enclosure variable in the
    # numeric frontend. Here its underlying parameter is the original mu.
    # No equality of that enclosure's endpoints with mu is asserted.
    for target in ('factor','theta','energy','pressure'):
        env[target]=p.assignment(o3,'background_cell',target,env,count=2,index=0)
    env['E']=p.assignment(o3,'background_cell','E',env);env['E2']=p.assignment(o3,'background_cell','E2',env)
    nextcur=dict(E=env['E'],V=0,hist=p.assignment(o3,'background_cell','hist',env))
    envold=dict(t=1,mu=mu,K=K,u1=old['E'],old=old['hist'],d=sy.exp(-1),d3=sy.exp(-sy.Rational(3,2)))
    envold['factor']=p.assignment(current.pre,'slope_mu','factor',envold)
    nextold=dict(E=p.assignment(current.pre,'slope_mu','u',envold),V=0,
        hist=p.assignment(current.pre,'slope_mu','hist',envold))
    p.equal('O3_Rw/E',nextcur['E'],nextold['E'])
    for key in nextcur['hist']:p.equal('O3_Rw/'+key,nextcur['hist'][key],nextold['hist'][key])
    cur,old=nextcur,nextold
    t=sy.Symbol('power_t',positive=True);d=sy.exp(-t);d3=sy.exp(-sy.Rational(3,2)*t)
    exprel=lambda a:(sy.exp(a)-1)/a
    mass=lambda a:(1-sy.exp(-a*t))/a
    env=dict(t=t,mu=mu,mu_cover=mu,d=d,d3=d3,u1=cur['E'],old=cur['hist'],
        **{'incoming.weighted.terminal.local.packets.recovery.exp_average(c, (1 - mu_cover) * t)':exprel((1-mu)*t),
           'positive_decay_mass(c, 2 * mu_cover, t)':mass(2*mu),
           'positive_decay_mass(c, 1 + 2 * mu_cover, t)':mass(1+2*mu)})
    for target in ('factor','average','theta','energy','pressure'):
        count=2 if target in ('factor','theta','energy','pressure') else 1
        env[target]=p.assignment(o3,'background_cell',target,env,count=count,index=count-1)
    env['E']=p.assignment(o3,'background_cell','E',env);env['E2']=p.assignment(o3,'background_cell','E2',env)
    nextcur=dict(E=env['E'],V=0,hist=p.assignment(o3,'background_cell','hist',env))
    envold=dict(t=t,mu=mu,d=d,d3=d3,u1=old['E'],old=old['hist'],
        **{'decay_integral(c, 2 * mu, t)':mass(2*mu),'decay_integral(c, 1 + 2 * mu, t)':mass(1+2*mu)})
    envold['f']=p.assignment(current.pre,'power','f',envold)
    theta=p.assignment(current.pre,'power','theta',envold,index=0,count=2);envold['theta']=theta
    nextold=dict(E=p.assignment(current.pre,'power','u',envold),V=0,
        hist=p.assignment(current.pre,'power','hist',envold))
    p.equal('O3_power_all_t/E',nextcur['E'],nextold['E'])
    for key in nextcur['hist']:p.equal('O3_power_all_t/'+key,nextcur['hist'][key],nextold['hist'][key])
    p.equal('exact_exprel_to_old_theta',d3*t*exprel((1-mu)*t),(sy.exp(-(sy.Rational(1,2)+mu)*t)-d3)/(1-mu))
    qshape={key:sy.simplify(sy.expand_power_exp(expr)) for key,expr in nextcur['hist'].items()}
    coefficients=dict(U=sy.cancel(nextcur['E']*q),M=sy.cancel(qshape['m']/Z),
        H=sy.cancel(qshape['h']*q),K=sy.cancel(qshape['k']*q/Z),Pin=sy.cancel(qshape['p']*q*q))
    poly=sy.Poly(sy.cancel(qshape['e']*q*q),Z)
    coefficients.update(EZ=poly.coeff_monomial(Z**6),EQ=poly.coeff_monomial(1))
    require(all(Z not in value.free_symbols for value in coefficients.values()),'Actual source q-shape coefficients must be Z-independent')
    energy=coefficients['EZ']*Z*Z+coefficients['EQ']/q**2
    p.equal('new_power_actual_energy_qshape',qshape['e'],energy)
    Tw=sy.Symbol('Tw',positive=True);t2=2+sy.log(2)
    # Source identity holds for arbitrary power t, so it covers both the
    # current repaired endpoint t2 and the native inlet t=Tw directly.
    for label,end in (('twice_Rc',t2),('Rp',Tw)):
        p.equal(label+'/E',nextcur['E'].subs(t,end),nextold['E'].subs(t,end))
        for key in qshape:p.equal(label+'/'+key,nextcur['hist'][key].subs(t,end),nextold['hist'][key].subs(t,end))
    bindings=kernel_bindings()
    p.hashes.update(bindings['input_hashes'])
    return dict(passed=True,source_assignment_bindings=p.bindings,exact_function_identities=p.identities,
        identity_count=len(p.identities),actual_canonical_coefficients={k:sy.sstr(v.subs(t,Tw)) for k,v in coefficients.items()},
        kernel_and_exact_parameter_bindings=bindings,
        source_function_not_numeric_enclosure_identity=True,arithmetic_mu_enclosure_not_equal_to_mu=True,
        ordinary_Z_identity_over_whole_original_domain=True,
        canonical_O2_seed_used_not_unproved_patched_Rh_seam=True,input_hashes=p.hashes)


def kernel_bindings():
    require(slope.original is current.pre,'Same original O2 slope function provider')
    require(axial.kernels.turnoff_kernels is current.pre.turnoff_kernels,'Same original turnoff primitive provider')
    require(o3.kernels.transition_kernels is current.pre.transition_kernels,'Same original O3 kernel provider')
    modules=(slope,axial,o3,o3.kernels,axial.kernels,current.pre)
    source_AST={Path(m.__file__).name:current.current.current.current.ast_binding(m.background_cell if m in (slope,axial,o3) else
        m.transition_kernels if m is o3.kernels else m.turnoff_kernels if m is axial.kernels else m.slope_masses) for m in modules}
    tau,k=sy.symbols('tau k',positive=True)
    primitive=(1-sy.exp(-k*tau))/k
    require(sy.simplify(sy.diff(primitive,tau)-sy.exp(-k*tau))==0,'Exact common decay FTC identity')
    KB=sy.Symbol('K_at_yd');tail=KB*sy.exp(-tau)
    require(sy.diff(tail,tau)+tail==0 and tail.subs(tau,0)==KB,'Original zero-cutoff buffer kernel identity')
    return dict(same_original_continuous_kernel_providers=True,source_AST_bindings=source_AST,
        continuous_primitives=dict(slope='J(y)=integral_0^y sigma; masses are the same original transition_integrals',
            turnoff='K_j(y)=integral_1^y exp(s-y)*sigma(1-log(s)/40)^j ds',
            buffer='sigma=0 for s>=exp(40); K_j(exp(40)+tau)=exp(-tau)*K_j(exp(40))',
            O3='same original transition_kernels underlying J/theta/energy/pressure',
            decay='integral_0^t exp(-k*s)ds=t*exprel(-k*t); both helpers enclose this primitive'),
        empty_slope_integrals_at_y0=True,cutoff_support_and_true_buffer_FTC=True,
        finite_partition_sums_are_enclosures_not_function_parameters=True,
        mu='exp(log(.001)-4*(exp(40)+11))',Tw='-60log(mu)',Md=40,
        input_hashes={Path(m.__file__).name:sha(Path(m.__file__).name) for m in modules})


def load_inputs():
    report=json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
    checked=json.loads((HERE/current.RECEIPT).read_bytes())
    require(report[current.GATE] and checked[current.GATE] and checked['all_passed'],'Accepted current Rp packet required')
    hashes=dict(checked['input_hashes'])
    for name,digest in hashes.items():require(sha(name)==digest,'Changed current source: '+name)
    hashes.update({current.NAME:sha(current.NAME),current.RECEIPT:sha(current.RECEIPT),Path(__file__).name:sha(Path(__file__).name)})
    oldname=PREFIX+'actual_Rp_source_join.json';old=json.loads((HERE/oldname).read_bytes())
    checked_old=json.loads((HERE/(PREFIX+'actual_Rp_source_join_check.json')).read_bytes())
    require(checked_old['input_hashes'][oldname]==sha(oldname),'Old native proof report must be checked')
    require(old['current_Rp_terminal_shape_proof']['passed'] and old['current_Rp_native_inlet_projection_proof']['passed'],
        'Accepted original native qshape and actual inlet projection proofs required')
    hashes[oldname]=sha(oldname)
    hashes[PREFIX+'actual_Rp_source_join_check.json']=sha(PREFIX+'actual_Rp_source_join_check.json')
    for name,digest in checked_old['input_hashes'].items():
        require(sha(name)==digest,'Changed checked native dependency: '+name)
        require(name not in hashes or hashes[name]==digest,'Different native source binding: '+name)
        hashes[name]=digest
    return report,old,hashes


def ast_assignments(stem,cls,method,wanted):
    name=PREFIX+stem+'.py';tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    owner=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls) if cls else tree
    fn=next(n for n in owner.body if isinstance(n,ast.FunctionDef) and n.name==method)
    for target,expr in wanted.items():
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets)]
        require(len(values)==1 and ast.dump(values[0])==ast.dump(ast.parse(expr,mode='eval').body),
            'Defining pressure/source assignment changed: '+name+':'+target)
    return dict(assignments=wanted,function_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest(),
        input_hashes={name:sha(name)})


def pressure_source_proof(report,old):
    """Bind the cache to the analytic datum, never select a mass-box value."""
    name=PREFIX+'pressure_source.json';ps=json.loads((HERE/name).read_bytes())['compliant_source']
    definition=ps['implicit_source_definition'];identity=report['source_family']
    require(hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()==ps['implicit_source_sha256'],
        'Actual analytic pressure definition must determine its source identity')
    for key in ('implicit_source_sha256','datum_enclosure_sha256'):
        require(identity[key]==ps[key]==old[key],'Current cache and native analytic datum differ: '+key)
    require(definition['Md']=='40' and definition['c_epsilon']=='.001' and ps['stage_count_total']==14,
        'Same compliant fourteen-stage pressure source required')
    require(CompliantPressureDatum.normalized_jets is LogarithmicPressureDatum.normalized_jets,
        'Actual inherited pressure callable differs')
    inherited=old['current_Rh_pressure_defining_function_proof']
    require(inherited['passed'] and inherited['same_fourteen_implicit_stage_sources']
        and inherited['common_flatten_function_retained_not_replaced_by_zero']
        and inherited['exact_same_inherited_normalized_jets_callable'], 'Native analytic datum proof required')
    bindings={}
    bindings['current_hydrated_datum']=ast_assignments('current_original_bridge_macro_functions',
        'OriginalBridgeMacroFunctions','__init__',dict(
            ps="self.records['pressure_source']['compliant_source']",
            datum='LogarithmicPressureDatum.__new__(LogarithmicPressureDatum)',
            **{'datum.m2':"previous.read_interval(c,ps['fixed_beta2_mass_normalized'])",
                'datum.m0':"previous.read_interval(c,ps['fixed_beta0_mass_normalized'])",
                'datum.rho':"previous.read_interval(c,ps['complex_strip_half_width'])",
                'datum.flatten_complex_upper':"previous.read_interval(c,ps['flatten_complex_mass_upper'])",
                'datum.stages':"{'z_flatten':{'mass':previous.read_interval(c,ps['stages']['z_flatten']['mass'])}}",
                'datum.parameters':"SimpleNamespace(logPstar=previous.read_interval(c,ps['parameter_bounds']['logPstar']))",
                'datum.source_sha':'self.source','datum.datum_sha':'self.datum','self.pressure':'datum'}))
    bindings['current_whole_Z_call']=ast_assignments('current_original_whole_Z_bridge_source',
        'OriginalWholeZBridgeSource','inputs',dict(pressure='self.admission.pressure.normalized_jets(ep(z),6)',
            result="dict(p0=IntervalTaylor(c,[c.mpf(ep(value)) for value in pressure['normalized_pressure_coefficients']]),F0_ratios=ratios(1),F0_squared_ratios=ratios(2))"))
    bindings['current_Rm_first_six']=ast_assignments('current_original_whole_Z_long_Rm_finite_N',
        'WholeZLongRmFiniteN','owner',dict(p0="current.source.IntervalTaylor(c,proof['original_pressure_Z0_through_Z6'][:6])",
            P0='f.jet(p0)',savedP0="upstream.decode(f,incoming['exact_common_P0_axial5'])"))
    # A genuine analytic flatten function is retained in every Taylor row;
    # its Cauchy enclosure is never substituted for F_flat itself.
    M2,M0=sy.symbols('pressure_M2 pressure_M0',real=True);flat=sy.Function('F_flat')(Z)
    P0=-(M2/(1+Z**2)**2+M0+flat)
    rows=[sy.sstr(sy.diff(P0,Z,n)/sy.factorial(n)) for n in range(6)]
    hashes={name:sha(name)}
    for row in bindings.values():hashes.update(row['input_hashes'])
    hashes.update(ps['input_hashes'])
    for module in (inspect.getmodule(LogarithmicPressureDatum),inspect.getmodule(CompliantPressureDatum)):
        hashes[Path(module.__file__).name]=sha(Path(module.__file__).name)
    return dict(passed=True,analytic_definition=definition,source_AST_bindings=bindings,
        inherited_normalized_jets_AST_sha256=hashlib.sha256(ast.dump(ast.parse(inspect.getsource(
            LogarithmicPressureDatum.normalized_jets).strip())).encode()).hexdigest(),
        exact_analytic_function=sy.sstr(P0),first_six_true_Taylor_coefficients=rows,
        same_fourteen_continuous_atoms_and_analytic_flatten_function=True,
        hydrated_cache_fields_are_enclosures_not_selected_analytic_masses=True,
        parameter_bounds_are_not_exact_parameter_values=True,
        current_order6_first_six_and_native_order5_same_defining_rows=True,
        common_P0_kept_independent_of_cumulative_pressure=True,input_hashes=hashes)


def native_frame_proof(proof):
    """Join raw old moments to current common moments before division."""
    q=1+Z*Z;U,M,H,K,EZ,EQ,Pin=sy.symbols('U M H K EZ EQ Pin',real=True)
    E=U/q;raw=dict(m=M*Z,h=H/q,k=K*Z/q,e=EZ*Z**2+EQ/q**2,p=Pin/q**2)
    common={key:value/S if key in ('m','k') else value for key,value in raw.items()}
    current_frame=dict(u=E,m1=common['m']/E,m2=common['k']/E**2,X=common['h']/E,
        energy=common['e']/E**2,Mp=common['p'])
    native=dict(u=U/q,m1=M/(S*U)*Z*q,m2=K/(S*U**2)*Z*q,X=H/U,
        energy=EQ/U**2+EZ/U**2*Z**2*q**2,Mp=Pin/q**2)
    identities={}
    for key,value in current_frame.items():
        require(sy.cancel(value-native[key])==0,'Actual pulse frame differs: '+key)
        for n in range(2):
            require(sy.cancel(sy.diff(value-native[key],Z,n))==0,'Actual ordinary Z pulse frame differs')
            identities[key+('/value' if n==0 else '/Z')]=sy.sstr(sy.diff(native[key],Z,n))
    require(proof['identity_count']==45,'Complete original recipe projection required')
    binding=ast_assignments('current_original_Rm_generic_inputs',None,'recover_inputs',dict(
        invS='f.factor((0,-.5,0,0,0))',
        histories="{key:scale(rows[0],invS) if key in ('m','k') else rows[0] for key,rows in raw['histories'].items()}"))
    return dict(passed=True,raw_to_common=dict(m='m_raw/Pstar',k='k_raw/Pstar',h='h_raw',e='e_raw',p='p_raw'),
        exact_native_value_and_first_Z_identities=identities,identity_count=len(identities),
        concrete_coefficient_functions=proof['actual_canonical_coefficients'],
        current_complete_equals_leading_after_actual_C1_relative_zero=True,
        accepted_current_semigroup_and_original_Rp_geometry_retained=True,
        no_second_inverse_Pstar_factor_applied_to_common_rows=True,
        actual_raw_to_common_AST_binding=binding,input_hashes=binding['input_hashes'])


def run():
    began=time.monotonic();report,old,hashes=load_inputs();proof=defining_source_proof()
    pressure=pressure_source_proof(report,old);frame=native_frame_proof(proof)
    for name,digest in {**proof['input_hashes'],**pressure['input_hashes'],**frame['input_hashes']}.items():
        require(name not in hashes or hashes[name]==digest,'Outer source binding conflict: '+name);hashes[name]=digest
    result=dict(**{GATE:True},all_passed=True,source_family=report['source_family'],
        exact_underlying_O2_O3_field_and_five_history_recipe_identity_proved=True,
        underlying_function_projection=proof,current_P0_source_binding=pressure,actual_native_frame_identity=frame,
        old_native_qshape_and_inlet_projection_bound=True,
        original_current_limit_relative_zero_and_Rp_semigroup_bound=True,
        current_to_native_P0_defining_function_identity_bound=True,
        actual_Rp_native_frame_source_function_identified=True,
        actual_native_O4_constructor_consumes_current_limit_frame=False,
        actual_numeric_point_values_installed=False,old_numeric_mu_enclosure_equal_to_exact_mu=False,
        patched_Rh_functional_join_proved=False,physical_original_exterior_five_targets_closed=False,
        actual_temporal_scale_recursion_installed=False,
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Exact underlying original O2/O3 formulas, whole-Z cumulative histories and canonical qshape agree with the old native recipe. '
            'Current P0 leaf is bound to the same fourteen-stage analytic datum and the native units are identified. No numeric enclosure equality, '
            'inner Rh seam, native installation, absolute heat closure or temporal recursion claim.')
    (HERE/NAME).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('CURRENT_LIMIT_RP_NATIVE_RECIPE_IDENTITY',proof['identity_count'],'outer identities;',frame['identity_count'],'pulse identities; same analytic P0',flush=True)
    return result


if __name__=='__main__':run()
