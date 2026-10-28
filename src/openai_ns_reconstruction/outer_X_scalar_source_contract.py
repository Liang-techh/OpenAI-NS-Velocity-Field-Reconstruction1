"""Check six actual source branches and their radial scalar dependencies.

This is a bounded AST/source-role audit; it does not run any original source
constructor, inspect graph edges, evaluate integrals or compare sampled data.
"""
import ast
import base64
import hashlib
import json
import sympy as s

OUTER='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py'
OUTER_SHA='4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3'
REPORT='experiments/root_st073/lei_ren_part1_paper_compliant_current_outer_X_common_factor.json'
REPORT_SHA='5b3867f7876e2cde17050a3b989cad395b651197d7dd71d377654bd0768c25e7'
RESULT='tmp/outer_X_common_factor_v1/result.json'
RESULT_SHA='b4e0d08ad407ea774546ab56e50f1db252ff088d8dac83ab97647542c6b599c3'
CHARTS=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')
PREV=("self.leading_Rh['fields']['E']",)+tuple("self.background_endpoints[%r]['fields']['E']"%c for c in CHARTS[:-1])
FACTOR=("g.unary('exp',g.mul(g.constant('1/10'),g.add(x,g.constant(5))))",
        "g.unary('exp',g.sub(g.mul(g.constant('1/10'),x),g.mul(g.constant('3/5'),J)))",
        "g.unary('exp',g.mul(g.constant('-1/2'),distance))",
        "g.unary('exp',g.mul(g.constant('-1/2'),x))",
        "g.unary('exp',g.sub(g.mul(g.constant('-1/2'),x),g.mul(self.mu,J)))",
        "g.unary('exp',g.neg(g.mul(g.add(g.constant('1/2'),self.mu),x)))")


def exact_expr(text):return ast.dump(ast.parse(text,mode='eval').body,include_attributes=False)


def checked_raw(root,path,expected):
    raw=(root/path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('input hash changed: '+path)
    return raw


def find_assignment(statements,target):
    found=[n.value for n in statements if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==target for t in n.targets)]
    if len(found)!=1:raise ValueError('expected one '+target+' assignment')
    return found[0]


def recover(root):
    source=checked_raw(root,OUTER,OUTER_SHA).decode('utf-8').replace('\r\n','\n')
    tree=ast.parse(source);cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentOriginalTrueOuterC6C5Functions')
    init=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    loops=[n for n in ast.walk(init) if isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id=='name' and isinstance(n.iter,ast.Name) and n.iter.id=='CHARTS']
    if len(loops)!=1:raise ValueError('outer chart loop is not unique')
    loop=loops[0]; branch=next(n for n in loop.body if isinstance(n,ast.If)); bodies={}
    while isinstance(branch,ast.If):
        test=branch.test
        if not (isinstance(test,ast.Compare) and isinstance(test.left,ast.Name) and test.left.id=='name' and len(test.ops)==1 and isinstance(test.ops[0],ast.Eq) and len(test.comparators)==1 and isinstance(test.comparators[0],ast.Constant)):
            raise ValueError('unexpected chart branch selector')
        name=test.comparators[0].value
        if name in bodies:raise ValueError('duplicate chart branch')
        bodies[name]=branch.body
        branch=branch.orelse[0] if len(branch.orelse)==1 and isinstance(branch.orelse[0],ast.If) else branch.orelse
    bodies['O3_power']=branch
    if set(bodies)!=set(CHARTS):raise ValueError('chart branch names differ')
    rows=[]
    for i,name in enumerate(CHARTS):
        expr=find_assignment(bodies[name],'E')
        if not (isinstance(expr,ast.Call) and ast.unparse(expr.func)=='a.scale' and len(expr.args)==2):raise ValueError('not a scale: '+name)
        if ast.dump(expr.args[0],include_attributes=False)!=exact_expr(PREV[i]) or ast.dump(expr.args[1],include_attributes=False)!=exact_expr(FACTOR[i]):
            raise ValueError('source branch predecessor or factor differs: '+name)
        rows.append({'chart':name,'predecessor':PREV[i],'factor':FACTOR[i]})
    integral="controls.integral(g,sigma(g,t),g.nodes[t.node]['name'],g.zero,x,exact_original_sigma_integral=True)"
    for name in ('O2_slope','O3_transition'):
        if ast.dump(find_assignment(bodies[name],'J'),include_attributes=False)!=exact_expr(integral):raise ValueError('sigma integral binding differs: '+name)
    distance="g.sub(g.unary('exp',g.mul(g.constant(40),x)),g.one)"
    if ast.dump(find_assignment(bodies['O2_axial'],'distance'),include_attributes=False)!=exact_expr(distance):raise ValueError('axial distance differs')
    anchors=["x=g.symbol('true_outer_'+name);offset,jac,domain=geometry(g,x,name)",
             "self.mu=g.unary('exp',self.logmu)",
             "self.logmu=g.sub(g.unary('log',g.constant('1/1000')),g.mul(g.constant(4),g.add(g.unary('exp',g.constant(40)),g.constant(11))))",
             "jac=g.mul(c(40),g.unary('exp',g.mul(c(40),x))) if chart=='O2_axial' else g.one",
             "return dict(m=V,h=E,k=a.mul(E,V)",
             "body=a.scale(a.at(density[key],x,t),g.mul(measure,kernel))",
             "memory[key]=decay;out[key]=a.add(a.scale(incoming[key],decay),local[key])"]
    if any(anchor not in source for anchor in anchors):raise ValueError('radial/source binding anchor missing')
    geometry=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='geometry')
    if any((isinstance(n,ast.Name) and n.id=='Z') or (isinstance(n,ast.Constant) and n.value=='Z') for n in ast.walk(geometry)):
        raise ValueError('explicit Z in geometry syntax')
    report=json.loads(checked_raw(root,REPORT,REPORT_SHA));data=base64.b64decode(report['artifacts'][RESULT]['base64'],validate=True)
    if hashlib.sha256(data).hexdigest()!=RESULT_SHA:raise ValueError('prior result hash differs')
    prior=json.loads(data)
    if prior.get('status')!='returned' or not prior['result']['checks_passed']:raise ValueError('prior factor model not returned')
    roles=prior['result']['saved_chart_roles']
    if [r['chart'] for r in roles]!=list(CHARTS):raise ValueError('saved chart roles differ')
    x,z=s.symbols('x z',real=True)
    checks={'source_exact_branch_count':len(rows)==6,
            'axial_jacobian':s.simplify(s.diff(6+s.exp(40*x),x)-40*s.exp(40*x))==0,
            'seed_denominator_positive':s.ask(s.Q.positive(1+z*z)) is True,
            'seed_exp_positive':s.ask(s.Q.positive(s.exp(-s.Rational(1,2)))) is True}
    return {'schema':'outer_X_scalar_source_contract_v1','status':'source_branch_scalar_contract_under_primitive_semantics',
            'checks':checks,'checks_passed':all(checks.values()),'branches':rows,'saved_chart_roles':roles,
            'source_input_pins':{'outer_source_sha256':OUTER_SHA,'prior_result_sha256':RESULT_SHA},
            'syntax_findings':{'two_sigma_integrals_are_radial':True,'mu_is_exp_of_fixed_constant_recipe':True,
                               'geometry_has_no_explicit_Z':True,'axial_physical_jacobian_is_40_exp_40x':True,
                               'h_density_is_E_and_Duhamel_keeps_inlet_memory':True},
            'conditions':['source graph primitives implement the declared exact scalar/function semantics',
                'sigma is the same real bounded unary source and radial integration is capture-free',
                'saved leading Rh source identity and same-graph chart endpoint composition hold',
                'scalar factors are evaluated at finite real chart coordinates and parameters',
                'terminal relative E/h corrections vanish as functions on the full Z domain'],
            'conditional_conclusion':'each leading E scale is a positive exponential of a finite real radial scalar; source-model h/E factorization follows under the declared Duhamel semantics',
            'actual_graph_to_function_semantics_proved':False,'actual_terminal_function_identity_proved':False,
            'actual_compact_P4_P5_zero_proved':False,'actual_X4_X5_zero_proved':False,
            'source_values_computed':False,'integrals_evaluated':False,'old_errors_removed':False,'runtime_model_metadata':None}
