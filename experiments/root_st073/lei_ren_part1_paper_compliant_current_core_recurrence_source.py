"""Exact production radial-scaling and leading-axis source equations.

This proves equations, not equality of two nonlinear enclosure algorithms.
The current symmetric S jets, finite residual/tail admission and core/first
functional join remain separate obligations. No temporal recursion claim.
"""
import ast
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import sympy as s
from lei_ren_part1_paper_functional_core_step import advance_one, initial_rows
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
RECEIPT=PREFIX+'compliant_current_core_recurrence_source_check.json'
CTX=SimpleNamespace(mpf=s.sympify)


def sha(name):return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def function(module,name):
    tree=ast.parse((HERE/(PREFIX+module+'.py')).read_bytes())
    rows=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name]
    if len(rows)!=1:raise ValueError('Unique production function required: '+module+'.'+name)
    return rows[0]


def binding(module,name,target,expression,allow_other_assignments=False):
    fn=function(module,name)
    normalized_target=ast.unparse(ast.parse(target,mode='eval').body)
    values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)==normalized_target for t in n.targets)]
    if not values:
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.keyword) and n.arg==target]
    expected=ast.dump(ast.parse(expression,mode='eval').body)
    matching=[v for v in values if ast.dump(v)==expected]
    if len(matching)!=1 or (len(values)!=1 and not allow_other_assignments):
        raise ValueError('Production expression changed: '+module+'.'+name+': '+target)
    return matching[0]


def full_production_AST_scaling():
    original=copy.deepcopy(function('functional_core_step','advance_one'))
    actual=copy.deepcopy(function('logarithmic_core_step','advance_scaled_one'))
    replacements=[
        ("add(rhs_a,scale(mul(W[i],a[n-i]),n-i+1),mul(H[i],a_derivative))",
         "add(rhs_a,scale(mul(W[i],a[n-i]),eps*(n-i+1)),mul(H[i],add(scale(diff(a[n-i]),eps),mul(ell,a[n-i]))))"),
        ("add(rhs_u,scale(mul(W[i],u[n-i]),n-i),mul(H[i],diff(u[n-i])))",
         "add(rhs_u,scale(mul(W[i],u[n-i]),eps*(n-i)),scale(mul(H[i],diff(u[n-i])),eps))"),
        ("add(rhs_a,scale(add(a[n],scale(mul(z,aa_u),-2)),delta_iv/2))",
         "add(rhs_a,scale(add(a[n],scale(mul(z,aa_u),-2)),eps*delta_iv/2))"),
        ("add(rhs_u,scale(add(u[n],scale(mul(z,uu),-2)),(1+delta_iv)/2),mul(d,diff(pressure[n])),scale(mul(z,pressure[n]),-2*(1+delta_iv)))",
         "add(rhs_u,scale(add(u[n],scale(mul(z,uu),-2)),eps*(1+delta_iv)/2),mul(d,diff(pressure[n])),scale(mul(z,pressure[n]),-2*(1+delta_iv)))")]
    table={ast.dump(ast.parse(a,mode='eval').body):ast.parse(b,mode='eval').body for a,b in replacements}
    seen={key:0 for key in table};removed=0
    derivative=ast.dump(ast.parse('add(diff(a[n-i]),mul(ell,a[n-i]))',mode='eval').body)
    class Normalize(ast.NodeTransformer):
        def visit_Assign(self,node):
            nonlocal removed
            if len(node.targets)==1 and ast.unparse(node.targets[0])=='a_derivative':
                if ast.dump(node.value)!=derivative:raise ValueError('Original amplitude derivative changed')
                removed+=1;return None
            key=ast.dump(node.value)
            if key in table:seen[key]+=1;node.value=copy.deepcopy(table[key])
            return node
    original=Normalize().visit(original)
    if removed!=1 or any(v!=1 for v in seen.values()):raise ValueError('All four scaling transformations required exactly once')
    original.name='advance_scaled_one'
    original.args.args.append(ast.arg(arg='epsilon',annotation=ast.Name(id='Any',ctx=ast.Load())))
    guard=ast.parse("eps = _as_interval(ctx, epsilon)\nif not (eps > 0 and eps <= 1):\n    raise ValueError('epsilon must lie in (0,1]')").body
    positions=[i for i,n in enumerate(original.body) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)=='convert' for t in n.targets)]
    if len(positions)!=1:raise ValueError('Original coercion position changed')
    original.body[positions[0]+1:positions[0]+1]=guard
    for fn in (original,actual):
        if not isinstance(fn.body[0],ast.Expr) or not isinstance(fn.body[0].value,ast.Constant):
            raise ValueError('Production function docstring required')
        fn.body=fn.body[1:]
    if ast.dump(original)!=ast.dump(actual):
        raise ValueError('Scaled production AST differs beyond derived normalization edits')
    return dict(entire_production_function_AST_equal_after_four_derived_scaling_edits=True,
        one_original_amplitude_derivative_inlined=True,original_helpers_validations_loops_and_all_other_terms_unchanged=True,
        positive_epsilon_domain_guard_retained=True,passed=True)


def all_radial_indices_homogeneity():
    n,i=s.symbols('n i',integer=True,nonnegative=True)
    # Lambda is Z-independent. Taylor differentiation and convolution retain
    # radial weight; the constant W0/H0 terms have weight zero at i=0.
    terms=[('A: W_i A_(n-i)',i+n-i,n+1,1),
        ('A: H_i dZ A_(n-i)',i+n-i,n+1,1),
        ('A: H_i ell A_(n-i)',i+1+n-i,n+1,0),
        ('A: delta A_n',n,n+1,1),('A: delta Z sum U_i A_(n-i)',i+n-i,n+1,1),
        ('Uz: W_i U_(n-i)',i+n-i,n+1,1),('Uz: H_i dZ U_(n-i)',i+n-i,n+1,1),
        ('Uz: affine U_n',n,n+1,1),('Uz: Z sum U_i U_(n-i)',i+n-i,n+1,1),
        ('Uz: d dZ P_n',n+1,n+1,0),('Uz: Z P_n',n+1,n+1,0),
        ('Uz: S Z sum A_i A_(n-1-i), n>=1',2+i+n-1-i,n+1,0),
        ('P: S sum A_i A_(n-i)',2+i+n-i,n+2,0)]
    rows=[]
    for label,weight,normalization,eps_power in terms:
        if s.simplify(weight-normalization+eps_power)!=0:raise ArithmeticError('Incorrect homogeneous source weight')
        rows.append(dict(term=label,original_Lambda_weight=str(s.expand(weight)),
            output_Lambda_normalization=str(normalization),resulting_epsilon_power=eps_power))
    return dict(all_thirteen_additive_term_families=rows,arbitrary_integer_n_and_convolution_index=True,
        arbitrary_axial_Taylor_length=True,Lambda_independent_of_Z=True,
        original_L_axial_constant_nonzero_required=True,
        W0_H0_constant_terms_have_zero_radial_weight=True,
        inverse_L_and_radial_denominators_have_zero_Lambda_weight=True,
        entire_AST_binding_required_before_homogeneity_claim=True,passed=True)


def independent_formal_production_jets():
    # lambda>1: SymPy proves the actual unmodified runtime guard. The exact
    # rational identities also extend to lambda=1. No production is patched.
    lam_minus_one=s.Symbol('Lambda_minus_one',positive=True)
    eps=1/(1+lam_minus_one);z,delta=s.symbols('Z delta')
    counts={}
    for n,K in ((0,4),(1,3),(2,3),(3,2)):
        size=n+K+1
        vector=lambda tag,length:[s.Symbol(tag+'_'+str(k)) for k in range(length)]
        scaled_fixed={key:vector(tag,size) for key,tag in
            (('ell_Z_taylor','ell'),('S_Z_taylor','S'),('U0_Z_taylor','U0'),('P0_Z_taylor','P0'))}
        fixed={key:[v/(eps if key in ('ell_Z_taylor','P0_Z_taylor') else eps**2 if key=='S_Z_taylor' else 1)
            for v in values] for key,values in scaled_fixed.items()}
        scaled={label:[vector(label+str(j),size-j) for j in range(n+1)] for label in ('A','Uz','P')}
        unscaled={label:[[v/eps**(j+(label=='P')) for v in row] for j,row in enumerate(rows)]
            for label,rows in scaled.items()}
        advance_one(CTX,fixed,unscaled,n,z,delta)
        advance_scaled_one(CTX,scaled_fixed,scaled,n,z,delta,eps)
        count=0
        for label in ('A','Uz','P'):
            for left,right in zip(scaled[label][-1],unscaled[label][-1]):
                if s.cancel(left-eps**(n+1+(label=='P'))*right)!=0:
                    raise ArithmeticError('Independent production formal jet differs: '+label+' n='+str(n))
                count+=1
        counts[str(n)]=count
    return dict(independent_formal_axial_coefficients=True,unmodified_production_functions_executed=True,
        original_Cauchy_Taylor_convention=True,new_row_scalar_identities_by_radial_index=counts,
        total_independent_formal_scalar_identities=sum(counts.values()),
        finite_fixture_is_not_the_all_n_proof=True,passed=True)


def leading_axis_source_equations():
    z,j,delta=s.symbols('Z j delta',real=True);sigma,eps=s.symbols('sigma epsilon',positive=True)
    one=s.Integer(1);self=SimpleNamespace(j=j,delta=delta,sigma=sigma)
    env={'z':z,'one':one,'self':self};core='compliant_core_physical_field'
    for target,expr in [('u','4*z+self.j'),('d','one-z*z'),('L','one-(z*z)*self.delta'),
        ('H','z*((1-self.delta)/2)+d*u')]:
        node=binding(core,'axis_inputs',target,expr);env[target]=eval(compile(ast.Expression(node),'<original core axis>','eval'),{},env)
    u,d,L,H=[env[key] for key in ('u','d','L','H')]
    P,Pz=s.symbols('physical_P0 physical_P0_Z',real=True);env.update(physicalP=P,pz=Pz)
    g=binding(core,'axis_inputs','g','-((1-2*z*u)*u*((1+self.delta)/2)+4*H+d*pz-z*physicalP*(2*(1+self.delta)))')
    env['g']=eval(compile(ast.Expression(g),'<original core pressure drive>','eval'),{},env)
    slope=binding(core,'axis_inputs','slope','-g/(2*L)')
    slope=eval(compile(ast.Expression(slope),'<original core axis slope>','eval'),{},env)
    binding(core,'axis_inputs','beta','((z*self.j+3-self.delta/2)/L).truncate(5)')
    binding(core,'axis_inputs','chi','h2/denominator',allow_other_assignments=True)
    binding(core,'axis_inputs','gradient','L*H/denominator')
    binding(core,'axis_inputs','ell','[-self.Lambda*gradient[k] for k in range(5)]')
    binding(core,'normalized_jets','phi[index]','-(source[\'chi\'][k]+self.epsilon*source[\'beta\'][k])*math.factorial(k)/4',allow_other_assignments=True)
    binding(core,'normalized_jets','model',"(source['slope'][k]*math.factorial(k))*(rho if i==0 else 1 if i==1 else 0)")
    binding(core,'normalized_jets','psi[index]','model',allow_other_assignments=True)
    binding(core,'normalized_jets','Vgrid[gridkey(i,k)]','base+self.epsilon*psi[gridkey(i,k)]')
    branches=[node for node in ast.walk(function(core,'normalized_jets')) if isinstance(node,ast.If)
        and ast.dump(node.test)==ast.dump(ast.parse('axis and i==1',mode='eval').body)]
    if len(branches)!=1 or not any(isinstance(node,ast.Assign)
        and any(ast.unparse(target)=='psi[index]' for target in node.targets)
        and ast.unparse(node.value)=='model' for node in branches[0].body):
        raise ValueError('Exact axis radial-one slope override changed')
    chi=H**2/(H**2+sigma**2);gradient=L*H/(H**2+sigma**2)
    beta=(z*j+3-delta/2)/L;S=s.Symbol('formal_exact_scaled_S',positive=True)
    # Positive symbolic epsilon with an actual domain-valid parametrization.
    t=s.Symbol('Lambda_minus_one',positive=True);e=1/(1+t)
    fixed={'ell_Z_taylor':[-gradient,0,0],'S_Z_taylor':[S,0,0],
        'U0_Z_taylor':[u,4,0],'P0_Z_taylor':[e*P,e*Pz,0]}
    rows=initial_rows(CTX,fixed,z,1,required_depth=1)
    advance_scaled_one(CTX,fixed,rows,0,z,delta,e)
    expected={'A':-(chi+e*beta)/4,'Uz':e*slope,'P':S}
    for label,value in expected.items():
        if s.cancel(rows[label][1][0]-value)!=0:raise ArithmeticError('Leading exact source equation differs: '+label)
    return dict(original_axis_u_d_L_H_g_slope_and_normalized_field_AST_bound=True,
        Atilde1='-(chi+epsilon*beta)/4',Utilde1='epsilon*(-g/(2*L))',
        Ptilde1='formal_exact_scaled_S; exact epsilon^2*F0^2 source binding NOT admitted',
        A_and_U_equations_identical_as_Z_functions=True,
        finite_Z_derivatives_follow_where_original_analytic_denominators_nonzero=True,
        actual_symmetric_S_seed_not_promoted_to_exact_source=True,passed=True)


def pressure_and_seed_bindings():
    binding('compliant_core_transfer','run','lam','ctx.exp(logLambda)')
    binding('compliant_core_transfer','run','eps','ctx.exp(-logLambda)')
    binding('compliant_core_physical_field','__init__','self.epsilon',"read(major,'epsilon')")
    binding('compliant_core_physical_field','__init__','self.datum',"CompliantPressureDatum('40',160)")
    binding('compliant_core_physical_field','axis_inputs','datum','self.datum.normalized_jets(endpoints(c.mpf(Zvalue)),6)')
    binding('compliant_core_physical_field','axis_inputs','pressure',"IntervalTaylor(c,[c.mpf(endpoints(v)) for v in datum['normalized_pressure_coefficients']])")
    binding('compliant_core_physical_field','axis_inputs','physicalP','pressure*c.exp(2*self.logP)')
    binding('compliant_core_physical_field','axis_inputs','pz','IntervalTaylor(c,[(k+1)*physicalP[k+1] for k in range(6)])')
    seed='compliant_core_coefficient_rebuild'
    binding(seed,'seed','u','4*z+self.core.j')
    binding(seed,'seed','L','one-z*z*self.core.delta')
    binding(seed,'seed','H','z*((1-self.core.delta)/2)+(one-z*z)*u')
    binding(seed,'seed','gradient','L*H/IntervalTaylor(c,coefficients)')
    binding(seed,'seed','coefficients[0]','H[0]**2+self.core.sigma**2')
    binding(seed,'seed','pressure','self.core.datum.normalized_jets(z0,order)')
    binding(seed,'seed','scale','c.exp(2*self.core.logP-self.core.logLambda)')
    binding(seed,'seed','fixed',"dict(ell_Z_taylor=[-v for v in gradient.coefficients],S_Z_taylor=S,U0_Z_taylor=list(u.coefficients),P0_Z_taylor=[scale*v for v in pressure['normalized_pressure_coefficients']])")
    l,p=s.symbols('logLambda logP',real=True)
    if s.simplify(s.exp(-l)*s.exp(l)-1)!=0 or s.simplify(s.exp(2*p-l)-s.exp(-l)*s.exp(p)**2)!=0:
        raise ArithmeticError('Exact seed parameter units differ')
    return dict(Lambda_times_core_epsilon_exactly_one=True,
        P0_step_equals_core_epsilon_times_Pstar_squared_times_normalized_pressure=True,
        seed_ell_is_already_scaled_negative_gradient=True,no_second_epsilon_applied_to_ell=True,
        same_core_datum_normalized_jets_method_bound=True,
        core_radial_epsilon_distinct_from_pressure_datum_parameter_epsilon=True,
        pressure_primitive_4C_equality_not_proved_here=True,passed=True)


def run():
    parent_names=[PREFIX+'compliant_current_actual_bridge_mixed_C4_check.json',
        PREFIX+'compliant_core_coefficient_rebuild_check.json']
    parents=[json.loads((HERE/name).read_bytes()) for name in parent_names]
    if not all(row['all_passed'] for row in parents):raise ValueError('Accepted current core coefficient parents required')
    identity_keys=('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')
    if any(parents[0][key]!=parents[1][key] for key in identity_keys):
        raise ValueError('Current bridge and fresh core coefficient source identities differ')
    hashes={}
    for name,row in zip(parent_names,parents):
        for dependency,digest in row['input_hashes'].items():
            if dependency in hashes and hashes[dependency]!=digest:raise ValueError('Conflicting accepted source hashes')
            if sha(dependency)!=digest:raise ValueError('Changed accepted source: '+dependency)
            hashes[dependency]=digest
        hashes[name]=sha(name)
    for name in [Path(__file__).name,PREFIX+'functional_core_step.py',PREFIX+'logarithmic_core_step.py',
        PREFIX+'compliant_core_coefficient_rebuild.py',PREFIX+'compliant_core_physical_field.py',PREFIX+'compliant_core_transfer.py']:
        hashes[name]=sha(name)
    result=dict(**{key:parents[0][key] for key in identity_keys},
        entire_production_AST_scaling=full_production_AST_scaling(),
        all_radial_indices_homogeneity=all_radial_indices_homogeneity(),
        independent_formal_production_jets=independent_formal_production_jets(),
        leading_axis_source_equations=leading_axis_source_equations(),
        pressure_and_seed_bindings=pressure_and_seed_bindings(),
        current_core_radial_production_scaling_certified=True,
        current_core_leading_axis_A_U_source_equations_certified=True,
        current_core_pressure_seed_scaling_certified=True,
        exact_scaled_S_source_jet_admitted=False,finite_rows_and_tails_bound_to_same_nonlinear_fixed_point=False,
        current_core_bridge_functional_mixed4_join_certified=False,all_current_bridge_functional_interfaces_certified=False,
        actual_point_moment_history_recovered=False,full_cartesian_vector_derivatives_certified=False,
        global_completed_tensor_admissibility=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
        remaining_dependency='Exact scaled S jet and selected complex-tube Cauchy bounds -> finite residual/tail common fixed-point provenance -> original core/first production ODE join',
        all_scoped_checks_passed=True,all_passed=True,input_hashes=hashes)
    (HERE/RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS full production radial scaling and leading axis A/Uz source equations; exact S/tails/core-first remain open',flush=True)
    return result


if __name__=='__main__':run()
