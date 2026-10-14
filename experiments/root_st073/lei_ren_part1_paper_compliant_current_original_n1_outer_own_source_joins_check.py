"""Actual reserved radial calculus and pending own seam identities."""
import json
from pathlib import Path
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins as current


class SourceAlgebra:
    def __init__(self,g,bindings):self.g,self.bindings,self.cache=g,bindings,{};self.flat_branches=set()
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.bindings[i]
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i];op=row['operation'];at=self.at
        if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op in ('current_original_source_parameter','original_source_parameter'):
            assert row['Z_independent']
            value=sy.Symbol('fixed_source_parameter_'+str(i),real=True)
        elif op=='bound_variable':value=sy.Symbol(row['name'],real=True)
        elif op=='sum':value=sy.Add(*(at(n) for n in row['arguments']))
        elif op=='product':value=sy.Mul(*(at(n) for n in row['arguments']))
        elif op=='negative':value=-at(row['argument'])
        elif op=='positive_quotient':value=at(row['numerator'])/at(row['denominator'])
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            value=at(row['expression']).subs(at(row['variable']),at(row['value']))
        elif op=='compact_raw_beta':
            x=at(row['argument']);value=sy.exp(-1/(1-x*x)) # interior formula; the flat branch is checked separately
        elif op=='original_lazy_flat_branch':
            assert row['original_beta_exact_flat_derivative'] and row['active_body_not_evaluated_on_flat_branch']
            assert row['eta']==1 and row['flat_value']==0
            predicate=self.g.nodes[row['flat_predicate']];assert predicate['operation']=='logical_or'
            expected={('<=',row['Delta'],self.g.constant(-1).node),('>=',row['Delta'],1)}
            comparisons=[]
            for j in predicate['arguments']:
                compare=self.g.nodes[j];assert compare['operation']=='exact_real_comparison'
                comparisons.append((compare['operator'],compare['left'],compare['right']))
            assert len(comparisons)==2 and set(comparisons)==expected
            self.flat_branches.add(i);value=at(row['active_body'])
        elif op=='analytic_unary':
            x=at(row['argument']);name=row['name']
            if name in ('exp','log','sqrt'):value=getattr(sy,name)(x)
            elif name=='compact_raw_beta_derivative':
                assert row['exact_original_derivative_not_its_global_cap']
                variable=sy.Symbol('raw_beta_argument',real=True)
                value=sy.diff(sy.exp(-1/(1-variable*variable)),variable).subs(variable,x)
            else:raise AssertionError('Unbound radial unary '+name)
        else:raise AssertionError('Unbound actual radial source '+op)
        self.cache[i]=value;return value


def bind(bindings,jet,body,Z):
    for n,ref in enumerate(current.rows(jet)):
        value=sy.diff(body,Z,n)
        if ref.node==0:assert value==0
        elif ref.node==1:assert value==1
        else:bindings[ref.node]=value


def radial_source_proof(owner):
    owner.assert_graph();g=owner.g
    Z=sy.Symbol('Z',real=True);Y=sy.Symbol(g.nodes[owner.reserved_phase.node]['name'],real=True)
    de=sy.Symbol('delta',positive=True)
    bindings={g.symbol('Z').node:Z,owner.delta.node:de}
    for key,jet in owner.coefficients.items():bind(bindings,jet,sy.Function('selected_'+key)(Z),Z)
    # Every scalar normalization integral in the actual physical bump
    # source has fixed support and no reserved query phase/radius.
    roots=[ref.node for key in ('Uz1','Utheta1','F1') for ref in current.rows(owner.reserved_fields[key])]
    closure=set().union(*(owner.operand_closure(i) for i in roots));normalizers=[]
    coefficient_roots={ref.node for jet in owner.coefficients.values() for ref in current.rows(jet)}
    stack=list(roots);seen=set()
    while stack:
        i=stack.pop()
        if i in seen or i in coefficient_roots:continue
        seen.add(i);row=g.nodes[i]
        if row['operation']=='definite_integral':
            assert owner.reserved_phase.node not in owner.operand_closure(i)
            assert owner.reserved_radius.node not in owner.operand_closure(i)
            bindings[i]=sy.Symbol('actual_normalizer_'+str(i),positive=True);normalizers.append(i)
            continue
        stack.extend(current.source_operands(row))
    alg=SourceAlgebra(g,bindings);velocity=[]
    for key in ('Uz1','Utheta1','F1'):
        for n,(original,derived) in enumerate(zip(current.rows(owner.reserved_fields[key]),current.rows(owner.reserved_velocity_y[key]))):
            assert sy.simplify(alg.at(derived)-sy.diff(alg.at(original),Y))==0,(key,n)
            velocity.append(key+'_y_Z'+str(n))
    b=owner.inputs[current.intermediate.previous.NAME]
    M=sy.Function('actual_own_mass')(Y,Z);bind(bindings,owner.jet(b['cumulative']['m11']),M,Z)
    alg=SourceAlgebra(g,bindings);U=alg.at(owner.reserved_fields['Uz1'].value);R=alg.at(owner.reserved_radius)
    Q=(2*Z*U-(1+de)*Z*M/R-(1-Z*Z)*sy.diff(M,Z)/R)/(1-de*Z*Z)
    recovery=[]
    for key,expected in (('Q1_y',sy.diff(Q,Y)),('V1_y',sy.diff(R*Q,Y))):
        for n,ref in enumerate(current.rows(owner.reserved_fields[key])):
            want=sy.diff(expected,Z,n)
            substitutions={sy.diff(M,Y,Z,k):R*sy.diff(U,Z,k) for k in range(4)}
            want=want.subs(substitutions,simultaneous=True)
            assert sy.simplify(alg.at(ref)-want)==0,(key,n)
            recovery.append(key+'_Z'+str(n))
    return dict(passed=True,actual_reserved_velocity_radial_C3_rows=velocity,
        actual_positive_order_Qy_Vy_C2_rows=recovery,
        actual_coefficient_phase_independent_rows=len(owner.coefficient_phase_independence),
        fixed_actual_bump_normalization_integrals=sorted(set(normalizers)),
        actual_mass_FTC_M_y_equals_R_U_consumed=True,
        exact_bump_derivative_source_not_global_caps=True,
        actual_outside_support_lazy_flat_branches_checked=len(alg.flat_branches),
        beta_boundary_derivatives_zero_by_original_smooth_compact_extension=True,
        own_one_plus_delta_recovery=True,seam_identities_and_radial_pressure_matching_not_claimed=True)


def save_construction_check(owner):
    """Save this exact partial proof without admitting any seam gate."""
    raw=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw['outer_own_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    proof=radial_source_proof(owner)
    receipt=dict(construction_only=True,source_radial_calculus_passed=proof['passed'],
        full_seam_acceptance_passed=False,source_family=owner.family,
        checks={'actual_reserved_radial_source':proof},**dict.fromkeys(current.GATES+current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        remaining_actual_seams=['2Rc','Rp','Rv'],radial_pressure_FTC_joins_not_admitted=True)
    owner.assert_graph()
    (current.HERE/current.CONSTRUCTION_RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('ACTUAL_RESERVED_RADIAL_SOURCE_CONSTRUCTION_CHECKED_SEAMS_OPEN',flush=True)
    return receipt
