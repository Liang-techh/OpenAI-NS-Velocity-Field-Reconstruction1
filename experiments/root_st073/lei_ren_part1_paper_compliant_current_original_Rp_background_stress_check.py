"""Independent Cartesian NS source assembly and signed 800-digit products."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from types import MethodType

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_background_stress as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError,KeyError):return True
    raise AssertionError('Invalid original background source admitted')


def graph_expression(graph,node,bindings,memo):
    node=node.node if hasattr(node,'node') else node
    if node in bindings:return bindings[node]
    if node in memo:return memo[node]
    data=graph.nodes[node];op=data['operation'];rec=lambda value:graph_expression(graph,value,bindings,memo)
    if op=='exact_rational':value=s.Rational(data['numerator'],data['denominator'])
    elif op=='sum':value=s.Add(*(rec(v) for v in data['arguments']))
    elif op=='product':value=s.Mul(*(rec(v) for v in data['arguments']))
    elif op=='negative':value=-rec(data['argument'])
    elif op=='positive_quotient':value=rec(data['numerator'])/rec(data['denominator'])
    elif op=='exact_operator_integer_power':value=rec(data['base'])**data['exponent']
    elif op=='analytic_unary' and data['name'] in ('sin','cos','exp','log'):value=getattr(s,data['name'])(rec(data['argument']))
    else:raise ValueError('Original exact physical operator required')
    memo[node]=value;return value


def raw_physical_expression(raw,graph,bindings,lam,memo):
    value=0
    for term in raw.terms:
        row=term.source_row;k,n=row.derivative
        primitive=s.diff(current.F[row.name],current.Y,k,current.Z,n)
        coefficient=graph_expression(graph,term.operator_function,bindings,memo)
        beta=graph_expression(graph,term.lambda_exponent,bindings,memo)
        radial=s.Rational(term.radial_power.numerator,term.radial_power.denominator)
        value+=coefficient*current.R**radial*lam**beta*primitive
    return value


def actual_cartesian_NS(owner,field,coordinates):
    record=owner.differential._require(field);lam=s.Symbol('lambda',positive=True)
    bindings={coordinates['Z'].node:current.Z,coordinates['cosine'].node:current.CS,
        coordinates['sine'].node:current.SN,owner.before.delta.node:current.DELTA,owner.before.mu.node:current.MU}
    memo={};raw=lambda component,index:raw_physical_expression(record['entries'][component,index]['raw'],owner.graph,bindings,lam,memo)
    components=('ux','uy','uz');axes=((1,0,0),(0,1,0),(0,0,1))
    velocities=[raw(component,(0,0,0)) for component in components]
    recovered=(2*current.Z*current.R*current.F['Uz']-(1-current.DELTA)*current.Z*current.F['Mz']
        -(1-current.Z**2)*s.diff(current.F['Mz'],current.Z))/((1-current.DELTA*current.Z**2)*s.sqrt(2*current.R))
    angle=s.Symbol('angle',real=True);proofs={}
    for i,(component,name) in enumerate(zip(components,('x','y','z'))):
        ns=raw(component,('t',))+sum(velocities[j]*raw(component,axes[j]) for j in range(3))
        ns+=raw('p',axes[i])-sum(raw(component,tuple(2*v for v in axis)) for axis in axes)
        proposed=sum(lam**s.sympify(beta)*expr for beta,expr in owner.operators['Cartesian_momentum_residual'][name])
        difference=s.expand_power_exp(ns-proposed)
        difference=current.reduce_source(difference.subs(current.F['Ur'],recovered).doit())
        difference=s.expand(s.powsimp(difference,force=True)).subs({current.CS:s.cos(angle),current.SN:s.sin(angle)})
        if s.trigsimp(s.cancel(difference))!=0:raise AssertionError('Actual Cartesian NS source operator differs: '+name)
        proofs[name]=True
    return proofs


def independent_tensor_completion():
    r,z,phi=s.symbols('r z phi',real=True);tt,tz=s.Function('Ttheta')(r,z),s.Function('Tz')(r,z)
    c,sn=s.cos(phi),s.sin(phi)
    basis=s.Matrix([[c,-sn,0],[sn,c,0],[0,0,1]])
    tensor=basis*s.Matrix([[0,tt,tz],[tt,r*s.diff(tz,z),0],[tz,0,0]])*basis.T
    derivatives=(lambda v:c*s.diff(v,r)-sn/r*s.diff(v,phi),lambda v:sn*s.diff(v,r)+c/r*s.diff(v,phi),lambda v:s.diff(v,z))
    div=s.Matrix([sum(derivatives[j](tensor[i,j]) for j in range(3)) for i in range(3)])
    expected=basis*s.Matrix([0,s.diff(tt,r)+2*tt/r,s.diff(tz,r)+tz/r])
    assert all(s.trigsimp(s.expand(v))==0 for v in div-expected)
    return dict(symmetric_Cartesian_completed_tensor=True,radial_divergence_zero=True,
        theta_radial_divergence_multiplier_2=True,axial_radial_divergence_multiplier_1=True)


def numeric(expr,c,bindings):
    if expr in bindings:return bindings[expr]
    if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
    if expr.is_Add:return sum((numeric(v,c,bindings) for v in expr.args),c.mpf(0))
    if expr.is_Mul:
        value=c.mpf(1)
        for part in expr.args:value*=numeric(part,c,bindings)
        return value
    if expr.is_Pow and expr.args[1].is_Integer:return numeric(expr.args[0],c,bindings)**int(expr.args[1])
    if expr.is_Pow and expr.args[0].is_Rational and expr.args[0].is_positive and expr.args[1].is_Rational:
        return c.exp(numeric(expr.args[1],c,bindings)*c.ln(numeric(expr.args[0],c,bindings)))
    raise ValueError('Independent directed coefficient required')


@source_precision
def run(owner,values,fields):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record and raw['original_operator_definitions']==owner.definitions
    assert raw['actual_radial_recovery_source_binding']==owner.radial_law
    assert all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name]==digest for name,digest in owner.before.hashes.items())
    tensor=independent_tensor_completion();observations={};products_checked=bounds_checked=0
    c=MPIntervalContext();c.dps=800
    for name,value in values.items():
        report=owner.report(value);field=fields[name];record=owner.differential._require(field)
        _,_,request=owner.product.amplitude.before._validate_delivery(record['delivery'])
        actual_NS=actual_cartesian_NS(owner,field,request.forward_coordinates)
        assert all(actual_NS.values())
        source=owner.product.before.source(record['delivery'])
        owner._validate_source(source)
        refined=owner.product.before.refined.source(record['delivery'])
        assert source['original_factorized_values']['P0'] is refined['original_factorized_values']['P0']
        assert source['original_factorized_values']['Mp'] is refined['original_factorized_values']['Mp']
        assert source['original_factorized_values']['pressure'] is not source['original_factorized_values']['Mp']
        assert 'source_equivalent_remaining_pressure_tail' in source
        issued_reader=owner.product.reader(record['delivery'])
        reader=current.correlated.CorrelatedGraphBounds(current.correlated.locator.CancelledGraphBounds(
            owner.graph,c,dict(issued_reader.bindings),issued_reader.allowed_exponentials),dict(issued_reader.aliases))
        angle=c.mpf(request.theta.numerator)/request.theta.denominator
        bindings={current.Z:c.mpf(request.Z.numerator)/request.Z.denominator,
            current.DELTA:c.mpf(issued_reader.bindings[owner.before.delta.node]),
            current.MU:c.mpf(issued_reader.bindings[owner.before.mu.node]),current.CS:c.cos(angle),current.SN:c.sin(angle)}
        signs={};zeros=unresolved=0
        for section,outputs in owner.operators.items():
            assert set(report[section])==set(outputs)
            signs[section]={}
            for component,result in report[section].items():
                groups={}
                for item in result['actual_source_product_ledger']:
                    coefficient=numeric(s.sympify(item['exact_operator_expression']),c,bindings)
                    units={'R':Fraction(item['radial_power']['numerator'],item['radial_power']['denominator'])}
                    product=c.mpf(1)
                    for atom in item['actual_shared_source_product']:
                        label=atom['source'];k,n=atom['ordinary_derivative'];power=atom['power']
                        source_row=source['log_radius_mixed_rows'][label]['y%d_Z%d'%(k,n)]
                        assert source_row.derivative==(k,n) and source_row.powers[-1]==0 and source_row.coefficients.order==0
                        assert atom['actual_source_units']==source_row.source_units
                        assert atom['actual_source_powers']==[current.signed.rational_record(q) for q in source_row.powers]
                        product*=c.mpf(source_row.coefficients[0])**power
                        for unit,q in zip(source_row.source_units,source_row.powers):units[unit]=units.get(unit,Fraction(0))+power*q
                    assert current.axial.differential.product.contains(item['signed_product_coefficient'],coefficient*product)
                    log=current.box.pulse.radius.FunctionRef(owner.graph,item['exact_product_log_scale'])
                    key=(tuple(sorted(reader.polynomial(log).items())),tuple(sorted(units.items())),item['lambda_exponent'])
                    group=groups.setdefault(key,dict(log=log,coefficient=c.mpf(0)))
                    group['coefficient']+=coefficient*product;products_checked+=1
                if result['exact_zero_enclosure']:
                    assert not groups and current.ends(result['ordinary_numeric_enclosure'])==(0,0)
                    zeros+=1;signs[section][component]=0;continue
                reference=current.box.pulse.radius.FunctionRef(owner.graph,result['exact_reference_log_scale_function'])
                total=c.mpf(0);failed=False
                for group in groups.values():
                    difference=owner.graph.sub(group['log'],reference);polynomial=reader.polynomial(difference)
                    if not polynomial:ratio=c.mpf(1)
                    else:
                        logratio=reader.at(difference);lo,hi=current.ends(logratio)
                        if hi>1000:failed=True;continue
                        if hi<=-1000:ratio=c.mpf([0,current.ends(c.exp(-1000))[1]])
                        elif lo<-1000:ratio=c.mpf([0,current.ends(c.exp(c.mpf(hi)))[1]])
                        else:ratio=c.exp(logratio)
                    total+=group['coefficient']*ratio
                if failed:
                    assert result['common_scale_coefficient_enclosure'] is None
                    signs[section][component]='ratio_unresolved';unresolved+=1
                else:
                    assert current.axial.differential.product.contains(result['common_scale_coefficient_enclosure'],total),(section,component)
                    signs[section][component]=result['signed_log_value']['sign'];bounds_checked+=1
                assert not result['ordinary_numeric_materialized']
                assert result['nonlinear_shared_source_products_collected_before_interval_bounds']
        assert zeros==5 and report['convention']=='R_B=-div(T_B)+E_B' and report['original_viscosity']==1
        assert report['completed_diagonal_r_partial_z_Tz_retained'] and report['actual_independent_P0_and_remaining_absolute_pressure_retained']
        assert report['regional_remainder_not_claimed_flat'] and all(not report[key] for key in current.OPEN)
        observations[name]=dict(actual_Cartesian_NS_source_operator_identities=actual_NS,output_count=30,
            exact_zero_outputs=zeros,unresolved_scale_ratio_outputs=unresolved,signs=signs,
            full_incoming_P0_and_pressure_retained=True,cone_flat_global_and_recursion_open=True)
    value=next(iter(values.values()));entry=owner._fields[id(value)]
    invalid=dict(copied_background_field=rejected(lambda:owner.report(copy.copy(value))),
        invalid_target=rejected(lambda:owner.evaluate(next(iter(fields.values())),'2')))
    pulse=owner.before.refined.pulse;saved_radial=pulse.radial
    try:
        pulse.radial=MethodType(lambda self,*args:None,pulse)
        invalid['substituted_original_radial_recovery']=rejected(lambda:owner.report(value))
    finally:pulse.radial=saved_radial
    invalid['missing_actual_radial_packet']=rejected(lambda:owner._validate_source({}))
    node=next(iter(entry[5]));saved=copy.deepcopy(owner.graph.nodes[node])
    try:
        owner.graph.nodes[node]['changed_current_background_operator']=True
        invalid['changed_defining_source_DAG']=rejected(lambda:owner.report(value))
    finally:owner.graph.nodes[node]=saved
    saved=owner.operators['cylindrical_stress']['theta_theta']
    try:
        owner.operators['cylindrical_stress']['theta_theta']=[]
        invalid['removed_completed_tensor_diagonal']=rejected(owner.assert_graph)
    finally:owner.operators['cylindrical_stress']['theta_theta']=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,original_operator_definitions=owner.definitions,
        actual_radial_recovery_source_binding=owner.radial_law,
        independent_completed_Cartesian_tensor_theorem=tensor,actual_background_stress_observations=observations,
        independently_checked_shared_source_products=products_checked,independently_checked_signed_output_bounds=bounds_checked,
        independent_800_digit_shared_source_product_bounds=True,actual_Cartesian_NS_assembled_from_original_derivative_rows=True,
        original_full_moments_independent_P0_and_absolute_pressure_retained=True,
        no_global_cone_flat_remainder_or_recursion_claim=True,invalid_inputs_rejected=invalid,input_hashes=hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_BACKGROUND_STRESS full original tensor, remainder and Cartesian NS source',flush=True)
    return receipt
