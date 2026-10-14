"""Actual reserved radial calculus and pending own seam identities."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import time
import sympy as sy
if __name__=='__main__':sys.setrecursionlimit(10000)
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
        else:raise AssertionError(('Unbound actual radial source',op,i))
        self.cache[i]=value;return value


def bind(bindings,jet,body,Z):
    for n,ref in enumerate(current.rows(jet)):
        value=sy.diff(body,Z,n)
        if ref.node==0:assert value==0
        elif ref.node==1:assert value==1
        else:bindings[ref.node]=value


class EndpointSourceAlgebra(SourceAlgebra):
    """Evaluate actual endpoint recipes before entering a flat bump body."""
    def __init__(self,g,bindings,evidence=None):
        super().__init__(g,bindings)
        self.evidence=evidence if evidence is not None else dict(flat_bumps=set(),empty_integrals=set())
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation']
        if op=='function_substitution':
            assert row['Z_independent_substitution']
            variable=self.g.nodes[row['variable']];assert variable['operation']=='bound_variable'
            symbol=sy.Symbol(variable['name'],real=True);value=self.at(row['value'])
            bindings={k:v.subs(symbol,value) if hasattr(v,'subs') else v for k,v in self.bindings.items()}
            bindings[row['variable']]=value
            result=EndpointSourceAlgebra(self.g,bindings,self.evidence).at(row['expression'])
        elif op=='definite_integral':
            assert row.get('endpoint_independent_of_Z'),('Unbound actual endpoint integral',i)
            assert sy.simplify(self.at(row['upper'])-self.at(row['lower']))==0,'Nonempty unbound actual primitive'
            self.evidence['empty_integrals'].add(i);result=sy.Integer(0)
        elif op=='compact_raw_beta' or (op=='analytic_unary' and row['name']=='compact_raw_beta_derivative'):
            x=sy.simplify(self.at(row['argument']))
            assert x.is_real is True
            if (x>=1) is sy.true or (x<=-1) is sy.true:
                self.evidence['flat_bumps'].add(i);result=sy.Integer(0)
            else:result=super().at(i)
        elif op=='analytic_unary' and row['name']=='exprel':
            x=self.at(row['argument']);result=sy.Integer(1) if x==0 else (sy.exp(x)-1)/x
        elif op=='original_lazy_flat_branch':
            # The branch predicate is checked by SourceAlgebra's radial proof.
            # At this endpoint inspect Delta before touching rational active terms.
            assert row['original_beta_exact_flat_derivative'] and row['active_body_not_evaluated_on_flat_branch']
            assert row['eta']==1 and row['flat_value']==0
            predicate=self.g.nodes[row['flat_predicate']];assert predicate['operation']=='logical_or'
            expected={('<=',row['Delta'],self.g.constant(-1).node),('>=',row['Delta'],1)}
            comparisons=[]
            for j in predicate['arguments']:
                compare=self.g.nodes[j];assert compare['operation']=='exact_real_comparison'
                comparisons.append((compare['operator'],compare['left'],compare['right']))
            assert len(comparisons)==2 and set(comparisons)==expected
            x=sy.simplify(self.at(row['Delta']))
            if (x>=1) is sy.true or (x<=-1) is sy.true:
                self.evidence['flat_bumps'].add(i);result=sy.Integer(0)
            else:result=super().at(i)
        else:return super().at(i)
        self.cache[i]=result;return result


def input_rows(bindings,ids,label):
    """Keep arbitrary actual incoming rows; structural constants stay fixed."""
    for n,i in enumerate(ids):
        if i not in (0,1):bindings.setdefault(i,sy.Symbol('actual_'+label+'_Z'+str(n),real=True))


def quiet_start_context(owner):
    """Source inputs for the compact exit and actual empty quiet primitive."""
    owner.assert_graph();g=owner.g;parent=current.previous
    c=owner.inlet.inputs[parent.compact.NAME];q=owner.inlet.inputs[parent.quiet.NAME]
    terminal=c['terminal']
    assert q['incoming']==terminal and q['incoming_integrals']==terminal['integrals']
    assert q['incoming_pressure']==terminal['pressure'] and q['initial']==c['initial']
    assert q['R2']==terminal['radius']
    xnode=g.nodes[terminal['modified_fields']['m'][0]]
    assert xnode['operation']=='function_substitution' and xnode['value']==g.constant(2).node
    x=xnode['variable'];X=sy.Symbol(g.nodes[x]['name'],real=True);Z=sy.Symbol('Z',real=True)
    bindings={g.symbol('Z').node:Z,x:X,q['s']:sy.Integer(0),c['mu']:sy.Symbol('mu',positive=True),
        c['N']:sy.Symbol('N',positive=True),c['delta']:sy.Symbol('delta',positive=True),
        c['Pstar']:sy.Symbol('Pstar',positive=True),c['Rc']:sy.Symbol('Rc',positive=True)}
    # Source inputs only: the selected coefficients, incoming amplitude and
    # actual complete mass. Do not bind the two boundaries or pressure outputs.
    for n,ids in enumerate(c['h']):bind(bindings,owner.jet(ids),sy.Function('compact_control_'+str(n))(Z),Z)
    bind(bindings,owner.jet(c['A']),sy.Function('actual_Rc_amplitude')(Z),Z)
    bind(bindings,owner.jet(c['modified_fields']['m']),sy.Function('actual_compact_mass')(X,Z),Z)
    for key,ids in c['integrals'].items():bind(bindings,owner.jet(ids),sy.Function('actual_compact_'+key)(X,Z),Z)
    for key,ids in c['initial'].items():input_rows(bindings,ids,'initial_'+key)
    normalizers=[];stack=[i for group in c['profiles']['radial_y_rows'] for key in ('V','delta_V') for i in group[key]];seen=set()
    while stack:
        i=stack.pop()
        if i in seen or i in bindings:continue
        seen.add(i);row=g.nodes[i]
        if row['operation']=='definite_integral':
            closure=owner.operand_closure(i)
            assert x not in closure and q['s'] not in closure and g.symbol('Z').node not in closure
            bindings[i]=sy.Symbol('actual_compact_normalizer_'+str(i),positive=True);normalizers.append(i)
        else:stack.extend(current.source_operands(row))
    return g,c,q,bindings,sorted(set(normalizers))


def quiet_start_source_proof(owner):
    """Actual empty primitive plus the compact/quiet finite boundary identity."""
    g,c,q,bindings,normalizers=quiet_start_context(owner);terminal=c['terminal']
    alg=EndpointSourceAlgebra(g,bindings);flat=[]
    for label,ids in (('V',terminal['modified_fields']['V']),('V_y',terminal['radial']['V_y']),
        ('V_yy',terminal['radial']['V_yy'])):
        for n,i in enumerate(ids):
            assert sy.simplify(alg.at(i))==0,(label,n)
            flat.append(label+'_Z'+str(n))
    empty=[]
    for key,ids in q['pressure_local'].items():
        for n,i in enumerate(ids):
            row=g.nodes[i]
            assert row['operation']=='definite_integral' and row['lower']==0 and row['upper']==q['s']
            assert row['physical_radius_power']==(1 if key=='S0' else 2)
            assert row['physical_Jacobian_applied_once'] and row['ordinary_slow_Z_derivative_order']==n
            assert alg.at(i)==0;empty.append(key+'_Z'+str(n))
    joined=[]
    for group in ('integrals','boundary','pressure'):
        for key,ids in q['inlet'][group].items():
            for n,(left,right) in enumerate(zip(ids,terminal[group][key])):
                assert sy.simplify(alg.at(left)-alg.at(right))==0,(group,key,n)
                joined.append(group+'_'+key+'_Z'+str(n))
    assert sy.simplify(alg.at(q['inlet']['radius'])-alg.at(terminal['radius']))==0
    increment=[]
    for key in ('m13','m15'):
        for n,(left,right) in enumerate(zip(q['inlet']['pressure'][key],q['incoming_pressure'][key])):
            assert sy.simplify(alg.at(left)-alg.at(right))==0,(key,n)
            increment.append(key+'_Z'+str(n))
    assert len(alg.evidence['empty_integrals'])==8 and alg.evidence['flat_bumps']
    return dict(passed=True,actual_compact_flat_velocity_radial_rows=flat,
        actual_empty_quiet_primitive_rows=empty,actual_source_boundary_integral_pressure_join_rows=joined,
        actual_quiet_pressure_increment_zero_rows=increment,
        actual_quiet_and_compact_radius_identity=True,
        nonzero_incoming_integrals_and_initial_finite_boundary_terms_retained=True,
        actual_pressure_outputs_and_boundary_outputs_not_bound=True,
        fixed_actual_compact_normalization_integrals=sorted(set(normalizers)),
        flat_source_nodes=sorted(alg.evidence['flat_bumps']),empty_source_integral_nodes=sorted(alg.evidence['empty_integrals']),
        own_full_tuple_and_radial_pressure_matching_not_yet_claimed=True)


def source_assignments(owner,stem,cls,function,expected):
    """Pin exact production assignments without replaying numerical owners."""
    name=current.PREFIX+stem+'.py';assert owner.hashes[name]==current.sha(name)
    tree=ast.parse((current.HERE/name).read_text(encoding='utf-8'))
    if cls:tree=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
    if function:tree=next(n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name==function)
    assignments={}
    def collect(target,value):
        if isinstance(target,(ast.Tuple,ast.List)):
            if not isinstance(value,(ast.Tuple,ast.List)):return
            assert len(target.elts)==len(value.elts)
            for t,v in zip(target.elts,value.elts):collect(t,v)
        else:assignments[ast.unparse(target)]=ast.unparse(value)
    for n in ast.walk(tree):
        if isinstance(n,ast.Assign):
            for target in n.targets:collect(target,n.value)
    checked={}
    for target,value in expected.items():
        canonical=ast.unparse(ast.parse(value,mode='eval').body)
        assert assignments.get(target)==canonical,(name,target,assignments.get(target),canonical)
        checked[target]=canonical
    return dict(module=name,module_sha256=current.sha(name),class_name=cls,function=function,
        AST_sha256=hashlib.sha256(ast.dump(tree).encode()).hexdigest(),actual_assignments=checked)


def selected_logC_source_projection(owner,old_node,true_node):
    """Resolve two parameter leaves to the same actual selected datum."""
    g=owner.g;old,new=g.nodes[old_node],g.nodes[true_node]
    assert old['operation']==new['operation']=='current_original_source_parameter'
    assert old['name']=='logC' and old['definition']=='same selected_logCstar and CurrentLongRadiusPhase exact singleton tuple'
    assert old['quantity_is_exact_original_parameter_not_range_endpoint'] and old['Z_independent']
    assert old['source_family']==new['source_family']==owner.family
    assert owner.hashes[old['current_source_receipt']]==old['current_source_receipt_sha256']==current.sha(old['current_source_receipt'])
    assert new['name']=='true_core_logC' and new['source_path']=='common.core.logC'
    assert new['exact_fixed_source_parameter'] and new['parameter_enclosure_not_selected_value'] and new['Z_independent']
    provider=new['source_provider'];assert provider['function']=='CurrentOriginalN1CoreFixedPointFunctions.attach_oracle'
    assert owner.hashes[provider['module']]==provider['module_sha256']==current.sha(provider['module'])
    receipt_name=current.PREFIX+'current_original_n1_core_fixed_point_functions_check.json'
    receipt=current.parent.linear.read(receipt_name)
    assert owner.hashes[receipt_name]==current.sha(receipt_name) and receipt['all_passed'] and receipt['source_family']==owner.family
    parameter=receipt['checks']['actual_oracle']['exact_parameter_and_native_P0_identity']
    scalar=parameter['actual_source_scalar_and_unique_anchor_identity']
    assert parameter['passed'] and scalar['logC_matches_original_selected_source'] and scalar['logC_leaf_uses_this_exact_source']
    specs=[
        ('current_original_whole_Z_long_Rm_finite_N','WholeZLongRmFiniteN','__init__',
            {'params':'json.loads((HERE/kernels.PARAMS).read_bytes())','self.logC':"read(c,params['selected_logCstar'])"}),
        ('current_original_reshape_terminal_kernels',None,None,{'PARAMS':"PREFIX+'physical_norm_family.json'"}),
        ('current_original_whole_Z_long_Rm_finite_N','WholeZLongRmFiniteN','owner',
            {'ref':"background.ActualReferenceRestoreFunctions(f,z,source.core.delta,terminal['actual_terminal_normalized_six_history_shapes'],fields['V'],P0,self.T,self.logC,self.provider)"}),
        ('current_original_reference_restore_functions','ActualReferenceRestoreFunctions','__init__',{'self.logC':'c.mpf(logC)'}),
        ('current_original_whole_Z_long_Rm_signed_functions','CurrentLongRadiusPhase','__init__',{'self.logC':'long_owner.reference.logC'}),
        ('current_core_common_fixed_point','CurrentCoreCommonFixedPoint','__init__',{'self.core':'self.source.core'}),
        ('core_physical_field','CompliantCorePhysicalField','__init__',
            {'norm':"self.records['physical_norm_family']",'self.logC':"read(norm,'selected_logCstar')"})]
    projections=[source_assignments(owner,*spec) for spec in specs]
    datum_name=current.PREFIX+'physical_norm_family.json';datum=current.parent.linear.read(datum_name)
    assert owner.hashes[datum_name]==current.sha(datum_name) and 'selected_logCstar' in datum
    return dict(passed=True,old_parameter_node=old_node,true_parameter_node=true_node,
        canonical_datum=datum_name,canonical_datum_sha256=current.sha(datum_name),
        canonical_datum_key='selected_logCstar',actual_source_assignments=projections,
        accepted_true_core_parameter_identity_consumed=True,
        same_selected_scalar_recipe_not_interval_overlap_or_new_common_radius=True)


def radius_source_proof(owner):
    """Actual old/current Rc, 2Rc, Rp, Rv recipes on one parameter source."""
    owner.assert_graph();g=owner.g;c=owner.inputs[current.intermediate.previous.previous.NAME]
    b=owner.inputs[current.intermediate.previous.NAME]
    row=g.nodes[c['repair_velocity']['physical_radius']]
    assert row['operation']=='analytic_unary' and row['name']=='exp'
    log_row=g.nodes[row['argument']];assert log_row['operation']=='sum' and len(log_row['arguments'])==2
    phase=c['repair_velocity']['log_R_over_Rc'];assert phase in log_row['arguments']
    old_logRc=next(i for i in log_row['arguments'] if i!=phase)
    old_leaves=[i for i in owner.operand_closure(old_logRc)
        if g.nodes[i]['operation']=='current_original_source_parameter' and g.nodes[i]['name']=='logC']
    assert len(old_leaves)==1
    true_logC=owner.inlet.seam.core.parameters['logC'].node
    projection=selected_logC_source_projection(owner,old_leaves[0],true_logC)
    mu_node=owner.inlet.geometry['mu'].node;mu_row=g.nodes[mu_node]
    assert mu_row['operation']=='analytic_unary' and mu_row['name']=='exp'
    log_mu=mu_row['argument'];plain=SourceAlgebra(g,{})
    LP=sy.exp(40)+11
    assert sy.simplify(plain.at(log_mu)+4*LP+sy.log(1000))==0
    C=sy.Symbol('actual_selected_logCstar',real=True);mu=sy.Symbol('mu',positive=True)
    alg=SourceAlgebra(g,{old_leaves[0]:C,true_logC:C,mu_node:mu,log_mu:sy.log(mu)})
    expected_logRc=sy.log(110)+10*C+11*LP+3
    def logzero(value):return sy.simplify(sy.expand_log(value,force=True))==0
    assert logzero(alg.at(old_logRc)-expected_logRc)
    current_Rc=owner.inlet.inputs[current.previous.outer.NAME]['endpoints']['O3_power']['radius']
    assert logzero(sy.log(alg.at(current_Rc))-expected_logRc)
    expected={'2Rc':expected_logRc+sy.log(2),'Rp':expected_logRc-2-60*sy.log(mu),
        'Rv':expected_logRc-2-60*sy.log(mu)+13/mu}
    checked=[]
    for seam,packet in owner.joins.items():
        for side in ('left','right'):
            node=packet[side]['radius'].node
            assert logzero(sy.log(alg.at(node))-expected[seam]),(seam,side,node)
            checked.append(dict(seam=seam,side=side,actual_radius_node=node))
    assert logzero(alg.at(b['tail_origins']['mass_log_radius'])-expected['2Rc'])
    assert logzero(alg.at(b['tail_origins']['pressure_log_radius'])-expected['Rv'])
    return dict(passed=True,selected_parameter_source_projection=projection,
        actual_old_logRc_node=old_logRc,actual_current_Rc_node=current_Rc,
        actual_six_seam_radius_sources=checked,
        original_mu_recipe='exp(-4*(exp(40)+11))/1000',
        actual_logRc='log(110)+10*logCstar+11*logPstar+3',
        actual_logRp='logRc-2-60*log(mu)',actual_logRv='logRp+13/mu',
        mass_and_pressure_tail_origins_matched=True,
        independent_radii_not_bound_to_a_manufactured_common_symbol=True)


def own_2Rc_source_join_proof(owner,radius_proof,quiet_proof):
    """Compare both saved own tuples from their actual common primitives."""
    assert radius_proof['passed'] and quiet_proof['passed']
    g,c,q,bindings,_=quiet_start_context(owner)
    projection=radius_proof['selected_parameter_source_projection']
    C=sy.Symbol('actual_selected_logCstar',real=True)
    bindings[projection['old_parameter_node']]=bindings[projection['true_parameter_node']]=C
    # Keep both full radius recipes. The only shared scalar substitution is
    # their independently proved canonical selected-datum projection.
    bindings.pop(c['Rc'])
    b=owner.inputs[current.intermediate.previous.NAME]
    Y=sy.Symbol(g.nodes[owner.reserved_phase.node]['name'],real=True);Z=sy.Symbol('Z',real=True)
    for key,jet in owner.coefficients.items():bind(bindings,jet,sy.Function('selected_'+key)(Z),Z)
    for key,ids in b['cumulative'].items():bind(bindings,owner.jet(ids),sy.Function('actual_own_'+key)(Y,Z),Z)
    bind(bindings,owner.jet(b['patch_bodies']['leading_Utheta0']),sy.Function('actual_leading_Utheta0')(Y,Z),Z)
    roots=[ref.node for key in ('Uz1','Utheta1','F1','Q1_y','V1_y') for ref in current.rows(owner.reserved_fields[key])]
    roots.extend(b['total_pressure_y'])
    normalizers=[];seen=set();stack=list(roots)
    while stack:
        i=stack.pop()
        if i in seen or i in bindings:continue
        seen.add(i);row=g.nodes[i]
        if row['operation']=='definite_integral':
            closure=owner.operand_closure(i)
            assert owner.reserved_phase.node not in closure and g.symbol('Z').node not in closure
            bindings[i]=sy.Symbol('actual_reserved_normalizer_'+str(i),positive=True);normalizers.append(i)
        else:stack.extend(current.source_operands(row))
    alg=EndpointSourceAlgebra(g,bindings);packet=owner.joins['2Rc'];checked=[]
    for group in ('fields','cumulative'):
        assert set(packet['left'][group])==set(packet['right'][group])
        for key,left in packet['left'][group].items():
            right=packet['right'][group][key]
            assert len(current.rows(left))==len(current.rows(right))
            for n,(l,r) in enumerate(zip(current.rows(left),current.rows(right))):
                try:assert sy.simplify(alg.at(l)-alg.at(r))==0,(group,key,n)
                except AssertionError as error:raise AssertionError((group,key,n,error.args)) from error
                checked.append(group+'_'+key+'_Z'+str(n))
    assert sy.simplify(alg.at(packet['left']['radius'])-alg.at(packet['right']['radius']))==0
    assert owner.charts['reserved']['physical_log_radius_jacobian']==g.one
    assert owner.charts['quiet']['physical_log_radius_jacobian']==g.one
    radial=[]
    for n,(left,right) in enumerate(zip(current.rows(owner.pressure_radial_traces['2Rc']['left']),
        current.rows(owner.pressure_radial_traces['2Rc']['right']))):
        assert sy.simplify(alg.at(left)-alg.at(right))==0,('actual_radial_pressure',n)
        radial.append('P1_y_Z'+str(n))
    return dict(passed=True,actual_own_2Rc_field_and_cumulative_trace_rows=checked,
        actual_2Rc_radial_pressure_trace_rows=radial,
        actual_both_radius_recipes_retained=True,
        both_actual_physical_log_radius_Jacobians_equal_one=True,
        arbitrary_nonzero_own_mass_and_pressure_primitives_retained=True,
        actual_flat_reserved_velocity_derivatives_and_compact_leading_velocity_consumed=True,
        quiet_zero_start_signed_pressure_identity_consumed=True,
        fixed_leading_N1p_and_own_one_plus_delta_not_confused=True,
        full_own_trace_outputs_not_bound_to_common_values=True,
        fixed_actual_reserved_normalization_integrals=sorted(set(normalizers)),
        actual_flat_source_nodes=sorted(alg.evidence['flat_bumps']),
        remaining_actual_seams=['Rp','Rv'],full_seam_acceptance_not_claimed=True)


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
            integral_closure=owner.operand_closure(i)
            assert owner.reserved_phase.node not in integral_closure
            assert owner.reserved_radius.node not in integral_closure and g.symbol('Z').node not in integral_closure
            bindings[i]=sy.Symbol('actual_normalizer_'+str(i),positive=True);normalizers.append(i)
            continue
        stack.extend(current.source_operands(row))
    alg=SourceAlgebra(g,bindings);velocity=[]
    for key in ('Uz1','Utheta1','F1'):
        for n,(original,derived) in enumerate(zip(current.rows(owner.reserved_fields[key]),current.rows(owner.reserved_velocity_y[key]))):
            assert sy.simplify(alg.at(derived)-sy.diff(alg.at(original),Y))==0,(key,n)
            velocity.append(key+'_y_Z'+str(n))
    b=owner.inputs[current.intermediate.previous.NAME]
    # Derive the mass FTC from the saved primitive itself before using it
    # in a differentiated recovery formula. The inlet is phase independent.
    mass_FTC=[]
    for n,(primitive,field_row) in enumerate(zip(b['patch_histories']['m11'],b['fields']['M1_y'])):
        row=g.nodes[primitive]
        assert row['operation']=='definite_integral' and row['lower']==0 and row['upper']==owner.reserved_phase.node
        assert row['endpoint_independent_of_Z'] and row['ordinary_slow_Z_derivative_order']==n
        assert row['actual_signed_repaired_n1_history'] and row['physical_measure']=='dR=R ds'
        assert owner.reserved_phase.node not in owner.operand_closure(b['incoming']['m11'][n])
        incoming=b['incoming']['m11'][n]
        prefix=sy.Integer(incoming) if incoming in (0,1) else sy.Symbol('actual_mass_inlet_'+str(n),real=True)
        history=sy.Symbol('actual_mass_primitive_'+str(n),real=True)
        cumulative=SourceAlgebra(g,{incoming:prefix,primitive:history})
        assert sy.simplify(cumulative.at(b['cumulative']['m11'][n])-prefix-history)==0,('actual_mass_inlet_plus_primitive',n)
        t=sy.Symbol(row['variable'],real=True)
        expected=alg.at(owner.reserved_radius)*alg.at(current.rows(owner.reserved_fields['Uz1'])[n])
        assert sy.simplify(alg.at(row['integrand']).subs(t,Y)-expected)==0,('mass_primitive_FTC',n)
        assert sy.simplify(alg.at(field_row)-expected)==0,('actual_mass_y_source_row',n)
        mass_FTC.append('M1_y_Z'+str(n))
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
        actual_mass_FTC_source_primitive_and_field_rows=mass_FTC,
        exact_bump_derivative_source_not_global_caps=True,
        actual_outside_support_lazy_flat_branches_checked=len(alg.flat_branches),
        beta_boundary_derivatives_zero_by_original_smooth_compact_extension=True,
        own_one_plus_delta_recovery=True,seam_identities_and_radial_pressure_matching_not_claimed=True)


def save_construction_check(owner):
    """Save this exact partial proof without admitting any seam gate."""
    began=time.monotonic()
    raw=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw['outer_own_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    proof=radial_source_proof(owner);radii=radius_source_proof(owner);quiet=quiet_start_source_proof(owner)
    seam=own_2Rc_source_join_proof(owner,radii,quiet)
    assert len(seam['actual_own_2Rc_field_and_cumulative_trace_rows'])==68
    assert len(seam['actual_2Rc_radial_pressure_trace_rows'])==4
    receipt=dict(construction_only=True,source_radial_calculus_passed=proof['passed'],
        full_seam_acceptance_passed=False,source_family=owner.family,
        checks={'actual_reserved_radial_source':proof,'actual_old_current_radius_sources':radii,
            'actual_quiet_start_source':quiet,'actual_own_2Rc_source_join':seam},
        actual_own_2Rc_source_join_passed=True,actual_2Rc_radial_pressure_C3_join_passed=True,
        actual_own_Rp_source_join_passed=False,actual_own_Rv_source_join_passed=False,
        **dict.fromkeys(current.GATES+current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        remaining_actual_seams=['Rp','Rv'],radial_pressure_FTC_joins_not_admitted=True,
        radial_pressure_all_seams_not_admitted=True,execution_seconds=time.monotonic()-began)
    owner.assert_graph()
    (current.HERE/current.CONSTRUCTION_RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('ACTUAL_2RC_SOURCE_JOIN_68_ROWS_AND_PRESSURE_Y_4_ROWS_CHECKED_RP_RV_OPEN',flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--construction',action='store_true',help='Restore source only and check the currently constructed seams.')
    args=parser.parse_args()
    if not args.construction:parser.error('--construction is required; full seam acceptance remains open')
    sys.setrecursionlimit(10000)
    save_construction_check(current.CurrentOriginalN1OuterOwnSourceJoins(require_checked=False))
