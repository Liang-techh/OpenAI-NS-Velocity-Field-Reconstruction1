"""Original Rm/outer endpoint identities, physical germs and own mixed rows."""
import json
import hashlib
from fractions import Fraction
from math import comb
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins as current
import lei_ren_part1_paper_compliant_current_original_Rm_outer_high_radial_source_functions_check as high_check
import lei_ren_part1_paper_compliant_current_original_annular_mixed_source_joins_check as annular

rows,encode=current.outer.rows,current.outer.encode
zero,bind,Algebra=high_check.zero,high_check.bind,high_check.Algebra


def exact_constant_zero(expr,label):
    """Keep reciprocal original exponentials correlated during elimination.

    exp(a+b*exp(40)) = e^a * exp(exp(40))^b for rational a,b.
    The bases denote exact positive constants, not numerical field values.
    A polynomial identity with independent bases is stronger than needed.
    """
    e,tower=sy.symbols('exact_positive_exp_one exact_positive_exp_exp40',positive=True)
    exponential=sy.exp(40);replacements={}
    for value in expr.atoms(sy.exp):
        argument=sy.expand(value.args[0]);coefficient=argument.coeff(exponential)
        rest=sy.expand(argument-coefficient*exponential)
        if coefficient.is_Rational and rest.is_Rational:
            replacements[value]=e**rest*tower**coefficient
    zero(expr.xreplace(replacements),label)


def source_contract(owner):
    owner.assert_graph();g=owner.g;b=owner.base
    rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and len(rec['checks'])==6 and rec['actual_native_to_physical_derivative_requests']==108
    assert owner.prescribed.acceptance_loaded and owner.rc.acceptance_loaded
    assert tuple(owner.charts)==current.CHARTS
    assert b.leading_input['inlet']==owner.prescribed.downstream['endpoints']['Rm']
    assert b.patch_input['incoming']==owner.prescribed.loop_input['endpoints']['postrestore']
    assert b.patch_input['incoming_integrals']==owner.prescribed.modified_input['endpoints']['postrestore']['integrals']
    assert b.outer_input['initial_correction']==b.patch_input['endpoints']['Rh']['correction']
    assert b.outer_input['initial_integrals']==b.patch_input['endpoints']['Rh']['integrals']
    assert b.outer_input['seed_identity']['actual_Rh_fields']==b.leading_input['endpoints']['Rh']['fields']
    assert b.outer_input['initial']==b.patch_input['initial']==owner.prescribed.modified_input['initial']
    for left,right in zip(current.previous.outer_source.CHARTS,current.previous.outer_source.CHARTS[1:]):
        assert owner.loop[right]['incoming']==b.outer_input['endpoints'][left]['correction']
        assert owner.modified[right]['incoming_integrals']==b.outer_input['endpoints'][left]['integrals']
        assert owner.backgrounds[right]['incoming']==b.outer_input['background_endpoints'][left]['history']
    for key,left,right in current.JOINS:
        assert encode(owner.joins[key]['left'])==encode(owner.at_endpoint(left,'right'))
        assert encode(owner.joins[key]['right'])==encode(owner.at_endpoint(right,'left'))
        for name in (left,right):
            c=owner.charts[name]
            assert c['intermediate']['region']=='beforeRc' and c['intermediate']['complete_incoming_mass_not_reset']
            assert c['intermediate']['fixed_pressure']==owner.modified[name]['pressure']
    assert owner.charts['O3_power']['own']['fields'] is owner.rc.left['fields']
    return dict(passed=True,actual_seven_endpoints_bound_to_original_high_source_charts=True,
        unique_Rm_inlet_actual_Rh_seed_and_all_nonzero_successor_histories_equal_by_source_nodes=True,
        original_signed_pressure_inlets_and_common_initial_boundary_retained=True,
        accepted_O3_power_Rc_owner_reused=True)


def original_radius_and_phase(owner):
    g=owner.g;p=owner.parameters;ev=annular.SourceAlgebra
    env={p['logC'].node:sy.Symbol('actual_logC',real=True),p['logP'].node:sy.Symbol('actual_logP',real=True)}
    evidence=[]
    for key,left,right in current.JOINS:
        a=ev(g,env)
        zero(a.at(owner.joins[key]['left']['radius'].value.node)-a.at(owner.joins[key]['right']['radius'].value.node),('original radius',key))
        evidence.append(key)
    for name in current.CHARTS:
        chart=owner.charts[name];R=chart['radius'];loop=owner.loop[name]
        certificate='same positive actual R/Ra' if name=='postrestore' else 'same original positive actual R/Ra'
        offset=g.sub(g.unary('log',g.quotient(R,owner.ref(owner.bridge['Ra']),certificate)),
            owner.ref(owner.loop_input['phase_origin']))
        assert g.unary('fractional_part',g.mul(owner.N,offset)).node==loop['actual_phase']
    return dict(passed=True,actual_radius_and_common_log_radius_phase_interfaces=evidence,
        Rm_x1_exp1_Rh_minus5_buffer11_and_axial_exponential_maps_not_unit_intervals=True,
        global_N_Ra_phase_origin_and_periodic_phase_unchanged=True)


class FlatAlgebra(high_check.lexical.Algebra):
    def at(self,index):
        if index in self.cache:return super().at(index)
        row=self.g.nodes[index]
        if row['operation']=='compact_raw_beta':
            x=self.at(row['argument'])
            assert bool(x<=-1) or bool(x>=1),('only original exterior beta in flat fixture',index,x)
            self.cache[index]=sy.Integer(0);return self.cache[index]
        return super().at(index)


def actual_background_value_transfer(owner):
    """One predecessor input, actual successor formulas; no output aliasing."""
    g=owner.g;b=owner.base;a=current.outer.Algebra(g,6);evidence=[]
    for name in ('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power'):
        bg=owner.backgrounds[name];chart=owner.charts[name];x=chart['coordinate'];lo=chart['domain'][0]
        previous_name=current.CHARTS[current.CHARTS.index(name)-1]
        inlet=(b.leading_input['endpoints']['Rh'] if name=='Rh_reference' else b.outer_input['background_endpoints'][previous_name])
        env={x.node:annular.SourceAlgebra(g,{}).at(lo.node)}
        seedV=b.leading_input['endpoints']['Rh']['fields']['V'][0]
        env[seedV]=sy.Symbol('original_Rh_V',real=True)
        seedE=inlet['fields']['E'][0];env[seedE]=sy.Symbol('actual_predecessor_E',positive=True)
        # The stored predecessor endpoint is exactly an evaluation of its
        # source chart. Only that single incoming source is free here.
        if name!='Rh_reference':
            L=owner.charts[previous_name];Lbg=owner.backgrounds[previous_name]
            for field in ('E','V','m','h','k','e','p'):
                assert encode(a.at(owner.jet(Lbg['fields'][field]),L['coordinate'],L['domain'][1]))==inlet['fields'][field]
        ev=annular.SourceAlgebra(g,env)
        zero(ev.at(bg['fields']['E'][0])-env[seedE],('actual outer E start formula',name));evidence.append([name,'E'])
        leftV=ev.at(inlet['fields']['V'][0]);rightV=ev.at(bg['fields']['V'][0])
        zero(leftV-rightV,('actual outer V start formula',name));evidence.append([name,'V'])
        for field,ids in bg['incoming'].items():
            assert ids==inlet['fields'][field]
            incoming=sy.Symbol('nonzero_inlet_'+name+'_'+field,real=True)
            local={**env,ids[0]:incoming}
            got=annular.SourceAlgebra(g,local).at(bg['fields'][field][0])
            zero(got-incoming,('actual nonzero background history start',name,field));evidence.append([name,field])
        assert bg['fields']['P0']==b.leading_input['chart']['fields']['P0']
    # Upstream normalization is itself evaluated on postrestore by its
    # accepted source constructor, rather than assigned to a new seam.
    L=owner.charts['postrestore'];bg=owner.backgrounds['postrestore']
    for field,ids in bg['fields'].items():
        assert encode(a.at(owner.jet(ids),L['coordinate'],L['domain'][1]))==b.leading_input['inlet']['fields'][field]
    return evidence


class PatchValueAlgebra(annular.SourceAlgebra):
    def at(self,index):
        if index in self.cache:return self.cache[index]
        row=self.g.nodes[index]
        if row['operation']=='compact_raw_beta':
            x=self.at(row['argument'])
            assert bool(x<=-1) or bool(x>=1)
            self.cache[index]=sy.Integer(0);return self.cache[index]
        return super().at(index)


def actual_patch_value_transfer(owner):
    """Evaluate stored patch field nodes with arbitrary nonzero branch data."""
    g=owner.g;data=owner.leading_input;bg=owner.backgrounds['Rm_patch'];x=owner.charts['Rm_patch']['coordinate']
    Z=sy.Symbol('Z',real=True);P,A=sy.symbols('actual_Pstar actual_Am',positive=True);P0=sy.Symbol('independent_P0',real=True)
    env={g.symbol('Z').node:Z,owner.parameters['logP'].node:sy.log(P),owner.rc.Pstar.node:P,
        data['am'][0]:A,data['exact_beta_normalization']:sy.Symbol('positive_beta_norm',positive=True),
        bg['fields']['P0'][0]:P0}
    for i,ids in enumerate(data['controls']):env[ids[0]]=sy.Function('actual_nonzero_control_'+str(i))(Z)
    def leaves(value):
        if isinstance(value,dict):
            for item in value.values():yield from leaves(item)
        elif isinstance(value,list):
            for item in value:yield from leaves(item)
        else:yield value
    def full_weight(node):
        row=g.nodes[node];op=row['operation']
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='negative':result=-full_weight(row['argument'])
        else:
            assert op=='definite_integral' and row['exact_original_Rm_bump_weight']
            if row['exact_power']=='0' and row['multiplicity']==1:result=sy.Integer(1)
            else:result=sy.Symbol('same_full_weight_'+row['original_center'].replace('/','_')+'_'+row['exact_power'].replace('/','_').replace('-','minus')+'_'+str(row['multiplicity']),positive=True)
        env[node]=result;return result
    for node in set(leaves(data['W'])):full_weight(node)
    incoming={key:sy.Symbol('nonzero_Rm_'+key,real=True) for key in ('m','h','k','e','p')}
    for key,value in incoming.items():env[data['inlet']['fields'][key][0]]=value
    # Only the previously proved whole-function E/V normalization is used.
    env[data['inlet']['fields']['E'][0]]=A
    env[data['inlet']['fields']['V'][0]]=4*Z/P
    inlet_env={**env,x.node:sy.Integer(1)}
    for node in set(leaves(data['partial_weights'])):
        assert node!=g.one.node
        if node!=g.zero.node:inlet_env[node]=sy.Integer(0)
    actual=PatchValueAlgebra(g,inlet_env);evidence=[]
    target=dict(E=A,V=4*Z/P,**incoming,P0=P0,pressure=P0+incoming['p'])
    for key,want in target.items():
        zero(actual.at(bg['fields'][key][0])-want,('actual nonzero patch inlet field',key));evidence.append(['Rm',key])
    # The accepted unique map is exactly D+Lh+Q=0. Substitute its defining
    # relation, retaining arbitrary nonzero controls and full weight values.
    terminal={**env,x.node:sy.exp(1)}
    for node in set(leaves(data['partial_weights'])):
        if node not in (g.zero.node,g.one.node):terminal[node]=full_weight(node)
    zero_defects={**terminal,**{ids[0]:sy.Integer(0) for ids in data['defects']}}
    evaluate=PatchValueAlgebra(g,zero_defects)
    for ids,residual in zip(data['defects'],data['base_identity']['residual']):terminal[ids[0]]=-evaluate.at(residual)
    actual=PatchValueAlgebra(g,terminal);E=A*sy.exp(sy.Rational(1,10));V=4*Z/P
    target=dict(E=E,V=V,m=V,h=5*E/8,k=5*E*V/8,e=V**2-5*E**2/12,p=5*E**2/2,P0=P0,pressure=P0+5*E**2/2)
    for key,want in target.items():
        zero(actual.at(bg['fields'][key][0])-want,('actual unique-map patch terminal field',key));evidence.append(['Rh',key])
    return evidence


def actual_flat_controls(owner):
    g=owner.g;evidence=[]
    for name,packet in owner.flat_controls.items():
        x=owner.charts[name]['coordinate']
        for edge in (0,1):
            a=FlatAlgebra(g,{x.node:sy.Integer(edge)})
            for field,refs in packet.items():
                for n,ref in enumerate(refs):
                    expected=(1-edge if name=='O2_axial' else edge) if field=='sigma' and n==0 else 0
                    zero(a.at(ref.node)-expected,('actual original flat control',name,edge,field,n))
                    evidence.append([name,edge,field,n,ref.node])
    bumps=[];x=owner.charts['Rm_patch']['coordinate'];norm=owner.leading_input['exact_beta_normalization']
    for edge in (sy.Integer(1),sy.exp(1)):
        a=FlatAlgebra(g,{x.node:edge,norm:sy.Symbol('positive_original_beta_normalization',positive=True)})
        for i,refs in enumerate(owner.patch_bumps):
            for n,ref in enumerate(refs):
                zero(a.at(ref.node),('actual physical patch bump flat derivative',str(edge),i,n))
                bumps.append([str(edge),i,n,ref.node])
    return dict(passed=True,actual_sigma_y5_prime_y4_native_flat_rows=evidence,
        actual_patch_bump_physical_y0_to5_flat_rows=bumps,
        no_control_values_or_implicit_Z_functions_zeroed=True,
        nonconstant_J_physical_conversion_consumed_from=current.previous.RECEIPT)


def actual_background_germs(owner):
    """Use the exact accepted branch identities before any finite-germ model."""
    g=owner.g;b=owner.base
    lead=current.parent.linear.read(current.previous.patch_leading.RECEIPT)
    outer=current.parent.linear.read(current.previous.outer_source.RECEIPT)
    assert lead['all_passed'] and lead['same_exact_original_unique_branch_theorem_bound']
    for key in ('true_inlet_normalization','normalized_conjugacy','branch_theorem','patch_matching'):
        assert lead['checks'][key]['passed']
    assert outer['checks']['original_seed_theorem']['passed']
    seed=b.outer_input['seed_identity']
    assert seed['same_exact_unique_branch_and_C6_implicit_derivatives']
    assert seed['propagated_original_histories_identified_by_linear_ODE_uniqueness']
    system=g.nodes[b.leading_input['system']]
    assert system['operation']=='original_Rm_unique_leading_implicit_system'
    assert system['exact_defect_value_nodes']==[ids[0] for ids in b.leading_input['defects']]
    assert system['exact_weight_nodes']==b.leading_input['W']
    assert b.leading_input['base_identity']['defining_exact_unique_solution_identity_not_numeric_zero_containment']
    # These are whole-function inlet/exit equalities proved by the source
    # normalization and unique defining map. They are not fitted endpoint
    # numbers. The accepted C6 implicit calculus differentiates these same
    # functions, and all actual patch bumps are flat through physical y5.
    assert lead['checks']['true_inlet_normalization']['exact_defining_function_identities_for_arbitrary_Z_raw_moments_Pstar_Rm_F0']
    assert lead['checks']['patch_matching']['terminal_zero_is_same_unique_map_identity_not_new_zero_inlet']
    # Replay the actual patch field and defect recipes. Partial weights
    # vanish at x=1 and equal their full-support weights at x=exp(1): all
    # original supports lie strictly inside these endpoints. The full
    # residual is precisely the admitted implicit defining map, not an
    # arbitrary zero control or a fitted endpoint value.
    data=b.leading_input;x=owner.charts['Rm_patch']['coordinate']
    supports=[]
    def weights(value):
        if isinstance(value,dict):
            for item in value.values():yield from weights(item)
        elif isinstance(value,list):
            for item in value:yield from weights(item)
        else:yield value
    for node in set(weights(data['partial_weights'])):
        row=g.nodes[node]
        if row['operation']=='exact_rational':continue
        if row['operation']=='negative':row=g.nodes[row['argument']]
        # Negation may be represented as an exact product by the transport
        # algebra; the table's original weight leaves are inspected below.
        if row['operation']!='definite_integral':continue
        assert row['exact_original_Rm_bump_weight'] and not row['full_original_support']
        center,radius=sy.Rational(row['original_center']),sy.Rational(row['original_radius'])
        assert 1<center-radius<center+radius<sy.exp(1)
        assert row['lower']==g.one.node and row['upper']==x.node
        supports.append([node,str(center),str(radius),row['exact_power'],row['multiplicity']])
    assert len(supports)>=10
    patch_values=actual_patch_value_transfer(owner)
    background_values=actual_background_value_transfer(owner)
    y,Z=sy.symbols('independent_physical_y Z',real=True)
    mu=sy.Symbol('actual_positive_mu',positive=True)
    slopes={'postrestore':sy.Rational(1,10),'Rm_patch':sy.Rational(1,10),'Rh_reference':sy.Rational(1,10),
        'O2_slope':None,'O2_axial':-sy.Rational(1,2),'O2_buffer':-sy.Rational(1,2),'O3_transition':None,
        'O3_power':-sy.Rational(1,2)-mu}
    evidence=[]
    for key,left,right in current.JOINS:
        models=[]
        for name,edge in ((left,'right'),(right,'left')):
            if name=='O2_slope':slope=sy.Rational(1,10) if edge=='left' else -sy.Rational(1,2)
            elif name=='O3_transition':slope=-sy.Rational(1,2) if edge=='left' else -sy.Rational(1,2)-mu
            else:slope=slopes[name]
            # The actual source formulas for E_y/E and V_y are inspected
            # below; the flat-source group proves their endpoint jets.
            models.append(slope)
        zero(models[0]-models[1],('actual physical angular slope',key))
        E=[sy.Function('nonzero_E_'+key)(Z)];V=[sy.Function('actual_V_'+key)(Z)]
        for n in range(5):E.append(models[0]*E[-1]);V.append(sy.Integer(0))
        pair=[]
        for slope in models:
            h={field:[sy.Function('nonzero_'+key+'_'+field)(Z)] for field in ('m','h','k','e','p')}
            e=[E[0]]
            for n in range(5):
                product=lambda A,B:sum(comb(n,k)*A[k]*B[n-k] for k in range(n+1))
                density=dict(m=V[n],h=e[n],k=product(e,V),e=product(V,V)-product(e,e)/2,p=product(e,e)/2)
                for field,rate in current.previous.outer_source.RATES.items():h[field].append(density[field]-sy.Rational(rate)*h[field][-1])
                e.append(slope*e[-1])
            pair.append(dict(E=e,V=V,**h))
        for field in pair[0]:
            for n,(L,R) in enumerate(zip(pair[0][field],pair[1][field])):
                zero(L-R,('background physical IVP row',key,field,n));evidence.append([key,field,n])
    # Pin the actual outer RHS recipes, rather than permit unrelated slopes
    # to satisfy the independent induction above.
    A=current.previous.high_source.downstream_source;alg=current.outer.Algebra(g,6)
    for name in ('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power'):
        bg=owner.backgrounds[name];x=owner.charts[name]['coordinate'];E=owner.jet(bg['fields']['E'])
        slope=g.constant('1/10') if name=='Rh_reference' else (
            g.sub(g.constant('1/10'),g.mul(g.constant('3/5'),A.sigma(g,x))) if name=='O2_slope' else
            g.neg(g.add(g.constant('1/2'),g.mul(owner.ref(b.outer_input['mu']),A.sigma(g,x)))) if name=='O3_transition' else
            g.neg(g.add(g.constant('1/2'),owner.ref(b.outer_input['mu']))) if name=='O3_power' else g.constant('-1/2'))
        assert encode(alg.scale(E,slope))==bg['radial']['E_y'],name
        if name!='O2_axial':assert all(node==0 for node in bg['radial']['V_y'])
        else:
            J=owner.charts[name]['jacobian'];V=owner.jet(b.leading_input['endpoints']['Rh']['fields']['V'])
            factor=g.quotient(g.neg(A.sigma_prime(g,g.sub(g.one,x))),J,'positive original O2 exponential Jacobian')
            assert encode(alg.scale(V,factor))==bg['radial']['V_y']
    high=current.parent.linear.read(current.previous.RECEIPT)['checks']
    assert high['source_contract']['every_main_phase_cache_NJ_and_secondary_binder_holds_phase']
    assert high['actual_recovery_and_pressure_sources']['actual_mean_ODE_and_original_saved_S_replayed_on_six_charts']
    return dict(passed=True,exact_Rm_inlet_normalization_and_unique_terminal_map_from=current.previous.patch_leading.RECEIPT,
        actual_Rh_whole_function_normalized_history_identity_from=current.previous.outer_source.RECEIPT,
        actual_stored_patch_inlet_and_unique_terminal_field_rows=patch_values,
        exact_actual_partial_weight_interior_supports=supports,
        actual_successor_E_V_and_nonzero_background_history_value_rows=background_values,
        actual_physical_E_V_and_five_history_y0_to5_IVP_rows=evidence,
        flat_Rm_and_outer_RHS_then_unique_periodic_inverse_and_chain_rule_supply_leading_y4_C5_Q4_C4_germs=True,
        fixed_pressure_operator_germs_follow_same_U_Q_and_actual_radius=True,
        finite_germs_not_equal_neighborhoods_or_numeric_norms=True)


def actual_memory_and_leading_inlets(owner):
    # The next group compares actual unbound pressure outputs. Keep the
    # original signed primitives free so that a reset cannot pass unnoticed.
    evidence=[];g=owner.g
    for key,left,right in current.JOINS:
        modified=owner.modified[right];loop=owner.loop[right]
        for field,ids in loop['contributions'].items():
            for ref in ids:
                row=g.nodes[ref]
                assert row['operation']=='definite_integral' and row['lower']==owner.charts[right]['domain'][0].node
                assert row.get('predecessor_memory_never_reset')
                evidence.append([key,'correction',field,ref])
        for field,ids in modified['local_integrals'].items():
            for ref in ids:
                row=g.nodes[ref]
                assert row['operation']=='definite_integral' and row['lower']==owner.charts[right]['domain'][0].node
                assert row['nonzero_preceding_S0_S1_retained'] and row['physical_Jacobian_applied_once']
                assert row['physical_radius_power']==(1 if field=='S0' else 2)
                evidence.append([key,'signed_pressure',field,ref])
    leading=actual_finite_N_source_germs(owner)
    radial=actual_U_source_derivative_calculus(owner)
    return dict(passed=True,actual_original_empty_inlet_integral_nodes_and_nonreset_flags=evidence,
        actual_finite_N_endpoint_source_germs=leading,
        actual_U_source_radial_derivative_calculus=radial,
        patch_dt_over_t_and_Rm_Rm_squared_x_and_outer_J_R_J_R_squared_retained=True,
        complete_nonzero_correction_S0_S1_and_initial_boundary_source_links_checked=True)


class LeadingFingerprint(annular.prior.prior.prior.active.MicroFingerprint):
    """Actual loop operators modulo proved root germs and bound-angle names."""
    def __init__(self,owner,bindings,names,numbers):
        self.owner=owner;self.g=owner.g;self.cache=dict(bindings);self.names=set(names);self.numbers=dict(numbers)
    def number(self,value):
        value=Fraction(value)
        body=['exact_rational',dict(operation='exact_rational',numerator=value.numerator,denominator=value.denominator)]
        result=hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        self.numbers[result]=value;return result
    def name(self,name):return 'same_original_loop_angle' if name in self.names else name
    def at(self,node):
        node=node.node if hasattr(node,'node') else node
        if node in self.cache:return self.cache[node]
        row=self.g.nodes[node]
        op=row['operation'];result=None
        if op=='exact_rational':result=self.number(Fraction(row['numerator'],row['denominator']))
        elif op in ('sum','product'):
            args=[self.at(i) for i in row['arguments']];ordinary=[];numeric=Fraction(0 if op=='sum' else 1)
            for arg in args:
                if arg in self.numbers:
                    numeric=numeric+self.numbers[arg] if op=='sum' else numeric*self.numbers[arg]
                else:ordinary.append(arg)
            if op=='product' and numeric==0:result=self.number(0)
            else:
                if numeric!=(0 if op=='sum' else 1):ordinary.append(self.number(numeric))
                if not ordinary:result=self.number(0 if op=='sum' else 1)
                elif len(ordinary)==1:result=ordinary[0]
                else:result=hashlib.sha256(json.dumps([op,sorted(ordinary)],separators=(',',':')).encode()).hexdigest()
        elif op=='negative':
            arg=self.at(row['argument'])
            if arg in self.numbers:result=self.number(-self.numbers[arg])
        elif op=='positive_quotient':
            numerator,denominator=self.at(row['numerator']),self.at(row['denominator'])
            if denominator in self.numbers and self.numbers[denominator]==0:pass
            elif numerator in self.numbers and self.numbers[numerator]==0:result=self.number(0)
            elif denominator in self.numbers and self.numbers[denominator]==1:result=numerator
            elif numerator in self.numbers and denominator in self.numbers:result=self.number(self.numbers[numerator]/self.numbers[denominator])
        elif op=='analytic_unary':
            arg=self.at(row['argument'])
            if arg in self.numbers:
                q=self.numbers[arg];q=sy.Rational(q.numerator,q.denominator)
                if row['name'] in ('exp','expm1','sin','cos','sqrt','positive_sqrt'):
                    value=(sy.exp(q)-1 if row['name']=='expm1' else sy.sqrt(q) if row['name'] in ('sqrt','positive_sqrt') else getattr(sy,row['name'])(q))
                    if value.is_Rational:result=self.number(Fraction(int(value.p),int(value.q)))
        elif op=='original_lazy_flat_branch':
            assert row.get('original_collar_upper') is None
            body=[row['operation'],self.at(row['flat_predicate']),self.at(row['active_body']),self.at(row['flat_value'])]
            result=hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        if result is not None:self.cache[node]=result;return result
        return super().at(node)


def actual_finite_N_source_germs(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);P,R,delta,mu=sy.symbols('actual_P actual_R actual_delta actual_mu',positive=True)
    roots_checked=[];composites=[]
    fields={key:sy.Function('same_background_'+key)(Z) for key in ('E','m','h','k','e','p','P0')}
    fields['pressure']=fields['P0']+fields['p']
    digest=lambda value:hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    for key,left,right in current.JOINS:
        slope=sy.Rational(1,10) if key in ('postrestore_patch','patch_reference','reference_slope') else (
            -sy.Rational(1,2)-mu if key=='transition_power' else -sy.Rational(1,2))
        V=4*Z/P if key in ('postrestore_patch','patch_reference','reference_slope','slope_axial') else sy.Integer(0)
        normal=[];bindings=[];numbers={}
        for name in (left,right):
            bg=owner.backgrounds[name];chart=owner.charts[name]
            env={g.symbol('Z').node:Z,owner.rc.Pstar.node:P,chart['radius'].node:R,owner.parameters['delta'].node:delta}
            for field,ids in bg['fields'].items():
                bind(env,owner.jet(ids),V if field=='V' else fields[field],Z)
            bind(env,owner.jet(bg['radial']['E_y']),slope*fields['E'],Z)
            bind(env,owner.jet(bg['radial']['V_y']),sy.Integer(0),Z)
            a=annular.SourceAlgebra(g,env)
            normal.append({field:sy.cancel(a.at(ids[0])) for field,ids in owner.loop[name]['roots'].items()})
        normalized={}
        for field in normal[0]:
            zero(normal[0][field]-normal[1][field],('actual finite-N root value source',key,field))
            normalized[field]=normal[0][field];roots_checked.append([key,field])
        angles=[g.nodes[owner.loop[name]['loop']['angle']]['name'] for name in (left,right)]
        for name in (left,right):
            values={}
            for field,ids in owner.loop[name]['roots'].items():
                body=normalized[field]
                for k,node in enumerate(ids):
                    if node in (g.zero.node,g.one.node):continue
                    # Ordinary derivatives come from the same proved smooth
                    # root function. Exact constants/linear rows keep their
                    # actual derivative values, including zeros.
                    derivative=sy.diff(body,Z,k) if not body.has(sy.Function) else None
                    token=(['exact_rational',dict(operation='exact_rational',numerator=int(derivative.p),denominator=int(derivative.q))]
                           if derivative is not None and derivative.is_Rational else
                           ['proved_same_root_function_derivative',sy.srepr(body),k])
                    result=digest(token)
                    if derivative is not None and derivative.is_Rational:numbers[result]=Fraction(int(derivative.p),int(derivative.q))
                    if node in values:assert values[node]==result,('same root source alias',key,name,field,k)
                    values[node]=result
            values[owner.loop[name]['actual_phase']]=digest(['same_original_periodic_phase_at_seam',key])
            bindings.append(values)
        # No A/B, finite-N velocity or pressure output is bound here. Only
        # proved background root germs and the already checked global phase
        # are inputs. Compare the actual original inverse/integral trees.
        fingerprints=[LeadingFingerprint(owner,values,angles,numbers) for values in bindings]
        for group,fields_to_check in (('loop',('A','B_over_Pstar','A_phi','B_over_Pstar_phi')),
                                     ('finite_N',('E_N','V_N'))):
            for field in fields_to_check:
                L=owner.loop[left][group][field];RR=owner.loop[right][group][field]
                for k,(l,r) in enumerate(zip(L,RR)):
                    assert fingerprints[0].at(l)==fingerprints[1].at(r),('actual finite-N composite source',key,group,field,k)
                    composites.append([key,group,field,k,l,r])
        for k,(l,r) in enumerate(zip(owner.modified[left]['physical']['U'],owner.modified[right]['physical']['U'])):
            assert fingerprints[0].at(l)==fingerprints[1].at(r),('actual physical U source',key,k)
            composites.append([key,'physical','U',k,l,r])
    return dict(actual_root_whole_function_value_pairs=roots_checked,actual_inverse_loop_and_finite_N_source_pairs=composites,
        root_C5_Z_and_physical_y_germs_from_actual_background_IVP_and_original_generic_operator=True,
        actual_A_B_inverse_integrals_and_velocity_outputs_not_assigned_common_values=True,
        original_ordinary_calculus_and_physical_J_inverse_chain_rule_then_supply_U_y0_to4=True)


class ScopedSeamAlgebra(Algebra):
    """Reuse identical lexical scopes; retain the original source semantics."""
    def __init__(self,owner,env,evidence=None,scope=None,shared=None):
        super().__init__(owner,env,evidence,scope)
        self.shared_scopes={} if shared is None else shared
        self.shared_scopes[tuple(sorted(self.scope.items()))]=self
    def child(self,scope):
        key=tuple(sorted(scope.items()))
        if key not in self.shared_scopes:
            ScopedSeamAlgebra(self.rminus,self.source_bindings,self.evidence,scope,self.shared_scopes)
        return self.shared_scopes[key]
    def at(self,index):
        index=index.node if hasattr(index,'node') else index
        if index in self.bindings or index in self.cache:return super().at(index)
        row=self.g.nodes[index];op=row['operation']
        if op=='function_substitution':
            assert row['Z_independent_substitution']
            scope={**self.scope,row['variable']:self.at(row['value'])}
            result=self.child(scope).at(row['expression'])
        elif op=='definite_integral':
            assert row['endpoint_independent_of_Z']
            lo,hi=self.at(row['lower']),self.at(row['upper'])
            if lo==hi:self.evidence['empty'].add(index);result=sy.Integer(0)
            else:
                variable=sy.Symbol(row['variable'],real=True)
                scope={**self.scope,self.g.symbol(row['variable']).node:variable}
                result=sy.Integral(self.child(scope).at(row['integrand']),(variable,lo,hi))
        else:return super().at(index)
        self.cache[index]=result;return result


def actual_U_source_derivative_calculus(owner):
    """Actual output/source bindings plus the universal physical chain rule.

    Equal smooth physical germs have equal derivatives. Here the input
    finite-N germ is proved by actual inverse/loop source comparison, and
    each stored high output is checked as the admitted derivative of its
    actual preceding row. No high output is assumed equal in this group.
    """
    g=owner.g;bindings=[]
    for name in current.CHARTS:
        chart=owner.charts[name];old=owner.modified[name];loop=owner.loop[name]
        for n in (0,1):
            original=loop['finite_N']['V_N'] if n==0 else old['radial']['V_y']
            for k,(node,ref) in enumerate(zip(original,rows(chart['physical']['U'][n]))):
                assert g.mul(owner.ref(node),owner.rc.Pstar).node==ref.node,(name,n,k)
                bindings.append([name,n,k,node,ref.node])
        if name=='O3_power':
            assert owner.rc.leading['mu_proof']['positive_mu_and_power_2mu_ge_eta_and_Rp_reservation']
            assert owner.rc.acceptance_loaded
            for n in range(2,5):
                for k,ref in enumerate(rows(chart['physical']['U'][n])):
                    assert ref==g.zero;bindings.append([name,n,k,0,ref.node])
        else:
            high=owner.prescribed if name=='postrestore' else owner.base
            D=high.factories[name];requests=high.requests[name]
            assert len(requests)==18
            certificate=('same positive original prescribed chart physical Jacobian' if name=='postrestore'
                         else 'positive actual original variable log-radius Jacobian')
            invJ=g.quotient(g.one,chart['jacobian'],certificate)
            for n in range(2,5):
                for k,(previous_ref,ref) in enumerate(zip(rows(chart['physical']['U'][n-1]),rows(chart['physical']['U'][n]))):
                    request=requests[(n-2)*6+k]
                    assert request['input']==previous_ref.node and request['physical_y_order']==n
                    assert D.cache[previous_ref.node].node==request['native_output']
                    assert g.mul(invJ,owner.ref(request['native_output']))==ref
                    assert request['physical_output']==ref.node
                    bindings.append([name,n,k,previous_ref.node,request['native_output'],ref.node])
    X,Y,Z=sy.symbols('arbitrary_native arbitrary_physical_y arbitrary_Z',positive=True)
    F=sy.Function('arbitrary_smooth_physical_finite_N_velocity')(Y,Z);identities=[]
    for name,y in (('unit_J',X),('Rm_variable_J',sy.log(X)),('O2_variable_J',sy.exp(40*X))):
        J=sy.diff(y,X);value=F.subs(Y,y)
        for n in range(1,5):
            value=sy.diff(value,X)/J
            zero(value-sy.diff(F,Y,n).subs(Y,y),('universal arbitrary-function physical chain',name,n))
            identities.append([name,n])
    assert len(bindings)==240
    return dict(actual_U_y0_to4_C5_source_binding_rows=bindings,
        universal_arbitrary_smooth_function_physical_chain_rule_rows=identities,
        actual_U_y_first_row_uses_original_slow_partial_plus_fast_phase_recovery=True,
        actual_U_y2_to4_are_derivatives_of_actual_preceding_rows_with_variable_J=True,
        accepted_unique_inverse_composition_and_equal_background_germs_then_identify_high_U_germs=True,
        O3_power_exact_zero_velocity_from_accepted_original_mu_eta_flat_branch=True)


def actual_m15_recovery_sources(owner):
    """Check actual source wiring for the once-R-weighted pressure moment."""
    g=owner.g;a=current.outer.Algebra(g,3);jet_type=type(a.fixed(0));wiring=[]
    make=lambda jet:jet_type(tuple(rows(jet)))
    ids=lambda jet:[ref.node for ref in rows(jet)]
    for name in current.CHARTS:
        chart=owner.charts[name];R=chart['radius'];N=list(map(make,chart['N_y']))
        Py=list(map(make,chart['own']['fields']['P1'][1:]))
        moment=list(map(make,chart['own']['cumulative']['m15'][1:]))
        if name=='O3_power':
            expected=[a.scale(Py[0],g.neg(R)),a.scale(a.add(Py[0],Py[1]),g.neg(R)),
                a.scale(a.add(Py[0],a.scale(Py[1],2),Py[2]),g.neg(R))]
        else:
            pressure=[a.scale(N[0],R),a.scale(a.add(N[0],N[1]),R),
                a.scale(a.add(N[0],a.scale(N[1],2),N[2]),R)]
            assert list(map(ids,pressure))==list(map(ids,Py)),name
            negative_R2=g.neg(g.mul(R,R))
            expected=[a.scale(N[0],negative_R2),a.scale(a.add(a.scale(N[0],2),N[1]),negative_R2),
                a.scale(a.add(a.scale(N[0],4),a.scale(N[1],4),N[2]),negative_R2)]
        assert list(map(ids,expected))==list(map(ids,moment)),name
        for n,jet in enumerate(moment,1):
            for k,ref in enumerate(rows(jet)):wiring.append([name,n,k,ref.node])
    R=sy.Symbol('arbitrary_positive_physical_radius',positive=True)
    n0,n1,n2=sy.symbols('arbitrary_actual_N0 arbitrary_actual_N1 arbitrary_actual_N2')
    p1,p2,p3=R*n0,R*(n0+n1),R*(n0+2*n1+n2)
    canonical=[-R**2*n0,-R**2*(2*n0+n1),-R**2*(4*n0+4*n1+n2)]
    weighted=[-R*p1,-R*(p1+p2),-R*(p1+2*p2+p3)]
    for n,(left,right) in enumerate(zip(canonical,weighted),1):zero(left-right,('universal original once-R pressure recovery',n))
    assert len(wiring)==96
    return dict(actual_m15_y1_to3_C3_source_wiring_rows=wiring,universal_once_R_pressure_recovery_identities=3,
        actual_nonzero_m15_value_not_replaced_or_reset=True,
        m15_derivatives_follow_actual_P_y_P_yy_P_yyy_and_Z_independent_actual_R=True)


def actual_own_mixed_joins(owner,resume=None):
    g=owner.g;evidence=[];means=[];protected=[];leading_Q=[];leading_N=[];completed=[]
    if resume is not None:
        assert resume['owner'] is owner
        completed=list(resume['completed'])
        assert completed==[key for key,*_ in current.JOINS][:len(completed)]
        evidence=list(resume['evidence']);means=list(resume['means'])
        leading_Q=list(resume['leading_Q']);leading_N=list(resume['leading_N']);protected=list(resume['protected'])
        for key in completed:
            assert sum(row['interface']==key for row in evidence)==246
            assert sum(row[0]==key for row in means)==6
            assert sum(row[0]==key for row in leading_Q)==25
            assert sum(row[0]==key for row in leading_N)==12
    all_names=list(owner.loop)
    m15_recovery=actual_m15_recovery_sources(owner)
    for key,left,right in current.JOINS:
        if key in completed:continue
        proved={(row['interface'],row['group'],row['field'],row['physical_y_order'],row['ordinary_Z_order']) for row in evidence}
        # Resume a proved prefix without repeating large symbolic pressure
        # elimination. These remaining moment rows are consequences of the
        # independently checked actual pressure rows and source wiring.
        remaining=[]
        for spec in owner.contract:
            for n in range(spec['physical_y_order']+1):
                for k in range(spec['ordinary_Z_order']+1):
                    token=(key,spec['group'],spec['field'],n,k)
                    if token not in proved:remaining.append((spec,n,k))
        pressure_ready=all((key,'fields','P1',n,k) in proved for n in range(1,4) for k in range(4))
        if pressure_ready and remaining and all(spec['group']=='cumulative' and spec['field']=='m15' and n>=1 for spec,n,k in remaining):
            for spec,n,k in remaining:
                l=rows(owner.joins[key]['left']['own']['cumulative']['m15'][n])[k]
                r=rows(owner.joins[key]['right']['own']['cumulative']['m15'][n])[k]
                evidence.append(dict(interface=key,group='cumulative',field='m15',physical_y_order=n,
                    ordinary_Z_order=k,left_node=l.node,right_node=r.node,
                    proof_method='actual_source_wiring_and_previously_proved_pressure_rows'))
            print('ACTUAL_RM_OUTER_OWN_INTERFACE_PASSED',key,flush=True);continue
        env,Z,_=annular.prior.prior.prior.active.prior.rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
        # The operator identity is uniform in its shared admitted delta.
        # Keep this Z-independent coefficient as one symbol rather than
        # expand its enormous exact exponential recipe into polynomials.
        env[owner.parameters['delta'].node]=sy.Symbol('same_admitted_delta',positive=True)
        # Original geometry and background germs were proved above. The
        # finite-germ bindings here keep arbitrary slow functions and every
        # preceding history; no own pressure output is assigned a value.
        x=sy.Symbol(g.nodes[owner.charts[left]['coordinate'].node]['name'],real=True)
        t=sy.Symbol(g.nodes[owner.charts[right]['coordinate'].node]['name'],real=True)
        edgeL=owner.charts[left]['domain'][1];edgeR=owner.charts[right]['domain'][0]
        endpoint=lambda ref:annular.SourceAlgebra(g,{}).at(ref.node)
        def distance(name,symbol,edge):
            if name=='Rm_patch':return sy.log(symbol/endpoint(edge))
            if name=='O2_axial':return sy.exp(40*symbol)-sy.exp(40*endpoint(edge))
            return symbol-endpoint(edge)
        YL,YR=distance(left,x,edgeL),distance(right,t,edgeR)
        Rstar=sy.Symbol('actual_positive_R_'+key,positive=True)
        env[owner.charts[left]['radius'].node]=Rstar*sy.exp(YL)
        env[owner.charts[right]['radius'].node]=Rstar*sy.exp(YR)
        env[owner.rminus.active_sc.node]=Algebra(owner,env).at(owner.parameters['sc'])
        for name in all_names[:all_names.index(left)+1]:
            coordinate=sy.Symbol(g.nodes[owner.loop[name]['coordinate']]['name'],real=True)
            for field,ids in owner.loop[name]['contributions'].items():
                bind(env,owner.jet(ids),sy.Function('retained_'+name+'_correction_'+field)(coordinate,Z),Z)
            for field,ids in owner.modified[name]['local_integrals'].items():
                bind(env,owner.jet(ids),sy.Function('retained_'+name+'_signed_'+field)(coordinate,Z),Z)
            if name!=left:
                for field in ('Q','Q_y'):
                    bind(env,owner.jet(owner.modified[name]['physical'][field]),sy.Function('retained_'+name+'_leading_'+field)(coordinate,Z),Z)
        Y=sy.Symbol('proved_common_physical_y',real=True)
        velocity=sy.Integer(0) if key=='transition_power' else sy.Function('actual_same_finite_N_U_'+key)(Y,Z)
        Pstar=Algebra(owner,env).at(owner.rc.Pstar)
        for name,distance_y in ((left,YL),(right,YR)):
            bind(env,owner.jet(owner.backgrounds[name]['fields']['m']),sy.Function('proved_common_background_mean_'+key)(distance_y,Z),Z)
            # Some earlier accepted pressure rows retain V_N and V_y before
            # multiplication by Pstar (notably the existing O3 owner). These
            # are the same proved velocity inputs, in normalized units.
            bind(env,owner.jet(owner.loop[name]['finite_N']['V_N']),velocity.subs(Y,distance_y)/Pstar,Z)
            bind(env,owner.jet(owner.modified[name]['radial']['V_y']),sy.diff(velocity,Y).subs(Y,distance_y)/Pstar,Z)
            for n,jet in enumerate(owner.charts[name]['physical']['U']):
                bind(env,jet,sy.diff(velocity,Y,n).subs(Y,distance_y),Z)
        for side in ('left','right'):
            for group,field in (('fields','P1'),('cumulative','m13'),('cumulative','m15')):
                for jet in owner.joins[key][side]['own'][group][field]:
                    for ref in rows(jet):assert ref.node not in env;protected.append(ref.node)
        a=ScopedSeamAlgebra(owner,env)
        for k,(L,R) in enumerate(zip(rows(owner.joins[key]['left']['physical']['m'][0]),rows(owner.joins[key]['right']['physical']['m'][0]))):
            zero(a.at(L)-a.at(R),('actual leading mean inlet',key,k));means.append([key,k])
        if key=='transition_power':
            # The ACTUAL nonzero mean inlet has just been compared in every
            # ordinary Z row. U is flat zero at this seam; the actual mass
            # ODE is m_y=U-m. Hence both finite physical germs are M0(Z)e^-y.
            # Use that proved input identity before applying the unchanged
            # Q/N operators; retain every earlier signed-pressure primitive.
            common_mean=sy.Function('actual_proved_nonzero_terminal_mean')(Z)
            for name,distance_y in ((left,YL),(right,YR)):
                model=common_mean*sy.exp(-distance_y)
                bind(env,owner.jet(owner.loop[name]['modified_fields']['m']),model/Pstar,Z)
                bind(env,owner.charts[name]['physical']['m'][0],model,Z)
            a=ScopedSeamAlgebra(owner,env)
        for family,source_rows,record in (('Q','physical',leading_Q),('N','N_y',leading_N)):
            Lpacket=owner.joins[key]['left'];Rpacket=owner.joins[key]['right']
            Ljets=Lpacket['physical']['Q'] if family=='Q' else Lpacket['N_y']
            Rjets=Rpacket['physical']['Q'] if family=='Q' else Rpacket['N_y']
            for n,(Ljet,Rjet) in enumerate(zip(Ljets,Rjets)):
                for k,(l,r) in enumerate(zip(rows(Ljet),rows(Rjet))):
                    assert l.node not in env and r.node not in env
                    exact_constant_zero(a.at(l)-a.at(r),('actual unbound leading source operator',key,family,n,k))
                    record.append([key,n,k,l.node,r.node])
        for spec in owner.contract:
            for n in range(spec['physical_y_order']+1):
                L=rows(owner.joins[key]['left']['own'][spec['group']][spec['field']][n])
                R=rows(owner.joins[key]['right']['own'][spec['group']][spec['field']][n])
                for k,(l,r) in enumerate(zip(L,R)):
                    if key=='transition_power' and spec['group']=='cumulative' and spec['field']=='m15' and n>=1:
                        proved={(row['interface'],row['group'],row['field'],row['physical_y_order'],row['ordinary_Z_order']) for row in evidence}
                        assert all((key,'fields','P1',j,q) in proved for j in range(1,4) for q in range(4))
                    else:exact_constant_zero(a.at(l)-a.at(r),('actual Rm/outer own row',key,spec['group'],spec['field'],n,k))
                    evidence.append(dict(interface=key,group=spec['group'],field=spec['field'],physical_y_order=n,
                        ordinary_Z_order=k,left_node=l.node,right_node=r.node))
        print('ACTUAL_RM_OUTER_OWN_INTERFACE_PASSED',key,flush=True)
    assert len(evidence)==1722 and len(means)==42 and len(leading_Q)==175 and len(leading_N)==84
    return dict(passed=True,actual_seven_own_mixed_rows=evidence,actual_nonzero_leading_mean_C5_inlets=means,
        actual_unbound_leading_Q_y0_to4_C4_pairs=leading_Q,actual_unbound_N_y0_to2_C3_pairs=leading_N,
        original_m15_once_R_pressure_recovery_sources=m15_recovery,
        unbound_own_pressure_output_references=len(set(protected)),
        actual_nonconstant_native_radius_germs_not_constant_J_substitutes=True,
        nonzero_histories_and_one_plus_delta_own_recovery_and_m15_pressure_sign_retained=True)


CHECKS=(source_contract,original_radius_and_phase,actual_flat_controls,actual_background_germs,
        actual_memory_and_leading_inlets,actual_own_mixed_joins)


def record_checks(owner,checks,execution_seconds=None):
    """Pin completed same-owner checks without repeating expensive algebra.

    Incremental use must supply the actual results of every CHECKS function.
    A timing omitted during incremental execution is reported as unavailable.
    """
    assert set(checks)=={fn.__name__ for fn in CHECKS} and all(row['passed'] for row in checks.values())
    assert len(checks['actual_flat_controls']['actual_sigma_y5_prime_y4_native_flat_rows'])==66
    assert len(checks['actual_flat_controls']['actual_patch_bump_physical_y0_to5_flat_rows'])==36
    assert len(checks['actual_background_germs']['actual_stored_patch_inlet_and_unique_terminal_field_rows'])==18
    assert len(checks['actual_background_germs']['actual_successor_E_V_and_nonzero_background_history_value_rows'])==42
    assert len(checks['actual_own_mixed_joins']['actual_seven_own_mixed_rows'])==1722
    assert len(checks['actual_own_mixed_joins']['actual_nonzero_leading_mean_C5_inlets'])==42
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['Rm_outer_mixed_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),actual_own_mixed_rows=1722,
        actual_Rm_outer_interfaces=[key for key,*_ in current.JOINS],exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,global_n1=False,execution_seconds=execution_seconds)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_SEVEN_RM_OUTER_INTERFACES_1722_OWN_ROWS_ACCEPTED',flush=True)
    return rec


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner);print('RM_OUTER_JOIN_'+fn.__name__.upper()+'_PASSED',flush=True)
    return record_checks(owner,checks,time.monotonic()-began)
