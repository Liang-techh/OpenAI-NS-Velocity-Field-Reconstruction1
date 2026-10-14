"""Exact physical change of variables before transported terminal cancellation."""
import copy
from fractions import Fraction
import json
from pathlib import Path
from types import SimpleNamespace
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities as current
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions_check as linear_check


class BetaIntegralAlgebra:
    """Interpret actual beta/exprel integral bodies, never numerical weights."""
    def __init__(self,g,bindings,ell):
        self.g,self.bindings,self.ell=g,bindings,ell
        self.beta=sy.Function('original_beta');self.J=sy.Function('original_beta_exponential_integral')
    def integrate(self,expression,t):
        if expression==0:return sy.Integer(0)
        # exprel has its continuous entire extension at t=0. Cancelling its
        # removable denominator below is an identity of the integrand.
        quotient=sy.cancel(expression/self.beta(t));out=sy.Integer(0)
        for term in sy.expand(quotient).as_ordered_terms():
            term=sy.powsimp(term,force=True)
            scalar,dependent=term.as_independent(t,as_Add=False)
            if dependent==1:slope=sy.Integer(0);offset=sy.Integer(0)
            elif dependent.func==sy.exp:
                exponent=sy.expand(dependent.args[0]);slope=sy.diff(exponent,t)
                offset=sy.simplify(exponent-slope*t)
                if slope.has(t) or offset.has(t):raise AssertionError('Nonlinear beta exponential body')
            else:raise AssertionError('Unproved beta moment body '+str(dependent))
            if scalar.has(t) or scalar.has(self.beta):raise AssertionError('Uncancelled beta/source factor')
            out+=scalar*sy.exp(offset)*self.J(sy.simplify(slope/self.ell))
        return sy.factor(out)
    def at(self,node,local=None):
        if local is not None and node in local:return local[node]
        if node in self.bindings:return self.bindings[node]
        row=self.g.nodes[node];op=row['operation'];at=lambda i:self.at(i,local)
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='sum':return sy.Add(*(at(i) for i in row['arguments']))
        if op=='product':return sy.Mul(*(at(i) for i in row['arguments']))
        if op=='negative':return -at(row['argument'])
        if op=='positive_quotient':return sy.cancel(at(row['numerator'])/at(row['denominator']))
        if op=='analytic_unary':
            q=at(row['argument']);name=row['name']
            if name=='exprel':return sy.Integer(1) if q==0 else (sy.exp(q)-1)/q
            if name in ('exp','log','sqrt'):return getattr(sy,name)(q)
        if op=='compact_raw_beta':
            assert row['definition']=='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise'
            assert row['defining_module']==current.PREFIX+'outer_pulse_map.py'
            assert row['defining_module_sha256']==current.sha(row['defining_module'])
            return self.beta(sy.simplify(at(row['argument'])))
        if op=='definite_integral':
            if at(row['lower'])!=-1 or at(row['upper'])!=1:
                raise AssertionError('Canonical beta integral requires original [-1,1] domain')
            var=next(i for i,r in enumerate(self.g.nodes) if r==dict(operation='bound_variable',name=row['variable']))
            t=sy.Symbol('normalized_beta_t',real=True)
            return self.integrate(self.at(row['integrand'],{**(local or {}),var:t}),t)
        raise AssertionError('Unbound exact physical integration node '+op)


def actual_integral_response(owner):
    """The 25 actual basis densities reduce to the same physical matrix."""
    f=owner.full;g=owner.g;W=f.linear.functions;bg=f.linear.background
    mu=sy.Symbol('mu',positive=True);Rc=sy.Symbol('Rc',positive=True)
    ell=sy.Symbol('ell',positive=True);centers=sy.symbols('c0:5',positive=True)
    A=sy.Symbol('Astar',positive=True);t=sy.Symbol('translated_beta_t',real=True)
    basic={W['mu'].node:mu,bg['logRc']:sy.log(Rc),W['ell'].node:ell,
        bg['Astar_ordinary_Z0_Z3'][0]:A,bg['Astar_ordinary_Z0_Z3'][1]:sy.Integer(0),
        **{q.node:c for q,c in zip(W['centers'],centers)}}
    alg=BetaIntegralAlgebra(g,basic,ell)
    assert alg.at(W['normal'].node)==alg.J(sy.Rational(1,2))
    assert sy.expand_log(alg.at(W['log2Rc'].node)-sy.log(2*Rc),force=True)==0
    # Verify the actual rational geometry separately from the general
    # translated-integral theorem; no arbitrary centers are selected.
    L=sy.Symbol('log_two',positive=True);geometry=BetaIntegralAlgebra(g,{g.unary('log',g.constant(2)).node:L},L/400)
    assert geometry.at(W['ell'].node)==L/400
    for j,q in zip(current.own.full.linear.CENTERS,W['centers']):assert geometry.at(q.node)==j*L/40
    supports=current.own.full.linear.exact_supports()
    assert supports==f.reserved_background['exact_new_and_old_supports']
    assert all(0<Fraction(a)<Fraction(b)<1 for a,b in supports['new_log_supports_over_log2'])
    assert list(f.coefficients)==list(current.own.full.COEFFICIENTS)
    checked=[];zero_columns=0
    for col,j in enumerate((*current.own.full.linear.AXIAL,*current.own.full.linear.ANGULAR)):
        bindings={**basic,owner.own.s.node:centers[j]+ell*t}
        for i,pair in enumerate(f.coefficients.values()):
            bindings[pair.value.node]=sy.Integer(i==col);bindings[pair.Z.node]=sy.Integer(0)
        basis=BetaIntegralAlgebra(g,bindings,ell)
        for key,pair in owner.own.patch_bodies['per_ds'].items():
            # This is the actual dR=R ds density and actual normalized
            # physical compact bump, evaluated at s=c_j+ell*t.
            translated=sy.expand_power_base(sy.powsimp(basis.at(pair.value.node)*ell,force=True),force=True)
            integral=basis.integrate(translated,t)
            matrix=basis.at(f.response_matrix[key][col].value.node)
            difference=sy.expand_power_base(sy.powsimp(integral-matrix,force=True),force=True)
            assert sy.simplify(sy.expand(difference))==0,(key,col,str(difference))
            if integral==0:zero_columns+=1
            checked.append(dict(moment=key,column=col,physical_bump_index=j,
                exact_integral=str(integral),physical_response=str(matrix),passed=True))
    # Ordinary Z1 follows the actual C1 bodies and Z-independent endpoints,
    # not differentiation of interval endpoints or numerical quadrature.
    upstream=current.parent.linear.read(current.own.RECEIPT)['checks']['densities']
    assert upstream['passed'] and len(upstream['actual_selected_density_graph_rows'])==10
    for key,pair in owner.own.patch_histories.items():
        for order,ref in enumerate((pair.value,pair.Z)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['endpoint_independent_of_Z']
            assert row['ordinary_slow_Z_derivative_order']==order
            assert row['integrand']==(owner.own.patch_bodies['bodies'][key].value,
                owner.own.patch_bodies['bodies'][key].Z)[order].node
    # Check the matrix's ordinary amplitude product rule independently.
    az=sy.Symbol('Astar_Z',real=True);matrix_alg=BetaIntegralAlgebra(g,{**basic,bg['Astar_ordinary_Z0_Z3'][1]:az},ell)
    derivative_rows=[]
    for key,row in f.response_matrix.items():
        for col,pair in enumerate(row):
            value=matrix_alg.at(pair.value.node);derivative=matrix_alg.at(pair.Z.node)
            assert sy.simplify(derivative-sy.diff(value,A)*az)==0,(key,col,'Z1')
            derivative_rows.append(key+'_column'+str(col)+'_Z1')
    return dict(passed=True,actual_translated_physical_basis_integrals=checked,
        basis_rows=len(checked),inactive_zero_columns=zero_columns,
        actual_C1_matrix_amplitude_derivative_rows=derivative_rows,
        all_five_value_and_Z1_patch_integral_response_identities=True,
        actual_beta_normalizer_and_support_domain_bound=True,
        divided_axial_exprel_integrand_reduced_exactly=True,
        removable_t_zero_singularity_resolved_by_entire_continuation=True,
        no_numeric_weights_or_enclosure_centers_consumed=True)


class EndpointAlgebra:
    """Expand transported end functions using only established integral equalities."""
    def __init__(self,g,bindings,s):self.g,self.bindings,self.s=g,bindings,s
    def at(self,node):
        if node in self.bindings:return self.bindings[node]
        row=self.g.nodes[node];op=row['operation'];at=self.at
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='sum':return sy.Add(*(at(i) for i in row['arguments']))
        if op=='product':return sy.Mul(*(at(i) for i in row['arguments']))
        if op=='negative':return -at(row['argument'])
        if op=='function_substitution':
            assert row['variable']==self.s.node and row['Z_independent_substitution']
            assert sy.simplify(BetaIntegralAlgebra(self.g,{},1).at(row['value'])-sy.log(2))==0
            return at(row['expression'])
        raise AssertionError('Unbound actual terminal partition node '+op)


def transported_partition(owner):
    g=owner.g;f=owner.full;a=current.parent.C1Algebra(g);bindings={};zero=g.zero.node
    bind=lambda ref:sy.Integer(0) if ref.node==zero else bindings.setdefault(ref.node,sy.Symbol('actual_source_'+str(ref.node)))
    pieces=(f.segments['inner'],f.lower.functions,*f.preRp.pieces.values(),f.pulse.functions)
    for piece in pieces:
        for pair in piece.values():
            bind(pair.value);bind(pair.Z)
    for key,pair in owner.own.base_pressure_prefix.items():
        old=f.preRp.pieces['compact'][key]
        bindings[pair.value.node]=bind(old.value);bindings[pair.Z.node]=bind(old.Z)
    for key,pair in owner.own.patch_histories.items():
        for ref,reduced in zip((pair.value,pair.Z),(owner.reduced_patch[key].value,owner.reduced_patch[key].Z)):
            symbol=sy.Symbol('proved_patch_integral_'+key+'_'+str(ref.node))
            bindings[ref.node]=symbol;bindings[reduced.node]=symbol
    endpoint=EndpointAlgebra(g,bindings,owner.own.s);checked=[]
    for key in MOMENTS:
        for label,actual,reduced in (('value',owner.terminal_sources[key].value,owner.reduced_terminal[key].value),
            ('Z1',owner.terminal_sources[key].Z,owner.reduced_terminal[key].Z)):
            assert sy.expand(endpoint.at(actual.node)-endpoint.at(reduced.node))==0,(key,label)
            checked.append(key+'_'+label)
    complete,_=current.own.pressure_prefix_in_band(g,f,owner.own.s,g.constant(2))
    assert complete==f.preRp.pieces['compact']
    return dict(passed=True,actual_transported_terminal_rows=checked,
        inner_lower_active_compact_quiet_pulse_counted_once=True,
        own_pressure_datum_not_regauged=True,m15_negative_weight_not_negated_again=True,
        actual_complete_compact_prefix_identity_consumed=True,
        patch_integral_response_identity_required_before_this_reduction=True)


def tail_origin_identity(owner):
    co=owner.own.coefficients;numeric=co.full;native=current.own.coefficients.full.pulse.native
    frame=numeric.pulse.native.selected.inlet.exact.frame
    importer=current.own.coefficients.full.true.parent.GraphImport(owner.g,frame.graph)
    params=frame.bridge.leading.parameters
    parameter=lambda key:importer(current.source.FunctionRef(frame.graph,params[key].node))
    mu=sy.Symbol('original_mu',positive=True);LP=sy.exp(40)+11
    bindings={parameter('mu').node:mu,parameter('logP').node:LP,
        parameter('logC').node:sy.Symbol('original_logCstar',positive=True)}
    for key in ('T','hbB','hbS','sc'):bindings[parameter(key).node]=sy.Symbol('original_'+key,positive=True)
    L=native.scalar_expression(owner.g,owner.tail_origins['mass_log_radius'],bindings)
    RV=native.scalar_expression(owner.g,owner.tail_origins['pressure_log_radius'],bindings)
    desired=13/mu-2-60*sy.log(mu)-sy.log(2)
    assert sy.simplify(RV-L-desired)==0
    actual_mu=native.scalar_expression(owner.g,owner.full.linear.functions['mu'],{})
    assert sy.simplify(actual_mu-sy.exp(-4*LP)/1000)==0
    # LP>0 gives 0<mu<1/1000<1/2; log(2)>1/2 gives
    # log(Rp/(2Rc)) >= 59*log(2)-2 > 55/2 > 0.
    return dict(passed=True,actual_Rv_over_2Rc_log_identity=str(desired),
        actual_original_mu_recipe_positive_and_less_than_half=True,
        Rp_over_2Rc_log_strict_lower_bound='55/2',Rv_strictly_beyond_2Rc=True,
        scalar_source_graph_and_phase_origin_bound=True)


def compact_own_tail(owner):
    """Actual beta recipes imply all own compact velocity C1 rows vanish."""
    g=owner.g;f=owner.full;W=f.linear.functions
    L=sy.Symbol('log_two',positive=True);u=sy.Symbol('tail_distance',nonnegative=True)
    Rc=sy.Symbol('Rc',positive=True)
    bindings={f.linear.background['logRc']:sy.log(Rc),owner.own.s.node:L*(1+u),
        W['ell'].node:L/400,g.unary('log',g.constant(2)).node:L}
    for j,center in zip(current.own.full.linear.CENTERS,W['centers']):bindings[center.node]=j*L/40
    for i,pair in enumerate(f.coefficients.values()):
        bindings[pair.value.node]=sy.Symbol('actual_coefficient_'+str(i))
        bindings[pair.Z.node]=sy.Symbol('actual_coefficient_Z_'+str(i))
    class TailBetaAlgebra(BetaIntegralAlgebra):
        def at(self,node,local=None):
            if self.g.nodes[node]['operation']=='compact_raw_beta':
                value=super().at(node,local)  # Validate the same actual beta source recipe.
                arg=sy.simplify(value.args[0])
                if sy.simplify(arg-1).is_nonnegative:return sy.Integer(0)
                return value
            return super().at(node,local)
    alg=TailBetaAlgebra(g,bindings,L/400);checked=[]
    for key in ('Uz1','Utheta1','F1'):
        pair=owner.own.patch[key]
        for order,ref in enumerate((pair.value,pair.Z)):
            assert sy.simplify(alg.at(ref.node))==0,(key,order)
            checked.append(key+'_Z'+str(order))
    return dict(passed=True,actual_compact_velocity_C1_rows=checked,
        physical_radius_domain='R >= 2Rc; s=log(2)*(1+nonnegative_distance)',
        least_actual_normalized_beta_argument=240,actual_beta_zero_branch_consumed=True,
        arbitrary_actual_coefficient_values_and_Z1_rows_retained=True)


def tail_recovery(owner):
    g=owner.g;f=owner.full;a=current.parent.C1Algebra(g)
    assert owner.terminal_sources==owner.own.terminal
    bindings={g.zero.node:sy.Integer(0)}
    for key in ('m11','m13'):
        pair=owner.terminal_sources[key];bindings[pair.value.node]=sy.Integer(0);bindings[pair.Z.node]=sy.Integer(0)
    delta=f.cutoff.inner.parameters['delta'];bindings[delta.node]=sy.Symbol('delta',positive=True)
    bindings[g.symbol('Z').node]=sy.Symbol('Z',real=True)
    bindings[owner.mass_coordinate.node]=sy.Symbol('tail_log_radius',nonnegative=True)
    bindings[owner.tail_origins['mass_log_radius'].node]=sy.Symbol('log2Rc',real=True)
    alg=BetaIntegralAlgebra(g,bindings,1);checked=[]
    for key,pair in {**owner.mass_tail,**owner.pressure_tail}.items():
        refs=(pair.value,pair.Z) if type(pair) is current.source.C1Function else (pair,)
        for order,ref in enumerate(refs):
            assert sy.simplify(alg.at(ref.node))==0,(key,order)
            checked.append(key+'_row'+str(order))
    support=f.beyond_Rv_support
    assert support['endpoint']['angular_and_pressure_not_zeroed']
    assert set(support['future_charts'])==set(current.own.full.support.FUTURE)
    assert all(row['exact_zero_unweighted_and_R_weighted_forcing_integrals']==[0,0]
        for row in support['future_charts'].values())
    return dict(passed=True,actual_mass_mean_K_QV_and_pressure_tail_rows=checked,
        zero_after_actual_terminal_integral_proof_only=True,
        compact_n1_velocities_zero_after_2Rc=True,
        original_future_fixed_N1p_zero_and_n1_F1_zero_imply_constant_P1=True,
        pressure_terminal_origin_is_Rv_not_2Rc=True,
        leading_angular_and_heat_pressure_not_zeroed=True,interior_QV_Z1_not_claimed=True,
        radius=tail_origin_identity(owner),actual_compact_own_velocity_tail=compact_own_tail(owner))


MOMENTS=current.MOMENTS


def run(owner):
    began=time.monotonic();raw=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw['terminal_identity_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks=dict(actual_integrals=actual_integral_response(owner))
    print('ACTUAL_N1_TERMINAL_BUMP_INTEGRALS_REDUCED_TO_PHYSICAL_MATRIX',flush=True)
    checks['transported_partition']=transported_partition(owner)
    proxy=SimpleNamespace(g=owner.g,linear=owner.full.linear,moments=owner.full.moments,
        linear_response_residuals=owner.reduced_terminal)
    checks['selected_inverse']=linear_check.universal_linear_identity(proxy)
    print('TEN_TRUE_TRANSPORTED_N1_TERMINAL_FUNCTION_ZEROS_PROVED',flush=True)
    checks['tail_recovery']=tail_recovery(owner)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_terminal_integral_functions_zero_not_only_linear_response=True,
        global_n1_derivatives_and_joins_complete=False,original_numerical_full_field_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_N1_TERMINAL_INTEGRAL_AND_TAIL_IDENTITIES_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the live true terminal graph')
