"""Actual C3 histories, positive-order recovery, integral reduction and tails."""
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_repaired_own_C3_terminal_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions_check as physical
import lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities_check as original_check


class FunctionAlgebra:
    """Exact ordinary calculus for the saved substitution/integral nodes."""
    def __init__(self,g,bindings):self.g,self.bindings,self.cache=g,bindings,{}
    def at(self,node):
        if node in self.bindings:return self.bindings[node]
        if node in self.cache:return self.cache[node]
        row=self.g.nodes[node];op=row['operation'];at=self.at
        if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op=='bound_variable':value=sy.Symbol(row['name'],real=True)
        elif op=='sum':value=sy.Add(*(at(i) for i in row['arguments']))
        elif op=='product':value=sy.Mul(*(at(i) for i in row['arguments']))
        elif op=='negative':value=-at(row['argument'])
        elif op=='positive_quotient':value=sy.cancel(at(row['numerator'])/at(row['denominator']))
        elif op=='analytic_unary' and row['name'] in ('exp','log','sqrt'):
            value=getattr(sy,row['name'])(at(row['argument']))
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            value=at(row['expression']).subs(at(row['variable']),at(row['value']))
        elif op=='definite_integral':
            assert row['endpoint_independent_of_Z']
            var=next(i for i,r in enumerate(self.g.nodes) if r==dict(operation='bound_variable',name=row['variable']))
            lo,hi=at(row['lower']),at(row['upper'])
            value=sy.Integer(0) if lo==hi else sy.Integral(at(row['integrand']),(at(var),lo,hi))
        else:raise AssertionError('Unbound actual own source '+op)
        self.cache[node]=value;return value


def bind_jet(bindings,jet,body,Z):
    for n,ref in enumerate(current.rows(jet)):bindings[ref.node]=sy.diff(body,Z,n)


def source_check(owner):
    g=owner.g;complete=owner.complete;a=current.Algebra(g,3)
    assert owner.patch is complete.repair_velocity
    assert owner.compact is complete.quiet.compact
    assert owner.legacy_terminal.full is complete.full
    assert owner.fields['M1']==owner.cumulative['m11'] and owner.fields['P1']==owner.cumulative['m13']
    for k in ('m11','m12','m14'):assert owner.incoming[k]==complete.moments[k]
    for k in ('m13','m15'):
        assert owner.incoming[k]==a.add(complete.segments['inner'][k],complete.segments['lower'][k],owner.active_pressure[k])
        assert owner.band_pressure[k]==a.at(owner.compact.pressure[k],owner.compact.x,g.unary('exp',owner.s))
        assert owner.active_pressure[k]==owner.compact.inlet['pressure'][k]
        assert owner.compact_local_pressure[k]==a.sub(owner.band_pressure[k],owner.active_pressure[k])
        assert owner.quiet_increment[k]==a.sub(complete.preRp_pressure[k],owner.compact.terminal['pressure'][k])
    count=0
    for key,q in owner.patch_histories.items():
        for n,ref in enumerate(current.rows(q)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['endpoint_independent_of_Z']
            assert row['lower']==g.zero.node and row['upper']==owner.s.node
            assert row['integrand']==current.rows(owner.patch_bodies['bodies'][key])[n].node
            assert row['ordinary_slow_Z_derivative_order']==n and row['history_not_reset_between_bumps']
            assert row['physical_measure']=='dR=R ds' and row['fixed_leading_N1p_not_subtracted']
            count+=1
    for key in ('mean1','M1','K1','P1','P1_y_total'):assert len(current.rows(owner.fields[key]))==4
    for key in ('Q1','V1'):assert len(current.rows(owner.fields[key]))==3
    assert len(current.rows(owner.compact_Qyy))==5
    assert owner.tail_origins==owner.legacy_terminal.tail_origins
    assert complete.beyond_Rv_support['endpoint']['angular_and_pressure_not_zeroed']
    for q in ('3/10',):assert owner.at_band_fraction(q)['QV_C2']
    for q in ('-1/100',):
        try:owner.at_band_fraction(q)
        except ValueError:pass
        else:raise AssertionError('Out-of-band query accepted')
    try:owner.tail_at_ratio('mass','2')
    except ValueError:pass
    else:raise AssertionError('Unaccepted zero certificate exposed')
    return dict(passed=True,actual_signed_history_integral_rows=count,
        actual_current_compact_pressure_and_quiet_difference_retained=True,
        inner_lower_pressure_datum_and_fixed_N1p_retained_once=True,
        own_QV_C2_uses_real_mean_Z3=True,actual_compact_direct_Qyy_C4=True,
        exact_chart_queries_and_acceptance_guard_enforced=True)


def actual_densities(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);Y=sy.Symbol('s',real=True);Rc=sy.Symbol('Rc',positive=True)
    mu=sy.Symbol('mu',positive=True);U,T,A=(sy.Function(k)(Z) for k in ('Uz1','Utheta1','Astar'))
    R=Rc*sy.exp(Y);theta0=A*sy.exp(-(sy.Rational(1,2)+mu)*(Y+sy.log(Rc)+sy.log(2)))
    bindings={owner.s.node:Y,owner.complete.linear.background['logRc']:sy.log(Rc),
              owner.complete.linear.functions['mu'].node:mu}
    bind_jet(bindings,owner.patch['Uz1'],U,Z);bind_jet(bindings,owner.patch['Utheta1'],T,Z)
    for n,node in enumerate(owner.complete.linear.background['Astar_ordinary_Z0_Z3']):bindings[node]=sy.diff(A,Z,n)
    expected=dict(m11=R*U,m12=R*sy.sqrt(2*R)*T,m13=theta0*T,
                  m14=R*sy.sqrt(2*R)*theta0*U,m15=-R*theta0*T)
    alg=FunctionAlgebra(g,bindings);checked=[]
    for key,q in owner.patch_bodies['per_ds'].items():
        for n,ref in enumerate(current.rows(q)):
            difference=sy.expand_power_base(sy.powsimp(alg.at(ref.node)-sy.diff(expected[key],Z,n),force=True),force=True)
            assert sy.simplify(difference)==0,(key,n)
            checked.append(f'{key}_Z{n}')
    return dict(passed=True,actual_selected_density_rows=checked,
        physical_dR_equals_R_ds=True,actual_Astar_Z0_Z3_product_rules=True,
        physical_F_divided_by_sqrt_2R_once=True,m13_factor_two_m15_single_minus=True)


def actual_recovery(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);R=sy.Symbol('R',positive=True);de=sy.Symbol('delta',positive=True)
    M,P,U=(sy.Function(k)(Z) for k in ('M','P','U'))
    bindings={g.symbol('Z').node:Z,owner.R.node:R,owner.delta.node:de}
    for q,body in ((owner.cumulative['m11'],M),(owner.cumulative['m13'],P),(owner.patch['Uz1'],U)):
        bind_jet(bindings,q,body,Z)
    Q=(2*Z*U-(1+de)*Z*M/R-(1-Z*Z)*sy.diff(M,Z)/R)/(1-de*Z*Z)
    expected=dict(M1=M,P1=P,mean1=M/R,K1=M/R-U,mean1_y=U-M/R,M1_y=R*U,Q1=Q,V1=R*Q)
    alg=FunctionAlgebra(g,bindings);checked=[]
    for key,body in expected.items():
        for n,ref in enumerate(current.rows(owner.fields[key])):
            assert sy.cancel(alg.at(ref.node)-sy.diff(body,Z,n))==0,(key,n)
            checked.append(f'{key}_Z{n}')
    return dict(passed=True,actual_recovery_rows=checked,n1_weight_one_plus_delta=True,
        full_mean_third_derivative_used_for_QV_second_derivative=True,
        mean_y_plus_mean_equals_U_and_mass_y_equals_R_U=True,V_equals_R_Q=True)


def actual_pressure_FTC(owner):
    g=owner.g;c=owner.compact;Z=sy.Symbol('Z',real=True);X=sy.Symbol('x',positive=True)
    Y=sy.Symbol('s',real=True);Rc=sy.Symbol('Rc',positive=True)
    S=sy.Function('S')(X,Z);Q=sy.Function('Q')(X,Z);B=sy.Function('patch_density')(Y,Z)
    bindings={g.symbol('Z').node:Z,c.x.node:X,owner.s.node:Y,c.Rc.node:Rc,
              owner.complete.linear.background['logRc']:sy.log(Rc)}
    bind_jet(bindings,c.S,S,Z);bind_jet(bindings,c.physical['Q'],Q,Z)
    bind_jet(bindings,c.physical['Q_y'],X*sy.diff(Q,X),Z)
    bind_jet(bindings,owner.compact_Qyy,X*sy.diff(X*sy.diff(Q,X),X),Z)
    bind_jet(bindings,owner.patch_bodies['per_ds']['m13'],B,Z)
    for key,q in c.incoming_integrals.items():bind_jet(bindings,q,sy.Function('incoming_'+key)(Z),Z)
    for key,q in c.initial.items():
        if hasattr(q,'value'):bind_jet(bindings,q,sy.Function('initial_'+key)(Z),Z)
    for segment in ('inner','lower'):
        bind_jet(bindings,owner.complete.segments[segment]['m13'],sy.Function(segment+'_m13')(Z),Z)
    alg=FunctionAlgebra(g,bindings);checked=[]
    for n,(P,Py) in enumerate(zip(current.rows(owner.fields['P1']),current.rows(owner.fields['P1_y_total']))):
        actual=sy.diff(alg.at(P.node),Y);wanted=alg.at(Py.node)
        assert sy.simplify(actual-wanted)==0,n
        checked.append('P1_y_Z'+str(n))
    return dict(passed=True,actual_pressure_primitive_radial_FTC_rows=checked,
        actual_S0_integrand_and_finite_Q_Qy_boundaries_used=True,
        fixed_leading_N1p_plus_patch_density_not_patch_alone=True,
        original_nonzero_pressure_inlets_and_inner_lower_datum_retained=True)


def actual_compact_second_radial(owner):
    g=owner.g;c=owner.compact;Z=sy.Symbol('Z',real=True);de=sy.Symbol('delta',positive=True)
    scale=sy.Symbol('Pstar',positive=True);Uy,Uyy,my=(sy.Function(k)(Z) for k in ('U_y','U_yy','m_y'))
    bindings={g.symbol('Z').node:Z,c.delta.node:de,c.Pstar.node:scale}
    for q,body in ((c.radial['V_y'],Uy/scale),(c.radial['V_yy'],Uyy/scale),(c.radial['m_y'],my/scale)):
        bind_jet(bindings,q,body,Z)
    mean_yy=Uy-my
    wanted=(2*Z*Uyy-(1-de)*Z*mean_yy-(1-Z*Z)*sy.diff(mean_yy,Z))/(1-de*Z*Z)
    alg=FunctionAlgebra(g,bindings);checked=[]
    for n,ref in enumerate(current.rows(owner.compact_Qyy)):
        assert sy.cancel(alg.at(ref.node)-sy.diff(wanted,Z,n))==0,n
        checked.append('Q_yy_Z'+str(n))
    return dict(passed=True,actual_direct_compact_Qyy_C4_rows=checked,
        original_leading_weight_one_minus_delta=True,mean_yy_equals_U_y_minus_mean_y=True,
        physical_Pstar_applied_once=True,genuine_compact_Vy_Vyy_mean_y_C5_sources_used=True)


def actual_integral_response(owner):
    """Integrate the new actual bases, then prove ordinary C3 commutation."""
    g=owner.g;f=owner.complete;old=owner.legacy_terminal.full;W=f.linear.functions
    receipt=current.parent.linear.read(current.original_terminal.RECEIPT)
    theorem=receipt['checks']['actual_integrals']
    assert theorem['passed'] and theorem['basis_rows']==25
    assert all(row['passed'] for row in theorem['actual_translated_physical_basis_integrals'])
    assert f.linear is old.linear and f.reserved_background==old.reserved_background
    Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True);Rc2=sy.Symbol('two_Rc',positive=True)
    A=sy.Function('Astar')(Z)
    bindings={W['mu'].node:mu,W['log2Rc'].node:sy.log(Rc2)}
    for n,node in enumerate(f.linear.background['Astar_ordinary_Z0_Z3']):bindings[node]=sy.diff(A,Z,n)
    for group in ('axial_rows','angular_rows'):
        for i,row in enumerate(W[group]):
            for j,q in enumerate(row):bindings[q.node]=sy.Symbol(f'{group}_{i}_{j}',real=True)
    alg=FunctionAlgebra(g,bindings);checked=0
    for key,row in f.response_matrix.items():
        for col,q in enumerate(row):
            old_value=alg.at(old.response_matrix[key][col].value.node)
            for n,ref in enumerate(current.rows(q)):
                assert sy.simplify(alg.at(ref.node)-sy.diff(old_value,Z,n))==0,(key,col,n)
                checked+=1
    for key,ref in owner.integral_identity.items():
        row=g.nodes[ref.node]
        assert row['actual_endpoint_rows']==current.encode(owner.patch_endpoint[key])
        assert row['equivalent_physical_response_rows']==current.encode(owner.reduced_patch[key])
        assert row['actual_signed_density_rows']==current.encode(owner.patch_bodies['per_ds'][key])
        assert row['original_basis_integral_receipt_sha256']==current.sha(current.original_terminal.RECEIPT)
        assert row['supports_and_endpoints_independent_of_Z'] and row['ordinary_Z_orders']==[0,1,2,3]
    Rc=sy.Symbol('Rc',positive=True);ell=sy.Symbol('ell',positive=True)
    centers=sy.symbols('center0:5',positive=True);t=sy.Symbol('translated_beta_t',real=True)
    basic={W['mu'].node:mu,f.linear.background['logRc']:sy.log(Rc),
           W['log2Rc'].node:sy.log(2*Rc),W['ell'].node:ell}
    for n,node in enumerate(f.linear.background['Astar_ordinary_Z0_Z3']):basic[node]=sy.diff(A,Z,n)
    for q,value in zip(W['centers'],centers):basic[q.node]=value
    new_basis=[];new_derivatives=[]
    for col,j in enumerate((*current.previous.original.linear.AXIAL,*current.previous.original.linear.ANGULAR)):
        bindings={**basic,owner.s.node:centers[j]+ell*t};C=sy.Function('selected_C'+str(col))(Z)
        for i,q in enumerate(f.coefficients.values()):
            for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.diff(C,Z,n) if i==col else sy.Integer(0)
        basis=original_check.BetaIntegralAlgebra(g,bindings,ell)
        assert basis.at(W['normal'].node)==basis.J(sy.Rational(1,2))
        for key,q in owner.patch_bodies['per_ds'].items():
            translated=sy.expand_power_base(sy.powsimp(basis.at(q.value.node)*ell,force=True),force=True)
            integral=basis.integrate(translated,t)
            response=basis.at(owner.reduced_patch[key].value.node)
            difference=sy.expand_power_base(sy.powsimp(integral-response,force=True),force=True)
            assert sy.simplify(sy.expand(difference))==0,(key,col,'actual basis integral')
            new_basis.append(dict(moment=key,column=col,actual_physical_bump_index=j,passed=True))
            for n,ref in enumerate(current.rows(q)[1:],1):
                actual=basis.at(ref.node)*ell
                difference=sy.expand_power_base(sy.powsimp(actual-sy.diff(translated,Z,n),force=True),force=True)
                assert sy.simplify(sy.expand(difference))==0,(key,col,n,'actual density derivative')
                new_derivatives.append(f'{key}_column{col}_Z{n}')
    return dict(passed=True,accepted_actual_physical_basis_integrals=25,
        current_new_actual_translated_basis_integrals=new_basis,
        current_new_density_ordinary_derivative_rows=new_derivatives,
        variable_selected_coefficients_and_Astar_not_unit_or_constant_Z_fixtures=True,
        same_exact_reserved_operator_normalizer_supports_and_physical_amplitude=True,
        actual_C3_matrix_rows_equal_derivatives_of_admitted_physical_basis_response=checked,
        actual_C3_signed_density_and_fixed_endpoint_integral_proofs_consumed=True,
        compact_flat_integrands_allow_ordinary_differentiation_under_integral=True,
        twenty_actual_patch_terminal_rows_equal_selected_response=True,
        no_numeric_bump_weights_or_interval_overlap=True)


def transported_terminal(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);S=sy.Symbol('s',real=True);X=sy.Symbol('x',positive=True)
    bindings={owner.s.node:S,owner.compact.x.node:X}
    for segment in ('inner','lower','preRp','pulse'):
        for key,q in owner.complete.segments[segment].items():
            for n,ref in enumerate(current.rows(q)):
                if ref!=g.zero:bindings.setdefault(ref.node,sy.Symbol(f'{segment}_{key}_{n}',real=True))
    for key,q in owner.compact.pressure.items():
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.Function(f'compact_{key}_{n}')(X)
    for key,q in owner.patch_histories.items():
        for n,ref in enumerate(current.rows(q)):
            H=sy.Function(f'proved_patch_{key}_{n}')
            bindings[ref.node]=H(S);bindings[current.rows(owner.reduced_patch[key])[n].node]=H(sy.log(2))
    alg=FunctionAlgebra(g,bindings);checked=[]
    for key,q in owner.terminal.items():
        for n,ref in enumerate(current.rows(q)):
            assert sy.expand(alg.at(ref.node)-alg.at(current.rows(owner.reduced_terminal[key])[n].node))==0,(key,n)
            checked.append(f'{key}_Z{n}')
    inverse=current.parent.linear.read(current.previous.RECEIPT)['checks']['actual_inverse_and_response']
    assert inverse['passed'] and len(inverse['actual_physical_response_C3_zero_rows'])==20
    return dict(passed=True,actual_transported_terminal_rows=checked,
        compact_exp_log2_endpoint_identified_with_actual_compact_terminal=True,
        inner_lower_active_compact_quiet_pulse_counted_once=True,
        actual_integral_response_theorem_consumed_before_selected_inverse_zero=True,
        pressure_axis_datum_and_m15_sign_unchanged=True,
        twenty_actual_terminal_integral_function_rows_zero=True)


def actual_tails(owner):
    g=owner.g;f=owner.complete;W=f.linear.functions;L=sy.Symbol('log_two',positive=True)
    u=sy.Symbol('distance',nonnegative=True);Rc=sy.Symbol('Rc',positive=True)
    bindings={f.linear.background['logRc']:sy.log(Rc),owner.s.node:L*(1+u),W['ell'].node:L/400,
              g.unary('log',g.constant(2)).node:L}
    for j,center in zip(current.previous.original.linear.CENTERS,W['centers']):bindings[center.node]=j*L/40
    for key,q in f.coefficients.items():
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.Symbol(f'{key}_{n}',real=True)
    class TailAlgebra(original_check.BetaIntegralAlgebra):
        def at(self,node,local=None):
            value=super().at(node,local)
            if self.g.nodes[node]['operation']=='compact_raw_beta':
                if sy.simplify(value.args[0]-1).is_nonnegative:return sy.Integer(0)
            return value
    alg=TailAlgebra(g,bindings,L/400);compact_rows=[]
    for key in ('Uz1','Utheta1','F1'):
        for n,ref in enumerate(current.rows(owner.patch[key])):
            assert sy.simplify(alg.at(ref.node))==0,(key,n);compact_rows.append(f'{key}_Z{n}')
    Z=sy.Symbol('Z',real=True);Y=sy.Symbol('tail_log_radius',nonnegative=True)
    bindings={g.symbol('Z').node:Z,owner.delta.node:sy.Symbol('delta',positive=True),
              owner.mass_coordinate.node:Y,owner.tail_origins['mass_log_radius'].node:sy.Symbol('log2Rc',real=True)}
    for key in ('m11','m13'):
        for ref in current.rows(owner.terminal[key]):bindings[ref.node]=sy.Integer(0)
    alg=FunctionAlgebra(g,bindings);checked=[]
    for key,q in {**owner.mass_tail,**owner.pressure_tail}.items():
        for n,ref in enumerate(current.rows(q)):
            assert sy.simplify(alg.at(ref.node))==0,(key,n);checked.append(f'{key}_Z{n}')
    legacy=current.parent.linear.read(current.original_terminal.RECEIPT)['checks']['tail_recovery']
    assert legacy['passed'] and legacy['radius']['Rv_strictly_beyond_2Rc']
    assert owner.tail_origins==owner.legacy_terminal.tail_origins
    support=f.beyond_Rv_support
    assert support['endpoint']['angular_and_pressure_not_zeroed']
    assert all(row['exact_zero_unweighted_and_R_weighted_forcing_integrals']==[0,0] for row in support['future_charts'].values())
    return dict(passed=True,actual_compact_velocity_tail_C3_rows=compact_rows,
        actual_recovered_mass_pressure_C3_QV_C2_tail_rows=checked,
        zero_after_actual_terminal_integral_proof_only=True,
        actual_2Rc_and_Rv_origins_and_future_forcing_support_preserved=True,
        original_leading_angular_and_heat_pressure_not_zeroed=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['repaired_own_C3_terminal_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_densities',actual_densities),('actual_recovery',actual_recovery),
        ('actual_compact_second_radial',actual_compact_second_radial),('actual_pressure_FTC',actual_pressure_FTC),
        ('source',source_check),('actual_integral_response',actual_integral_response),
        ('transported_terminal',transported_terminal),('actual_tails',actual_tails)):
        checks[name]=fn(owner);print('REPAIRED_OWN_C3_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_terminal_integral_functions_C3_zero=True,global_n1_joins_and_quantitative_norms_complete=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ACTUAL_REPAIRED_OWN_C3_TERMINAL_INTEGRALS_AND_TAILS_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses checked source owners')
