"""Actual full-bump terminal integrals and the genuine compact branch.

No terminal field or residual is bound to zero. Actual integrands are
projected using translated defining bump moments, then the current branch
equation and its actual IFT rows prove the resulting residuals vanish.
"""
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source as current


def prove(owner,Expression,zero):
    c=owner.inputs[current.compact.NAME];g=owner.g;ref=owner.ref;W=c['W']
    Z=sy.Symbol('Z',real=True);mu,N,J0=sy.symbols('mu N J0',positive=True)
    ell=sy.log(2)/40;centers=[sy.log(2)*sy.Rational(s) for s in ('1/5','1/2','4/5')]
    H=sy.symbols('Htheta Henergy Hpressure',positive=True);Ha=sy.Symbol('Hax',positive=True)
    Gram=sy.symbols('Gcross Genergy Gpressure',positive=True)
    base={g.symbol('Z').node:Z,c['mu']:mu,c['N']:N,W['normal']:J0}
    assert ell>0 and all(0<ci-ell and ci+ell<sy.log(2) for ci in centers)
    assert all(abs(a-b)>2*ell for i,a in enumerate(centers) for b in centers[i+1:])
    actual=current.compact.compact.exact_weights(g,ref(c['mu']))
    assert current.previous.encode(actual)==W
    canonical_bumps=current.compact.bump_rows(g,actual,ref(c['background']['phase']))
    assert current.previous.encode(canonical_bumps)==c['profiles']['bumps']
    assert c['bindings']['original_weights']==current.parent.active.current.ast_binding(current.compact.compact.exact_weights)
    class WeightExpression(Expression):
        def at(self,i):
            i=i.node if hasattr(i,'node') else i
            if i in self.bindings:return self.bindings[i]
            row=self.g.nodes[i]
            if row['operation']=='compact_raw_beta':return sy.Function('raw_beta')(super().at(row['argument']))
            return super().at(i)
    w=WeightExpression(g,base);normal=g.nodes[W['normal']]
    v=sy.Symbol(normal['variable'],real=True);beta=sy.Function('raw_beta')(v)
    def integral(i,wanted):
        row=g.nodes[i]
        assert row['operation']=='definite_integral' and row['exact_function_integral'] and row['numerical_value_not_installed']
        assert row['variable']==normal['variable'] and w.at(row['lower'])==-1 and w.at(row['upper'])==1
        assert zero(w.at(row['integrand'])-wanted),i
    integral(W['normal'],beta)
    powers=(sy.Rational(1,2),-sy.Rational(1,2)-mu,-sy.Rational(3,2)-mu)
    for i,p,h in zip(W['H'],powers,H):
        row=g.nodes[i];assert row['operation']=='positive_quotient' and row['denominator']==W['normal']
        integral(row['numerator'],beta*sy.exp(p*ell*v));base[row['numerator']]=J0*h
    row=g.nodes[W['Hax']];assert row['denominator']==W['normal']
    integral(row['numerator'],beta*sy.exp(-mu*ell*v));base[row['numerator']]=J0*Ha
    # D_i is exactly (translated negative-mu moment - 1)/mu.
    for i,ci in zip(W['D'],(centers[0],centers[2])):
        row=g.nodes[i];assert row['denominator']==W['normal']
        integral(row['numerator'],beta*(sy.exp(-mu*(ci+ell*v))-1)/mu)
        base[row['numerator']]=J0*(sy.exp(-mu*ci)*Ha-1)/mu
    def closure(roots):
        seen=set();stack=list(roots)
        while stack:
            i=stack.pop()
            if i in seen:continue
            seen.add(i);row=g.nodes[i]
            if row['operation'] not in current.SCHEMA:raise ValueError('Unknown weight operand')
            for key in current.SCHEMA[row['operation']]:
                value=row.get(key)
                if value is not None:stack.extend(value if isinstance(value,list) else [value])
        return seen
    roots=W['cross']+W['energy']+W['pressure'];found=set()
    for i in closure(roots):
        row=g.nodes[i]
        if row['operation']!='definite_integral' or i==W['normal']:continue
        body=w.at(row['integrand'])
        for j,p in enumerate((-sy.Rational(1,2),-1,-2)):
            if zero(body-beta**2*sy.exp(p*ell*v)):
                integral(i,beta**2*sy.exp(p*ell*v));base[i]=ell*J0**2*Gram[j];found.add(j);break
        else:raise AssertionError('Unrecognized actual quadratic bump primitive')
    assert found=={0,1,2}
    weights=WeightExpression(g,base)
    assert zero(weights.at(W['ell'])-ell)
    for i,ci in zip(W['centers'],centers):assert zero(weights.at(i)-ci)
    B=sy.Matrix([[weights.at(i) for i in row] for row in W['B']])
    for j,ci in enumerate(centers):
        assert zero(B[2,j+2]-H[0]*sy.exp(ci/2))
        assert zero(B[3,j+2]+H[1]*sy.exp(powers[1]*ci))
        assert zero(B[4,j+2]-H[2]*sy.exp(powers[2]*ci))
        assert zero(weights.at(W['energy'][j])-Gram[1]*sy.exp(-ci))
        assert zero(weights.at(W['pressure'][j])-Gram[2]*sy.exp(-2*ci))
    for i,j in enumerate((0,2)):assert zero(weights.at(W['cross'][i])-Gram[0]*sy.exp(-centers[j]/2))
    # Validate the actual branch definition, not a certificate's zero_rows.
    system=g.nodes[c['system']];transfer=c['branch_transfer'];admission=current.parent.linear.read(transfer['original_admission_receipt'])
    geometry=current.compact.previous.original_geometry
    assert transfer['original_admission_receipt']==geometry.RECEIPT
    original_admission=current.parent.linear.read(geometry.NAME)
    for name in (geometry.NAME,geometry.RECEIPT):assert owner.hashes[name]==current.sha(name)
    assert admission['all_passed'] and admission[geometry.GATE] and admission['source_family']==owner.family
    assert original_admission['source_family']==owner.family
    original_proof=original_admission['actual_whole_Z_frequency_connection']
    assert transfer['original_C1_limit_proof']==original_proof['C1_limit']
    assert transfer['selected_integer']==original_proof['selected_integer']
    assert transfer['same_original_B_Q_N_and_exact_bump_integrals']
    assert transfer['original_C1_contraction_not_new_C5_numeric_norm_admission']
    assert transfer['original_C1_limit_proof']['exact_C1_limit_exists_unique_for_this_repair_only_integer']
    assert transfer['original_C1_limit_proof']['contraction_at_most']=='1/2'
    assert transfer['selected_integer']['J']==g.nodes[c['N']]['J']
    assert system['operation']=='original_compact_unique_Banach_limit_system'
    assert system['actual_selected_N']==c['N'] and system['original_weights']==W and system['original_branch_transfer']==transfer
    assert system['source_family']==owner.family and system['actual_target_values']==[ids[0] for ids in c['targets'].values()]
    for j,ids in enumerate(c['h']):
        row=g.nodes[ids[0]]
        assert row['operation']=='original_compact_unique_Banach_limit_component' and row['system']==c['system']
        assert row['component']==j and row['mathematical_function_not_numeric_interval'] and row['source_family']==owner.family
        assert row['Z_variable']==g.symbol('Z').node
    generic=[ref(i) for i in system['variables']]
    Wlive=actual;Q=current.compact.compact.quadratic(g,Wlive,generic)
    vector=[g.add(ref(ids[0]),g.quotient(q,ref(c['N']),'same original exact positive integer N')) for ids,q in zip(c['targets'].values(),Q)]
    actual_map=[g.neg(q) for q in current.compact.compact.inverse_action(g,Wlive,vector)]
    assert [q.node for q in actual_map]==system['exact_map']
    scalar_h=sy.symbols('h0:5');scalar_d=sy.symbols('d0:5');scalar=dict(base)
    for i,h in zip(system['variables'],scalar_h):scalar[i]=h
    for ids,h in zip(c['h'],scalar_h):scalar[ids[0]]=h
    for ids,d in zip(c['targets'].values(),scalar_d):scalar[ids[0]]=d
    calc=WeightExpression(g,scalar);T=sy.Matrix([calc.at(i) for i in system['exact_map']])
    F=sy.Matrix([calc.at(ids[0]) for ids in c['residual'].values()])
    for value in F-B*(sy.Matrix(scalar_h)-T):assert zero(value)
    # All higher actual residual rows have J*h_n + known_rhs. The actual
    # stored h_n is -J^-1*known_rhs; J is the same admitted invertible one.
    jet_bind=dict(base);hs=[[sy.Symbol(f'h{j}_{n}') for n in range(6)] for j in range(5)]
    ds=[[sy.Symbol(f'd{j}_{n}') for n in range(6)] for j in range(5)]
    for ids,values in zip(c['h'],hs):jet_bind.update(zip(ids,values))
    for ids,values in zip(c['targets'].values(),ds):jet_bind.update(zip(ids,values))
    jet=WeightExpression(g,jet_bind);JF=sy.Matrix([[jet.at(i) for i in row] for row in c['implicit']['Jacobian']])
    F0=sy.Matrix([jet.at(ids[0]) for ids in c['residual'].values()])
    for i in range(5):
        for j in range(5):assert zero(JF[i,j]-sy.diff(F0[i],hs[j][0]))
    inverse,det=current.compact.leading.inverse_matrix(g,[[ref(i) for i in row] for row in c['implicit']['Jacobian']])
    assert current.previous.encode(inverse)==c['implicit']['inverse'] and det.node==c['implicit']['determinant']
    higher=[]
    for n,step in enumerate(c['implicit']['recurrence'],1):
        assert step['order']==n and step['same_full_Jacobian'] and step['only_unknown_highest_row_temporarily_zero']
        for j,ids in enumerate(c['h']):
            solved=g.neg(g.add(*(g.mul(inverse[j][k],ref(step['known_rhs'][k])) for k in range(5))))
            assert solved.node==ids[n]==step['solved_ordinary_row'][j]
        for j,ids in enumerate(c['residual'].values()):
            assert zero(jet.at(ids[n])-sum(JF[j,k]*hs[k][n] for k in range(5))-jet.at(step['known_rhs'][j]))
            higher.append(f'residual_{j}_Z{n}')
    # Actual partial integrands, projected at x=2 by full supports. The
    # physical beta_i(t) already includes 1/(ell*J0*t); its moments have
    # exactly the translated weights verified above.
    A=sy.Function('A')(Z);h=[sy.Function('control_'+str(j))(Z) for j in range(5)]
    incoming={key:sy.Function('incoming_'+key)(Z) for key in current.outer.RATES}
    physical=dict(base)
    def bind(ids,body):
        for n,i in enumerate(ids):
            if i==0:assert sy.diff(body,Z,n)==0
            else:physical[i]=sy.diff(body,Z,n)
    bind(c['A'],A)
    for ids,body in zip(c['h'],h):bind(ids,body)
    for key,ids in c['incoming_correction'].items():bind(ids,incoming[key])
    betas=sy.symbols('b0:3')
    for triplet,beta_i in zip(c['profiles']['bumps'],betas):physical[triplet[0]]=beta_i
    variable=c['background']['phase'];physical[variable]=sy.Integer(2)
    # Rates m/e coincide, so choose the chart's actual key via its integral
    # node identity rather than guessing from the recovery rate.
    integral_keys={i:key for key,ids in c['contributions'].items() for i in ids}
    class Terminal(WeightExpression):
        def __init__(self,g,bindings,phase_value=sy.Integer(2)):
            super().__init__(g,bindings);self.phase_value=phase_value
        def at(self,i):
            i=i.node if hasattr(i,'node') else i
            if i in self.bindings:return self.bindings[i]
            row=self.g.nodes[i]
            if row['operation']=='function_substitution':
                assert row['Z_independent_substitution'] and row['variable']==variable
                value=self.at(row['value'])
                if value==2:assert row['value']==g.constant(2).node
                elif value==1:
                    assert row['value']==g.one.node
                    assert row['expression']==g.unary('log',ref(variable)).node
                else:
                    assert value==self.phase_value and isinstance(value,sy.Symbol)
                    assert row['value']==g.symbol(str(value)).node
                b=dict(self.bindings);b[row['variable']]=value
                return Terminal(self.g,b,self.phase_value).at(row['expression'])
            if row['operation']=='definite_integral':
                assert i in integral_keys and row['finite_N_correction'] and row['physical_Jacobian_applied_once']
                assert self.at(row['lower'])==1 and self.at(row['upper'])==2
                key=integral_keys[i];rate=sy.Rational(str(current.outer.RATES[key]))
                assert sy.Rational(row['original_rate'])==rate
                t=sy.Symbol(row['variable'],positive=True);b=dict(self.bindings)
                b[g.symbol(row['variable']).node]=t
                body=Terminal(g,b,t).at(row['integrand']);poly=sy.Poly(sy.expand(body),*betas)
                out=sy.Integer(0)
                for powers_beta,coefficient in poly.terms():
                    degree=sum(powers_beta);assert degree in (1,2)
                    if degree==1:
                        j=powers_beta.index(1);p=dict(m=0,h=sy.Rational(1,2),k=-mu,e=-sy.Rational(1,2)-mu,p=-sy.Rational(3,2)-mu)[key]
                        moment=(sy.Integer(1) if key=='m' else weights.at(W['B'][2][j+2]) if key=='h'
                            else 1+mu*weights.at(W['D'][0 if j==0 else 1]) if key=='k' and j in (0,2)
                            else -weights.at(W['B'][3][j+2]) if key=='e' else weights.at(W['B'][4][j+2]) if key=='p' else None)
                        assert moment is not None
                    else:
                        p=dict(k=sy.Rational(1,2),e=0,p=-1)[key]
                        if 2 in powers_beta:
                            j=powers_beta.index(2)
                            moment=weights.at(W['cross'][0 if j==0 else 1]) if key=='k' and j in (0,2) else weights.at(W['energy'][j]) if key=='e' else weights.at(W['pressure'][j]) if key=='p' else None
                            assert moment is not None
                        else:
                            indices=[j for j,power in enumerate(powers_beta) if power==1]
                            assert len(indices)==2 and all(power in (0,1) for power in powers_beta)
                            assert abs(centers[indices[0]]-centers[indices[1]])>2*ell
                            moment=sy.Integer(0) # the pinned canonical bump supports are disjoint
                    independent=sy.cancel(sy.powsimp(sy.expand_power_exp((coefficient/t**p).rewrite(sy.exp)),force=True))
                    assert not independent.has(t),(key,powers_beta,'wrong physical moment measure')
                    out+=independent*moment
                return out
            return super().at(i)
    terminal=Terminal(g,physical);projected=[]
    for key,ids in c['terminal']['correction'].items():
        for n,(i,rhs) in enumerate(zip(ids,c['terminal_rhs'][key])):
            assert zero(terminal.at(i)-terminal.at(rhs)),(key,n,'actual integral-to-residual')
            projected.append(f'{key}_Z{n}')
    return dict(passed=True,actual_integral_to_residual_C5_rows=projected,
        actual_branch_value_equation_rows=5,actual_implicit_residual_derivative_rows=higher,
        current_exact_map_is_same_actual_B_inverse_target_and_quadratic=True,
        actual_full_bump_kernel_normalization_translation_and_disjoint_supports=True,
        all_nine_actual_bump_profile_rows_pinned_to_canonical_center_width_and_argument=True,
        actual_compact_substitutions_limited_to_inlet_log_terminal_and_integral_variable=True,
        actual_unique_branch_canonical_Z_and_original_admission_hashes_pinned=True,
        actual_relative_terminal_functions_zero_from_integrals_and_unique_branch=True,
        absolute_histories_and_interior_relative_histories_not_zeroed=True)
