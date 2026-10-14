"""Actual compact/own radial calculus and source joins, without numeric probes."""
from pathlib import Path
import json
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_compact_higher_radial_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_pulse_higher_radial_pressure_C3_functions_check as high_check

outer_check=high_check.outer_check
zero,bind=high_check.zero,high_check.bind


class CompactAlgebra(high_check.HigherAlgebra):
    def __init__(self,owner,bindings):
        super().__init__(owner.high,bindings);self.compact_owner=owner
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation']
        if op=='compact_raw_beta':
            result=high_check.radial_check.Beta(sy.expand(self.at(row['argument'])))
        elif op=='analytic_unary' and row['name']=='exprel':
            x=self.at(row['argument']);result=sy.Integer(1) if x==0 else (sy.exp(x)-1)/x
        elif op=='original_lazy_flat_branch':
            assert row['original_beta_exact_flat_derivative'] and row['active_body_not_evaluated_on_flat_branch']
            assert row['eta']==1 and row['flat_value']==0
            x=sy.simplify(self.at(row['Delta']))
            if (x<=-1) is sy.true or (x>=1) is sy.true:
                self.flat_nodes.add(i);result=sy.Integer(0)
            else:result=self.at(row['active_body'])
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            var=self.g.nodes[row['variable']];assert var['operation']=='bound_variable'
            symbol=sy.Symbol(var['name'],real=True);value=self.at(row['value'])
            env={key:body.subs(symbol,value) if hasattr(body,'subs') else body for key,body in self.bindings.items()}
            env[row['variable']]=value
            child=CompactAlgebra(self.compact_owner,env);result=child.at(row['expression'])
            self.flat_nodes.update(child.flat_nodes)
        else:return super().at(i)
        self.cache[i]=result;return result


def contract(owner):
    owner.assert_graph();pins=high_check.ordinary_source_contract(owner.high)
    for module,receipt in ((current.compact,current.compact.RECEIPT),(current.outer,current.outer.CONSTRUCTION_RECEIPT)):
        for name in (module.NAME,receipt,Path(module.__file__).name,Path(module.__file__).stem+'_check.py'):
            assert owner.hashes[name]==current.sha(name),name
    for rows in (owner.profiles,owner.histories):
        assert len(rows)==5 and all(len(current.outer.rows(jet))==6 for row in rows for jet in row.values())
    for key,jet in owner.physical.items():assert len(current.outer.rows(jet))==(5 if key.startswith('Q') else 6)
    for key,rows in owner.own_reserved.items():
        assert len(rows)==5 and all(len(current.outer.rows(jet))==(3 if key in ('Q1','V1') else 4) for jet in rows)
    assert all(len(current.outer.rows(jet))==4 for jet in owner.higher_density.values())
    for chain in list(owner.compact_beta_chains)+list(owner.own_beta_chains.values()):
        for ref in chain[1:3]:
            row=owner.g.nodes[ref.node];assert row['operation']=='original_lazy_flat_branch'
            predicate=owner.g.nodes[row['flat_predicate']];assert predicate['operation']=='logical_or'
            actual={(owner.g.nodes[i]['operator'],owner.g.nodes[i]['left'],owner.g.nodes[i]['right']) for i in predicate['arguments']}
            assert actual=={('<=',row['Delta'],owner.g.constant(-1).node),('>=',row['Delta'],owner.g.one.node)}
    return pins


def profile_context(owner):
    g=owner.g;c=owner.c;Z=sy.Symbol('Z',real=True)
    X=sy.Symbol(g.nodes[owner.x.node]['name'],real=True)
    mu=sy.Symbol('mu',positive=True);P=sy.Symbol('actual_Pstar',positive=True)
    N=sy.Symbol('actual_selected_N',positive=True,integer=True);Rc=sy.Symbol('actual_compact_Rc',positive=True)
    delta=sy.Symbol('delta',positive=True);norm=sy.Symbol('actual_compact_beta_normalization',positive=True)
    bindings={g.symbol('Z').node:Z,owner.x.node:X,c['mu']:mu,c['Pstar']:P,c['N']:N,c['Rc']:Rc,
        c['delta']:delta,c['W']['normal']:norm}
    bind(bindings,owner.jet(c['A']),sy.Function('actual_compact_A')(Z),Z)
    for j,ids in enumerate(c['h']):bind(bindings,owner.jet(ids),sy.Function('actual_compact_control_'+str(j))(Z),Z)
    return dict(bindings=bindings,Z=Z,X=X,mu=mu,P=P,N=N,Rc=Rc,delta=delta,norm=norm)


def actual_profile_calculus(owner):
    pins=contract(owner);rec=high_check.admitted_receipt(owner,current.compact)
    assert rec['checks']['actual_source']['original_compact_log_transport_measure_and_relative_only_terminal_theorem']
    context=profile_context(owner);Z=context['Z'];Y=sy.Symbol('physical_compact_y',real=True)
    bindings=dict(context['bindings']);bindings[owner.x.node]=sy.exp(Y)
    alg=CompactAlgebra(owner,bindings);scalar=[];bumps=[];profiles=[];derived=[];flat=[]
    chains=list(owner.compact_beta_chains)+list(owner.own_beta_chains.values())
    for i in sorted({ref.node for chain in chains for ref in chain[3:]}):
        row=owner.g.nodes[i];n=row['derivative_order'];u=sy.Symbol('actual_beta_argument',real=True)
        source=CompactAlgebra(owner,{row['argument']:u})
        zero(source.at(i)-sy.diff(high_check.radial_check.Beta(u),u,n),('actual compact beta ordinary body',i,n))
        scalar.append(i)
        for boundary in (-2,-1,1,2):
            guard=CompactAlgebra(owner,{row['argument']:sy.Integer(boundary)})
            assert guard.at(i)==0 and i in guard.flat_nodes
    # Prove the low scalar chain bodies once as well. Then consume their
    # exact derivative relations in the logarithmic Leibniz rule. Repeated
    # expansion of the same rational beta bodies adds no source evidence.
    for chain in chains:
        for n,ref in enumerate(chain[1:3],1):
            row=owner.g.nodes[ref.node];u=sy.Symbol('actual_low_beta_argument',real=True)
            zero(CompactAlgebra(owner,{row['Delta']:u}).at(ref)-sy.diff(high_check.radial_check.Beta(u),u,n),
                ('actual low beta chain body',ref.node,n))
    symbolic_beta=sy.Function('actual_proved_beta');u=sy.Symbol('actual_beta_argument',real=True)
    scalar_bindings=dict(bindings)
    for chain in owner.compact_beta_chains:
        argument=sy.expand(alg.at(owner.g.nodes[chain[0].node]['argument']))
        for n,ref in enumerate(chain):scalar_bindings[ref.node]=sy.diff(symbolic_beta(u),u,n).subs(u,argument)
    scalar_alg=CompactAlgebra(owner,scalar_bindings)
    for j,packet in enumerate(owner.bumps):
        base=scalar_alg.at(packet[0])
        for n in (3,4):
            zero(scalar_alg.at(packet[n])-sy.diff(base,Y,n),('actual compact log-bump source',j,n));bumps.append('b'+str(j)+'_y'+str(n))
    profile_bindings=dict(bindings)
    for j,packet in enumerate(owner.bumps):
        B=sy.Function('actual_proved_compact_bump_'+str(j))(Y)
        for n,ref in enumerate(packet):profile_bindings[ref.node]=sy.diff(B,Y,n)
    profile_alg=CompactAlgebra(owner,profile_bindings)
    for n in (3,4):
        for key,jet in owner.profiles[n].items():
            base=profile_alg.at(owner.profiles[0][key].value)
            zero(profile_alg.at(jet.value)-sy.diff(base,Y,n),('actual compact profile source',n,key))
            profiles.append(key+'_y'+str(n)+'_Z0')
            derived.extend(key+'_y'+str(n)+'_Z'+str(k) for k in range(1,6))
    for value in (1,2):
        env=dict(context['bindings']);env[owner.x.node]=sy.Integer(value);at=CompactAlgebra(owner,env)
        for n,row in enumerate(owner.profiles):
            for key in ('F','G','delta_E','delta_V','V'):
                zero(at.at(row[key].value),('actual compact support endpoint',value,n,key));flat.append(str(value)+'_'+key+'_y'+str(n)+'_Z0')
    return dict(passed=True,actual_beta3_beta4_interior_body_nodes=scalar,actual_bump_y3_y4_rows=bumps,
        actual_whole_compact_profile_y3_y4_rows=profiles,higher_profile_Z_rows_by_ordinary_commutation=derived,
        actual_flat_endpoint_whole_profile_rows=flat,ordinary_jet_algebra_source_bindings=pins,
        actual_original_log_width_centers_normalization_and_selected_N_retained=True,
        beta_chain_relations_consumed_only_after_actual_body_proof=True,
        profile_bump_inputs_bound_only_after_actual_radial_source_proof=True,
        exact_guarded_ordinary_scalar_derivatives_not_Taylor_coefficients=True)


def actual_history_and_Q_calculus(owner):
    pins=contract(owner);g=owner.g;Z=sy.Symbol('Z',real=True);Y=sy.Symbol('physical_y',real=True)
    E,V=sy.Function('actual_E')(Y,Z),sy.Function('actual_V')(Y,Z)
    H={key:sy.Function('actual_'+key+'_history')(Y,Z) for key in current.compact.RATES}
    bindings={g.symbol('Z').node:Z}
    for j,row in enumerate(owner.profiles):
        for key,body in (('E',E),('V',V)):bind(bindings,row[key],sy.diff(body,Y,j),Z)
    for key,body in H.items():bind(bindings,owner.histories[0][key],body,Z)
    rhs=dict(m=V-H['m'],h=E-sy.Rational(3,2)*H['h'],k=E*V-sy.Rational(3,2)*H['k'],
        e=V**2-E**2/2-H['e'],p=E**2/2)
    expected=[H];replace={sy.diff(H[key],Y):body for key,body in rhs.items()}
    alg=CompactAlgebra(owner,bindings);actual=[];derived=[]
    for j in range(1,5):
        expected.append({key:sy.diff(body,Y).xreplace(replace) for key,body in expected[-1].items()})
        for key,jet in owner.histories[j].items():
            zero(alg.at(jet.value)-expected[j][key],('actual complete modified history ODE',key,j))
            actual.append(key+'_y'+str(j)+'_Z0');derived.extend(key+'_y'+str(j)+'_Z'+str(n) for n in range(1,6))
    P=sy.Symbol('actual_Pstar',positive=True);delta=sy.Symbol('delta',positive=True)
    M=sy.Function('actual_normalized_mean')(Y,Z);U=M+sy.diff(M,Y)
    bindings={g.symbol('Z').node:Z,owner.Pstar.node:P,owner.delta.node:delta}
    for j,row in enumerate(owner.profiles):bind(bindings,row['V'],sy.diff(U,Y,j),Z)
    bind(bindings,owner.histories[0]['m'],M,Z)
    alg=CompactAlgebra(owner,bindings)
    Q=P*(2*Z*U-(1-delta)*Z*M-(1-Z**2)*sy.diff(M,Z))/(1-delta*Z**2)
    qrows=[];qderived=[]
    for j,suffix in enumerate(('','_y','_yy','_yyy','_yyyy')):
        for prefix,body in (('U',P*U),('m',P*M),('Q',Q)):
            key=prefix+suffix
            zero(alg.at(owner.physical[key].value)-sy.diff(body,Y,j),('actual compact physical radial source',key))
            qrows.append(key+'_Z0');qderived.extend(key+'_Z'+str(n) for n in range(1,len(current.outer.rows(owner.physical[key]))))
    return dict(passed=True,actual_whole_complete_history_ODE_rows=actual,higher_history_Z_rows_by_ordinary_commutation=derived,
        actual_whole_physical_U_mass_Q_rows=qrows,higher_physical_Z_rows_by_ordinary_commutation=qderived,
        complete_original_histories_are_inputs_not_zero_inlets=True,original_pressure_rate_zero_kept=True,
        actual_Q_outputs_not_bound=True,Pstar_scaled_once=True,leading_one_minus_delta=True,
        ordinary_jet_algebra_source_bindings=pins)


def actual_compact_pressure_calculus(owner):
    pins=contract(owner);g=owner.g;Z=sy.Symbol('Z',real=True);Y=sy.Symbol('physical_y',real=True)
    delta=sy.Symbol('delta',positive=True);Rc=sy.Symbol('actual_Rc',positive=True)
    context=profile_context(owner);env=dict(context['bindings']);env[owner.x.node]=sy.exp(Y)
    R=CompactAlgebra(owner,env).at(owner.R)
    for j in (1,2):zero(sy.diff(R,Y,j)-R,('actual compact radius',j))
    U,Q=sy.Function('actual_U')(Y,Z),sy.Function('actual_Q')(Y,Z)
    bindings={g.symbol('Z').node:Z,owner.R.node:R,owner.delta.node:delta}
    for j,suffix in enumerate(('','_y','_yy','_yyy','_yyyy')):
        bind(bindings,owner.physical['Q'+suffix],sy.diff(Q,Y,j),Z)
        if j<3:bind(bindings,owner.physical['U'+suffix],sy.diff(U,Y,j),Z)
    alg=CompactAlgebra(owner,bindings);S=alg.at(owner.density['S'].value)
    Py=alg.at(owner.density['P1_y'].value);N=alg.at(owner.density['N1p'].value)
    expected=dict(S_y=sy.diff(S,Y),S_yy=sy.diff(S,Y,2),P1_yy=sy.diff(Py,Y),P1_yyy=sy.diff(Py,Y,2),
        N1p_y=sy.diff(N,Y),N1p_yy=sy.diff(N,Y,2),m15_yy=sy.diff(-R*Py,Y),m15_yyy=sy.diff(-R*Py,Y,2))
    actual=[];derived=[]
    for key,jet in owner.higher_density.items():
        zero(alg.at(jet.value)-expected[key],('actual compact pressure radial body',key))
        actual.append(key+'_Z0');derived.extend(key+'_Z'+str(j) for j in range(1,4))
    return dict(passed=True,actual_radius_y_and_yy_identity=True,actual_whole_pressure_radial_rows=actual,
        higher_pressure_Z_rows_by_ordinary_commutation=derived,actual_pressure_outputs_not_bound=True,
        same_original_signed_pressure_FTC_consumed=True,ordinary_jet_algebra_source_bindings=pins)


def own_context(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);Y=sy.Symbol(g.nodes[owner.phase.node]['name'],real=True)
    bindings=outer_check.radius_context(owner.own_outer,owner.high.joins.radius_proof)
    bindings.update({g.symbol('Z').node:Z,owner.phase.node:Y,owner.mu.node:sy.Symbol('mu',positive=True)})
    assert owner.coefficient_phase_independence==owner.own_outer.coefficient_phase_independence
    assert len(owner.coefficient_phase_independence)==20
    for key,jet in owner.own_outer.coefficients.items():bind(bindings,jet,sy.Function('actual_selected_'+key)(Z),Z)
    roots=[ref.node for key in ('Uz1','Utheta1','F1') for ref in current.outer.rows(owner.own_reserved[key][0])]
    normalizers=[];seen=set();stack=roots[:]
    while stack:
        i=stack.pop()
        if i in seen or i in bindings:continue
        seen.add(i);row=g.nodes[i]
        if row['operation']=='definite_integral':
            closure=owner.own_outer.operand_closure(i)
            assert owner.phase.node not in closure and g.symbol('Z').node not in closure
            bindings[i]=sy.Symbol('actual_own_fixed_normalizer_'+str(i),positive=True);normalizers.append(i)
        else:stack.extend(current.outer.source_operands(row))
    return dict(bindings=bindings,Z=Z,Y=Y,normalizers=sorted(set(normalizers)))


def actual_own_velocity_calculus(owner):
    pins=contract(owner);context=own_context(owner);Y=context['Y'];Z=context['Z'];bindings=context['bindings']
    # Consume the pinned scalar derivative source proved in the profile
    # group and accepted parent. Keep source guards for endpoint checks.
    original_bindings=dict(bindings);source=CompactAlgebra(owner,bindings)
    u=sy.Symbol('actual_own_beta_argument',real=True);beta=sy.Function('actual_proved_own_beta')
    for chain in owner.own_beta_chains.values():
        argument=sy.expand(source.at(owner.g.nodes[chain[0].node]['argument']))
        for n,ref in enumerate(chain):bindings[ref.node]=sy.diff(beta(u),u,n).subs(u,argument)
    alg=CompactAlgebra(owner,bindings);actual=[];derived=[];flat=[]
    R=alg.at(owner.own_R)
    for j in range(1,5):zero(sy.diff(R,Y,j)-R,('actual own radius derivative',j))
    for key in ('Uz1','Utheta1','F1'):
        base=alg.at(owner.own_reserved[key][0].value)
        for j in range(1,5):
            zero(alg.at(owner.own_reserved[key][j].value)-sy.diff(base,Y,j),('actual own selected bump radial source',key,j))
            actual.append(key+'_y'+str(j)+'_Z0');derived.extend(key+'_y'+str(j)+'_Z'+str(n) for n in range(1,4))
    zero(alg.at(owner.own_reserved['F1'][0].value)-alg.at(owner.own_reserved['Utheta1'][0].value)/sy.sqrt(2*R),
        'actual physical F1 versus Utheta1')
    for endpoint in (sy.Integer(0),sy.log(2)):
        env=dict(original_bindings);env[owner.phase.node]=endpoint;at=CompactAlgebra(owner,env)
        for key in ('Uz1','Utheta1','F1'):
            for j,jet in enumerate(owner.own_reserved[key]):
                zero(at.at(jet.value),('actual own bump support endpoint',endpoint,key,j));flat.append(str(endpoint)+'_'+key+'_y'+str(j)+'_Z0')
    return dict(passed=True,actual_whole_own_velocity_radial_rows=actual,higher_own_velocity_Z_rows_by_ordinary_commutation=derived,
        actual_flat_own_support_endpoint_rows=flat,actual_R_y_to_y4_equals_R=True,
        actual_F1_Utheta1_sqrt_2R_prefactor_retained=True,actual_selected_coefficient_phase_independence_rows=20,
        actual_fixed_normalizer_nodes=context['normalizers'],ordinary_jet_algebra_source_bindings=pins)


def actual_own_recovery_calculus(owner):
    pins=contract(owner);g=owner.g;Z=sy.Symbol('Z',real=True);Y=sy.Symbol('physical_own_y',real=True)
    R0=sy.Symbol('actual_own_origin_radius',positive=True);R=R0*sy.exp(Y);delta=sy.Symbol('delta',positive=True)
    M=sy.Function('actual_own_complete_mass')(Y,Z);U=sy.diff(M,Y)/R;mean=M/R
    bindings={g.symbol('Z').node:Z,owner.own_R.node:R,owner.delta.node:delta}
    assert current.outer.rows(owner.own_reserved['M1'][0])==current.outer.rows(owner.own_outer.reserved_fields['M1'])
    bind(bindings,owner.own_reserved['M1'][0],M,Z)
    for j,jet in enumerate(owner.own_reserved['Uz1']):bind(bindings,jet,sy.diff(U,Y,j),Z)
    alg=CompactAlgebra(owner,bindings);actual=[];derived=[]
    Q=(2*Z*U-(1+delta)*Z*mean-(1-Z**2)*sy.diff(mean,Z))/(1-delta*Z**2)
    for key,body in (('M1',M),('mean1',mean),('K1',mean-U),('Q1',Q),('V1',R*Q)):
        for j,jet in enumerate(owner.own_reserved[key]):
            if key=='M1' and j==0:continue # complete moment is the retained source input
            zero(alg.at(jet.value)-sy.diff(body,Y,j),('actual own radial recovery source',key,j))
            actual.append(key+'_y'+str(j)+'_Z0');derived.extend(key+'_y'+str(j)+'_Z'+str(n) for n in range(1,len(current.outer.rows(jet))))
    return dict(passed=True,actual_whole_own_mass_mean_K_Q_V_radial_rows=actual,
        higher_own_recovery_Z_rows_by_ordinary_commutation=derived,actual_own_Q_outputs_not_bound=True,
        source_complete_mass_and_actual_mass_FTC_retained=True,own_one_plus_delta=True,
        ordinary_jet_algebra_source_bindings=pins)


def mass_ode_reduce(expression,M,V,Y,Z):
    for _ in range(6):
        substitutions={}
        for derivative in expression.atoms(sy.Derivative):
            if derivative.expr!=M:continue
            counts=dict(derivative.variable_count);ny=counts.get(Y,0)
            if ny:substitutions[derivative]=sy.diff(V-M,Y,ny-1,Z,counts.get(Z,0))
        if not substitutions:return expression
        expression=expression.xreplace(substitutions)
    raise AssertionError('Unreduced actual mass ODE derivative')


def theta0_context(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);Y=sy.Symbol(g.nodes[owner.phase.node]['name'],real=True)
    theta=owner.theta0_radial[0]
    closure=owner.own_outer.operand_closure(theta.value.node)
    powers=[i for i in closure if g.nodes[i]['operation']=='analytic_unary' and g.nodes[i]['name']=='exp'
        and outer_check.free_query_dependency(owner.own_outer,i,owner.phase.node)]
    assert len(powers)==1;power=powers[0]
    bindings=outer_check.radius_context(owner.own_outer,owner.high.joins.radius_proof)
    bindings.update({g.symbol('Z').node:Z,owner.phase.node:Y,owner.mu.node:sy.Symbol('mu',positive=True)})
    amplitude=[]
    for n,ref in enumerate(current.outer.rows(theta)):
        row=g.nodes[ref.node];assert row['operation']=='product' and power in row['arguments']
        factors=[i for i in row['arguments'] if i!=power]
        assert factors and not any(outer_check.free_query_dependency(owner.own_outer,i,owner.phase.node) for i in factors)
        # One actual slow input coefficient per ordinary row. The shared
        # radial power is identified from the actual producer before this.
        coefficient=g.mul(*(current.source.FunctionRef(g,i) for i in factors))
        bindings[coefficient.node]=sy.diff(sy.Function('actual_Astar')(Z),Z,n);amplitude.append(coefficient.node)
    alg=CompactAlgebra(owner,bindings);real_power=alg.at(power);mu=bindings[owner.mu.node]
    zero(sy.diff(real_power,Y)+(sy.Rational(1,2)+mu)*real_power,'actual theta0 radial power source')
    return dict(bindings=bindings,Z=Z,Y=Y,actual_power_node=power,actual_amplitude_nodes=amplitude)


def actual_total_pressure_source(owner):
    pins=contract(owner);g=owner.g;repaired_mod=current.outer.intermediate.previous
    rec=high_check.admitted_receipt(owner,repaired_mod)
    assert rec[repaired_mod.GATES[-1]]
    context=profile_context(owner);Z=context['Z'];Y=sy.Symbol(g.nodes[owner.phase.node]['name'],real=True)
    V=sy.Function('actual_compact_V')(Y,Z);M=sy.Function('actual_compact_m')(Y,Z)
    bindings=dict(context['bindings']);bindings.update({owner.x.node:sy.exp(Y),owner.phase.node:Y})
    for j,row in enumerate(owner.profiles):bind(bindings,row['V'],sy.diff(V,Y,j),Z)
    bind(bindings,owner.histories[0]['m'],M,Z)
    # These admitted low-order compact source recipes are retained rather
    # than binding either old or new pressure output to a common symbol.
    bind(bindings,owner.jet(owner.c['radial']['m_y']),V-M,Z)
    bind(bindings,owner.jet(owner.c['radial']['m_yy']),sy.diff(V,Y)-V+M,Z)
    alg=CompactAlgebra(owner,bindings)
    old=alg.at(owner.repaired['fixed_pressure_y'][0]);new=alg.at(owner.density['P1_y'].value)
    zero(new-old,'actual compact-to-own fixed pressure-y source projection')
    actual=['fixed_P1_y_Z0'];derived=['fixed_P1_y_Z'+str(n) for n in range(1,4)]
    for n,key in ((1,'P1_yy'),(2,'P1_yyy')):
        want=mass_ode_reduce(sy.diff(old,Y,n),M,V,Y,Z)
        zero(alg.at(owner.higher_density[key].value)-want,('actual higher fixed pressure source',key))
        actual.append(key+'_Z0');derived.extend(key+'_Z'+str(j) for j in range(1,4))
    fixed_bindings=dict(bindings)
    th=theta0_context(owner);bindings=dict(th['bindings']);Y=th['Y'];Z=th['Z']
    Theta=sy.Function('actual_own_Utheta')(Y,Z)
    for j,jet in enumerate(owner.own_reserved['Utheta1']):bind(bindings,jet,sy.diff(Theta,Y,j),Z)
    alg=CompactAlgebra(owner,bindings);base=alg.at(owner.patch_pressure_density[0].value)
    theta0=alg.at(owner.theta0_radial[0].value)
    zero(base-theta0*Theta,'actual patch pressure density theta0 Utheta1')
    patch=[]
    for j,jet in enumerate(owner.patch_pressure_density):
        zero(alg.at(jet.value)-sy.diff(base,Y,j),('actual own pressure patch radial source',j));patch.append('patch_m13_density_y'+str(j)+'_Z0')
    # The patch-only negative m15 density is compared with patch m13.
    # The full signed derivative also includes the leading fixed pressure.
    geometry=outer_check.radius_context(owner.own_outer,owner.high.joins.radius_proof)
    geometry.update({owner.x.node:sy.exp(Y),owner.phase.node:Y})
    geometry_alg=CompactAlgebra(owner,geometry)
    zero(geometry_alg.at(owner.R)-geometry_alg.at(owner.own_R),'actual whole reserved/compact physical radius source')
    R=context['Rc']*sy.exp(Y)
    full=dict(fixed_bindings);full.update(th['bindings']);full[owner.own_R.node]=R
    for j,jet in enumerate(owner.own_reserved['Utheta1']):bind(full,jet,sy.diff(Theta,Y,j),Z)
    alg=CompactAlgebra(owner,full)
    patch_m15=alg.at(owner.repaired['patch_bodies']['per_ds']['m15'][0])
    zero(patch_m15+R*alg.at(owner.patch_pressure_density[0].value),'actual negative patch m15 sign')
    # Check the actual total pressure expressions, keeping both inputs.
    # Neither old/new P1_y nor new P1_yy/P1_yyy outputs are bound.
    total=alg.at(owner.own_pressure['P1_y'].value);total_rows=[]
    zero(total-alg.at(owner.repaired['fixed_pressure_y'][0])-alg.at(owner.patch_pressure_density[0].value),
        'actual retained leading-plus-patch pressure FTC')
    for key,n in (('P1_yy',1),('P1_yyy',2)):
        zero(alg.at(owner.own_pressure[key].value)-mass_ode_reduce(sy.diff(total,Y,n),M,V,Y,Z),
            ('actual own total pressure radial body',key));total_rows.append(key+'_Z0')
    weighted=[]
    for key,n in (('m15_yy',1),('m15_yyy',2)):
        zero(alg.at(owner.own_m15_higher[key].value)-mass_ode_reduce(sy.diff(-R*total,Y,n),M,V,Y,Z),
            ('actual full weighted pressure derivative',key));weighted.append(key+'_Z0')
    return dict(passed=True,actual_whole_fixed_pressure_source_rows=actual,higher_fixed_pressure_Z_rows_by_ordinary_commutation=derived,
        actual_whole_patch_pressure_density_rows=patch,actual_total_own_pressure_radial_rows=total_rows,
        actual_weighted_total_pressure_rows=weighted,
        theta0_actual_power_node=th['actual_power_node'],theta0_actual_query_independent_amplitude_nodes=th['actual_amplitude_nodes'],
        old_and_new_fixed_pressure_outputs_not_bound=True,
        actual_total_pressure_outputs_not_bound=True,actual_whole_compact_reserved_radius_source_identity=True,
        patch_m15_negative_sign_and_full_leading_plus_patch_FTC_retained=True,
        full_nonzero_P1_and_m15_value_anchors_retained=True,
        admitted_actual_total_pressure_FTC_receipt=repaired_mod.RECEIPT,
        ordinary_jet_algebra_source_bindings=pins)


def actual_outside_own_radial_source(owner):
    pins=contract(owner);g=owner.g;Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True)
    Rp=sy.Symbol('actual_Rp',positive=True);M=sy.Function('actual_retained_own_mass')(Z)
    actual=[];derived=[];independent=[]
    for name,fields in owner.own_outside.items():
        if name=='quiet':
            data=owner.own_outer.charts['quiet'];coordinate=data['coordinate'];jac=sy.Integer(1)
        else:
            data=owner.high.charts[name];coordinate=data['coordinate'];jac=1 if name=='end' else mu
        X=sy.Symbol(g.nodes[coordinate.node]['name'],real=True)
        bindings={g.symbol('Z').node:Z,owner.mu.node:mu,owner.delta.node:sy.Symbol('delta',positive=True),
            owner.high.inlet.geometry['Rp'].node:Rp,coordinate.node:X}
        for ref in current.outer.rows(fields['M1'][0]):
            assert not outer_check.free_query_dependency(owner.own_outer,ref.node,coordinate.node);independent.append(name+'_M1_'+str(ref.node))
        bind(bindings,fields['M1'][0],M,Z)
        alg=CompactAlgebra(owner,bindings)
        for key,rows in fields.items():
            base=alg.at(rows[0].value)
            for j,jet in enumerate(rows[1:],1):
                zero(alg.at(jet.value)-jac**j*sy.diff(base,X,j),('actual flat outside own radial source',name,key,j))
                actual.append(name+'_'+key+'_y'+str(j)+'_Z0')
                derived.extend(name+'_'+key+'_y'+str(j)+'_Z'+str(n) for n in range(1,len(current.outer.rows(jet))))
    return dict(passed=True,actual_whole_outside_own_y1_to_y4_source_rows=actual,
        higher_outside_own_Z_rows_by_ordinary_commutation=derived,
        actual_saved_M1_radial_query_independence_rows=independent,
        actual_nonzero_saved_mass_and_existing_pressure_values_retained=True,
        actual_physical_main_gap_mu_end_quiet_one_Jacobian=True,
        ordinary_jet_algebra_source_bindings=pins)


def actual_compact_pressure_FTC_source(owner):
    pins=contract(owner);g=owner.g;c=owner.c;context=profile_context(owner)
    integral_sources={key:current.parent.active.current.ast_binding(getattr(current.outer.Algebra,key))
        for key in ('at','integral')}
    integral_sources['signed_pressure_geometry']=current.parent.active.current.ast_binding(current.compact.previous.pressure_geometry)
    for pin in integral_sources.values():assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    Z,X=context['Z'],context['X'];bindings=dict(context['bindings'])
    V=sy.Function('actual_compact_V')(X,Z);M=sy.Function('actual_compact_m')(X,Z)
    bind(bindings,owner.profiles[0]['V'],V,Z)
    bind(bindings,owner.profiles[1]['V'],sy.Function('actual_compact_Vy')(X,Z),Z)
    bind(bindings,owner.histories[0]['m'],M,Z)
    bind(bindings,owner.jet(c['radial']['m_y']),V-M,Z)
    alg=CompactAlgebra(owner,bindings);S=alg.at(c['S'][0]);R=alg.at(owner.R)
    recipes=[]
    for key,power in (('S0',1),('S1',2)):
        for n,i in enumerate(c['pressure_local'][key]):
            row=g.nodes[i]
            assert row['operation']=='definite_integral' and row['lower']==g.one.node and row['upper']==owner.x.node
            assert row['physical_radius_power']==power and row['physical_Jacobian_applied_once']
            assert row['nonzero_preceding_S0_S1_retained'] and row['endpoint_independent_of_Z']
            assert row['ordinary_slow_Z_derivative_order']==n
            if n==0:
                t=sy.Symbol(row['variable'],real=True)
                zero(alg.at(row['integrand'])-(R**power*S/X).subs(X,t),('actual compact pressure integral body',key))
            recipes.append(key+'_Z'+str(n))
    Y=sy.Symbol('physical_pressure_y',real=True);V=sy.Function('actual_V')(Y,Z);M=sy.Function('actual_m')(Y,Z)
    bindings=dict(context['bindings']);bindings[owner.x.node]=sy.exp(Y)
    for j,row in enumerate(owner.profiles):bind(bindings,row['V'],sy.diff(V,Y,j),Z)
    bind(bindings,owner.histories[0]['m'],M,Z)
    bind(bindings,owner.jet(c['radial']['m_y']),V-M,Z)
    bind(bindings,owner.jet(c['radial']['m_yy']),sy.diff(V,Y)-V+M,Z)
    J={key:sy.Function('actual_same_partial_'+key)(Y,Z) for key in ('S0','S1')}
    for key,body in J.items():bind(bindings,owner.jet(c['pressure_local'][key]),body,Z)
    retained=[]
    for label,packet in (('incoming',c['incoming_integrals']),('initial',c['initial'])):
        for key,ids in packet.items():
            for i in ids:assert not outer_check.free_query_dependency(owner.own_outer,i,owner.x.node)
            bind(bindings,owner.jet(ids),sy.Function('actual_retained_'+label+'_'+key)(Z),Z);retained.append(label+'_'+key)
    alg=CompactAlgebra(owner,bindings);S=alg.at(c['S'][0]);R=alg.at(owner.R);Py=alg.at(owner.density['P1_y'].value)
    ftc={sy.diff(J['S0'],Y):R*S,sy.diff(J['S1'],Y):R**2*S}
    actual=[];derived=[]
    for key,want in (('m13',Py),('m15',-R*Py)):
        pressure=alg.at(c['pressure'][key][0])
        derivative=mass_ode_reduce(sy.diff(pressure,Y).xreplace(ftc),M,V,Y,Z)
        zero(derivative-want,('actual retained compact signed pressure FTC',key));actual.append(key+'_y_Z0')
        derived.extend(key+'_y_Z'+str(n) for n in range(1,4))
    return dict(passed=True,actual_same_defining_pressure_integral_recipe_rows=recipes,
        actual_retained_signed_pressure_FTC_whole_rows=actual,higher_pressure_FTC_Z_rows_by_ordinary_commutation=derived,
        retained_query_independent_nonzero_input_packets=retained,
        actual_input_integrals_bound_only_after_same_body_measure_and_endpoint_proof=True,
        actual_pressure_outputs_not_bound=True,original_measure_dx_over_x_applied_once=True,
        ordinary_integral_constructor_and_substitution_source_bindings=integral_sources,
        higher_integral_rows_from_same_typed_body_and_named_binder=True,
        ordinary_jet_algebra_source_bindings=pins)


def actual_2Rc_source_joins(owner):
    pins=contract(owner);g=owner.g;c=owner.c;q=owner.high.quiet['source']
    construction=current.parent.linear.read(current.outer.CONSTRUCTION_RECEIPT)
    assert construction['all_three_ordinary_axial_trace_identities_passed']
    assert construction['source_family']==owner.family
    assert construction['ordinary_own_field_and_cumulative_trace_rows']==204
    assert q['incoming']==c['terminal'] and q['R2']==c['terminal']['radius']
    assert not outer_check.free_query_dependency(owner.own_outer,owner.own_R.node,g.symbol('Z').node)
    context=profile_context(owner);Z=context['Z'];X=context['X'];env=dict(context['bindings'])
    env.update({owner.x.node:sy.Integer(2),q['s']:sy.Integer(0)})
    for key,jet in owner.histories[0].items():bind(env,jet,sy.Function('actual_nonzero_terminal_'+key)(X,Z),Z)
    alg=CompactAlgebra(owner,env);whole=[];derived=[]
    for j in range(5):
        for key in ('E','V'):
            left=owner.compact_endpoint['profiles'][j][key]
            right=owner.jet(q['flow']['profiles_y0_to_y4'][j][key])
            zero(alg.at(left.value)-alg.at(right.value),('actual compact/quiet radial profile',key,j))
            whole.append(key+'_y'+str(j)+'_Z0');derived.extend(key+'_y'+str(j)+'_Z'+str(n) for n in range(1,6))
        for key in current.compact.RATES:
            left=owner.compact_endpoint['histories'][j][key]
            right=owner.jet(q['flow']['histories_y0_to_y4'][j][key])
            zero(alg.at(left.value)-alg.at(right.value),('actual compact/quiet complete history',key,j))
            whole.append(key+'_history_y'+str(j)+'_Z0');derived.extend(key+'_history_y'+str(j)+'_Z'+str(n) for n in range(1,6))
    flat_env=own_context(owner)['bindings'];flat_env.update(theta0_context(owner)['bindings'])
    flat_env[owner.phase.node]=sy.log(2);flat_alg=CompactAlgebra(owner,flat_env);flat_patch=[]
    for j,jet in enumerate(owner.patch_pressure_density):
        zero(flat_alg.at(jet.value),('actual2Rc patch pressure density flat source',j));flat_patch.append('patch_m13_density_y'+str(j)+'_Z0')
    # Actual native linear-Z background mass remains nonzero. Use the
    # pinned quiet theorem, with both physical endpoint radii preserved.
    native=high_check.admitted_receipt(owner,current.rp.quiet)
    assert native['checks']['native_mass_identity']['actual_compact_derived_K_equals_original_Rp_Pstar_kappa']
    R2=sy.Symbol('actual_R2',positive=True);kappa=sy.Symbol('actual_normalized_terminal_mass',nonzero=True,real=True)
    env=dict(context['bindings']);env.update({owner.x.node:sy.Integer(2),owner.phase.node:sy.log(2),
        c['Rc']:R2/2,q['R2']:R2,q['s']:sy.Integer(0),owner.own_R.node:R2})
    for j,row in enumerate(owner.profiles):bind(env,row['V'],sy.Integer(0),Z)
    bind(env,owner.histories[0]['m'],kappa*Z,Z)
    bind(env,owner.jet(c['radial']['m_y']),-kappa*Z,Z)
    bind(env,owner.jet(c['radial']['m_yy']),kappa*Z,Z)
    alg=CompactAlgebra(owner,env);leading=[]
    for key in ('P1_yy','P1_yyy','N1p_y','N1p_yy','m15_yy','m15_yyy'):
        right=current.outer.Algebra(g,3).at(owner.high.quiet['higher_density'][key],owner.high.quiet['coordinate'],g.zero)
        for n,(left_ref,right_ref) in enumerate(zip(current.outer.rows(owner.compact_endpoint['higher_density'][key]),current.outer.rows(right))):
            zero(alg.at(left_ref)-alg.at(right_ref),('actual compact/quiet higher pressure',key,n));leading.append(key+'_Z'+str(n))
    # Both actual own mass traces are the accepted nonzero source function
    # at2Rc. Only endpoint input values are abstracted, never an interior
    # cumulative history or a pressure output.
    M=sy.Function('actual_own_endpoint_mass')(Z)
    for key in ('Uz1','Utheta1','F1'):
        for jet in owner.own_reserved[key]:bind(env,jet,sy.Integer(0),Z)
    bind(env,owner.own_reserved['M1'][0],M,Z)
    bind(env,owner.own_outside['quiet']['M1'][0],M,Z)
    for jet in (owner.own_reserved['M1'][0],owner.own_outside['quiet']['M1'][0]):
        assert len(current.outer.rows(jet))==4
    alg=CompactAlgebra(owner,env);own=[]
    for key,left_rows in owner.own_2Rc_traces['left'].items():
        right_rows=owner.own_2Rc_traces['right'][key]
        for j,(left,right) in enumerate(zip(left_rows,right_rows)):
            for n,(lr,rr) in enumerate(zip(current.outer.rows(left),current.outer.rows(right))):
                zero(alg.at(lr)-alg.at(rr),('actual own2Rc radial source join',key,j,n));own.append(key+'_y'+str(j)+'_Z'+str(n))
    bind(env,owner.theta0_radial[0],sy.Function('actual_nonzero_theta0_trace')(Z),Z)
    alg=CompactAlgebra(owner,env);pressure=[]
    for packet in (owner.pressure_2Rc_traces,owner.m15_2Rc_traces):
        for key,left in packet['left'].items():
            right=packet['right'][key]
            for n,(lr,rr) in enumerate(zip(current.outer.rows(left),current.outer.rows(right))):
                zero(alg.at(lr)-alg.at(rr),('actual total own2Rc higher pressure join',key,n));pressure.append(key+'_Z'+str(n))
    return dict(passed=True,actual_whole_compact_quiet_profile_history_y0_to_y4_rows=whole,
        higher_compact_quiet_Z_rows_by_ordinary_commutation=derived,
        actual_leading_2Rc_higher_pressure_C3_rows=leading,actual_own_2Rc_y0_to_y4_C3_QV_C2_rows=own,
        actual_total_own_2Rc_pressure_and_m15_C3_rows=pressure,
        actual_patch_pressure_density_y0_y1_y2_flat_whole_rows=flat_patch,
        actual_low_trace_and_radius_source_receipt=current.outer.CONSTRUCTION_RECEIPT,
        actual_native_nonzero_linear_mass_theorem=current.rp.quiet.RECEIPT,
        actual_original_low_pressure_values_and_nonzero_incoming_histories_retained=True,
        actual_compact_and_reserved_velocity_flat_source_lemmas_consumed_before_input_abstraction=True,
        actual_leading_one_minus_delta_and_own_one_plus_delta_kept_separate=True,
        own_radius_Z_independence_checked=True,pressure_outputs_not_bound=True,
        ordinary_jet_algebra_source_bindings=pins)


CHECKS={'actual_profile_calculus','actual_history_and_Q_calculus','actual_compact_pressure_calculus',
    'actual_own_velocity_calculus','actual_own_recovery_calculus','actual_compact_pressure_FTC_source',
    'actual_total_pressure_source','actual_outside_own_radial_source','actual_2Rc_source_joins'}


def save(owner,checks,execution_seconds):
    contract(owner);assert set(checks)==CHECKS and all(proof['passed'] for proof in checks.values())
    raw=current.parent.linear.read(current.NAME)
    assert raw['compact_higher_radial_pressure_C3_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    counts=checks['actual_2Rc_source_joins']
    assert len(counts['actual_leading_2Rc_higher_pressure_C3_rows'])==24
    assert len(counts['actual_own_2Rc_y0_to_y4_C3_QV_C2_rows'])==150
    assert len(counts['actual_total_own_2Rc_pressure_and_m15_C3_rows'])==16
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        leading_2Rc_higher_pressure_C3_rows=24,own_2Rc_y4_C3_QV_C2_rows=150,
        total_own_2Rc_higher_pressure_C3_rows=16,exact_graph_nodes=len(raw['exact_graph_nodes']),
        numerical_oracle_available=False,execution_seconds=execution_seconds)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ACTUAL_COMPACT_HIGHER_RADIAL_AND_OWN_2RC_SOURCE_JOINS_CHECKED',flush=True)
    return rec


def run(owner):
    began=time.monotonic();checks={}
    for fn in (actual_profile_calculus,actual_history_and_Q_calculus,actual_compact_pressure_calculus,
        actual_own_velocity_calculus,actual_own_recovery_calculus,actual_compact_pressure_FTC_source,
        actual_total_pressure_source,actual_outside_own_radial_source,actual_2Rc_source_joins):
        checks[fn.__name__]=fn(owner);print('COMPACT_HIGHER_'+fn.__name__.upper()+'_PASSED',flush=True)
    return save(owner,checks,time.monotonic()-began)


if __name__=='__main__':
    import sys
    sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
    run(current.CurrentOriginalCompactHigherRadialPressureC3Functions(require_checked=False))
