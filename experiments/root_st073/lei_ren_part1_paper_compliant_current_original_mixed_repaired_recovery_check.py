"""Independent compact analytic, mixed recovery, physical-rate and cap checks."""
import gzip
import json
import math
from pathlib import Path
import time
from types import FunctionType,SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_mixed_repaired_recovery as current
from lei_ren_part1_paper_compliant_current_original_C2_repaired_recovery_check import symbolic


def bind_jet(bindings,q,expression,z,order=2):
    for j,key in enumerate(('value','Z','ZZ','ZZZ')[:order+1]):
        bindings[getattr(q,key).node]=s.diff(expression,z,j)


def compact_profiles(field):
    g=field.graph;functions=field.functions;profiles=functions['actual_C3_profiles_y4']
    t,y,z=s.symbols('t y Z');ell,J=s.symbols('ell J',positive=True);centers=s.symbols('c0:3')
    beta=s.exp(-1/(1-t*t));polynomials=current.phase.flat_source.BETA_POLYNOMIALS
    for n,row in enumerate(polynomials):
        polynomial=sum(v*t**k for k,v in enumerate(row))
        assert s.cancel(s.diff(beta,t,n)/beta-polynomial/(1-t*t)**(2*n))==0,n
        assert current.phase.flat_source.BETA_CONSTANTS[n]==sum(abs(v) for v in row)
    assert current.phase.flat_source.BETA_CONSTANTS[3:5]==(88,1096)
    built=field.lower.band.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter'];W=built['repair_weights']
    x=field.lower.band.functions['original_band_variable'];B=s.Function('B')
    bindings={x.node:s.exp(y),g.unary('log',x).node:y,W['ell'].node:ell,W['normal'].node:J}
    bsymbols=s.symbols('beta0:5');bump_count=0
    for i,center in enumerate(W['centers']):
        bindings[center.node]=centers[i]
        for k,q in enumerate(profiles['raw_beta_ordinary_rows'][i]):
            node=g.nodes[q.node]
            assert node['ordinary_beta_derivative_order']==k and node['polynomial_coefficients']==polynomials[k]
            assert node['Taylor_coefficient_factorial']==math.factorial(k)
            assert node['ordinary_function_not_Taylor_coefficient'] and node['lazy_outside_support']
            assert node['all_support_endpoint_jets_zero'] and not node['derivative_of_range_endpoint']
            assert node['original_derivative_source']==current.current.ast_binding(current.phase.flat_source.beta_jets)
            assert node['original_polynomial_source']==current.current.ast_binding(current.phase.flat_source.beta_polynomials)
            bindings[q.node]=bsymbols[k]
        arg=(y-centers[i])/ell
        for j in (3,4):
            expected=s.diff(s.exp(-y)*B(arg)/(ell*J),y,j)
            replacements={q:bsymbols[sum(v for _,v in q.expr.variable_count)] for q in expected.atoms(s.Subs)}
            expected=expected.xreplace(replacements).replace(lambda expr:expr.func==B,lambda expr:bsymbols[0])
            got=symbolic(g,profiles['bump_ordinary_y_rows'][j][i].node,bindings)
            assert s.cancel(got-expected)==0,(i,j)
            bump_count+=1
    # Verify all new profile Z rows with independent arbitrary amplitude/control functions.
    A=s.Function('A')(z);controls=[s.Function('h'+str(i))(z) for i in range(5)]
    amplitude=field.lower.band.control.ranges.target.functions['actual_terminal_amplitude']
    bind_jet(bindings,amplitude,A,z,3)
    for q,expr in zip(field.lower.band.control.control_functions(),controls):bind_jet(bindings,q,expr,z,3)
    N,alpha,power=s.symbols('N alpha power',positive=True)
    bindings[built['N'].node]=N;bindings[profiles['original_alpha'].node]=alpha;bindings[profiles['original_power'].node]=power
    count=0
    for j in (3,4):
        b=s.symbols('b'+str(j)+'_0:3')
        for q,expr in zip(profiles['bump_ordinary_y_rows'][j],b):bindings[q.node]=expr
        F=(b[0]*controls[2]+b[1]*controls[3]+b[2]*controls[4])/N
        G=(b[0]*controls[0]+b[2]*controls[1])/N
        original=A*power*(-alpha)**j
        expected=dict(F=F,G=G,original_E=original,delta_E=A*F,delta_V=A*G,E=original+A*F,V=A*G)
        for key,q in profiles['radial_y_rows'][j].items():
            for order,row in enumerate(current.target.rows(q)):
                assert s.cancel(symbolic(g,row.node,bindings)-s.diff(expected[key],z,order))==0,(j,key,order)
                count+=1
    # exp(-1/w) times every finite inverse power tends to zero at w=0.
    w=s.Symbol('w',positive=True)
    assert s.limit(s.exp(-1/w)/w**8,w,0,dir='+')==0
    return dict(independent_ordinary_beta_derivative_orders=list(range(5)),
        independent_y3_y4_exponential_width_normalization_checks=bump_count,
        independent_new_profile_Z0_Z1_Z2_Z3_rows=count,
        compact_endpoint_flatness_checked=True,ordinary_beta_not_Taylor_rows=True)


def signed_histories(field):
    y,z=s.symbols('y Z');E=s.Function('E')(y,z);V=s.Function('V')(y,z)
    M={key:s.Function(key)(y,z) for key in current.current.RATES}
    density=dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
    rates={key:s.Rational(str(rate)) for key,rate in current.current.RATES.items()}
    rules={s.diff(M[key],y):density[key]-rates[key]*M[key] for key in M}
    expected=[M]
    for j in range(4):expected.append({key:s.diff(expr,y).subs(rules) for key,expr in expected[-1].items()})
    bindings={};f=field.functions;g=field.graph
    for j,radial in enumerate(f['actual_C3_profiles_y4']['radial_y_rows']):
        for key,expr in (('E',E),('V',V)):bind_jet(bindings,radial[key],s.diff(expr,y,j),z,3)
    for j in (0,1):
        for key,q in f['actual_C3_complete_history_y_rows'][j].items():bind_jet(bindings,q,expected[j][key],z,3)
    count=0
    for j in (2,3,4):
        for key,q in f['actual_C3_complete_history_y_rows'][j].items():
            for order,row in enumerate(current.target.rows(q)):
                assert s.expand(symbolic(g,row.node,bindings)-s.diff(expected[j][key],z,order))==0,(j,key,order)
                count+=1
    for q in f['actual_C3_independent_P0_y_rows'][1:]:assert current.target.rows(q)==(g.zero,)*4
    return dict(independent_repeated_signed_own_rate_FTC_Z3_rows=count,
        original_nonzero_complete_histories_retained=True,independent_P0_positive_y_exact_zero=True)


def polynomials():
    y,z=s.symbols('y Z');E=2+y/3+y*y/7+y**3/11+y**4/13+z/5+y*z*z/17+z**3/19
    V=1-y/5+y*y/9-y**3/15+y**4/21+z/7-y*z*z/23+z**3/29
    H={key:(i+1)/s.Integer(3)+y/(i+4)+y*y/(i+9)+y**3/(i+13)+y**4/(i+17)
        +z/(i+5)+y*z*z/(i+19)+z**3/(i+23) for i,key in enumerate(current.current.RATES)}
    P0=s.Rational(9,10)+3*z/10-z*z/5+z**3/31
    return y,z,E,V,H,P0,s.Rational(1,5),7*s.exp(y),s.Integer(3)


def independent_expected(y,z,E,V,H,P0,delta,R,S):
    d=1-z*z;L=1-delta*z*z;P=P0+H['p']
    m,h,k,e=(H[key] for key in ('m','h','k','e'))
    transport=(1-delta)*z*m+d*s.diff(m,z);Q=(2*z*V-transport)/L
    sectors=dict(theta_linear=(-E+(1-delta/2)*h-(1-delta)*z*s.diff(h,z)/2)/L,
        theta_quadratic=((2*delta-1)*z*k-d*s.diff(k,z)+E*transport)/L,
        axial_linear=(-V+(1-delta)*(m-z*s.diff(m,z))/2)/L,
        axial_quadratic=(V*transport+2*delta*z*e-d*s.diff(e,z)+2*(1+delta)*z*P-d*s.diff(P,z))/L)
    C=E-2*s.diff(E,y);B=2*s.diff(V,y)
    prefactor=S*s.sqrt(R/2);shear_prefactor=S/s.sqrt(2*R)
    derivatives=lambda expr,n:[s.diff(expr,y,j) for j in range(n)]
    return dict(common_radial_Q_y_rows=derivatives(Q,5),common_absolute_pressure_y_rows=derivatives(P,5),
        common_absolute_pressure_y_Z_rows=derivatives(s.diff(P,z),5),
        full_signed_inertial_sector_y_rows={key:derivatives(expr,4) for key,expr in sectors.items()},
        actual_generic_C_B_y_rows=dict(C=derivatives(C,4),B=derivatives(B,4)),
        actual_generic_inertial_y_rows={key:derivatives(R*(sectors[key+'_linear']+S*sectors[key+'_quadratic']),4)
            for key in ('theta','axial')},
        original_cylindrical_velocity_pressure_y_rows={key:derivatives(expr,5)
            for key,expr in dict(Utheta=S*E,Uz=S*V,Ur=prefactor*Q,Pi=S*S*P).items()},
        original_physical_signed_stress_y_rows={key:derivatives(expr,4) for key,expr in
            {**{key:prefactor*(S if 'quadratic' in key else 1)*expr for key,expr in sectors.items()},
                'shear_theta':-shear_prefactor*C,'shear_axial':shear_prefactor*B}.items()})


def walk_pairs(got,expected):
    if isinstance(got,dict):
        assert got.keys()==expected.keys()
        for key in got:yield from walk_pairs(got[key],expected[key])
    elif isinstance(got,list):
        assert len(got)==len(expected)
        for a,b in zip(got,expected):yield from walk_pairs(a,b)
    else:yield got,expected


def recovery_calculus(field):
    y,z,E,V,H,P0,delta,R,S=polynomials();f=field.functions;bindings={};g=field.graph
    for j,radial in enumerate(f['actual_C2_profile_inputs']):
        for key,expr in (('E',E),('V',V)):bind_jet(bindings,radial[key],s.diff(expr,y,j),z)
    for j,rows in enumerate(f['actual_C2_history_inputs']):
        for key,q in rows.items():bind_jet(bindings,q,s.diff(H[key],y,j),z)
    for j,rows in enumerate(f['actual_C2_history_Z_inputs']):
        for key,q in rows.items():bind_jet(bindings,q,s.diff(H[key],y,j,z),z)
    P=f['actual_C3_independent_P0_y_rows'][0]
    bind_jet(bindings,P,P0,z,3)
    parameters=f['actual_source_parameters']
    bindings[parameters['Z'].value.node]=z;bindings[parameters['delta'].value.node]=delta
    bindings[parameters['R'].value.node]=R;bindings[parameters['Pstar'].value.node]=S
    expected=independent_expected(y,z,E,V,H,P0,delta,R,S);count=0
    for q,expr in walk_pairs(f['recovered_mixed_C2_functions'],expected):
        for j,key in enumerate(('value','Z','ZZ')):
            got=symbolic(g,getattr(q,key).node,bindings)
            assert s.cancel(got-s.diff(expr,z,j))==0,(count,j)
            count+=1
    assert count==273
    return dict(independent_mixed_ordinary_y_Z_calculus_rows=count,
        generic_radius_rate_one_and_physical_half_rates_checked=True,
        independent_signed_meridional_pressure_terms_checked=True)


def directed_fixture(field):
    y,z,E,V,H,P0,delta,R,S=polynomials();c=field.c;bd=current.lower.ranges.Bounds(c)
    point={y:s.Rational(1,7),z:s.Rational(1,4)}
    def evaluate(expr):
        if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
        if expr.func==s.Add:return sum((evaluate(q) for q in expr.args),c.mpf(0))
        if expr.func==s.Mul:return math.prod(evaluate(q) for q in expr.args)
        if expr.func==s.Pow:return c.power(evaluate(expr.args[0]),evaluate(expr.args[1]))
        if expr.func==s.exp:return c.exp(evaluate(expr.args[0]))
        raise ValueError('Independent exact fixture expression '+str(expr))
    scalar=lambda expr:evaluate(expr.subs(point))
    jet=lambda expr:current.lower.ranges.JetBound(*(bd.constant(scalar(s.diff(expr,z,j))) for j in range(3)))
    radial=[dict(E=jet(s.diff(E,y,j)),V=jet(s.diff(V,y,j))) for j in range(5)]
    M=[{key:jet(s.diff(expr,y,j)) for key,expr in H.items()} for j in range(5)]
    MZ=[{key:jet(s.diff(expr,y,j,z)) for key,expr in H.items()} for j in range(5)]
    parameters=dict(Z=jet(z),d=jet(1-z*z),L=jet(1-delta*z*z),delta=jet(delta),R=jet(R),Pstar=jet(S),
        Pstar_sqrt_R_over_2=jet(S*s.sqrt(R/2)),Pstar_over_sqrt_2R=jet(S/s.sqrt(2*R)))
    caps=current.mixed_recover(current.lower.MagnitudeAlgebra(c,{'L':c.ln(c.mpf('.8'))}),radial,M,MZ,jet(P0),jet(s.diff(P0,z)),parameters)
    expected=independent_expected(y,z,E,V,H,P0,delta,R,S);count=0
    for q,expr in walk_pairs(caps,expected):
        for j,key in enumerate(('value','Z','ZZ')):
            value=scalar(s.diff(expr,z,j));cap=getattr(q,key)
            lo,hi=current.ep(value);absolute=max(abs(lo),abs(hi))
            upper=c.mpf(0) if cap.log is None else c.exp(cap.log)
            assert absolute<=current.ep(upper)[1],(count,j)
            count+=1
    assert count==273
    return dict(independent_directed_polynomial_mixed_rows=count,physical_rates_and_positive_L_quotient_caps_checked=True,
        actual_field_values_not_materialized=True)


def unchanged_generic_replay():
    """Replay the original method body with symbolic scalar arithmetic only."""
    y,z,E,V,H,P0,delta,R,S=polynomials()
    method=current.generic.GenericMomentRecovery.field_rows
    environment=dict(method.__globals__)
    environment['check_jets']=lambda *args:None
    environment['axial_derivative']=lambda q:s.diff(q,z)
    replay=FunctionType(method.__code__,environment,method.__name__,method.__defaults__,method.__closure__)
    ctx=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    owner=SimpleNamespace(ctx=ctx,P0=P0,own=lambda:H,field=lambda **kwargs:{})
    erows=[s.diff(E,y,j) for j in range(5)];vrows=[s.diff(V,y,j) for j in range(5)]
    result=replay(owner,Z=z,delta=delta,E_rows=erows,V_rows=vrows)
    class SymbolicAlgebra:
        fixed=lambda self,k:s.Rational(str(k))
        add=lambda self,*q:s.Add(*q)
        neg=lambda self,q:-q
        mul=lambda self,a,b:a*b
        scale=lambda self,q,k:q*s.Rational(str(k))
        divide=lambda self,n,d,name:n/d
    histories=[{key:rows[j] for key,rows in result['own_normalized_history_ordinary_y_rows'].items()} for j in range(5)]
    shifted=[{key:s.diff(q,z) for key,q in row.items()} for row in histories]
    expected=current.mixed_recover(SymbolicAlgebra(),[dict(E=e,V=v) for e,v in zip(erows,vrows)],
        histories,shifted,P0,s.diff(P0,z),dict(Z=z,d=1-z*z,L=1-delta*z*z,delta=delta,R=R,Pstar=S,
            Pstar_sqrt_R_over_2=S*s.sqrt(R/2),Pstar_over_sqrt_2R=S/s.sqrt(2*R)))
    baseprefactor=S*s.sqrt(R/2);shearprefactor=S/s.sqrt(2*R);count=0
    for key,rows in result['physical_velocity_pressure_ordinary_y_rows'].items():
        name={'theta':'Utheta','axial':'Uz','radial':'Ur','pressure':'Pi'}[key]
        factor={'theta':S,'axial':S,'radial':baseprefactor,'pressure':S*S}[key]
        for got,want in zip(rows,expected['original_cylindrical_velocity_pressure_y_rows'][name]):
            assert s.cancel(factor*got-want)==0,(key,count);count+=1
    for key,rows in result['full_signed_stress_ordinary_y_rows'].items():
        name=key.removeprefix('inertial_') if key.startswith('inertial_') else key
        factor=shearprefactor if key.startswith('shear_') else baseprefactor*(S if 'quadratic' in key else 1)
        for got,want in zip(rows,expected['original_physical_signed_stress_y_rows'][name]):
            assert s.cancel(factor*got-want)==0,(key,count);count+=1
    assert replay.__code__ is method.__code__ and count==44
    return dict(unchanged_generic_method_body_replayed=True,
        independent_generic_physical_velocity_and_stress_rows=count,
        only_symbolic_jet_checks_and_axial_derivative_adapters_used=True)


def aliases_and_ranges(field,raw):
    f=field.functions;prior=field.lower.band.functions;old=field.lower.functions['recovered_C2_functions'];count=0
    for new,previous in zip(f['actual_C3_profiles_y4']['radial_y_rows'][:3],field.lower.band.velocity_functions()):
        for key in new:assert new[key] is previous[key];count+=4
    for new,previous in zip(f['actual_C3_complete_history_y_rows'][:2],(prior['complete_histories'],prior['complete_history_y'])):
        for key in new:assert new[key] is previous[key];count+=4
    assert f['actual_C3_independent_P0_y_rows'][0] is prior['original_leading_power']['independent_P0']
    assert count==124
    recovered=f['recovered_mixed_C2_functions'];alias_count=0
    pairs=[(recovered['common_radial_Q_y_rows'][0],old['common_radial_Q']),
        (recovered['common_radial_Q_y_rows'][1],old['common_radial_Q_y']),
        (recovered['common_absolute_pressure_y_rows'][0],old['common_absolute_pressure']),
        (recovered['common_absolute_pressure_y_Z_rows'][0],old['common_absolute_pressure_Z'])]
    pairs.extend((rows[0],old['full_signed_inertial_sectors'][key]) for key,rows in recovered['full_signed_inertial_sector_y_rows'].items())
    pairs.extend((rows[0],old['actual_generic_numerators'][key]) for key,rows in recovered['actual_generic_C_B_y_rows'].items())
    pairs.extend((rows[0],old['actual_generic_numerators']['inertial_'+key]) for key,rows in recovered['actual_generic_inertial_y_rows'].items())
    pairs.extend((rows[0],old['original_cylindrical_velocity_pressure'][key]) for key,rows in recovered['original_cylindrical_velocity_pressure_y_rows'].items())
    pairs.append((recovered['original_cylindrical_velocity_pressure_y_rows']['Ur'][1],old['original_cylindrical_velocity_pressure']['Ur_y']))
    for new,previous in pairs:assert new is previous;alias_count+=3
    assert alias_count==51
    input_aliases=0
    for i,rows in enumerate(f['actual_C3_complete_history_y_rows']):
        for key,q in rows.items():
            low=f['actual_C2_history_inputs'][i][key];shift=f['actual_C2_history_Z_inputs'][i][key]
            assert (low.value,low.Z,low.ZZ)==current.target.rows(q)[:3]
            assert (shift.value,shift.Z,shift.ZZ)==current.target.rows(q)[1:]
            input_aliases+=6
    assert input_aliases==150
    encoded=lambda q:json.loads(json.dumps(current.band.control.encoded(current.ranges.record(q))))
    cells=0
    for saved in raw['actual_four_Z_mixed_repaired_recovery_ranges']:
        ends=tuple(saved['exact_Z_cell']);assert encoded(field.quantitative_range(ends))==saved
        assert len(saved['actual_profile_y0_to_y4_Z3_bounds'])==5
        assert len(saved['actual_complete_history_y0_to_y4_Z3_bounds'])==5
        assert saved['range_caps_not_function_values']
        source=field.lower.parameter_source(ends)
        assert current.ep(source['actual_Rc'].coefficient)[0]>0
        cells+=1
    assert cells==4
    return dict(accepted_C3_profile_history_alias_rows=count,accepted_C2_recovered_alias_rows=alias_count,
        ordinary_C3_to_C2_history_consumption_alias_rows=input_aliases,actual_whole_band_Z_cells=cells,
        independent_P0_and_actual_source_positive_radius_retained=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_mixed_repaired_recovery_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentMixedRepairedRecovery(require_checked=False)
        assert field.lower.acceptance_loaded and not field.acceptance_loaded
        assert field.graph.nodes[:len(field.prefix)]==field.prefix and field.graph.nodes==raw['exact_graph_nodes']
        encoded=lambda q:json.loads(json.dumps(current.band.control.encoded(current.ranges.record(q))))
        assert encoded(field.functions)==raw['actual_mixed_repaired_recovery_functions']
        assert raw['source_bindings']==current.source_bindings() and raw['source_family']==field.identity
        compact=compact_profiles(field);history=signed_histories(field)
        calculus=recovery_calculus(field);fixture=directed_fixture(field);aliases=aliases_and_ranges(field,raw)
        generic=unchanged_generic_replay()
        assert not raw['current_numeric_point_field_oracle_installed'] and not raw['actual_generic_stress_cone_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            independent_compact_profile_calculus=compact,independent_signed_history_calculus=history,
            independent_mixed_recovery_calculus=calculus,independent_directed_fixture=fixture,actual_aliases_and_ranges=aliases,
            unchanged_generic_recovery_equivalence=generic,
            original_source_coordinates_before_time_map=True,current_numeric_point_field_oracle_installed=False,
            actual_generic_stress_cone_installed=False,physical_time_Cartesian_global_Rh_heat_and_temporal_recursion_installed=False,
            **dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.band.control.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual mixed repair y4/Z3 input and y4/Z2 velocity, y3/Z2 physical stress checks passed',flush=True)
    return result


if __name__=='__main__':run()
