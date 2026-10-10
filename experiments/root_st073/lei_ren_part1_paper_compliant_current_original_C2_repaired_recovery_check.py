"""Independent signed recovery, C3 input, radius and positivity checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C2_repaired_recovery as current


def symbolic(g,node,bindings):
    if node in bindings:return bindings[node]
    n=g.nodes[node];op=n['operation'];at=lambda i:symbolic(g,i,bindings)
    if op=='exact_rational':return s.Rational(n['numerator'],n['denominator'])
    if op=='sum':return s.Add(*(at(i) for i in n['arguments']))
    if op=='product':return s.Mul(*(at(i) for i in n['arguments']))
    if op=='negative':return -at(n['argument'])
    if op=='positive_quotient':return at(n['numerator'])/at(n['denominator'])
    if op=='analytic_unary':
        v=at(n['argument'])
        if n['name']=='exp':return s.exp(v)
        if n['name']=='log':return s.log(v)
        if n['name']=='positive_sqrt':return s.sqrt(v)
    raise ValueError('Unbound repaired recovery node '+str(n))


def independent_expected(z,delta,R,S,profiles,hist,P0,my):
    E,V,Ey,Vy=(profiles[key] for key in ('E','V','E_y','V_y'))
    m,h,k,e,p=(hist[key] for key in ('m','h','k','e','p'))
    d=1-z*z;L=1-delta*z*z;P=P0+p
    transport=(1-delta)*z*m+d*s.diff(m,z)
    Q=(2*z*V-transport)/L
    Qy=(2*z*Vy-(1-delta)*z*my-d*s.diff(my,z))/L
    sectors=dict(theta_linear=(-E+(1-delta/2)*h-(1-delta)*z*s.diff(h,z)/2)/L,
        theta_quadratic=((2*delta-1)*z*k-d*s.diff(k,z)+E*transport)/L,
        axial_linear=(-V+(1-delta)*(m-z*s.diff(m,z))/2)/L,
        axial_quadratic=(V*transport+2*delta*z*e-d*s.diff(e,z)+2*(1+delta)*z*P-d*s.diff(P,z))/L)
    C=E-2*Ey;B=2*Vy;It=R*(sectors['theta_linear']+S*sectors['theta_quadratic'])
    Iz=R*(sectors['axial_linear']+S*sectors['axial_quadratic'])
    den=C*E;kap=C*C+B*B;excess=kap-2*den;H=C*It-B*Iz;D=H-kap;J=C*Iz+B*It
    return dict(common_E=E,common_V=V,common_radial_Q=Q,common_radial_Q_y=Qy,
        common_absolute_pressure=P,common_absolute_pressure_Z=s.diff(P,z),full_signed_inertial_sectors=sectors,
        actual_generic_numerators=dict(E=E,C=C,B=B,inertial_theta=It,inertial_axial=Iz,
            positive_denominator=den,kappa=kap,kappa_minus2=excess,H0_minus2=H-2*den,D=D,J=J),
        full_signed_generic_quotients=dict(a=C/E,b=B/E,t0=-B/C,p1=It/E,p2=Iz/E,
            kappa=kap/den,D=D/den,J=J/den,Delta=kap/den-2),
        full_signed_generic_discriminant=(2*D*D*den-excess*J*J)/den**3,
        original_cylindrical_velocity_pressure=dict(Utheta=S*E,Uz=S*V,Ur=S*s.sqrt(R/2)*Q,
            Ur_y=S*s.sqrt(R/2)*(Qy+Q/2),Pi=S*S*P))


def independent_functions(field):
    g=field.band.control.ranges.phase.built['graph'];f=field.functions
    z=s.Symbol('Z');delta=s.Symbol('delta',positive=True);R,S=s.symbols('R S',positive=True)
    profiles={key:s.Function(key)(z) for key in ('E','V','E_y','V_y')}
    hist={key:s.Function(key)(z) for key in current.current.RATES};P0=s.Function('P0')(z)
    bindings={};count=0
    def bind(q,expr):
        for order,key in enumerate(('value','Z','ZZ')):bindings[getattr(q,key).node]=s.diff(expr,z,order)
    def check(q,expr,label):
        nonlocal count
        for order,key in enumerate(('value','Z','ZZ')):
            got=symbolic(g,getattr(q,key).node,bindings)
            assert s.cancel(s.together(got-s.diff(expr,z,order)))==0,(label,order)
            count+=1
    inputs=f['actual_C3_repaired_inputs'];params=f['actual_source_parameters']
    for key,expr in profiles.items():bind(inputs[key],expr)
    for key,expr in hist.items():
        bind(inputs['histories'][key],expr);bind(inputs['history_Z'][key],s.diff(expr,z))
    bind(inputs['P0'],P0);bind(inputs['P0_Z'],s.diff(P0,z))
    bind(inputs['history_y']['m'],profiles['V']-hist['m'])
    bind(inputs['history_y_Z']['m'],s.diff(profiles['V']-hist['m'],z))
    bindings[params['Z'].value.node]=z;bindings[params['delta'].value.node]=delta
    bindings[params['R'].value.node]=R;bindings[params['Pstar'].value.node]=S
    expected=independent_expected(z,delta,R,S,profiles,hist,P0,profiles['V']-hist['m'])
    recovered=f['recovered_C2_functions']
    for key in ('common_E','common_V','common_radial_Q','common_radial_Q_y',
            'common_absolute_pressure','common_absolute_pressure_Z'):
        check(recovered[key],expected[key],key)
    for key,q in recovered['full_signed_inertial_sectors'].items():check(q,expected['full_signed_inertial_sectors'][key],key)
    Q,Qy=expected['common_radial_Q'],expected['common_radial_Q_y'];d=1-z*z;L=1-delta*z*z
    assert s.cancel(L*(Q+Qy)-((1+delta)*z*profiles['V']-d*s.diff(profiles['V'],z)+2*z*profiles['V_y']))==0
    # Separate verified sectors from algebraic composition, avoiding an
    # unnecessary fully expanded quartic symbolic rational expression.
    sectors={key:s.Function(key)(z) for key in recovered['full_signed_inertial_sectors']}
    for key,expr in sectors.items():bind(recovered['full_signed_inertial_sectors'][key],expr)
    E,V,Ey,Vy=(profiles[key] for key in ('E','V','E_y','V_y'))
    C=E-2*Ey;B=2*Vy;It=R*(sectors['theta_linear']+S*sectors['theta_quadratic'])
    Iz=R*(sectors['axial_linear']+S*sectors['axial_quadratic'])
    den=C*E;kap=C*C+B*B;H=C*It-B*Iz
    nums=dict(E=E,C=C,B=B,inertial_theta=It,inertial_axial=Iz,positive_denominator=den,
        kappa=kap,kappa_minus2=kap-2*den,H0_minus2=H-2*den,D=H-kap,J=C*Iz+B*It)
    for key,q in recovered['actual_generic_numerators'].items():check(q,nums[key],key)
    abstract={key:s.Function('N_'+key)(z) for key in nums}
    for key,expr in abstract.items():bind(recovered['actual_generic_numerators'][key],expr)
    n=abstract
    quotient=dict(a=n['C']/n['E'],b=n['B']/n['E'],t0=-n['B']/n['C'],
        p1=n['inertial_theta']/n['E'],p2=n['inertial_axial']/n['E'],
        kappa=n['kappa']/n['positive_denominator'],D=n['D']/n['positive_denominator'],
        J=n['J']/n['positive_denominator'],Delta=n['kappa']/n['positive_denominator']-2)
    for key,q in recovered['full_signed_generic_quotients'].items():check(q,quotient[key],key)
    disc=(2*n['D']**2*n['positive_denominator']-n['kappa_minus2']*n['J']**2)/n['positive_denominator']**3
    check(recovered['full_signed_generic_discriminant'],disc,'discriminant')
    # Bind the verified common pressure/Q and the E/V alias used by the
    # physical prefactors; derivative rows remain independent arbitrary jets.
    radial,radial_y,pressure=[s.Function(key)(z) for key in ('Q','Q_y','P')]
    bind(recovered['common_radial_Q'],radial);bind(recovered['common_radial_Q_y'],radial_y)
    bind(recovered['common_absolute_pressure'],pressure)
    physical=dict(Utheta=S*n['E'],Uz=S*V,Ur=S*s.sqrt(R/2)*radial,
        Ur_y=S*s.sqrt(R/2)*(radial_y+radial/2),Pi=S*S*pressure)
    for key,q in recovered['original_cylindrical_velocity_pressure'].items():check(q,physical[key],key)
    assert count==108
    return dict(independent_signed_recovery_and_quotient_rows=count,
        full_pressure_linear_quadratic_meridional_sectors_checked=True,
        exact_own_moment_radial_divergence_identity_checked=True,
        physical_radial_first_y_half_power_checked=True)


def independent_directed_fixture(field):
    c=field.c;bd=current.ranges.Bounds(c);z=s.Symbol('Z');at=s.Rational(1,4)
    delta=s.Rational(1,5);R=s.Integer(7);S=s.Integer(3)
    profiles=dict(E=2+z*z/10,V=1-z/5+z*z/7,E_y=-s.Rational(2,5)+z*z/30,V_y=z*z/10+s.Rational(1,50))
    hist=dict(m=1+z+z*z,h=2-z+z**3/5,k=3+2*z/5-z*z/5,e=7*z*z/5-s.Rational(1,2),p=-z/5+z**3+s.Rational(3,10))
    P0=9/s.Integer(10)+3*z/10-z*z/5;my=profiles['V']-hist['m']
    scalar=lambda expr:c.mpf(str(s.N(expr.subs(z,at),100)))
    def jet(expr):return current.ranges.JetBound(*(bd.constant(scalar(s.diff(expr,z,j))) for j in range(3)))
    inputs={**{key:jet(expr) for key,expr in profiles.items()},
        'histories':{key:jet(expr) for key,expr in hist.items()},
        'history_Z':{key:jet(s.diff(expr,z)) for key,expr in hist.items()},
        'history_y':{'m':jet(my)},'history_y_Z':{'m':jet(s.diff(my,z))},'P0':jet(P0),'P0_Z':jet(s.diff(P0,z))}
    parameters=dict(Z=jet(z),d=jet(1-z*z),L=jet(1-delta*z*z),delta=jet(delta),R=jet(R),Pstar=jet(S),
        Pstar_sqrt_R_over_2=jet(S*s.sqrt(R/2)))
    lowers=dict(L=c.ln(c.mpf('.8')),E=c.mpf(0),C=c.mpf(0),den=c.mpf(0),den3=c.mpf(0))
    caps=current.recover(current.MagnitudeAlgebra(c,lowers),inputs,parameters)
    expected=independent_expected(z,delta,R,S,profiles,hist,P0,my);count=0
    def check(a,b):
        nonlocal count
        if isinstance(a,current.ranges.JetBound):
            for j,key in enumerate(('value','Z','ZZ')):
                value=scalar(s.diff(b,z,j));cap=getattr(a,key)
                upper=c.mpf(0) if cap.log is None else c.exp(cap.log)
                assert ep_abs(value)[1]<=current.ep(upper)[0],(j,b)
                count+=1
        else:
            for key,v in a.items():check(v,b[key])
    def ep_abs(q):
        lo,hi=current.ep(q);return min(abs(lo),abs(hi)),max(abs(lo),abs(hi))
    check(caps,expected);assert count==108
    return dict(independent_directed_finite_fixture_rows=count,
        independent_polynomial_jets_and_positive_quotient_arithmetic_checked=True,
        fixture_only_no_actual_field_values_materialized=True)


def actual_inputs_and_sources(field,raw):
    c=field.c;bf=field.band.functions;inputs=field.functions['actual_C3_repaired_inputs']
    alias_count=0
    def aliases(q,prior,shift=False):
        nonlocal alias_count
        want=current.target.rows(prior)[1:4] if shift else current.target.rows(prior)[:3]
        assert (q.value,q.Z,q.ZZ)==want
        alias_count+=3
    for name,y,key in (('E',0,'E'),('V',0,'V'),('E_y',1,'E'),('V_y',1,'V')):
        aliases(inputs[name],field.band.velocity_functions()[y][key])
    for name,prior,shift in (('histories',bf['complete_histories'],False),('history_Z',bf['complete_histories'],True),
            ('history_y',bf['complete_history_y'],False),('history_y_Z',bf['complete_history_y'],True)):
        for key,q in inputs[name].items():aliases(q,prior[key],shift)
    P0=bf['original_leading_power']['independent_P0']
    aliases(inputs['P0'],P0);aliases(inputs['P0_Z'],P0,True)
    assert alias_count==78
    counts=dict(actual_Z_cells=0,actual_parameter_radius_scale_and_delta_bindings=0)
    encoded=lambda q:json.loads(json.dumps(current.band.control.encoded(current.ranges.record(q))))
    for saved in raw['actual_four_Z_C2_repaired_recovery_ranges']:
        ends=tuple(saved['exact_Z_cell']);assert encoded(field.quantitative_range(ends))==saved
        source=field.parameter_source(ends);op=field.band.control.ranges.phase.outer.owner.owner(ends);f=op.flow
        Rc=op.Rm_factor*f.factor((0,0,0,0,0),c.exp(40)+20)
        assert current.current.current.previous.equivalent_rows([Rc],[source['actual_Rc']])
        assert source['actual_Pstar'] is op.Pstar and source['actual_Rm_factor'] is op.Rm_factor
        assert source['actual_delta'] is op.reference.delta and source['original_source_P0_is_same_object']
        logP=c.exp(40)+11;de=c.exp(-4*logP-30)
        assert current.ep(de)[1]<current.ep(c.mpf('1e-200'))[0]
        identity=source['actual_similarity_parameter_identity']
        assert identity['exact_live_delta_tuple_equals_branch'] and identity['same_source_log_basis_and_delta_context']
        assert op.reference.delta.ctx is c and op.reference.logP.ctx is c
        assert op.reference.delta._mpi_==c.exp(-4*op.reference.logP-30)._mpi_
        assert op.reference.logP._mpi_==(f.logs[1]/2)._mpi_
        assert current.current.current.previous.equivalent_rows([f.factor((0,.5,0,0,0))],[op.Pstar])
        proof=field.proofs[ends]
        assert current.ep(proof['actual_meridional_L_lower'])[0]>0
        for q in proof['actual_selected_N_F_Fy_log_caps']:
            assert current.ep(q)[1]<current.ep(c.ln(c.mpf(1)/8))[0]
        assert current.ep(proof['actual_alpha'])[1]<current.ep(c.mpf(2)/3)[0]
        assert not saved['actual_generic_cone_signed_margins_not_admitted'] is False
        counts['actual_Z_cells']+=1;counts['actual_parameter_radius_scale_and_delta_bindings']+=4
    assert counts['actual_Z_cells']==4
    g=field.band.control.ranges.phase.built['graph'];functions=field.functions
    Rc_node=g.nodes[functions['actual_Rc'].node]
    assert Rc_node['quantity']=='actual_source_radius' and Rc_node['Z_independent']
    assert Rc_node['coordinate']==g.constant(2).node and Rc_node['quantity_path']==['original_closed_O3_background','actual_source_radius']
    x=functions['original_band_variable'];assert functions['actual_R']==g.mul(functions['actual_Rc'],x)
    built=field.band.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    symbols={key:s.Symbol(key,real=True) for key in built['parameters']}
    bindings={q.node:symbols[key] for key,q in built['parameters'].items()};bindings[x.node]=s.Symbol('x',positive=True)
    P,C,B,sc=[symbols[key] for key in ('logP','logC','hbB','sc')]
    Rm=s.log(25)+1000+4*P+s.log(s.Rational(11,10))-B*sc/2+10*(C+P)-6
    absolute=symbolic(g,functions['original_absolute_radius_offset'].node,bindings)
    relative=symbolic(g,functions['source_relative_Rm_radius_offset'].node,bindings)
    # logP is the original exact exp(40)+11, not an independent parameter.
    assert s.simplify(s.expand_log((absolute-Rm-relative).subs(P,s.exp(40)+11),force=True))==0
    assert functions['original_Rm_factor_not_counted_twice']
    return dict(accepted_C3_input_derivative_aliases=alias_count,**counts,
        exact_source_radius_without_repeated_Rm_checked=True,
        actual_independent_pressure_and_same_Z_source_checked=True,
        formal_absolute_and_relative_radius_identity_checked=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C2_repaired_recovery_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC2RepairedRecovery(require_checked=False)
        assert field.band.acceptance_loaded and not field.acceptance_loaded
        g=field.band.control.ranges.phase.built['graph']
        assert g.nodes[:len(field.prefix)]==field.prefix and g.nodes==raw['exact_graph_nodes']
        encoded=lambda q:json.loads(json.dumps(current.band.control.encoded(current.ranges.record(q))))
        assert encoded(field.functions)==raw['actual_C2_repaired_recovery_functions']
        assert raw['source_bindings']==current.source_bindings() and raw['source_family']==field.identity
        functions=independent_functions(field);fixture=independent_directed_fixture(field)
        inputs=actual_inputs_and_sources(field,raw)
        # Elementary factor proof, independent of the magnitude provider.
        assert s.Rational(1,2)-s.Rational(1,8)==s.Rational(3,8)>s.Rational(1,4)
        assert 2*s.Rational(1,2)-3*s.Rational(1,8)==s.Rational(5,8)>s.Rational(1,2)
        assert s.Rational(2)**s.Rational(-2,3)>s.Rational(1,2)
        assert not raw['current_numeric_point_field_oracle_installed'] and not raw['actual_generic_stress_cone_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            actual_C2_repaired_radial_inertial_functions_installed=True,
            actual_repaired_positive_L_E_C_margins_installed=True,
            independent_signed_function_calculus=functions,independent_directed_arithmetic=fixture,
            actual_input_radius_and_parameter_checks=inputs,
            whole_actual_selected_N_band_positive_margin_theorem_checked=True,
            current_numeric_point_field_oracle_installed=False,actual_generic_stress_cone_installed=False,
            physical_time_Cartesian_dispatcher_global_Rh_heat_and_temporal_recursion_installed=False,
            **dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.band.control.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual repaired radial/inertial C2 signed calculus, source radius and positive margin checks passed',flush=True)
    return result


if __name__=='__main__':run()
