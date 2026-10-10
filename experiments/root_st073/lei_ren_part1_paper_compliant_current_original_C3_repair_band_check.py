"""Independent mixed bump/profile calculus and actual C3 band closure checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_repair_band as current
import lei_ren_part1_paper_compliant_current_original_C2_density_transport_check as calculus

phase=current.phase;control=current.control


def independent_bump_calculus():
    r=s.Symbol('r');beta=s.exp(-1/(1-r*r))
    assert s.simplify(s.diff(beta,r,2)*(1-r*r)**4/beta-(6*r**4-2))==0
    y,center,ell,J0=s.symbols('y center ell J0',real=True)
    arg=(y-center)/ell;b=beta.subs(r,arg)*s.exp(-y)/(ell*J0)
    first=s.diff(beta,r).subs(r,arg)*s.exp(-y)/(ell**2*J0)-b
    second=s.diff(beta,r,2).subs(r,arg)*s.exp(-y)/(ell**3*J0) \
        -2*s.diff(beta,r).subs(r,arg)*s.exp(-y)/(ell**2*J0)+b
    assert s.simplify(s.diff(b,y)-first)==0 and s.simplify(s.diff(b,y,2)-second)==0
    w=s.Symbol('w',positive=True)
    for n,C,wmax,cap in ((1,2,s.Rational(1,2),8*s.exp(-2)),(2,8,s.Rational(1,4),2048*s.exp(-4))):
        envelope=C*s.exp(-1/w)*w**(-2*n)
        assert s.simplify(s.diff(s.log(envelope),w)-(1-2*n*w)/w**2)==0
        assert envelope.subs(w,wmax)==cap
    assert phase.flat_source.BETA_POLYNOMIALS[2]==[-2,0,0,0,6]
    for power in range(5):assert s.limit(s.exp(-1/w)*w**(-power),w,0,dir='+')==0
    return dict(original_beta_second_polynomial_checked=True,
        actual_log_radius_first_and_second_bump_derivatives_checked=True,
        original_beta_global_derivative_caps=['e^-1','8e^-2','2048e^-4'],
        original_support_endpoint_y0_y1_y2_jets_flat=True)


def independent_profiles(field):
    g=field.control.ranges.phase.built['graph'];built=field.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    z=s.Symbol('Z');A=s.Function('Am')(z);h=[s.Function('h'+str(i))(z) for i in range(5)]
    alpha,P,N=s.symbols('alpha power N',positive=True);bindings={built['N'].node:N}
    amplitude=field.control.ranges.target.functions['actual_terminal_amplitude']
    for order,key in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(amplitude,key).node]=s.diff(A,z,order)
    for q,f in zip(field.control.control_functions(),h):
        for order,key in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(q,key).node]=s.diff(f,z,order)
    count=0
    for profiles in (field.functions['profiles_at_band_x'],field.functions['profiles_at_integration_t']):
        local=dict(bindings);local[profiles['original_power'].node]=P;local[profiles['original_alpha'].node]=alpha
        for j,(bumps,row) in enumerate(zip((profiles['bumps'],profiles['bump_y'],profiles['bump_yy']),profiles['radial_y_rows'])):
            bb=s.symbols('b%d_0:3'%j)
            local.update({q.node:v for q,v in zip(bumps,bb)})
            F=sum(bb[i]*h[i+2] for i in range(3))/N;G=(bb[0]*h[0]+bb[2]*h[1])/N
            original=(-alpha)**j*P*A
            expected=dict(F=F,G=G,original_E=original,delta_E=A*F,delta_V=A*G,E=original+A*F,V=A*G)
            for key,q in row.items():
                for order,handle in enumerate((q.value,q.Z,q.ZZ,q.ZZZ)):
                    got=calculus.symbolic(g,handle.node,local)
                    assert s.simplify(got-s.diff(expected[key],z,order))==0,(j,key,order)
                    count+=1
    return dict(actual_profile_and_mixed_y_Z_rows_checked=count,
        radial_orders=[0,1,2],axial_orders=[0,1,2,3],
        all_actual_amplitude_cross_terms_and_limit_third_rows_retained=True)


def independent_density(field):
    g=field.control.ranges.phase.built['graph'];z=s.Symbol('Z');bindings={}
    row=field.functions['profiles_at_integration_t']['radial_y_rows'][0]
    E,dE,dV=[s.Function(k)(z) for k in ('E','delta_E','delta_V')]
    for key,f in zip(('original_E','delta_E','delta_V'),(E,dE,dV)):
        for order,name in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(row[key],name).node]=s.diff(f,z,order)
    expected=dict(m=dV,h=dE,k=(E+dE)*dV,e=dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    for key,q in field.functions['signed_density'].items():
        for order,handle in enumerate((q.value,q.Z,q.ZZ,q.ZZZ)):
            assert s.simplify(calculus.symbolic(g,handle.node,bindings)-s.diff(expected[key],z,order))==0
    return dict(independent_signed_band_density_rows=20,energy_and_pressure_signs_and_cubic_cross_terms_checked=True)


def independent_endpoints(field):
    z=s.Symbol('Z');A=s.Function('Am')(z);N,mu=s.symbols('N mu',positive=True)
    g=field.control.ranges.phase.built['graph'];built=field.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    bindings={built['N'].node:N,built['parameters']['mu'].node:mu}
    bindings[g.unary('log',g.constant(2)).node]=s.log(2)
    amplitude=field.control.ranges.target.functions['actual_terminal_amplitude']
    for order,key in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(amplitude,key).node]=s.diff(A,z,order)
    R={row:s.Function('R'+str(i))(z) for i,row in enumerate(current.controls.ROWS)}
    for row,q in field.functions['control_residual_C3'].items():
        for order,key in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(q,key).node]=s.diff(R[row],z,order)
    expected=dict(m=A*R['M']/(2*N),h=A*R['I']/(2**s.Rational(3,2)*N),
        e=A*A*R['S']/(2*N),p=A*A*R['Cp']/N,
        k=A*A*(R['M']+mu*R[current.controls.ROWS[1]])/(2**s.Rational(3,2)*N))
    for key,q in field.functions['relative_terminal_residual_rhs'].items():
        for order,handle in enumerate((q.value,q.Z,q.ZZ,q.ZZZ)):
            assert s.simplify(calculus.symbolic(g,handle.node,bindings)-s.diff(expected[key],z,order))==0,(key,order)
    # Third-row FTC/inlet, with fixed native x endpoint and own recovery rate.
    x,t=s.symbols('x t',positive=True);seed=s.Function('incoming')(z);D=s.Function('density')
    for rate in current.current.RATES.values():
        lam=s.Rational(rate.numerator,rate.denominator)
        H=x**(-lam)*(seed+s.Integral(t**(lam-1)*D(t,z),(t,1,x)))
        HZZ=s.diff(H,z,3)
        assert s.simplify(HZZ.subs(x,1)-s.diff(seed,z,3))==0
        assert s.simplify(x*s.diff(HZZ,x)+lam*HZZ-s.diff(D(x,z),z,3))==0
    return dict(independent_relative_terminal_C3_rows=20,
        independent_own_rate_third_inlet_and_FTC_identities=5,
        identity_scope='relative correction history at x=2; not absolute exterior closure')


def independent_leading_power(field):
    g=field.control.ranges.phase.built['graph'];built=field.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    z,T=s.symbols('Z T');mu=s.Symbol('mu',positive=True);A=s.Function('Am')(z)
    bindings={built['parameters']['mu'].node:mu}
    x=field.functions['original_band_variable'];bindings[g.unary('log',x).node]=T
    amplitude=field.control.ranges.target.functions['actual_terminal_amplitude']
    for order,key in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(amplitude,key).node]=s.diff(A,z,order)
    leading=field.functions['original_leading_power'];seeds={key:s.Function(key+'0')(z) for key in current.current.RATES}
    for key,q in leading['actual_endpoint_seeds'].items():
        for order,row in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(q,row).node]=s.diff(seeds[key],z,order)
    expected=dict(m=seeds['m']*s.exp(-T),k=seeds['k']*s.exp(-3*T/2),
        h=s.exp(-3*T/2)*(seeds['h']+A*(s.exp((1-mu)*T)-1)/(1-mu)),
        e=s.exp(-T)*(seeds['e']-A*A*(1-s.exp(-2*mu*T))/(4*mu)),
        p=seeds['p']+A*A*(1-s.exp(-(1+2*mu)*T))/(2*(1+2*mu)))
    for key,q in leading['leading_histories'].items():
        for order,handle in enumerate((q.value,q.Z,q.ZZ,q.ZZZ)):
            assert s.simplify(calculus.symbolic(g,handle.node,bindings)-s.diff(expected[key],z,order))==0,(key,order)
    return dict(independent_original_power_semigroup_C3_rows=20,
        original_endpoint_histories_and_amplitude_cross_terms_retained=True)


def independent_complete_history_y(field):
    g=field.control.ranges.phase.built['graph'];z=s.Symbol('Z');bindings={}
    E,V=[s.Function(key)(z) for key in ('E','V')]
    profiles=field.velocity_functions()[0]
    hist={key:s.Function('complete_'+key)(z) for key in current.current.RATES}
    for key,f in [('E',E),('V',V)]:
        for order,row in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(profiles[key],row).node]=s.diff(f,z,order)
    for key,f in hist.items():
        for order,row in enumerate(('value','Z','ZZ','ZZZ')):bindings[getattr(field.history_functions()[key],row).node]=s.diff(f,z,order)
    expected=dict(m=-hist['m']+V,h=-s.Rational(3,2)*hist['h']+E,
        k=-s.Rational(3,2)*hist['k']+E*V,e=-hist['e']+V*V-E*E/2,p=E*E/2)
    for key,q in field.functions['complete_history_y'].items():
        assert s.simplify(calculus.symbolic(g,q.ZZZ.node,bindings)-s.diff(expected[key],z,3))==0,key
    return dict(independent_complete_history_first_y_third_Z_rows=5,
        original_complete_own_rate_and_signed_cubic_FTC_checked=True)


def preserve_C2_handles(field):
    old=field.lower.functions;new=field.functions;count=0
    def same(a,b):
        nonlocal count
        assert current.target.rows(a)[:3]==(b.value,b.Z,b.ZZ)
        count+=3
    for name in ('profiles_at_band_x','profiles_at_integration_t'):
        for row,prior in zip(new[name]['radial_y_rows'],old[name]['radial_y_rows']):
            for key,q in row.items():same(q,prior[key])
    for name in ('signed_density','partial_correction_contributions','partial_correction_histories',
            'complete_histories','complete_history_y','relative_terminal_actual','relative_terminal_residual_rhs'):
        for key,q in new[name].items():same(q,old[name][key])
    for key,q in new['actual_incoming_correction_C3'].items():same(q,old['actual_incoming_correction_C2'][key])
    for key,q in new['control_residual_C3'].items():same(q,old['control_residual_C2'][key])
    for name in ('actual_endpoint_seeds','leading_histories'):
        for key,q in new['original_leading_power'][name].items():same(q,old['original_leading_power'][name][key])
    same(new['original_leading_power']['independent_P0'],old['original_leading_power']['independent_P0'])
    return dict(accepted_C2_function_handle_aliases_checked=count,
        existing_second_functions_and_independent_P0_exactly_retained=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C3_repair_band_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC3RepairBand(require_checked=False)
        assert not field.acceptance_loaded and field.control.acceptance_loaded and field.lower.acceptance_loaded
        g=field.control.ranges.phase.built['graph'];functions=field.functions
        encoded=lambda q:json.loads(json.dumps(control.encoded(current.ranges.record(q))))
        assert g.nodes[:len(field.prefix)]==field.prefix and g.nodes==raw['exact_graph_nodes']
        assert raw['original_C3_control_and_C2_band_graph_prefix_length']==len(field.prefix)
        assert encoded(functions)==raw['actual_C3_repair_band_functions'] and raw['source_family']==field.identity
        assert raw['actual_selected_repair_integer']==field.control.bounds['actual_same_repair_integer']
        assert field.control.lower is field.lower.control
        assert field.control.ranges.phase is field.lower.control.ranges.phase
        aliases=preserve_C2_handles(field)
        profiles=independent_profiles(field);density=independent_density(field)
        endpoints=independent_endpoints(field);leading=independent_leading_power(field)
        dy=independent_complete_history_y(field)
        assert raw['source_bindings']==current.source_bindings()
        built=field.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
        for key,q in functions['partial_correction_histories'].items():
            assert q.ZZZ!=g.zero and functions['actual_incoming_correction_C3'][key].ZZZ!=g.zero
            n=g.nodes[functions['relative_terminal_zero_certificates'][key].node]
            assert n['relative_not_absolute_exterior'] and n['source_family']==field.identity
            assert n['equivalent_residual_rows']==control.encoded(functions['relative_terminal_residual_rhs'][key])
            assert n['expression_rows']==control.encoded(functions['relative_terminal_actual'][key])
            assert n['actual_third_zero_theorem_nodes']==current.current.encode_graph(field.control.functions['actual_third_zero_theorem_nodes'])
            assert n['actual_C2_relative_zero_theorem']==field.lower.functions['relative_terminal_zero_certificates'][key].node
            t=g.symbol(functions['partial_integration_variable']);x=functions['original_band_variable'];rate=current.current.RATES[key]
            kernel=g.unary('exp',g.mul(g.constant(rate-1),g.unary('log',t)))
            decay=g.unary('exp',g.neg(g.mul(g.constant(rate),g.unary('log',x))))
            term=functions['partial_correction_contributions'][key].ZZZ
            arguments=g.nodes[term.node]['arguments'] if g.nodes[term.node]['operation']=='product' else [term.node]
            integral_nodes=[g.nodes[i] for i in arguments if g.nodes[i]['operation']=='definite_integral']
            assert len(integral_nodes)==1 and decay.node in arguments
            integral=integral_nodes[0]
            assert integral['integrand']==g.mul(kernel,functions['signed_density'][key].ZZZ).node
            assert integral['lower']==g.one.node and integral['upper']==x.node
            assert integral['ordinary_slow_Z_derivative_order']==3 and integral['original_physical_Jacobian_applied_once']
            assert integral['endpoint_independent_of_Z'] and integral['original_recovery_rate']==str(rate)
            substitution=g.nodes[functions['relative_terminal_actual'][key].ZZZ.node]
            assert substitution['expression']==q.ZZZ.node and substitution['value']==g.constant(2).node
        P0=functions['original_leading_power']['independent_P0'];n=g.nodes[P0.ZZZ.node]
        assert n['quantity']=='P0' and n['Z_order']==3 and n['Taylor_coefficient_factorial']==6
        assert not n['derivative_of_range_endpoint']
        target=json.loads(gzip.decompress((current.HERE/current.ranges.NAME).read_bytes()))
        counts=dict(actual_Z_cells=0,actual_endpoint_history_rows=0,actual_endpoint_P0_rows=0,
            original_incoming_third_memories_retained=5,relative_third_terminal_identities=5)
        for saved,parent_row in zip(raw['actual_four_Z_C3_band_ranges'],target['actual_four_Z_C3_target_transports']):
            ends=tuple(saved['exact_Z_cell']);replay=field.quantitative_range(ends,parent_row)
            assert encoded(replay)==saved
            endpoint=field.endpoint_source(ends);op=field.control.ranges.phase.outer.owner.owner(ends);f=op.flow
            inlet=field.source.query(ends,(1,1));actual=inlet['raw']['raw_current_radius_y_derivative_axial_coefficients']
            for key,rows in actual['histories'].items():
                normalized=[q*f.factor((0,-.5,0,0,0)) for q in rows[0]] if key in ('m','k') else rows[0]
                want=[normalized[0],normalized[1],2*normalized[2],6*normalized[3]]
                assert current.current.current.previous.equivalent_rows(want,endpoint['ordinary_leading_history_C3'][key])
                counts['actual_endpoint_history_rows']+=4
            for q,want in zip(endpoint['ordinary_independent_P0_C3'],(op.P0[0],op.P0[1],2*op.P0[2],6*op.P0[3])):
                assert q.ctx is field.c and q.scale.bases is f.logs and q.ledger is f.ledger
                assert current.current.current.previous.equivalent_rows([q],[want])
                counts['actual_endpoint_P0_rows']+=1
            whole=field.leading_range(ends)
            assert saved['actual_whole_band_leading_source']['same_original_P0']
            assert all(q.zero for rows in whole['actual_whole_band_V_y_C3'] for q in rows)
            for key,rows in whole['actual_whole_reserved_band_leading_history_C3'].items():
                assert current.current.current.previous.equivalent_rows(rows[:3],
                    field.lower.leading_range(ends)['actual_whole_reserved_band_leading_history_C2'][key])
            assert saved['range_caps_not_function_values'] and saved['selected_repair_integer_unchanged']
            for key,rate in current.current.RATES.items():
                r=field.c.mpf(rate.numerator)/rate.denominator
                want=field.c.ln(2) if not rate else -field.c.expm1(-r*field.c.ln(2))/r
                assert encoded(want)==saved['actual_own_rate_positive_masses'][key]
                assert len(saved['actual_complete_history_C3_bounds'][key])==4
                assert len(saved['actual_complete_history_y_C3_bounds'][key])==4
            counts['actual_Z_cells']+=1
        assert counts['actual_Z_cells']==4 and counts['actual_endpoint_history_rows']==80 and counts['actual_endpoint_P0_rows']==16
        assert functions['next_generic_recovery_available_axial_order']==2
        assert raw['radial_inertial_recovery_C2_installed'] is False
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            C3_repair_band_fields_and_histories_installed=True,
            independent_mixed_profile_calculus=profiles,independent_signed_density=density,
            independent_relative_endpoints=endpoints,independent_original_power_semigroup=leading,
            independent_complete_history_y=dy,accepted_C2_aliases=aliases,actual_replay_counts=counts,
            actual_C2_function_handles_and_independent_P0_retained=True,
            radial_inertial_recovery_C2_installed=False,C2_Rh_join_and_absolute_exterior_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(control.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual C3 repaired-band mixed velocity, signed histories, relative terminal and range checks passed',flush=True)
    return result


if __name__=='__main__':run()
