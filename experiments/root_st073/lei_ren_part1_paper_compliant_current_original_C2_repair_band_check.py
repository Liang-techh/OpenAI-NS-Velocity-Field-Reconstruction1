"""Independent mixed bump/profile calculus and actual C2 band closure checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C2_repair_band as current
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
    g=field.control.ranges.phase.built['graph'];built=field.control.functions['exact_C1_limit_adapter']
    z=s.Symbol('Z');A=s.Function('Am')(z);h=[s.Function('h'+str(i))(z) for i in range(5)]
    alpha,P,N=s.symbols('alpha power N',positive=True);bindings={built['N'].node:N}
    amplitude=field.control.ranges.target.functions['actual_terminal_amplitude']
    for order,key in enumerate(('value','Z','ZZ')):bindings[getattr(amplitude,key).node]=s.diff(A,z,order)
    for q,f in zip(field.control.control_functions(),h):
        for order,key in enumerate(('value','Z','ZZ')):bindings[getattr(q,key).node]=s.diff(f,z,order)
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
                for order,handle in enumerate((q.value,q.Z,q.ZZ)):
                    got=calculus.symbolic(g,handle.node,local)
                    assert s.simplify(got-s.diff(expected[key],z,order))==0,(j,key,order)
                    count+=1
    return dict(actual_profile_and_mixed_y_Z_rows_checked=count,
        radial_orders=[0,1,2],axial_orders=[0,1,2],
        both_actual_amplitude_cross_terms_and_limit_second_rows_retained=True)


def independent_density(field):
    g=field.control.ranges.phase.built['graph'];z=s.Symbol('Z');bindings={}
    row=field.functions['profiles_at_integration_t']['radial_y_rows'][0]
    E,dE,dV=[s.Function(k)(z) for k in ('E','delta_E','delta_V')]
    for key,f in zip(('original_E','delta_E','delta_V'),(E,dE,dV)):
        for order,name in enumerate(('value','Z','ZZ')):bindings[getattr(row[key],name).node]=s.diff(f,z,order)
    expected=dict(m=dV,h=dE,k=(E+dE)*dV,e=dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    for key,q in field.functions['signed_density'].items():
        for order,handle in enumerate((q.value,q.Z,q.ZZ)):
            assert s.simplify(calculus.symbolic(g,handle.node,bindings)-s.diff(expected[key],z,order))==0
    return dict(independent_signed_band_density_rows=15,energy_and_pressure_signs_and_second_cross_terms_checked=True)


def independent_endpoints(field):
    z=s.Symbol('Z');A=s.Function('Am')(z);N,mu=s.symbols('N mu',positive=True)
    g=field.control.ranges.phase.built['graph'];built=field.control.functions['exact_C1_limit_adapter']
    bindings={built['N'].node:N,built['parameters']['mu'].node:mu}
    bindings[g.unary('log',g.constant(2)).node]=s.log(2)
    amplitude=field.control.ranges.target.functions['actual_terminal_amplitude']
    for order,key in enumerate(('value','Z','ZZ')):bindings[getattr(amplitude,key).node]=s.diff(A,z,order)
    R={row:s.Function('R'+str(i))(z) for i,row in enumerate(current.controls.ROWS)}
    for row,q in field.functions['control_residual_C2'].items():
        for order,key in enumerate(('value','Z','ZZ')):bindings[getattr(q,key).node]=s.diff(R[row],z,order)
    expected=dict(m=A*R['M']/(2*N),h=A*R['I']/(2**s.Rational(3,2)*N),
        e=A*A*R['S']/(2*N),p=A*A*R['Cp']/N,
        k=A*A*(R['M']+mu*R[current.controls.ROWS[1]])/(2**s.Rational(3,2)*N))
    for key,q in field.functions['relative_terminal_residual_rhs'].items():
        for order,handle in enumerate((q.value,q.Z,q.ZZ)):
            assert s.simplify(calculus.symbolic(g,handle.node,bindings)-s.diff(expected[key],z,order))==0,(key,order)
    # Second-row FTC/inlet, with fixed native x endpoint and own recovery rate.
    x,t=s.symbols('x t',positive=True);seed=s.Function('incoming')(z);D=s.Function('density')
    for rate in current.current.RATES.values():
        lam=s.Rational(rate.numerator,rate.denominator)
        H=x**(-lam)*(seed+s.Integral(t**(lam-1)*D(t,z),(t,1,x)))
        HZZ=s.diff(H,z,2)
        assert s.simplify(HZZ.subs(x,1)-s.diff(seed,z,2))==0
        assert s.simplify(x*s.diff(HZZ,x)+lam*HZZ-s.diff(D(x,z),z,2))==0
    return dict(independent_relative_terminal_C2_rows=15,
        independent_own_rate_second_inlet_and_FTC_identities=5,
        identity_scope='relative correction history at x=2; not absolute exterior closure')


def independent_leading_power(field):
    g=field.control.ranges.phase.built['graph'];built=field.control.functions['exact_C1_limit_adapter']
    z,T=s.symbols('Z T');mu=s.Symbol('mu',positive=True);A=s.Function('Am')(z)
    bindings={built['parameters']['mu'].node:mu}
    x=field.functions['original_band_variable'];bindings[g.unary('log',x).node]=T
    amplitude=field.control.ranges.target.functions['actual_terminal_amplitude']
    for order,key in enumerate(('value','Z','ZZ')):bindings[getattr(amplitude,key).node]=s.diff(A,z,order)
    leading=field.functions['original_leading_power'];seeds={key:s.Function(key+'0')(z) for key in current.current.RATES}
    for key,q in leading['actual_endpoint_seeds'].items():
        for order,row in enumerate(('value','Z','ZZ')):bindings[getattr(q,row).node]=s.diff(seeds[key],z,order)
    expected=dict(m=seeds['m']*s.exp(-T),k=seeds['k']*s.exp(-3*T/2),
        h=s.exp(-3*T/2)*(seeds['h']+A*(s.exp((1-mu)*T)-1)/(1-mu)),
        e=s.exp(-T)*(seeds['e']-A*A*(1-s.exp(-2*mu*T))/(4*mu)),
        p=seeds['p']+A*A*(1-s.exp(-(1+2*mu)*T))/(2*(1+2*mu)))
    for key,q in leading['leading_histories'].items():
        for order,handle in enumerate((q.value,q.Z,q.ZZ)):
            assert s.simplify(calculus.symbolic(g,handle.node,bindings)-s.diff(expected[key],z,order))==0,(key,order)
    return dict(independent_original_power_semigroup_C2_rows=15,
        original_endpoint_histories_and_amplitude_cross_terms_retained=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C2_repair_band_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC2RepairBand(require_checked=False)
        assert not field.acceptance_loaded and field.control.acceptance_loaded
        g=field.control.ranges.phase.built['graph'];functions=field.functions
        encoded=lambda q:json.loads(json.dumps(control.encoded(current.ranges.record(q))))
        parent=json.loads(gzip.decompress((current.HERE/control.NAME).read_bytes()))
        assert g.nodes[:len(parent['exact_graph_nodes'])]==parent['exact_graph_nodes']
        assert g.nodes[:len(field.prefix)]==field.prefix and g.nodes==raw['exact_graph_nodes']
        assert encoded(functions)==raw['actual_C2_repair_band_functions'] and raw['source_family']==field.identity
        assert raw['actual_selected_repair_integer']==field.control.bounds['actual_same_repair_integer']
        profiles=independent_profiles(field);density=independent_density(field)
        bump=independent_bump_calculus();endpoints=independent_endpoints(field);leading=independent_leading_power(field)
        bindings=dict(profiles_C2=current.current.ast_binding(current.profiles_C2),
            signed_density_C2=current.current.ast_binding(current.signed_density_C2),
            leading_C2=current.current.ast_binding(current.leading_C2),build=current.current.ast_binding(current.build),
            actual_endpoint_projection=current.current.ast_binding(current.CurrentC2RepairBand.endpoint_source),
            actual_reserved_source_query=current.current.ast_binding(current.CurrentC2RepairBand.leading_range),
            quantitative_range=current.current.ast_binding(current.CurrentC2RepairBand.quantitative_range),
            original_C1_profiles=current.current.ast_binding(current.limit.band.profiles_at),
            original_C1_signed_density=current.current.ast_binding(current.limit.band.density_pairs),
            original_C1_power_background=current.current.ast_binding(current.limit.original_power_graph),
            original_reserved_power_source=current.limit.BACKGROUND_BINDING)
        assert raw['source_bindings']==bindings
        built=field.control.functions['exact_C1_limit_adapter'];old=functions['original_C1_background']
        for key,q in functions['partial_correction_histories'].items():
            prior=built['partial_band_histories'][key];assert q.value==prior.value and q.Z==prior.Z
            assert q.ZZ!=g.zero and functions['actual_incoming_correction_C2'][key].ZZ!=g.zero
            n=g.nodes[functions['relative_terminal_zero_certificates'][key].node]
            assert n['relative_not_absolute_exterior'] and n['source_family']==field.identity
            assert n['equivalent_residual_rows']==control.encoded(functions['relative_terminal_residual_rhs'][key])
            assert n['actual_second_zero_theorem_nodes']==current.current.encode_graph(field.control.functions['actual_second_zero_theorem_nodes'])
            t=g.symbol(functions['partial_integration_variable']);x=functions['original_band_variable']
            rate=current.current.RATES[key]
            kernel=g.unary('exp',g.mul(g.constant(rate-1),g.unary('log',t)))
            decay=g.unary('exp',g.neg(g.mul(g.constant(rate),g.unary('log',x))))
            term=functions['partial_correction_contributions'][key].ZZ
            arguments=g.nodes[term.node]['arguments'] if g.nodes[term.node]['operation']=='product' else [term.node]
            integral_nodes=[g.nodes[i] for i in arguments if g.nodes[i]['operation']=='definite_integral']
            assert len(integral_nodes)==1 and decay.node in arguments
            integral=integral_nodes[0]
            assert integral['integrand']==g.mul(kernel,functions['signed_density'][key].ZZ).node
            assert integral['lower']==g.one.node and integral['upper']==x.node
            assert integral['ordinary_slow_Z_derivative_order']==2 and integral['original_physical_Jacobian_applied_once']
        for key,q in functions['complete_histories'].items():
            assert q.value==old['complete_histories'][key].value and q.Z==old['complete_histories'][key].Z
        P0=functions['original_leading_power']['independent_P0'];prior=old['independent_original_P0']
        assert P0.value==prior.value and P0.Z==prior.Z
        assert g.nodes[P0.ZZ.node]['quantity']=='P0' and g.nodes[P0.ZZ.node]['Z_order']==2
        assert g.nodes[P0.ZZ.node]['Taylor_coefficient_factorial']==2
        target=json.loads(gzip.decompress((current.HERE/current.ranges.NAME).read_bytes()))
        counts=dict(actual_Z_cells=0,actual_endpoint_history_rows=0,actual_endpoint_P0_rows=0,
            original_incoming_second_memories_retained=5,relative_second_terminal_identities=5)
        for saved,parent_row in zip(raw['actual_four_Z_C2_band_ranges'],target['actual_four_Z_C2_target_transports']):
            ends=tuple(saved['exact_Z_cell']);replay=field.quantitative_range(ends,parent_row)
            assert encoded(replay)==saved
            endpoint=field.endpoint_source(ends);op=field.control.ranges.phase.outer.owner.owner(ends);f=op.flow
            inlet=field.source.query(ends,(1,1));actual=inlet['raw']['raw_current_radius_y_derivative_axial_coefficients']
            for key,rows in actual['histories'].items():
                normalized=[q*f.factor((0,-.5,0,0,0)) for q in rows[0]] if key in ('m','k') else rows[0]
                want=[normalized[0],normalized[1],2*normalized[2]]
                assert current.current.current.previous.equivalent_rows(want,endpoint['ordinary_leading_history_C2'][key])
                counts['actual_endpoint_history_rows']+=3
            for q,want in zip(endpoint['ordinary_independent_P0_C2'],(op.P0[0],op.P0[1],2*op.P0[2])):
                assert q.ctx is field.c and q.scale.bases is f.logs and q.ledger is f.ledger
                assert current.current.current.previous.equivalent_rows([q],[want])
                counts['actual_endpoint_P0_rows']+=1
            assert saved['actual_whole_band_leading_source']['same_original_P0']
            assert all(q.zero for rows in field.leading_range(ends)['actual_whole_band_V_y_C2'] for q in rows)
            assert saved['range_caps_not_function_values'] and saved['selected_repair_integer_unchanged']
            counts['actual_Z_cells']+=1
        assert counts['actual_Z_cells']==4 and counts['actual_endpoint_history_rows']==60 and counts['actual_endpoint_P0_rows']==12
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            C2_repair_band_fields_and_histories_installed=True,
            independent_mixed_profile_calculus=profiles,independent_signed_density=density,
            independent_original_beta=bump,independent_relative_endpoints=endpoints,
            independent_original_power_semigroup=leading,actual_replay_counts=counts,
            actual_C1_function_handles_and_independent_P0_retained=True,
            C2_Rh_join_and_absolute_exterior_installed=False,current_numeric_point_field_oracle_installed=False,
            **dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(control.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual C2 repaired-band mixed velocity, signed histories, relative terminal and range checks passed',flush=True)
    return result


if __name__=='__main__':run()
