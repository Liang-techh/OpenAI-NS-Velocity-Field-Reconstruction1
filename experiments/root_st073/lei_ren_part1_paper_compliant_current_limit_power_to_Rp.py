"""Same repaired original source continued through the quiet power to Rp.

Complete histories and P0 are transported, never reset. The exact pulse
input frame is exported; identifying every old native inlet constant with
this new frame remains a separate source-function bridge requirement.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_limit_band_physical_recovery as current
import lei_ren_part1_paper_compliant_actual_Rp_source_join as native
import lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 as pre
import lei_ren_part1_paper_compliant_power_inlet_C4 as inlet

source,packets=current.source,current.packets
HERE,PREFIX,sha,require=current.HERE,current.PREFIX,current.sha,current.require
NAME=PREFIX+'current_limit_power_to_Rp.json.gz'
RECEIPT=PREFIX+'current_limit_power_to_Rp_check.json'
GATE='current_repaired_original_C1_absolute_history_quiet_power_to_Rp_and_pulse_frame_defined'


def load_inputs():
    report=json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
    checked=json.loads((HERE/current.RECEIPT).read_bytes())
    require(report[current.GATE] and checked[current.GATE] and checked['all_passed'],
        'Accepted current repaired physical source required')
    hashes=dict(checked['input_hashes'])
    old_name=PREFIX+'actual_Rp_source_join_check.json'
    old=json.loads((HERE/old_name).read_bytes())
    require(old['all_passed'] and old['current_Rp_external_pulse_join_certified'],
        'Accepted original native Rp recipe required')
    for key,value in report['source_family'].items():
        require(old[key]==value,'Same original native source family required: '+key)
    for name,digest in old['input_hashes'].items():
        require(name not in hashes or hashes[name]==digest,'Current/native source hash conflict: '+name)
        hashes[name]=digest
    for module in (current,native,pre,inlet):
        name=Path(module.__file__).name;hashes[name]=sha(name)
    hashes.update({current.NAME:sha(current.NAME),current.RECEIPT:sha(current.RECEIPT),old_name:sha(old_name)})
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,digest in hashes.items():require(sha(name)==digest,'Changed prerequisite: '+name)
    return report,old,hashes


def restore(report):
    g=source.FunctionTransportGraph();nodes=report['exact_graph_nodes']
    require(nodes[:2]==g.nodes,'Exact same original source graph required')
    g.nodes=[dict(row) for row in nodes]
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    require(len(g.keys)==len(g.nodes),'Duplicate accepted graph expressions')
    raw=report['exact_physical_recovery'];ref=lambda q:source.FunctionRef(g,q)
    pair=lambda q:source.C1Function(*map(ref,q))
    return dict(graph=g,source_family=report['source_family'],source_graph_prefix_length=len(nodes),
        parameters={k:ref(v) for k,v in raw['parameters'].items()},
        S=ref(raw['S']),delta=ref(raw['delta']),Rc=ref(raw['Rc']),N=ref(raw['N']),
        own={k:pair(v) for k,v in raw['own'].items()},fields={k:pair(v) for k,v in raw['fields'].items()},
        original_own={k:pair(v) for k,v in raw['original_leading'].items()},
        original_fields={k:pair(v) for k,v in raw['original_fields'].items()},P0=pair(raw['P0']),
        relative_terminal_zero_certificates=raw['relative_terminal_zero_certificates'],
        original_leaf_binding=raw['original_leaf_binding'],
        physical_coordinates={k:ref(v) for k,v in raw['physical_coordinate_functions'].items() if k!='domain'},
        final_heat_reference=raw['final_heat_reference'])


def positive_c1_quotient(g,left,right,certificate):
    return source.C1Function(g.quotient(left.value,right.value,certificate),
        g.quotient(g.sub(g.mul(left.Z,right.value),g.mul(left.value,right.Z)),
            g.mul(right.value,right.value),certificate))


def power_flow(g,A,own,mu,s):
    """Exact original all-five semigroup, with true ordinary-Z pairs."""
    alpha=g.add(g.constant('1/2'),mu)
    decay=lambda r:g.unary('exp',g.neg(g.mul(g.constant(r),s)))
    mass=lambda k:g.mul(s,g.unary('exprel',g.neg(g.mul(k,s))))
    E=g.c1scale(g.unary('exp',g.neg(g.mul(alpha,s))),A)
    AA=g.c1mul(A,A);zero=source.C1Function(g.zero,g.zero)
    theta=g.mul(decay(Fraction(3,2)),s,g.unary('exprel',g.mul(g.sub(g.one,mu),s)))
    hist=dict(m=g.c1scale(decay(1),own['m']),
        h=g.c1add(g.c1scale(decay(Fraction(3,2)),own['h']),g.c1scale(theta,A)),
        k=g.c1scale(decay(Fraction(3,2)),own['k']),
        e=g.c1add(g.c1scale(decay(1),own['e']),g.c1scale(g.mul(g.constant('-1/2'),decay(1),mass(g.mul(g.constant(2),mu))),AA)),
        p=g.c1add(own['p'],g.c1scale(g.mul(g.constant('1/2'),mass(g.add(g.one,g.mul(g.constant(2),mu)))),AA)))
    EE=g.c1mul(E,E)
    dy=dict(m=g.c1scale(g.constant(-1),hist['m']),h=g.c1add(E,g.c1scale(g.constant('-3/2'),hist['h'])),
        k=g.c1scale(g.constant('-3/2'),hist['k']),
        e=g.c1add(g.c1scale(g.constant(-1),hist['e']),g.c1scale(g.constant('-1/2'),EE)),
        p=g.c1scale(g.constant('1/2'),EE))
    return dict(fields=dict(E=E,V=zero,E_y=g.c1scale(g.neg(alpha),E),V_y=zero),own=hist,own_y=dy)


def pulse_input(g,E,own,P0):
    """Common-unit histories -> actual native Utheta-based pulse units.

    m,k already include 1/S in the accepted current common frame. Dividing
    them by S again would erase the original axial/joint histories.
    """
    E2=g.c1mul(E,E);positive='Same original strictly positive power swirl at Rp'
    return dict(u=E,m1=positive_c1_quotient(g,own['m'],E,positive),
        m2=positive_c1_quotient(g,own['k'],E2,positive),
        X=positive_c1_quotient(g,own['h'],E,positive),
        energy=positive_c1_quotient(g,own['e'],E2,positive),Mp=own['p'],P0=P0,
        pressure=g.c1add(P0,own['p']),
        definitions=dict(u='Utheta/Pstar',m1='Mz/(R*Utheta)',
            m2='J/(sqrt(2)*R^(3/2)*Utheta^2)',X='I/(sqrt(2)*R^(3/2)*Utheta)',
            energy='S_energy/(R*Utheta^2)',Mp='Cp/Pstar^2',P0='Pi_axis/Pstar^2'),
        common_m_k_already_divided_by_S=True,independent_original_P0_preserved=True,
        original_C1_quotient_rule_includes_amplitude_Z=True)


def build(report):
    built=restore(report);g=built['graph'];x=g.symbol('repair_x');Z=g.symbol('Z');s=g.symbol('post_2Rc_power_s')
    mu=built['parameters']['mu'];logmu=g.unary('log',mu);Tw=g.mul(g.constant(-60),logmu)
    length=g.sub(g.sub(Tw,g.constant(2)),g.unary('log',g.constant(2)))
    at2=lambda pair:source.C1Function(*[current.substitution(g,getattr(pair,row),[x],[g.constant(2)]) for row in ('value','Z')])
    own2={k:at2(v) for k,v in built['own'].items()};A2=at2(built['fields']['E'])
    flow=power_flow(g,A2,own2,mu,s)
    original2={k:at2(v) for k,v in built['original_own'].items()}
    original_flow=power_flow(g,at2(built['original_fields']['E']),original2,mu,s)
    endpoint=lambda pair:source.C1Function(*[current.substitution(g,getattr(pair,row),[s],[length]) for row in ('value','Z')])
    terminal={k:endpoint(v) for k,v in flow['own'].items()};Ep=endpoint(flow['fields']['E'])
    R2=g.mul(g.constant(2),built['Rc']);R=g.mul(R2,g.unary('exp',s));Rp=g.mul(R2,g.unary('exp',length))
    state=current.native_recovery(g,flow['fields'],flow['own'],flow['own_y'],built['P0'],Z,built['delta'])
    paper=current.paper_moments(g,flow['own'],built['P0'],R,built['S'])
    native_frame=pulse_input(g,Ep,terminal,built['P0'])
    physical=built['physical_coordinates'];lam=physical['lambda_scale'];r=physical['r']
    sp=g.unary('log',g.quotient(physical['x_band'],g.constant(2),'positive physical radius and same Rc'))
    at=lambda q:current.substitution(g,q,[Z,s],[physical['Z'],sp])
    nu=g.symbol('constant_viscosity');px=g.symbol('physical_x');py=g.symbol('physical_y')
    power=lambda q:g.unary('exp',g.mul(q,g.unary('log',lam)))
    Ur=g.mul(g.unary('sqrt',nu),power(g.constant(-1)),built['S'],
        g.unary('sqrt',g.mul(g.constant('1/2'),physical['R'])),at(state['Q']))
    Ut=g.mul(g.unary('sqrt',nu),power(g.neg(g.add(g.one,built['delta']))),built['S'],at(state['E'].value))
    cosine=g.quotient(px,r,'quiet power r>0');sine=g.quotient(py,r,'quiet power r>0')
    cartesian=dict(u=g.sub(g.mul(cosine,Ur),g.mul(sine,Ut)),
        v=g.add(g.mul(sine,Ur),g.mul(cosine,Ut)),w=g.zero,
        p=g.mul(nu,built['S'],built['S'],power(g.neg(g.mul(g.constant(2),g.add(g.one,built['delta'])))),at(state['pressure'].value)))
    physical_domain=g.node('source_region_guard',tau=physical['tau'].node,nu=nu.node,r=r.node,
        native_coordinate=sp.node,length=length.node,definition='tau>0,nu>0,r>0,0<=post_2Rc_power_s<=length',
        source_family=built['source_family'],outside_quiet_power_requires_other_actual_region=True)
    continuity=g.node('original_quiet_power_C1_relative_zero_propagation_theorem',
        inlet_certificate_nodes=built['relative_terminal_zero_certificates'],
        complete_inlet={k:g.ids([q.value,q.Z]) for k,q in own2.items()},
        leading_inlet={k:g.ids([q.value,q.Z]) for k,q in original2.items()},
        complete_E=g.ids([A2.value,A2.Z]),leading_E=g.ids([at2(built['original_fields']['E']).value,at2(built['original_fields']['E']).Z]),
        exact_own_rates=['1','3/2','3/2','1','0'],same_profiles_on_whole_quiet_power=True,
        correction_equations='D_y=-rate*D; D(2Rc)=D_Z(2Rc)=0; density differences exactly zero',
        conclusion='All five complete/leading histories and first Z rows agree on [2Rc,Rp]',
        pressure_constant_not_reset=True,does_not_set_absolute_histories_to_zero=True)
    domain=g.node('original_power_continuation_domain',coordinate=s.node,length=length.node,
        definition='0<=s<=Tw-2-log(2); s=log(R/(2Rc)); same original Z in[-1,1]',
        geometry='Rw=Rc*exp(-2); Rp=Rw*exp(Tw)',
        positive_width_binding=current.current.current.ast_binding(current.current.outer.o3.source_mu_admission),
        original_positive_mu_and_quiet_power_source=True,
        derivative_coordinate='ordinary y=logR; ds/dy=1, ds/dZ=0')
    built.update(Tw=Tw,length=length,R2=R2,R=R,Rp=Rp,s=s,domain=domain,A2=A2,incoming_at_2Rc=own2,
        continuation=flow,leading_continuation=original_flow,relative_zero_propagation=continuity,
        recovery=state,paper_moments=paper,terminal_Rp_histories=terminal,terminal_Rp_E=Ep,
        terminal_Rp_paper=current.paper_moments(g,terminal,built['P0'],Rp,built['S']),
        pulse_input_frame=native_frame,
        cartesian_velocity_pressure=cartesian,physical_domain=physical_domain,
        native_bridge_requirements=dict(original_Rp_power_recipe_bound=True,
            actual_Rp_native_frame_source_function_identified=False,
            exact_native_constructor_consumes_current_limit_frame=False,
            required_common_Rw_parent_identity='Identify Rc recipe original O3 parent with native SharedOuterBuffer.slope_mu(Z,1) as functions',
            required_constant_identities='u=U/q; m1=C1*Z*q; m2=C2*Z*q; X=Xp; energy=C0+C_E*Z^2*q^2; Mp=Pin/q^2',
            q='1+Z^2',same_family_and_range_overlap_are_insufficient=True,
            pressure_and_all_five_history_zero_reset_forbidden=True),
        source_recipe_bindings=dict(original_background=current.current.current.ast_binding(current.current.original.background_cell),
            original_pre_pulse_power=current.current.current.ast_binding(pre.CompliantPrePulseMixedC4.power),
            canonical_incoming=current.current.current.ast_binding(inlet.CompliantPowerInletC4.incoming),
            canonical_shape_proof=current.current.current.ast_binding(native.terminal_shape_proof)))
    return built


def symbolic_checks():
    s,t,mu=sy.symbols('s t mu',positive=True);z=sy.Symbol('Z',real=True)
    A=sy.Function('A')(z);old={k:sy.Function(k+'0')(z) for k in ('m','h','k','e','p')}
    def mass(k,s):return (1-sy.exp(-k*s))/k
    def flow(A,H,s):
        E=A*sy.exp(-(sy.Rational(1,2)+mu)*s)
        return E,dict(m=H['m']*sy.exp(-s),h=H['h']*sy.exp(-sy.Rational(3,2)*s)+A*(sy.exp(-(sy.Rational(1,2)+mu)*s)-sy.exp(-sy.Rational(3,2)*s))/(1-mu),
            k=H['k']*sy.exp(-sy.Rational(3,2)*s),e=sy.exp(-s)*(H['e']-A*A*mass(2*mu,s)/2),p=H['p']+A*A*mass(1+2*mu,s)/2)
    E,H=flow(A,old,s);E2,H2=flow(E,H,t);Eall,Hall=flow(A,old,s+t)
    drivers=dict(m=0,h=E,k=0,e=-E*E/2,p=E*E/2)
    rates=dict(m=1,h=sy.Rational(3,2),k=sy.Rational(3,2),e=1,p=0)
    for key in H:
        require(sy.simplify(sy.diff(H[key],s)+rates[key]*H[key]-drivers[key])==0,'Actual power moment ODE: '+key)
        require(sy.simplify(H[key].subs(s,0)-old[key])==0,'Actual complete power inlet: '+key)
        diff=sy.simplify(sy.expand_power_exp(H2[key]-Hall[key]))
        require(diff==0 and sy.diff(diff,z)==0,'Same-parent original semigroup split: '+key)
    require(sy.simplify(sy.expand_power_exp(E2-Eall))==0,'Same-parent original swirl split')
    # Canonical denominators use the actual C1 E, not an amplitude cap.
    S,R=sy.symbols('S R',positive=True);M,I,J,En,Cp,P0=sy.symbols('M I J En Cp P0')
    norm=dict(m=M/(R*S),h=I/(sy.sqrt(2)*R**sy.Rational(3,2)*S),
        k=J/(sy.sqrt(2)*R**sy.Rational(3,2)*S*S),e=En/(R*S*S),p=Cp/(S*S))
    actual=dict(m1=norm['m']/E,m2=norm['k']/E**2,X=norm['h']/E,energy=norm['e']/E**2,Mp=norm['p'],P0=P0)
    Ut=S*E;expected=dict(m1=M/(R*Ut),m2=J/(sy.sqrt(2)*R**sy.Rational(3,2)*Ut**2),
        X=I/(sy.sqrt(2)*R**sy.Rational(3,2)*Ut),energy=En/(R*Ut**2),Mp=Cp/(S*S),P0=P0)
    for key in expected:require(sy.simplify(actual[key]-expected[key])==0,'Original native pulse units: '+key)
    delta=sy.symbols('D_m D_h D_k D_e D_p')
    require(all(sy.diff(d*sy.exp(-rates[k]*s),s)+rates[k]*d*sy.exp(-rates[k]*s)==0 for k,d in zip(rates,delta)),
        'Zero relative corrections must propagate without erasing absolute state')
    Tw,Rc=sy.symbols('Tw Rc',positive=True)
    require(sy.simplify(sy.expand_power_exp(2*Rc*sy.exp(Tw-2-sy.log(2))-Rc*sy.exp(Tw-2)))==0,'Rp radius frame mismatch')
    return dict(all_five_complete_power_ODEs_and_inlets=True,same_parent_all_five_C1_semigroup=True,
        independent_native_pulse_unit_conversion=True,true_amplitude_Z_quotient_rule_required=True,
        homogeneous_relative_zero_propagation=True,absolute_history_zero_not_inferred=True,
        exact_Rp_Rw_Rc_radius_identity=True,original_native_parent_function_identification_still_required=True)


def run():
    began=time.monotonic();report,old,hashes=load_inputs();built=build(report)
    result=dict(**{GATE:True},source_family=report['source_family'],exact_graph_nodes=built['graph'].nodes,
        exact_power_to_Rp=current.current.current.encode_graph(built),symbolic_checks=symbolic_checks(),
        complete_absolute_C1_power_continuation_defined=True,actual_Rp_pulse_input_units_defined=True,
        current_original_source_recipe_and_native_Rp_receipt_bound=True,
        actual_Rp_native_frame_source_function_identified=False,
        actual_native_O4_constructor_consumes_current_limit_frame=False,
        actual_preheat_pressure_frame_identified=False,physical_original_exterior_five_targets_closed=False,
        actual_numeric_point_source_oracle_installed=False,actual_numeric_controls_evaluated=False,
        actual_global_frequency_admitted=False,full_recovered_velocity_pressure_and_heat_joins_admitted=False,
        higher_Z_velocity_jets_and_stress_cone_admitted=False,actual_temporal_scale_recursion_installed=False,
        **dict.fromkeys(current.current.outer.OPEN,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Same original repaired source absolute C1 histories and n=0 recovery continue over 2Rc..Rp; '
            'native pulse input units are exact. Old native constructor frame identification, O4/preheat/exterior/heat, '
            'actual numerics, high-Z/global stress/temporal recursion and full NS remain open.')
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(packets.encode(result),separators=(',',':'))+'\n').encode(),mtime=0))
    print('CURRENT_LIMIT_POWER_TO_RP',len(built['graph'].nodes),'nodes; all five absolute C1 histories and exact pulse units',flush=True)
    return result


if __name__=='__main__':run()
