"""Current 17-chart C1 limit on the actual reserved original power band.

The exact source-bound Banach functions are mathematical function handles,
not numerical controls or range endpoints. Relative terminal cancellation
does not close the absolute exterior/heat/pressure problem.
"""
from fractions import Fraction
import ast
import copy
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection as current
import lei_ren_part1_paper_compliant_current_native_repair_band_cone_budget as band

source,controls,packets=current.source,current.controls,current.packets
outer=current.current
original=outer.o3.original
HERE,PREFIX,sha,ep,require=current.HERE,current.PREFIX,current.sha,current.ep,current.require
NAME=PREFIX+'current_original_whole_Z_limit_repair_band.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_limit_repair_band_check.json'
GATE='current_original_17_chart_source_bound_C1_limit_actual_power_band_relative_five_moments_closed'


def band_fraction(value):
    if not isinstance(value,tuple) or len(value)!=2 or any(type(q) is not int for q in value) or value[1]<=0:
        raise ValueError('Exact rational repair x coordinate required')
    q=Fraction(*value)
    if not 1<=q<=2:raise ValueError('Reserved original power band x in[1,2] required')
    return q


def extended_background():
    """Only adapt the coordinate prefix; retain original background arithmetic.

    The historical power[0,2] guard stays strict. This separate API uses
    x in[1,2] and t=2+log(x), inside the same original Tw reservation.
    """
    fn=ast.parse(textwrap.dedent(inspect.getsource(original.background_cell))).body[0]
    split=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='key' for t in n.targets))
    suffix=copy.deepcopy(fn.body[split:])
    prefix=ast.parse('''lo,hi=band_fraction(left),band_fraction(right)
if hi<lo:raise ValueError('Ordered reserved original power x cell required')
if chart!='power':raise ValueError('Only reserved original power continuation is extended')
f,c=op.flow,op.c
cv=lambda q:c.mpf(q.numerator)/q.denominator
t=c.mpf([ep(2+c.ln(cv(lo)))[0],ep(2+c.ln(cv(hi)))[1]])
''').body
    binding=dict(original_background=current.ast_binding(original.background_cell),
        unchanged_original_background_suffix_AST_sha256=hashlib.sha256(ast.dump(ast.Module(
            body=suffix,type_ignores=[])).encode()).hexdigest(),
        original_suffix_first_assignment='key',only_coordinate_prefix_adapted=True,
        original_domain_guard_unchanged=True,new_domain='x in[1,2]; t=2+log(x)',
        original_leading_histories_P0_radial_and_Z_arithmetic_unchanged=True,
        arithmetic_mu_cover_remains_an_enclosure_not_a_function_parameter=True)
    fn.name='reserved_original_power_background';fn.body=prefix+suffix
    namespace={**vars(original),'band_fraction':band_fraction}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original power background on reserved x band>','exec'),namespace)
    return namespace[fn.name],binding


BACKGROUND,BACKGROUND_BINDING=extended_background()


class ActualReservedPowerSource:
    """Actual leading-source range access, separate from exact graph values."""
    def __init__(self,owner):
        if type(owner) is not outer.WholeZAllNOuterRcFunctions:
            raise TypeError('Same checked current all-N outer source owner required')
        self.owner,self.family=owner,owner.identity
        logmu,Tw,margin=outer.o3.source_mu_admission(owner.c,owner.owner.source.eta_log)
        require(logmu._mpi_==owner.logmu._mpi_,'Same original formal mu required')
        self.reservation=dict(actual_source_log_mu=logmu,actual_Tw=Tw,
            actual_quiet_power_log_margin=margin,
            reserved_band_upper_log_offset=2+owner.c.ln(2),
            original_Tw_strictly_exceeds_band=True)

    def query(self,ends,left,right=None):
        right=left if right is None else right
        lo,hi=band_fraction(left),band_fraction(right)
        if tuple(ends) not in outer.CELLS or hi<lo:raise ValueError('Admitted Z and ordered x cell required')
        op=self.owner.owner.owner(tuple(ends))
        parent=self.owner.leading_packet(tuple(ends),'O3_transition',(1,1))
        result=BACKGROUND(op,parent,'power',left,right,self.owner.logmu)
        raw=result['raw']['raw_current_radius_y_derivative_axial_coefficients']
        require(result['raw']['original_P0_normalized_axial5'] is op.P0,'Original independent P0 identity lost')
        require(all(q.zero for rows in raw['velocity']['axial'] for q in rows),'Original power V=Vy=0 required')
        result.update(source_family=self.family,fixed_leading_input_N0=self.owner.N0,
            exact_Z_cell=list(ends),exact_x_cell=[list(left),list(right)],
            exact_common_P0_axial5=op.P0,actual_phase_Z_exact_zero=True,
            same_original_transition_parent_retained=True,source_binding=BACKGROUND_BINDING,
            original_reservation=self.reservation,finite_N_density_and_support_not_called=True)
        return result


def load_current():
    report=json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
    checked=json.loads((HERE/current.RECEIPT).read_bytes())
    require(report[current.GATE] and checked[current.GATE] and checked['all_passed'],
        'Accepted current 17-chart control connection required')
    require(report['source_family']==checked['source_family'],'Same current source family required')
    hashes=dict(checked['input_hashes'])
    for name,digest in hashes.items():require(sha(name)==digest,'Changed accepted source: '+name)
    hashes[current.NAME]=sha(current.NAME);hashes[current.RECEIPT]=sha(current.RECEIPT)
    for module in (current,band,controls,source,original):
        name=Path(module.__file__).name;hashes[name]=sha(name)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    return report,hashes


def restore_current(report):
    """Restore exact handles without replaying source integrations or owners."""
    g=source.FunctionTransportGraph();nodes=report['exact_function_graph_nodes']
    require(nodes[:2]==g.nodes,'Exact original zero/one graph prefix required')
    g.nodes=[dict(row) for row in nodes]
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    require(len(g.keys)==len(g.nodes),'Duplicate current expression identities')
    raw=report['exact_current_integral_and_control_graph']
    ref=lambda i:source.FunctionRef(g,i)
    pair=lambda q:source.C1Function(*map(ref,q))
    return dict(graph=g,N=ref(raw['N']),source_family=report['source_family'],
        parameters={k:ref(v) for k,v in raw['parameters'].items()},amplitude=pair(raw['amplitude']),
        N_scaled_targets={k:pair(v) for k,v in raw['N_scaled_targets'].items()},
        history={k:pair(v) for k,v in raw['history'].items()},
        normalized_history={k:pair(v) for k,v in raw['normalized_history'].items()},
        original_Rc_offset=ref(raw['windows'][-1]['right_offset']),
        original_graph_prefix_length=len(nodes),band_variable='repair_x')


def control_equations(built,h,W):
    g,N=built['graph'],built['N'];d=[built['N_scaled_targets'][key] for key in controls.ROWS]
    overN=lambda q:g.quotient(q,N,'same selected exact positive integer N')
    Q=controls.quadratic(g,W,[q.value for q in h])
    QZ=controls.quadratic_action(g,W,[q.value for q in h],[q.Z for q in h])
    next_h=controls.c1_vectors(g,
        [g.neg(q) for q in controls.inverse_action(g,W,[g.add(q.value,overN(r)) for q,r in zip(d,Q)])],
        [g.neg(q) for q in controls.inverse_action(g,W,[g.add(q.Z,overN(r)) for q,r in zip(d,QZ)])])
    residual=controls.c1_vectors(g,
        [g.add(q.value,g.add(*(g.mul(b,v.value) for b,v in zip(row,h))),overN(r)) for q,row,r in zip(d,W['B'],Q)],
        [g.add(q.Z,g.add(*(g.mul(b,v.Z) for b,v in zip(row,h))),overN(r)) for q,row,r in zip(d,W['B'],QZ)])
    columns=[controls.quadratic_action(g,W,[q.value for q in h],[g.one if j==i else g.zero for j in range(5)]) for i in range(5)]
    jac=[[g.add(W['B'][i][j],overN(columns[j][i])) for j in range(5)] for i in range(5)]
    return dict(map=next_h,residual=residual,implicit_Z_matrix=jac,implicit_Z_rhs=[g.neg(q.Z) for q in d])


def exact_limit_adapter(built,report):
    g=built['graph'];W=controls.exact_weights(g,built['parameters']['mu'])
    generic=[source.C1Function(g.symbol('repair_control_'+str(i)),g.symbol('repair_control_Z_'+str(i))) for i in range(5)]
    template=control_equations(built,generic,W)
    proof=report['actual_whole_Z_frequency_connection']
    require(proof['C1_limit']['exact_C1_limit_exists_unique_for_this_repair_only_integer'],
        'Actual source-bound C1 limit certificate required')
    vector=g.node('source_bound_C1_Banach_limit',
        variable_pairs=current.encode_graph(generic),map_pairs=current.encode_graph(template['map']),
        initial_pairs=[[g.zero.node,g.zero.node]]*5,
        exact_target_pairs=current.encode_graph(built['N_scaled_targets']),shared_N=built['N'].node,
        source_family=built['source_family'],global_Z_domain=[-1,1],
        source_certificate=current.RECEIPT,source_certificate_sha256=sha(current.RECEIPT),
        source_connection_report=current.NAME,source_connection_report_sha256=sha(current.NAME),
        selected_integer=proof['selected_integer'],C1_limit_proof=proof['C1_limit'],
        definition='unique C1 limit of this vector map from h0=0 in the certified original ball',
        finite_iterate_not_used_as_limit=True,numerical_value_not_installed=True)
    values=[g.node('C1_Banach_limit_component',limit_vector=vector.node,
        component=i,Z_order=0,genuine_C1_derivative=False) for i in range(5)]
    provisional=control_equations(built,[source.C1Function(q,g.zero) for q in values],W)
    linear_system=g.node('implicit_C1_Z_linear_system',
        matrix=current.encode_graph(provisional['implicit_Z_matrix']),
        rhs=current.encode_graph(provisional['implicit_Z_rhs']),
        limit_vector=vector.node,definition='(B+DQ(h*)/N)*h*_Z=-d_Z',
        invertibility_proof='B^-1*(B+DQ/N)=I+B^-1*DQ/N; norm of latter perturbation <=1/2 from same C1 contraction',
        C1_derivative_proof='differentiate actual C1 fixed point; Z-independent B,mu,N and exact target Z pairs',
        source_certificate_sha256=sha(current.RECEIPT),numerical_value_not_installed=True)
    h=[source.C1Function(q,g.node('C1_Banach_limit_component',limit_vector=vector.node,
        component=i,Z_order=1,genuine_C1_derivative=True,implicit_Z_linear_system=linear_system.node))
        for i,q in enumerate(values)]
    equations=control_equations(built,h,W)
    zeros=[g.node('proved_C1_zero_identity',expression_pair=current.encode_graph(q),
        theorem='Banach fixed point for actual same source target and selected repair integer',
        limit_vector=vector.node,zero_pair=[g.zero.node,g.zero.node],row=row)
        for row,q in zip(controls.ROWS,equations['residual'])]
    # The legacy pure band helper selects the last sequence element. Supply
    # exactly the h* functions, retain their actual residual expressions,
    # and label this adapter explicitly so no finite iterate is promoted.
    out=dict(**built,repair_weights=W,finite_picard_sequence=[h],
        control_residual=equations['residual'],implicit_Z_matrix=equations['implicit_Z_matrix'],
        implicit_Z_rhs=equations['implicit_Z_rhs'],exact_limit_vector=vector,exact_limit_controls=h,
        implicit_Z_linear_system=linear_system,
        generic_control_variables=generic,generic_control_map=template['map'],
        residual_zero_theorem_nodes=zeros,band_control_binding='exact source-bound C1 limit, not finite Picard')
    out=band.exact_partial_band_graph(out)
    out['exact_limit_terminal_residual_identities']=out.pop('finite_terminal_residual_identities')
    out.pop('finite_picard_sequence')
    return out


def profiles_at_limit(built,x):
    """Private compatibility adapter; exported controls are only exact h*."""
    return band.profiles_at(dict(**built,finite_picard_sequence=[built['exact_limit_controls']]),x)


def leading_history_binding(report):
    recipe=copy.deepcopy(report['exact_source_recipe_manifest']['O3_power'])
    recovery=outer.o3.reference.long
    recipe['quantity_paths'].update({**{('history_'+key):
        'original_generic_source.common_own_five_histories_axial5.'+key+'[Z_order]'
        for key in current.RATES},'P0':'original_generic_source.common_original_P0_axial5[Z_order]'})
    recipe['actual_original_frontend_binding']=outer.FRONTEND_BINDINGS['o3']
    recipe['original_recovery_binding']=dict(compiler=current.ast_binding(recovery.compile_recovery),
        original_recovery=current.ast_binding(recovery.generic.recover_inputs),
        unchanged_signed_recovery_binding=recovery.RECOVERY_BINDING,
        actual_histories='m,k: raw value rows / original S; h,e,p: unchanged raw value rows',
        normalization='invS=f.factor((0,-.5,0,0,0)); same original fixed source bases',
        P0='packet.original_P0_normalized_axial5 is op.P0; returned common_original_P0_axial5 is op.P0')
    recipe['exact_source_call']='WholeZAllNOuterRcFunctions.leading_packet(ends, O3_power, (2,1))'
    recipe['exact_Rc_coordinate_identity']='Rw*exp(2)=Rc; reserved x=1 gives t=2+log1=2'
    recipe['source_family']=report['source_family']
    recipe['fixed_leading_input_N0']=report['fixed_leading_input_N0']
    return recipe


def original_power_graph(built,report):
    """Continue every leading history and its Z jet; P0 stays independent."""
    g=built['graph'];mu=built['parameters']['mu'];x=g.symbol(built['band_variable'])
    t=g.unary('log',x);alpha=g.add(g.constant('1/2'),mu)
    decay=lambda r:g.unary('exp',g.neg(g.mul(g.constant(r),t)))
    mass=lambda k:g.mul(t,g.unary('exprel',g.neg(g.mul(k,t))))
    A=built['amplitude'];AA=g.c1mul(A,A);zero=source.C1Function(g.zero,g.zero)
    recipe=leading_history_binding(report)
    def rc_leaf(quantity,order):
        return g.node('current_original_leading_function_recipe',source_family=built['source_family'],
            native_chart='O3_power',recipe=recipe,quantity=quantity,Z_order=order,
            coordinate=g.constant(2).node,Z_variable=g.symbol('Z').node,
            defining_quantity_not_a_range_value=True,
            quantity_path=(['original_generic_source','common_own_five_histories_axial5',quantity[8:]]
                if quantity.startswith('history_') else ['original_generic_source','common_original_P0_axial5']))
    old={key:source.C1Function(rc_leaf('history_'+key,0),rc_leaf('history_'+key,1)) for key in current.RATES}
    P0=source.C1Function(rc_leaf('P0',0),rc_leaf('P0',1))
    E=g.c1mul(A,source.C1Function(g.unary('exp',g.neg(g.mul(alpha,t))),g.zero))
    theta=g.mul(decay(Fraction(3,2)),t,g.unary('exprel',g.mul(g.sub(g.one,mu),t)))
    energy=g.mul(decay(1),mass(g.mul(g.constant(2),mu)))
    pressure=mass(g.add(g.one,g.mul(g.constant(2),mu)))
    history=dict(m=g.c1scale(decay(1),old['m']),
        h=g.c1add(g.c1scale(decay(Fraction(3,2)),old['h']),g.c1scale(theta,A)),
        k=g.c1scale(decay(Fraction(3,2)),old['k']),
        e=g.c1add(g.c1scale(decay(1),old['e']),g.c1scale(g.mul(g.constant('-1/2'),energy),AA)),
        p=g.c1add(old['p'],g.c1scale(g.mul(g.constant('1/2'),pressure),AA)))
    EE=g.c1mul(E,E)
    dy=dict(m=g.c1scale(g.constant(-1),history['m']),
        h=g.c1add(E,g.c1scale(g.constant('-3/2'),history['h'])),
        k=g.c1scale(g.constant('-3/2'),history['k']),
        e=g.c1add(g.c1scale(g.constant(-1),history['e']),g.c1scale(g.constant('-1/2'),EE)),
        p=g.c1scale(g.constant('1/2'),EE))
    fields=profiles_at_limit(built,x)
    W=built['repair_weights'];controls_h=built['exact_limit_controls']
    bump_y=[]
    for center,bump in zip(W['centers'],fields['bumps']):
        arg=g.quotient(g.sub(t,center),W['ell'],'original positive log bump width')
        beta_y=g.node('compact_raw_beta_derivative',argument=arg.node,
            definition='-2*t*exp(-1/(1-t^2))/(1-t^2)^2 for abs(t)<1; zero otherwise',
            defining_module=PREFIX+'outer_pulse_map.py',defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'),
            lazy_outside_support=True)
        bump_y.append(g.sub(g.quotient(beta_y,g.mul(W['ell'],W['ell'],W['normal'],x),
            'original positive ell^2 J0 x'),bump))
    overN=lambda q:g.quotient(q,built['N'],'same selected exact positive integer N')
    Fy=source.C1Function(*[overN(g.add(*(g.mul(b,getattr(controls_h[i+2],row)) for i,b in enumerate(bump_y)))) for row in ('value','Z')])
    Gy=source.C1Function(*[overN(g.add(g.mul(bump_y[0],getattr(controls_h[0],row)),
        g.mul(bump_y[2],getattr(controls_h[1],row)))) for row in ('value','Z')])
    fields.update(F_y=Fy,G_y=Gy,E_y=g.c1add(g.c1scale(g.neg(alpha),fields['original_E']),g.c1mul(A,Fy)),V_y=g.c1mul(A,Gy))
    complete={key:g.c1add(history[key],built['partial_band_histories'][key]) for key in current.RATES}
    Ecomplete,Vcomplete=fields['E'],fields['V'];EC2=g.c1mul(Ecomplete,Ecomplete)
    dy_complete=dict(m=g.c1add(g.c1scale(g.constant(-1),complete['m']),Vcomplete),
        h=g.c1add(g.c1scale(g.constant('-3/2'),complete['h']),Ecomplete),
        k=g.c1add(g.c1scale(g.constant('-3/2'),complete['k']),g.c1mul(Ecomplete,Vcomplete)),
        e=g.c1add(g.c1scale(g.constant(-1),complete['e']),g.c1add(g.c1mul(Vcomplete,Vcomplete),g.c1scale(g.constant('-1/2'),EC2))),
        p=g.c1scale(g.constant('1/2'),EC2))
    return dict(original_source_background_binding=BACKGROUND_BINDING,leading_history_P0_source_binding=recipe,
        original_Rc_leading_histories=old,independent_original_P0=P0,
        local_original_power_offset=g.add(g.constant(2),t),
        original_band_radius_offset=g.add(built['original_Rc_offset'],t),
        E=E,V=zero,E_y=g.c1scale(g.neg(alpha),E),V_y=zero,
        leading_histories=history,leading_history_y_Z_pairs=dy,
        corrected_fields=fields,complete_histories=complete,complete_history_y_Z_pairs=dy_complete,
        P0_not_reset_or_merged_into_p_history=True,
        exact_original_formal_mu_used_in_functions=True,
        source_range_enclosures_not_used_as_graph_values=True)


def symbolic_theorems():
    """Original ODE/semigroup and signed endpoint map, including genuine Z."""
    mu,t,s=sy.symbols('mu t s',positive=True);z=sy.symbols('Z',real=True)
    A=sy.Function('A')(z);old={key:sy.Function(key+'0')(z) for key in current.RATES}
    mass=lambda k,u:(1-sy.exp(-k*u))/k
    def flow(u,seed,amp):
        return dict(m=seed['m']*sy.exp(-u),
            h=seed['h']*sy.exp(-3*u/2)+amp*sy.exp(-3*u/2)*(sy.exp((1-mu)*u)-1)/(1-mu),
            k=seed['k']*sy.exp(-3*u/2),
            e=seed['e']*sy.exp(-u)-amp**2*sy.exp(-u)*mass(2*mu,u)/2,
            p=seed['p']+amp**2*mass(1+2*mu,u)/2)
    h=flow(t,old,A);E=A*sy.exp(-(sy.Rational(1,2)+mu)*t)
    densities=dict(m=0,h=E,k=0,e=-E**2/2,p=E**2/2)
    for key,rate in current.RATES.items():
        require(sy.simplify(sy.diff(h[key],t)+sy.Rational(rate.numerator,rate.denominator)*h[key]-densities[key])==0,
            'Original leading history ODE failed: '+key)
        require(sy.simplify(sy.diff(h[key].subs(t,0)-old[key],z))==0,'Original inlet Z jet failed')
    advanced=flow(s,h,A*sy.exp(-(sy.Rational(1,2)+mu)*t))
    for key in current.RATES:
        require(sy.simplify(advanced[key]-h[key].subs(t,t+s))==0,'Original power semigroup failed: '+key)
    # Five original bump integrals are independent formal exact symbols.
    N=sy.symbols('N',positive=True);alpha=sy.Rational(1,2)+mu
    a0,a2,e0,e1,e2=[sy.Function(k)(z) for k in ('a0','a2','e0','e1','e2')]
    D0,D2=sy.symbols('D0 D2');cross=sy.symbols('cross0 cross2')
    I=sy.symbols('I0:3');S=sy.symbols('S0:3');P=sy.symbols('P0:3')
    dM,dD,dI,dS,dP=[sy.Function(k)(z) for k in ('dM','dD','dI','dS','dP')]
    ee=(e0,e1,e2);joint=cross[0]*a0*e0+cross[1]*a2*e2
    Qe=sy.symbols('Qe0:3');Qp=sy.symbols('Qp0:3')
    energy_quadratic=Qe[0]*a0*a0+Qe[2]*a2*a2-sum(v*e*e for v,e in zip(Qe,ee))/2
    pressure_quadratic=sum(v*e*e for v,e in zip(Qp,ee))/2
    res=[dM+a0+a2,dD+D0*a0+D2*a2+joint/(mu*N),
        dI+sum(i*e for i,e in zip(I,ee)),
        dS+sum(v*e for v,e in zip(S,ee))+energy_quadratic/N,
        dP+sum(v*e for v,e in zip(P,ee))+pressure_quadratic/N]
    endpoint=dict(m=A*(dM+a0+a2)/(2*N),
        h=2**(-sy.Rational(3,2))*A*(dI+sum(v*e for v,e in zip(I,ee)))/N,
        k=2**(-sy.Rational(3,2))*A*A*(dM+mu*dD+(1+mu*D0)*a0+(1+mu*D2)*a2+joint/N)/N,
        e=(A*A*dS/N+A*A*sum(v*e for v,e in zip(S,ee))/N+A*A*energy_quadratic/N**2)/2,
        p=A*A*dP/N+A*A*sum(v*e for v,e in zip(P,ee))/N+A*A*pressure_quadratic/N**2)
    expected=dict(m=A*res[0]/(2*N),h=2**(-sy.Rational(3,2))*A*res[2]/N,
        k=2**(-sy.Rational(3,2))*A*A*(res[0]+mu*res[1])/N,
        e=A*A*res[3]/(2*N),p=A*A*res[4]/N)
    for key in current.RATES:
        require(sy.expand(endpoint[key]-expected[key])==0,'Signed endpoint identity failed: '+key)
        require(sy.expand(sy.diff(endpoint[key]-expected[key],z))==0,'Genuine endpoint Z identity failed: '+key)
    L=sy.log(2);supports=[(q*L-L/40,q*L+L/40) for q in (sy.Rational(1,5),sy.Rational(1,2),sy.Rational(4,5))]
    require(all(sy.simplify(a)>0 and sy.simplify(L-b)>0 for a,b in supports),'Compact supports reach band ends')
    require(all(sy.simplify(supports[i+1][0]-supports[i][1])>0 for i in range(2)),'Compact supports overlap')
    return dict(original_five_leading_history_ODEs=True,original_five_history_semigroup=True,
        genuine_inlet_and_terminal_Z_product_rules=True,signed_five_endpoint_identities=True,
        original_joint_row_reconstructed_as_res_M_plus_mu_res_D=True,
        original_pressure_rate_zero_memory_retained=True,
        exact_bump_log_supports=[['7*log2/40','9*log2/40'],['19*log2/40','21*log2/40'],['31*log2/40','33*log2/40']],
        correction_fields_vanish_on_endpoint_neighborhoods_with_all_radial_derivatives=True,
        limits_agree_with_leading_background_only_at_right_terminal_histories=True,
        incoming_correction_histories_not_reset_at_left=True,
        absolute_physical_exterior_or_heat_closure_not_implied=True)


def build(report):
    built=exact_limit_adapter(restore_current(report),report)
    built['original_power_background']=original_power_graph(built,report)
    g=built['graph'];x=g.symbol(built['band_variable']);at2=g.constant(2)
    endpoints={key:source.C1Function(*[g.node('function_substitution',expression=getattr(q,row).node,
        variable=x.node,value=at2.node,Z_independent_substitution=True) for row in ('value','Z')])
        for key,q in built['partial_band_histories'].items()}
    rhs=built['exact_limit_terminal_residual_identities']
    built['relative_terminal_Duhamel_functions']=endpoints
    built['relative_terminal_zero_certificates']={key:g.node('proved_C1_zero_identity',
        expression_pair=current.encode_graph(endpoints[key]),
        equivalent_residual_pair=current.encode_graph(rhs[key]),
        theorem='original signed density/bump integrals + own-rate Duhamel + exact current C1 limit',
        residual_zero_theorem_nodes=[q.node for q in built['residual_zero_theorem_nodes']],
        zero_pair=[g.zero.node,g.zero.node],relative_not_absolute_exterior=True) for key in current.RATES}
    return built


def run():
    began=time.monotonic();report,hashes=load_current();theorems=symbolic_theorems();built=build(report)
    result=dict(**{GATE:True},source_family=report['source_family'],
        fixed_leading_input_N0=report['fixed_leading_input_N0'],
        selected_repair_integer=report['actual_whole_Z_frequency_connection']['selected_integer'],
        exact_graph_nodes=built['graph'].nodes,exact_limit_repair_band=current.encode_graph(built),
        original_power_background_source_binding=BACKGROUND_BINDING,symbolic_theorems=theorems,
        actual_power_repair_band_source_extended=True,
        exact_source_bound_C1_limit_band_functions_defined=True,
        relative_terminal_five_moment_C0_Z_function_identities_proved=True,
        endpoint_correction_field_jets_vanish_by_original_flat_supports=True,
        current17_graph_and228_source_integrations_reused=True,
        actual_numeric_controls_evaluated=False,actual_global_frequency_admitted=False,
        physical_original_exterior_five_targets_closed=False,
        full_recovered_velocity_pressure_and_heat_joins_admitted=False,
        actual_temporal_scale_recursion_installed=False,**dict.fromkeys(outer.OPEN,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Exact current17 C1 Banach limit, original reserved Rc..2Rc power/leading histories/P0, '
            'compact-bump fields and relative five-moment terminal identities in C1(Z). '
            'No numeric field/control oracle, absolute exterior/heat, global frequency or temporal recursion admission.')
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(packets.encode(result),separators=(',',':'))+'\n').encode(),mtime=0))
    print('CURRENT_SOURCE_LIMIT_REPAIR_BAND',len(built['graph'].nodes),'nodes; relative terminal C0/Z identities',flush=True)
    return result


if __name__=='__main__':run()
