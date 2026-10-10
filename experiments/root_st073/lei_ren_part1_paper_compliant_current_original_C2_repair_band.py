"""Actual same-limit C2 repaired power band, histories and relative closure.

The original reserved x in[1,2] is used, with t=2+log(x). Signed function
handles and directed source/range witnesses remain separate throughout.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C2_limit_controls as control

phase,current,limit=control.phase,control.current,control.limit
ranges,source,controls=control.ranges,control.source,control.controls
HERE,PREFIX,sha,ep=control.HERE,control.PREFIX,control.sha,control.ep
NAME=PREFIX+'current_original_C2_repair_band.json.gz'
RECEIPT=PREFIX+'current_original_C2_repair_band_check.json'
GATES=('current_original_same_limit_repair_band_C2_velocity_histories_installed',
    'current_original_same_limit_repair_band_C2_relative_terminal_identities_installed',
    'current_original_same_limit_repair_band_C2_directed_magnitude_ranges_installed')
OPEN=control.OPEN


def preserve(pair,second):return phase.C2Function(pair.value,pair.Z,second)


def profiles_C2(field,x):
    built=field.control.functions['exact_C1_limit_adapter'];g=built['graph'];W=built['repair_weights']
    alg=phase.C2Algebra(g);old=limit.profiles_at_limit(built,x)
    h=field.control.control_functions();A=field.control.ranges.target.functions['actual_terminal_amplitude']
    logx=g.unary('log',x);alpha=g.add(g.constant('1/2'),W['mu'])
    power=g.unary('exp',g.neg(g.mul(alpha,logx)))
    overN=lambda v:g.quotient(v,built['N'],'one common original positive integer N')
    bump_y=[];bump_yy=[]
    for center,bump in zip(W['centers'],old['bumps']):
        arg=g.quotient(g.sub(logx,center),W['ell'],'original positive log bump width')
        prime=g.node('compact_raw_beta_derivative',argument=arg.node,
            definition='-2*t*exp(-1/(1-t^2))/(1-t^2)^2 for abs(t)<1; zero otherwise',
            defining_module=PREFIX+'outer_pulse_map.py',defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'),
            lazy_outside_support=True)
        second=g.node('compact_raw_beta_second_derivative',argument=arg.node,
            definition='(6*t^4-2)*exp(-1/(1-t^2))/(1-t^2)^4 for abs(t)<1; zero otherwise',
            original_derivative_source=current.ast_binding(phase.flat_source.beta_jets),
            original_polynomial_source=current.ast_binding(phase.flat_source.beta_polynomials),
            ordinary_beta_derivative_order=2,lazy_outside_support=True,
            derivative_of_range_endpoint=False)
        prime_term=g.quotient(prime,g.mul(W['ell'],W['ell'],W['normal'],x),
            'original positive ell^2 J0 x')
        second_term=g.quotient(second,g.mul(W['ell'],W['ell'],W['ell'],W['normal'],x),
            'original positive ell^3 J0 x')
        bump_y.append(g.sub(prime_term,bump))
        bump_yy.append(g.add(second_term,g.neg(g.mul(g.constant(2),prime_term)),bump))
    def correction(bumps,indices):
        return phase.C2Function(*[overN(g.add(*(g.mul(bumps[j],getattr(h[i],row)) for j,i in indices)))
            for row in ('value','Z','ZZ')])
    rows=[]
    for j,bumps in enumerate((old['bumps'],bump_y,bump_yy)):
        F=correction(bumps,((0,2),(1,3),(2,4)))
        G=correction(bumps,((0,0),(2,1)))
        original=alg.mul(A,alg.fixed(g.mul(power,g.constant((-1)**j),*[alpha]*j)))
        dE,dV=alg.mul(A,F),alg.mul(A,G)
        row=dict(F=F,G=G,original_E=original,delta_E=dE,delta_V=dV,E=alg.add(original,dE),V=dV)
        if j==0:
            row={key:preserve(old[key],q.ZZ) for key,q in row.items()}
        rows.append(row)
    return dict(radial_y_rows=rows,bumps=old['bumps'],bump_y=bump_y,bump_yy=bump_yy,
        original_power=power,original_alpha=alpha,ordinary_radial_coordinate='y=log(x)',
        compact_support_all_endpoint_jets_zero=True,amplitude_Z2_not_dropped=True)


def signed_density_C2(g,E,dE,dV):
    alg=phase.C2Algebra(g);EF,FF,VV=alg.mul(E,dE),alg.mul(dE,dE),alg.mul(dV,dV)
    exact=dict(m=dV,h=dE,k=alg.mul(alg.add(E,dE),dV),
        e=alg.add(VV,alg.neg(EF),alg.neg(alg.scale(FF,'1/2'))),
        p=alg.add(EF,alg.scale(FF,'1/2')))
    old=limit.band.density_pairs(g,source.C1Function(E.value,E.Z),
        source.C1Function(dE.value,dE.Z),source.C1Function(dV.value,dV.Z))
    return {key:preserve(old[key],q.ZZ) for key,q in exact.items()}


def leading_C2(field,old_background,x):
    built=field.control.functions['exact_C1_limit_adapter'];g=built['graph'];alg=phase.C2Algebra(g)
    A=field.control.ranges.target.functions['actual_terminal_amplitude'];AA=alg.mul(A,A)
    mu=built['parameters']['mu'];t=g.unary('log',x)
    decay=lambda r:g.unary('exp',g.neg(g.mul(g.constant(r),t)))
    mass=lambda k:g.mul(t,g.unary('exprel',g.neg(g.mul(k,t))))
    def second_leaf(pair):
        node=dict(g.nodes[pair.value.node]);node['Z_order']=2
        node.update(ordinary_slow_Z_derivative_order=2,Taylor_coefficient_factorial=2,
            source_projection_binding=current.ast_binding(CurrentC2RepairBand.endpoint_source),
            source_projection='ordinary ZZ=2*same normalized endpoint Taylor coefficient[2]',
            derivative_of_range_endpoint=False)
        return preserve(pair,g.node(node.pop('operation'),**node))
    seeds={key:second_leaf(q) for key,q in old_background['original_Rc_leading_histories'].items()}
    P0=second_leaf(old_background['independent_original_P0'])
    theta=g.mul(decay('3/2'),t,g.unary('exprel',g.mul(g.sub(g.one,mu),t)))
    energy=g.mul(decay(1),mass(g.mul(g.constant(2),mu)))
    pressure=mass(g.add(g.one,g.mul(g.constant(2),mu)))
    histories=dict(m=alg.scale(seeds['m'],decay(1)),
        h=alg.add(alg.scale(seeds['h'],decay('3/2')),alg.scale(A,theta)),
        k=alg.scale(seeds['k'],decay('3/2')),
        e=alg.add(alg.scale(seeds['e'],decay(1)),alg.scale(AA,g.mul(g.constant('-1/2'),energy))),
        p=alg.add(seeds['p'],alg.scale(AA,g.mul(g.constant('1/2'),pressure))))
    histories={key:preserve(old_background['leading_histories'][key],q.ZZ) for key,q in histories.items()}
    return dict(actual_endpoint_seeds=seeds,leading_histories=histories,independent_P0=P0,
        original_Rc_source_coordinate=2,reserved_power_offset=g.add(g.constant(2),t),
        original_radius_offset=old_background['original_band_radius_offset'],
        m_k_S_normalization_already_in_endpoint_source=True,P0_stays_separate=True)


def build(field):
    built=field.control.functions['exact_C1_limit_adapter'];g=built['graph'];alg=phase.C2Algebra(g)
    old_background=limit.original_power_graph(built,field.control.ranges.phase.report)
    x=g.symbol(built['band_variable']);variable=built['partial_integration_variable'];t=g.symbol(variable)
    at_x=profiles_C2(field,x);at_t=profiles_C2(field,t)
    # Retain accepted existing first-y fields from the same original source.
    prior=old_background['corrected_fields'];radial=at_x['radial_y_rows'][1]
    for key in ('F','G','E','V'):radial[key]=preserve(prior[key+'_y'],radial[key].ZZ)
    density=signed_density_C2(g,*[at_t['radial_y_rows'][0][key] for key in ('original_E','delta_E','delta_V')])
    invN=g.quotient(g.one,built['N'],'same selected exact positive integer N')
    incoming={key:preserve(built['history'][key],g.mul(invN,q.ZZ))
        for key,q in field.control.ranges.target.functions['terminal_N_scaled_histories'].items()}
    partial={};contributions={}
    for key,rate in current.RATES.items():
        kernel=g.unary('exp',g.mul(g.constant(rate-1),g.unary('log',t)))
        decay=g.unary('exp',g.neg(g.mul(g.constant(rate),g.unary('log',x))))
        integral=controls.integral(g,g.mul(kernel,density[key].ZZ),variable,g.one,x,
            measure='dt; exact t^(rate-1) original recovery weight',original_recovery_rate=str(rate),
            partial_endpoint='repair_x in[1,2]',endpoint_independent_of_Z=True,
            ordinary_slow_Z_derivative_order=2,original_physical_Jacobian_applied_once=True)
        contributions[key]=preserve(built['partial_band_contributions'][key],g.mul(decay,integral))
        partial[key]=preserve(built['partial_band_histories'][key],g.mul(decay,g.add(incoming[key].ZZ,integral)))
    leading=leading_C2(field,old_background,x)
    complete={key:alg.add(leading['leading_histories'][key],partial[key]) for key in current.RATES}
    complete={key:preserve(old_background['complete_histories'][key],q.ZZ) for key,q in complete.items()}
    E,V=at_x['radial_y_rows'][0]['E'],at_x['radial_y_rows'][0]['V'];EE=alg.mul(E,E)
    dy=dict(m=alg.add(alg.scale(complete['m'],-1),V),
        h=alg.add(alg.scale(complete['h'],'-3/2'),E),
        k=alg.add(alg.scale(complete['k'],'-3/2'),alg.mul(E,V)),
        e=alg.add(alg.scale(complete['e'],-1),alg.mul(V,V),alg.neg(alg.scale(EE,'1/2'))),
        p=alg.scale(EE,'1/2'))
    dy={key:preserve(old_background['complete_history_y_Z_pairs'][key],q.ZZ) for key,q in dy.items()}
    R={row:preserve(pair,second) for row,pair,second in zip(controls.ROWS,built['control_residual'],
        field.control.functions['actual_second_control_residual'])}
    A=field.control.ranges.target.functions['actual_terminal_amplitude'];AA=alg.mul(A,A)
    endpoint_rhs={}
    for key,row,amplitude in (('m','M',A),('h','I',A),('e','S',AA),('p','Cp',AA)):
        decay=g.unary('exp',g.neg(g.mul(g.constant(current.RATES[key]),g.unary('log',g.constant(2)))))
        exact=alg.scale(alg.mul(amplitude,R[row]),g.mul(decay,invN))
        endpoint_rhs[key]=preserve(built['exact_limit_terminal_residual_identities'][key],exact.ZZ)
    joint=alg.add(R['M'],alg.scale(R[controls.ROWS[1]],built['parameters']['mu']))
    decay=g.unary('exp',g.neg(g.mul(g.constant(current.RATES['k']),g.unary('log',g.constant(2)))))
    exact=alg.scale(alg.mul(AA,joint),g.mul(decay,invN))
    endpoint_rhs['k']=preserve(built['exact_limit_terminal_residual_identities']['k'],exact.ZZ)
    endpoint_actual={key:phase.C2Function(*[g.node('function_substitution',expression=getattr(q,row).node,
        variable=x.node,value=g.constant(2).node,Z_independent_substitution=True) for row in ('value','Z','ZZ')])
        for key,q in partial.items()}
    certificates={key:g.node('proved_C2_relative_terminal_zero_identity',
        expression_rows=control.encoded(endpoint_actual[key]),equivalent_residual_rows=control.encoded(endpoint_rhs[key]),
        actual_C1_zero_theorem_nodes=current.encode_graph(built['residual_zero_theorem_nodes']),
        actual_second_zero_theorem_nodes=current.encode_graph(field.control.functions['actual_second_zero_theorem_nodes']),
        theorem='same original signed bump integrals and own-rate endpoint identity, differentiated twice on the same C2 limit',
        zero_rows=[g.zero.node]*3,relative_not_absolute_exterior=True,source_family=field.identity)
        for key in current.RATES}
    return dict(profiles_at_band_x=at_x,profiles_at_integration_t=at_t,signed_density=density,
        actual_incoming_correction_C2=incoming,partial_correction_contributions=contributions,
        partial_correction_histories=partial,original_leading_power=leading,
        complete_histories=complete,complete_history_y=dy,control_residual_C2=R,
        relative_terminal_actual=endpoint_actual,relative_terminal_residual_rhs=endpoint_rhs,
        relative_terminal_zero_certificates=certificates,
        original_C1_background=old_background,original_band_variable=x,
        partial_integration_variable=variable,actual_reserved_band=[1,2],
        radial_y_orders=[0,1,2],axial_Z_orders=[0,1,2],
        existing_C1_profiles_histories_and_P0_retained=True,
        absolute_exterior_pressure_heat_and_Rh_C2_join_not_admitted=True)


class CurrentC2RepairBand:
    def __init__(self,control_field=None,owner=None,require_checked=True):
        self.control=control_field if control_field is not None else control.CurrentC2LimitControls(owner=owner)
        if not self.control.acceptance_loaded:raise ValueError('Accepted genuine C2 controls required')
        self.identity=self.control.identity;self.c=self.control.ranges.c;self.hashes=dict(self.control.hashes)
        self.source=limit.ActualReservedPowerSource(self.control.ranges.phase.outer)
        self.cache={};self.prefix=[dict(n) for n in self.control.ranges.phase.built['graph'].nodes]
        for name in (control.NAME,control.RECEIPT,PREFIX+'flat_pulse_derivatives_check.json',Path(__file__).name):
            self.hashes[name]=sha(name)
        beta=json.loads((HERE/(PREFIX+'flat_pulse_derivatives_check.json')).read_bytes())
        if not beta['all_passed'] or not beta['analytic_original_beta_and_flat_envelope_checks']['beta_direct_derivative_2']:
            raise ValueError('Checked original beta second derivative required')
        self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C2 repair band required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual C2 band source '+name)
            self.acceptance_loaded=True

    def endpoint_source(self,ends):
        """Same already-normalized endpoint histories, ordinary Z2 once."""
        ends=tuple(ends);phase_field=self.control.ranges.phase
        if ends not in current.current.CELLS:raise ValueError('Admitted actual Z cell required')
        key=('endpoint',ends)
        if key not in self.cache:
            packet=phase_field.leading_source(ends,'O3_power',(2,1))
            op=phase_field.outer.owner.owner(ends);f=op.flow;generic=packet['original_generic_source']
            if generic['common_original_P0_axial5'] is not op.P0:raise ValueError('Same actual independent P0 required')
            def project(row):
                if len(row)<3 or any(q.ctx is not self.c or q.scale.bases is not f.logs or q.ledger is not f.ledger for q in row):
                    raise ValueError('Actual Taylor rows in same context/basis/ledger required')
                return [row[0],row[1],2*row[2]]
            self.cache[key]=dict(source_family=self.identity,exact_Z_cell=ends,
                ordinary_leading_history_C2={key:project(row) for key,row in generic['common_own_five_histories_axial5'].items()},
                ordinary_independent_P0_C2=project(op.P0),
                ordinary_terminal_amplitude_C2=project(generic['common_velocity_E_axial5']),
                original_Rc_source_coordinate=2,reserved_band_left_coordinate=1,
                m_k_S_normalization_already_applied=True,second_Taylor_factorial_applied_once=True,
                range_witness_not_function_value=True)
        return self.cache[key]

    def leading_range(self,ends):
        """Independent actual full reserved band, with m/k source units once."""
        ends=tuple(ends);key=('whole_band',ends)
        if key not in self.cache:
            packet=self.source.query(ends,(1,1),(2,1));op=self.control.ranges.phase.outer.owner.owner(ends);f=op.flow
            raw=packet['raw']['raw_current_radius_y_derivative_axial_coefficients'];invS=f.factor((0,-.5,0,0,0))
            project=lambda row:[row[0],row[1],2*row[2]]
            histories={key:project([q*invS for q in rows[0]]) if key in ('m','k') else project(rows[0])
                for key,rows in raw['histories'].items()}
            self.cache[key]=dict(actual_whole_reserved_band_leading_history_C2=histories,
                actual_whole_band_E_y_C2=[project(row) for row in raw['velocity']['theta']],
                actual_whole_band_V_y_C2=[project([q*invS for q in row]) for row in raw['velocity']['axial']],
                actual_independent_P0_C2=project(op.P0),source_family=self.identity,
                exact_Z_cell=ends,actual_band_x=[1,2],same_original_P0=packet['exact_common_P0_axial5'] is op.P0,
                actual_power_offset='t in[2,2+log2]',original_reservation=self.source.reservation,
                source_binding=limit.BACKGROUND_BINDING,m_k_S_normalization_applied_once=True,
                second_Taylor_factorial_applied_once=True,range_witness_not_function_value=True)
        return self.cache[key]

    def quantitative_range(self,ends,parent_row):
        c=self.c;bd=ranges.Bounds(c);read=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(read(row) for row in rows))
        source=self.endpoint_source(ends);whole=self.leading_range(ends)
        A=ranges.JetBound(*(bd.row(q) for q in source['ordinary_terminal_amplitude_C2']))
        bounds=self.control.bounds;rho=read(bounds['original_C1_ball_radius_upper']);H2=read(bounds['actual_limit_and_iterate_ZZ_upper'])
        H=ranges.JetBound(rho,rho,H2);epsilon=read(bounds['actual_inverse_epsilon_upper'])
        freq=self.control.ranges.phase.report['actual_whole_Z_frequency_connection'];weights=freq['fresh_exact_integral_weights']
        ell=current.packets.interval(c,weights['radius']);normal=current.packets.interval(c,weights['raw_normalization'])
        if ep(ell)[0]<=0 or ep(normal)[0]<=0:raise ValueError('Actual positive original ell/J0 required')
        beta=[c.exp(-1),8*c.exp(-2),2048*c.exp(-4)]
        shape=[bd.constant(beta[0]/(ell*normal)),
            bd.sum(bd.constant(beta[1]/(ell**2*normal)),bd.constant(beta[0]/(ell*normal))),
            bd.sum(bd.constant(beta[2]/(ell**3*normal)),bd.constant(2*beta[1]/(ell**2*normal)),bd.constant(beta[0]/(ell*normal)))]
        mu=c.exp(self.control.ranges.phase.outer.logmu);alpha=c.mpf('.5')+mu
        if ep(alpha)[1]>=ep(c.mpf(2)/3)[0]:raise ValueError('Original power alpha<2/3 required')
        radial=[]
        for j,cap in enumerate(shape):
            F=bd.scaled(H,bd.scale(cap*epsilon,3));G=bd.scaled(H,bd.scale(cap*epsilon,2))
            original=bd.scaled(A,bd.constant(alpha**j));dE,dV=bd.product(A,F),bd.product(A,G)
            radial.append(dict(F=F,G=G,original_E=original,delta_E=dE,delta_V=dV,E=bd.add(original,dE),V=dV))
        row=radial[0];E,dE,dV=[row[key] for key in ('original_E','delta_E','delta_V')]
        EF,FF,VV=bd.product(E,dE),bd.product(dE,dE),bd.product(dV,dV)
        density=dict(m=dV,h=dE,k=bd.product(bd.add(E,dE),dV),
            e=bd.add(VV,EF,bd.scaled(FF,bd.constant(c.mpf('.5')))),
            p=bd.add(EF,bd.scaled(FF,bd.constant(c.mpf('.5')))))
        incoming={key:bd.scaled(jet(rows),epsilon) for key,rows in parent_row['actual_normalized_N_scaled_histories_C2'].items()}
        histories={};masses={};complete={}
        for key,rate in current.RATES.items():
            r=c.mpf(rate.numerator)/rate.denominator
            mass=c.ln(2) if not rate else -c.expm1(-r*c.ln(2))/r
            if ep(mass)[0]<=0:raise ValueError('Original own positive band mass required')
            histories[key]=bd.add(incoming[key],bd.scaled(density[key],bd.constant(mass)));masses[key]=mass
            leading=ranges.JetBound(*(bd.row(q) for q in whole['actual_whole_reserved_band_leading_history_C2'][key]))
            complete[key]=bd.add(leading,histories[key])
        return dict(source_family=self.identity,exact_Z_cell=ends,band_x=[1,2],
            actual_endpoint_source=control.encoded(source),actual_whole_band_leading_source=control.encoded(whole),
            actual_selected_N_profile_y0_y1_y2_C2_bounds=ranges.record(radial),
            actual_signed_density_magnitude_C2_bounds=ranges.record(density),
            actual_incoming_correction_C2_bounds=ranges.record(incoming),
            actual_own_rate_positive_masses=masses,actual_partial_correction_C2_bounds=ranges.record(histories),
            actual_complete_history_C2_bounds=ranges.record(complete),
            original_beta_derivative_global_caps=beta,original_bump_y0_y1_y2_caps=ranges.record(shape),
            source_derivative_units_and_independent_P0_retained=True,selected_repair_integer_unchanged=True,
            range_caps_not_function_values=True,absolute_exterior_heat_Rh_C2_and_physical_oracle_not_admitted=True)

    def velocity_functions(self):return self.functions['profiles_at_band_x']['radial_y_rows']
    def history_functions(self):return self.functions['complete_histories']


def run(control_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC2RepairBand(control_field=control_field,owner=owner,require_checked=False)
        parent=json.loads(gzip.decompress((HERE/ranges.NAME).read_bytes()))
        rows=[field.quantitative_range(tuple(row['exact_Z_cell']),row) for row in parent['actual_four_Z_C2_target_transports']]
        bindings=dict(profiles_C2=current.ast_binding(profiles_C2),signed_density_C2=current.ast_binding(signed_density_C2),
            leading_C2=current.ast_binding(leading_C2),build=current.ast_binding(build),
            actual_endpoint_projection=current.ast_binding(CurrentC2RepairBand.endpoint_source),
            actual_reserved_source_query=current.ast_binding(CurrentC2RepairBand.leading_range),
            quantitative_range=current.ast_binding(CurrentC2RepairBand.quantitative_range),
            original_C1_profiles=current.ast_binding(limit.band.profiles_at),
            original_C1_signed_density=current.ast_binding(limit.band.density_pairs),
            original_C1_power_background=current.ast_binding(limit.original_power_graph),
            original_reserved_power_source=limit.BACKGROUND_BINDING)
        report=dict(candidate_actual_C2_repair_band_constructed=True,source_family=field.identity,
            original_C2_control_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.control.ranges.phase.built['graph'].nodes,
            actual_C2_repair_band_functions=control.encoded(field.functions),actual_four_Z_C2_band_ranges=ranges.record(rows),
            source_bindings=bindings,actual_selected_repair_integer=field.control.bounds['actual_same_repair_integer'],
            C2_Rh_join_and_absolute_exterior_installed=False,current_numeric_point_field_oracle_installed=False,
            **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(control.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual same-limit C2 repaired-band velocities, histories, relative endpoints and directed ranges constructed',flush=True)
    return field


if __name__=='__main__':run()
