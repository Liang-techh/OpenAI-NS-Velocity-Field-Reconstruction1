"""Same original repaired band: genuine C3 velocity and signed histories.

Every accepted C2 function handle is retained. Ordinary third source rows,
signed integral functions, relative terminal identities, and directed magnitude
ranges are separate. The reserved band is x in [1,2], t=2+log(x).
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C3_limit_controls as control
import lei_ren_part1_paper_compliant_current_original_C2_repair_band as lower

phase,current,limit=lower.phase,lower.current,lower.limit
ranges,target,controls=control.ranges,control.target,control.controls
HERE,PREFIX,sha,ep=control.HERE,control.PREFIX,control.sha,control.ep
NAME=PREFIX+'current_original_C3_repair_band.json.gz'
RECEIPT=PREFIX+'current_original_C3_repair_band_check.json'
GATES=('current_original_same_limit_repair_band_C3_velocity_histories_installed',
    'current_original_same_limit_repair_band_C3_relative_terminal_identities_installed',
    'current_original_same_limit_repair_band_C3_directed_magnitude_ranges_installed')
OPEN=control.OPEN


def profiles_C3(field,prior):
    built=field.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    g=built['graph'];alg=target.C3Algebra(g);N=built['N']
    h=field.control.control_functions();A=field.control.ranges.target.functions['actual_terminal_amplitude']
    rows=[]
    for j,(bumps,old) in enumerate(zip((prior['bumps'],prior['bump_y'],prior['bump_yy']),prior['radial_y_rows'])):
        def correction(key,indices):
            third=g.quotient(g.add(*(g.mul(bumps[b],h[i].ZZZ) for b,i in indices)),N,
                'one common original positive integer N')
            return target.preserve(old[key],third)
        F=correction('F',((0,2),(1,3),(2,4)));G=correction('G',((0,0),(2,1)))
        fixed=g.mul(prior['original_power'],g.constant((-1)**j),*[prior['original_alpha']]*j)
        original=alg.scale(A,fixed);dE,dV=alg.mul(A,F),alg.mul(A,G)
        exact=dict(F=F,G=G,original_E=original,delta_E=dE,delta_V=dV,E=alg.add(original,dE),V=dV)
        rows.append({key:target.preserve(old[key],q.ZZZ) for key,q in exact.items()})
    return {**prior,'radial_y_rows':rows,'amplitude_Z3_not_dropped':True,
        'accepted_C2_value_Z_ZZ_handles_retained':True}


def signed_density_C3(g,E,dE,dV,prior):
    alg=target.C3Algebra(g);EF,FF,VV=alg.mul(E,dE),alg.mul(dE,dE),alg.mul(dV,dV)
    exact=dict(m=dV,h=dE,k=alg.mul(alg.add(E,dE),dV),
        e=alg.add(VV,alg.neg(EF),alg.neg(alg.scale(FF,'1/2'))),
        p=alg.add(EF,alg.scale(FF,'1/2')))
    return {key:target.preserve(prior[key],q.ZZZ) for key,q in exact.items()}


def leading_C3(field,prior,x):
    built=field.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    g=built['graph'];alg=target.C3Algebra(g)
    A=field.control.ranges.target.functions['actual_terminal_amplitude'];AA=alg.mul(A,A)
    mu=built['parameters']['mu'];t=g.unary('log',x)
    decay=lambda r:g.unary('exp',g.neg(g.mul(g.constant(r),t)))
    mass=lambda k:g.mul(t,g.unary('exprel',g.neg(g.mul(k,t))))
    def third_leaf(q):
        node=dict(g.nodes[q.value.node]);node['Z_order']=3
        node.update(ordinary_slow_Z_derivative_order=3,Taylor_coefficient_factorial=6,
            source_projection_binding=current.ast_binding(CurrentC3RepairBand.endpoint_source),
            source_projection='ordinary ZZZ=6*same normalized endpoint Taylor coefficient[3]',
            derivative_of_range_endpoint=False)
        return target.preserve(q,g.node(node.pop('operation'),**node))
    seeds={key:third_leaf(q) for key,q in prior['actual_endpoint_seeds'].items()}
    P0=third_leaf(prior['independent_P0'])
    theta=g.mul(decay('3/2'),t,g.unary('exprel',g.mul(g.sub(g.one,mu),t)))
    energy=g.mul(decay(1),mass(g.mul(g.constant(2),mu)))
    pressure=mass(g.add(g.one,g.mul(g.constant(2),mu)))
    exact=dict(m=alg.scale(seeds['m'],decay(1)),
        h=alg.add(alg.scale(seeds['h'],decay('3/2')),alg.scale(A,theta)),
        k=alg.scale(seeds['k'],decay('3/2')),
        e=alg.add(alg.scale(seeds['e'],decay(1)),alg.scale(AA,g.mul(g.constant('-1/2'),energy))),
        p=alg.add(seeds['p'],alg.scale(AA,g.mul(g.constant('1/2'),pressure))))
    histories={key:target.preserve(prior['leading_histories'][key],q.ZZZ) for key,q in exact.items()}
    return {**prior,'actual_endpoint_seeds':seeds,'leading_histories':histories,'independent_P0':P0,
        'accepted_C2_value_Z_ZZ_handles_retained':True}


def build(field):
    parent=field.control.functions['exact_C2_limit_functions'];built=parent['exact_C1_limit_adapter']
    g=built['graph'];alg=target.C3Algebra(g);prior=field.lower.functions
    x=prior['original_band_variable'];variable=prior['partial_integration_variable'];t=g.symbol(variable)
    at_x=profiles_C3(field,prior['profiles_at_band_x']);at_t=profiles_C3(field,prior['profiles_at_integration_t'])
    density=signed_density_C3(g,*[at_t['radial_y_rows'][0][key] for key in ('original_E','delta_E','delta_V')],
        prior['signed_density'])
    invN=g.quotient(g.one,built['N'],'same selected exact positive integer N')
    incoming={key:target.preserve(prior['actual_incoming_correction_C2'][key],g.mul(invN,q.ZZZ))
        for key,q in field.control.ranges.target.functions['terminal_N_scaled_histories'].items()}
    partial={};contributions={}
    for key,rate in current.RATES.items():
        kernel=g.unary('exp',g.mul(g.constant(rate-1),g.unary('log',t)))
        decay=g.unary('exp',g.neg(g.mul(g.constant(rate),g.unary('log',x))))
        integral=controls.integral(g,g.mul(kernel,density[key].ZZZ),variable,g.one,x,
            measure='dt; exact t^(rate-1) original recovery weight',original_recovery_rate=str(rate),
            partial_endpoint='repair_x in[1,2]',endpoint_independent_of_Z=True,
            ordinary_slow_Z_derivative_order=3,original_physical_Jacobian_applied_once=True)
        contributions[key]=target.preserve(prior['partial_correction_contributions'][key],g.mul(decay,integral))
        partial[key]=target.preserve(prior['partial_correction_histories'][key],
            g.mul(decay,g.add(incoming[key].ZZZ,integral)))
    leading=leading_C3(field,prior['original_leading_power'],x)
    complete={key:target.preserve(prior['complete_histories'][key],
        alg.add(leading['leading_histories'][key],partial[key]).ZZZ) for key in current.RATES}
    E,V=at_x['radial_y_rows'][0]['E'],at_x['radial_y_rows'][0]['V'];EE=alg.mul(E,E)
    dy=dict(m=alg.add(alg.scale(complete['m'],-1),V),
        h=alg.add(alg.scale(complete['h'],'-3/2'),E),
        k=alg.add(alg.scale(complete['k'],'-3/2'),alg.mul(E,V)),
        e=alg.add(alg.scale(complete['e'],-1),alg.mul(V,V),alg.neg(alg.scale(EE,'1/2'))),
        p=alg.scale(EE,'1/2'))
    dy={key:target.preserve(prior['complete_history_y'][key],q.ZZZ) for key,q in dy.items()}
    R={row:target.preserve(prior['control_residual_C2'][row],third)
        for row,third in zip(controls.ROWS,field.control.functions['actual_third_control_residual'])}
    A=field.control.ranges.target.functions['actual_terminal_amplitude'];AA=alg.mul(A,A)
    rhs={}
    for key,row,amplitude in (('m','M',A),('h','I',A),('e','S',AA),('p','Cp',AA)):
        decay=g.unary('exp',g.neg(g.mul(g.constant(current.RATES[key]),g.unary('log',g.constant(2)))))
        exact=alg.scale(alg.mul(amplitude,R[row]),g.mul(decay,invN))
        rhs[key]=target.preserve(prior['relative_terminal_residual_rhs'][key],exact.ZZZ)
    joint=alg.add(R['M'],alg.scale(R[controls.ROWS[1]],built['parameters']['mu']))
    decay=g.unary('exp',g.neg(g.mul(g.constant(current.RATES['k']),g.unary('log',g.constant(2)))))
    exact=alg.scale(alg.mul(AA,joint),g.mul(decay,invN))
    rhs['k']=target.preserve(prior['relative_terminal_residual_rhs']['k'],exact.ZZZ)
    actual={key:target.preserve(prior['relative_terminal_actual'][key],
        g.node('function_substitution',expression=q.ZZZ.node,variable=x.node,value=g.constant(2).node,
            Z_independent_substitution=True)) for key,q in partial.items()}
    certificates={key:g.node('proved_C3_relative_terminal_zero_identity',
        expression_rows=control.encoded(actual[key]),equivalent_residual_rows=control.encoded(rhs[key]),
        actual_C2_relative_zero_theorem=prior['relative_terminal_zero_certificates'][key].node,
        actual_third_zero_theorem_nodes=current.encode_graph(field.control.functions['actual_third_zero_theorem_nodes']),
        theorem='same original signed bump integrals and own-rate endpoint identity differentiated thrice on the same C3 limit',
        zero_rows=[g.zero.node]*4,relative_not_absolute_exterior=True,source_family=field.identity)
        for key in current.RATES}
    return dict(profiles_at_band_x=at_x,profiles_at_integration_t=at_t,signed_density=density,
        actual_incoming_correction_C3=incoming,partial_correction_contributions=contributions,
        partial_correction_histories=partial,original_leading_power=leading,
        complete_histories=complete,complete_history_y=dy,control_residual_C3=R,
        relative_terminal_actual=actual,relative_terminal_residual_rhs=rhs,
        relative_terminal_zero_certificates=certificates,original_band_variable=x,
        partial_integration_variable=variable,actual_reserved_band=[1,2],
        radial_y_orders=[0,1,2],axial_Z_orders=[0,1,2,3],
        existing_C2_profiles_histories_and_P0_retained=True,
        next_generic_recovery_available_axial_order=2,
        radial_inertial_recovery_C2_not_yet_installed=True,
        absolute_exterior_pressure_heat_and_Rh_C2_join_not_admitted=True)


class CurrentC3RepairBand:
    def __init__(self,control_field=None,owner=None,require_checked=True):
        self.control=control_field if control_field is not None else control.CurrentC3LimitControls(owner=owner)
        if not self.control.acceptance_loaded:raise ValueError('Accepted genuine C3 controls required')
        self.lower=lower.CurrentC2RepairBand(control_field=self.control.lower)
        self.identity=self.control.identity;self.c=self.control.ranges.c
        if self.lower.identity!=self.identity:raise ValueError('Same C2/C3 band source required')
        self.source=self.lower.source;self.cache={};self.hashes={**self.lower.hashes,**self.control.hashes}
        for name in (control.NAME,control.RECEIPT,lower.NAME,lower.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.prefix=[dict(n) for n in self.control.ranges.phase.built['graph'].nodes]
        self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C3 repair band required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual C3 band source '+name)
            self.acceptance_loaded=True

    def project(self,row,op):
        f=op.flow
        if len(row)<4 or any(q.ctx is not self.c or q.scale.bases is not f.logs or q.ledger is not f.ledger for q in row):
            raise ValueError('Actual Taylor rows in same context/basis/ledger required')
        return [row[0],row[1],2*row[2],6*row[3]]

    def endpoint_source(self,ends):
        """Already-normalized actual endpoint; third factorial once."""
        ends=tuple(ends);key=('endpoint',ends);phase_field=self.control.ranges.phase
        if ends not in current.current.CELLS:raise ValueError('Admitted actual Z cell required')
        if key not in self.cache:
            packet=phase_field.leading_source(ends,'O3_power',(2,1));generic=packet['original_generic_source']
            op=phase_field.outer.owner.owner(ends)
            if packet['exact_common_P0_axial5'] is not op.P0:raise ValueError('Same actual independent P0 required')
            if not current.current.previous.equivalent_rows(generic['common_original_P0_axial5'],op.P0):
                raise ValueError('Same complete normalized P0 rows required')
            self.cache[key]=dict(source_family=self.identity,exact_Z_cell=ends,
                ordinary_leading_history_C3={key:self.project(row,op) for key,row in generic['common_own_five_histories_axial5'].items()},
                ordinary_independent_P0_C3=self.project(op.P0,op),
                ordinary_terminal_amplitude_C3=self.project(generic['common_velocity_E_axial5'],op),
                original_Rc_source_coordinate=2,reserved_band_left_coordinate=1,
                m_k_S_normalization_already_applied=True,third_Taylor_factorial_applied_once=True,
                range_witness_not_function_value=True)
        return self.cache[key]

    def leading_range(self,ends):
        """Independent whole reserved band: raw V/m/k S-normalized once."""
        ends=tuple(ends);key=('whole_band',ends)
        if ends not in current.current.CELLS:raise ValueError('Admitted actual Z cell required')
        if key not in self.cache:
            packet=self.source.query(ends,(1,1),(2,1));op=self.control.ranges.phase.outer.owner.owner(ends);f=op.flow
            if packet['exact_common_P0_axial5'] is not op.P0:raise ValueError('Same actual full-band P0 required')
            raw=packet['raw']['raw_current_radius_y_derivative_axial_coefficients'];invS=f.factor((0,-.5,0,0,0))
            histories={key:self.project([q*invS for q in rows[0]],op) if key in ('m','k') else self.project(rows[0],op)
                for key,rows in raw['histories'].items()}
            self.cache[key]=dict(actual_whole_reserved_band_leading_history_C3=histories,
                actual_whole_band_E_y_C3=[self.project(row,op) for row in raw['velocity']['theta']],
                actual_whole_band_V_y_C3=[self.project([q*invS for q in row],op) for row in raw['velocity']['axial']],
                actual_independent_P0_C3=self.project(op.P0,op),source_family=self.identity,
                exact_Z_cell=ends,actual_band_x=[1,2],same_original_P0=True,
                actual_power_offset='t in[2,2+log2]',original_reservation=self.source.reservation,
                source_binding=limit.BACKGROUND_BINDING,m_k_S_normalization_applied_once=True,
                third_Taylor_factorial_applied_once=True,range_witness_not_function_value=True)
        return self.cache[key]

    def quantitative_range(self,ends,parent_row):
        c=self.c;bd=ranges.Bounds(c)
        read=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(read(row) for row in rows))
        source=self.endpoint_source(ends);whole=self.leading_range(ends)
        A=ranges.JetBound(*(bd.row(q) for q in source['ordinary_terminal_amplitude_C3']))
        bounds=self.control.bounds;rho=read(bounds['original_C1_ball_radius_upper'])
        H=ranges.JetBound(rho,rho,read(bounds['actual_limit_and_iterate_ZZ_upper']),read(bounds['actual_limit_and_iterate_ZZZ_upper']))
        epsilon=read(bounds['actual_inverse_epsilon_upper'])
        weights=self.control.ranges.phase.report['actual_whole_Z_frequency_connection']['fresh_exact_integral_weights']
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
        incoming={key:bd.scaled(jet(rows),epsilon) for key,rows in parent_row['actual_normalized_N_scaled_histories_C3'].items()}
        histories={};masses={};complete={}
        for key,rate in current.RATES.items():
            r=c.mpf(rate.numerator)/rate.denominator
            mass=c.ln(2) if not rate else -c.expm1(-r*c.ln(2))/r
            if ep(mass)[0]<=0:raise ValueError('Original own positive band mass required')
            histories[key]=bd.add(incoming[key],bd.scaled(density[key],bd.constant(mass)));masses[key]=mass
            leading=ranges.JetBound(*(bd.row(q) for q in whole['actual_whole_reserved_band_leading_history_C3'][key]))
            complete[key]=bd.add(leading,histories[key])
        row=radial[0];E,V=row['E'],row['V'];EE=bd.product(E,E)
        dy=dict(m=bd.add(complete['m'],V),h=bd.add(bd.scaled(complete['h'],bd.constant(c.mpf('1.5'))),E),
            k=bd.add(bd.scaled(complete['k'],bd.constant(c.mpf('1.5'))),bd.product(E,V)),
            e=bd.add(complete['e'],bd.product(V,V),bd.scaled(EE,bd.constant(c.mpf('.5')))),
            p=bd.scaled(EE,bd.constant(c.mpf('.5'))))
        return dict(source_family=self.identity,exact_Z_cell=ends,band_x=[1,2],
            actual_endpoint_source=control.encoded(source),actual_whole_band_leading_source=control.encoded(whole),
            actual_selected_N_profile_y0_y1_y2_C3_bounds=ranges.record(radial),
            actual_signed_density_magnitude_C3_bounds=ranges.record(density),
            actual_incoming_correction_C3_bounds=ranges.record(incoming),actual_own_rate_positive_masses=masses,
            actual_partial_correction_C3_bounds=ranges.record(histories),actual_complete_history_C3_bounds=ranges.record(complete),
            actual_complete_history_y_C3_bounds=ranges.record(dy),
            original_beta_derivative_global_caps=beta,original_bump_y0_y1_y2_caps=ranges.record(shape),
            source_derivative_units_and_independent_P0_retained=True,selected_repair_integer_unchanged=True,
            range_caps_not_function_values=True,absolute_exterior_heat_Rh_C2_and_physical_oracle_not_admitted=True)

    def velocity_functions(self):return self.functions['profiles_at_band_x']['radial_y_rows']
    def history_functions(self):return self.functions['complete_histories']


def source_bindings():
    return dict(profiles_C3=current.ast_binding(profiles_C3),signed_density_C3=current.ast_binding(signed_density_C3),
        leading_C3=current.ast_binding(leading_C3),build=current.ast_binding(build),
        actual_endpoint_projection=current.ast_binding(CurrentC3RepairBand.endpoint_source),
        actual_ordinary_third_projection=current.ast_binding(CurrentC3RepairBand.project),
        actual_reserved_source_query=current.ast_binding(CurrentC3RepairBand.leading_range),
        quantitative_range=current.ast_binding(CurrentC3RepairBand.quantitative_range),
        accepted_C2_band_receipt=lower.RECEIPT,accepted_C2_band_receipt_sha256=sha(lower.RECEIPT),
        original_reserved_power_source=limit.BACKGROUND_BINDING)


def run(control_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC3RepairBand(control_field=control_field,owner=owner,require_checked=False)
        parent=json.loads(gzip.decompress((HERE/ranges.NAME).read_bytes()))
        rows=[field.quantitative_range(tuple(row['exact_Z_cell']),row) for row in parent['actual_four_Z_C3_target_transports']]
        report=dict(candidate_actual_C3_repair_band_constructed=True,source_family=field.identity,
            original_C3_control_and_C2_band_graph_prefix_length=len(field.prefix),
            exact_graph_nodes=field.control.ranges.phase.built['graph'].nodes,
            actual_C3_repair_band_functions=control.encoded(field.functions),actual_four_Z_C3_band_ranges=ranges.record(rows),
            source_bindings=source_bindings(),actual_selected_repair_integer=field.control.bounds['actual_same_repair_integer'],
            radial_inertial_recovery_C2_installed=False,C2_Rh_join_and_absolute_exterior_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(GATES+OPEN,False),
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(control.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Same-limit C3 repaired-band velocity, signed histories, relative endpoints and directed ranges constructed',flush=True)
    return field


if __name__=='__main__':run()
