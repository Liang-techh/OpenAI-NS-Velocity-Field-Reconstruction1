"""Actual repaired-band radial/inertial C2 functions and directed ranges.

One Z derivative is consumed from the accepted C3 histories and independent
P0. Positive L/E/C margins are proved on the same selected repair limit.
Physical radius is the source Rc*x, without duplicating the Rm factor.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C3_repair_band as band
import lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs as original
import lei_ren_part1_paper_compliant_current_original_reference_restore_functions as reference
import lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_finite_N as long_source

phase,current,target=band.phase,band.current,band.target
ranges=band.ranges.lower
HERE,PREFIX,sha,ep=band.HERE,band.PREFIX,band.sha,band.ep
NAME=PREFIX+'current_original_C2_repaired_recovery.json.gz'
RECEIPT=PREFIX+'current_original_C2_repaired_recovery_check.json'
GATES=('current_original_repaired_band_C2_radial_and_full_inertial_functions_installed',
    'current_original_repaired_band_positive_L_E_C_margins_installed',
    'current_original_repaired_band_C2_recovery_magnitude_ranges_installed')
OPEN=band.OPEN


def lower_rows(q):return phase.C2Function(q.value,q.Z,q.ZZ)
def axial_derivative(q):return phase.C2Function(q.Z,q.ZZ,q.ZZZ)


class SignedAlgebra(phase.C2Algebra):
    def divide(self,n,d,name):return self.div(n,d)


class MagnitudeAlgebra:
    """Absolute majorants for the same signed formulas, not field values."""
    def __init__(self,c,lowers):self.c=c;self.bd=ranges.Bounds(c);self.lowers=lowers
    def fixed(self,k):
        if isinstance(k,str):k=Fraction(k)
        if isinstance(k,Fraction):k=self.c.mpf(k.numerator)/k.denominator
        return self.bd.fixed(self.bd.constant(k))
    def add(self,*q):return self.bd.add(*q)
    def neg(self,q):return q
    def scale(self,q,k):return self.bd.product(q,self.fixed(k))
    def mul(self,a,b):return self.bd.product(a,b)
    def divide(self,n,d,name):return self.bd.quotient(n,d,self.lowers[name])


def recover(alg,inputs,parameters):
    """Unchanged n=0 signed moment recovery, with ordinary Z2 algebra."""
    E,V,Ey,Vy=(inputs[key] for key in ('E','V','E_y','V_y'))
    m,h,k,e,p=(inputs['histories'][key] for key in ('m','h','k','e','p'))
    dm,dh,dk,de,dp=(inputs['history_Z'][key] for key in ('m','h','k','e','p'))
    z,d,L,delta,R,S=(parameters[key] for key in ('Z','d','L','delta','R','Pstar'))
    one=alg.fixed(1);minus_delta=alg.add(one,alg.neg(delta))
    pressure=alg.add(inputs['P0'],p);pressure_Z=alg.add(inputs['P0_Z'],dp)
    transport=alg.add(alg.mul(alg.mul(m,z),minus_delta),alg.mul(dm,d))
    divL=lambda q:alg.divide(q,L,'L')
    Q=divL(alg.add(alg.scale(alg.mul(V,z),2),alg.neg(transport)))
    itl=divL(alg.add(alg.neg(E),alg.mul(h,alg.add(one,alg.neg(alg.scale(delta,'1/2')))),
        alg.neg(alg.scale(alg.mul(alg.mul(dh,z),minus_delta),'1/2'))))
    itq=divL(alg.add(alg.mul(alg.mul(k,z),alg.add(alg.scale(delta,2),alg.neg(one))),
        alg.neg(alg.mul(dk,d)),alg.mul(E,transport)))
    izl=divL(alg.add(alg.neg(V),alg.scale(alg.mul(alg.add(m,alg.neg(alg.mul(dm,z))),minus_delta),'1/2')))
    izq=divL(alg.add(alg.mul(V,transport),alg.scale(alg.mul(alg.mul(e,z),delta),2),
        alg.neg(alg.mul(de,d)),alg.scale(alg.mul(alg.mul(pressure,z),alg.add(one,delta)),2),
        alg.neg(alg.mul(pressure_Z,d))))
    my=inputs['history_y']['m'];my_Z=inputs['history_y_Z']['m']
    transport_y=alg.add(alg.mul(alg.mul(my,z),minus_delta),alg.mul(my_Z,d))
    Qy=divL(alg.add(alg.scale(alg.mul(Vy,z),2),alg.neg(transport_y)))
    C=alg.add(E,alg.neg(alg.scale(Ey,2)));B=alg.scale(Vy,2)
    It=alg.mul(R,alg.add(itl,alg.mul(S,itq)));Iz=alg.mul(R,alg.add(izl,alg.mul(S,izq)))
    den=alg.mul(C,E);kap=alg.add(alg.mul(C,C),alg.mul(B,B))
    excess=alg.add(kap,alg.neg(alg.scale(den,2)))
    H=alg.add(alg.mul(C,It),alg.neg(alg.mul(B,Iz)))
    D=alg.add(H,alg.neg(kap));J=alg.add(alg.mul(C,Iz),alg.mul(B,It))
    quotients={key:alg.divide(n,q,lower) for key,n,q,lower in (
        ('a',C,E,'E'),('b',B,E,'E'),('t0',alg.neg(B),C,'C'),
        ('p1',It,E,'E'),('p2',Iz,E,'E'),('kappa',kap,den,'den'),('D',D,den,'den'),('J',J,den,'den'))}
    quotients['Delta']=alg.add(quotients['kappa'],alg.fixed(-2))
    den3=alg.mul(alg.mul(den,den),den)
    discriminant=alg.divide(alg.add(alg.scale(alg.mul(alg.mul(D,D),den),2),
        alg.neg(alg.mul(excess,alg.mul(J,J)))),den3,'den3')
    prefactor=parameters['Pstar_sqrt_R_over_2'];S2=alg.mul(S,S)
    cylindrical=dict(Utheta=alg.mul(S,E),Uz=alg.mul(S,V),Ur=alg.mul(prefactor,Q),
        Ur_y=alg.mul(prefactor,alg.add(Qy,alg.scale(Q,'1/2'))),Pi=alg.mul(S2,pressure))
    return dict(common_E=E,common_V=V,common_radial_Q=Q,common_radial_Q_y=Qy,
        common_absolute_pressure=pressure,common_absolute_pressure_Z=pressure_Z,
        full_signed_inertial_sectors=dict(theta_linear=itl,theta_quadratic=itq,axial_linear=izl,axial_quadratic=izq),
        actual_generic_numerators=dict(E=E,C=C,B=B,inertial_theta=It,inertial_axial=Iz,
            positive_denominator=den,kappa=kap,kappa_minus2=excess,
            H0_minus2=alg.add(H,alg.neg(alg.scale(den,2))),D=D,J=J),
        full_signed_generic_quotients=quotients,full_signed_generic_discriminant=discriminant,
        original_cylindrical_velocity_pressure=cylindrical)


def denominator_proof(field,ends,parent):
    c=field.c;bd=ranges.Bounds(c);op=field.band.control.ranges.phase.outer.owner.owner(ends)
    endpoint=field.band.endpoint_source(ends);A=endpoint['ordinary_terminal_amplitude_C3'][0]
    logA=ranges.positive_lower(A);delta=op.reference.delta
    if ep(delta)[0]<=0 or ep(delta)[1]>=ep(c.mpf('.5'))[0]:raise ValueError('Same original 0<delta<1/2 required')
    profiles=parent['actual_selected_N_profile_y0_y1_y2_C3_bounds']
    caps=[current.packets.interval(c,profiles[j]['F'][0]['log_absolute_upper']) for j in (0,1)]
    threshold=c.ln(c.mpf(1)/8)
    if any(ep(q)[1]>=ep(threshold)[0] for q in caps):raise ValueError('Same selected repair F/Fy<1/8 required')
    mu=c.exp(field.band.control.ranges.phase.outer.logmu);alpha=c.mpf('.5')+mu
    if ep(alpha)[0]<ep(c.mpf('.5'))[0] or ep(alpha)[1]>=ep(c.mpf(2)/3)[0]:raise ValueError('Original 1/2<=alpha<2/3 required')
    # Factor the same actual positive amplitude: E=A*(power+F),
    # C=A*((1+2alpha)*power+F-2Fy); power>1/2 on x in[1,2].
    logE=c.mpf(ep(logA-c.ln(4))[0]);logC=c.mpf(ep(logA-c.ln(2))[0])
    Llower=1-delta;logL=c.mpf(ep(c.ln(Llower))[0]);logden=c.mpf(ep(logE+logC)[0])
    return dict(source_family=field.identity,exact_Z_cell=ends,band_x=[1,2],
        actual_positive_terminal_amplitude_log_lower=logA,actual_delta=delta,
        actual_meridional_L_lower=Llower,actual_L_log_lower=logL,
        actual_E_log_lower=logE,actual_C_log_lower=logC,actual_C_times_E_log_lower=logden,
        actual_selected_N_F_Fy_log_caps=caps,actual_perturbation_threshold_log=threshold,
        actual_alpha=alpha,strict_relative_margins=dict(E='E/A>3/8>1/4',C='C/A>5/8>1/2'),
        proof='Same A>0 factors both functions; power>=2^-alpha>1/2, |F|,|Fy|<1/8, 1+2alpha>=2.',
        same_selected_repair_integer=field.band.control.bounds['actual_same_repair_integer'],
        source_range_is_not_function_value=True,stress_cone_not_admitted=True)


def build(field):
    bf=field.band.functions;g=field.band.control.ranges.phase.built['graph'];alg=SignedAlgebra(g)
    built=field.band.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    profiles=field.band.velocity_functions();hist=field.band.history_functions();hy=bf['complete_history_y']
    P0=bf['original_leading_power']['independent_P0'];z=g.symbol('Z');x=bf['original_band_variable']
    Z=phase.C2Function(z,g.one,g.zero);logP=built['parameters']['logP']
    delta=alg.fixed(g.unary('exp',g.sub(g.neg(g.mul(g.constant(4),logP)),g.constant(30))))
    d=alg.add(alg.fixed(1),alg.neg(alg.mul(Z,Z)));L=alg.add(alg.fixed(1),alg.neg(alg.mul(delta,alg.mul(Z,Z))))
    recipe=field.band.lower.functions['original_C1_background']['leading_history_P0_source_binding']
    Rc=g.node('current_original_leading_function_recipe',source_family=field.identity,
        native_chart='O3_power',recipe=recipe,quantity='actual_source_radius',Z_order=0,
        coordinate=g.constant(2).node,Z_variable=z.node,Z_independent=True,
        quantity_path=['original_closed_O3_background','actual_source_radius'],
        source_projection_binding=current.ast_binding(CurrentC2RepairedRecovery.parameter_source),
        exact_definition='Rc=op.Rm_factor*exp(exp40+20); actual source coordinate2',
        defining_quantity_not_a_range_value=True,source_background_binding=band.limit.BACKGROUND_BINDING)
    R=g.mul(Rc,x);S=g.unary('exp',logP)
    prefactor=g.mul(S,g.unary('positive_sqrt',g.mul(g.constant('1/2'),R)))
    inputs=dict(E=lower_rows(profiles[0]['E']),V=lower_rows(profiles[0]['V']),
        E_y=lower_rows(profiles[1]['E']),V_y=lower_rows(profiles[1]['V']),
        histories={key:lower_rows(q) for key,q in hist.items()},history_Z={key:axial_derivative(q) for key,q in hist.items()},
        history_y={key:lower_rows(q) for key,q in hy.items()},history_y_Z={key:axial_derivative(q) for key,q in hy.items()},
        P0=lower_rows(P0),P0_Z=axial_derivative(P0))
    parameters=dict(Z=Z,d=d,L=L,delta=delta,R=alg.fixed(R),Pstar=alg.fixed(S),Pstar_sqrt_R_over_2=alg.fixed(prefactor))
    functions=recover(alg,inputs,parameters)
    proofs={str(ends):g.node('proved_actual_repaired_band_positive_denominators',
        actual_positive_proof=band.control.encoded(proof),expressions=dict(L=L.value.node,
            E=functions['common_E'].value.node,C=functions['actual_generic_numerators']['C'].value.node),
        actual_checked_C3_band_source=band.RECEIPT,actual_checked_C3_band_sha256=sha(band.RECEIPT),
        theorem_binding=current.ast_binding(denominator_proof),only_same_actual_selected_limit=True)
        for ends,proof in field.proofs.items()}
    return dict(actual_C3_repaired_inputs=inputs,actual_source_parameters=parameters,
        recovered_C2_functions=functions,positive_denominator_certificates=proofs,
        actual_Rc=Rc,actual_R=R,original_absolute_radius_offset=bf['original_leading_power']['original_radius_offset'],
        source_relative_Rm_radius_offset=g.add(logP,g.constant(9),g.unary('log',x)),
        original_Rm_factor_not_counted_twice=True,original_band_variable=x,
        axial_orders=[0,1,2],radial_Q_y_orders=[0,1],
        ordinary_derivative_shift_consumes_actual_C3_row=True,independent_P0_and_full_signed_pressure_retained=True,
        same_original_S_normalization_no_second_division=True,
        original_cylindrical_prefactors_before_physical_time_map=True,
        full_radial_y4_stress_C3_Rh_heat_global_and_temporal_recursion_not_admitted=True)


class CurrentC2RepairedRecovery:
    def __init__(self,band_field=None,owner=None,require_checked=True):
        self.band=band_field if band_field is not None else band.CurrentC3RepairBand(owner=owner)
        if not self.band.acceptance_loaded:raise ValueError('Accepted actual C3 repaired band required')
        self.c=self.band.c;self.identity=self.band.identity;self.hashes=dict(self.band.hashes);self.cache={}
        for name in (band.NAME,band.RECEIPT,Path(original.__file__).name,
                PREFIX+'current_generic_shear_moment_recovery.py',PREFIX+'current_generic_shear_inputs.py',
                PREFIX+'pressure_source.json',Path(__file__).name):self.hashes[name]=sha(name)
        pressure=json.loads((HERE/(PREFIX+'pressure_source.json')).read_bytes())['compliant_source']
        if pressure['implicit_source_sha256']!=self.identity['implicit_source_sha256']:raise ValueError('Same actual similarity source required')
        definition=pressure['implicit_source_definition']
        if (definition['Md']!='40' or definition['logPstar']!='exp(Md)+11'
                or definition['delta']!='min(1e-200,exp(-4logPstar-30))'):
            raise ValueError('Original Md40 logP/delta definition required')
        self.parameter_definition=definition
        report=json.loads(gzip.decompress((HERE/band.NAME).read_bytes()))
        self.rows={tuple(row['exact_Z_cell']):row for row in report['actual_four_Z_C3_band_ranges']}
        if set(self.rows)!=set(current.current.CELLS):raise ValueError('All four same actual Z cells required')
        self.parameter_identities={ends:self.similarity_parameter_identity(ends) for ends in self.rows}
        self.proofs={ends:denominator_proof(self,ends,row) for ends,row in self.rows.items()}
        self.prefix=[dict(n) for n in self.band.control.ranges.phase.built['graph'].nodes]
        self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C2 repaired recovery required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed repaired recovery source '+name)
            self.acceptance_loaded=True

    def similarity_parameter_identity(self,ends):
        """Bind the live delta to the exact source expression and min branch."""
        op=self.band.control.ranges.phase.outer.owner.owner(ends);c=self.c;ref=op.reference
        if not isinstance(ref,reference.ActualReferenceRestoreFunctions) or ref.c is not c:
            raise ValueError('Same actual reference parameter provider required')
        if ref.delta.ctx is not c or ref.logP.ctx is not c:
            raise ValueError('Same live delta/logP context required')
        if ref.logP._mpi_!=(op.flow.logs[1]/2)._mpi_:
            raise ValueError('Same canonical Pstar log basis required')
        branch=c.exp(-4*ref.logP-30)
        if ep(branch)[1]>=ep(c.mpf('1e-200'))[0] or ref.delta._mpi_!=branch._mpi_:
            raise ValueError('Live original delta must be the exact exponential min branch')
        return dict(source_family=self.identity,exact_Z_cell=ends,
            original_implicit_parameter_definition=self.parameter_definition,
            canonical_source_logP=ref.logP,canonical_source_delta=ref.delta,
            actual_delta_exponential_branch=branch,exact_live_delta_tuple_equals_branch=True,
            same_source_log_basis_and_delta_context=True,
            formal_graph_logP='exp(40)+11',formal_graph_delta='exp(-4*logP-30)',
            identity='Checked Md40 source definition, inherited original parameter provider, and exact live min-branch tuple agree.',
            actual_reference_provider=current.ast_binding(reference.ActualReferenceRestoreFunctions.__init__),
            actual_source_to_reference_provider=current.ast_binding(long_source.WholeZLongRmFiniteN.owner),
            source_interval_is_not_selected_function_parameter=True)

    def parameter_source(self,ends):
        ends=tuple(ends)
        if ends not in self.rows:raise ValueError('Admitted actual Z cell required')
        if ends not in self.cache:
            op=self.band.control.ranges.phase.outer.owner.owner(ends);f=op.flow;c=self.c
            inlet=self.band.source.query(ends,(1,1));whole=self.band.source.query(ends,(1,1),(2,1))
            if inlet['exact_common_P0_axial5'] is not op.P0 or whole['exact_common_P0_axial5'] is not op.P0:
                raise ValueError('Same actual independent source P0 required')
            z=f.jet(op.reference.z)
            if not current.current.previous.equivalent_rows(z,op.zrows):raise ValueError('Same actual axial variable required')
            Rc=inlet['actual_source_radius'];R=whole['actual_source_radius']
            for q in (Rc,R,op.Pstar,op.Rm_factor,*z):
                if q.ctx is not c or q.scale.bases is not f.logs or q.ledger is not f.ledger:
                    raise ValueError('Same source context/basis/ledger required')
            self.cache[ends]=dict(source_family=self.identity,exact_Z_cell=ends,
                actual_Rc=Rc,actual_whole_band_R=R,actual_Pstar=op.Pstar,actual_Rm_factor=op.Rm_factor,
                actual_delta=op.reference.delta,actual_axial_Taylor_rows=z,
                actual_similarity_parameter_identity=self.parameter_identities[ends],
                exact_radius_definition='R=op.Rm_factor*exp(exp40+20+log(x)); Rc=R(x=1)',
                exact_Rm_relative_log_offset=c.exp(40)+20,actual_source_geometry=whole['raw']['geometry'],
                original_source_P0_is_same_object=True,radius_and_scale_Z_independent=True)
        return self.cache[ends]

    def quantitative_range(self,ends):
        parent=self.rows[ends];proof=self.proofs[ends];c=self.c;bd=ranges.Bounds(c)
        read=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(read(row) for row in rows[:3]))
        shifted=lambda rows:ranges.JetBound(*(read(row) for row in rows[1:4]))
        profiles=parent['actual_selected_N_profile_y0_y1_y2_C3_bounds']
        histories=parent['actual_complete_history_C3_bounds'];hy=parent['actual_complete_history_y_C3_bounds']
        P0=self.band.endpoint_source(ends)['ordinary_independent_P0_C3'];P=[bd.row(q) for q in P0]
        inputs=dict(E=jet(profiles[0]['E']),V=jet(profiles[0]['V']),E_y=jet(profiles[1]['E']),V_y=jet(profiles[1]['V']),
            histories={key:jet(row) for key,row in histories.items()},history_Z={key:shifted(row) for key,row in histories.items()},
            history_y={key:jet(row) for key,row in hy.items()},history_y_Z={key:shifted(row) for key,row in hy.items()},
            P0=ranges.JetBound(*P[:3]),P0_Z=ranges.JetBound(*P[1:4]))
        source=self.parameter_source(ends);delta=source['actual_delta'];zero=bd.zero;one=bd.one
        R,S=bd.row(source['actual_whole_band_R']),bd.row(source['actual_Pstar'])
        prefactor=ranges.LogUpper(c,S.log+(R.log-c.ln(2))/2)
        parameters=dict(Z=ranges.JetBound(one,one,zero),d=ranges.JetBound(one,bd.constant(2),bd.constant(2)),
            L=ranges.JetBound(one,bd.constant(2*delta),bd.constant(2*delta)),delta=bd.fixed(bd.constant(delta)),
            R=bd.fixed(R),Pstar=bd.fixed(S),Pstar_sqrt_R_over_2=bd.fixed(prefactor))
        lowers=dict(L=proof['actual_L_log_lower'],E=proof['actual_E_log_lower'],C=proof['actual_C_log_lower'],
            den=proof['actual_C_times_E_log_lower'],den3=3*proof['actual_C_times_E_log_lower'])
        caps=recover(MagnitudeAlgebra(c,lowers),inputs,parameters)
        return dict(source_family=self.identity,exact_Z_cell=ends,band_x=[1,2],
            actual_positive_denominator_proof=proof,actual_parameter_source=band.control.encoded(source),
            actual_repaired_recovery_C2_magnitude_bounds=ranges.record(caps),
            source_cylindrical_prefactors_before_time_map=True,
            actual_generic_cone_signed_margins_not_admitted=True,current_numeric_point_field_oracle_not_installed=True)

    def radial_functions(self):
        q=self.functions['recovered_C2_functions']
        return q['common_radial_Q'],q['common_radial_Q_y']
    def cylindrical_functions(self):return self.functions['recovered_C2_functions']['original_cylindrical_velocity_pressure']


def source_bindings():
    return dict(actual_repaired_recovery=current.ast_binding(recover),ordinary_C3_derivative_shift=current.ast_binding(axial_derivative),
        positive_denominator_proof=current.ast_binding(denominator_proof),build=current.ast_binding(build),
        similarity_parameter_identity=current.ast_binding(CurrentC2RepairedRecovery.similarity_parameter_identity),
        actual_parameter_source=current.ast_binding(CurrentC2RepairedRecovery.parameter_source),
        quantitative_range=current.ast_binding(CurrentC2RepairedRecovery.quantitative_range),
        unchanged_original_recovery=original.source_bindings(),original_recover_inputs=current.ast_binding(original.recover_inputs),
        same_original_reserved_source=band.limit.BACKGROUND_BINDING)


def run(band_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC2RepairedRecovery(band_field=band_field,owner=owner,require_checked=False)
        rows=[field.quantitative_range(ends) for ends in current.current.CELLS]
        report=dict(candidate_actual_C2_repaired_recovery_constructed=True,source_family=field.identity,
            original_C3_band_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.band.control.ranges.phase.built['graph'].nodes,
            actual_C2_repaired_recovery_functions=band.control.encoded(field.functions),actual_four_Z_C2_repaired_recovery_ranges=ranges.record(rows),
            source_bindings=source_bindings(),**dict.fromkeys(GATES+OPEN,False),
            current_numeric_point_field_oracle_installed=False,actual_generic_stress_cone_installed=False,
            physical_time_Cartesian_dispatcher_global_Rh_heat_and_temporal_recursion_installed=False,
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(band.control.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual repaired radial/inertial C2 functions, positive denominators and whole-band bounds constructed',flush=True)
    return field


if __name__=='__main__':run()
