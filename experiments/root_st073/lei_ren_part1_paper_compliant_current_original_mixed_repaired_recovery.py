"""Same repaired limit: ordinary y4/Z3 inputs and y4/Z2 recovery.

Higher compact-bump rows are signed analytic functions, not range samples.
Complete moment rows follow the original own-rate FTC. Physical radial,
inertial and shear prefactors are differentiated once in source coordinates.
"""
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C2_repaired_recovery as lower
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as generic

band,phase,current,target=lower.band,lower.phase,lower.current,lower.target
ranges=band.ranges
HERE,PREFIX,sha,ep=lower.HERE,lower.PREFIX,lower.sha,lower.ep
NAME=PREFIX+'current_original_mixed_repaired_recovery.json.gz'
RECEIPT=PREFIX+'current_original_mixed_repaired_recovery_check.json'
GATES=('current_original_repaired_band_y4_Z3_profile_history_functions_installed',
    'current_original_repaired_band_y4_Z2_velocity_and_y3_Z2_stress_recovery_installed',
    'current_original_repaired_band_mixed_directed_magnitude_ranges_installed')
OPEN=lower.OPEN


def product_row(alg,a,b,j):
    return alg.add(*(alg.scale(alg.mul(a[k],b[j-k]),math.comb(j,k)) for k in range(j+1)))


def shifted_rows(alg,rows,rate):
    from fractions import Fraction
    rate=Fraction(str(rate))
    return [alg.add(*(alg.scale(rows[k],math.comb(j,k)*rate**(j-k)) for k in range(j+1)))
        for j in range(len(rows))]


def profiles_y4(field):
    prior=field.lower.band.functions['profiles_at_band_x']
    built=field.lower.band.control.functions['exact_C2_limit_functions']['exact_C1_limit_adapter']
    g=built['graph'];alg=target.C3Algebra(g);W=built['repair_weights'];x=field.lower.band.functions['original_band_variable']
    logx=g.unary('log',x);ell=W['ell'];normal=W['normal']
    bumps=[prior['bumps'],prior['bump_y'],prior['bump_yy']];raw=[]
    for center in W['centers']:
        arg=g.quotient(g.sub(logx,center),ell,'original positive log bump width')
        raw.append([g.node('original_compact_beta_ordinary_derivative',argument=arg.node,
            ordinary_beta_derivative_order=n,polynomial_coefficients=phase.flat_source.BETA_POLYNOMIALS[n],
            definition='exp(-1/(1-t^2))*P_n(t)/(1-t^2)^(2*n) for abs(t)<1; zero otherwise',
            original_derivative_source=current.ast_binding(phase.flat_source.beta_jets),
            original_polynomial_source=current.ast_binding(phase.flat_source.beta_polynomials),
            Taylor_coefficient_factorial=math.factorial(n),ordinary_function_not_Taylor_coefficient=True,
            lazy_outside_support=True,all_support_endpoint_jets_zero=True,derivative_of_range_endpoint=False)
            for n in range(5)])
    for j in (3,4):
        row=[]
        for derivatives in raw:
            terms=[g.quotient(g.mul(g.constant(math.comb(j,k)*(-1)**(j-k)),derivatives[k]),
                g.mul(*([ell]*(k+1)),normal,x),'original positive ell^(k+1) J0 x') for k in range(j+1)]
            row.append(g.add(*terms))
        bumps.append(row)
    rows=list(prior['radial_y_rows']);h=field.lower.band.control.control_functions()
    A=field.lower.band.control.ranges.target.functions['actual_terminal_amplitude']
    for j in (3,4):
        def correction(indices):
            return target.C3Function(*[g.quotient(g.add(*(g.mul(bumps[j][b],getattr(h[i],key)) for b,i in indices)),
                built['N'],'same selected exact positive integer N') for key in ('value','Z','ZZ','ZZZ')])
        F,G=correction(((0,2),(1,3),(2,4))),correction(((0,0),(2,1)))
        fixed=g.mul(prior['original_power'],g.constant((-1)**j),*[prior['original_alpha']]*j)
        original=alg.scale(A,fixed);dE,dV=alg.mul(A,F),alg.mul(A,G)
        rows.append(dict(F=F,G=G,original_E=original,delta_E=dE,delta_V=dV,E=alg.add(original,dE),V=dV))
    return dict(radial_y_rows=rows,bump_ordinary_y_rows=bumps,raw_beta_ordinary_rows=raw,
        original_power=prior['original_power'],original_alpha=prior['original_alpha'],
        accepted_y0_y1_y2_C3_handles_retained=True,ordinary_radial_coordinate='y=log(x)',
        no_extra_Taylor_factorial=True,endpoint_flatness_from_original_compact_beta=True)


def history_y4(alg,E,V,histories,first_y=None):
    rows=[histories]
    for j in range(4):
        EE=product_row(alg,E,E,j)
        row=dict(m=alg.add(V[j],alg.neg(rows[j]['m'])),
            h=alg.add(E[j],alg.neg(alg.scale(rows[j]['h'],'3/2'))),
            k=alg.add(product_row(alg,E,V,j),alg.neg(alg.scale(rows[j]['k'],'3/2'))),
            e=alg.add(product_row(alg,V,V,j),alg.neg(alg.scale(EE,'1/2')),alg.neg(rows[j]['e'])),
            p=alg.scale(EE,'1/2'))
        rows.append(first_y if j==0 and first_y is not None else row)
    return rows


def mixed_recover(alg,profiles,histories,history_Z,P0,P0_Z,parameters):
    """Original signed recovery with ordinary y-binomial products and Z2 jets."""
    E,V=([row[key] for row in profiles] for key in ('E','V'))
    z,d,L,delta,R,S=(parameters[key] for key in ('Z','d','L','delta','R','Pstar'))
    one=alg.fixed(1);minus_delta=alg.add(one,alg.neg(delta));divL=lambda q:alg.divide(q,L,'L')
    P=[alg.add(P0,histories[0]['p'])]+[row['p'] for row in histories[1:]]
    PZ=[alg.add(P0_Z,history_Z[0]['p'])]+[row['p'] for row in history_Z[1:]]
    transport=[alg.add(alg.mul(alg.mul(row['m'],z),minus_delta),alg.mul(derivative['m'],d))
        for row,derivative in zip(histories,history_Z)]
    Q=[divL(alg.add(alg.scale(alg.mul(v,z),2),alg.neg(t))) for v,t in zip(V,transport)]
    sectors={key:[] for key in ('theta_linear','theta_quadratic','axial_linear','axial_quadratic')}
    for j in range(4):
        m,h,k,e=(histories[j][key] for key in ('m','h','k','e'))
        dm,dh,dk,de=(history_Z[j][key] for key in ('m','h','k','e'))
        sectors['theta_linear'].append(divL(alg.add(alg.neg(E[j]),
            alg.mul(h,alg.add(one,alg.neg(alg.scale(delta,'1/2')))),
            alg.neg(alg.scale(alg.mul(alg.mul(dh,z),minus_delta),'1/2')))))
        sectors['theta_quadratic'].append(divL(alg.add(alg.mul(alg.mul(k,z),alg.add(alg.scale(delta,2),alg.neg(one))),
            alg.neg(alg.mul(dk,d)),product_row(alg,E,transport,j))))
        sectors['axial_linear'].append(divL(alg.add(alg.neg(V[j]),
            alg.scale(alg.mul(alg.add(m,alg.neg(alg.mul(dm,z))),minus_delta),'1/2'))))
        sectors['axial_quadratic'].append(divL(alg.add(product_row(alg,V,transport,j),
            alg.scale(alg.mul(alg.mul(e,z),delta),2),alg.neg(alg.mul(de,d)),
            alg.scale(alg.mul(alg.mul(P[j],z),alg.add(one,delta)),2),alg.neg(alg.mul(PZ[j],d)))))
    C=[alg.add(E[j],alg.neg(alg.scale(E[j+1],2))) for j in range(4)]
    B=[alg.scale(V[j+1],2) for j in range(4)]
    inertial={}
    for name,linear,quadratic in (('theta','theta_linear','theta_quadratic'),('axial','axial_linear','axial_quadratic')):
        combined=[alg.add(a,alg.mul(S,b)) for a,b in zip(sectors[linear],sectors[quadratic])]
        inertial[name]=[alg.mul(R,q) for q in shifted_rows(alg,combined,1)]
    prefactor=parameters['Pstar_sqrt_R_over_2'];shear_prefactor=parameters['Pstar_over_sqrt_2R']
    physical=dict(Utheta=[alg.mul(S,q) for q in E],Uz=[alg.mul(S,q) for q in V],
        Ur=[alg.mul(prefactor,q) for q in shifted_rows(alg,Q,'1/2')],Pi=[alg.mul(alg.mul(S,S),q) for q in P])
    stress={key:[alg.mul(alg.mul(prefactor,S) if 'quadratic' in key else prefactor,q)
        for q in shifted_rows(alg,rows,'1/2')] for key,rows in sectors.items()}
    stress['shear_theta']=[alg.mul(shear_prefactor,q) for q in shifted_rows(alg,[alg.neg(q) for q in C],'-1/2')]
    stress['shear_axial']=[alg.mul(shear_prefactor,q) for q in shifted_rows(alg,B,'-1/2')]
    return dict(common_radial_Q_y_rows=Q,common_absolute_pressure_y_rows=P,common_absolute_pressure_y_Z_rows=PZ,
        full_signed_inertial_sector_y_rows=sectors,actual_generic_C_B_y_rows=dict(C=C,B=B),
        actual_generic_inertial_y_rows=inertial,original_cylindrical_velocity_pressure_y_rows=physical,
        original_physical_signed_stress_y_rows=stress)


def build(field):
    previous=field.lower.functions;g=field.graph;alg=target.C3Algebra(g)
    profiles=profiles_y4(field);radial=profiles['radial_y_rows']
    histories=history_y4(alg,[row['E'] for row in radial],[row['V'] for row in radial],
        field.lower.band.history_functions(),field.lower.band.functions['complete_history_y'])
    P0=field.lower.band.functions['original_leading_power']['independent_P0'];zero=alg.fixed(0)
    if (len(radial)!=5 or len(histories)!=5 or not isinstance(P0,target.C3Function)
            or any(not isinstance(q,target.C3Function) for row in radial+histories for q in row.values())):
        raise ValueError('Five ordinary y rows of genuine C3 profiles and histories plus C3 P0 required')
    inputs=[{key:lower.lower_rows(q) for key,q in row.items()} for row in radial]
    M=[{key:lower.lower_rows(q) for key,q in row.items()} for row in histories]
    MZ=[{key:lower.axial_derivative(q) for key,q in row.items()} for row in histories]
    params=dict(previous['actual_source_parameters']);a=lower.SignedAlgebra(g)
    params['Pstar_over_sqrt_2R']=a.fixed(g.quotient(params['Pstar_sqrt_R_over_2'].value,params['R'].value,
        'original positive source radius R=Rc*x'))
    recovered=mixed_recover(a,inputs,M,MZ,lower.lower_rows(P0),lower.axial_derivative(P0),params)
    # Exact historical aliases are kept despite algebraically equal alternate parenthesization.
    old=previous['recovered_C2_functions']
    recovered['common_radial_Q_y_rows'][:2]=[old['common_radial_Q'],old['common_radial_Q_y']]
    recovered['common_absolute_pressure_y_rows'][0]=old['common_absolute_pressure']
    recovered['common_absolute_pressure_y_Z_rows'][0]=old['common_absolute_pressure_Z']
    for key in recovered['full_signed_inertial_sector_y_rows']:
        recovered['full_signed_inertial_sector_y_rows'][key][0]=old['full_signed_inertial_sectors'][key]
    for key in ('C','B'):recovered['actual_generic_C_B_y_rows'][key][0]=old['actual_generic_numerators'][key]
    for key in ('theta','axial'):
        recovered['actual_generic_inertial_y_rows'][key][0]=old['actual_generic_numerators']['inertial_'+key]
    for key in recovered['original_cylindrical_velocity_pressure_y_rows']:
        recovered['original_cylindrical_velocity_pressure_y_rows'][key][0]=old['original_cylindrical_velocity_pressure'][key]
    recovered['original_cylindrical_velocity_pressure_y_rows']['Ur'][1]=old['original_cylindrical_velocity_pressure']['Ur_y']
    return dict(actual_C3_profiles_y4=profiles,actual_C3_complete_history_y_rows=histories,
        actual_C3_independent_P0_y_rows=[P0,zero,zero,zero,zero],actual_C2_profile_inputs=inputs,
        actual_C2_history_inputs=M,actual_C2_history_Z_inputs=MZ,actual_source_parameters=params,
        recovered_mixed_C2_functions=recovered,profile_and_history_y_orders=list(range(5)),
        input_axial_orders=[0,1,2,3],recovered_axial_orders=[0,1,2],radial_velocity_y_orders=list(range(5)),
        inertial_and_shear_y_orders=list(range(4)),same_original_selected_repair_integer=True,
        accepted_C3_lower_profiles_histories_and_C2_recovery_handles_retained=True,
        physical_radial_inertial_half_power_and_shear_negative_half_power_applied_once=True,
        generic_inertial_R_prefactor_rate_one_applied_once=True,
        original_source_coordinates_before_time_and_Cartesian_mapping=True,
        absolute_Rh_heat_global_cone_numeric_oracle_and_temporal_recursion_not_admitted=True)


class CurrentMixedRepairedRecovery:
    def __init__(self,recovery_field=None,owner=None,require_checked=True):
        self.lower=recovery_field if recovery_field is not None else lower.CurrentC2RepairedRecovery(owner=owner)
        if not self.lower.acceptance_loaded:raise ValueError('Accepted actual C2 repaired recovery required')
        self.c,self.identity=self.lower.c,self.lower.identity
        self.graph=self.lower.band.control.ranges.phase.built['graph'];self.hashes=dict(self.lower.hashes)
        for name in (lower.NAME,lower.RECEIPT,Path(phase.flat_source.__file__).name,Path(generic.__file__).name,Path(__file__).name):
            self.hashes[name]=sha(name)
        self.prefix=[dict(n) for n in self.graph.nodes];self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual mixed repair recovery required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed mixed recovery source '+name)
            self.acceptance_loaded=True

    def quantitative_range(self,ends):
        c=self.c;bd=ranges.Bounds(c);prior=self.lower.rows[ends]
        read=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(read(row) for row in rows))
        radial=[{key:jet(row) for key,row in p.items()} for p in prior['actual_selected_N_profile_y0_y1_y2_C3_bounds']]
        A=ranges.JetBound(*(bd.row(q) for q in self.lower.band.endpoint_source(ends)['ordinary_terminal_amplitude_C3']))
        control=self.lower.band.control;bounds=control.bounds
        rho=read(bounds['original_C1_ball_radius_upper'])
        H=ranges.JetBound(rho,rho,read(bounds['actual_limit_and_iterate_ZZ_upper']),read(bounds['actual_limit_and_iterate_ZZZ_upper']))
        epsilon=read(bounds['actual_inverse_epsilon_upper'])
        W=control.ranges.phase.report['actual_whole_Z_frequency_connection']['fresh_exact_integral_weights']
        ell=current.packets.interval(c,W['radius']);normal=current.packets.interval(c,W['raw_normalization'])
        if ep(ell)[0]<=0 or ep(normal)[0]<=0:raise ValueError('Original positive ell/J0 required')
        beta=[c.exp(-1)]+[phase.flat_source.BETA_CONSTANTS[n]*(2*n)**(2*n)*c.exp(-2*n) for n in range(1,5)]
        shape=[bd.sum(*(bd.constant(math.comb(j,k)*beta[k]/(ell**(k+1)*normal)) for k in range(j+1))) for j in (3,4)]
        alpha=c.mpf('.5')+c.exp(control.ranges.phase.outer.logmu)
        for j,cap in zip((3,4),shape):
            F,G=(bd.scaled(H,bd.scale(cap*epsilon,n)) for n in (3,2))
            original=bd.scaled(A,bd.constant(alpha**j));dE,dV=bd.product(A,F),bd.product(A,G)
            radial.append(dict(F=F,G=G,original_E=original,delta_E=dE,delta_V=dV,E=bd.add(original,dE),V=dV))
        class C3Magnitude:
            fixed=lambda self,k:bd.fixed(bd.constant(k))
            add=lambda self,*q:bd.add(*q)
            neg=lambda self,q:q
            mul=lambda self,a,b:bd.product(a,b)
            def scale(self,q,k):
                from fractions import Fraction
                k=Fraction(str(k));return bd.scaled(q,bd.constant(c.mpf(k.numerator)/k.denominator))
        M=history_y4(C3Magnitude(),[row['E'] for row in radial],[row['V'] for row in radial],
            {key:jet(q) for key,q in prior['actual_complete_history_C3_bounds'].items()},
            {key:jet(q) for key,q in prior['actual_complete_history_y_C3_bounds'].items()})
        br=lower.ranges.Bounds(c)
        trim=lambda q:lower.ranges.JetBound(q.value,q.Z,q.ZZ)
        shift=lambda q:lower.ranges.JetBound(q.Z,q.ZZ,q.ZZZ)
        P0=ranges.JetBound(*(bd.row(q) for q in self.lower.band.endpoint_source(ends)['ordinary_independent_P0_C3']))
        source=self.lower.parameter_source(ends);de=source['actual_delta'];zero,one=br.zero,br.one
        R,S=br.row(source['actual_whole_band_R']),br.row(source['actual_Pstar'])
        pre=lower.ranges.LogUpper(c,S.log+(R.log-c.ln(2))/2)
        # The inverse sqrt uses the actual positive radius lower, never the upper radius cap.
        logRmin=lower.ranges.positive_lower(source['actual_Rc'])
        shear=lower.ranges.LogUpper(c,S.log-(logRmin+c.ln(2))/2)
        parameters=dict(Z=lower.ranges.JetBound(one,one,zero),d=lower.ranges.JetBound(one,br.constant(2),br.constant(2)),
            L=lower.ranges.JetBound(one,br.constant(2*de),br.constant(2*de)),delta=br.fixed(br.constant(de)),
            R=br.fixed(R),Pstar=br.fixed(S),Pstar_sqrt_R_over_2=br.fixed(pre),Pstar_over_sqrt_2R=br.fixed(shear))
        caps=mixed_recover(lower.MagnitudeAlgebra(c,{'L':self.lower.proofs[ends]['actual_L_log_lower']}),
            [{key:trim(q) for key,q in row.items()} for row in radial],
            [{key:trim(q) for key,q in row.items()} for row in M],
            [{key:shift(q) for key,q in row.items()} for row in M],trim(P0),shift(P0),parameters)
        return dict(source_family=self.identity,exact_Z_cell=ends,band_x=[1,2],
            actual_profile_y0_to_y4_Z3_bounds=ranges.record(radial),actual_complete_history_y0_to_y4_Z3_bounds=ranges.record(M),
            original_beta_ordinary_derivative_global_caps=beta,original_bump_y3_y4_caps=ranges.record(shape),
            actual_mixed_recovery_Z2_magnitude_bounds=lower.ranges.record(caps),
            actual_source_radius_log_lower=logRmin,original_source_coordinates_before_time_map=True,
            range_caps_not_function_values=True,absolute_heat_cone_global_field_and_temporal_recursion_not_admitted=True)

    def velocity_functions(self):return self.functions['recovered_mixed_C2_functions']['original_cylindrical_velocity_pressure_y_rows']
    def stress_functions(self):return self.functions['recovered_mixed_C2_functions']['original_physical_signed_stress_y_rows']


def source_bindings():
    return dict(profiles_y4=current.ast_binding(profiles_y4),history_y4=current.ast_binding(history_y4),
        mixed_recover=current.ast_binding(mixed_recover),ordinary_y_product=current.ast_binding(product_row),
        physical_prefactor_shifts=current.ast_binding(shifted_rows),build=current.ast_binding(build),
        quantitative_range=current.ast_binding(CurrentMixedRepairedRecovery.quantitative_range),
        original_beta_derivative=current.ast_binding(phase.flat_source.beta_jets),
        original_beta_polynomial=current.ast_binding(phase.flat_source.beta_polynomials),
        original_beta_tail=current.ast_binding(phase.flat_source.beta_tail_bound),
        unchanged_generic_mixed_recovery=current.ast_binding(generic.GenericMomentRecovery.field_rows),
        unchanged_generic_physical_shift=current.ast_binding(generic.shifted_rows),
        accepted_C2_recovery_receipt=lower.RECEIPT,accepted_C2_recovery_receipt_sha256=sha(lower.RECEIPT))


def run(recovery_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentMixedRepairedRecovery(recovery_field=recovery_field,owner=owner,require_checked=False)
        rows=[field.quantitative_range(ends) for ends in current.current.CELLS]
        report=dict(candidate_actual_mixed_repaired_recovery_constructed=True,source_family=field.identity,
            accepted_C2_recovery_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.graph.nodes,
            actual_mixed_repaired_recovery_functions=band.control.encoded(field.functions),
            actual_four_Z_mixed_repaired_recovery_ranges=ranges.record(rows),source_bindings=source_bindings(),
            **dict.fromkeys(GATES+OPEN,False),current_numeric_point_field_oracle_installed=False,
            actual_generic_stress_cone_installed=False,physical_time_Cartesian_global_Rh_heat_and_temporal_recursion_installed=False,
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(band.control.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual repaired y4/Z3 profiles and moments, y4/Z2 velocity and y3/Z2 physical stress constructed',flush=True)
    return field


if __name__=='__main__':run()
