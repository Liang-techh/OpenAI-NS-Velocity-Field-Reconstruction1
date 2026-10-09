"""Genuine current Rsh/reference/restoration/Rm source cells and finite-N drivers.

Original full finite background kernels and the nonlinear candidate source
are retained separately. The real R110 finite-N correction is still an
unsupplied affine argument. No earlier-owner N1024 inlet is transplanted.
"""
import ast
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_long_reshape_finite_N as long
import lei_ren_part1_paper_compliant_current_original_reference_restore_functions as background

fields,base,ep=long.fields,long.base,long.ep
parameters,primitives,phase,bounds=long.parameters,long.primitives,long.phase,long.bounds
HERE,PREFIX,sha=long.HERE,long.PREFIX,long.sha
NAME=PREFIX+'current_original_reference_restore_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_reference_restore_finite_N_check.json'
GATE='current_original_Rsh_Rm_full_reference_restoration_source_and_finite_N_local_drivers_installed'
RATES,PARTITION,C0,Z=long.RATES,long.PARTITION,long.C0,long.Z
CHARTS=('reference','restoration','postrestore')


def serialized(value):
    if isinstance(value,long.IntervalTaylor):return list(value.coefficients)
    if isinstance(value,dict):return {('y%d_Z%d'%key if isinstance(key,tuple) else key):serialized(row) for key,row in value.items()}
    if isinstance(value,(tuple,list)):return [serialized(row) for row in value]
    return fields.serialized(value)


def scalar_hull(c,left,right):
    return c.mpf([min(ep(left)[0],ep(right)[0]),max(ep(left)[1],ep(right)[1])])


def source_bindings():
    """Pin both source recipes; exact coordinate identities are proved separately."""
    original=background.source_bindings();tree=ast.parse(Path(background.__file__).read_text(encoding='utf8'))
    owner=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='ActualReferenceRestoreFunctions')
    recipes={
        'reference':dict(g='self.gap*t',
            logu="-logq+((self.T/10-self.logC-self.logP)*(1-t)-c.mpf('.8')*t)",
            logR='c.ln(110)+self.T*(1-t)+(self.logref-8)*t'),
        'restoration':dict(alpha='[1-cutoff[0]]+[-cutoff[k]*math.factorial(k) for k in range(1,5)]',
            V="self.V if q==0 else self.Vref if q==1 else f.add(self.Vref,f.scale(self.E,alpha[0]))",
            offset='-8+t',logu='-logarithm(1+self.z*self.z)+offset/10',logR='c.ln(110)+self.logref+offset'),
        'postrestore':dict(s='x+7',logu='-logarithm(1+self.z*self.z)+x/10',logR='c.ln(110)+self.logref+x')}
    assignments=[]
    for name,wanted in recipes.items():
        fn=next(node for node in owner.body if getattr(node,'name',None)==name)
        for variable,expression in wanted.items():
            target=ast.dump(ast.parse(expression,mode='eval').body)
            actual=[node.value for node in ast.walk(fn) if isinstance(node,ast.Assign)
                and any(isinstance(item,ast.Name) and item.id==variable for item in node.targets)]
            if len(actual)!=1 or ast.dump(actual[0])!=target:
                raise ValueError('Original current reference/restoration source changed: '+name+'.'+variable)
            assignments.append(name+'.'+variable+'='+expression)
    own=ast.parse(Path(__file__).read_text(encoding='utf8'))
    fn=next(node for node in own.body if getattr(node,'name',None)=='background_cell')
    wanted=("g=op.gap*t", "logR=c.ln(110)+op.T+g",
        "logu=-background.logarithm(1+op.z*op.z)+(op.T/10-op.logC-op.logP)*(1-t)-c.mpf('.8')*t",
        "radius=f.factor((0,5,0,0,0),10*op.logC+c.ln(110)-8-op.gap*(1-t))",
        "K=background.original.restoration_kernels(c,t,cells=128)",
        "V=f.add(op.Vref,f.scale(op.E,alpha[0]))", "Vy=f.scale(op.E,alpha[1])",
        "logu=-background.logarithm(1+op.z*op.z)+offset/10",
        "radius=f.factor((0,5,0,0,0),10*op.logC+c.ln(110)+offset)")
    for statement in wanted:
        target=ast.dump(ast.parse(statement).body[0])
        if sum(ast.dump(node)==target for node in ast.walk(fn) if isinstance(node,ast.Assign))!=1:
            raise ValueError('Current interval source identity changed: '+statement)
    return dict(original_full_centered_and_kernel_source_bindings=original,
        original_current_wrapper_assignment_bindings=assignments,new_interval_assignment_bindings=wanted,
        coordinate_equivalence_requires_independent_symbolic_receipt=True,
        exact_physical_E_y_equals_E_over10_source_identity=True,
        axial_velocity_forcing_is_V_Rsh_minus4Z_not_shear_amplitude=True,
        full_original_restoration_kernel_function_sha256=sha(Path(background.original.__file__).name))


def shapes_from_centered(op,centered):
    f,c=op.flow,op.c;zero=[f.scalar(0)]*6
    unit=lambda value:[f.scalar(value)]+zero[1:]
    H=f.add(centered['angular_error'],unit(c.mpf(5)/8))
    return dict(theta=H,theta_z=f.add(f.scale(f.multiply(op.zrows,H),4),centered['mixed_error']),
        mean=f.add(op.Vref,centered['mean_error']),
        axial=f.add(f.scale(f.multiply(op.zrows,op.zrows),16),
            f.scale(f.multiply(op.zrows,centered['mean_error']),8),centered['axial_square']),
        swirl=f.add(centered['swirl_error'],unit(c.mpf(5)/6)),pressure=f.add(centered['pressure_error'],unit(5)))


def background_cell(op,chart,t):
    """Interval function of the original full source, never endpoint sampling."""
    f,c=op.flow,op.c;t=c.mpf(t)
    if chart not in CHARTS or ep(t)[0]<0 or ep(t)[1]>1:
        raise ValueError('Original reference/restoration/postrestore phase cell in[0,1] required')
    zero=[f.scalar(0)]*6;K=None
    if chart=='reference':
        g=op.gap*t;decays={name:op.scalar_exp(-c.mpf(rate)*g) for name,rate in background.RATES.items()}
        initial=op.centered
        centered=dict(mean_error=f.add(op.E,f.scale(f.add(initial['mean_error'],f.scale(op.E,-1)),decays['mean_error'])),
            angular_error=f.scale(initial['angular_error'],decays['angular_error']),
            mixed_error=f.add(f.scale(initial['mixed_error'],decays['mixed_error']),
                f.scale(op.E,(f.scalar(1)-decays['mixed_error'])*(c.mpf(5)/8))),
            axial_square=f.add(op.E2,f.scale(f.add(initial['axial_square'],f.scale(op.E2,-1)),decays['axial_square'])),
            swirl_error=f.scale(initial['swirl_error'],decays['swirl_error']),
            pressure_error=f.scale(initial['pressure_error'],decays['pressure_error']))
        V,Vy=op.V,zero;alpha=[c.mpf(1),c.mpf(0)]
        logu=-background.logarithm(1+op.z*op.z)+(op.T/10-op.logC-op.logP)*(1-t)-c.mpf('.8')*t
        radius=f.factor((0,5,0,0,0),10*op.logC+c.ln(110)-8-op.gap*(1-t))
        logR=c.ln(110)+op.T+g;length=op.gap
    else:
        initial=(op.reference((1,1))['actual_centered_histories'] if chart=='restoration'
                 else op.restoration((1,1))['actual_centered_histories'])
        decays={name:op.scalar_exp(-c.mpf(rate)*t) for name,rate in background.RATES.items()}
        centered={name:f.scale(row,decays[name]) for name,row in initial.items()}
        if chart=='restoration':
            # Original pure kernel accepts whole endpoint intervals and retains
            # the complete uncertain segment; no old profile owner is created.
            K=background.original.restoration_kernels(c,t,cells=128)
            centered['mean_error']=f.add(centered['mean_error'],f.scale(op.E,K['mean']))
            centered['mixed_error']=f.add(centered['mixed_error'],f.scale(op.E,K['mixed']))
            centered['axial_square']=f.add(centered['axial_square'],f.scale(op.E2,K['square']))
            cutoff=long.sigma_jets(c,t)
            alpha=[long.intersection(c,1-cutoff[0],c.mpf([0,1])),
                   long.intersection(c,-cutoff[1],c.mpf([-8,0]))]
            V=f.add(op.Vref,f.scale(op.E,alpha[0]));Vy=f.scale(op.E,alpha[1]);offset=-8+t
        else:
            alpha=[c.mpf(0),c.mpf(0)];V,Vy=op.Vref,zero;offset=-7+t
        logu=-background.logarithm(1+op.z*op.z)+offset/10
        radius=f.factor((0,5,0,0,0),10*op.logC+c.ln(110)+offset)
        logR=c.ln(110)+op.logref+offset;length=c.mpf(1)
    shapes=shapes_from_centered(op,centered)
    if chart=='reference' and ep(t)==(0,0):shapes=op.initial
    return dict(chart=chart,phase=t,original_window_length=length,actual_centered_histories=centered,
        actual_normalized_six_history_shapes=shapes,actual_velocity_V=V,actual_velocity_V_y=Vy,
        actual_centered_velocity_forcing=op.E,actual_alpha_C0_y=alpha,
        original_log_E_axial5=logu,actual_physical_radius=radius,actual_log_radius=logR,
        actual_incoming_background_decays=decays,full_original_restoration_kernels=K,
        original_P0_normalized_axial5=op.P0,source_rows_are_original_interval_functions=True,
        original_full_source_window_not_shortened=True,radial_rows_use_physical_log_radius=True,
        phase_radius_and_window_length_Z_independent=True)


def recover_cell(op,source):
    f,c=op.flow,op.c;shapes=source['actual_normalized_six_history_shapes'];logu=source['original_log_E_axial5']
    relative=long.IntervalTaylor(c,[c.mpf(0)]+list(logu.coefficients[1:])).exp()
    E=f.scale(f.jet(relative),f.factor((0,0,0,0,0),logu[0]));E2=f.multiply(E,E)
    V,Vy=source['actual_velocity_V'],source['actual_velocity_V_y'];invP2=f.factor((0,-1,0,0,0))
    hist=dict(m=shapes['mean'],h=f.multiply(E,shapes['theta']),k=f.multiply(E,shapes['theta_z']),
        e=f.add(f.scale(shapes['axial'],invP2),f.scale(f.multiply(E2,shapes['swirl']),-c.mpf('.5'))),
        p=f.scale(f.multiply(E2,shapes['pressure']),c.mpf('.5')))
    V2=f.multiply(V,V)
    dy=dict(m=f.add(V,f.scale(hist['m'],-1)),h=f.add(E,f.scale(hist['h'],-c.mpf('1.5'))),
        k=f.add(f.multiply(E,V),f.scale(hist['k'],-c.mpf('1.5'))),
        e=f.add(f.scale(V2,invP2),f.scale(hist['e'],-1),f.scale(E2,-c.mpf('.5'))),p=f.scale(E2,c.mpf('.5')))
    raw=dict(histories={name:[row,dy[name]] for name,row in hist.items()},
        velocity=dict(theta=[E,f.scale(E,c.mpf('.1'))],axial=[V,Vy]))
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=raw,original_P0_normalized_axial5=op.P0,
        geometry=dict(chart=source['chart'],phase=source['phase'],actual_physical_radius=source['actual_physical_radius'],
            actual_log_radius=source['actual_log_radius'],original_window_length=source['original_window_length'],
            phase_and_radius_Z_independent=True))
    proxy=SimpleNamespace(flow=f,c=c,reference=op,zrows=op.zrows,P0=op.P0,
        Pstar=f.factor((0,.5,0,0,0)),source_radius=source['actual_physical_radius'],correlated_C=f.scale(E,c.mpf('.8')))
    recovered=long.RECOVER(proxy,packet)
    recovered.pop('source_frame_conditional_on_same_actual_Rm_inlet')
    recovered['source_frame_conditional_on_same_current_Rsh_background']=True
    recovered['exact_E_y_equals_E_over10_and_a_four_fifths']=True
    return proxy,packet,recovered


def general_q_C0_Z(f,a,Delta,eta_log):
    """Original smooth q, with an explicit active/flat source-domain union."""
    ap=parameters.positive_source(f,a[0],'actual_reference_shear_a')
    exact=parameters.shear_q_jets(f,a,Delta,eta_log,ap)
    if exact['analytic_axial_jets_available']:
        qr={C0:exact['q_axial_coefficients'][0],Z:exact['q_axial_coefficients'][1]}
        return qr,dict(original_q_source=exact,used_original_exact_axial_branch=True,
            q_is_exact_zero=qr[C0].zero,q_Z_is_exact_zero=qr[Z].zero,
            q_globally_positive=False,selected_compatible_field_claimed=False)
    loop=exact['original_C0_lazy_cutoff'];eta=f.factor((0,0,0,0,0),eta_log)
    # This positivity is ONLY on Delta<eta. The other source domain is exact
    # flat zero. No global gamma positivity or selected source value follows.
    gamma=eta*2-Delta[0]
    if ep(gamma.coefficient)[0]<=0:gamma=gamma.positive_intersection(eta_log)
    ratio=gamma.positive_divide(a[0]*2,ap['source_log_lower']+f.c.ln(2))
    root=parameters.original.nonnegative_sqrt(ratio)
    dz=primitives.absolute(Delta[1]);az=primitives.absolute(a[1])
    derivative=primitives.absolute(root)*(dz.positive_divide(eta,eta_log)*(f.c.mpf(17)/2)
        +az.positive_divide(a[0]*2,ap['source_log_lower']+f.c.ln(2)))
    qr={C0:loop['q'],Z:bounds.symmetric(f,derivative)}
    return qr,dict(original_q_source=loop,used_original_exact_axial_branch=False,
        active_enclosure_subdomain='same source predicate Delta<eta',
        flat_enclosure_subdomain='same source predicate Delta>=eta',
        mixed_branch_union=loop['branch']=='requires_source_box_refinement',
        active_root_cover_only=root,active_q_Z_cover_only=derivative,
        active_gamma_positive_intersection_scope_only=True,active_gamma_lower_is_eta=True,
        original_global_sigma_prime_bound=8,original_flat_sigma_and_first_derivative_zero=True,
        exact_flat_q_and_q_Z_zero_only_on_flat_subdomain=True,
        derivative_of_selected_cap_not_used=True,all_domain_active=loop['branch']=='active',
        q_globally_positive=False,selected_compatible_field_claimed=False,
        original_q_Z_majorant_formula='root*(17*abs(Delta_Z)/(2*eta)+abs(a_Z)/(2*a))')


def quotient_cell(proxy,recovered,eta_log):
    f,c=proxy.flow,proxy.c;n=recovered['actual_generic_source_numerators'];E=n['E']
    pe=parameters.positive_source(f,E[0],'actual_reference_E')
    a=[f.scalar(c.mpf('.8'))]+[f.scalar(0)]*5
    b=parameters.quotient(f,n['B'],E,pe);t0=parameters.scale(b,-c.mpf(5)/4)
    b2=parameters.multiply(f,b,b);b2[0]=parameters.original.square(b[0])
    kappa=parameters.add(f,a,parameters.scale(b2,c.mpf(5)/4))
    Delta=parameters.add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*5)
    qr,qproof=general_q_C0_Z(f,a,Delta,eta_log)
    sectors=recovered['full_signed_inertial_sectors_axial4'];inertial={}
    for name,left,right in (('p1','theta_linear','theta_quadratic'),('p2','axial_linear','axial_quadratic')):
        inertial[name]=parameters.quotient(f,parameters.add(f,sectors[left],parameters.scale(sectors[right],proxy.Pstar)),E,pe)
    p2=parameters.scale(inertial['p2'],proxy.source_radius)
    roots={name:{C0:rows[0],Z:rows[1]} for name,rows in dict(a=a,t0=t0,E=E,p2=p2).items()}
    return dict(q=qr[C0],roots=roots),qr,dict(source_positive_E=pe,
        actual_correlated_a_axial5=a,actual_nonzero_b_axial5=b,actual_t0_axial5=t0,
        actual_kappa_axial5=kappa,actual_Delta_axial5=Delta,original_general_q_C0_Z=qproof,
        full_signed_inertial_before_R_axial4=inertial,full_signed_p2_axial4=p2,
        one_actual_radius_factor=proxy.source_radius,b_Z_and_all_square_product_terms_retained=True,
        q_parameters_are_definitions_not_current_owner_cone_admission=True)


def density_cell(op,source,eta_log,dstar_log,N):
    proxy,raw,recovered=recover_cell(op,source)
    roots,qr,quotients=quotient_cell(proxy,recovered,eta_log)
    got=primitives.all_u_primitive_bounds(op.flow,roots,qr,dstar_log,op.c.mpf([0,1]))
    if qr[C0].zero and qr[Z].zero:
        # Original inverse: q=0 implies Phi=psi/(2pi), A=0, B=0 even
        # when t0!=0. Flat smoothness also gives exact first Z zeros.
        got['values']={key:op.flow.scalar(0) for key in got['values']}
        got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
    E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
    density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
    return dict(original_background_interval_source=source,original_raw_current_source=raw,
        original_generic_source=recovered,original_full_source_quotients=quotients,
        original_roots=roots['roots'],original_q_rows=qr,original_primitive_values=got['values'],
        original_primitive_proof=got['record'],original_signed_five_density_C0_Z=density,
        actual_radius_phase_definition='frac(N*(logR-logRa-hb*s_c/2))',
        actual_radius_phase_cover=op.c.mpf([0,1]),actual_radius_phase_Z_exact_zero=True,
        actual_phase_full_period_cover_not_selected_inverse=True,
        finite_N_incoming_correction_supplied=False,**dict.fromkeys(fields.previous.OPEN,False))


class OriginalReferenceRestoreFiniteN:
    mode='current_original_full_reference_restoration_postrestore_source_and_local_finite_N_Duhamel'
    def __init__(self,dps=500,require_checked=True):
        self.long=long.OriginalLongReshapeFiniteN(dps);self.parameters=self.long.parameters
        self.reference=self.long.reference_wrapper;self.c=self.long.c;self.family=self.long.family
        self.hashes=dict(self.long.hashes);self.cache={}
        self.saved_long=json.loads(gzip.decompress((HERE/long.NAME).read_bytes()))
        if self.saved_long['source_family']!=self.family or not self.saved_long[long.GATE]:
            raise ValueError('Same current accepted long-window source required')
        self.bindings=source_bindings()
        for module in (background,background.original,primitives,phase,bounds):
            fields.previous.bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted current reference/restoration finite-N receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def query(self,label,chart,left,right=None,N=257):
        N=phase.candidate_N(N);lo=long.fraction(left);hi=lo if right is None else long.fraction(right)
        if hi<lo:raise ValueError('Ordered actual reference/restoration phase cell required')
        key=(label,chart,left,right,N)
        if key in self.cache:return self.cache[key]
        op=self.reference.owner(label);f,c=op.flow,op.c
        with mp.workdps(c.dps+40):
            t=scalar_hull(c,c.mpf(lo.numerator)/lo.denominator,c.mpf(hi.numerator)/hi.denominator)
            source=background_cell(op,chart,t)
            packet=density_cell(op,source,self.parameters.eta_log,self.parameters.dstar_log,N)
            packet.update(source_family=self.family,source_frame=label,candidate_N=N,chart=chart,
                exact_phase_left=left,exact_phase_right=left if right is None else right,
                exact_common_P0_axial5=op.P0,source_function_hash=sha(Path(background.__file__).name),
                current_Rsh_source_hash=sha(background.previous.NAME),
                current_source_owner_basis_ledger_and_P0_preserved=True)
        self.cache[key]=packet;return packet

    def contribution(self,label,N=257):
        N=phase.candidate_N(N)
        if self.saved_long['candidate_N']!=N:
            raise ValueError('Accepted long driver candidate N must match; no cross-N transplant')
        op=self.reference.owner(label);f,c=op.flow,op.c
        windows={};local={name:[f.scalar(0),f.scalar(0)] for name in RATES}
        with mp.workdps(c.dps+40):
            for chart in CHARTS:
                length=op.gap if chart=='reference' else c.mpf(1)
                total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
                for left,right in zip(PARTITION,PARTITION[1:]):
                    source=self.query(label,chart,left,right,N);l,r=long.fraction(left),long.fraction(right)
                    width=length*c.mpf((r-l).numerator)/(r-l).denominator
                    suffix=length*c.mpf((1-r).numerator)/(1-r).denominator;rows={};weights={}
                    for name,rate in RATES.items():
                        mass,decay,tail=long.kernel_weight(f,width,suffix,rate)
                        density=source['original_signed_five_density_C0_Z']
                        pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                              for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):total[name][n]+=row
                        rows[name]=pair;weights[name]=dict(full_positive_mass=mass,cell_memory_decay=decay,
                            downstream_suffix_decay=tail,own_rate=str(rate),physical_log_measure_once=True)
                    cells.append(dict(source=source,actual_log_radius_width=width,
                        actual_log_radius_suffix=suffix,signed_cell_driver_C0_Z=rows,original_own_rate_weights=weights))
                memory={name:f.scalar(1) if not rate else f.factor((0,0,0,0,0),-length*c.mpf(rate.numerator)/rate.denominator)
                        for name,rate in RATES.items()}
                for name in RATES:
                    local[name]=[local[name][n]*memory[name]+total[name][n] for n in range(2)]
                windows[chart]=dict(actual_full_source_cells=cells,original_full_log_radius_length=length,
                    actual_local_signed_driver_C0_Z=total,actual_nonzero_incoming_memory=memory)
            accepted=self.saved_long['frames'][label]
            if accepted['candidate_N']!=N:raise ValueError('Accepted long driver candidate N must match; no cross-N transplant')
            # Derived display/log-majorant fields depend on evaluation context;
            # bind the actual coefficient and formal scale source tuples.
            expression=lambda row:{key:row[key] for key in ('coefficient_interval','formal_positive_scale','exact_zero')}
            if ([expression(row) for row in base.encoded(fields.serialized(op.P0))]
                    !=[expression(row) for row in accepted['exact_common_P0_axial5']]):
                raise ValueError('Same exact canonical P0 source expression required when composing current source windows')
            restore=lambda row:background.endpoint.restore_row(f,row)
            long_driver={name:[restore(row) for row in rows] for name,rows in accepted['actual_local_finite_N_driver_C0_Z'].items()}
            downstream_memory={name:f.scalar(1) if not rate else f.factor((0,0,0,0,0),-(op.gap+2)*c.mpf(rate.numerator)/rate.denominator)
                               for name,rate in RATES.items()}
            complete={name:[long_driver[name][n]*downstream_memory[name]+local[name][n] for n in range(2)] for name in RATES}
            inlet_memory={name:f.scalar(1) if not rate else f.factor((0,0,0,0,0),-(op.T+op.gap+2)*c.mpf(rate.numerator)/rate.denominator)
                          for name,rate in RATES.items()}
        return dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
            actual_source_windows=windows,actual_Rsh_Rm_local_driver_C0_Z=local,
            actual_R110_Rm_local_driver_C0_Z=complete,actual_R110_Rm_incoming_memory=inlet_memory,
            accepted_long_driver_rehydrated_in_same_live_source_algebra=True,
            exact_affine_boundary_formula='deltaH(Rm)=exp(-lambda*(T+gap+2))*deltaH(R110)+all_local_source_drivers',
            genuine_current_finite_N_R110_correction_still_unsupplied=True,
            actual_finite_N_Rm_incoming_correction_supplied=False,no_saved_N1024_inlet_transplanted=True,
            all_downstream_source_windows_integrated_not_homogeneous_only=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalReferenceRestoreFiniteN(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Current original full reference/restoration/postrestore finite-N drivers',label,flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_interval_source_bindings=owner.bindings,
        frames=serialized(frames),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
