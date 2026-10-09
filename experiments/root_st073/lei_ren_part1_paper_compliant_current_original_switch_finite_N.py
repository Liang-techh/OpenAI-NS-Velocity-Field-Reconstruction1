"""Current original R100/R110 switch source cells and finite-N correction drivers.

The microscopic physical width stays factored. Cumulative background
histories keep their actual inlet memory. Real finite-N R100 correction
is still a separate unsupplied affine argument, never a background row.
"""
import ast
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_reference_restore_finite_N as downstream
import lei_ren_part1_paper_compliant_current_original_first_switch_functions as first
import lei_ren_part1_paper_compliant_current_original_second_switch_R110 as second

long=downstream.long;parameters,primitives,phase,bounds=long.parameters,long.primitives,long.phase,long.bounds
fields,base,ep=long.fields,long.base,long.ep
HERE,PREFIX,sha=long.HERE,long.PREFIX,long.sha
NAME=PREFIX+'current_original_switch_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_switch_finite_N_check.json'
GATE='current_original_R100_R110_whole_switch_source_and_finite_N_local_drivers_installed'
RATES,PARTITION,C0,Z=long.RATES,long.PARTITION,long.C0,long.Z
CHARTS=('first_switch','second_switch','post_power')


def source_bindings():
    original=second.source_bindings();denominator=first.fixed_comparison_denominator_binding()
    from lei_ren_part1_paper_compliant_first_switch_leading import switch_control_source_bridge
    switch_controls=switch_control_source_bridge()
    tree=ast.parse(Path(first.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    statement="phi_phase=f.scale(f.multiply(Dend,at['phi']),-h*h*c.mpf('.5'))"
    target=ast.dump(ast.parse(statement).body[0])
    if sum(ast.dump(n)==target for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
        raise ValueError('Original first-switch phi phase ODE changed')
    recipes=dict(prefix_cell=("partial=h*ds*average",
        "a=f.scale(D,h)",
        "a=f.add(f.scale(D,h*complement),[f.scalar(c.mpf('.8')*sig)]+[f.scalar(0)]*5)",
        "histories[name]=f.add(seed['actual_six_histories'][name],f.scale(seed['actual_six_histories'][name],dm),f.scale(source[name],partial))"),
        post_cell=("length=f.scalar(Y)-f.h*2",
        "powers[name]=f.scalar(c.exp(-c.mpf(p)*Y*t))*first_op.scalar_series(2*c.mpf(p),t)[0]",
        "a=[f.scalar(c.mpf('.8'))]+[f.scalar(0)]*5"),
        recover_source=("amplitude=rootR*F*c.sqrt(2)",
        "Ey=f.scale(dress(f.add(source['phi_y'],f.scale(phi,c.mpf('.5')))),amplitude)",
        "C=f.multiply(a,E)"))
    own=ast.parse(Path(__file__).read_text(encoding='utf8'))
    for name,texts in recipes.items():
        recipe=next(n for n in own.body if isinstance(n,ast.FunctionDef) and n.name==name)
        for text in texts:
            wanted=ast.dump(ast.parse(text).body[0])
            if sum(ast.dump(n)==wanted for n in ast.walk(recipe) if isinstance(n,ast.Assign))!=1:
                raise ValueError('Current whole-cell original source recipe changed: '+text)
    post=next(n for n in own.body if isinstance(n,ast.FunctionDef) and n.name=='post_cell')
    wanted=ast.dump(ast.parse('[f.scalar(0)]*6',mode='eval').body)
    if sum(n.arg=='V_y' and ast.dump(n.value)==wanted for n in ast.walk(post) if isinstance(n,ast.keyword))!=1:
        raise ValueError('Current original post-power exact zero V_y changed')
    source_path=HERE/(PREFIX+'inner_switch_profiles.py')
    original_tree=ast.parse(source_path.read_text(encoding='utf8'))
    phase_source=next(n for n in ast.walk(original_tree) if isinstance(n,ast.FunctionDef) and n.name=='phase')
    original_angular_source='-.5*hb^2*[integral_0^min(phase,1)Dbar(100exp(hb*t),Z)dt + integral_0^max(phase-1,0)(1-sigma(t))*Dbar(100exp(hb*(1+t)),Z)dt] -.4*hb*integral_0^max(phase-1,0)sigma(t)dt'
    original_axial_source='V100 - hb^2*integral_0^min(phase,1)(1-sigma(t))*(phi_actual/barphi)*drive(100exp(hb*t),Z)dt'
    for key,text in (('angular_source',original_angular_source),('exact_Uz_source',original_axial_source)):
        values=[n.value.value for n in ast.walk(phase_source) if isinstance(n,ast.keyword)
                and n.arg==key and isinstance(n.value,ast.Constant)]
        if values!=[text]:raise ValueError('Original cumulative switch defining source changed: '+key)
    return dict(original_second_switch_power_and_fixed_comparison_bindings=original,
        original_first_second_and_post_switch_control_bindings=switch_controls,
        current_original_post_power_V_y_exact_zero_binding='[f.scalar(0)]*6',
        original_fixed_comparison_denominator=denominator,first_phi_phase_ODE_binding=statement,
        original_physical_endpoint_binding=first.endpoint.source_bindings(),
        exact_conversion_requires_independent_symbolic_receipt=True,
        original_generic_signed_recovery_binding=long.RECOVERY_BINDING,current_whole_cell_recipe_assignments=recipes,
        original_switch_defining_angular_source=original_angular_source,
        original_switch_defining_axial_source=original_axial_source,
        phase_y_Jacobian='dy=hb*ds; ordinary phi_y=phi_s/hb',
        microscopic_width_is_formal_source_not_numerical_cap=True)


def positive_jet(f,row,role):return parameters.positive_source(f,row[0],role)


def canonical_expression(rows):
    return [{key:row[key] for key in ('coefficient_interval','formal_positive_scale','exact_zero')}
            for row in base.encoded(fields.serialized(rows))]


def guard_saved_source(family,label,N,P0,accepted):
    if (accepted['source_family']!=family or accepted['source_frame']!=label
            or accepted['candidate_N']!=N):raise ValueError('Same family, frame and finite N required for composition')
    expression=[{key:row[key] for key in ('coefficient_interval','formal_positive_scale','exact_zero')}
                for row in accepted['exact_common_P0_axial5']]
    if canonical_expression(P0)!=expression:
        raise ValueError('Same exact canonical analytic P0 source expression required for composition')


def bounded_exponent_cover(f,got,N):
    """Outward A range only; keep the factored source and genuine A_Z."""
    row=got['values']['A']
    if row.zero:return got
    lower,upper=ep(row.scale.evaluate())
    if lower < -1000 <= upper:
        log_upper=ep(row.record()['log_absolute_upper'])[1]
        if not mp.isfinite(log_upper) or abs(log_upper)>1000:
            raise ValueError('Bounded exponent cover must have a finite materializable log upper')
        if log_upper>ep(f.c.ln(f.c.mpf(N)))[0]:
            raise ValueError('Actual primitive exponent requires a larger candidate N')
        cap=f.c.exp(f.c.mpf(log_upper))
        got['record']['original_A_C0_before_outward_exponent_cover']=row.record()
        got['values']['A']=f.scalar(f.c.mpf([-ep(cap)[1],ep(cap)[1]]))
        got['record']['bounded_A_C0_cover_only_not_source_or_A_Z_replacement']=True
        got['record']['density_A_C0_outward_enclosure']=got['values']['A'].record()
    return got


def prefix_cell(first_op,second_op,chart,left,right):
    """True partial-cell ODE integrals from a current cumulative left endpoint."""
    f,c=first_op.flow,first_op.c;lo,hi=long.fraction(left),long.fraction(right)
    if chart not in ('first_switch','second_switch') or not 0<=lo<=hi<=1:raise ValueError('Ordered original switch cell in [0,1] required')
    l=c.mpf(lo.numerator)/lo.denominator;r=c.mpf(hi.numerator)/hi.denominator
    t=downstream.scalar_hull(c,l,r);ds=c.mpf([0,ep(r-l)[1]])
    sig=long.intersection(c,first.sigma_jets(c,t)[0],c.mpf([0,1]));complement=1-sig
    h=f.h;seed=(first_op if chart=='first_switch' else second_op).evaluate(left)
    if chart=='first_switch':
        angular=first_op.angular(t);delta,G,_=first_op.force_delta(t,angular)
        phi=angular['phi'];mass=c.mpf([0,ep(first_op.mass(lo,hi))[1]])
        V=f.add(seed['actual_fields']['V'],*[f.scale(f.add(first_op.force0[part],delta[part]),
            -h*h*first_op.scales[part]*mass) for part in fields.PARTS])
        Vy=f.add(*[f.scale(f.add(first_op.force0[part],delta[part]),
            -h*first_op.scales[part]*complement) for part in fields.PARTS])
        D=f.add(*[f.scale(row,first_op.scalar_series(1-j,t)[0]) for j,row in enumerate(first_op.D)])
        dp=positive_jet(f,D,'current_first_Dbar');a=f.scale(D,h)
        phi_y=f.scale(f.multiply(a,phi),-c.mpf('.5'));phase_offset=t
    else:
        mass=c.mpf([0,ep(first_op.mass(lo,hi))[1]])
        prefix=seed['complement_prefix_mass']+mass
        prefix=long.intersection(c,prefix,c.mpf([0,ep(t)[1]]))
        D=f.add(second_op.D0,second_op.Ddelta(t))
        WD=f.add(seed['complete_weighted_original_D_integral'],f.scale(D,mass))
        phi=second_op.angular(t,prefix,WD)['phi'];V=second_op.inlet_fields['V'];Vy=[f.scalar(0)]*6
        dp=positive_jet(f,D,'current_second_Dbar')
        a=f.add(f.scale(D,h*complement),[f.scalar(c.mpf('.8')*sig)]+[f.scalar(0)]*5)
        # The exact convex interpolation is positive even when an arithmetic
        # cover loses its microscopic lower end next to the .8 branch.
        lower=min(ep(dp['source_log_lower']+f.logs[0])[0],ep(c.ln(c.mpf('.8')))[0])
        if ep(a[0].coefficient)[0]<=0:a[0]=a[0].positive_intersection(c.mpf(lower))
        phi_y=f.scale(f.multiply(a,phi),-c.mpf('.5'));phase_offset=1+t
    source=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),
        B=f.multiply(phi,phi),C=f.multiply(phi,phi))
    histories={};weights={}
    for name,rate in first.moments.RATES.items():
        dm,err=first_op.scalar_series(-rate,ds,minus_one=True)
        average,mass_error=first_op.scalar_series(-rate,ds,average=True)
        partial=h*ds*average
        histories[name]=f.add(seed['actual_six_histories'][name],
            f.scale(seed['actual_six_histories'][name],dm),f.scale(source[name],partial))
        weights[name]=dict(actual_partial_positive_mass=partial,actual_memory_minus_one=dm,
            complete_decay_error=err,complete_mass_error=mass_error,cumulative_left_history_preserved=True)
    radial=first_op.scalar_series(1,phase_offset)[0]
    root_radial=first_op.scalar_series(c.mpf('.5'),phase_offset)[0]
    return dict(chart=chart,phase=t,fields=dict(phi=phi,V=V),histories=histories,
        phi_y=phi_y,V_y=Vy,correlated_a_axial5=a,actual_Dbar_source=D,actual_positive_Dbar_proof=dp,
        radius=f.scalar(100)*radial,root_radius=f.scalar(10)*root_radial,
        physical_log_radius_over100=h*phase_offset,window_length=h,
        actual_current_left_endpoint=left,current_left_source_packet=seed,
        true_partial_cell_history_weights=weights,original_sigma=sig,
        physical_phase_measure='dy=hb*ds',source_interval_functions_not_endpoint_hulls=True)


def post_cell(first_op,second_op,t):
    f,c=first_op.flow,first_op.c;t=c.mpf(t)
    if ep(t)[0]<0 or ep(t)[1]>1:raise ValueError('Original post-power phase in[0,1] required')
    Y=c.ln(c.mpf(11)/10);length=f.scalar(Y)-f.h*2
    parameters.positive_source(f,length,'actual_R2_R110_positive_length')
    seed=second_op.evaluate((1,1));phi2=seed['actual_fields']['phi'];V=seed['actual_fields']['V']
    powers={}
    for name,p in (('one',1),('two',2),('angular',c.mpf('.4')),('squared',c.mpf('.8'))):
        powers[name]=f.scalar(c.exp(-c.mpf(p)*Y*t))*first_op.scalar_series(2*c.mpf(p),t)[0]
    def positive(value):
        ordinary=f.ordinary_cover(value);lo,hi=ep(ordinary)
        if hi<0:raise ValueError('Original theta kernel must enclose a nonnegative function')
        return f.scalar(c.mpf([max(mp.mpf(0),lo),hi]))
    angular=positive((powers['angular']-powers['two'])*(c.mpf(5)/4))
    swirl=positive((powers['squared']-powers['two'])*(c.mpf(5)/6))
    pressure=positive((powers['squared']-powers['one'])*5);axial=positive(1-powers['one'])
    initial=seed['actual_six_histories'];phi=f.scale(phi2,powers['angular'])
    histories=dict(H=f.add(f.scale(initial['H'],powers['two']),f.scale(phi2,angular)),
        M=f.add(f.scale(initial['M'],powers['one']),f.scale(V,axial)),
        K=f.add(f.scale(initial['K'],powers['two']),f.scale(f.multiply(phi2,V),angular)),
        A=f.add(f.scale(initial['A'],powers['one']),f.scale(f.multiply(V,V),axial)),
        B=f.add(f.scale(initial['B'],powers['two']),f.scale(f.multiply(phi2,phi2),swirl)),
        C=f.add(f.scale(initial['C'],powers['one']),f.scale(f.multiply(phi2,phi2),pressure)))
    radius=f.scalar(100*c.exp(Y*t))*first_op.scalar_series(2,1-t)[0]
    root_radius=f.scalar(10*c.exp(Y*t/2))*first_op.scalar_series(1,1-t)[0]
    if ep(t)==(1,1):radius,root_radius=f.scalar(110),f.scalar(c.sqrt(110))
    a=[f.scalar(c.mpf('.8'))]+[f.scalar(0)]*5
    return dict(chart='post_power',phase=t,fields=dict(phi=phi,V=V),histories=histories,
        phi_y=f.scale(phi,-c.mpf('.4')),V_y=[f.scalar(0)]*6,correlated_a_axial5=a,
        radius=radius,root_radius=root_radius,physical_log_radius_over100=f.h*2+length*t,
        window_length=length,current_second_exit_seed=seed,original_theta_powers=powers,
        original_positive_theta_kernels=dict(angular=angular,swirl=swirl,pressure=pressure,axial=axial),
        physical_phase_measure='dy=(log(1.1)-2hb)*dt',
        source_interval_functions_not_endpoint_hulls=True,nonzero_micro_radius_terms_retained=True)


def recover_source(reference,inlet,source):
    f,c=reference.flow,reference.c
    ratio=inlet['original_F0_derivative_ratios_ordinary'];ratio2=inlet['original_F0_squared_derivative_ratios_ordinary']
    phi,V=source['fields']['phi'],source['fields']['V'];s=source['histories'];R,rootR=source['radius'],source['root_radius']
    F=f.factor((0,-.5,.5,0,0));F2=f.factor((0,-1,1,0,0));S2inv=f.factor((0,-1,0,0,0))
    dress=lambda row:f.multiply(row,ratio)
    amplitude=rootR*F*c.sqrt(2);E=f.scale(dress(phi),amplitude)
    histories=dict(m=s['M'],h=f.scale(dress(s['H']),rootR*F*(1/c.sqrt(2))),
        k=f.scale(dress(s['K']),rootR*F*(1/c.sqrt(2))),
        e=f.add(f.scale(s['A'],S2inv),f.scale(f.multiply(s['B'],ratio2),-R*F2)),
        p=f.scale(f.multiply(s['C'],ratio2),R*F2))
    E2=f.multiply(E,E);V2=f.multiply(V,V)
    dy=dict(m=f.add(V,f.scale(histories['m'],-1)),h=f.add(E,f.scale(histories['h'],-c.mpf('1.5'))),
        k=f.add(f.multiply(E,V),f.scale(histories['k'],-c.mpf('1.5'))),
        e=f.add(f.scale(V2,S2inv),f.scale(histories['e'],-1),f.scale(E2,-c.mpf('.5'))),p=f.scale(E2,c.mpf('.5')))
    Ey=f.scale(dress(f.add(source['phi_y'],f.scale(phi,c.mpf('.5')))),amplitude)
    a=source['correlated_a_axial5'];C=f.multiply(a,E)
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=dict(
        histories={name:[row,dy[name]] for name,row in histories.items()},
        velocity=dict(theta=[E,Ey],axial=[V,source['V_y']])),original_P0_normalized_axial5=reference.P0,
        geometry=dict(chart=source['chart'],phase=source['phase'],actual_physical_radius=R,
            actual_log_radius_over100=source['physical_log_radius_over100'],phase_and_radius_Z_independent=True))
    proxy=SimpleNamespace(flow=f,c=c,reference=reference,zrows=reference.zrows,P0=reference.P0,
        Pstar=f.factor((0,.5,0,0,0)),source_radius=R,correlated_C=C)
    recovered=long.RECOVER(proxy,packet);recovered.pop('source_frame_conditional_on_same_actual_Rm_inlet')
    recovered['source_frame_conditional_on_same_current_R100_background']=True
    return proxy,packet,recovered,a


def general_quotients(proxy,recovered,a,eta_log):
    f,c=proxy.flow,proxy.c;n=recovered['actual_generic_source_numerators'];E=n['E']
    pe=positive_jet(f,E,'actual_switch_E');pa=positive_jet(f,a,'actual_switch_a')
    b=parameters.quotient(f,n['B'],E,pe);t0=parameters.quotient(f,parameters.scale(b,-1),a,pa)
    square=parameters.multiply(f,b,b);square[0]=parameters.original.square(b[0])
    kappa=parameters.add(f,a,parameters.quotient(f,square,a,pa))
    Delta=parameters.add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*5)
    loop=parameters.original.q_enclosure(a[0],Delta[0],eta_log,pa['source_log_lower'])
    q=loop['q'];eta=f.factor((0,0,0,0,0),eta_log);scope={}
    if loop['branch']=='flat':qz=f.scalar(0);scope=dict(exact_flat_q_and_q_Z_zero=True)
    elif ep(Delta[0].coefficient)[1]<=0:
        gamma=eta*2-Delta[0];pg=parameters.positive_source(f,gamma,'actual_negative_Delta_gamma')
        qz=q*(parameters.scale(Delta,-1)[1].positive_divide(gamma,pg['source_log_lower'])
            -a[1].positive_divide(a[0],pa['source_log_lower']))*c.mpf('.5')
        scope=dict(original_sigma_identically_one_on_negative_Delta=True,
            signed_q_Z_from_true_root_identity_not_cap_derivative=True)
    else:
        gamma=eta*2-Delta[0]
        if ep(gamma.coefficient)[0]<=0:gamma=gamma.positive_intersection(eta_log)
        root=parameters.original.nonnegative_sqrt(gamma.positive_divide(a[0]*2,pa['source_log_lower']+c.ln(2)))
        majorant=primitives.absolute(root)*(primitives.absolute(Delta[1]).positive_divide(eta,eta_log)*(c.mpf(17)/2)
            +primitives.absolute(a[1]).positive_divide(a[0]*2,pa['source_log_lower']+c.ln(2)))
        qz=bounds.symmetric(f,majorant)
        scope=dict(active_enclosure_subdomain='same source Delta<eta',flat_enclosure_subdomain='same source Delta>=eta',
            mixed_branch_union=loop['branch']=='requires_source_box_refinement',active_gamma_lower_only_on_active_subdomain=True,
            active_root_cover_only=root,active_q_Z_cover_only=majorant,q_globally_positive=False,
            original_flat_sigma_and_first_derivative_zero=True,derivative_of_cap_not_used=True)
    sectors=recovered['full_signed_inertial_sectors_axial4'];inertial={}
    for name,l,r in (('p1','theta_linear','theta_quadratic'),('p2','axial_linear','axial_quadratic')):
        inertial[name]=parameters.quotient(f,parameters.add(f,sectors[l],parameters.scale(sectors[r],proxy.Pstar)),E,pe)
    p2=parameters.scale(inertial['p2'],proxy.source_radius)
    roots={name:{C0:row[0],Z:row[1]} for name,row in dict(a=a,t0=t0,E=E,p2=p2).items()};qr={C0:q,Z:qz}
    return dict(q=q,roots=roots),qr,dict(actual_positive_E=pe,actual_positive_a=pa,
        actual_a_axial5=a,actual_b_axial5=b,actual_t0_axial5=t0,actual_kappa_axial5=kappa,actual_Delta_axial5=Delta,
        original_q_C0=loop,original_q_Z_scope=scope,full_signed_inertial_before_R_axial4=inertial,
        full_signed_p2_axial4=p2,one_actual_radius_factor=proxy.source_radius,
        source_parameters_not_new_owner_cone_or_global_N_admission=True)


def own_weights(first_op,chart,left,right,rate):
    f,c=first_op.flow,first_op.c;l,r=long.fraction(left),long.fraction(right)
    width=c.mpf((r-l).numerator)/(r-l).denominator;suffix=c.mpf((1-r).numerator)/(1-r).denominator
    if ep(width)[0]<=0:raise ValueError('Positive actual source cell width required')
    if chart in ('first_switch','second_switch'):
        w=f.h*width
        if not rate:return w,f.scalar(1),f.scalar(1)
        lam=c.mpf(rate.numerator)/rate.denominator
        mass=w*first_op.scalar_series(-lam,width,average=True)[0]
        return mass,first_op.scalar_series(-lam,width)[0],first_op.scalar_series(-lam,suffix)[0]
    Y=c.ln(c.mpf(11)/10);L=f.scalar(Y)-f.h*2;w=L*width
    if not rate:return w,f.scalar(1),f.scalar(1)
    lam=c.mpf(rate.numerator)/rate.denominator
    mass=w*long.exp_average(c,-lam*f.ordinary_cover(w))
    decay=f.scalar(c.exp(-lam*Y*width))*first_op.scalar_series(2*lam,width)[0]
    tail=f.scalar(c.exp(-lam*Y*suffix))*first_op.scalar_series(2*lam,suffix)[0]
    return mass,decay,tail


class OriginalSwitchFiniteN:
    mode='current_original_complete_R100_R110_switch_source_and_finite_N_local_affine_drivers'
    def __init__(self,dps=500,require_checked=True):
        self.downstream=downstream.OriginalReferenceRestoreFiniteN(dps)
        self.second=self.downstream.long.long_wrapper.upstream;self.first=self.second.upstream;self.r100=self.first.upstream
        self.c=self.downstream.c;self.family=self.downstream.family;self.hashes=dict(self.downstream.hashes);self.cache={}
        self.saved_downstream=json.loads(gzip.decompress((HERE/downstream.NAME).read_bytes()))
        if self.saved_downstream['source_family']!=self.family or not self.saved_downstream[downstream.GATE]:
            raise ValueError('Same current downstream driver owner required')
        self.bindings=source_bindings()
        for name,digest in self.bindings['original_first_second_and_post_switch_control_bindings']['input_hashes'].items():
            fields.previous.bind(self.hashes,name,digest)
        for module in (first,second):fields.previous.bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual switch finite-N receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def query(self,label,chart,left,right=None,N=257):
        N=phase.candidate_N(N);right=left if right is None else right;key=(label,chart,left,right,N)
        if key in self.cache:return self.cache[key]
        first_op,second_op=self.first.owner(label),self.second.owner(label)
        ref=self.downstream.reference.owner(label);f,c=ref.flow,ref.c
        _,_,inlet=self.r100.owner(label)
        if first_op.flow is not f or second_op.flow is not f:raise ValueError('One live current source flow required')
        if canonical_expression(inlet['pressure_axis_over_Pstar_squared_axial5'])!=canonical_expression(ref.P0):
            raise ValueError('Same exact original R100 and reference analytic P0 source required')
        with mp.workdps(c.dps+40):
            l,r=long.fraction(left),long.fraction(right)
            if r<l:raise ValueError('Ordered source phase required')
            t=downstream.scalar_hull(c,c.mpf(l.numerator)/l.denominator,c.mpf(r.numerator)/r.denominator)
            source=(post_cell(first_op,second_op,t) if chart=='post_power'
                    else prefix_cell(first_op,second_op,chart,left,right))
            proxy,raw,recovered,a=recover_source(ref,inlet,source)
            roots,qr,proof=general_quotients(proxy,recovered,a,self.downstream.parameters.eta_log)
            got=primitives.all_u_primitive_bounds(f,roots,qr,self.downstream.parameters.dstar_log,c.mpf([0,1]))
            got=bounded_exponent_cover(f,got,N)
            if qr[C0].zero and qr[Z].zero:
                got['values']={key:f.scalar(0) for key in got['values']}
                got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
            E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
            density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,chart=chart,
                exact_common_P0_axial5=ref.P0,actual_background_source=source,original_raw_source=raw,
                original_generic_source=recovered,original_full_source_quotients=proof,
                original_roots=roots['roots'],original_q_C0_Z=qr,original_primitive_values=got['values'],
                original_primitive_proof=got['record'],original_signed_five_density_C0_Z=density,
                actual_phase_definition='frac(N*(log100+physical_log_radius_over100-logRa-hb*s_c/2))',
                actual_phase_Z_exact_zero=True,actual_phase_full_period_cover=c.mpf([0,1]),
                real_finite_N_R100_boundary_correction_supplied=False,**dict.fromkeys(fields.previous.OPEN,False))
        self.cache[key]=result;return result

    def contribution(self,label,N=257):
        N=phase.candidate_N(N)
        if N!=self.saved_downstream['candidate_N']:raise ValueError('Same current candidate N required for driver composition')
        op=self.first.owner(label);ref=self.downstream.reference.owner(label);f,c=ref.flow,ref.c
        local={name:[f.scalar(0),f.scalar(0)] for name in RATES};windows={}
        with mp.workdps(c.dps+40):
            for chart in CHARTS:
                total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
                for left,right in zip(PARTITION,PARTITION[1:]):
                    source=self.query(label,chart,left,right,N);rows={};weights={}
                    for name,rate in RATES.items():
                        mass,decay,tail=own_weights(op,chart,left,right,rate)
                        positive_jet(f,[mass],'actual_positive_cell_Duhamel_mass')
                        density=source['original_signed_five_density_C0_Z']
                        pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                            for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):total[name][n]+=row
                        rows[name]=pair;weights[name]=dict(positive_full_mass=mass,incoming_decay=decay,
                            downstream_suffix_decay=tail,own_rate=str(rate),true_physical_measure_once=True)
                    cells.append(dict(source=source,signed_cell_driver_C0_Z=rows,own_rate_weights=weights))
                memory={name:own_weights(op,chart,(0,1),(1,1),rate)[1] for name,rate in RATES.items()}
                for name in RATES:local[name]=[local[name][n]*memory[name]+total[name][n] for n in range(2)]
                windows[chart]=dict(actual_full_source_cells=cells,actual_local_signed_driver_C0_Z=total,
                    retained_actual_incoming_memory=memory)
                print('Current original switch source window',label,chart,flush=True)
            accepted=self.saved_downstream['frames'][label]
            guard_saved_source(self.family,label,N,ref.P0,accepted)
            restore=lambda row:first.endpoint.restore_row(f,row)
            child={name:[restore(row) for row in rows] for name,rows in accepted['actual_R110_Rm_local_driver_C0_Z'].items()}
            suffix={name:restore(row) for name,row in accepted['actual_R110_Rm_incoming_memory'].items()}
            complete={name:[local[name][n]*suffix[name]+child[name][n] for n in range(2)] for name in RATES}
            memory={name:suffix[name]*f.scalar(c.exp(-c.mpf(rate.numerator)/rate.denominator*c.ln(c.mpf(11)/10)))
                    for name,rate in RATES.items()}
        return dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=ref.P0,
            actual_source_windows=windows,actual_R100_R110_local_driver_C0_Z=local,
            actual_R100_Rm_local_driver_C0_Z=complete,retained_R100_Rm_incoming_memory=memory,
            actual_R100_R110_log_radius_length=c.ln(c.mpf(11)/10),
            exact_affine_boundary_formula='deltaH(Rm)=memory(R100,Rm)*deltaH(R100)+all_local_source_drivers(R100,Rm)',
            genuine_finite_N_R100_correction_still_unsupplied=True,actual_finite_N_Rm_incoming_correction_supplied=False,
            accepted_downstream_driver_rehydrated_in_same_live_source_algebra=True,
            source_owned_local_windows_not_real_inlet_or_five_moment_closure=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalSwitchFiniteN(require_checked=False);frames={}
    for label in ('0','.5'):frames[label]=owner.contribution(label)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_bindings=owner.bindings,frames=downstream.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
