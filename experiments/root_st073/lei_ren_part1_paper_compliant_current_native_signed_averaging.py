"""Original zero-mean density split and endpoint-retaining C0/Z averaging.

The source functions remain the original loop. Bounds of the leading slow
mixed jets do not define their values. All-N integral covers are rebuilt
from the actual zero inlet; existing tighter direct covers remain valid.
"""
import gzip
import json
import math
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets as current
import lei_ren_part1_paper_compliant_current_generic_loop_function_sources as functions
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets,ep,prior=current.packets,current.ep,current.prior
native,serial,repair=current.native,current.serial,current.repair
LogUpper=repair.LogUpper;ZERO=(0,0);DY=(1,0);DZ=(0,1);DYZ=(1,1)
ORDERS=(ZERO,DY,DZ,DYZ);MIN_N=current.MIN_N;RATES=current.RATES
NAME=PREFIX+'current_native_signed_averaging.json'
RECEIPT=PREFIX+'current_native_signed_averaging_check.json'
VIEWS=PREFIX+'current_native_signed_averaging_views.json.gz'
GATE='current_original_native_zero_mean_C0_Z_Nminus2_source_integral_covers_executed'


class SlowJet:
    """Directed absolute caps for ordinary fixed-phi 1 by 1 slow jets."""
    def __init__(self,c,rows):
        if set(rows)!=set(ORDERS) or any(v.ctx is not c for v in rows.values()):
            raise ValueError('Same-context C0,y,Z,yZ source caps required')
        self.ctx,self.rows=c,rows
    def __getitem__(self,k):return self.rows[k]
    @classmethod
    def constant(cls,c,v):return cls(c,{k:LogUpper.constant(c,v if k==ZERO else 0) for k in ORDERS})
    def __add__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Common slow-jet context required')
        return SlowJet(self.ctx,{k:LogUpper.add(self.ctx,[self[k],other[k]]) for k in ORDERS})
    def __mul__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Common slow-jet context required')
        c=self.ctx;rows={}
        for j,k in ORDERS:
            rows[(j,k)]=LogUpper.add(c,[self[(i,ell)]*other[(j-i,k-ell)]
                for i in range(j+1) for ell in range(k+1)])
        return SlowJet(c,rows)
    def record(self):return {'y%d_Z%d'%k:v.record() for k,v in self.rows.items()}


def read_cap(c,row):
    if row['exact_zero']!=(row['log_absolute_upper'] is None):raise ValueError('Exact-zero cap encoding differs')
    return LogUpper(c,None if row['exact_zero'] else packets.interval(c,row['log_absolute_upper']))


def minimum(*caps):
    if any(v.log is None for v in caps):return LogUpper(caps[0].ctx,None)
    return min(caps,key=lambda v:ep(v.log)[1])


def leading_caps(c,raw,loop,primitive):
    """Normalized E,V,A and loop B/Pstar; no extra width or S conversion."""
    read=lambda rows:SlowJet(c,{k:read_cap(c,rows['y%d_Z%d'%k]) for k in ORDERS})
    E,V=read(raw['E']),read(raw['V']);A=read(loop['slow_phase_held_log_bounds']['A'])
    B=read(loop['slow_phase_held_log_bounds']['B'])
    # The new periodic C0/Z certificates intersect only absolute bounds;
    # their signed source functions are the same original loop functions.
    for key,jet,names in (('A',A,('A','A_Z')),('B_over_Pstar',B,('B_over_Pstar','B_Z_over_Pstar'))):
        for order,name in zip((ZERO,DZ),names):
            jet.rows[order]=minimum(jet[order],read_cap(c,current.magnitude(primitive[name]).record()))
    two=SlowJet.constant(c,2)
    f=dict(m=B,h=E*A,k=E*(V*A+B),e=two*V*B+E*E*A,p=E*E*A)
    # e uses absolute caps here; the signed formula is 2VB-E^2A.
    return dict(E=E,V=V,A=A,B=B,leading=f)


def remainder_caps(jets):
    """Exact exponential second remainder, uniformly for N>=160."""
    E,V,A,B=(jets[k] for k in ('E','V','A','B'));c=E.ctx
    C=lambda value:LogUpper.constant(c,value);add=lambda *v:LogUpper.add(c,v)
    expcap=LogUpper(c,c.mpf(159)/MIN_N)
    # |A|<=159 is an original analytic theorem, not a selected A value.
    R2=expcap*C('.5');R3=expcap*C(c.mpf(1)/6)
    E0,EZ,A0,AZ,B0,BZ,V0,VZ=(E[ZERO],E[DZ],A[ZERO],A[DZ],B[ZERO],B[DZ],V[ZERO],V[DZ])
    F=E0*A0*expcap;FZ=add(EZ*A0*expcap,E0*AZ*expcap)
    RE=E0*A0*A0*R2
    REZ=add(EZ*A0*A0*R2,E0*A0*AZ*R2*C(2),E0*A0*A0*AZ*R3*C(c.mpf(1)/MIN_N))
    zero=C(0);half=C('.5')
    values=dict(m=zero,h=RE,k=add(V0*RE,F*B0),e=add(E0*RE,B0*B0,F*F*half),p=add(E0*RE,F*F*half))
    Z=dict(m=zero,h=REZ,k=add(VZ*RE,V0*REZ,FZ*B0,F*BZ),
        e=add(EZ*RE,E0*REZ,B0*BZ*C(2),F*FZ),p=add(EZ*RE,E0*REZ,F*FZ))
    return dict(values=values,Z_derivatives=Z,record=dict(
        delta_E='E*A/N+R_E/N^2; R_E=E*A^2*R2(A/N)',delta_V='B_over_Pstar/N',
        R2='integral_0^1 (1-s)*exp(s*A/N) ds',R3='integral_0^1 s*(1-s)*exp(s*A/N) ds',
        R_E_Z='E_Z*A^2*R2+2*E*A*A_Z*R2+E*A^2*A_Z*R3/N',
        F_N='E*A*integral_0^1 exp(s*A/N) ds',
        remainder_density=dict(m='0',h='R_E',k='V*R_E+F_N*B',e='-E*R_E+B^2-F_N^2/2',p='E*R_E+F_N^2/2'),
        nonlinear_F_N_not_assumed_zero_mean=True,N_dependent_remainder_coefficient_functions=True,
        uniform_N_lower=MIN_N,original_A_cap=159,all_signed_cross_and_quadratic_terms_retained=True))


def exact_theorem():
    E,V,A,B,R,F,N=sy.symbols('E V A B R_E F_N N',nonzero=True)
    # F_N=E*A+R_E/N follows from the exact exponential remainder.
    dE=E*A/N+R/N**2;dV=B/N
    delta=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    f=dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E*E*A,p=E*E*A)
    rem=dict(m=0,h=R,k=V*R+F*B,e=-E*R+B*B-F*F/2,p=E*R+F*F/2)
    for key in RATES:
        if sy.expand(delta[key]-f[key]/N-sy.sympify(rem[key]).subs(F,E*A+R/N)/N**2)!=0:
            raise ArithmeticError('Exact original zero-mean split failed: '+key)
    y,Z,phi=sy.symbols('y Z phi');lam=sy.Symbol('lambda',nonnegative=True)
    G=sy.Function('G')(y,Z,phi);K=sy.Function('K')(y)
    # D_y G(y,Z,Ny)=G_y_slow+N*G_phi; K_y=lambda*K.
    for function in (G,sy.diff(G,Z)):
        total=lam*K*function+K*(sy.diff(function,y)+N*sy.diff(function,phi))
        if sy.expand(K*sy.diff(function,phi)/N-(total-K*(sy.diff(function,y)+lam*function))/N**2)!=0:
            raise ArithmeticError('Original C0/Z IBP sign identity failed')
    primitive_theorem=original.exact_theorem()
    return dict(passed=True,original_primitive_parity_theorem=primitive_theorem,
        signed_leading_densities=dict(m='B',h='E*A',k='E*(V*A+B)',e='2*V*B-E^2*A',p='E^2*A'),
        source_parity='A(1-phi)=-A(phi), B(1-phi)=-B(phi), with original t0=-b/a',
        uniform_fractional_phi_mean_of_leading_and_fixed_phase_slow_jets_exact_zero=True,
        primitive='G_j(y,Z,phi)=integral_0^phi f1_j(y,Z,s)ds; G_j(0)=G_j(1)=0',
        primitive_derivative_caps='|G_alpha|<=sup_phi|D_slow^alpha f1|/2 for alpha=0,y,Z,yZ',
        exact_five_signed_linear_second_remainder_identities=True,
        source_IBP='integral_a^b K*f1(y,Z,frac(N*y+phi0))/N dy = ([K*G]_a^b-integral_a^b K*(G_y_slow+lambda*G)dy)/N^2',
        Z_IBP_uses_G_yZ_slow=True,global_phase_y_Z_exact_zero=True,
        each_cell_endpoint_terms_retained=True,seam_cancellation_not_claimed=True,
        actual_zero_inlet_required_for_whole_route_Nminus2=True,incoming_pressure_rate0_memory_retained=True,
        slow_y_not_total_y=True,higher_jets_or_actual_point_inverse_not_installed=True)


def mixed_source_graph(view):
    """Exact source definitions of mixed inverse/A/B and leading 1x1 jets.

    Formal partial derivative instructions differentiate the original free
    angle graph. They are function operators, never derivatives of caps.
    """
    g=functions.FunctionGraph(view['original_signed_input_graph'])
    g.nodes=list(view['function_graph_nodes']);g.keys={json.dumps(v,sort_keys=True):i for i,v in enumerate(g.nodes)}
    one,two=g.one,g.integer(2);pi=g.node('mathematical_pi');twopi=g.mul(two,pi)
    half=g.div(one,two,'exact_positive_integer2');inv2pi=g.div(one,twopi,'exact_positive_2pi')
    roots=view['roots'];inverse=roots['inverse_angle'];Phi=roots['Phi'];lam=roots['lambda_at_free_angle'];T1=roots['T1_at_free_angle']
    t=roots['direction_at_free_angle'];at=lambda body:g.at_inverse(body,inverse)
    def derivative(body,variables):return g.node('formal_original_function_partial_derivative',body=body,
        variables=variables,other_free_variables_fixed=True,global_eta_dstar_N_and_Pstar_fixed=True,
        differentiates_function_graph_not_saved_cover=True)
    ilam=g.div(one,at(lam),'same_original_monotone_inverse_lambda_positive')
    psiy=g.neg(g.mul(at(derivative(Phi,['y'])),ilam));psiZ=g.neg(g.mul(at(derivative(Phi,['Z'])),ilam))
    psiyZ=g.neg(g.mul(ilam,g.add(at(derivative(Phi,['y','Z'])),
        g.mul(at(derivative(lam,['y'])),psiZ),g.mul(at(derivative(lam,['Z'])),psiy),
        g.mul(at(derivative(Phi,['psi','psi'])),psiy,psiZ))))
    phi=g.node('function_variable',name='phi',domain='periodic_modulo1');chi=g.sub(phi,g.mul(inverse,inv2pi))
    a,b,E=(g.root(key) for key in ('a','b','E'))
    ay,aZ,ayZ=(g.root('a',k) for k in (DY,DZ,DYZ));byZ=g.root('b',DYZ)
    T=at(T1);Ty=g.add(at(derivative(T1,['y'])),g.mul(at(t),psiy));TZ=g.add(at(derivative(T1,['Z'])),g.mul(at(t),psiZ))
    TyZ=g.add(at(derivative(T1,['y','Z'])),g.mul(at(derivative(t,['y'])),psiZ),
        g.mul(at(derivative(t,['Z'])),psiy),g.mul(at(derivative(t,['psi'])),psiy,psiZ),g.mul(at(t),psiyZ))
    M=g.sub(g.neg(g.mul(a,T,inv2pi)),g.mul(b,phi))
    My=g.sub(g.neg(g.mul(g.add(g.mul(ay,T),g.mul(a,Ty)),inv2pi)),g.mul(g.root('b',DY),phi))
    MZ=g.sub(g.neg(g.mul(g.add(g.mul(aZ,T),g.mul(a,TZ)),inv2pi)),g.mul(g.root('b',DZ),phi))
    MyZ=g.sub(g.neg(g.mul(g.add(g.mul(ayZ,T),g.mul(ay,TZ),g.mul(aZ,Ty),g.mul(a,TyZ)),inv2pi)),g.mul(byZ,phi))
    AYZ=g.mul(half,g.sub(g.mul(ayZ,chi),g.mul(g.add(g.mul(ay,psiZ),g.mul(aZ,psiy),g.mul(a,psiyZ)),inv2pi)))
    BYZ=g.mul(half,g.add(g.mul(g.root('E',DYZ),M),g.mul(g.root('E',DY),MZ),
        g.mul(g.root('E',DZ),My),g.mul(E,MyZ)))
    delta=roots['Delta'];eta=roots['eta']
    AYZ=g.flat(AYZ,delta,eta);BYZ=g.flat(BYZ,delta,eta)
    V={}
    for j,k in ORDERS:
        V[(j,k)]=g.node('original_normalized_packet_ordinary_velocity_derivative',component='axial',
            ordinary_y_order=j,ordinary_Z_order=k,packet_source=packets.encode(view['original_signed_input_graph']['source_provenance']),
            source_definition="ordinary_axial_coefficient(packet.velocity['axial'][j],k)",
            source_unit='Uz/Pstar',no_additional_Pstar_width_or_radial_half_shift=True)
    ej={k:g.root('E',k) for k in ORDERS}
    aj={ZERO:roots['A'],DY:roots['A_y_slow'],DZ:roots['A_Z_slow'],DYZ:AYZ}
    bj={ZERO:roots['B_over_Pstar'],DY:roots['B_y_slow'],DZ:roots['B_Z_slow'],DYZ:BYZ}
    def add(*js):return {k:g.add(*(j[k] for j in js)) for k in ORDERS}
    def neg(j):return {k:g.neg(v) for k,v in j.items()}
    def mul(l,r):return {(j,k):g.add(*(g.mul(l[(i,ell)],r[(j-i,k-ell)])
        for i in range(j+1) for ell in range(k+1))) for j,k in ORDERS}
    twice=lambda j:{k:g.mul(two,v) for k,v in j.items()}
    EA=mul(ej,aj);EEA=mul(ej,EA)
    leading=dict(m=bj,h=EA,k=mul(ej,add(mul(V,aj),bj)),e=add(twice(mul(V,bj)),neg(EEA)),p=EEA)
    return dict(chart=view['chart'],source_family=view['source_family'],function_graph_nodes=g.nodes,
        original_signed_input_graph=view['original_signed_input_graph'],
        roots=dict(original_A=roots['A'],original_B_over_Pstar=roots['B_over_Pstar'],inverse_yZ_slow=psiyZ,
            A_yZ_slow=AYZ,B_over_Pstar_yZ_slow=BYZ),
        leading_signed_density_slow_roots={key:{'y%d_Z%d'%k:v for k,v in rows.items()} for key,rows in leading.items()},
        common_phase_binding=view['common_phase_binding'],partial_derivative_nodes_have_function_semantics=True,
        actual_point_inverse_or_mixed_jet_evaluator_installed=False,
        derivative_caps_not_function_values=True,high_mixed_derivatives_not_installed=True)


def symmetric_cap(cap,coords):
    if cap.log is None:return coords.scalar(0)
    return prior.ScaledEnclosure(prior.FormalScale(coords.bases,offset=cap.log),coords.ctx.mpf((-1,1)),coords.ledger)


def ibp_contribution(jets,remainder,factors,key,Z=False):
    c=jets.ctx;half=LogUpper.constant(c,'.5');order=DZ if Z else ZERO;slow=DYZ if Z else DY
    G=jets[order]*half;Gy=jets[slow]*half
    decay=current.magnitude(factors['decay']);mass=current.magnitude(factors['mass'])
    endpoints=LogUpper.add(c,[LogUpper.constant(c,1),decay])*G
    rate=RATES[key];rate=c.mpf(rate.numerator)/rate.denominator
    integral=mass*LogUpper.add(c,[Gy,G*LogUpper.constant(c,rate)])
    rest=mass*remainder
    return LogUpper.add(c,[endpoints,integral,rest]),dict(endpoint_cap=endpoints.record(),
        slow_and_kernel_derivative_integral_cap=integral.record(),quadratic_remainder_integral_cap=rest.record(),
        primitive_cap=G.record(),primitive_slow_cap=Gy.record(),all_each_cell_endpoint_terms_kept=True,
        true_log_radius_mass_used_once=True)


class NativeSignedAveraging:
    def __init__(self,owner):
        if type(owner) is not current.NativeRcParameterTargets:raise ValueError('Same original Rc all-N target owner required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.coordinates=owner.coordinates;self.service=owner.service
        checked=json.loads((HERE/current.RECEIPT).read_bytes())
        if not checked['all_passed'] or not checked[current.GATE] or checked['source_family']!=self.family:
            raise ValueError('Checked same original all-N native target family required')
        self.service.bind_hashes(checked['input_hashes']);self.service.bind_hashes({current.RECEIPT:sha(current.RECEIPT),current.NAME:sha(current.NAME)})
        for module in (functions,functions.signed.current):
            checked=json.loads((HERE/module.RECEIPT).read_bytes());data=json.loads((HERE/module.NAME).read_bytes())
            if not checked['all_passed'] or not checked[module.GATE] or data['source_family']!=self.family:
                raise ValueError('Same original checked loop graph/slow jet caps required')
            self.service.bind_hashes(checked['input_hashes']);self.service.bind_hashes({module.NAME:sha(module.NAME),module.RECEIPT:sha(module.RECEIPT)})
        self.loop=data['current_actual_loop_jet_log_bounds_by_chart']
        self.raw=json.loads((HERE/(PREFIX+'current_generic_shear_source_bounds.json')).read_bytes())['current_original_source_log_bound_charts']
        self.raw.update(json.loads((HERE/(PREFIX+'current_generic_shear_O3_sources.json')).read_bytes())['original_O3_quotient_log_norms'])
        self.graph_views=json.loads(gzip.decompress((HERE/functions.VIEWS).read_bytes()))
        self.graphs={chart:mixed_source_graph(v) for chart,v in self.graph_views.items()}
        if set(self.graphs)!=set(self.loop):raise ValueError('Same17 original graph/slow cap charts required')
        for chart,row in self.loop.items():
            if row['source_family']!=self.family or row['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar':
                raise ValueError('Same normalized original loop B units required')
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name),functions.VIEWS:sha(functions.VIEWS),
            PREFIX+'current_generic_shear_loop.py':sha(PREFIX+'current_generic_shear_loop.py')})
        self.theorem=exact_theorem()

    @native.inlet.source_precision
    def route(self,Z=(-1,1)):
        c=self.ctx;coords=self.coordinates;direct=self.owner.route(Z)
        zero=lambda:{key:coords.scalar(0) for key in RATES}
        running,runningZ=zero(),zero();cells=[]
        for cell in direct['cells']:
            label=cell['record']['label'];chart=cell['record']['chart'];primitive=cell['primitives'];details={}
            flat=primitive is None or all(v.zero for v in primitive['values'].values())
            if flat:values,jets=zero(),zero();source=dict(original_flat_support_exact_zero=True)
            else:
                bounded=leading_caps(c,self.raw[chart]['ordinary_mixed_source_log_norms'],self.loop[chart],primitive['values'])
                remainder=remainder_caps(bounded);values={};jets={}
                for key in RATES:
                    cap,proof=ibp_contribution(bounded['leading'][key],remainder['values'][key],cell['factors'][key],key)
                    capZ,proofZ=ibp_contribution(bounded['leading'][key],remainder['Z_derivatives'][key],cell['factors'][key],key,True)
                    values[key]=symmetric_cap(cap,coords);jets[key]=symmetric_cap(capZ,coords)
                    details[key]=dict(C0=proof,Z=proofZ)
                source=dict(original_flat_support_exact_zero=False,
                    native_original_source_provenance=cell['record']['original_signed_source_C1']['source_provenance'],
                    source_graph_chart=chart,source_graph_views=VIEWS,
                    normalized_leading_density_slow_caps={key:row.record() for key,row in bounded['leading'].items()},
                    normalized_remainder_C0_caps={key:row.record() for key,row in remainder['values'].items()},
                    normalized_remainder_Z_caps={key:row.record() for key,row in remainder['Z_derivatives'].items()},
                    exact_remainder_contract=remainder['record'],same_original_chart_wide_slow_caps_cover_native_cell=True,
                    old_coarse_y_yZ_caps_retained_and_not_promoted_to_tight_new_bounds=True,
                    same_original_parameter_uniform_C0_Z_bounds_intersected=True)
            incoming,incomingZ=running,runningZ
            running={key:cell['factors'][key]['decay']*incoming[key]+values[key] for key in RATES}
            runningZ={key:cell['factors'][key]['decay']*incomingZ[key]+jets[key] for key in RATES}
            record=dict(label=label,chart=chart,original_geometry=cell['geometry']['record'],original_signed_source=source,
                per_density_IBP_bounds=details,
                true_width_source_Nminus2_C0_covers={key:v.record() for key,v in values.items()},
                true_width_source_Nminus2_Z_covers={key:v.record() for key,v in jets.items()},
                inherited_Nminus2_C0_covers={key:v.record() for key,v in incoming.items()},
                inherited_Nminus2_Z_covers={key:v.record() for key,v in incomingZ.items()},
                right_Nminus2_C0_covers={key:v.record() for key,v in running.items()},
                right_Nminus2_Z_covers={key:v.record() for key,v in runningZ.items()},
                each_cell_endpoint_terms_retained=True,incoming_memory_not_reset=True)
            cells.append(dict(record=record,values=values,Z_derivatives=jets,incoming=incoming,incoming_Z=incomingZ,
                cumulative=running,cumulative_Z=runningZ,factors=cell['factors'],flat=flat))
        values={key:{-1:coords.scalar(0),-2:value} for key,value in running.items()}
        jets={key:{-1:coords.scalar(0),-2:value} for key,value in runningZ.items()}
        targets=current.target_rows(values,jets,direct['amplitude'],direct['amplitude_Z'],direct['logA'],direct['mu'],direct['logmu'])
        caps=current.target_caps(targets)
        W=self.owner.repair_manifest['fresh_exact_weight_definitions_and_enclosures'];matrix=self.owner.repair_manifest['fresh_divided_linear_inverse_and_enclosures']
        conditions=repair.contraction_log_conditions(c,caps,matrix,W,direct['logmu'],c.ln(MIN_N))
        # These two certificates bound the same source functions. Keep both:
        # the N^-2 estimate need not be sharper at N=160 or N=1024.
        comparison={}
        for key in repair.ROWS:
            new=read_cap(c,caps['transformed_N_scaled_target_C1_caps'][key])
            old=read_cap(c,direct['caps']['transformed_N_scaled_target_C1_caps'][key])
            coefficient=LogUpper.add(c,[current.magnitude(targets['values'][key][-2]),current.magnitude(targets['Z_derivatives'][key][-2])])
            threshold=None if coefficient.log is None or old.log is None else c.mpf(max(ep(c.ln(MIN_N))[1],
                ep(coefficient.log-c.mpf(ep(old.log)[1]))[1]))
            comparison[key]=dict(averaging_N_scaled_C1_cap_at_floor=new.record(),prior_direct_uniform_N_scaled_C1_cap=old.record(),
                averaging_certificate_sharper_at_floor=(new.log is None or (old.log is not None and ep(new.log)[1]<ep(old.log)[0])),
                sufficient_log_N_to_beat_prior_uniform_certificate=threshold,
                actual_N_scaled_target_C1_bound='min(prior_uniform_cap, averaging_coefficient/N)',
                comparison_is_between_upper_bound_certificates_not_true_errors=True)
        record=dict(source_family=self.family,Z_box=c.mpf(Z),uniform_all_integer_N_lower=MIN_N,
            cells=[row['record'] for row in cells],original_true_cell_count=24,original_chart_count=17,
            same_original_direct_target_record=direct['record'],
            actual_Rc_C0_Nminus2_covers={key:v.record() for key,v in running.items()},
            actual_Rc_Z_Nminus2_covers={key:v.record() for key,v in runningZ.items()},
            actual_Rc_target_C0_orders=current.records(targets['values']),actual_Rc_target_Z_orders=current.records(targets['Z_derivatives']),
            actual_uniform_N_scaled_repair_C1_caps=caps,conditional_repair_log_conditions=conditions,
            same_source_direct_and_averaged_bound_comparison=comparison,
            original_zero_inlet_and_initial_flat_collar_retained=True,
            every_source_increment_rebuilt_not_rescaled_from_old_N1024=True,
            whole_route_Nminus2_cover_not_N_independent_coefficient_identity=True,
            cell_endpoint_cancellation_or_selected_phase_not_assumed=True,
            same_original_P0_and_P0_Z_unchanged=True,
            tightened_first_bridge_yZ_caps_or_actual_controls_installed=False,
            one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,values=values,Z_derivatives=jets,targets=targets,direct=direct)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        target_owner=current.NativeRcParameterTargets(current.preceding.NativeRcC1Histories(
            current.preceding.preceding.NativeO2C1Histories(current.preceding.preceding.make_middle_owner(bridge))))
        owner=NativeSignedAveraging(target_owner);regions={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            regions[name]=owner.route(Z)['record'];print('Original signed source averaging:',name,flush=True)
    (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(packets.encode(owner.graphs),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    owner.service.bind_hashes({VIEWS:sha(VIEWS)})
    result=dict(source_family=owner.family,**{GATE:True},native_Z_query_count=2,all_integer_N_lower=MIN_N,
        actual_original_native_signed_averaging_records=regions,exact_source_averaging_theorem=owner.theorem,
        original_mixed_source_graph_views=VIEWS,source_graph_chart_count=17,
        native_source_generated_integrals_have_C0_Z_Nminus2_covers=True,
        mixed_source_roots_are_formal_exact_function_definitions_not_point_evaluators=True,
        actual_control_functions_or_terminal_closure_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='All24 original cells/17 charts and two Z domains: exact zero-mean linear five-density split, source-linked fixed-phi yZ roots, endpoint-retaining IBP and quadratic remainder give uniform C0/Z N^-2 integral covers from the actual zero inlet. Old coarse slow y/yZ caps remain; direct certificates are retained. No actual controls, terminal identities, common N/global cone/higher jets/energy/recursion/full NS.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
