"""Current original R110-to-Rsh source, finite-N densities and Duhamel.

The full finite backward kernels enclose the genuine leading source in
the same live algebra as the Rm reference owner. Candidate density drivers
are integrated separately from homogeneous memory. An actual finite-N
R110 boundary correction remains required, not inferred from background
histories or a saved earlier-owner N1024 route.
"""
import ast
import copy
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_positive_quotients_q as parameters
import lei_ren_part1_paper_compliant_current_original_Rm_conditioned_phase_density as phase
import lei_ren_part1_paper_compliant_current_original_Rm_all_u_density_integrals as primitives
import lei_ren_part1_paper_compliant_current_original_Rm_weighted_averaging as bounds
import lei_ren_part1_paper_compliant_long_reshape_profiles as original
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import exp_average

generic=parameters.previous
fields,base,ep=generic.fields,generic.base,generic.ep
HERE,PREFIX,sha=generic.HERE,generic.PREFIX,generic.sha
NAME=PREFIX+'current_original_long_reshape_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_long_reshape_finite_N_check.json'
GATE='current_original_R110_Rsh_source_owned_finite_N_density_and_local_Duhamel_installed'
RATES=bounds.RATES
PARTITION=((0,1),(1,4),(1,2),(3,4),(1,1))
C0,Z=(0,0),(0,1)


def serialized(value):
    if isinstance(value,dict):return {('y%d_Z%d'%key if isinstance(key,tuple) else key):serialized(row) for key,row in value.items()}
    if isinstance(value,(tuple,list)):return [serialized(row) for row in value]
    return fields.serialized(value)


def compile_recovery():
    """Keep original signed recovery; change only explicit source geometry/C."""
    tree=ast.parse(Path(generic.__file__).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if getattr(n,'name',None)=='recover_inputs'))
    fn.name='recover_long_inputs';nodes=[];skipping=False;replacements=[]
    for node in fn.body:
        name=node.targets[0].id if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) else None
        if name=='geometry':
            nodes.append(ast.parse("geometry=packet['geometry']").body[0]);skipping=True
        elif skipping:
            if name=='R':
                assert ast.dump(node.value)==ast.dump(ast.parse('op.Rm_factor*x',mode='eval').body)
                nodes.append(ast.parse('R=op.source_radius').body[0]);skipping=False;replacements.append('explicit physical R=110*exp(y)')
        elif name=='C':
            assert ast.dump(node.value)==ast.dump(ast.parse('add(E,scale(Ey,-2))',mode='eval').body)
            nodes.append(ast.parse('C=op.correlated_C').body[0]);replacements.append('same C=E*(.8+2B*sigma_prime/T) identity')
        else:nodes.append(node)
    assert not skipping and len(replacements)==2
    fn.body=nodes;module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))
    scope=dict(ep=ep,fields=fields)
    exec(compile(module,'<original generic recovery with explicit long radius>','exec'),scope)
    return scope[fn.name],dict(original_recovery_source=Path(generic.__file__).name,
        signed_inertial_pressure_meridional_assignments_unchanged=True,geometry_and_correlated_C_changes=replacements)


RECOVER,RECOVERY_BINDING=compile_recovery()


def fraction(value):
    if type(value) is not tuple or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact rational original reshape phase required')
    q=Fraction(*value)
    if not 0<=q<=1:raise ValueError('Original reshape phase in[0,1] required')
    return q


def bounded_mass(f,width,rate,tail):
    """Finite positive exponential mass, stable under microscopic refinement."""
    if ep(rate*width)[1]<=1:
        return width*exp_average(f.c,-rate*width)
    return (1-f.ordinary_cover(tail))/rate


def finite_kernel(f,B,T,y,kind):
    """Full original finite integral with ordinary axial Bell enclosures."""
    c=f.c;k,m,rmin,rmax={'theta':('1.6',1,'1.55','1.65'),
        'pressure':('.2',2,'.1','.3'),'swirl':('1.2',2,'1.1','1.3')}[kind]
    k,rmin,rmax=c.mpf(k),c.mpf(rmin),c.mpf(rmax)
    if ep(y)[0]<0 or ep(T)[0]<=0:raise ValueError('Original finite positive T and y>=0 required')
    beta=[c.mpf(max(abs(v) for v in ep(row))) for row in B.coefficients]
    ratio=8*m*beta[0]/T
    if ep(ratio)[1]>ep(k-rmin)[0] or ep(ratio)[1]>ep(rmax-k)[0]:
        raise ValueError('Original finite kernel rates must follow from actual B/T')
    yl,yh=(c.mpf(v) for v in ep(y))
    lower_tail=f.factor((0,0,0,0,0),-rmax*yl);upper_tail=f.factor((0,0,0,0,0),-rmin*yh)
    if ep(y)==(0,0):rows=[f.scalar(0)]*6
    else:
        # Tiny tails remain explicit formal sources. Directed ordinary covers
        # are used only to enclose this bounded kernel, never select its value.
        lower=ep(bounded_mass(f,yl,rmax,lower_tail))[0]
        upper=ep(bounded_mass(f,yh,rmin,upper_tail))[1]
        coefficients=[c.mpf([max(mp.mpf(0),lower),max(mp.mpf(0),upper)])];polys=[[c.mpf(1)]]
        for n in range(1,6):
            polynomial=[c.mpf(0)]*(n+1)
            for j in range(1,n+1):
                for power,term in enumerate(polys[n-j]):polynomial[power+1]+=beta[j]*(8*m/T)*j/n*term
            polys.append(polynomial);cap=c.mpf(0)
            for power in range(1,n+1):
                finite=yh**(power+1)/(power+1);complete=c.mpf(math.factorial(power))/rmin**(power+1)
                cap+=polynomial[power]*c.mpf(min(ep(finite)[1],ep(complete)[1]))
            coefficients.append(c.mpf([-ep(cap)[1],ep(cap)[1]]))
        rows=f.jet(IntervalTaylor(c,coefficients))
    return rows,dict(original_full_finite_definition='integral_0^y exp(-k*t+m*B*(sigma(y/T)-sigma((y-t)/T)))dt',
        k=k,m=m,rate_min=rmin,rate_max=rmax,actual_B_rate_ratio=ratio,
        full_finite_lower_tail=lower_tail,full_finite_upper_tail=upper_tail,
        genuine_B_axial_Bell_terms_through5=True,full_source_T_not_shortened=True,
        source_kernel_enclosures_not_selected_values=True)


def source_packet(long,reference,physical_radius,reshape_phase,y):
    """Actual background plus raw first y/axial5 rows in current units."""
    f,c=long.flow,long.c
    if reference.flow is not f or reference.T._mpi_!=long.T._mpi_ or reference.logC._mpi_!=long.logC._mpi_:
        raise ValueError('Same current long/reference flow, frozen T and logC required')
    for n,row in enumerate(reference.P0):
        if row.record()!=f.scalar(long.p0[n]).record():raise ValueError('Canonical current analytic P0 tuple identity required')
    z=reference.z;B=long.B;T=long.T;cutoff=sigma_jets(c,reshape_phase)
    sig=intersection(c,cutoff[0],c.mpf([0,1]));ds=intersection(c,cutoff[1],c.mpf([0,8]))
    logu=B*(1-sig)-original.logarithm(1+z*z)+(y/10-long.logC-f.logs[1]/2)
    ratio=IntervalTaylor(c,[c.mpf(0)]+list(logu.coefficients[1:])).exp()
    E=f.scale(f.jet(ratio),f.factor((0,0,0,0,0),logu[0]));slope=IntervalTaylor.constant(c,'.1',5)-B*(ds/T)
    a=IntervalTaylor.constant(c,'.8',5)+B*(2*ds/T)
    a=IntervalTaylor(c,[intersection(c,a[0],c.mpf(['.72','.88']))]+list(a.coefficients[1:]))
    arows=f.jet(a);correlated_C=f.multiply(arows,E);V=long.V;zero=[f.scalar(0)]*6
    shapes={};kernels={};inherited={}
    for name,kind,m,k in (('theta','theta',1,'1.6'),('theta_z','theta',1,'1.6'),
                         ('pressure','pressure',2,'.2'),('swirl','swirl',2,'1.2')):
        if kind not in kernels:kernels[kind]=finite_kernel(f,B,T,y,kind)
        bell=IntervalTaylor(c,[c.mpf(0)]+[row*m*sig for row in B.coefficients[1:]]).exp()
        decay=f.scale(f.jet(bell),f.factor((0,0,0,0,0),B[0]*(m*sig)-c.mpf(k)*y))
        inherited[name]=f.multiply(long.inlets[name],decay)
        source=kernels[kind][0] if name!='theta_z' else f.multiply(V,kernels[kind][0])
        shapes[name]=f.add(inherited[name],source)
    dmean=f.factor((0,0,0,0,0),-y);V2=f.multiply(V,V)
    shapes['mean']=f.add(V,f.scale(f.add(long.inlets['mean'],f.scale(V,-1)),dmean))
    shapes['axial']=f.add(V2,f.scale(f.add(long.inlets['axial'],f.scale(V2,-1)),dmean))
    E2=f.multiply(E,E);invP2=f.factor((0,-1,0,0,0))
    histories=dict(m=shapes['mean'],h=f.multiply(E,shapes['theta']),k=f.multiply(E,shapes['theta_z']),
        e=f.add(f.scale(shapes['axial'],invP2),f.scale(f.multiply(E2,shapes['swirl']),-c.mpf('.5'))),
        p=f.scale(f.multiply(E2,shapes['pressure']),c.mpf('.5')))
    # Raw m/k/V are in physical axial units; original recovery divides each
    # by Pstar once. Exact original five ODEs supply ordinary first y rows.
    derivatives=dict(m=f.add(V,f.scale(histories['m'],-1)),
        h=f.add(E,f.scale(histories['h'],-c.mpf('1.5'))),
        k=f.add(f.multiply(E,V),f.scale(histories['k'],-c.mpf('1.5'))),
        e=f.add(f.scale(V2,invP2),f.scale(histories['e'],-1),f.scale(E2,-c.mpf('.5'))),
        p=f.scale(E2,c.mpf('.5')))
    raw=dict(histories={name:[row,derivatives[name]] for name,row in histories.items()},
        velocity=dict(theta=[E,f.multiply(E,f.jet(slope))],axial=[V,zero]))
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=raw,original_P0_normalized_axial5=reference.P0,
        geometry=dict(chart='current_original_long_reshape',reshape_phase=reshape_phase,log_R_over110=y,
            original_frozen_T=T,actual_physical_radius=physical_radius,phase_and_radius_Z_independent=True))
    op=SimpleNamespace(flow=f,c=c,reference=reference,zrows=reference.zrows,P0=reference.P0,
        Pstar=f.factor((0,.5,0,0,0)),source_radius=physical_radius,correlated_C=correlated_C)
    recovered=RECOVER(op,packet)
    recovered.pop('source_frame_conditional_on_same_actual_Rm_inlet')
    recovered['source_frame_conditional_on_same_current_R110_background']=True
    return op,packet,recovered,arows,dict(actual_normalized_six_history_shapes=shapes,
        actual_nonzero_inherited_angular_history_rows=inherited,actual_full_finite_kernel_proofs={name:value[1] for name,value in kernels.items()},
        exact_mean_axial_incoming_decay=dmean,actual_E_log_source_axial5=list(logu.coefficients),
        actual_E_log_y_slope_axial5=list(slope.coefficients),actual_correlated_shear_a_axial5=arows,
        exact_V_y_and_V_yZ_zero=True,one_Pstar_conversion_for_raw_m_k_V=True,
        original_raw_first_y_histories_from_exact_five_ODEs=True,exact_common_P0_axial5=reference.P0)


def source_quotients(op,recovered,a,eta_log):
    f=op.flow;n=recovered['actual_generic_source_numerators'];E=n['E']
    proof=parameters.positive_source(f,E[0],'actual_long_E');ap=parameters.positive_source(f,a[0],'actual_long_a')
    Delta=parameters.add(f,a,[f.scalar(-2)]+[f.scalar(0)]*5)
    q=parameters.shear_q_jets(f,a,Delta,eta_log,ap)
    if not q['analytic_axial_jets_available'] or ep(Delta[0].coefficient)[1]>=0:
        raise ValueError('Original long active negative-Delta branch required')
    sectors=recovered['full_signed_inertial_sectors_axial4'];quotients={}
    for name,left,right in (('p1','theta_linear','theta_quadratic'),('p2','axial_linear','axial_quadratic')):
        before=parameters.add(f,sectors[left],parameters.scale(sectors[right],op.Pstar))
        quotients[name]=parameters.quotient(f,before,E,proof)
    p2=parameters.scale(quotients['p2'],op.source_radius)
    roots={name:{C0:rows[0],Z:rows[1]} for name,rows in dict(a=a,t0=[f.scalar(0)]*6,E=E,p2=p2).items()}
    qr={C0:q['q_axial_coefficients'][0],Z:q['q_axial_coefficients'][1]}
    return dict(q=qr[C0],roots=roots),qr,dict(full_original_shear_a_axial5=a,Delta_axial5=Delta,
        exact_b_and_t0_zero=True,original_active_q_proof=q,full_signed_inertial_before_R_axial4=quotients,
        exact_full_p2_axial4=p2,one_actual_positive_radius_factor=op.source_radius,
        source_positive_E=proof,source_positive_a=ap,selected_parameters_not_new_cone_admission=True)


def kernel_weight(f,width,suffix,rate):
    c=f.c;lam=c.mpf(rate.numerator)/rate.denominator
    if ep(width)[0]<=0 or ep(suffix)[0]<0:raise ValueError('Positive original cell width and nonnegative suffix required')
    if not rate:return f.scalar(width),f.scalar(1),f.scalar(1)
    decay=f.factor((0,0,0,0,0),-lam*width);downstream=f.factor((0,0,0,0,0),-lam*suffix)
    mass=(f.scalar(width*exp_average(c,-lam*width)) if ep(lam*width)[1]<=1
          else (f.scalar(1)-decay)*(c.mpf(1)/lam))
    return mass,decay,downstream


class OriginalLongReshapeFiniteN:
    mode='current_original_full_R110_Rsh_source_finite_N_signed_density_C0_Z_and_affine_local_Duhamel'
    def __init__(self,dps=500,require_checked=True):
        self.parameters=parameters.OriginalRmPositiveQuotientsQ(dps);self.generic=self.parameters.upstream
        self.c=self.parameters.c;self.family=self.parameters.family;self.hashes=dict(self.parameters.hashes);self.cache={}
        self.reference_wrapper=self.generic.upstream.upstream.upstream
        self.long_wrapper=self.reference_wrapper.upstream
        for name in (Path(__file__).name,Path(original.__file__).name,Path(primitives.__file__).name,
                     Path(phase.__file__).name,Path(bounds.__file__).name,
                     PREFIX+'current_O3_independent_repair_operator.py'):fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted current long-reshape finite-N receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def inlet_binding(self,label,long):
        packet=self.long_wrapper.saved['post_power_packets'][label][-1]
        geometry=packet['post_power_function_evaluation']['geometry']
        if (packet['source_family']!=self.family or packet['source_frame']!=label
                or geometry['exact_fixed_radius']!=[110,1]):
            raise ValueError('Same current source family, frame and fixed R110 geometry required')
        stored=packet['original_anchored_log_shape_B_axial5']
        if any(fields.previous.read_interval(self.c,row)._mpi_!=value._mpi_
               for row,value in zip(stored,long.B.coefficients)) or len(stored)!=6:
            raise ValueError('Same current R110 anchored B source tuple required')
        if long.T._mpi_!=self.long_wrapper.T._mpi_ or long.logC._mpi_!=self.long_wrapper.logC._mpi_:
            raise ValueError('Same current source frozen T and logC required')
        return dict(source_family=self.family,source_frame=label,exact_fixed_inlet_radius=[110,1],
            current_R110_source_receipt=PREFIX+'current_original_second_switch_R110.json',
            same_actual_anchored_B_axial5=True,same_frozen_T_and_logC=True,
            leading_source_binding_not_finite_N_incoming_correction=True)

    def query(self,label,left,right=None,N=257):
        N=phase.candidate_N(N);lo=fraction(left);hi=lo if right is None else fraction(right)
        if hi<lo:raise ValueError('Ordered original reshape phase cell required')
        key=(label,left,right,N)
        if key in self.cache:return self.cache[key]
        long=self.long_wrapper.owner(label);reference=self.reference_wrapper.owner(label);f,c=long.flow,long.c
        inlet=self.inlet_binding(label,long)
        with mp.workdps(c.dps+40):
            rho=c.mpf([ep(c.mpf(lo.numerator)/lo.denominator)[0],ep(c.mpf(hi.numerator)/hi.denominator)[1]])
            y=long.T*rho;radius=f.factor((0,0,0,0,0),c.ln(110)+y)
            op,packet,recovered,a,source_record=source_packet(long,reference,radius,rho,y)
            source,qr,quotient_record=source_quotients(op,recovered,a,self.parameters.eta_log)
            got=primitives.all_u_primitive_bounds(f,source,qr,self.parameters.dstar_log,c.mpf([0,1]))
            density=phase.densities.density_Z_kernels(source['roots']['E'][C0],source['roots']['E'][Z],
                recovered['common_velocity_V_axial5'][0],recovered['common_velocity_V_axial5'][1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,source_geometry=packet['geometry'],
                exact_common_P0_axial5=reference.P0,original_raw_current_source=packet,original_generic_source=recovered,
                current_R110_source_binding=inlet,
                original_background_function_record=source_record,original_source_quotient_record=quotient_record,
                original_roots=source['roots'],original_q_rows=qr,original_primitive_values=got['values'],
                original_primitive_proof=got['record'],original_signed_five_density_C0_Z=density,
                actual_radius_phase_definition='frac(N*(log(110)+y-logRa-hb*s_c/2))',
                actual_radius_phase_cover=c.mpf([0,1]),actual_radius_phase_Z_exact_zero=True,
                full_period_cover_encloses_actual_phase_not_selected_inverse=True,
                current_original_source_owner_and_background_histories_retained=True,
                genuine_finite_N_R110_boundary_correction_supplied=False,whole_axis_provider_installed=False,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[key]=result;return result

    def contribution(self,label,N=257):
        N=phase.candidate_N(N);long=self.long_wrapper.owner(label);reference=self.reference_wrapper.owner(label)
        f,c=long.flow,long.c;total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
        with mp.workdps(c.dps+40):
            for left,right in zip(PARTITION,PARTITION[1:]):
                source=self.query(label,left,right,N);l,r=fraction(left),fraction(right)
                width=long.T*(c.mpf((r-l).numerator)/(r-l).denominator);suffix=long.T*(c.mpf((1-r).numerator)/(1-r).denominator)
                rows={};weights={}
                for name,rate in RATES.items():
                    mass,decay,downstream=kernel_weight(f,width,suffix,rate)
                    density=source['original_signed_five_density_C0_Z']
                    pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*downstream)
                          for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):total[name][n]+=row
                    rows[name]=pair;weights[name]=dict(exact_positive_mass=mass,cell_incoming_decay=decay,
                        downstream_suffix_decay=downstream,own_rate=str(rate),physical_log_radius_measure_once=True,
                        original_full_finite_cell_not_shortened=True,nonlinear_source_driver_not_zeroed=True)
                cells.append(dict(source=source,original_log_width=width,original_suffix_to_Rsh=suffix,
                    actual_cell_signed_driver_C0_Z_enclosures=rows,original_own_rate_weights=weights))
            memory={name:f.scalar(1) if not rate else f.factor((0,0,0,0,0),-long.T*c.mpf(rate.numerator)/rate.denominator)
                    for name,rate in RATES.items()}
        return dict(source_family=self.family,source_frame=label,candidate_N=N,original_frozen_T=long.T,
            exact_common_P0_axial5=reference.P0,exact_partition=PARTITION,actual_full_window_source_cells=cells,
            actual_local_finite_N_driver_C0_Z=total,actual_nonzero_incoming_memory_decays=memory,
            exact_affine_boundary_formula='deltaH(Rsh)=exp(-lambda*T)*deltaH(R110)+local_signed_driver',
            actual_R110_correction_is_still_unsupplied=True,leading_R110_histories_not_finite_N_boundary=True,
            no_saved_N1024_boundary_transplanted=True,actual_Rm_incoming_correction_supplied=False,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalLongReshapeFiniteN(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Current original long-reshape finite-N source and signed Duhamel',label,flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_generic_recovery_binding=RECOVERY_BINDING,frames=serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
