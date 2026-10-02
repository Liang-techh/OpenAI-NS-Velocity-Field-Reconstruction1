"""Connect actual Rsh/Rm histories to the original implicit five-bump patch.

Axial5 derivatives belong to one smooth implicit coefficient family. Exact
terminal closure is an identity of that map, not a replacement of defects
by fitted zero data. Radial mixed4/full-field/stress lift remain separate.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_reference_restore_profiles import (
    CompliantReferenceRestoreProfiles,IntervalTaylor,scale_history)
from lei_ren_part1_paper_compliant_five_moment_repair import SharedFiveMomentRepair,pack
from lei_ren_part1_paper_interval_five_bump_inverse import (
    certify,nonlinear_jacobian,identity_error,weighted_norm,dot,mag)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    square,derivative,coefficient_lists)
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def nonlinear_series(c,W,h,invAm2):
    zero=IntervalTaylor.constant(c,0,h[0].order)
    return [zero,h[0]*h[2]*W['fg'][0]+h[1]*h[4]*W['fg'][1],zero,
        (square(h[0])*W['gg'][0]+square(h[1])*W['gg'][1])*invAm2
        -sum((square(h[i+2])*W['ff'][i] for i in range(3)),zero)/2,
        sum((square(h[i+2])*W['ff_over_x'][i] for i in range(3)),zero)/2]


def implicit_axial_jets(c,W,d,invAm2,scales,tightening=12):
    """C0/C1 existence plus each higher coefficient via the SAME Jacobian.

    For degree n, insert h_n=0 into the quadratic Taylor map. Its n-th
    coefficient contains every known lower-order convolution and explicit
    invAm2 derivative. The missing term is J(h0)*h_n, so a rigorously
    bounded inverse recovers it without independent point fitting.
    """
    if any(row.order!=5 for row in d) or invAm2.order!=5:
        raise ValueError('Actual defect and amplitude jets through5 required')
    initial=certify(c,W,[row.truncate(1) for row in d],invAm2.truncate(1),
                    scales=scales,tightening=tightening)
    if not initial['certified']:raise ArithmeticError('Actual defect C1 inverse failed')
    L=W['L'];mid=mp.matrix([[(endpoints(v)[0]+endpoints(v)[1])/2 for v in row] for row in L])
    ri=mp.inverse(mid);R=[[c.mpf(ri[i,j]) for j in range(5)] for i in range(5)]
    h0=[row[0] for row in initial['controls']]
    Jq=nonlinear_jacobian(c,W,h0,invAm2[0])
    J=[[L[i][j]+Jq[i][j] for j in range(5)] for i in range(5)]
    K=identity_error(c,R,J);q=weighted_norm(c,K,scales)
    if endpoints(q)[1]>=1:raise ArithmeticError('Actual higher-jet Jacobian inverse failed')
    coefficients=[list(row.coefficients) for row in initial['controls']];proofs=[]
    for n in range(2,6):
        known=[IntervalTaylor(c,row+[c.mpf(0)]) for row in coefficients]
        Q=nonlinear_series(c,W,known,invAm2.truncate(n))
        rhs=[d[i][n]+Q[i][n] for i in range(5)]
        b=[-dot(row,rhs) for row in R]
        bnorm=max((mag(c,b[i])/scales[i] for i in range(5)),key=lambda v:endpoints(v)[1])
        radius=bnorm/(1-q)
        hn=[c.mpf([-endpoints(radius*si)[1],endpoints(radius*si)[1]]) for si in scales]
        for _ in range(tightening):
            hn=[intersection(c,a,b[i]+dot(K[i],hn)) for i,a in enumerate(hn)]
        residual=[dot(J[i],hn)+rhs[i] for i in range(5)]
        if any(not endpoints(value)[0]<=0<=endpoints(value)[1] for value in residual):
            raise ArithmeticError('Implicit higher-jet diagnostic excludes zero')
        proofs.append(dict(order=n,known_nonlinear_coefficient=[row[n] for row in Q],
            source_rhs=rhs,preconditioned_rhs=b,inverse_contraction=q,
            weighted_radius=radius,linear_equation_residual=residual))
        for i in range(5):coefficients[i].append(hn[i])
    controls=[IntervalTaylor(c,row) for row in coefficients]
    Q=nonlinear_series(c,W,controls,invAm2)
    map_jets=[sum((controls[j]*L[i][j] for j in range(5)),controls[0]*0)+Q[i]+d[i] for i in range(5)]
    if any(not endpoints(v)[0]<=0<=endpoints(v)[1] for row in map_jets for v in row.coefficients):
        raise ArithmeticError('Actual implicit map Taylor diagnostic excludes zero')
    return dict(initial_C1_inverse=initial,controls=controls,preconditioner=R,
        actual_Jacobian=J,preconditioned_error=K,higher_inverse_contraction=q,
        higher_order_proofs=proofs,implicit_map_axial5_residual_enclosures=map_jets,
        residual_zero_containment_only_diagnostic=True,
        source_function_closure_proof='unique uniform implicit solution of Lh+Q(h,Am^-2)+d=0; its Taylor coefficients use the same invertible Jacobian',
        smooth_implicit_coefficient_family_through_axial5=True)


class CompliantActualMomentPatch:
    def __init__(self,require_checked=True):
        self.reference=CompliantReferenceRestoreProfiles();self.core=self.reference.core
        self.ctx=c=self.reference.ctx;self.family=self.reference.family;self.source=self.reference.source
        self.repair=SharedFiveMomentRepair();self.hashes=dict(self.reference.hashes);self.cache={}
        for part in ('reference_restore_profiles_check','five_moment_repair_check'):
            name=PREFIX+part+'.json';raw=json.loads((HERE/name).read_bytes())
            for path,digest in raw['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Actual patch dependency changed: '+path)
            accepted=(raw.get('all_passed',False) if part=='reference_restore_profiles_check' else
                raw.get('actual_implicit_functional_five_moment_closure_independently_checked',False)
                and raw.get('corrected_partial_moments_and_same_pressure_independently_checked',False)
                and raw.get('residual_zero_containment_not_used_as_closure_proof',False))
            if not accepted:raise ValueError('Accepted actual patch prerequisites required: '+part)
            if part=='reference_restore_profiles_check' and not (
                raw['actual_source_binding']['exact_actual_core_bridge_switch_Rsh_axial_source_chain']
                and raw['actual_source_binding']['same_original_cutoff_complement_identity']
                and raw['original_axial_restoration_reaches_exact_4Z']):
                raise ValueError('Actual axial source and identical restoration cutoff must be bound')
            if raw.get('actual_five_defect_family_sha256',self.family)!=self.family:
                raise ValueError('Actual patch family mismatch')
            self.hashes.update(raw['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        if self.repair.admit['actual_five_defect_family_sha256']!=self.family or self.repair.datum.source_sha!=self.source:
            raise ValueError('Actual histories and original five-bump source differ')
        self.hashes.update(self.repair.hashes)
        convert=lambda v:c.mpf(endpoints(v))
        self.W={name:([[convert(v) for v in row] for row in values] if name=='L'
                      else [convert(v) for v in values]) for name,values in self.repair.W.items()}
        self.scales=[convert(v) for v in self.repair.scales]
        self.oldtail=read_interval(c,self.repair.admit['actual_tail_C1_positive_cap'])
        self.oldrho=read_interval(c,self.repair.admit['actual_v1_minus_4Z_minus_j_C1_upper'])
        self.e=read_interval(c,self.repair.admit['complete_actual_functional_e_upper'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        if require_checked:
            name=PREFIX+'actual_moment_patch_check.json';check=json.loads((HERE/name).read_bytes())
            for path,digest in check['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Actual patch acceptance changed: '+path)
            if not (check['all_passed'] and check['actual_five_defect_family_sha256']==self.family
                    and check['source_transport_and_implicit_structure']['passed']
                    and check['original_restoration_integrals_directly_source_bound']
                    and check['independent_coefficient_fixture']['passed']):
                raise ValueError('Independent actual source transport and implicit axial5 checks required')
            self.hashes.update(check['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()

    def actual_data(self,Z):
        c=self.ctx;Z=c.mpf(Z);inp=self.reference.inputs(Z);z=inp['z'];E=inp['E']
        initial=inp['source_centered_Rsh'];gap=self.reference.loggap-8;proofs=[]
        decay=lambda row,rate:scale_history(c,row,rate,gap,self.reference.logcap,proofs)
        invAm2=square(1+square(z))*c.exp(c.mpf('1.2')-2*self.core.logP)
        # These are exact source differences, including the negative E/E^2
        # terms from the actual incoming history. No tail is set to zero.
        tails=[decay(initial['mean_error']-E,1)*c.exp(-2),
            decay(initial['mixed_error']-E/c.mpf('1.6'),'1.6')*c.exp(c.mpf('-3.2')),
            decay(initial['angular_error'],'1.6')*c.exp(c.mpf('-3.2')),
            (decay(initial['axial_square']-square(E),1)*invAm2)*c.exp(-2)
                -decay(initial['swirl_error'],'1.2')*(c.exp(c.mpf('-2.4'))/2),
            decay(initial['pressure_error'],'.2')*(c.exp(c.mpf('-.4'))/2)]
        admitted=self.reference.defect_admission
        kernels={name:read_interval(c,row) for name,row in admitted['exact_signed_kernel_enclosures'].items()}
        zero=IntervalTaylor.constant(c,0,5)
        dominant=[E*kernels['d1'],E*kernels['d2'],zero,(square(E)*invAm2)*kernels['d4'],zero]
        exact=[a+b for a,b in zip(dominant,tails)]
        transported=self.reference.defects(Z)['actual_normalized_five_defect_axial5_coefficients']
        # Both formulas represent the identical actual source functions.
        # The independent checker proves their transport algebra, so their
        # enclosure intersection sharpens those functions, not chosen data.
        d=[IntervalTaylor(c,[intersection(c,a,b) for a,b in zip(row.coefficients,other)])
           for row,other in zip(exact,transported)]
        tail_C1=[mag(c,row[0])+mag(c,row[1]) for row in tails]
        if any(endpoints(v)[1]>endpoints(self.oldtail)[0] for v in tail_C1):
            raise ArithmeticError('Actual transported tail not in original C1 admission')
        if endpoints(inp['rho'])[1]>endpoints(self.oldrho)[1]:
            raise ArithmeticError('Actual axial C2 theorem exceeds admitted C1 source class')
        return dict(Z=Z,E=E,rho=inp['rho'],source_centered_Rsh=initial,
            source_log_gap=gap,actual_dominant_kernels=kernels,actual_dominant_functions=dominant,
            actual_tail_functions=tails,actual_tail_C1_norm_upper=tail_C1,
            admitted_tail_C1_cap=self.oldtail,actual_defects=d,invAm2=invAm2,
            retained_source_decay_proofs=proofs,
            exact_tail_definitions=['e^-2*e^-gap*(em_Rsh-E)',
                'e^-3.2*e^-1.6gap*(ek_Rsh-E/1.6)', 'e^-3.2*e^-1.6gap*eh_Rsh',
                'e^-2*Am^-2*e^-gap*(aa_Rsh-E^2)-.5*e^-2.4*e^-1.2gap*eb_Rsh',
                '.5*e^-.4*e^-.2gap*ep_Rsh'],
            correlated_actual_C1_source_admission_retained=True,
            arbitrary_expanded_Taylor_box_not_assigned_Banach_norm=True)

    def coefficients(self,Z):
        Z=self.ctx.mpf(Z);key=Z._mpi_
        if key not in self.cache:
            data=self.actual_data(Z)
            inverse=implicit_axial_jets(self.ctx,self.W,data['actual_defects'],data['invAm2'],self.scales)
            self.cache[key]=(inverse,data)
        return self.cache[key]

    def evaluate(self,x,Z,partial_cells=64):
        """Actual patch velocities/moments through Z5, mean-recovered Ur Z4."""
        c=self.ctx;x=c.mpf(x);Z=c.mpf(Z)
        if endpoints(x)[0]<1 or endpoints(x)[1]>endpoints(c.exp(1))[1]:
            raise ValueError('Original patch 1<=R/Rm<=e required')
        inverse,data=self.coefficients(Z);h=inverse['controls'];d=data['actual_defects']
        z=IntervalTaylor.variable(c,Z,5);am=(1+square(z)).reciprocal()*c.exp(c.mpf('-.6'))
        zero=IntervalTaylor.constant(c,0,5);f=zero;fx=zero;g=zero;gx=zero
        centers=(c.mpf(5)/4,c.mpf(3)/2,c.mpf(7)/4)
        values=[self.repair.beta(c.mpf(endpoints(x-cc))) for cc in centers]
        for k in range(3):f+=h[k+2]*c.mpf(endpoints(values[k][0]));fx+=h[k+2]*c.mpf(endpoints(values[k][1]))
        for i,k in enumerate((0,2)):g+=h[i]*c.mpf(endpoints(values[k][0]));gx+=h[i]*c.mpf(endpoints(values[k][1]))
        H=f+x**c.mpf('.1');u=am*H;uy=am*(fx*x+x**c.mpf('.1')/10)
        V=z*4+g;Vy=gx*x;cache={}
        def w(k,p,m=1):
            key=(k,p,m)
            if key not in cache:cache[key]=c.mpf(endpoints(self.repair.partial_weight(c.mpf(endpoints(x)),('5/4','3/2','7/4')[k],p,m,partial_cells)))
            return cache[key]
        change=[sum((h[i]*w(k,'0') for i,k in enumerate((0,2))),zero),
            sum((h[i]*w(k,'.6')+h[i]*h[k+2]*w(k,'.5',2) for i,k in enumerate((0,2))),zero),
            sum((h[k+2]*w(k,'.5') for k in range(3)),zero),
            sum((square(h[i])*data['invAm2']*w(k,'0',2) for i,k in enumerate((0,2))),zero)
                -sum((h[k+2]*w(k,'.1')+square(h[k+2])*(w(k,'0',2)/2) for k in range(3)),zero),
            sum((h[k+2]*w(k,'-.9')+square(h[k+2])*(w(k,'-1',2)/2) for k in range(3)),zero)]
        defects=[a+b for a,b in zip(d,change)];unrefined=defects
        terminal=endpoints(x)[0]>=mp.mpf(71)/40
        if terminal:
            if any(not endpoints(v)[0]<=0<=endpoints(v)[1] for row in defects for v in row.coefficients):
                raise ArithmeticError('Terminal implicit-map diagnostic excludes closure')
            # The exact unique h solves the five functional equations. All
            # partial integrals are now full weights of that identical map.
            defects=[zero]*5
        mass=z*4+defects[0]/x
        theta=defects[2]+x**c.mpf('1.6')*c.mpf('.625')
        mixed=(z*theta)*4+defects[1]
        energy=defects[3]-x**c.mpf('1.2')*c.mpf(5)/12
        pressure_moment=defects[4]+x**c.mpf('.2')*c.mpf('2.5')
        p0=self.reference.inputs(Z)['original_axis_pressure'];P=p0+square(am)*pressure_moment
        dz=1-square(z);L=1-square(z)*self.core.delta
        Ur=(2*z*V-(z*mass)*(1-self.core.delta)-dz*derivative(mass))/L
        if endpoints(H[0])[0]<=0 or endpoints(u[0])[0]<=0:raise ArithmeticError('Actual patch positive swirl lost')
        return dict(x=x,Z=Z,Utheta_over_Pstar_axial5=list(u.coefficients),
            Uz_axial5=list(V.coefficients),Utheta_y_over_Pstar_axial5=list(uy.coefficients),Uz_y_axial5=list(Vy.coefficients),
            Ur_over_sqrt_R_over_2_axial4=list(Ur.coefficients),Mz_over_R_axial5=list(mass.coefficients),
            Mtheta_over_sqrt2_Rm_1p5_Am_axial5=list(theta.coefficients),
            Mtheta_z_over_sqrt2_Rm_1p5_Am_axial5=list(mixed.coefficients),
            Mztheta_minus_8ZMz_plus_16Z2R_over_RmAm2_axial5=list(energy.coefficients),
            Mp_over_Am2_axial5=list(pressure_moment.coefficients),P_over_Pstar2_axial5=list(P.coefficients),
            original_P0_axial5=list(p0.coefficients),actual_normalized_five_defects=defects,
            actual_unreferenced_terminal_diagnostic=unrefined,actual_partial_weights={str(k):v for k,v in cache.items()},
            actual_coefficient_functions=h,actual_Rm_history_retained=True,axis_pressure_not_changed=True,
            exact_terminal_closure_from_same_implicit_map=terminal,
            terminal_refinement_is_functional_identity_not_moment_reset=terminal,
            radial_mixed4_certified=False,full_cartesian_vector_derivatives_certified=False,temporal_recursion=False)

    def report(self):
        c=self.ctx;inverse,data=self.coefficients([-1,1])
        t=c.mpf(endpoints(self.repair.t));CSt=t*self.repair.CS
        bw=CSt*(256*self.repair.KN)
        kappa=c.mpf('.9')+(16*CSt)**2/c.mpf('.7')
        margin=8*(c.mpf('.7')-bw)-c.mpf('1.8')
        if endpoints(CSt)[1]>=mp.mpf('.01') or endpoints(bw)[1]>=mp.mpf('.1') or endpoints(kappa)[1]>=1 or endpoints(margin)[0]<=0:
            raise ArithmeticError('Actual connected patch relaxed cone theorem failed')
        samples=[self.evaluate(x,z) for x,z in ((1,[-1,1]),('1.25','0'),('1.5','.5'),('1.75','0'),(2,[-1,1]))]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            reference_join_family_sha256=self.repair.admit['reference_join_family_sha256'],
            actual_connected_source_data=data,actual_implicit_coefficient_inverse_axial5=inverse,
            actual_patch_samples=samples,whole_actual_patch=self.evaluate(c.mpf([1,endpoints(c.exp(1))[1]]),[-1,1]),
            actual_terminal_Rh=self.evaluate(c.exp(1),[-1,1]),
            coefficient_order=['c1','c2','xi1','xi2','xi3'],coefficient_axial_order=5,
            inherited_correlated_C1_e_upper=self.e,inherited_coefficient_C1_t_upper=t,
            CS_times_C1_coefficient_bound=CSt,local_relaxed_bw_upper=bw,local_relaxed_kappa_upper=kappa,
            local_relaxed_cone_lower_margin=margin,
            actual_five_moment_patch_connected=True,actual_patch_coefficient_axial5_enclosures_available=True,
            actual_patch_velocity_pressure_moment_axial5_enclosures_available=True,
            actual_patch_radial_recovery_axial4_available=True,
            actual_five_functional_terminal_identities_connected=True,
            actual_patch_local_relaxed_cone_theorem_bound=True,
            point_coefficient_functions_newly_recomputed=False,original_pressure_datum_retained=True,
            radial_mixed4_certified=False,all_radial_endpoint_joins_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            full_NS_background_completed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    # Bootstrap writes a new receipt; the independent checker accepts it.
    # Ordinary callers require that accepted receipt by default.
    with mp.workdps(280):result=CompliantActualMomentPatch(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual transported five defects connected to implicit patch and axial5 coefficients',flush=True)
    return result


if __name__=='__main__':run()
