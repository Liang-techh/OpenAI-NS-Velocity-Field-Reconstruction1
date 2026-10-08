"""Original whole near-midplane mixed source and regular C1 phase integration.

Only Rh_reference[-5,0]: a=4/5, b=t0=0 and q/nu have zero slow
derivatives. O2 parameter variation is not covered. Original four source
rows, smooth regular Fourier jets and complete implicit products precede
the unchanged actual-N own-rate integration program.
"""
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
from types import MappingProxyType,SimpleNamespace
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_near_midplane_integrals as near
import lei_ren_part1_paper_compliant_current_original_reference_mixed_C1_integrals as mixed

whole=near.whole;points=near.points;base=near.base;prior=near.prior;ep=near.ep
HERE,PREFIX,sha=near.HERE,near.PREFIX,near.sha
C0,Y,Z,YZ=mixed.ORDERS;MixedJet=mixed.MixedJet
NAME=PREFIX+'current_original_reference_near_midplane_mixed_C1_integrals.json.gz'
RECEIPT=PREFIX+'current_original_reference_near_midplane_mixed_C1_integrals_check.json'
GATE='original_whole_near_midplane_yZ_source_regular_primitives_C1_integrals_enclosed'


def regular_fourier_tail_bounds(c,R,M=48):
    """Bounds for actual analytic tails and their first two r partials."""
    R=c.mpf(R)
    if ep(R)[0]<0 or ep(R)[1]>c.mpf('.25') or type(M) is not int or M<4:
        raise ValueError('Original bounded regular radius and M>=4 required')
    g=1-R
    S1=lambda n:R**n*((n+1)/g+R/g**2)
    S2=lambda n:R**n*((n+1)*(n+2)/g+(2*n+3)*R/g**2+R*(1+R)/g**3)
    W=(R**M/((M+1)*g),R**(M-1)/g,S1(M-2))
    # T2/q²=2psi + sum c_k*sin(kpsi)/k,
    # c_1=4r; c_k=2(k-1)r^(k-2)+(6-2k)r^k for k>=2.
    # For omitted k>M>=4, both coefficient magnitudes/k are <=2.
    T=(2*(R**(M-1)+R**(M+1))/g,
       2*(S1(M-2)+S1(M)),
       2*(S2(M-3)+S2(M-1)))
    return dict(W1=W,T2_over_q_squared=T)


def tail_jet(a,rr,bounds):
    c=a.ctx
    def symmetric(bound):
        upper=ep(bound)[1]
        return a.scalar(c.mpf((-upper,upper)))
    b0,b1,b2=map(symmetric,bounds)
    return MixedJet(a,{C0:b0,Y:b1*rr[Y],Z:b1*rr[Z],
        YZ:a.add(b1*rr[YZ],b2*rr[Y]*rr[Z])})


def finite_offset_anchor(a,value):
    """Enclose a finite offset ONCE; retain every original formal power.

    Repeated Fourier jet sums must not repeatedly subtract two independent
    copies of the same radial offset hull. This is an enclosing arithmetic
    change, not a selected source value or a derivative of that hull.
    """
    if value.zero:return a.scalar(0)
    if value.scale.bases is not a.bases or value.ledger is not a.ledger:
        raise ValueError('Same original atlas/ledger required')
    if max(abs(v) for v in ep(value.scale.offset))>1000:
        raise ArithmeticError('Collect original large powers before finite offset anchoring')
    return prior.ScaledEnclosure(prior.FormalScale(a.bases,value.scale.powers),
        value.coefficient*value.bounded_exp(value.scale.offset),a.ledger)


def regular_fixed_and_implicit_mixed(a,kernel,roots,coordinate):
    """Whole regular branch, including both signs and r=0; no r inverse."""
    c=a.ctx;const=lambda value:MixedJet.constant(a,value)
    if kernel.geometry!='small_r_series':
        raise ValueError('Original regular reference branch required')
    for name in ('a','b','t0'):
        if any(not roots[name][key].zero for key in (Y,Z,YZ)):
            raise ValueError('Constant original reference parameters required; O2 is separate')
    if not roots['b'][C0].zero or not roots['t0'][C0].zero:
        raise ValueError('Original reference b=t0=0 required')
    x=c.mpf(coordinate)
    if not 0<=ep(x)[0]<=ep(x)[1]<=1:
        raise ValueError('Closed original psi fraction required')
    q=kernel.q;h=kernel.hinv;u=a.scalar(whole.conditioned.bounded(kernel.u))
    raw_ui={key:(roots['p2'][key]*q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate())
        for key in (Y,Z,YZ)}
    ui={key:finite_offset_anchor(a,value) for key,value in raw_ui.items()}
    h3=h*h*h;h5=h3*h*h
    rr=MixedJet(a,{C0:a.scalar(kernel.r),Y:h3*ui[Y],Z:h3*ui[Z],
        YZ:a.add(h3*ui[YZ],-u*h5*ui[Y]*ui[Z]*3)})
    hh=MixedJet(a,{C0:h,Y:-u*h3*ui[Y],Z:-u*h3*ui[Z],
        YZ:a.add(-u*h3*ui[YZ],a.add(base.current.square(u)*h5*3,-h3)*ui[Y]*ui[Z])})
    R=max(abs(v) for v in ep(kernel.r));tails=regular_fourier_tail_bounds(c,R)
    psi=2*c.pi*x;M=48;powers=[const(1)]
    for k in range(1,M+1):powers.append(powers[-1]*rr)
    W1=const(0);P2=const(2*psi)
    for k in range(1,M+1):
        sine=c.sin(k*psi)/k
        W1=W1+powers[k-1]*sine
        coefficient=powers[1]*4 if k==1 else powers[k-2]*(2*(k-1))+powers[k]*(6-2*k)
        P2=P2+coefficient*sine
    W1=W1+tail_jet(a,rr,tails['W1'])
    P2=P2+tail_jet(a,rr,tails['T2_over_q_squared'])
    T1=const(q)*hh*W1*2
    T2=const(base.current.square(q))*P2
    # The original rational direction and its actual fixed-psi slow partials.
    Dpsi=const(1)-rr*(2*c.cos(psi))+rr*rr
    t=const(q)*hh*(const(c.cos(psi))-rr)*Dpsi.reciprocal()*2
    invDpsi=Dpsi.reciprocal()[C0]
    tpsi=-q*h3*invDpsi*invDpsi*(2*c.sin(psi))
    D=a.add(a.scalar(1),base.current.square(t[C0]))
    invD=a.scalar(1).positive_divide(D,0)
    cross=a.add(t[Y]*T2[Z],t[Z]*T2[Y])
    psiYZ=a.sum((-T2[YZ]*invD,t[C0]*cross*invD*invD*2,
        -t[C0]*tpsi*T2[Y]*T2[Z]*invD*invD*invD*2))
    totalY=a.add(T1[Y],-t[C0]*T2[Y]*invD)
    totalZ=a.add(T1[Z],-t[C0]*T2[Z]*invD)
    t2=base.current.square(t[C0])
    totalYZ=a.sum((T1[YZ],-t[C0]*T2[YZ]*invD,
        a.add(t2,-a.scalar(1))*cross*invD*invD,
        a.add(a.scalar(1),-t2)*tpsi*T2[Y]*T2[Z]*invD*invD*invD))
    total=MixedJet(a,{C0:T1[C0],Y:totalY,Z:totalZ,YZ:totalYZ})
    Bmixed=-(roots['a']*roots['E']*total)*(c.mpf(1)/(4*c.pi))
    AYZ=-roots['a'][C0]*psiYZ*(c.mpf(1)/(4*c.pi))
    primitive=kernel.primitives(x,'psi');first={}
    for order in (Y,Z):
        directional={name:{(0,0):row[C0],(0,1):row[order]} for name,row in roots.items()}
        values,proof=points.slow.slow_values(kernel,directional,x,'psi')
        first[order]=dict(A=values['A_Z_slow'],B=values['B_Z_slow'],proof=proof)
    # Exact original periodic and halfperiod traces are zero for every
    # parameter. The full neighborhood is never declared exact-midplane.
    symmetry=ep(x)[0]==ep(x)[1] and ep(x)[0] in (0,mp.mpf('.5'),1)
    if symmetry:AYZ=Bmixed_zero=a.scalar(0)
    else:Bmixed_zero=Bmixed[YZ]
    A=MixedJet(a,{C0:primitive['A'],Y:first[Y]['A'],Z:first[Z]['A'],YZ:AYZ})
    B=MixedJet(a,{C0:primitive['B_over_Pstar'],Y:first[Y]['B'],Z:first[Z]['B'],YZ:Bmixed_zero})
    return dict(A=A,B=B,record=dict(
        branch='original_regular_Fourier_mixed_through_zero',
        original_fixed_angle_T1=T1.record(),original_fixed_angle_T2=T2.record(),
        actual_original_u_mixed_rows={str(key):value.record() for key,value in ui.items()},
        original_u_mixed_rows_before_finite_offset_anchor={str(key):value.record() for key,value in raw_ui.items()},
        finite_radial_source_offset_hulls_enclosed_once_before_Fourier_products=True,
        finite_offset_anchor_retains_all_original_parameter_powers=True,
        original_r_hinv_mixed_jets=dict(r=rr.record(),hinv=hh.record()),
        original_fixed_angle_direction=t.record(),original_t_psi=tpsi.record(),
        original_total_T1_mixed=total.record(),original_inverse_psi_yZ=psiYZ.record(),
        original_A_and_B_mixed=dict(A=A.record(),B=B.record()),
        original_Fourier_tail_bounds=tails,M=M,whole_regular_r_magnitude_upper=R,
        original_regular_polynomial_identity='T2/q^2=2psi+sum c_k*sin(kpsi)/k; c1=4r; ck=2(k-1)r^(k-2)+(6-2k)r^k',
        tail_derivatives_are_actual_analytic_tail_partials_not_bound_derivatives=True,
        original_reference_a_b_t0_q_nu_slow_derivatives_exact_zero=True,
        both_fixed_angle_cross_terms_and_implicit_curvature_retained=True,
        denominator_one_plus_t_squared_lower=1,
        no_division_by_r_or_p2=True,symmetry_trace_exact_zero=symmetry,
        source_enclosures_not_saved_cap_or_midpoint_values=True))


@dataclass(frozen=True)
class OriginalNearMidplaneMixedFrame:
    owner:object
    family:object
    left:object
    right:object
    zeta_bounds:tuple
    roots:object
    query:object
    record:dict


class OriginalReferenceNearMidplaneMixed(mixed.OriginalReferenceMixedC1):
    def __init__(self,*,zeta_lower='-1/1000000',zeta_upper='1/1000000'):
        self.owner=near.OriginalReferenceNearMidplane(zeta_lower=zeta_lower,zeta_upper=zeta_upper)
        self.parent=SimpleNamespace(owner=self.owner)
        self.atlas=a=self.owner.atlas;self.ctx=c=a.ctx;self.family=self.owner.family;self.Z=None
        self.frames={};self.cache={};self.source_cache={};self.primitive_cache={}
        self.hashes=dict(self.owner.hashes)
        for module in (near,mixed):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted original near-source and mixed integration dependencies required')
            for name,digest in {**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Original near mixed dependency changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Near mixed source families differ')
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.templates={};self.terms={};self.compiler_proof={}
        radial=self.owner.templates['inputs'][1:5]
        z=self.owner.z;Q=self.owner.Qsymbol;alpha=self.owner.alpha
        P0,P0Z,P0ZZ=self.owner.templates['pressure_symbols']
        baseline={P0:-alpha/Q**2,P0Z:4*alpha*z/Q**3,P0ZZ:alpha*(4-20*z*z)/Q**4}
        errors=(s.Rational(5,2)/Q**2,10*z/Q**3,22/Q**2)
        for (name,order),rows in self.owner.templates['rows'].items():
            derivative_rows=[];compiled=[];proofs=[]
            for powers,expr in rows:
                derivative=powers[0]*expr+sum(s.diff(expr,x)*rate*x for x,rate in
                    zip(radial,(s.Rational(1,10),s.Rational(1,10),s.Rational(1,5),s.Rational(1,5)),strict=True))
                derivative_rows.append((powers,derivative))
                pieces=[(None,derivative.xreplace(baseline))]
                pieces.extend((j,s.diff(derivative,P)*factor) for j,(P,factor) in
                    enumerate(zip(self.owner.templates['pressure_symbols'],errors,strict=True)))
                for error_order,value in pieces:
                    if value==0:continue
                    terms,proof=self.owner.Z_polynomial(value)
                    actual_powers=powers if error_order is None else (powers[0],powers[1]-1,powers[2],powers[3])
                    for k,coefficient in terms:
                        fn=s.lambdify(self.owner.arguments,coefficient,modules=[{'mpf':c.mpf},'mpmath'])
                        compiled.append(dict(powers=actual_powers,k=k,expr=coefficient,fn=fn,error_order=error_order))
                    proofs.append(dict(original_factor_powers=powers,pressure_error_order=error_order,**proof))
            self.templates[(name,order)]=tuple(derivative_rows)
            self.terms[(name,order)]=tuple(compiled);self.compiler_proof[str((name,order))]=proofs

    def source_frame(self,left,right):
        with mp.workdps(self.ctx.dps+40):
            left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
            if not left<right:raise ValueError('Strict original reference cell required')
            if (left,right) in self.source_cache:return self.source_cache[(left,right)]
            original=self.owner.source_frame(left,right)
            a=self.atlas;c=self.ctx;y=c.mpf((ep(a.rational(left))[0],ep(a.rational(right))[1]))
            zl,zh=a.zeta_bounds;zeta=c.mpf((ep(a.rational(zl))[0],ep(a.rational(zh))[1]))
            f=c.exp(y/10);values=(f,c.mpf(5)/8*f,c.mpf(5)/12*f*f,c.mpf(5)/2*f*f,
                a.copy_interval(self.owner.inputs.alpha_enclosure),a.Q)
            derivatives={};records=[]
            for key,terms in self.terms.items():
                value=a.scalar(0);parts=[]
                for row in terms:
                    coefficient=c.mpf(row['fn'](*values))
                    if row['error_order'] is not None:
                        coefficient=coefficient*c.exp(c.mpf(3)/5)*c.mpf((-1,1))
                    term=a.zterm(row['powers'],coefficient,coordinate=y,zeta=zeta,Z_power=row['k'])
                    value=a.add(value,term)
                    parts.append(dict(original_factor_powers=row['powers'],original_Z_power=row['k'],
                        pressure_remainder_order=row['error_order'],finite_coefficient_enclosure=coefficient,
                        exact_collected_original_derivative_term=term.record()))
                derivatives[key]=value;records.append(dict(original_input=key,ordinary_y_derivative_terms=parts))
            roots={name:MixedJet(a,{C0:original.roots[name][C0],Z:original.roots[name][Z],
                Y:derivatives[(name,0)],YZ:derivatives[(name,1)]}) for name in ('E','V','b','p1','p2')}
            for name in ('a','t0'):roots[name]=MixedJet.constant(a,original.roots[name][C0])
            frame=OriginalNearMidplaneMixedFrame(self,self.family,left,right,a.zeta_bounds,
                MappingProxyType(roots),dict(kernel=original.kernel),
                dict(source_family=self.family,exact_reference_cell=[str(left),str(right)],
                    exact_zeta_interval=[str(v) for v in a.zeta_bounds],
                    original_C0_Z_source_frame=original.record,
                    original_y_yZ_source_terms=records,actual_root_jets={name:row.record() for name,row in roots.items()},
                    original_native_y_and_yZ_template_compiler=self.compiler_proof,
                    original_reference_rates=['1/10','1/10','1/5','1/5'],
                    native_radius_factor_y_derivative_retained=True,
                    original_pressure_datum_y_and_yZ_exact_zero=True,
                    original_L_Z_not_applied_twice=True,ordinary_Z_not_zeta_derivative=True,
                    source_and_pressure_jet_hulls_not_exact_arithmetic_correlation=True))
            self.frames[id(frame)]=frame;self.source_cache[(left,right)]=frame;return frame

    def primitive(self,frame,coordinate):
        if type(frame) is not OriginalNearMidplaneMixedFrame or self.frames.get(id(frame)) is not frame or frame.owner is not self:
            raise ValueError('Issued original near-midplane mixed source frame required')
        with mp.workdps(self.ctx.dps+40):
            x=self.ctx.mpf(coordinate);key=(id(frame),ep(x))
            if key not in self.primitive_cache:
                self.primitive_cache[key]=regular_fixed_and_implicit_mixed(self.atlas,frame.query['kernel'],frame.roots,x)
            return self.primitive_cache[key]

    def integrate(self,*,count,N):
        with mp.workdps(self.ctx.dps+40):
            report=super().integrate(count=count,N=N);report.pop('original_Z_exact')
            report['actual_mixed_source_and_C1_averaging_installed_on_this_fixed_Z_reference_window']=False
            report.update(exact_zeta_interval=[str(v) for v in self.atlas.zeta_bounds],
                physical_Z_map='Z=zeta/(Pstar^11*Cstar^10)',
                actual_mixed_source_and_C1_averaging_installed_on_whole_near_midplane_window=True,
                original_reference_constant_q_nu_required=True,O2_nonconstant_parameter_extension_installed=False,
                native_large_ordinary_Z_and_yZ_factors_retained=True,
                pressure_and_L_interval_hulls_not_exact_joint_function_correlation=True,
                original_yZ_source_not_derivative_of_envelope=True)
            return report


def run():
    begin=time.monotonic();owner=OriginalReferenceNearMidplaneMixed();levels=[]
    for count,N in ((4,160),(16,160),(16,16384)):
        levels.append(owner.integrate(count=count,N=N))
        print('Original whole near-midplane mixed/C1 integration:',count,'cells, N',N,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,
        actual_original_near_midplane_mixed_C1_levels=levels,
        exact_zeta_interval=[str(v) for v in owner.atlas.zeta_bounds],
        original_y_yZ_polynomial_compiler=owner.compiler_proof,
        original_native_C0_y_Z_yZ_source_and_regular_primitives_installed=True,
        no_signed_r_inverse_used_through_zero=True,actual_original_N_C0_Z_averaging_installed=True,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(near.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Actual original Rh_reference[-5,0] whole signed zeta[-1e-6,1e-6] C0/y/Z/yZ sources, smooth regular M48 Fourier mixed jets, fixed-true-phase chain rules, full E/V products and original finite-N own-rate C0/Z phase averaged/direct contribution bounds. Ordinary Z and yZ native factors remain. Separate source/pressure/L hulls, not exact joint arithmetic correlation. Not O2 varying q/nu, full Z atlas, terminal/all-route controls/global N, recursion or corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()
