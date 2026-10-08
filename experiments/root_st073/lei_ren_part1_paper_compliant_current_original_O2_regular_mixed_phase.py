"""Actual varying-q O2 slope regular mixed phase and complete A/B products.

The regular branch is proved from the original u, never from its derivatives.
This first source-backed application covers the exact axis across all O2 y.
Signed/near-axis domain coverage and changed five integrals remain separate.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_mixed_source_jets as source

base,prior,ep=source.base,source.prior,source.ep
regular=source.axial.regular
MixedJet=source.MixedJet;C0,Y,Z,YZ=source.C0,source.Y,source.Z,source.YZ
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_original_O2_regular_mixed_phase.json.gz'
RECEIPT=PREFIX+'current_original_O2_regular_mixed_phase_check.json'
GATE='original_O2_varying_q_regular_fixed_phi_mixed_AB_source_enclosed'


def regular_fixed_phi_mixed(a,kernel,roots,coordinate,*,phase):
    """Fixed-angle integrals first, then fixed true phase inverse derivatives."""
    c=a.ctx;const=lambda value:MixedJet.constant(a,value)
    if kernel.geometry!='small_r_series':raise ValueError('Proved original regular source geometry required')
    if kernel.q is not roots['q'][C0]:raise ValueError('Original q root and phase kernel differ')
    for name in ('b','t0'):
        if any(not roots[name][key].zero for key in (C0,Y,Z,YZ)):
            raise ValueError('Only original O2 slope b=t0=0 supported')
    for name in ('a','q','nu'):
        if any(not roots[name][key].zero for key in (Z,YZ)):
            raise ValueError('Original O2 slope parameter Z identities required')
    x,phi=c.mpf(coordinate),c.mpf(phase)
    if any(not 0<=ep(v)[0]<=ep(v)[1]<=1 for v in (x,phi)):
        raise ValueError('Closed source angle and fixed true phase fractions required')
    q=roots['q'];q0=q[C0];h=kernel.hinv
    u0=a.scalar(source.positive.bounded(kernel.u))
    raw_u=roots['p2']*q*const(a.scalar(1).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate()))
    # Enclose finite radial offsets once before long Fourier arithmetic;
    # all P/C/L/q powers and actual derivatives remain.
    ui={key:regular.finite_offset_anchor(a,raw_u[key]) for key in (Y,Z,YZ)}
    h3=h*h*h;h5=h3*h*h
    rr=MixedJet(a,{C0:a.scalar(kernel.r),Y:h3*ui[Y],Z:h3*ui[Z],
        YZ:a.add(h3*ui[YZ],-u0*h5*ui[Y]*ui[Z]*3)})
    hh=MixedJet(a,{C0:h,Y:-u0*h3*ui[Y],Z:-u0*h3*ui[Z],
        YZ:a.add(-u0*h3*ui[YZ],a.add(base.current.square(u0)*h5*3,-h3)*ui[Y]*ui[Z])})
    R=max(abs(v) for v in ep(kernel.r));tails=regular.regular_fourier_tail_bounds(c,R)
    psi0=2*c.pi*x;M=48;powers=[const(1)]
    for k in range(1,M+1):powers.append(powers[-1]*rr)
    W1=const(0);P2=const(2*psi0)
    for k in range(1,M+1):
        sine=c.sin(k*psi0)/k
        W1=W1+powers[k-1]*sine
        coefficient=powers[1]*4 if k==1 else powers[k-2]*(2*(k-1))+powers[k]*(6-2*k)
        P2=P2+coefficient*sine
    W1=W1+regular.tail_jet(a,rr,tails['W1'])
    P2=P2+regular.tail_jet(a,rr,tails['T2_over_q_squared'])
    T1=q*hh*W1*2;T2=q*q*P2
    Dpsi=const(1)-rr*(2*c.cos(psi0))+rr*rr
    t=q*hh*(const(c.cos(psi0))-rr)*Dpsi.reciprocal()*2
    invDpsi=Dpsi.reciprocal()[C0]
    tpsi=-q0*h3*invDpsi*invDpsi*(2*c.sin(psi0))
    q2=base.current.square(q0)
    # Exact source identity 2*q*q_y=nu_y/2 prevents an independent
    # microscopic q^-1 factor from entering the implicit numerator.
    q2y=roots['nu'][Y]*(c.mpf(1)/2)
    FY=a.add(q2*P2[Y],q2y*a.add(P2[C0],-a.scalar(4*c.pi*phi)))
    FZ=q2*P2[Z]
    FYZ=a.add(q2*P2[YZ],q2y*P2[Z])
    D=a.add(a.scalar(1),base.current.square(t[C0]));invD=a.scalar(1).positive_divide(D,0)
    py=-FY*invD;pz=-FZ*invD
    pyz=-a.sum((FYZ,t[C0]*t[Y]*pz*2,t[C0]*t[Z]*py*2,
        t[C0]*tpsi*py*pz*2))*invD
    psi=MixedJet(a,{C0:a.scalar(psi0),Y:py,Z:pz,YZ:pyz})
    total=MixedJet(a,{C0:T1[C0],Y:a.add(T1[Y],t[C0]*py),Z:a.add(T1[Z],t[C0]*pz),
        YZ:a.sum((T1[YZ],t[Y]*pz,t[Z]*py,tpsi*py*pz,t[C0]*pyz))})
    A=roots['a']*(const(phi)-psi*(c.mpf(1)/(2*c.pi)))*(c.mpf(1)/2)
    B=-(roots['a']*roots['E']*total)*(c.mpf(1)/(4*c.pi))
    symmetry=ep(x)[0]==ep(x)[1] and ep(x)[0] in (0,mp.mpf('.5'),1) and ep(phi)==ep(x)
    if symmetry:
        A=const(0);B=const(0)
    return dict(A=A,B=B,record=dict(
        branch='original_O2_regular_varying_q_fixed_true_phi',fixed_true_phase_fraction=phi,
        source_angle_fraction=x,phase_is_held_fixed_not_differentiated=True,
        actual_original_u_mixed_rows=raw_u.record(),original_r_hinv_mixed_jets=dict(r=rr.record(),hinv=hh.record()),
        original_fixed_angle_T1=T1.record(),original_fixed_angle_T2=T2.record(),
        original_fixed_angle_direction=t.record(),original_t_psi=tpsi.record(),
        actual_F_y_Z_yZ={str(key):value.record() for key,value in ((Y,FY),(Z,FZ),(YZ,FYZ))},
        F_definition='psi+T2-2*pi*phi*nu; phi fixed',
        exact_source_cancellation='2*q*q_y=nu_y/2; F_y=q^2*P2_y+(nu_y/2)*(P2-4*pi*phi)',
        actual_nu_y_retained=True,q_y_in_original_direction_and_T1_T2_retained=True,
        original_inverse_psi_mixed=psi.record(),original_total_T1_mixed=total.record(),
        complete_original_A_B_mixed=dict(A=A.record(),B=B.record()),
        all_a_E_first_and_mixed_product_terms_retained=True,
        original_Fourier_tail_bounds=tails,M=M,whole_regular_r_magnitude_upper=R,
        tail_derivatives_are_actual_analytic_partials_not_bound_derivatives=True,
        denominator_one_plus_t_squared_lower=1,no_division_by_r_or_p2=True,
        exact_common_symmetry_trace=symmetry,
        actual_changed_five_integrals_installed=False,
        no_constant_reference_q_a_nu_derivative_assumption=True))


class OriginalO2RegularMixed:
    def __init__(self):
        self.source=source.OriginalO2MixedSources();self.ctx=self.source.ctx;self.family=self.source.family
        receipt=json.loads((HERE/source.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(source.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted actual original O2 mixed source required')
        identity=receipt.get('actual_varying_q_nu_calculus',{})
        if not identity.get('passed') or identity.get('exact_q_and_phase_normalization_nu_identities')!=4:
            raise ValueError('Accepted original q/nu identity and derivative calculus required')
        self.parameter_identity_receipt=source.RECEIPT
        self.hashes=dict(self.source.hashes)
        for name,digest in {**receipt['input_hashes'],source.RECEIPT:sha(source.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original O2 mixed phase dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('O2 mixed phase hash closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.kernels={}

    def kernel(self,frame):
        record=self.source.describe(frame)
        if record['exact_Z_range']!=['0','0']:
            raise ValueError('This issued original O2 regular adapter currently supports exact axis only')
        if id(frame) in self.kernels:return self.kernels[id(frame)]
        roots=frame.roots;a=roots['q'].atlas
        with mp.workdps(a.ctx.dps+40):
            logd=a.copy_interval(self.source.owner.scales.logs['d_star'])
            d=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=logd),1,a.ledger)
            u=(roots['p2'][C0]*roots['q'][C0]).positive_divide(d,logd)
            directional={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()}
            kernel=source.positive.PositiveLogQPhase(dict(q=roots['q'][C0],roots=directional,original_u_source=u),logd)
            if kernel.geometry!='small_r_series':
                raise ValueError('Actual original source requires signed geometry or domain refinement')
            # Exact original nu formula, with the original positive eta kept.
            kernel.nu=roots['nu'][C0]
            kernel.n0=kernel.normalized(kernel.scalar(1),1,True)
            kernel.nt=kernel.normalized(kernel.scalar(0),1,True)
            kernel.nq=kernel.normalized(base.current.square(kernel.q),a.ctx.mpf('.5'),True)
            kernel.ntq=kernel.normalized(kernel.scalar(0),1/(2*a.ctx.sqrt(2)),False)
            self.kernels[id(frame)]=kernel;return kernel

    def whole_phase(self,frame):
        kernel=self.kernel(frame);a=frame.roots['q'].atlas
        with mp.workdps(a.ctx.dps+40):
            result=regular_fixed_phi_mixed(a,kernel,frame.roots,a.ctx.mpf((0,1)),phase=a.ctx.mpf((0,1)))
            return dict(source=self.source.describe(frame),actual_whole_phase_mixed=result['record'],
                original_source_and_true_phase_domains_continuous=True,
                conservative_whole_phase_inverse_graph_outer_enclosure=True,
                independent_angle_phase_rectangle_not_exact_compatible_inverse_graph=True,
                accepted_original_q_nu_identity_receipt=self.parameter_identity_receipt,
                scope='Exact original O2 axis source y cell, whole fixed true phi[0,1], ordinary y/Z/yZ phase and A/B jets; no neighboring-Z coverage claim.')

    def evaluate(self,frame,phase,*,bits=20):
        kernel=self.kernel(frame);a=frame.roots['q'].atlas
        with mp.workdps(a.ctx.dps+40):
            selected=kernel.evaluate(phase,bits=bits)['selected_inverse']
            if selected['chart']!='psi':raise ValueError('Original regular psi inverse required')
            result=regular_fixed_phi_mixed(a,kernel,frame.roots,selected['coordinate_interval'],phase=phase)
            return dict(original_source=self.source.describe(frame),selected_original_inverse=selected,
                actual_fixed_phi_mixed=result['record'],source_inverse_not_selected_field=True,
                accepted_original_q_nu_identity_receipt=self.parameter_identity_receipt)


def run():
    began=time.monotonic();owner=OriginalO2RegularMixed();count=min(owner.source.parent.parent.levels)
    whole=[];selected=[]
    for index in range(count):
        frame=owner.source.source_frame(count,index,Z_lower=0,Z_upper=0)
        whole.append(owner.whole_phase(frame))
        if index in (0,count//2,count-1):
            for phase in ('.137','.663'):
                selected.append(owner.evaluate(frame,phase,bits=20))
    report=dict(**{GATE:True},source_family=owner.family,exact_source_domain=dict(y=['0','1'],Z=['0','0'],phi=['0','1']),
        whole_original_axis_source_cells=whole,actual_true_phase_inverse_queries=selected,
        actual_varying_q_nu_and_full_A_B_mixed_on_regular_axis_installed=True,
        whole_neighboring_Z_or_signed_domain_coverage=False,actual_changed_five_integrals_installed=False,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original O2 slope regular fixed-phi inverse and complete A/B mixed jets, all y at exact Z0; positive original endpoint q, actual nu_y and full products. Not nonzero-Z coverage, signed branch, changed five integrals, terminal, all-route/global N or full reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual varying-q regular O2 axis mixed phase:',count,'whole source/phase cells and',len(selected),'inverse queries',flush=True)
    return report


if __name__=='__main__':run()
