"""Actual varying-q fixed-phi mixed A/B on the full O2 predicate union.

The angle/phase rectangles are outer covers of the unique true inverse
graph. Curvature cancellation is used on that graph, never on arbitrary
independent angle/phase pairs. Every derivative is an actual source row.
"""
import gzip
import json
from pathlib import Path
from types import FunctionType,SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_predicate_source_frames as frames

source=frames.source;regular=frames.regular;curve=frames.curve
base,prior,ep=frames.base,frames.prior,frames.ep
MixedJet=frames.MixedJet;C0,Y,Z,YZ=frames.C0,frames.Y,frames.Z,frames.YZ
HERE,PREFIX,sha=frames.HERE,frames.PREFIX,frames.sha
NAME=PREFIX+'current_original_O2_predicate_mixed_phase.json.gz'
RECEIPT=PREFIX+'current_original_O2_predicate_mixed_phase_check.json'
GATE='original_O2_full_predicate_union_varying_q_fixed_phi_mixed_AB_enclosed'


def bounded_O2_offset_anchor(a,value):
    """Keep native powers; arbitrary negative offsets are safe small tails.

    The reference helper's absolute1000 cutoff rejected the last two O2
    cells, where the actual flat derivative has a negative log offset.
    This does not materialize a positive large scale or floor original q.
    """
    if value.zero:return a.scalar(0)
    if value.scale.bases is not a.bases or value.ledger is not a.ledger:
        raise ValueError('Same original O2 atlas and ledger required')
    if ep(value.scale.offset)[1]>1000:
        raise ArithmeticError('Collect positive large powers before O2 offset anchoring')
    return prior.ScaledEnclosure(prior.FormalScale(a.bases,value.scale.powers),
        value.coefficient*value.bounded_exp(value.scale.offset),a.ledger)


# Reuse the accepted regular mathematical code unchanged. Only the finite
# offset arithmetic callback is widened to safely enclose negative logs.
_regular_function=regular.regular_fixed_phi_mixed
_regular_env=dict(_regular_function.__globals__,regular=SimpleNamespace(**dict(
    vars(regular.regular),finite_offset_anchor=bounded_O2_offset_anchor)))
regular_fixed_phi_mixed=FunctionType(_regular_function.__code__,_regular_env,
    _regular_function.__name__,_regular_function.__defaults__,_regular_function.__closure__)
regular_fixed_phi_mixed.__kwdefaults__=_regular_function.__kwdefaults__


def symmetric_cap(a,value):
    """An absolute formal bound becomes a signed range with powers intact."""
    lo,hi=ep(value.coefficient)
    if lo<0:raise ValueError('Positive absolute source bound required')
    return prior.ScaledEnclosure(value.scale,a.ctx.mpf((-hi,hi)),a.ledger)


def signed_fixed_phi_mixed(a,kernel,roots,u,gamma_g,gamma_q,coordinate,*,phase):
    c=a.ctx;const=lambda value:MixedJet.constant(a,value)
    if type(kernel) is not frames.O2PredicatePhase or kernel.geometry!='signed_Mobius':
        raise ValueError('Proved original O2 signed source geometry required')
    if kernel.q is not roots['q'][C0] or u.atlas is not a:
        raise ValueError('Same original q/u roots and phase kernel required')
    for name in ('b','t0'):
        if any(not roots[name][key].zero for key in (C0,Y,Z,YZ)):
            raise ValueError('Original O2 slope b=t0=0 required')
    for name in ('a','q','nu'):
        if any(not roots[name][key].zero for key in (Z,YZ)):
            raise ValueError('Original O2 slope transverse parameter identities required')
    x,phi=c.mpf(coordinate),c.mpf(phase)
    if any(not 0<=ep(v)[0]<=ep(v)[1]<=1 for v in (x,phi)):
        raise ValueError('Closed angle and fixed true phase fractions required')
    q=roots['q'];q0=q[C0];h=kernel.hinv;s0=kernel.s_source;r0=a.scalar(kernel.r)
    # u*h=r exactly, before bounding. The gamma_g cover encloses the
    # derivative of the same original carrier, not a selected C1 function.
    KY=r0*a.add(gamma_g,gamma_q)
    KZ=bounded_O2_offset_anchor(a,u[Z]*h)
    KYZ=bounded_O2_offset_anchor(a,u[YZ]*h)
    rr=MixedJet(a,{C0:r0,Y:KY*s0,Z:KZ*s0,YZ:s0*a.add(KYZ,-r0*KY*KZ*3)})
    ss=MixedJet(a,{C0:s0,Y:-r0*KY*s0*2,Z:-r0*KZ*s0*2,
        YZ:-s0*a.add(r0*KYZ,a.add(s0,-base.current.square(r0)*3)*KY*KZ)*2})
    hh=MixedJet(a,{C0:h,Y:-r0*KY*h,Z:-r0*KZ*h,
        YZ:h*a.add(a.add(base.current.square(r0)*3,-a.scalar(1))*KY*KZ,-r0*KYZ)})
    psif,chif=kernel.angles(x,'E');psi0=2*c.pi*psif;chi0=2*c.pi*chif
    sn=a.scalar(c.sin(chi0));cs=a.scalar(c.cos(chi0));psi_fixed=const(psi0)
    chi=MixedJet(a,{C0:a.scalar(chi0),Y:KY*sn*2,Z:KZ*sn*2,
        YZ:sn*a.sum((KYZ,-r0*KY*KZ,KY*KZ*cs*2))*2})
    sine,cosine=chi.sincos();inv_r=rr.reciprocal()
    T1=q*hh*inv_r*(chi-psi_fixed)
    T2=q*q*inv_r*inv_r*((const(2)-ss*3)*chi+ss*psi_fixed+rr*sine*2)
    # Actual q_y and nu_y are retained. At fixed phi, the implicit
    # equation is psi+T2-2*pi*phi*nu=0, and nu_Z=nu_yZ=0.
    FY=a.add(T2[Y],-roots['nu'][Y]*(2*c.pi*phi));FZ=T2[Z];FYZ=T2[YZ]
    factor=a.add(cosine[C0]*2,-r0)
    betaY=a.add(KY*factor,gamma_q);betaZ=KZ*factor
    gammaY=-q0*KY*h*2;gammaZ=-q0*KZ*h*2
    R0=a.scalar(c.mpf((0,1)));Rt=a.scalar(c.mpf(['-.5','.5']))
    Rtt=a.scalar(c.mpf(['0','.5']));R2t=a.scalar(c.mpf((-1,1)))
    G2=a.scalar(c.mpf((-1,ep(c.mpf(1)/8)[1])));G2t=Rt
    certificate=curve.positive_q_weighted_curvature(a,q0,r_lower=3/c.sqrt(265))
    caps={name:symmetric_cap(a,value) for name,value in certificate['bounds'].items()}
    Bnu=a.scalar(c.mpf((-ep(4*c.pi)[1],ep(4*c.pi)[1])))
    # On the true inverse graph: FY=KY*J+(q_y/q)*Bnu,
    # FZ=KZ*J and Bnu=4*pi*phi-2*psi. Both curvature pieces remain.
    inverse_curvature=a.add(caps['inverse_J_squared']*KY*KZ,caps['inverse_J']*gamma_q*KZ*Bnu)
    primitive_curvature=a.add(caps['primitive_J_squared']*KY*KZ,caps['primitive_J']*gamma_q*KZ*Bnu)
    py=-R0*FY;pz=-R0*FZ
    pyz=a.sum((-R0*FYZ,a.add(betaY*Rtt,gammaY*R2t)*FZ,
        a.add(betaZ*Rtt,gammaZ*R2t)*FY,-inverse_curvature))
    total=MixedJet(a,{C0:T1[C0],Y:a.add(T1[Y],-Rt*FY),Z:a.add(T1[Z],-Rt*FZ),
        YZ:a.sum((T1[YZ],-Rt*FYZ,a.add(betaY*G2t,gammaY*G2)*FZ,
            a.add(betaZ*G2t,gammaZ*G2)*FY,primitive_curvature))})
    psi=MixedJet(a,{C0:a.scalar(psi0),Y:py,Z:pz,YZ:pyz})
    A=roots['a']*(const(phi)-psi*(c.mpf(1)/(2*c.pi)))*(c.mpf(1)/2)
    B=-(roots['a']*roots['E']*total)*(c.mpf(1)/(4*c.pi))
    symmetry=ep(x)[0]==ep(x)[1] and ep(x)[0] in (0,mp.mpf('.5'),1) and ep(phi)==ep(x)
    if symmetry:A=const(0);B=const(0)
    return dict(A=A,B=B,record=dict(branch='original_O2_signed_varying_q_fixed_true_phi',
        fixed_true_phase_fraction=phi,source_Mobius_fraction=x,source_angle_fraction=psif,
        original_fixed_angle_T1=T1.record(),original_fixed_angle_T2=T2.record(),
        actual_original_u_mixed_rows=u.record(),actual_source_K_rows={str(key):value.record() for key,value in ((Y,KY),(Z,KZ),(YZ,KYZ))},
        original_r_s_hinv_mixed_jets=dict(r=rr.record(),s=ss.record(),hinv=hh.record()),
        fixed_angle_chi_mixed=chi.record(),actual_F_y_Z_yZ={str(key):value.record() for key,value in ((Y,FY),(Z,FZ),(YZ,FYZ))},
        F_definition='psi+T2-2*pi*phi*nu; phi fixed',
        source_correlated_K_y='K_y=r*(g_y/g+q_y/q), from p2=Z*g and u=p2*q/dstar',
        gamma_g_is_same_source_ratio_enclosure_not_selected_C1_function=True,
        actual_nu_y_retained=True,q_y_in_original_direction_and_T1_T2_retained=True,
        original_inverse_psi_mixed=psi.record(),original_total_T1_mixed=total.record(),
        complete_original_A_B_mixed=dict(A=A.record(),B=B.record()),
        all_a_E_first_and_mixed_product_terms_retained=True,
        phase_is_held_fixed_not_differentiated=True,no_constant_reference_q_a_nu_derivative_assumption=True,
        actual_q_y_direction_beta_and_extra_curvature_retained=True,
        correlated_weighted_curvature=certificate['record'],
        extra_variable_q_curvature=dict(inverse=inverse_curvature.record(),primitive=primitive_curvature.record(),
            exact_source_identity='FY=KY*J+(q_y/q)*(4*pi*phi-2*psi); FZ=KZ*J',
            source_Bnu_absolute_bound=4*c.pi,valid_on_true_inverse_graph_only=True),
        rational_function_ranges=dict(R0=[0,1],Rt=['-1/2','1/2'],Rtt=[0,'1/2'],R2t=[-1,1],G2=[-1,'1/8'],G2t=['-1/2','1/2']),
        positive_original_q_and_hinv_not_replaced_by_flat_limit=True,
        inverse_graph_enclosure_not_arbitrary_independent_angle_phase_identity=True,
        exact_common_symmetry_trace=symmetry,actual_changed_five_integrals_installed=False))


class OriginalO2PredicateMixed:
    def __init__(self):
        self.source=frames.OriginalO2PredicateSources();self.ctx=self.source.ctx;self.family=self.source.family
        receipt=json.loads((HERE/frames.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(frames.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted actual original O2 full predicate source interface required')
        self.hashes=dict(self.source.hashes)
        for name,digest in {**receipt['input_hashes'],frames.RECEIPT:sha(frames.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original O2 predicate phase dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('O2 predicate phase hash closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def primitive(self,frame,coordinate,*,phase):
        # This API accepts only issued conditional source frames. Its result
        # carries the defining predicate; consumers cannot apply it to the
        # entire outer physical rectangle as if the branch held everywhere.
        record=self.source.describe(frame);a=frame.roots['q'].atlas
        with mp.workdps(a.ctx.dps+40):
            if frame.branch=='regular':
                result=regular_fixed_phi_mixed(a,frame.kernel,frame.roots,coordinate,phase=phase)
            else:
                result=signed_fixed_phi_mixed(a,frame.kernel,frame.roots,frame.u,frame.gamma_g,frame.gamma_q,coordinate,phase=phase)
            result['record'].update(defining_original_source_predicate=record['actual_original_source_predicate'],
                conditional_domain_not_entire_outer_rectangle=True,
                negative_flat_source_log_offsets_supported_without_q_floor=True,
                same_original_source_q_u_pressure_histories_and_ordinary_derivatives=True)
            return result

    def whole_phase(self,frame):
        result=self.primitive(frame,self.ctx.mpf((0,1)),phase=self.ctx.mpf((0,1)))
        return dict(source=self.source.describe(frame),actual_whole_phase_mixed=result['record'],
            conservative_whole_true_inverse_graph_outer_enclosure=True,
            independent_angle_phase_rectangle_not_exact_compatible_inverse_graph=True,
            named_original_predicate_required_for_all_downstream_consumers=True)

    def evaluate(self,frame,phase,*,bits=20):
        self.source.describe(frame);kernel=frame.kernel;c=self.ctx
        with mp.workdps(c.dps+40):
            inverse=kernel.evaluate(phase,bits=bits)
            if inverse['status']!='enclosed':raise ArithmeticError('Original conditional phase inverse not enclosed')
            selected=inverse['selected_inverse'];coordinate=selected['coordinate_interval']
            if frame.branch=='regular':
                if selected['chart']!='psi':raise ValueError('Regular original psi inverse required')
            elif selected['chart']!='E':
                coordinate=kernel.angles(c.mpf(coordinate),selected['chart'])[1]
            result=self.primitive(frame,coordinate,phase=phase)
            return dict(original_source=self.source.describe(frame),selected_original_inverse=selected,
                actual_backend_coordinate_fraction=coordinate,actual_backend_chart='psi' if frame.branch=='regular' else 'E',
                actual_fixed_phi_mixed=result['record'],source_inverse_not_selected_field=True)


def run():
    began=time.monotonic();owner=OriginalO2PredicateMixed();count=min(owner.source.source.parent.parent.levels)
    whole=[];selected=[]
    for index in range(count):
        row={}
        for branch in ('regular','positive','negative'):
            frame=owner.source.frame(count,index,branch=branch);row[branch]=owner.whole_phase(frame)
            if index in (0,count//2,count-1):
                for phase in ('.137','.663'):selected.append(owner.evaluate(frame,phase,bits=20))
        whole.append(row)
        if (index+1)%16==0:print('Actual full O2 predicate mixed phase:',index+1,'/',count,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,exact_outer_source_domain=dict(y=['0','1'],Z=['-1','1'],phi=['0','1']),
        whole_original_O2_predicate_mixed_cells=whole,actual_true_phase_inverse_queries=selected,
        regular_positive_negative_original_source_predicate_union_installed=True,
        actual_varying_q_nu_full_A_B_mixed_installed_on_each_named_predicate=True,
        outer_rectangles_not_claimed_entirely_regular_or_signed=True,
        actual_changed_five_integrals_installed=False,all_17_chart_or_24_cell_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual full original O2 regular/signed conditional fixed-phi inverse and A/B C0/y/Z/yZ; named predicates cover y[0,1], Z[-1,1]. Original varying q/a/nu, carrier-correlated y, transverse source jets, all products and positive-q mixed curvature with extra q_y retained. Conservative inverse graph enclosures, not selected fields, integrals, global N or full reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Full original O2 predicate mixed phase:',count*3,'source/phase frames and',len(selected),'inverse queries',flush=True)
    return report


if __name__=='__main__':run()
