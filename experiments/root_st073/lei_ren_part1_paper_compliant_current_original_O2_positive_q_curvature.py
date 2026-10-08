"""Original signed mixed curvature theorem for every positive O2 q.

Keep q formal in the primitive J^2 bound. The source predicate |u|>=3/16
implies |r|>=3/sqrt(265); that lower is a theorem, not a selected r value.
This module supplies derivative bounds; it does not install signed phases.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_mixed_source_jets as source

base,prior,ep=source.base,source.prior,source.ep
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_original_O2_positive_q_curvature.json.gz'
RECEIPT=PREFIX+'current_original_O2_positive_q_curvature_check.json'
GATE='original_O2_all_positive_q_signed_weighted_mixed_curvature_bound'


def positive_q_weighted_curvature(a,q,*,r_lower):
    c=a.ctx
    if type(q) is not prior.ScaledEnclosure or q.scale.bases is not a.bases or q.ledger is not a.ledger:
        raise ValueError('Same original source atlas and q ledger required')
    if q.zero or ep(q.coefficient)[0]<=0:
        raise ValueError('Strict positive original q factor required; a rounded finite lower may be zero')
    if any(not mp.isfinite(v) for v in ep(q.scale.evaluate())):
        raise ValueError('Finite original positive logq source required')
    rf=c.mpf(r_lower)
    if not 0<ep(rf)[0]<=ep(rf)[1]<=1:raise ValueError('Proved positive signed |r| lower required')
    qf=source.positive.bounded(q)
    if ep(qf)[1]>1:raise ValueError('Original O2 q<=1 required')
    M=c.mpf(max(ep(qf)[1],mp.mpf('.5')));rl=c.mpf(ep(rf)[0]);D0=2*c.pi/rl
    C25=2*(20*c.power(5,c.mpf(5)/4)+c.power(5,c.mpf(3)/4)*D0*D0)
    C2=2*(20*c.power(5,c.mpf(3)/2)+c.power(5,c.mpf(3)/4)*D0*D0)
    C3=5*c.sqrt(20)+c.power(5,c.mpf(3)/4)*D0
    bounds=dict(inverse_J_squared=a.scalar(16*M**5/rl**4*C25),
        primitive_J_squared=q*a.scalar(8*M**4/rl**4*C2),
        inverse_J=a.scalar(8*M**3/rl**2*C3),primitive_J=a.scalar(4*M**3/rl**2*C3))
    return dict(bounds=bounds,record=dict(
        original_q=q.record(),original_q_bounded_conversion=qf,strict_signed_r_lower=rf,
        source_upper_M_not_selected_parameter=M,C25=C25,C2=C2,C3=C3,
        bounds={name:value.record() for name,value in bounds.items()},
        definitions=['h=hinv>0; s=h^2=1-r^2','n=r+coschi; k=abs(n)/h; t=2*q*n/h',
            'D=1+t^2; Dchi=s+2*r*n','t_psi=-2*q*sinchi*Dchi/h^3',
            'J=2*q^2/r^2*(sinchi*(Dchi-3*s)+s*(chi-psi)/r)'],
        positive_q_denominator='m=min(1,2*q)>0; D>=m^2*(1+k^2); D>=1',
        exact_weighted_q_identities=['q^5/m^5=max(q,1/2)^5',
            'q^5/m^4=q*max(q,1/2)^4','q^3/m^3=max(q,1/2)^3'],
        bound_definitions=dict(inverse_J_squared='abs(2*t*t_psi*J^2/D^3)',
            primitive_J_squared='abs((1-t^2)*t_psi*J^2/D^3)',
            inverse_J='abs(2*t*t_psi*J/D^3)',primitive_J='abs((1-t^2)*t_psi*J/D^3)'),
        additional_q_y_transport='F_y=K_y*J+(q_y/q)*(4*pi*phi-2*psi), F_Z=K_Z*J on true inverse graph; |4*pi*phi-2*psi|<=4*pi',
        extra_geometry_proof='|sinchi|*Dchi*P/h^3 <= (h+2*k)^2*(4*h+2*k)+sqrt(h)*(h+2*k)^(3/2)*D0 <= C3*(1+k^2)^(3/2); use D^(5/2),D^2 >= D^(3/2) >= m^3*(1+k^2)^(3/2)',
        accepted_reference_geometry_receipt=source.axial.mixed.RECEIPT,
        reference_q_at_least_half_hypothesis_replaced_by_exact_positive_q_weights=True,
        original_q_factor_retained_in_primitive_J_squared=True,
        no_native_q_lower_floor_or_exact_flat_replacement=True,
        constant_q_and_varying_q_curvature_parts_available=True,
        actual_signed_mixed_phase_primitives_installed=False))


class OriginalO2PositiveQCurvature:
    def __init__(self):
        self.source=source.OriginalO2MixedSources();self.ctx=self.source.ctx;self.family=self.source.family
        accepted=json.loads((HERE/source.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(source.GATE) or accepted['source_family']!=self.family:
            raise ValueError('Accepted original whole O2 mixed source required')
        self.hashes=dict(self.source.hashes)
        for name,digest in {**accepted['input_hashes'],source.RECEIPT:sha(source.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original positive-q curvature dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original positive-q closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def cell(self,count,index):
        frame=self.source.source_frame(count,index,Z_lower=-1,Z_upper=1);a=frame.roots['q'].atlas;c=a.ctx
        with mp.workdps(c.dps+40):
            r_lower=3/c.sqrt(265)
            proof=positive_q_weighted_curvature(a,frame.roots['q'][source.C0],r_lower=r_lower)
            return dict(source_family=self.family,exact_y_cell=[str(frame.left),str(frame.right)],
                exact_outer_Z_window=['-1','1'],conditional_source_predicate='abs(original u=p2*q/dstar)>=3/16',
                predicate_implies_signed_r_lower='abs(r)=abs(u)/sqrt(1+u^2)>=3/sqrt(265)',
                source_q_is_original_positive_function_not_point_value=True,
                actual_conditional_weighted_curvature=proof['record'],
                source_frame_pressure_and_derivative_hash_closure_retained=True,
                not_claimed_full_rectangle_signed_domain=True)


def run():
    began=time.monotonic();owner=OriginalO2PositiveQCurvature();count=min(owner.source.parent.parent.levels)
    rows=[owner.cell(count,index) for index in range(count)]
    report=dict(**{GATE:True},source_family=owner.family,whole_original_O2_positive_q_conditional_bounds=rows,
        original_q_retained_positive_through_endpoint_one=True,
        actual_signed_mixed_phase_primitives_installed=False,regular_signed_domain_coverage_installed=False,
        actual_changed_five_integrals_installed=False,all_17_chart_or_24_cell_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Original O2 signed weighted mixed curvature for every positive q, under actual |u|>=3/16 predicate, all y[0,1] and outer Z[-1,1]. q^5/m^5,q^5/m^4 and q^3/m^3 weights, original q in primitive bound and additional varying-q part. Conditional derivative theorem only, not installed signed inverse/AB, domain union, five integrals or full reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Original all-positive-q signed weighted curvature:',count,'continuous source cells',flush=True)
    return report


if __name__=='__main__':run()
