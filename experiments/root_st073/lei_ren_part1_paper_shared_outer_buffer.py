"""Same-family O.3 slope-mu and pure-power buffer through the pulse inlet.

All five accumulated primitives are transported from the actual F31 Rd
packet. Only exact relative offsets are used. No axial moment correction,
flatten correction or heat match is claimed by this initial outer segment.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_shared_outer_initial import SharedOuterInitial,stable_sigma
from lei_ren_part1_paper_shared_five_moment_repair import pack
from lei_ren_part1_paper_interval_outer_slope_field import fraction_box
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def transition_kernels(c,t,mu,cells=256):
    """Closed cells with exact exponential weights, enclosing J(t)."""
    t=c.mpf(t)
    if endpoints(t)[0]<0 or endpoints(t)[1]>1 or not isinstance(cells,int) or cells<1:
        raise ValueError('t in [0,1] and positive integer cells required')
    J=c.mpf(0);K=[c.mpf(0),c.mpf(0),c.mpf(0)]
    for i in range(cells):
        a=t*i/cells;b=t*(i+1)/cells;da=b-a
        sa=stable_sigma(c,a)[0];sb=stable_sigma(c,b)[0]
        sj=c.mpf([endpoints(sa)[0],endpoints(sb)[1]])
        nextJ=J+da*sj;jcell=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
        weights=(c.exp(b)-c.exp(a),da,c.exp(-a)-c.exp(-b))
        for n,power in enumerate((1,2,2)):
            K[n]+=weights[n]*c.exp(-mu*power*jcell)
        J=nextJ
    if endpoints(t)==(mp.mpf(1),mp.mpf(1)):J=c.mpf('.5')
    return dict(J=J,theta=K[0],energy=K[1],pressure=K[2],cells=cells)


def decay_integral(c,k,t):
    """Positive integral with a cancellation-free small-exponent bound."""
    if endpoints(k)[0]<=0 or endpoints(t)[0]<0:raise ValueError('k>0,t>=0 required')
    if endpoints(t)==(mp.mpf(0),mp.mpf(0)):return c.mpf(0)
    if endpoints(k*t)[0]>mp.mpf('.001'):
        value=(1-c.exp(-k*t))/k
        return c.mpf([max(mp.mpf(0),endpoints(value)[0]),endpoints(value)[1]])
    return c.mpf([max(mp.mpf(0),endpoints(t*c.exp(-k*t))[0]),endpoints(t)[1]])


class SharedOuterBuffer:
    def __init__(self):
        self.initial=SharedOuterInitial();self.ctx=c=self.initial.ctx;self.params=self.initial.params
        self.hashes=dict(self.initial.hashes)
        for n in ('shared_outer_initial','shared_outer_initial_check'):
            name=PREFIX+n+'.json';record=json.loads((HERE/name).read_bytes())
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Same-family buffer input changed: '+source)
                self.hashes[source]=digest
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def slope_mu(self,Z,offset,cells=256):
        c=self.ctx;t=fraction_box(c,offset);z,qi=self.initial.coordinate(Z);mu=self.params.mu
        inlet=self.initial.axial(Z,buffer_offset=11,cells=cells)
        get=lambda key:IntervalTaylor(c,inlet[key])
        kernels=transition_kernels(c,t,mu,cells)
        f=c.exp(-t/2-mu*kernels['J']);sig=stable_sigma(c,t)[0]
        u1=get('Utheta_over_Pstar');u=u1*f;slope=-c.mpf('.5')-mu*sig
        decay=c.exp(-t);d3=c.exp(-c.mpf('1.5')*t)
        m=get('Mz_over_R')*decay
        h=get('Mtheta_over_sqrt2_R_3half_Pstar')*d3+u1*(d3*kernels['theta'])
        k=get('Mtheta_z_over_sqrt2_R_3half_Pstar')*d3
        e=get('Mztheta_over_R_Pstar_squared')*decay-u1*u1*(decay*kernels['energy']/2)
        p=get('Mp_over_Pstar_squared')+u1*u1*(kernels['pressure']/2)
        packet=self.initial.packet(Z,None,u,u*slope,z*0,z*0,m,h,k,e,p,'O.3 same-family slope-mu transition',
            dict(stage_coordinate=dict(origin='d',local_offset=t),transition_kernels=kernels,
                 accumulated_axial_moment_not_reset_after_Uz_zero=True))
        return packet

    def power(self,Z,phase,cells=256):
        c=self.ctx;phase=fraction_box(c,phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Buffer phase in [0,1] required')
        t=self.params.Tw*phase;z,qi=self.initial.coordinate(Z);mu=self.params.mu
        inlet=self.slope_mu(Z,1,cells);get=lambda key:IntervalTaylor(c,inlet[key])
        slope=-c.mpf('.5')-mu;f=c.exp(slope*t);decay=c.exp(-t);d3=c.exp(-c.mpf('1.5')*t)
        u1=get('Utheta_over_Pstar');u=u1*f
        if endpoints(t)==(mp.mpf(0),mp.mpf(0)):theta_kernel=c.mpf(0)
        else:
            theta_kernel=(f-d3)/(1-mu)
            theta_kernel=c.mpf([max(mp.mpf(0),endpoints(theta_kernel)[0]),endpoints(theta_kernel)[1]])
        m=get('Mz_over_R')*decay
        h=get('Mtheta_over_sqrt2_R_3half_Pstar')*d3+u1*theta_kernel
        k=get('Mtheta_z_over_sqrt2_R_3half_Pstar')*d3
        e=get('Mztheta_over_R_Pstar_squared')*decay-u1*u1*(decay*decay_integral(c,2*mu,t)/2)
        p=get('Mp_over_Pstar_squared')+u1*u1*(decay_integral(c,1+2*mu,t)/2)
        return self.initial.packet(Z,None,u,u*slope,z*0,z*0,m,h,k,e,p,'O.3 same-family pure-power buffer',
            dict(stage_coordinate=dict(origin='w',local_offset=t,phase=phase,total_length=self.params.Tw),
                 accumulated_axial_moment_not_reset_after_Uz_zero=True,
                 pulse_inlet=bool(endpoints(phase)==(mp.mpf(1),mp.mpf(1)))))

    def report(self):
        with mp.workdps(210):
            samples=[self.slope_mu('.5',t) for t in ('0','.5','1')]
            samples += [self.power('.5',t) for t in ('0','.5','1')]
            return dict(actual_five_defect_family_sha256=self.initial.family,
                implicit_source_sha256=self.initial.datum.source_sha,datum_enclosure_sha256=self.initial.datum.datum_sha,
                same_family_O3_callable_through_pulse_inlet=True,samples=samples,
                all_five_moments_transported=True,axis_pressure_unchanged=True,
                short_offsets_and_buffer_length_not_subtracted_from_absolute_radii=True,
                decay_integrals_not_divided_by_mu_after_cancellation=True,
                axial_pulse_and_moment_repairs_installed=False,heat_exterior_matched=False,
                whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    field=SharedOuterBuffer();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Same-family O.3 slope-mu and power buffer through pulse inlet callable; 6 examples',flush=True)
    return result


if __name__=='__main__':run()
