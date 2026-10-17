"""New source-data/dressing/finite-width/scale controls for swirl recovery."""
from pathlib import Path
import sys
import json
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from openai_ns_reconstruction.original_swirl_source_point import (
    relative_amplitude_from_dressed_moment,relative_amplitude_from_gradient,finite_width_swirl_comparison,scaled_atom_enclosure,
)
from openai_ns_reconstruction.original_signed_C1_source import main_coefficient


def run_fixture():
    c=mp.mp.clone();c.dps=110
    iv=MPIntervalContext();iv.dps=180
    controls=[];differences=[]
    def check(name,condition):
        if not condition:raise AssertionError(name)
        controls.append({'name':name,'passed':True})
    def inside(value,box):
        a,b=[c.make_mpf(v) for v in box._mpi_];return a<=value<=b
    def box(value):return iv.mpf(value)+iv.mpf([-1,1])*iv.mpf('1e-100')
    def exact(x):return {'source_enclosure_exact_mpf_tuples':[list(c.mpf(x)._mpf_)]*2}
    for z in map(c.mpf,('-.3','0','.3')):
        raw=lambda x:2+c.cos(x)**2
        relative=lambda x:c.exp((x-z)/5+(c.sin(x)-c.sin(z))/7)
        dressed=lambda x:raw(x)*relative(x)
        r=[c.diff(raw,z,k) for k in range(4)];d=[c.diff(dressed,z,k) for k in range(4)]
        point=relative_amplitude_from_dressed_moment(d,r,normalized_base_proved=True)
        truth=relative_amplitude_from_dressed_moment(list(map(box,d)),list(map(box,r)),normalized_base_proved=True)
        G=lambda x:-(x/5+c.sin(x)/7)
        gradient=[c.diff(G,z,k+1) for k in range(3)]
        bell=relative_amplitude_from_gradient(gradient,c.mpf('.5'),source_proved=True)
        bell_truth=relative_amplitude_from_gradient(list(map(box,gradient)),iv.mpf('.5'),source_proved=True)
        for k in range(4):
            expected=c.diff(relative,z,k);error=abs(point[k]-expected);differences.append(error)
            check(f'relative_ratio_{z}_Z{k}',error<c.mpf('1e-95'))
            check(f'relative_ratio_enclosure_{z}_Z{k}',inside(expected,truth[k]))
            error=abs(bell[k]-expected);differences.append(error)
            check(f'anchored_gradient_Bell_{z}_Z{k}',error<c.mpf('1e-95'))
            check(f'anchored_gradient_Bell_enclosure_{z}_Z{k}',inside(expected,bell_truth[k]))
        amplitude=lambda x:c.exp(x/5)
        bare=lambda x:c.exp(x)+c.sin(2*x)
        field=lambda x:amplitude(x)*bare(x)
        normalized=[c.diff(field,z,k)/amplitude(z) for k in range(3)]
        delta=c.mpf('.02');mass=c.mpf('.3')
        component=main_coefficient(z,delta,normalized,mass,mass_Z_independent=True)
        full=lambda x:((1-delta)*x*field(x)+(1-x*x)*c.diff(field,x))/(1-delta*x*x)*mass
        for k in range(2):
            error=abs(component[k]*amplitude(z)-c.diff(full,z,k));differences.append(error)
            check(f'frozen_base_C1_{z}_Z{k}',error<c.mpf('1e-95'))
    h=c.mpf('1e-20');y=2*h;w=iv.mpf('.01')
    domain={'joint_core_domain':{'Cauchy_Z_radius':exact(1),'Cauchy_rho_radius':'.5',
             'Xh_norm_majorants':{'Phi':exact(2)},'core_value_factor':4},
            'denominator_bounds':{'Phi_core_modulus_lower':exact('.5')}}
    p=lambda z:c.mpf('.8')+c.sin(z)/50
    rate=lambda z:c.exp(z)/100
    b0=lambda z:c.mpf('.2')+c.cos(z)/100
    c0=lambda z:c.mpf('.3')+c.sin(z)/100
    phi=lambda z:p(z)*c.exp(y*rate(z))
    bm=lambda z:c.exp(-2*y)*b0(z)+p(z)**2*(c.exp(2*rate(z)*y)-c.exp(-2*y))/(2+2*rate(z))
    cm=lambda z:c.exp(-y)*c0(z)+p(z)**2*(c.exp(2*rate(z)*y)-c.exp(-y))/(1+2*rate(z))
    for z in map(c.mpf,('-.3','0','.3')):
        core={n:[iv.mpf(c.diff(fn,z,k)) for k in range(4)] for n,fn in (('Phi',p),('B',b0),('C',c0))}
        truth,proof=finite_width_swirl_comparison(c,iv,core,domain,w,iv.mpf(h),iv.log(iv.mpf(h)),
                                                 source_proved=True,error_cap='1e-6')
        for n,fn,initial in (('Phi',phi,p),('B',bm,b0),('C',cm,c0)):
            for k in range(4):
                actual=c.diff(fn,z,k)
                check(f'finite_width_{z}_{n}_Z{k}',inside(actual,truth[n][k]))
                bound=proof['ordinary_error_proofs'][n][k]['hb_coefficient_bound_exact_mpf_tuples'][1]
                check(f'original_coefficient_{z}_{n}_Z{k}',abs(actual-c.diff(initial,z,k))<=h*c.make_mpf(tuple(bound)))
    coefficient=[iv.mpf([2,3]),iv.mpf([-5,-4]),iv.mpf([0,1])]
    point,truth,proof=scaled_atom_enclosure(c,iv,coefficient,iv.mpf(-100),source_proved=True,error_cap='1e-30')
    for k,value in enumerate((c.mpf('2.5'),c.mpf('-4.5'),c.mpf('.5'))):
        actual=c.exp(-100)*value
        check(f'nonzero_actual_amplitude_{k}',actual!=0 and inside(actual,truth[k]) and point[k]==0)
    check('actual_amplitude_not_zeroed',proof['actual_positive_amplitude_is_zero'] is False and proof['underflow_used_to_define_field'] is False)
    for name,fn in (
        ('missing_ratio_Z3',lambda:relative_amplitude_from_dressed_moment([c.mpf(1)]*3,[c.mpf(1)]*4,normalized_base_proved=True)),
        ('unproved_ratio_base',lambda:relative_amplitude_from_dressed_moment([c.mpf(1)]*4,[c.mpf(1)]*4,normalized_base_proved=False)),
        ('nonpositive_ratio_denominator',lambda:relative_amplitude_from_dressed_moment([c.mpf(1)]*4,[c.mpf(-1)]*4,normalized_base_proved=True)),
        ('missing_G_gradient',lambda:relative_amplitude_from_gradient([c.mpf(1)]*2,c.mpf(1),source_proved=True)),
        ('unproved_G_gradient',lambda:relative_amplitude_from_gradient([c.mpf(1)]*3,c.mpf(1),source_proved=False)),
        ('unproved_scale',lambda:scaled_atom_enclosure(c,iv,coefficient,iv.mpf(-100),source_proved=False)),
        ('inadmissible_scale_cap',lambda:scaled_atom_enclosure(c,iv,coefficient,iv.mpf(1),source_proved=True)),
        ('unproved_comparison',lambda:finite_width_swirl_comparison(c,iv,core,domain,w,iv.mpf(h),iv.mpf(-100),source_proved=False)),
        ('negative_inverse_margin',lambda:finite_width_swirl_comparison(c,iv,core,domain,iv.mpf(1),iv.mpf(h),iv.log(iv.mpf(h)),source_proved=True))):
        try:fn()
        except ValueError:check(name,True)
        else:check(name,False)
    return {'passed':True,'control_count':len(controls),'controls':controls,
            'maximum_supporting_difference':c.nstr(max(differences),25),
            'independent_oracle':'mp.diff of explicit nonconstant relative amplitude; exact exponential Phi/B/C histories; finite nonzero absolute amplitude',
            'original_inverse_or_NS':False,'accepted_ancestor_replayed':False}


if __name__=='__main__':print(json.dumps(run_fixture(),indent=2))
