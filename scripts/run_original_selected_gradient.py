"""Check selected cut/own physical values, gradient, vorticity and width sectors."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_selected_gradient import OriginalSelectedGradient
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_selected_gradient.json'


def independent_benchmark():
    """Differentiate a finite explicit physical field with a log-radial cutoff.

    A cubic transition is used only for this independent algebra benchmark.
    Its own primitive is integrated exactly before Q recovery. Finite hb and
    phase_in make both inverse-width derivative terms measurable.
    """
    z,r,q,La=sy.symbols('Z rho q Lambda',positive=True)
    hb=sy.Rational(1,7);a=sy.Rational(1,11);de=sy.Rational(1,17)
    d=1-z*z;L=1-de*z*z;lam=sy.exp(-1)/sy.sqrt(d);F=sy.exp(-z*z-2*z)
    arg=4*q-2;chi=1-3*arg**2+2*arg**3
    mass=sy.integrate(chi,q)
    raw_u=(1+z)/r;u=raw_u*chi
    mean_u=(1+z*z+hb*a*(1+z)*mass)/r
    psi=1+z*r+z*z*r*r/5;mean=sy.integrate(psi,r)/r
    U0=4*z+sy.Rational(1,19)
    Q=lambda v,m,b:(2*z*v-b*z*m-d*sy.diff(m,z))/L
    raw_f=1+z*r/11;f=raw_f*chi
    profiles=dict(U0=U0,Phi=1+z*r/13+r*r/17,Psi=psi,
        Qaxis=((1+de)*z*U0-4*d)/L,Qpsi=Q(psi,mean,1-de),
        f=f,u=u,q=Q(u,mean_u,1+de))
    mr=(u-mean_u)/r
    ur0=sy.diff(raw_u,r)*chi;urs=raw_u*sy.diff(chi,q)/r
    qr0=(2*z*ur0-(1+de)*z*mr-d*sy.diff(mr,z))/L
    qrs=2*z*urs/L
    regular={name:sy.diff(expr,r) for name,expr in profiles.items()}
    regular.update(f=sy.diff(raw_f,r)*chi,u=ur0,q=qr0)
    singular=dict(f=raw_f*sy.diff(chi,q)/r,u=urs,q=qrs)
    q_of_r=sy.log(r/4)/(hb*a)
    actual={name:expr.subs(q,q_of_r) for name,expr in profiles.items()}
    physical_r=lam*sy.sqrt(2*r/La);physical_z=z*lam**(1-de)
    ur=lam**-1*sy.sqrt(r/2)*(La**sy.Rational(-1,2)*actual['Qaxis']
        +La**sy.Rational(-3,2)*(actual['Qpsi']+lam**(2*de)*actual['q']))
    ut=F*lam**-1*sy.sqrt(2*r)*(La**sy.Rational(-1,2)*lam**(-de)*actual['Phi']
        +La**sy.Rational(1,2)*lam**de*actual['f'])
    uz=lam**-1*(lam**(-de)*U0+La**-1*(lam**(-de)*psi+lam**de*actual['u']))
    Dr=lambda v:sy.diff(v,r)/sy.diff(physical_r,r)
    Dz=lambda v:(sy.diff(v,z)-sy.diff(v,r)*sy.diff(physical_r,z)/sy.diff(physical_r,r))/sy.diff(physical_z,z)
    J=sy.Matrix([[Dr(ur),-ut/physical_r,Dz(ur)],
        [Dr(ut),ur/physical_r,Dz(ut)],[Dr(uz),0,Dz(uz)]])
    angle=sy.Rational(1,4);R=sy.Matrix([[sy.cos(angle),-sy.sin(angle),0],
        [sy.sin(angle),sy.cos(angle),0],[0,0,1]])
    omega=sy.Matrix([-Dz(ut),Dz(ur)-Dr(uz),Dr(ut)+ut/physical_r])
    assert sy.simplify(L*r*singular['q']-2*r*z*singular['u'])==0
    c=mp.mp.clone();c.dps=95;iv=MPIntervalContext();iv.dps=135
    amp=object.__new__(OriginalCoreLogAmplitude);amp.ctx=c;amp.iv=iv;amp.precision=95
    amp.lp=(c.log(17)-30)/4;amp.lp_iv=(iv.ln(17)-30)/4
    amp.ll=970+c.log(17);amp.ll_iv=970+iv.ln(17)
    amp.delta=c.mpf(1)/17;amp.delta_iv=iv.mpf(1)/17
    amp.eps=c.exp(-amp.ll);amp.eps_iv=iv.exp(-amp.ll_iv)
    amp.selected_logC_tuple=[0,0,0,0]
    amp.global_amplitude=lambda Z:dict(log_F0_relative_to_anchor=amp.summarize(
        -c.mpf(Z)**2-2*c.mpf(Z),-iv.mpf(Z)**2-2*iv.mpf(Z)))
    fake=object.__new__(OriginalSelectedGradient)
    fake.ctx=c;fake.iv=iv;fake.amplitude=amp;fake.family={'manufactured_benchmark':True}
    fake._selected_geometry=lambda fraction:amp.summarize(c.mpf(fraction)/154,iv.mpf(fraction)/154)
    fake._log_width=lambda:amp.summarize(-c.log(7),-iv.ln(7))
    fake._phase_log=lambda:amp.summarize(-c.log(11),-iv.ln(11))
    at={z:sy.Rational(1,5),r:4*sy.exp(hb*a*sy.Rational(5,8)),q:sy.Rational(5,8),La:17*sy.exp(970)}
    def number(expr):return str(expr.subs(at).evalf(150))
    data={'cutoff_packet':dict(chart='manufactured_cutoff',geometry=dict(
        exact_sc_mpf_tuple=None,exact_phase_in_mpf_tuple=None),normalized_own_primitive={},cutoff={})}
    for directed,label in ((False,'point'),(True,'truth')):
        ctx=iv if directed else c;evaluate=lambda expr:ctx.mpf(number(expr))
        zz=ctx.mpf('.2');delta=amp.delta_iv if directed else amp.delta
        eps=amp.eps_iv if directed else amp.eps
        data[label]=dict(profiles={name:(evaluate(expr),evaluate(sy.diff(expr,z)),evaluate(regular[name]))
            for name,expr in profiles.items()},radial_singular={name:evaluate(expr) for name,expr in singular.items()},
            rho=evaluate(r),g=eps*(2*zz+2),delta=delta,d=1-zz*zz,L=1-delta*zz*zz)
    fake._native=lambda Z,s,lower=False:data
    packet=fake.evaluate('.2','.625',log_tau='-2',azimuth='.25')
    def resolve(refs):
        total=c.mpf(0)
        for ref in refs:
            term=packet['term_pool'][ref['term']];exponent=c.mpf(0)
            for row in term['log_scale_terms']:
                atom=packet['log_atoms'][row['atom']]
                value=-c.make_mpf(tuple(atom['negative_of_exact_mpf_tuple'])) if 'negative_of_exact_mpf_tuple' in atom else c.mpf(atom['approximation'])
                exponent+=c.mpf(row['multiplier'])*value
            assert abs(exponent)<100000,'Only the manufactured field may be materialized'
            total+=c.mpf(ref['factor']['approximation'])*c.mpf(term['coefficient']['approximation'])*c.exp(exponent)
        return total
    checks=0
    def check(refs,expr):
        nonlocal checks
        got=resolve(refs);expected=c.mpf(number(expr))
        assert abs(got-expected)<=max(abs(expected),c.mpf(1))*c.mpf('1e-65'),(got,expected)
        checks+=1
    for i in range(3):
        for j in range(3):
            check(packet['gradient_cylindrical'][i][j],J[i,j])
            check(packet['gradient_Cartesian'][i][j],(R*J*R.T)[i,j])
    for i,name in enumerate(('r','theta','z')):check(packet['vorticity_cylindrical'][name],omega[i])
    for i,name in enumerate(('x','y','z')):check(packet['vorticity_Cartesian'][name],(R*omega)[i])
    for name,expr in (('ur',ur),('utheta',ut),('uz',uz)):check(packet['cylindrical_velocity'][name],expr)
    for i,name in enumerate(('u','v','w')):check(packet['Cartesian_velocity'][name],(R*sy.Matrix([ur,ut,uz]))[i])
    return checks


def overlap(a,b):
    c=mp.mp.clone();c.dps=180
    diff=abs(c.mpf(a['approximation'])-c.mpf(b['approximation']))
    error=c.mpf(a['absolute_error_upper'])+c.mpf(b['absolute_error_upper'])
    assert diff<=error,(diff,error)


def reference(field,stem):
    path=DIRECTORY/('lei_ren_part1_paper_'+stem+'.json')
    report=json.loads(path.read_text(encoding='utf-8'))
    assert report['checks']['passed'] and report['source_family']==field.family
    for name,digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        field.code_hashes[name]=digest
    field.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    return report


def run():
    start=time.perf_counter();benchmark=independent_benchmark()
    print('Independent finite-width cutoff physical differentiation passes',flush=True)
    field=OriginalSelectedGradient();accepted=reference(field,'compliant_current_original_n1_cutoff_numeric')
    previous=reference(field,'compliant_current_original_first_collar_gradient')
    samples=[];source_checks=divergence_checks=flat_checks=width_checks=0
    requests=[('.3',q,False) for q in ('0','.5','.625','.75','1')]
    requests += [(Z,'.625',False) for Z in ('0','-.3')]
    requests += [('.3',q,True) for q in ('1','1.5','2')]
    requests += [(Z,'2',True) for Z in ('0','-.3')]
    for Z,fraction,lower in requests:
        evaluate=field.lower_collar if lower else field.evaluate
        packet=evaluate(Z,fraction,log_tau='-2',azimuth='.25')
        assert len(packet['term_pool'])==43
        for rows in packet['normalized_meridional_divergence_width_sectors'].values():
            for row in rows.values():
                lo,hi=_bounds(field.ctx,field._restore(row));assert lo<=0<=hi,(Z,fraction,lo,hi)
                divergence_checks+=1
        source=next(p for p in accepted['lower_collar_samples' if lower else 'selected_inner_samples']
                    if p['Z']==Z and p['coordinate']==fraction)
        for name in ('f','u','q'):
            row=packet['native_source_rows'][name]
            for key,k in (('value',0),('Z',1)):
                overlap(row[key],source['normalized_original_n1'][name][k]);source_checks+=1
            for sector in ('0,0','-1,-1'):
                overlap(row['rho_width_sectors'][sector],source['native_first_radial_width_sectors'][name+'_rho'][sector][0]);source_checks+=1
        for name,term in packet['term_pool'].items():
            logs={r['atom']:r['multiplier'] for r in term['log_scale_terms']}
            is_width='_width:' in name
            if is_width:
                assert field.ctx.mpf(logs['original_log_hb'])==-1 and field.ctx.mpf(logs['log_phase_in'])==-1;width_checks+=1
            else:assert 'original_log_hb' not in logs and 'log_phase_in' not in logs
            if (lower or field.ctx.mpf(fraction)>=field.ctx.mpf('.75')) and term['source'] in ('f','u'):
                lo,hi=_bounds(field.ctx,field._restore(term['coefficient']));assert lo==hi==0;flat_checks+=1
            if (lower or not field.ctx.mpf('.5')<field.ctx.mpf(fraction)<field.ctx.mpf('.75')) and is_width:
                lo,hi=_bounds(field.ctx,field._restore(term['coefficient']));assert lo==hi==0;flat_checks+=1
        samples.append(packet)
        print('Selected physical field: Z='+Z+' fraction='+fraction+' lower='+str(lower),flush=True)
    by=lambda q,lower=False:next(p for p in samples if p['native']['Z']=='.3' and p['native']['fraction']==q and ('Rin_to_rminus' in p['native']['chart'])==lower)
    seam_checks=0
    core=next(p for p in previous['core_boundary_samples'] if p['native']['Z']=='.3')
    inner=by('0')
    for name,row in core['term_pool'].items():
        overlap(row['coefficient'],inner['term_pool'][name]['coefficient']);seam_checks+=1
        assert row['log_scale_terms']==inner['term_pool'][name]['log_scale_terms']
    for name,row in by('1')['term_pool'].items():
        other=by('1',True)['term_pool'][name]
        overlap(row['coefficient'],other['coefficient']);seam_checks+=1
        assert row['log_scale_terms']==other['log_scale_terms']
    assert by('1')['log_atoms']==by('1',True)['log_atoms']
    assert by('1')['gradient_Cartesian']==by('1',True)['gradient_Cartesian']
    assert by('1')['own_history']['primitive']==by('2',True)['own_history']['primitive']
    for name in ('q',):
        lo,hi=_bounds(field.ctx,field._restore(by('2',True)['native_source_rows'][name]['value']))
        assert hi<0 or lo>0
    comparisons=0;coarse=OriginalSelectedGradient(precision=60,amplitude_precision=90)
    for fraction,lower in (('.625',False),('2',True)):
        evaluate=coarse.lower_collar if lower else coarse.evaluate
        old=evaluate('.3',fraction,log_tau='-2',azimuth='.25');new=by(fraction,lower)
        for name,row in new['term_pool'].items():
            overlap(row['coefficient'],old['term_pool'][name]['coefficient']);comparisons+=1
            assert row['log_scale_terms']==old['term_pool'][name]['log_scale_terms']
        for name,row in new['log_atoms'].items():
            if 'approximation' in row:overlap(row,old['log_atoms'][name]);comparisons+=1
            elif name=='original_log_hb':
                x0,x1=_bounds(field.ctx,field._restore(row));y0,y1=_bounds(field.ctx,coarse._restore(old['log_atoms'][name]))
                assert max(x0,y0)<=min(x1,y1);comparisons+=1
            else:assert row==old['log_atoms'][name]
    scaled=field.lower_collar('.3','2',log_tau=None,scaled_time='1e18',azimuth='1.7')
    rejected=0
    for lower,Z,q,kwargs in ((False,'1','.625',{}),(False,'.3','1.001',{}),(False,'.3','-.1',{}),
            (True,'.3','.9',{}),(True,'.3','2.001',{}),(False,'.3','.625',{'log_tau':'1'}),
            (False,'.3','.625',{'scaled_time':'1'}),(False,'.3','.625',{'azimuth':'nan'})):
        try:(field.lower_collar if lower else field.evaluate)(Z,q,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError(('Unsupported selected physical request admitted',Z,q,lower,kwargs))
    implementation=dict(field.code_hashes)
    implementation[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,checks=dict(passed=True,
        independently_differentiated_finite_width_physical_value_gradient_vorticity_checks=benchmark,
        actual_selected_physical_packets=len(samples),accepted_cut_own_native_value_Z_radial_overlaps=source_checks,
        normalized_divergence_width_sector_enclosures_containing_zero=divergence_checks,
        exact_zero_postsupport_or_nontransition_slots=flat_checks,
        inverse_hb_phase_in_scale_terms_correct=width_checks,core_and_Rin_term_seam_overlaps=seam_checks,
        Rin_physical_gradient_reference_structure_equal=True,terminal_own_primitive_retained=True,
        lower_q_nonzero_after_f_and_u_cut=True,cross_precision_coefficient_and_log_overlaps=comparisons,
        source_scaled_time_supported=True,rejected_unsupported_requests=rejected),
        samples=samples,source_scaled_time_samples=[scaled],input_hashes=field.hashes,
        implementation_hashes=implementation,
        limitations=dict(inner_fraction_range=['0','1'],lower_fraction_range=['1','2'],stops_at_rminus=True,
            physical_derivative_order=1,ordinary_float_gradient=False,useful_combined_physical_relative_error=False,
            physical_second_derivatives=False,complete_Cartesian_C2_join=False,global_matched_field=False,
            global_moment_repair=False,heat_stress_layers=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-start,3))


def inspect(path):
    report=json.loads(path.read_text(encoding='utf-8'));assert report['checks']['passed']
    for name,digest in report['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    report=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds']),indent=2),flush=True)
