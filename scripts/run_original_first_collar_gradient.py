"""Evaluate original physical first gradient, vorticity and core-boundary seam."""
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
from openai_ns_reconstruction.original_first_collar_gradient import OriginalFirstCollarGradient
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_first_collar_gradient.json'


def independent_benchmark():
    """Explicit finite manufactured field; differentiate the physical map directly.

    This exercises production assembly without using the original source
    generator or its derivative formulas as the expected result.
    """
    z,r,La=sy.symbols('Z rho Lambda',positive=True)
    de=sy.Rational(1,17);d=1-z*z;L=1-de*z*z
    lam=sy.exp(-1)/sy.sqrt(d);F=sy.exp(-z*z-2*z)
    U0=4*z+sy.Rational(1,19)
    psi=1+z*r+z*z*r*r/5;mean=sy.integrate(psi,r)/r
    u=r*(1+z)+z*z*r*r/7;mu=sy.integrate(u,r)/r
    Q=lambda v,m,b:(2*z*v-b*z*m-d*sy.diff(m,z))/L
    profiles=dict(U0=U0,Phi=1+z*r/13+r*r/17,Psi=psi,
        Qaxis=((1+de)*z*U0-4*d)/L,Qpsi=Q(psi,mean,1-de),
        f=r*(1+z*r/11),u=u,q=Q(u,mu,1+de))
    physical_r=lam*sy.sqrt(2*r/La);physical_z=z*lam**(1-de)
    assert sy.simplify(sy.diff(physical_z,z)-lam**(1-de)*L/d)==0
    assert sy.simplify(sy.diff(physical_r,r)-physical_r/(2*r))==0
    assert sy.simplify(-sy.diff(physical_r,z)/sy.diff(physical_r,r)+2*r*z/d)==0
    ur=lam**-1*sy.sqrt(r/2)*(La**sy.Rational(-1,2)*profiles['Qaxis']
        +La**sy.Rational(-3,2)*(profiles['Qpsi']+lam**(2*de)*profiles['q']))
    ut=F*lam**-1*sy.sqrt(2*r)*(La**sy.Rational(-1,2)*lam**(-de)*profiles['Phi']
        +La**sy.Rational(1,2)*lam**de*profiles['f'])
    uz=lam**-1*(lam**(-de)*U0+La**-1*(lam**(-de)*psi+lam**de*u))
    Dr=lambda v:sy.diff(v,r)/sy.diff(physical_r,r)
    # Differentiate at fixed physical radius before dividing by dz/dZ.
    Dz=lambda v:(sy.diff(v,z)-sy.diff(v,r)*sy.diff(physical_r,z)/sy.diff(physical_r,r))/sy.diff(physical_z,z)
    J=sy.Matrix([[Dr(ur),-ut/physical_r,Dz(ur)],
                 [Dr(ut),ur/physical_r,Dz(ut)],[Dr(uz),0,Dz(uz)]])
    angle=sy.Rational(1,4);R=sy.Matrix([[sy.cos(angle),-sy.sin(angle),0],[sy.sin(angle),sy.cos(angle),0],[0,0,1]])
    omega=sy.Matrix([-Dz(ut),Dz(ur)-Dr(uz),Dr(ut)+ut/physical_r])
    c=mp.mp.clone();c.dps=95;iv=MPIntervalContext();iv.dps=135
    amp=object.__new__(OriginalCoreLogAmplitude);amp.ctx=c;amp.iv=iv;amp.precision=95
    amp.lp=(c.log(17)-30)/4;amp.lp_iv=(iv.ln(17)-30)/4
    amp.ll=970+c.log(17);amp.ll_iv=970+iv.ln(17)
    amp.delta=c.mpf(1)/17;amp.delta_iv=iv.mpf(1)/17
    amp.eps=c.exp(-amp.ll);amp.eps_iv=iv.exp(-amp.ll_iv)
    amp.selected_logC_tuple=[0,0,0,0]
    fake=object.__new__(OriginalFirstCollarGradient)
    fake.ctx=c;fake.iv=iv;fake.amplitude=amp;fake.family={'manufactured_benchmark':True}
    fake._geometry_log=lambda s,m:amp.summarize(c.mpf(0),iv.mpf(0))
    amp.global_amplitude=lambda Z:dict(log_F0_relative_to_anchor=amp.summarize(-c.mpf(Z)**2-2*c.mpf(Z),-iv.mpf(Z)**2-2*iv.mpf(Z)))
    at={z:sy.Rational(1,5),r:sy.Integer(4),La:17*sy.exp(970)}
    def number(expr):return str(expr.subs(at).evalf(150))
    data={}
    for directed,label in ((False,'point'),(True,'truth')):
        ctx=iv if directed else c
        evaluate=lambda expr:ctx.mpf(number(expr))
        prof={name:(evaluate(expr),evaluate(sy.diff(expr,z)),evaluate(sy.diff(expr,r))) for name,expr in profiles.items()}
        zz=ctx.mpf('.2');delta=amp.delta_iv if directed else amp.delta
        eps=amp.eps_iv if directed else amp.eps
        data[label]=dict(profiles=prof,rho=ctx.mpf(4),g=eps*(2*zz+2),delta=delta,d=1-zz*zz,L=1-delta*zz*zz)
    fake._native=lambda Z,s,core_boundary=False:data
    packet=fake.evaluate('.2','0',log_tau='-2',azimuth='.25')
    def resolve(refs):
        total=c.mpf(0)
        for ref in refs:
            term=packet['term_pool'][ref['term']];exponent=c.mpf(0)
            for row in term['log_scale_terms']:
                atom=packet['log_atoms'][row['atom']]
                value=-c.make_mpf(tuple(atom['negative_of_exact_mpf_tuple'])) if 'negative_of_exact_mpf_tuple' in atom else c.mpf(atom['approximation'])
                exponent+=c.mpf(row['multiplier'])*value
            assert abs(exponent)<100000,'Manufactured benchmark only; do not materialize original amplitudes'
            total+=c.mpf(ref['factor']['approximation'])*c.mpf(term['coefficient']['approximation'])*c.exp(exponent)
        return total
    count=0
    def check(refs,expr):
        nonlocal count
        actual=resolve(refs);expected=c.mpf(number(expr))
        difference=abs(actual-expected)
        assert difference<=max(abs(expected),c.mpf(1))*c.mpf('1e-70'),(actual,expected,difference)
        count+=1
    for i in range(3):
        for j in range(3):
            check(packet['gradient_cylindrical'][i][j],J[i,j])
            check(packet['gradient_Cartesian'][i][j],(R*J*R.T)[i,j])
    for i,name in enumerate(('r','theta','z')):check(packet['vorticity_cylindrical'][name],omega[i])
    for i,name in enumerate(('x','y','z')):check(packet['vorticity_Cartesian'][name],(R*omega)[i])
    for name,expr in (('ur',ur),('utheta',ut),('uz',uz)):check(packet['cylindrical_velocity'][name],expr)
    return 3,count


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
    start=time.perf_counter();chain,benchmark=independent_benchmark()
    print('Independent physical-coordinate manufactured gradient/vorticity benchmark passes',flush=True)
    field=OriginalFirstCollarGradient()
    native_reference=reference(field,'compliant_current_original_first_collar_axial2')
    radial_reference=reference(field,'compliant_current_original_first_collar_derivatives')
    samples=[];divergence_checks=native_checks=radial_checks=0
    for Z in ('-.3','0','.3'):
        for phase in ('0','.25','.5'):
            packet=field.evaluate(Z,phase,log_tau='-2',azimuth='.25')
            for row in packet['normalized_meridional_divergence_enclosures'].values():
                lo,hi=_bounds(field.ctx,field._restore(row));assert lo<=0<=hi,(Z,phase,lo,hi)
                divergence_checks+=1
            old=next(p for p in native_reference['samples'] if p['Z']==Z and p['phase']==phase)
            for name in ('f','u','q'):
                for key,k in (('value',0),('Z',1)):
                    overlap(packet['native_source_rows'][name][key],old['normalized_original_n1'][name][k]);native_checks+=1
                overlap(packet['native_source_rows'][name]['rho'],old['native_first_radial_derivatives'][name+'_rho'][0]);native_checks+=1
            old=next(p for p in radial_reference['samples'] if p['Z']==Z and p['phase']==phase)
            for name in ('Phi','Psi'):
                row=old['leading_radial_derivatives'][name+'_rho']['original_hb_powers']['0'][0]
                overlap(packet['native_source_rows'][name]['rho'],row);radial_checks+=1
            samples.append(packet)
        print('Original physical gradient/vorticity: Z='+Z,flush=True)
    samples.append(field.evaluate('.17','.37',log_tau='-3',azimuth='3.2'))
    for row in samples[-1]['normalized_meridional_divergence_enclosures'].values():
        lo,hi=_bounds(field.ctx,field._restore(row));assert lo<=0<=hi;divergence_checks+=1
    seam_checks=0;core_packets=[]
    for Z in ('-.3','0','.3'):
        core=field.core_endpoint(Z,log_tau='-2',azimuth='.25')
        collar=next(p for p in samples if p['native']['Z']==Z and p['native']['phase']=='0')
        assert set(core['term_pool'])==set(collar['term_pool'])
        for name,row in core['term_pool'].items():
            after=collar['term_pool'][name]
            overlap(row['coefficient'],after['coefficient']);seam_checks+=1
            assert row['log_scale_terms']==after['log_scale_terms']
        for key in ('gradient_cylindrical','gradient_Cartesian','vorticity_cylindrical','vorticity_Cartesian'):
            assert core[key]==collar[key],key
        core_packets.append(core)
    print('Core-boundary physical first-derivative and vorticity seam passes',flush=True)
    coarse=OriginalFirstCollarGradient(precision=60,amplitude_precision=90);comparisons=0
    for Z in ('-.3','0','.3'):
        old=coarse.evaluate(Z,'.5',log_tau='-2',azimuth='.25')
        new=next(p for p in samples if p['native']['Z']==Z and p['native']['phase']=='.5')
        for name,row in new['term_pool'].items():
            overlap(row['coefficient'],old['term_pool'][name]['coefficient']);comparisons+=1
            assert row['log_scale_terms']==old['term_pool'][name]['log_scale_terms']
        for name,row in new['log_atoms'].items():
            if 'approximation' in row:overlap(row,old['log_atoms'][name]);comparisons+=1
            else:assert row==old['log_atoms'][name]
        print('Physical gradient coefficient/log error comparison: Z='+Z,flush=True)
    scaled=[]
    for Z in ('-.3','0','.3'):
        later=field.evaluate(Z,'.5',log_tau=None,scaled_time='1e18',azimuth='1.7')
        assert later['axial_vorticity_order_comparison']['available']
        assert later['axial_vorticity_order_comparison']['classification']=='n1_smaller'
        scaled.append(later)
    rejected=0
    for Z,s,kwargs in (('1','.5',{}),('.3','.5001',{}),('.3','-.1',{}),
                       ('.3','.5',dict(log_tau='1')),('.3','.5',dict(scaled_time='1')),
                       ('.3','.5',dict(azimuth='nan'))):
        try:field.evaluate(Z,s,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError(('Unsupported gradient request admitted',Z,s,kwargs))
    implementation=dict(field.code_hashes)
    implementation[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,checks=dict(passed=True,
        independent_similarity_chain_identities=chain,
        independently_differentiated_manufactured_physical_value_gradient_vorticity_checks=benchmark,
        original_collar_packets=len(samples),source_scaled_time_packets=len(scaled),
        normalized_divergence_enclosures_containing_zero=divergence_checks,
        accepted_full_n1_native_value_Z_rho_overlaps=native_checks,
        accepted_leading_first_radial_overlaps=radial_checks,
        core_boundary_physical_term_coefficient_seam_overlaps=seam_checks,
        core_boundary_gradient_and_vorticity_reference_structure_equal=True,
        cross_precision_coefficient_and_log_error_comparisons=comparisons,
        original_F0_and_delta_derivative_scales_separate=True,
        axial_vorticity_n1_smaller_at_tested_source_scaled_times=True,
        rejected_unsupported_requests=rejected),
        samples=samples,core_boundary_samples=core_packets,source_scaled_time_samples=scaled,
        input_hashes=field.hashes,implementation_hashes=implementation,
        limitations=dict(native_phase_required=True,admitted_first_phase_max='.5',nonzero_radius=True,
            physical_spatial_derivative_order=1,ordinary_float_gradient=False,
            uniform_useful_combined_physical_relative_error=False,
            physical_second_derivatives=False,complete_Cartesian_C2_join=False,
            global_matched_field=False,global_moment_repair=False,heat_stress_layers=False,
            temporal_n_ge_2_recursion=False,full_corrected_NS=False),
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
