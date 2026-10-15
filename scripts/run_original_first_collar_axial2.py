"""Evaluate original first-half collar axial Z2, q_Z and native radial rows."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_first_collar_axial2 import OriginalFirstCollarAxial2
from openai_ns_reconstruction.original_first_collar_transport import Laurent

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_first_collar_axial2.json'

def exact_source_checks():
    """Compare production rows with independently differentiated raw sources."""
    sys.path.insert(0, str(DIRECTORY))
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    ctx = SimpleNamespace(mpf=sy.sympify)
    z, r = sy.symbols('z r', real=True)
    de, e, P2 = sy.Rational(1,17), sy.Rational(1,3), sy.Integer(5)
    a = 1+z*r/7+r*r/11+z**3*r/13
    u = 2+z/3+r*z*z/5+r*r/19
    mean = sy.integrate(u,r)/r
    d,L = 1-z*z,1-de*z*z
    q = (2*z*u-(1-de)*z*mean-d*sy.diff(mean,z))/L
    g = e*(2+2*z)
    T = lambda weight,value:(weight*z*value+d*sy.diff(value,z)-2*z*r*sy.diff(value,r))/L
    flow = (q+(1-de)*z*sy.diff(q,z)/2+r*sy.diff(q,r))/L
    flow += q*(q/2+r*sy.diff(q,r))+u*T(-2,q)
    np = -flow/(2*P2)+(2*sy.diff(q,r)+r*sy.diff(q,r,2))/(e*P2)
    F0 = sy.exp(-z*z-2*z)
    expected = dict(p=2*sy.sqrt(r)*np,
        fx=-2*e*e*T(-3,T(-2-de,F0*a))/F0,
        ux=-2*T(-2,T(-1-de,u))-4*z*r*e*P2*np/L)
    at = {z:sy.Rational(1,5),r:sy.Rational(1,2)}
    def jet(expr, point=at):
        return J(ctx,[sy.diff(expr,z,k).subs(point)/math.factorial(k) for k in range(6)])
    def radial(expr, point=at):
        return [Laurent({0:jet(sy.diff(expr,r,i),point)}) for i in range(3)]
    owner = SimpleNamespace(leading=SimpleNamespace(_derivative=lambda v:
        J(v.ctx,[(k+1)*v[k+1] for k in range(v.order)])))
    v = dict(z=jet(z),rho=jet(r),x=jet(sy.sqrt(r)),de=de,e=e,P2=P2,
        Phi=radial(a),U=radial(u),Q=radial(q),g=jet(g))
    result = OriginalFirstCollarAxial2._forcing(owner,v)
    count = 0
    for name, expr in expected.items():
        for k in (0,1,2):
            assert sy.simplify(result[name].rows[0][k]-sy.diff(expr,z,k).subs(at)/math.factorial(k)) == 0,(name,k)
            count += 1
    # Independent singular radial rows and the physical dx/hb measure.
    at = {z:sy.Rational(1,5),r:sy.Integer(4)}
    sp = sy.Rational(7,3)
    phir, psir = 1+z+z*z, 2-z+z**3
    # Pass the new point explicitly; defaults retain the regular test point.
    v = dict(z=jet(z,at),rho=jet(r,at),x=jet(sy.sqrt(r),at),g=jet(g,at),
        Phi=radial(a,at),U=radial(u,at),Q=radial(q,at),de=de,e=e,P2=P2)
    v['Phi'][2] += Laurent({-1:jet(-sp*phir/r,at)})
    v['U'][2] += Laurent({-1:jet(-sp*e*psir/r,at)})
    v['Q'][2] += Laurent({-1:jet(-2*z*sp*e*psir/(r*L),at)})
    singular = OriginalFirstCollarAxial2._forcing(owner,v)
    expected = dict(p=-8*z*psir/(P2*L),fx=32*e*e*z*z*phir/L**2,
                    ux=64*e*z*z*psir/L**2)
    for name,expr in expected.items():
        measured = singular[name].rows[-1]*v['x']/2/sp
        for k in (0,1,2):
            assert sy.simplify(measured[k]-sy.diff(expr,z,k).subs(at)/math.factorial(k)) == 0,(name,k,'singular')
            count += 1
    return count


def exact_recovery_checks():
    z,r=sy.symbols('z r',positive=True);de=sy.Rational(1,17)
    u=1+z*r+r*r*z*z/3+z**3/5
    mean=sy.integrate(u,r)/r;k=mean-u;ux=2*sy.sqrt(r)*sy.diff(u,r)
    q=(2*z*u-(1+de)*z*mean-(1-z*z)*sy.diff(mean,z))/(1-de*z*z)
    qr=(z*ux/sy.sqrt(r)+(1+de)*z*k/r+(1-z*z)*sy.diff(k,z)/r)/(1-de*z*z)
    assert sy.simplify(sy.diff(mean,r)+k/r)==0
    for order in (0,1):
        assert sy.simplify(sy.diff(qr,z,order)-sy.diff(q,r,z,order))==0
    return 3


def compare(current,previous):
    high=mp.mp.clone();high.dps=180
    difference=abs(high.mpf(current['approximation'])-high.mpf(previous['approximation']))
    error=high.mpf(current['absolute_error_upper'])+high.mpf(previous['absolute_error_upper'])
    assert difference<=error,(difference,error)


def run():
    start=time.perf_counter();exact=exact_source_checks();recovery=exact_recovery_checks()
    print('Independent Z0/Z1/Z2 source, singular measure and radial recovery identities pass',flush=True)
    field=OriginalFirstCollarAxial2(precision=85)
    samples=[];zero_checks=seam_checks=0
    for Z in ('-.3','0','.3'):
        for phase in ('0','.25','.5'):
            packet=field.evaluate(Z,phase)
            if phase=='0':
                for rows in packet['factored_first_collar_increments'].values():
                    for row in rows['coefficient']:
                        assert row['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]]
                        zero_checks+=1
                for name,rows in packet['normalized_original_n1'].items():
                    for a,b in zip(rows,packet['incoming_core'][name]):
                        compare(a,b);seam_checks+=1
            samples.append(packet)
        print('Original collar full Z2 and q_Z: Z='+Z,flush=True)
    samples.append(field.evaluate('.17','.37'))
    previous=json.loads((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_first_collar_transport.json').read_text(encoding='utf-8'))
    legacy_checks=0
    for old in previous['samples']:
        new=next(p for p in samples if p['Z']==old['Z'] and p['phase']==old['phase'])
        for name,rows in old['normalized_original_n1'].items():
            current=new['normalized_original_n1'][name]
            pairs=zip(current,rows) if isinstance(rows,list) else [(current[0],rows)]
            for a,b in pairs:compare(a,b);legacy_checks+=1
    coarse=OriginalFirstCollarAxial2(precision=60);comparisons=0
    for Z in ('-.3','0','.3'):
        old=coarse.evaluate(Z,'.5');new=next(p for p in samples if p['Z']==Z and p['phase']=='.5')
        for section in ('normalized_original_n1','native_first_radial_derivatives'):
            for name,rows in new[section].items():
                for a,b in zip(rows,old[section][name]):compare(a,b);comparisons+=1
        for name in field.fields:
            for a,b in zip(new['factored_first_collar_increments'][name]['coefficient'],old['factored_first_collar_increments'][name]['coefficient']):
                compare(a,b);comparisons+=1
        print('Complete Z2 error comparison: Z='+Z,flush=True)
    example=next(p for p in samples if p['Z']=='.3' and p['phase']=='.5');accuracy={}
    for name in ('p','fx','ux'):
        rows=example['factored_first_collar_increments'][name]['coefficient']
        for k,row in enumerate(rows):
            value=field.ctx.mpf(row['approximation']);error=field.ctx.mpf(row['absolute_error_upper'])
            assert value!=0,(name,k)
            relative=error/abs(value)
            assert relative<field.ctx.mpf('1e-12'),(name,k,relative)
            accuracy[name+str(k)]=field.ctx.nstr(relative,35)
    rejected=0
    for Z,s in (('.3','-.1'),('.3','.5001'),('1','.5'),('.3','nan')):
        try:field.evaluate(Z,s)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError(('Unsupported Z2 collar request admitted',Z,s))
    implementation=dict(field.code_hashes)
    implementation[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    proofs=list({json.dumps(p,sort_keys=True):p for p in field.leading.scale_proofs}.values())
    return dict(source_family=field.family,checks=dict(passed=True,
        independent_raw_source_and_singular_measure_Z0_Z1_Z2_identities=exact,
        independent_radial_recovery_and_mixed_identities=recovery,
        full_original_first_half_collar_packets=len(samples),
        exact_zero_phase_increment_coefficients=zero_checks,incoming_core_Z2_seam_comparisons=seam_checks,
        accepted_collar_Z1_interface_comparisons=legacy_checks,
        cross_precision_complete_Z2_error_comparisons=comparisons,
        nonzero_factored_jump_Taylor_relative_errors=accuracy,rejected_unsupported_requests=rejected,
        regular_remainder_uses_bounded_W_ZZZ=True,numeric_W_ZZZ_not_fabricated=True,
        genuine_leading_axial_order=5,distinct_original_width_log_proofs=len(proofs)),
        samples=samples,original_width_log_proofs=proofs,input_hashes=field.hashes,implementation_hashes=implementation,
        coefficient_convention='D_Z^k/k!: ordinary second derivative and its error are twice order2',
        limitations=dict(admitted_first_phase_max='.5',full_n1_state_Taylor_order=2,mean_Taylor_order=2,q_Taylor_order=1,
            q_rho_Taylor_order=1,q_ZZ_available=False,complete_Cartesian_C2_join=False,
            remaining_first_second_lower_transport=False,physical_mixed_derivatives=False,
            global_moment_repair=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False),
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
