"""Compute original finite-width n1 values and Z1 on the admitted collar."""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from openai_ns_reconstruction.original_first_collar_transport import (
    OriginalFirstCollarTransport, Laurent)
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY = ROOT/'experiments/root_st073'
NAME = 'lei_ren_part1_paper_compliant_current_original_first_collar_transport.json'


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
        return J(ctx,[sy.diff(expr,z,k).subs(point)/math.factorial(k) for k in range(5)])
    def radial(expr, point=at):
        return [Laurent({0:jet(sy.diff(expr,r,i),point)}) for i in range(3)]
    owner = SimpleNamespace(leading=SimpleNamespace(_derivative=lambda v:
        J(v.ctx,[(k+1)*v[k+1] for k in range(v.order)])))
    v = dict(z=jet(z),rho=jet(r),x=jet(sy.sqrt(r)),de=de,e=e,P2=P2,
        Phi=radial(a),U=radial(u),Q=radial(q),g=jet(g))
    result = OriginalFirstCollarTransport._forcing(owner,v)
    count = 0
    for name, expr in expected.items():
        for k in (0,1):
            assert sy.simplify(result[name].rows[0][k]-sy.diff(expr,z,k).subs(at)) == 0,(name,k)
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
    singular = OriginalFirstCollarTransport._forcing(owner,v)
    expected = dict(p=-8*z*psir/(P2*L),fx=32*e*e*z*z*phir/L**2,
                    ux=64*e*z*z*psir/L**2)
    for name,expr in expected.items():
        measured = singular[name].rows[-1]*v['x']/2/sp
        for k in (0,1):
            assert sy.simplify(measured[k]-sy.diff(expr,z,k).subs(at)) == 0,(name,k,'singular')
            count += 1
    return count


def compare(c, current, previous):
    # Subtract serialized decimals above BOTH output precisions: otherwise
    # cancellation of O(1) values can consume every significant error digit.
    high = mp.mp.clone();high.dps = max(160,c.dps+80)
    difference = abs(high.mpf(current['approximation'])-high.mpf(previous['approximation']))
    error = high.mpf(current['absolute_error_upper'])+high.mpf(previous['absolute_error_upper'])
    assert difference <= error,(difference,error)


def run():
    started = time.perf_counter()
    exact = exact_source_checks()
    print('Raw differential-source and inverse-width measure identities pass',flush=True)
    field = OriginalFirstCollarTransport(precision=85)
    samples, zero_checks = [],0
    for Z in ('-.3','0','.3'):
        for phase in ('0','.25','.5'):
            packet = field.evaluate(Z,phase)
            if phase == '0':
                for name,rows in packet['factored_first_collar_increments'].items():
                    for row in rows['coefficient']:
                        assert row['source_enclosure_exact_mpf_tuples'] == [[0,0,0,0],[0,0,0,0]],name
                        zero_checks += 1
                for name in field.fields:
                    for before,after in zip(packet['incoming_core'][name],packet['normalized_original_n1'][name]):
                        compare(field.ctx,before,after)
            samples.append(packet)
        print('Original full n1 collar values/Z1: Z='+Z,flush=True)
    coarse = OriginalFirstCollarTransport(precision=60)
    comparisons = 0
    for Z in ('-.3','0','.3'):
        old = coarse.evaluate(Z,'.5')
        new = next(p for p in samples if p['Z']==Z and p['phase']=='.5')
        for name,rows in new['normalized_original_n1'].items():
            previous = old['normalized_original_n1'][name]
            pairs = zip(rows,previous) if isinstance(rows,list) else [(rows,previous)]
            for a,b in pairs:compare(field.ctx,a,b);comparisons += 1
        for name in field.fields:
            for a,b in zip(new['factored_first_collar_increments'][name]['coefficient'],
                           old['factored_first_collar_increments'][name]['coefficient']):
                compare(field.ctx,a,b);comparisons += 1
    example = next(p for p in samples if p['Z']=='.3' and p['phase']=='.5')
    accuracy = {}
    for name in ('p','fx','ux'):
        row = example['factored_first_collar_increments'][name]['coefficient'][0]
        point = field.ctx.mpf(row['approximation'])
        relative = field.ctx.mpf(row['absolute_error_upper'])/abs(point)
        assert relative < field.ctx.mpf('1e-15'),(name,relative)
        accuracy[name] = field.ctx.nstr(relative,35)
    assert field.ctx.mpf(example['factored_first_collar_increments']['p']['coefficient'][0]['approximation']) < -1
    rejected = 0
    for Z,phase in (('.3','-.1'),('.3','.5001'),('1','.5'),('.3','nan')):
        try:field.evaluate(Z,phase)
        except (ValueError,ArithmeticError):rejected += 1
        else:raise AssertionError(('Unsupported request admitted',Z,phase))
    code = dict(field.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    proofs = list({json.dumps(p,sort_keys=True):p for p in field.leading.scale_proofs}.values())
    return dict(source_family=field.family,
        checks=dict(passed=True,independent_raw_source_and_singular_measure_identities=exact,
            full_n1_original_first_half_collar_packets=len(samples),
            exact_zero_phase_increment_coefficients=zero_checks,
            cross_precision_complete_error_comparisons=comparisons,
            factored_nonzero_jump_relative_errors=accuracy,rejected_unsupported_requests=rejected,
            actual_core_history_not_reintegrated=True,full_feedback_bounded=True,
            centered_singular_width_variation_included=True,
            distinct_original_width_log_proofs=len(proofs)),
        samples=samples,original_width_log_proofs=proofs,input_hashes=field.hashes,
        implementation_hashes=code,
        remainder_method=dict(local_integral_only=True,dx_over_hb='sqrt(rho)/2',
            principal_pressure='-8 Z Psi_rho sigma / (Pstar^2 L)',
            principal_fx='epsilon^2 * 32 Z^2 Phi_rho sigma / L^2',
            principal_ux='epsilon Pstar^2 * 64 Z^2 (Psi_rho/Pstar^2) sigma / L^2',
            Psi_definition='(Uz_leading-(4Z+j))/epsilon',
            remaining_regular_integral='original hb times row-specific regular ODE/source norm',
            singular_remainder='original hb times centered coefficient-variation norm'),
        limitations=dict(original_finite_width_n1_value_and_Z1_enclosures=True,
            admitted_first_phase_max='.5',exact_original_point_oracle=False,
            remaining_first_second_lower_transport=False,physical_collar_velocity=False,
            global_moment_repair=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-started,3))


def inspect(path):
    result = json.loads(path.read_text(encoding='utf-8'))
    assert result['checks']['passed']
    for name,digest in result['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in result['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args = parser.parse_args();path = DIRECTORY/NAME
    report = inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report['checks'],indent=2),flush=True)
