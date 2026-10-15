"""Compute physical leading+n1 core velocity packets and directed inverses."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_core_velocity import OriginalCoreVelocity
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_core_velocity.json'


def exact_streamfunction():
    R,Z,lam,de,beta=s.symbols('R Z lambda delta beta',positive=True)
    d,L=1-Z*Z,1-de*Z*Z
    M=s.Function('M')(R,Z)
    mr,mz=s.diff(M,R),s.diff(M,Z)
    # At fixed physical r,t: R_Z=-2ZR/d, lambda_Z=lambda*Z/d,
    # z_Z=lambda^(1-delta)*L/d. The independent mass is arbitrary.
    radial_numerator=2*Z*R*mr-beta*Z*M-d*mz
    radial_div=lam**(beta-3+de)*s.diff(radial_numerator,R)/L
    axial_div=lam**(beta-3+de)*((beta-2)*Z*mr+d*s.diff(mr,Z)-2*Z*R*s.diff(mr,R))/L
    assert s.simplify(radial_div+axial_div)==0
    comparisons=1
    for order,weight in ((0,1-de),(1,1+de)):
        uz_weight=s.simplify((beta-2).subs(beta,weight))
        ur_weight=s.simplify((beta-2+de).subs(beta,weight))
        assert uz_weight == (-1-de if order==0 else -1+de)
        assert ur_weight == (-1 if order==0 else -1+2*de)
        comparisons+=2
    eps,rho=s.symbols('epsilon rho',positive=True)
    assert s.simplify(s.sqrt(eps*rho/2)*eps-s.sqrt(rho/2)*eps**s.Rational(3,2))==0
    assert s.simplify(s.sqrt(2*eps*rho)/eps-s.sqrt(2*rho/eps))==0
    comparisons+=2
    # Every axisymmetric swirl cancels independently from Cartesian div.
    x,y=s.symbols('x y');T=s.Function('T')((x*x+y*y)/2,Z)
    assert s.simplify(s.diff(-T*y,x)+s.diff(T*x,y))==0
    return comparisons+1


def packet_ratio(field,packet):
    """Derive the order ratio from returned terms after common-log cancellation."""
    c,iv=field.ctx,field.iv
    leading,n1=packet['cylindrical']['utheta']
    key=lambda row:json.dumps(row,sort_keys=True)
    old=Counter(key(row) for row in leading['log_scale_terms'])
    new=Counter(key(row) for row in n1['log_scale_terms'])
    common=old&new;old-=common;new-=common
    point=c.log(abs(c.mpf(n1['coefficient']['approximation'])/c.mpf(leading['coefficient']['approximation'])))
    truth=iv.ln(abs(field._restore(n1['coefficient']))/field._restore(leading['coefficient']))
    for sign,rows in ((1,new),(-1,old)):
        for text,count in rows.items():
            row=json.loads(text)
            assert 'approximation' in row,'Huge shared amplitude logs must cancel exactly'
            point+=sign*count*c.mpf(row['approximation'])
            truth+=sign*count*field._restore(row)
    saved=packet['swirl_order_comparison']['log_absolute_ratio']
    saved_box=field._restore(saved)
    a,b=_bounds(c,truth);lo,hi=_bounds(c,saved_box)
    assert max(a,lo)<=min(b,hi)
    error=c.mpf(saved['absolute_error_upper'])+max(abs(v) for v in _bounds(c,truth-iv.mpf(c.nstr(point,field.amplitude.precision))))
    assert abs(point-c.mpf(saved['approximation']))<=error
    return dict(independent_log_ratio=c.nstr(point,40),huge_common_logs_cancelled=True,
        classification=packet['swirl_order_comparison']['classification'])


def physical_request(Z,rho,angle,logtau):
    """Generate exact decimal test requests, not an assumed exact round trip."""
    c=mp.mp.clone();c.dps=240
    z,r,p,lt=map(c.mpf,(Z,rho,angle,logtau))
    lp=c.exp(40)+11;ll=4*lp+1000;de=c.exp(-4*lp-30)
    ell=(lt-c.log(1-z*z))/2
    radius=c.mpf(0) if not r else c.exp(ell+(c.log(2*r)-ll)/2)
    return tuple(c.nstr(v,230) for v in (radius*c.cos(p),radius*c.sin(p),z*c.exp((1-de)*ell)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args()
    if args.inspect_report:
        r=json.loads((DIRECTORY/NAME).read_text(encoding='utf-8'))
        for name,digest in r['implementation_hashes'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        for name,digest in r['input_hashes'].items():
            assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
        assert r['checks']['passed']
        print(json.dumps(dict(checks=r['checks'],archive_only=True),indent=2));return
    started=time.monotonic()
    exact=exact_streamfunction()
    field=OriginalCoreVelocity()
    native,scaled,ratios=[],[],[]
    for Z in ('-0.3','0','0.3'):
        for rho in ('0','2','4'):
            packet=field.evaluate(Z,rho,log_tau='-2',azimuth='0.25')
            native.append(packet)
            assert packet['coefficient_orders_included']==[0,1]
            assert len(packet['Cartesian']['u'])==len(packet['Cartesian']['v'])==5
            assert len(packet['Cartesian']['w'])==len(packet['streamfunction'])==3
            if rho=='0':
                assert all(field.ctx.mpf(row['coefficient']['approximation'])==0 for key in ('u','v') for row in packet['Cartesian'][key])
                assert field.ctx.mpf(packet['Cartesian']['w'][-1]['coefficient']['approximation'])==0
            else:
                ratios.append(packet_ratio(field,packet))
                later=field.evaluate(Z,rho,log_tau=None,scaled_time='1e18',azimuth='0.25')
                scaled.append(later);ratios.append(packet_ratio(field,later))
                assert packet['swirl_order_comparison']['classification']=='n1_larger'
                assert later['swirl_order_comparison']['classification']=='n1_smaller'
                assert any(row['name']=='leading_axial_correction' for row in later['Cartesian']['w'])
        print('Computed physical core/time packets at Z='+Z,flush=True)
    inverse=[]
    for Z,rho,angle in (('0','0','0'),('0.3','0','0'),('-0.3','2','1.7'),('0.17','1.5','3.2')):
        request=physical_request(Z,rho,angle,'-2')
        packet=field.cartesian(*request,log_tau='-2')
        assert packet['physical_inverse']['source_input_boxes_contain_entire_inverse']
        assert packet['physical_inverse']['no_uncontrolled_interpolation_from_native_point']
        inverse.append(packet)
        print(f'Computed directed Cartesian inverse near Z={Z}, rho={rho}',flush=True)
    physical_time=[]
    c=mp.mp.clone();c.dps=240
    request=physical_request('0.17','1.5','0.4',c.nstr(c.log(c.mpf('0.5')),230))
    for values in (('0','0','0'),request):
        packet=field.at_time(*values,'0.5')
        assert packet['time']['exact_t']=='0.5'
        assert packet['physical_inverse']['source_time_boxes_contain_true_log_tau']
        physical_time.append(packet)
    rejected=0
    for call in (lambda:field.evaluate('1','2'),lambda:field.evaluate('0.3','4.1'),
                 lambda:field.evaluate('0.3','2',log_tau='1'),
                 lambda:field.evaluate('0.3','2',scaled_time='1'),
                 lambda:field.cartesian('1','0','0',log_tau='-2'),
                 lambda:field.at_time('0','0','0','1')):
        try:call()
        except ValueError:rejected+=1
        else:raise AssertionError('Unsupported request accepted')
    code=dict(field.code_hashes)
    for path in (Path(__file__),ROOT/'src/openai_ns_reconstruction/original_core_velocity.py'):
        code[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(source_family=field.family,checks=dict(passed=True,
        exact_streamfunction_weight_divergence_and_swirl_identities=exact,
        native_physical_velocity_packets=len(native),scaled_time_packets=len(scaled),
        independent_returned_term_log_ratio_comparisons=len(ratios),
        directed_Cartesian_inverse_packets=len(inverse),unsupported_requests_rejected=rejected,
        exact_physical_time_packets=len(physical_time),
        exact_axis_included=True,offgrid_native_source_computed=True,
        leading_corrections_and_temporal_orders_remain_factored=True,
        ordinary_time_and_late_time_swirl_dominance_distinguished=True,
        original_absolute_F0_materialized=False),
        native_samples=native,scaled_time_samples=scaled,ratio_comparisons=ratios,
        Cartesian_samples=inverse,physical_time_samples=physical_time,input_hashes=field.hashes,implementation_hashes=code,
        limitations=dict(core_only=True,physical_velocity_is_signed_exponential_sum=True,
            ordinary_float_vector=False,useful_combined_physical_relative_error=False,
            full_temporal_series_admitted=False,global_matched_field=False,
            actual_collar_values_and_histories_evaluated=False,
            genuine_n_ge_2_temporal_recursion=False,full_corrected_NS_solution=False),
        execution_seconds=round(time.monotonic()-started,3))
    if args.write_report:(DIRECTORY/NAME).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds']),indent=2),flush=True)


if __name__=='__main__':main()
