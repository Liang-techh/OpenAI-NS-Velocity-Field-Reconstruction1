"""Original nonzero-mu quiet continuation from actual2Rc toRp."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_quiet_scalar import OriginalQuietScalar
from openai_ns_reconstruction.original_quiet_numeric import OriginalQuietVelocity
from run_original_second_micro_numeric import overlap

NAME = PREFIX+'compliant_current_original_quiet_numeric.json'


def independent_kernels():
    c = mp.mp.clone(); c.dps = 115
    iv = MPIntervalContext(); iv.dps = 85
    moments = corrections = 0
    for mu in ('.003', '.02'):
        f = OriginalQuietScalar(c, iv, c.mpf(mu), iv.mpf(mu))
        for text in ('0', '1e-120', '1e-8', '.2', '2', '8', '30'):
            s, si = c.mpf(text), iv.mpf(text)
            values = f.evaluate(si, iv)
            for n in range(1, 5):
                for forward in (False, True):
                    kernel = (lambda v: c.exp(-s*(1-v))*v**n) if forward else (lambda v: c.exp(-s*v)*v**n)
                    truth = s**(n+1)*c.quad(kernel, [0, c.mpf('.5'), 1])
                    got = values[('K' if forward else 'B')+str(n)]
                    lo, hi = _bounds(c, got)
                    assert lo <= truth <= hi, (mu, text, n, forward)
                    assert lo >= 0; moments += 1
            muv = c.mpf(mu)
            exprel=lambda z: c.expm1(z)/z if z else c.mpf(1)
            truth = dict(E_mu_divided=c.expm1(-muv*s)/muv,
                h_mu_divided=-s*s*c.quad(lambda v:c.exp(-s*(1-v))*v*exprel(-muv*s*v), [0,1]),
                e_mass_mu_divided=-2*s*s*c.quad(lambda v:v*exprel(-2*muv*s*v), [0,1]),
                p_mass_mu_divided=-2*s*s*c.quad(lambda v:c.exp(-s*v)*v*exprel(-2*muv*s*v), [0,1]))
            for name, value in truth.items():
                lo, hi = _bounds(c, values[name])
                assert lo <= value <= hi, (mu, text, name)
                corrections += 1
    # Long positive moments are checked independently through standard
    # incomplete-gamma/hypergeometric defining integrals at higher precision.
    f = OriginalQuietScalar(c, iv, c.mpf('.001'), iv.mpf('.001'))
    for text in ('1000', '1e20'):
        s = c.mpf(text)
        for n in range(1, 5):
            B = c.gammainc(n+1, 0, s)
            K = s**(n+1)*c.hyp1f1(n+1, n+2, s)*c.exp(-s)/(n+1)
            for forward, truth in ((False,B),(True,K)):
                lo, hi = _bounds(c, f.moments(iv.mpf(text),n,iv,forward=forward))
                assert lo <= truth <= hi, (text,n,forward)
                moments += 1
    return dict(independent_positive_own_rate_moments=moments,
                independent_nonzero_mu_divided_integrals=corrections,
                independent_precision=115, directed_precision=85)


def run():
    start = time.perf_counter(); independent = independent_kernels()
    print('Independent short/long positive moments and nonzero-mu divided integrals pass', flush=True)
    f = OriginalQuietVelocity(); c, iv = f.ctx, f.iv
    samples = []; seams = pressure = memory = nonzero = physical = scalar = frame = geometry = 0
    for Z in ('-.3', '0', '.3'):
        for q in ('0', '1e-160', '1e-20', '.25', '.5', '1'):
            p = f.evaluate(Z,q,log_tau='-2',azimuth='.25'); samples.append(p)
            inlet = f.numeric.inlets[Z]
            assert p['selected_integer'] == inlet['selected_integer']; memory += 1
            assert all(p['contract'].values()); memory += 1
            assert p['independent_P0_over_Pstar_squared'] == inlet['independent_P0_over_Pstar_squared']; memory += 1
            for name in ('m','h','k','e','p'):
                for row in p['finite_N_history_errors'][name]:
                    assert _bounds(c,f._restore(row)) == (0,0); memory += 1
            for n in range(1,5):
                for prefix in ('K','B'):
                    assert _bounds(c,f._restore(p['normalized_mu_profiles'][prefix+str(n)]))[0] >= 0; scalar += 1
            assert _bounds(c,f._restore(p['original_mu']['mu']))[0] > 0; nonzero += 1
            if q != '0':
                for name, sign in (('E',-1),('h',-1),('e',1),('p',-1)):
                    lo, hi = _bounds(c,f._restore(p['quiet_mu_divided_fields'][name][0]))
                    assert (hi<0 if sign<0 else lo>0),(Z,q,name); nonzero += 1
            else:
                for name in p['generic_leading_fields']:
                    for a,b in zip(p['generic_leading_fields'][name],inlet['generic_leading_fields'][name]):
                        overlap(a,b); seams += 1
                for source in p['quiet_mu_divided_fields'].values():
                    for row in source: assert _bounds(c,f._restore(row)) == (0,0); seams += 1
            terms = p['reserved_n1_pressure']['P1_C0_C1_terms']
            original = inlet['reserved_n1_pressure']['P1_C0_C1_terms']
            assert terms[:-1] == original; pressure += 1
            increment = terms[-1]
            assert increment['name'] == 'quiet_fixed_leading'; pressure += 1
            assert _bounds(c,f._restore(increment['coefficient_C0_C1'][1])) == (0,0); pressure += 1
            lo,hi = _bounds(c,f._restore(increment['coefficient_C0_C1'][0]))
            assert (lo==hi==0 if q=='0' else lo>0); pressure += 1
            assert len(p['physical']['pressure']) == 9; physical += 1
            assert p['reserved_n1_pressure']['limits']['beyond_2Rc']; pressure += 1
            assert p['physical']['coefficient_orders_included'] == [0,1]; physical += 1
            assert p['physical']['pressure_orders_included'] == [0,1]; physical += 1
            assert p['physical']['positive_mu_physical_terms_separate']; physical += 1
            assert p['physical']['leading_pressure_includes_same_P0_once']; physical += 1
            for name in ('ur','utheta','uz'):
                n1 = [v for v in p['physical']['cylindrical'][name] if v['temporal_order']==1]
                assert n1
                for v in n1: assert _bounds(c,f._restore(v['coefficient'])) == (0,0); physical += 1
            if Z!='0':
                assert _bounds(c,f._restore(p['generic_leading_fields']['m'][0])) != (0,0); memory += 1
            for name in ('h','e'):
                assert _bounds(c,f._restore(p['generic_leading_fields'][name][0])) != (0,0); memory += 1
            if q=='1':
                packet=p['native_pulse_input_frame']
                jets={n:f.numeric._jet(rows,True) for n,rows in packet.items()}
                for name,expected in (('m',jets['m1']*jets['u']),('k',jets['m2']*jets['u']*jets['u']),
                        ('h',jets['X']*jets['u']),('e',jets['energy']*jets['u']*jets['u']),
                        ('pressure',jets['P0']+jets['Mp'])):
                    actual=f.numeric._jet(p['generic_leading_fields'][name],True)
                    for a,b in zip(expected.coefficients,actual.coefficients):
                        lo,hi=_bounds(c,a-b);assert lo<=0<=hi;frame += 1
                terms=p['coordinate']['global_y_log_terms']
                finite=f._restore(terms[-2])+f._restore(terms[-1])
                exact=-60*f.numeric.outer.logmui-2
                lo,hi=_bounds(c,finite-exact);assert lo<=0<=hi;geometry += 1
        print('Actual quiet field reaches Rp with retained P1: Z='+Z,flush=True)
    scaled = 0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'1',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        normal=next(a for a in samples if a['Z']==Z and a['native_quiet_fraction']=='1')
        for a,b in zip(p['physical']['pressure'],normal['physical']['pressure']):
            overlap(a['coefficient'],b['coefficient']);scaled += 1
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled += 1
    rejected=0
    for Z,q in (('1','.5'),('nan','.5'),('.3','-.01'),('.3','1.01'),('.3','nan')):
        try:f.numeric.evaluate(Z,q)
        except ValueError:rejected += 1
        else:raise AssertionError('Invalid quiet request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,**independent,actual_original_quiet_physical_packets=len(samples),
            actual_2Rc_numerical_seam_conditions=seams,nonzero_mu_sector_conditions=nonzero,
            positive_quiet_scalar_conditions=scalar,nonzero_history_and_same_N_conditions=memory,
            incoming_and_integrated_P1_conditions=pressure,physical_orders_units_and_zero_n1_velocity_conditions=physical,
            Rp_native_pulse_frame_unit_conditions=frame,original_Rp_radius_recipe_conditions=geometry,
            scaled_time_conditions=scaled,rejected_requests=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(through_original_Rp=True,physical_velocity_pressure_orders=[0,1],leading_Z2=True,P1_ordinary_Z1=True,
            quiet_positive_mu_kept_separately=True,resolved_signed_P1_centers=False,relative_pressure_cancellation_accuracy=False,
            physical_gradients=False,global_n1_joins=False,Rp_selected_pulse_numerical_join=False,common_global_N=False,
            exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
