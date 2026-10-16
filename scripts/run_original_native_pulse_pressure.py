"""Focused defining-density and actual-source repaired pulse pressure checks."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_native_pulse_pressure import OriginalNativePulsePressureNumeric, OriginalNativePulsePressureVelocity
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from lei_ren_part1_paper_compliant_axial_pulse_field import gp, decay_integral, raw_beta
from run_original_second_micro_numeric import overlap

NAME = PREFIX+'compliant_current_original_native_pulse_pressure.json'


def independent_pressure():
    """Synthetic finite-mu check against the original local N1p density.

    The reference integrates Q_y/Q_yy directly, without radial integration
    by parts, endpoint tensors, the new W/WT reduction, or its cell sums.
    Two independent Gauss orders check reference convergence.
    """
    c = mp.mp.clone(); c.dps = 48
    iv = MPIntervalContext(); iv.dps = 85
    e = object.__new__(OriginalNativePulsePressureNumeric)
    e.ctx, e.iv, e.J, e.cells = c, iv, J, 128
    e.mu, e.delta, e.normal, e.ell = iv.mpf('.02'), iv.mpf('.001'), iv.mpf('.7'), iv.mpf('.15')
    e.a, e.b = iv.mpf('.5')-e.mu, iv.mpf('.5')+e.mu
    e._kernels, e._bumps, e._quadratics = {}, {}, {}
    e.gp, e.raw_beta = gp, raw_beta
    def data(Z):
        z = J.variable(iv, iv.mpf(Z), 4)
        C = [-c.mpf('.7')+z*c.mpf('.1'), c.mpf('.5')+z*c.mpf('.2')+z*z*c.mpf('.1')]
        return dict(selected_scaled_end_coefficient_Taylor=C), None, c.mpf('.9')/(1+z*z), None, None
    e.native = SimpleNamespace(_data=data, decay_integral=decay_integral)
    e.term = lambda packet, name, value, kind, extra=(): dict(name=name, jet=value, kind=kind, extra=extra)
    nodes = {}
    def rule(n):
        if n not in nodes: nodes[n] = c.gauss_quadrature(n, 'legendre')
        return nodes[n]
    def quad(f, left, right, n):
        if right <= left: return c.mpf(0)
        xs, ws = rule(n); h = (right-left)/2; m = (right+left)/2
        return h*sum((ws[k]*f(m+h*xs[k]) for k in range(n)), c.mpf(0))
    def beta(v, center):
        ell = c.mpf('.15'); r = (v-center)/ell
        if abs(r) >= 1: return c.mpf(0), c.mpf(0), c.mpf(0)
        f = c.exp(-1/(1-r*r))/(ell*c.mpf('.7'))
        d = -2*r/(ell*(1-r*r)**2)
        dd = -2/(ell**2*(1-r*r)**2)-8*r*r/(ell**2*(1-r*r)**3)
        return f, f*d, f*(d*d+dd)
    def reference(s, Z, n):
        a, b, de = c.mpf('.48'), c.mpf('.52'), c.mpf('.001')
        support_cuts = [c.mpf(center)+c.mpf('.15')*j/2 for center in (-3,-1) for j in range(-2,3)]
        cuts = sorted(set([s, c.mpf(0)]+[v for v in support_cuts if s < v < 0]))
        totals = [c.mpf(0), c.mpf(0)]
        xs, ws = rule(n)
        for lo, hi in zip(cuts, cuts[1:]):
            for k in range(n):
                v, weight = (lo+hi)/2+(hi-lo)*xs[k]/2, (hi-lo)*ws[k]/2
                beta_rows = [beta(v, center) for center in (-3, -1)]
                H = []
                for center in (-3,-1):
                    lower,upper=max(v,c.mpf(center)-c.mpf('.15')),c.mpf(center)+c.mpf('.15')
                    inner=sorted(set([lower,upper]+[c.mpf(center)+c.mpf('.15')*j/2 for j in (-1,0,1) if lower<c.mpf(center)+c.mpf('.15')*j/2<upper]))
                    H.append(-sum((quad(lambda w:c.exp(a*(w-v))*beta(w,center)[0],l,r,n) for l,r in zip(inner,inner[1:])),c.mpf(0)) if lower<upper else c.mpf(0))
                def density(z):
                    C = lambda w: (-c.mpf('.7')+c.mpf('.1')*w, c.mpf('.5')+c.mpf('.2')*w+c.mpf('.1')*w*w)
                    theta = lambda w: c.mpf('.9')/(1+w*w)*c.exp(-b*v)
                    U = lambda w: theta(w)*sum(C(w)[j]*beta_rows[j][0] for j in range(2))
                    Uy = lambda w: theta(w)*sum(C(w)[j]*(beta_rows[j][1]-b*beta_rows[j][0]) for j in range(2))
                    Uyy = lambda w: theta(w)*sum(C(w)[j]*(beta_rows[j][2]-2*b*beta_rows[j][1]+b*b*beta_rows[j][0]) for j in range(2))
                    mean = lambda w: theta(w)*sum(C(w)[j]*H[j] for j in range(2))
                    op = lambda f, m, w: (2*w*f(w)-(1-de)*w*m(w)-(1-w*w)*c.diff(m,w))/(1-de*w*w)
                    Q = lambda w: op(U, mean, w)
                    Qy = lambda w: op(Uy, lambda q: U(q)-mean(q), w)
                    Qyy = lambda w: op(Uyy, lambda q: Uy(q)-U(q)+mean(q), w)
                    q, qy, qz = Q(z), Qy(z), c.diff(Q,z)
                    L = 1-de*z*z
                    S = (q+(1-de)*z*qz/2+qy)/L+q*(q/2+qy)+U(z)*(-2*z*q+(1-z*z)*qz-2*z*qy)/L
                    return -c.exp(v)*S/2+qy+Qyy(z)
                totals[0] -= weight*density(Z)
                totals[1] -= weight*c.diff(density,Z)
        return totals
    checked = 0; convergence = []
    with mp.workdps(115):
        for text in ('-4','-3','-1'):
            s, z = c.mpf(text), c.mpf('.2')
            terms = e.end('.2', text, {})[0]
            jet = terms[0]['jet']+terms[1]['jet']+terms[2]['jet']*iv.exp(-e.b*iv.mpf(text))
            coarse, fine = reference(s,z,64), reference(s,z,96)
            for row in (0,1):
                discrepancy = abs(coarse[row]-fine[row]); convergence.append(c.nstr(discrepancy,12))
                assert discrepancy < c.mpf('2e-8'), (text,row,'reference quadrature convergence',discrepancy)
                lo,hi = _bounds(c,jet[row]); assert lo <= fine[row] <= hi,(text,row,lo,fine[row],hi)
                checked += 1
        # Independently integrate the new first-time moment, not its producer's cells.
        def sigma(x):
            if x <= 0: return c.mpf(0)
            if x >= 1: return c.mpf(1)
            return 1/(1+c.exp(1/x**2-1/(1-x)**2))
        def shape(x):
            if x <= 0 or x >= 11: return c.mpf(0)
            primitive = x-c.mpf('.01') if x >= c.mpf('.02') else quad(lambda w: sigma(50*w),0,x,72)
            return primitive*sigma(11-x)
        for text in ('1','25','550'):
            t = c.mpf(text); xi = c.mpf('.02')*t
            got = e.kernel_pair(iv.mpf(c.nstr(xi,48)),iv.mpf(text))
            cuts = sorted(set([c.mpf(0),t]+[(xi-x)/c.mpf('.02') for x in map(c.mpf,('0','.01','.02','10','10.5','11')) if 0 < (xi-x)/c.mpf('.02') < t]))
            for row in (0,1):
                truth = sum((quad(lambda v: v**row*c.exp(-c.mpf('.48')*v)*shape(xi-c.mpf('.02')*v),a,b,72) for a,b in zip(cuts,cuts[1:])),c.mpf(0))
                lo,hi = _bounds(c,got[row]);assert lo <= truth <= hi,(text,row,lo,truth,hi)
                checked += 1
    return dict(independent_defining_N1p_C0_C1=6, independent_positive_main_kernel_moments=6,
                total=checked, reference_Gauss_orders=[64,96], reference_convergence=convergence)


def run():
    start = time.perf_counter(); independent = independent_pressure()
    print('Independent defining N1p C0/C1 and first-moment kernels pass',flush=True)
    f = OriginalNativePulsePressureVelocity();c,iv = f.ctx,f.iv
    samples = [];conditions=seams=zeros=gap=units=0
    points = [('entrance','0'),('entrance','1'),('main','1'),('main','10.5'),('main','11'),
              ('gap','11'),('gap','12'),('gap_end','-1000'),('gap_end','-4'),
              ('end','-4'),('end','-3'),('end','-2'),('end','-1'),('end','-.5'),('end','0')]
    for Z in ('-.3','0','.3'):
        for chart,x in points:
            p=f.evaluate(Z,chart,x,log_tau='-2');samples.append(p);q=p['reserved_n1_pressure']
            assert p['physical']['pressure_orders_included']==[0,1];conditions+=1
            assert p['physical']['full_n1_pressure_included'];conditions+=1
            assert q['full_incoming_P1_at_Rp_retained']==p['incoming_P1_at_Rp'];conditions+=1
            assert q['terminal_pressure_value_zero_certified'];conditions+=1
            assert q['partial_quadratic_errors_not_full_signed_caps'];conditions+=1
            for source,physical in zip(q['P1_C0_C1_terms'],p['physical']['pressure'][1:]):
                assert len(source['coefficient_C0_C1'])==2;conditions+=1
                assert physical['temporal_order']==1;units+=1
                assert physical['log_scale_terms'][2:]==source['log_scale_terms'];units+=1
                assert physical['coefficient']==source['coefficient_C0_C1'][0];units+=1
            if chart=='entrance' and x=='0':
                assert q['P1_C0_C1_terms']==p['incoming_P1_at_Rp']['P1_C0_C1_terms'];seams+=1
            if chart=='end' and x in ('-.5','0'):
                for term in q['P1_C0_C1_terms']:
                    for row in term['coefficient_C0_C1']:
                        assert _bounds(c,f._restore(row))==(0,0);zeros+=1
            if chart in ('gap','gap_end'):
                assert q['representation']=='backward_terminal_and_gap_history';gap+=1
        main=next(p for p in samples if p['Z']==Z and p['native_pulse_chart']=='main' and p['native_pulse_coordinate']=='11')
        other=next(p for p in samples if p['Z']==Z and p['native_pulse_chart']=='gap' and p['native_pulse_coordinate']=='11')
        assert main['reserved_n1_pressure']==other['reserved_n1_pressure'];seams+=1
        e=next(p for p in samples if p['Z']==Z and p['native_pulse_chart']=='end' and p['native_pulse_coordinate']=='-4')
        g=next(p for p in samples if p['Z']==Z and p['native_pulse_chart']=='gap_end' and p['native_pulse_coordinate']=='-4')
        for a,b in zip(e['reserved_n1_pressure']['P1_C0_C1_terms'],g['reserved_n1_pressure']['P1_C0_C1_terms']):
            assert a==b;seams+=1
        print('Actual repaired pulse P1 reaches Rv: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'end','-3',log_tau=None,scaled_time='1e18');samples.append(p)
        regular=next(v for v in samples if v['Z']==Z and v['native_pulse_chart']=='end' and v['native_pulse_coordinate']=='-3')
        assert p['reserved_n1_pressure']==regular['reserved_n1_pressure'];scaled+=1
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent=independent,actual_pressure_packets=len(samples),source_and_C1_conditions=conditions,
            actual_Rp_main_gap_end_source_seams=seams,source_certified_terminal_and_empty_support_zeros=zeros,
            nonzero_gap_history_representations=gap,physical_lambda_minus2_units=units,scaled_time_conditions=scaled),
        execution_seconds=round(time.perf_counter()-start,3),limits=dict(original_selected_pulse_P1_C1=True,
            actual_incoming_and_terminal_theorem_both_consumed=True,physical_velocity_pressure_orders=[0,1],
            signed_point_centers_resolved=False,relative_pressure_cancellation_accuracy=False,m15_numeric=False,
            global_physical_gradients=False,common_global_N=False,exact_heat_exterior=False,stress_cone=False,
            genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-report',action='store_true');p.add_argument('--inspect-report',action='store_true')
    args=p.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
