"""Focused independent compact/quiet weighted m15 transport evidence."""
import argparse
import ast
import hashlib
import inspect as inspectlib
import json
from pathlib import Path
import sys
import textwrap
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'experiments/root_st073')]
DIRECTORY = ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_preRp_m15 import OriginalPreRpM15Numeric
from openai_ns_reconstruction.original_compact_n1_pressure import OriginalReservedN1PressureNumeric
from openai_ns_reconstruction.original_reserved_n1_numeric import OriginalReservedN1Numeric
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets

NAME = PREFIX+'compliant_current_original_preRp_m15_numeric.json'


def independent_units_and_baseline():
    import sympy as sy
    import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as original
    lc, lp, lu = sy.symbols('logRc logP logUp', real=True)
    lm = -4*lp-sy.log(1000)
    ctx = SimpleNamespace(mpf=lambda v: sy.Rational(str(v)), log=sy.log, ln=sy.log,
                          exp=sy.exp, sqrt=sy.sqrt, pi=sy.pi)
    namespace = dict(original.pressure_unit_log.__globals__)
    namespace['interval'] = lambda c, v: v
    tree = ast.parse(textwrap.dedent(inspectlib.getsource(original.pressure_unit_log)))
    exec(compile(tree, '<original p=2 unit>', 'exec'), namespace)
    original_unit = namespace['pressure_unit_log'](ctx, dict(logP=lp, logUp=lu, logRc=lc), 2)
    background = SimpleNamespace(patch=SimpleNamespace(leading=SimpleNamespace(
        defect=SimpleNamespace(logP=lp), reference=SimpleNamespace(_summary=lambda p, v: v))),
        outer=SimpleNamespace(logmui=lm))
    reserved = object.__new__(OriginalReservedN1Numeric)
    reserved.iv=ctx; reserved.ctx=ctx; reserved.background=background; reserved.logUp=lu
    reserved._source_row=lambda v: v; reserved._mid=lambda v: v
    reader = object.__new__(OriginalReservedN1PressureNumeric)
    reader.iv=ctx; reader.background=background; reader.reserved=reserved
    reader._row=lambda v: v; reader._rows=lambda v: v
    f = object.__new__(OriginalPreRpM15Numeric)
    f.iv=ctx; f.reader=reader; f.reserved=reserved; f.logk=sy.Symbol('logk', real=True)
    packet = dict(coordinate=dict(global_y_log_terms=[lc-sy.log(110), 0, 0, 0, sy.Symbol('local_s')]))
    actual = sum(f.fixed_unit_atoms(packet))
    assert sy.simplify(sy.expand(actual-original_unit)) == 0
    sigma = -lp/2+sy.Rational(1, 2)
    own_unit = sum(reserved.beta_unit_atoms(packet))+lp+sigma+lc
    assert sy.simplify(sy.expand(own_unit-actual)) == 0
    assert sy.simplify(sum(f.leading_unit_atoms(packet))-(2*lc+2*f.logk-sy.log(4))) == 0
    R, Rc, K, L = sy.symbols('R Rc K L', positive=True)
    q = -K/R
    # Integrate the original R-weighted radial-free source independently.
    boundary_linear = R**2*q-Rc**2*q.subs(R, Rc)
    integral_linear = sy.integrate(R*q, (R, Rc, R))
    boundary_square = R**2*q*q-Rc**2*q.subs(R, Rc)**2
    integral_square = sy.integrate(R*q*q, (R, Rc, R))
    S1 = (boundary_linear-integral_linear)/L+(boundary_square-integral_square)/2
    viscous = R*(R*sy.diff(q, R))-Rc*(R*sy.diff(q, R)).subs(R, Rc)
    primitive = sy.simplify(S1/2-viscous)
    assert sy.simplify(primitive+K*K*sy.log(R/Rc)/4) == 0
    assert sy.simplify(sy.diff(primitive, R)+R*q*q/4) == 0
    assert sy.simplify(viscous) == 0
    return dict(actual_p2_unit_recipe=True, actual_Utheta1_to_p2_unit_cancellation=True,
        actual_leading_K_squared_unit=True, original_S1_weighted_baseline=True,
        direct_weighted_density_derivative=True, both_baseline_RQ_y_endpoints_cancel=True)


def independent_physical_swirl():
    c = mp.mp.clone(); c.dps=65
    iv = MPIntervalContext(); iv.dps=85
    f = object.__new__(OriginalPreRpM15Numeric)
    f.ctx=c; f.iv=iv; f.J=J; f.cells=256; f.ell=iv.log(2)/400; f.beta_jets=beta_jets
    f._pressure_integrals={}
    normalp = c.quad(lambda t: c.exp(-1/(1-t*t)+c.log(2)*t/800)
                     if abs(t)<1 else c.mpf(0), [-1, 0, 1])
    normal = iv.mpf(normalp)
    # A high-precision independent normalizer is widened beyond its rounding error.
    nlo,nhi=_bounds(c,normal); normal=iv.mpf([nlo-c.mpf('1e-60'),nhi+c.mpf('1e-60')])
    centers=[iv.log(2)*j/40 for j in (12,13,14,15,16)]
    f.reserved=SimpleNamespace(centers=centers, normal=normal, _decode=lambda v:v)
    f.reader=SimpleNamespace(_hull=lambda vs:iv.mpf([min(_bounds(c,v)[0] for v in vs),max(_bounds(c,v)[1] for v in vs)]))
    count=prefixes=negative=0
    Rc, Pstar, Sigma, unit = map(c.mpf, ('3','5','7','11'))
    for mu in ('.001','.2','.4'):
        f.mu=iv.mpf(mu)
        def arc(z):return c.exp(-c.mpf('1.7')-c.mpf('2.5')*c.mpf(mu))/(1+z*z)
        def bs(z):return [1+z, -c.mpf('.3')+z*z, c.mpf('.2')-z**3]
        z=c.mpf('.23'); a=arc(z); az=c.diff(arc,z)
        bv=bs(z); bz=[c.diff(lambda zz:bs(zz)[j],z) for j in range(3)]
        cell=dict(actual_A_rc_C0_C1_coefficient_enclosures=[iv.mpf(a),iv.mpf(az)],
            coefficients={name:dict(normalized_C0_C1_enclosures=[iv.mpf(bv[j]),iv.mpf(bz[j])])
                          for j,name in enumerate(('beta1','beta2','beta3'))})
        for fraction in ('.31','.35','.39'):
            s=c.log(2)*c.mpf(fraction); got, factors=f.own_swirl([cell],iv.mpf(s))
            moments=[]
            for j in (1,2,3):
                center=c.log(2)*(12+j)/40; ell=c.log(2)/400
                left,right=Rc*c.exp(center-ell),min(Rc*c.exp(s),Rc*c.exp(center+ell))
                def density(r):
                    u=c.log(r/Rc); t=(u-center)/ell
                    if abs(t)>=1:return c.mpf(0)
                    shape=c.exp(-1/(1-t*t))/(ell*c.exp(center/2)*normalp)
                    theta0=Pstar*Sigma*c.exp(-(c.mpf('.5')+c.mpf(mu))*u)
                    theta1=unit/(Pstar*Sigma*Rc)*shape
                    return theta0*theta1/unit
                moment=c.mpf(0) if right<=left else c.quad(density,[left,(left+right)/2,right])
                moments.append(moment)
                lo,hi=_bounds(c,factors[j-1]); assert lo<=moment<=hi,(mu,fraction,j,'physical prefix')
                prefixes+=1
            def reference(zz):return -arc(zz)*sum(b*m for b,m in zip(bs(zz),moments))
            truth=[reference(z),c.diff(reference,z)]
            for k in (0,1):
                lo,hi=_bounds(c,got[k]); assert lo<=truth[k]<=hi,(mu,fraction,k,lo,truth[k],hi)
                count+=1
            if fraction=='.39':assert _bounds(c,got[0])[1]<0;negative+=1
    return dict(direct_physical_dR_swirl_C0_C1_comparisons=count,
                independent_positive_weighted_prefixes=prefixes,nonzero_negative_complete_swirl_cases=negative,
                finite_mu_values=['.001','.2','.4'])


def run():
    start=time.perf_counter()
    algebra=independent_units_and_baseline(); direct=independent_physical_swirl()
    print('Independent original p=2 units, S1 baseline and signed physical swirl C0/C1 pass',flush=True)
    f=OriginalPreRpM15Numeric();c=f.ctx;iv=f.iv
    def restore(row):return f.reserved._source_row(row)
    samples=[];source=finite=history=zeros=seams=0
    points=[('compact','1'),('compact',c.nstr(c.exp(c.log(2)*c.mpf('.31')),95))]
    for j in (13,14,15):
        points.append(('compact',c.nstr(c.exp(c.log(2)*c.mpf(j)/40),95)))
        points.append(('compact',c.nstr(c.exp(c.log(2)*(c.mpf(j)/40+c.mpf(1)/800)),95)))
    points += [('compact','1.9'),('compact','2'),('quiet','0'),('quiet','.5'),('quiet','1')]
    for Z in ('-.3','0','.3'):
        for chart,x in points:
            p=f.evaluate(Z,chart,x); samples.append(p)
            assert all(p['contract'].values());source+=len(p['contract'])
            assert p['selected_integer']==f.quiet.inlets[Z]['selected_integer'];source+=1
            assert not p['limits']['independent_signed_numeric_seam_residual'];source+=1
            terms={t['name']:t for t in p['m15_C0_C1_terms']}
            assert set(terms)=={'inner','lower','active','compact_finite_N','compact_exact_leading','own_swirl'}|({'quiet_exact_leading'} if chart=='quiet' else set());source+=1
            for order in p['actual_partial_finite_N_normalized_log_bounds']:
                for row in order:assert _bounds(c,restore(row))[1]<-1000;finite+=1
            for name in ('compact_exact_leading','quiet_exact_leading'):
                if name not in terms:continue
                row=terms[name]['coefficient_C0_C1'];assert _bounds(c,restore(row[1]))==(0,0);history+=1
                s=(iv.log(iv.mpf(x)) if chart=='compact' else iv.log(2)) if name=='compact_exact_leading' else f.quiet.lengthi*iv.mpf(x)
                lo,hi=_bounds(c,restore(row[0]));ll,hh=_bounds(c,-s)
                assert max(lo,ll)<=min(hi,hh);history+=1
                if _bounds(c,s)[0]>0:assert hi<0;history+=1
            if chart=='compact' and x=='1':
                for name in ('compact_finite_N','compact_exact_leading','own_swirl'):
                    for row in terms[name]['coefficient_C0_C1']:assert _bounds(c,restore(row))==(0,0);zeros+=1
            for row in p['swirl_partial_integral_factors']:
                assert _bounds(c,restore(row))[0]>=0;source+=1
            assert p['incoming_true_function']==f.reader.binding['incoming']['m15'];source+=1
        compact=f.evaluate(Z,'compact','2');quiet=f.evaluate(Z,'quiet','0')
        assert compact['m15_C0_C1_terms']==quiet['m15_C0_C1_terms'][:-1];seams+=1
        assert all(_bounds(c,restore(v))==(0,0) for v in quiet['m15_C0_C1_terms'][-1]['coefficient_C0_C1']);seams+=1
        rp=f.evaluate(Z,'quiet','1')
        assert rp['source_partition_identity']['actual_Rp_source_equals_negative_complete_pulse'];seams+=1
        print('Independent repaired weighted histories reach Rp: Z='+Z,flush=True)
    rejected=0
    for args in [('1','compact','1'),('nan','quiet','0'),('.3','compact','0'),('.3','compact','2.1'),('.3','quiet','-1'),('.3','quiet','1.1'),('.3','unknown','0')]:
        try:f.evaluate(*args)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original preRp m15 query admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_source_algebra=algebra,independent_physical_swirl=direct,
            actual_compact_quiet_m15_packets=len(samples),source_and_scope_conditions=source,
            partial_finite_N_bounds=finite,exact_nonzero_K_squared_history_conditions=history,
            inlet_zero_prefix_conditions=zeros,compact_quiet_seams_and_source_partition_conditions=seams,
            rejected_requests=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(independent_preRp_m15_C0_C1_enclosure_transport=True,preRp_pulse_signed_source_identity=True,
            resolved_signed_point_centers=False,independent_signed_numeric_seam_residual=False,
            relative_cancellation_accuracy=False,global_five_moment_numeric_join=False,common_global_N=False,
            exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    data=json.loads(path.read_text(encoding='utf-8'));assert data['checks']['passed']
    for name,digest in data['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in data['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    data=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
