"""Original compact background, supported kernels and source-relative closure."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_compact_numeric import OriginalCompactNumeric,supported_mass
from openai_ns_reconstruction.original_compact_background_velocity import OriginalCompactBackgroundVelocity
from openai_ns_reconstruction.original_o3_scalar import OriginalO3Scalar
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_compact_numeric.json'


def independent_kernels():
    c=mp.mp.clone();c.dps=70;iv=MPIntervalContext();iv.dps=90
    raw=[('.2','.3'),('.6','.8')];supports=[(iv.mpf(a),iv.mpf(b)) for a,b in raw];count=0
    for s in ('0','.15','.25','.5','.7','.9','1'):
        x=c.mpf(s)
        for rate in ('0','1','1.5'):
            r=c.mpf(rate)
            for direction in ('forward','backward'):
                got,_=supported_mass(iv,iv.mpf(s),iv.mpf(1),rate,direction,supports);truth=c.mpf(0)
                for a,b in raw:
                    left,right=(c.mpf(a),min(x,c.mpf(b))) if direction=='forward' else (max(x,c.mpf(a)),c.mpf(b))
                    if right>left:truth+=c.quad(lambda u:c.exp(-r*(x-u)),[left,right])
                lo,hi=_bounds(c,got);assert lo-c.mpf('1e-60')<=truth<=hi+c.mpf('1e-60'),(s,rate,direction);count+=1
    # A time box must cover every true clipped mass, including support crossings.
    for lo,hi in (('.19','.21'),('.29','.31'),('.79','.81')):
        for direction in ('forward','backward'):
            got,_=supported_mass(iv,iv.mpf([lo,hi]),iv.mpf(1),'1',direction,supports)
            for point in (lo,hi):
                truth,_=supported_mass(iv,iv.mpf(point),iv.mpf(1),'1',direction,supports)
                a,b=_bounds(c,got);l,h=_bounds(c,truth);assert a<=l+c.mpf('1e-60') and b>=h-c.mpf('1e-60');count+=1
    return count


def independent_mu_continuation():
    c=mp.mp.clone();c.dps=70;iv=MPIntervalContext();iv.dps=90;mu=c.mpf('.001')
    f=object.__new__(OriginalCompactNumeric);f.ctx,f.iv=c,iv
    f.outer=SimpleNamespace(scalar=OriginalO3Scalar(c,iv,mu,iv.mpf('.001')),mu=mu,mui=iv.mpf('.001'))
    restore=lambda row:iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])
    f.patch=SimpleNamespace(leading=SimpleNamespace(reference=SimpleNamespace(previous=SimpleNamespace(_restore=restore))))
    initial=dict(J=c.mpf('2.5'),U1=c.expm1(-mu*c.mpf('2.5'))/mu,U2=c.expm1(-2*mu*c.mpf('2.5'))/mu,C1=c.mpf('-.3'),C0=c.mpf('-.2'),Cm=c.mpf('-.1'))
    def row(v):
        t=iv.mpf(v)
        return dict(approximation=c.nstr(v,70),source_enclosure_exact_mpf_tuples=[list(q) for q in t._mpi_])
    inlet=dict(normalized_mu_profiles={n:row(v) for n,v in initial.items()});count=0
    for s in ('0','.1','.69'):
        x=c.mpf(s);got=f._scalar(inlet,x,iv.mpf(s))
        for n,alpha,factor in (('C1',1,1),('C0',0,2),('Cm',-1,2)):
            truth=initial[n]+c.exp(3*alpha)*c.quad(lambda u:c.exp(alpha*u)*c.expm1(-factor*mu*(c.mpf('2.5')+u))/mu,[0,x])
            lo,hi=_bounds(c,got['directed'][n]);assert lo-c.mpf('1e-60')<=truth<=hi+c.mpf('1e-60'),(s,n);count+=1
    return count


def run():
    start=time.perf_counter();kernels=independent_kernels();mu=independent_mu_continuation()
    print('Independent supported forward/backward kernels and finite-mu continuation pass',flush=True)
    f=OriginalCompactBackgroundVelocity();samples=[];seams=pressure=history=terminal=support=physical=0
    old=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_outer_o3_numeric.json')).read_text(encoding='utf-8'))
    for Z in ('-.3','0','.3','.17','.23'):
        for x in ('1','1.13','1.15','1.4','1.42','1.74','1.9','2'):
            p=f.evaluate(Z,x,log_tau='-2',azimuth='.25');samples.append(p);inlet=f.numeric._inlets[Z]
            assert p['selected_integer']==inlet['selected_integer'];history+=1
            assert not p['limits']['original_control_signed_values_resolved'];history+=1
            assert p['independent_P0_over_Pstar_squared']==inlet['independent_P0_over_Pstar_squared'];pressure+=1
            for k in range(3):
                v=f._restore(p['generic_leading_fields']['pressure'][k])-f._restore(p['generic_leading_fields']['p'][k])-f._restore(p['independent_P0_over_Pstar_squared'][k])
                lo,hi=_bounds(f.ctx,v);assert lo<=0<=hi;pressure+=1
            assert p['physical']['coefficient_orders_included']==[0] and not p['physical']['compact_repaired_n1_field_included'];physical+=1
            assert p['physical']['preRc_n1_homogeneous_history_not_extrapolated'];physical+=1
            assert _bounds(f.ctx,f._restore(p['original_mu']['mu']))[0]>0;physical+=1
            if x=='1':
                for group in ('generic_modified_fields','physical_similarity','finite_N_history_errors'):
                    for name,rows in p[group].items():
                        for a,b in zip(rows,inlet[group][name]):overlap(a,b);seams+=1
                previous=next(a for a in old['samples'] if a['Z']==Z and a['native_O3_chart']=='O3_power' and a['native_O3']=='2' and a['physical']['time']['kind']=='log_tau')
                for name,terms in p['physical']['cylindrical'].items():
                    overlap(terms[0]['coefficient'],previous['physical']['cylindrical'][name][0]['coefficient']);physical+=1
            if x in ('1.9','2'):
                assert p['compact_source_cover']['active_bumps']==[];support+=1
                for n,rows in p['finite_N_history_errors'].items():
                    for row in rows:assert _bounds(f.ctx,f._restore(row))==(0,0);terminal+=1
                    assert p['relative_terminal_source_identities'][n]['relative_only'];terminal+=1
                # Full fields retain their nonzero leading memories and P0.
                for n in ('E','h','p'):
                    assert _bounds(f.ctx,f._restore(p['generic_leading_fields'][n][0]))[0]>0;terminal+=1
            if x=='1.15':
                assert p['compact_source_cover']['active_bumps']==[0];support+=1
                for rows in p['finite_N_local_errors'].values():
                    assert any(_bounds(f.ctx,f._restore(row))[0]<_bounds(f.ctx,f._restore(row))[1] for row in rows);support+=1
                for n in ('h','p'):
                    assert any(_bounds(f.ctx,f._restore(row))[0]<_bounds(f.ctx,f._restore(row))[1] for row in p['finite_N_history_errors'][n]);history+=1
        print('Actual original compact background and correlated histories: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'2',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        original=next(a for a in samples if a['Z']==Z and a['native_compact_x']=='2')
        for n,terms in p['physical']['cylindrical'].items():overlap(terms[0]['coefficient'],original['physical']['cylindrical'][n][0]['coefficient']);scaled+=1
    rejected=0
    for x,kwargs in [('0.9',{}),('2.1',{}),('nan',{}),('1.5',dict(log_tau='1')),('1.5',dict(azimuth='nan'))]:
        try:f.evaluate('.3',x,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid compact request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_supported_kernel_enclosures=kernels,independent_finite_mu_continuation_values=mu,
            actual_original_compact_background_packets=len(samples),Rc_numerical_seam_coefficients=seams,
            original_N_and_nonzero_history_conditions=history,P0_unit_and_once_only_checks=pressure,
            source_identity_relative_zero_and_nonzero_full_field_conditions=terminal,support_conditions=support,
            physical_background_only_and_Rc_seam_conditions=physical,scaled_time_checks=scaled,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),limits=dict(through_original_2Rc_background=True,physical_temporal_orders=[0],
            resolved_signed_finite_N_centers=False,original_control_signed_values_resolved=False,compact_n1_numeric_fields=False,
            relative_terminal_zero_from_checked_source_identity=True,independent_signed_control_residual_evaluation=False,
            total_radial_or_fast_phase_derivatives=False,full_n1_pressure=False,reduced_m13_m15_integrals_continued=False,
            relative_pressure_cancellation_accuracy=False,admissible_stress_cone_margin=False,global_N_admission=False,
            global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
