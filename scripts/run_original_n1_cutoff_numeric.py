"""Evaluate actual selected n1 cutoff and Rin-to-rminus own histories."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_n1_cutoff_numeric import OriginalN1CutoffNumeric
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric.json'


def exact_checks():
    r,z,h,a,de=sy.symbols('rho Z hb phase delta',positive=True)
    sigma=sy.Function('sigma');f=sy.Function('f')(r,z);u=sy.Function('u')(r,z)
    arg=4*sy.log(r/4)/(h*a)-2;chi=1-sigma(arg)
    t=sy.Symbol('t');prime=sy.diff(sigma(t),t).subs(t,arg)
    assert sy.simplify(sy.diff(chi*f,r)-chi*sy.diff(f,r)+4*prime*f/(h*a*r))==0
    M=sy.Function('M')(r,z);mean=M/r;U=chi*u;d=1-z*z;L=1-de*z*z
    Q=(2*z*U-(1+de)*z*mean-d*sy.diff(mean,z))/L
    mr=(U-mean)/r
    expected=(2*z*sy.diff(U,r)-(1+de)*z*mr-d*sy.diff(mr,z))/L
    relation={sy.diff(M,r):U,sy.diff(M,r,z):sy.diff(U,z)}
    assert sy.simplify((sy.diff(Q,r)-expected).subs(relation).doit())==0
    div=L*(Q+r*sy.diff(Q,r))+d*sy.diff(U,z)-2*r*z*sy.diff(U,r)+(-1+de)*z*U
    assert sy.simplify(div.subs(relation).doit())==0
    S,eps,P2,A,F=sy.symbols('S epsilon Pstar_squared Phi f')
    assert sy.simplify((2*S*A*F/(eps**3*P2)).subs(S,eps**2*sy.Symbol('F0')**2)-2*sy.Symbol('F0')**2*A*F/(eps*P2))==0
    return 4


def compare(a,b):
    c=mp.mp.clone();c.dps=180
    diff=abs(c.mpf(a['approximation'])-c.mpf(b['approximation']))
    error=c.mpf(a['absolute_error_upper'])+c.mpf(b['absolute_error_upper'])
    assert diff<=error,(diff,error)


def all_rows(packet):
    for rows in packet['normalized_original_n1'].values():yield from rows
    for rows in packet['normalized_own_primitive'].values():
        if isinstance(rows,list):yield from rows
    for sectors in packet['native_first_radial_width_sectors'].values():
        for rows in sectors.values():yield from rows


def run():
    start=time.perf_counter();identities=exact_checks()
    field=OriginalN1CutoffNumeric();samples=[];lower=[]
    flat=seams=terminal=0
    for Z in ('-.3','0','.3'):
        for q in ('0','.5','.625','.75','1'):
            p=field.evaluate(Z,q)
            if q in ('.75','1'):
                for name in ('f','u'):
                    for row in p['normalized_original_n1'][name]:
                        assert row['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]];flat+=1
                    for rows in p['native_first_radial_width_sectors'][name+'_rho'].values():
                        for row in rows:
                            assert row['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]];flat+=1
            if q=='0':
                core=field._core[Z]['normalized_original_n1']
                for name,rows in p['normalized_original_n1'].items():
                    for new,old in zip(rows,core[name]):compare(new,old);seams+=1
                for row in p['normalized_own_primitive']['increment_per_original_hb_phase_in']:
                    assert row['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]]
            if q=='.625':assert p['cutoff']['value']['approximation']=='0.5'
            samples.append(p)
        at_zero=next(p for p in samples if p['Z']==Z and p['coordinate']=='.75')
        at_exit=next(p for p in samples if p['Z']==Z and p['coordinate']=='1')
        for q in ('1','1.5','2'):
            p=field.lower_collar(Z,q)
            assert p['normalized_own_primitive']==at_zero['normalized_own_primitive']
            assert p['cutoff_control_mass']['approximation']=='0.625';terminal+=1
            assert p['pressure_history']['leading_N1p_retained_once']
            if q=='1':
                for a,b in zip(all_rows(p),all_rows(at_exit)):compare(a,b);seams+=1
            lower.append(p)
        print('Actual cut n1 and retained lower history: Z='+Z,flush=True)
    example=next(p for p in lower if p['Z']=='.3' and p['coordinate']=='2')
    assert field.ctx.mpf(example['normalized_original_n1']['mean_u'][0]['approximation'])!=0
    assert field.ctx.mpf(example['normalized_original_n1']['q'][0]['approximation'])!=0
    transition=next(p for p in samples if p['Z']=='.3' and p['coordinate']=='.625')
    assert field.ctx.mpf(transition['native_first_radial_width_sectors']['u_rho']['-1,-1'][0]['approximation'])!=0
    coarse=OriginalN1CutoffNumeric(precision=60,control_cells=64);comparisons=0
    for chart,coordinate in (('inner','.625'),('lower','2')):
        old=coarse.evaluate('.3',coordinate) if chart=='inner' else coarse.lower_collar('.3',coordinate)
        new=transition if chart=='inner' else example
        for a,b in zip(all_rows(new),all_rows(old)):compare(a,b);comparisons+=1
    rejected=0
    for Z,q,chart in (('1','.5','inner'),('.3','-.1','inner'),('.3','1.01','inner'),
                      ('.3','.9','lower'),('.3','2.01','lower'),('.3','nan','inner')):
        try:field.evaluate(Z,q) if chart=='inner' else field.lower_collar(Z,q)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError(('Unsupported selected-cutoff query admitted',Z,q,chart))
    implementation=dict(field.code_hashes)
    implementation[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,checks=dict(passed=True,independent_cutoff_chain_own_recovery_and_pressure_scale_identities=identities,
        selected_inner_packets=len(samples),selected_Rin_to_rminus_packets=len(lower),
        exact_flat_zero_velocity_and_radial_coefficients=flat,core_and_inner_lower_seam_error_comparisons=seams,
        terminal_primitive_same_source_packets=terminal,cross_precision_coefficient_comparisons=comparisons,
        mean_and_Q_nonzero_after_velocity_zero=True,transition_inverse_hb_phase_in_derivative_nonzero=True,
        exact_original_sc_singleton_retained=True,original_raw_domain_not_extended=True,
        leading_pressure_forcing_not_cut_or_double_added=True,rejected_unsupported_requests=rejected),
        selected_inner_samples=samples,lower_collar_samples=lower,input_hashes=field.hashes,
        implementation_hashes=implementation,
        limitations=dict(actual_selected_inner_and_Rin_to_rminus_only=True,
            beyond_rminus_requires_actual_modified_leading_charts=True,
            physical_velocity_assembled=False,complete_five_moment_repair=False,
            exterior_pressure_m13_m15_global_tails=False,global_field=False,
            complete_Cartesian_C2_join=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-start,3))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
