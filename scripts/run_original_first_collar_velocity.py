"""Compute original first-collar physical velocity and inlet seam evidence."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_first_collar_velocity import OriginalFirstCollarVelocity
from openai_ns_reconstruction.original_core_velocity import OriginalCoreVelocity
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_first_collar_velocity.json'


def exact_physical_identities():
    R,Z,lam,de,beta,rho,h,s=sy.symbols('R Z lambda delta beta rho h s',positive=True)
    mass=sy.Function('M')(R,Z);d,L=1-Z*Z,1-de*Z*Z
    mr,mz=sy.diff(mass,R),sy.diff(mass,Z)
    numerator=2*Z*R*mr-beta*Z*mass-d*mz
    divr=sy.diff(numerator,R)
    divz=(beta-2)*Z*mr+d*sy.diff(mr,Z)-2*Z*R*sy.diff(mr,R)
    assert sy.simplify(divr+divz)==0
    count=1
    for weight in (1-de,1+de):
        assert sy.simplify((divr+divz).subs(beta,weight))==0;count+=1
    for base,power in ((sy.sqrt(rho/2),sy.Rational(1,2)),
                       (sy.sqrt(2*rho),sy.Rational(1,2)),
                       (rho,1)):
        assert sy.simplify(base.subs(rho,4*sy.exp(h*s))/base.subs(rho,4)-sy.exp(power*h*s))==0
        count+=1
    ur,ut,angle=sy.symbols('ur utheta theta',real=True)
    u=sy.cos(angle)*ur-sy.sin(angle)*ut
    v=sy.sin(angle)*ur+sy.cos(angle)*ut
    assert sy.trigsimp(u*u+v*v-ur*ur-ut*ut)==0
    return count+1


def overlap(field,a,b):
    lo,hi=_bounds(field.ctx,field._restore(a));left,right=_bounds(field.ctx,field._restore(b))
    assert max(lo,left)<=min(hi,right)


def structural_geometry_checks(packet):
    count=0
    for name,rows in packet['cylindrical'].items():
        for row in rows:
            geometry=[q for q in row['log_scale_terms'] if q.get('kind')=='original_hb_linear_log']
            assert len(geometry)==(0 if name=='uz' else 1)
            if geometry:assert geometry[0]['multiplier']=='0.5'
            widths=[q for q in row['log_scale_terms'] if q.get('kind')=='original_log_hb']
            assert len(widths)==(1 if row['name'].startswith('leading_') and row['name'].endswith('_collar_increment') else 0)
            count+=1
    for row in packet['streamfunction']:
        geometry=[q for q in row['log_scale_terms'] if q.get('kind')=='original_hb_linear_log']
        assert len(geometry)==1 and geometry[0]['multiplier']=='1.0';count+=1
    assert not any(q.get('kind')=='original_hb_linear_log' for q in packet['coordinates']['z']['log_scale_terms'])
    assert packet['native']['rho_definition']=='4*exp(hb*phase)'
    assert packet['native']['original_log_hb']['source_value_selected'] is False
    assert packet['direct_Cartesian_inverse_available'] is False
    return count


def returned_ratio(field,packet):
    """Recover the ratio from emitted physical terms, cancelling common logs."""
    c,iv=field.ctx,field.iv
    rows={row['name']:row for row in packet['cylindrical']['utheta']}
    old,new=rows['leading_swirl'],rows['n1_swirl']
    phi=field._restore(old['coefficient'])
    increment=field._restore(rows['leading_swirl_collar_increment']['coefficient'])
    if _bounds(c,abs(increment))[1]:
        hlog=field._restore(packet['native']['original_log_hb'])
        bound=hlog+iv.ln(abs(increment))
        assert _bounds(c,bound)[1]<_bounds(c,iv.ln(iv.mpf('1e-1000')))[0]
        phi+=iv.mpf(['-1e-1000','1e-1000'])
    f=field._restore(new['coefficient'])+field._restore(rows['n1_swirl_collar_increment']['coefficient'])
    truth=iv.ln(abs(f)/phi)
    key=lambda row:json.dumps(row,sort_keys=True)
    a,b=Counter(map(key,old['log_scale_terms'])),Counter(map(key,new['log_scale_terms']))
    common=a&b;a-=common;b-=common
    for sign,terms in ((-1,a),(1,b)):
        for text,n in terms.items():
            row=json.loads(text)
            assert 'approximation' in row,'Original huge amplitude and width geometry must cancel'
            truth+=sign*n*field._restore(row)
    saved=field._restore(packet['swirl_order_comparison']['log_absolute_ratio'])
    lo,hi=_bounds(c,truth);left,right=_bounds(c,saved)
    assert max(lo,left)<=min(hi,right)
    return dict(Z=packet['native']['Z'],phase=packet['native']['phase'],
        classification=packet['swirl_order_comparison']['classification'],
        original_F0_and_width_geometry_cancelled_from_returned_terms=True)


def run():
    started=time.perf_counter();identities=exact_physical_identities()
    field=OriginalFirstCollarVelocity()
    native,scaled,ratios=[],[],[]
    geometry=0
    for Z in ('-.3','0','.3'):
        for phase in ('0','.25','.5'):
            packet=field.evaluate(Z,phase,log_tau='-2',azimuth='.25')
            assert {name:len(rows) for name,rows in packet['Cartesian'].items()}==dict(u=9,v=9,w=5)
            geometry+=structural_geometry_checks(packet)
            native.append(packet);ratios.append(returned_ratio(field,packet))
        later=field.evaluate(Z,'.5',log_tau=None,scaled_time='1e18',azimuth='1.7')
        assert later['swirl_order_comparison']['classification']=='n1_smaller'
        scaled.append(later);ratios.append(returned_ratio(field,later))
        print('Original physical first-collar values/time: Z='+Z,flush=True)
    offgrid=field.evaluate('.17','.37',log_tau='-3',azimuth='3.2')
    geometry+=structural_geometry_checks(offgrid)
    ratios.append(returned_ratio(field,offgrid));native.append(offgrid)
    # Independent accepted core entry point: compare the entire phase-zero seam.
    core=OriginalCoreVelocity(precision=85)
    seams=0
    for Z in ('-.3','0','.3'):
        old=core.evaluate(Z,'4',log_tau='-2',azimuth='.25')
        new=next(q for q in native if q['native']['Z']==Z and q['native']['phase']=='0')
        groups=[(old['Cartesian'][name],new['Cartesian'][name]) for name in ('u','v','w')]
        groups.append((old['streamfunction'],new['streamfunction']))
        for before,after in groups:
            current={row['name']:row for row in after}
            for row in before:
                overlap(field,row['coefficient'],current[row['name']]['coefficient']);seams+=1
            for row in after:
                if row['name'].endswith('_collar_increment'):
                    assert row['coefficient']['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]]
        for name in ('r','x','y','z'):
            overlap(field,old['coordinates'][name]['coefficient'],new['coordinates'][name]['coefficient']);seams+=1
        for row in new['cylindrical']['ur']+new['cylindrical']['utheta']+new['streamfunction']:
            for term in row['log_scale_terms']:
                if term.get('kind')=='original_hb_linear_log':
                    assert term['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]]
    rejected=0
    for args in (('1','.5',{}),('.3','.50001',{}),('.3','-.1',{}),
                 ('.3','.5',dict(log_tau='1')),('.3','.5',dict(scaled_time='1')),
                 ('.3','.5',dict(azimuth='nan'))):
        try:field.evaluate(args[0],args[1],**args[2])
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError(('Unsupported collar request admitted',args))
    code=dict(field.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,
        checks=dict(passed=True,independent_physical_geometry_streamfunction_and_rotation_identities=identities,
            native_physical_collar_packets=len(native),source_scaled_time_packets=len(scaled),
            physical_term_geometry_and_width_factor_checks=geometry,
            independent_returned_term_order_ratios=len(ratios),
            independently_evaluated_core_seam_coefficient_overlaps=seams,
            original_zero_phase_increments_and_geometry_exact=True,
            offgrid_original_phase_computed=True,unsupported_requests_rejected=rejected,
            actual_positive_width_not_chosen_from_its_enclosure=True,
            original_hb_not_applied_twice_to_n1_increments=True),
        native_samples=native,scaled_time_samples=scaled,ratio_checks=ratios,
        input_hashes=field.hashes,implementation_hashes=code,
        limitations=dict(native_phase_required=True,first_phase_max='.5',
            signed_exponential_sum_with_source_log_enclosures=True,ordinary_float_vector=False,
            direct_Cartesian_collar_inverse=False,useful_combined_physical_relative_error=False,
            numerical_spatial_derivatives_or_join_order_one_certified=False,
            global_matched_field=False,global_moment_repair=False,heat_stress_layers=False,
            temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-started,3))


def inspect(path):
    report=json.loads(path.read_text(encoding='utf-8'));assert report['checks']['passed']
    for name,digest in report['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    report=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds']),indent=2),flush=True)
