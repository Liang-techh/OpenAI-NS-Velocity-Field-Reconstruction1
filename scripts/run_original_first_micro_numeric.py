"""Compute original selected-N first-micro velocity enclosures and prefixes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_first_micro_numeric import OriginalFirstMicroNumeric
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_first_micro_numeric.json'


def independent_prefix_and_exponential_checks():
    c=mp.mp.clone();c.dps=90;iv=MPIntervalContext();iv.dps=130
    fake=object.__new__(OriginalFirstMicroNumeric)
    fake.transport=SimpleNamespace(ctx=c,iv=iv,phase_iv=iv.mpf(1)/32)
    fake.leading=SimpleNamespace(core=SimpleNamespace(lp_iv=iv.mpf(0),ll_iv=iv.mpf(0)),logh=-iv.ln(7))
    fake.logN=iv.ln(17)
    def cap(v):
        if not v:return dict(exact_zero=True)
        q=iv.ln(v)
        return dict(exact_zero=False,log_absolute_upper=dict(lower_exact_mpf_tuple=list(q._mpi_[0]),upper_exact_mpf_tuple=list(q._mpi_[1])))
    caprows=[cap(v) for v in (3,3,2,0)]
    cells=[]
    for left,right in ((['sc_half'],[1,4]),([1,4],[1,2]),([1,2],[3,4]),([3,4],[1,1])):
        cells.append(dict(exact_left=left,exact_right=right,
            actual_source_and_primitive_range=dict(A=caprows,B=caprows),
            actual_density_ranges=dict(N_scaled={k:caprows for k in ('m','h','k','e','p')})))
    fake.ranges=[dict(exact_Z_cell=['0','.5'],actual_17_chart_C3_transports=[dict(chart='first_micro',cells=cells)])]
    checks=0
    rates=dict(m=c.mpf(1),h=c.mpf('1.5'),k=c.mpf('1.5'),e=c.mpf(1),p=c.mpf(0))
    for text in ('.125','.25','.37','.5','.88','1'):
        s=iv.mpf(text);local,prefix,frames=fake._covers('.3',s)
        A,B,H=fake._correction_logs(local,prefix,False)
        for name,rate in rates.items():
            for k,bound in enumerate((c.mpf('1.39'),c.mpf('1.6'),c.mpf(2))):
                exact=c.quad(lambda t:c.exp(-rate*(c.mpf(text)-t)/7)*bound*c.cos(17*(t-c.mpf(1)/16)/7)/119,
                             [c.mpf(1)/16,c.mpf(text)])
                upper=c.exp(_bounds(c,H[name][k])[1])
                assert abs(exact)<=upper,(text,name,k,exact,upper)
                checks+=1
        total=sum((q[1] for q in prefix),iv.mpf(0))
        assert _bounds(c,total-(iv.mpf(text)-iv.mpf(1)/16))[1]>=0
        checks+=1
    # Direct differentiation of exp(A(Z)/N)-1 at Z=.3, independent of
    # the production normalized Taylor recurrence.
    z=c.mpf('.3');N=17;f=lambda x:c.exp((1+x+x*x)/N)-1
    x=[iv.mpf('1.39')/N,iv.mpf('1.6')/N,iv.mpf(1)/N]
    bounds=fake._exponential_increment_bounds(x)
    for k in range(3):
        exact=abs(c.diff(f,z,k)/mp.factorial(k))
        assert exact<=_bounds(c,bounds[k])[1];checks+=1
    return checks


def overlap(a,b):
    c=mp.mp.clone();c.dps=180
    assert abs(c.mpf(a['approximation'])-c.mpf(b['approximation']))<=c.mpf(a['absolute_error_upper'])+c.mpf(b['absolute_error_upper'])


def run():
    began=time.perf_counter();bench=independent_prefix_and_exponential_checks()
    print('Independent partial-prefix and exponential calculus passes',flush=True)
    field=OriginalFirstMicroNumeric();samples=[];native_checks=flat=seams=precision=0
    reference=json.loads((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_selected_gradient.json').read_text(encoding='utf-8'))
    assert reference['checks']['passed'] and reference['source_family']==field.family
    field.hashes['lei_ren_part1_paper_compliant_current_original_selected_gradient.json']=hashlib.sha256((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_selected_gradient.json').read_bytes()).hexdigest()
    for Z in ('-.3','0','.3'):
        for phase in ('.25','.625','1'):
            packet=field.evaluate(Z,phase,log_tau='-2',azimuth='.25')
            assert packet['selected_integer']['recipe']['J']==field.selection['J']
            for row in packet['small_error_proofs']:assert row['numerical_absolute_error_cap']=='1e-1000';native_checks+=1
            physical=packet['physical'];assert physical['actual_modified_leading_errors_included']
            assert not physical['physical_gradient_or_vorticity_available']
            assert physical['n1_own_history']['f_and_u_exact_zero']
            assert physical['shared_streamfunction_meridional_recovery']
            assert not packet['limits']['resolved_fast_phase_or_loop_point']
            for row in physical['cylindrical']['utheta']:
                if row['temporal_order']==1:
                    assert field.ctx.mpf(row['coefficient']['absolute_error_upper'])==0;flat+=1
            for name in ('delta_Psi','delta_Mpsi'):
                log=packet['normalized_third_derivative_log_error_bounds'][name]
                lo=mp.make_mpf(tuple(log['log_absolute_upper']['lower_exact_mpf_tuple']))
                assert lo>0,(phase,name)
            samples.append(packet)
        print('Actual selected-N first-micro enclosures: Z='+Z,flush=True)
    for multiple in ('.5','.75','1'):
        packet=field.selected('.3',multiple,log_tau='-2',azimuth='.25')
        if multiple in ('.5','.75'):
            assert packet['source_cover']['exact_flat_collar']
            for rows in packet['modified_error_contributions'].values():
                for row in rows:
                    lo,hi=_bounds(field.ctx,field._restore(row));assert lo==hi==0;flat+=1
        samples.append(packet)
    left=next(p for p in samples if p['coordinate']==dict(kind='selected_sc_multiple',multiple='.5',exact_sc_mpf_tuple=field.transport.sc_tuple))
    old=next(p for p in reference['samples'] if p['native']['Z']=='.3' and p['native']['fraction']=='2')
    mapping=dict(leading_radial_axis='ur:Qaxis',leading_radial_correction='ur:Qpsi',n1_radial='ur:q',
        leading_swirl='utheta:Phi',n1_swirl='utheta:f',leading_axis='uz:U0',leading_axial_correction='uz:Psi',n1_axial='uz:u')
    for key in ('ur','utheta','uz'):
        for row in left['physical']['cylindrical'][key]:
            overlap(row['coefficient'],old['term_pool'][mapping[row['name']]]['coefficient']);seams+=1
    actual=next(p for p in samples if p['Z']=='.3' and p['coordinate']==dict(kind='exact_decimal_phase',value='1'))
    primitive=left['physical']['n1_own_history']['actual_terminal_primitive']
    assert actual['physical']['n1_own_history']['actual_terminal_primitive']==primitive
    # Independent inspection of the physical leading value term weights is
    # supplied by the accepted core adapter; compare its original source seam.
    coarse=OriginalFirstMicroNumeric(precision=60,amplitude_precision=90)
    old=coarse.evaluate('.3','1',log_tau='-2',azimuth='.25')
    for name,rows in actual['normalized_modified_leading'].items():
        for a,b in zip(rows,old['normalized_modified_leading'][name]):overlap(a,b);precision+=1
    for component,rows in actual['physical']['cylindrical'].items():
        for a,b in zip(rows,old['physical']['cylindrical'][component]):overlap(a['coefficient'],b['coefficient']);precision+=1
    rejected=0
    for selected,Z,s,kwargs in ((False,'1','.25',{}),(False,'.3','0',{}),(False,'.3','1.1',{}),
            (True,'.3','.4',{}),(True,'.3','1.1',{}),(False,'.3','1',{'azimuth':'nan'}),
            (False,'.3','1',{'log_tau':'1'})):
        try:(field.selected if selected else field.evaluate)(Z,s,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Unsupported first-micro request admitted')
    implementation=dict(field.code_hashes)
    implementation[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,checks=dict(passed=True,
        independent_partial_prefix_and_exponential_comparisons=bench,actual_modified_first_micro_packets=len(samples),
        normalized_selected_N_small_error_proofs=native_checks,exact_flat_error_or_n1_swirl_checks=flat,
        original_rminus_physical_velocity_coefficient_seams=seams,terminal_n1_primitive_retained=True,
        cross_precision_native_and_physical_coefficient_overlaps=precision,selected_original_N_not_changed=True,
        third_derivative_caps_retained_as_unresolved_positive_logs=True,rejected_unsupported_requests=rejected),
        samples=samples,input_hashes=field.hashes,implementation_hashes=implementation,
        limitations=dict(full_fast_phase_point_resolution=False,leading_Phi_Psi_mean_Z_order=2,Q_order=1,
            numerical_U_Z3_and_QZZ=False,physical_gradient_and_vorticity_beyond_rminus=False,
            exact_signed_correction_integrals_resolved=False,global_N_admission=False,global_field=False,
            global_five_moment_repair=False,genuine_n_ge_2=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-began,3))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for name,digest in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    r=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=r['checks'],execution_seconds=r['execution_seconds']),indent=2),flush=True)
