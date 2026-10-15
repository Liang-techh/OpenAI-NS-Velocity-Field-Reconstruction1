"""Evaluate original-N macro values, true geometry and nonzero own history."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_macro_numeric import OriginalMacroNumeric
from openai_ns_reconstruction.original_core_numeric import _bounds
from run_original_second_micro_numeric import overlap

NAME='lei_ren_part1_paper_compliant_current_original_macro_numeric.json'


def independent_checks():
    c=mp.mp.clone();c.dps=90;iv=MPIntervalContext();iv.dps=130
    module=importlib.import_module('lei_ren_part1_paper_compliant_actual_bridge_integrals')
    checks=0
    # Ordinary quadrature of every modal kernel, including resonances.
    for text in ('.3','1.7'):
        L=c.mpf(text);r=c.mpf(3);r1=r*c.exp(L)
        for power in (1,2):
            a=module.radius_kernels(iv,iv.mpf(3),iv.mpf(c.nstr(r1,90)),iv.mpf(text),power)
            b=module.volterra_primitive_kernels(iv,iv.mpf(3),iv.mpf(c.nstr(r1,90)),iv.mpf(text),power,1)
            for j in range(3):
                value=c.quad(lambda y:r**power*c.exp((power-j)*y),[0,L])
                center=sum(_bounds(c,a[j]))/2
                assert abs(value-center)<c.mpf('1e-80');checks+=1
                # Swap the Volterra integrals; no closed-form kernel is used.
                value=c.quad(lambda y:r**power*c.exp((power-j)*y)*(1-c.exp(-(L-y))),[0,L])
                center=sum(_bounds(c,b[j]))/2
                assert abs(value-center)<c.mpf('1e-80');checks+=1
    # Nonzero incoming mean, a changing velocity, and its Z derivatives.
    for L in (c.mpf('.3'),c.mpf('1.7')):
        z=c.mpf('.31');initial=lambda z:2+z+z*z;core=lambda z:3-z
        amplitude=lambda z:1+z+z*z
        exact=lambda z:c.exp(-L)*initial(z)+c.quad(lambda y:c.exp(-(L-y))*(core(z)+amplitude(z)*(1-c.exp(-y))),[0,L])
        baseline=lambda z:c.exp(-L)*initial(z)+(1-c.exp(-L))*core(z)
        for k in range(3):
            error=abs(c.diff(exact,z,k)-c.diff(baseline,z,k))
            bound=(1-c.exp(-L))*abs(c.diff(amplitude,z,k))
            assert 0<error<=bound;checks+=1
    fake=object.__new__(OriginalMacroNumeric)
    fake.transport=SimpleNamespace(ctx=c,iv=iv,phase_iv=iv.mpf(1)/32)
    h=iv.mpf(1)/101
    fake.leading=SimpleNamespace(core=SimpleNamespace(lp_iv=iv.mpf(0),ll_iv=iv.mpf(0)),
        logh=-iv.ln(101),hcap=h,_width_error=lambda x:h*abs(x))
    fake.Y=iv.ln(25);fake.logN=iv.ln(17)
    def cap(v):
        t=iv.ln(v)
        return dict(exact_zero=False,log_absolute_upper=dict(lower_exact_mpf_tuple=list(t._mpi_[0]),upper_exact_mpf_tuple=list(t._mpi_[1])))
    rows=[cap(3)]*4
    def cell(left,right):
        return dict(exact_left=left,exact_right=right,actual_source_and_primitive_range=dict(A=rows,B=rows),
            actual_density_ranges=dict(N_scaled={name:rows for name in ('m','h','k','e','p')}))
    fake.ranges=[dict(exact_Z_cell=['0','.5'],actual_17_chart_C3_transports=[
        dict(chart='first_micro',cells=[cell(['sc_half'],[1,1])]),
        dict(chart='second_micro',cells=[cell([1,1],[2,1])]),
        dict(chart='frozen_macro',cells=[cell([0,1],[1,2]),cell([1,2],[1,1])])])]
    for q in ('0','.2','1'):
        qi=iv.mpf(q);qp=c.mpf(q)
        rho,theta,mass,base,tiny=fake._geometry(qi)
        actualrho=4*c.exp(c.ln(25)*qp+2*(1-qp)/101)
        decay=c.exp(-(c.ln(25)-c.mpf(2)/101)*qp)
        for box,value in ((rho,actualrho),(theta,decay),(mass,1-decay)):
            lo,hi=_bounds(c,box);assert lo-c.mpf('1e-80')<=value<=hi+c.mpf('1e-80');checks+=1
        local,prefix,_=fake._covers('.3',qi)
        _,_,hist=fake._correction_logs(local,prefix,False)
        macro=[p for p in prefix if p[3]=='frozen_macro']
        width=sum((p[1] for p in macro),iv.mpf(0))
        assert abs(_bounds(c,width)[1]-qp)<c.mpf('1e-80');checks+=1
        assert all(p[2]._mpi_==iv.ln(fake.Y)._mpi_ for p in macro);checks+=1
        # Constant signed source: macro measure is L*dq, not hb*dq.
        for name,rate in (('m',1),('h',c.mpf('1.5')),('k',c.mpf('1.5')),('e',1),('p',0)):
            length=(c.ln(25)-c.mpf(2)/101)*qp
            start=c.mpf(1)/16;micro=(2-start)/101
            expected=c.mpf(3)/17*(c.quad(lambda y:c.exp(-rate*(micro+length-y)),[0,micro])
                +c.quad(lambda y:c.exp(-rate*(length-y)),[0,length]))
            assert 0<expected<=c.exp(_bounds(c,hist[name][0])[1]);checks+=1
    return checks


def run():
    began=time.perf_counter();independent=independent_checks()
    print('Independent true-radius, Volterra, mean and Jacobian checks pass',flush=True)
    f=OriginalMacroNumeric();samples=[];proofs=0
    for Z in ('-.3','0','.3'):
        for q in ('0','1e-18','.5','1'):
            p=f.evaluate(Z,q,log_tau='-2',azimuth='.25');samples.append(p)
            assert [len(p['normalized_modified_leading'][n]) for n in ('Phi','Psi','Mpsi','Qpsi')]==[4,4,3,2]
            assert p['source_cover']['cumulative_history_charts']==['first_micro','second_micro','frozen_macro']
            assert p['physical']['n1_own_history']['current_rho_not_rho4_used']
            assert not p['physical']['physical_gradient_or_vorticity_available']
            assert not p['limits']['global_field'] and not p['limits']['actual_mean_Z3_numeric']
            for row in p['small_error_proofs']:assert row['numerical_absolute_error_cap']=='1e-1000';proofs+=1
        print('Original-N true macro geometry/own history: Z='+Z,flush=True)
    samples.append(f.evaluate('.17','.37',log_tau='-3',azimuth='1.1'))
    reference=json.loads((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_second_micro_numeric.json').read_text(encoding='utf-8'))
    seams=histories=geometry=0
    for Z in ('-.3','0','.3'):
        previous=next(p for p in reference['samples'] if p['Z']==Z and p['coordinate']['value']=='2')
        inlet=next(p for p in samples if p['Z']==Z and p['coordinate']['value']=='0')
        for name,rows in previous['normalized_modified_leading'].items():
            for a,b in zip(rows,inlet['normalized_modified_leading'][name]):overlap(a,b);seams+=1
        for name,rows in previous['physical']['cylindrical'].items():
            for a,b in zip(rows,inlet['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);seams+=1
        # Same nonzero incoming original history, allowing a changed sum
        # majorant but requiring overlapping directed logarithmic bounds.
        for name,rows in previous['original_five_correction_history_bounds'].items():
            for a,b in zip(rows,inlet['original_five_correction_history_bounds'][name]):
                assert a==b;histories+=1
        narrow=next(p for p in samples if p['Z']==Z and p['coordinate']['value']=='1e-18')
        assert narrow['normalized_modified_leading']['Mpsi'][0]['approximation']!=inlet['normalized_modified_leading']['Mpsi'][0]['approximation'];histories+=1
        for p in [p for p in samples if p['Z']==Z]:
            assert p['physical']['n1_own_history']['actual_terminal_primitive']==inlet['physical']['n1_own_history']['actual_terminal_primitive'];histories+=1
            n=f._numeric(Z,p['coordinate']['value'],False)
            rho=n[3];row=p['physical']['n1_own_history']['mean_u']
            truth=f.transport.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])
            primitive=p['physical']['n1_own_history']['actual_terminal_primitive']['actual_primitive'][0]
            pbox=f.transport.iv.mpf([mp.make_mpf(tuple(v)) for v in primitive['source_enclosure_exact_mpf_tuples']])
            lo,hi=_bounds(f.transport.ctx,truth*rho-pbox);assert lo<=0<=hi;geometry+=1
        exit=next(p for p in samples if p['Z']==Z and p['coordinate']['value']=='1')
        assert exit['physical']['native']['R100_endpoint_exact'];geometry+=1
        for a,b in zip(exit['normalized_modified_leading']['Mpsi'],exit['normalized_modified_leading']['Psi']):overlap(a,b);histories+=1
    coarse=OriginalMacroNumeric(precision=60,amplitude_precision=90)
    old=coarse.evaluate('.3','1e-18',log_tau='-2',azimuth='.25')
    new=next(p for p in samples if p['Z']=='.3' and p['coordinate']['value']=='1e-18')
    comparisons=0
    for name,rows in new['normalized_modified_leading'].items():
        for a,b in zip(rows,old['normalized_modified_leading'][name]):overlap(a,b);comparisons+=1
    for name,rows in new['physical']['cylindrical'].items():
        for a,b in zip(rows,old['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);comparisons+=1
    for a,b in zip(new['physical']['coordinates']['r']['log_scale_terms'],old['physical']['coordinates']['r']['log_scale_terms']):
        if 'approximation' in a and 'approximation' in b:overlap(a,b);comparisons+=1
    rejected=0
    for Z,q,kwargs in (('1','.5',{}),('.3','-.01',{}),('.3','1.01',{}),('.3','nan',{}),('.3','1',{'log_tau':'1'}),('.3','1',{'azimuth':'nan'})):
        try:f.evaluate(Z,q,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Unsupported macro request admitted')
    try:f.selected('.3','.5')
    except ValueError:rejected+=1
    else:raise AssertionError('Selected-sc coordinate admitted as macro')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,
        independent_modal_Volterra_mean_geometry_and_measure_checks=independent,
        actual_modified_macro_packets=len(samples),selected_N_small_error_proofs=proofs,
        second_micro_macro_native_physical_coefficient_seams=seams,
        nonzero_incoming_histories_mean_relaxation_and_terminal_memory_checks=histories,
        actual_rho_primitive_and_R100_geometry_checks=geometry,
        cross_precision_coefficients_and_log_geometry_comparisons=comparisons,
        rejected_unsupported_requests=rejected),samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        limitations=dict(Phi_Psi_Z_order=3,own_mean_Z_order=2,Qpsi_Z_order=1,
            actual_comparison_frozen_but_fields_evolve=True,total_radial_derivatives=False,
            resolved_fast_phase=False,resolved_signed_correction_integrals=False,
            mean_Z3_and_QZZ_numeric=False,later_charts=False,global_N_admission=False,
            global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False),
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
