"""Focused finite-long logarithm, kernel, R110 seam and history checks."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_long_numeric import OriginalLongNumeric,gradient_jet
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_long_numeric.json'


def independent_checks():
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    from lei_ren_part1_paper_compliant_long_reshape_profiles import backward_kernel,normalized_amplitude_decay
    c=mp.mp.clone();c.dps=65;iv=MPIntervalContext();iv.dps=100
    f=object.__new__(OriginalLongNumeric);f.ctx=c;f.iv=iv;f.J=J
    from types import SimpleNamespace
    f.leading=SimpleNamespace(_abs=lambda x:abs(x),_symmetric=lambda x:iv.mpf([-_bounds(c,x)[1],_bounds(c,x)[1]]))
    def cap(log,proof,label):
        if log is None:return iv.mpf(0)
        assert _bounds(c,log)[1]<_bounds(c,iv.ln(iv.mpf('1e-1000')))[0]
        proof.append(label);return iv.mpf('1e-1000')
    f.owner=SimpleNamespace(_cap=cap)
    def contains(value,truth):
        lo,hi=_bounds(c,truth)
        assert lo-c.mpf('1e-55')*max(1,abs(value))<=value<=hi+c.mpf('1e-55')*max(1,abs(value)),(value,lo,hi)
    count=0
    def sigma(q):
        if q<=0:return c.mpf(0)
        if q>=1:return c.mpf(1)
        odds=1/(1-q)**2-1/q**2
        return 1/(1+c.exp(-odds)) if odds>=0 else c.exp(odds)/(1+c.exp(odds))
    z=c.mpf('.3');T=c.mpf(7);Bi=J.variable(iv,iv.mpf('.3'),2)
    Bi=(1+Bi+Bi*Bi)/100
    B=(1+z+z*z)/100;B1=(1+2*z)/100;B2=c.mpf(1)/100
    initial=lambda Z:2+Z+Z*Z
    initialj=J(iv,[iv.mpf(initial(z)),iv.mpf(1+2*z),iv.mpf(1)])
    for qvalue in ('.2','.65'):
        q=c.mpf(qvalue);y=T*q;s=sigma(q)
        for rate,m in (('1.6',1),('1.2',2),('.2',2)):
            k=c.mpf(rate);radius=8*m*abs(B)/T
            kernel=backward_kernel(iv,Bi,iv.mpf(7),iv.mpf(str(y)),m,
                iv.mpf(str(k-radius)),iv.mpf(str(k+radius)),[])
            for order in range(3):
                def integrand(t):
                    diff=s-sigma(q-t/T);value=c.exp(-k*t+m*B*diff)
                    return value*(1 if order==0 else m*B1*diff if order==1 else
                                  m*B2*diff+(m*B1*diff)**2/2)
                actual=c.quad(integrand,[0,y/2,y]);contains(actual,kernel[order]);count+=1
            logj=Bi*iv.mpf(str(m*s))-iv.mpf(str(k*y))
            memory=normalized_amplitude_decay(iv,logj,iv.mpf(str(k-radius)),
                iv.mpf(str(k+radius)),iv.mpf(str(y)),[],initial=initialj)
            for order in range(3):
                exact=c.diff(lambda Z:initial(Z)*c.exp(-k*y+m*(1+Z+Z*Z)*s/100),z,order)/math.factorial(order)
                contains(exact,memory[order]);count+=1
    # Test the actual log-profile operator independently of a copied formula.
    from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm
    pj=J(c,[B,B1,B2]);tj=J(iv,[iv.mpf(B),iv.mpf(B1),iv.mpf(B2)])
    for qraw in ('0','.2','.65','1'):
        q=c.mpf(qraw);sbp,sbt,cbp,cbt=f._cutoff_times(pj,tj,q,iv.mpf(qraw),[])
        for order in range(3):
            contains(pj[order]*sigma(q),sbt[order]);contains(pj[order]*(1-sigma(q)),cbt[order]);count+=2
    # A flat source is multiplied before its enormous coefficient is capped.
    huge=J(iv,[iv.exp(10000)]*3);proof=[]
    for qraw in ('1e-20','0.99999999999999999999'):
        _,a,_,b=f._cutoff_times(J(c,[c.exp(10000)]*3),huge,c.mpf(qraw),iv.mpf(qraw),proof)
        small=a if c.mpf(qraw)<c.mpf('.5') else b
        assert all(_bounds(c,abs(v))[1]<=c.mpf('1e-999') for v in small.coefficients);count+=1
    assert len(proof)==6;count+=1
    for zraw in ('-.3','0','.17'):
        z=c.mpf(zraw);j=c.mpf('0.001');delta=c.mpf('.0001');width=j/500
        result=gradient_jet(J,iv,iv.mpf(zraw),iv.mpf(str(j)),iv.mpf(str(delta)),iv.mpf(str(width)),2)
        def rational(Z):
            H=-4*Z**3-j*Z**2+(9-delta)*Z/2+j
            return (1-delta*Z*Z)*H/(H*H+width*width)
        for k in range(3):contains(c.diff(rational,z,k)/math.factorial(k),result[k]);count+=1
        phi=J(iv,[iv.mpf(1+z*z),iv.mpf(2*z),iv.mpf(1)])
        logphi=logarithm(phi)
        for k in range(3):contains(c.diff(lambda Z:c.ln(1+Z*Z),z,k)/math.factorial(k),logphi[k]);count+=1
    return count


def run():
    began=time.perf_counter();independent=independent_checks()
    print('Independent finite kernels, incoming memories, cutoff and analytic G jets pass',flush=True)
    f=OriginalLongNumeric();samples=[];proofs=0
    for Z in ('-.3','0','.3'):
        for q in ('0','.5','1'):samples.append(f.evaluate(Z,q))
        print('Actual original finite long: Z='+Z,flush=True)
    samples.extend(f.at_log_radius('.3',y) for y in ('1','10'))
    reference=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_r110_numeric.json')).read_text(encoding='utf-8'))
    seams=history=terminal=finite=measure=0;c,iv=f.ctx,f.iv;J=f.J
    pick=lambda Z,q:next(p for p in samples if p['Z']==Z and p['coordinate']['kind']=='exact_decimal_q' and p['coordinate']['value']==q)
    for Z in ('-.3','0','.3'):
        old=next(p for p in reference['samples'] if p['Z']==Z and p['source_cover']['local_source_chart']=='post_power' and p['coordinate']['value']=='1')
        start,end=pick(Z,'0'),pick(Z,'1');inp=f.inputs(Z)
        for n in ('Psi','Mpsi','Qpsi'):
            for a,b in zip(old['normalized_modified_leading'][n],start['normalized_modified_leading'][n]):overlap(a,b);seams+=1
        # LogPhi_N is joined against log of the actual modified inlet jet.
        rows=old['normalized_modified_leading']['Phi']
        jt=J(iv,[f._restore(row) for row in rows]);jp=J(c,[c.mpf(row['approximation']) for row in rows])
        for a,b in zip(f.leading._jet_summary(f.logarithm(jp),f.logarithm(jt)),start['normalized_modified_leading']['logPhi']):overlap(a,b);seams+=1
        for name,raw,power in (('H_over_Phi','H',1),('Kpsi_over_Phi','Kpsi',1),('B_over_Phi2','B',2),('C_over_Phi2','C',2),('Mpsi','Mpsi',0),('Apsi','Apsi',0)):
            p,t=inp;ep=p['six'][raw].truncate(2);et=t['six'][raw].truncate(2)
            if power:
                ep=ep/(p['values']['Phi'].truncate(2)**power);et=et/(t['values']['Phi'].truncate(2)**power)
            for a,b in zip(f.leading._jet_summary(ep,et),start['unmodified_six_current_amplitude_units'][name]):overlap(a,b);seams+=1
        assert start['original_five_correction_history_bounds']==old['original_five_correction_history_bounds'];history+=20
        assert start['positive_order_n1_own_history']['actual_terminal_primitive']==old['physical']['n1_own_history']['actual_terminal_primitive'];history+=1
        # All G and inlet derivatives must cancel from the terminal swirl.
        z=J.variable(iv,iv.mpf(Z),3);zp=J.variable(c,c.mpf(Z),3)
        expected=-f.logarithm(1+z*z)+f.T/10-iv.mpf(f.owner.amplitude.logC_iv)
        nominal=-f.logarithm(1+zp*zp)+f.Tp/10-c.mpf(f.owner.amplitude.logC)
        for a,b in zip(f.leading._jet_summary(nominal,expected),end['physical_similarity_log_Utheta']):overlap(a,b);terminal+=1
        for p in (start,end):
            assert p['original_T_singleton_exact_mpf_tuple']==f.T_tuple
            assert p['selected_integer']['recipe']==f.owner.integer;history+=2
    for p in samples:
        assert [len(p['normalized_modified_leading'][n]) for n in ('logPhi','Psi','Mpsi','Qpsi')]==[4,4,3,2]
        assert len(p['unmodified_six_current_amplitude_units'])==6
        assert not p['limits']['global_field'] and not p['limits']['genuine_temporal_n_ge_2']
        assert p['six_history_contract']['all_nonzero_actual_R110_inlets_retained'];history+=1
        for row in p['small_error_proofs']:
            if isinstance(row,dict) and row.get('numerical_absolute_error_cap')=='1e-1000':proofs+=1
        for rows in p['normalized_modified_leading'].values():
            for row in rows:
                assert c.isfinite(c.mpf(row['approximation'])) and c.isfinite(c.mpf(row['absolute_error_upper']));finite+=1
        q=f._restore(p['coordinate']['q']);y=f._restore(p['coordinate']['y'])
        yl,yh=_bounds(c,y-f.T*q);assert yl<=0<=yh;measure+=1
    # Native finite-y prefixes must keep T*dq and cannot reset to equilibrium.
    for p in samples[-2:]:
        assert p['coordinate']['kind']=='exact_decimal_R110_log_offset'
        assert p['unmodified_six_current_amplitude_units']['C_over_Phi2'][0]['approximation']!=pick('.3','1')['unmodified_six_current_amplitude_units']['C_over_Phi2'][0]['approximation'];finite+=1
    coarse=OriginalLongNumeric(precision=60,amplitude_precision=90)
    old=coarse.evaluate('.3','.5');new=pick('.3','.5');comparisons=0
    for group in ('normalized_modified_leading','unmodified_six_current_amplitude_units','finite_backward_kernels'):
        for n,rows in new[group].items():
            for a,b in zip(rows,old[group][n]):overlap(a,b);comparisons+=1
    for a,b in zip(new['physical_similarity_log_Utheta'],old['physical_similarity_log_Utheta']):overlap(a,b);comparisons+=1
    rejected=0
    for Z,q in (('1','.5'),('.3','-.01'),('.3','1.01'),('.3','nan'),('.3','.')):
        try:f.evaluate(Z,q)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Unsupported original long request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,independent_finite_kernel_memory_cutoff_G_checks=independent,
        actual_original_long_packets=len(samples),selected_N_and_flat_tail_small_error_proofs=proofs,
        R110_field_and_six_history_seam_coefficients=seams,incoming_primitive_integer_and_history_checks=history,
        exact_terminal_G_and_inlet_cancellations=terminal,finite_log_coefficients_and_nonconstant_histories=finite,
        exact_T_native_measure_relations=measure,cross_precision_coefficients=comparisons,rejected_requests=rejected),
        samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        limitations=samples[0]['limits'],execution_seconds=round(time.perf_counter()-began,3))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for name,digest in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],execution_seconds=result['execution_seconds']),indent=2),flush=True)
