"""Independent tanh-sinh check of saved ordinary/peak rate-two beta kernels.

Reads accepted coefficient enclosures without replaying the reconstruction.
"""
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds

DIRECTORY=ROOT/'experiments/root_st073'


def run():
    start=time.perf_counter();c=mp.mp.clone();c.dps=140;iv=MPIntervalContext();iv.dps=160
    path=DIRECTORY/(PREFIX+'compliant_current_original_axial_m15_numeric.json')
    raw=path.read_bytes();data=json.loads(raw);assert data['checks']['passed']
    checks=[]
    exact=lambda v:c.make_mpf(tuple(v['exact_mpf_tuple']))
    for value in ('.5','1'):
        for power,Z,mode,termname in ((1,'0','beta_rate2','G_m_minus_Ib'),(2,'-.3','beta2_rate2','G_b2')):
            sample=next(p for p in data['samples'] if p['Z']==Z and p['native_coordinate']==value)
            proof=sample['positive_kernel_proofs'][mode]
            anchor,origin,scale,left,right=(exact(proof[n]) for n in
                ('normalized_log_anchor','q_origin','q_scale','integration_left','integration_right'))
            U=c.exp(40*c.mpf(value));r0=1-c.mpf(value)
            def integrand(z):
                q=origin+scale*z;r=r0-c.log1p(-q/U)/40
                odds=1/(1-r)**2-1/r**2
                logbeta=-c.log1p(c.exp(-odds)) if odds>0 else odds-c.log1p(c.exp(odds))
                return c.exp(power*logbeta-2*q-anchor)
            split=sorted({left,right}|{c.mpf(v) for v in ('-16','-8','-4','-1','0','1','4','8','16') if left<c.mpf(v)<right})
            direct=c.quad(integrand,split)
            term=next(t for t in sample['m15_C0_C1_terms'][176:] if t['name']=='O2_axial_'+termname)
            tuples=term['coefficient_C0_C1'][0]['source_enclosure_exact_mpf_tuples']
            interval=iv.mpf([c.make_mpf(tuple(v)) for v in tuples])
            if power==1:factor=iv.mpf(-2)
            else:
                z=iv.mpf(Z);delta=iv.exp(-4*(iv.exp(40)+11)-30);A=1-delta*z*z
                P=z**4*(delta*delta+delta)-z*z*(2*delta+3)+3
                factor=16*z*z*P/A**3
            lo,hi=_bounds(c,interval/factor)
            assert lo<=direct<=hi,(value,power,c.nstr(direct,20),c.nstr(lo,20),c.nstr(hi,20))
            checks.append(dict(s=value,power=power,normalized_direct_value=c.nstr(direct,45),
                saved_interval_relative_width=c.nstr((hi-lo)/lo,12),independent_quadrature_inside_saved_enclosure=True))
            print('Independent saved rate-two kernel: s='+value+', beta power='+str(power),flush=True)
    source='scripts/verify_original_axial_kernels.py'
    return dict(checks=dict(passed=True,independent_ordinary_and_peak_integrals=len(checks)),samples=checks,
        report_file=path.name,report_sha256=hashlib.sha256(raw).hexdigest(),
        implementation_hashes={source:hashlib.sha256((ROOT/source).read_bytes()).hexdigest()},
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(saved_kernel_validation_only=True,ancestors_replayed=False,new_global_physical_claim=False))


if __name__=='__main__':
    data=run();path=DIRECTORY/(PREFIX+'compliant_current_original_axial_kernels_check.json')
    path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
