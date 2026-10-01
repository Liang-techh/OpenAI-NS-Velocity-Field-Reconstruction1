"""Independent known nonlinear inverse and implicit derivative fixture."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_five_bump_inverse import certify
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    c=MPIntervalContext();c.dps=90
    with mp.workdps(140):
        centers=[mp.mpf('1.25'),mp.mpf('1.5'),mp.mpf('1.75')]
        L=mp.matrix(5,5)
        L[0,0]=L[0,1]=1
        L[1,0]=centers[0]**mp.mpf('.6');L[1,1]=centers[2]**mp.mpf('.6')
        for i,x in enumerate(centers):
            L[2,i+2]=mp.sqrt(x);L[3,i+2]=-x**mp.mpf('.1');L[4,i+2]=x**mp.mpf('-.9')
        fg=[mp.mpf(30),mp.mpf(40)];gg=[mp.mpf(31),mp.mpf(32)]
        ff=[mp.mpf(33),mp.mpf(34),mp.mpf(35)];fx=[mp.mpf(36),mp.mpf(37),mp.mpf(38)]
        W=dict(L=[[c.mpf(L[i,j]) for j in range(5)] for i in range(5)],
            fg=list(map(c.mpf,fg)),gg=list(map(c.mpf,gg)),ff=list(map(c.mpf,ff)),ff_over_x=list(map(c.mpf,fx)))
        def known(z):
            return [-(1+z)*mp.mpf('1e-15'),(2-z)*mp.mpf('1e-15'),
                (1+z)*mp.mpf('1e-34'),-(2+z)*mp.mpf('1e-34'),(3-z)*mp.mpf('1e-34')]
        def amplitude(z):return mp.mpf('1e-11')+z*mp.mpf('1e-12')
        # Independent MP formula builds data from a prescribed exact solution.
        def data(z):
            h=known(z)
            q=mp.matrix([0,fg[0]*h[0]*h[2]+fg[1]*h[1]*h[4],0,
                amplitude(z)*(gg[0]*h[0]**2+gg[1]*h[1]**2)-sum(ff[i]*h[i+2]**2 for i in range(3))/2,
                sum(fx[i]*h[i+2]**2 for i in range(3))/2])
            return -(L*mp.matrix(h)+q)
        z=mp.mpf('.5')
        val=data(z);tangent=[mp.diff(lambda zz:data(zz)[i],z) for i in range(5)]
        d=[IntervalTaylor(c,[c.mpf(val[i]),c.mpf(tangent[i])]) for i in range(5)]
        inv=IntervalTaylor(c,[c.mpf(amplitude(z)),c.mpf('1e-12')])
        cert=certify(c,W,d,inv)
        if not cert['certified']:raise AssertionError('known map contraction failed')
        count=0
        for i,jet in enumerate(cert['controls']):
            for k,ref in enumerate((known(z)[i],mp.diff(lambda zz:known(zz)[i],z))):
                lo,hi=endpoints(jet[k])
                if not lo<=ref<=hi:raise AssertionError(('known implicit coefficient excluded',i,k))
                count+=1
        rejected=certify(c,W,[IntervalTaylor(c,[c.mpf('.1'),c.mpf(0)])]+d[1:],inv)
        if rejected['certified'] or rejected['self_map_strictly_inside']:
            raise AssertionError('outside-box problem accepted')
    here=Path(__file__).parent
    report=dict(passed=True,fixture_only=True,known_solution_value_derivative_coefficients_contained=count,
        outside_box_data_rejected=True,input_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest()
            for name in ('lei_ren_part1_paper_interval_five_bump_inverse.py',Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Known nonlinear inverse:10 value/derivative coefficients contained; outside-box data rejected',flush=True)
    return report

if __name__=='__main__':run()
