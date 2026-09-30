"""Positive range bounds for startup, plateau and cutoff pulse energy.

Startup primitive bounds use monotone sigma Riemann sums; no nested numerical
quadrature is treated as exact. The plateau has an analytic antiderivative.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_basis_enclosure import ContinuousBasisEnclosure


class ContinuousPulseEnergyEnclosure(ContinuousBasisEnclosure):
    def sigma_node(self,index,panels):
        v=self.iv
        if index==0:return v.mpf(0)
        if index==panels:return v.mpf(1)
        # Symmetry avoids exponentially large logistic numerator values.
        if 2*index>panels:return 1-self.sigma_node(panels-index,panels)
        q=v.mpf(index)/panels
        tiny=v.exp(-1/q**2+1/(1-q)**2)
        return tiny/(1+tiny)

    def atoms(self,*,panels=4096):
        if not isinstance(panels,int) or panels<16:
            raise ValueError('Integer panels >=16 required')
        v=self.iv;step=v.mpf(1)/panels
        sigma=[self.sigma_node(j,panels) for j in range(panels+1)]
        startup=v.mpf(0);primitive=v.mpf(0);cutoff=v.mpf(0)
        for j in range(panels):
            left=v.mpf(j)/panels;right=v.mpf(j+1)/panels
            next_primitive=primitive+step*self._hull(sigma[j],sigma[j+1])
            # gp(xi)=int_0^(50xi) sigma(q)dq / 50 on startup.
            gp=self._hull(primitive/50,next_primitive/50)
            xi=self._hull(left/50,right/50)
            startup+=step/50*v.exp(-2*xi)*gp**2
            primitive=next_primitive
            xi=self._hull(10+left,10+right)
            # sigma(11-xi) is decreasing along this panel.
            switch=self._hull(sigma[panels-j-1],sigma[panels-j])
            cutoff+=step*v.exp(-2*xi)*(xi-v.mpf('.01'))**2*switch**2
        def antiderivative(x):
            y=x-v.mpf('.01')
            return -v.exp(-2*x)*(y**2/2+y/2+v.mpf('.25'))
        plateau=antiderivative(v.mpf(10))-antiderivative(v.mpf('.02'))
        return dict(startup=startup,plateau=plateau,cutoff=cutoff,
                    total=startup+plateau+cutoff)

    def report(self,*,panels=4096,nominal=None):
        atoms=self.atoms(panels=panels)
        result=dict(panels=panels,interval_precision=self.precision,
            parts={k:self.describe(x) for k,x in atoms.items()},
            pulse_energy_quadrature_enclosed=True,
            parameter_scope='Exact decimal startup .02, shift .01 and cutoff [10,11] of ContinuousAxialPulse',
            finite_energy_certified=False,full_axial_closure_enclosed=False)
        with mp.workdps(self.precision+20):
            lo,hi=self.bounds(atoms['total'])
            result['relative_width']=mp.nstr((hi-lo)/lo,30)
            if nominal is not None:result['nominal_contained']=bool(lo<=nominal<=hi)
        return result


def run():
    source=json.loads(Path(__file__).with_name('lei_ren_part1_paper_continuous_incoming_outer.json').read_text())
    with mp.workdps(120):nominal=mp.mpf(source['continuous_solve']['K_p'])
    enclosure=ContinuousPulseEnergyEnclosure(precision=80)
    reports=[]
    for panels in (1024,4096):
        report=enclosure.report(panels=panels,nominal=nominal)
        if not report['nominal_contained']:raise AssertionError('Nominal pulse energy outside bounds')
        reports.append(report)
    result=dict(method='Monotone startup primitive bounds, analytic plateau, cutoff interval rectangles',reports=reports)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([dict(panels=r['panels'],relative_width=r['relative_width'],nominal_contained=r['nominal_contained']) for r in reports]),flush=True)
    return result


if __name__=='__main__':run()
