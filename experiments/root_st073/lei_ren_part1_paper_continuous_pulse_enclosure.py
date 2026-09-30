"""Outward range-integration bounds for the complete continuous pulse rows.

The finite centered band is enclosed by interval rectangles. All omitted
support is positive and retained in the upper bound. Decimal input parameters
are treated as declared real numbers, not as certified source parameters.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_basis_enclosure import ContinuousBasisEnclosure


class ContinuousPulseEnclosure(ContinuousBasisEnclosure):
    def atoms(self,mu,row,*,panels=4096,band=48):
        if row not in (1,2) or not isinstance(panels,int) or panels<16:
            raise ValueError('Rows 1 or 2 and integer panels >=16 required')
        v=self.iv;mu=v.mpf(mu);band=v.mpf(band)
        if self.bounds(mu)[0]<=0 or self.bounds(band)[0]<=0:
            raise ValueError('Positive mu and band required')
        lam=v.mpf('.5')-row*mu;k=lam/mu
        if self.bounds(lam)[0]<=0:raise ValueError('Positive lambda required')
        u0=v.exp(v.log(2/k)/3);L=1/u0**2;w=u0/v.sqrt(6)
        if self.bounds(band*w)[1]>=self.bounds(v.mpf('.25'))[0]:
            raise ValueError('Centered band is too wide')
        if self.bounds(u0*(1+band*w))[1]>=self.bounds(v.mpf('.5'))[0]:
            raise ValueError('Band must remain inside u < .5')
        total=v.mpf(0);step=2*band/panels
        for j in range(panels):
            left=-band+2*band*j/panels
            right=-band+2*band*(j+1)/panels
            x=self._hull(left,right)
            scale=1+x*w;u=u0*scale
            # Interval square preserves nonnegativity on panels crossing zero.
            difference=x**2*(2*scale+1)/(6*scale**2)
            smooth=1/(1-u)**2
            f=(v.mpf('10.99')-u)*v.exp(smooth-difference)/(1+v.exp(-1/u**2+smooth))
            total+=step*f
        logpref=-2*k-3*L+2*v.log(u0)-v.log(6)/2-v.log(mu)
        core=v.exp(logpref)*total
        phases=[band**2*(2*(1+s*band*w)+1)/(6*(1+s*band*w)**2) for s in (-1,1)]
        # Sum of endpoint exponentials bounds their maximum, avoiding a
        # rounded ordinary-mpf minimum in the tail certificate.
        tails=[v.mpf('5.5')/mu*v.exp(4-2*k-3*L)*sum(v.exp(-p) for p in phases),
               v.mpf('5.5')/mu*v.exp(-v.mpf('2.5')*k),
               11/lam*v.exp(-3*k)]
        omitted=sum(tails)
        full=self._hull(core,core+omitted)
        return dict(centered=total,core=core,omitted=omitted,tails=tails,full=full)

    def report(self,mu,row,*,panels=4096,band=48,nominal=None):
        atoms=self.atoms(mu,row,panels=panels,band=band)
        result=dict(mu=str(mu),row=row,panels=panels,band=band,
            interval_precision=self.precision,
            centered_integral=self.describe(atoms['centered']),
            full_integral=self.describe(atoms['full']),
            log_full_integral=self.describe(self.iv.log(atoms['full'])),
            positive_omitted_upper_bound=self.describe(atoms['omitted']),
            omitted_piece_bounds=[self.describe(t) for t in atoms['tails']],
            pulse_quadrature_enclosed=True,arithmetic_rounding_enclosed=True,
            parameter_scope='Declared real decimal mu; inherited source uncertainty is not enclosed',
            full_axial_closure_enclosed=False,finite_energy_certified=False,
            incoming_uncertainty_enclosed=False,scale_recursion_certified=False)
        with mp.workdps(self.precision+30):
            lo,hi=self.bounds(atoms['centered'])
            result['centered_relative_width']=mp.nstr((hi-lo)/lo,24)
            if nominal is not None:
                saved=nominal.full_row(mu,row,band=band)
                center=mp.mpf(saved['centered_integral'])
                result['nominal_center_contained']=bool(lo<=center<=hi)
        return result


def run():
    from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
    source=json.loads(Path(__file__).with_name('lei_ren_part1_paper_continuous_incoming_outer.json').read_text())
    mu=source['continuous_solve']['input_mu']
    enclosure=ContinuousPulseEnclosure(precision=80)
    nominal=ContinuousAxialPulse(precision=100)
    reports=[]
    for panels in (1024,4096):
        for row in (1,2):
            print(f'pulse enclosure row={row} panels={panels}',flush=True)
            report=enclosure.report(mu,row,panels=panels,nominal=nominal)
            if not report['nominal_center_contained']:
                raise AssertionError('Nominal centered integral outside range enclosure')
            reports.append(report)
    result=dict(method='Outward interval rectangles plus positive omitted-support bounds',reports=reports)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([dict(row=r['row'],panels=r['panels'],relative_width=r['centered_relative_width']) for r in reports]),flush=True)
    return result


if __name__=='__main__':run()
