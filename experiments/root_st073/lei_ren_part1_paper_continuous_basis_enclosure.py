"""Outward-rounded midpoint enclosures for the continuous end-bump atoms.

The quadrature error is bounded analytically, separately from interval
arithmetic rounding. Exact dyadic interval endpoints are serialized alongside
approximate decimal displays. Pulse/input/energy-target errors are not covered.
"""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


class ContinuousBasisEnclosure:
    def __init__(self,*,precision=60,ell='.15'):
        if precision<40:raise ValueError('At least 40 interval digits required')
        self.precision=int(precision)
        self.iv=MPIntervalContext()
        self.iv.dps=self.precision
        self.ell=self.iv.mpf(ell)
        if self.bounds(self.ell)[0]<=0:raise ValueError('Positive ell required')
        v=self.iv
        self.D0=v.exp(-v.mpf(1))
        self.D1=8*v.exp(-v.mpf(2))
        self.D2=1024*v.exp(-v.mpf(4))+216*v.exp(-v.mpf(3))+8*v.exp(-v.mpf(2))

    def bounds(self,value):
        with mp.workdps(self.precision+20):
            return tuple(mp.mpf(endpoint) for endpoint in value._mpi_)

    def _hull(self,lower,upper):
        return self.iv.make_mpf((lower._mpi_[0],upper._mpi_[1]))

    def raw_integral(self,k,*,panels=4096,squared=False):
        """Enclose int[-1,1] exp(k*x)*phi(x)**power dx, power=1 or 2.

        Midpoints are made from exact integers/rational interval operations.
        With a=1/(1-x^2), |phi'|<=D1 and |phi''|<=D2. The composite midpoint
        remainder is at most 8*sup|f''|/(24*n^2). All constants are intervals.
        """
        if not isinstance(panels,int) or panels<16:
            raise ValueError('Integer panel count >=16 required')
        v=self.iv;k=v.mpf(k);n=v.mpf(panels)
        total=v.mpf(0)
        for index in range(panels):
            # This interval contains the exact rational midpoint, even when
            # its denominator is not a power of two.
            x=-v.mpf(1)+v.mpf(2*index+1)/n
            phi=v.exp(-v.mpf(1)/(1-x*x))
            total+=v.exp(k*x)*(phi*phi if squared else phi)
        midpoint=2*total/n
        if squared:
            d0=self.D0*self.D0
            d1=2*self.D0*self.D1
            d2=2*self.D1*self.D1+2*self.D0*self.D2
        else:
            d0,d1,d2=self.D0,self.D1,self.D2
        derivative_bound=v.exp(abs(k))*(d2+2*abs(k)*d1+k*k*d0)
        quadrature_bound=derivative_bound/(3*n*n)
        enclosure=midpoint+v.mpf([-1,1])*quadrature_bound
        return dict(enclosure=enclosure,midpoint=midpoint,
                    quadrature_bound=quadrature_bound,
                    derivative_bound=derivative_bound)

    def basis_atoms(self,mu,*,panels=4096):
        v=self.iv;mu=v.mpf(mu)
        lo,hi=self.bounds(mu)
        if lo<=0 or hi>self.bounds(v.mpf(1)/60)[0]:
            raise ValueError('Require positive source mu <=1/60')
        normalizer=self.raw_integral(0,panels=panels)
        I0=normalizer['enclosure']
        if self.bounds(I0)[0]<=0:
            raise ArithmeticError('Normalizer enclosure is not strictly positive')
        rows=[];B=[]
        for index in (1,2):
            lam=v.mpf('.5')-index*mu
            raw=self.raw_integral(lam*self.ell,panels=panels)
            atom=raw['enclosure']/I0
            B.append(atom)
            rows.append(dict(lambda_value=lam,raw=raw,B=atom,
                matrix=[v.exp(center*lam)*atom for center in (-3,-1)]))
        # sinh(mu) lies between mu and mu*exp(mu) for positive mu. This
        # avoids subtracting nearly equal exponentials or matrix products.
        sinh=self._hull(mu,mu*v.exp(mu))
        determinant=-2*v.exp(-2+6*mu)*sinh*B[0]*B[1]
        energy_raw=self.raw_integral(-2*mu*self.ell,panels=panels,squared=True)
        gram=energy_raw['enclosure']/(self.ell*I0*I0)
        return dict(normalizer=normalizer,rows=rows,determinant=determinant,
                    energy_raw=energy_raw,energy_gram=gram)

    def describe(self,value):
        """Export exact binary endpoints; decimal strings are only displays."""
        def endpoint(t):
            sign,mantissa,exponent,bits=t
            with mp.workdps(self.precision+20):
                display=mp.nstr(mp.mpf(t),40)
            return dict(sign=int(sign),mantissa=str(mantissa),
                        binary_exponent=int(exponent),decimal_display=display)
        lo,hi=self.bounds(value)
        with mp.workdps(self.precision+20):
            width_display=mp.nstr(hi-lo,40)
        return dict(lower=endpoint(value._mpi_[0]),upper=endpoint(value._mpi_[1]),
                    width_display=width_display,
                    exact_endpoint_definition='(-1)^sign * mantissa * 2^binary_exponent')

    def report(self,mu,*,panels=4096,nominal=None):
        atoms=self.basis_atoms(mu,panels=panels)
        def integral_report(atom):
            return dict(enclosure=self.describe(atom['enclosure']),
                midpoint=self.describe(atom['midpoint']),
                analytic_quadrature_bound=self.describe(atom['quadrature_bound']))
        report=dict(mu=str(mu),panels=panels,interval_precision=self.precision,
            normalizer=integral_report(atoms['normalizer']),
            weighted_rows=[dict(B=self.describe(row['B']),
                raw=integral_report(row['raw']),
                matrix=[self.describe(entry) for entry in row['matrix']])
                for row in atoms['rows']],
            determinant=self.describe(atoms['determinant']),
            determinant_strictly_negative=self.bounds(atoms['determinant'])[1]<0,
            energy_gram=self.describe(atoms['energy_gram']),
            energy_raw=integral_report(atoms['energy_raw']),
            basis_quadrature_enclosed=True,arithmetic_rounding_enclosed=True,
            parameter_scope='Declared real decimal mu and ell, not an enclosure of inherited source parameter uncertainty',
            full_axial_closure_enclosed=False,pulse_quadrature_enclosed=False,
            incoming_uncertainty_enclosed=False,finite_energy_certified=False,
            scale_recursion_certified=False)
        if nominal is not None:
            with mp.workdps(max(self.precision+20,nominal.precision)):
                mu_mp=mp.mpf(mu)
                matrix,det=nominal.matrix(mu_mp)
                def contains(interval,value):
                    lo,hi=self.bounds(interval)
                    return bool(lo<=value<=hi)
                report['nominal_containment']=dict(
                    normalizer=contains(atoms['normalizer']['enclosure'],nominal.normalizer),
                    matrix=all(contains(row['matrix'][j],matrix[i][j])
                        for i,row in enumerate(atoms['rows']) for j in range(2)),
                    determinant=contains(atoms['determinant'],det),
                    energy_gram=contains(atoms['energy_gram'],nominal.energy_gram(mu_mp)))
        return report


def run():
    from lei_ren_part1_paper_continuous_axial_basis import ContinuousAxialBump
    folder=Path(__file__).parent
    source=json.loads((folder/'lei_ren_part1_paper_continuous_incoming_outer.json').read_text())
    mu=source['continuous_solve']['input_mu']
    encloser=ContinuousBasisEnclosure(precision=60)
    nominal=ContinuousAxialBump(precision=100)
    reports=[]
    for panels in (1024,4096):
        print('basis enclosure panels='+str(panels),flush=True)
        report=encloser.report(mu,panels=panels,nominal=nominal)
        if not report['determinant_strictly_negative'] or not all(report['nominal_containment'].values()):
            raise AssertionError('Nominal atoms failed independent interval containment')
        reports.append(report)
    result=dict(method='Outward interval midpoint rule plus analytic global second-derivative remainder',
                reports=reports,global_finite_energy_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(containment=reports[-1]['nominal_containment'],
                         determinant_strictly_negative=True)),flush=True)
    return result


if __name__=='__main__':run()
