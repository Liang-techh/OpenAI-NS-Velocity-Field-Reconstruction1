"""One heat-kernel functional with separately retained tiny deficit and jets.

H(x)=E[(1+x V)^(-h)], V~Gamma(1+h,1). Small-x jets use
one quadratic Taylor polynomial, with analytic truncation bounds. Bounds do
not enclose floating arithmetic; non-small-x quadrature remains uncertified.
This provider is not yet installed in the shared field or its moment targets.
"""
import json
from pathlib import Path
import mpmath as mp


class ContinuousHeatKernel:
    def __init__(self,h,*,precision=443,quadrature_precision=80):
        self.precision=int(precision)
        self.quadrature_precision=int(quadrature_precision)
        with mp.workdps(self.precision):
            self.h=mp.mpf(str(h))
            if self.h<=0:raise ValueError('Require h>0')

    def jet(self,xi):
        with mp.workdps(self.precision):
            x=mp.mpf(str(xi));h=self.h
            if x<0 or not mp.isfinite(x):raise ValueError('Require finite xi>=0')
            c1=h*(1+h);c2=h*(h+1)**2*(h+2)
            c3=h*(h+1)**2*(h+2)**2*(h+3)
            if x<=mp.mpf('1e-8'):
                deficit=c1*x-c2*x*x/2
                prime=-c1+c2*x;second=c2
                bounds=[c3*x**3/6,c3*x*x/2,c3*x]
                method='shared_quadratic_Taylor_with_analytic_truncation_bounds'
                prime_correction=c2*x
            else:
                with mp.workdps(self.quadrature_precision):
                    # expm1 retains H-1 rather than subtracting two full values.
                    gamma=mp.gamma(1+h)
                    def weight(v):return mp.exp(-v)*v**h/gamma
                    # Normalize by h before quadrature so absolute stopping
                    # tolerances cannot discard the whole tiny integral.
                    deficit=h*mp.quad(lambda v:weight(v)*(-mp.expm1(-h*mp.log1p(x*v))/h),[0,1,4,16,mp.inf])
                    prime=-h*mp.quad(lambda v:weight(v)*v/(1+x*v)**(h+1),[0,1,4,16,mp.inf])
                    second=h*(h+1)*mp.quad(lambda v:weight(v)*v*v/(1+x*v)**(h+2),[0,1,4,16,mp.inf])
                prime_correction=prime+c1
                bounds=None;method='shared_positive_Gamma_integrals'
            if deficit<0 or deficit>=1:raise ArithmeticError('Heat deficit outside [0,1)')
            return dict(H=1-deficit,H_minus_one=-deficit,deficit=deficit,
                logH=mp.log1p(-deficit),H_prime=prime,H_second=second,
                H_prime_at_zero=-c1,H_prime_correction=prime_correction,
                truncation_absolute_bounds=bounds,method=method,
                arithmetic_error_enclosed=False,quadrature_error_enclosed=False,
                installed_in_shared_field=False,finite_energy_certified=False)


def run():
    with mp.workdps(100):
        engine=ContinuousHeatKernel('5e-201')
        rows=[]
        for x in ('0','1e-12','.2','1'):
            row=engine.jet(x);xx=mp.mpf(x)
            if not(row['H_prime']<0 and row['H_second']>0):
                raise ArithmeticError('Heat monotonicity/convexity lost')
            if xx>mp.mpf('1e-8'):
                step=mp.mpf('1e-5')
                neighbors={i:engine.jet(xx+i*step) for i in (-2,-1,1,2)}
                dminus=sum(w*neighbors[i]['H_minus_one'] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
                dprime=sum(w*neighbors[i]['H_prime'] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
                errors=[abs(dminus/row['H_prime']-1),abs(dprime/row['H_second']-1)]
                if max(errors)>mp.mpf('1e-12'):raise ArithmeticError('Heat functional jet mismatch')
            else:errors=[]
            rows.append(dict(xi=x,deficit=mp.nstr(row['deficit'],40),
                H_prime=mp.nstr(row['H_prime'],40),H_second=mp.nstr(row['H_second'],40),
                derivative_relative_errors=[mp.nstr(e,40) for e in errors],method=row['method']))
        tiny=engine.jet(mp.exp(-mp.mpf('1e6')))
        if not(tiny['deficit']>0 and tiny['logH']<0 and tiny['H_prime_correction']>0):
            raise ArithmeticError('Tiny heat terms were lost')
        report=dict(samples=rows,tiny_log_deficit=mp.nstr(mp.log(tiny['deficit']),40),
            tiny_logH_retained=True,tiny_prime_correction_retained=True,
            common_functional=True,analytic_Taylor_truncation_bounds=True,
            arithmetic_error_enclosed=False,quadrature_error_enclosed=False,
            installed_in_shared_field=False,complete_heat_moments=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':run()
