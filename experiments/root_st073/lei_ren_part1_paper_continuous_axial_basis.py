"""Shared continuous bump values, jets and weighted primitives.

Definition: beta_ell(s)=exp(-1/(1-(s/ell)^2))/(ell*I0) on |s|<ell,
I0=integral[-1,1] exp(-1/(1-x^2)) dx. All full/partial integrals use
this same definition. MP quadrature evaluations have no certified enclosure;
this layer alone does not close the pulse, incoming mass or physical tail.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp


class ContinuousAxialBump:
    def __init__(self,*,precision=160,ell='.15'):
        self.precision=precision
        with mp.workdps(precision):
            self.ell=mp.mpf(ell)
            if self.ell<=0:raise ValueError('Positive source support required')
            self.normalizer=self._integrate_raw(mp.mpf(0),mp.mpf(1))

    @staticmethod
    def raw(x):
        return mp.exp(-1/(1-x*x)) if abs(x)<1 else mp.mpf(0)

    def _integrate_raw(self,lam,upper):
        if upper<=-1:return mp.mpf(0)
        upper=min(upper,mp.mpf(1))
        points=[mp.mpf(-1)]
        for p in ('-.5','0','.5'):
            if mp.mpf(p)<upper:points.append(mp.mpf(p))
        points.append(upper)
        if upper<mp.mpf('-.5') and abs(lam*self.ell)<=1:
            # A flat endpoint integral may be far below mp.quad's absolute
            # stopping scale. Factor out the monotone endpoint maximum so
            # quadrature sees an O(1) integrand and refines relative shape.
            log_endpoint=lam*self.ell*upper-1/(1-upper*upper)
            def scaled(x):
                if abs(x)>=1:return mp.mpf(0)
                return mp.exp(lam*self.ell*x-1/(1-x*x)-log_endpoint)
            return mp.exp(log_endpoint)*mp.quad(scaled,points)
        return mp.quad(lambda x:mp.exp(lam*self.ell*x)*self.raw(x),points)

    def values(self,s):
        with mp.workdps(self.precision):
            s=mp.mpf(s); x=s/self.ell
            if abs(x)>=1:return dict(beta=mp.mpf(0),beta_s=mp.mpf(0),beta_ss=mp.mpf(0))
            beta=self.raw(x)/(self.ell*self.normalizer)
            dlog=-2*x/(self.ell*(1-x*x)**2)
            ddlog=-2*(1+3*x*x)/(self.ell**2*(1-x*x)**3)
            return dict(beta=beta,beta_s=beta*dlog,beta_ss=beta*(dlog*dlog+ddlog))

    @lru_cache(maxsize=32)
    def full(self,lam_key):
        with mp.workdps(self.precision):
            lam=mp.mpf(lam_key)
            # Identical continuous integral atom at lambda=0, not a
            # separately measured total being overwritten with one.
            raw=self.normalizer if lam==0 else self._integrate_raw(lam,mp.mpf(1))
            return raw/self.normalizer

    def primitive(self,lam,s):
        with mp.workdps(self.precision):
            lam=mp.mpf(lam);s=mp.mpf(s)
            if s<=-self.ell:return mp.mpf(0)
            if s>=self.ell:return self.full(mp.nstr(lam,self.precision))
            return self._integrate_raw(lam,s/self.ell)/self.normalizer

    def primitive_jet(self,lam,s):
        with mp.workdps(self.precision):
            lam=mp.mpf(lam);s=mp.mpf(s);v=self.values(s)
            return dict(value=self.primitive(lam,s),
                derivative=mp.exp(lam*s)*v['beta'],
                second_derivative=mp.exp(lam*s)*(lam*v['beta']+v['beta_s']))

    def tail(self,lam,s):
        """Remaining weighted mass, integrated directly by reflection.

        beta is even, so int_s^ell exp(lam*t) beta(t) dt equals
        primitive(-lam,-s). This preserves a small positive endpoint tail
        without subtracting two nearly equal complete/partial integrals.
        """
        with mp.workdps(self.precision):
            return self.primitive(-mp.mpf(lam),-mp.mpf(s))

    def matrix(self,mu):
        with mp.workdps(self.precision):
            mu=mp.mpf(mu)
            if not 0<mu<=mp.mpf(1)/60:raise ValueError('Source mu outside range')
            lambdas=[mp.mpf('.5')-i*mu for i in (1,2)]
            B=[self.full(mp.nstr(l,self.precision)) for l in lambdas]
            matrix=[[mp.exp(center*l)*b for center in (-3,-1)]
                for l,b in zip(lambdas,B)]
            # Positive B and mu establish invertibility of this continuous
            # basis by the analytic identity, without subtracting close rows.
            determinant=-2*mp.exp(-2+6*mu)*mp.sinh(mu)*B[0]*B[1]
            return matrix,determinant

    def energy_gram(self,mu):
        with mp.workdps(self.precision):
            mu=mp.mpf(mu)
            value=mp.quad(lambda x:mp.exp(-2*mu*self.ell*x)*self.raw(x)**2,
                [-1,-.5,0,.5,1])/(self.ell*self.normalizer**2)
            return value


def run():
    print('integrating shared continuous axial bump basis',flush=True)
    source=json.loads(Path(__file__).with_name('lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    basis=ContinuousAxialBump(precision=160)
    higher=ContinuousAxialBump(precision=200)
    with mp.workdps(200):
        mu=mp.mpf(source['axial']['input_mu'])
        lam=mp.mpf('.5')-mu;s=mp.mpf('-.0375')
        n=lambda x:mp.nstr(x,60)
        matrix,det=basis.matrix(mu)
        high_matrix,high_det=higher.matrix(mu)
        refinement=max(abs(x/y-1) for row,other in zip(matrix,high_matrix) for x,y in zip(row,other))
        old=[[mp.mpf(x) for x in row] for row in source['axial']['linear_matrix']]
        shift=max(abs(x/y-1) for row,other in zip(matrix,old) for x,y in zip(row,other))
        jet=basis.primitive_jet(lam,s); rows=[]
        for h in (mp.mpf('1e-5'),mp.mpf('5e-6')):
            values={i:basis.primitive(lam,s+i*h) for i in (-2,-1,1,2)}
            derivative=(values[-2]-8*values[-1]+8*values[1]-values[2])/(12*h)
            rows.append(dict(step=n(h),relative_derivative_error=n(abs(derivative/jet['derivative']-1))))
        report=dict(precision=160,refined_precision=200,mu=n(mu),
            continuous_normalizer=n(basis.normalizer),
            full_mass_from_identical_atom=n(basis.full('0')),
            continuous_matrix=[[n(x) for x in row] for row in matrix],
            analytic_determinant=n(det),
            determinant_relative_refinement=n(abs(det/high_det-1)),
            matrix_relative_refinement=n(refinement),
            matrix_relative_change_from_float_raw_canonical_basis=n(shift),
            energy_gram=n(basis.energy_gram(mu)),primitive_derivative_replay=rows,
            endpoint_values={side:{k:n(v) for k,v in basis.values(side*basis.ell).items()} for side in (-1,1)},
            pointwise_and_primitive_definition_shared=True,
            numerical_quadrature_enclosure_certified=False,
            installed_in_global_profile=False,continuous_mean_closure_certified=False,
            finite_energy_certified=False,
            scope='Continuous bump sublayer only; incoming/core primitive, pulse full/partial integrals and Z jets still require shared provenance.')
        if refinement>mp.mpf('1e-130'):raise ArithmeticError('Continuous basis precision refinement failed')
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({k:report[k] for k in ('matrix_relative_refinement','matrix_relative_change_from_float_raw_canonical_basis','primitive_derivative_replay')}),flush=True)
        return report


if __name__=='__main__':run()
