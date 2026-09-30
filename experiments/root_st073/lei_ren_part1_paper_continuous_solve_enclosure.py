"""Conditional interval coefficient solve over continuous integral bounds.

Incoming rows, target and pulse energy are supplied intervals, or explicitly
declared stored dyadic values. Only basis/pulse integration is certified here.
The materialized runtime coefficients are never overwritten by this report.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_basis_enclosure import ContinuousBasisEnclosure
from lei_ren_part1_paper_continuous_pulse_enclosure import ContinuousPulseEnclosure


class ContinuousSolveEnclosure(ContinuousBasisEnclosure):
    def declared_value(self,value):
        """A supplied iv interval, decimal real, or exact stored MP dyadic.

        A stored MP number is not a certificate of its originating integral.
        Its exact dyadic value can nevertheless be a conditional parameter.
        """
        if hasattr(value,'_mpi_'):return self.iv.make_mpf(value._mpi_)
        if isinstance(value,mp.mpf):
            return self.iv.make_mpf((value._mpf_,value._mpf_))
        return self.iv.mpf(value)

    def solve_atoms(self,mu,base,target,Kp,*,panels=4096):
        v=self.iv;mu=v.mpf(mu)
        basis=self.basis_atoms(mu,panels=panels)
        pulse=ContinuousPulseEnclosure(precision=self.precision)
        pulse_atoms=[pulse.atoms(mu,i,panels=panels) for i in (1,2)]
        p=[self.declared_value(atom['full']) for atom in pulse_atoms]
        matrix=[row['matrix'] for row in basis['rows']]
        det=basis['determinant']
        if self.bounds(det)[1]>=0:raise ArithmeticError('Determinant not strictly negative')
        base=[self.declared_value(x) for x in base]
        target=self.declared_value(target);Kp=self.declared_value(Kp)
        if self.bounds(Kp)[0]<=0:raise ValueError('Pulse energy must be positive')
        def affine(rhs):
            return [(-rhs[0]*matrix[1][1]+matrix[0][1]*rhs[1])/det,
                    (-matrix[0][0]*rhs[1]+rhs[0]*matrix[1][0])/det]
        c0=affine(base);cp=affine(p)
        K=[v.exp(-26+i*mu)*basis['energy_gram'] for i in (6,2)]
        A=Kp+mu*sum(k*x**2 for k,x in zip(K,cp))
        B=2*mu*sum(k*x*y for k,x,y in zip(K,c0,cp))
        C=mu*sum(k*x**2 for k,x in zip(K,c0))-target
        if self.bounds(A)[0]<=0 or self.bounds(C)[1]>=0:
            raise ArithmeticError('Bounds do not establish a unique positive energy branch')
        D=B**2-4*A*C
        if self.bounds(D)[0]<=0:raise ArithmeticError('Discriminant not strictly positive')
        # Both are identities for the unique positive root. Select an
        # interval form whose denominator is positive; intersect when both
        # can be evaluated, retaining the common guaranteed root enclosure.
        direct=(-B+v.sqrt(D))/(2*A)
        denominator=B+v.sqrt(D)
        amplitude=direct
        if self.bounds(denominator)[0]>0:
            rational=-2*C/denominator
            # Endpoint comparison must be performed at enough MP precision.
            with mp.workdps(self.precision+20):
                lower=max((direct._mpi_[0],rational._mpi_[0]),key=mp.mpf)
                upper=min((direct._mpi_[1],rational._mpi_[1]),key=mp.mpf)
                if mp.mpf(lower)>mp.mpf(upper):raise ArithmeticError('Root enclosures disjoint')
            amplitude=v.make_mpf((lower,upper))
        if self.bounds(amplitude)[0]<=0:raise ArithmeticError('Positive root enclosure unresolved')
        coefficients=[x+amplitude*y for x,y in zip(c0,cp)]
        return dict(basis=basis,pulse=p,matrix=matrix,base=base,target=target,Kp=Kp,
            K=K,c0=c0,cp=cp,quadratic=[A,B,C],discriminant=D,
            amplitude=amplitude,coefficients=coefficients)

    def report(self,runtime,*,panels=4096):
        atoms=self.solve_atoms(runtime.solve_receipt['input_mu'],runtime.base,
            runtime.target,runtime.Kp,panels=panels)
        def contains(interval,point):
            lo,hi=self.bounds(interval)
            return bool(lo<=point<=hi)
        def relative_width(interval):
            with mp.workdps(self.precision+20):
                lo,hi=self.bounds(interval)
                scale=min(abs(lo),abs(hi)) if lo*hi>0 else max(abs(lo),abs(hi))
                return mp.nstr((hi-lo)/scale,30) if scale else None
        # Fixed stored coefficients are tested against true integral ranges.
        # Containment of zero means compatibility, never proof of equality.
        stored_a=self.declared_value(runtime.a)
        stored_c=[self.declared_value(c) for c in runtime.c]
        balances=[atoms['base'][i]+stored_a*atoms['pulse'][i]+
            sum(atoms['matrix'][i][j]*stored_c[j] for j in range(2)) for i in range(2)]
        return dict(panels=panels,interval_precision=self.precision,
            conditional_unique_positive_branch=True,
            amplitude=self.describe(atoms['amplitude']),
            coefficients=[self.describe(c) for c in atoms['coefficients']],
            relative_widths=dict(amplitude=relative_width(atoms['amplitude']),
                coefficients=[relative_width(c) for c in atoms['coefficients']]),
            nominal_containment=dict(amplitude=contains(atoms['amplitude'],runtime.a),
                coefficients=[contains(x,y) for x,y in zip(atoms['coefficients'],runtime.c)]),
            quadratic=[self.describe(x) for x in atoms['quadratic']],
            fixed_materialized_row_balances=[self.describe(x) for x in balances],
            fixed_materialized_row_zero_compatible=[contains(x,0) for x in balances],
            input_parameter_scope='Exact stored dyadic incoming rows, energy target and Kp; their source uncertainty is NOT enclosed. Mu is the declared decimal receipt value.',
            basis_and_pulse_quadrature_enclosed=True,
            inherited_input_uncertainty_enclosed=False,pulse_energy_uncertainty_enclosed=False,
            materialized_coefficients_replaced=False,global_mean_closed=False,
            finite_energy_certified=False,scale_recursion_certified=False)


def run():
    from lei_ren_part1_paper_continuous_axial_runtime import _fixture,SharedContinuousAxialRuntime
    source,receipt=_fixture();runtime=SharedContinuousAxialRuntime(source,receipt)
    encloser=ContinuousSolveEnclosure(precision=80)
    reports=[]
    for panels in (1024,4096):
        print(f'conditional coefficient enclosure panels={panels}',flush=True)
        report=encloser.report(runtime,panels=panels)
        if not report['nominal_containment']['amplitude'] or not all(report['nominal_containment']['coefficients']):
            raise AssertionError('Nominal solve outside conditional interval solution')
        reports.append(report)
    result=dict(reports=reports,full_axial_closure_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([dict(panels=r['panels'],widths=r['relative_widths'],containment=r['nominal_containment']) for r in reports]),flush=True)
    return result


if __name__=='__main__':run()
