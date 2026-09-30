"""Callable same-source background with incremental five-bump inverse."""
import operator
import mpmath as mp
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_five_bump_inverse import iterate_inverse
from lei_ren_part1_paper_five_bump_field import evaluate_correction,join


class FiveBumpInverseField:
    """Resolve the five-bump map from each common source Z evaluation.

The source remains responsible for its axial domain, pressure datum, bounds,
and remainders. No sample or cached point is extrapolated to another Z.
"""
    def __init__(self,reference,moment_map,*,steps=10):
        if not hasattr(reference,'defect_data') or not hasattr(reference,'evaluate_x'):
            raise TypeError('reference requires defect_data(Z) and evaluate_x(x,Z)')
        self.reference=reference;self.moment_map=moment_map
        self.steps=operator.index(steps)
        if self.steps<0:raise ValueError('steps must be nonnegative')
        self._cache={}

    def inverse_data(self,Z):
        with mp.workdps(self.moment_map.precision):
            z=mp.mpf(Z);key=mp.nstr(z,self.moment_map.precision)
            if key not in self._cache:
                source=self.reference.defect_data(z)
                amplitude=AxialDual(source['Am'],source['Am_Z'])
                result=iterate_inverse(self.moment_map,[source['d_dual'][i] for i in range(1,6)],
                                       amplitude,steps=self.steps)
                self._cache[key]=(amplitude,result)
            return self._cache[key]

    def evaluate_x(self,x,Z):
        with mp.workdps(self.moment_map.precision):
            baseline=self.reference.evaluate_x(x,Z)
            amplitude,inverse=self.inverse_data(Z)
            correction=evaluate_correction(self.moment_map,inverse['h'],amplitude,
                       self.reference.Rm,Z,x,delta=self.reference.delta)
            field=join(baseline,correction,delta=self.reference.delta)
            field['inverse_increments']=inverse['increments']
            field['inverse_terminal_residual']=inverse['terminal_residual']
            field['inverse_metadata']=inverse['metadata']
            field['finite_inverse_updates']=self.steps
            field['infinite_response_error_enclosed']=False
            return field

    def evaluate_phase(self,t,Z):
        with mp.workdps(self.moment_map.precision):
            return self.evaluate_x(mp.exp(mp.mpf(t)-2),Z)
