"""Pressure-component propagation through the Section 9.23 comparison ODE.

This adapter preserves a finite pressure-parameter Taylor jet through the
same core log-slope, axial, moment and pressure equations. It is the auxiliary
comparison, not the prescribed exit bridge, reshape, or global installer.
"""
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_component_pressure_core import PressurePolynomial
from lei_ren_part1_paper_pressure_parameter_jet import PressureParameterJet
from lei_ren_part1_paper_exit_comparison import Section923Comparison


class _ComponentCoreView:
    def __init__(self,core,order):
        self.core=core;self.order=order;self.precision=core.precision
        self.Lambda=core.axis.Lambda;self.delta=core.axis.delta
    def convert(self,value):
        if isinstance(value,PressurePolynomial):return PressureParameterJet(value,order=self.order)
        if isinstance(value,dict):return {k:self.convert(v) for k,v in value.items()}
        if isinstance(value,list):return [self.convert(v) for v in value]
        return value
    @lru_cache(maxsize=32)
    def coefficients(self,Z):
        original=self.core.coefficients(Z)
        result=dict(original)
        for name in ('F','Uz','P'):result[name]=self.convert(original[name])
        return result
    def evaluate(self,R,Z):return self.convert(self.core.evaluate(R,Z))


class ComponentExitComparison(Section923Comparison):
    def __init__(self,bundle,*,pressure_order=9,h_b='.005',transition_steps=16):
        if int(pressure_order)!=pressure_order or pressure_order<1:
            raise ValueError('Positive integer pressure parameter order required')
        core=bundle.get('component_pressure_core')
        if core is None:raise ValueError('Bundle requires a complete component pressure core')
        self.pressure_order=int(pressure_order)
        revised=dict(bundle)
        revised.update(core=_ComponentCoreView(core,self.pressure_order),precision=core.precision,
            scalar_converter=lambda x:PressureParameterJet(x,order=self.pressure_order))
        super().__init__(revised,h_b=h_b,transition_steps=transition_steps)
    def metadata(self):
        result=super().metadata()
        result.update(pressure_source='complete H=1 preheat pressure components',
            pressure_parameter_order=self.pressure_order,pressure_parameter_order_truncated=True,
            pressure_parameter_remainder_enclosed=False,transition_ODE_error_enclosed=False,
            tiny_collar_width_increment_error_enclosed=False,
            prescribed_exit_bridge_installed=False,reshape_installed=False,global_field_installed=False)
        return result


def fixture():
    """Resolved-tail scalar comparisons at zero and nonzero pressure input."""
    from types import SimpleNamespace
    from lei_ren_part1_paper_core_recursion import core_coefficients
    from lei_ren_part1_paper_core_adapter import CorePolynomial
    from lei_ren_part1_paper_component_pressure_core import ComponentPressureCore
    with mp.workdps(100):
        f=[mp.mpf(2),mp.mpf('.3')]+[mp.mpf(0)]*3
        u=[mp.mpf('.4'),mp.mpf('.2')]+[mp.mpf(0)]*3
        p=[mp.mpf(-1),mp.mpf('.2')]+[mp.mpf(0)]*3
        tail=[mp.mpf('.0001'),mp.mpf('.0003'),mp.mpf('-.0002')]+[mp.mpf(0)]*2
        kw=dict(F0_Z_taylor=f,U0_Z_taylor=u,radial_degree=3,precision=100)
        formal=core_coefficients('.3','.01',P0_Z_taylor=[PressurePolynomial({0:a,1:b}) for a,b in zip(p,tail)],scalar_converter=PressurePolynomial,**kw)
        axis=SimpleNamespace(precision=100,Lambda=mp.mpf(10),delta=mp.mpf('.01'))
        component=ComponentPressureCore(axis,SimpleNamespace(precision=100),3)
        component.coefficients=lambda Z:formal
        base=dict(axis=axis,Lambda=axis.Lambda,precision=100,radial_degree=3,component_pressure_core=component)
        comparison=ComponentExitComparison(base,pressure_order=9,transition_steps=8)
        maximum=mp.mpf(0)
        for parameter in (0,1):
            scalar=core_coefficients('.3','.01',P0_Z_taylor=[a+parameter*b for a,b in zip(p,tail)],**kw)
            scalar_core=CorePolynomial(lambda Z:scalar,Lambda=10,delta='.01',precision=100)
            scalar_comparison=Section923Comparison(dict(base,core=scalar_core),transition_steps=8)
            for y in ('0','.0075','.01','.02'):
                a=comparison.evaluate(y,'.3');b=scalar_comparison.evaluate(y,'.3')
                for name in ('F','Uz','P','PZ','D','E','T_theta','T_z'):
                    error=abs(a[name].evaluate(parameter)-b[name])/max(1,abs(b[name]))
                    maximum=max(maximum,error)
                for name in a['moments']:
                    error=abs(a['moments'][name].evaluate(parameter)-b['moments'][name])/max(1,abs(b['moments'][name]))
                    maximum=max(maximum,error)
            ca=comparison.frozen_driver_coefficients('.3');cb=scalar_comparison.frozen_driver_coefficients('.3')
            for name in ('D','E'):
                for key,value in ca[name].items():
                    maximum=max(maximum,abs(value.evaluate(parameter)-cb[name][key])/max(1,abs(cb[name][key])))
        assert maximum<mp.mpf('1e-28'),maximum
        return dict(maximum_scaled_difference=mp.nstr(maximum,40),pressure_parameter_order=9,
            pressure_parameter_remainder_enclosed=False,scalar_comparison_parameters=[0,1])


if __name__=='__main__':print(fixture())
