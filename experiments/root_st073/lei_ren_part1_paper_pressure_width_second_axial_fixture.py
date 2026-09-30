"""Independent axial stencils through comparison, bridge, moments and Ur."""
import json
from math import comb
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_component_pressure_core import PressurePolynomial
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_pressure_width_comparison import PressureWidthComparison
from lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_second_axial_comparison import SecondAxialPressureWidthComparison
from lei_ren_part1_paper_pressure_width_second_axial_bridge import SecondAxialPressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_second_axial_continuation import SecondAxialPressureWidthExitContinuation
from lei_ren_part1_paper_pressure_width_axial_comparison import AxialPressureWidthComparison
from lei_ren_part1_paper_pressure_width_axial_bridge import AxialPressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_continuation import PressureWidthExitContinuation
from lei_ren_part1_paper_pressure_width_switches import PressureWidthExitSwitches
from lei_ren_part1_paper_pressure_width_second_axial_switches import SecondAxialPressureWidthExitSwitches


def run():
    with mp.workdps(120):
        z=mp.mpf('.3');h=mp.mpf('1e-8');width=mp.mpf('1e-4')
        axis=SimpleNamespace(precision=120,Lambda=mp.mpf(10),delta=mp.mpf('.01'))
        formal=core_coefficients(z,axis.delta,
            F0_Z_taylor=[mp.mpf(v) for v in ('2','.3','.02','.001','.0004','.0001','.00001','.000001')],
            U0_Z_taylor=[mp.mpf(v) for v in ('.4','.2','.01','.001','.0004','.0001','.00001','.000001')],
            P0_Z_taylor=[PressurePolynomial({0:mp.mpf(v)}) for v in ('-1','.2','.03','.002','.0003','.0001','.00001','.000001')],
            radial_degree=3,precision=120,scalar_converter=PressurePolynomial)
        def bundle(offset):
            shifted=dict(formal)
            for name in ('F','Uz','P'):
                shifted[name]=[[sum((comb(k,j)*offset**(k-j)*row[k] for k in range(j,len(row))),PressurePolynomial(0))
                                for j in range(len(row))] for row in formal[name]]
            component=SimpleNamespace(axis=axis,precision=120,coefficients=lambda _Z:shifted)
            return dict(axis=axis,precision=120,component_pressure_core=component)
        comparison=SecondAxialPressureWidthComparison(bundle(0),h_b=width,pressure_order=2,width_order=1)
        bridge=SecondAxialPressureWidthExitBridge(comparison,steps=8)
        out_comparison=comparison.evaluate('1.5',z)
        out_bridge=bridge.evaluate('.8',z)
        continuation=SecondAxialPressureWidthExitContinuation(bridge)
        out_continuation=continuation.evaluate_R('.8',z)
        switches=SecondAxialPressureWidthExitSwitches(continuation,steps=8)
        out_switches=switches.evaluate_R(110,z)
        samples={}
        for k in (-2,-1,0,1,2):
            scalar=PressureWidthComparison(bundle(k*h),h_b=width,pressure_order=2,width_order=1)
            scalar_bridge=PressureWidthExitBridge(scalar,steps=8)
            first=AxialPressureWidthComparison(bundle(k*h),h_b=width,pressure_order=2,width_order=1)
            first_bridge=AxialPressureWidthExitBridge(first,steps=8)
            first_continuation=PressureWidthExitContinuation(first_bridge)
            first_switches=PressureWidthExitSwitches(first_continuation,steps=8)
            samples[k]=(scalar.evaluate('1.5',z+k*h),scalar_bridge.evaluate('.8',z+k*h),
                        first_continuation.evaluate_R('.8',z+k*h),first_switches.evaluate_R(110,z+k*h))
        def nominal(v):return v.evaluate(pressure=0,width=1)
        errors={}
        def check(name,second,values):
            oracle=(-values[2]+16*values[1]-30*values[0]+16*values[-1]-values[-2])/(12*h*h)
            errors[name]=abs(nominal(second)-oracle)/max(1,abs(oracle))
        for name in ('F','Uz','P','I_theta','I_z','D','E'):
            check('comparison_'+name,out_comparison[name].second,
                  {k:nominal(v[0][name]) for k,v in samples.items()})
        for name in ('F','Uz','P'):
            check('bridge_'+name,out_bridge[name+'_ZZ'],
                  {k:nominal(v[1][name]) for k,v in samples.items()})
        for name,second in out_bridge['moments_ZZ'].items():
            check('moment_'+name,second,{k:nominal(v[1]['moments'][name]) for k,v in samples.items()})
        for name in ('F','Uz','P'):
            check('continuation_'+name,out_continuation[name+'_ZZ'],
                  {k:nominal(v[2][name]) for k,v in samples.items()})
        for name,second in out_continuation['moments_ZZ'].items():
            check('continuation_moment_'+name,second,
                  {k:nominal(v[2]['moments'][name]) for k,v in samples.items()})
        for name,second in out_continuation['raw_quadratic_integrals_ZZ'].items():
            check('continuation_raw_'+name,second,
                  {k:nominal(v[2]['raw_quadratic_integrals'][name]) for k,v in samples.items()})
        continuation_ur_z=(nominal(samples[-2][2]['Ur'])-8*nominal(samples[-1][2]['Ur'])
                          +8*nominal(samples[1][2]['Ur'])-nominal(samples[2][2]['Ur']))/(12*h)
        errors['continuation_Ur_Z']=abs(nominal(out_continuation['Ur_Z'])-continuation_ur_z)/max(1,abs(continuation_ur_z))
        for name in ('F','Uz','P'):
            check('switches_'+name,out_switches[name+'_ZZ'],
                  {k:nominal(v[3][name]) for k,v in samples.items()})
        for name,second in out_switches['moments_ZZ'].items():
            check('switches_moment_'+name,second,{k:nominal(v[3]['moments'][name]) for k,v in samples.items()})
        for name,second in out_switches['raw_quadratic_integrals_ZZ'].items():
            check('switches_raw_'+name,second,{k:nominal(v[3]['raw_quadratic_integrals'][name]) for k,v in samples.items()})
        switch_ur_z=(nominal(samples[-2][3]['Ur'])-8*nominal(samples[-1][3]['Ur'])
                    +8*nominal(samples[1][3]['Ur'])-nominal(samples[2][3]['Ur']))/(12*h)
        errors['switches_Ur_Z']=abs(nominal(out_switches['Ur_Z'])-switch_ur_z)/max(1,abs(switch_ur_z))
        # Independent radial identity uses only scalar bridge moments and scalar
        # comparison finite-Z derivative, without second-jet propagation.
        def ur(k):
            raw=samples[k][1];zz=z+k*h;R=raw['R']
            mz=raw['moments']['z'];uz=raw['Uz']
            # Obtain Mz_Z by a scalar five-point first stencil with local shift.
            eps=mp.mpf('1e-16');local={}
            for j in (-2,-1,1,2):
                c=PressureWidthComparison(bundle(k*h+j*eps),h_b=width,pressure_order=2,width_order=1)
                b=PressureWidthExitBridge(c,steps=8)
                local[j]=b.evaluate('.8',zz+j*eps)['moments']['z']
            mz_z=(local[-2]-8*local[-1]+8*local[1]-local[2])/(12*eps)
            return nominal((2*zz*R*uz-(1-axis.delta)*zz*mz-(1-zz*zz)*mz_z)/((1-axis.delta*zz*zz)*(2*R).sqrt()))
        urvals={k:ur(k) for k in (-2,-1,1,2)}
        ur_z=(urvals[-2]-8*urvals[-1]+8*urvals[1]-urvals[2])/(12*h)
        errors['Ur_Z']=abs(nominal(out_bridge['Ur_Z'])-ur_z)/max(1,abs(ur_z))
        maximum=max(errors.values())
        print('max second-Z / radial-first-Z stencil error',mp.nstr(maximum,30),flush=True)
        assert maximum<mp.mpf('1e-24'),errors
        report=dict(maximum_scaled_error=mp.nstr(maximum,60),errors={k:mp.nstr(v,60) for k,v in errors.items()},
                    second_stencil_step='1e-8',radial_inner_stencil_step='1e-16',
                    oracle='shifted scalar finite core and scalar RK replay',
                    core_Z_depth_required=3,radial_second_Z_available=False,
                    actual_source_certified=False,uniform_Z_certified=False,ODE_error_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':run()
