"""Uniform pressure-driven perturbation of the degree18 finite core.

Shared axis data are enclosed analytically, and paired interval arithmetic
retains exactly zero differences for common uncertain nominal factors.
This does not bound infinite radial-series or RK continuation errors.
"""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_interval_difference import IntervalDifference
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_component_pressure_core import evaluate_component_core_coefficients
from lei_ren_part1_paper_pressure_core_interval_propagation import core_moment_second_derivatives
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def q_coefficients(iv,beta,z,degree):
    beta=iv.mpf(beta);q=1+z**2
    c=[iv.exp(-beta*iv.ln(q))]
    if degree:c.append(-2*beta*z*c[0]/q)
    for n in range(1,degree):c.append(-(2*z*(n+beta)*c[n]+(n-1+2*beta)*c[n-1])/(q*(n+1)))
    return c


def run(degree=18,depth=2,*,profile=None,pressure_receipt=None,output=None,profile_origin=None):
    precision=473;iv=MPIntervalContext();iv.dps=precision;length=degree+depth+1
    receipt_path=Path(pressure_receipt) if pressure_receipt else Path(__file__).with_name('lei_ren_part1_paper_uniform_pressure_high_derivatives.json')
    receipt=json.loads(receipt_path.read_text())
    assert length-1<=receipt['maximum_derivative_order']
    profile=source_profile() if profile is None else profile;datum=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(precision):
        j=mp.mpf('1e-14');lam=mp.mpf('1e36');logC=mp.mpf('5e151');delta=mp.mpf('1e-200')
        a=iv.mpf(receipt['radius']);z=iv.mpf([-endpoints(a)[1],endpoints(a)[1]])
        axis=uniform_axis_jets(iv,radius=a,j=j,Lambda=lam,logC=logC,delta=delta,length=length)
        finite=datum.taylor_components(0,center=0)['components'];pressure=[iv.mpf(0) for _ in range(length)]
        for row in finite.values():
            atoms=row['atoms'] if row['kind']=='q_power_atoms' else [dict(beta=row.get('beta',0),atom_at_Z0=row['value_at_Z0'])]
            for atom in atoms:
                c=q_coefficients(iv,atom['beta'],z,length-1)
                for k in range(length):pressure[k]-=iv.mpf(atom['atom_at_Z0'])*c[k]
        physical=iv.exp(iv.mpf(28));p0=[]
        for p,row in zip(pressure,receipt['derivatives']):
            size=mp.make_mpf(tuple(row['normalized_Taylor_coefficient_error_upper']['exact_mpf_tuple']))
            p0.append(IntervalDifference(iv,physical*p,physical*iv.mpf([-size,size])))
        converter=IntervalDifference.converter(iv)
        core=core_coefficients(z,delta,F0_Z_taylor=axis['F0'],U0_Z_taylor=axis['U0'],P0_Z_taylor=p0,
            radial_degree=degree,precision=precision,scalar_converter=converter)
        r=converter(iv.mpf(4)/iv.mpf(lam));zp=converter(z);dp=converter(delta)
        def sqrt_constant(value):
            lo,hi=endpoints(value.difference)
            if lo!=0 or hi!=0:raise ValueError('Square root only supplied for a common fixed radius')
            return converter(iv.sqrt(value.nominal))
        fields=evaluate_component_core_coefficients(core,r,zp,dp,square_root=sqrt_constant)
        second=core_moment_second_derivatives(core,r)
        fields['Uz_ZZ']=sum(2*row[2]*r**n for n,row in enumerate(core['Uz']))
        fields['F_ZZ']=sum(2*row[2]*r**n for n,row in enumerate(core['F']))
        fields['P_ZZ']=2*core['P'][0][2]+second['p']
        d=1-zp*zp;L=1-dp*zp*zp
        N=2*zp*r*fields['Uz']-(1-dp)*zp*fields['moments']['z']-d*fields['moments_Z']['z']
        NZ=2*r*fields['Uz']+2*zp*r*fields['Uz_Z']-(1-dp)*(fields['moments']['z']+zp*fields['moments_Z']['z'])+2*zp*fields['moments_Z']['z']-d*second['z']
        fields['Ur_Z']=(NZ+2*dp*zp*N/L)/(L*sqrt_constant(2*r))
        def error(value):
            lo,hi=endpoints(value.difference)
            return dict(error_interval=value.difference,absolute_error_upper=max(abs(lo),abs(hi)))
        field_errors={name:error(fields[name]) for name in ('F','F_Z','F_ZZ','Uz','Uz_Z','Uz_ZZ','Ur','Ur_Z','P','P_Z','P_ZZ')}
        moment_errors={name:dict(value=error(fields['moments'][name]),first_Z=error(fields['moments_Z'][name]),second_Z=error(second[name])) for name in second}
        out=dict(pressure_error_receipt=receipt_path.name,profile_origin=profile_origin or 'source_profile helper',radius=receipt['radius'],radial_degree=degree,Z_jet_depth=depth,precision=precision,
            R='4/Lambda',field_pressure_error_bounds=field_errors,five_core_moment_pressure_error_bounds=moment_errors,
            axis_G_absolute_upper=axis['G_absolute_upper'],axis_F0_interval=axis['F0_interval'],
            uniform_pressure_perturbation_enclosed=True,uniform_axis_jets_enclosed=True,
            shared_axis_data_difference_exactly_zero=True,parameter_scope='stored MP axis constants and stored schedule pressure datum',
            original_parameter_errors_enclosed=False,radial_series_remainder_enclosed=False,
            core_RK_error_enclosed=False,full_five_defect_interval_closure=False,temporal_recursion=False)
        (Path(output) if output else Path(__file__).with_suffix('.json')).write_text(json.dumps(encode(out),indent=2)+'\n')
        for name in ('Uz','Uz_Z','Uz_ZZ','Ur','Ur_Z'):
            print(name,'uniform pressure-driven finite core error',mp.nstr(field_errors[name]['absolute_error_upper'],16),flush=True)


if __name__=='__main__':run()
