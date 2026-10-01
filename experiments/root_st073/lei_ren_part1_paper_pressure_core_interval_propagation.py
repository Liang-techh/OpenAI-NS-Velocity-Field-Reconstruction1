"""Pressure-only interval propagation through the actual finite radial core.

A pointwise coefficient calculation at Z=.3 with uniform pressure-error
inputs. Axis jets and delta are held at their stored MP values. No radial
series remainder, axis-data error, RK continuation or moment closure follows.
"""
import json, math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets
from lei_ren_part1_paper_core_ra_experiment import _axis_F0_taylor,_axis_u0_taylor
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_component_pressure_core import evaluate_component_core_coefficients
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def q_coefficients(iv,beta,z,degree):
    beta=iv.mpf(beta);z=iv.mpf(z);q=1+z*z
    values=[iv.exp(-beta*iv.ln(q))]
    if degree:
        values.append(-2*beta*z*values[0]/q)
    for n in range(1,degree):
        values.append(-(2*z*(n+beta)*values[n]+(n-1+2*beta)*values[n-1])/(q*(n+1)))
    return values


def core_moment_second_derivatives(coefficients,r):
    """Exact second axial derivatives of the five finite polynomial integrals."""
    f=coefficients['F'];u=coefficients['Uz'];zero=f[0][0]*0
    result={name:zero for name in ('theta','z','theta_z','z_theta','p')}
    def product_second(a,b):return 2*(a[2]*b[0]+a[1]*b[1]+a[0]*b[2])
    for n,row in enumerate(f):result['theta']+=4*row[2]*r**(n+2)/(n+2)
    for n,row in enumerate(u):result['z']+=2*row[2]*r**(n+1)/(n+1)
    for i,fi in enumerate(f):
        for j,uj in enumerate(u):
            n=i+j+2;result['theta_z']+=2*product_second(fi,uj)*r**n/n
        for j,fj in enumerate(f):
            n=i+j+1;value=product_second(fi,fj)
            result['p']+=value*r**n/n
            result['z_theta']-=value*r**(n+1)/(n+1)
    for i,ui in enumerate(u):
        for j,uj in enumerate(u):
            n=i+j+1;result['z_theta']+=product_second(ui,uj)*r**n/n
    return result


def run(radial_degree=18,jet_depth=2):
    precision=473;iv=MPIntervalContext();iv.dps=precision
    report=json.loads(Path(__file__).with_name('lei_ren_part1_paper_uniform_pressure_high_derivatives.json').read_text())
    length=radial_degree+jet_depth+1
    assert length-1<=report['maximum_derivative_order']
    profile=source_profile();datum=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(precision):
        center=mp.mpf('.3');axis=RegularCoreAxisJets(j='1e-14',Lambda='1e36',logC='5e151',delta='1e-200',precision=precision)
        finite=datum.taylor_components(0,center=0)['components']
        pressure=[iv.mpf(0) for _ in range(length)]
        for stage,row in finite.items():
            atoms=row['atoms'] if row['kind']=='q_power_atoms' else [dict(beta=row.get('beta',0),atom_at_Z0=row['value_at_Z0'])]
            for atom in atoms:
                coeff=q_coefficients(iv,atom['beta'],center,length-1)
                for k in range(length):pressure[k]-=iv.mpf(atom['atom_at_Z0'])*coeff[k]
        physical=iv.exp(iv.mpf(28));errors=[]
        for k,row in enumerate(report['derivatives'][:length]):
            size=iv.mpf(mp.make_mpf(tuple(row['normalized_Taylor_coefficient_error_upper']['exact_mpf_tuple'])))
            upper=endpoints(size)[1];errors.append(iv.mpf([-upper,upper]))
        uncertain=[physical*(p+e) for p,e in zip(pressure,errors)]
        nominal=[physical*p for p in pressure]
        f0=_axis_F0_taylor(axis,center,length);u0=_axis_u0_taylor(center,j=axis.j,length=length)
        kwargs=dict(F0_Z_taylor=f0,U0_Z_taylor=u0,radial_degree=radial_degree,precision=precision,scalar_converter=iv.mpf)
        true=core_coefficients(center,axis.delta,P0_Z_taylor=uncertain,**kwargs)
        fixed=core_coefficients(center,axis.delta,P0_Z_taylor=nominal,**kwargs)
        radius=iv.mpf('4')/iv.mpf('1e36');rows={};field={}
        for key in ('F','Uz','P'):
            parts=[];sums=[iv.mpf(0) for _ in range(jet_depth+1)]
            for n,(trow,frow) in enumerate(zip(true[key],fixed[key])):
                errors_n=[]
                for k in range(min(jet_depth+1,len(trow),len(frow))):
                    diff=trow[k]-frow[k];lo,hi=endpoints(diff);size=iv.mpf(max(abs(lo),abs(hi)))
                    weighted=math.factorial(k)*size*radius**n
                    sums[k]+=weighted
                    errors_n.append(dict(axial_order=k,coefficient_error_interval=diff,
                        error_contribution_at_R4_over_Lambda=weighted))
                parts.append(dict(radial_order=n,axial_errors=errors_n))
            rows[key]=parts;field[key]=dict(zip(('value','first_Z','second_Z'),[endpoints(v)[1] for v in sums]))
        true_fields=evaluate_component_core_coefficients(true,radius,iv.mpf(center),iv.mpf(axis.delta),square_root=iv.sqrt)
        fixed_fields=evaluate_component_core_coefficients(fixed,radius,iv.mpf(center),iv.mpf(axis.delta),square_root=iv.sqrt)
        true_fields['moments_ZZ']=core_moment_second_derivatives(true,radius)
        fixed_fields['moments_ZZ']=core_moment_second_derivatives(fixed,radius)
        moment_errors={}
        for name in ('theta','z','theta_z','z_theta','p'):
            moment_errors[name]={}
            for label,key in (('value','moments'),('first_Z','moments_Z'),('second_Z','moments_ZZ')):
                diff=true_fields[key][name]-fixed_fields[key][name]
                lo,hi=endpoints(diff)
                moment_errors[name][label]=dict(error_interval=diff,absolute_error_upper=max(abs(lo),abs(hi)))
        out=dict(Z=center,radial_degree=radial_degree,Z_jet_depth=jet_depth,precision=precision,
            field_coefficient_error_upper_bounds_at_R4_over_Lambda=field,coefficient_errors=rows,
            five_core_moment_pressure_error_bounds=moment_errors,core_moment_second_Z_error_enclosed=True,
            pressure_error_propagated=True,pressure_datum='same stored finite masses and exact q powers plus uniform analytic-integral error',
            axis_data_held_fixed=True,pointwise_core_only=True,uniform_core_error_enclosed=False,
            radial_series_remainder_enclosed=False,core_RK_error_enclosed=False,
            original_parameter_errors_enclosed=False,five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(out),indent=2)+'\n')
        for key in field:print(key,'pressure-driven finite-core value error upper',mp.nstr(field[key]['value'],16),flush=True)


if __name__=='__main__':run()
