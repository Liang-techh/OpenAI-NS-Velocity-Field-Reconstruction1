"""Global real-Z finite-core C3 contributions and axial exit tests.

This bounds degree18 finite polynomials; it does not certify the infinite
core or its reciprocal. The coherent analytic preheat datum is unchanged.
"""
import json,math,hashlib
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_uniform_core_pressure_propagation import q_coefficients
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_interval_difference import IntervalDifference
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def absolute_upper(value):
    if isinstance(value,IntervalDifference):value=value.value
    lo,hi=endpoints(value);return max(abs(lo),abs(hi))


def finite_norms(core,r,ctx):
    """C3 upper sums in x=R/Ra in [0,1], Z in [-1,1]."""
    out={}
    # The physical finite-model pressure restores the integral of the whole
    # degree18 F polynomial, not just the radial pressure recursion prefix.
    f=core['F'];pressure_rows=[core['P'][0]]
    for n in range(2*len(f)-1):
        row=[]
        for k in range(4):
            value=f[0][0]*0
            for i in range(max(0,n-len(f)+1),min(n,len(f)-1)+1):
                for a in range(k+1):value+=f[i][a]*f[n-i][k-a]
            row.append(value/(n+1))
        pressure_rows.append(row)
    for name in ('F','Uz','P'):
        rows=pressure_rows if name=='P' else core[name]
        terms={}
        for i in range(4):
            for k in range(4-i):
                total=ctx.mpf(0)
                for n,row in enumerate(rows):
                    if n<i:continue
                    factor=math.factorial(n)//math.factorial(n-i)*math.factorial(k)
                    total+=factor*ctx.mpf(absolute_upper(row[k]))*r**n
                terms[f'x{i}:Z{k}']=total
        norm=sum(terms.values(),ctx.mpf(0))
        out[name]=dict(mixed_derivative_upper_bounds=terms,C3_sum_upper=norm,log_C3_sum_upper=ctx.log(norm))
    return out


def run():
    base=Path(__file__).parent;precision=473;ctx=MPIntervalContext();ctx.dps=precision
    error_path=base/'lei_ren_part1_paper_global_pressure_high_derivatives.json'
    error=json.loads(error_path.read_text());profile,alignment=accepted_profile()
    datum=ContinuousPreheatPressure(profile,quadrature_order=192)
    with mp.workdps(precision+40):
        degree=18;length=degree+4;z=ctx.mpf([-1,1]);lam=ctx.mpf('1e36');dt=ctx.mpf('1e-200');j=ctx.mpf('1e-14')
        axis=uniform_axis_jets(ctx,radius=1,j=j,Lambda=lam,logC='5e151',delta=dt,length=length)
        print('global axis jets ready',flush=True)
        components=datum.taylor_components(0,center=0)['components']
        pressure=[ctx.mpf(0) for _ in range(length)]
        for name,row in components.items():
            expected=mp.make_mpf(tuple(alignment['stages'][name]['accepted_mass']['exact_mpf_tuple']))
            if row['value_at_Z0']!=expected:raise AssertionError('Accepted datum changed: '+name)
            atoms=row['atoms'] if row['kind']=='q_power_atoms' else [dict(beta=row.get('beta',0),atom_at_Z0=row['value_at_Z0'])]
            for atom in atoms:
                values=q_coefficients(ctx,atom['beta'],z,length-1)
                for k in range(length):pressure[k]-=ctx.mpf(atom['atom_at_Z0'])*values[k]
        physical=ctx.exp(ctx.mpf(28));p0=[]
        for nominal,row in zip(pressure,error['derivatives']):
            size=mp.make_mpf(tuple(row['normalized_Taylor_coefficient_error_upper']['exact_mpf_tuple']))
            p0.append(IntervalDifference(ctx,physical*nominal,physical*ctx.mpf([-size,size])))
        convert=IntervalDifference.converter(ctx)
        core=core_coefficients(z,'1e-200',F0_Z_taylor=axis['F0'],U0_Z_taylor=axis['U0'],P0_Z_taylor=p0,
                               radial_degree=degree,precision=precision,scalar_converter=convert)
        print('global finite recurrence ready',flush=True)
        r=ctx.mpf(4)/lam;norms=finite_norms(core,r,ctx)
        # Subtract the known affine axis analytically, retaining correlated Z.
        # This avoids [-4,4]-[-4,4] interval cancellation in the exit tests.
        deviations={};moments={}
        for k in range(3):
            leading=convert(j if k==0 else 0)
            velocity=leading+sum((math.factorial(k)*row[k]*r**n for n,row in enumerate(core['Uz']) if n),convert(0))
            moment=leading+sum((math.factorial(k)*row[k]*r**n/(n+1) for n,row in enumerate(core['Uz']) if n),convert(0))
            deviations[str(k)]=dict(interval=velocity.value,absolute_upper=absolute_upper(velocity))
            moments[str(k)]=dict(interval=moment.value,absolute_upper=absolute_upper(moment))
        total=sum((ctx.mpf(deviations[str(k)]['absolute_upper'])+ctx.mpf(moments[str(k)]['absolute_upper']) for k in range(3)),ctx.mpf(0))
        from lei_ren_part1_paper_collar_interval_receipt_reader import read_inlet
        inlet,inlet_source=read_inlet();checks=0
        z0=inlet['z'].value.evaluate(width=0).value;r0=inlet['R'].value.evaluate(width=0).value
        for k,slot in enumerate(('value','tangent','second')):
            affine=4*z0 if k==0 else ctx.mpf(4 if k==1 else 0)
            U=getattr(inlet['Uz'],slot).evaluate(width=0).value-affine
            M=getattr(inlet['moments']['z'],slot).evaluate(width=0).value/r0-affine
            for value,box in ((U,deviations[str(k)]['interval']),(M,moments[str(k)]['interval'])):
                lo,hi=endpoints(value);bl,bh=endpoints(box)
                if not bl<=lo<=hi<=bh:raise AssertionError('Independent actual inlet outside global axial bound')
                checks+=1
            for name in ('F','Uz','P'):
                value=getattr(inlet[name],slot).evaluate(width=0).value
                if absolute_upper(value)>endpoints(norms[name]['mixed_derivative_upper_bounds'][f'x0:Z{k}'])[1]:
                    raise AssertionError('Independent actual inlet outside global field norm')
                checks+=1
        # Datum C3 includes the axis pressure before radial integration.
        P0norm=sum((ctx.mpf(absolute_upper(p0[k]))*math.factorial(k) for k in range(4)),ctx.mpf(0))
        report=dict(radius='1',radial_degree=degree,Z_jet_depth=3,precision=precision,
            normalized_radial_domain=['0','1'],radial_coordinate='x=R/R_a, R_a=4/Lambda',
            accepted_schedule_sha256=alignment['accepted_schedule']['sha256'],
            pressure_error_receipt=error_path.name,pressure_error_sha256=hashlib.sha256(error_path.read_bytes()).hexdigest(),
            core_C3_contributions=norms,axis_pressure_C3_upper=P0norm,
            exit_Uz_minus_4Z_derivative_bounds=deviations,exit_mz_over_Ra_minus_4Z_derivative_bounds=moments,
            combined_exit_C2_deviation_upper=total,
            independent_actual_inlet_receipt=inlet_source,independent_actual_inlet_checks_passed=checks,
            exact_finite_F_squared_pressure_integral_restored=True,
            degree18_finite_core_enclosed=True,whole_real_axial_domain_enclosed=True,
            pressure_integral_error_included=True,reciprocal_core_F_bound_certified=False,
            full_K_bound_certified=False,frozen_profile_cone_bounds_certified=False,
            infinite_radial_remainder_enclosed=False,source_parameter_derivation_certified=False,
            actual_finite_width_cone_certified=False,temporal_recursion=False)
        (base/(Path(__file__).stem+'.json')).write_text(json.dumps(encode(report),indent=2)+'\n')
        print('whole-Z finite exit C2 deviation',mp.nstr(endpoints(total)[1],18),flush=True)
        for name,data in norms.items():print(name,'log C3 contribution',mp.nstr(endpoints(data['log_C3_sum_upper'])[1],18),flush=True)


if __name__=='__main__':run()
