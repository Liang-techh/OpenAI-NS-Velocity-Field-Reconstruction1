"""Accepted preheat-pressure complex modulus from true radial mass bounds.

Includes all14 stages, preserving P0 and physical Pstar^2 normalization.
Claims are relative to the accepted stored schedule, not parameter derivation.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    base=Path(__file__).parent
    names=['lei_ren_part1_paper_analytic_radial_tail.json',
           'lei_ren_part1_paper_coherent_uniform_fixed_beta_error.json',
           'lei_ren_part1_paper_coherent_uniform_preheat_error_budget.json']
    tube,fixed,complete=[json.loads((base/n).read_text()) for n in names]
    if not fixed['accepted_source_alignment_certified'] or not complete['accepted_source_alignment_certified']:
        raise AssertionError('Accepted schedule alignment missing')
    if fixed['accepted_schedule_sha256']!=complete['accepted_schedule_sha256']:
        raise AssertionError('Pressure source schedules differ')
    if len(fixed['stages'])!=13 or not complete['all_pressure_stages_included']:
        raise AssertionError('Pressure stages incomplete')
    ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        def interval(row):
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        def scalar(row):return ctx.mpf(mp.make_mpf(tuple(row['exact_mpf_tuple'])))
        eta=interval(tube['complex_tube_radius']);h=interval(tube['Xh_parameter'])
        qlower=interval(tube['pressure_q_modulus_lower'])
        if endpoints(qlower)[0]<=0:raise AssertionError('Pressure branch not protected')
        masses={};total=ctx.mpf(0)
        for stage,row in fixed['stages'].items():
            if row['beta'] not in (0,2):raise AssertionError('Unexpected fixed pressure exponent')
            mass=interval(row['true_mass_interval'])
            if endpoints(mass)[0]<0:raise AssertionError('Negative radial mass')
            mass_upper=ctx.mpf([0,endpoints(mass)[1]])
            modulus=mass_upper*qlower**(-row['beta'])
            masses[stage]=dict(beta=row['beta'],true_mass_upper=mass_upper,
                               normalized_complex_modulus_upper=modulus)
            total+=modulus
        flatten=scalar(complete['flatten_true_mass_upper'])
        # Flatten is a positive radial mixture of q^(-beta), beta in[0,2].
        # Its existing true-mass bound is radial and independent of Z.
        flatten_modulus=flatten*qlower**-2
        masses['z_flatten']=dict(beta_range=[0,2],true_mass_upper=flatten,
                                normalized_complex_modulus_upper=flatten_modulus)
        total+=flatten_modulus
        pressure=ctx.exp(28)*total
        # P0 holomorphic on the full capsule. At any point of the half
        # capsule a disk of radius eta/2 stays in the full capsule.
        pressure_derivative_cauchy=2*pressure/eta
        # Differentiate the same positive mass representation instead of
        # paying eta^-1. For beta in[0,2], |d_Z q^-beta| is bounded by
        # 4|Z| q_lower^-1 times the modulus majorant.
        pressure_derivative=4*(1+eta)*pressure/qlower
        a=1+eta/2;dt=ctx.mpf('1e-200');j=ctx.mpf('1e-14')
        Uupper=4*a+j;Hupper=interval(tube['H0_modulus_upper'])
        Llower=interval(tube['L_modulus_lower'])
        gupper=(1+dt)/2*(1+2*a*Uupper)*Uupper+4*Hupper
        gupper+=2*(1+dt)*a*pressure+(1+a*a)*pressure_derivative
        g_over_L=gupper/Llower
        ratio=h/(eta/2)
        if endpoints(ratio)[1]>mp.mpf('.4'):raise AssertionError('Cauchy weights too large')
        weight=ctx.mpf([1,max(mp.mpf(1),endpoints(4*ratio)[1])])
        gnorm=g_over_L*weight
        report=dict(input_hashes={n:hashlib.sha256((base/n).read_bytes()).hexdigest() for n in names},
             accepted_schedule_sha256=fixed['accepted_schedule_sha256'],precision=160,
             stage_count=len(masses),stages=masses,physical_pressure_normalization='Pstar^2=exp(28)',
             normalized_complex_pressure_modulus_upper=total,
             physical_P0_complex_modulus_upper=pressure,
             P0_first_derivative_modulus_upper_on_half_capsule=pressure_derivative,
             P0_first_derivative_Cauchy_alternative_upper=pressure_derivative_cauchy,
             pressure_derivative_from_same_radial_mass_representation=True,
             full_capsule_radius=eta,half_capsule_radius=eta/2,Xh_parameter=h,
             g_over_L_modulus_upper_on_half_capsule=g_over_L,
             g_over_L_Xh_norm_upper=gnorm,Psi_model_Xh_norm_upper=40*gnorm,
             pressure_complex_modulus_enclosed_relative_to_accepted_schedule=True,
             all_true_pressure_stages_included=True,
             fixed_pressure_real_derivative_errors_used_as_complex_errors=False,
             physical_P0_replaced=False,
             original_parameter_errors_enclosed=False,adapter_evaluation_roundoff_enclosed=False,
             core_source_error_enclosed=False,five_defect_interval_closure=False,
             nonlinear_Kstar_certified=False,
             actual_Lambda_contraction_certified=False,infinite_nonlinear_core_remainder_enclosed=False,
             temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('P0 complex modulus upper',mp.nstr(endpoints(pressure)[1],22),
              'Psi model norm upper',mp.nstr(endpoints(40*gnorm)[1],22),flush=True)
        return report


if __name__=='__main__':run()
