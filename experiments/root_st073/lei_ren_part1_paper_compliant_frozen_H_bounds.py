"""Whole-axis full frozen-profile test9.14 via high/low chi partition.

Same preheat pressure, analytic core and nonzero implicit amplitude. This
proves a frozen inertial direction test, not the actual exit cone/stress lift.
"""
# Recomputed for the distinct compliant pressure source.
# Formula origin: lei_ren_part1_paper_shared_frozen_H_bounds.py; legacy source/receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum as LogarithmicPressureDatum
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    names=['lei_ren_part1_paper_compliant_frozen_angular_bounds.json',
           'lei_ren_part1_paper_compliant_core_uniform_bounds.json',
           'lei_ren_part1_paper_compliant_core_transfer.json',
           'lei_ren_part1_paper_shared_analytic_tube.json']
    angular,uniform,major,tube=[json.loads((HERE/n).read_bytes()) for n in names]
    for r in (angular,uniform,major,tube):
        if r['analytic_core_family_sha256']!=major['analytic_core_family_sha256']:
            raise ValueError('Frozen H core family mismatch')
        for n,d in r['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:
                raise ValueError('Frozen H dependency changed: '+n)
    if not angular['frozen_D_positive_full_axis_Ra_R110_certified']:
        raise ValueError('Positive frozen D required for AM-GM')
    datum=LogarithmicPressureDatum('40');c=datum.ctx
    if datum.source_sha!=major['implicit_source_sha256'] or datum.datum_sha!=major['datum_enclosure_sha256']:
        raise ValueError('Preheat pressure source mismatch')
    with mp.workdps(c.dps+40):
        get=lambda r,n:restore_value(c,r[n])
        lo=lambda x:c.mpf(endpoints(x)[0])
        hi=lambda x:c.mpf(endpoints(x)[1])
        eps=get(major,'epsilon');lam=get(major,'Lambda');j=get(major,'required_j')
        dt=hi(get(tube,'delta'));Lmin=lo(1-dt);sigma=lo(get(tube,'sigma'))
        h=get(tube,'Xh_parameter');radius=c.mpf('4.1')
        phiMin=lo(get(uniform,'actual_Phi_lower'));phiMax=hi(get(uniform,'actual_Phi_upper'))
        corr=get(major,'scaled_map_size_upper')
        psiNorm=get(major,'fresh_Psi_model_Xh_norm_upper')+corr
        psi=hi(psiNorm*embedding(c,h,radius,0,0));psiZ=hi(psiNorm*embedding(c,h,radius,0,1))
        deltaPhi=hi(corr*embedding(c,h,radius,0,0))
        deltaPhiR=hi(corr*embedding(c,h,radius,1,0))
        # High chi: q=2chi in [1.98,2], inside paper8.28 [1.9,2.1].
        # Its p4 majorant is decreasing and p4(1.9)<-.02.
        qa=c.mpf('1.9');qb=c.mpf('2.1')
        p4=lambda q:3-5*q/2+7*q*q/12-q**3/16+11*q**4/2880
        p4_derivative_upper=-c.mpf('2.5')+7*qb/6-3*qa**2/16+11*qb**3/720
        slopePolynomial=p4(qa)+3*deltaPhi+8*deltaPhiR
        if endpoints(p4_derivative_upper)[1]>=0 or endpoints(slopePolynomial)[1]>=0:
            raise ArithmeticError('High chi entry D>=3 not proved')
        square=uniform['signed_square_completion']
        M=hi(get(square,'model_chi_coefficient'))
        B=hi(get(square,'linear_sqrt_chi_B'))
        other=hi(get(square,'other_correction_upper'))
        S0high=lo(3-dt/2-j-eps*(psi+psiZ)+c.mpf('.99')*lam*Lmin-M-B-other)
        muBound=hi(get(angular,'normalized_entry_moment_error_bound'))
        highCoefficient=lo(4*eps*(S0high-muBound)/2)
        if endpoints(highCoefficient)[0]<=mp.mpf('1.9'):
            raise ArithmeticError('High chi frozen source coefficient too small')
        Dhigh=2*c.sqrt(c.mpf('1.9')*c.mpf('1.1'))

        # Low chi confines the axis to a small negative-Z interval.
        lowHfrac=hi(c.sqrt(99)/500)
        coarseCoefficient=lo(c.mpf('.5')-dt/2)
        if not (endpoints(lowHfrac)[1]<mp.mpf('.02') and
                endpoints(j)[1]<mp.mpf('.01') and
                endpoints(3*coarseCoefficient-1)[0]>endpoints(lowHfrac)[1]):
            raise ArithmeticError('Low chi coarse Z localization failed')
        ratioLower=lo((-1-c.mpf('.02'))/(c.mpf('4.5')-dt/2-36*j*j))
        ratioUpper=hi((c.mpf('.02')-1+9*j*j)/c.mpf('4.5'))
        if not (endpoints(ratioLower)[0]>-mp.mpf(1)/3 and
                endpoints(ratioUpper)[1]<-mp.mpf(1)/5):
            raise ArithmeticError('Low chi Z/j in [-1/3,-1/5] not proved')
        # |4Z+j|<=j/3 in this interval, throughout core and frozen branch.
        Vmax=hi(j/3+eps*psi);VZmax=hi(4+eps*psiZ)
        if not (endpoints(Vmax)[1]<endpoints(j)[0] and endpoints(VZmax)[1]<5):
            raise ArithmeticError('Low chi velocity bounds fail')
        Wlow=hi(1+j*Vmax/3+VZmax)
        averageLinear=hi((Vmax+j*VZmax/3)/2)
        if not (endpoints(Wlow)[1]<7 and endpoints(averageLinear)[1]<endpoints(2*j)[0]):
            raise ArithmeticError('Low chi cumulative average bounds fail')
        # Physical F0 stays implicit. Only its genuine complex guard upper
        # exp(-2logLambda-1000) is materialized, not the amplitude itself.
        F0upper=c.exp(-2*get(major,'logLambda')-1000)
        fullPhiZ=next(r for r in uniform['derivative_bounds'] if
                     r['scaled_radial_order']==0 and r['axial_order']==1)
        phiZ=hi(get(fullPhiZ,'Phi_absolute_bound'))
        Fsquare=hi(F0upper**2*phiMax**2)
        FsquareZ=hi(2*F0upper**2*(lam*phiMax**2/sigma+phiMax*phiZ))
        energyBound=hi(j*j+55*Fsquare)
        smallEnergyTerm=hi(2*dt*j*energyBound/3)
        if endpoints(smallEnergyTerm)[1]>=endpoints(j)[0]:
            raise ArithmeticError('Low chi energy term exceeds j')
        # All positive preheat atoms have exponent theta in [0,2]. ForZ<0,
        # both terms of P_operator P0 are nonnegative; six beta2 atoms
        # have m2>=2.5 and (1+Z^2)^2<=4, giving >=.25*j*Pstar^2.
        if endpoints(datum.m2)[0]<mp.mpf('2.5') or any(
                endpoints(row['mass'])[0]<0 for row in datum.stages.values()):
            raise ValueError('Positive same-source pressure atoms required')
        P2=c.exp(2*get(major,'logPstar'))
        pressureFloor=c.mpf('.25')*j*P2
        Nlow=lo(pressureFloor-20*j-440*Fsquare-165*FsquareZ)
        Ntarget=hi(c.mpf('.125')*j*P2)
        if endpoints(Nlow)[0]<=endpoints(Ntarget)[1]:
            raise ArithmeticError('Low chi frozen axial stress lower bound failed')
        fUpper=hi(F0upper*phiMax)
        Elow=Nlow*c.sqrt(2*eps)/fUpper
        Hlow=2*Elow
        threshold=c.mpf('2.04')
        if min(endpoints(Dhigh)[0],endpoints(Hlow)[0])<endpoints(threshold)[1]:
            raise ArithmeticError('Frozen H test9.14 unresolved')
        result=dict(
            analytic_core_family_sha256=major['analytic_core_family_sha256'],
            implicit_source_sha256=datum.source_sha,datum_enclosure_sha256=datum.datum_sha,
            real_axial_domain=['-1','1'],physical_radial_domain=['Ra=4/Lambda','110'],
            chi_partition=['chi>=.99','chi<=.99'],gamma='.01',Hf_threshold=threshold,
            high_chi=dict(entry_slope_polynomial_upper=slopePolynomial,
                polynomial_derivative_upper=p4_derivative_upper,S0_lower=S0high,
                frozen_source_coefficient_lower=highCoefficient,Df_lower=Dhigh,
                proof='Df>=3/x+1.9(x-1/x)>=2sqrt(1.9*1.1), x=R/Ra'),
            low_chi=dict(H0_over_j_absolute_upper=lowHfrac,
                Z_over_j_lower=ratioLower,Z_over_j_upper=ratioUpper,
                V_absolute_bound=Vmax,V_Z_absolute_bound=VZmax,
                Mz_over_R_absolute_bound=Vmax,Mz_over_R_Z_absolute_bound=VZmax,
                cumulative_average_proof='integral averages of core/frozen V and V_Z obey the same positive sup bounds',
                W_absolute_bound=Wlow,axial_linear_average_term_bound=averageLinear,
                energy_average_absolute_bound=energyBound,
                energy_average_Z_absolute_bound=10*j+55*FsquareZ,
                P_operator_pressure_increment_absolute_bound=110*(4*Fsquare+FsquareZ),
                F0_actual_upper=F0upper,F_squared_upper=Fsquare,F_squared_Z_upper=FsquareZ,
                P_operator_preheat_pressure_lower=pressureFloor,frozen_N_lower=Nlow,
                Ef_absolute_lower=Elow,Hf_lower=Hlow,
                proof='positive theta0..2 pressure atoms; bound all cumulative moments; Hf>=2|Ef| by Df>0'),
            frozen_D_positive_full_axis_Ra_R110_certified=True,
            frozen_D_at_least4_full_axis_R100_R110_certified=True,
            frozen_H_at_least2_plus4gamma_full_axis_Ra_R110_certified=True,
            full_frozen_profile_test_9_14_certified=True,
            finite_coefficient_state_used=False,actual_exit_cone_certified=False,
            admissible_stress_lift_constructed=False,full_Section9_parameter_admission=False,
            full_physical_C3_K_norms_certified=False,temporal_recursion=False,
            input_hashes={**angular['input_hashes'],**datum.input_hashes,
                **{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names+[Path(__file__).name]}})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Whole-axis full frozen9.14 PASS; high chi H lower',mp.nstr(endpoints(Dhigh)[0],14),
              'low chi H lower log',mp.nstr(mp.log(endpoints(Hlow)[0]),14))
        return result


if __name__=='__main__':run()
