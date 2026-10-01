"""Actual same-family functional five defects, not fitted point targets.

The early and swirl tails remain exact unknown smooth functions with proved
positive C1 enclosures. Their enormous logarithmic bounds are never
exponentiated or replaced by zero. The dominant signed axial contributions
are computed by directed integrals of the existing restoration cutoff.
"""
# Recomputed for the distinct compliant pressure source.
# Formula origin: lei_ren_part1_paper_shared_five_defect_admission.py; legacy source/receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_physical_norm_family import LogBounds
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def restoration_integrals(c,cells=2048):
    """Closed-cell rectangle enclosures; no midpoint cutoff substitution."""
    if not isinstance(cells,int) or cells<1:
        raise ValueError('Positive cell count required')
    total={n:c.mpf(0) for n in ('linear_1','linear_8_over_5','square_1')}
    for i in range(cells):
        left=c.mpf(i)/cells;right=c.mpf(i+1)/cells
        t=c.mpf([endpoints(left)[0],endpoints(right)[1]])
        phi=alpha_box(c,t+1,c.mpf(1))
        total['linear_1']+=c.exp(t)*phi/cells
        total['linear_8_over_5']+=c.exp(t*c.mpf('1.6'))*phi/cells
        total['square_1']+=c.exp(t)*phi**2/cells
    return total


def signed_jets(c,z,j,rho,logP,kernels,tail_cap):
    """Taylor enclosures of the ACTUAL defect functions at any real Z box.

    Each true function is defined by the same axis primitives and source
    field. Independent boxes enclose its tiny unresolved tails; they do not
    choose different functions, pressure data, or defect targets.
    """
    zl,zh=endpoints(z)
    if zl<-1 or zh>1:
        raise ValueError('Real axial domain[-1,1] required')
    rh=endpoints(rho)[1];th=endpoints(tail_cap)[1]
    E=IntervalTaylor(c,[j+c.mpf([-rh,rh]),c.mpf([-rh,rh])])
    q=IntervalTaylor(c,[1+z**2,2*z])
    invAm2=q*q*c.exp(c.mpf('1.2')-2*logP)
    tail=lambda:IntervalTaylor(c,[c.mpf([-th,th]),c.mpf([-th,th])])
    d=[E*kernels['d1']+tail(),E*kernels['d2']+tail(),tail(),
       invAm2*E*E*kernels['d4']+tail(),tail()]
    return d,invAm2


def run(cells=2048):
    names=['compliant_reference_join_bounds','compliant_reference_join_check',
           'compliant_K1_ledger','compliant_global_exit_certificate','compliant_physical_norm_family',
           'shared_bump_constants','compliant_core_transfer','shared_analytic_tube']
    records={n:json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in names}
    hashes={}
    for record in records.values():
        for n,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
                raise ValueError('Actual defect dependency changed: '+n)
            hashes[n]=digest
    join=records['compliant_reference_join_bounds'];exit=records['compliant_global_exit_certificate']
    ledger=records['compliant_K1_ledger'];norms=records['compliant_physical_norm_family']
    fixed=records['shared_bump_constants'];major=records['compliant_core_transfer']
    if (not join['same_family_R110_Rh_relaxed_cone_analytically_certified']
            or not join['actual_moments_continuously_inherited']
            or records['compliant_reference_join_check']['reference_join_family_sha256']!=join['reference_join_family_sha256']
            or join['admitted_inner_parameter_family_sha256']!=ledger['admitted_inner_parameter_family_sha256']
            or exit['admitted_inner_parameter_family_sha256']!=ledger['admitted_inner_parameter_family_sha256']
            or norms['uniform_Cstar_family_sha256']!=join['uniform_Cstar_family_sha256']):
        raise ValueError('Actual source join or independent identities missing')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        b=LogBounds(c);get=lambda r,n:read_interval(c,r[n])
        hi=lambda x:c.mpf(endpoints(x)[1]);lo=lambda x:c.mpf(endpoints(x)[0])
        A=get(norms,'A_upper');logC=get(norms,'selected_logCstar');logP=get(major,'logPstar')
        j=get(ledger,'required_j');rho=hi(get(ledger,'rho_core_C2_bound')+get(exit,'rho_bridge_C2_upper'))
        eta=hi(j+rho);e_star=get(fixed,'e_star')
        R0=c.mpf(110);lam=A+1
        # Conservative explicit C1 coefficients of the five centered early
        # integrals. The reference pressure primitive is integrated exactly
        # even though Fref^2 is singular but integrable at the axis.
        early=dict(z=2*R0*eta,
            mixed=4*R0**2*(lam+2)*eta,
            theta=4*(1+lam)*R0**2+(4*c.sqrt(2)/3)*R0**c.mpf('1.5'),
            energy_axial=3*R0*eta**2,
            energy_swirl=8*(1+2*lam)*R0**2+c.mpf(25)*R0/12,
            pressure=16*R0*(1+2*lam)+c.mpf('12.5'))
        # Exact defining scales, enclosed only when forming a norm bound.
        logRm=c.ln(R0)-6+10*(logC+logP)
        logBeta=c.ln(R0)-logRm
        logAlpha=400*A-10*(logC+logP)+6
        if endpoints(A)[0]<10 or endpoints(logAlpha)[1]>=-2:
            raise ArithmeticError('Required0<beta<alpha<exp(-2) topology failed')
        d0log={
            '1':b.hi(c.ln(early['z'])-logRm),
            '2':b.hi(c.ln(early['mixed'])-logC+c.ln(4)+c.mpf('.6')
                        -c.ln(2)/2-c.mpf('1.5')*logRm-logP),
            '3':b.hi(c.ln(early['theta'])-logC+c.ln(4)+c.mpf('.6')
                        -c.ln(2)/2-c.mpf('1.5')*logRm-logP),
            '4':b.add(b.hi(c.ln(early['energy_axial'])+c.ln(12)+c.mpf('1.2')-logRm-2*logP),
                      b.hi(c.ln(early['energy_swirl'])-2*logC+c.ln(12)+c.mpf('1.2')-logRm-2*logP)),
            '5':b.hi(c.ln(early['pressure'])-2*logC+c.ln(12)+c.mpf('1.2')-2*logP)}
        e0log=b.add(*d0log.values())
        alpha_term_log=b.hi(c.ln(16)+logAlpha/5)
        # Stronger than the sufficient thirds budget. These positive caps
        # make the complete norm numerically usable without exp(logCstar).
        e0cap=e_star/1000000;alphacap=e_star/1000000
        if (endpoints(e0log)[1]>=endpoints(c.ln(e0cap))[0]
                or endpoints(alpha_term_log)[1]>=endpoints(c.ln(alphacap))[0]
                or endpoints(eta)[1]>=endpoints(e_star/100)[0]):
            raise ArithmeticError('Complete five-defect preparation failed')
        invP2=c.exp(-2*logP)  # compact mpf exponent; unlike exp(Abar), finite to represent
        e_paper=hi(e0cap+alphacap+3*eta+40*eta**2*invP2)
        if endpoints(e_paper)[1]>=endpoints(e_star)[0]:
            raise ArithmeticError('Actual complete10.19 bound exceeds e_star')

        # Exact dominant signed contributions and the remaining actual tails.
        ints=restoration_integrals(c,cells)
        kernels=dict(d1=c.exp(-2)*(1+ints['linear_1']),
                     d2=c.exp(c.mpf('-3.2'))*(c.mpf('.625')+ints['linear_8_over_5']),
                     d4=c.exp(-2)*(1+ints['square_1']))
        tails={
            '1':b.hi(c.ln(2*R0*eta)-logRm),
            '2':b.add(d0log['2'],b.hi(c.ln(c.mpf('.625')*eta)+c.mpf('1.6')*logBeta),
                      b.hi(c.ln(2*eta)+c.mpf('1.6')*logAlpha)),
            '3':b.add(d0log['3'],b.hi(c.ln(2)+c.mpf('1.6')*logAlpha)),
            '4':b.add(b.hi(c.ln(2*R0*eta**2)+c.ln(12)+c.mpf('1.2')-logRm-2*logP),
                      b.hi(c.ln(early['energy_swirl'])-2*logC+c.ln(12)+c.mpf('1.2')-logRm-2*logP),
                      b.hi(c.mpf('1.2')*logAlpha)),
            '5':b.add(d0log['5'],b.hi(c.ln(10)+logAlpha/5))}
        logTailCap=-800*c.ln(10)-2*logP
        tail_cap=c.exp(logTailCap)
        if any(endpoints(v)[1]>=endpoints(logTailCap)[0] for v in tails.values()):
            raise ArithmeticError('Actual signed functional tail cap failed')
        e_signed=hi((kernels['d1']+kernels['d2'])*eta
                    +12*c.exp(c.mpf('1.2'))*invP2*kernels['d4']*eta**2+5*tail_cap)
        e_complete=c.mpf(min(endpoints(e_paper)[1],endpoints(e_signed)[1]))
        if endpoints(e_complete)[1]>=endpoints(e_star)[0]:
            raise ArithmeticError('Signed actual functional defect norm not admitted')
        definition=dict(reference_join_family_sha256=join['reference_join_family_sha256'],
            target='the five actual axis moment differences atRh=exp(-5)Rref; same P0',
            row_order=['Delta_z/Rm','(Delta_theta_z-4ZDelta_theta)/(sqrt2 Rm^1.5 Am)',
                       'Delta_theta/(sqrt2 Rm^1.5 Am)','(Delta_ztheta-8ZDelta_z)/(Rm Am^2)','Delta_p/Am^2'],
            Am='exp(-.6)*Pstar/(1+Z^2)',Rm='exp(-6)*110*(Cstar Pstar)^10',
            exact_dominant_functions=['k1*(v1-4Z)','k2*(v1-4Z)','0 dominant part; actual tail retained',
                                      'k4*Am^-2*(v1-4Z)^2','0 dominant part; actual tail retained'],
            kernels='exact integrals of the fixed flat restoration, enclosed by closed interval cells',
            tails='the exact differences between actual primitives and dominant contributions, not independent chosen targets',
            norm='sum_j (sup[-1,1]|d_j|+sup[-1,1]|d_j_Z|)',
            source_core_and_pressure_changed=False)
        sha=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        d,invAm2=signed_jets(c,c.mpf([-1,1]),j,rho,logP,kernels,tail_cap)
        signed_positive_floors={str(i+1):lo(d[i][0]) for i in (0,1,3)}
        if min(endpoints(v)[0] for v in signed_positive_floors.values())<=0:
            raise ArithmeticError('Signed actual d1,d2,d4 positivity lost to retained tails')
        result=dict(actual_five_defect_family_sha256=sha,definition=definition,
            reference_join_family_sha256=join['reference_join_family_sha256'],
            admitted_inner_parameter_family_sha256=ledger['admitted_inner_parameter_family_sha256'],
            uniform_Cstar_family_sha256=norms['uniform_Cstar_family_sha256'],
            implicit_source_sha256=major['implicit_source_sha256'],
            datum_enclosure_sha256=major['datum_enclosure_sha256'],
            early_centered_integral_coefficient_bounds=early,
            early_centered_integral_Cstar_powers=dict(z=0,mixed=-1,theta=-1,energy_axial=0,energy_swirl=-2,pressure=-2),
            inherited_d0_C1_log_upper_bounds=d0log,inherited_e0_C1_log_upper=e0log,
            log_alpha_enclosure=logAlpha,shape16_alpha_one_fifth_log_upper=alpha_term_log,
            inherited_e0_positive_cap=e0cap,shape16_alpha_one_fifth_positive_cap=alphacap,
            actual_v1_minus_4Z_minus_j_C1_upper=rho,actual_axial_eta_C1_upper=eta,
            required_j=j,logPstar=logP,inverse_Pstar_squared=invP2,
            e_star=e_star,complete_paper_sufficient_e_upper=e_paper,
            directed_restoration_integrals=ints,restoration_cells=cells,exact_signed_kernel_enclosures=kernels,
            actual_tail_C1_log_upper_bounds=tails,actual_tail_C1_positive_cap=tail_cap,
            actual_tail_cap_definition='10^-800/Pstar^2; positive and checked against all exact log bounds',
            actual_source_C1_norm_not_independent_Taylor_box_norm=True,
            independent_box_norm_can_be_larger_than_source_C1_bound=True,
            signed_actual_functional_e_upper=e_signed,complete_actual_functional_e_upper=e_complete,
            signed_d1_d2_d4_value_lower_bounds=signed_positive_floors,
            topology_beta_alpha_restore_certified=True,
            whole_axis_signed_actual_defect_Taylor_enclosures={str(i+1):list(row.coefficients) for i,row in enumerate(d)},
            whole_axis_inverse_Am_squared_Taylor_enclosure=list(invAm2.coefficients),
            complete_10_19_actual_functional_defect_test_certified=True,
            directed_fixed_cutoff_integrals_not_midpoint_fit=True,unknown_actual_tails_not_set_to_zero=True,
            actual_defect_functions_defined_by_source_primitives=True,source_core_pressure_or_Cstar_changed=False,
            five_moment_correction_constructed=False,finite_actual_defect_evaluator_complete=False,
            heat_exterior_matched=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
            proof=dict(
                early='C1 product norm, F<=4/Cstar, |logF_Z|<=Abar+1; integrate each centered integrand before applying inverse normalizers; ||q^-1||C1<=2,||q^-2||C1<=5,||q||C1=4,||q^2||C1=12 forq=1+Z^2; reference pressure primitive=5/(2Cstar^2q^2), reference energy half-integral norm<=25R0/(12Cstar^2)',
                complete='paper10.17-10.19 with T=400Abar, B C2<=2Abar; actual e0 andalpha log gates pluseta<e_star/100 imply the complete bound, without changing Cstar',
                signed='on the reference branch E=E1 untilxz=exp(-2), thenE=E1*(1-sigma(t)) untilexp(-1); integrate exact scalar kernels. Earlier field and reshaping differences are the five retained tails',
                tails='bound exact early centered integrals and all shape corrections by10.18; tail1 also subtracts110E1/Rm andtail4 subtracts110E1^2/(Rm Am^2); retained log bounds imply the positive10^-800/Pstar^2 cap',
                C1='E1-j has the inherited summed C2 boundrho; each actual tail has a C1 norm bound, so its ordinary value/derivative coefficients lie in the displayed boxes; actual smooth functions remain coupled by the same source',
                boxes='the true source satisfies sup|E1-j|+sup|E1_Z|<=rho and each actual tail has summed C1 norm<=tail_cap; independent Taylor boxes enlarge this correlated class and may have normj+2rho or2tail_cap. The analytic Banach e bound applies to the true source functions, not every arbitrary function in the enlarged independent boxes. The separate directed inverse encloses the enlarged boxes.',
                scope='signed function enclosures plus sufficient norm admission; no exact finite point defects, repaired exterior, ortemporal coefficients claimed'),
            input_hashes={**hashes,**{PREFIX+n+'.json':hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest() for n in names},
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Actual whole-axis five-defect10.19 admission PASS; Cstar unchanged',flush=True)
        print('Signed functional e upper=',mp.nstr(endpoints(e_complete)[1],14),
              'e/e_star<=',mp.nstr(endpoints(e_complete/lo(e_star))[1],14),flush=True)
        return result


if __name__=='__main__':
    run()
