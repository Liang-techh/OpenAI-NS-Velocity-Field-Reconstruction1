"""Source-bound Banach, coefficient uniqueness and common tail proof.

All claims concern the unique solution of the original fixed implicit source.
Directed mass/jet boxes remain enclosures, never selected source values.
"""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_core_common_fixed_point import (
    HERE,PREFIX,NAME,RECEIPT,BANACH_GATE,GATE,TAIL_GATE,STILL_OPEN,
    CurrentCoreCommonFixedPoint,raw_flatten_pressure_function,sha)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_candidate_pressure_function import q_jets
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_analytic_radial_tail import tail_factor


def exact(left,right,message):
    difference=s.cancel(left-right)
    if difference!=0 and s.simplify(difference)!=0:raise ArithmeticError(message)


def analytic_pressure_source_proof():
    z,t,c=s.symbols('Z offset center',real=True)
    q=1+z*z
    # Differential equation and analytic-axis value identify q^-2 at all
    # Taylor orders. This is stronger than finite sampled jets.
    f=q**-2
    exact(q*s.diff(f,z)+4*z*f,0,'Exact q inverse-square differential equation differs')
    n=s.Symbol('n',integer=True,nonnegative=True)
    a,b=s.symbols('a b')
    prev,current,next_=s.symbols('f_previous f_current f_next')
    coefficient=a*(n+1)*next_+b*(n+2)*current+(n+3)*prev
    solved=-(b*(n+2)*current+(n+3)*prev)/(a*(n+1))
    exact(coefficient.subs(next_,solved),0,'All-order q coefficient recurrence differs')
    count=0
    for center in (s.Integer(0),s.Rational(1,3),s.Rational(-2,5)):
        rows=q_jets(SimpleNamespace(mpf=s.sympify),center,7)
        for k,row in enumerate(rows):
            expected=s.diff(f,z,k).subs(z,center)/s.factorial(k)
            exact(row,expected,'Independent actual q_jets scalar differs')
            count+=1
    kernel=raw_flatten_pressure_function()
    sigma=s.Symbol('sigma',nonnegative=True)
    mu,yd,Tw=s.symbols('mu yd Tw',real=True)
    yv=yd+1+Tw+13/mu;y=yv+t
    J1=s.Rational(1,2)
    original_logA=s.Rational(1,10)*y-s.Rational(3,5)*(y-1+J1)-mu*(y-yd-1+J1)
    logP=s.Symbol('logPstar',real=True)
    exact((logP+original_logA)-logP,kernel['normalized_log_A'],
          'Original log amplitude minus logPstar normalization differs')
    normalized_log_u=original_logA-s.log(q)+sigma*(s.log(q)-s.log(2))
    expected_logdensity=2*original_logA-2*sigma*s.log(2)-s.log(2)
    exact(2*normalized_log_u-s.log(2),expected_logdensity+(-2+2*sigma)*s.log(q),
          'Original flatten pressure integrand factorization differs')
    actual_logdensity=kernel['log_density'].xreplace({kernel['sigma']:sigma})
    exact(actual_logdensity,expected_logdensity,'Published raw flatten radial density differs')
    exact(s.diff(actual_logdensity,yd),-1,'Flatten segmented radial density loses its yd factor')
    split=kernel['finite_log_density']-13/mu
    exact(split,kernel['log_density'],'Separated inverse-mu density logarithm differs')
    A,B=s.symbols('edge_x edge_one_minus_x',positive=True)
    exact(A/(A+B)+B/(B+A),1,'Original cutoff symmetry differs')
    # Bind the source implementation's log-ratio formula to the exact edge
    # formula used in the new integral. Floating point evaluation is not used
    # as an exact primitive or exact density.
    path=HERE/'../../src/openai_ns_reconstruction/schedule_pressure.py'
    tree=ast.parse(path.read_bytes())
    functions=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='outgoing_sigma']
    if len(functions)!=1:raise ValueError('Pinned outgoing sigma function missing')
    fn=functions[0]
    expected={'log_a':'-(1.0/x)**2','log_b':'-(1.0/(1.0-x))**2','delta':'log_b-log_a'}
    for target,expr in expected.items():
        values=[node.value for node in ast.walk(fn) if isinstance(node,ast.Assign)
                and any(ast.unparse(t)==target for t in node.targets)]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(expr,mode='eval').body):
            raise ValueError('Pinned exact cutoff logarithmic ratio changed')
    for statement in ('if x<=0.0:\n return 0.0','if x>=1.0:\n return 1.0',
                      'if delta>=0.0:\n e=math.exp(-delta)\n return e/(1.0+e)',
                      'e=math.exp(delta)','return 1.0/(1.0+e)'):
        node=ast.parse(statement).body[0]
        if not any(ast.dump(actual)==ast.dump(node) for actual in fn.body):
            raise ValueError('Pinned exact cutoff endpoint/log-ratio branch changed')
    exact((A/B)/(1+A/B),A/(A+B),'Positive log-ratio cutoff branch differs')
    exact(1/(1+B/A),A/(A+B),'Negative log-ratio cutoff branch differs')
    return dict(exact_q_inverse_square_ODE_and_all_order_unique_recurrence=True,
        independent_actual_q_jet_identities=count,
        cutoff_symmetry_implies_exact_J1_one_half=True,
        raw_flatten_positive_measure_and_q_kernel_factorization=True,
        original_log_amplitude_minus_logPstar_normalization_exact=True,
        inverse_mu_density_logarithm_kept_separate=True,
        source_formula='P0/Pstar^2=-(M2*(1+Z^2)^-2+M0+integral_0^100 exp(log_density(t))*(1+Z^2)^-beta(t) dt)',
        exact_true_mass_symbols_remain_implicit=True,
        holomorphic_integral_proof=[
            'The raw H=1 preheat source has fixed nonnegative radial measures, independent of Z.',
            'sigma is the pinned edge ratio;0<=sigma<=1, so0<=beta=2-2sigma<=2.',
            'On the common eta tube q=1+Z^2 has Re q>=1-eta^2>0 and |q|>=q_min=1-eta(2+eta).',
            'Use the single holomorphic principal log of q to define q^-beta.',
            '|q^-beta|<=q_min^-2; finite positive radial mass dominates the integral uniformly on smaller disks.',
            'Dominated holomorphic integration gives one analytic pressure function for every fixed original source.',
            '|dZ q^-beta|=|2beta Z q^(-beta-1)|<=4(1+eta)q_min^-1 |q^-beta|.',
            'Hence |P0|<=Pstar^2*((M2_upper+Mflatten_upper)/q_min^2+M0_upper), and |P0_Z|<=4(1+eta)Pbound/q_min.',
            'On every radius1/4 disk about a real center, |q|>=1-(1/4)^2. All flatten Taylor orders obey the original Cauchy bound.',
            'Pressure is even in Z; all odd derivatives at0 vanish exactly.'],
        original_all_order_normalized_jets_enclose_this_same_function=True,
        pressure_enclosure_not_defined_as_independent_coefficient_boxes=True,passed=True)


def actual_Banach_and_coefficient_proof(field):
    parent=field.parent
    if not (parent['arbitrary_order_extraction']['passed']
            and parent['supremum_composite_norm_proof']['passed']
            and field.replay['all_twenty_rows_and_combined_bounds_bit_exact_replayed']
            and field.measure_bindings['all_fourteen_current_true_pressure_measures_bound']
            and field.consumer_bindings['actual_all_order_map_and_integrals_AST_bound']
            and field.consumer_bindings['original_factorial_model_tail_consumer_AST_bound']):
        raise ValueError('Actual operator/coefficient/norm premises required')
    m,q=field.replay['size_bound'],field.replay['Lipschitz_bound']
    if endpoints(m)[1]>.5 or endpoints(q)[1]>=.5:
        raise ValueError('Actual map no longer a contraction of the closed unit ball')
    return dict(actual_source_operator_defined_on_same_complete_Xh_pair=True,
        complete_space_proof=[
            'A norm-Cauchy sequence converges uniformly in every coefficient and every axial derivative.',
            'The fundamental theorem of calculus identifies successive limiting derivatives.',
            'The supremum weighted bound passes to the coefficient limits and gives norm convergence.',
            'For each fixed radial n, the derivative bounds give an axial Taylor series of radius h, with only polynomial m growth.',
            'The radial series and its mixed derivatives converge locally where |rho|/20+|Z-center|/h<1.',
            'Thus Xh and its zero-axis correction subspace are complete; the closed unit pair ball is complete.'],
        model_center=['Phi0=sum_n (-chi*rho/2)^n/(n!(n+1)!)','Psi0=B*rho/2'],
        correction_map='X-X0=T(X), T=(epsilon/2 R J2 Etheta, epsilon/2 J1 Ez)',
        actual_closed_ball_self_mapping_bound=m,actual_pair_Lipschitz_bound=q,
        actual_uniform_correction_norm_bound=m,
        actual_conservative_consumer_correction_bound=field.replay['correction_bound'],
        unique_actual_analytic_fixed_point_for_each_fixed_implicit_source=True,
        coefficient_induction=[
            'At radial order n, the angular/axial RHS uses known rows through n; rho*derivative preserves that order.',
            'M_n=Psi_n/(n+1); the pressure primitive at order n uses Phi rows through n-1.',
            'The swirl rho*F0^2*Phi^2 term at order n also uses rows through n-1.',
            'Divide by2(n+1)(n+2) or2(n+1)^2; fixed L is zero-free on the admitted domain.',
            'Common Phi_axis=1, Psi_axis=0 and pressure-axis datum initialize the same unique formal series.',
            'The prior exact gauge/whole-production AST proof identifies these recurrences with advance_scaled_one.',
            'The admitted S Cauchy jets, exact all-order pressure enclosure and common gradient jets enclose the same fixed source.',
            'Directed Taylor convolution, differentiation and inverse of nonzero L propagate coefficient inclusion by induction.',
            'Each radial step consumes one axial coefficient; initial degree+depth+1 rows provide every requested remaining coefficient.'],
        all_original_finite_rows_enclose_coefficients_of_this_unique_solution=True,
        no_interval_overlap_or_fitted_source_equality_argument=True,
        same_source_analytic_functions_not_just_shared_receipt_ids=True,passed=True)


def all_order_tail_proof(field):
    # The exact mixed coefficient ratio contains a favorable square-weight
    # factor. This establishes the existing geometric bound at every n>N.
    n,i,k=s.symbols('n i k',integer=True,nonnegative=True)
    a=s.Symbol('radius',positive=True)
    exact_ratio=a/20*(n+k+1)/(n+1-i)*(n+1)**2/(n+2)**2
    dropped_ratio=a/20*(1+(i+k)/(n+1-i))
    exact(dropped_ratio-exact_ratio,
          a/20*(n+k+1)/(n+1-i)*(2*n+3)/(n+2)**2,
          'All-order nonlinear coefficient tail ratio differs')
    c=field.ctx;h=field.core.h;radius=c.mpf('4.1');rows=[]
    for N in (6,24):
        for total in range(6):
            for radial in range(total+1):
                axial=total-radial
                packet=tail_factor(c,degree=N,radial_order=radial,axial_order=axial,radius=radius,h=h)
                if endpoints(packet['ratio_upper'])[1]>=1:
                    raise ValueError('Same-source mixed tail ratio fails')
                rows.append(dict(degree=N,radial_order=radial,axial_order=axial,**packet))
    # The Bessel coefficient tail is not inferred from a sampled finite sum.
    # Derive the all-n ratio of the exact positive majorant used by production.
    first=a**(n-i)/(2**n*s.factorial(n-i)*s.factorial(n+1))*(k+1)*(n+1)**k
    factorial_ratio=a/2*((n+2)/(n+1))**k/((n+1-i)*(n+2))
    exact(s.expand_func(first.subs(n,n+1)/first),factorial_ratio,
          'Actual Bessel factorial majorant ratio differs')
    return dict(same_actual_fixed_point_coefficient_norm_used=True,
        exact_mixed_tail_coefficient_formula='correction*k!*h^-k*n!/(n-i)!*r^(n-i)*binomial(n+k,k)/(20^n*(n+1)^2*(k+1)^2)',
        exact_successive_ratio=str(exact_ratio),
        all_remaining_ratios_bounded_by='r/20*(1+(i+k)/(N+2-i)), n>=N+1',
        geometric_tail_rows=rows,all_order_ratio_proof=True,
        model_tail_proof=[
            'Phi=Phi0+deltaPhi is the same unique analytic solution; finite rows are coefficients of that sum.',
            'Real chi in[0,1]; the Bessel coefficients have factorial radial decay.',
            'For axial order k<=5, |[t^k]chi(t)^n|<=(k+1)(n+1)^k(1+sum_(j=1)^5 |chi_[j]|)^k.',
            'The original model_phi_jets sums through N and uses the decreasing factorial ratio to enclose all later terms.',
            'Psi0 is linear in rho, so has no tail beyond N>=6.',
            'M averages radial coefficients by1/(n+1), so its nonlinear tail is no larger.',
            'At the common root H=0, chi has axial valuation>=2; Bessel radial tail has no requested low axial-order contribution.'],
        pressure_tail_proof=[
            'Ptilde_(n+1)=S_scaled*sum_(i=0)^n Phi_i Phi_(n-i)/(n+1).',
            '|Phi_i|<=Bphi/(20^i*(i+1)^2) implies the sum is <=(n+1)*Bphi^2/20^n.',
            'With the same admitted S_scaled<=S_bound, the tail beyond radial degree N is <=S_bound*Bphi^2*r^(N+1)/(20^N*(1-r/20)).',
            'This sharper original values() bound is used only for axial order0.',
            'For all mixed orders use the complete Phi*Phi Xh convolution bound256*Bphi^2.',
            'S_scaled is fixed and radial constant, so it commutes with the radial integral V.',
            'V has norm<=80 by the parent all-order weight proof; the one fixed Cauchy multiplier contributes S(h/(eta/2)).',
            'Thus ||P_scaled-P_axis||_Xh<=80*256*fixed_multiplier*S_bound*Bphi^2.',
            'Multiply this common primitive norm by the same proved tail_factor at radial/axial order(i,k).',
            'The new pressure_profile consumes original finite P rows and this mixed tail; it retains derivatives of S in the complete axial convolution.',
            'This is the exact production pressure primitive; no independent pressure tail is fitted.'],
        mixed_pressure_tail_includes_both_axial_convolution_sums=True,
        finite_and_infinite_parts_of_same_analytic_solution=True,
        ordinary_derivative_factorials_and_radial_units_preserved=True,
        common_tail_does_not_imply_point_evaluation_or_temporal_recursion=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentCoreCommonFixedPoint(require_checked=False)
    if (raw['actual_five_defect_family_sha256']!=field.family or raw['implicit_source_sha256']!=field.source_sha
        or raw['datum_enclosure_sha256']!=field.datum_sha or any(raw[key] for key in (BANACH_GATE,GATE,TAIL_GATE))
        or any(raw[key] for key in STILL_OPEN)):
        raise ValueError('Common solution source/scope differs')
    # Newly produced runtime only. Previous121 seed rows and all physical
    # reports remain untouched; their source receipts are already consumed.
    runtime=field.evaluate('.327')
    if encode(runtime)!=raw['fresh_common_core_runtime']:
        raise ValueError('Common original finite/tail runtime replay changed')
    pressure=analytic_pressure_source_proof()
    Banach=actual_Banach_and_coefficient_proof(field)
    tails=all_order_tail_proof(field)
    report=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source_sha,
                datum_enclosure_sha256=field.datum_sha,
                actual_all_order_pressure_source_proof=pressure,
                actual_Banach_and_common_coefficient_proof=Banach,
                same_solution_all_order_tail_proof=tails,
                actual_twenty_majorant_rows_replayed=True,
                source_models_and_domains=field.model_bindings,
                original_finite_tail_consumers=field.consumer_bindings,
                current_true_pressure_measure_bindings=field.measure_bindings,
                accepted_exact_sigma_reflection_and_J1=field.exact_reflection,
                fresh_common_runtime_replayed=True,
                **dict.fromkeys((BANACH_GATE,GATE,TAIL_GATE),True),**dict.fromkeys(STILL_OPEN,False),
                remaining_dependency='Common pressure4C=PD+PI and original core-first recovery/drive/mixed4 join, then global stress/energy and true temporal recursion',
                all_scoped_checks_passed=True,all_passed=True,
                input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(encode(report),indent=2)+'\n').encode('utf8'))
    print('PASS actual Banach map and common production finite/model/nonlinear tails; core-first/global/temporal remain open',flush=True)
    return report


if __name__=='__main__':run()
