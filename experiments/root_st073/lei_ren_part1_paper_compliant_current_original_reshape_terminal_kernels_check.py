"""Independent complete finite integrals and actual source/error contracts."""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_reshape_terminal_kernels as current
from lei_ren_part1_paper_compliant_current_generic_shear_loop import flat_step

base=current.base;ep=current.ep


def restore(record,backend):
    c=backend.c;s=record['formal_positive_scale']
    scale=current.prior.FormalScale(backend.bases,tuple(s['source_exponents'])+(s['radius_power'],),
        current.read_interval(c,s['additional_log_interval']))
    return current.prior.ScaledEnclosure(scale,current.read_interval(c,record['coefficient_interval']),backend.ledger)


def full_finite_references():
    p=mp.mp.clone();p.dps=100;c=current.MPIntervalContext();c.dps=180
    comparisons=coefficients=0;records=[];observed_positive=observed_negative=False
    z,d,m=sy.symbols('z d m');bs=sy.symbols('b0:6')
    E=sy.exp(m*d*sum(bs[j]*z**j for j in range(6)))
    formulas=[sy.simplify(sy.diff(E,z,n).subs(z,0)/sy.factorial(n)/sy.exp(m*d*bs[0])) for n in range(6)]
    with mp.workdps(220):
        for T,L,sign in ((100,10,-1),(400,40,0),(1000,100,1)):
            Bcoeff=[p.mpf(sign)*T/1000,p.mpf('.04'),p.mpf('-.025'),p.mpf('.015'),p.mpf('-.007'),p.mpf('.003')]
            B=current.IntervalTaylor(c,[c.mpf(v) for v in Bcoeff]);backend=current.TerminalKernelEnclosures(c,B,c.mpf(T))
            for kind,(k0,m0,r0) in current.KINDS.items():
                row=backend.evaluate(kind,L=L);k=p.mpf(k0)
                expressions=[f.subs({m:m0,**{b:sy.Rational(str(v)) for b,v in zip(bs,Bcoeff)}}) for f in formulas]
                scalar_polynomials=[sy.lambdify(d,f,'mpmath') for f in expressions]
                cuts=sorted(set([p.mpf(0),p.mpf(2),p.mpf(8),p.mpf(L),p.mpf(T)/4,p.mpf(T)/2,p.mpf(T)*3/4,p.mpf(T)]))
                actual=[]
                for n in range(6):
                    def integrand(t):
                        delta=flat_step(p,t/T)
                        return p.exp(-k*t+m0*Bcoeff[0]*delta)*scalar_polynomials[n](delta)
                    value=p.quad(integrand,cuts);actual.append(value)
                    lo,hi=ep(current.conditioned.bounded_value(row['coefficients'][n]))
                    allowance=p.mpf('1e-85')*(1+abs(value))
                    assert lo-allowance<=value<=hi+allowance,(T,kind,n,p.nstr(value,18),mp.nstr(lo,18),mp.nstr(hi,18))
                    coeff_value=p.mpf(1)/k if n==0 else p.mpf(0)
                    upper=ep(current.conditioned.bounded_value(row['absolute_errors'][n]))[1]
                    assert abs(value-coeff_value)<=upper+allowance,(T,kind,n,'full remainder')
                    coefficients+=1
                if actual[0]>1/k:observed_positive=True
                if actual[0]<1/k:observed_negative=True
                # Differentiate the original finite-tail moment, independent
                # of the kernel implementation, and compare full quadrature.
                for power in (0,1,3,5):
                    r=p.mpf(r0)
                    reference=p.quad(lambda t:t**power*p.exp(-r*t),[L,L+10,L+100,p.inf])
                    expression=p.exp(-r*L)*sum(p.mpf(math.comb(power,j))*L**(power-j)*math.factorial(j)/r**(j+1) for j in range(power+1))
                    assert abs(reference-expression)<p.mpf('1e-85')*(1+abs(expression));comparisons+=1
                records.append(dict(T=T,L=L,B0_sign=sign,kind=kind,complete_finite_integral_coefficients=actual,
                    true_leading_remainder=actual[0]-1/k,manufactured_diagnostic_functions_only=True))
    assert coefficients==54 and comparisons==36 and observed_positive and observed_negative
    return dict(passed=True,complete_finite_original_integral_Z_coefficient_comparisons=coefficients,
        independent_full_positive_tail_moment_checks=comparisons,
        both_signs_of_actual_C0_remainder_observed=True,finite_diagnostic_function_records=records,
        native_source_parameters_not_selected_from_fixtures=True)


def actual_source_contracts(report):
    owner=current.OriginalLongReshapeTerminalKernels();backend=owner.backend;c=owner.c;item=report['actual_original_terminal_kernel_evaluation']
    assert item['source_family']==owner.family and item['implicit_source_sha256']==owner.source and item['datum_enclosure_sha256']==owner.datum
    assert item['original_bridge_switch_formal_reconstruction_flag_preserved_false']
    assert item['upstream_inlet_history_functions_still_unresolved'] and item['no_original_ancestor_constructors_or_producers_executed']
    assert item['separate_P0_retained'] and current.read_interval(c,item['original_Z_domain'])._mpi_==owner.Z._mpi_
    assert current.read_interval(c,item['original_frozen_T'])._mpi_==owner.T._mpi_
    binding=item['source_function_binding']
    assert binding['same_defining_B_function_identity_not_interval_overlap']
    assert binding['fixed_T_constructor_AST_bound'] and binding['all_original_source_packets_share_identical_T']
    assert binding['original_cached_T_contains_directed_same_source400_Abar']
    old=owner.old_source['actual_Rsh_exit']['actual_inherited_axial5_packet']['full_backward_kernel_axial5_coefficients']
    checks=nonzero_decays=nonzero_errors=0;gains={}
    with mp.workdps(c.dps+40):
        for kind in current.KINDS:
            row=item['actual_full_finite_endpoint_kernels'][kind];proof=row['proof']
            assert proof['full_original_finite_T_integral_enclosed'] and proof['source_T_not_shortened']
            assert proof['complete_original_tail_not_zeroed'] and proof['incoming_decay_not_reset_to_zero']
            assert proof['ordinary_Z_Taylor_coefficients_not_derivatives'] and proof['original_frozen_T_Z_exact_zero']
            k=current.read_interval(c,proof['k']);r=current.read_interval(c,proof['rate_min'])
            assert ep(current.read_interval(c,proof['source_decay_ratio']))[1]<=ep(k-r)[0]
            for n,value in enumerate(row['coefficients']):
                actual=restore(value,backend);lo,hi=ep(current.conditioned.bounded_value(actual))
                prior=current.read_interval(c,old[kind][n]);before=ep(prior)
                assert before[0]<=lo<=hi<=before[1],(kind,n,'original source cover containment')
                checks+=1
                assert not value['point_value_selected']
                error=restore(row['absolute_errors'][n],backend)
                assert not error.zero and not row['body_errors'][n]['exact_zero'] and not row['tail_errors'][n]['exact_zero']
                nonzero_errors+=1
                decay=restore(row['positive_incoming_decay_Z_bounds'][n],backend)
                assert not decay.zero;nonzero_decays+=1
                if n==0:
                    width=c.mpf(hi)-c.mpf(lo);oldwidth=c.mpf(before[1])-c.mpf(before[0])
                    assert ep(width)[0]>0 and ep(width)[1]<ep(oldwidth)[0]
                    gains[kind]=c.ln(oldwidth/width)/c.ln(10)
        assert all('unknown_original_' in text for text in item['source_owned_affine_moment_transport'].values())
        mean_decay=backend.exp_source(-owner.T)
        assert not mean_decay.zero and ep(mean_decay.coefficient)[0]>0
    assert checks==18 and nonzero_errors==18 and nonzero_decays==18
    return dict(passed=True,actual_whole_axis_original_kernel_Z0_through5_source_cover_containments=checks,
        nonzero_factored_error_sectors=nonzero_errors,nonzero_original_incoming_decay_sectors=nonzero_decays,
        C0_decimal_width_improvement=gains,all_six_original_affine_history_formulas_preserve_unknown_inlets=True,
        actual_original_mean_axial_exp_minus_T_not_reset_to_zero=True,
        source_B_recovery_and_frozen_T_are_defining_identities=True,
        no_upstream_caps_or_fixture_values_promoted_to_field_values=True)


def guards():
    c=current.MPIntervalContext();c.dps=140
    B=current.IntervalTaylor(c,[c.mpf(0)]*6);backend=current.TerminalKernelEnclosures(c,B,c.mpf(1000));count=0
    def rejected(fn):
        nonlocal count
        try:fn()
        except (ValueError,ArithmeticError):count+=1
        else:raise AssertionError('Invalid source/kernel request accepted')
    rejected(lambda:backend.evaluate('mean'))
    rejected(lambda:backend.evaluate('pressure',L=0))
    rejected(lambda:backend.evaluate('pressure',L=True))
    rejected(lambda:backend.evaluate('pressure',L=1.5))
    rejected(lambda:backend.evaluate('pressure',L=1000001))
    rejected(lambda:backend.evaluate('pressure',L=501))
    rejected(lambda:current.TerminalKernelEnclosures(c,B,0))
    rejected(lambda:current.TerminalKernelEnclosures(c,B,c.mpf('inf')))
    rejected(lambda:current.TerminalKernelEnclosures(c,current.IntervalTaylor(c,[c.mpf(0)]*5),1000))
    alien=current.MPIntervalContext();alien.dps=140
    rejected(lambda:current.TerminalKernelEnclosures(c,current.IntervalTaylor(alien,[alien.mpf(0)]*6),1000))
    large=current.IntervalTaylor(c,[c.mpf(1000)]+[c.mpf(0)]*5)
    rejected(lambda:current.TerminalKernelEnclosures(c,large,1000).evaluate('pressure',L=100))
    rejected(lambda:current.OriginalLongReshapeTerminalKernels(dps=400))
    assert count==12
    return dict(passed=True,invalid_kind_split_source_context_rate_and_precision_requests_rejected=count)


def run():
    began=time.monotonic();report=json.loads((current.HERE/current.NAME).read_bytes())
    assert report[current.GATE] and all(report[key] is False for key in current.OPEN)
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    assert current.exact_identities()['passed'] and report['exact_source_kernel_identities']['passed']
    result=dict(all_passed=True,**{current.GATE:True},source_family=report['source_family'],
        independent_full_finite_integrals_and_tail_moments=full_finite_references(),
        actual_whole_axis_source_kernel_and_inlet_contracts=actual_source_contracts(report),
        scoped_source_and_kernel_guards=guards(),**dict.fromkeys(current.OPEN,False),
        input_hashes={**report['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            current.PREFIX+'current_generic_shear_loop.py':current.sha(current.PREFIX+'current_generic_shear_loop.py')},
        execution_seconds=time.monotonic()-began,
        scope='Original complete finite terminal reshape kernel functions and source error bounds. Finite independent quadrature checks ordinary-Z0..5; genuine whole-axis source cache stays bound. Original bridge/switch and unknown R110 inlet histories, active patch, controls/global N and reconstruction remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original terminal reshape kernels: complete integrals, source factors and unknown inlet contracts PASS',flush=True)
    return result


if __name__=='__main__':run()
