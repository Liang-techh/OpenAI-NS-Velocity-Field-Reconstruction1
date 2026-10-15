"""Check numerical axial jets and own-history Q/pressure on the first collar."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from openai_ns_reconstruction.original_first_collar_jets import OriginalFirstCollarJets, bessel_jet_tail
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY = ROOT/'experiments/root_st073'
NAME = 'lei_ren_part1_paper_compliant_current_original_first_collar_jets.json'


def core_recovery_identity():
    z, e, de, j = sy.symbols('Z epsilon delta j')
    psi, mean = sy.Function('Psi')(z), sy.Function('Mpsi')(z)
    U0, d, L = 4*z+j, 1-z*z, 1-de*z*z
    V, M = U0+e*psi, U0+e*mean
    raw = (2*z*V-(1-de)*z*M-d*sy.diff(M,z))/L
    axis = ((1+de)*z*U0-4*d)/L
    correction = (2*z*psi-(1-de)*z*mean-d*sy.diff(mean,z))/L
    difference = sy.simplify(raw-axis-e*correction)
    for k in range(5):
        assert sy.diff(difference,z,k) == 0
    return 5


def bessel_checks(field):
    """Independent Bessel special-function derivatives vs interval tails."""
    c, iv, J = mp.mp.clone(), field.iv, field.Jet
    c.dps = 160
    frac = field.core.j_fraction
    j, de = c.mpf(frac.numerator)/frac.denominator, c.exp(-4*(c.exp(40)+11)-30)
    count = 0
    for Z in ('0', '.3'):
        def chi_exact(z):
            H = (1-de)*z/2+(1-z*z)*(4*z+j)
            return H*H/(H*H+(j/500)**2)
        z = J.variable(iv, iv.mpf(Z), 5)
        H = z*((1-field.core.delta_iv)/2)+(1-z*z)*(4*z+field.core.j_iv)
        h2 = J(iv, [H[0]**2]+list((H*H).coefficients[1:]))
        chi = h2/(h2+(field.core.j_iv/500)**2)
        lo, hi = _bounds(field.ctx, chi[0])
        chi = J(iv, [iv.mpf([max(0,lo),min(1,hi)])]+list(chi.coefficients[1:]))
        for i in range(3):
            def exact(z):
                a = chi_exact(z)
                q = 2*a
                return (-a/2)**i*c.besselj(i+1,2*c.sqrt(q))/q**(c.mpf(i+1)/2)
            reference = c.taylor(exact, c.mpf(Z), 5)
            for N in (8,16,24):
                finite = J.constant(iv,0,5)
                for n in range(i,N+1):
                    finite += chi**n*((-1)**n*iv.mpf(4)**(n-i)
                        /(iv.mpf(2)**n*math.factorial(n-i)*math.factorial(n+1)))
                error = bessel_jet_tail(iv,chi,N,radial_order=i)
                for k in range(6):
                    box = finite[k]+field._symmetric(error[k])
                    other = iv.mpf(c.nstr(reference[k],145))
                    a,b = _bounds(field.ctx,box);lo,hi = _bounds(field.ctx,other)
                    assert max(a,lo) <= min(b,hi), ('Bessel derivative',Z,i,N,k)
                    count += 1
    return count


def overlap(field, p, q):
    c, iv = field.ctx, field.iv
    def restore(v):
        return iv.mpf([mp.make_mpf(t) for t in v['source_enclosure_exact_mpf_tuples']])
    a,b = _bounds(c,restore(p));lo,hi = _bounds(c,restore(q))
    assert max(a,lo) <= min(b,hi)


def run():
    start = time.perf_counter()
    field = OriginalFirstCollarJets(precision=90)
    identities, tail_count = core_recovery_identity(), bessel_checks(field)
    samples, comparisons, scalar_matches, joins = [],0,0,0
    for Z in ('-.3','0','.3'):
        for phase in ('0','.5','1'):
            fine = field.evaluate(Z,phase,degree=24,order=5)
            coarse = field.evaluate(Z,phase,degree=18,order=5)
            scalar = super(type(field),field).evaluate(Z,phase,degree=24)
            pairs = []
            for name, rows in fine['normalized_fields_and_histories'].items():
                for kind in ('inlet','increment_per_original_hb'):
                    pairs.extend(zip(rows[kind],coarse['normalized_fields_and_histories'][name][kind]))
                for kind in ('inlet','increment_per_original_hb'):
                    overlap(field,rows[kind][0],scalar['normalized_fields_and_own_histories'][name][kind])
                    scalar_matches += 1
            for name, rows in fine['radial_recovery'].items():
                pairs.extend(zip(rows,coarse['radial_recovery'][name]))
            for name in ('P0_over_Pstar_squared','C_inlet_dressed_by_relative_F0_squared','C_increment_per_hb_dressed_by_relative_F0_squared'):
                pairs.extend(zip(fine['pressure'][name],coarse['pressure'][name]))
            for p,q in pairs:
                error = field.ctx.mpf(p['absolute_error_upper'])+field.ctx.mpf(q['absolute_error_upper'])
                assert abs(field.ctx.mpf(p['approximation'])-field.ctx.mpf(q['approximation'])) <= error
                comparisons += 1
            if phase == '0':
                for rows in fine['normalized_fields_and_histories'].values():
                    for row in rows['increment_per_original_hb']:
                        assert row['source_enclosure_exact_mpf_tuples'] == [[0,0,0,0],[0,0,0,0]]
                legacy = field._prepared[Z]['inputs']['p0']
                for k,row in enumerate(fine['pressure']['P0_over_Pstar_squared']):
                    lo,hi = _bounds(field.ctx,field._transfer(legacy[k]))
                    a,b = [field.ctx.make_mpf(t) for t in row['source_enclosure_exact_mpf_tuples']]
                    assert max(a,lo) <= min(b,hi), ('P0 source join',Z,k)
                    joins += 1
                legacy_C = field._prepared[Z]['data']['moments']['C']
                legacy_ratios = field._prepared[Z]['inputs']['F0_squared_ratios']
                for k,row in enumerate(fine['pressure']['C_inlet_dressed_by_relative_F0_squared']):
                    body = sum((legacy_C[l]*legacy_ratios[k-l]/math.factorial(k-l) for l in range(k+1)),legacy_C[0]*0)
                    lo,hi = _bounds(field.ctx,field._transfer(body))
                    a,b = [field.ctx.make_mpf(t) for t in row['source_enclosure_exact_mpf_tuples']]
                    assert max(a,lo) <= min(b,hi),('Original F0 squared pressure dress',Z,k)
                    joins += 1
            samples.append(fine)
        print('Computed actual first-collar Z5 jets, Q Z4 and pressure Z5: Z='+Z,flush=True)
    example = next(p for p in samples if p['Z']=='.3' and p['phase']=='.5')
    accuracy = {}
    for name, row in (('Qpsi_increment_value',example['radial_recovery']['Qpsi_increment_per_hb'][0]),
                      ('Phi_increment_Z5_coefficient',example['normalized_fields_and_histories']['Phi']['increment_per_original_hb'][5])):
        point = field.ctx.mpf(row['approximation'])
        relative = field.ctx.mpf(row['absolute_error_upper'])/abs(point)
        assert relative < field.ctx.mpf('1e-4'), (name,relative)
        accuracy[name] = field.ctx.nstr(relative,40)
    rejected = 0
    for args in (('.3','-.1',5),('1','.5',5),('.3','.5',6),('.3','.5',0)):
        try:
            field.evaluate(args[0],args[1],order=args[2])
        except (ValueError,ArithmeticError):
            rejected += 1
        else:
            raise AssertionError(('unsupported input accepted',args))
    code = dict(field.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    # Preserve distinct log proofs without repeating identical constant
    # denominator/product records at every phase and refinement level.
    proofs = list({json.dumps(p,sort_keys=True):p for p in field.scale_proofs}.values())
    return dict(source_family=field.family,
        checks=dict(passed=True,exact_Q_own_mean_split_derivative_identities=identities,
            independent_Bessel_special_function_mixed_tail_checks=tail_count,
            original_axial5_first_collar_packets=len(samples),Q_axial_order=4,pressure_axial_order=5,
            directed_degree18_degree24_comparisons=comparisons,
            independent_accepted_scalar_enclosure_overlaps=scalar_matches,
            compatible_original_pressure_source_coefficient_joins=joins,
            actual_zero_phase_increment_jets_exact=True, sampled_increment_relative_errors=accuracy,
            distinct_actual_width_log_proofs=len(proofs), rejected_unsupported_inputs=rejected,
            F0_gauge_canceled_before_width_difference_bounds=True,
            original_S_epsilon_squared_F0_squared_retained=True),
        samples=samples,normalized_source_width_log_proofs=proofs,input_hashes=field.hashes,
        implementation_hashes=code,
        limitations=dict(original_first_collar_leading_Z5_Q4_pressure5_available=True,
            absolute_F0_or_hb_materialized=False,finite_width_derivatives_enclosed=True,
            exact_original_point_oracle=False,first_second_derivative_join=False,
            n1_collar_numerical_transport=False,physical_Cartesian_collar_velocity=False,
            global_moment_repair=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-start,3))


def inspect(path):
    report = json.loads(path.read_text(encoding='utf-8'))
    assert report['checks']['passed']
    for name,digest in report['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result['checks'],indent=2),flush=True)
