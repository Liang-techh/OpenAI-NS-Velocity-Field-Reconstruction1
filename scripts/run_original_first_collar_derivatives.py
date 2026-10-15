"""Evaluate actual radial derivatives and check the original core source seam."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_first_collar_derivatives import OriginalFirstCollarDerivatives
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_first_collar_derivatives.json'


def exact_identities():
    x,r,h,sp=sy.symbols('x rho hb sigma_prime',positive=True)
    F=sy.Function('F')(x);rhs=sy.Symbol('rhs')
    radial=lambda expr:sy.diff(expr,x)/(2*x)
    count=0
    for damping,head in ((3,4),(1,2)):
        second=radial(radial(F)).subs(sy.diff(F,x,2),rhs-damping*sy.diff(F,x)/x)
        assert sy.simplify(second-(rhs-head*sy.diff(F,x)/x)/(4*x*x))==0;count+=1
    u,k=sy.Function('u')(r),sy.Function('k')(r)
    k1=-sy.diff(u,r)-k/r
    k2=sy.diff(k1,r).subs(sy.diff(k,r),k1)
    assert sy.simplify(k2-(-sy.diff(u,r,2)+sy.diff(u,r)/r+2*k/r**2))==0;count+=1
    z,de,ur,kv,kz=sy.symbols('Z delta ur k kZ');d,L=1-z*z,1-de*z*z
    raw=((1-de)*z*ur-(1+de)*z*(-ur-kv/r)-d*(-kz/r))/L
    expected=(2*z*ur+(1+de)*z*kv/r+d*kz/r)/L
    assert sy.simplify(raw-expected)==0;count+=1
    assert sy.simplify(-(1-h)*sp/(h*r)-(-sp/(h*r)+sp/r))==0;count+=1
    phi,pr,prr=sy.symbols('Phi Phi_rho Phi_rhorho',nonzero=True)
    A,Ar=pr/phi,prr/phi-(pr/phi)**2
    assert sy.simplify(phi*(A*A+Ar)-prr)==0;count+=1
    # Own-mean identity, including its differentiated k relation.
    assert sy.simplify(sy.diff(u,r)+k1+k/r)==0;count+=1
    return count


def compare(a,b):
    c=mp.mp.clone();c.dps=180
    diff=abs(c.mpf(a['approximation'])-c.mpf(b['approximation']))
    error=c.mpf(a['absolute_error_upper'])+c.mpf(b['absolute_error_upper'])
    assert diff<=error,(diff,error)


def overlaps(field,a,b):
    left,right=_bounds(field.ctx,a);lo,hi=_bounds(field.ctx,b)
    assert max(left,lo)<=min(right,hi)


def seam_checks(field,Z):
    packet,data=field._prepare(Z,'0')
    true=data['truth'];leading,_,forcing,_=true['derivatives']
    core=field.transport.n1;iv=core.iv
    source,raw=core._forcing_polynomial(Z,18,True)
    errors,_=core._source_tail(raw,Z,'4',18)
    count=0
    for name,target in (('f',forcing['fx']),('u',forcing['ux']),
                        ('p',forcing['p']/(2*true['inputs']['x']))):
        value=sum((row[0]*iv.mpf(4)**n for n,row in enumerate(source[name].rows)),iv.mpf(0))
        bound=_bounds(core.ctx,errors[name][0])[1]
        enclosure=field.leading._transfer(value+iv.mpf([-bound,bound]))
        overlaps(field,target.rows[0][0],enclosure);count+=1
        assert target.rows[-1][0]._mpi_[0]==target.rows[-1][0]._mpi_[1]==(0,0,0,0)
    inlet=field.leading.inlet_jets(Z,degree=24,order=4)['truth']['radial']
    for name in ('Phi','Psi'):
        for i,suffix in ((1,'_rho'),(2,'_rhorho')):
            for k in (0,1):
                overlaps(field,leading[name+suffix].rows[0][k],inlet[name][i][k]);count+=1
    return count


def run():
    started=time.perf_counter();identities=exact_identities()
    field=OriginalFirstCollarDerivatives(precision=85)
    samples=[];seams=zeros=0
    for Z in ('-.3','0','.3'):
        for phase in ('0','.25','.5'):
            packet=field.evaluate(Z,phase)
            if phase=='0':
                seams+=seam_checks(field,Z)
                for group in ('leading_radial_derivatives','n1_radial_derivatives'):
                    for rows in packet[group].values():
                        for row in rows['original_hb_powers'].get('-1',[]):
                            assert row['source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]]
                            zeros+=1
            samples.append(packet)
        print('Actual first-collar radial derivatives and source seam: Z='+Z,flush=True)
    samples.append(field.evaluate('.17','.37'))
    coarse=OriginalFirstCollarDerivatives(precision=60);comparisons=0
    for Z in ('-.3','0','.3'):
        old=coarse.evaluate(Z,'.5')
        new=next(p for p in samples if p['Z']==Z and p['phase']=='.5')
        for group in ('leading_radial_derivatives','n1_radial_derivatives'):
            for name,rows in new[group].items():
                previous=old[group][name]['original_hb_powers']
                for power,coefficients in rows['original_hb_powers'].items():
                    for a,b in zip(coefficients,previous[power]):
                        compare(a,b);comparisons+=1
    example=next(p for p in samples if p['Z']=='.3' and p['phase']=='.5')
    accuracy={}
    for group,name in (('leading_radial_derivatives','Phi_rhorho'),
                       ('n1_radial_derivatives','f_rhorho'),
                       ('n1_radial_derivatives','u_rhorho'),
                       ('n1_radial_derivatives','p_rho')):
        row=example[group][name]['original_hb_powers']['-1'][0]
        point=field.ctx.mpf(row['approximation'])
        relative=field.ctx.mpf(row['absolute_error_upper'])/abs(point)
        assert relative<field.ctx.mpf('1e-15'),(name,relative)
        accuracy[name]=field.ctx.nstr(relative,35)
    assert field.ctx.mpf(example['n1_radial_derivatives']['p_rho']['original_hb_powers']['-1'][0]['approximation'])<-10
    rejected=0
    for Z,phase in (('1','.5'),('.3','.5001'),('.3','-.1'),('.3','nan')):
        try:field.evaluate(Z,phase)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError(('Unsupported derivative request admitted',Z,phase))
    code=dict(field.code_hashes)
    code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    proofs=list({json.dumps(p,sort_keys=True):p for p in field.leading.scale_proofs}.values())
    return dict(source_family=field.family,
        checks=dict(passed=True,exact_chain_rule_kinematic_recovery_and_width_split_identities=identities,
            original_radial_derivative_packets=len(samples),
            independent_core_known_source_and_leading_radial_seam_comparisons=seams,
            exactly_zero_inverse_width_seam_coefficients=zeros,
            cross_precision_complete_coefficient_error_comparisons=comparisons,
            sample_inverse_width_coefficient_relative_errors=accuracy,
            rejected_unsupported_requests=rejected,all_width_powers_kept_factored=True,
            second_radial_values_use_full_local_equation=True,
            finite_radial_unknown_polynomial_not_differentiated=True,
            distinct_original_width_log_proofs=len(proofs)),
        samples=samples,original_width_log_proofs=proofs,
        input_hashes=field.hashes,implementation_hashes=code,
        limitations=dict(native_radial_derivatives_only=True,first_phase_max='.5',
            n1_second_rho_Z_order=0,n1_q_rho_Z_order=0,n1_pressure_rho_Z_order=0,
            complete_Cartesian_C2_velocity_join=False,pressure_second_rho=False,
            original_hb_inverse_not_materialized=True,remaining_collars=False,
            global_moment_repair=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-started,3))


def inspect(path):
    report=json.loads(path.read_text(encoding='utf-8'));assert report['checks']['passed']
    for name,digest in report['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    report=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds']),indent=2),flush=True)
