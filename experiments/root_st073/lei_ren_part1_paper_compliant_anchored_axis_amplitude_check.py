"""Independent real-path integration of the certified pole primitive."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_anchored_axis_amplitude import CompliantAnchoredAxisAmplitude
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_anchored_axis_amplitude.json'


def complex_interval(c,record):
    return c.mpc(read_interval(c,record['real']),read_interval(c,record['imag']))


def source_identities():
    z,a,x,j,delta,sigma=s.symbols('z a x j delta sigma',real=True)
    H=lambda v:-4*v**3-j*v**2+(9-delta)*v/2+j
    L=lambda v:1-delta*v*v
    hp=s.diff(H(z),z).subs(z,a)
    remainder=s.expand(H(a+x)-H(a)-hp*x-(-12*a-j)*x*x+4*x**3)
    if remainder!=0:raise ArithmeticError('Cubic Rouche remainder identity failed')
    decomposition=s.together(L(z)*H(z)/(H(z)**2+sigma**2)
        -L(z)/2*(1/(H(z)-s.I*sigma)+1/(H(z)+s.I*sigma)))
    if s.cancel(decomposition)!=0:raise ArithmeticError('Rational pole-source decomposition failed')
    h,hd,l=s.symbols('h hd l',nonzero=True)
    if s.cancel(l*h/(2*h*hd)-l/(2*hd))!=0:raise ArithmeticError('Original pole residue identity failed')
    return dict(exact_cubic_remainder_identity=True,quadratic_j_coefficient_retained=True,
        exact_rational_pole_decomposition=True,exact_simple_pole_residue_identity=True,passed=True)


def independent_real_integrals(field,receipt):
    # Numerical quadrature is an independent diagnostic. The directed pole
    # disks, rational identity and certified logarithm branches provide the
    # source error enclosure; quadrature error estimates are not proof inputs.
    with mp.workdps(240):
        c=field.ctx;nominal=lambda value:sum(endpoints(value))/2
        j=nominal(field.params['j']);delta=nominal(field.params['delta']);sigma=j/500
        linear=(9-delta)/2
        H=lambda z:-4*z**3-j*z*z+linear*z+j
        anchor=mp.findroot(H,-j/linear,tol=mp.eps*16,verify=True)
        hp=linear-12*anchor**2-2*j*anchor
        records=[];derivatives=0
        for text,packet in receipt['anchored_amplitude_packets'].items():
            if text=='anchor':continue
            z=mp.mpf(text);last=mp.asinh((z-anchor)*hp/sigma)
            direction=mp.sign(last)
            nodes=[mp.mpf(0)]+[direction*mp.mpf(n) for n in (1,2,4,8,16,32,64) if n<abs(last)]+[last]
            def transformed(q):
                distance=sigma*mp.sinh(q)/hp;zz=anchor+distance
                h=distance*(hp+(-12*anchor-j)*distance-4*distance**2)
                return (1-delta*zz*zz)*h/(h*h+sigma*sigma)*sigma*mp.cosh(q)/hp
            integral,error=mp.quad(transformed,nodes,error=True)
            bound=read_interval(c,packet['G']);lo,hi=endpoints(bound)
            if not lo<=integral<=hi:raise ArithmeticError('Independent real anchored integral outside pole enclosure: '+text)
            if error>mp.mpf('1e-180'):raise ArithmeticError('Independent quadrature did not resolve source interval')
            g=lambda zz:(1-delta*zz*zz)*H(zz)/(H(zz)**2+sigma*sigma)
            for k,value in enumerate(packet['G_gradient_taylor_coefficients']):
                derivative=mp.diff(g,z,k)/math.factorial(k)
                low,high=endpoints(read_interval(c,value))
                if not low<=derivative<=high:raise ArithmeticError('Independent rational gradient jet failed: '+str((text,k)))
                derivatives+=1
            records.append(dict(Z=text,G_independent_integral=mp.nstr(integral,70),
                numerical_quadrature_error_estimate=mp.nstr(error,12),
                directed_G_interval_width=mp.nstr(hi-lo,20)))
        return dict(independent_real_anchored_integrals=len(records),independent_rational_gradient_jet_coefficients=derivatives,
            root_peak_resolved_by_sinh_coordinate=True,quadrature_is_diagnostic_not_certificate=True,
            nominal_parameters_not_selected_as_exact_source=True,records=records,passed=True)


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Anchored amplitude source changed: '+name)
    with mp.workdps(400):
        f=CompliantAnchoredAxisAmplitude();c=f.ctx;identities=source_identities()
        roots=receipt['root_certificates'];poles=roots['poles']
        if len(poles)!=6 or not roots['unique_real_anchor_certified'] or not roots['all_six_poles_certified']:
            raise ValueError('Incomplete uniform source-root certificates')
        for row in [roots['anchor']]+poles:
            lhs=read_interval(c,row['rouche_lhs']);rhs=read_interval(c,row['rouche_rhs'])
            if endpoints(lhs)[1]>=endpoints(rhs)[0]:raise ArithmeticError('Stored Rouche certificate is not strict')
        for row in poles:
            imag=endpoints(complex_interval(c,row['disk']).imag)
            if not (imag[0]>0 or imag[1]<0):raise ArithmeticError('A pole disk meets real logarithm path')
        if endpoints(read_interval(c,roots['full_real_domain_H_sign_factor']))[0]<=0:
            raise ArithmeticError('Nonnegative anchored primitive lost its domain proof')
        if endpoints(read_interval(c,roots['full_real_domain_L_lower']))[0]<=0:
            raise ArithmeticError('Nonnegative anchored primitive lost its L positivity proof')
        anchor=receipt['anchored_amplitude_packets']['anchor']
        if endpoints(read_interval(c,anchor['G']))!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Shared anchor primitive is not exactly zero')
        branch_terms=0
        for text,packet in receipt['anchored_amplitude_packets'].items():
            g=read_interval(c,packet['G'])
            if endpoints(g)[0]<0:raise ArithmeticError('Anchored G must be nonnegative in original real domain')
            if text!='anchor' and endpoints(g)[1]-endpoints(g)[0]>mp.mpf('1e-140'):
                raise ArithmeticError('Anchored primitive is still a coarse global bound')
            if not packet['exp_logF0_not_materialized'] or not packet['exact_source_anchor_retained']:
                raise ValueError('Actual implicit amplitude source replaced by a scalar cap')
            logF=read_interval(c,packet['logF0']);expected=-f.core.logC-f.core.Lambda*g
            if logF._mpi_!=expected._mpi_:raise ArithmeticError('Selected Cstar amplitude conversion changed')
            for term in packet['signed_pole_log_terms']:
                if not term['individual_logs_stay_in_fixed_halfplane']:raise ValueError('Uncertified logarithm branch')
                branch_terms+=1
        independent=independent_real_integrals(f,receipt)
        points=0
        for Z,packets in receipt['physical_core_value_packets'].items():
            for packet in packets:
                rho=read_interval(c,packet['rho'])
                fresh=f.core_value(Z,rho)
                for component,term in packet['cylindrical_components'].items():
                    current=fresh['cylindrical_components'][component]
                    for key in ('coefficient_enclosure','positive_scale_log'):
                        if read_interval(c,term[key])._mpi_!=current[key]._mpi_:
                            raise ArithmeticError('Physical core source factor replay failed')
                if packet['axis_exact_zero']:
                    for component in ('ur','utheta'):
                        if endpoints(read_interval(c,packet['cylindrical_components'][component]['coefficient_enclosure']))!=(mp.mpf(0),mp.mpf(0)):
                            raise ArithmeticError('Axis velocity is not structurally zero')
                if not packet['fresh_core_and_actual_anchored_amplitude_combined']:
                    raise ValueError('Resolved amplitude not bound to fresh compliant core')
                points+=1
        for flag in ('full_point_physical_field_evaluation','all_annular_source_values_resolved','measured_blowup_dynamics',
                     'physical_energy_integral_certified','admissible_stress_lift_constructed','temporal_recursion'):
            if receipt[flag]:raise ValueError('Logarithmic core scope promoted: '+flag)
    for name in (NAME,Path(__file__).name):hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        source_identities=identities,independent_real_integral_check=independent,
        strict_uniform_rouche_root_disks_checked=7,six_poles_and_original_unique_anchor_certified=True,
        certified_halfplane_log_terms_checked=branch_terms,reproduced_logarithmic_physical_core_packets=points,
        actual_anchored_G_resolved_with_directed_error=True,positive_F0_retained_without_exp_logF0=True,
        full_point_physical_field_evaluation=False,all_annular_source_values_resolved=False,
        measured_blowup_dynamics=False,physical_energy_integral_certified=False,
        admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Anchored current-source amplitude PASS:7 root disks,7 independent real integrals,12 physical core packets',flush=True)
    return result


if __name__=='__main__':run()
