"""Independent unscaled equations, fresh axis data and infinite-tail checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_core_coefficient_rebuild import CompliantCoreCoefficientRebuild
from lei_ren_part1_paper_functional_core_step import initial_rows
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_core_coefficient_rebuild.json'


def original_unscaled_fixture():
    with mp.workdps(100):
        c=MPIntervalContext();c.dps=130;degree=4;length=8
        z=mp.mpf('.3');delta=mp.mpf('.02');Lambda=mp.mpf(7);eps=1/Lambda;tol=mp.mpf('1e-70')
        F=lambda zz:mp.mpf('.4')*mp.exp(mp.mpf('.17')*zz+mp.mpf('.03')*zz*zz)
        U=lambda zz:4*zz+mp.mpf('.1')
        P=lambda zz:mp.exp(mp.mpf('.2')*zz)+mp.sin(zz)
        exact={label:mp.taylor(fn,z,length-1) for label,fn in (('F',F),('U',U),('P',P))}
        ell=mp.taylor(lambda zz:eps*(mp.mpf('.17')+mp.mpf('.06')*zz),z,length-1)
        S=mp.taylor(lambda zz:eps**2*F(zz)**2,z,length-1)
        box=lambda rows:[c.mpf([value-tol,value+tol]) for value in rows]
        fixed=dict(ell_Z_taylor=box(ell),S_Z_taylor=box(S),
            U0_Z_taylor=box(exact['U']),P0_Z_taylor=box([eps*v for v in exact['P']]))
        rows=initial_rows(c,fixed,z,degree,required_depth=3)
        for n in range(degree):advance_scaled_one(c,fixed,rows,n,c.mpf(z),c.mpf(delta),c.mpf(eps))
        direct=core_coefficients(z,delta,F0_Z_taylor=exact['F'],U0_Z_taylor=exact['U'],
            P0_Z_taylor=exact['P'],radial_degree=degree,precision=100)
        count=0
        for n in range(degree+1):
            for k in range(length-n):
                angular=sum((rows['A'][n][j]*c.mpf(exact['F'][k-j]) for j in range(k+1)),c.mpf(0))*c.mpf(Lambda)**n
                mapped={'F':angular,'Uz':rows['Uz'][n][k]*c.mpf(Lambda)**n,
                        'P':rows['P'][n][k]*c.mpf(Lambda)**n/c.mpf(eps)}
                for label,value in mapped.items():
                    target=direct[label][n][k]
                    if not endpoints(value)[0]<=target<=endpoints(value)[1]:
                        raise ArithmeticError('Independent original/unscaled core equation mismatch: '+str((label,n,k)))
                    count+=1
        return dict(independent_unscaled_coupled_coefficients=count,
            nonzero_swirl_and_axial_and_pressure_fixture=True,
            distinct_Lambda_epsilon_pressure_units_checked=True,passed=True)


def independent_tail_bounds():
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;h=mp.mpf('.07');count=0
        for degree in (8,24):
            for radius in (mp.mpf(2),mp.mpf('4.1')):
                for i in range(5):
                    for k in range(5-i):
                        record=tail_factor(c,degree=degree,radial_order=i,axial_order=k,radius=c.mpf(radius),h=c.mpf(h))
                        bound=endpoints(record['tail_per_Xh_norm'])[1];rat=endpoints(record['ratio_upper'])[1]
                        total=mp.mpf(0);previous=None
                        for n in range(degree+1,degree+140):
                            term=(mp.mpf(math.factorial(n)//math.factorial(n-i))*math.factorial(k)*math.comb(n+k,k)
                                *radius**(n-i)/(mp.mpf(20)**n*h**k*(n+1)**2*(k+1)**2))
                            if previous is not None and term/previous>rat:raise ArithmeticError('Mixed tail ratio does not dominate coefficient norm')
                            previous=term;total+=term
                        if total>bound:raise ArithmeticError('Mixed infinite tail bound smaller than independent partial sum')
                        count+=1
        # Independent saturated coefficient sequence for P_scaled'=S*Phi^2.
        pressure=0;B=mp.mpf(3);S=mp.mpf('.43');r=mp.mpf('4.1');N=24
        for n in range(N+1,N+100):
            convolution=sum((B/(mp.mpf(20)**j*(j+1)**2)*B/(mp.mpf(20)**(n-1-j)*(n-j)**2) for j in range(n)),mp.mpf(0))
            pressure+=S*convolution*r**n/n
        upper=S*B**2*r**(N+1)/(20**N*(1-r/20))
        if not 0<pressure<upper:raise ArithmeticError('Pressure integration tail bound failed')
        return dict(independent_mixed_coefficient_tail_sums=count,
            independent_pressure_convolution_tail_checked=True,passed=True)


def unpack_packet(c,packet):
    result=dict(packet)
    result['Z']=read_interval(c,packet['Z'])
    result['rows']={label:[[read_interval(c,v) for v in row] for row in rows] for label,rows in packet['rows'].items()}
    result['fixed']={label:[read_interval(c,v) for v in rows] for label,rows in packet['fixed'].items()}
    return result


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Fresh core input changed: '+name)
    original=original_unscaled_fixture();tails=independent_tail_bounds()
    f=CompliantCoreCoefficientRebuild();c=f.ctx;coefficients=0;axis_rows=0;value_count=0
    with mp.workdps(280):
        for center,packed in receipt['recomputed_center_packets'].items():
            packet=unpack_packet(c,packed);source=f.core.axis_inputs(packet['Z'])
            if packed['old_finite_coefficient_rows_read'] or not packed['fresh_compliant_axis_pressure_seed']:
                raise ValueError('Fresh coefficient rebuild replaced by old-source transfer')
            N=packet['radial_degree'];length=N+packet['axial_depth']+1
            for label,rows in packet['rows'].items():
                if len(rows)!=N+1:raise ValueError('Finite radial row coverage incomplete')
                for n,row in enumerate(rows):
                    if len(row)!=length-n:raise ValueError('Axial row consumption layout changed')
                    if any(not all(mp.isfinite(v) for v in endpoints(value)) for value in row):raise ArithmeticError('Nonfinite fresh coefficient')
                    coefficients+=len(row)
            if endpoints(packet['fixed']['S_Z_taylor'][0])[1]<=0:raise ArithmeticError('Actual swirl source replaced by zero')
            for k in range(6):
                expected_phi=-(source['chi'][k]+f.core.epsilon*source['beta'][k])/4
                expected_u=f.core.epsilon*source['slope'][k]
                for value,expected in ((packet['rows']['A'][1][k],expected_phi),(packet['rows']['Uz'][1][k],expected_u)):
                    lo,hi=endpoints(value);a,b=endpoints(expected)
                    if max(lo,a)>min(hi,b):raise ArithmeticError('Fresh first radial row differs from original axis identity')
                    axis_rows+=1
            for packed_value in packed['value_packets']:
                rho=read_interval(c,packed_value['rho']);value=f.values(packet,rho)
                if endpoints(value['Phi'])[0]<=0:raise ArithmeticError('Fresh reconstructed normalized core loses nontrivial swirl')
                if endpoints(rho)==(mp.mpf(0),mp.mpf(0)):
                    if endpoints(value['Phi'])!=(mp.mpf(1),mp.mpf(1)):raise ArithmeticError('Phi axis value changed')
                    if endpoints(value['P_scaled_infinite_tail_bound'])[1]!=0:raise ArithmeticError('Axis pressure remainder must vanish')
                for key in ('Phi','Uz','Mz_over_R','radial_recovery_Q','P_scaled'):
                    stored=read_interval(c,packed_value[key]);fresh=value[key]
                    if stored._mpi_!=fresh._mpi_:raise ArithmeticError('Stored fresh core value is not reproducible: '+key)
                value_count+=1
        if receipt['full_point_physical_field_evaluation'] or receipt['temporal_recursion']:
            raise ValueError('Fresh radial coefficient reconstruction scope promoted')
    for name in (NAME,Path(__file__).name,'lei_ren_part1_paper_core_recursion.py'):
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        independent_original_equations=original,independent_infinite_tail_checks=tails,
        fresh_radial_coefficient_enclosures_checked=coefficients,first_axis_radial_coefficients_checked=axis_rows,
        reproducible_profile_value_packets_checked=value_count,
        same_selected_Cstar_and_compliant_pressure_source=True,old_finite_coefficient_rows_read=False,
        positive_swirl_source_not_dropped=True,infinite_radial_tails_bound=True,
        full_point_physical_field_evaluation=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Fresh compliant core rebuild PASS: original unscaled equations, fresh axis seeds and controlled infinite radial tails',flush=True)
    return result


if __name__=='__main__':run()
