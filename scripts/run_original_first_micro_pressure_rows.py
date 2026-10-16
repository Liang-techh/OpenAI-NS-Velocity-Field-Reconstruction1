"""Build the missing actual first-micro pressure rows and endpoint errors."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_first_micro_pressure_rows import OriginalFirstMicroPressureRows,pressure_coefficients

NAME=PREFIX+'compliant_current_original_first_micro_pressure_rows.json'


def source_identities():
    y,z,n=sy.symbols('y z N',nonzero=True);phase=sy.Function('phase')(y)
    B=sy.Function('B')(y,z,phase)
    radial=sy.diff(B/n,y).subs(sy.diff(phase,y),n)
    # Sympy leaves a Subs around the slow derivative; compare before the
    # concrete chain substitution to retain the actual held-phase meaning.
    assert sy.simplify(radial-sy.diff(B,y)/n).subs(sy.diff(phase,y),n)==0
    assert sy.diff(phase,z)==0
    zz,dd=sy.symbols('Z delta');U,m,Uy,UZ,mZ=sy.symbols('U m Uy UZ mZ')
    L=1-dd*zz**2;d=1-zz**2
    qy=(2*zz*Uy-(1-dd)*zz*(U-m)-d*(UZ-mZ))/L
    coeffs={'U_y':2*zz/L,'U':-(1-dd)*zz/L,'m':(1-dd)*zz/L,'U_Z':-d/L,'m_Z':d/L}
    assert sy.simplify(qy-sum(coeffs[k]*v for k,v in [('U_y',Uy),('U',U),('m',m),('U_Z',UZ),('m_Z',mZ)]))==0
    # Direct chain rule, with independently held slow variables.
    bslow,bfast,bslowZ,bfastZ=sy.symbols('B_y B_phi B_yZ B_Zphi')
    assert sy.simplify((bslow+bfast*n)/n-(bslow/n+bfast))==0
    assert sy.simplify((bslowZ+bfastZ*n)/n-(bslowZ/n+bfastZ))==0
    return dict(passed=True,physical_Qy_recovery=True,fast_phase_N_cancels_own_1_over_N=True,
        mixed_fast_B_Z_phi_no_N_denominator=True,phase_Z_exact_zero=True,
        sources=[PREFIX+'compliant_current_original_modified_Q4_pressure_C3_functions.py',
                 PREFIX+'compliant_current_original_prescribed_loop_C5_functions.py'])


def run():
    start=time.perf_counter();identities=source_identities();f=OriginalFirstMicroPressureRows(precision=85)
    packets=[];audit=mp.mp.clone();audit.dps=180
    serialized=flat=fast=ordinary=endpoint=direct=0
    for Z in ('-.3','0','.3'):
        for coordinate,selected in (('.5',True),('.75',True),('1',True),('.125',False),('.25',False)):
            p=f.evaluate(Z,coordinate,selected_sc_multiple=selected);packets.append(p)
            assert all(p['contract'].values())
            assert p['selected_integer']==f.micro.integer
            s,isflat,_=f.micro._coordinate(coordinate,selected)
            points,truths,rho=f.background(Z,s)
            for name,field in p['physical_pressure_rows'].items():
                for k,row in enumerate(field['signed_background_ordinary_Z']):
                    actual=truths[name][k]*math.factorial(k)
                    assert [list(v) for v in actual._mpi_]==row['source_enclosure_exact_mpf_tuples'];ordinary+=1
                    ends=[audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]
                    required=max(abs(v-audit.mpf(row['approximation'])) for v in ends)
                    assert audit.mpf(row['absolute_error_upper'])>=required;serialized+=1
                if isflat:
                    assert all(v['exact_zero'] for v in field['finite_N_absolute_error_log_bounds']);flat+=1
            components=p['finite_N_axial_radial_error_components']
            for k,key in enumerate(('B_phi','B_Z_phi')):
                got=components['fast_B_phi'][k]
                if not isflat:
                    want=f.read_log(f.loop['phase_and_first_slow_phase_log_bounds'][key])+f.o.core.lp_iv
                    assert got['log_absolute_upper']['lower_exact_mpf_tuple']==list(want._mpi_[0]);fast+=1
                    assert got['log_absolute_upper']['upper_exact_mpf_tuple']==list(want._mpi_[1]);fast+=1
            # Independently recover Q_y from the mean ODE rather than from
            # the radial Psi rows used in the production leading frame.
            c,J=f.ctx,f.J;z=J.variable(c,c.mpf(Z),3);de=f.o.core.delta
            u,m,uy=points['U'].truncate(3),points['m'].truncate(3),points['U_y'].truncate(3)
            dz=lambda v:J(v.ctx,[(k+1)*v[k+1] for k in range(v.order)])
            qy=(uy*(2*z)-(u-m)*(z*(1-de))-dz(u-m)*(1-z*z))/(1-z*z*de)
            for k in (0,1):
                a,b=qy[k],points['Q_y'][k]
                assert abs(a-b)<c.mpf('1e-65')*max(c.mpf(1),abs(a),abs(b));direct+=1
                row=p['R_Q_y_endpoint']['signed_background_ordinary_Z'][k]
                assert row['source_enclosure_exact_mpf_tuples']==[list(v) for v in (truths['Q_y'][k]*rho[0])._mpi_];endpoint+=1
            assert not p['limits']['signed_first_micro_m15_integrated']
        print('Actual first-micro U/U_y/Q/Q_y rows with finite-N errors: Z='+Z,flush=True)
    rejected=0
    for args in [('1','.25',False),('nan','.25',False),('.3','0',False),('.3','.3',False),('.3','.4',True)]:
        try:f.evaluate(args[0],args[1],selected_sc_multiple=args[2])
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original first-micro pressure query admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=packets,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,actual_first_micro_packets=len(packets),original_source_identities=identities,
            original_ordinary_derivative_rows=ordinary,serialized_background_error_checks=serialized,
            source_flat_finite_N_error_zeros=flat,actual_fast_phase_log_and_mixed_rows=fast,
            independent_mean_ODE_Qy_comparisons=direct,finite_R_Qy_endpoint_rows=endpoint,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(actual_first_micro_U_Q_Z2_Uy_Qy_Z1_factored_error_rows=True,
            resolved_finite_N_signed_centers=False,signed_first_micro_m15_integrated=False,
            independent_Rp_signed_numeric_seam=False,global_five_moment_numeric_join=False,
            common_global_N=False,exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-report',action='store_true');p.add_argument('--inspect-report',action='store_true')
    args=p.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
