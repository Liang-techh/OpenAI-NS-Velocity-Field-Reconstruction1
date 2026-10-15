"""Evaluate the original selected-N first switch with unequal source controls."""
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_first_switch_numeric import OriginalFirstSwitchNumeric
from openai_ns_reconstruction.original_core_numeric import _bounds
from run_original_second_micro_numeric import overlap

NAME='lei_ren_part1_paper_compliant_current_original_first_switch_numeric.json'


def independent_source_checks():
    c=mp.mp.clone();c.dps=65;iv=MPIntervalContext();iv.dps=110
    module=importlib.import_module('lei_ren_part1_paper_compliant_actual_bridge_integrals');J=module.IntervalTaylor
    b=object.__new__(module.CompliantActualBridgeIntegrals)
    h=iv.mpf(1)/101;hp=c.mpf(1)/101;eps=c.mpf('.5')
    b.ctx=iv;b.cap=h;b.weight=iv.mpf('.1');b.logh=iv.ln(h)
    b.pressure_log=iv.ln(11);b.swirl_log=iv.ln(13);b.pulse_cache={}
    z=c.mpf('.3')
    functions={
        'D_over_R':[lambda z:2+z,lambda z:-1+z*z,lambda z:1-z+z**3],
        'drive_hydro':[lambda z:1+z,lambda z:1-z*z,lambda z:z**3],
        'drive_pressure':[lambda z:2-z,lambda z:-z,lambda z:1+z*z],
        'drive_swirl':[lambda z:-1-z*z,lambda z:z,lambda z:1-z**3]}
    factor=lambda name:c.mpf('1e-4')
    def jet(fn):return J(iv,[iv.mpf(c.nstr(c.diff(fn,z,k)/math.factorial(k),65)) for k in range(4)])
    modes={name:[jet(lambda z,fn=fn,name=name:fn(z)*factor(name)) for fn in rows] for name,rows in functions.items()}
    p=dict(data=dict(phi0=J.constant(iv,1,3)),R0=iv.mpf(2),modes=modes,
        quotient=jet(lambda z:1+z/5),angular_prefix_norm_per_hb=iv.mpf(0))
    leading=SimpleNamespace(iv=iv,bridge=b,_prepared={'.3':p},core=SimpleNamespace(epsilon_iv=iv.mpf('.5')))
    leading._transfer=lambda x:iv.mpf([mp.make_mpf(v) for v in x._mpi_])
    leading._abs=lambda x:abs(x)
    leading._width_error=lambda x:h*abs(x)
    leading._scaled=lambda row,term,divide_eps=False:row*iv.exp(term['additional_positive_scale_log_enclosure']+(term['width_power']-1)*b.logh)/(iv.mpf('.5') if divide_eps else 1)
    fake=object.__new__(OriginalFirstSwitchNumeric);fake.leading=leading
    fake.transport=SimpleNamespace(ctx=c,iv=iv);fake._switch_increment_cache={};fake.width_identity={'finite_parameter_fixture_only':True};fake.Y=iv.ln(50)
    def sigma(x):
        if x<=0:return c.mpf(0)
        if x>=1:return c.mpf(1)
        odds=1/(1-x)**2-1/x**2
        if odds>1000:return c.mpf(1)
        if odds<-1000:return c.mpf(0)
        return 1/(1+c.exp(-odds))
    def primitive_D(z,t):
        total=c.mpf(0)
        for j,fn in enumerate(functions['D_over_R']):
            a=1-j
            integral=t if a==0 else c.expm1(a*hp*t)/(a*hp)
            total+=fn(z)*factor('D_over_R')*c.mpf(2)**j*100**(1-j)*integral
        return total
    def ell(z,t):return -hp*hp*primitive_D(z,t)/2
    def drive(z,t):
        total=c.mpf(0)
        for name,power,scale in (('hydro',1,1),('pressure',1,11),('swirl',2,13)):
            for j,fn in enumerate(functions['drive_'+name]):
                total+=scale*fn(z)*factor(name)*c.mpf(2)**j*100**(power-j)*c.exp((power-j)*hp*t)
        return total
    checks=0
    for phase in ('.4','1'):
        s=c.mpf(phase);Eell,Epsi,evidence=fake._switch_increments('.3',iv.mpf(phase))
        # Independent scalar source quadrature and genuine Z derivatives
        # test production positive-modal bounds at representable finite h.
        delta_phi=lambda z:c.expm1(ell(z,s))
        delta_psi=lambda z:-hp*hp/eps*c.quad(lambda t:(1-sigma(t))*(1+z/5)*c.exp(ell(z,t))*drive(z,t),[0,s/2,s])
        for k in range(4):
            phi_value=abs(c.diff(delta_phi,z,k)/math.factorial(k))
            psi_value=abs(c.diff(delta_psi,z,k)/math.factorial(k))
            phi_bound=_bounds(c,Eell*iv.exp(Eell)/b.weight**k)[1]
            psi_bound=_bounds(c,Epsi/b.weight**k)[1]
            assert phi_value<=phi_bound and psi_value<=psi_bound;checks+=2
        mass=c.quad(lambda t:1-sigma(t),[0,s/2,s])
        box=iv.mpf([mp.make_mpf(tuple(v)) for v in evidence['axial_prefix_mass_exact_mpf_tuples']])
        lo,hi=_bounds(c,box);assert lo-c.mpf('1e-55')<=mass<=hi+c.mpf('1e-55');checks+=1
        assert abs(mass-s)>c.mpf('1e-5');checks+=1
        rho,base,tiny=fake._switch_geometry(iv.mpf(phase))
        lo,hi=_bounds(c,rho);actual=100/eps*c.exp(hp*s)
        assert lo-c.mpf('1e-55')<=actual<=hi+c.mpf('1e-55');checks+=1
    # A nonzero incoming mean is not reset, even when velocity is constant.
    for s in (c.mpf('.4'),c.mpf(1)):
        V=c.mpf(3);M0=c.mpf(2)
        mean=c.exp(-hp*s)*M0+hp*c.quad(lambda t:c.exp(-hp*(s-t))*V,[0,s])
        assert 0<abs(mean-M0)<=hp*s*abs(V-M0);checks+=1
    return checks


def run():
    began=time.perf_counter();independent=independent_source_checks()
    print('Independent unequal controls, full source, Z3 and mean-history checks pass',flush=True)
    f=OriginalFirstSwitchNumeric();samples=[];proofs=0
    for Z in ('-.3','0','.3'):
        for s in ('0','.5','1'):
            p=f.evaluate(Z,s,log_tau='-2',azimuth='.25');samples.append(p)
            assert [len(p['normalized_modified_leading'][n]) for n in ('Phi','Psi','Mpsi','Qpsi')]==[4,4,3,2]
            control=p['background_control_contract']
            assert control['angular_weight']=='1' and control['axial_weight']=='1-sigma(s)'
            assert control['original_hb_power_angular']==control['original_hb_power_axial']==2
            assert p['source_cover']['local_source_chart']=='first_switch'
            assert p['physical']['n1_own_history']['current_rho_not_rho4_used']
            assert not p['physical']['physical_gradient_or_vorticity_available']
            assert not p['limits']['global_field'] and not p['limits']['actual_mean_Z3_numeric']
            for row in p['small_error_proofs']:assert row['numerical_absolute_error_cap']=='1e-1000';proofs+=1
        print('Actual original-N first-switch values and own history: Z='+Z,flush=True)
    samples.append(f.evaluate('.17','.37',log_tau='-3',azimuth='1.1'))
    reference=json.loads((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_macro_numeric.json').read_text(encoding='utf-8'))
    seams=histories=geometry=0
    for Z in ('-.3','0','.3'):
        previous=next(p for p in reference['samples'] if p['Z']==Z and p['coordinate']['value']=='1')
        inlet=next(p for p in samples if p['Z']==Z and p['coordinate']['value']=='0')
        for name,rows in previous['normalized_modified_leading'].items():
            for a,b in zip(rows,inlet['normalized_modified_leading'][name]):overlap(a,b);seams+=1
        for name,rows in previous['physical']['cylindrical'].items():
            for a,b in zip(rows,inlet['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);seams+=1
        assert previous['original_five_correction_history_bounds']==inlet['original_five_correction_history_bounds'];histories+=20
        assert inlet['physical']['native']['R100_inlet_exact'];geometry+=1
        for p in [p for p in samples if p['Z']==Z]:
            own=p['physical']['n1_own_history'];assert own['actual_terminal_primitive']==previous['physical']['n1_own_history']['actual_terminal_primitive'];histories+=1
            rho=f._numeric(Z,p['coordinate']['value'],False)[3]
            mean=f.transport.iv.mpf([mp.make_mpf(tuple(v)) for v in own['mean_u']['source_enclosure_exact_mpf_tuples']])
            row=own['actual_terminal_primitive']['actual_primitive'][0]
            primitive=f.transport.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])
            lo,hi=_bounds(f.transport.ctx,mean*rho-primitive);assert lo<=0<=hi;geometry+=1
        exit=next(p for p in samples if p['Z']==Z and p['coordinate']['value']=='1')
        assert exit['physical']['native']['R1_endpoint_exact'];geometry+=1
    coarse=OriginalFirstSwitchNumeric(precision=60,amplitude_precision=90)
    old=coarse.evaluate('.3','1',log_tau='-2',azimuth='.25')
    new=next(p for p in samples if p['Z']=='.3' and p['coordinate']['value']=='1')
    comparisons=0
    for name,rows in new['normalized_modified_leading'].items():
        for a,b in zip(rows,old['normalized_modified_leading'][name]):overlap(a,b);comparisons+=1
    for name,rows in new['physical']['cylindrical'].items():
        for a,b in zip(rows,old['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);comparisons+=1
    for a,b in zip(new['physical']['coordinates']['r']['log_scale_terms'],old['physical']['coordinates']['r']['log_scale_terms']):
        if 'approximation' in a and 'approximation' in b:overlap(a,b);comparisons+=1
    rejected=0
    for Z,s,kwargs in (('1','.5',{}),('.3','-.01',{}),('.3','1.01',{}),('.3','nan',{}),('.3','1',{'log_tau':'1'}),('.3','1',{'azimuth':'nan'})):
        try:f.evaluate(Z,s,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Unsupported first-switch request admitted')
    try:f.selected('.3','.5')
    except ValueError:rejected+=1
    else:raise AssertionError('Selected-sc coordinate admitted as first switch')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,
        independent_unequal_control_full_source_Z3_geometry_and_mean_checks=independent,
        actual_original_first_switch_packets=len(samples),selected_N_small_error_proofs=proofs,
        R100_native_physical_coefficient_seams=seams,all_incoming_correction_and_terminal_primitive_memory_checks=histories,
        actual_rho_and_endpoint_geometry_checks=geometry,cross_precision_coefficients_and_log_geometry_comparisons=comparisons,
        rejected_unsupported_requests=rejected),samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        limitations=dict(Phi_Psi_Z_order=3,own_mean_Z_order=2,Qpsi_Z_order=1,
            total_radial_derivatives=False,resolved_fast_phase=False,resolved_signed_correction_integrals=False,
            mean_Z3_and_QZZ_numeric=False,second_switch_and_R110=False,later_charts=False,global_N_admission=False,
            global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False),execution_seconds=round(time.perf_counter()-began,3))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for name,digest in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    r=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=r['checks'],execution_seconds=r['execution_seconds']),indent=2),flush=True)
