"""Actual reference, axial restoration and postrestore values to unpatched Rm."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_reference_numeric import OriginalReferenceNumeric
from openai_ns_reconstruction.original_reference_velocity import OriginalReferenceVelocity
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_reference_velocity.json'


def independent_raw_ODE_checks():
    """Integrate raw H/M/K/A/B/C equations, independently of sector formulas."""
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    import lei_ren_part1_paper_compliant_reference_restore_profiles as source
    from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
    c=mp.mp.clone();c.dps=55;iv=MPIntervalContext();iv.dps=100
    f=object.__new__(OriginalReferenceNumeric);f.ctx=c;f.iv=iv;f.J=J;f.source=source;f.sigma_jets=sigma_jets
    f._kernels={};f.lengthp=c.mpf('1.2');f.length=iv.mpf('1.2')
    f.leading=SimpleNamespace(_abs=lambda v:abs(v),_symmetric=lambda v:iv.mpf([-_bounds(c,v)[1],_bounds(c,v)[1]]))
    f.owner=SimpleNamespace(_cap=lambda *args:(_ for _ in ()).throw(AssertionError('No synthetic large-gap cap expected')))
    z=J.variable(c,c.mpf('.3'),2);zi=J.variable(iv,iv.mpf('.3'),2)
    def data(ctx,Z):
        return dict(H_over_Phi=1+Z+Z*Z,Kpsi_over_Phi=2-Z+Z*Z,
            B_over_Phi2=3+Z,C_over_Phi2=4-Z*Z,Mpsi=1-Z/2,Apsi=2+Z*Z)
    pp,tt=data(c,z),data(iv,zi);vp=2+z*z;vt=2+zi*zi
    inp=dict(six={n:(pp[n],tt[n]) for n in pp},Psi=(J(c,list(vp.coefficients)+[c.mpf(0)]),J(iv,list(vt.coefficients)+[iv.mpf(0)])))
    f.inputs=lambda Z:inp
    phi0=1+z*z/5;j=c.mpf('.03');eps=c.mpf('.2');u0=4*z+j
    initial=dict(H=pp['H_over_Phi']*phi0,M=u0+pp['Mpsi']*eps,
        K=u0*pp['H_over_Phi']*phi0+pp['Kpsi_over_Phi']*phi0*eps,
        A=u0*u0+2*u0*pp['Mpsi']*eps+pp['Apsi']*eps*eps,
        B=pp['B_over_Phi2']*phi0*phi0,C=pp['C_over_Phi2']*phi0*phi0)
    def alpha(s):
        if s<=f.lengthp:return c.mpf(1)
        t=s-f.lengthp
        if t>=1:return c.mpf(0)
        odds=1/(1-t)**2-1/t**2
        return c.mpf(0) if odds>10000 else c.mpf(1) if odds<-10000 else 1/(1+c.exp(odds))
    checks=0
    cases=[('reference','.6'),('restoration','.37'),('restoration','1'),('postrestore','.4')]
    for chart,t in cases:
        tp=c.mpf(t)
        if chart=='reference':six,b,v=f._reference_six('.3',tp,iv.mpf(t),[]);total=tp
        elif chart=='restoration':six,b,v=f._restoration_six('.3',t,[]);total=f.lengthp+tp
        else:six,b,v=f._post_six('.3',tp,iv.mpf(t),[]);total=f.lengthp+1+tp
        def raw(ctx,Z,index):
            take=lambda n:six[n][index];base=lambda n:b[n][index]
            jj,ee=ctx.mpf('.03'),ctx.mpf('.2')
            F=(1+Z*Z/5)*ctx.exp(-ctx.mpf('.4')*ctx.mpf(str(total)))
            H=take('H_over_Phi')*F
            return dict(H=H,M=4*Z+base('gamma1')*jj+take('Mpsi')*ee,
                K=4*Z*H+base('Kalpha_over_Phi')*F*jj+take('Kpsi_over_Phi')*F*ee,
                A=16*Z*Z+8*Z*jj*base('gamma1')+base('gamma2')*(jj*jj)+8*Z*take('Mpsi')*ee+
                    base('Dpsi')*(2*jj*ee)+take('Apsi')*(ee*ee),
                B=take('B_over_Phi2')*F*F,C=take('C_over_Phi2')*F*F)
        computed=raw(iv,zi,1)
        breaks=[c.mpf(0)]+[s for s in (f.lengthp,f.lengthp+1) if s<total]+[total]
        for name,rate in dict(H=2,M=1,K=2,A=1,B=2,C=1).items():
            for k in range(3):
                def integrand(s):
                    F=phi0*c.exp(-c.mpf('.4')*s);a=alpha(s);V=4*z+j*a+vp*(eps*a)
                    target=dict(H=F*2,M=V,K=F*V*2,A=V*V,B=F*F,C=F*F)[name]
                    return c.exp(-rate*(total-s))*target[k]
                actual=initial[name][k]*c.exp(-rate*total)+c.quad(integrand,breaks)
                lo,hi=_bounds(c,computed[name][k]);tol=c.mpf('1e-42')*max(1,abs(actual))
                assert lo-tol<=actual<=hi+tol,(chart,t,name,k,actual,lo,hi);checks+=1
    return checks


def physical_seam(field,left,right):
    count=0;c,iv=field.ctx,field.iv
    def terms(a,b):
        # All sides now use the same primitive/leading coefficient units.
        # Signed interval overlap remains meaningful at a zero crossing;
        # log(abs(coefficient)) would incorrectly require a resolved sign.
        overlap(a['coefficient'],b['coefficient'])
        x=sum((field._restore(log) for log in a['log_scale_terms']),iv.mpf(0))
        y=sum((field._restore(log) for log in b['log_scale_terms']),iv.mpf(0))
        lo,hi=_bounds(c,x);ll,hh=_bounds(c,y)
        assert max(lo,ll)<=min(hi,hh),a['name']
        return 2
    for name in ('ur','utheta','uz'):
        for a,b in zip(left['cylindrical'][name],right['cylindrical'][name]):
            count+=terms(a,b)
    for a,b in zip(left['streamfunction'],right['streamfunction']):
        count+=terms(a,b)
    assert left['n1_own_history']['actual_terminal_primitive']==right['n1_own_history']['actual_terminal_primitive'];count+=1
    return count


def run():
    started=time.perf_counter();independent=independent_raw_ODE_checks()
    print('Independent raw six-history ODEs and separated j/epsilon transport pass',flush=True)
    f=OriginalReferenceVelocity();samples=[];seams=fields=history=geometry=0
    pick=lambda Z,chart,t:next(p for p in samples if p['Z']==Z and p['chart']==chart and p['coordinate']['value']==t and p['coordinate']['kind']!='exact_reference_gap')
    reference=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_long_velocity.json')).read_text(encoding='utf-8'))
    for Z in ('-.3','0','.3'):
        for method in ('reference','restoration','postrestore'):
            for t in ('0','.5','1'):samples.append(getattr(f,method)(Z,t,log_tau='-2',azimuth='.25'))
        print('Actual reference/restoration/postrestore: Z='+Z,flush=True)
        old=next(p for p in reference['samples'] if p['Z']==Z and p['coordinate']['kind']=='exact_decimal_q' and p['coordinate']['value']=='1')
        pairs=((old,pick(Z,'reference','0')),(pick(Z,'reference','1'),pick(Z,'restoration','0')),
               (pick(Z,'restoration','1'),pick(Z,'postrestore','0')))
        for a,b in pairs:
            for n,rows in a['normalized_modified_leading'].items():
                for aa,bb in zip(rows,b['normalized_modified_leading'][n]):overlap(aa,bb);fields+=1
            for n,rows in a['unmodified_six_current_amplitude_units'].items():
                for aa,bb in zip(rows,b['unmodified_six_current_amplitude_units'][n]):overlap(aa,bb);fields+=1
            assert a['original_five_correction_history_bounds']==b['original_five_correction_history_bounds'];history+=1
            seams+=physical_seam(f,a['physical'],b['physical'])
        for p in [p for p in samples if p['Z']==Z]:
            assert p['baseline_contract']['j_and_epsilon_sectors_never_merged'];history+=1
            assert p['physical']['original_selected_integer']==old['physical']['original_selected_integer'];history+=1
            mass=next(r for r in p['physical']['streamfunction'] if r['name']=='n1_mass')
            assert not any(l.get('kind')=='original_reference_geometry_log' for l in mass['log_scale_terms']);geometry+=1
        for chart,t,offset in (('restoration','0',-8),('restoration','1',-7),('postrestore','1',-6)):
            p=pick(Z,chart,t);row=p['coordinate']['y_log_terms'][-1]
            lo,hi=_bounds(f.ctx,f._restore(row));assert lo<=offset<=hi;geometry+=1
            # Pstar cancels before comparing the finite swirl shape.
            value=f._restore(p['physical_similarity_log_Utheta'][0])-f.amplitude.lp_iv
            exact=f.iv.mpf(offset)/10-f.iv.ln(1+f.iv.mpf(Z)**2)
            lo,hi=_bounds(f.ctx,value);ll,hh=_bounds(f.ctx,exact);assert max(lo,ll)<=min(hi,hh);geometry+=1
        end=pick(Z,'postrestore','1');assert f.ctx.mpf(end['physical']['baseline_coefficients']['U']['approximation'])==4*f.ctx.mpf(Z)
        assert end['normalized_modified_leading']['Psi'][0]['approximation']=='0.0';geometry+=2
    samples.append(f.at_reference_gap('.3','1',log_tau='-3',azimuth='1.1'))
    samples.append(f.at_reference_gap('.3','10000',log_tau='-3',azimuth='1.1'))
    samples.append(f.postrestore('.3','1',log_tau=None,scaled_time='1e18',azimuth='1.7'))
    coarse=OriginalReferenceVelocity(precision=60,amplitude_precision=90)
    old=coarse.restoration('.3','.5',log_tau='-2',azimuth='.25');new=pick('.3','restoration','.5');comparisons=0
    for group in ('normalized_modified_leading','unmodified_six_current_amplitude_units'):
        for name,rows in new[group].items():
            for a,b in zip(rows,old[group][name]):overlap(a,b);comparisons+=1
    for name in ('alpha','gamma1','gamma2'):overlap(new['baseline_axial_histories'][name],old['baseline_axial_histories'][name]);comparisons+=1
    for name,rows in new['physical']['cylindrical'].items():
        for a,b in zip(rows,old['physical']['cylindrical'][name]):
            overlap(a['coefficient'],b['coefficient']);comparisons+=1
            for aa,bb in zip(a['log_scale_terms'],b['log_scale_terms']):overlap(aa,bb);comparisons+=1
    rejected=0
    for method,t,kwargs in (('reference','-.1',{}),('restoration','1.1',{}),('postrestore','nan',{}),
                             ('restoration','.5',dict(log_tau='1')),('postrestore','.5',dict(azimuth='nan'))):
        try:getattr(f,method)('.3',t,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Invalid reference request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,independent_raw_six_ODE_coefficients=independent,
        original_reference_restoration_postrestore_packets=len(samples),all_three_field_six_history_seam_coefficients=fields,
        all_three_physical_term_seams=seams,selected_N_and_nonzero_history_checks=history,
        exact_small_offset_swirl_target_and_primitive_cancellation_checks=geometry,cross_precision_comparisons=comparisons,rejected_requests=rejected),
        samples=samples,input_hashes=f.hashes,implementation_hashes=code,execution_seconds=round(time.perf_counter()-started,3),
        limits=dict(native_coordinates_required=True,through_unpatched_Rm=True,direct_physical_inverse=False,
            radial_or_fast_phase_derivatives=False,merged_modified_six=False,pressure_exported=False,actual_patch_connected=False,
            global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for n,d in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==d,n
    for n,d in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==d,n
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],execution_seconds=result['execution_seconds']),indent=2),flush=True)
