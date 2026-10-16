"""Physical selected-N Rm/Rh values, with independent streamfunction recovery."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
from openai_ns_reconstruction.original_Rm_patch_velocity import OriginalRmPatchVelocity
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_Rm_patch_velocity.json'


def independent_streamfunction_checks():
    """Materialize moderate test data and differentiate at fixed physical r/z.

    This checks the native-to-physical map and n1 powers without differentiating
    the wrapper's velocity formulas or using unrepresentable original radii.
    """
    c=mp.mp.clone();c.dps=80;iv=MPIntervalContext();iv.dps=100
    f=object.__new__(OriginalRmPatchVelocity)
    a=object.__new__(OriginalCoreLogAmplitude);a.ctx=c;a.iv=iv;a.precision=70
    a.delta=c.mpf('.11');a.delta_iv=iv.mpf('.11')
    a.lp=c.ln(3);a.lp_iv=iv.ln(3);a.ll=c.ln(7);a.ll_iv=iv.ln(7)
    f.ctx=c;f.iv=iv;f.amplitude=a;f.J=J;f.family={'manufactured_test_only':True};f.n1_support={}
    logC=c.mpf('.01');Rm=110*c.exp(10*logC+10*a.lp-6)
    def profiles(ctx,z,R):
        C=1+z+z*z/3
        return 4*z+z*z*R/10**6,4*z+z*z*R/(2*10**6)+C/R,(1+z*z)*R**ctx.mpf('.1')/50
    def primitive(z):return (1+z*z)*c.mpf('.03')-z*c.mpf('.02')
    def point(row):return c.mpf(row['approximation'])
    def value(term):return point(term['coefficient'])*c.exp(sum((point(q) for q in term['log_scale_terms']),c.mpf(0)))
    count=0;worst=c.mpf(0)
    for Z in ('-.3','0','.3'):
        z=c.mpf(Z);R=Rm*c.mpf('1.4')
        zp=J.variable(c,z,2);zi=J.variable(iv,iv.mpf(Z),2)
        U,M,E=profiles(c,zp,R);Ui,Mi,Ei=profiles(iv,zi,iv.mpf(c.nstr(R,80)))
        def Q(ctx,Z,U,M,delta):
            derivative=J(ctx,[(k+1)*M[k+1] for k in range(M.order)])
            return (2*Z*U-Z*M*(1-delta)-(1-Z*Z)*derivative)/(1-Z*Z*delta)
        Qp,Qi=Q(c,zp,U,M,a.delta),Q(iv,zi,Ui,Mi,a.delta_iv)
        jp,ji=primitive(zp),(1+zi*zi)*iv.mpf('.03')-zi*iv.mpf('.02')
        row=lambda p,t:f._summary((p,t))
        packet=dict(source_family=f.family,native_x='1.4',selected_integer={'manufactured':True},limits={},
            physical_similarity={n:[row(p[0],t[0])] for n,p,t in
                [('Uz',U,Ui),('Mz_over_R',M,Mi),('Utheta_over_Pstar',E,Ei),('Ur_over_sqrt_R_over2',Qp,Qi)]},
            generic_modified_fields={'pressure':[row(c.mpf(2),iv.mpf(2))]},
            positive_order_n1_own_history=dict(f_and_u_exact_zero=True,
                actual_terminal_primitive=dict(exact_same_terminal_history_used=True,
                    actual_primitive=[row(jp[k],ji[k]) for k in range(3)])),
            coordinate={'global_y_log_terms':[row(10*logC,iv.mpf('.1')),row(10*a.lp,10*a.lp_iv),
                row(c.mpf(-6),iv.mpf(-6)),row(c.ln(c.mpf('1.4')),iv.ln(iv.mpf('1.4')))]})
        for scaled in (False,True):
            lt=-2*c.mpf('1.2')/a.delta if scaled else c.mpf(-2)
            lam=lambda zz:c.exp(lt/2)/c.sqrt(1-zz*zz)
            rr=lam(z)*c.sqrt(2*R);zz=z*lam(z)**(1-a.delta)
            def psi(rphysical,zphysical,order):
                native=c.findroot(lambda v:v*lam(v)**(1-a.delta)-zphysical,(z-c.mpf('.01'),z+c.mpf('.01')))
                ll=lam(native);r=rphysical*rphysical/(2*ll*ll)
                return ll**(1-a.delta)*r*profiles(c,native,r)[1] if order==0 else primitive(native)*ll**(1+a.delta)/49
            ur=[-c.diff(lambda v:psi(rr,v,k),zz)/rr for k in (0,1)]
            uz=[c.diff(lambda v:psi(v,zz,k),rr)/rr for k in (0,1)]
            ut=3*profiles(c,z,R)[2]*lam(z)**(-1-a.delta)
            for angle in ('0','.7'):
                out=f._assemble(packet,Z,None if scaled else '-2','1.2' if scaled else None,angle)['physical']
                expected=dict(ur=ur,uz=uz,utheta=[ut,c.mpf(0)])
                pairs=[(sum(value(t) for t in out['cylindrical'][n]),sum(v)) for n,v in expected.items()]
                theta=c.mpf(angle);urall=sum(ur);uzall=sum(uz)
                pairs.extend((sum(value(t) for t in out['Cartesian'][n]),v) for n,v in
                    [('u',urall*c.cos(theta)-ut*c.sin(theta)),('v',urall*c.sin(theta)+ut*c.cos(theta)),('w',uzall)])
                pairs.extend((value(out['coordinates'][n]),v) for n,v in
                    [('r',rr),('x',rr*c.cos(theta)),('y',rr*c.sin(theta)),('z',zz)])
                pairs.extend((value(t),psi(rr,zz,k)) for k,t in enumerate(out['streamfunction']))
                for got,exact in pairs:
                    error=abs(got-exact)/max(1,abs(exact));worst=max(worst,error)
                    assert error<c.mpf('1e-60'),(Z,scaled,angle,error);count+=1
    return dict(comparisons=count,max_scaled_error=c.nstr(worst,12),
        independent_fixed_physical_coordinate_streamfunction_differentiation=True)


def inlet_physical_seam(f,old,new):
    """Compare sums in common factored units, including signed zero crossings.

    Identical huge source atoms cancel before interval arithmetic. All remaining
    logarithms are ordinary finite logs, never an exponent of original logC.
    """
    c,iv=f.ctx,f.iv
    def reduced(term,common):
        left=list(term['log_scale_terms']);right=list(common)
        for q in left[:]:
            key=q['source_enclosure_exact_mpf_tuples']
            match=next((r for r in right if r['source_enclosure_exact_mpf_tuples']==key),None)
            if match is not None:left.remove(q);right.remove(match)
        assert all(abs(c.mpf(q['approximation']))<c.mpf('1e30') for q in left+right),'Huge geometry must cancel exactly'
        shift=sum((f._restore(q) for q in left),iv.mpf(0))-sum((f._restore(q) for q in right),iv.mpf(0))
        return f._restore(term['coefficient'])*iv.exp(shift)
    count=0
    for left,right in [(old['cylindrical'][n],new['cylindrical'][n]) for n in ('ur','utheta','uz')]+[(old['streamfunction'],new['streamfunction'])]:
        for order in (0,1):
            aa=[t for t in left if t['temporal_order']==order];bb=[t for t in right if t['temporal_order']==order]
            common=bb[0]['log_scale_terms']
            x=sum((reduced(t,common) for t in aa),iv.mpf(0));y=sum((reduced(t,common) for t in bb),iv.mpf(0))
            lo,hi=_bounds(c,x);ll,hh=_bounds(c,y)
            assert max(lo,ll)<=min(hi,hh),(aa[0]['name'],c.nstr(lo,20),c.nstr(ll,20));count+=1
    return count


def run():
    started=time.perf_counter();independent=independent_streamfunction_checks()
    print('Independent fixed-physical-coordinate streamfunction and Cartesian recovery pass',flush=True)
    f=OriginalRmPatchVelocity();samples=[];geometry=history=seams=scaled_checks=0
    reference_path=DIRECTORY/(PREFIX+'compliant_current_original_reference_velocity.json')
    reference=json.loads(reference_path.read_text(encoding='utf-8'))
    assert reference['checks']['passed'] and reference['source_family']==f.family
    f.hashes[reference_path.name]=hashlib.sha256(reference_path.read_bytes()).hexdigest()
    for name,digest in reference['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        f.code_hashes[name]=digest
    for Z,x in list(f._packets):
        p=f.evaluate(Z,x,log_tau='-2',azimuth='.25');samples.append(p);v=p['physical']
        assert v['original_selected_integer']==p['selected_integer'];history+=1
        assert v['n1_own_history']['actual_terminal_primitive']==p['positive_order_n1_own_history']['actual_terminal_primitive'];history+=1
        assert not v['n1_own_history']['constant_primitive_support_admission']['global_high_order_n1_field_joins_admitted'];history+=1
        mass=v['streamfunction'][1];radial=v['cylindrical']['ur'][1]
        assert not any(q.get('kind')=='original_Rm_geometry_y_log' for q in mass['log_scale_terms']);geometry+=1
        assert all(q['multiplier']=='-.5' for q in radial['log_scale_terms'] if q.get('kind')=='original_Rm_geometry_y_log');geometry+=1
        assert v['cylindrical']['uz'][1]['coefficient']['approximation']=='0.0';geometry+=1
        assert v['cylindrical']['utheta'][1]['coefficient']['approximation']=='0.0';geometry+=1
        if x=='1' and Z in ('-.3','0','.3'):
            old=next(a for a in reference['samples'] if a['Z']==Z and a['chart']=='postrestore' and a['coordinate']['value']=='1')
            seams+=inlet_physical_seam(f,old['physical'],v)
        print('Actual physical patch: Z='+Z+', x='+x[:18],flush=True)
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'Rh',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        for name,terms in p['physical']['cylindrical'].items():
            ordinary=next(a for a in samples if a['Z']==Z and a['native_x']=='Rh' and a['physical']['time']['kind']=='log_tau')
            for a,b in zip(terms,ordinary['physical']['cylindrical'][name]):
                overlap(a['coefficient'],b['coefficient']);scaled_checks+=1
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled_checks+=1
    rejected=0
    for kwargs in (dict(log_tau='1'),dict(azimuth='nan'),dict(log_tau='-2',scaled_time='1')):
        try:f.evaluate('.3','1',**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid physical request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,independent_streamfunction=independent,
        actual_original_finite_N_physical_packets=len(samples),reference_to_Rm_physical_sum_seams=seams,
        original_N_history_and_local_support_checks=history,geometry_and_zero_local_n1_checks=geometry,
        scaled_time_coefficient_and_cancellation_checks=scaled_checks,rejected_requests=rejected),
        samples=samples,input_hashes=f.hashes,implementation_hashes=code,execution_seconds=round(time.perf_counter()-started,3),
        limits=dict(native_coordinates_required=True,through_original_Rh=True,physical_Cartesian_values=True,
            finite_N_values_are_point_plus_directed_error=True,resolved_signed_finite_N_centers=False,
            leading_pressure_only=True,full_n1_pressure=False,physical_gradients=False,direct_physical_inverse=False,
            global_n1_high_order_joins=False,global_N_admission=False,global_field=False,
            genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    result=json.loads(path.read_text(encoding='utf-8'));assert result['checks']['passed']
    for name,digest in result['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in result['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],execution_seconds=result['execution_seconds']),indent=2),flush=True)
