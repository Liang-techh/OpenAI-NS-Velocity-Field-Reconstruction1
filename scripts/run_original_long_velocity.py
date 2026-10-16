"""Assemble and check physical finite-long Cartesian leading+n1 values."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_long_velocity import OriginalLongVelocity
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_long_velocity.json'


def independent_physical_checks():
    """Differentiate an independent physical streamfunction at fixed r,z,t."""
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    f=object.__new__(OriginalLongVelocity);c=mp.mp.clone();c.dps=65
    iv=MPIntervalContext();iv.dps=110;f.ctx=c;f.iv=iv;f.J=J;f.family={'manufactured':True}
    a=SimpleNamespace(ctx=c,iv=iv,precision=65,delta=c.mpf('.07'),delta_iv=iv.mpf('.07'),
        ll=c.ln(3),ll_iv=iv.ln(3),j=c.mpf('.01'),j_iv=iv.mpf('.01'))
    a.summarize=lambda p,t:OriginalCoreLogAmplitude.summarize(a,p,t);f.amplitude=a
    summary=lambda value:f._summary((value,iv.mpf(value)))
    count=0
    def materialize(row):
        return c.mpf(row['coefficient']['approximation'])*c.exp(sum((c.mpf(x['approximation']) for x in row['log_scale_terms']),c.mpf(0)))
    M=lambda R,Z:c.mpf('.2')+Z*Z+R*Z/100
    Psi=lambda R,Z:c.mpf('.2')+Z*Z+R*Z/50
    primitive=lambda Z:c.mpf('.4')+Z**3
    for Z,y in ((c.mpf('.3'),c.mpf('.7')),(c.mpf('-.2'),c.mpf('1.3'))):
        R=110*c.exp(y);delta=a.delta;L=1-delta*Z*Z;d=1-Z*Z
        qpsi=(2*Z*Psi(R,Z)-(1-delta)*Z*M(R,Z)-d*c.diff(lambda zz:M(R,zz),Z))/L
        prim=[summary(c.diff(primitive,Z,k)/c.factorial(k)) for k in range(3)]
        logswirl=c.ln(c.mpf('.7'))+c.ln(220)/2+y/10-c.mpf('.03')*Z*Z
        packet=dict(source_family=f.family,Z=str(Z),chart='long_reshape',
            normalized_modified_leading={n:[summary(v)] for n,v in (('Psi',Psi(R,Z)),('Mpsi',M(R,Z)),('Qpsi',qpsi))},
            coordinate=dict(kind='exact_decimal_R110_log_offset',value=str(y),y=summary(y)),
            original_T_singleton_exact_mpf_tuple=list(c.mpf(7)._mpf_),
            physical_similarity_log_Utheta=[summary(logswirl)],
            positive_order_n1_own_history=dict(actual_terminal_primitive=dict(actual_primitive=prim)),
            selected_integer=dict(manufactured=True),limits={})
        for lt,S in (('-2',None),(None,'.11')):
            output=f._assemble(packet,str(Z),lt,S,'.7')['physical']
            tau=c.exp(c.mpf(lt)) if S is None else c.exp(-2*c.mpf(S)/delta)
            lam=c.sqrt(tau/(1-Z*Z));r=lam*c.sqrt(2*R);zz=Z*lam**(1-delta)
            def psi_physical(rr,zphys):
                def zmap(az):return az*(tau/(1-az*az))**((1-delta)/2)-zphys
                az=c.findroot(zmap,Z);ll=c.sqrt(tau/(1-az*az));ar=rr*rr/(2*ll*ll)
                return ll**(1-delta)*ar*(4*az+a.j+M(ar,az)/3)+ll**(1+delta)*primitive(az)/9
            ur=-c.diff(lambda zphys:psi_physical(r,zphys),zz)/r
            uz=c.diff(lambda rr:psi_physical(rr,zz),r)/r
            ut=c.exp(logswirl)*lam**(-1-delta)
            expected=dict(ur=ur,uz=uz,utheta=ut)
            for name,value in expected.items():
                got=sum((materialize(row) for row in output['cylindrical'][name]),c.mpf(0))
                assert abs(got-value)<c.mpf('1e-50')*max(1,abs(value)),(name,got,value);count+=1
            expectedcart=dict(u=ur*c.cos(c.mpf('.7'))-ut*c.sin(c.mpf('.7')),
                v=ur*c.sin(c.mpf('.7'))+ut*c.cos(c.mpf('.7')),w=uz)
            for name,value in expectedcart.items():
                got=sum((materialize(row) for row in output['Cartesian'][name]),c.mpf(0))
                assert abs(got-value)<c.mpf('1e-50')*max(1,abs(value));count+=1
            got=sum((materialize(row) for row in output['streamfunction']),c.mpf(0))
            assert abs(got-psi_physical(r,zz))<c.mpf('1e-50')*max(1,abs(got));count+=1
            for name,value in (('r',r),('z',zz),('x',r*c.cos(c.mpf('.7'))),('y',r*c.sin(c.mpf('.7')))):
                assert abs(materialize(output['coordinates'][name])-value)<c.mpf('1e-50')*max(1,abs(value));count+=1
    return count


def structure_checks(packet):
    p=packet['physical'];count=0
    assert {k:len(v) for k,v in p['Cartesian'].items()}==dict(u=5,v=5,w=3);count+=1
    assert p['n1_own_history']['y_cancellation_before_enclosure'];count+=1
    for group,name,multiplier in ((p['cylindrical']['ur'],'n1_radial','-.5'),
                                  (p['streamfunction'],'n1_mass',None)):
        row=next(x for x in group if x['name']==name)
        logs=[x for x in row['log_scale_terms'] if x.get('kind')=='original_finite_long_y_log']
        assert [x['multiplier'] for x in logs]==([] if multiplier is None else [multiplier]);count+=1
    swirl=p['cylindrical']['utheta'][0]
    assert [x['kind'] for x in swirl['log_scale_terms'] if 'kind' in x]==['original_long_similarity_swirl_log'];count+=1
    for group in (p['cylindrical']['uz'],p['cylindrical']['utheta']):
        row=group[-1];lo,hi=_bounds(mp.mp,mp.iv.mpf([mp.make_mpf(tuple(v)) for v in row['coefficient']['source_enclosure_exact_mpf_tuples']]))
        assert lo==hi==0;count+=1
    assert p['shared_streamfunction_meridional_recovery'] and not p['physical_gradient_or_vorticity_available'];count+=1
    assert packet['limits']['full_physical_vector_assembled'] and not packet['limits']['global_field'];count+=1
    return count


def normalized_log(field,row,*,coefficient_divisor=0,remove_anchor=False):
    """Compare term magnitudes after analytic coefficient/unit conversion."""
    iv=field.iv;co=field._restore(row['coefficient'])
    if _bounds(field.ctx,abs(co))[1]==0:return None
    lo,hi=_bounds(field.ctx,co)
    if lo<=0<=hi:raise AssertionError('Seam coefficient sign unresolved')
    out=iv.ln(abs(co))-coefficient_divisor
    for log in row['log_scale_terms']:
        if 'negative_of_exact_mpf_tuple' in log:
            if not remove_anchor:out-=iv.mpf(mp.make_mpf(tuple(log['negative_of_exact_mpf_tuple'])))
        else:out+=field._restore(log)
    return out


def run():
    start=time.perf_counter();identities=independent_physical_checks()
    print('Independent fixed-time physical streamfunction, Cartesian and scaled-time checks pass',flush=True)
    f=OriginalLongVelocity();samples=[];structure=0
    for Z in ('-.3','0','.3'):
        for q in ('0','.5','1'):
            packet=f.evaluate(Z,q,log_tau='-2',azimuth='.25');structure+=structure_checks(packet);samples.append(packet)
        print('Actual physical long-chart values: Z='+Z,flush=True)
    samples.append(f.at_log_radius('.3','1',log_tau='-3',azimuth='1.1'))
    samples.append(f.evaluate('.3','.5',log_tau=None,scaled_time='1e18',azimuth='1.7'))
    reference=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_r110_numeric.json')).read_text(encoding='utf-8'))
    seams=histories=0;iv,c=f.iv,f.ctx
    for Z in ('-.3','0','.3'):
        old=next(p['physical'] for p in reference['samples'] if p['Z']==Z and p['source_cover']['local_source_chart']=='post_power' and p['coordinate']['value']=='1')
        new=next(p['physical'] for p in samples if p['Z']==Z and p['coordinate']['value']=='0')
        for group in ('ur','utheta','uz'):
            for a,b in zip(old['cylindrical'][group],new['cylindrical'][group]):
                # Both sides use the same signed physical term. Unit changes
                # of n1 coefficient are offset by the emitted log changes.
                if group=='utheta' and a['name']=='leading_swirl':
                    left=normalized_log(f,a,remove_anchor=True)
                    right=normalized_log(f,b)+iv.mpf(f.amplitude.logC_iv)
                else:left,right=normalized_log(f,a),normalized_log(f,b)
                if left is None or right is None:assert left is right
                else:
                    al,ah=_bounds(c,left);bl,bh=_bounds(c,right);assert max(al,bl)<=min(ah,bh),(Z,group,a['name']);seams+=1
        for a,b in zip(old['streamfunction'],new['streamfunction']):
            left,right=normalized_log(f,a),normalized_log(f,b)
            al,ah=_bounds(c,left);bl,bh=_bounds(c,right);assert max(al,bl)<=min(ah,bh);seams+=1
        for name in ('r','x','y','z'):
            left,right=normalized_log(f,old['coordinates'][name]),normalized_log(f,new['coordinates'][name])
            if left is None or right is None:assert left is right
            else:
                al,ah=_bounds(c,left);bl,bh=_bounds(c,right);assert max(al,bl)<=min(ah,bh);seams+=1
        assert old['n1_own_history']['actual_terminal_primitive']==new['n1_own_history']['actual_terminal_primitive'];histories+=1
    for p in samples:
        assert p['physical']['original_selected_integer']==p['selected_integer'];histories+=1
        assert p['physical']['native']['coordinate']==p['coordinate'];histories+=1
    coarse=OriginalLongVelocity(precision=60,amplitude_precision=90)
    before=coarse.evaluate('.3','.5',log_tau='-2',azimuth='.25')['physical']
    after=next(p['physical'] for p in samples if p['Z']=='.3' and p['coordinate']['value']=='.5')
    comparisons=0
    for group in ('cylindrical','Cartesian'):
        for key,rows in after[group].items():
            for a,b in zip(rows,before[group][key]):
                overlap(a['coefficient'],b['coefficient']);comparisons+=1
                assert len(a['log_scale_terms'])==len(b['log_scale_terms'])
                for aa,bb in zip(a['log_scale_terms'],b['log_scale_terms']):overlap(aa,bb);comparisons+=1
    rejected=0
    for kwargs in (dict(log_tau='1'),dict(log_tau='nan'),dict(azimuth='nan'),dict(scaled_time='1'),dict(log_tau=None,scaled_time='-1')):
        try:f.evaluate('.3','.5',**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Invalid physical time/angle admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,independent_fixed_time_streamfunction_Cartesian_scaled_time_checks=identities,
        actual_original_long_physical_packets=len(samples),physical_geometry_and_analytic_cancellation_checks=structure,
        R110_physical_term_log_seams=seams,original_history_integer_coordinate_checks=histories,
        cross_precision_physical_coefficient_and_log_comparisons=comparisons,rejected_physical_requests=rejected),
        samples=samples,input_hashes=f.hashes,implementation_hashes=code,execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(native_Z_q_or_log_radius_required=True,physical_values_available=True,physical_gradient=False,
            direct_Cartesian_inverse=False,later_charts=False,global_N_admission=False,global_field=False,
            useful_combined_physical_relative_error=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for name,digest in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],execution_seconds=result['execution_seconds']),indent=2),flush=True)
