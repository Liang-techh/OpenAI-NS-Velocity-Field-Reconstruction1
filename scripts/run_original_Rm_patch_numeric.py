"""Original leading implicit controls and independent five raw patch integrals."""
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
from openai_ns_reconstruction.original_Rm_patch_numeric import OriginalRmPatchNumeric,_Bump,implicit_Z2
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_Rm_patch_numeric.json'


def fixture():
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    from lei_ren_part1_paper_interval_five_bump_inverse import weights
    c=mp.mp.clone();c.dps=65;iv=MPIntervalContext();iv.dps=110
    f=object.__new__(OriginalRmPatchNumeric);f.ctx=c;f.iv=iv;f.J=J;f.precision=65
    raw=json.loads((DIRECTORY/(PREFIX+'bump_integral_enclosures_check.json')).read_text(encoding='utf-8'))
    constants=json.loads((DIRECTORY/(PREFIX+'shared_bump_constants.json')).read_text(encoding='utf-8'))
    f.W=weights(iv,raw);f.bump=_Bump(iv,raw,constants)
    f.normp=c.quad(lambda s:c.exp(-1/(1-s*s)) if abs(s)<1 else c.mpf(0),[-1,0,1])
    f._weights={};f._point_weights={};f._packets={}
    def restore(row):return iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])
    def rows(p,t):
        return [dict(approximation=c.nstr(a,65),source_enclosure_exact_mpf_tuples=[list(v) for v in b._mpi_]) for a,b in zip(p.coefficients,t.coefficients)]
    f.reference=SimpleNamespace(previous=SimpleNamespace(_restore=restore),_rows=rows,
        leading=SimpleNamespace(core=SimpleNamespace(delta=c.mpf('.0001'),delta_iv=iv.mpf('.0001'))))
    f.family={};f.defect=SimpleNamespace(logP=iv.ln(3),logPp=c.ln(3))
    return f


def independent_known_implicit_family():
    from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
    f=fixture();c,iv,J=f.ctx,f.iv,f.J;z=J.variable(iv,iv.mpf('.17'),2)
    W={n:v if n in ('L','gg') else [q*iv.mpf('1e-12') for q in v] for n,v in f.W.items()}
    alpha=(1+z*z)*(1+z*z)*iv.mpf('1e-12')
    exact=[(1+z*(i+1)+z*z*(i-2))*iv.mpf(str((i+1)*1e-5)) for i in range(5)]
    Q=nonlinear_series(iv,W,exact,alpha)
    d=[-Q[i]-sum((exact[j]*W['L'][i][j] for j in range(5)),exact[0]*0) for i in range(5)]
    with mp.workdps(130):control,proof=implicit_Z2(iv,J,W,d,alpha)
    count=0
    for a,b in zip(exact,control):
        for aa,bb in zip(a.coefficients,b.coefficients):
            lo,hi=_bounds(c,aa);ll,hh=_bounds(c,bb);assert max(lo,ll)<=min(hi,hh);count+=1
    assert _bounds(c,proof['first']['contraction_bound'])[1]<1;count+=1
    assert _bounds(c,proof['second_row_contraction'])[1]<1;count+=1
    return count


def independent_raw_patch_integrals():
    """Integrate actual five densities in x with their original own rates."""
    from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
    f=fixture();c,iv,J=f.ctx,f.iv,f.J;Z='.23';p,t=[],[]
    Wp=f._W(c.mpf(2),False,True)
    for directed,ctx,W in ((False,c,Wp),(True,iv,f.W)):
        z=J.variable(ctx,ctx.mpf(Z),2);am=(1+z*z).reciprocal()*ctx.exp(ctx.mpf('-.6'))
        h=[(1+z*(i+1)+z*z*(i-2))*ctx.mpf(str((i+1)*.00001)) for i in range(5)]
        alpha=(am*am).reciprocal()/9;Q=nonlinear_series(ctx,W,h,alpha)
        d=[-Q[i]-sum((h[j]*W['L'][i][j] for j in range(5)),h[0]*0) for i in range(5)]
        (t if directed else p).extend((h,d,-1-z*z))
    rows=f.reference._rows
    inlet=dict(independent_P0_over_Pstar_squared=rows(p[2],t[2]),selected_integer='manufactured fixture',
        finite_N_signed_unknown_history_errors={},positive_order_n1_own_history={})
    data=dict(point=p[0],truth=t[0],normalized_point=p[0],normalized_truth=t[0],
        defects_point=p[1],defects_truth=t[1],inlet=inlet,proof={})
    f.coefficients=lambda Z:data
    z=J.variable(c,c.mpf(Z),2);am=(1+z*z).reciprocal()*c.exp(c.mpf('-.6'));h,d=p[:2]
    M0=z*4+d[0]
    initial=dict(m=M0/3,h=am*(d[2]+c.mpf(5)/8),
        k=am*(z*(d[2]+c.mpf(5)/8)*4+d[1])/3,
        e=(z*M0*8-z*z*16)/9+am*am*(d[3]-c.mpf(5)/12),p=am*am*(d[4]+c.mpf(5)/2))
    def beta(s,center):
        u=(s-center)/c.mpf('.025')
        return c.mpf(0) if abs(u)>=1 else c.exp(-1/(1-u*u))/(c.mpf('.025')*f.normp)
    def densities(s):
        bs=[beta(s,c.mpf(q.numerator)/q.denominator) for q in f.centers]
        U=z*4+h[0]*bs[0]+h[1]*bs[2]
        E=am*(s**c.mpf('.1')+sum((h[i+2]*bs[i] for i in range(3)),h[0]*0))
        return dict(m=U/3,h=E,k=E*U/3,e=U*U/9-E*E/2,p=E*E/2)
    count=0
    for x in ('1.5','1.8'):
        result=f.evaluate(Z,x);xx=c.mpf(x);actual={}
        breaks=sorted(set([c.mpf(1),xx]+[c.mpf(v) for v in ('1.225','1.25','1.275','1.475','1.5','1.525','1.725','1.75','1.775') if c.mpf(v)<xx]))
        for name,rate in dict(m='1',h='1.5',k='1.5',e='1',p='0').items():
            r=c.mpf(rate);values=[]
            for k in range(3):
                v=initial[name][k]/xx**r+c.quad(lambda s:densities(s)[name][k]*s**(r-1),breaks)/xx**r
                values.append(v)
                lo,hi=_bounds(c,f._restore(result['generic_leading_fields'][name][k]))
                tol=c.mpf('1e-40')*max(1,abs(v));assert lo-tol<=v<=hi+tol,(x,name,k,v,lo,hi);count+=1
            actual[name]=J(c,values)
        M=actual['m']*3;Mz=J(c,[M[1],2*M[2]])
        bs=[beta(xx,c.mpf(q.numerator)/q.denominator) for q in f.centers]
        U=z*4+h[0]*bs[0]+h[1]*bs[2];delta=c.mpf('.0001')
        Q=((z*U*2-z*M*(1-delta)-(1-z*z)*Mz)/(1-z*z*delta)).truncate(1)
        for k,row in enumerate(result['physical_similarity']['Ur_over_sqrt_R_over2']):
            lo,hi=_bounds(c,f._restore(row));assert lo-c.mpf('1e-40')<=Q[k]<=hi+c.mpf('1e-40');count+=1
    return count


def run():
    start=time.perf_counter();manufactured=independent_known_implicit_family();raw=independent_raw_patch_integrals()
    print('Known coupled implicit Z2 family and independent raw five-moment patch integrals pass',flush=True)
    f=OriginalRmPatchNumeric();samples=[];closure=seams=memory=pressure=0;contractions=[]
    for Z in ('-.3','0','.3','.17','.23'):
        for x in ('1','Rh'):samples.append(f.evaluate(Z,x))
        if Z=='.3':samples.extend(f.evaluate(Z,x) for x in ('1.25','1.5','1.75','1.8'))
        print('Actual original leading implicit patch: Z='+Z,flush=True)
    for p in samples:
        c,iv=f.ctx,f.iv;first=p['implicit_certificate']['first']
        def upper(r):return c.mpf(mp.make_mpf(tuple(r['upper_exact_mpf_tuple'])))
        assert first['certified'] and first['self_map_strictly_inside'];contractions.append(c.nstr(upper(first['contraction_bound']),20))
        assert upper(first['contraction_bound'])<1 and upper(p['implicit_certificate']['second_row_contraction'])<1
        assert p['contract']['finite_N_not_used_to_resolve_leading_controls']
        assert not p['limits']['modified_patch_fields'];memory+=2
        if p['contract']['terminal_leading_five_function_identities']:
            for row in p['actual_partial_defects']:
                for v in row:assert _bounds(c,f._restore(v))==(0,0);closure+=1
        if p['native_x']=='1':
            before=f._inlet(p['Z'])
            for n,rows in before['generic_own_moment_enclosures'].items():
                # The old frame includes finite-N error; the leading branch
                # has the same unmodified source inside that enclosure.
                for a,b in zip(rows,p['generic_leading_fields'][n]):overlap(a,b);seams+=1
            assert p['finite_N_incoming_history_errors']==before['finite_N_signed_unknown_history_errors'];memory+=1
            assert p['positive_order_n1_own_history']==before['positive_order_n1_own_history'];memory+=1
        for k,row in enumerate(p['generic_leading_fields']['pressure']):
            axis=f._restore(p['independent_P0_over_Pstar_squared'][k]);inc=f._restore(p['generic_leading_fields']['p'][k])
            lo,hi=_bounds(c,axis+inc);ll,hh=_bounds(c,f._restore(row));assert max(lo,ll)<=min(hi,hh);pressure+=1
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,manufactured_coupled_implicit_Z2_coefficients_and_contractions=manufactured,
            independent_raw_five_moment_and_radial_recovery_coefficients=raw,actual_original_leading_patch_packets=len(samples),
            leading_terminal_five_identity_coefficients=closure,actual_Rm_inlet_generic_coefficients=seams,
            leading_vs_finite_N_history_checks=memory,same_P0_once_only_coefficients=pressure,
            maximum_actual_self_map_contraction=max(contractions,key=lambda x:f.ctx.mpf(x))),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(whole_axis_new_numeric_cover=False,finite_N_patch_transport=False,modified_patch=False,
            physical_Cartesian_assembly=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for n,d in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==d,n
    for n,d in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==d,n
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    r=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=r['checks'],execution_seconds=r['execution_seconds']),indent=2),flush=True)
