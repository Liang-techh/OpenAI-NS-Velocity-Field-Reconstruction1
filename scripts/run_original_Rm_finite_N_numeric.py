"""Actual selected-N patch transport and independent signed-density checks."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from decimal import Decimal,localcontext

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_Rm_finite_N_numeric import OriginalRmFiniteNNumeric,transport_step,canonical_defect_increment,RATES
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_Rm_finite_N_numeric.json'


def independent_signed_transport():
    """Integrate actual finite-N signed density differences in physical dx/x."""
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    c=mp.mp.clone();c.dps=60;iv=MPIntervalContext();iv.dps=105
    zp=J.variable(c,c.mpf('.23'),2);zt=J.variable(iv,iv.mpf('.23'),2)
    def source(ctx,z,y):
        N=ctx.mpf(13) # Manufactured fixture only; never the actual source N.
        E=(1+z*z)*ctx.mpf('.3')*(1+y/5)
        V=(z*4+(1+z)*ctx.mpf('.03')*y)/ctx.mpf('3.1')
        a=(1-z)*ctx.mpf('.1')*(1+y/3);b=(z*z-1)*ctx.mpf('.04')*(1-y/4)
        en=E*(a/N).exp();vn=V+b/N
        return dict(m=vn-V,h=en-E,k=en*vn-E*V,
            e=vn*vn-en*en/2-V*V+E*E/2,p=en*en/2-E*E/2)
    incomingp={n:(1+zp*zp/5)*c.mpf(str((i+1)*.017)) for i,n in enumerate(RATES)}
    incomingt={n:(1+zt*zt/5)*iv.mpf(str((i+1)*.017)) for i,n in enumerate(RATES)}
    history=dict(incomingt);checks=0
    bounds=('0','.2','.7','1')
    for a,b in zip(bounds,bounds[1:]):
        y=iv.mpf([a,b]);source_box=source(iv,zt,y)
        # The adapter consumes an unsigned range of an unknown signed driver.
        for name,rate in RATES.items():
            signed=source_box[name]
            cap=J(iv,[iv.mpf([-max(abs(q) for q in _bounds(c,v)),max(abs(q) for q in _bounds(c,v))])
                      for v in signed.coefficients])
            width=iv.ln(iv.exp(iv.mpf(b))/iv.exp(iv.mpf(a)))
            history[name],decay,mass=transport_step(iv,history[name],cap,width,rate)
            if rate=='0':
                assert _bounds(c,decay)==(1,1);checks+=1
            if b not in ('.7','1'):continue
            xp=c.exp(c.mpf(b));r=c.mpf(rate)
            for k in range(3):
                actual=incomingp[name][k]/xp**r+c.quad(lambda x:source(c,zp,c.ln(x))[name][k]*x**(r-1),[1,xp])/xp**r
                lo,hi=_bounds(c,history[name][k]);tol=c.mpf('1e-42')*max(1,abs(actual))
                assert lo-tol<=actual<=hi+tol,(b,name,k,actual,lo,hi);checks+=1
    # With no driver, retain the actual nonzero inlet with exactly one own
    # decay. This catches an extra1/N or a silent patch-seam restart.
    for name,rate in RATES.items():
        out,_,_=transport_step(iv,incomingt[name],J.constant(iv,0,2),iv.mpf('.4'),rate)
        expected=incomingp[name]*c.exp(-c.mpf(rate)*c.mpf('.4'))
        for k in range(3):
            lo,hi=_bounds(c,out[k]);tol=c.mpf('1e-48')
            assert lo-tol<=expected[k]<=hi+tol;checks+=1
    return checks


def independent_canonical_defects():
    """Invert raw generic m/e units before taking their signed differences."""
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    c=mp.mp.clone();c.dps=65;z=J.variable(c,c.mpf('.23'),2);P=c.mpf('3.1')
    am=(1+z*z).reciprocal()*c.exp(c.mpf('-.6'));old={};inc={}
    for i,n in enumerate(RATES):
        old[n]=(1+z*(i+1)+z*z)*c.mpf(str((i+1)*.07))
        inc[n]=(1-z*(i+1)+z*z*(i-2))*c.mpf(str((i+1)*-.003))
    def defects(fields,x):
        M=fields['m']*P
        return [(M-z*4)*x,
            (fields['k']*P-z*fields['h']*4)*x**c.mpf('1.5')/am,
            fields['h']*x**c.mpf('1.5')/am-x**c.mpf('1.6')*c.mpf(5)/8,
            (fields['e']-(z*M*8-z*z*16)/(P*P))*x/(am*am)+x**c.mpf('1.2')*c.mpf(5)/12,
            fields['p']/(am*am)-x**c.mpf('.2')*c.mpf(5)/2]
    checks=0
    for x in (c.mpf('1.25'),c.mpf('1.5'),c.exp(1)):
        before=defects(old,x);after=defects({n:old[n]+inc[n] for n in RATES},x)
        change=canonical_defect_increment(c,J,'.23',x,am,P,1/P,inc)
        for a,b,d in zip(before,after,change):
            for k in range(3):assert abs((b-a)[k]-d[k])<c.mpf('1e-48');checks+=1
    return checks


def run():
    start=time.perf_counter();independent=independent_signed_transport();canonical=independent_canonical_defects()
    print('Independent signed finite-N density, true dx/x and nonzero-inlet memory checks pass',flush=True)
    f=OriginalRmFiniteNNumeric();samples=[]
    for Z in ('-.3','0','.3','.17','.23'):
        for x in ('1','Rh'):samples.append(f.evaluate(Z,x))
        if Z=='.3':samples.extend(f.evaluate(Z,x) for x in ('1.25','1.5','1.75','1.8'))
        print('Actual original finite-N patch: Z='+Z,flush=True)
    samples.append(f.evaluate('.23','1.4'))
    with localcontext() as dc:
        dc.prec=180
        near=[format(Decimal('1.225')+sign*Decimal('1e-140'),'f') for sign in (-1,1)]
    for x,count in zip(near,(1,2)):
        p=f.evaluate('.3',x);assert len(p['physical_transport_steps'])==count;samples.append(p)
    seams=pressure=radial=memory=measures=terminal=0
    for p in samples:
        c,iv=f.ctx,f.iv;old=f._leading[p['Z'],p['native_x']]
        assert p['leading_controls_Z2']==old['controls_Z2'];memory+=1
        assert p['selected_integer']==old['selected_integer'];memory+=1
        assert p['positive_order_n1_own_history']==old['positive_order_n1_own_history'];memory+=1
        assert p['contract']['finite_incoming_already_contains_one_over_N']
        assert p['contract']['new_density_divided_by_N_once']
        assert not p['limits']['modified_terminal_five_moment_closure']
        assert p['source_cover']['unique_native_cells']==11
        total=iv.mpf(0)
        for step in p['physical_transport_steps']:
            total+=f.owner._decode(step['physical_log_width'])
            weight=step['own_weights']['p'];assert _bounds(c,f.owner._decode(weight['decay']))==(1,1);measures+=1
        expected=f.leading._restore(p['coordinate']['log_x'])
        lo,hi=_bounds(c,total);ll,hh=_bounds(c,expected);assert max(lo,ll)<=min(hi,hh);measures+=1
        if p['native_x']=='Rh':
            assert len(p['physical_transport_steps'])==11
            for rows in p['actual_canonical_partial_defect_enclosures']:
                lo,hi=_bounds(c,f.leading._restore(rows[0]));assert lo<0<hi;terminal+=1
        if p['native_x']=='1':
            assert not p['physical_transport_steps']
            before=f.leading._inlet(p['Z'])
            for name,rows in before['generic_own_moment_enclosures'].items():
                for a,b in zip(rows,p['generic_modified_fields'][name]):overlap(a,b);seams+=1
            for n,rows in p['finite_N_history_errors'].items():
                for a,b in zip(rows,p['finite_N_nonzero_incoming_errors'][n]):overlap(a,b);memory+=1
        for k,row in enumerate(p['generic_modified_fields']['pressure']):
            axis=f.leading._restore(p['independent_P0_over_Pstar_squared'][k]);inc=f.leading._restore(p['generic_modified_fields']['p'][k])
            lo,hi=_bounds(c,axis+inc);ll,hh=_bounds(c,f.leading._restore(row));assert max(lo,ll)<=min(hi,hh);pressure+=1
        # Recover radial velocity independently from the complete generic
        # fields, including the full modified own mass derivative.
        J=f.J;z=J.variable(iv,iv.mpf(p['Z']),2)
        v=f._jet(p['generic_modified_fields']['V'],True);m=f._jet(p['generic_modified_fields']['m'],True)
        mz=J(iv,[m[1],2*m[2]]);delta=iv.mpf(f.leading.reference.leading.core.delta_iv)
        Q=((z*v*2-z*m*(1-delta)-(1-z*z)*mz)/(1-z*z*delta)*iv.exp(f.leading.defect.logP)).truncate(1)
        for k,row in enumerate(p['physical_similarity']['Ur_over_sqrt_R_over2']):
            lo,hi=_bounds(c,Q[k]);ll,hh=_bounds(c,f.leading._restore(row));assert max(lo,ll)<=min(hi,hh);radial+=1
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_signed_density_dx_over_x_and_memory_coefficients=independent,
            independent_raw_canonical_defect_difference_coefficients=canonical,exact_close_cell_boundary_classifications=2,
            actual_original_finite_N_patch_packets=len(samples),actual_Rm_inlet_own_coefficients=seams,
            fixed_leading_N_primitive_and_incoming_history_checks=memory,
            unique_native_measure_and_pressure_memory_checks=measures,
            independent_generic_radial_recovery_coefficients=radial,same_P0_once_only_coefficients=pressure,
            terminal_finite_N_defects_retained_as_nonzero_error_intervals=terminal),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(signed_finite_N_centers_resolved=False,physical_Cartesian_assembly=False,
            fast_phase_or_total_native_y_rows=False,modified_terminal_five_moment_closure=False,
            global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
