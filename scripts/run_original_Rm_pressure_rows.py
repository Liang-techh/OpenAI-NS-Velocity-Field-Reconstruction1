"""Evaluate original active-patch pressure and genuine C3 finite-N memories."""
import argparse
import ast
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_Rm_pressure_rows import OriginalRmPressureRows
from openai_ns_reconstruction.original_Rm_finite_N_numeric import RATES,transport_step
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_Rm_pressure_rows_numeric.json'
SOURCE=PREFIX+'compliant_current_original_true_Rm_patch_C6_functions.py'


def independent_source():
    tree=ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('nonlinear','patch_fields')]
    class Graph:
        zero,one=sy.Integer(0),sy.Integer(1)
        symbol=staticmethod(sy.Symbol)
        constant=staticmethod(sy.Rational)
        mul=staticmethod(lambda *v:sy.prod(v))
        quotient=staticmethod(lambda a,b,reason:a/b)
    class Algebra:
        def __init__(self,g,order):pass
        fixed=staticmethod(sy.sympify)
        trim=staticmethod(sy.sympify)
        add=staticmethod(lambda *v:sum(v))
        mul=staticmethod(lambda *v:sy.prod(v))
        sub=staticmethod(lambda a,b:a-b)
        scale=staticmethod(lambda a,b:a*sy.sympify(b))
    scope=dict(Algebra=Algebra,FunctionJet=lambda values:values[0],Fraction=Fraction,
               power=lambda g,x,p:x**sy.Rational(p.numerator,p.denominator))
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(DIRECTORY/SOURCE),'exec'),scope)
    z=sy.Symbol('Z');x=sy.Rational(63,50);delta=sy.Rational(1,17);P=sy.Integer(17)
    h=[sy.Rational(i+1,1000)*(1+z+sy.Rational(i+2,7)*z*z+z**3/11) for i in range(5)]
    defects=[sy.Rational(i+1,300)*(1-z+z*z/5-z**3/13) for i in range(5)]
    W=dict(L=[[sy.Rational(i+j+1,23) for j in range(5)] for i in range(5)],
           fg=[sy.Rational(1,7)]*2,gg=[sy.Rational(1,5)]*2,
           ff=[sy.Rational(1,11)]*3,ff_over_x=[sy.Rational(1,13)]*3)
    # The true first row contains only the two axial mass integrals.
    W['L'][0][2:]=[sy.Integer(0)]*3
    c=mp.mp.clone();c.dps=100
    norm=c.quad(lambda s:c.exp(-1/(1-s*s)) if abs(s)<1 else c.mpf(0),[-1,0,1])
    gammas=[]
    for center in (sy.Rational(5,4),sy.Rational(3,2),sy.Rational(7,4)):
        u=40*(x-center)
        if abs(u)>=1:gammas.append((sy.Integer(0),sy.Integer(0)))
        else:
            value=40*sy.exp(-1/(1-u*u))/sy.Float(str(norm),100)
            gammas.append((value,-80*u*value/(1-u*u)**2))
    am=sy.exp(-sy.Rational(3,5))/(1+z*z)
    raw,radial,proof=scope['patch_fields'](Graph(),x,am,P,0,defects,h,1/(P*P*am*am),W,gammas)
    U,M,Uy=raw['V']*P,raw['m']*P,radial['V_y']*P
    A,D=1-delta*z*z,1-z*z
    Q=(2*z*U-(1-delta)*z*M-D*sy.diff(M,z))/A
    Qy=(2*z*Uy-(1-delta)*z*(U-M)-D*sy.diff(U-M,z))/A
    # Independent native derivative: x*d/dx, with dW_mass/dx=beta.
    xx=sy.Symbol('x',positive=True);w0,w2=sy.Function('w0')(xx),sy.Function('w2')(xx)
    b0,b2=sy.Function('b0')(xx),sy.Function('b2')(xx)
    uu=4*z+h[0]*b0+h[1]*b2;mm=4*z+(defects[0]+h[0]*w0+h[1]*w2)/xx
    replacements={sy.diff(w0,xx):b0,sy.diff(w2,xx):b2}
    assert sy.expand((xx*sy.diff(mm,xx)-(uu-mm)).subs(replacements))==0
    qq=(2*z*uu-(1-delta)*z*mm-D*sy.diff(mm,z))/A
    qy=(2*z*xx*sy.diff(uu,xx)-(1-delta)*z*(uu-mm)-D*sy.diff(uu-mm,z))/A
    assert sy.cancel((xx*sy.diff(qq,xx)-qy).subs(replacements))==0
    def jet(expr,zz):
        return J(c,[c.mpf(str(sy.N(sy.diff(expr,z,k).subs(z,zz),100)))/math.factorial(k) for k in range(4)])
    fake=object.__new__(OriginalRmPressureRows);fake.J=J;fake.ctx=c
    fake.o=SimpleNamespace(core=SimpleNamespace(delta=c.mpf(1)/17))
    Wp={name:[[c.mpf(str(v.p))/int(v.q) for v in row] for row in values] if name=='L'
        else [c.mpf(str(v.p))/int(v.q) for v in values] for name,values in W.items()}
    fake.patch=SimpleNamespace(defect=SimpleNamespace(logPp=c.log(17)),
        centers=(Fraction(5,4),Fraction(3,2),Fraction(7,4)),normp=norm,_W=lambda *args:Wp)
    exact=dict(U=U,m=M,Q=Q,U_y=Uy,Q_y=Qy,E=raw['E'],E_y=radial['E_y'])
    count=0
    for zz in (sy.Rational(-3,10),sy.Integer(0),sy.Rational(3,10)):
        data=dict(point=[jet(v,zz) for v in h],source=dict(point=dict(defects=[jet(v,zz) for v in defects])))
        fake.controls=SimpleNamespace(coefficients=lambda Z:data)
        got=fake._fields(str(c.mpf(str(zz.p))/int(zz.q)),c,c.mpf(63)/50,False,False)
        for name,expr in exact.items():
            for k in range(got[name].order+1):
                want=c.mpf(str(sy.N(sy.diff(expr,z,k).subs(z,zz),100)))/math.factorial(k)
                assert abs(got[name][k]-want)<c.mpf('1e-85')*max(1,abs(want)),(name,k)
                count+=1
    # Independent finite-rate integrals, including the pressure rate-zero case.
    duhamel=0
    for rate in RATES.values():
        r=c.mpf(rate);width=c.mpf('.37');incoming=c.mpf('.23');driver=c.mpf('.41')
        got,decay,mass=transport_step(c,incoming,driver,width,rate)
        want=incoming*c.exp(-r*width)+c.quad(lambda t:driver*c.exp(-r*(width-t)),[0,width])
        assert abs(got-want)<c.mpf('1e-95');duhamel+=1
    return dict(passed=True,original_patch_fields_source_executed=True,
        independent_finite_pressure_Taylor_comparisons=count,native_mean_ODE_identity=True,
        independent_native_Qy_derivative_identity=True,independent_Duhamel_integrals=duhamel,source=SOURCE)


def run():
    start=time.perf_counter();fixture=independent_source()
    print('Original patch source / independent native pressure derivative / Duhamel fixtures pass',flush=True)
    f=OriginalRmPressureRows(precision=85);audit=mp.mp.clone();audit.dps=190;samples=[]
    counts=dict(serialized_pressure_rows=0,serialized_error_enclosures=0,original_inlet_seams=0,
        full_C3_history_rows=0,positive_native_weights=0,ordinary_factorial_conditions=0,
        prefix_uniform_inclusions=0,leading_terminal_identities=0,nonzero_terminal_finite_N_rows=0,
        original_lower_pressure_overlaps=0,slow_fast_unit_conditions=0,exact_inlet_history_replays=0,
        merged_Z_width_conditions=0,full_weight_residual_diagnostics=0)
    def bounds(row):return tuple(audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])
    def highlog(row):return audit.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
    oldpath=DIRECTORY/(PREFIX+'compliant_current_original_Rm_patch_numeric.json')
    old=json.loads(oldpath.read_bytes());assert old['checks']['passed'] and old['source_family']==f.family
    f.hashes[oldpath.name]=hashlib.sha256(oldpath.read_bytes()).hexdigest()
    oldrows={(r['Z'],r['native_x']):r for r in old['samples']}
    for Z in ('-.3','0','.3'):
        for x in ('1','1.25','1.26','1.4','1.5','1.75','1.775','Rh'):
            out=f.evaluate(Z,x);samples.append(out);point,truth,geo=f.background(Z,x)
            assert f.evaluate(Z,x) is out and all(out['contract'].values())
            for name,rows in out['physical_pressure_rows'].items():
                for k,value in enumerate(rows['signed_background_ordinary_Z']):
                    actual=truth[name][k]*math.factorial(k)
                    assert value['source_enclosure_exact_mpf_tuples']==[list(v) for v in actual._mpi_]
                    lo,hi=bounds(value);p=audit.mpf(value['approximation'])
                    assert max(abs(lo-p),abs(hi-p))<=audit.mpf(value['absolute_error_upper'])
                    counts['serialized_pressure_rows']+=1;counts['ordinary_factorial_conditions']+=1
                assert len(rows['finite_N_absolute_error_log_bounds'])==len(rows['signed_background_ordinary_Z'])
                for error in rows['finite_N_absolute_error_log_bounds']:
                    assert not error['exact_zero'];counts['serialized_error_enclosures']+=1
            for name,rows in out['finite_N_history_error_logs'].items():
                assert len(rows)==4 and not any(row['exact_zero'] for row in rows)
                counts['full_C3_history_rows']+=4
            hist=f.history_logs(Z,x);width=sum((step['actual_log_width'] for step in hist['steps']),f.iv.mpf(0))
            # Directed sums enclose exactly log(x), regardless of overlapping Z frames.
            wl,wh=_bounds(audit,width);target=audit.mpf(1) if x=='Rh' else audit.log(audit.mpf(x))
            assert wl<=target<=wh,(Z,x,'width');counts['merged_Z_width_conditions']+=1
            seen=set()
            for step in hist['steps']:
                cell=(str(step['exact_left']),str(step['exact_right']));assert cell not in seen;seen.add(cell)
                for name,values in step['weights'].items():
                    assert _bounds(audit,values['positive_mass'])[0]>=0
                    if name=='p':assert values['positive_mass']==step['actual_log_width']
                    counts['positive_native_weights']+=1
            if x=='1':
                assert not hist['steps'] and out['finite_N_history_error_logs']==out['finite_N_incoming_history_error_logs']
                counts['exact_inlet_history_replays']+=1
                ref=f.inlet_rows.evaluate(Z,'1')['physical_pressure_rows']
                for name,values in out['physical_pressure_rows'].items():
                    for a,b in zip(values['signed_background_ordinary_Z'],ref[name]['signed_background_ordinary_Z']):
                        assert overlap(bounds(a),bounds(b)),(Z,name,'inlet');counts['original_inlet_seams']+=1
            uniform=f.uniform_background(Z,x)
            for name in ('U','m','Q','U_y','Q_y'):
                for k in range(truth[name].order+1):
                    ul,uh=_bounds(audit,uniform[name][k]);pl,ph=_bounds(audit,truth[name][k])
                    # Point-query intervals may be marginally wider by final rounding;
                    # require full overlap and inclusion of the independent point proxy.
                    assert max(ul,pl)<=min(uh,ph),(Z,x,name,k,'uniform overlap')
                    pv=audit.mpf(point[name][k]);assert ul<=pv<=uh,(Z,x,name,k,'uniform point')
                    counts['prefix_uniform_inclusions']+=1
            if x in ('1.775','Rh'):
                assert all(v==0 for row in truth['partial_defects'] for v in row.coefficients)
                for row in truth['unrefined_partial_defects']:
                    for value in row.coefficients:
                        lo,hi=_bounds(audit,value);assert lo<=0<=hi;counts['full_weight_residual_diagnostics']+=1
                assert truth['U_y'][0]==truth['Q_y'][0]==0
                assert all(not row['exact_zero'] for row in out['physical_pressure_rows']['Q']['finite_N_absolute_error_log_bounds'])
                counts['leading_terminal_identities']+=1;counts['nonzero_terminal_finite_N_rows']+=3
            if (Z,x) in oldrows:
                oldp=oldrows[Z,x]['physical_similarity']
                for name,oldname in (('U','Uz'),('m','Mz_over_R'),('Q','Ur_over_sqrt_R_over2'),('U_y','Uz_y')):
                    for k,row in enumerate(oldp[oldname][:len(out['physical_pressure_rows'][name]['signed_background_ordinary_Z'])]):
                        target=tuple(v/math.factorial(k) for v in bounds(out['physical_pressure_rows'][name]['signed_background_ordinary_Z'][k]))
                        assert overlap(target,bounds(row)),(Z,x,name,k,'old lower');counts['original_lower_pressure_overlaps']+=1
            units=out['finite_N_axial_radial_error_components']
            for k,key in enumerate(('B_phi','B_Z_phi')):
                source=f.read_log(f.loop['phase_and_first_slow_phase_log_bounds'][key])+f.o.core.lp_iv
                assert highlog(units['fast_B_phi'][k])==_bounds(audit,source)[1]
                counts['slow_fast_unit_conditions']+=1
            assert not out['limits']['active_patch_signed_m15_integrated']
        print('Original Rm patch pressure / genuine C3 finite-N histories: Z='+Z,flush=True)
    rejected=0
    for Z,x in (('1','1'),('nan','1'),('0','0'),('0','2.8'),('0','nan')):
        try:f.evaluate(Z,x)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original patch query admitted')
    code=dict(f.code_hashes);name='scripts/run_original_Rm_pressure_rows.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_source_and_native_pressure=fixture,actual_original_packets=len(samples),
                    **counts,rejected_invalid_queries=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(actual_patch_UQ_Z2_mean_Z3_UyQy_Z1_with_complete_absolute_errors=True,
            all_five_finite_N_patch_histories_C3_transported=True,active_patch_signed_m15_integrated=False,
            resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            new_whole_axis_C3_certificate=False,global_five_moment_numeric_join=False,common_global_N=False,
            exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    data=json.loads(path.read_text(encoding='utf-8'));assert data['checks']['passed']
    for base,group in ((DIRECTORY,data['input_hashes']),(ROOT,data['implementation_hashes'])):
        for name,digest in group.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;data=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
