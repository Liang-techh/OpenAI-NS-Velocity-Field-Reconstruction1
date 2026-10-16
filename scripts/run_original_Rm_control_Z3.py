"""Execute genuine private Rm third-source rows and implicit controls."""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_Rm_control_Z3 import OriginalRmControlZ3,implicit_Z3
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_Rm_control_Z3_numeric.json'
SOURCE=PREFIX+'compliant_current_original_true_Rm_patch_C6_functions.py'


def independent_inverse():
    tree=ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='nonlinear')
    class Algebra:
        fixed=staticmethod(sy.sympify)
        add=staticmethod(lambda *v:sum(v))
        mul=staticmethod(lambda *v:sy.prod(v))
        sub=staticmethod(lambda a,b:a-b)
        scale=staticmethod(lambda value,factor:value*sy.Rational(factor))
    scope={};exec(compile(ast.Module(body=[node],type_ignores=[]),str(DIRECTORY/SOURCE),'exec'),scope)
    z=sy.Symbol('Z');h=[sy.Rational(i+1,10**6)*(1+z+(i+1)*z*z+z**3) for i in range(5)]
    alpha=(1+z*z)**2/10**6
    W=dict(L=[[1,1,0,0,0],[1,2,0,0,0],[0,0,1,1,1],[0,0,1,2,3],[0,0,1,4,9]],
           fg=[sy.Rational(1,10**6),sy.Rational(2,10**6)],gg=[1,1],
           ff=[sy.Rational(1,10**6)]*3,ff_over_x=[sy.Rational(i,10**6) for i in (1,2,3)])
    Q=scope['nonlinear'](Algebra(),W,h,alpha)
    defects=[-sum(W['L'][i][j]*h[j] for j in range(5))-Q[i] for i in range(5)]
    iv=MPIntervalContext();iv.dps=115;audit=mp.mp.clone();audit.dps=190
    def jet(expr,zz):
        values=[]
        for k in range(4):
            value=sy.cancel(sy.diff(expr,z,k).subs(z,zz)/math.factorial(k))
            values.append(iv.mpf(str(value.p))/int(value.q))
        return J(iv,values)
    cv=lambda v:iv.mpf(str(sy.Rational(v).p))/int(sy.Rational(v).q)
    Wi={name:[[cv(v) for v in row] for row in values] if name=='L' else [cv(v) for v in values] for name,values in W.items()}
    conditions=0
    for zz in (sy.Rational(-1,3),sy.Rational(0),sy.Rational(2,5)):
        result,proof=implicit_Z3(iv,J,Wi,[jet(d,zz) for d in defects],jet(alpha,zz),precision=160)
        assert _bounds(audit,proof['third_inverse_contraction'])[1]<1
        for i,row in enumerate(result):
            for k in range(4):
                exact=sy.cancel(sy.diff(h[i],z,k).subs(z,zz)/math.factorial(k))
                target=audit.mpf(str(exact.p))/int(exact.q);lo,hi=_bounds(audit,row[k])
                assert lo<=target<=hi,(i,k,str(lo),str(target),str(hi));conditions+=1
    return dict(passed=True,original_nonlinear_source_executed=True,
        independently_known_control_Taylor_rows_enclosed=conditions,known_cubic_controls_at_axial_points=3,
        third_coefficient_factorial_applied_once=True,source=SOURCE)


def run():
    start=time.perf_counter();fixture=independent_inverse()
    print('Original nonlinear map with independent known cubic-control inverse passes',flush=True)
    f=OriginalRmControlZ3(precision=85);audit=mp.mp.clone();audit.dps=190;samples=[]
    counts=dict(serialized_coefficient_errors=0,source_axial5_overlap_conditions=0,
        control_axial5_overlap_conditions=0,retained_lower_defect_rows=0,source_six_Z3_rows=0,
        third_defect_nonzero_enclosures=0,third_control_nonzero_enclosures=0,original_mean_pressure_seams=0,
        strict_third_inverse_conditions=0,angular_coefficient_before_decay_conditions=0)
    pressure_path=DIRECTORY/(PREFIX+'compliant_current_original_reference_m15_numeric.json')
    pressure=json.loads(pressure_path.read_bytes())
    assert pressure['checks']['passed'] and pressure['source_family']==f.family
    for base,group,target in ((DIRECTORY,pressure['input_hashes'],f.hashes),(ROOT,pressure['implementation_hashes'],f.code_hashes)):
        for name,digest in group.items():
            assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
            assert target.get(name,digest)==digest,name
            target[name]=digest
    oracle={(p['Z'],p['coordinate']['value']):p for p in pressure['pressure_row_samples']
        if p['source_cover']['generic_derivative_chart']=='restore_buffer'}
    def row_interval(row):return [audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]
    def archived(row):return [audit.make_mpf(tuple(row[k])) for k in ('lower_exact_mpf_tuple','upper_exact_mpf_tuple')]
    for Z in ('-.3','0','.3'):
        out=f.evaluate(Z);samples.append(out);data=f.coefficients(Z)
        assert f.evaluate(Z) is out and all(out['contract'].values())
        q=_bounds(audit,data['proof']['third_inverse_contraction'])[1];assert q<1;counts['strict_third_inverse_conditions']+=1
        cells=[p for p in f.archive['source_cells'] if audit.mpf(p['exact_Z_cell'][0])<=audit.mpf(Z)<=audit.mpf(p['exact_Z_cell'][1])]
        assert cells
        for group in ('normalized_defects_Z3','normalized_controls_Z3','physical_controls_Z3'):
            for jet in out[group]:
                assert len(jet)==4
                for value in jet:
                    lo,hi=row_interval(value);point=audit.mpf(value['approximation'])
                    assert max(abs(lo-point),abs(hi-point))<=audit.mpf(value['absolute_error_upper'])
                    counts['serialized_coefficient_errors']+=1
        for group in ('normalized_defects_Z3','normalized_controls_Z3'):
            for i,row in enumerate(out[group]):
                lo,hi=row_interval(row[3]);assert lo!=hi or lo!=0
                counts['third_defect_nonzero_enclosures' if group=='normalized_defects_Z3' else 'third_control_nonzero_enclosures']+=1
                for cell in cells:
                    source=cell['coefficient_solution']
                    jets=source['normalized_five_defect_axial5'] if group=='normalized_defects_Z3' else source['normalized_implicit_inverse_axial5']['controls']
                    for k in range(4):
                        a,b=row_interval(row[k]);x,y=archived(jets[i]['coefficients'][k]);assert max(a,x)<=min(b,y),(Z,group,i,k)
                        counts['source_axial5_overlap_conditions' if group=='normalized_defects_Z3' else 'control_axial5_overlap_conditions']+=1
        old=f.patch._inlet(Z)
        for i,row in enumerate(old['signed_unmodified_defects']):
            for k in range(3):
                a,b=row_interval(row[k]);x,y=_bounds(audit,data['source']['truth']['defects'][i][k])
                assert max(a,x)<=min(b,y),(Z,i,k);counts['retained_lower_defect_rows']+=1
        for rows in out['genuine_six_Rm_Z3'].values():assert len(rows)==4;counts['source_six_Z3_rows']+=4
        ref=oracle[Z,'1']['physical_pressure_rows']['m']['signed_background_ordinary_Z']
        for k,row in enumerate(out['physical_mean_at_Rm_Z3']):
            a,b=row_interval(row);x,y=row_interval(ref[k]);a*=math.factorial(k);b*=math.factorial(k)
            assert max(a,x)<=min(b,y),(Z,'mean',k);counts['original_mean_pressure_seams']+=1
        for proof in out['source_decay_proofs']:
            if proof.get('quantity','').startswith('genuine_'):
                assert proof['bound_is_not_signed_field_value'];counts['angular_coefficient_before_decay_conditions']+=1
        print('Original Rm genuine third source and same-Jacobian controls: Z='+Z,flush=True)
    rejected=0
    for Z in ('1','nan'):
        try:f.evaluate(Z)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid axial source admitted')
    hashes=dict(f.hashes);hashes[pressure_path.name]=hashlib.sha256(pressure_path.read_bytes()).hexdigest()
    code=dict(f.code_hashes);name='scripts/run_original_Rm_control_Z3.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_known_inverse=fixture,actual_original_packets=len(samples),**counts,
                    archive_is_independent_enclosure_not_source_values=True,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(original_leading_Rm_source_and_implicit_controls_Z3_with_absolute_errors=True,
            Rm_patch_pressure_rows_Z2=False,finite_N_patch_C3_transported=False,Rm_Rh_signed_m15_integrated=False,
            new_whole_axis_Z3_certificate=False,control_C6_numeric_certificate=False,
            global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
            stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for base,group in ((DIRECTORY,d['input_hashes']),(ROOT,d['implementation_hashes'])):
        for name,digest in group.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;data=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
