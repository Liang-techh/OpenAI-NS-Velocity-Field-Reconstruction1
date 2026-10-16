"""Execute changing-source restoration and postrestore m15 through Rm."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_restore_m15 import OriginalRestoreM15,pressure_polynomials
from openai_ns_reconstruction.original_restore_kernels import OriginalRestoreKernels,powers,MONOMIALS
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_restore_m15_numeric.json'
SOURCE=PREFIX+'compliant_current_original_n1_preRp_moment_functions.py'


def original_identity():
    z,de=sy.symbols('Z delta');A,D=1-de*z*z,1-z*z
    tree=ast.parse((DIRECTORY/SOURCE).read_text(encoding='utf-8'))
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='pressure_functions')
    class Algebra:
        def __init__(self,g):pass
        fixed=staticmethod(sy.sympify)
        add=staticmethod(lambda *v:sum(v))
        neg=staticmethod(lambda v:-v)
        sub=staticmethod(lambda a,b:a-b)
        mul=staticmethod(lambda *v:sy.prod(v))
        scale=staticmethod(lambda a,b:sy.Rational(a)*b)
        div=staticmethod(lambda a,b:a/b)
    scope=dict(C1Algebra=Algebra,source=SimpleNamespace(C1Function=lambda v,d:v))
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(DIRECTORY/SOURCE),'exec'),scope)
    names=('q','q_Z','q_squared','U_q_Z','U_squared','U_U_Z')
    integrals={f'I{p}_{n}':sy.Symbol(f'I{p}_{n}') for p in (0,1) for n in names}
    enames=('R_q','R_squared_q','R_q_squared','R_squared_q_squared','R_U_squared',
            'R_squared_U_squared','q_plus_q_y','R_q_y')
    low={n:sy.Symbol('low_'+n) for n in enames};up={n:sy.Symbol('up_'+n) for n in enames}
    value=scope['pressure_functions'](SimpleNamespace(symbol=lambda name:z,one=1),de,integrals,low,up)['m15']
    bound=lambda name:up[name]-low[name]
    expected=(bound('R_squared_q')/(2*A)+bound('R_squared_q_squared')/4
        -z*z*bound('R_squared_U_squared')/A**2-integrals['I1_q']/(2*A)
        +z*(1-de)*integrals['I1_q_Z']/(4*A)-integrals['I1_q_squared']/4
        +D*integrals['I1_U_q_Z']/(2*A)+z*z*(1-de)*integrals['I1_U_squared']/A**2
        +z*D*integrals['I1_U_U_Z']/A**2-bound('R_q_y'))
    assert sy.cancel(value-expected)==0
    return dict(passed=True,original_reduced_S1_identity_executed=True,
        existing_F_G_coefficients_already_include_half=True,source=SOURCE)


def independent_fixture():
    c=mp.mp.clone();c.dps=100;iv=MPIntervalContext();iv.dps=160
    z,t,eps=sy.symbols('Z t epsilon');de=sy.Rational(1,17);j=sy.Rational(7,5)
    psi=1+z+z*z/3+z**3/7;memory=sy.Rational(2,3)-z/2+z*z/5-z**3/11-psi
    alpha=1-3*t*t+2*t**3;gamma=2*t**3-9*t*t+18*t-17+18*sy.exp(-t)
    U=4*z+(j+eps*psi)*alpha;M=4*z+(j+eps*psi)*gamma+eps*memory*sy.exp(-t)
    A,D=1-de*z*z,1-z*z;Q=(2*z*U-(1-de)*z*M-D*sy.diff(M,z))/A
    F=Q/(2*A)+Q*Q/4-z*z*U*U/A**2
    G=(-Q/(2*A)+z*(1-de)*sy.diff(Q,z)/(4*A)-Q*Q/4+D*U*sy.diff(Q,z)/(2*A)
       +z*z*(1-de)*U*U/A**2+z*D*U*sy.diff(U,z)/A**2)
    assert sy.simplify(sy.diff(gamma,t)-(alpha-gamma))==0
    assert sy.simplify(sy.diff(M,t)-(U-M))==0
    physical=sy.exp(2*t)*(2*F+G+sy.diff(F,t))
    direct=[sy.lambdify((t,z,eps),sy.diff(physical,z,k),'mpmath') for k in (0,1)]
    ends=[sy.lambdify((t,z,eps),sy.diff(sy.exp(2*t)*F,z,k),'mpmath') for k in (0,1)]
    scalar=OriginalRestoreKernels(SimpleNamespace(ctx=c,iv=iv),cells=64)
    scalar.alpha=lambda ctx,x:1-3*x*x+2*x*x*x
    coefficients=integrals=inclusions=gamma_checks=0
    # Independent scalar targets need more precision than the160-digit
    # directed enclosure, especially the exact constant monomial.
    with mp.workdps(190):
        nodes,weights=mp.gauss_quadrature(64,'legendre')
        al=lambda x:1-3*x*x+2*x**3
        ga=lambda x:2*x**3-9*x*x+18*x-17+18*mp.exp(-x)
        for end in (mp.mpf('.1'),mp.mpf('.7'),mp.mpf(1)):
            exact={}
            for code in MONOMIALS:
                i,jj,k=powers(code)
                exact[code]=mp.quad(lambda x:mp.exp((2-k)*x)*al(x)**i*ga(x)**jj,[0,end/2,end])
            ks=scalar.evaluate(c.mpf(end),iv.mpf(str(end)))
            for code,value in exact.items():
                lo,hi=_bounds(c,ks['truth'][code]);assert lo<=value<=hi,(str(end),code,str(value),str(lo),str(hi));inclusions+=1
            gl,gh=_bounds(c,ks['gamma_truth']);assert gl<=ga(end)<=gh;gamma_checks+=1
            for zz in ('-.3','0','.3'):
                def jet(expr):return J(c,[c.mpf(str(sy.N(sy.diff(expr,z,k).subs(z,sy.Rational(zz)),100)))/c.factorial(k) for k in range(4)])
                zj=J.variable(c,c.mpf(zz),3);pj,dj=jet(psi),jet(memory)
                poly=pressure_polynomials(J,c,zz,c.mpf(1)/17,
                    {(0,0):4*zj,(0,1):J.constant(c,c.mpf(7)/5,3),(1,1):pj},
                    {(0,0):4*zj,(0,3):J.constant(c,c.mpf(7)/5,3),(1,3):pj,(1,9):dj})
                ee=mp.mpf('.23')
                for k in (0,1):
                    gpart=sum((c.mpf(row[k])*ee**d*exact[code] for (d,code),row in poly['G'].items()),mp.mpf(0))
                    result=ends[k](end,mp.mpf(zz),ee)-ends[k](0,mp.mpf(zz),ee)+gpart
                    want=sum((w*direct[k]((n+1)*end/2,mp.mpf(zz),ee)*end/2 for n,w in zip(nodes,weights)),mp.mpf(0))
                    assert abs(result-want)<mp.mpf('1e-65')*max(1,abs(want));integrals+=1
                coefficients+=len(poly['G'])
    return dict(passed=True,changing_source_physical_value_Z1_integrals=integrals,
        independent_scalar_monomial_inclusions=inclusions,independent_gamma_inclusions=gamma_checks,
        scalar_polynomial_coefficients_used=coefficients,exact_varying_U_mean_ODE=True)


def run():
    start=time.perf_counter();identity,fixture=original_identity(),independent_fixture()
    print('Original reduced-S1 source and independent changing-U integral fixture pass',flush=True)
    restore=OriginalRestoreM15(precision=85,cells=128);post=restore.next_field();samples=[]
    audit=mp.mp.clone();audit.dps=190
    counts=dict(exact_inlet_replays=0,serialized_new_coefficient_errors=0,exact_empty_rows=0,
        lower_finite_N_endpoint_conditions=0,actual_changing_RQy_rows=0,atomic_offset_conditions=0,
        scalar_kernel_conditions=0,retained_reference_decay_proofs=0)
    def ends(row):return [audit.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]
    def check(field,Z,coordinate):
        p=field.evaluate(Z,coordinate);samples.append(p);n=124 if field.chart=='restoration' else 140
        assert len(p['m15_C0_C1_terms'])==n+16
        assert p['m15_C0_C1_terms'][:n]==p['incoming_m15_terms']==field.inlet(Z)['m15_C0_C1_terms']
        assert p['actual_selected_integer']==field.inlet(Z)['actual_selected_integer']
        assert all(p['contract'].values()) and field.evaluate(Z,coordinate) is p;counts['exact_inlet_replays']+=1
        for term in p['m15_C0_C1_terms'][n:]:
            for value in term['coefficient_C0_C1']:
                lo,hi=ends(value);v=audit.mpf(value['approximation'])
                assert max(abs(lo-v),abs(hi-v))<=audit.mpf(value['absolute_error_upper'])
                counts['serialized_new_coefficient_errors']+=1
                if p['exact_empty_increment']:
                    assert lo==hi==v==0;counts['exact_empty_rows']+=1
            if term['name'].endswith('_error'):continue
            atoms=term['log_scale_terms'];power=1 if 'RQy_' in term['name'] else 2
            assert audit.mpf(atoms[3]['approximation'])==power*(-8 if field.chart=='restoration' else -7)
            counts['atomic_offset_conditions']+=1
            if 'RQy_' in term['name'] and not p['exact_empty_increment']:
                counts['actual_changing_RQy_rows']+=1
        if not p['exact_empty_increment']:
            lower=p['actual_lower_finite_N_endpoint_budget'];radius_hi=ends(lower['log_radius'])[1]
            for name,power in (('endpoint',2),('viscous',1)):
                for k in (0,1):
                    v=lower[name][k]
                    if v['exact_zero']:continue
                    want=audit.make_mpf(tuple(v['log_absolute_upper']['upper_exact_mpf_tuple']))+power*radius_hi
                    got=audit.make_mpf(tuple(p['finite_N_error_logs_in_physical_units'][name][k]['log_absolute_upper']['upper_exact_mpf_tuple']))
                    assert got>=want or abs(got-want)<=abs(got)*audit.mpf('1e-150')
                    counts['lower_finite_N_endpoint_conditions']+=1
            ks=p['scalar_kernel_receipt'];assert all(ks['contract'].values())
            counts['scalar_kernel_conditions']+=len(ks['monomials'])
        counts['retained_reference_decay_proofs']+=len(p['reference_mean_memory_proofs'])
        assert not p['limits']['genuine_temporal_n_ge_2']
        return p
    for Z in ('-.3','0','.3'):
        for field in (restore,post):
            for q in ('0','.25','.75','1'):check(field,Z,q)
        print('Changing-source140-term restoration and156-term postrestore m15: Z='+Z,flush=True)
    # Independent axial point executes the uncached complete moment inlet.
    for field in (restore,post):check(field,'.17','.3')
    original_scalar_checks=0
    for coordinate in ('.25','.75','1'):
        t=restore.ctx.mpf(coordinate);q=restore.iv.mpf(coordinate);ks=restore.kernels.evaluate(t,q)
        kp,kt=restore.rows.reference.kernels(coordinate)['mean']
        exact=restore.ctx.exp(-t)+kp;lo,hi=_bounds(audit,ks['gamma_truth'])
        assert lo<=audit.mpf(exact)<=hi;original_scalar_checks+=1
    # Refine only a scalar quadrature, without rerunning accepted ancestors.
    fine=OriginalRestoreKernels(restore.rows.reference,cells=256)
    coarse=restore.kernels.evaluate(restore.ctx.mpf('.75'),restore.iv.mpf('.75'))
    refined=fine.evaluate(restore.ctx.mpf('.75'),restore.iv.mpf('.75'))
    contractions=0
    for code in MONOMIALS:
        a,b=_bounds(audit,coarse['truth'][code]);x,y=_bounds(audit,refined['truth'][code])
        assert y-x<=b-a or abs((y-x)-(b-a))<audit.mpf('1e-150');contractions+=1
    # Directed decimal endpoint segment is included even when point rounds.
    long_decimal='0.'+'3'*100
    uncertain=restore.kernels.evaluate(restore.ctx.mpf(long_decimal),restore.iv.mpf(long_decimal))
    assert uncertain['uncertain_endpoint_segment']
    rejected=0
    for field,Z,q in [(restore,'1','.5'),(restore,'nan','.5'),(restore,'.3','-.1'),(post,'.3','1.1'),(post,'.3','g:.1')]:
        try:field.evaluate(Z,q)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid restoration moment query admitted')
    hashes=dict(restore.hashes);hashes.update(post.hashes);code=dict(restore.code_hashes);code.update(post.code_hashes)
    name='scripts/run_original_restore_m15.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=restore.family,samples=samples,input_hashes=hashes,implementation_hashes=code,
        checks=dict(passed=True,original_pressure_identity=identity,independent_changing_source_fixture=fixture,
            actual_packets=len(samples),**counts,original_flat_source_gamma_inclusions=original_scalar_checks,
            refined_scalar_width_contractions=contractions,uncertain_endpoint_segment_retained=True,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(restoration_postrestore_cumulative_signed_m15_with_complete_absolute_errors=True,
            cumulative_pressure_moment_to_Rm=True,resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,common_global_N=False,
            exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
