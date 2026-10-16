"""Signed original core m15: source identity, direct density and error evidence."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_n1_core_m15_numeric import OriginalN1CoreM15Numeric,signed_prefix_rows

NAME=PREFIX+'compliant_current_original_n1_core_m15_numeric.json'


def original_density_identity():
    """Execute the original five-density body, then check physical IBP."""
    import sympy as sy
    r,z,eps,P2=sy.symbols('rho Z epsilon Pstar_squared',positive=True)
    u0=1+z+r*z*z/7
    u=r*(2-z+z*z)/5+r*r*(1+z)/11
    p=-r*(1+z*z)/3+r**3*(2-z)/13
    leadingF,F1=sy.symbols('leadingF F1')
    original_path=DIRECTORY/(PREFIX+'compliant_current_original_n1_inner_cutoff_functions.py')
    tree=ast.parse(original_path.read_text(encoding='utf-8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='signed_rho_integrands')
    algebra=SimpleNamespace(scale=lambda a,b:a*b,mul=lambda a,b:a*b,
        add=lambda *v:sy.Add(*v),neg=lambda a:-a)
    env=dict(c3=SimpleNamespace(C3Algebra=lambda g:algebra))
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<original signed five rho densities>','exec'),env)
    source=env['signed_rho_integrands'](SimpleNamespace(constant=lambda v:sy.Rational(str(v)),mul=lambda *v:sy.Mul(*v)),
        dict(F_rho_rows=[leadingF],U_rho_rows=[u0]),dict(F1=[F1],Uz1=[eps*u]),
        dict(rho=r,N1p=P2*sy.diff(p,r)-2*leadingF*F1),eps)['m15']
    primitive=eps**2*2*sy.integrate(u0*u,(r,0,r))+eps**2*P2*(sy.integrate(p,(r,0,r))-r*p)
    count=0
    for k in (0,1):
        assert sy.simplify(sy.diff(primitive,r,z,k)-sy.diff(source,z,k))==0;count+=1
    assert sy.simplify(primitive.subs(r,0))==0;count+=1
    # A nonzero constant pressure gauge must cancel from this weighted moment.
    gauge=sy.Symbol('constant_pressure_gauge')
    assert sy.simplify(sy.integrate(p+gauge,(r,0,r))-r*(p+gauge)-(sy.integrate(p,(r,0,r))-r*p))==0;count+=1
    return dict(actual_original_density_C0_C1_and_IBP_checks=count,
        source_SHA256=hashlib.sha256(original_path.read_bytes()).hexdigest(),
        original_epsilon_and_Pstar_squared_weights=True,finite_pressure_endpoint_required=True)


def direct_polynomial_density(f,Z,degree,kernel):
    """Independent Gauss quadrature of 2Uu and -rho*p_rho, including Z1."""
    c=f.ctx;r=c.mpf(4)
    total,_,_,_,_=f.n1._prepare(Z,degree,kernel,1)
    _,data=f.n1._forcing_polynomial(Z,degree,False)
    def value(rows,x,k,derivative=False):
        if derivative:return c.fsum(n*row[k]*x**(n-1) for n,row in enumerate(rows) if n)
        return c.fsum(row[k]*x**n for n,row in enumerate(rows))
    results=[]
    for order in (64,96):
        nodes,weights=c.gauss_quadrature(order,'legendre')
        values=[[c.mpf(0)]*2 for _ in range(2)]
        for node,weight in zip(nodes,weights):
            x=(node+1)*r/2;w=weight*r/2
            U=[value(data['rows']['Uz'],x,k) for k in (0,1)]
            u=[value(total['u'].rows,x,k) for k in (0,1)]
            pr=[value(total['p'].rows,x,k,True) for k in (0,1)]
            values[0][0]+=2*U[0]*u[0]*w
            values[0][1]+=2*(U[1]*u[0]+U[0]*u[1])*w
            for k in (0,1):values[1][k]+=-x*pr[k]*w
        results.append(values)
    got=f.evaluate(Z,'4',degree=degree,kernel_degree=kernel)
    checks=0;convergence=[]
    for sector in range(2):
        for k in (0,1):
            coarse,fine=results[0][sector][k],results[1][sector][k]
            error=abs(fine-coarse)/(1+abs(fine));convergence.append(c.nstr(error,10))
            assert error<c.mpf('1e-60'),(Z,sector,k,'reference convergence',error)
            row=got['m15_core_C0_C1_terms'][sector]['coefficient_C0_C1'][k]
            approximation=c.mpf(row['approximation'])
            assert abs(approximation-fine)/(1+abs(fine))<c.mpf('1e-60'),(Z,sector,k,'finite integral')
            lo,hi=_bounds(c,f.iv.mpf([c.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']]))
            # A 75-digit numerical reference is not an exact real value;
            # the directed production interval can be sharper than its
            # roundoff. Keep this comparison tolerance separate from the
            # original source/Volterra error, which is never enlarged.
            comparison_roundoff=c.mpf('1e-60')*(1+abs(fine))
            # The nominal feedback iterate need not lie inside the true
            # principal-origin interval. Its distance is already charged
            # in the published complete approximation error.
            declared_error=c.mpf(row['absolute_error_upper'])
            assert lo-declared_error-comparison_roundoff<=fine<=hi+declared_error+comparison_roundoff,(Z,sector,k,'source comparison')
            checks+=1
    return dict(direct_density_C0_C1_comparisons=checks,reference_Gauss_orders=[64,96],
        relative_reference_convergence=convergence,pressure_reference_uses_minus_rho_p_rho=True,
        numerical_reference_comparison_tolerance='1e-60*(1+abs(reference)); not added to production error')


def independent_error_majorants(f):
    """Integrate perturbed exact polynomials to exercise error integration."""
    c,iv=f.ctx,f.iv
    # Signed errors with independent Z slopes. Uniform source/feedback
    # contributions are tested at their finite rho powers, including endpoint.
    r=iv.mpf(4);U=[[iv.mpf(1),iv.mpf('.2')],[iv.mpf('.3'),iv.mpf('-.1')]]
    u=[[iv.mpf(0)]*2,[iv.mpf('-.4'),iv.mpf('.5')]]
    p=[[iv.mpf(0)]*2,[iv.mpf('-.7'),iv.mpf('.1')]]
    du=[[iv.mpf(0)]*2,[iv.mpf('.03'),iv.mpf('-.04')]]
    dp=[[iv.mpf(0)]*2,[iv.mpf('.05'),iv.mpf('-.06')]]
    dU=[[iv.mpf('.02'),iv.mpf('-.01')],[iv.mpf(0)]*2]
    add=lambda a,b:[[a[i][k]+b[i][k] for k in (0,1)] for i in range(len(a))]
    base=signed_prefix_rows(U,u,p,r,iv);changed=signed_prefix_rows(add(U,dU),add(u,du),add(p,dp),r,iv)
    # Independent absolute product/endpoint inequalities, without sampling.
    SU=[iv.mpf('2.2'),iv.mpf('.6')];Bu=[iv.mpf('.4'),iv.mpf('.5')]
    EU=[iv.mpf('.02'),iv.mpf('.01')];eu=[iv.mpf('.03'),iv.mpf('.04')];ep=[iv.mpf('.05'),iv.mpf('.06')]
    bound=[r*r*(SU[0]*eu[0]+Bu[0]*EU[0]+EU[0]*eu[0]),
           r*r*(SU[1]*eu[0]+Bu[0]*EU[1]+EU[1]*eu[0]+SU[0]*eu[1]+Bu[1]*EU[0]+EU[0]*eu[1])]
    count=0
    for k in (0,1):
        assert max(abs(t) for t in _bounds(c,changed[0][k]-base[0][k]))<=_bounds(c,bound[k])[1];count+=1
        assert max(abs(t) for t in _bounds(c,changed[1][k]-base[1][k]))<=_bounds(c,iv.mpf('1.5')*r*r*ep[k])[1];count+=1
    return count


def run():
    start=time.perf_counter();identity=original_density_identity()
    f=OriginalN1CoreM15Numeric(precision=75,quadrature_nodes=128)
    audit=mp.mp.clone();audit.dps=150
    error_checks=independent_error_majorants(f)
    print('Original signed density/physical weights and integrated error identities pass',flush=True)
    samples=[];direct=[];source=zero=signs=overlaps=0;relative=[]
    for Z in ('-.3','0','.3'):
        for degree,kernel in ((12,16),(16,24)):
            for rho in ('0','1','2','4'):
                p=f.evaluate(Z,rho,degree=degree,kernel_degree=kernel);samples.append(p)
                assert all(p['contract'].values());source+=len(p['contract'])
                assert not p['limits']['complete_signed_inner_m15'];source+=1
                for term in p['m15_core_C0_C1_terms']:
                    for row in term['coefficient_C0_C1']:
                        # Preserve all directed endpoint bits when checking
                        # a 75-digit decimal's much smaller roundoff error.
                        lo,hi=[audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]
                        a=audit.mpf(row['approximation']);e=audit.mpf(row['absolute_error_upper'])
                        assert max(abs(lo-a),abs(hi-a))<=e,(Z,degree,rho,term['name']);source+=1
                        if rho=='0':assert (lo,hi)==(0,0) and a==0;zero+=1
                if degree==16 and rho!='0':
                    row=p['m15_core_C0_C1_terms'][1]['coefficient_C0_C1'][0]
                    lo,hi=_bounds(f.ctx,f.iv.mpf([f.ctx.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]))
                    assert lo>0;signs+=1
                    rel=f.ctx.mpf(row['relative_error_upper']);assert rel<f.ctx.mpf('1e-18')
                    relative.append(f.ctx.nstr(rel,12))
                if degree==16:
                    coarse=f.evaluate(Z,rho,degree=12,kernel_degree=16)
                    for a,b in zip(coarse['m15_core_C0_C1_terms'],p['m15_core_C0_C1_terms']):
                        for x,y in zip(a['coefficient_C0_C1'],b['coefficient_C0_C1']):
                            pairs=[tuple(audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']) for row in (x,y)]
                            assert max(q[0] for q in pairs)<=min(q[1] for q in pairs);overlaps+=1
            print('Signed core m15 evaluated: Z='+Z+', degree='+str(degree),flush=True)
        direct.append(direct_polynomial_density(f,Z,16,24))
    rejected=0
    for args in [('1','1'),('nan','1'),('.3','-1'),('.3','4.1')]:
        try:f.evaluate(*args)
        except ValueError:rejected+=1
        else:raise AssertionError('Outside original core admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,original_density_and_physical_weights=identity,
            perturbed_polynomial_error_majorants=error_checks,direct_signed_density=direct,
            actual_core_C0_C1_packets=len(samples),source_error_and_scope_conditions=source,
            exact_axis_moment_zeros=zero,positive_pressure_moment_sectors=signs,
            pressure_sector_C0_relative_errors=relative,source_degree_refinement_enclosure_comparisons=overlaps,
            rejected_requests=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(signed_core_m15_C0_C1_approximation_and_complete_error=True,
            complete_signed_inner_m15=False,inner_cutoff_first_collar=False,lower_collar=False,
            active_compact_signed_pressure=False,independent_preRp_signed_numeric_seam=False,
            global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
            stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
