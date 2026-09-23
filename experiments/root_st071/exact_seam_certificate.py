"""Exact rational Bernstein signs for the STORED finite polynomial seam only.
IEEE-754 coefficients are interpreted as their exact binary rational values.
This certifies neither a solved PDE nor an admissible full stress cone.
"""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import math,json,hashlib,time,argparse
import numpy as np
ROOT=Path(__file__).resolve().parent

def add(a,b):
    out=[Q(0)]*max(len(a),len(b))
    for i,x in enumerate(a):out[i]+=x
    for i,x in enumerate(b):out[i]+=x
    return out

def mul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def times(a,c):return [x*c for x in a]
def cheb_to_power(co):
    ts=[[Q(1)],[Q(0),Q(1)]]
    for k in range(2,len(co)):ts.append(add([Q(0)]+times(ts[-1],2),times(ts[-2],-1)))
    out=[Q(0)]*len(co)
    for x,t in zip(co,ts):out=add(out,times(t,x))
    return out

def shifted_power(p,a,b):
    n=len(p)-1;step=b-a;qa=[Q(0)]*(n+1)
    ap=[a**j for j in range(n+1)];sp=[step**j for j in range(n+1)]
    for i in range(n+1):qa[i]=sp[i]*sum((p[j]*math.comb(j,i)*ap[j-i] for j in range(i,n+1)),Q(0))
    return qa

def bernstein(p,a,b):
    power=shifted_power(p,a,b);n=len(power)-1
    return [sum((power[i]*Q(math.comb(k,i),math.comb(n,i)) for i in range(k+1)),Q(0)) for k in range(n+1)]

def evaluate(p,x):
    out=Q(0)
    for c in p[::-1]:out=out*x+c
    return out

def seam_polynomials(F,U):
    # T_j(1)=1 and T_j'(1)=j^2; d/dX=8 d/dxi.
    fc=[sum((Q(float(F[i,j])) for i in range(F.shape[0])),Q(0)) for j in range(F.shape[1])]
    fx=[sum((8*i*i*Q(float(F[i,j])) for i in range(F.shape[0])),Q(0)) for j in range(F.shape[1])]
    ux=[sum((8*i*i*Q(float(U[i,j])) for i in range(U.shape[0])),Q(0)) for j in range(U.shape[1])]
    f,dx,du=[cheb_to_power(c) for c in [fc,fx,ux]]
    # G=F^2 a (v-2), with a=-Fx/(2F), b=Ux/(sqrt(2)F).
    G=add(add(times(mul(dx,dx),Q(1,4)),times(mul(du,du),Q(1,2))),mul(dx,f))
    A=add(times(dx,Q(-1,2)),times(f,Q(-1,20)))
    return dict(F=f,A=A,G=G),dict(F=f,Fx=dx,Ux=du)

def pack(q):return str(q.numerator)+'/'+str(q.denominator)
def run(model,subdivisions=48):
    start=time.monotonic();model=Path(model)
    with np.load(model,allow_pickle=False) as d:polys,jets=seam_polynomials(d['F_cheb'],d['U_cheb'])
    cert=dict(schema='exact_rational_bernstein_seam_v1',model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),coefficient_semantics='Exact binary rational values of stored IEEE coefficients, no coefficient uncertainty',variable='s=2eta in [-1,1], seam X=1/4',source_necessary_condition='F>0, a>1/20, v>2',subdivisions=subdivisions,polynomials={},scope='Only three sign inequalities for this stored polynomial trace. No leading equation, full stress, annulus, PDE, global energy or scale-recursion certification.')
    for name,p in polys.items():
        rows=[]
        for k in range(subdivisions):
            a=Q(-1)+Q(2*k,subdivisions);b=Q(-1)+Q(2*(k+1),subdivisions);bs=bernstein(p,a,b)
            rows.append(dict(left=pack(a),right=pack(b),coefficients=[pack(q) for q in bs]))
        vals=[Q(q) for row in rows for q in row['coefficients']];low=min(vals)
        cert['polynomials'][name]=dict(degree=len(p)-1,power_coefficients=[pack(q) for q in p],segments=rows,min_coefficient=pack(low),min_coefficient_float=float(low),all_strictly_positive=all(q>0 for q in vals))
        print('EXACT',name,float(low),all(q>0 for q in vals),flush=True)
    cert['sign_certificate_pass']=all(x['all_strictly_positive'] for x in cert['polynomials'].values());cert['elapsed']=time.monotonic()-start;return cert

def witness(model,eta=Q(3429,10000)):
    with np.load(model,allow_pickle=False) as d:ps,jets=seam_polynomials(d['F_cheb'],d['U_cheb'])
    v={k:evaluate(p,2*eta) for k,p in ps.items()}
    return dict(model_sha256=hashlib.sha256(Path(model).read_bytes()).hexdigest(),eta=pack(eta),values={k:pack(q) for k,q in v.items()},values_float={k:float(q) for k,q in v.items()},necessary_condition_fails=bool(v['F']>0 and v['A']>0 and v['G']<0),scope='An exact rational point witness, not an all-domain PDE statement')

def main():
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,default=ROOT/'data/ST071-FP.npz');p.add_argument('--out',type=Path,default=ROOT/'evidence/exact_seam_certificate.json');a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    c=run(a.model);a.out.write_text(json.dumps(c,indent=2)+'\n')
    w=witness(ROOT/'data/ST071-RN.npz');(ROOT/'evidence/exact_RN_counterexample.json').write_text(json.dumps(w,indent=2)+'\n')
if __name__=='__main__':main()
