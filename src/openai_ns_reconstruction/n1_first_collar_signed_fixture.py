"""Independent symbolic and direct-integral fixture; synthetic inputs only."""
from __future__ import annotations


def run(matrix, green, c):
    import sympy as S
    import mpmath as mp
    z=S.symbols('z'); r=S.Rational(9,4); x=S.Rational(3,2); de=S.Rational(1,7); ep=S.Rational(1,13); z0=S.Rational(1,5)
    F=2+z+2*z*z+z**3;Fr=1-2*z+z*z;FZ=S.diff(F,z)
    U=3-z+z*z-z**3;Ur=2+z-2*z*z;UZ=S.diff(U,z);Q=1+z+3*z*z
    d=1-z*z;L=1-de*z*z;H=z*(1-de)/2+d*U;E=z*U-S.Rational(1,2);Cf=F+r*Fr
    Zf=((-2-de)*z*F+d*FZ-2*z*r*Fr)/L
    Zw=((-1-de)*z*U+d*UZ-2*z*r*Ur)/L
    B=Q+(1-2*z*U)/L
    expected={'C0':{'4,1':4*ep*x*F,'5,1':2*ep*(Q+E*(-2+de)/L),
       '5,2':ep*(2*Zf-2*(-1+de)*z*Cf/L),'5,3':-2*ep*(1+de)*z*Cf/L,
       '5,5':ep*x*B,'6,1':-8*ep**2*z*r*F/L,
       '6,2':2*ep*(E*(-1+de)/L+Zw-Ur*z*(-1+de)*r/L),
       '6,3':-2*ep*(1+de)*z*r*Ur/L,'6,4':-4*ep*z/L,'6,6':ep*x*B},
     'B1':{'5,1':2*ep*H/L,'5,2':-2*ep*d*Cf/L,'5,3':-2*ep*d*Cf/L,
       '6,2':2*ep*(H-Ur*d*r)/L,'6,3':-2*ep*d*r*Ur/L,'6,4':2*ep*d/L}}
    def cv(q):
        q=S.Rational(q);return c.mpf(int(q.p))/int(q.q)
    def jet(q):return tuple(cv(S.diff(q,z,k).subs(z,z0)) for k in range(3))
    actual=matrix.evaluate_raw_matrices(cv(r),cv(x),cv(de),cv(ep),cv(z0),*(jet(q) for q in (F,Fr,FZ,U,Ur,UZ,Q)))
    checks={}
    for name,entries in expected.items():
        for key,expr in entries.items():
            for k in range(3):
                ref=cv(S.diff(expr,z,k).subs(z,z0));got=actual[name][key][k]
                checks['matrix_'+name+'_'+key+'_'+str(k)]=bool(got.a<=ref.a and got.b>=ref.b)
    hb=c.mpf(1)/4;sc=c.mpf(1)/2048;s=sc/4
    core0={j:tuple(c.mpf((j+1)*(i+1))/7 for i in range(6)) for j in range(4)}
    core1=tuple(c.mpf(i-2)/9 for i in range(6))
    forcing={j:{str(row):c.mpf((j+1)*(row-5))/11 for row in (4,5,6)} for j in range(4)}
    const={name:{k:{key:c.mpf(((i+2)*(k+1))*(-1 if i%2 else 1))/37
         for i,key in enumerate(sorted(keys))} for k in range(3)}
         for name,keys in (('C0',green.C0_KEYS),('B1',green.B1_KEYS))}
    result=green.enclose(c,hb=hb,sc=sc,phase=s,core_W0_Z=core0,core_W1_Z2=core1,
         forcing_Z=lambda v:forcing,matrix_Z=lambda v:const,inner_cells=4,outer_cells=5,source_family='synthetic')
    inlet=green.enclose(c,hb=hb,sc=sc,phase=c.mpf(0),core_W0_Z=core0,core_W1_Z2=core1,
         forcing_Z=lambda v:forcing,matrix_Z=lambda v:const,inner_cells=1,outer_cells=1,source_family='synthetic')
    mp.mp.dps=80
    mh=mp.mpf(1)/4;ms=mp.mpf(1)/8192
    def transport(v,a):
        e=mp.exp(-mh*v);f,u,k,p,fx,ux=a
        return (f+(1-e)*fx,u+mh*v*ux,e*k-(1-e)*ux,p,mp.exp(-3*mh*v/2)*fx,mp.exp(-mh*v/2)*ux)
    def word(v,j):
        X=2*mp.exp(mh*v/2)
        A=(X**2-4)/4-(X**4-16)/(8*X**2)
        B=X**2/4-2*mp.log(X/2)-1
        mass=(A,B,-A,X-2,(X**4-16)/(4*X**3),(X**2-4)/(2*X))
        g={row:mp.mpf((j+1)*(int(row)-5))/11 for row in ('4','5','6')}
        initial=transport(v,[mp.mpf((j+1)*(i+1))/7 for i in range(6)])
        return [initial[i]+mass[i]*g[row] for i,row in enumerate(green.ROWS)]
    specs=(('C0',2,0,1),('C0',1,1,2),('C0',0,2,1),('B1',2,1,1),('B1',1,2,2),('B1',0,3,1))
    def source(v):
        w={j:word(v,j) for j in range(4)};out={str(k):mp.mpf(0) for k in (4,5,6)}
        for name,mo,wo,mul in specs:
            keys=green.C0_KEYS if name=='C0' else green.B1_KEYS
            for i,key in enumerate(sorted(keys)):
                row,col=key.split(',');val=mp.mpf((i+2)*(mo+1)*(-1 if i%2 else 1))/37
                out[row]+=mul*val*w[wo][int(col)-1]
        return out
    def integrand(v,i):
        u=2*mp.exp(mh*v/2);dist=mh*(ms-v);ratio=mp.exp(-dist/2)
        K=u*(1-mp.exp(-dist))/2;ks=(K,u*dist/2,-K,mp.mpf(1),ratio**3,ratio)
        return source(v)[green.ROWS[i]]*ks[i]*mh*u/2
    initial=transport(ms,[mp.mpf(i-2)/9 for i in range(6)])
    for i,name in enumerate(green.FIELDS):
        value=initial[i]+mp.quad(lambda v:integrand(v,i),[0,ms])
        got=result['finite_W1_Z2'][name];reference=c.mpf(mp.nstr(value,70))
        checks['signed_green_'+name]=bool(got.a<=reference.a and got.b>=reference.b)
        ref=c.mpf(i-2)/9;at=inlet['finite_W1_Z2'][name]
        checks['inlet_'+name]=bool(at.a<=ref.a and at.b>=ref.b)
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,
       'matrix_reference':'independent symbolic ordinary Z differentiation',
       'green_reference':'independent high-precision direct integral with analytic W0 source masses',
       'synthetic_only':True}
