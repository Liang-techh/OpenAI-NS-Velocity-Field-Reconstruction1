"""Local nonlinear radial Taylor recursion for the exact Section 8.2 core.

Inputs are Taylor coefficients in Z of the actual axis data and pressure.
This radial recursion is not the time-scale recursion of later sections.
No global analytic-contraction or matching claim follows from a finite jet.
"""
import mpmath as mp


def pressure_taylor(pressure,Z,*,degree=8,spacing='.02',precision=160):
    """Local polynomial jets of the same axis pressure; explicit FD input.

Return coefficients, not a certification of the underlying pressure
derivatives. Call at another spacing to measure interpolation sensitivity.
"""
    with mp.workdps(precision):
        z=mp.mpf(str(Z)); h=mp.mpf(str(spacing)); count=degree+1
        if degree%2 or degree<2 or h<=0 or abs(z)+degree*h/2>=1:
            raise ValueError('Even degree and an interior positive-spacing stencil required')
        nodes=[mp.mpf(i-degree//2) for i in range(count)]
        matrix=mp.matrix([[t**k for k in range(count)] for t in nodes])
        values=mp.matrix([mp.mpf(pressure(z+h*t)) for t in nodes])
        coefficients=mp.lu_solve(matrix,values)
        return [coefficients[k]/h**k for k in range(count)]


def core_coefficients(Z,delta,*,F0_Z_taylor,U0_Z_taylor,P0_Z_taylor,
                      radial_degree=3,precision=160):
    """Return radial coefficients, each represented by a local Z polynomial.

The input coefficient k equals derivative_k/k!, not derivative_k.
At least radial_degree+2 initial Z coefficients are required.
"""
    with mp.workdps(precision):
        degree=int(radial_degree)
        if degree<1 or degree!=radial_degree:
            raise ValueError('Positive integer radial_degree required')
        K=min(len(F0_Z_taylor),len(U0_Z_taylor),len(P0_Z_taylor))
        initial_count=K
        if K<degree+2:raise ValueError('Insufficient axis Z Taylor coefficients')
        def pad(v):return list(map(mp.mpf,v))+[mp.mpf(0)]*(K-len(v))
        def const(v):return pad([v])
        def add(*args):return [sum(a[k] for a in args) for k in range(K)]
        def scale(a,c):return [v*c for v in a]
        def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(K)]
        def diff(a):return [(k+1)*a[k+1] if k+1<len(a) else mp.mpf(0) for k in range(K)]
        def inverse(a):
            q=[1/a[0]]
            for k in range(1,K):q.append(-sum(a[i]*q[k-i] for i in range(1,k+1))/a[0])
            return q
        z=pad([mp.mpf(str(Z)),1]); dt=mp.mpf(str(delta)); one=const(1)
        d=add(one,scale(mul(z,z),-1)); L=add(one,scale(mul(z,z),-dt))
        invL=inverse(L)
        f=[pad(F0_Z_taylor)]; u=[pad(U0_Z_taylor)]; pressure=[pad(P0_Z_taylor)]
        for n in range(degree):
            # One Z derivative is consumed by each radial recurrence.
            # Retain only jets that can still affect a trusted coefficient.
            K=initial_count-n-1
            W=[]; H=[]
            for i in range(n+1):
                wi=scale(add(scale(mul(z,u[i]),1-dt),mul(d,diff(u[i]))),-mp.mpf(1)/(i+1))
                if i==0:wi=add(one,wi)
                W.append(wi)
                hi=mul(d,u[i])
                if i==0:hi=add(hi,scale(z,(1-dt)/2))
                H.append(hi)
            uf=const(0); uu=const(0); rhsf=const(0); rhsu=const(0)
            for i in range(n+1):
                uf=add(uf,mul(u[i],f[n-i])); uu=add(uu,mul(u[i],u[n-i]))
                rhsf=add(rhsf,scale(mul(W[i],f[n-i]),n-i+1),mul(H[i],diff(f[n-i])))
                rhsu=add(rhsu,scale(mul(W[i],u[n-i]),n-i),mul(H[i],diff(u[n-i])))
            rhsf=add(rhsf,scale(add(f[n],scale(mul(z,uf),-2)),dt/2))
            rhsu=add(rhsu,scale(add(u[n],scale(mul(z,uu),-2)),(1+dt)/2),
                     mul(d,diff(pressure[n])),scale(mul(z,pressure[n]),-2*(1+dt)))
            if n:
                ff=const(0)
                for i in range(n):ff=add(ff,mul(f[i],f[n-1-i]))
                rhsu=add(rhsu,scale(mul(z,ff),-2))
            f.append(scale(mul(invL,rhsf),mp.mpf(1)/(2*(n+1)*(n+2))))
            u.append(scale(mul(invL,rhsu),mp.mpf(1)/(2*(n+1)**2)))
            ff=const(0)
            for i in range(n+1):ff=add(ff,mul(f[i],f[n-i]))
            pressure.append(scale(ff,mp.mpf(1)/(n+1)))
        return {'F':f,'Uz':u,'P':pressure,'Z':mp.mpf(str(Z)),
                'delta':dt,'radial_degree':degree,'initial_Z_degree':initial_count-1,'precision':precision,
                'scope':'Local exact-equation radial jets; no temporal recursion or global matching.'}


def evaluate_core_jets(coefficients,R):
    """Center-Z values and derivatives of one finite radial polynomial."""
    with mp.workdps(max(80,coefficients['precision'],mp.mp.dps)):
        r=mp.mpf(str(R)); result={}
        for name in ('F','Uz','P'):
            rows=coefficients[name]
            result[name]=sum(row[0]*r**n for n,row in enumerate(rows))
            result[name+'_Z']=sum(row[1]*r**n for n,row in enumerate(rows))
            result[name+'_R']=sum(n*row[0]*r**(n-1) for n,row in enumerate(rows) if n)
            result[name+'_RR']=sum(n*(n-1)*row[0]*r**(n-2) for n,row in enumerate(rows) if n>1)
        result['Mz_over_R']=sum(row[0]*r**n/(n+1) for n,row in enumerate(coefficients['Uz']))
        result['Mz_Z_over_R']=sum(row[1]*r**n/(n+1) for n,row in enumerate(coefficients['Uz']))
        return result


def core_equation_defects(coefficients,R):
    with mp.workdps(max(80,coefficients['precision'],mp.mp.dps)):
        a=evaluate_core_jets(coefficients,R); r=mp.mpf(str(R))
        z=coefficients['Z']; dt=coefficients['delta']; d=1-z*z; L=1-dt*z*z
        W=1-(1-dt)*z*a['Mz_over_R']-d*a['Mz_Z_over_R']
        H=(1-dt)*z/2+d*a['Uz']
        lhsf=2*L*(r*a['F_RR']+2*a['F_R'])
        rhsf=W*(a['F']+r*a['F_R'])+dt/2*(1-2*z*a['Uz'])*a['F']+H*a['F_Z']
        lhsu=2*L*(r*a['Uz_RR']+a['Uz_R'])
        rhsu=W*r*a['Uz_R']+(1+dt)/2*(1-2*z*a['Uz'])*a['Uz']+H*a['Uz_Z']
        rhsu+=d*a['P_Z']-2*(1+dt)*z*a['P']-2*z*r*a['F']**2
        return {'angular':lhsf-rhsf,'axial':lhsu-rhsu,
                'pressure':a['P_R']-a['F']**2}
