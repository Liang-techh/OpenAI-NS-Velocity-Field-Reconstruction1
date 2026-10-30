"""Synthetic symbolic differentiation of the low-order core matrix."""
def run(leaf,c):
    import sympy as s
    z=s.symbols('z');z0=s.Rational(1,5);rho=s.Integer(1);delta=s.Rational(1,3);lam=s.Rational(7,4);j=s.Rational(5,3)
    D=1-z*z;L=1-delta*z*z;H0=z*(1-delta)/2+D*(4*z+j);ell=-lam*L*H0/(H0**2+(j/500)**2)
    ell0=ell.subs(z,z0);ell1=s.diff(ell,z).subs(z,z0);h=z-z0;relative=1+ell0*h+(ell0**2+ell1)*h**2/2
    phi=(1+2*z+3*z*z,2-z+z*z,1+z);u=(2+z+2*z*z,3-z,1-z);mean=(4+2*z+z*z,2+z,3)
    def grid(polys):
        return {f'rho{i}_Z{k}':exact(s.diff(p,z,k).subs(z,z0)) for i,p in enumerate(polys) for k in range(6-i)}
    def exact(x):
        x=s.cancel(x);return c.mpf(int(s.numer(x)))/int(s.denom(x))
    x=s.Integer(1);f=phi[0]*relative;fr=phi[1]*relative;fz=s.diff(f,z);w=u[0];wr=u[1];wz=s.diff(w,z);m=mean[0];mz=s.diff(m,z)
    q=(2*z*w-(1-delta)*z*m-D*mz)/L;H=z*(1-delta)/2+D*w;B=q+(1-2*z*w)/L;E=z*w-s.Rational(1,2)
    Cf=f+rho*fr;Zf=((-2-delta)*z*f+D*fz-2*z*rho*fr)/L;Zw=((-1-delta)*z*w+D*wz-2*z*rho*wr)/L
    c0={'4,1':(1,-1,4*x*f),'5,1':(0,-1,2*(q+E*(-2+delta)/L)),'5,2':(1,-1,2*Zf-2*(-1+delta)*z*Cf/L),'5,3':(1,-1,-2*(1+delta)*z*Cf/L),'5,5':(0,-1,B*x),'6,1':(1,-2,-8*z*rho*f/L),'6,2':(0,-1,2*(E*(-1+delta)/L+Zw-wr*z*(-1+delta)*rho/L)),'6,3':(0,-1,-2*(1+delta)*z*rho*wr/L),'6,4':(0,-1,-4*z/L),'6,6':(0,-1,B*x)}
    b1={'5,1':(0,-1,2*H/L),'5,2':(1,-1,-2*D*Cf/L),'5,3':(1,-1,-2*D*Cf/L),'6,2':(0,-1,2*(H-wr*D*rho)/L),'6,3':(0,-1,-2*D*rho*wr/L),'6,4':(0,-1,2*D/L)}
    got=leaf.build(c,exact(rho),exact(z0),exact(delta),exact(lam),exact(j),grid(phi),grid(u),grid(mean));checks={}
    for name,expected in (('C0',c0),('B1',b1)):
        for order in (0,1):
            checks[f'{name}_{order}_support']=set(got[name][str(order)])==set(expected)
            for key,(fp,lp,expr) in expected.items():
                term=got[name][str(order)][key][0];ref=exact(s.diff(expr,z,order).subs(z,z0));actual=term['coefficient']
                checks[f'{name}_{order}_{key}']=term['F0_power']==fp and term['Lambda_power']==lp and actual.a<=ref.a and actual.b>=ref.b
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,'reference':'symbolic ordinary derivative of independent rational polynomial source and anchored F0 second jet'}
