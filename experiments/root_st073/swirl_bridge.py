"""C2 radial swirl bridge; kinematic component only, not a solved NS transition."""
import sys,json,numpy as np
from radial_continuation import ROOT,FullRadialField
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from heat_exterior import heat_factor

def quintic(left,right,width):
    a=np.zeros(6);a[:3]=[left[0],width*left[1],width**2*left[2]/2]
    a[3:]=np.linalg.solve([[1,1,1],[3,4,5],[6,12,20]],np.array([right[0]-sum(a[:3]),width*right[1]-a[1]-2*a[2],width**2*right[2]-2*a[2]]))
    return a

def traces(f,eta,tau,c,X=3/64,outer_ratio=2.):
    q=tau/(1-eta**2);rho=np.sqrt(2*q*X);ri=np.sqrt(f.nu)*rho;ro=outer_ratio*ri
    b=f.coefficients(float(eta),float(q))[1,:,0];B=np.polynomial.polynomial.polyval(X,b);Bs=np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(b))/(2*q);Bss=np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(b,2))/(2*q)**2
    left=np.array([np.sqrt(f.nu)*rho*B,B+2*rho**2*Bs,(6*rho*Bs+4*rho**3*Bss)/np.sqrt(f.nu)],float)
    rr=ro/np.sqrt(f.nu);s=rr**2/2;A=.5+f.h;Z=2*tau/s
    H,H1,H2=[float(heat_factor(Z,f.h,j)) for j in range(3)]
    K=c*s**(-A)*H;Ks=c*s**(-A-1)*(-A*H-Z*H1);Kss=c*s**(-A-2)*(A*(A+1)*H+2*(A+1)*Z*H1+Z**2*H2)
    right=np.array([np.sqrt(f.nu)*K,rr*Ks,(Ks+rr**2*Kss)/np.sqrt(f.nu)])
    return ri,ro,left,right

def run():
    f=FullRadialField.load(ROOT/'radial_continuation/candidate.json');c=json.loads((ROOT/'heat_join/screen.json').read_text())['heat_amplitude'];rows=[]
    for k in (0,3,6):
        for eta in (-.4,0,.4):
            ri,ro,left,right=traces(f,eta,.5*2**(-k),c);w=ro-ri;a=quintic(left,right,w)
            errors=[]
            for j in (0,1,2):
                coeff=np.polynomial.polynomial.polyder(a,j)/w**j
                errors.append([float(np.polynomial.polynomial.polyval(y,coeff)-target[j]) for y,target in [(0,left),(1,right)]])
            values=np.polynomial.polynomial.polyval(np.linspace(0,1,101),a)
            rows.append(dict(k=k,eta=eta,inner_radius=ri,outer_radius=ro,coefficients_in_unit_radial_coordinate=a.tolist(),endpoint_derivative_errors=errors,swirl_min=float(values.min()),swirl_max=float(values.max())))
    report=dict(rows=rows,heat_amplitude=c,scope='Autonomous quintic interpolation between corrected inner swirl and paper-inspired heat exterior. C2 in radius, matches tangential viscous traction. Poloidal transport, pressure, time evolution residual, axial closure and stress realization NOT solved. No PDE acceptance.',pde_validated=False,global_field_ready=False)
    out=ROOT/'swirl_bridge';out.mkdir(exist_ok=True);(out/'coefficients.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print('max endpoint errors by derivative',np.max(np.abs([r['endpoint_derivative_errors'] for r in rows]),axis=(0,2)).tolist());print('min swirl',min(r['swirl_min'] for r in rows))
if __name__=='__main__':run()
