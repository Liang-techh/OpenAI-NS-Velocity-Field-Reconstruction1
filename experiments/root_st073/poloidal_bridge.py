"""C3 streamfunction bridge from corrected ST073 core to zero poloidal exterior."""
import json,math,numpy as np
from radial_continuation import ROOT,FullRadialField

def septic(left,width):
    a=np.zeros(8)
    for j in range(4):a[j]=left[j]*width**j/math.factorial(j)
    matrix=np.array([[math.factorial(i)/math.factorial(i-j) for i in range(4,8)] for j in range(4)],float)
    rhs=[-np.polynomial.polynomial.polyval(1,np.polynomial.polynomial.polyder(a[:4],j)) for j in range(4)]
    a[4:]=np.linalg.solve(matrix,rhs)
    return a

def coefficients(f,eta,tau):
    X=3/64;q=tau/(1-eta*eta);rho=np.sqrt(2*q*X);r=np.sqrt(f.nu)*rho
    c=f.coefficients(float(eta),float(q))[2,:,0]
    C=np.polynomial.polynomial.polyval(X,c);Cs=np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(c))/(2*q);Css=np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(c,2))/(2*q)**2
    psi=f.nu**1.5*q*np.polynomial.polynomial.polyval(X,np.r_[0,c/np.arange(1,len(c)+1)])
    left=np.array([psi,r*np.sqrt(f.nu)*C,np.sqrt(f.nu)*C+r*2*rho*Cs,6*rho*Cs+4*rho**3*Css],float)
    return r,2*r,left,septic(left,r)

def run():
    f=FullRadialField.load(ROOT/'radial_continuation/candidate.json');rows=[]
    for k in (0,3,6):
        for eta in (-.4,0,.4):
            ri,ro,left,a=coefficients(f,eta,.5*2**(-k));w=ro-ri
            errors=[]
            for j in range(4):
                d=np.polynomial.polynomial.polyder(a,j)/w**j
                errors.append([float(np.polynomial.polynomial.polyval(0,d)-left[j]),float(np.polynomial.polynomial.polyval(1,d))])
            y=np.linspace(0,1,201);r=ri+w*y;uz=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(a))/w/r
            rows.append(dict(k=k,eta=eta,inner_radius=ri,outer_radius=ro,streamfunction_coefficients=a.tolist(),endpoint_jet_errors=errors,inner_axial_flux=float(2*np.pi*left[0]),annulus_axial_flux=float(-2*np.pi*left[0]),axial_velocity_min=float(uz.min()),axial_velocity_max=float(uz.max())))
    out=ROOT/'poloidal_bridge';out.mkdir(exist_ok=True)
    report=dict(rows=rows,definition='u_r=-psi_z/r, u_z=psi_r/r. Septic matches inner psi radial jets through order3 and zero outer jets. Derivatives must include z/t motion of interpolation endpoints.',scope='Autonomous kinematic streamfunction construction, not paper moment restoration or dynamic correction. Forces zero. No full momentum acceptance or axial global closure. Zero outer psi forces annular return flux exactly opposite the core cross-section flux.',pde_validated=False,global_field_ready=False)
    (out/'coefficients.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows[-1],indent=2))
if __name__=='__main__':run()
