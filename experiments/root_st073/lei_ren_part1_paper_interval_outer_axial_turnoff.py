"""Parametric O.2 axial turnoff with inherited physical C1 moments.

No accepted Md=0.5 pressure is silently relabeled as a new schedule. This
adapter takes an explicit inlet packet and its pressure-schedule Md. A new
production datum must be regenerated before using a different Md with it.
"""
from fractions import Fraction
from lei_ren_part1_paper_interval_outer_slope_field import fraction_box
from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative
from lei_ren_part1_paper_interval_repaired_reference_field import stress_values
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def cutoff_integrals(c,Md,phase,cells=256):
    """Full directed B and B^2 masses through y=exp(Md*phase).

    K_j=int_1^y exp(s)*B(log(s)/Md)^j ds. Substitute s=exp(Md*q)
    and integrate over q in [0,phase]; positive interval rectangles include
    every point, including flat endpoints. No omitted-tail approximation.
    """
    md=Fraction(Md);q=Fraction(phase)
    if md<=0 or not 0<=q<=1 or not isinstance(cells,int) or cells<1:
        raise ValueError('positive Md, phase in [0,1], positive integer cells required')
    mdv=fraction_box(c,md);dq=fraction_box(c,q/cells)
    mass=c.mpf(0);square=c.mpf(0)
    for i in range(cells):
        left=fraction_box(c,q*i/cells);right=fraction_box(c,q*(i+1)/cells)
        sigleft=sigma_value_derivative(c,left)[0]
        sigright=sigma_value_derivative(c,right)[0]
        cutoff=c.mpf([endpoints(1-sigright)[0],endpoints(1-sigleft)[1]])
        s=c.exp(mdv*c.mpf([endpoints(left)[0],endpoints(right)[1]]))
        weight=mdv*s*c.exp(s)*dq
        mass+=weight*cutoff;square+=weight*cutoff**2
    return dict(B_mass=mass,B_squared_mass=square)


class IntervalOuterAxialTurnoff:
    def __init__(self,c,z,delta,inlet,Md,pressure_schedule_Md):
        md=Fraction(Md)
        if md<=1:raise ValueError('paper O.2 requires Md>1')
        if md!=Fraction(pressure_schedule_Md):
            raise ValueError('regenerate pressure datum for the chosen Md before attaching this layer')
        recorded=inlet.get('pressure_schedule_Md')
        if recorded is None or Fraction(recorded)!=md:
            raise ValueError('inlet packet must record the same pressure-schedule Md')
        self.c=c;self.z=z;self.delta=delta;self.inlet=inlet;self.Md=md
        # C1 radial input must be the end of the angular slope transition.
        if inlet.get('y')!='1' or inlet.get('stage')!='O.2 angular slope transition':
            raise ValueError('O.2 slope endpoint packet required')
        self.R1=inlet['R'];self.u1=inlet['Utheta'];self.P0=inlet['P0']
        self.m0=inlet['physical_moments']

    def evaluate_phase(self,phase,cells=256):
        """phase=log(y)/Md, y=log(R/Rref), 0<=phase<=1."""
        q=Fraction(phase)
        if not 0<=q<=1:raise ValueError('phase in [0,1] required')
        c=self.c;md=fraction_box(c,self.Md);qq=fraction_box(c,q)
        y=c.exp(md*qq);t=y-1
        sig,ds=sigma_value_derivative(c,qq)
        B=1-sig;By=-ds/(md*y)
        K=cutoff_integrals(c,self.Md,q,cells)
        return self._packet(y,t,B,By,K,q)

    def evaluate_buffer(self,offset,cells=256):
        """Zero-Uz buffer after the cutoff, y=exp(Md)+offset, 0<=offset<=11.

        offset=11 is R_d, where the following -1/2-mu turn starts.
        The retained Mz and mixed moment remain nonzero.
        """
        q=Fraction(offset)
        if not 0<=q<=11:raise ValueError('buffer offset in [0,11] required')
        c=self.c;y=c.exp(fraction_box(c,self.Md))+fraction_box(c,q)
        K=cutoff_integrals(c,self.Md,1,cells)
        return self._packet(y,y-1,c.mpf(0),c.mpf(0),K,None,str(q))

    def _packet(self,y,t,B,By,K,phase,buffer_offset=None):
        c=self.c;z=self.z;R=self.R1*c.exp(t)
        u=self.u1*c.exp(-t/2);uy=-u/2;V=z*(4*B);Vy=z*(4*By)
        dtheta=self.u1*(c.sqrt(2)*self.R1**c.mpf('1.5')*(c.exp(t)-1))
        radial_mass=self.R1*c.exp(-1)*K['B_mass']
        radial_square=self.R1*c.exp(-1)*K['B_squared_mass']
        mixed_scale=self.u1*(c.sqrt(2)*self.R1**c.mpf('1.5')*c.exp(-1)*K['B_mass'])
        swirl=self.u1*self.u1*(self.R1*t/2)
        pressure=self.u1*self.u1*((1-c.exp(-t))/2)
        m=dict(z=self.m0['z']+z*(4*radial_mass),
            theta=self.m0['theta']+dtheta,
            theta_z=self.m0['theta_z']+z*mixed_scale*4,
            z_theta=self.m0['z_theta']+z*z*(16*radial_square)-swirl,
            p=self.m0['p']+pressure)
        P=self.P0+m['p'];root=c.sqrt(2*R);F=u/root
        st=c.mpf(-2);sz=2*Vy[0]/u[0]
        stresses=stress_values(c,z,self.delta,R,u,uy,V,Vy,m,P)
        Itheta=stresses['stress']['I_theta'];Iz=stresses['stress']['I_z']
        tt=Itheta/F[0]+st;tz=Iz/F[0]+sz
        stresses['stress'].update(S_theta=F[0]*st,S_z=F[0]*sz,
            T_theta=Itheta+F[0]*st,T_z=Iz+F[0]*sz)
        stresses['normalized_stress']=dict(S_theta_over_F=st,S_z_over_F=sz,
            T_theta_over_F=tt,T_z_over_F=tz)
        stresses['cone']=cone(c,st,sz,tt,tz)
        zv=z[0];L=1-self.delta*zv*zv;mz=m['z']
        Ur=(2*zv*R*V[0]-(1-self.delta)*zv*mz[0]-(1-zv*zv)*mz[1])/(L*root)
        Ur_R=((1+self.delta)*zv*V[0]+2*zv*Vy[0]-(1-zv*zv)*V[1])/(L*root)-Ur/(2*R)
        return dict(y=y,phase=str(phase) if phase is not None else None,
            buffer_offset=buffer_offset,stage='O.2 axial turnoff / retained-moment buffer',
            R=R,F=F,Utheta=u,Utheta_y=uy,Uz=V,Uz_y=Vy,P=P,P0=self.P0,
            physical_moments=m,cutoff=B,cutoff_y=By,cutoff_integrals=K,
            Ur_value=Ur,Ur_R_value=Ur_R,
            moment_radial_rhs=dict(z=V,theta=u*root,theta_z=u*V*root,
                z_theta=V*V-u*u/2,p=u*u/(2*R)),
            all_five_moments_transported=True,axis_pressure_datum_unchanged=True,
            Md=str(self.Md),declared_pressure_schedule_Md_matches=True,
            inlet_pressure_is_analytic_preheat_datum_certified=False,
            whole_outer_cone_certified=False,Ur_Z_available=False,
            exact_heat_exterior_installed=False,temporal_recursion=False,**stresses)
