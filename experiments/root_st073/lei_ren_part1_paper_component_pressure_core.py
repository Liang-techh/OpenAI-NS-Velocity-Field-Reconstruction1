"""Finite radial core recursion retaining the complete preheat tail as a formal input.

P0=Pprefix+lambda*Ppost is propagated polynomially through the same nonlinear
radial recurrence. lambda=1 is the requested datum. Different lambda powers
are never added during construction; no pressure-perturbation order is cut.
This is a finite radial/Z jet, not an analytic or global field certificate.
"""
from functools import lru_cache
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_core_recursion import core_coefficients, evaluate_core_jets
from lei_ren_part1_paper_core_ra_experiment import _axis_F0_taylor, _axis_u0_taylor


class PressurePolynomial:
    """Polynomial in a formal pressure-tail parameter with MP coefficients."""
    def __init__(self, value=0):
        if isinstance(value, PressurePolynomial):
            self.atoms=dict(value.atoms)
        elif isinstance(value, dict):
            self.atoms={int(k):mp.mpf(v) for k,v in value.items() if v!=0}
            if any(k<0 for k in self.atoms):raise ValueError('Nonnegative powers required')
        else:
            v=mp.mpf(value)
            self.atoms={0:v} if v else {}

    def __add__(self, other):
        other=PressurePolynomial(other); atoms=dict(self.atoms)
        for k,v in other.atoms.items():atoms[k]=atoms.get(k,mp.mpf(0))+v
        return PressurePolynomial(atoms)
    __radd__=__add__
    def __neg__(self):return PressurePolynomial({k:-v for k,v in self.atoms.items()})
    def __sub__(self, other):return self+-PressurePolynomial(other)
    def __rsub__(self, other):return PressurePolynomial(other)+-self
    def __mul__(self, other):
        other=PressurePolynomial(other); atoms={}
        for k,v in self.atoms.items():
            for j,w in other.atoms.items():atoms[k+j]=atoms.get(k+j,mp.mpf(0))+v*w
        return PressurePolynomial(atoms)
    __rmul__=__mul__
    def __truediv__(self, other):
        other=PressurePolynomial(other)
        if any(k!=0 for k in other.atoms):raise ValueError('Only tail-independent divisors are supported')
        v=other.atoms.get(0,mp.mpf(0))
        if not v:raise ZeroDivisionError()
        return PressurePolynomial({k:a/v for k,a in self.atoms.items()})
    def __rtruediv__(self, other):return PressurePolynomial(other)/self
    def __pow__(self, exponent):
        if int(exponent)!=exponent or exponent<0:raise ValueError('Nonnegative integer power required')
        result=PressurePolynomial(1); base=self; n=int(exponent)
        while n:
            if n%2:result=result*base
            base=base*base;n//=2
        return result
    def evaluate(self, parameter=1):
        return mp.fsum(v*mp.mpf(parameter)**k for k,v in self.atoms.items())
    def component(self, power):return self.atoms.get(int(power),mp.mpf(0))


def build_component_coefficients(axis, pressure_datum, center, radial_degree=18):
    """Propagate all powers of the supplied complete post-Rv datum."""
    precision=max(axis.precision,pressure_datum.precision)
    with mp.workdps(precision):
        degree=int(radial_degree);length=degree+2
        z=mp.mpf(str(center))
        jets=pressure_datum.taylor_components(degree=length-1,center=z)
        physical=mp.exp(2*pressure_datum.log_pstar)
        prefix=[physical*v for v in jets['dominant_pressure_coefficients']]
        post=[physical*v for v in jets['post_Rv_pressure_coefficients']]
        p0=[PressurePolynomial({0:a,1:b}) for a,b in zip(prefix,post)]
        result=core_coefficients(z,axis.delta,
            F0_Z_taylor=_axis_F0_taylor(axis,z,length),
            U0_Z_taylor=_axis_u0_taylor(z,j=axis.j,length=length),P0_Z_taylor=p0,
            radial_degree=degree,precision=precision,scalar_converter=PressurePolynomial)
        result.update(pressure_datum_components=jets,
            complete_preheat_input_included=True,pressure_parameter_order_truncated=False,
            component_arithmetic_error_enclosed=False,post_stage_aggregation_error_enclosed=False,
            pressure_datum_quadrature_enclosed=False,global_field_installed=False)
        return result


def evaluate_component_core_coefficients(coefficients,r,z,delta):
    """Shared core/moment equations for scalar or component-valued radii."""
    jets=evaluate_core_jets(coefficients,r,radial_converter=lambda value:value)
    root=(2*r).sqrt() if hasattr(r,"sqrt") else mp.sqrt(2*r)
    f=coefficients['F'];u=coefficients['Uz']
    moments={k:coefficients['F'][0][0]*0 for k in ('theta','z','theta_z','z_theta','p')}
    moments_Z=dict(moments)
    for n,row in enumerate(f):
        moments['theta']+=2*row[0]*r**(n+2)/(n+2)
        moments_Z['theta']+=2*row[1]*r**(n+2)/(n+2)
    for n,row in enumerate(u):
        moments['z']+=row[0]*r**(n+1)/(n+1)
        moments_Z['z']+=row[1]*r**(n+1)/(n+1)
    for i,fi in enumerate(f):
        for j,uj in enumerate(u):
            n=i+j+2
            moments['theta_z']+=2*fi[0]*uj[0]*r**n/n
            moments_Z['theta_z']+=2*(fi[1]*uj[0]+fi[0]*uj[1])*r**n/n
        for j,fj in enumerate(f):
            n=i+j+1;value=fi[0]*fj[0];tangent=fi[1]*fj[0]+fi[0]*fj[1]
            moments['p']+=value*r**n/n
            moments_Z['p']+=tangent*r**n/n
            moments['z_theta']-=value*r**(n+1)/(n+1)
            moments_Z['z_theta']-=tangent*r**(n+1)/(n+1)
    for i,ui in enumerate(u):
        for j,uj in enumerate(u):
            n=i+j+1
            moments['z_theta']+=ui[0]*uj[0]*r**n/n
            moments_Z['z_theta']+=(ui[1]*uj[0]+ui[0]*uj[1])*r**n/n
    # Restore the exact F-polynomial integral, not a truncated pressure jet.
    jets['P']=coefficients['P'][0][0]+moments['p']
    jets['P_Z']=coefficients['P'][0][1]+moments_Z['p']
    jets['P_R']=jets['F']**2
    jets['Utheta']=root*jets['F']
    jets['Ur']=((2*z*r*jets['Uz']-(1-delta)*z*moments['z']
        -(1-z*z)*moments_Z['z'])/((1-delta*z*z)*root)) if r else coefficients["F"][0][0]*0
    # Differentiate the radial-flux polynomial independently by powers.
    flux_R=coefficients['F'][0][0]*0
    for n,row in enumerate(u):
        flux_R+=(((2*(n+1)-(1-delta))*z*row[0]
            -(1-z*z)*row[1])*r**n/(1-delta*z*z))
    jets['divergence_numerator']=((1-delta*z*z)*flux_R+(1-z*z)*jets['Uz_Z']
        - (1+delta)*z*jets['Uz']-2*z*r*jets['Uz_R'])
    jets['moments']=moments;jets['moments_Z']=moments_Z
    return jets


class ComponentPressureCore:
    """Local component jets for the complete preheat datum; no scalar installer."""
    def __init__(self,axis,pressure_datum,radial_degree=18):
        self.axis=axis;self.pressure_datum=pressure_datum;self.degree=int(radial_degree)
        self.precision=max(axis.precision,pressure_datum.precision)
    @lru_cache(maxsize=32)
    def coefficients(self, Z):
        return build_component_coefficients(self.axis,self.pressure_datum,Z,self.degree)
    def evaluate(self,R,Z):
        with mp.workdps(self.precision):
            r=mp.mpf(str(R));z=mp.mpf(str(Z))
            if not 0<=r<=mp.mpf('4.1')/self.axis.Lambda or abs(z)>=1:
                raise ValueError('Local core query outside its domain')
            coefficients=self.coefficients(Z)
            return evaluate_component_core_coefficients(coefficients,r,z,self.axis.delta)


    def physical_chart(self,R,Z,logq,*,nu='.01',phi=0):
        """Cartesian velocity components of this same local component core."""
        with mp.workdps(self.precision):
            r=mp.mpf(str(R));z=mp.mpf(str(Z));q=mp.exp(mp.mpf(str(logq)))
            viscosity=mp.mpf(str(nu));angle=mp.mpf(str(phi))
            if viscosity<=0:raise ValueError('Positive viscosity required')
            jets=self.evaluate(R,Z)
            radius=mp.sqrt(2*viscosity*q*r)
            axial=mp.sqrt(viscosity)*q**((1-self.axis.delta)/2)*z
            radial_scale=mp.sqrt(viscosity/q)
            scale=mp.sqrt(viscosity)*q**(-(1+self.axis.delta)/2)
            ur=radial_scale*jets['Ur'];ut=scale*jets['Utheta'];uz=scale*jets['Uz']
            return dict(xyz=(radius*mp.cos(angle),radius*mp.sin(angle),axial),
                uvw=(ur*mp.cos(angle)-ut*mp.sin(angle),ur*mp.sin(angle)+ut*mp.cos(angle),uz),
                pressure=viscosity*q**(-1-self.axis.delta)*jets['P'],tau=q*(1-z*z),
                scope='Component-valued local core only; no outer installation or global energy claim.')


def fixture():
    """Independent scalar replay of the original recurrence at resolved scales."""
    with mp.workdps(100):
        f=[mp.mpf(2),mp.mpf('.3')]+[mp.mpf(0)]*4
        u=[mp.mpf('.4'),mp.mpf('.2')]+[mp.mpf(0)]*4
        p=[mp.mpf(-1),mp.mpf('.2')]+[mp.mpf(0)]*4
        tail=[mp.mpf('.01'),mp.mpf('.03'),mp.mpf('-.02')]+[mp.mpf(0)]*3
        formal=core_coefficients('.3','.01',F0_Z_taylor=f,U0_Z_taylor=u,
            P0_Z_taylor=[PressurePolynomial({0:a,1:b}) for a,b in zip(p,tail)],
            radial_degree=4,precision=100,scalar_converter=PressurePolynomial)
        from lei_ren_part1_paper_core_adapter import CorePolynomial
        component=ComponentPressureCore(SimpleNamespace(precision=100,Lambda=mp.mpf(10),delta=mp.mpf('.01')),SimpleNamespace(precision=100),4)
        component.coefficients=lambda Z:formal
        errors=[]
        for parameter in (0,1,-1,mp.mpf('.5')):
            scalar=core_coefficients('.3','.01',F0_Z_taylor=f,U0_Z_taylor=u,
                P0_Z_taylor=[a+parameter*b for a,b in zip(p,tail)],radial_degree=4,precision=100)
            for name in ('F','Uz','P'):
                for formal_row,scalar_row in zip(formal[name],scalar[name]):
                    for a,b in zip(formal_row,scalar_row):errors.append(abs(a.evaluate(parameter)-b)/max(1,abs(b)))
            scalar_core=CorePolynomial(lambda Z:scalar,Lambda=10,delta='.01',precision=100)
            a=component.evaluate('.1','.3');b=scalar_core.evaluate('.1','.3')
            for name,value in a.items():
                if isinstance(value,PressurePolynomial) and name in b:
                    errors.append(abs(value.evaluate(parameter)-b[name])/max(1,abs(b[name])))
            for name in a['moments']:
                errors.append(abs(a['moments'][name].evaluate(parameter)-b['moments'][name])/max(1,abs(b['moments'][name])))
            if parameter==1:
                ca=component.physical_chart('.1','.3','-4');cb=scalar_core.physical_chart('.1','.3','-4')
                for x,y in zip(ca['uvw'],cb['uvw']):errors.append(abs(x.evaluate()-y)/max(1,abs(y)))
        maximum=max(errors)
        assert maximum<mp.mpf('1e-90'),maximum
        return {'maximum_scaled_difference':mp.nstr(maximum,30),'parameters':[0,1,-1,.5]}


if __name__=='__main__':
    print(fixture())
