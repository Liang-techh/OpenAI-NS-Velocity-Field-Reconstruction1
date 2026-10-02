"""Actual five-bump patch mixed x/Z and logR/Z derivatives through4.

The same implicit axial5 coefficient family and actual initial histories
are used. Exact Rm remains a formal common positive radial prefactor.
This stage is a local spatial profile certificate, not the whole field.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_actual_moment_patch import CompliantActualMomentPatch,IntervalTaylor
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square,derivative
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def power_derivatives(c,x,p,order=4):
    x=c.mpf(x);p=c.mpf(p);out=[];coefficient=c.mpf(1)
    for k in range(order+1):
        out.append(coefficient*x**(p-k));coefficient*=p-k
    return out


def product_derivatives(a,b,k):
    return sum((a[j]*b[k-j]*math.comb(k,j) for j in range(k+1)),a[0]*0)


def stirling_second(n,k):
    if n==0:return int(k==0)
    if k<=0 or k>n:return 0
    return k*stirling_second(n-1,k)+stirling_second(n-1,k-1)


def log_radial_derivatives(c,x,rows):
    return [rows[0]]+[sum((rows[j]*(stirling_second(k,j)*x**j) for j in range(1,k+1)),rows[0]*0)
        for k in range(1,len(rows))]


def patch_mixed(c,x,Z,delta,am,invAm2,invP2,H,V,g,initial,p0):
    """Differentiate actual physical primitive RHSs with all prefactors.

    H,V,g contain ordinary x derivatives with axial5 Taylor jets. No
    finite-difference derivative of an enclosure or capped shape is used.
    """
    if any(len(rows)!=5 or any(row.order!=5 for row in rows) for rows in (H,V,g)):
        raise ValueError('Ordinary x orders0..4 with axial5 required')
    z=IntervalTaylor.variable(c,c.mpf(Z),5);zero=z*0;one=zero+1
    sqrtx=power_derivatives(c,x,'.5');inverses=power_derivatives(c,x,-1)
    U=[am*row for row in H]
    G=[initial['mass']*x]+V[:4]
    mass=[sum((G[j]*inverses[k-j]*math.comb(k,j) for j in range(k+1)),zero) for k in range(5)]
    theta=[initial['theta']];mixed=[initial['mixed']];energy=[initial['energy']];pressure=[initial['pressure']]
    HV=[product_derivatives(H,V,k) for k in range(4)]
    H2=[product_derivatives(H,H,k) for k in range(4)]
    g2=[product_derivatives(g,g,k) for k in range(4)]
    for k in range(4):
        theta.append(sum((H[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),zero))
        mixed.append(sum((HV[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),zero))
        energy.append(g2[k]*invAm2-H2[k]/2)
        pressure.append(sum((H2[j]*(inverses[k-j]*math.comb(k,j)/2) for j in range(k+1)),zero))
    L=1-square(z)*delta;d=1-square(z)
    Q=[(2*z*V[k]-(z*mass[k])*(1-delta)-d*derivative(mass[k]))/L for k in range(5)]
    Ur=[sum((Q[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),Q[0]*0) for k in range(5)]
    P=[p0+square(am)*pressure[0]]+[square(am)*row for row in pressure[1:]]
    xder=[c.mpf(x),c.mpf(1),c.mpf(0),c.mpf(0),c.mpf(0)]
    physical_primitives=dict(Mz_over_Rm=G,
        Mtheta_over_sqrt2_Rm_1p5_Pstar=[am*row for row in theta],
        Mtheta_z_over_sqrt2_Rm_1p5_Pstar=[am*row for row in mixed],
        Mztheta_over_Rm_Pstar2=[((z*G[k])*8-square(z)*(16*xder[k]))*invP2+square(am)*energy[k] for k in range(5)],
        Mp_over_Pstar2=[square(am)*row for row in pressure])
    physical=dict(Utheta_over_Pstar=U,Uz=V,Ur_over_sqrt_Rm_over_2=Ur,P_over_Pstar2=P)
    grid=lambda rows:{'x'+str(k)+'_Z'+str(n):row[n]*math.factorial(n)
        for k,row in enumerate(rows) for n in range(5-k)}
    ygrid=lambda rows:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n)
        for k,row in enumerate(log_radial_derivatives(c,x,rows)) for n in range(5-k)}
    return dict(physical_velocity_pressure_x_derivative_Taylor={n:[row.truncate(4-k) for k,row in enumerate(rows)] for n,rows in physical.items()},
        physical_velocity_pressure_y_derivative_Taylor={n:[row.truncate(4-k) for k,row in enumerate(log_radial_derivatives(c,x,rows))] for n,rows in physical.items()},
        physical_velocity_pressure_x_Z_mixed4={n:grid(rows) for n,rows in physical.items()},
        physical_velocity_pressure_y_Z_mixed4={n:ygrid(rows) for n,rows in physical.items()},
        physical_five_primitive_x_Z_mixed4={n:grid(rows) for n,rows in physical_primitives.items()},
        actual_normalized_primitive_x_derivative_axial5=dict(mean=mass,theta=theta,mixed=mixed,energy=energy,pressure=pressure),
        actual_Q_x_derivative_axial4=Q,
        radial_prefactors_differentiated_before_mixed_grid=True)


class CompliantActualPatchMixedC4:
    def __init__(self):
        self.patch=CompliantActualMomentPatch();self.ctx=c=self.patch.ctx
        self.family=self.patch.family;self.source=self.patch.source;self.hashes=dict(self.patch.hashes)
        name=PREFIX+'flat_pulse_derivatives_check.json';flat=json.loads((HERE/name).read_bytes())
        if not (flat['all_passed'] and flat['original_radial_shape_derivatives_C4_available']
                and flat['actual_five_defect_family_sha256']==self.family):
            raise ValueError('Original compact beta derivative certificate required')
        for path,digest in flat['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Patch beta source changed: '+path)
        self.hashes.update(flat['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.radius=c.mpf(1)/40;self.N=c.mpf(endpoints(self.patch.repair.normalization))
        self.invP2=c.exp(-2*self.patch.core.logP)
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def gamma(self,x,center):
        c=self.ctx;argument=(c.mpf(x)-c.mpf(center))/self.radius
        rows=beta_jets(c,argument)
        derivatives=[rows[k]*(math.factorial(k)/(self.radius**(k+1)*self.N)) for k in range(5)]
        old=self.patch.repair.beta(c.mpf(endpoints(c.mpf(x)-c.mpf(center))))
        derivatives[0]=intersection(c,derivatives[0],c.mpf(endpoints(old[0])))
        derivatives[1]=intersection(c,derivatives[1],c.mpf(endpoints(old[1])))
        return derivatives

    def evaluate(self,x,Z):
        c=self.ctx;x=c.mpf(x);Z=c.mpf(Z);parent=self.patch.evaluate(x,Z)
        inverse,data=self.patch.coefficients(Z);h=inverse['controls'];z=IntervalTaylor.variable(c,Z,5)
        zero=z*0;one=zero+1;am=(1+square(z)).reciprocal()*c.exp(c.mpf('-.6'))
        H=[one*value for value in power_derivatives(c,x,'.1')];g=[zero for _ in range(5)]
        gamma=[self.gamma(x,c.mpf(i)/4) for i in (5,6,7)]
        for k in range(5):
            for i in range(3):H[k]+=h[i+2]*gamma[i][k]
            g[k]=h[0]*gamma[0][k]+h[1]*gamma[2][k]
        V=[z*4+g[0]]+g[1:]
        jet=lambda name:IntervalTaylor(c,parent[name])
        # Zeroth histories are the accepted ACTUAL partial integrals.
        initial=dict(mass=jet('Mz_over_R_axial5'),theta=jet('Mtheta_over_sqrt2_Rm_1p5_Am_axial5'),
            mixed=jet('Mtheta_z_over_sqrt2_Rm_1p5_Am_axial5'),
            energy=jet('Mztheta_minus_8ZMz_plus_16Z2R_over_RmAm2_axial5'),pressure=jet('Mp_over_Am2_axial5'))
        p0=jet('original_P0_axial5')
        mixed=patch_mixed(c,x,Z,self.patch.core.delta,am,data['invAm2'],self.invP2,H,V,g,initial,p0)
        return dict(x=x,Z=Z,**mixed,actual_gamma_x_derivatives=gamma,
            actual_H_x_derivative_axial5=H,actual_V_x_derivative_axial5=V,
            actual_inherited_patch_packet=parent,
            same_actual_coefficient_family_and_P0_retained=True,ordinary_derivatives_not_Taylor_radial_coefficients=True,
            exact_formal_prefactors=dict(Utheta='Pstar',Uz='1',Ur='sqrt(Rm/2)',pressure='Pstar^2',
                radial_x_to_R='partial_R^k=Rm^-k*partial_x^k',log_radius='partial_y=x*partial_x; y=log(R)'),
            actual_patch_all_mixed_derivatives_total_order_le4_available=True,
            original_beta_support_flat_joins_certified=True,
            patch_Rm_and_Rh_functional_profile_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)

    def report(self):
        c=self.ctx;edges=[c.mpf(i)/40 for i in (49,51,59,61,69,71)]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            whole_actual_patch=self.evaluate(c.mpf([1,endpoints(c.exp(1))[1]]),[-1,1]),
            actual_Rm_inlet=self.evaluate(1,[-1,1]),actual_Rh_exit=self.evaluate(c.exp(1),[-1,1]),
            source_bump_edge_packets=[self.evaluate(edge,[-1,1]) for edge in edges],
            interior_packets=[self.evaluate(x,z) for x,z in (('1.25','0'),('1.5','.5'),('1.75','0'))],
            exact_original_radius=self.radius,original_normalization=self.N,
            original_beta_flat_join_proof='beta derivatives0..4 vanish at |argument|=1 by the accepted flat tail majorants; original normalized supports retained',
            actual_endpoint_join_proof='Rm is below all supports: field and all derivatives equal unpatched reference with inherited defects; after71/40 all partial moments equal the SAME implicit map full weights and terminal defects vanish as Z functions, so Rh matches exact reference through mixed4',
            actual_patch_all_mixed_derivatives_total_order_le4_available=True,
            patch_Rm_and_Rh_functional_profile_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantActualPatchMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual patch physical mixed x/Z and logR/Z derivatives through4 generated',flush=True)
    return result


if __name__=='__main__':run()
