"""Coefficientwise actual core atoms for the original annular inlet.

Integrates fresh coupled Phi/Uz radial rows, including differentiated
product tails. The bridge's raw axial V=Uz is distinct from the pressure
primitive int Phi^2. No covers or midpoint parameters define the atoms.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_rooted_core_field import CompliantRootedCoreField
from lei_ren_part1_paper_compliant_core_physical_field import symmetric,intersection
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def finite_atom_coefficients(c,phi,uz,order,radius=4):
    """Ordinary axial Taylor coefficients, integrated in scaled rho."""
    r=c.mpf(radius)
    integral=lambda n,weight:r**(n+weight+1)/(n+weight+1)
    def linear(rows,weight,normalization):
        return [sum((row[k]*integral(n,weight) for n,row in enumerate(rows)),c.mpf(0))/normalization
                for k in range(order+1)]
    def product(left,right,weight,normalization):
        return [sum((left[n][ell]*right[m][k-ell]*integral(n+m,weight)
                     for ell in range(k+1) for n in range(len(left)) for m in range(len(right))),c.mpf(0))/normalization
                for k in range(order+1)]
    return dict(H=linear(phi,1,8),M=linear(uz,0,4),K=product(phi,uz,1,8),
                A=product(uz,uz,0,4),B=product(phi,phi,1,16),C=product(phi,phi,0,4))


def product_tail(c,left_norm,right_norm,left_tail,right_tail,k):
    # Taylor coefficient convolution: derivative factorials are absent.
    return sum((left_norm[ell]*right_tail[k-ell]+left_tail[ell]*right_norm[k-ell]
                +left_tail[ell]*right_tail[k-ell] for ell in range(k+1)),c.mpf(0))


def bessel_radial_tail_coefficients(c,chi,degree,order,radius,radial_order,axial_weight):
    """Factorial radial-model tail using a scaled axial coefficient norm."""
    if order>chi.order or axial_weight<=0:raise ValueError('Finite positive axial Taylor scale required')
    lo,hi=endpoints(chi[0])
    if lo<0 or hi>1:raise ArithmeticError('Real chi0 must lie in[0,1]')
    variation=sum((abs(chi[k])*axial_weight**k for k in range(1,order+1)),c.mpf(0))
    N=int(degree);n=N+1;i=int(radial_order);r=c.mpf(radius);tau=c.mpf(axial_weight)
    tails=[]
    for k in range(order+1):
        first=(r**(n-i)*(k+1)*(n+1)**k*(1+variation)**k
               /(tau**k*2**n*math.factorial(n-i)*math.factorial(n+1)))
        ratio=r/2*(c.mpf(n+2)/(n+1))**k/((n+1-i)*(n+2))
        if endpoints(ratio)[1]>=1:raise ArithmeticError('Scaled axial factorial radial tail does not contract')
        tails.append(first/(1-ratio))
    return tails


class CompliantCoreIntegralAtoms:
    def __init__(self):
        self.field=CompliantRootedCoreField();self.rebuild=self.field.rebuild
        self.core=self.field.core;self.ctx=self.field.ctx;self.cache={}
        tau_model=self.core.h/(100*max(endpoints(self.rebuild.phi_norm)[1],mp.mpf(1)))
        tau_micro=endpoints(self.core.sigma)[0]/1000
        self.axial_weight=self.ctx.mpf(min(endpoints(tau_model)[0],tau_micro))
        self.hashes=dict(self.field.hashes);name=PREFIX+'rooted_core_field_check.json'
        receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or receipt['implicit_source_sha256']!=self.core.source:
            raise ValueError('Accepted same-source core required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Core changed: '+path)
        self.hashes.update(receipt['input_hashes'])
        for path in (name,Path(__file__).name):self.hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()

    def model_tail_coefficients(self,Z,degree,order):
        """Original Bessel-model radial tail, through axial Taylor order6.

        For real Z, 0<=chi0<=1. Let M=sum |chi_k|, k1..order.
        Coefficient k of chi^n has modulus at most
        (k+1)*(n+1)^k*(1+M)^k. The factorial radial series then contracts.
        """
        c=self.ctx;z=IntervalTaylor.variable(c,Z,order);one=IntervalTaylor.constant(c,1,order)
        H=z*((1-self.core.delta)/2)+(one-z*z)*(4*z+self.core.j)
        square=list((H*H).coefficients);square[0]=H[0]**2
        h2=IntervalTaylor(c,square);denominator=h2+self.core.sigma**2
        chi=h2/denominator
        return bessel_radial_tail_coefficients(c,chi,degree,order,4,0,self.axial_weight)

    def atoms_from_packet(self,packet,order,shared_root=False):
        c=self.ctx;N=packet['radial_degree'];rows=packet['rows'];r=c.mpf(4)
        phi=[row[:order+1] for row in rows['A'][:N+1]]
        uz=[row[:order+1] for row in rows['Uz'][:N+1]]
        finite=finite_atom_coefficients(c,phi,uz,order)
        pnorm=[sum((abs(row[k])*r**n for n,row in enumerate(phi)),c.mpf(0)) for k in range(order+1)]
        unorm=[sum((abs(row[k])*r**n for n,row in enumerate(uz)),c.mpf(0)) for k in range(order+1)]
        if shared_root:
            if 2*(N+1)<=order:raise ValueError('Root model-tail valuation insufficient')
            model=[c.mpf(0)]*(order+1)
        else:model=self.model_tail_coefficients(packet['Z'],N,order)
        ep=[];eu=[]
        for k in range(order+1):
            factor=tail_factor(c,degree=N,radial_order=0,axial_order=k,radius=r,h=self.core.h)['tail_per_Xh_norm']
            correction=self.core.correction*factor/math.factorial(k)
            ep.append(model[k]+correction);eu.append(self.core.epsilon*correction)
        errors=dict(H=ep,M=eu,
            K=[product_tail(c,pnorm,unorm,ep,eu,k) for k in range(order+1)],
            A=[product_tail(c,unorm,unorm,eu,eu,k) for k in range(order+1)],
            B=[product_tail(c,pnorm,pnorm,ep,ep,k)/2 for k in range(order+1)],
            C=[product_tail(c,pnorm,pnorm,ep,ep,k) for k in range(order+1)])
        atoms={name:[value+symmetric(c,errors[name][k]) for k,value in enumerate(values)] for name,values in finite.items()}
        # Positive independent same-source Phi bounds sharpen only the value.
        floor=c.mpf(endpoints(self.core.phi_floor)[0]);ceiling=c.mpf(endpoints(self.core.phi_ceiling)[1])
        for name,lo,hi in (('H',floor,ceiling),('B',floor**2/2,ceiling**2/2),('C',floor**2,ceiling**2)):
            atoms[name][0]=intersection(c,atoms[name][0],c.mpf([endpoints(lo)[0],endpoints(hi)[1]]))
        endpoint={}
        for name,polynomial,tail in (('phi',phi,ep),('v',uz,eu)):
            endpoint[name]=[sum((row[k]*r**n for n,row in enumerate(polynomial)),c.mpf(0))+symmetric(c,tail[k])
                            for k in range(order+1)]
        return dict(Z=packet['Z'],radial_domain=[0,4],axial_order=order,radial_degree=N,
            axial_coefficient_units='Taylor coefficient k = ordinary Z derivative / k!',
            finite_integrated_coefficients=finite,directed_integrated_tail_coefficient_bounds=errors,
            actual_core_atom_axial_coefficients=atoms,
            actual_core_atom_ordinary_axial_derivatives={name:[value*math.factorial(k) for k,value in enumerate(values)]
                for name,values in atoms.items()},
            core_exit_profile_axial_coefficients=dict(**endpoint,mean=atoms['M']),
            pressure_primitive_at_exit_axial_coefficients=[4*value for value in atoms['C']],
            pressure_identity='V_pressure(4,Z)=int_0^4 Phi^2 drho=4C; P-P0=epsilon F0^2 V_pressure=R F0^2 C',
            raw_bridge_axial_V_is_Uz_not_pressure_primitive=True,
            finite_profile_absolute_coefficient_bounds=dict(Phi=pnorm,Uz=unorm),
            radial_tail_coefficient_bounds=dict(Phi=ep,Uz=eu,Phi_model=model),
            exact_shared_H_root_used_before_enclosure=shared_root,
            fresh_rows_radial_depth=N,old_finite_rows_read=False,
            cover_endpoints_used_as_atom_values=False,point_parameter_representatives_selected=False,
            same_coupled_rows_and_nonlinear_tails=True,all_annular_source_values_resolved=False,
            temporal_recursion=False)

    def core_atoms(self,Z,z_order=6,degree=24):
        if not isinstance(z_order,int) or not 0<=z_order<=6 or not isinstance(degree,int) or degree<6:
            raise ValueError('Axial order0..6 and radial degree>=6 required')
        z=self.ctx.mpf(Z);key=(z._mpi_,z_order,degree)
        if key not in self.cache:
            packet=self.rebuild.rebuild(z,degree,depth=max(1,z_order))
            self.cache[key]=self.atoms_from_packet(packet,z_order)
        return self.cache[key]

    def root_atoms(self,z_order=6,degree=24):
        if not isinstance(z_order,int) or not 0<=z_order<=6 or not isinstance(degree,int) or degree<6:
            raise ValueError('Axial order0..6 and radial degree>=6 required')
        key=('exact_root',z_order,degree)
        if key not in self.cache:
            packet=self.field.build_root_rows(degree,depth=max(1,z_order))
            self.cache[key]=self.atoms_from_packet(packet,z_order,shared_root=True)
        return self.cache[key]

    def report(self):
        points={}
        for z in ('.3','.5','-.5','0'):
            points[z]=self.core_atoms(z);print('Actual coefficientwise core atoms: Z='+z+', axial6',flush=True)
        points['exact_shared_root']=self.root_atoms()
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,core_atom_packets=points,
            actual_core_atoms_coefficientwise_integrated=True,actual_core_atoms_axial6_available=True,
            same_nonlinear_tail_product_derivatives_included=True,
            source_domain='Same analytic core, Z in[-1,1], rho integration in[0,4]; API accepts source intervals',
            original_signed_bridge_integrals_resolved=False,all_annular_source_values_resolved=False,
            full_point_physical_field_evaluation=False,measured_blowup_dynamics=False,
            physical_energy_integral_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
            input_hashes=self.hashes)


def run():
    with mp.workdps(400):result=CompliantCoreIntegralAtoms().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
