"""Actual-core five exit primitives and explicit frozen comparison (9.12).

This is the auxiliary UNSMOOTHED frozen field, not the prescribed-shear bridge.
Positive F0 stays formal. The axial inertial direction is supplied as the
drive sqrt(R/2)*F*E, cancelling its unrepresentable inverse-F0 factor.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_core_physical_field import CompliantCorePhysicalField,gridkey
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
MTH,MZ,MTHZ,MZT,MP='angular','axial','angular_axial','axial_quadratic','pressure'


def derivative(jet):
    return IntervalTaylor(jet.ctx,[(k+1)*jet[k+1] for k in range(jet.order)])


def square(jet):
    values=list((jet*jet).coefficients)
    values[0]=jet[0]**2
    return IntervalTaylor(jet.ctx,values)


def axial_jet(c,grid,i=0,order=5):
    return IntervalTaylor(c,[grid[gridkey(i,k)]/math.factorial(k) for k in range(order+1)])


def dress(jet,ratios):
    """True Z derivative jet divided by the amplitude at its base point."""
    c=jet.ctx
    return jet*IntervalTaylor(c,[ratios[k]/math.factorial(k) for k in range(jet.order+1)])


def frozen_algebra(c,theta,Z,delta,exit):
    """Exact source formulas in theta=r/R; all inputs are ordinary axial jets.

    Normalizations: Mtheta/(F0 R^2), Mz/R, Mtheta_z/(F0 R^2),
    Mztheta/R=A-R F0^2 B, Mp/(R F0^2). F=F0*phi.
    """
    theta=c.mpf(theta); z=IntervalTaylor.variable(c,Z,5)
    z2=square(z); d=1-z2; L=1-z2*delta
    phi,v,mean=(exit[n] for n in ('phi','v','mean'))
    th2=theta**2
    moments={MTH:exit['H']*th2+phi*(1-th2),
        MZ:mean*theta+v*(1-theta),
        MTHZ:exit['K']*th2+phi*v*(1-th2),
        MZT:dict(axial=exit['A']*theta+square(v)*(1-theta),
                 swirl=exit['B']*th2+square(phi)*((1-th2)/2)),
        MP:exit['C']*theta+square(phi)*(1-theta)}
    m=moments[MZ]; a=moments[MZT]['axial']
    h=dress(moments[MTH],exit['F0_ratios']); k=dress(moments[MTHZ],exit['F0_ratios'])
    f=dress(phi,exit['F0_ratios'])
    b=dress(moments[MZT]['swirl'],exit['F0_squared_ratios'])
    p=dress(moments[MP],exit['F0_squared_ratios']); p0=exit['p0']
    W=1-(z*m)*(1-delta)-d*derivative(m)
    angular=h*(1-delta/2)-(z*derivative(h))*((1-delta)/2)-d*derivative(k)+(z*k)*(2*delta-1)
    D_over_R=(-W+angular/(2*f))/L
    drive_hydro=(-W*v+(m-z*derivative(m))*((1-delta)/2)
                 +(z*a)*(2*delta)-d*derivative(a))/(2*L)
    drive_pressure=((z*p0)*(2*(1+delta))-d*derivative(p0))/(2*L)
    drive_swirl=(-(z*b)*(2*delta)+d*derivative(b)
                 +(z*p)*(2*(1+delta))-d*derivative(p))/(2*L)
    Q=(2*z*v-(z*m)*(1-delta)-d*derivative(m))/L
    qv=((z*v)*(1+delta)-d*derivative(v))/L
    diff=mean-v
    qdiff=-((z*diff)*(1-delta)+d*derivative(diff))/L
    return dict(moments=moments,Q=Q,Q_constant=qv,Q_theta=qdiff,
                phi=phi,v=v,F=dress(phi,exit['F0_ratios']),P0=p0,PI=p,
                PI_constant=dress(square(phi),exit['F0_squared_ratios']),
                frozen_D_over_R=D_over_R,
                axial_direction_drive_parts=dict(hydro=drive_hydro,pressure=drive_pressure,swirl=drive_swirl))


def profile_grids(c,data,theta):
    """True y/Z derivatives divided by each prefactor at the base point."""
    grids={n:{} for n in ('Ur','Utheta','Uz','P0','PI')}
    for i in range(5):
        ur=data['Q_constant']*(c.mpf('.5')**i)+data['Q_theta']*(theta*c.mpf('-.5')**i)
        ut=data['F']*(c.mpf('.5')**i)
        vv=data['v'] if i==0 else data['v']*0
        pp=data['P0'] if i==0 else data['P0']*0
        pi=data['PI'] if i==0 else data['PI_constant']
        for k in range(5-i):
            index='y'+str(i)+'_Z'+str(k)
            for name,jet in zip(grids,(ur,ut,vv,pp,pi)):
                grids[name][index]=jet[k]*math.factorial(k)
    return grids


class CompliantFrozenComparisonField:
    def __init__(self):
        self.core=CompliantCorePhysicalField();self.ctx=c=self.core.ctx
        self.family=self.core.family;self.source=self.core.source;self.delta=self.core.delta
        self.hashes=dict(self.core.hashes)
        name=PREFIX+'core_physical_field_check.json'; check=json.loads((HERE/name).read_bytes())
        if not (check['all_passed'] and check['whole_Z_analytic_core_profile_enclosures_through_order5_available']
                and check['actual_five_defect_family_sha256']==self.family and check['implicit_source_sha256']==self.source):
            raise ValueError('Accepted same-source nonlinear analytic core required')
        for path,digest in check['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Core comparison prerequisite changed: '+path)
        self.hashes.update(check['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        self.r=4*self.core.epsilon
        self.cache={}

    def exit_primitives(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=Z._mpi_
        if key in self.cache:return self.cache[key]
        at=self.core.normalized_jets(4,Z);cover=self.core.normalized_jets([0,4],Z)
        phi=axial_jet(c,at['Phi']);v=axial_jet(c,at['Uz']);mean=axial_jet(c,at['Mz_over_R'])
        pc=axial_jet(c,cover['Phi']);vc=axial_jet(c,cover['Uz'])
        # Integral enclosures use exact positive masses 8 and 4. The
        # integrals stay source definitions, not midpoint quadrature.
        # H=1/8 int rho Phi; K=1/8 int rho Phi Uz;
        # A=1/4 int Uz^2; B=1/16 int rho Phi^2; C=1/4 int Phi^2.
        inputs=dict(phi=phi,v=v,mean=mean,H=pc,K=pc*vc,A=square(vc),B=square(pc)/2,C=square(pc),
            p0=at['source']['pressure'].truncate(5),
            F0_ratios=at['source']['relative_F0_derivatives'],
            F0_squared_ratios=at['source']['relative_F0_squared_derivatives'])
        self.cache[key]=inputs
        return inputs

    def packet(self,Z,theta):
        c=self.ctx;Z=c.mpf(Z);theta=c.mpf(theta)
        floor=endpoints(self.r/110)[0]
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(theta)[0]<floor or endpoints(theta)[1]>1:
            raise ValueError('Frozen comparison requires Z in[-1,1], r/110<=theta=r/R<=1')
        exit=self.exit_primitives(Z);data=frozen_algebra(c,theta,Z,self.delta,exit)
        grids=profile_grids(c,data,theta)
        def coeffs(value):
            if isinstance(value,IntervalTaylor):return list(value.coefficients)
            return {name:coeffs(jet) for name,jet in value.items()}
        core_true={MTH:dress(exit['H'],exit['F0_ratios']),MZ:exit['mean'],
            MTHZ:dress(exit['K'],exit['F0_ratios']),
            MZT:dict(axial=exit['A'],swirl=dress(exit['B'],exit['F0_squared_ratios'])),
            MP:dress(exit['C'],exit['F0_squared_ratios'])}
        frozen_true={MTH:dress(data['moments'][MTH],exit['F0_ratios']),MZ:data['moments'][MZ],
            MTHZ:dress(data['moments'][MTHZ],exit['F0_ratios']),
            MZT:dict(axial=data['moments'][MZT]['axial'],swirl=dress(data['moments'][MZT]['swirl'],exit['F0_squared_ratios'])),
            MP:dress(data['moments'][MP],exit['F0_squared_ratios'])}
        return dict(Z=Z,theta_r_over_R=theta,normalization_radius_r=self.r,
            core_exit_normalized_integral_shape_axial_coefficients=coeffs({k:v for k,v in exit.items() if isinstance(v,IntervalTaylor)}),
            core_exit_true_moment_derivatives_divided_by_basepoint_prefactors=coeffs(core_true),
            frozen_moment_shape_axial_coefficients=coeffs(data['moments']),
            frozen_true_moment_derivatives_divided_by_basepoint_prefactors=coeffs(frozen_true),
            moment_prefactors=dict(angular='F0*R^2',axial='R',angular_axial='F0*R^2',
                axial_quadratic='R*A-R^2*F0^2*B (separate axial/swirl parts)',pressure='R*F0^2'),
            mixed_profile_derivatives_total_order_le4=grids,
            profile_prefactors=dict(Ur='sqrt(R/2)',Utheta='sqrt(2R)*F0',Uz='1',P0='Pstar^2',PI='R*F0^2'),
            frozen_D_over_R_coefficients=list(data['frozen_D_over_R'].coefficients),
            axial_direction_drive_coefficients=coeffs(data['axial_direction_drive_parts']),
            axial_drive_definition='sqrt(R/2)*F*Ef = R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl',
            F0_log_enclosure=c.mpf([endpoints(-self.core.logC-self.core.Lambda*self.core.Gbar)[0],endpoints(-self.core.logC)[1]]),
            amplitude_is_exact_positive_formal_source=True,
            actual_unsmoothed_frozen_comparison=True,actual_smooth_comparison=False,actual_prescribed_shear_bridge=False)

    def report(self):
        c=self.ctx;low=endpoints(self.r/110)[0]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            source_equations='paper(9.12),(9.13); Ra=4/Lambda; all five core primitives retained; P=P0+Mp',
            core_exit_integral_definitions=dict(H='(1/8)*integral_0^4 rho*Phi drho',
                K='(1/8)*integral_0^4 rho*Phi*Uz drho',A='(1/4)*integral_0^4 Uz^2 drho',
                B='(1/16)*integral_0^4 rho*Phi^2 drho',C='(1/4)*integral_0^4 Phi^2 drho'),
            radial_integrals_are_enclosures_not_point_recomputed_coefficients=True,
            actual_core_exit_five_primitive_axial5_enclosures_available=True,
            frozen_comparison_all_five_primitives_available=True,
            frozen_comparison_profile_mixed4_available=True,
            frozen_inertial_direction_axial4_enclosures_available=True,
            whole_domain=self.packet([-1,1],[low,1]),
            exit_domain=self.packet([-1,1],1),
            terminal_domain=self.packet([-1,1],self.r/c.mpf([100,110])),
            samples=[self.packet(z,theta) for z in ('0','.5') for theta in ('1','.5')],
            frozen_direction_used_as_actual_bridge=False,
            actual_smooth_comparison_installed=False,actual_prescribed_shear_bridge_installed=False,
            core_inner_annulus_interfaces_certified=False,whole_outer_cone_certified=False,
            full_cartesian_vector_derivatives_certified=False,temporal_recursion=False,
            next_dependency='smooth comparison(9.23) with its own moments -> actual prescribed-shear bridge(9.26), actual moment recovery and interfaces',
            input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantFrozenComparisonField().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Same-source core exit five primitives and unsmoothed frozen comparison/direction generated',flush=True)
    return result


if __name__=='__main__':run()
