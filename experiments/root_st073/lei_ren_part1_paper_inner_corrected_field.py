"""Actual five-bump field adapter with explicit unresolved defect intervals."""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_restore import AxialRestore, defect_receipt
from lei_ren_part1_paper_inner_moment_map import MomentMap
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_exit_field import cone_receipt


class CorrectedInnerField:
    def __init__(self,provider,*,order=128,map_precision=160):
        self.provider=provider; self.core=provider.core; self.precision=provider.precision
        self.r=provider.reshape.switches.r
        self.moment_map=MomentMap(precision=map_precision,order=order)

    @lru_cache(maxsize=16)
    def coefficients(self,Zkey):
        with mp.workdps(self.precision):
            z=mp.mpf(Zkey)
            data=defect_receipt(self.provider,z,digits=self.precision)
            unresolved=set(data['cancellation_unresolved_entries'])
            original=[mp.mpf(v) for v in data['centered_defects_from_finite_subtraction']]
            originalZ=[mp.mpf(v) for v in data['centered_Z_derivatives_from_finite_subtraction']]
            bounds=[mp.mpf(v) for v in data['conditional_centered_value_bounds']]
            C1bounds=[mp.mpf(v) for v in data['conditional_centered_C1_bounds_using_sampled_data']]
            # Symmetric intervals have representative zero; this choice is
            # explicit and NEVER an assertion that the physical defect is zero.
            d=mp.matrix([0 if i+1 in unresolved else v for i,v in enumerate(original)])
            dZ=mp.matrix([0 if i+1 in unresolved else v for i,v in enumerate(originalZ)])
            radii=[bounds[i] if i+1 in unresolved else mp.mpf(0) for i in range(5)]
            radiiZ=[C1bounds[i] if i+1 in unresolved else mp.mpf(0) for i in range(5)]
            entry=self.provider.evaluate_phase(2,z)
            Rm=entry['R']; Am=mp.sqrt(2*Rm)*entry['F']; zeta=-2*z/(1+z*z)
            solution=self.moment_map.solve(d,Am)
            h=solution['h']; hZ=self.moment_map.solve_Z(h,dZ,Am,zeta*Am)
            CA=self.moment_map.inverse_l1_norm
            CQ=self.moment_map.quadratic_l1_bound(Am)
            ball=2*CA*(sum(abs(v) for v in d)+sum(radii))
            factor=2*CA*CQ*ball
            if factor>=1:
                raise ArithmeticError('Conditional pointwise contraction envelope failed')
            uncertainty=CA*sum(radii)/(1-factor)
            inverse_bound=CA/(1-factor)
            coeffZ_uncertainty=inverse_bound*(sum(radiiZ)
                +2*CQ*uncertainty*sum(abs(v) for v in hZ)
                +2*abs(zeta)*CQ*uncertainty*(2*sum(abs(v) for v in h)+uncertainty))
            return dict(data=data,d=d,dZ=dZ,h=h,hZ=hZ,Rm=Rm,Am=Am,zeta=zeta,
                unresolved=sorted(unresolved),input_radii=radii,input_Z_radii=radiiZ,
                solution=solution,conditional_coeff_l1_uncertainty=uncertainty,
                conditional_coeff_Z_l1_uncertainty=coeffZ_uncertainty,
                contraction_factor=factor,contraction_ball=ball)

    def evaluate_x(self,x,Z):
        """x=R/Rm in [1,e], including partial supports and terminal join."""
        with mp.workdps(self.precision):
            x=mp.mpf(x); z=mp.mpf(Z)
            if not 1<=x<=mp.e or not abs(z)<1:
                raise ValueError('Require 1<=R/Rm<=e and |Z|<1')
            c=self.coefficients(mp.nstr(z,self.provider.reshape.precision))
            v=self.provider.evaluate_phase(2+mp.log(x),z)
            h=c['h']; hZ=c['hZ']; Am=c['Am']; AmZ=c['zeta']*Am; Rm=c['Rm']
            bumps=self.moment_map.bump_values(x)
            gamma=bumps['gamma']; gx=bumps['gamma_x']
            f=sum(h[i+2]*gamma[i] for i in range(3))
            fx=sum(h[i+2]*gx[i] for i in range(3))
            fZ=sum(hZ[i+2]*gamma[i] for i in range(3))
            g=h[0]*gamma[0]+h[1]*gamma[2]
            gZ=hZ[0]*gamma[0]+hZ[1]*gamma[2]
            gy=x*(h[0]*gx[0]+h[1]*gx[2])
            partial_x=min(x,mp.mpf(2))
            q=self.moment_map.partial_integrals(partial_x,h,Am)
            qZ=self.moment_map.partial_Z(partial_x,h,hZ,Am,AmZ)
            thetaScale=mp.sqrt(2)*Rm**mp.mpf('1.5')*Am
            energyScale=Rm*Am*Am; zeta=c['zeta']
            mass=Rm*q[0]; massZ=Rm*qZ[0]
            theta=thetaScale*q[2]; thetaZ=thetaScale*(qZ[2]+zeta*q[2])
            mixed=thetaScale*q[1]+4*z*theta
            mixedZ=thetaScale*(qZ[1]+zeta*q[1])+4*theta+4*z*thetaZ
            energy=energyScale*q[3]+8*z*mass
            energyZ=energyScale*(qZ[3]+2*zeta*q[3])+8*mass+8*z*massZ
            pressure=Am*Am*q[4]; pressureZ=Am*Am*(qZ[4]+2*zeta*q[4])
            delta=dict(theta=theta,z=mass,theta_z=mixed,z_theta=energy,p=pressure)
            deltaZ=dict(theta=thetaZ,z=massZ,theta_z=mixedZ,z_theta=energyZ,p=pressureZ)
            m={k:v['moments'][k]+delta[k] for k in delta}
            mz={k:v['momentsZ'][k]+deltaZ[k] for k in deltaZ}
            u=Am*(x**mp.mpf('.1')+f); uy=Am*(x**mp.mpf('.1')/10+x*fx)
            uZ=AmZ*(x**mp.mpf('.1')+f)+Am*fZ
            V=4*z+g; VZ=4+gZ; R=Rm*x
            F=u/mp.sqrt(2*R); FZ=uZ/mp.sqrt(2*R)
            a=1-2*uy/u; b=2*gy/u
            P=v['P']+pressure; PZ=v['PZ']+pressureZ
            stress=evaluate_mp_stress(v['logR'],z,self.core.delta,
                Utheta=u,Uz=V,Utheta_y=uy,Utheta_Z=uZ,Uz_y=gy,Uz_Z=VZ,
                moments=m,moments_Z=mz,P=P,P_Z=PZ,shear_theta=-a*F,
                shear_z=mp.sqrt(2*R)*gy/R,precision=self.precision)
            # Keep the separate quadratic primitives for later energy/exterior work.
            g2=self.moment_map.axial_square_increment(partial_x,h)
            g2Z=self.moment_map.axial_square_increment_Z(partial_x,h,hZ)
            axial_change=8*z*mass+Rm*g2
            axialZ_change=8*mass+8*z*massZ+Rm*g2Z
            raw=v['raw_quadratic_integrals']
            raw=dict(axial=raw['axial']+axial_change,axial_Z=raw['axial_Z']+axialZ_change,
                swirl=raw['swirl']+axial_change-energy,
                swirl_Z=raw['swirl_Z']+axialZ_change-energyZ)
            return dict(R=R,logR=v['logR'],Z=z,x=x,F=F,FZ=FZ,Uz=V,UZ=VZ,
                Ur=stress['U_r'],P=P,PZ=PZ,moments=m,momentsZ=mz,stress=stress,
                raw_quadratic_integrals=raw,a=a,b=b,g_y=uy/u-mp.mpf('.5'),u_y=gy,
                F_R=F*(uy/u-mp.mpf('.5'))/R,Uz_R=gy/R,
                relative_swirl_correction=f/x**mp.mpf('.1'),axial_correction=g,
                partial_centered_change=q,partial_centered_change_Z=qZ,
                unresolved_input_entries=c['unresolved'],
                conditional_coeff_l1_uncertainty=c['conditional_coeff_l1_uncertainty'],
                conditional_coeff_Z_l1_uncertainty=c['conditional_coeff_Z_l1_uncertainty'],
                source_constants_certified=False,global_moment_repair_certified=False,
                finite_energy_certified=False,region='five_bump_inner_correction')

    def evaluate(self,y,Z):
        """Exit-to-Rh profile at y=log(R/Ra); physical wrapper adds the core."""
        with mp.workdps(self.precision):
            logR=mp.log(self.r)+mp.mpf(y); z=mp.mpf(Z)
            offset=logR-mp.log(110); reshape=self.provider.reshape
            if offset<=0:
                return reshape.switches.evaluate(y,z)
            if offset<=reshape.reference_end:
                return reshape.evaluate_log_offset(offset,z)
            t=offset-reshape.reference_end
            if t<=2:return self.provider.evaluate_phase(t,z)
            if t>3:raise ValueError('Complete inner construction currently ends at Rh')
            return self.evaluate_x(mp.exp(t-2),z)

    def physical_field(self,*,nu='.01',T=0):
        """Cartesian velocity callable on the constructed inner domain."""
        from lei_ren_part1_paper_exit_field import LocalCoreExitField
        return LocalCoreExitField(self,nu=nu,T=T)


def build_candidate(*,continuous_pressure=False,pressure_order=192):
    from lei_ren_part1_paper_core_adapter import build_source_core
    from lei_ren_part1_paper_exit_comparison import Section923Comparison
    from lei_ren_part1_paper_exit_tangents import ExitTangents
    from lei_ren_part1_paper_exit_continuation import ExitContinuation
    from lei_ren_part1_paper_exit_switches import ExitSwitches
    from lei_ren_part1_paper_long_reshape import LongReshape
    with mp.workdps(260):
        bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
            logC='5e151',logPstar='14',delta='1e-200',
            continuous_pressure=continuous_pressure,pressure_order=pressure_order)
        hb=mp.exp(-100-100*mp.mpf('1e152'))
        comparison=Section923Comparison(bundle,h_b=hb,transition_steps=16)
        switches=ExitSwitches(ExitContinuation(ExitTangents(comparison,epsilon=hb,
            steps=16,derivative_step='1e-45')),steps=16)
        return AxialRestore(LongReshape(switches,order=32),order=64)


def run():
    provider=build_candidate(); field=CorrectedInnerField(provider)
    path=Path(__file__).with_suffix('.json')
    with mp.workdps(field.precision):
        z=mp.mpf('.3'); c=field.coefficients(mp.nstr(z,provider.reshape.precision))
        n=lambda v:mp.nstr(v,60)
        receipt=dict(Z='.3',coefficients=[n(v) for v in c['h']],
            coefficients_Z=[n(v) for v in c['hZ']],
            unresolved_input_entries=c['unresolved'],
            unresolved_input_radii=[n(v) for v in c['input_radii']],
            input_defects=c['data'],conditional_coeff_l1_uncertainty=n(c['conditional_coeff_l1_uncertainty']),
            conditional_coeff_Z_l1_uncertainty=n(c['conditional_coeff_Z_l1_uncertainty']),
            conditional_contraction_factor=n(c['contraction_factor']),
            source_constants_certified=False,finite_energy_certified=False,
            scale_recursion_established=False,global_moment_repair_certified=False,rows=[])
        audit=MomentMap(precision=160,order=64)
        audit_solution=audit.solve(c['d'],c['Am'])
        receipt['coefficient_64_128_refinement']=[n(c['h'][i]-audit_solution['h'][i]) for i in range(5)]
        receipt['quadrature_refinement_is_not_certified_error_bound']=True
        receipt['moment_map_order']=128
        receipt['inverse_matrix_l1_norm']=n(field.moment_map.inverse_l1_norm)
        receipt['quadratic_l1_bound']=n(field.moment_map.quadratic_l1_bound(c['Am']))
        residual=c['d']+field.moment_map.evaluate(c['h'],c['Am'])
        receipt['representative_terminal_centered_residual']=[n(v) for v in residual]
        for x in ('1','1.2375','1.25','1.2625','1.4875','1.5','1.7375','1.75','1.7625','2',mp.e):
            print('computing bump '+n(mp.mpf(x)),flush=True)
            v=field.evaluate_x(x,z)
            row=dict(x=n(mp.mpf(x)),a=n(v['a']),b=n(v['b']),cone=cone_receipt(v),
                relative_swirl_correction=n(v['relative_swirl_correction']),
                axial_correction=n(v['axial_correction']))
            if x=='1':
                old=provider.evaluate_phase(2,z)
                receipt['entry_matching']={k:n(abs(v[k]-old[k])/abs(old[k]) if old[k] else abs(v[k]))
                    for k in ('F','FZ','Uz','UZ','Ur','P','PZ','g_y','u_y')}
                receipt['entry_moment_matching']={kind:{k:n(abs(v[kind][k]-old[kind][k])/abs(old[kind][k])
                    if old[kind][k] else abs(v[kind][k])) for k in v[kind]}
                    for kind in ('moments','momentsZ')}
            receipt['rows'].append(row)
            path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
            print(json.dumps(row),flush=True)
        # A higher quadrature order evaluates the same solved coefficients,
        # separating the algebraic residual from quadrature consistency.
        higher=MomentMap(precision=160,order=192)
        receipt['terminal_residual_at_order192_same_coefficients']=[n(v)
            for v in c['d']+higher.evaluate(c['h'],c['Am'])]
        x=mp.mpf('1.2375'); center=field.evaluate_x(x,z); divergence=[]
        for h in (mp.mpf('1e-5'),mp.mpf('5e-6')):
            values={i:field.evaluate_x(x*mp.exp(i*h),z)['Ur'] for i in (-2,-1,1,2)}
            derivative=(values[-2]-8*values[-1]+8*values[1]-values[2])/(12*h)
            R=center['R']; dt=field.core.delta; L=1-dt*z*z
            terms=[mp.sqrt(2*R)*derivative/R,center['Ur']/mp.sqrt(2*R),
                ((1-z*z)*center['UZ']-2*z*center['u_y']-(1+dt)*z*center['Uz'])/L]
            divergence.append(dict(h=n(h),mapped_q_divergence=n(sum(terms)),
                relative_cancellation=n(abs(sum(terms))/sum(abs(v) for v in terms))))
        receipt['bump_flank_radial_divergence']=divergence
        from lei_ren_part1_paper_exit_field import physical_chart
        chart=physical_chart(center,z,'-4',delta=field.core.delta,phi='.4',precision=field.precision)
        physical=field.physical_field()
        back=physical.evaluate(*chart['xyz'],-chart['tau'])
        receipt['physical_callable_roundtrip']=dict(region=back['region'],
            Z_error=n(back['Z']-z),R_relative_error=n((back['R']-center['R'])/center['R']),
            cylindrical_velocity_relative_errors=[n(abs(a-b)/abs(b) if b else abs(a))
                for a,b in zip(back['cylindrical_velocity'],chart['cylindrical_velocity'])],
            scope='inner domain through Rh; no heat exterior or finite-energy certificate')
        path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
        return receipt


if __name__=='__main__':run()
