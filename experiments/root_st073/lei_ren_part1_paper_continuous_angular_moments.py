"""Actual seeded angular and swirl-energy primitives through R_tail.

Forward normalized propagation uses the installed continuous angular
schedule. Signed bump moments and actual inner offsets are retained as
separate terms, so adding a tiny correction to a large baseline does not
silently certify cancellation. Heat and the other three moments remain open.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from numpy.polynomial.legendre import leggauss

from lei_ren_part1_paper_axial_correction import signed_log
from lei_ren_part1_paper_joined_outer import _mp


class ContinuousAngularMoments:
    def __init__(self, profile, *, order=96, mp_nodes=False):
        self.profile=profile
        self.schedule=profile.schedule
        self.precision=profile.precision
        self.order=int(order)
        if self.order<16:
            raise ValueError('At least 16 quadrature nodes required')
        self.correction=profile.angular_correction_provider
        with mp.workdps(self.precision):
            self.mu=_mp(str(self.schedule.mu))
            self.delta=_mp(str(self.schedule.delta))
            self.MP_nodes=bool(mp_nodes)
            if self.MP_nodes:
                nodes,weights=mp.gauss_quadrature(self.order,'legendre')
                self.nodes=[((x+1)/2,w/2) for x,w in zip(nodes,weights)]
            else:
                nodes,weights=leggauss(self.order)
                self.nodes=[((_mp(str(float(x)))+1)/2,_mp(str(float(w)))/2)
                            for x,w in zip(nodes,weights)]
        self.stages=[
            ('slope_transition_ref',None),('axial_turnoff','-.5'),
            ('slope_transition_mu',None),('power_buffer','pulse'),
            ('pulse_reserved','pulse'),('z_flatten',None),
            ('power_buffer_rel','pulse'),('steep_transition_in',None),
            ('steep_power','-1.5'),('steep_transition_out',None),
            ('waiting','waiting')]
        from lei_ren_part1_paper_continuous_heat_moments import ContinuousHeatMoments
        self.heat_provider=ContinuousHeatMoments(self)

    def _radius(self,y):
        return self.profile.log_at(self.schedule.logRref,str(y))

    def _preheat_y(self,logR):
        # Validate against the recorded absolute checkpoint before forming
        # a relative MP coordinate. Subtracting O(1e152) log radii can round
        # the exact endpoint a few last-place bits above its stored offset.
        if self.profile.offset(logR,self.schedule.logR_tail)>0:
            raise ValueError('Continuous preheat moments stop at R_tail')
        y=_mp(str(self.profile.offset(logR,self.schedule.logRref)))
        return min(y,_mp(str(self.schedule.y_tail)))

    @lru_cache(maxsize=4096)
    def _base(self,y_key):
        """Z-independent normalization, evaluated at the identical log R."""
        row=self.schedule.at_log_radius(self._radius(y_key),0)
        return _mp(row['log_angular_amplitude'])

    def _ratio_jet(self,y,z):
        """U_base(Z)/U_base(0) and its derivative before heat."""
        if y<=_mp(str(self.schedule.y_v)):
            q=mp.mpf(0)
        elif y>=_mp(str(self.schedule.y_f)):
            q=mp.mpf(1)
        else:
            from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
            q=ContinuousAxialPulse.sigma_pair((y-_mp(str(self.schedule.y_v)))/_mp(str(self.schedule.Tf)))[0]
        L=mp.log1p(z*z)
        G=mp.exp(-(1-q)*L)
        return G,-(1-q)*2*z/(1+z*z)

    def _transfer(self,state,left,end,z,slope):
        """Integrate the normalized forward ODE without huge-radius sampling."""
        length=end-left
        if not length:return list(state)
        if slope is not None:
            rate_s=(-mp.mpf('.5')-self.mu if slope=='pulse' else
                    -(1+self.delta)/2 if slope=='waiting' else mp.mpf(slope))
            G,LZ=self._ratio_jet(left,z)
            result=[]
            for i,(a,b) in enumerate(((mp.mpf('1.5'),1),(mp.mpf(1),2))):
                rate=a+b*rate_s
                damp=mp.exp(-rate*length)
                atom=-mp.expm1(-rate*length)/rate if rate else length
                result.append(state[i]*damp+G**b*atom)
            for i,(a,b) in enumerate(((mp.mpf('1.5'),1),(mp.mpf(1),2))):
                rate=a+b*rate_s
                damp=mp.exp(-rate*length)
                atom=-mp.expm1(-rate*length)/rate if rate else length
                result.append(state[i+2]*damp+b*G**b*LZ*atom)
            return result
        Eleft=self._base(mp.nstr(left,self.precision))
        Eend=self._base(mp.nstr(end,self.precision))
        result=[]
        for i,(a,b) in enumerate(((mp.mpf('1.5'),1),(mp.mpf(1),2))):
            damp=mp.exp(-a*length-b*(Eend-Eleft))
            value=state[i]*damp
            derivative=state[i+2]*damp
            for unit,weight in self.nodes:
                t=length*unit
                y=left+t
                Et=self._base(mp.nstr(y,self.precision))
                G,LZ=self._ratio_jet(y,z)
                atom=length*weight*mp.exp(-a*(length-t)-b*(Eend-Et))*G**b
                value+=atom
                derivative+=b*LZ*atom
            result.extend((value,derivative))
        return [result[0],result[2],result[1],result[3]]

    @lru_cache(maxsize=256)
    def _baseline(self,y_key,z_key):
        with mp.workdps(self.precision):
            y=_mp(y_key);z=_mp(z_key)
            G,LZ=self._ratio_jet(mp.mpf(0),z)
            state=[mp.mpf(5)/8*G,mp.mpf(5)/6*G**2,
                   mp.mpf(5)/8*G*LZ,mp.mpf(5)/3*G**2*LZ]
            if y<=0:return tuple(state)
            bounds=self.schedule._make_stage_bounds()
            for name,slope in self.stages:
                left,right=map(lambda v:_mp(str(v)),bounds[name])
                if y<=left:break
                end=min(y,right)
                state=self._transfer(state,left,end,z,slope)
                if end==y:break
            return tuple(state)

    @lru_cache(maxsize=256)
    def angular_difference_from_axis(self,y_key,z_key):
        """Propagate Xtheta(Z)-Xtheta(0) without subtracting two O(1) means.

        This retained difference is needed after imposing the coherent Z=0
        preheat waiting equation; it is not a replacement for the actual mean.
        """
        with mp.workdps(self.precision):
            y=_mp(y_key);z=_mp(z_key);logq=mp.log1p(z*z)
            def difference(at):
                v=_mp(str(self.schedule.y_v));f=_mp(str(self.schedule.y_f))
                if at<=v:flatten=mp.mpf(0)
                elif at>=f:flatten=mp.mpf(1)
                else:
                    from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
                    flatten=ContinuousAxialPulse.sigma_pair((at-v)/(f-v))[0]
                return mp.expm1(-(1-flatten)*logq)
            state=mp.mpf(5)/8*difference(mp.mpf(0))
            bounds=self.schedule._make_stage_bounds()
            for name,slope in self.stages:
                left,right=(_mp(str(v)) for v in bounds[name])
                if y<=left:break
                end=min(y,right);length=end-left
                if slope is not None:
                    rate_s=(-mp.mpf('.5')-self.mu if slope=='pulse' else
                            -(1+self.delta)/2 if slope=='waiting' else mp.mpf(slope))
                    rate=mp.mpf('1.5')+rate_s
                    atom=-mp.expm1(-rate*length)/rate if rate else length
                    state=state*mp.exp(-rate*length)+difference(left)*atom
                else:
                    El=self._base(mp.nstr(left,self.precision))
                    Ee=self._base(mp.nstr(end,self.precision))
                    state*=mp.exp(-mp.mpf('1.5')*length-(Ee-El))
                    for unit,weight in self.nodes:
                        t=length*unit;at=left+t
                        Et=self._base(mp.nstr(at,self.precision))
                        state+=length*weight*mp.exp(-mp.mpf('1.5')*(length-t)-(Ee-Et))*difference(at)
                if end==y:break
            return state

    def _bump_moments(self,logR,z):
        """Separate normalized corrections and their Z derivatives."""
        t=_mp(str(self.profile.offset(logR,self.schedule.logR_rel)))
        if t<=mp.mpf('-3.15'):
            return [mp.mpf(0)]*4
        coeff=self.correction.coefficients(z)
        lam_theta=1-self.mu;lam_energy=-2*self.mu
        theta=energy=thetaZ=energyZ=mp.mpf(0)
        for key,center in (('d1',-3),('d2',-1)):
            d=_mp(coeff[key]['arbitrary_exponent_value'])
            dz=_mp(coeff[key+'_Z']['arbitrary_exponent_value'])
            theta_atom=mp.exp(lam_theta*center)*self.correction.weighted_atom(
                lam_theta,upper=t-center)
            linear=mp.exp(lam_energy*center)*self.correction.weighted_atom(
                lam_energy,upper=t-center)
            square=mp.exp(lam_energy*center)*self.correction.weighted_atom(
                lam_energy,power=2,upper=t-center)
            theta+=d*theta_atom;thetaZ+=dz*theta_atom
            energy+=2*d*linear+d*d*square
            energyZ+=2*dz*linear+2*d*dz*square
        Erel=self._base(str(self.schedule.y_rel))
        y=self._preheat_y(logR)
        E=self._base(mp.nstr(y,self.precision))
        theta_factor=mp.exp(-mp.mpf('1.5')*t-(E-Erel))
        energy_factor=mp.exp(-t-2*(E-Erel))
        return [theta*theta_factor,energy*energy_factor,
                thetaZ*theta_factor,energyZ*energy_factor]

    @lru_cache(maxsize=128)
    def _inner_offsets(self,z_key):
        """Both primitives use the actual corrected inner field at Rh."""
        with mp.workdps(self.precision):
            source=self.profile.seed_source
            z=_mp(z_key)
            value=source.inner.evaluate_x(mp.e,z)
            radius=mp.nstr(source.logRh,self.precision)
            y=_mp(str(self.profile.offset(radius,self.schedule.logRref)))
            E=self._base(mp.nstr(y,self.precision))
            base=self._baseline(mp.nstr(y,self.precision),z_key)
            theta_scale=mp.sqrt(2)*mp.exp(mp.mpf('1.5')*source.logRh+E)
            energy_scale=mp.exp(source.logRh+2*E)
            raw=value['raw_quadratic_integrals']
            return [value['moments']['theta']-theta_scale*base[0],
                    2*raw['swirl']-energy_scale*base[1],
                    value['momentsZ']['theta']-theta_scale*base[2],
                    2*raw['swirl_Z']-energy_scale*base[3]]

    def moments_jet(self,logR,Z,*,include_inner=True):
        with mp.workdps(self.precision):
            if self.profile.offset(logR,self.schedule.logR_tail)>0:
                return self.heat_provider.moments_jet(logR,Z,include_inner=include_inner)
            z=_mp(Z)
            if abs(z)>=1:raise ValueError('Actual inner seed requires |Z|<1')
            y=self._preheat_y(logR)
            z_key=mp.nstr(z,self.precision)
            baseline=self._baseline(mp.nstr(y,self.precision),z_key)
            correction=self._bump_moments(logR,z)
            E=self._base(mp.nstr(y,self.precision))
            rlog=_mp(logR)
            scales=[mp.sqrt(2)*mp.exp(mp.mpf('1.5')*rlog+E),
                    mp.exp(rlog+2*E)]
            offsets=self._inner_offsets(z_key) if include_inner else [mp.mpf(0)]*4
            boundary=_mp(mp.nstr(self.profile.seed_source.logRh,self.precision))
            if include_inner and rlog<boundary:
                raise ValueError('Actual seeded outer primitives start at Rh')
            values=[scales[i%2]*(baseline[i]+correction[i])+offsets[i] for i in range(4)]
            return dict(theta=values[0],swirl_energy=values[1],theta_Z=values[2],
                swirl_energy_Z=values[3],baseline_normalized=list(baseline),
                bump_normalized=correction,inner_offsets=offsets,
                inner_offsets_reapplied_once=include_inner,
                inner_seed_uncertainty_enclosed=False,
                bump_atoms_shared=True,quadrature_order=self.order,
                quadrature_error_enclosed=False,complete_five_moments=False,
                heat_supported=False,finite_energy_certified=False,
                scale_recursion_certified=False)


def run():
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    from lei_ren_part1_paper_joined_outer import build_joined_field
    print('building shared continuous field',flush=True)
    source=build_joined_field()
    folder=Path(__file__).parent
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    engine=field.outer.angular_moment_provider
    with mp.workdps(field.precision):
        z=mp.mpf('.3');checks=[]
        for label,radius in (('Rh',mp.nstr(source.logRh,field.precision)),
                              ('flatten_mid',field.outer.log_at(source.schedule.logR_v,'50')),
                              ('first_bump',field.outer.log_at(source.schedule.logR_rel,'-2.96')),
                              ('Rtail',str(source.schedule.logR_tail))):
            row=engine.moments_jet(radius,z)
            if label=='Rh':
                inner=source.inner.evaluate_x(mp.e,z)
                seed_errors=[abs(row['theta']/inner['moments']['theta']-1),
                        abs(row['swirl_energy']/(2*inner['raw_quadratic_integrals']['swirl'])-1),
                        abs(row['theta_Z']/inner['momentsZ']['theta']-1),
                        abs(row['swirl_energy_Z']/(2*inner['raw_quadratic_integrals']['swirl_Z'])-1)]
                if max(seed_errors)>mp.mpf('1e-70'):raise ArithmeticError('Actual inner moment seed mismatch')
            checks.append(dict(label=label,theta=signed_log(row['theta'],field.precision),
                theta_Z=signed_log(row['theta_Z'],field.precision),
                swirl_energy=signed_log(row['swirl_energy'],field.precision),
                swirl_energy_Z=signed_log(row['swirl_energy_Z'],field.precision),
                bump_normalized=[signed_log(v,field.precision) for v in row['bump_normalized']]))
        # A local derivative checks the actual velocity integrands in the
        # flattening region. Divide out large scale factors before comparing.
        center=field.outer.log_at(source.schedule.logR_v,'50');h=mp.mpf('1e-4')
        neighbors={i:engine.moments_jet(field.outer.log_at(center,i*h),z,include_inner=False)
                   for i in (-2,-1,1,2)}
        row=engine.moments_jet(center,z,include_inner=False)
        y=_mp(str(field.outer.offset(center,source.schedule.logRref)))
        E=engine._base(mp.nstr(y,field.precision));Rlog=_mp(str(center))
        velocity=field.outer.values_with_jets(center,z)
        targets=[mp.sqrt(2)*mp.exp(mp.mpf('1.5')*Rlog)*velocity['Utheta'],
                 mp.exp(Rlog)*velocity['Utheta']**2,
                 mp.sqrt(2)*mp.exp(mp.mpf('1.5')*Rlog)*velocity['Utheta_Z'],
                 2*mp.exp(Rlog)*velocity['Utheta']*velocity['Utheta_Z']]
        errors=[]
        for key,target in zip(('theta','swirl_energy','theta_Z','swirl_energy_Z'),targets):
            derivative=sum(w*neighbors[i][key] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*h)
            errors.append(mp.nstr(abs(derivative/target-1),40))
        if max(map(mp.mpf,errors))>mp.mpf('1e-8'):
            raise ArithmeticError('Angular moment derivative differs from actual velocity integrand')
        # Test the separate tiny correction primitive. Testing total moments
        # here would hide it under the much larger baseline.
        center=field.outer.log_at(source.schedule.logR_rel,'-2.96');h=mp.mpf('1e-5')
        bump=engine._bump_moments(center,z)
        neighbors={i:engine._bump_moments(field.outer.log_at(center,i*h),z) for i in (-2,-1,1,2)}
        point=engine.correction.value_jet(mp.mpf('-2.96'),z)
        targets=[point['h'],2*point['h']+point['h']**2,
                 point['h_Z'],2*point['h_Z']*(1+point['h'])]
        bump_errors=[]
        for j,target in enumerate(targets):
            derivative=sum(w*neighbors[i][j] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*h)
            rate=1-engine.mu if j%2==0 else -2*engine.mu
            bump_errors.append(mp.nstr(abs((derivative+rate*bump[j])/target-1),40))
        if max(map(mp.mpf,bump_errors))>mp.mpf('1e-10'):
            raise ArithmeticError('Tiny bump primitive differs from the point velocity correction')
        report=dict(samples=checks,Rh_seed_relative_errors=[mp.nstr(v,40) for v in seed_errors],
            flattening_integrand_relative_errors=errors,
            separate_bump_integrand_relative_errors=bump_errors,
            actual_inner_seed_installed=True,continuous_bump_atoms_shared=True,
            quadrature_error_enclosed=False,complete_five_moments=False,
            heat_supported=False,finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'flattening_integrand_relative_errors':errors,'samples':len(checks)}),flush=True)
    return report


if __name__=='__main__':run()
