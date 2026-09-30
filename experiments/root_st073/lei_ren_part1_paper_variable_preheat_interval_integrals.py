"""Convergent interval subdivision for variable preheat radial stages.

Uses a single monotone switch-primitive grid, interval log-amplitudes, and
positive interval density integration. Scope is stored schedule parameters.
"""
import mpmath as mp
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


class VariablePreheatIntervalIntegrals:
    def __init__(self,enclosures):
        self.e=enclosures;self.iv=enclosures.iv;self.schedule=enclosures.schedule

    def _primitive(self,argument,grid,panels):
        iv=self.iv
        def at(x):
            if x<=0:return iv.mpf(0)
            if x>=1:return iv.mpf(x)-iv.mpf('.5')
            # The index is a convenience only. The directed Lipschitz shift
            # below remains valid even if scalar index rounding picks a
            # neighboring node.
            k=max(0,min(panels,int(x*panels)))
            shift=iv.mpf(x)-iv.mpf(k)/panels
            lo,hi=endpoints(shift)
            correction=iv.mpf([min(mp.mpf(0),lo),max(mp.mpf(0),hi)])
            value=grid[k]+correction
            vl,vh=endpoints(value)
            return iv.mpf([max(mp.mpf(0),vl),min(mp.mpf('.5'),vh)])
        lo,hi=endpoints(argument)
        lower=at(lo);upper=at(hi)
        return iv.mpf([endpoints(lower)[0],endpoints(upper)[1]])

    def _log_amplitude(self,y,grid,panels):
        iv=self.iv;s=self.schedule;c=self.e.scalar
        mu=c(s.mu);dt=c(s.delta)
        Js=[self._primitive(arg,grid,panels) for arg in
            (y,y-c(s.y_d),y-c(s.y_rel),y-c(s.y_rel)-1-c(s.Ts))]
        return iv.mpf('.1')*y-iv.mpf('.6')*Js[0]-mu*Js[1]-(1-mu)*Js[2]+(1-dt/2)*Js[3]

    def integrate_stage(self,stage,*,panels=128):
        if isinstance(panels,bool) or not isinstance(panels,int) or panels<2:
            raise ValueError('Require integer panels>=2')
        allowed=('slope_transition_ref','slope_transition_mu','steep_transition_in','steep_transition_out','z_flatten')
        if stage not in allowed:raise ValueError('Only variable preheat stages supported')
        # Reuse the paper parameter gate and stage ordering contract.
        self.e.stage_pressure_upper(stage)
        iv=self.iv;s=self.schedule;c=self.e.scalar
        with mp.workdps(self.e.precision+30):
            grid=self.e.primitive_grid_enclosures(panels)
            left,right=s._stage_bounds[stage]
            L=c(right)-c(left);total=iv.mpf(0)
            for i in range(panels):
                u=iv.mpf([endpoints(iv.mpf(i)/panels)[0],endpoints(iv.mpf(i+1)/panels)[1]])
                y=c(left)+L*u
                ell=self._log_amplitude(y,grid,panels)
                if stage=='z_flatten':
                    theta=1-self.e.sigma_interval((y-c(s.y_v))/c(s.Tf))
                    factor=iv.exp((2*theta-3)*iv.ln(2))
                else:factor=iv.mpf('.5' if stage in ('slope_transition_ref','slope_transition_mu') else '.125')
                total+=L/panels*factor*iv.exp(2*ell)
            lo,hi=endpoints(total)
            return dict(stage=stage,panels=panels,normalized_mass_at_Z0_interval=total,
                interval_width=hi-lo,relative_width=(hi-lo)/lo if lo>0 else mp.inf,
                primitive_grid_panels=panels,directed_interval_arithmetic=True,
                scope='positive radial integral at Z=0 for stored schedule parameters',
                original_parameter_errors_enclosed=False,
                uniform_axial_derivative_quadrature_errors_enclosed=False,
                five_defect_interval_closure=False)

    def integrate_negative_stage(self,stage):
        """Analytic exponential lower/upper integrals from stage slope bounds.

        Only stages with a constant angular density prefactor are accepted.
        This avoids subdivision of the extremely long negative-slope stages.
        Bounds remain valid when the slope is not exactly constant.
        """
        allowed=('axial_turnoff','power_buffer','pulse_reserved','power_buffer_rel',
                 'steep_power','waiting')
        if stage not in allowed:raise ValueError('Require a supported negative-slope stage')
        iv=self.iv
        record=self.e.stage_pressure_upper(stage)
        L=record['length'];ell=record['left_log_amplitude'];slo,shi=record['slope_bounds']
        if endpoints(shi)[1]>=0:raise ValueError('Negative upper slope required')
        factor=iv.mpf('.5' if stage in ('axial_turnoff','power_buffer','pulse_reserved') else '.125')
        if self.schedule._stage_bounds[stage][0]==self.schedule._stage_bounds[stage][1]:
            mass=iv.mpf(0)
        else:
            low=factor*iv.exp(2*ell)*(1-iv.exp(2*slo*L))/(-2*slo)
            high=factor*iv.exp(2*ell)*(1-iv.exp(2*shi*L))/(-2*shi)
            mass=iv.mpf([max(mp.mpf(0),endpoints(low)[0]),endpoints(high)[1]])
        return dict(stage=stage,normalized_mass_at_Z0_interval=mass,
            directed_interval_arithmetic=True,radial_quadrature_used=False,
            source_scope='stored parameters and declared negative slope envelope',
            uniform_axial_derivative_quadrature_errors_enclosed=False)


def finite_mass_error_interval(integral_result,finite_mass,enclosures):
    """Signed true-minus-finite mass interval, not an entire field error."""
    return integral_result['normalized_mass_at_Z0_interval']-enclosures.iv.mpf(finite_mass)
