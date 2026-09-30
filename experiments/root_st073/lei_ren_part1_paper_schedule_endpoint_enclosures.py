"""Directed interval endpoint and positive preheat stage bounds.

Scope: the continuous formula with the supplied stored Decimal parameters.
Does not certify how those parameters were calculated from transcendental
paper inputs. The default endpoint enclosures use no primitive quadrature;
``primitive_refined`` and ``primitive_grid_enclosures`` provide optional
directed Darboux enclosures for the monotone switch primitive.
"""
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


def endpoints(value):
    return tuple(mp.make_mpf(v) for v in value._mpi_)


class ScheduleEndpointEnclosures:
    def __init__(self,schedule,*,precision=None):
        self.schedule=schedule
        self.precision=max(schedule.decimal_precision+30,precision or 0)
        self.iv=MPIntervalContext();self.iv.dps=self.precision
        # Keys contain only exact scalar endpoint tuples and the panel count.
        # Interval objects themselves are deliberately not used as keys.
        self._primitive_scalar_cache={}

    def scalar(self,value):return self.iv.mpf(str(value))

    def _as_interval(self,value):
        return value if hasattr(value,'_mpi_') else self.scalar(value)

    @staticmethod
    def _panel_count(panels):
        try:n=int(panels)
        except (TypeError,ValueError,OverflowError):raise ValueError('panels must be a positive integer')
        if n<1 or n!=panels:raise ValueError('panels must be a positive integer')
        return n

    def _sigma_point(self,x):
        """Evaluate the known increasing logistic switch at one scalar endpoint."""
        iv=self.iv;one=iv.mpf(1)
        if x<=0:return iv.mpf(0)
        if x>=1:return one
        xv=iv.mpf(x)
        phase=one/(xv*xv)-one/((one-xv)*(one-xv))
        plo,phi=endpoints(phase)
        # Use the reciprocal logistic branch on the positive side so that
        # tiny switches do not require forming exp(large positive phase).
        if plo>=0:
            tiny=iv.exp(-phase)
            return tiny/(one+tiny)
        if phi<=0:
            tiny=iv.exp(phase)
            return one/(one+tiny)
        return one/(one+iv.exp(phase))

    def sigma_interval(self,value):
        """Return a directed hull for the increasing switch ``sigma``.

        The phase is ``1/x**2 - 1/(1-x)**2`` on ``0 < x < 1``.  Exact
        endpoint values are used outside that interval, and the monotonic
        endpoint hull handles a genuine input interval without sampling.
        """
        value=self._as_interval(value);lo,hi=endpoints(value)
        if lo>hi:raise ValueError('invalid interval with lower endpoint above upper endpoint')
        if not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('sigma interval requires finite endpoints')
        if hi<=0:return self.iv.mpf(0)
        if lo>=1:return self.iv.mpf(1)
        left=mp.mpf(0) if lo<0 else lo
        right=mp.mpf(1) if hi>1 else hi
        lower=self._sigma_point(left);upper=self._sigma_point(right)
        slo,_=endpoints(lower);_,shi=endpoints(upper)
        return self.iv.mpf([slo,shi])

    def _primitive_darboux_scalar(self,x,panels):
        """Darboux enclosure for J(x) with scalar endpoint x in [0,1]."""
        iv=self.iv;n=self._panel_count(panels)
        if x<=0:return iv.mpf(0)
        if x>=1:return iv.mpf('.5')
        xv=iv.mpf(x)
        h=xv/iv.mpf(n)
        lower_sum=iv.mpf(0);upper_sum=iv.mpf(0)
        for k in range(n):
            left=h*k;right=h*(k+1)
            lower_sum += self.sigma_interval(left)
            upper_sum += self.sigma_interval(right)
        lower=h*lower_sum;upper=h*upper_sum
        return iv.mpf([endpoints(lower)[0],endpoints(upper)[1]])

    def _primitive_scalar(self,x,panels):
        """Cached scalar-endpoint primitive enclosure; interval callers do not cache."""
        iv=self.iv;n=self._panel_count(panels);xv=iv.mpf(x)
        xlo,xhi=endpoints(xv);key=(xlo._mpf_,xhi._mpf_,n)
        cached=self._primitive_scalar_cache.get(key)
        if cached is not None:return cached
        if x<=0:
            result=iv.mpf(0)
        elif x>=1:
            result=xv-iv.mpf('.5')
        else:
            result=self._primitive_darboux_scalar(x,n)
        self._primitive_scalar_cache[key]=result
        return result

    def primitive_refined(self,value,panels=128):
        """Directed Darboux enclosure of monotone J at a scalar or interval."""
        n=self._panel_count(panels);value=self._as_interval(value)
        lo,hi=endpoints(value)
        if lo>hi:raise ValueError('invalid interval with lower endpoint above upper endpoint')
        if not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('primitive interval requires finite endpoints')
        if hi<=0:return self.iv.mpf(0)
        if lo>=1:return self.iv.mpf(value)-self.iv.mpf('.5')
        if lo<=0:lower=self.iv.mpf(0)
        else:lower=self._primitive_scalar(lo,n)
        if hi>=1:upper=self.iv.mpf(value)-self.iv.mpf('.5')
        else:upper=self._primitive_scalar(hi,n)
        lower_endpoint=endpoints(lower)[0];upper_endpoint=endpoints(upper)[1]
        return self.iv.mpf([lower_endpoint,upper_endpoint])

    def primitive_grid_enclosures(self,panels=128):
        """Return O(N) directed prefix enclosures for J(i/N), i=0,...,N.

        A single uniform sigma grid supplies lower and upper Darboux prefix
        sums.  The terminal prefix is replaced by the exact symmetry value
        J(1)=1/2.
        """
        iv=self.iv;n=self._panel_count(panels)
        h=iv.mpf(1)/iv.mpf(n)
        grid=[h*i for i in range(n+1)]
        sigma=[self.sigma_interval(x) for x in grid]
        result=[iv.mpf(0)]
        lower_sum=iv.mpf(0);upper_sum=iv.mpf(0)
        for i in range(1,n+1):
            lower_sum += sigma[i-1]
            upper_sum += sigma[i]
            lower=h*lower_sum;upper=h*upper_sum
            if i==n:
                result.append(iv.mpf('.5'))
            else:
                result.append(iv.mpf([endpoints(lower)[0],endpoints(upper)[1]]))
        return result

    def primitive(self,value):
        """Monotone J interval, exact outer branches and conservative interior.

        sigma(x)+sigma(1-x)=1 gives integral_0^1 sigma=1/2.
        For a partial argument, 0<=J(x)<=x; no quadrature is used.
        """
        lo,hi=endpoints(value);iv=self.iv
        def bound(x):
            v=iv.mpf(x)
            if x<=0:return iv.mpf(0),iv.mpf(0)
            if x>=1:return v-iv.mpf('0.5'),v-iv.mpf('0.5')
            # Lipschitz continuity from J(0)=0 and J(1)=1/2.
            return (iv.mpf(0) if x<=mp.mpf('.5') else v-iv.mpf('.5')),(v if x<=mp.mpf('.5') else iv.mpf('.5'))
        lower,_=bound(lo);_,upper=bound(hi)
        return iv.mpf([endpoints(lower)[0],endpoints(upper)[1]])

    def log_amplitude_ratio(self,y):
        iv=self.iv;s=self.schedule;y=self.scalar(y);mu=self.scalar(s.mu);delta=self.scalar(s.delta)
        arguments=[y,y-self.scalar(s.y_d),y-self.scalar(s.y_rel),
                   y-self.scalar(s.y_rel)-1-self.scalar(s.Ts)]
        Js=[self.primitive(v) for v in arguments]
        value=iv.mpf('.1')*y-iv.mpf('.6')*Js[0]-mu*Js[1]-(1-mu)*Js[2]+(1-delta/2)*Js[3]
        return dict(interval=value,primitive_arguments=arguments,
                    primitive_intervals=Js,primitive_quadrature_used=False)

    def stage_pressure_upper(self,stage,*,axial_radius='.8'):
        """Positive preheat mass/derivative upper bounds before heat collar.

        Prefactors bound the normalized density. Constant slopes use an
        infinite exponential integral cap to avoid huge stage-length bounds.
        These are conservative integral bounds, not Gauss error estimates.
        """
        iv=self.iv;s=self.schedule;mu=self.scalar(s.mu);dt=self.scalar(s.delta)
        slopes={
            'slope_transition_ref':(iv.mpf('-.5'),iv.mpf('.1')),
            'axial_turnoff':(iv.mpf('-.5'),iv.mpf('-.5')),
            'slope_transition_mu':(-iv.mpf('.5')-mu,-iv.mpf('.5')),
            'power_buffer':(-iv.mpf('.5')-mu,-iv.mpf('.5')),
            'pulse_reserved':(-iv.mpf('.5')-mu,-iv.mpf('.5')-mu),
            'z_flatten':(-iv.mpf('.5')-mu,-iv.mpf('.5')-mu),
            'power_buffer_rel':(-iv.mpf('.5')-mu,-iv.mpf('.5')-mu),
            'steep_transition_in':(-iv.mpf('1.5'),-iv.mpf('.5')-mu),
            # Rounded stored y_s/y_q can straddle a switching boundary.
            # Use the negative envelope rather than assert an exact slope.
            'steep_power':(-iv.mpf('1.5'),-iv.mpf('.5')-mu),
            'steep_transition_out':(-iv.mpf('1.5'),-(1+dt)/2),
            'waiting':(-iv.mpf('1.5'),-(1+dt)/2)}
        if stage not in slopes:raise ValueError('Only finite pre-collar stages supported')
        if not 0<endpoints(mu)[0]<=endpoints(mu)[1]<=endpoints(iv.mpf(1)/60)[0]:
            raise ValueError('Require the paper gate 0<mu<=1/60')
        if not 0<endpoints(dt)[0]<=endpoints(dt)[1]<1:raise ValueError('Require 0<delta<1')
        left,right=s._stage_bounds[stage];L=iv.mpf(0) if left==right else self.scalar(right)-self.scalar(left)
        if endpoints(L)[0]<0:raise ValueError('Unresolved or negative stage length')
        ell=self.log_amplitude_ratio(left);slo,shi=slopes[stage]
        # For negative slopes the left endpoint is the amplitude supremum.
        # For a sign-changing stage, ell(left)+positive slope*L is valid.
        slope_upper=endpoints(shi)[1]
        sup=ell['interval']+iv.mpf(max(mp.mpf(0),slope_upper))*L
        factor=iv.mpf('.5') if stage in ('slope_transition_ref','axial_turnoff','slope_transition_mu','power_buffer','pulse_reserved','z_flatten') else iv.mpf('.125')
        mass=factor*L*iv.exp(2*sup)
        cap=None
        if slope_upper<0:
            cap=factor*iv.exp(2*ell['interval'])/(-2*shi)
        upper=endpoints(mass)[1]
        if cap is not None:upper=min(upper,endpoints(cap)[1])
        a=self.scalar(axial_radius)
        if not 0<=endpoints(a)[0]<=endpoints(a)[1]<1:raise ValueError('Require radius<1')
        beta=iv.mpf(2 if stage in ('slope_transition_ref','axial_turnoff','slope_transition_mu','power_buffer','pulse_reserved','z_flatten') else 0)
        factors=[iv.mpf(1),2*beta*a,2*beta+4*beta*(beta+1)*a*a,
                 12*beta*(beta+1)*a+8*beta*(beta+1)*(beta+2)*a**3]
        bounds=[endpoints(iv.mpf(upper)*f)[1] for f in factors]
        return dict(stage=stage,length=L,left_log_amplitude=ell['interval'],slope_bounds=(slo,shi),
            mass_upper=upper,normalized_derivative_upper_bounds=bounds,
            constant_or_negative_slope_cap_used=cap is not None,
            relative_to_stored_schedule_parameters=True,endpoint_primitive_quadrature_used=False,
            directed_interval_arithmetic=True,transcendental_parameter_errors_enclosed=False,
            heat_collar_included=False,core_propagation_errors_enclosed=False)

    def terminal_preheat_bounds(self):
        """Enclose H=1 heat collar and exact power tail of the shared datum.

        Uses the stored shared log_c_inf, not a substituted normalization.
        K0 lies in [1-epsilon,1], is independent of Z, and equals1 for t>=3.
        The independent endpoint normalization route is only a diagnostic.
        """
        iv=self.iv;s=self.schedule;eps=self.scalar(s.epsilon);lam=1+self.scalar(s.delta)
        if not 0<endpoints(eps)[0]<=endpoints(eps)[1]<1:raise ValueError('Require 0<epsilon<1')
        if s.y_b-s.y_tail!=3:raise ValueError('Require exact stored heat collar length3')
        # Preserve the actual stored source heat scale through cancellation
        # of large logarithms with enough directed interval precision.
        logscale=2*(self.scalar(s._log_c_inf)-self.scalar(s.logPstar))-iv.ln(2)-lam*self.scalar(s.logR_tail)
        scale=iv.exp(logscale);tail_factor=iv.exp(-3*lam)/lam
        collar_factor=(1-iv.exp(-3*lam))/lam
        lower=endpoints(scale*(1-eps)**2*collar_factor)[0]
        upper=endpoints(scale*collar_factor)[1]
        collar=iv.mpf([lower,upper]);exterior=scale*tail_factor
        combined=collar+exterior
        ell=self.log_amplitude_ratio(s.y_tail)['interval']
        expected_logc=self.scalar(s.logPstar)+ell+lam*self.scalar(s.logR_tail)/2-iv.ln(2*(1-eps))
        return dict(stored_heat_log_scale=logscale,stored_heat_scale=scale,
            heat_collar_mass_interval=collar,exterior_mass_interval=exterior,
            combined_terminal_mass_interval=combined,
            endpoint_normalization_log_residual=self.scalar(s._log_c_inf)-expected_logc,
            normalization_reset=False,collar_K0_bounds=iv.mpf([endpoints(1-eps)[0],1]),
            normalized_derivative_upper_bounds=[endpoints(combined)[1],mp.mpf(0),mp.mpf(0),mp.mpf(0)],
            H_replaced_by_one=True,Z_independent=True,
            relative_to_stored_schedule_parameters=True,directed_interval_arithmetic=True,
            endpoint_primitive_quadrature_used=False,terminal_integral_quadrature_used=False,
            actual_heat_velocity_field_certified=False,finite_energy_certified=False,
            transcendental_parameter_errors_enclosed=False)
