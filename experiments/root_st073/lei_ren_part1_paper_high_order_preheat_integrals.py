"""Taylor-remainder interval quadrature for the first two preheat transitions."""
import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,integrate_symmetric
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_variable_preheat_interval_integrals import VariablePreheatIntervalIntegrals


def switch_taylor(ctx,value,order):
    """Derivative coefficients of sigma at a point or entire interior cell."""
    x=IntervalTaylor.variable(ctx,value,order)
    phase=x**(-2)-(1-x)**(-2)
    lo,hi=endpoints(phase.coefficients[0])
    if lo>=0:
        tiny=(-phase).exp();return tiny/(1+tiny)
    return 1/(1+phase.exp())


class HighOrderPreheatIntegrals:
    def __init__(self,enclosures):
        self.e=enclosures;self.iv=enclosures.iv;self._cache={}

    def primitive_grid(self,panels=128,order=12):
        """Integrate sigma with analytic edge bounds and interior Taylor panels."""
        if not isinstance(panels,int) or panels<4 or not isinstance(order,int) or order<2:
            raise ValueError('Require panels>=4 and order>=2')
        key=(panels,order)
        if key in self._cache:return self._cache[key]
        iv=self.iv;h=iv.mpf(1)/panels;r=h/2
        edge=h*self.e.sigma_interval(h)
        total=iv.mpf([0,endpoints(edge)[1]]);grid=[iv.mpf(0),total]
        for i in range(1,panels-1):
            center=(iv.mpf(i)+iv.mpf('.5'))/panels
            cell=iv.mpf([endpoints(iv.mpf(i)/panels)[0],endpoints(iv.mpf(i+1)/panels)[1]])
            value=integrate_symmetric(switch_taylor(iv,center,order),switch_taylor(iv,cell,order),r)
            lo,hi=endpoints(value)
            value=iv.mpf([max(mp.mpf(0),lo),min(endpoints(h)[1],hi)])
            total+=value;grid.append(total)
        # J(1)=1/2 is exact; intermediate prefixes retain all panel errors.
        grid.append(iv.mpf('.5'));self._cache[key]=grid
        return grid

    def _node_primitive(self,argument,index,grid,panels):
        iv=self.iv;shift=argument-iv.mpf(index)/panels
        lo,hi=endpoints(shift)
        value=grid[index]+iv.mpf([min(mp.mpf(0),lo),max(mp.mpf(0),hi)])
        vl,vh=endpoints(value)
        return iv.mpf([max(mp.mpf(0),vl),min(mp.mpf('.5'),vh)])

    def integrate_stage(self,stage,*,panels=64,order=12):
        if stage not in ('slope_transition_ref','slope_transition_mu'):
            raise ValueError('Only first two constant-beta transitions supported')
        iv=self.iv;s=self.e.schedule;c=self.e.scalar
        self.e.stage_pressure_upper(stage)
        with mp.workdps(self.e.precision+40):
            grid=self.primitive_grid(2*panels,order)
            left,right=s._stage_bounds[stage];L=c(right)-c(left)
            if stage=='slope_transition_ref':base=iv.mpf(0);slope=iv.mpf('.1');coefficient=iv.mpf('.6')
            else:base=self.e.log_amplitude_ratio(left)['interval'];slope=iv.mpf('-.5');coefficient=c(s.mu)
            total=iv.mpf(0);r=iv.mpf(1)/(2*panels)
            helper=VariablePreheatIntervalIntegrals(self.e)
            for i in range(panels):
                uc=(iv.mpf(i)+iv.mpf('.5'))/panels;xc=L*uc
                cell=iv.mpf([endpoints(L*iv.mpf(i)/panels)[0],endpoints(L*iv.mpf(i+1)/panels)[1]])
                Jc=self._node_primitive(xc,2*i+1,grid,2*panels)
                Jcell=helper._primitive(cell,grid,2*panels)
                def density(value,J):
                    # Coefficients are with respect to normalized panel u.
                    sigma=switch_taylor(iv,value,order-1)
                    coefficients=[base+slope*value-coefficient*J]
                    coefficients.extend((-coefficient*sigma.coefficients[k-1]*L**k/k) for k in range(1,order+1))
                    coefficients[1]+=slope*L
                    return (2*IntervalTaylor(iv,coefficients)).exp()*iv.mpf('.5')*L
                if i==0 or i==panels-1:
                    # Derivatives of sigma at the flat endpoints are not
                    # evaluated by singular formulas. On these tiny edges,
                    # compare to the exact outer primitive branch instead.
                    yl=L*iv.mpf(i)/panels;yr=L*iv.mpf(i+1)/panels
                    y=iv.mpf([endpoints(yl)[0],endpoints(yr)[1]])
                    if i==0:
                        correction_upper=endpoints(y*self.e.sigma_interval(y))[1]
                        edge_base=base;edge_slope=slope
                    else:
                        d=1-y;small=self.e.sigma_interval(d)
                        correction_upper=max(mp.mpf(0),endpoints(d*small)[1])
                        edge_base=base+coefficient/2;edge_slope=slope-coefficient
                    edge=iv.mpf('.5')*iv.exp(2*edge_base)*(iv.exp(2*edge_slope*yr)-iv.exp(2*edge_slope*yl))/(2*edge_slope)
                    multiplier=iv.mpf([endpoints(iv.exp(-2*coefficient*iv.mpf(correction_upper)))[0],1])
                    total+=edge*multiplier
                else:
                    total+=integrate_symmetric(density(xc,Jc),density(cell,Jcell),r)
            lo,hi=endpoints(total)
            return dict(stage=stage,panels=panels,Taylor_order=order,
                normalized_mass_at_Z0_interval=total,interval_width=hi-lo,
                primitive_method='Taylor remainder with analytic flat edge bounds',
                radial_edge_method='exact exponential integral plus positive flat correction bound',
                source_scope='stored schedule parameters',directed_interval_arithmetic=True,
                full_pressure_error_enclosed=False,five_defect_interval_closure=False)
