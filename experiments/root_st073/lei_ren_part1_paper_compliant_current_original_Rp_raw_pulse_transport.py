"""Six current active pulse histories and the common closed exterior caller.

The exact radial decay is separated from directed axial coefficients. Active
axial/linear histories and independent P0 survive; no post-Rv zero theorem is
applied inside the pulse. The full 15-chart caller selects closed heat output.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_same_repair_heat_closure as closed
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

raw=closed.raw_source;radius=raw.radius
HERE,PREFIX,sha=closed.HERE,closed.PREFIX,closed.sha
NAME=PREFIX+'current_original_Rp_raw_pulse_transport.json.gz'
RECEIPT=PREFIX+'current_original_Rp_raw_pulse_transport_check.json'
GATES=('current_original_Rp_six_pulse_absolute_scale_functions_installed',
    'current_original_Rp_six_pulse_factorized_raw_histories_callable',
    'current_original_Rp_six_pulse_first_native_history_derivatives_installed',
    'current_original_Rp_fifteen_chart_raw_caller_uses_closed_heat')
OPEN=closed.OPEN
PULSE=tuple(radius.post.selected.PULSE_CHARTS)
CHARTS=PULSE+raw.CHARTS
UNITS=('R','pulse_F','Pstar')
POWERS=dict(Mz=(1,1,1),Mtheta=(Fraction(3,2),1,1),
    Mtheta_z=(Fraction(3,2),2,2),Mztheta=(1,2,2),Mp=(0,0,2),
    P0=(0,0,2),pressure=(0,0,2),Utheta=(0,1,1),Uz=(0,1,1),Ur=(Fraction(1,2),1,1))
VIEWS=(('entrance','pulse_entrance','7/2'),('main','pulse_main',1),
    ('exit','pulse_exit','21/2'),('gap','pulse_gap','23/2'),
    ('gap_end','pulse_gap_end',-5),('active_end','pulse_end',-3),
    ('flatten','flatten',37),('power','outer_power','1/2'),
    ('angular','outer_angular',-3),('entry','steep_entry','1/2'),
    ('steep','steep_power','1/2'),('exit_steep','steep_exit','1/2'),
    ('waiting','waiting','1/2'),('collar','heat_collar',1),('exterior','heat_exterior',4))


@dataclass(frozen=True)
class FactorizedPulseSourceTaylor:
    """R^r Fpulse^e Pstar^p times directed ordinary axial Taylor rows.

    Fpulse=exp(-(.5+mu)*log(R/Rp)) is an exact Z-independent
    function. Fourth axial order is retained for recovered radial velocity.
    """
    name: str
    powers: tuple
    log_scale_parts: tuple
    coefficients: IntervalTaylor

    def in_exact_units(self,powers):
        if tuple(map(Fraction,powers))!=self.powers:raise ValueError('Identical explicit pulse source units required')
        return self.coefficients

    @source_precision
    def report(self):
        c=self.coefficients.ctx
        widths=[c.mpf(radius.post.selected.inlet.endpoints(v)[1])-c.mpf(radius.post.selected.inlet.endpoints(v)[0])
            for v in self.coefficients.coefficients]
        return dict(name=self.name,scale_units=UNITS,
            exact_scale_powers=[dict(numerator=v.numerator,denominator=v.denominator) for v in self.powers],
            positive_scale_log_parts={key:value.node for key,value in self.log_scale_parts},
            signed_scaled_axial_Taylor_enclosures=self.coefficients,
            directed_scaled_coefficient_width_bounds=widths,retained_axial_order=self.coefficients.order,
            pulse_F_definition='exp(-(.5+mu)*log(R/Rp)); exact source, not an exponential cap',
            coefficient_convention='dZ^j/j!; exact radial scale is Z independent',
            source_enclosure_not_a_polynomial_field=True,numerical_absolute_point_value_installed=False)


class CurrentOriginalRpRawPulseTransport:
    @source_precision
    def __init__(self,heat=None,require_checked=True):
        self.closed=heat if heat is not None else closed.CurrentOriginalRpSameRepairHeatClosure()
        if type(self.closed) is not closed.CurrentOriginalRpSameRepairHeatClosure or not self.closed.acceptance_loaded:
            raise ValueError('Accepted current same-repair closed heat source required')
        self.closed.assert_graph()
        # Only the expression view is new. All 15 current providers and the
        # accepted selection/repair remain the exact same live objects.
        self.raw=raw.CurrentOriginalRpRawHistoryTransport(self.closed.before.before)
        self.radius=self.raw.radius;self.graph=self.raw.graph;self.post=self.raw.post
        self.selected=self.post.before;self.ctx=self.post.ctx;self.family_record=self.post.family_record
        self.hashes=dict(self.closed.hashes)
        radius.post.selected.inlet.add_hashes(self.hashes,self.raw.hashes)
        for name in (closed.NAME,closed.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.call_trace=[];self.acceptance_loaded=False
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current pulse transport receipt/scope differs')
            if receipt['source_family']!=self.family_record:raise ValueError('Current pulse family differs')
            radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        b=self.selected;c=self.closed
        graph=dict(accepted_current_selected_source=b.acceptance_loaded and all(b.assert_graph().values()),
            same_postpulse_and_closed_heat_provider=self.post is c.post,
            same_selected_pulse_on_all_six_routes=all(self.post.provider(chart) is b.pulse for chart in PULSE),
            same_unique_repair_and_future=b.seed.exact.repair is b.future.repair is c.exact.repair
                and b.future is c.exact.future,
            same_independent_analytic_P0=b.datum is c.exact.repair.angular.initial.datum,
            independent_expression_graph_same_absolute_origin=self.graph is not c.before.graph
                and self.radius.logRp.node==b.current_Rp_physical_log_radius.node,
            unchanged_current_beta_and_C5_context=b.pulse.flat.ctx is b.fifth.ctx is self.ctx,
            closed_heat_radius_and_amplitude_source_binding=c.absolute_binding['passed'],
            complete_fifteen_chart_registry=set(CHARTS)==set(self.post.registry) and len(CHARTS)==15)
        if not all(graph.values()):raise ValueError('Pulse/closed exterior provider graph differs: '+str(graph))
        return graph

    def scale(self,chart,coordinate,powers):
        """Combine R and pulse decay before astronomical source arithmetic."""
        r,e,p=map(Fraction,powers);g=self.graph;mu=self.radius.functions['mu']
        value=radius.exact_coordinate(coordinate);v=g.constant(value);weight=r-e/2
        inv=g.zero;finite=g.zero;local=g.zero
        if chart=='pulse_entrance':
            local=g.mul(g.sub(g.constant(weight),g.mul(g.constant(e),mu)),v)
        elif chart in ('pulse_main','pulse_exit','pulse_gap'):
            inv=g.quotient(g.mul(g.constant(weight),v),mu,'accepted current positive mu')
            finite=g.mul(g.constant(-e),v)
        elif chart in ('pulse_gap_end','pulse_end'):
            inv=g.quotient(g.constant(13*weight),mu,'accepted current positive mu')
            finite=g.constant(-13*e)
            local=g.mul(g.sub(g.constant(weight),g.mul(g.constant(e),mu)),v)
        else:raise ValueError('Actual current pulse source chart required')
        return dict(absolute_Rp_origin=g.mul(g.constant(r),self.radius.logRp),
            logPstar=g.mul(g.constant(p),self.radius.functions['logP']),
            combined_inverse_mu_pulse=inv,finite_pulse_decay=finite,local_native_offset=local)

    def factor(self,name,chart,coordinate,powers,coefficients):
        if type(coefficients) is not IntervalTaylor or coefficients.ctx is not self.ctx or coefficients.order not in (4,5):
            raise ValueError('Actual current same-context axial C4/C5 source rows required')
        return FactorizedPulseSourceTaylor(name,tuple(map(Fraction,powers)),
            tuple(self.scale(chart,coordinate,powers).items()),coefficients)

    def expression_owner(self,chart):
        if chart not in CHARTS:raise ValueError('Unknown current source chart')
        return self.raw

    @source_precision
    def pulse_view(self,chart,Z,coordinate):
        self.assert_graph();geometry=self.radius.geometry(chart,coordinate)
        if chart not in PULSE:raise ValueError('Actual pulse chart required')
        value=radius.exact_coordinate(coordinate);c=self.ctx;v=c.mpf(value.numerator)/value.denominator
        packet=self.selected.evaluate(chart,Z,v)['source_packet']
        u=packet['inlet_Utheta_over_Pstar_Taylor'];B=packet['Uz_over_Utheta']
        if u.order!=5 or B.order!=5:raise ValueError('True current inlet and pulse axial C5 rows required')
        jets=dict(Mz=u*packet['Mz_over_R_Utheta'],
            Mtheta=u*packet['Mtheta_over_sqrt2_R_3half_Utheta']*c.sqrt(2),
            Mtheta_z=u*u*packet['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']*c.sqrt(2),
            Mztheta=u*u*packet['Mztheta_over_R_Utheta_squared'],
            Mp=packet['pressure']['Mp_over_Pstar_squared'],P0=packet['pressure']['P0_over_Pstar_squared'],
            pressure=packet['pressure']['P_over_Pstar_squared'],Utheta=u,Uz=u*B,
            Ur=u.truncate(4)*packet['Ur_over_sqrt_R_over_2_Utheta']/c.sqrt(2))
        histories={name:self.factor(name,chart,coordinate,POWERS[name],value) for name,value in jets.items()}
        # These are the general cumulative equations, including nonzero Uz.
        jac=self.radius.local_step(chart,coordinate,coordinate)['native_to_log_radius_jacobian_bound']
        zero=IntervalTaylor.constant(c,0,5)
        densities=dict(Mz=u*B,Mtheta=u*c.sqrt(2),Mtheta_z=u*u*B*c.sqrt(2),
            Mztheta=u*u*(B*B-c.mpf('1/2')),Mp=u*u/2,P0=zero,pressure=u*u/2)
        derivative_powers={**POWERS,'Mp':(0,2,2),'pressure':(0,2,2)}
        derivatives={name:self.factor('d_native_'+name,chart,coordinate,derivative_powers[name],value*jac)
            for name,value in densities.items()}
        return dict(chart=chart,Z=packet['Z'],geometry=geometry,raw_histories=histories,
            similarity_velocity={name:histories[name] for name in ('Ur','Utheta','Uz')},
            first_native_radial_derivatives=derivatives,source_packet=packet,
            directed_native_Jacobian_bound=jac,exact_native_Jacobian_node=geometry['native_to_log_radius_jacobian'],
            original_P0_and_incoming_history_memory_retained=True,
            post_Rv_terminal_zeros_not_applied_in_active_pulse=True,
            no_radial_decay_cap_used_as_defining_amplitude=True,
            physical_time_Cartesian_velocity_installed=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        """One current caller: six pulse charts plus nine postpulse charts."""
        if chart not in CHARTS:raise ValueError('Unknown current pulse/postpulse chart')
        if chart in PULSE:
            view=self.pulse_view(chart,Z,coordinate);kind='active_current_selected_pulse'
        elif chart in ('heat_collar','heat_exterior'):
            view=dict(self.closed.evaluate(chart,Z,coordinate));kind='same_repair_closed_heat'
            # Rebase only exact scale expressions into this caller's graph.
            # Keep the proved closed coefficients, not the forward diagnostics.
            geometry=self.radius.geometry(chart,coordinate)
            for key in ('raw_histories','first_native_radial_derivatives'):
                view[key]={name:self.raw.factor(value.name,chart,geometry,value.powers,value.coefficients)
                    for name,value in view[key].items()}
            view['geometry']=view['current_absolute_source_geometry']=geometry
        else:
            view=dict(self.raw.evaluate(chart,Z,coordinate));kind='current_postpulse_raw'
        if chart not in PULSE:
            zero=IntervalTaylor.constant(self.ctx,0,5);geometry=view['geometry']
            factor_owner=self.expression_owner(chart)
            view['similarity_velocity']=dict(Utheta=view['raw_histories']['Utheta'],
                Uz=factor_owner.factor('Uz',chart,geometry,(0,1,0),zero),
                Ur=factor_owner.factor('Ur',chart,geometry,(Fraction(1,2),1,0),zero))
            view['chart']=chart
        self.call_trace.append(dict(chart=chart,actual_provider_kind=kind,
            same_selected_repair=True,closed_heat_provider_used=chart in ('heat_collar','heat_exterior')))
        view.update(source_family=self.family_record,actual_provider_kind=kind,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        return view


def view_report(view):
    omitted={'raw_histories','first_native_radial_derivatives','similarity_velocity','source_packet',
        'source_normalized_histories','current_source_packet','closed_normalized_histories'}
    return {**{key:value for key,value in view.items() if key not in omitted},
        'raw_histories':{key:value.report() for key,value in view['raw_histories'].items()},
        'similarity_velocity':{key:value.report() for key,value in view['similarity_velocity'].items()},
        'first_native_radial_derivatives':{key:value.report() for key,value in view['first_native_radial_derivatives'].items()}}


@source_precision
def run(heat=None):
    began=time.monotonic();owner=CurrentOriginalRpRawPulseTransport(heat,require_checked=False)
    Z='.479';views={}
    for name,chart,coordinate in VIEWS:
        views[name]=view_report(owner.evaluate(chart,Z,coordinate))
        print('Actual current pulse/closed exterior:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_pulse_transport_constructed=True,
        actual_same_current_provider_graph=owner.assert_graph(),fresh_Z=Z,
        actual_fifteen_current_views=views,actual_source_call_trace=owner.call_trace,
        original_selected_repair_and_analytic_P0_unchanged=True,
        pulse_scale_definition='R^r*Fpulse^e*Pstar^p; Fpulse=exp(-(.5+mu)*log(R/Rp))',
        source_scope='Six pulse raw histories/velocity components and first native transport; common 15-chart source caller with closed heat. Unrestricted numerical points and global/time Cartesian admission remain open',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner


if __name__=='__main__':run()
