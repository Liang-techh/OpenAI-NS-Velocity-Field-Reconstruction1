"""Current postpulse raw histories and their native radial transport.

Raw values are executable factorized source enclosures: a lazy positive
absolute scale times directed axial Taylor coefficients. No enormous R,
vanishing Ev0, interval endpoint or amplitude cap becomes a point value.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius as radius
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=radius.HERE,radius.PREFIX,radius.sha
NAME=PREFIX+'current_original_Rp_raw_history_transport.json.gz'
RECEIPT=PREFIX+'current_original_Rp_raw_history_transport_check.json'
GATES=('current_original_Rp_factorized_raw_postpulse_histories_callable',
       'current_original_Rp_raw_postpulse_first_native_radial_transport_installed')
OPEN=radius.OPEN
CHARTS=('flatten',)+tuple(radius.post.METHODS)
UNITS=('R','Ev0','Pstar')
POWERS=dict(Mz=(1,1,0),Mtheta=(Fraction(3,2),1,0),
    Mtheta_z=(Fraction(3,2),2,0),Mztheta=(1,2,0),Mp=(0,0,2),
    P0=(0,0,2),pressure=(0,0,2),Utheta=(0,1,0))
VIEWS=(('flatten_inlet','flatten',0),('flatten_active','flatten',37),
    ('flatten_exit','flatten',100),('power_inlet','outer_power',0),
    ('power_exit','outer_power',1),('angular_inlet','outer_angular',-4),
    ('angular_active','outer_angular',-3),('angular_exit','outer_angular',0),
    ('entry_inlet','steep_entry',0),('entry_exit','steep_entry',1),
    ('steep_inlet','steep_power',0),('steep_exit','steep_power',1),
    ('exit_inlet','steep_exit',0),('exit_exit','steep_exit',1),
    ('waiting_inlet','waiting',0),('waiting_exit','waiting',1),
    ('collar_inlet','heat_collar',0),('collar_exit','heat_collar',3),
    ('exterior_inlet','heat_exterior',3),('exterior_fresh','heat_exterior',4))


def packed(value):return radius.post.selected.inlet.encode(radius.post.selected.inlet.pack(value))


@dataclass(frozen=True)
class FactorizedSourceTaylor:
    """Finite axial derivative enclosures, multiplied by a true scale function.

    The scale is Z independent, so every coefficient has the same factor.
    Directed widths are error bounds in these units; no midpoint is chosen.
    """
    name: str
    powers: tuple
    log_scale_parts: tuple
    coefficients: IntervalTaylor

    def in_exact_units(self,powers):
        if tuple(Fraction(v) for v in powers)!=self.powers:
            raise ValueError('Explicit identical absolute scale units required')
        return self.coefficients

    @source_precision
    def report(self):
        c=self.coefficients.ctx
        widths=[]
        for value in self.coefficients.coefficients:
            lo,hi=radius.post.selected.inlet.endpoints(value)
            widths.append(c.mpf(hi)-c.mpf(lo))
        return dict(name=self.name,scale_units=UNITS,
            exact_scale_powers=[dict(numerator=v.numerator,denominator=v.denominator) for v in self.powers],
            positive_scale_log_parts={key:value.node for key,value in self.log_scale_parts},
            signed_scaled_axial_Taylor_enclosures=self.coefficients,
            directed_scaled_coefficient_width_bounds=widths,
            coefficient_convention='dZ^j/j!; exact scale is Z independent',
            source_enclosure_not_a_polynomial_field=True,
            numerical_absolute_point_value_installed=False)


class CurrentOriginalRpRawHistoryTransport:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else radius.CurrentOriginalRpSegmentedRadius()
        if type(self.before) is not radius.CurrentOriginalRpSegmentedRadius or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current absolute-radius owner required')
        self.before.assert_graph()
        # Build only an independent, cheap radius expression view. Existing
        # pulse/repair/future providers remain the same live objects.
        self.radius=radius.CurrentOriginalRpSegmentedRadius(self.before.before)
        self.post=self.before.before;self.ctx=self.post.ctx
        self.family_record=self.post.family_record;self.hashes=dict(self.before.hashes)
        self.graph=self.radius.graph;self.prefix=copy.deepcopy(self.graph.nodes)
        g=self.graph;frame=self.radius.frame.functions['actual_identified_Rp_native_frame_C3']
        self.U0=g.node('function_substitution',expression=radius.frame_source.target.rows(frame['u'])[0].node,
            variable=g.symbol('Z').node,value=g.zero.node,Z_independent_substitution=True)
        self.logU0=g.unary('log',self.U0)
        self.logP=self.radius.functions['logP'];self.mu=self.radius.functions['mu']
        self.logEv0_parts=dict(logPstar=self.logP,inlet_logU=self.logU0,
            pulse_decay=g.quotient(g.constant('-13/2'),self.mu,'accepted positive mu'),finite=g.constant(-13))
        self.call_trace=[];self.acceptance_loaded=False
        for name in (radius.NAME,radius.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[k] for k in GATES) or any(record[k] for k in OPEN):
                raise ValueError('Current raw-history receipt or scope differs')
            if record['source_family']!=self.family_record:raise ValueError('Raw history family differs')
            radius.post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        b=self.before;p=self.post;r=self.radius
        graph=dict(accepted_absolute_radius=b.acceptance_loaded and all(b.assert_graph().values()),
            same_actual_postpulse_owner=r.before is p is b.before,
            independent_radius_graph=r.graph is self.graph and self.graph is not b.graph,
            original_radius_prefix_unchanged=self.graph.nodes[:len(self.prefix)]==self.prefix,
            same_unique_selected_repair=p.heat.repair is p.before.future.repair is p.before.seed.exact.repair,
            exact_current_inlet_source=r.frame is b.frame,
            same_current_context=self.ctx is p.ctx is r.ctx,
            nine_postpulse_routes=set(CHARTS)=={'flatten'}|set(radius.post.METHODS))
        if not all(graph.values()):raise ValueError('Raw history source graph differs: '+str(graph))
        return graph

    def scale(self,chart,geometry,powers):
        """Cancel the common pulse terms algebraically, before any arithmetic."""
        r,e,p=map(Fraction,powers);g=self.graph
        pulse_weight=13*(r-e/2)
        return dict(absolute_Rp_origin=g.mul(g.constant(r),self.radius.logRp),
            local_radius_offset=g.mul(g.constant(r),radius.FunctionRef(g,geometry['offset'])),
            combined_inverse_mu_pulse=g.quotient(g.constant(pulse_weight),self.mu,'accepted positive mu'),
            logPstar=g.mul(g.constant(e+p),self.logP),
            inlet_logU=g.mul(g.constant(e),self.logU0),finite_amplitude=g.constant(-13*e))

    def factor(self,name,chart,geometry,powers,coefficients):
        if type(coefficients) is not IntervalTaylor or coefficients.ctx is not self.ctx or coefficients.order!=5:
            raise ValueError('Actual same-context axial C5 source enclosures required')
        powers=tuple(map(Fraction,powers))
        return FactorizedSourceTaylor(name,powers,tuple(self.scale(chart,geometry,powers).items()),coefficients)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.assert_graph()
        if chart not in CHARTS:raise ValueError('An actual postpulse history chart is required')
        # Validate exact source geometry before invoking expensive providers.
        geometry=self.radius.geometry(chart,coordinate)
        packet=self.post.evaluate(chart,Z,coordinate)['source_packet']
        hats=self.post.histories(chart,Z,coordinate,packet)
        theta=packet['theta_over_Ev0_Taylor'];theta2=theta*theta;c=self.ctx
        if radius.post.selected.inlet.endpoints(theta[0])[0]<=0:
            raise ArithmeticError('Positive current swirl source required')
        sqrt2=c.sqrt(c.mpf(2));zero=IntervalTaylor.constant(c,0,5)
        jets=dict(Mz=theta*hats['Mz_over_R_Utheta'],
            Mtheta=theta*hats['Mtheta_over_sqrt2_R_3half_Utheta']*sqrt2,
            Mtheta_z=theta2*hats['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']*sqrt2,
            Mztheta=theta2*hats['Mztheta_over_R_Utheta_squared'],
            Mp=hats['Mp_over_Pstar_squared'],P0=hats['P0_over_Pstar_squared'],
            pressure=hats['P_over_Pstar_squared'],Utheta=theta)
        raw={name:self.factor(name,chart,geometry,POWERS[name],value) for name,value in jets.items()}
        # These are the actual five cumulative-history equations in logR.
        # Uz=0 here follows the admitted current terminal theorem. Native
        # coordinate derivatives receive exactly one dlogR/dcoordinate.
        jac=self.radius.local_step(chart,coordinate,coordinate)['native_to_log_radius_jacobian_bound']
        native=dict(Mz=zero,Mtheta=theta*sqrt2*jac,Mtheta_z=zero,
            Mztheta=-theta2*jac/2,Mp=theta2*jac/2,P0=zero,pressure=theta2*jac/2)
        d_powers={**POWERS,'Mp':(0,2,0),'pressure':(0,2,0)}
        derivatives={name:self.factor('d_native_'+name,chart,geometry,d_powers[name],value)
            for name,value in native.items()}
        self.call_trace.append(dict(chart=chart,actual_current_source_called=True,
            exact_radius_and_Ev0_functions_used=True,all_five_histories_transported=True,
            one_native_coordinate_Jacobian_applied=True))
        return dict(chart=chart,Z=packet['Z'],geometry=geometry,
            raw_histories=raw,first_native_radial_derivatives=derivatives,
            source_normalized_histories=hats,current_source_packet=packet,
            exact_Ev0_log_parts={key:value.node for key,value in self.logEv0_parts.items()},
            exact_native_Jacobian_node=geometry['native_to_log_radius_jacobian'],
            directed_native_Jacobian_bound=jac,
            P0_is_independent_and_not_reset=True,heat_terminal_constants_not_eliminated=True,
            no_absolute_radius_amplitude_or_frequency_materialized=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


def view_report(view):
    omitted={'source_normalized_histories','current_source_packet','raw_histories','first_native_radial_derivatives'}
    return {**{key:value for key,value in view.items() if key not in omitted},
        'raw_histories':{key:value.report() for key,value in view['raw_histories'].items()},
        'first_native_radial_derivatives':{key:value.report() for key,value in view['first_native_radial_derivatives'].items()}}


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpRawHistoryTransport(before,require_checked=False)
    Z='.439';views={}
    for name,chart,coordinate in VIEWS:
        views[name]=view_report(owner.evaluate(chart,Z,coordinate))
        print('Actual current raw-history transport:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_raw_history_transport_constructed=True,
        actual_same_source_graph=owner.assert_graph(),fresh_Z=Z,actual_twenty_source_views=views,
        actual_source_call_trace=owner.call_trace,accepted_radius_graph_prefix_length=len(owner.prefix),
        appended_exact_function_nodes=owner.graph.nodes[len(owner.prefix):],
        exact_current_inlet_U0_node=owner.U0.node,exact_current_inlet_logU0_node=owner.logU0.node,
        retained_heat_terminal_constants=owner.post.heat_constants(Z),
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,
        scope='Factorized raw source histories and first native radial transport; no numeric absolute point field or terminal closure',
        execution_seconds=time.monotonic()-began)
    data=json.dumps(packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner


if __name__=='__main__':run()
