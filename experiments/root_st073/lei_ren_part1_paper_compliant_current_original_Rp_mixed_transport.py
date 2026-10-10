"""Current fifteen-chart mixed velocity, pressure and raw history transport.

All ordinary logR/Z derivatives through total order four retain explicit
source units. Native derivatives carry the exact positive Jacobian power
as a separate scale function. Amplitude/pressure caps are never derivative
definitions. This does not certify uniform seams or physical time points.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_common_unit_seam as seam
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

pulse=seam.pulse
HERE,PREFIX,sha=seam.HERE,seam.PREFIX,seam.sha
NAME=PREFIX+'current_original_Rp_mixed_transport.json.gz'
RECEIPT=PREFIX+'current_original_Rp_mixed_transport_check.json'
GATES=('current_original_Rp_fifteen_chart_raw_history_mixed_logR_Z4_installed',
    'current_original_Rp_fifteen_chart_velocity_pressure_mixed_logR_Z4_installed',
    'current_original_Rp_native_mixed_derivatives_exact_Jacobian_scale_installed')
OPEN=seam.OPEN
CHARTS=pulse.CHARTS
UT='Utheta_over_Pstar_without_common_theta_radial_factor'
UZ='Uz_over_Pstar_without_common_theta_radial_factor'
UR='Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor'


def product_row(left,right,k,order):
    result=left[0].truncate(order)*0
    for j in range(k+1):
        result+=left[j].truncate(order)*right[k-j].truncate(order)*math.comb(k,j)
    return result


def shifted_row(rows,rate,k,order):
    result=rows[0].truncate(order)*0
    for j in range(k+1):result+=rows[j].truncate(order)*(math.comb(k,j)*rate**(k-j))
    return result


@dataclass(frozen=True)
class FactorizedMixedSourceRow:
    name: str
    derivative: tuple
    source_units: tuple
    powers: tuple
    log_scale_parts: tuple
    coefficients: IntervalTaylor

    @source_precision
    def report(self):
        lo,hi=pulse.radius.post.selected.inlet.endpoints(self.coefficients[0])
        c=self.coefficients.ctx
        return dict(name=self.name,ordinary_derivative_index=self.derivative,
            scale_units=self.source_units,
            exact_scale_powers=[dict(numerator=v.numerator,denominator=v.denominator) for v in self.powers],
            positive_scale_log_parts={k:v.node for k,v in self.log_scale_parts},
            signed_scaled_derivative_enclosure=self.coefficients,
            directed_scaled_width_bound=c.mpf(hi)-c.mpf(lo),
            scalar_is_an_ordinary_derivative_not_a_Taylor_polynomial=True,
            numerical_absolute_point_value_installed=False)


class CurrentOriginalRpMixedTransport:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else seam.CurrentOriginalRpCommonUnitSeam()
        if type(self.before) is not seam.CurrentOriginalRpCommonUnitSeam or not self.before.acceptance_loaded:
            raise ValueError('Accepted current common-unit Rv and closed fifteen-chart caller required')
        self.before.assert_graph();self.owner=self.before.before
        self.graph=self.owner.graph;self.ctx=self.owner.ctx;self.family_record=self.owner.family_record
        self.hashes=dict(self.before.hashes)
        for name in (seam.NAME,seam.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self.call_trace=[]
        self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[k] for k in GATES) or any(record[k] for k in OPEN):
                raise ValueError('Current all-chart mixed receipt/scope differs')
            if record['source_family']!=self.family_record:raise ValueError('Current mixed family differs')
            pulse.radius.post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        o=self.owner;c=o.closed
        result=dict(accepted_Rv_same_source_common_unit_join=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            one_current_fifteen_chart_provider=self.owner is self.before.before and all(o.assert_graph().values()),
            one_current_output_expression_graph=self.graph is o.graph is o.raw.graph,
            closed_heat_shape_and_forward_source=c.angular.heat.ctx is self.ctx
                and c.angular.exact is c.exact and c.absolute_binding['passed'],
            same_current_repair_and_P0=c.exact.repair is o.selected.seed.exact.repair
                and o.selected.datum is c.exact.flatten.inlet.datum,
            original_selected_N_unchanged=True)
        if not all(result.values()):raise ValueError('Current mixed source graph differs: '+str(result))
        return result

    @source_precision
    def velocity_rows(self,chart,Z,coordinate,view):
        c=self.ctx
        if chart in pulse.PULSE:
            packet=view['source_packet']
            rows=packet['physical_velocity_and_pressure_y_derivative_Taylor']
            return {key:list(rows[label]) if key!='Ur' else [v/c.sqrt(2) for v in rows[label]]
                for key,label in (('Utheta',UT),('Uz',UZ),('Ur',UR))},dict(
                    source='actual_current_pulse_physical_radial_product_rows',
                    physical_radial_prefactors_already_differentiated=True)
        if chart not in ('heat_collar','heat_exterior'):
            packet=view['current_source_packet']
            if chart=='flatten':
                valid_units=packet['component_units']['velocity']=='U/Ev0; exact formal Ev0/Pstar kept separately'
            else:
                valid_units=packet['velocity_units']=='U/Ev0; exact formal Ev0/Pstar retained separately'
            if not valid_units:raise ValueError('Actual postpulse velocity rows must use admitted Ev0 units')
            rows=packet['physical_velocity_and_pressure_y_derivative_Taylor']
            return {key:list(rows[label]) if key!='Ur' else [v/c.sqrt(2) for v in rows[label]]
                for key,label in (('Utheta',UT),('Uz',UZ),('Ur',UR))},dict(
                    source='actual_current_postpulse_velocity_ODE_rows',
                    old_cap_scaled_pressure_derivatives_not_consumed=True)
        # Use precisely the closed owner used by the current field, including
        # its source-bound heat view. Do not mix in a second forward branch.
        close=self.owner.closed;heat=close.angular.heat;z=c.mpf(Z)
        value=pulse.radius.exact_coordinate(coordinate);t=c.mpf(value.numerator)/value.denominator
        shape=(pulse.closed.closure.angular.shape_radial5(heat,z,t) if chart=='heat_collar'
            else heat.local_Gamma(z,t))
        K=shape['K_rows'][:5]
        forward=close.angular.terminal_constants(z)['current_repaired_forward_terminal']
        if view['source_packet']['exact_pressure_scale'] is not forward['pressure_scale']:
            raise ValueError('Closed heat rows must retain the same source-bound forward amplitude datum')
        th=[(shifted_row(K,-heat.bh,k,5)*forward['theta_base']*c.exp(-heat.bh*t)).truncate(4-k)
            for k in range(5)]
        zero=IntervalTaylor.constant(c,0,5)
        return dict(Utheta=th,Uz=[zero.truncate(4-k) for k in range(5)],
            Ur=[zero.truncate(4-k) for k in range(5)]),dict(
                source='same_closed_current_full_Gamma_K_product_rows',
                full_infinite_Gamma_not_a_finite_series=True,
                closed_coefficients_and_forward_diagnostics_kept_separate=True)

    def histories(self,velocities):
        """General cumulative densities, including active nonzero Uz.

        Velocity rows already include derivatives of their radial scales.
        The shifts below differentiate only the remaining R or R^(3/2).
        P0 is independent; Mp/pressure radial rows use true theta^2 units.
        """
        c=self.ctx;E=velocities['Utheta'];V=velocities['Uz'];result={}
        for k in range(1,5):
            n=4-k
            uv=[product_row(E,V,j,n) for j in range(k)]
            uu=[product_row(E,E,j,n) for j in range(k)]
            vv=[product_row(V,V,j,n) for j in range(k)]
            signed=[v-u/2 for v,u in zip(vv,uu)]
            zero=E[0].truncate(n)*0
            result[k]=dict(Mz=shifted_row(V,1,k-1,n),
                Mtheta=shifted_row(E,c.mpf('3/2'),k-1,n)*c.sqrt(2),
                Mtheta_z=shifted_row(uv,c.mpf('3/2'),k-1,n)*c.sqrt(2),
                Mztheta=shifted_row(signed,1,k-1,n),
                Mp=uu[k-1]/2,P0=zero,pressure=uu[k-1]/2)
        return result

    def factor(self,name,chart,coordinate,geometry,powers,value,k,j,native):
        if type(value) is not IntervalTaylor or value.order!=0 or value.ctx is not self.ctx:
            raise ValueError('Same-context directed ordinary derivative enclosure required')
        powers=tuple(map(Fraction,powers));g=self.graph
        if chart in pulse.PULSE:
            units=pulse.UNITS;scale=self.owner.scale(chart,coordinate,powers)
        else:
            units=pulse.raw.UNITS;scale=self.owner.raw.scale(chart,geometry,powers)
        jac=Fraction(k if native else 0)
        J=pulse.radius.FunctionRef(g,geometry['native_to_log_radius_jacobian'])
        scale['native_coordinate_Jacobian']=g.mul(g.constant(jac),g.unary('log',J)) if jac else g.zero
        return FactorizedMixedSourceRow(name,(k,j),units+('native_Jacobian',),powers+(jac,),
            tuple(scale.items()),value)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.assert_graph()
        if chart not in CHARTS:raise ValueError('Actual current source chart required')
        view=self.owner.evaluate(chart,Z,coordinate)
        geometry=view['geometry'];base=dict(view['raw_histories']);base.update(view['similarity_velocity'])
        powers={key:value.powers for key,value in base.items()}
        dp={**powers,'Mp':(0,2,2 if chart in pulse.PULSE else 0),
            'pressure':(0,2,2 if chart in pulse.PULSE else 0)}
        velocities,evidence=self.velocity_rows(chart,Z,coordinate,view)
        for key,rows in velocities.items():
            if len(rows)!=5 or any(row.order<4-k or row.ctx is not self.ctx for k,row in enumerate(rows)):
                raise ValueError('Actual velocity radial/axial rows through total order four required: '+key)
        histories=self.histories(velocities);logrows={};nativerows={}
        for key,value in base.items():
            logrows[key]={};nativerows[key]={}
            for k in range(5):
                row=(value.coefficients if k==0 else velocities[key][k] if key in velocities else histories[k][key])
                pp=dp[key] if k else powers[key]
                for j in range(5-k):
                    coefficient=IntervalTaylor.constant(self.ctx,row[j]*math.factorial(j),0)
                    label='y'+str(k)+'_Z'+str(j)
                    logrows[key][label]=self.factor(key,chart,coordinate,geometry,pp,coefficient,k,j,False)
                    nativerows[key]['n'+str(k)+'_Z'+str(j)]=self.factor(key,chart,coordinate,geometry,pp,coefficient,k,j,True)
        self.call_trace.append(dict(chart=chart,coordinate=coordinate,actual_current_source_called=True,
            exact_native_Jacobian_is_scale_not_bound_multiplier=True,closed_heat=chart.startswith('heat_')))
        return dict(chart=chart,Z=view.get('Z',self.ctx.mpf(Z)),geometry=geometry,
            original_factorized_values=base,log_radius_mixed_rows=logrows,native_coordinate_mixed_rows=nativerows,
            velocity_radial_source_rows=velocities,velocity_source_evidence=evidence,
            five_general_cumulative_density_equations_used=True,
            active_axial_quadratic_terms_and_independent_P0_preserved=True,
            radial_velocity_source_order4_retained=True,
            derivative_scope='ordinary log similarity R and axial Z total order<=4; native radial rows multiply exact J^k',
            original_forward_source_packet=view.get('source_packet',view.get('current_source_packet')),
            absolute_point_values_and_uniform_global_seams_installed=False,
            physical_time_Cartesian_derivatives_installed=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


def view_report(view):
    skip=('original_factorized_values','log_radius_mixed_rows','native_coordinate_mixed_rows',
        'velocity_radial_source_rows','original_forward_source_packet')
    return {**{k:v for k,v in view.items() if k not in skip},
        'original_factorized_values':{k:v.report() for k,v in view['original_factorized_values'].items()},
        **{group:{name:{k:v.report() for k,v in rows.items()} for name,rows in view[group].items()}
            for group in ('log_radius_mixed_rows','native_coordinate_mixed_rows')}}


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpMixedTransport(before,require_checked=False);views={}
    for name,chart,coordinate in pulse.VIEWS:
        views[name]=view_report(owner.evaluate(chart,'.521',coordinate))
        print('Actual current raw mixed source:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_fifteen_chart_mixed_transport_constructed=True,
        actual_source_graph=owner.assert_graph(),fresh_Z='.521',actual_fifteen_chart_mixed_views=views,
        actual_source_call_trace=owner.call_trace,
        source_scope='Raw histories, velocities and pressure mixed logR/Z4 source rows on fifteen actual routes; exact native Jacobian powers. Uniform/global seams and physical-time point evaluation remain open',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner


if __name__=='__main__':run()
