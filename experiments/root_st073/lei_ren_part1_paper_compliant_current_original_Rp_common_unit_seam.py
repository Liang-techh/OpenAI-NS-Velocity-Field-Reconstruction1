"""Actual current Rv endpoint in common source units and mixed logR/Z rows.

Only pulse_end(s=0) and flatten(t=0) are joined here. Enclosure overlap is
diagnostic; the exact function join consumes the same-source endpoint proof.
No absolute amplitude, enormous radius or unrestricted physical point is
materialized. Every derivative retains its explicit positive source scale.
"""
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_raw_pulse_transport as pulse
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=pulse.HERE,pulse.PREFIX,pulse.sha
NAME=PREFIX+'current_original_Rp_common_unit_seam.json.gz'
RECEIPT=PREFIX+'current_original_Rp_common_unit_seam_check.json'
GATES=('current_original_Rp_Rv_common_source_unit_adapter_installed',
    'current_original_Rp_Rv_same_source_five_history_function_join_installed',
    'current_original_Rp_Rv_log_radius_mixed_rows_total_order4_installed')
OPEN=pulse.OPEN
POWERS={key:(Fraction(r),Fraction(e),Fraction(p-e)) for key,(r,e,p) in pulse.POWERS.items()}


class CurrentOriginalRpCommonUnitSeam:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else pulse.CurrentOriginalRpRawPulseTransport()
        if type(self.before) is not pulse.CurrentOriginalRpRawPulseTransport or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current pulse/closed heat caller required')
        self.before.assert_graph()
        self.raw=self.before.raw;self.graph=self.before.graph;self.ctx=self.before.ctx
        self.selected=self.before.selected;self.family_record=self.before.family_record
        self.binding=self.selected.flatten_binding
        if not self.binding['passed'] or self.binding['interval_overlap_used_as_join_proof']:
            raise ValueError('Current same-source defining endpoint theorem required')
        self.U0=self.raw.U0;self.logU0=self.raw.logU0
        # This scalar box encloses the admitted exact U0 function. It is used
        # only for diagnostic unit conversion, never as a selected amplitude.
        self.U0_bound=self.selected.constants['U']
        if pulse.radius.post.selected.inlet.endpoints(self.U0_bound)[0]<=0:
            raise ArithmeticError('Positive current U0 function bound required')
        self.hashes=dict(self.before.hashes)
        for name in (pulse.NAME,pulse.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        theorem=self.binding['canonical_endpoint_theorem']
        self.endpoint_receipt=json.loads((HERE/theorem).read_bytes())
        proof=self.endpoint_receipt
        if not (proof['all_passed'] and proof['pulse_end_flatten_full_moment_stress_pressure_functional_join_verified']
                and proof['exact_functional_join_identities']==155
                and proof['actual_local_inlet_and_datum_source_identities']==34
                and proof['actual_endpoint_source_identities']==24
                and not proof['source_caps_used_as_defining_field_values']
                and not proof['interval_overlap_used_as_join_proof']):
            raise ValueError('Accepted canonical endpoint function theorem required')
        if (proof['actual_five_defect_family_sha256'],proof['implicit_source_sha256'])!=(
                self.selected.family,self.selected.source):raise ValueError('Endpoint theorem family differs')
        pulse.radius.post.selected.inlet.add_hashes(self.hashes,proof['input_hashes'])
        self.hashes[theorem]=sha(theorem)
        self.acceptance_loaded=False;self.call_trace=[]
        self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[k] for k in GATES) or any(record[k] for k in OPEN):
                raise ValueError('Current common-unit endpoint receipt/scope differs')
            if record['source_family']!=self.family_record:raise ValueError('Current endpoint family differs')
            pulse.radius.post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        b=self.before;s=self.selected;f=s.flatten
        result=dict(accepted_fifteen_chart_closed_caller=b.acceptance_loaded and all(b.assert_graph().values()),
            actual_same_pulse_and_flatten=f.pulse is s.pulse is b.post.outer.pulse,
            same_current_future_energy=f.pulse.selection.future is s.pulse.selection.future
                and s.pulse.selection.future.high.energy.future.__self__ is s.fifth
                and s.fifth.fourth.energy.base is s.future,
            same_current_independent_P0=f.inlet.datum is s.datum,
            admitted_current_endpoint_function_binding=self.binding is s.flatten_binding and self.binding['passed'],
            exact_current_U0_and_common_graph=self.U0 is self.raw.U0 and self.graph is self.raw.graph,
            same_positive_U0_bound=self.U0_bound is s.constants['U'] is f.inlet.constants['U'],
            same_source_context=self.ctx is s.ctx is f.ctx)
        if not all(result.values()):raise ValueError('Current Rv endpoint source differs: '+str(result))
        return result

    def factor(self,name,powers,coefficients,geometry):
        if type(coefficients) is not IntervalTaylor or coefficients.ctx is not self.ctx or not 0<=coefficients.order<=5:
            raise ValueError('Actual current directed Taylor source rows required')
        powers=tuple(map(Fraction,powers))
        return pulse.raw.FactorizedSourceTaylor(name,powers,
            tuple(self.raw.scale('flatten',geometry,powers).items()),coefficients)

    def converted_terminal(self,view,geometry):
        if view['chart']!='pulse_end' or view['geometry']['exact_native_coordinate']!={'numerator':0,'denominator':1}:
            raise ValueError('Common pulse-to-post conversion is defined only at the exact Rv endpoint')
        result={}
        for key,value in view['raw_histories'].items():
            r,e,p=value.powers
            if e.denominator!=1:raise ValueError('Integer pulse amplitude degree required')
            result[key]=self.factor(key,(r,e,p-e),value.coefficients/self.U0_bound**int(e),geometry)
        return result

    def mixed(self,base,theta,geometry):
        """True d_logR^k dZ^j, k+j<=4, in cancellation-safe source units.

        At Rv both flat schedules have vanishing derivatives. The general
        cumulative equations and theta_y=-bp*theta give the rows below.
        Pressure radial derivatives use Ev0^2, not a tiny decay cap/Pstar^2.
        """
        c=self.ctx;mu=self.selected.pulse.mu;bp=c.mpf('1/2')+mu
        a=1-mu;prate=2*bp;zero=theta*0;rows={}
        for key,value in base.items():
            derivatives=[(value.powers,value.coefficients)]
            for k in range(1,5):
                powers=value.powers
                if key=='Utheta':jet=theta*(-bp)**k
                elif key=='Mtheta':jet=theta*c.sqrt(2)*a**(k-1)
                elif key=='Mztheta':jet=-theta*theta/2*(-2*mu)**(k-1)
                elif key in ('Mp','pressure'):
                    powers=(0,2,0);jet=theta*theta/2*(-prate)**(k-1)
                else:jet=zero
                derivatives.append((powers,jet))
            rows[key]={}
            for k,(powers,jet) in enumerate(derivatives):
                for j in range(5-k):
                    coefficient=IntervalTaylor.constant(c,jet[j]*math.factorial(j),0)
                    rows[key]['y'+str(k)+'_Z'+str(j)]=self.factor(
                        'd_logR_'+str(k)+'_dZ_'+str(j)+'_'+key,powers,coefficient,geometry)
        return rows

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();c=self.ctx
        left=self.before.evaluate('pulse_end',Z,0);right=self.before.evaluate('flatten',Z,0)
        geometry=right['geometry'];packet=left['source_packet']
        z=packet['Z'];q=IntervalTaylor(c,[1+z*z,2*z,1,0,0,0]);theta=q.reciprocal()
        # The accepted native inlet is u=U0/q. Reduce the exact U0 factor
        # before interval arithmetic; the divided original rows remain below.
        jets=dict(Mz=theta*packet['Mz_over_R_Utheta'],
            Mtheta=theta*packet['Mtheta_over_sqrt2_R_3half_Utheta']*c.sqrt(2),
            Mtheta_z=theta*theta*packet['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']*c.sqrt(2),
            Mztheta=theta*theta*packet['Mztheta_over_R_Utheta_squared'],
            Mp=packet['pressure']['Mp_over_Pstar_squared'],P0=packet['pressure']['P0_over_Pstar_squared'],
            pressure=packet['pressure']['P_over_Pstar_squared'],Utheta=theta,Uz=theta*packet['Uz_over_Utheta'],
            Ur=theta.truncate(4)*packet['Ur_over_sqrt_R_over_2_Utheta']/c.sqrt(2))
        reduced={key:self.factor(key,POWERS[key],value,geometry) for key,value in jets.items()}
        post=dict(right['raw_histories']);post.update(right['similarity_velocity'])
        converted=self.converted_terminal(left,geometry)
        mixed_left=self.mixed(reduced,theta,geometry)
        mixed_right=self.mixed(post,right['similarity_velocity']['Utheta'].coefficients,geometry)
        self.call_trace.append(dict(Z=z,actual_pulse_end_zero_and_flatten_zero_called=True,
            same_repair_and_P0=True,exact_amplitude_factor_cancelled_before_arithmetic=True))
        return dict(Z=z,geometry=geometry,common_scale_units=pulse.raw.UNITS,
            converted_pulse_terminal_enclosures=converted,
            reduced_pulse_terminal=reduced,actual_flatten_terminal=post,
            pulse_log_radius_mixed_rows=mixed_left,flatten_log_radius_mixed_rows=mixed_right,
            current_source_function_binding=self.binding,
            retained_canonical_function_join_identity_counts=dict(functional=155,inlet_datum=34,endpoint=24),
            current_exact_U0_node=self.U0.node,current_exact_logU0_node=self.logU0.node,
            independent_P0_and_nonzero_quadratic_pressure_memory_retained=True,
            ordinary_derivative_convention='d_logR^k dZ^j; k+j<=4; factorial applied only to axial Taylor coefficients',
            source_function_join_is_not_enclosure_overlap=True,
            pressure_radial_rows_use_exact_Ev0_squared=True,
            no_cap_absolute_amplitude_or_radius_used_as_defining_value=True,
            uniform_or_global_mixed_admission=False,physical_time_Cartesian_velocity_installed=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


def view_report(view):
    flat=('converted_pulse_terminal_enclosures','reduced_pulse_terminal','actual_flatten_terminal')
    mixed=('pulse_log_radius_mixed_rows','flatten_log_radius_mixed_rows')
    return {**{k:v for k,v in view.items() if k not in flat+mixed},
        **{k:{name:value.report() for name,value in view[k].items()} for k in flat},
        **{k:{name:{label:value.report() for label,value in rows.items()}
            for name,rows in view[k].items()} for k in mixed}}


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpCommonUnitSeam(before,require_checked=False)
    views={}
    for name,Z in (('fresh','.521'),('cell',['.520','.522'])):
        views[name]=view_report(owner.evaluate(Z));print('Actual Rv common-unit seam:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_common_unit_Rv_seam_constructed=True,
        actual_source_graph=owner.assert_graph(),actual_Rv_seam_views=views,actual_source_call_trace=owner.call_trace,
        exact_unit_conversion='(R,Fpulse,Pstar):(r,e,p) -> (R,Ev0,Pstar):(r,e,p-e), coefficient/U0^e at Rv',
        source_scope='Same-current Rv function join, directed common-unit histories and mixed logR/Z rows through total order4; no uniform/global or physical point admission',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner


if __name__=='__main__':run()
