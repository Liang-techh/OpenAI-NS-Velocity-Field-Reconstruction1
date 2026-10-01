"""Callable coherent new-source leading layers through the O.2 axial buffer.

Equations and dimensionless bump bounds are reused; physical source moments,
pressure, defects and controls all come from the new Md1.1 core chain.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_Md11_comparison_jets import Md11ComparisonJets
from lei_ren_part1_paper_interval_long_reshape_field import IntervalLongReshapeField,T_EXACT,RM_LOG_EXACT
from lei_ren_part1_paper_interval_axial_restore_field import IntervalAxialRestoreField
from lei_ren_part1_paper_interval_functional_defects import logjet
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_jet,_pack
from lei_ren_part1_paper_interval_repaired_reference_field import evaluate_reference,validate_receipt
from lei_ren_part1_paper_interval_partial_five_bump_map import IntervalPartialFiveBumpMap
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition,reference_moments,fraction_box
from lei_ren_part1_paper_interval_outer_axial_turnoff import IntervalOuterAxialTurnoff
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


class Md11LongReshapeField(IntervalLongReshapeField):
    def __init__(self,calc):
        self.calc=calc;c=calc.ctx
        self.source_name='lei_ren_part1_paper_Md11_exit_switch.json'
        self.defects_name='lei_ren_part1_paper_Md11_functional_defects.json'
        self.source=json.loads((HERE/self.source_name).read_bytes())
        self.defects=json.loads((HERE/self.defects_name).read_bytes())
        for raw in (self.source,self.defects):
            validate_receipt(raw,'new-source connecting data')
            if raw['state_sha256']!=calc.state_hash or raw['new_schedule_sha256']!=calc.new_schedule_sha256:
                raise ValueError('new-source connecting identity mismatch')
        with mp.workdps(calc.precision+60):
            self.z=calc.z.truncate(1);self.P0=calc.p0.truncate(1)
            self.B=restore_jet(c,self.defects['B'],1)
            if endpoints(self.B[0])[1]>=0:raise ValueError('negative reshape B required')
            self.g=restore_jet(c,self.defects['g'],1);self.V=self.z*4+self.g
            self.m0={k:restore_jet(c,v,1) for k,v in self.source['physical_moments'].items()}
            self.u0=restore_jet(c,self.source['F'],1)*c.sqrt(220)
            self.I={k:restore_jet(c,v['jet'],1) for k,v in self.defects['flat_kernel_bounds'].items()}
            self.log_reference=-logjet(1+self.z*self.z)-c.mpf('5e151')
            self.T=fraction_box(c,T_EXACT)


class Md11LeadingField:
    def __init__(self):
        self.calc=Md11ComparisonJets(4);c=self.calc.ctx
        self.reshape=Md11LongReshapeField(self.calc)
        self.connecting=IntervalAxialRestoreField(self.reshape)
        self.inverse_name='lei_ren_part1_paper_Md11_five_bump_inverse.json'
        inverse=json.loads((HERE/self.inverse_name).read_bytes())
        validate_receipt(inverse,'new-source inverse')
        if (inverse['state_sha256']!=self.calc.state_hash or not inverse['certified']
            or not inverse['uniform_implicit_C1_family_exists']):
            raise ValueError('new-source uniform implicit inverse required')
        with mp.workdps(self.calc.precision+60):
            self.z=self.calc.z.truncate(1);self.P0=self.calc.p0.truncate(1)
            self.h=[restore_jet(c,v,1) for v in inverse['controls']]
            self.d=[restore_jet(c,self.reshape.defects['rows'][str(i)],1) for i in range(1,6)]
            self.A=(1+self.z*self.z).reciprocal()*c.exp(14)
            self.Am=self.A*c.exp(c.mpf('-.6'))
            self.Rm=110*c.exp(fraction_box(c,RM_LOG_EXACT));self.Rref=self.Rm*c.exp(6)
            self.map=IntervalPartialFiveBumpMap(c)
            endpoint=self.evaluate_outer_y(1)
            self.turnoff=IntervalOuterAxialTurnoff(c,self.z,self.calc.delta,endpoint,'1.1','1.1')

    def _source(self,data):
        data.update(new_schedule_sha256=self.calc.new_schedule_sha256,state_sha256=self.calc.state_hash,
            pressure_schedule_Md='1.1',pressure_datum_kind='analytic_preheat_stored_schedule',
            local_axial_family=self.calc.center_family,whole_axis=False,
            whole_outer_cone_certified=False,original_parameter_errors_enclosed=False,
            exact_heat_exterior_installed=False,temporal_recursion=False)
        return data

    def evaluate_connecting_offset(self,y):
        return self._source(self.connecting.evaluate_log_offset(y))

    def evaluate_repair_x(self,x):
        with mp.workdps(self.calc.precision+60):
            data=evaluate_reference(self.calc.ctx,self.z,self.calc.delta,self.Rm,self.Am,
                self.P0,self.h,self.d,x,self.map)
            data.update(implicit_inverse_family_enclosed=True,controls_projected_to_midpoints=False,
                terminal_moment_identities_implied_by_uniform_inverse=Fraction(x)==2)
            return self._source(data)

    def evaluate_reference_x(self,x):
        with mp.workdps(self.calc.precision+60):
            c=self.calc.ctx;q=Fraction(x);xx=fraction_box(c,q)
            if q<2 or endpoints(c.ln(xx))[1]>6:raise ValueError('2<=x<=exp(6) required')
            R=self.Rm*xx;u=self.Am*xx**c.mpf('.1')
            from lei_ren_part1_paper_interval_long_reshape_field import physical_packet
            data=physical_packet(c,self.z,self.calc.delta,R,u,u/10,self.z*4,
                reference_moments(c,self.z,R,u),self.P0,c.mpf('-.8'))
            return self._source(data)

    def evaluate_outer_y(self,y,cells=256):
        with mp.workdps(self.calc.precision+60):
            return self._source(evaluate_transition(self.calc.ctx,self.z,self.calc.delta,
                self.Rref,self.A,self.P0,y,cells))

    def evaluate_axial_phase(self,phase,cells=256):
        with mp.workdps(self.calc.precision+60):
            return self._source(self.turnoff.evaluate_phase(phase,cells))

    def evaluate_axial_buffer(self,offset,cells=256):
        with mp.workdps(self.calc.precision+60):
            return self._source(self.turnoff.evaluate_buffer(offset,cells))


def run():
    field=Md11LeadingField();c=field.calc.ctx
    with mp.workdps(field.calc.precision+60):
        samples=dict(connecting_terminal=field.evaluate_connecting_offset(RM_LOG_EXACT),
            repaired_terminal=field.evaluate_repair_x(2),reference=field.evaluate_reference_x(10),
            slope_endpoint=field.evaluate_outer_y(1),axial_midpoint=field.evaluate_axial_phase('.5'),
            axial_terminal=field.evaluate_axial_phase(1),buffer_Rd=field.evaluate_axial_buffer(11))
        end=samples['axial_terminal'];buffer=samples['buffer_Rd']
        if any(v._mpi_!=c.mpf(0)._mpi_ for key in ('Uz','Uz_y') for v in end[key].coefficients):
            raise ValueError('new physical axial terminal not exactly flat')
        if any(end['physical_moments']['z'][k]._mpi_!=buffer['physical_moments']['z'][k]._mpi_ for k in (0,1)):
            raise ValueError('new accumulated axial mass not retained')
        if any(v._mpi_!=field.P0[k]._mpi_ for data in samples.values() for k,v in enumerate(data['P0'].coefficients)):
            raise ValueError('P0 changed between new-source layers')
        result=dict(new_schedule_sha256=field.calc.new_schedule_sha256,state_sha256=field.calc.state_hash,
            center_family=field.calc.center_family,samples=samples,
            coherent_new_source_field_through_Rd_installed=True,
            exact_terminal_axial_cutoff=True,retained_Mz_verified=True,P0_unchanged_across_layers=True,
            terminal_local_implicit_five_moment_identity=True,
            new_whole_radial_cone_certified=False,new_heat_exterior_installed=False,
            original_parameter_errors_enclosed=False,whole_axis=False,temporal_recursion=False)
        names=(Path(__file__).name,'lei_ren_part1_paper_Md11_transition_pipeline.json',field.inverse_name,
            field.reshape.source_name,field.reshape.defects_name,'lei_ren_part1_paper_Md11_comparison_jets.py',
            'lei_ren_part1_paper_interval_long_reshape_field.py','lei_ren_part1_paper_interval_axial_restore_field.py',
            'lei_ren_part1_paper_interval_repaired_reference_field.py','lei_ren_part1_paper_interval_partial_five_bump_map.py',
            'lei_ren_part1_paper_interval_outer_slope_field.py','lei_ren_part1_paper_interval_outer_axial_turnoff.py')
        result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
        print('Coherent new Md1.1 leading field throughRd installed; axial cutoff and retained moments verified',flush=True)
        return result

if __name__=='__main__':run()
