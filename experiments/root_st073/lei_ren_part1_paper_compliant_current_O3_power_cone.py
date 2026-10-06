"""Current whole O3 power cone from original correlated histories."""
import gzip
import json
from lei_ren_part1_paper_compliant_current_pulse_entrance_cone import (
    CurrentPulseEntranceCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN)
from lei_ren_part1_paper_compliant_current_O3_power_cone_operator import (
    DOMAIN,original_correlated_O2_O3_theorem,current_O3_power_source_theorem,
    whole_current_O3_power_bounds)

NAME=PREFIX+'current_O3_power_cone.json'
RECEIPT=PREFIX+'current_O3_power_cone_check.json'
VIEWS_NAME=PREFIX+'current_O3_power_cone_views.json.gz'
GATES=('current_whole_O3_power_signed_two_vector_cone_certified',
    'current_O3_power_actual_M_K_X_full_energy_pressure_correlation_bound',
    'current_O3_power_transition_and_entrance_function_joins_consumed')


class CurrentO3PowerCone:
    @source_precision
    def __init__(self,entrancecone=None,require_checked=True):
        self.entrancecone=entrancecone if entrancecone is not None else CurrentPulseEntranceCone()
        if type(self.entrancecone) is not CurrentPulseEntranceCone or not self.entrancecone.acceptance_loaded:
            raise ValueError('Checked current entrance cone required')
        self.registry=self.entrancecone.registry;self.ctx=self.entrancecone.ctx
        self.family=self.entrancecone.family;self.source=self.entrancecone.source;self.datum_sha=self.entrancecone.datum_sha
        self.assert_graph();self.theorem=original_correlated_O2_O3_theorem()
        self.bindings=current_O3_power_source_theorem(self.entrancecone)
        self.whole_view=self.registry.native('O3_power',(-1,1),DOMAIN,('-3','-1'),None,'1')
        pre=self.registry.owners['incoming'].physical.pre
        self.slope=pre.slope(0,1);self.initial=pre.power(0,0)
        self.bounds=whole_current_O3_power_bounds(self.entrancecone,self.whole_view,self.slope,self.initial)
        self.hashes={**self.entrancecone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_O3_theta_correlation','current_O3_power_cone_operator','current_O3_power_cone'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current O3 power cone receipt exceeds regional scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_O3_power_correlated_bounds']:
                raise ValueError('Actual current O3 power bounds differ from receipt')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.entrancecone.assert_graph()
        if self.registry is not self.entrancecone.registry or self.ctx is not self.entrancecone.ctx:
            raise ValueError('Same checked complete current source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();eligible=region=='O3_power';shear=None
        if eligible:
            view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity)
            c=self.ctx;mu=c.mpf(self.registry.owners['incoming'].pulse.mu)
            velocity=view['original_complete_view']['current_source_three_component_velocity_rows']
            if any(endpoints(value)!=(0,0) for row in velocity['axial'] for value in row.coefficients):
                raise ValueError('Actual O3 axial velocity/shear differs from canonical source')
            shear=dict(a=2+2*mu,bs=c.mpf(0),vs_minus2=2*mu,exact_a_minus2_source=2*mu,
                full_current_velocity_local_and_incoming=velocity,
                axial_velocity_and_shear_zero_by_actual_source_identity=True,
                full_incoming_radial_velocity_and_original_stress_retained=True,
                positive_viscosity_and_lambda_cancel_from_common_source_factor=True)
        else:view=self.entrancecone.native(region,Z,coordinate,log_tau,theta,viscosity)
        return dict(region=region,original_complete_current_source_view=view,current_source_shear=shear,
            regional_correlated_original_cone_bounds=self.bounds if eligible else None,
            strict_O3_power_two_vector_cone_admitted=self.acceptance_loaded and eligible,
            independent_signed_sector_boxes_not_subtracted_to_test_cone=True,
            cone_status='strict_by_current_whole_O3_power_source_bound' if eligible else 'inherited_region',
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:self.acceptance_loaded and eligible for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_original_correlated_O2_O3_theorem=self.theorem,
            current_actual_O3_power_full_source_composition=self.bindings,
            current_whole_O3_power_correlated_bounds=self.bounds,complete_current_signed_source_views=VIEWS_NAME,
            current_strict_nonzero_regions_including_inherited=15,current_exact_zero_exterior_regions=1,
            remaining_registry_regions_without_current_cone=17,
            scope='Whole actual O3 power phase[0,1], all Z[-1,1], finite positive lambda and constant nu. Actual shared O2 history and nonnegative transition integral correlated before enclosure; complete theta M/K/X, inverse-R shear, full energy and same signed absolute pressure retained. Transition/power and power/entrance function joins consumed. Fourteen strict downstream regions and exact zero heat exterior inherited separately. Variable O3 transition, O2/inner/global cones, completed global tensor/lift, actual waves, n-dependent recursion, temporal flatness, corrected NS and energy remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3PowerCone(require_checked=False)
    result=field.manifest();views=dict(whole_current_O3_power=field.whole_view,
        actual_O2_slope_endpoint=field.slope,actual_O3_power_phase0=field.initial)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole O3 power cone constructed; 15 nonzero regions, variable transition/recursion/waves open',flush=True)
    return result


if __name__=='__main__':run()
