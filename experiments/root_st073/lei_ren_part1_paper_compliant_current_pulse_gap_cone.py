"""Current whole gap/gap-end strict cone with real signed full histories."""
import gzip
import json
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone import (
    CurrentPulseMainExitCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN)
from lei_ren_part1_paper_compliant_current_pulse_gap_cone_operator import (
    DOMAINS,generic_original_gap_cone_theorem,current_gap_cone_source_theorem,whole_current_gap_bounds)

NAME=PREFIX+'current_pulse_gap_cone.json'
RECEIPT=PREFIX+'current_pulse_gap_cone_check.json'
VIEWS_NAME=PREFIX+'current_pulse_gap_cone_views.json.gz'
GATES=('current_whole_pulse_gap_signed_two_vector_cone_certified',
       'current_whole_reciprocal_gap_end_signed_two_vector_cone_certified',
       'current_gap_actual_zero_axial_shear_and_full_history_cone_reduction_available')


class CurrentPulseGapCone:
    @source_precision
    def __init__(self,maincone=None,require_checked=True):
        self.maincone=maincone if maincone is not None else CurrentPulseMainExitCone()
        if type(self.maincone) is not CurrentPulseMainExitCone or not self.maincone.acceptance_loaded:
            raise ValueError('Checked current main/exit cone and original source graph required')
        self.registry=self.maincone.registry;self.ctx=self.maincone.ctx
        self.family=self.maincone.family;self.source=self.maincone.source;self.datum_sha=self.maincone.datum_sha
        self.assert_graph();self.theorem=generic_original_gap_cone_theorem()
        self.bindings=current_gap_cone_source_theorem(self.maincone,self.theorem)
        self.whole_views={region:self.registry.native(region,(-1,1),limits,('-3','-1'),None,'1')
            for region,limits in DOMAINS.items()}
        self.bounds=whole_current_gap_bounds(self.maincone,self.whole_views)
        self.hashes={**self.maincone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_pulse_gap_cone_operator','current_pulse_gap_cone'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current gap cone receipt exceeds regional scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_gap_correlated_bounds']:
                raise ValueError('Current whole gap bounds differ from checked producer')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.maincone.assert_graph()
        if self.registry is not self.maincone.registry or self.ctx is not self.maincone.ctx:
            raise ValueError('Same current gap/main/exit source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();eligible=region in DOMAINS
        if eligible:
            view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity)
            velocity=view['original_complete_view']['current_source_three_component_velocity_rows']
            if any(endpoints(row[n])!=(0,0) for row in velocity['axial'] for n in range(6)):
                raise ValueError('Current native gap shear is not structurally zero')
            mu=self.ctx.mpf(self.registry.owners['gap'].pulse.mu)
            shear=dict(a=2+2*mu,bs=self.ctx.mpf(0),vs_minus2=2*mu,
                actual_zero_axial_velocity_jets=velocity['axial'],nonzero_radial_history_preserved=True)
        else:view=self.maincone.native(region,Z,coordinate,log_tau,theta,viscosity);shear=None
        return dict(region=region,original_complete_current_source_view=view,
            current_source_shear=shear,regional_correlated_original_cone_bounds=self.bounds if eligible else None,
            strict_gap_two_vector_cone_admitted=self.acceptance_loaded and eligible,
            cone_status='strict_by_whole_current_gap_function_bound' if eligible else 'inherited_region_scope',
            scope='Current whole gap/gap-end two-vector cone and complete signed tensor; global completed tensor and actual waves remain open.',
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:self.acceptance_loaded and eligible and (k!=GATES[0] or region=='pulse_gap')
                and (k!=GATES[1] or region=='pulse_gap_end') for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,original_generic_gap_signed_cone_theorem=self.theorem,
            current_actual_gap_source_and_three_attachment_theorem=self.bindings,
            current_whole_gap_correlated_bounds=self.bounds,original_complete_signed_whole_views=VIEWS_NAME,
            scope='Whole current gap xi[11,12] and reciprocal gap-end phase[0,1], all Z[-1,1], positive lambda and finite positive constant nu. Actual D0/D1/D2, pressure, full weights and nonzero radial histories remain. Main/exit/end regional admissions inherited separately. Global completed tensor, uniform edges, actual waves, coefficient recursion, flatness, NS and energy remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseGapCone(require_checked=False)
    result=field.manifest()
    payload=(json.dumps(encode(pack(field.whole_views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole gap/gap-end signed cone constructed; remaining region/wave/recursion scope open',flush=True)
    return result


if __name__=='__main__':run()
