"""Whole current main/exit strict cone from the original correlated source.

All original signed physical rows remain available. The analytic theorem
can resolve a cone even when the direct interval ratio loses correlation.
"""
import gzip
import json
from lei_ren_part1_paper_compliant_current_original_cone import (
    CurrentOriginalCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN,GATES as ORIGINAL_GATES)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_operator import (
    DOMAINS,generic_original_main_cone_theorem,current_main_cone_source_theorem,
    whole_current_main_exit_bounds)

NAME=PREFIX+'current_pulse_main_exit_cone.json'
RECEIPT=PREFIX+'current_pulse_main_exit_cone_check.json'
VIEWS_NAME=PREFIX+'current_pulse_main_exit_cone_views.json.gz'
GATES=('current_whole_pulse_main_signed_two_vector_cone_certified',
       'current_whole_pulse_exit_signed_two_vector_cone_certified',
       'current_pulse_main_exit_original_kernel_cone_reduction_available')


class CurrentPulseMainExitCone:
    @source_precision
    def __init__(self,cone=None,require_checked=True):
        self.cone=cone if cone is not None else CurrentOriginalCone()
        if type(self.cone) is not CurrentOriginalCone or not self.cone.acceptance_loaded:
            raise ValueError('Checked current original cone/source graph required')
        self.ctx=self.cone.ctx;self.family=self.cone.family;self.source=self.cone.source
        self.datum_sha=self.cone.datum_sha;self.registry=self.cone.registry
        self.assert_graph();self.theorem=generic_original_main_cone_theorem()
        self.bindings=current_main_cone_source_theorem(self.cone)
        self.whole_views={region:self.cone.native(
            region,(-1,1),limits,('-3','-1'),None,'1') for region,limits in DOMAINS.items()}
        self.bounds=whole_current_main_exit_bounds(self.cone,self.whole_views)
        self.hashes={**self.cone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_pulse_main_exit_cone_operator','current_pulse_main_exit_cone'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current main/exit cone receipt exceeds regional scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_main_exit_correlated_bounds']:
                raise ValueError('Current whole-domain bounds differ from checked producer')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def assert_graph(self):
        self.cone.assert_graph()
        if self.registry is not self.cone.registry or self.ctx is not self.cone.ctx:
            raise ValueError('Same original current cone/source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();view=self.cone.native(region,Z,coordinate,log_tau,theta,viscosity)
        eligible=region in DOMAINS
        return dict(region=region,original_complete_current_cone_view=view,
            regional_correlated_original_cone_bounds=self.bounds if eligible else None,
            cone_status='strict_by_whole_source_function_bound' if eligible else view.get('status',
                view.get('original_leading_two_vector_cone',{}).get('status')),
            strict_two_vector_cone_admitted=self.acceptance_loaded and eligible,
            direct_interval_margin_preserved_even_when_inconclusive=True,
            physical_positive_lambda_and_constant_viscosity_factor_cancels=True,
            scope='Whole current main/exit two-vector cone; full signed tensor preserved. No global completed tensor or actual wave lift.',
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:self.acceptance_loaded and eligible and (k!=GATES[0] or region=='pulse_main')
                and (k!=GATES[1] or region=='pulse_exit') for k in GATES},
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,original_generic_signed_main_exit_cone_theorem=self.theorem,
            current_actual_main_exit_source_composition=self.bindings,
            current_whole_main_exit_correlated_bounds=self.bounds,
            original_complete_signed_whole_views=VIEWS_NAME,whole_domain_source_function_proof=True,
            actual_homogeneous_covariance_pulses_still_unconstructed=True,
            scope='Whole current main xi[.02,10] and exit xi[10,11], all Z[-1,1], positive lambda and constant nu. Current pulse-end admission inherited separately. Other regions, uniform edges, global completed tensor, waves, recursion, flatness, NS and energy remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseMainExitCone(require_checked=False)
    result=field.manifest()
    payload=(json.dumps(encode(pack(field.whole_views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole current main/exit original signed cone constructed; global waves/recursion remain open',flush=True)
    return result


if __name__=='__main__':run()
