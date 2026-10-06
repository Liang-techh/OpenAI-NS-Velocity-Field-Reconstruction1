"""Checked whole entrance cone, preserving all current source histories."""
import gzip
import json
from lei_ren_part1_paper_compliant_current_angular_tail_cone import (
    CurrentAngularTailCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN)
from lei_ren_part1_paper_compliant_current_original_cone import normalized_signed_pulse,cone_margins
from lei_ren_part1_paper_compliant_current_pulse_entrance_cone_operator import (
    DOMAIN,entrance_memory_theorem,current_entrance_source_theorem,whole_current_entrance_bounds)

NAME=PREFIX+'current_pulse_entrance_cone.json'
RECEIPT=PREFIX+'current_pulse_entrance_cone_check.json'
VIEWS_NAME=PREFIX+'current_pulse_entrance_cone_views.json.gz'
GATES=('current_whole_pulse_entrance_signed_two_vector_cone_certified',
    'current_entrance_actual_pre_power_memory_and_full_histories_bound',
    'current_entrance_incoming_and_main_source_joins_consumed')


class CurrentPulseEntranceCone:
    @source_precision
    def __init__(self,tailcone=None,require_checked=True):
        self.tailcone=tailcone if tailcone is not None else CurrentAngularTailCone()
        if type(self.tailcone) is not CurrentAngularTailCone or not self.tailcone.acceptance_loaded:
            raise ValueError('Checked current angular/tail cone required')
        self.registry=self.tailcone.registry;self.ctx=self.tailcone.ctx
        self.family=self.tailcone.family;self.source=self.tailcone.source;self.datum_sha=self.tailcone.datum_sha
        self.assert_graph();self.theorem=entrance_memory_theorem()
        self.bindings=current_entrance_source_theorem(self.tailcone)
        self.whole_view=self.registry.native('pulse_entrance',(-1,1),DOMAIN,('-3','-1'),None,'1')
        self.preTw=self.registry.owners['incoming'].physical.pre.power(0,0)
        self.bounds=whole_current_entrance_bounds(self.tailcone,self.whole_view,self.preTw)
        self.hashes={**self.tailcone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_pulse_entrance_cone_operator','current_pulse_entrance_cone'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current entrance cone receipt exceeds regional scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_entrance_correlated_bounds']:
                raise ValueError('Actual current entrance bounds differ from receipt')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.tailcone.assert_graph()
        if self.registry is not self.tailcone.registry or self.ctx is not self.tailcone.ctx:
            raise ValueError('Same checked complete current source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();eligible=region=='pulse_entrance';shear=None;reduced=None;tested=None
        if eligible:
            view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity)
            c=self.ctx;owner=self.registry.owners['incoming'];mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta)
            velocity=view['original_complete_view']['current_source_three_component_velocity_rows']
            C=velocity['local']['theta'][0][0]
            if endpoints(C)[0]<=0:raise ValueError('Actual positive local theta coefficient required')
            bs=2*(velocity['local']['axial'][1][0]+velocity['incoming']['axial'][1][0])/C
            a=2+2*mu;km=2*mu+bs**2/a
            shear=dict(a=a,bs=bs,vs_minus2=km,exact_a_minus2_source=2*mu,
                full_current_velocity_local_and_incoming=velocity,full_axial_shear_retained=True,
                incoming_radial_velocity_and_stress_cross_terms_retained=True,
                positive_viscosity_and_lambda_cancel_from_common_source_factor=True)
            reduced=normalized_signed_pulse(c,view,mu,delta)
            tested=cone_margins(c,a,bs,reduced['signed_totals']['theta'],reduced['signed_totals']['axial'],km)
        else:view=self.tailcone.native(region,Z,coordinate,log_tau,theta,viscosity)
        return dict(region=region,original_complete_current_source_view=view,current_source_shear=shear,
            normalized_current_signed_stress=reduced,direct_interval_cone_diagnostic=tested,
            regional_correlated_original_cone_bounds=self.bounds if eligible else None,
            strict_entrance_two_vector_cone_admitted=self.acceptance_loaded and eligible,
            direct_interval_diagnostic_not_used_as_continuous_proof=True,
            cone_status='strict_by_current_whole_entrance_source_bound' if eligible else 'inherited_region',
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:self.acceptance_loaded and eligible for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_original_power_memory_theorem=self.theorem,
            current_actual_entrance_full_source_composition=self.bindings,
            current_whole_entrance_correlated_bounds=self.bounds,complete_current_signed_source_views=VIEWS_NAME,
            current_strict_nonzero_regions_including_inherited=14,current_exact_zero_exterior_regions=1,
            remaining_registry_regions_without_current_cone=18,
            scope='Whole current entrance xi[0,.02], all Z[-1,1], finite positive lambda and constant nu. Actual pre-Tw memory, current incoming moments/energy/absolute pressure, full signed stress and radial velocity retained. Incoming/entrance and entrance/main function joins consumed. Thirteen strict downstream regions and exact zero heat exterior inherited separately. O3/O2/core upstream cones, global completed tensor/lift, actual oscillatory fields, n-dependent recursion, temporal flatness, NS and energy remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseEntranceCone(require_checked=False)
    result=field.manifest();views=dict(whole_current_entrance=field.whole_view,actual_pre_Tw_power_source=field.preTw)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole entrance cone constructed; 14 nonzero regions, upstream recursion/waves open',flush=True)
    return result


if __name__=='__main__':run()
