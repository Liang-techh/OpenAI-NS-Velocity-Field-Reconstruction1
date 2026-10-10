"""Source-bound current Rp functions in the actual native interval inlet.

Reuse accepted directed enclosures of the now-identified defining functions.
This adapter supplies true six-row local inlet data, never point values.
The selected pulse and downstream physical assembly remain separate.
"""
import copy
import json
from pathlib import Path
import time
import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_Rp_native_constants as exact
import lei_ren_part1_paper_compliant_power_inlet_C4 as original
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=exact.HERE,exact.PREFIX,exact.sha
NAME=PREFIX+'current_original_Rp_interval_inlet.json'
RECEIPT=PREFIX+'current_original_Rp_interval_inlet_check.json'
GATES=('current_original_Rp_source_bound_interval_constant_view_installed',
       'current_original_Rp_local_native_power_interval_callback_installed')
PENDING=('native_interval_inlet_callback_installed',
         'selected_native_pulse_constructor_consumes_current_C3_frame',
         'current_numeric_point_field_oracle_installed')+exact.OPEN
REUSED=PREFIX+'power_inlet_C4_check.json'
FRAME_KEYS=frozenset(('u','m1','m2','X','energy','Mp','P0','pressure'))


def add_hashes(target,source):
    for name,digest in source.items():
        if sha(name)!=digest or (name in target and target[name]!=digest):
            raise ValueError('Different inlet defining source: '+name)
        target[name]=digest


class CurrentOriginalRpIntervalInlet(original.CompliantPowerInletC4):
    @source_precision
    def __init__(self,constants=None,native=None,require_checked=True):
        self.exact=constants if constants is not None else exact.CurrentOriginalRpNativeConstants()
        if not self.exact.acceptance_loaded:raise ValueError('Accepted current/native constant identities required')
        self.before=native if native is not None else original.CompliantPowerInletC4()
        if type(self.before) is not original.CompliantPowerInletC4:
            raise TypeError('Original typed source enclosure constructor required')
        family=self.exact.identity
        if (self.before.family,self.before.source,self.before.datum.datum_sha)!=tuple(family[key] for key in
                ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')):
            raise ValueError('Current Rp/native interval family differs')
        definition=self.exact.frame.reports[exact.frame.native.NAME]['current_P0_source_binding']['analytic_definition']
        if self.before.datum.definition!=definition:raise ValueError('Independent analytic P0 definition differs')
        proof=json.loads((HERE/REUSED).read_bytes())
        if not (proof['all_passed'] and proof['two_sided_O3_pulse_join_certified']
                and proof['exact_functional_production_and_join_identities'][
                    'single_sample_only_encloses_proved_Z_independent_Hp_Pin_constants']):
            raise ValueError('Original source-enclosure and canonical-shape proof required')
        self.family_record=family;self.family=self.before.family;self.source=self.before.source
        self.datum_sha=self.before.datum.datum_sha;self.hashes=dict(self.exact.hashes)
        add_hashes(self.hashes,self.before.hashes);add_hashes(self.hashes,proof['input_hashes'])
        for name in (REUSED,exact.NAME,exact.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.ctx=c=self.before.ctx
        # Preserve full directed boxes, their tails and their defining recipe.
        # The exact identity, not numerical overlap, permits their current use.
        self.constants=dict(self.before.constants)
        self.mu=self.before.mu;self.delta=self.before.delta
        self.rate=self.before.rate;self.prate=self.before.prate
        self.inlet_H=self.before.inlet_H;self.inlet_P=self.before.inlet_P;self.Xp=self.before.Xp
        self.datum=self.before.datum;self.radial=copy.copy(self.before.radial)
        self.interval_constants={key:self.constants[key] for key in exact.KEYS if key in self.constants}
        self.interval_constants.update(H=self.inlet_H,Pin=self.inlet_P,Xp=self.Xp)
        if set(self.interval_constants)!=exact.KEYS:raise ValueError('Complete current interval constant view required')
        if self.constants is self.before.constants or self.radial is self.before.radial:
            raise ValueError('Separate mutable inlet views required')
        self.calls=[];self._caller='direct interval callback';self.acceptance_loaded=False
        self.original_constant_ids={key:id(value) for key,value in self.before.constants.items()}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES):
                raise ValueError('Accepted current interval caller handshake required')
            add_hashes(self.hashes,receipt['input_hashes'])
            if receipt['source_family']!=family or any(receipt[key] for key in PENDING):
                raise ValueError('Current interval receipt family/scope differs')
            self.acceptance_loaded=True

    @source_precision
    def incoming(self,Z):
        result=super().incoming(Z)
        result['pressure']=result['P0']+result['Mp']
        if set(result)!=FRAME_KEYS or any(jet.order!=5 or jet.ctx is not self.ctx for jet in result.values()):
            raise ValueError('All eight true six-row interval inlet functions required')
        self.calls.append(dict(Z=endpoints(self.ctx.mpf(Z)),rows=48,
            caller=self._caller,source_receipt=exact.RECEIPT))
        return result

    @source_precision
    def current_power(self,Z,distance=0):
        before=len(self.calls)
        previous=self._caller;self._caller='CompliantPowerInletC4.power'
        try:packet=super().power(Z,distance)
        finally:self._caller=previous
        if len(self.calls)!=before+1:raise ValueError('Actual native power caller must invoke current callback')
        return packet

    def ownership(self):
        return dict(current_function_receipt_consumed=self.exact.acceptance_loaded,
            separate_constant_dictionary=self.constants is not self.before.constants,
            separate_radial_view=self.radial is not self.before.radial,
            original_interval_objects_retained=all(self.constants[key] is value for key,value in self.before.constants.items()),
            original_constant_objects_unmutated=self.original_constant_ids=={key:id(value) for key,value in self.before.constants.items()},
            independent_P0_owner_retained=self.datum is self.before.datum,
            actual_native_power_method_retained=self.power.__func__ is original.CompliantPowerInletC4.power,
            source_enclosures_used_without_point_selection=True)


@source_precision
def run(constants=None,native=None):
    began=time.monotonic();owner=CurrentOriginalRpIntervalInlet(constants,native,require_checked=False)
    views={}
    for name,Z in (('axis','0'),('fresh','.371'),('cell',['.370','.372'])):
        incoming=owner.incoming(Z)
        packet=owner.current_power(Z,0)
        views[name]=dict(Z=owner.ctx.mpf(Z),incoming=incoming,actual_native_power_packet=packet)
    result=dict(source_family=owner.family_record,candidate_current_interval_inlet_constructed=True,
        current_interval_constant_functions=owner.interval_constants,
        source_bound_six_row_views=views,actual_native_caller_invocations=owner.calls,
        actual_current_inlet_owner_graph=owner.ownership(),
        unchanged_td_Tw_B_squared_mass_and_far_tail_retained=True,
        exact_current_native_constant_identity_receipt=exact.RECEIPT,
        original_directed_source_enclosure_receipt=REUSED,
        **dict.fromkeys(GATES+PENDING,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8',newline='\n')
    print('Current interval Rp inlet: actual native power caller consumed six-row source functions',flush=True)
    return owner


if __name__=='__main__':run()
