"""Whole current flatten/outer-power strict cone, full histories and joins."""
import gzip
import json
from lei_ren_part1_paper_compliant_current_pulse_gap_cone import (
    CurrentPulseGapCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN)
from lei_ren_part1_paper_compliant_current_flatten_power_cone_operator import (
    DOMAINS,generic_flatten_power_theorem,current_flatten_power_source_theorem,
    whole_current_flatten_power_bounds,sigma_jets)

NAME=PREFIX+'current_flatten_power_cone.json'
RECEIPT=PREFIX+'current_flatten_power_cone_check.json'
VIEWS_NAME=PREFIX+'current_flatten_power_cone_views.json.gz'
GATES=('current_whole_flatten_signed_two_vector_cone_certified',
    'current_whole_outer_power_signed_two_vector_cone_certified',
    'current_flatten_power_full_history_shear_and_three_source_joins_available')


class CurrentFlattenPowerCone:
    @source_precision
    def __init__(self,gapcone=None,require_checked=True):
        self.gapcone=gapcone if gapcone is not None else CurrentPulseGapCone()
        if type(self.gapcone) is not CurrentPulseGapCone or not self.gapcone.acceptance_loaded:
            raise ValueError('Checked current gap cone and complete same source graph required')
        self.registry=self.gapcone.registry;self.ctx=self.gapcone.ctx
        self.family=self.gapcone.family;self.source=self.gapcone.source;self.datum_sha=self.gapcone.datum_sha
        self.assert_graph();self.theorem=generic_flatten_power_theorem()
        self.bindings=current_flatten_power_source_theorem(self.gapcone)
        self.whole_views={region:self.registry.native(region,(-1,1),domain,('-3','-1'),None,'1')
            for region,domain in DOMAINS.items()}
        self.bounds=whole_current_flatten_power_bounds(self.gapcone,self.whole_views)
        self.hashes={**self.gapcone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_flatten_power_cone_operator','current_flatten_power_cone'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current flatten/power cone receipt exceeds regional scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_flatten_power_correlated_bounds']:
                raise ValueError('Actual whole current flatten/power bounds differ from receipt')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.gapcone.assert_graph()
        if self.registry is not self.gapcone.registry or self.ctx is not self.gapcone.ctx:
            raise ValueError('Same current checked complete tensor source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();eligible=region in DOMAINS;shear=None
        if eligible:
            view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity)
            c=self.ctx;owner=self.registry.owners['flatten'];mu=c.mpf(owner.outer.mu)
            rho=c.ln((1+c.mpf(Z)**2)/2)
            slope=sigma_jets(c,c.mpf(coordinate)/100)[1]/100 if region=='flatten' else c.mpf(0)
            m=2*mu-2*rho*slope
            shear=dict(a=2+m,bs=c.mpf(0),vs_minus2=m,
                exact_axial_velocity_zero_from_current_terminal_and_full_FTC=True,
                source_inverse_radius_strictly_positive=True,source_caps_define_field_values=False)
        else:view=self.gapcone.native(region,Z,coordinate,log_tau,theta,viscosity)
        return dict(region=region,original_complete_current_source_view=view,current_source_shear=shear,
            regional_correlated_original_cone_bounds=self.bounds['regions'][region] if eligible else None,
            strict_flatten_power_two_vector_cone_admitted=self.acceptance_loaded and eligible,
            cone_status='strict_by_whole_current_correlated_history_bound' if eligible else 'inherited_region_scope',
            scope='Current whole flatten/power two-vector cone; full E/P/X/K and three joins. Remaining global regions, actual waves and coefficient recursion open.',
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:self.acceptance_loaded and eligible and (k!=GATES[0] or region=='flatten')
                and (k!=GATES[1] or region=='outer_power') for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,original_generic_flatten_power_cone_theorem=self.theorem,
            current_actual_flatten_power_source_and_join_theorem=self.bindings,
            current_whole_flatten_power_correlated_bounds=self.bounds,original_complete_signed_whole_views=VIEWS_NAME,
            scope='Whole current flatten offset[0,100] and outer_power phase[0,1], all Z[-1,1], finite positive lambda and constant nu. Full forward X, variable K_Z, complete remaining E/P and absolute pressure retained. Current pulse-end/flatten, flatten/power and power/angular functional tensor joins retained. Five preceding pulse cone admissions inherited separately. Global completed tensor, uniform edges, actual waves, coefficient recursion, flatness, NS and energy remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentFlattenPowerCone(require_checked=False)
    result=field.manifest()
    payload=(json.dumps(encode(pack(field.whole_views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole flatten/power strict cone constructed; remaining regions/waves/recursion open',flush=True)
    return result


if __name__=='__main__':run()
