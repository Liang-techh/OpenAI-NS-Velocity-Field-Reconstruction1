"""Checked current angular/tail cone, including the exact zero heat edge."""
import gzip
import json
from lei_ren_part1_paper_compliant_current_flatten_power_cone import (
    CurrentFlattenPowerCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN)
from lei_ren_part1_paper_compliant_current_angular_tail_cone_operator import (
    DOMAINS,generic_angular_tail_theorem,current_angular_tail_source_theorem,
    whole_current_angular_tail_bounds,sigma_jets)

NAME=PREFIX+'current_angular_tail_cone.json'
RECEIPT=PREFIX+'current_angular_tail_cone_check.json'
VIEWS_NAME=PREFIX+'current_angular_tail_cone_views.json.gz'
GATES=('current_whole_angular_entry_power_exit_waiting_signed_cone_certified',
    'current_whole_heat_collar_nonzero_signed_cone_certified',
    'current_heat_zero_edge_uniform_direction_and_unbounded_exterior_available',
    'current_angular_tail_full_histories_and_seven_source_joins_available')


class CurrentAngularTailCone:
    @source_precision
    def __init__(self,flattencone=None,require_checked=True):
        self.flattencone=flattencone if flattencone is not None else CurrentFlattenPowerCone()
        if type(self.flattencone) is not CurrentFlattenPowerCone or not self.flattencone.acceptance_loaded:
            raise ValueError('Checked current flatten/power cone required')
        self.registry=self.flattencone.registry;self.ctx=self.flattencone.ctx
        self.family=self.flattencone.family;self.source=self.flattencone.source;self.datum_sha=self.flattencone.datum_sha
        self.assert_graph();self.theorem=generic_angular_tail_theorem()
        self.bindings=current_angular_tail_source_theorem(self.flattencone)
        self.whole_views={region:self.registry.native(region,(-1,1),domain,('-3','-1'),None,'1')
            for region,domain in DOMAINS.items()}
        self.bounds=whole_current_angular_tail_bounds(self.flattencone,self.whole_views)
        self.unbounded_exterior=self.registry.owners['heat'].unbounded_exterior()
        self.hashes={**self.flattencone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_angular_tail_cone_operator','current_angular_tail_cone'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current angular/tail cone receipt exceeds scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_angular_tail_correlated_bounds']:
                raise ValueError('Actual current angular/tail bounds differ from receipt')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.flattencone.assert_graph()
        if self.registry is not self.flattencone.registry or self.ctx is not self.flattencone.ctx:
            raise ValueError('Same complete current checked source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();eligible=region in DOMAINS;shear=None;strict=False;zero=False;partition=None
        if eligible:
            view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity)
            c=self.ctx;h=self.registry.owners['heat'].heat;mu=h.mu;delta=h.delta;a=h.a;k=h.k
            x=c.mpf(coordinate);z=c.mpf(Z)
            if region=='outer_angular':
                native=self.registry.owners['angular'].outer.angular(z,x)
                F=native['swirl_factor_one_plus_h_Taylor'][0]
                m=2*mu-2*native['actual_angular_bump_y_derivatives'][1][0]/F
            elif region=='steep_entry':m=2*mu+2*(1-mu)*sigma_jets(c,x)[0]
            elif region=='steep_power':m=c.mpf(2)
            elif region=='steep_exit':m=delta+2*k*(1-sigma_jets(c,x)[0])
            elif region=='waiting':m=delta
            else:
                source=h.shape(z,x) if region=='heat_collar' else h.local_Gamma(z,x)
                m=delta-2*source['K_rows'][1][0]/source['K_rows'][0][0]
            strict=region!='heat_exterior' and (region!='heat_collar' or endpoints(x)[1]<3)
            zero=region=='heat_exterior' or (region=='heat_collar' and endpoints(x)==(3,3))
            if region=='heat_collar':
                lower,upper=endpoints(x);pieces=[]
                if lower<=1:pieces.append(dict(proof='sigmoid_part',query_intersection=c.mpf([lower,min(upper,1)])))
                if upper>1 and lower<3:pieces.append(dict(proof='heat_part',query_intersection=c.mpf([max(lower,1),upper]),strict_upper_t3_excluded=True))
                partition=dict(subdomain_proofs=pieces,same_original_source_at_t1=True,
                    strict_proof_composed_across_t1=True,zero_t3_separate=True)
            shear=dict(a=2+m,vs_minus2=m,bs=c.mpf(0),
                actual_axial_velocity_zero_from_current_terminal_and_full_FTC=True,
                exact_source_inverse_radius_positive=True,source_caps_define_field_values=False)
        else:view=self.flattencone.native(region,Z,coordinate,log_tau,theta,viscosity)
        return dict(region=region,original_complete_current_source_view=view,current_source_shear=shear,
            current_heat_collar_proof_partition=partition,
            regional_original_cone_bounds=self.bounds['regions'][region] if eligible else None,
            strict_angular_tail_two_vector_cone_admitted=self.acceptance_loaded and strict,
            regional_cone_theorem_strict_where_nonzero=self.acceptance_loaded and eligible and region!='heat_exterior',
            exact_heat_zero_tensor_admitted=self.acceptance_loaded and zero,
            heat_collar_boxes_containing_t3_not_strict=True,
            cone_status=('exact_zero_heat_source' if zero else 'strict_by_current_continuous_source_bound' if strict else 'scoped_theorem_or_inherited_region'),
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{GATES[0]:self.acceptance_loaded and region in ('outer_angular','steep_entry','steep_power','steep_exit','waiting'),
               GATES[1]:self.acceptance_loaded and region=='heat_collar',
               GATES[2]:self.acceptance_loaded and region in ('heat_collar','heat_exterior'),
               GATES[3]:self.acceptance_loaded and eligible},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,original_generic_angular_tail_cone_theorem=self.theorem,
            current_actual_angular_tail_source_theorem=self.bindings,
            current_whole_angular_tail_correlated_bounds=self.bounds,
            current_unbounded_exact_zero_exterior=self.unbounded_exterior,
            original_complete_signed_whole_views=VIEWS_NAME,
            current_strict_nonzero_regions_including_inherited=13,
            current_exact_zero_exterior_regions=1,remaining_registry_regions_without_current_cone=19,
            scope='Current whole angular s[-4,0], entry/power/exit/waiting phase[0,1], collar t[0,3) and exact zero edge/exterior t>=3, all Z[-1,1], finite positive lambda and constant nu. Original full A/E/P/K and pressure retained, seven adjacent plus four internal functional tensor joins. Seven preceding cone regions inherited separately. Strict cone excludes zero edge. Uniform heat edge direction is regional. Global completed tensor/lift, actual oscillatory fields, n-dependent recursion, temporal flatness, NS and energy remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentAngularTailCone(require_checked=False)
    result=field.manifest();payload=(json.dumps(encode(pack(field.whole_views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current angular/tail cone constructed; 13 nonzero regions plus exact zero exterior, global waves/recursion open',flush=True)
    return result


if __name__=='__main__':run()
