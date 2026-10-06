"""Checked variable transition direction and uninstalled shear-loop construction."""
import gzip
import json
import mpmath as mp
from lei_ren_part1_paper_compliant_current_O3_power_cone import (
    CurrentO3PowerCone,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,accepted,_verify_hashes,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_current_O3_transition_direction_operator import (
    DOMAIN,actual_variable_transition_theorem,current_variable_transition_bindings,whole_current_transition_bounds)

NAME=PREFIX+'current_O3_transition_direction.json'
RECEIPT=PREFIX+'current_O3_transition_direction_check.json'
VIEWS_NAME=PREFIX+'current_O3_transition_direction_views.json.gz'
GATES=('current_whole_O3_variable_transition_direction_bound_certified',
    'current_O3_transition_actual_full_energy_pressure_and_variable_jets_bound',
    'current_O3_transition_nonzero_zero_shear_endpoint_obstruction_certified',
    'current_O3_transition_mean_preserving_periodic_shear_loop_candidate_certified')
CLOSED_OPEN=('current_whole_O3_transition_signed_two_vector_cone_certified',
    'current_O3_transition_finite_N_modified_profiles_installed',
    'current_O3_transition_modified_five_moment_repair_certified')
OPEN=EARLIER_OPEN+CLOSED_OPEN


class CurrentO3TransitionDirection:
    @source_precision
    def __init__(self,powercone=None,require_checked=True):
        self.powercone=powercone if powercone is not None else CurrentO3PowerCone()
        if type(self.powercone) is not CurrentO3PowerCone or not self.powercone.acceptance_loaded:raise ValueError('Checked actual current O3 power cone required')
        self.registry=self.powercone.registry;self.ctx=self.powercone.ctx
        self.family=self.powercone.family;self.source=self.powercone.source;self.datum_sha=self.powercone.datum_sha
        self.assert_graph();self.theorem=actual_variable_transition_theorem();self.bindings=current_variable_transition_bindings(self.powercone)
        self.whole_view=self.registry.native('O3_slope_mu',(-1,1),DOMAIN,('-3','-1'),None,'1')
        self.inlet=self.registry.owners['o3'].physical.pre.axial(0,buffer_offset=11)
        self.bounds=whole_current_transition_bounds(self.powercone,self.whole_view,self.inlet)
        self.hashes={**self.powercone.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_O3_transition_direction_operator','current_O3_transition_direction'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Variable direction receipt exceeds scope')
            producer=json.loads((HERE/NAME).read_bytes())
            if encode(pack(self.bounds))!=producer['current_whole_variable_transition_bounds']:raise ValueError('Actual transition bounds differ')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.powercone.assert_graph()
        if self.registry is not self.powercone.registry or self.ctx is not self.powercone.ctx:raise ValueError('Same checked source graph required')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1',loop_phase=None):
        self.assert_graph();eligible=region=='O3_slope_mu';shear=None;loop=None
        if eligible:
            view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity);c=self.ctx
            v=c.mpf(coordinate);mu=c.mpf(self.registry.owners['o3'].pulse.mu);sig=sigma_jets(c,v)[0]
            shear=dict(a=2+2*mu*sig,bs=c.mpf(0),vs_minus2=2*mu*sig,
                exact_variable_shear_excess='2*mu*sigma(offset)',original_cutoff_value=sig,
                zero_shear_excess_at_exact_offset0=endpoints(v)==(0,0),
                strict_shear_positive_by_source_on_open_offset=endpoints(v)[0]>0)
            if loop_phase is not None:
                phi=c.mpf(loop_phase)
                if not all(mp.isfinite(x) for x in endpoints(phi)):raise ValueError('Finite periodic phase required')
                dA=mu*c.cos(4*c.pi*phi);b=2*c.sqrt(mu)*c.sin(2*c.pi*phi)
                velocity=view['original_complete_view']['current_source_three_component_velocity_rows']
                E0=velocity['theta'][0];lp=c.mpf(self.registry.owners['o3'].physical.logP)
                loop=dict(phase=phi,a=2+2*mu*sig+dA,a_minus2=2*mu*sig+dA,b=b,
                    periodic_A=-mu*c.sin(4*c.pi*phi)/(8*c.pi),
                    periodic_B_signed_coefficient=E0*(-c.sqrt(mu)*c.cos(2*c.pi*phi)/(2*c.pi)),
                    periodic_B_source_logPstar=lp,actual_original_E0_raw_function=E0,
                    exact_average_shear_preserved=True,
                    original_inviscid_vector_retained_with_changed_stress_target=True,
                    current_source_loop_bounds=self.bounds['source_correlated_periodic_shear_loop'],
                    finite_N_physical_profile_and_moment_repair_installed=False)
        else:view=self.powercone.native(region,Z,coordinate,log_tau,theta,viscosity)
        return dict(region=region,original_complete_current_source_view=view,current_source_variable_shear=shear,
            regional_variable_direction_bounds=self.bounds if eligible else None,current_source_periodic_shear_loop=loop,
            whole_closed_transition_cone_admitted=False,
            open_offset_strict_two_vector_cone_by_current_source=self.acceptance_loaded and eligible and endpoints(self.ctx.mpf(coordinate))[0]>0,
            cone_status=('strict_only_on_open_offset; offset0_requires_original_shear_repair' if eligible else 'inherited_region'),
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:self.acceptance_loaded and eligible for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_variable_transition_and_periodic_shear_theorem=self.theorem,
            current_actual_full_variable_source_bindings=self.bindings,current_whole_variable_transition_bounds=self.bounds,
            complete_current_signed_source_views=VIEWS_NAME,current_strict_nonzero_whole_regions_including_inherited=15,
            current_exact_zero_exterior_regions=1,remaining_registry_regions_without_current_whole_cone=17,
            scope='Complete actual variable transition direction/positive theta on closed offset[0,1], Z[-1,1]. Full variable logU jets, all histories, original absolute functional pressure and remainder retained. Strict cone on offset(0,1] only. At offset0 original vs-2=0 while stress is nonzero. A current-source mean-preserving periodic shear loop and primitives are constructed, using the same original inviscid vector and changed loop stress. Spatial support, finite uniform N profiles, independent five-moment restoration and affected joins are not installed. Whole-region count stays15; global cone/lift, actual waves, coefficient recursion, flatness, NS/energy open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3TransitionDirection(require_checked=False)
    result=field.manifest();views=dict(whole_current_variable_transition=field.whole_view,actual_current_O2_buffer_inlet=field.inlet)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual variable transition direction and shear-loop constructed; closed cone still needs finite N/moment repair',flush=True)
    return result


if __name__=='__main__':run()
