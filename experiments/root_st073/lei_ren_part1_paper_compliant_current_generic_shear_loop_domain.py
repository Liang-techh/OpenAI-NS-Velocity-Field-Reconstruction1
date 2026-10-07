"""Current original two-sided loop collars and later new-repair geometry.

Microscopic offsets stay formal even when adding them to logRa would lose
them. The two-sided original collar query does not rebuild an ancestor graph.
Whole modification-box input margins and actual loop scales remain separate.
"""
import ast
import json
from pathlib import Path
from types import SimpleNamespace
import lei_ren_part1_paper_compliant_current_generic_shear_O3_sources as source

packets=source.packets
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_generic_shear_loop_domain.json'
RECEIPT=PREFIX+'current_generic_shear_loop_domain_check.json'
GATE='current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'
OPEN=source.OPEN


def decode(c,value):
    if isinstance(value,dict):
        if 'lower' in value and 'upper' in value:return packets.interval(c,value)
        return {k:decode(c,v) for k,v in value.items()}
    if isinstance(value,list):return [decode(c,v) for v in value]
    return value


def two_sided_query(c,proof,family):
    """Replay the original query alone with the checked original collar proof."""
    asts=packets.recovery.numeric.transport.SourceAST()
    fn=asts.method('current_inner_exit_strict_collar','query')
    strict=[kw.value for node in ast.walk(fn) if isinstance(node,ast.Call) for kw in node.keywords
        if kw.arg=='strict_nonzero_stress_cone_certified_for_entire_query_box']
    if len(strict)!=1 or ast.dump(strict[0])!=ast.dump(ast.parse('lo>0',mode='eval').body):
        raise ValueError('Original whole-query strict collar rule changed')
    import lei_ren_part1_paper_compliant_current_inner_exit_strict_collar as original
    env=dict(vars(original))
    query=packets.compile_function(fn,env,'<original collar query only>')
    owner=SimpleNamespace(ctx=c,proof=proof,family=family)
    result=query(owner,('.25','.75'),(-1,1))
    if not result['strict_nonzero_stress_cone_certified_for_entire_query_box']:
        raise ArithmeticError('Two-sided original left collar unavailable')
    return result,dict(passed=True,original_query_AST_replayed_without_constructor=True,
        strict_fraction_lower_positive_rule_bound=True,input_hashes=asts.hashes)


class CurrentLoopDomain:
    def __init__(self,provider=None):
        self.provider=provider if provider is not None else source.CurrentO3Sources()
        self.service=self.provider.service;self.ctx=self.provider.ctx;self.family=self.provider.family
        for stem,gate in (('current_generic_shear_O3_sources',source.GATE),
            ('current_inner_exit_strict_collar','current_inner_exit_strict_collar_attached_to_current_source_graph_certified')):
            name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate):raise ValueError('Checked original domain prerequisite required')
            family=receipt.get('source_family') or {k:receipt[k] for k in packets.FAMILY_KEYS}
            if family!=self.family:raise ValueError('Original loop collar/source/pressure family differs')
            self.service.bind_hashes(receipt['input_hashes']);self.service.bind_hashes({name:sha(name)})
        original=PREFIX+'current_inner_exit_strict_collar.json'
        self.left=json.loads((HERE/original).read_bytes());self.service.bind_hashes({original:sha(original)})
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})

    def domain(self):
        c=self.ctx;read=lambda v:packets.interval(c,v);endpoints=packets.recovery.endpoints
        proof=decode(c,self.left['explicit_current_inner_exit_strict_collar'])
        query,theorem=two_sided_query(c,proof,self.family['actual_five_defect_family_sha256'])
        self.service.bind_hashes(theorem['input_hashes'])
        sc=proof['selected_first_phase_endpoint'];logh=proof['exact_source_width_log_enclosure']
        gamma=proof['positive_margins']['positive_gamma']
        reserve=self.provider.right_reservation();mu=self.provider.denominator('O3_power')['actual_positive_mu']
        norm=json.loads((HERE/(PREFIX+'physical_norm_family.json')).read_bytes())
        logC=read(norm['selected_logCstar']);logP=self.service.data['logP']
        # hb*s_c<log(4.1/4) is already proved for this exact width/source.
        offset_cap=c.ln(c.mpf('4.1')/4)/2
        gap_lower=c.ln(c.mpf(110)/4)+10*logC+15*logP+1002-offset_cap
        left_to_R110=c.ln(c.mpf(110)/4)+4*logP+1000-offset_cap
        for value in (sc,gamma,mu,logC,logP,gap_lower,left_to_R110):
            if endpoints(value)[0]<=0:raise ArithmeticError('Original loop radius/collar ordering unresolved')
        if endpoints(proof['positive_margins']['subcollar_inside_original_rho4point1_continuation'])[0]<=0:
            raise ArithmeticError('Original positive microscopic left-radius cap required')
        left_excess=c.mpf('1.9')*gamma;right_excess=2*mu
        logboundary=c.mpf(min(endpoints(c.ln(left_excess))[0],endpoints(c.ln(right_excess))[0]))
        Tw=reserve['log_radius_offsets_from_same_Rw']['Rp']
        if endpoints(Tw-c.mpf('1.25'))[0]<=0:raise ArithmeticError('Strict right collar must lie in original power chart')
        return dict(source_family=self.family,
            formal_radii=dict(Ra='4*epsilon_core',r_minus='Ra*exp(hb*s_c/2)',
                r_plus='Rw*exp(1)',Rc='Rw*exp(2)',twice_Rc='Rw*exp(2+log2)',Rb='later original exterior boundary, not selected here'),
            exact_log_Ra_without_lost_offset=c.ln(4)-4*logP-1000,
            exact_log_Rw=c.ln(110)+10*logC+11*logP+1,
            formal_left_log_radius_offset=dict(source='hb*s_c/2',log_positive_offset=logh+c.ln(sc)-c.ln(2),
                bound_above=offset_cap,not_added_to_huge_log_Ra=True,positive_width_not_materialized=True),
            positive_log_r_plus_over_r_minus_lower=c.mpf(endpoints(gap_lower)[0]),
            positive_log_R110_over_r_minus_lower=c.mpf(endpoints(left_to_R110)[0]),
            left_two_sided_original_source_query=query,left_query_source_theorem=theorem,
            left_center_fraction=c.mpf('.5'),left_collar_fraction_box=c.mpf(['.25','.75']),
            right_center_original_power_offset=c.mpf(1),right_collar_original_power_offset_box=c.mpf(['.75','1.25']),
            right_collar_mu=mu,right_collar_kappa_minus2=right_excess,
            left_collar_kappa_minus2_lower=left_excess,
            both_boundary_kappa_excess_log_lower=logboundary,
            allowed_eta_log_upper_for_automatic_constant_edges=logboundary-c.ln(2),
            boundary_log_eta_constraint_is_not_the_full_loop_scale_selection=True,
            original_Ra_to_r_minus_relaxed_admission=dict(
                receipt=PREFIX+'global_exit_certificate.json',domain='open Ra<R<=r_minus, subset of original Ra<R<=110',
                attached_same_source_theorem=self.left['exact_current_exit_source_attachment'],
                Ra_endpoint='separate exact stress-free core case; not a nonzero strict point'),
            original_r_plus_through_twice_Rc_admission_and_reservation=reserve,
            selected_geometry_Ra_r_minus_r_plus_Rc_twice_Rc_Rp_strict=True,
            original_both_loop_edge_collars_strict=True,
            independent_spatial_multiplier_or_taper_not_introduced=True,
            old_one_sided_support_metadata_not_used_as_two_sided_collar=True,
            complete_modification_box_H0_minus2_certified=False,
            original_outer_admission_twice_Rc_through_Rb_certified=False,
            complete_Section11_loop_domain_certified=False,whole_generic_scales_instantiated=False,
            **dict.fromkeys(OPEN,False))

    def run(self):
        domain=self.domain()
        result=dict(source_family=self.family,current_original_generic_loop_domain=domain,
            **{GATE:True},**dict.fromkeys(OPEN,False),source_graph_ancestor_constructors_called=False,
            scope='Actual original two-sided strict loop-edge collars, left radius ordering and right/new-repair reservation. Full modification-box H0/scales and original post-repair admission remain open.',
            input_hashes=self.service.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        print('Original generic loop two-sided collars and ordered new repair geometry PASS',flush=True)
        return result


def run():return CurrentLoopDomain().run()


if __name__=='__main__':run()
