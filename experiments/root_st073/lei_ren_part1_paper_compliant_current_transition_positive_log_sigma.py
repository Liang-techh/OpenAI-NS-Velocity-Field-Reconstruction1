"""Original positive logarithmic cutoff and genuine transition collar C1.

Refines the same sigma(t), Delta=2mu*sigma(t), a=2+Delta and b=0 source;
the selected eta, defining q, native radius phase and first-jet formulas stay
unchanged. No accepted source module or prior numerical archive is mutated.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_transition_right_half as previous

generic=previous.query_source;tail=generic.tail;rc=generic.rc;common=generic.common
first=rc.density.first;slow=first.slow;prior=slow.prior
HERE,PREFIX,sha=generic.HERE,generic.PREFIX,generic.sha;N=generic.N
ep,iv,encode=generic.ep,generic.iv,generic.encode;KEYS=generic.KEYS;ZERO=(0,0)
NAME=PREFIX+'current_transition_positive_log_sigma.json'
RECEIPT=PREFIX+'current_transition_positive_log_sigma_check.json'
GATE='original_positive_log_sigma_transition_left_collar_C1_executed'
COLLAR_LEFT=Fraction(1,2**128)


def original_active_support_proof(coordinates,mu,eta_log):
    """Source-derived support width; no parameter or field value selected.

On0<t<=1/2, sigma(t)>=exp(-1/t^2)/2. Thus Delta>=mu*exp(-1/t^2)
and q=0 whenever t>=D^(-1/2), D=log(mu)-log(eta)>0.
The whole parameter cover stays directed; its upper width bounds support.
"""
    c=coordinates.ctx;D=c.ln(mu)-eta_log
    if ep(D)[0]<=0:raise ArithmeticError('Original positive logarithmic mu/eta gap required')
    logwidth=-c.ln(D)/2
    width=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,offset=logwidth),1,coordinates.ledger)
    return dict(original_selected_mu=mu,original_selected_eta_log=eta_log,
        positive_original_log_mu_over_eta_cover=D,
        original_q_active_support_width=width.record(),
        original_support_theorem='0<t<=1/2: L(t)>=-1/t^2, 1+exp(L)<=2, Delta=2mu*sigma>=mu*exp(-1/t^2); t>=1/sqrt(log(mu)-log(eta)) implies Delta>=eta and original q/all slow jets=0',
        finite_entire_source_parameter_covers_retained=True,
        support_width_upper_bounds_the_original_possible_nonzero_q_domain=True,
        original_t_zero_active_endpoint_not_deleted=True,
        native_coordinate_support_width_not_vortex_core_or_project_completion_percentage=True,
        no_tiny_width_exponential_or_huge_integer_cycles_materialized=True)


def positive_log_sigma(template,t):
    """Same original logistic sigma and ordinary y0/y1/y2 source rows."""
    c=template.ctx;lo,hi=ep(c.mpf(t));scalar=template.scalar
    if lo<0 or hi>ep(c.mpf('.5'))[1]:raise ValueError('Original left cutoff source interval in[0,1/2] required')
    if lo==hi==0:
        values={order:scalar(0) for order in slow.ORDERS}
        return dict(values=values,record=dict(branch='exact_original_endpoint_zero',coordinate=c.mpf(t),
            original_sigma_and_all_ordinary_derivatives_exact_zero=True))
    if lo<=0:
        coefficients=prior.sigma_jets(c,c.mpf(t))
        values={order:scalar(coefficients[order[0]]*math.factorial(order[0]) if not order[1] else 0)
            for order in slow.ORDERS}
        return dict(values=values,record=dict(branch='original_flat_endpoint_tail_cover',coordinate=c.mpf(t),
            nonzero_lower_not_claimed_on_interval_containing_zero=True,
            original_ordinary_rows={str(k):v.record() for k,v in values.items()}))
    def odds(x):return 1/(1-x)**2-1/x**2
    L=c.mpf((ep(odds(c.mpf(lo)))[0],ep(odds(c.mpf(hi)))[1]))
    assert ep(L)[1]<=0
    bounded_exponential=template.bounded_exp(L)
    inv=1/(1+bounded_exponential)
    sigma=prior.ScaledEnclosure(prior.FormalScale(template.scale.bases,offset=L),inv,template.ledger)
    # 1-sigma = 1/(1+exp(L)); retain its defining factor, avoiding the
    # cancellation of a positive source lower bound in 1-sigma.
    complement=scalar(inv);x=c.mpf(t)
    L1=2/(1-x)**3+2/x**3;L2=6/(1-x)**4-6/x**4
    logistic_product=sigma*complement
    ordinary=(sigma,logistic_product*L1,logistic_product*(L2+(2*inv-1)*L1**2))
    values={order:ordinary[order[0]] if not order[1] else scalar(0) for order in slow.ORDERS}
    record=dict(branch='original_positive_logistic_source',coordinate=x,
        original_logistic_odds='L(t)=1/(1-t)^2-1/t^2; sigma=exp(L)/(1+exp(L))',
        monotone_original_odds_interval=L,positive_denominator_inverse=inv,
        bounded_exponential_used_only_in_positive_denominator=bounded_exponential,
        original_complement_factor=complement.record(),ordinary_odds_derivatives=dict(y1=L1,y2=L2),
        original_sigma_ordinary_y2_Z1_rows={str(k):v.record() for k,v in values.items()},
        positive_source_factor_retained_even_when_exponential_cannot_be_materialized=True,
        cutoff_source_function_not_changed=True,source_or_sigma_cap_point_not_selected=True,
        ordinary_derivatives_not_Taylor_coefficients=True,all_original_Z_sigma_derivatives_exact_zero=True)
    return dict(values=values,record=record)


class PositiveLogTransition:
    def __init__(self,owner):
        self.oracle=generic.GenericTailPhaseIntegrals(owner);self.owner=owner
        self.c=owner.ctx;self.coordinates=owner.coordinates;self.family=owner.family
        self.density=self.oracle.density_owner;self.first=self.density.owner
        self.signed=self.first.owner.owner.owner
        self.positive=self.signed.decode(self.signed.inventory['O3_slope_mu']['actual_positive_denominator_theorem'])
        self.mu=iv(self.c,self.positive['actual_positive_mu'])
        self.eta=iv(self.c,self.signed.scales['selected_positive_eta_log'])
        self.dstar=self.first.dstar

    def query(self,Z,left,right):
        geometry=self.oracle.geometry('O3_slope_mu',Fraction(left),Fraction(right))
        old=self.density.spatial_query('O3_slope_mu',Z,geometry['coordinate'],N)
        original=old['source'];roots=dict(original['roots']);base=roots['E'][ZERO]
        sigma=positive_log_sigma(base,geometry['coordinate'])
        delta={k:v*(2*self.mu) for k,v in sigma['values'].items()}
        assert all(v.zero for v in roots['b'].values())
        roots.update(kappa_minus2=delta,a={k:v+(2 if k==ZERO else 0) for k,v in delta.items()})
        roots['t0']={k:base.scalar(0) for k in slow.ORDERS}
        loop=slow.current.q_enclosure(roots['a'][ZERO],delta[ZERO],self.eta,self.positive['log_actual_a_positive_lower'])
        qjet=slow.original_q_jet(roots,self.eta,self.positive['log_actual_a_positive_lower'],loop)
        source_record=dict(original['record'],
            correlated_shear_and_signed_root_enclosures={name:{'y%d_Z%d'%order:v.record() for order,v in rows.items()} for name,rows in roots.items()},
            original_q_enclosure={k:v.record() if isinstance(v,prior.ScaledEnclosure) else v for k,v in loop.items()},
            same_original_positive_log_sigma_source_refinement=True,
            original_Delta_equals_two_mu_sigma_and_a_equals_two_plus_Delta=True)
        refined=dict(original,roots=roots,q=loop['q'],loop=loop,record=source_record)
        qrecord={k:v for k,v in qjet.items() if k!='rows'}
        qrecord.update(source_family=self.family,chart='O3_slope_mu',source_provenance=original['packet'].provenance,
            original_correlated_shear_and_q=source_record,
            current_refined_original_excess_ordinary_rows={'y%d_Z%d'%k:v.record() for k,v in delta.items()},
            current_refined_original_q_C0=loop['q'].record(),
            original_q_ordinary_slow_derivative_enclosures=None if qjet['rows'] is None else {'y%d_Z%d'%k:v.record() for k,v in qjet['rows'].items()},
            original_definition_sigma_eta_q_unchanged=True,
            original_selected_mu=self.mu,original_selected_eta_log=self.eta,
            original_actual_a_positive_lower_log=self.positive['log_actual_a_positive_lower'])
        # Same native normalized axial velocity and ordinary Z conversion.
        packet=original['packet']
        leaf=lambda k:prior.signed.expressions.RadiusPolynomial(packet.algebra,
            {0:prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)})
        V,VZ=[self.signed.leaf(leaf(k),base.scale.bases,base.ledger) for k in (0,1)]
        cells=[];values=[]
        for phi in old['geometry']['phase_boxes']:
            got=first.conditioned_first_jets(refined,qjet['rows'],self.dstar,phi)
            chain=None
            if got['values'] is not None:
                primitive=got['values'];chain=dict(A_y_total=primitive['A_y']+N*primitive['A_phi'],
                    B_y_total_over_Pstar=primitive['B_y_over_Pstar']+N*primitive['B_phi_over_Pstar'],
                    A_Z_total=primitive['A_Z'],B_Z_total_over_Pstar=primitive['B_Z_over_Pstar'])
                density=rc.density.density_Z_kernels(base,roots['E'][(0,1)],V,VZ,primitive,N)
                values.append(density)
            fr=dict(got['record'],derived_spatial_fractional_phase_box=phi,phase_is_independent_input_on_this_query=False,
                actual_spatial_fast_N_chain_installed=chain is not None,
                original_total_spatial_first_derivative_enclosures=None if chain is None else {k:v.record() for k,v in chain.items()})
            row=dict(status=fr['status'],original_phase_first_jet_source=fr)
            if got['values'] is not None:
                row.update(original_and_candidate_normalized_velocity_Z_enclosures={k:v.record() for k,v in density['velocities'].items()},
                    five_original_signed_density_C0_enclosures={k:v.record() for k,v in density['kernels'].items()},
                    five_original_signed_density_Z_derivative_enclosures={k:v.record() for k,v in density['Z_derivatives'].items()},
                    original_radius_phase_Z_derivative_exactly_zero=True,same_original_factor_basis_and_ledger=True,
                    source_native_width_or_Pstar_conversion_not_reapplied=True)
            cells.append(row)
        enclosed=len(values)==len(cells)
        native=dict(source_family=self.family,chart='O3_slope_mu',source_provenance=packet.provenance,candidate_N=N,
            actual_original_radius_phase=old['geometry']['record'],original_q_slow_jet_source=qrecord,
            actual_spatial_signed_density_Z_cells=cells,five_signed_density_Z_functions_installed=enclosed)
        record=dict(chart='O3_slope_mu',exact_native_endpoints=[str(Fraction(left)),str(Fraction(right))],candidate_N=N,
            source_family=self.family,exact_Z_range=list(Z),actual_true_geometry=geometry['record'],
            actual_original_phase_inverse_and_density_source=native,
            status='enclosed' if enclosed else 'requires_original_source_phase_subdivision',
            actual_N_log_radius_phase_not_independent_phase_samples=True,whole_period_primitive_caps_not_used_as_A_B=True,
            original_ordinary_Z_density_derivatives_not_cap_derivatives=True,
            same_original_positive_log_sigma_refinement=sigma['record'],
            original_unrefined_q_slow_source=old['record']['original_q_slow_jet_source'],
            original_native_source_log_bases=base.scale.bases)
        if not enclosed:return dict(record=record,values=None)
        kernels={k:rc.transfer.local.same_source_union([v['kernels'][k] for v in values]) for k in KEYS}
        jets={k:rc.transfer.local.same_source_union([v['Z_derivatives'][k] for v in values]) for k in KEYS}
        factors={k:rc.transfer.true_width_kernel(self.coordinates,geometry,rate) for k,rate in rc.RATES.items()}
        adds={k:self.coordinates.rebase(kernels[k],self.family)*factors[k]['mass'] for k in KEYS}
        derivs={k:self.coordinates.rebase(jets[k],self.family)*factors[k]['mass'] for k in KEYS}
        record.update(full_phase_union_signed_C0_density={k:v.record() for k,v in kernels.items()},
            full_phase_union_signed_Z_density={k:v.record() for k,v in jets.items()},
            true_positive_kernel_factors={k:dict(branch=f['branch'],mass=f['mass'].record(),decay=f['decay'].record()) for k,f in factors.items()},
            original_local_C0_integral={k:v.record() for k,v in adds.items()},original_local_Z_integral={k:v.record() for k,v in derivs.items()},
            true_log_radius_Jacobian_integrated_once=True,pressure_P0_not_added_to_density_or_integral=True)
        return dict(record=record,geometry=geometry,values=adds,Z_derivatives=derivs)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/previous.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,previous,accepted['source_family'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(tail.baseline.build_owner(bridge)))
        query=PositiveLogTransition(owner);archives=[];summaries=[]
        support=original_active_support_proof(owner.coordinates,query.mu,query.eta)
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            probes=[]
            for lo,hi in ((0,COLLAR_LEFT),(COLLAR_LEFT,Fraction(1,2))):
                got=query.query(Z,lo,hi);probes.append(got['record'])
                print('Original positive-log transition source:',tag,str(lo),str(hi),got['record']['status'],
                    got['record']['actual_original_phase_inverse_and_density_source']['original_q_slow_jet_source']['branch'],flush=True)
            collar=got;assert collar['values'] is not None,'Original positive-log collar must be enclosed'
            assert all(v.zero for v in (*collar['values'].values(),*collar['Z_derivatives'].values()))
            operator=rc.history.C1DuhamelOperator(owner.coordinates)
            rc.transfer.append_true_cell(operator,collar['geometry'],collar['values'],collar['Z_derivatives'],owner.family)
            old_archive=next(a for a in accepted['actual_original_transition_right_half_source_archives'] if a['exact_Z_range']==list(Z))
            right=json.loads(gzip.decompress((HERE/old_archive['filename']).read_bytes()))['original_entire_transition_right_half_source_query']
            g=right['actual_true_geometry'];geometry=dict(width=common.restore_common_source(g['positive_true_log_radius_width'],owner.coordinates),
                regular=iv(owner.ctx,g['regular_true_log_radius_width_cover']),scalar_cover=iv(owner.ctx,g['scalar_width_cover_used_only_for_directed_kernel_bounds']))
            add={k:common.restore_common_source(v,owner.coordinates) for k,v in right['original_local_C0_integral'].items()}
            jets={k:common.restore_common_source(v,owner.coordinates) for k,v in right['original_local_Z_integral'].items()}
            rc.transfer.append_true_cell(operator,geometry,add,jets,owner.family)
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                original_left_half_positive_log_sigma_queries=probes,
                accepted_right_half_source_archive=old_archive,
                original_combined_quiet_transition_window=[str(COLLAR_LEFT),'1'],exact_true_log_radius_width=str(1-COLLAR_LEFT),
                original_combined_quiet_transition_C1_operator=operator.record(),
                common_directed_coordinate_theorem=owner.coordinates.record(),
                unresolved_original_transition_prefix_not_skipped=True,
                original_mu_eta_sigma_q_and_source_family_unchanged=True,
                entire_transition_or_tail_numerical_integrals_admitted=False)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_transition_positive_log_sigma_'+tag+'.json.gz'
            (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
            archives.append(dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),
                uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()))
            summaries.append(dict(exact_Z_range=list(Z),query_statuses=[r['status'] for r in probes],
                combined_quiet_native_window=[str(COLLAR_LEFT),'1'],own_C0_Z_integral_increments_exact_zero=True))
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Source closures disagree: '+name)
            hashes[name]=digest
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            actual_original_positive_log_sigma_source_archives=archives,actual_transition_collar_summaries=summaries,
            original_source_derived_nonzero_q_support=support,
            same_original_source_logistic_sigma_and_ordinary_y2_Z1_refined=True,
            original_combined_quiet_transition_C1_window=[str(COLLAR_LEFT),'1'],
            original_transition_unresolved_prefix=['0',str(COLLAR_LEFT)],
            original_source_seed_evidence=seed,input_hashes=hashes,execution_seconds=time.monotonic()-began,
            complete_original_transition_or_tail_numerical_integrals_or_Rc_targets_admitted=False,
            **dict.fromkeys(common.current.FLAGS,False),
            scope='Same original positive-log sigma and ordinary y2/Z1 source on both strict-sign tiles. Actual original left collar[2^-128,1/2] flat C1 integral joined with accepted right half gives exact quiet operator[2^-128,1]. Original prefix/axial/targets/controls/global N/stress/recursion/full NS remain open.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
