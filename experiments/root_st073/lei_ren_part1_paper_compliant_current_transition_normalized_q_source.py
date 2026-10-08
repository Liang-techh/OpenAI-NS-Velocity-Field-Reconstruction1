"""Same-original microscopic O3 coordinates and correlated ordinary q jets.

xi=s*sqrt(D) resolves the bulk; k=1/s**2-D resolves the cutoff seam.
Whole original parameter covers and true nonzero widths are retained.
This is a source/q/geometry component, not a complete velocity integral.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_transition_positive_log_sigma as previous

rc,common,prior,slow=previous.rc,previous.common,previous.prior,previous.slow
first=previous.first;current=slow.current
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
ep,iv,encode=previous.ep,previous.iv,previous.encode
ZERO=(0,0);ORDERS=slow.ORDERS;N=previous.N
NAME=PREFIX+'current_transition_normalized_q_source.json'
RECEIPT=PREFIX+'current_transition_normalized_q_source_check.json'
GATE='original_transition_typed_microscopic_coordinates_and_correlated_q_y2_Z1_executed'


class OriginalTransitionCoordinates:
    """Typed endpoint expressions; scalar covers never define endpoints."""
    def __init__(self,coordinates,mu,eta_log):
        self.coordinates=coordinates;self.c=c=coordinates.ctx
        self.mu=c.mpf(mu);self.eta_log=c.mpf(eta_log)
        self.D=c.ln(self.mu)-self.eta_log
        if ep(self.D)[0]<=4:raise ValueError('Original D>4 support scale required')
        self.logD=c.ln(self.D)
        self.W=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,offset=-self.logD/2),1,coordinates.ledger)
        # Compact mpf exponents, not an enormous Fraction denominator or exp(-D).
        self.W_cover=1/c.sqrt(self.D)
        self.parameter_record=dict(original_mu=self.mu,original_eta_log=self.eta_log,
            original_D_expression='log(mu)-log(eta)',entire_original_D_cover=self.D,
            original_W_expression='(log(mu)-log(eta))**(-1/2)',
            positive_original_W=self.W.record(),directed_W_cover=self.W_cover,
            entire_parameter_covers_retained=True,parameter_midpoint_or_endpoint_selected=False,
            finite_mpf_exponent_tuple_not_huge_rational_denominator=True)

    def cv(self,q):
        q=Fraction(q);return self.c.mpf(q.numerator)/q.denominator

    def geometry(self,kind,left,right):
        c=self.c;lo,hi=Fraction(left),Fraction(right)
        if kind=='xi':
            if not 0<=lo<hi<=1:raise ValueError('Original ordered xi endpoints in[0,1]')
            a=self.W_cover*self.cv(lo);b=self.W_cover*self.cv(hi)
            width=self.W*self.cv(hi-lo)
            specs=[dict(kind='xi',exact_multiplier=str(q),source_expression='s=W*xi') for q in (lo,hi)]
            expression='s=xi/sqrt(D); dy/dxi=W; original dy/ds=1'
            x=c.mpf((ep(self.cv(lo))[0],ep(self.cv(hi))[1]))
            inverse_square=self.D/x**2 if lo else None
        elif kind=='k':
            if not 0<=hi<lo:raise ValueError('Original increasing s requires decreasing k>=0')
            a=1/c.sqrt(self.D+self.cv(lo));b=1/c.sqrt(self.D+self.cv(hi))
            # Exact rationalized endpoint difference, shared D collected first.
            ra=c.sqrt(1+self.cv(lo)/self.D);rb=c.sqrt(1+self.cv(hi)/self.D)
            width=prior.ScaledEnclosure(prior.FormalScale(self.coordinates.bases,offset=-3*self.logD/2),
                self.cv(lo-hi)/(ra*rb*(ra+rb)),self.coordinates.ledger)
            specs=[dict(kind='k',exact_offset=str(q),source_expression='s=(D+k)**(-1/2)') for q in (lo,hi)]
            expression='k=1/s**2-D; s=(D+k)**(-1/2); ds/dk=-(D+k)**(-3/2)/2'
            inverse_square=self.D+c.mpf((ep(self.cv(hi))[0],ep(self.cv(lo))[1]))
        else:raise ValueError('Original xi or k source coordinate required')
        coordinate=c.mpf((ep(a)[0],ep(b)[1]))
        if ep(coordinate)[0]<0 or ep(coordinate)[1]>ep(c.mpf('.5'))[0]:raise ValueError('Original left half required')
        if ep(width.coefficient)[0]<=0:raise ArithmeticError('Positive original typed width required')
        # Only exp(log(width)) is materialized: its exponent tuple is finite
        # and compact. exp(-D), with an unallocatable exponent, is never used.
        scalar_width=width.coefficient*c.exp(width.scale.evaluate())
        record=dict(chart='O3_slope_mu',source_family=self.coordinates.family,
            typed_coordinate_kind=kind,exact_normalized_endpoints=[str(lo),str(hi)],
            original_endpoint_specification=specs,original_native_coordinate_expression=expression,
            original_parameter_source=self.parameter_record,native_coordinate_box=coordinate,
            positive_true_log_radius_width=width.record(),regular_true_log_radius_width_cover=c.mpf(0),
            scalar_width_cover_used_only_for_directed_kernel_bounds=scalar_width,
            endpoint_difference_collected_symbolically_before_numeric_subtraction=True,
            ordinary_source_rows_are_y_derivatives_no_second_coordinate_Jacobian=True,
            integration_Jacobian_installed_once_in_true_width=True,width_and_endpoints_independent_of_Z=True,
            source_endpoints_not_scalar_parameter_selections=True)
        return dict(record=record,kind=kind,left=lo,right=hi,coordinate=coordinate,
            width=width,regular=c.mpf(0),scalar_cover=scalar_width,inverse_square=inverse_square)

    def normalized_excess(self,template,geometry):
        """Delta/eta and ordinary y jets, with original D cancellation."""
        c=self.c;scalar=template.scalar;s=geometry['coordinate']
        if ep(s)[0]<=0:
            return dict(status='requires_open_positive_prefix_source',rows=None,
                record=dict(exact_original_s_zero_endpoint_retained=True,
                    whole_interval_starting_at_zero_not_assumed_active=True))
        F=1/(1-s)**2
        if geometry['kind']=='xi':
            x=c.mpf((ep(self.cv(geometry['left']))[0],ep(self.cv(geometry['right']))[1]))
            residual=self.D*(1-1/x**2)
            formula='log(Delta/eta)=log(2)+D*(1-xi**(-2))+1/(1-s)**2-log(1+exp(L))'
        else:
            k=c.mpf((ep(self.cv(geometry['right']))[0],ep(self.cv(geometry['left']))[1]))
            residual=-k
            formula='log(Delta/eta)=log(2)+1/(1-s)**2-k-log(1+exp(L)); L=1/(1-s)**2-D-k'
        L=F-geometry['inverse_square']
        if ep(L)[1]>=0:raise ArithmeticError('Original negative logistic odds required')
        tail=template.bounded_exp(L);complement=1/(1+tail)
        G=c.ln(2)+residual+F
        value=prior.ScaledEnclosure(prior.FormalScale(template.scale.bases,offset=G),complement,template.ledger)
        L1=2/(1-s)**3+2/s**3;L2=6/(1-s)**4-6/s**4
        ordinary=(value,value*complement*L1,value*complement*(L2+(2*complement-1)*L1**2))
        rows={k:ordinary[k[0]] if not k[1] else scalar(0) for k in ORDERS}
        record=dict(status='enclosed',same_original_normalized_Delta_source=formula,
            original_logistic_odds_cover=L,shared_D_terms_cancelled_before_evaluation=True,
            normalized_excess_log_numerator_cover=G,
            positive_original_logistic_denominator_inverse=complement,
            nonzero_exponential_tail_not_set_to_zero=tail,
            ordinary_odds_derivatives=dict(y1=L1,y2=L2),
            original_normalized_Delta_over_eta_ordinary_rows={'y%d_Z%d'%k:v.record() for k,v in rows.items()},
            all_original_normalized_excess_Z_rows_exact_zero=True,
            no_xi_or_k_derivative_substituted_for_ordinary_y_derivative=True,
            original_sigma_mu_eta_and_q_functions_unchanged=True)
        return dict(status='enclosed',rows=rows,record=record)


def correlated_original_q(normalized,eta_log,log_a_lower):
    """Original lazy q and smooth seam, where the original root stays real."""
    v=normalized[ZERO];c=v.ctx;scalar=v.scalar
    eta=prior.ScaledEnclosure(prior.FormalScale(v.scale.bases,offset=eta_log),1,v.ledger)
    Delta={k:eta*row for k,row in normalized.items()}
    a={k:row+(2 if k==ZERO else 0) for k,row in Delta.items()}
    # A normalized interval may have a bounded upper log and an arbitrarily
    # negative lower log. Keep its directed tail; finite_interval would
    # incorrectly require both log endpoints to be scalar-materializable.
    d=c.mpf(0) if v.zero else v.coefficient*v.bounded_exp(v.scale.evaluate())
    dlo,dhi=ep(d)
    if dlo>=1:
        rows={k:scalar(0) for k in ORDERS}
        return dict(rows=rows,roots=dict(a=a,kappa_minus2=Delta),record=dict(status='enclosed',
            branch='flat',normalized_original_excess_cover=d,active_body_enclosure_evaluated=False,
            all_original_q_slow_jets_exact_zero=True,original_q_C0=rows[ZERO].record(),
            original_q_ordinary_y2_Z1_rows={'y%d_Z%d'%k:x.record() for k,x in rows.items()}))
    if dhi>=2:
        return dict(rows=None,roots=dict(a=a,kappa_minus2=Delta),record=dict(status='requires_positive_body_source_subdivision',
            branch='unresolved',normalized_original_excess_cover=d,active_body_enclosure_evaluated=False))
    argument={k:scalar(1-d) if k==ZERO else -row for k,row in normalized.items()}
    coefficients=prior.sigma_jets(c,1-d)
    cutoff=slow.compose_sigma(argument,[coefficients[n]*math.factorial(n) for n in range(4)])
    # Collect eta*(2-Delta/eta) before independent giant-log subtraction.
    gamma={k:eta*(2-d) if k==ZERO else -row for k,row in Delta.items()}
    ratio=current.quotient_jet(gamma,{k:row*2 for k,row in a.items()},log_a_lower+c.ln(2))
    log_a_upper=ep(a[ZERO].record()['log_absolute_upper'])[1]
    lower=(eta_log+c.ln(c.mpf(2)-c.mpf(dhi))-c.mpf(log_a_upper)-c.ln(2))/2
    root=slow.sqrt_jet(ratio,lower)
    rows=current.multiply_jet(cutoff,root)
    if any(not row.zero for k,row in rows.items() if k[1]):raise ArithmeticError('Original O3 q_Z must be exact zero')
    record=dict(status='enclosed',branch='active' if dhi<1 else 'smooth_cutoff_seam',
        normalized_original_excess_cover=d,original_sigma_argument_cover=1-d,
        original_positive_body_argument_cover=2-d,
        positive_body_theorem='Delta/eta<2 implies 2eta-Delta=eta*(2-Delta/eta)>0; original flat cutoff contributes exactly zero when Delta>=eta',
        original_active_gamma_cover=gamma[ZERO].record(),positive_original_root_lower_log=lower,
        original_q_C0=rows[ZERO].record(),active_body_enclosure_evaluated=True,
        original_q_ordinary_y2_Z1_rows={'y%d_Z%d'%k:x.record() for k,x in rows.items()},
        same_original_flat_cutoff_evaluated_through_seam=True,
        all_original_q_Z_rows_exact_zero=True,ordinary_derivatives_not_Taylor_or_normalized_coordinate_derivatives=True)
    return dict(rows=rows,roots=dict(a=a,kappa_minus2=Delta),record=record)


def typed_original_phase(binder,coordinates,geometry,Z,N,base=None):
    """Original affine O3 radius, keeping its microscopic offset separate."""
    c=coordinates.ctx
    if base is None:base=binder.query('O3_slope_mu',Z,0,N)
    if geometry['kind']=='xi':
        lo,hi=geometry['left'],geometry['right']
        coefficient=c.mpf((ep(c.mpf(lo.numerator)/lo.denominator)[0],ep(c.mpf(hi.numerator)/hi.denominator)[1]))
        micro=common.restore_common_source(geometry['record']['original_parameter_source']['positive_original_W'],coordinates)*coefficient
    else:
        D=geometry['record']['original_parameter_source']['entire_original_D_cover']
        k=c.mpf((ep(c.mpf(geometry['right'].numerator)/geometry['right'].denominator)[0],
                 ep(c.mpf(geometry['left'].numerator)/geometry['left'].denominator)[1]))
        micro=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,offset=-c.ln(D+k)/2),1,coordinates.ledger)
    zero_base=base['offset']
    origin=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,offset=zero_base.scale.evaluate()),zero_base.coefficient,coordinates.ledger)
    projection=[]
    if base['record']['periodic_projection']['full_period']:projection=[c.mpf((0,1))]
    else:
        for box in base['phase_boxes']:
            projection.extend(first.spatial.ordinary_mod_one(c,box+geometry['coordinate']*N)['boxes'])
    return dict(original_base_source_phase=base['record'],typed_source_geometry=geometry['record'],
        original_microscopic_s_component=micro.record(),original_microscopic_N_s_component=(micro*N).record(),
        actual_original_log_R_over_r_minus=(origin+micro).record(),actual_phase_boxes=projection,
        same_original_affine_radius_identity=binder.identity,
        actual_phase_is_N_times_original_s_not_N_times_xi_or_k=True,
        nonzero_microscopic_phase_component_retained=True,phase_independent_of_Z=True,
        phase_box_not_an_independently_selected_sample=True)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/previous.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,previous,accepted['source_family'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(previous.generic.tail.baseline.build_owner(bridge)))
        native=previous.PositiveLogTransition(owner)
        coord=OriginalTransitionCoordinates(owner.coordinates,native.mu,native.eta)
        base=owner.coordinates.scalar(1);rows=[]
        binder=owner.transfer.geometry.binder
        phase0=binder.query('O3_slope_mu',('.36','.38'),0,N)
        for label,kind,left,right in (
            ('open_bulk','xi','1/4','3/4'),('near_seam_active','k',4,2),
            ('seam_active_collar','k',2,'7/4'),('cutoff_seam','k','7/4','5/3'),
            ('seam_flat_collar','k','5/3',0)):
            geometry=coord.geometry(kind,left,right)
            normalized=coord.normalized_excess(base,geometry);assert normalized['rows'] is not None
            q=correlated_original_q(normalized['rows'],native.eta,native.positive['log_actual_a_positive_lower'])
            assert q['rows'] is not None
            masses={k:rc.transfer.true_width_kernel(owner.coordinates,geometry,rate) for k,rate in rc.RATES.items()}
            phase=typed_original_phase(binder,owner.coordinates,geometry,('.36','.38'),N,phase0)
            rows.append(dict(label=label,actual_original_geometry=geometry['record'],
                original_correlated_normalized_excess=normalized['record'],original_q_source=q['record'],
                original_true_width_mass_and_decay={k:dict(branch=x['branch'],mass=x['mass'].record(),decay=x['decay'].record()) for k,x in masses.items()},
                actual_original_typed_radius_phase=phase,
                original_amplitude_inverse_density_and_integral_not_evaluated=True))
            print('Original microscopic source:',label,q['record']['branch'],flush=True)
        endpoint=correlated_original_q({k:base.scalar(0) for k in ORDERS},native.eta,native.positive['log_actual_a_positive_lower'])
        payload=dict(source_family=owner.family,candidate_N=N,original_parameter_source=coord.parameter_record,
            common_directed_coordinate_theorem=owner.coordinates.record(),original_source_coordinate_queries=rows,
            exact_original_s_zero_active_q_endpoint=endpoint['record'],
            original_positive_a_lower_log=native.positive['log_actual_a_positive_lower'],
            source_q_Z_independent_on_original_axial_domain=True,
            original_bulk_prefix_before_xi_quarter_and_bulk_to_seam_not_integrated=True,
            full_original_prefix_amplitude_inverse_density_integral_admitted=False)
        raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
        filename=PREFIX+'current_transition_normalized_q_source_evidence.json.gz'
        (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Same original source closures required')
            hashes[name]=digest
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,original_parameter_source=coord.parameter_record,
            original_typed_microscopic_source_archive=dict(filename=filename,compressed_bytes=len(compressed),
                uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
            query_branches={r['label']:r['original_q_source']['branch'] for r in rows},
            exact_original_zero_endpoint_active_q_retained=True,
            original_smooth_cutoff_seam_y2_Z1_source_enclosed=True,
            genuine_original_nonzero_q_bulk_and_seam_source_enclosures_executed=True,
            actual_micro_width_and_original_N_s_phase_retained=True,
            full_prefix_or_tail_numerical_integrals_or_Rc_targets_admitted=False,
            native_amplitude_source_evaluator_requires_safe_original_transition_kernels=True,
            **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,source_seed_evidence=seed,
            execution_seconds=time.monotonic()-began,
            scope='Same original mu/eta/q. Typed xi bulk and k cutoff-seam coordinates retain whole parameter covers, positive true width and actual N*s phase; original q C0/y2/Z1 at active, smooth-seam, flat and exact-zero branches. Original E/p2 amplitude/inverse/density/full prefix integral and all later gates remain open.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
