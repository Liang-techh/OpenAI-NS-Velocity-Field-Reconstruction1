"""Actual original O2 midplane five own-rate contributions and C0 transport.

Accepted entire-cell source ranges are reused as ranges, never field points.
Both primitives feed the original graph. Positive exact Duhamel cell masses
transport signed densities, original nonzero histories and explicit incoming
defects. Z-functional/global matching and C1 stress recovery remain open.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_pressure_own_integral as pressure

density=pressure.density;base=pressure.base;HERE,PREFIX,sha=pressure.HERE,pressure.PREFIX,pressure.sha
ep=pressure.ep;RATES=density.recovery.RATES;UNITS=density.recovery.UNITS
NAME=PREFIX+'current_original_O2_five_own_integrals.json.gz'
RECEIPT=PREFIX+'current_original_O2_five_own_integrals_check.json'
GATE='actual_original_O2_midplane_all_five_own_rate_contributions_and_C0_affine_transport_enclosed'


def interval(c,record):return base.conditioned.packets.interval(c,record)


def restored_range(record,bases,ledger):
    """Restore the whole source range, including its formal positive scale."""
    c=bases[0].ctx;scale=record['formal_positive_scale']
    powers=tuple(scale['source_exponents'])+(scale['radius_power'],)
    result=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,powers,
        interval(c,scale['additional_log_interval'])),interval(c,record['coefficient_interval']),ledger)
    if result.zero!=record['exact_zero']:raise ValueError('Restored whole-cell range zero contract differs')
    return result


def original_inlet(c):
    return dict(m=c.mpf(0),h=c.mpf(5)/8,k=c.mpf(0),e=-c.mpf(5)/12,p=c.mpf(5)/2)


def masses(c,width,left,right):
    """Exact positive local and final-endpoint Duhamel masses."""
    local={};final={};decay={}
    for key,r in RATES.items():
        rate=c.mpf(r);decay[key]=c.exp(-rate*width)
        local[key]=width if r==0 else width*density.recovery.exp_average(c,-rate*width)
        final[key]=local[key]*c.exp(-rate*(1-right))
        if ep(final[key])[0]<=0:raise ArithmeticError('Positive original Duhamel mass lost')
    return local,final,decay


def apply_incoming(c,record,incoming,*,source_family,original_P0_datum_sha256):
    """C0 midplane values only; no default zero incoming five defects."""
    if source_family!=record['source_family'] or original_P0_datum_sha256!=source_family['datum_enclosure_sha256']:
        raise ValueError('Same original source family and P0 datum required')
    if set(incoming)!=set(RATES):raise ValueError('All five explicit incoming normalized defect enclosures required')
    result={}
    for key in RATES:
        value=c.mpf(incoming[key])
        if any(not mp.isfinite(x) for x in ep(value)):raise ValueError('Finite C0 inlet enclosures required')
        old=interval(c,record['original_five_histories_at_y1'][key])
        change=interval(c,record['five_own_rate_integral_contributions'][key])
        decay=interval(c,record['incoming_defect_decay_coefficients'][key])
        result[key]=old+change+decay*value
    return dict(source_family=source_family,original_Z_exact='0',own_five_histories_at_y1=result,
        incoming_defects_preserved=True,original_P0_datum_sha256=original_P0_datum_sha256,
        P0_not_added_to_or_reset_as_cumulative_pressure=True,C1_or_Z_functional_transport_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False)


class OriginalO2FiveOwnIntegrals:
    def __init__(self):
        admitted=json.loads((HERE/pressure.RECEIPT).read_bytes())
        if not admitted.get('all_passed') or not admitted.get(pressure.GATE):raise ValueError('Accepted genuine original pressure integral required')
        self.parent=pressure.OriginalO2PressureOwnIntegral();self.owner=self.parent.owner;self.family=self.parent.family
        if admitted['source_family']!=self.family:raise ValueError('Original five transport family differs')
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**admitted['input_hashes'],pressure.RECEIPT:sha(pressure.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original five-integral dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Pressure/five-integral source dependencies disagree')
            self.hashes[name]=digest
        compact=admitted['compressed_producer_report']['filename']
        content=gzip.decompress((HERE/compact).read_bytes())
        import hashlib
        if hashlib.sha256(content).hexdigest()!=admitted['compressed_producer_report']['lossless_original_json_sha256']:
            raise ValueError('Lossless source cell archive differs from accepted receipt')
        self.saved=json.loads(content);self.compact=compact
        if self.saved['source_family']!=self.family or not self.saved[pressure.GATE]:raise ValueError('Accepted entire-source-cell archive required')
        self.graph=self.parent.density.graph
        E=sy.Symbol('E',positive=True)
        assert density.recovery.history_densities(E,sy.Integer(0))==dict(m=0,h=E,k=0,e=-E*E/2,p=E*E/2)
        # These exact inlet constants come from the same defining O2 source.
        profile=self.owner.inputs.frame.owner.profiles;row=profile.radial(0);p=profile.ctx
        assert row['f']==1 and row['H']==p.mpf(5)/8 and row['D']==p.mpf(5)/12 and row['P']==p.mpf(5)/2
        self.inlet_binding=dict(passed=True,original_midplane_normalized_reference_inlet=['0','5/8','0','-5/12','5/2'],
            exact_zero_defining_J_and_masses_at_y0=True,original_nonzero_h_e_p_preserved=True,
            same_original_P0_separate_from_cumulative_p=True,source_midplane_history_density_identity=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def integrate_cached_level(self,refinement,bits=24):
        if type(bits) is not int or not 8<=bits<=256:raise ValueError('Explicit original inverse bits in[8,256] required')
        c=self.owner.ctx;count=refinement['ordered_source_cells'];N=refinement['explicit_candidate_N']
        if refinement['source_family']!=self.family or refinement['original_Z_exact']!='0' or refinement['exact_y_window']!=['0','1']:
            raise ValueError('Same original full-window midplane source-cell cache required')
        began=time.monotonic()
        with mp.workdps(c.dps+40):
            width=c.mpf(1)/count;logP=c.exp(40)+11;logdelta=-4*logP-30
            logC=c.mpf(mp.mp.make_mpf(self.owner.inputs.frame.selected_logCstar_mpf_tuple))
            logRref=c.ln(110)+10*(logC+logP)
            change={key:c.mpf(0) for key in RATES};old_integral={key:c.mpf(0) for key in RATES}
            prefix_change={key:c.mpf(0) for key in RATES};prefix_old=original_inlet(c)
            rows=[];prefixes=[];phase_pieces=0
            for i,saved in enumerate(refinement['whole_source_cells']):
                left,right=c.mpf(i)/count,c.mpf(i+1)/count
                if saved['exact_y_cell']!=[str(i)+'/'+str(count),str(i+1)+'/'+str(count)]:raise ValueError('Ordered source cache has a gap')
                if not saved['entire_cell_source_and_inverse_enclosure_not_point_quadrature'] or saved['source_caps_or_midpoints_selected_as_field_values']:
                    raise ValueError('Actual entire-cell source/error ranges required')
                ys=c.mpf((ep(left)[0],ep(right)[1]));bases=(logP,logdelta,c.mpf(0),c.mpf(0),logRref+ys)
                ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                    positive_function_root_intersections=0,directed_independent_log_rescalings=0)
                scalar=lambda value:base.prior.ScaledEnclosure(base.prior.FormalScale(bases),value,ledger)
                q=restored_range(saved['q_source_cell'],bases,ledger);eta=restored_range(saved['positive_eta_source'],bases,ledger)
                if q.zero or eta.zero:raise ValueError('A source range touching zero is not an exact flat branch')
                if ep(eta.scale.offset)!=ep(self.owner.scales.logs['eta']):raise ValueError('Selected positive original eta differs')
                E=scalar(interval(c,saved['f_source_cell']));a=scalar(interval(c,saved['a_source_cell']))
                roots={key:{(0,0):value} for key,value in dict(E=E,a=a,V=scalar(0),p2=scalar(0),t0=scalar(0)).items()}
                kernel=base.conditioned.ConditionedPhase(dict(q=q,roots=roots),self.owner.scales.logs['d_star'])
                pieces={key:[] for key in RATES}
                for phi_record in saved['true_common_N_phase_boxes']:
                    phi=interval(c,phi_record);inverse=kernel.evaluate(phi,bits=bits)
                    if inverse['status']!='enclosed':raise ArithmeticError('Original whole-cell inverse needs refinement')
                    psi=2*c.pi*inverse['selected_inverse']['coordinate_interval']
                    A=(a*base.current.square(q)).positive_divide(kernel.nu,0)*(c.sin(2*psi)/(4*c.pi))
                    B=-a*E*q*(c.sin(psi)/(2*c.pi))
                    interpreter=density.BoundDensityGraph(self.graph,dict(kernel=kernel,roots=roots,ledger=ledger),
                        dict(A=A,B_over_Pstar=B),None,N)
                    for key,value in interpreter.values()['densities'].items():pieces[key].append(value.finite_interval())
                    phase_pieces+=1
                hulls={key:c.mpf((min(ep(value)[0] for value in values),max(ep(value)[1] for value in values)))
                    for key,values in pieces.items()}
                local,final,decay=masses(c,width,left,right)
                # Original V=0 is exact here. All changed cross terms still
                # execute in the graph, including nonzero B/delta_V.
                old_rhs=dict(m=c.mpf(0),h=E.finite_interval(),k=c.mpf(0),
                    e=(-base.current.square(E)*c.mpf('.5')).finite_interval(),
                    p=(base.current.square(E)*c.mpf('.5')).finite_interval())
                weighted={}
                for key in RATES:
                    weighted[key]=hulls[key]*final[key];change[key]+=weighted[key]
                    old_integral[key]+=old_rhs[key]*final[key]
                    prefix_change[key]=prefix_change[key]*decay[key]+hulls[key]*local[key]
                    prefix_old[key]=prefix_old[key]*decay[key]+old_rhs[key]*local[key]
                rows.append(dict(original_source_cache_cell_index=i,exact_y_cell=saved['exact_y_cell'],
                    true_common_N_phase_piece_count=len(saved['true_common_N_phase_boxes']),
                    five_signed_density_whole_cell_hulls=hulls,positive_original_final_endpoint_masses=final,
                    five_signed_final_endpoint_contributions=weighted,
                    B_and_delta_V_not_zeroed_with_original_V=True,
                    source_ranges_not_selected_as_points=True))
                if (i+1)%(count//4)==0:
                    prefixes.append(dict(exact_y=str(i+1)+'/'+str(count),
                        original_own_histories=prefix_old.copy(),window_contributions=prefix_change.copy(),
                        physical_incoming_defects_not_set_to_zero=True))
                if (i+1)%max(64,count//8)==0:print('Actual O2 five own integrals:',count,i+1,flush=True)
            inlet=original_inlet(c);total_decay={key:c.exp(-c.mpf(r)) for key,r in RATES.items()}
            old_final={key:inlet[key]*total_decay[key]+old_integral[key] for key in RATES}
            zero_window_end={key:old_final[key]+change[key] for key in RATES}
        return dict(source_family=self.family,original_Z_exact='0',explicit_candidate_N=N,exact_y_window=['0','1'],
            accepted_whole_source_cache=dict(filename=self.compact,source_cell_level=count),
            ordered_source_cells=count,inverse_bits=bits,own_rates=RATES,normalized_own_units=UNITS,
            five_own_rate_integral_contributions=change,original_five_histories_at_y1=old_final,
            candidate_end_before_adding_incoming_defects=zero_window_end,
            incoming_defect_decay_coefficients=total_decay,
            exact_affine_own_history_rule='own_out_j=original_out_j+window_contribution_j+exp(-rate_j)*actual_same_family_incoming_defect_j',
            original_reference_inlet_binding=self.inlet_binding,
            original_P0_datum_sha256=self.family['datum_enclosure_sha256'],original_P0_not_reset=True,
            physical_Pstar_and_R_normalization_remain_formal=True,
            exact_positive_Duhamel_mass_not_an_extra_R_Jacobian=True,
            cached_ranges_are_entire_original_source_function_covers_not_field_points=True,
            source_phase_piece_count=phase_pieces,whole_source_cell_density_and_mass_records=rows,
            actual_cumulative_window_contribution_prefixes=prefixes,
            C0_midplane_five_transport_only=True,C1_or_Z_functional_transport_installed=False,
            numerical_original_source_point_or_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False,
            execution_seconds=time.monotonic()-began)


def run():
    began=time.monotonic();owner=OriginalO2FiveOwnIntegrals();records=[]
    for saved in owner.saved['actual_original_pressure_integral_refinements']:
        records.append(owner.integrate_cached_level(saved))
    result=dict(**{GATE:True},source_family=owner.family,
        mode='actual_original_O2_midplane_five_own_rate_C0_integrals_and_affine_history_transport',
        actual_original_five_own_integral_refinements=records,
        all_five_actual_original_midplane_own_rate_contributions_enclosed=True,
        original_nonzero_inlet_histories_and_explicit_incoming_defect_memory_preserved=True,
        no_original_defining_quadratures_or_pressure_cell_producer_reexecuted=True,
        C1_or_Z_functional_transport_installed=False,
        actual_changed_five_moment_integral_evaluated=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Five genuine original O2 C0 midplane contributions at candidate N7, complete source window and prescribed own rates, original inlet histories and explicit affine incoming-defect memory. Not Z-functional five closure, all-chart own transport, global N/controls, stress, recursion or full NS admission.')
    content=json.dumps(base.encoded(result),indent=2).encode('utf8')+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(content,compresslevel=9,mtime=0))
    return result


if __name__=='__main__':run()
