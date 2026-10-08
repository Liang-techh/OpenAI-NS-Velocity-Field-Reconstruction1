"""Actual original O2-slope joint source and inherited C1 obstruction.

Read saved 2048-cell source/inverse boxes, without rerunning source owners or
inverse solves. Restore the original native factor basis before evaluating
primitives, then collect k-a_eff*m before phase/source unions.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_Rc_joint_terminal_defects as terminal
import lei_ren_part1_paper_compliant_current_inner_reference_weighted_pressure as weighted

common=terminal.common;native=terminal.native;prior=terminal.prior;rc=terminal.rc
HERE,PREFIX,sha,encode,ep,iv=terminal.HERE,terminal.PREFIX,terminal.sha,terminal.encode,terminal.ep,terminal.iv
source=weighted.common;original=source.previous;N=terminal.N;KEYS=terminal.KEYS
NAME=PREFIX+'current_slope_joint_source.json';RECEIPT=PREFIX+'current_slope_joint_source_check.json'
GATE='original_2048_cell_slope_joint_source_and_pre_slope_C1_attribution_executed'


def native_factored_expm1(value):
    """Same integral identity; retain microscopic native log-|u| factors."""
    if value.zero:return value
    c=value.ctx;lo,hi=ep(original.bounded(value))
    if max(abs(lo),abs(hi))>1:raise ArithmeticError('Original source exponent must remain in [-1,1]')
    return value*c.mpf((ep(c.exp(c.mpf(min(lo,0))))[0],ep(c.exp(c.mpf(max(hi,0))))[1]))


def native_restore(record,bases,ledger):
    scale=record['formal_positive_scale'];c=bases[0].ctx
    value=prior.ScaledEnclosure(prior.FormalScale(bases,tuple(scale['source_exponents'])+(scale['radius_power'],),iv(c,scale['additional_log_interval'])),iv(c,record['coefficient_interval']),ledger)
    if value.zero!=record['exact_zero'] or record['point_value_selected']:raise ValueError('Saved original cover required')
    return value


def native_to_half(value,coordinates,logC,y):
    """Collect original logPstar, delta and R powers before conversion."""
    c=coordinates.ctx;p,d,ell,u,r=value.scale.powers
    logP=coordinates.logP_squared/2
    if value.scale.bases[0]._mpi_!=logP._mpi_:raise ValueError('Same original Pstar required')
    offset=value.scale.offset-30*d+ell*value.scale.bases[2]+u*value.scale.bases[3]+r*(c.ln(110)+10*logC+y)
    return prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,(0,0,0,2*(p-4*d+10*r),0),offset),value.coefficient,coordinates.ledger)


def original_kernel(coordinates,cell,phase,logC,dstar):
    c=coordinates.ctx;src=cell['source'];row=src['original_source_function_coefficient_record']
    assert row['original_basis_contract']['order']==['logPstar','logdelta','logL','zero','logR']
    assert row['original_basis_contract']['delta_and_R_Z_independent']
    z=c.mpf(src['exact_Z_range']);left,right=map(Fraction,src['exact_y_cell'])
    y=c.mpf((ep(c.mpf(left.numerator)/left.denominator)[0],ep(c.mpf(right.numerator)/right.denominator)[1]))
    logP=coordinates.logP_squared/2;logL=c.ln(iv(c,row['original_basis_contract']['positive_L_interval']))
    geom=src['conditioned_source_piece_geometry'][phase['source_piece_index']]
    urecord=geom['original_u'];assert urecord['formal_positive_scale']['source_exponents']==[0,0,0,1]
    assert not urecord['exact_zero'] and all(abs(x)==1 for x in ep(iv(c,urecord['coefficient_interval'])))
    logu=iv(c,urecord['log_absolute_upper'])
    bases=(logP,-4*logP-30,logL,logu,c.ln(110)+10*(logC+logP)+y)
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    roots={k:{(0,int(order[-1])):native_restore(v,bases,ledger) for order,v in rows.items()} for k,rows in row['native_source_root_enclosures'].items()}
    roots['p2'][(0,1)]=native_restore(phase['correlated_original_p2_Z_carrier']['correlated_same_original_source_root'],bases,ledger)
    q=native_restore(src['native_positive_q'],bases,ledger);u=native_restore(urecord,bases,ledger)
    kernel=original.PositiveLogQPhase(dict(roots=roots,q=q,original_u_source=u),dstar)
    assert kernel.geometry==geom['branch']
    return kernel,roots,y,z,bases


def original_joint_density(kernel,roots,y,z,phase,mu_interval):
    """a_eff=qi*exp(3/10-y/2-5mu/2), distance logPstar+3-y."""
    c=kernel.c;scalar=kernel.scalar
    inv=phase['original_inverse'];primitive=kernel.primitives(iv(c,inv['coordinate_interval']),inv['chart'])
    jets={k:native_restore(v,kernel.q.scale.bases,kernel.q.ledger) for k,v in phase['native_original_primitive_Z_enclosures'].items()}
    E,EZ=roots['E'][(0,0)],roots['E'][(0,1)];V,VZ=roots['V'][(0,0)],roots['V'][(0,1)]
    x=primitive['A']*(c.mpf(1)/N);xZ=jets['A_Z_slow']*(c.mpf(1)/N)
    increment=native_factored_expm1(x);exponential=c.exp(original.bounded(x))
    dE=E*increment;dEZ=EZ*increment+E*exponential*xZ
    dV=primitive['B_over_Pstar']*(c.mpf(1)/N);dVZ=jets['B_Z_slow']*(c.mpf(1)/N)
    base=scalar(c.exp(c.mpf('.3')-y/2)/(1+z*z))
    eps=rc.density.density.factored_expm1(scalar(mu_interval)*c.mpf('-2.5'))
    ae=base+base*eps;aeZ=ae*(-2*z/(1+z*z))
    # Keep the tiny mu term separately from the moderate slope mismatch.
    mismatch=E-base-base*eps;mismatchZ=EZ-base*(-2*z/(1+z*z))-base*eps*(-2*z/(1+z*z))
    C=V*dE+dE*dV+mismatch*dV
    CZ=VZ*dE+V*dEZ+dEZ*dV+dE*dVZ+mismatchZ*dV+mismatch*dVZ
    return dict(C=C,C_Z=CZ,record=dict(actual_original_A=primitive['A'].record(),actual_original_B_over_Pstar=primitive['B_over_Pstar'].record(),
        original_effective_amplitude=ae.record(),original_effective_amplitude_Z=aeZ.record(),
        original_nonzero_relative_mu_increment=eps.record(),original_E_minus_effective_amplitude=mismatch.record(),
        actual_joint_density=C.record(),actual_joint_density_Z=CZ.record(),
        actual_primitive_Z_uses_original_derivative_inverse_box=True,
        original_saved_inverse_coordinate_used_no_solver_rerun=True,nonlinear_joint_before_phase_or_source_hull=True))


def source_row(coordinates,cell,index,logC,dstar,mu_interval):
    rows=[];densities=[]
    for phase in cell['actual_source_phase_held_Z_records']:
        kernel,roots,y,z,bases=original_kernel(coordinates,cell,phase,logC,dstar)
        got=original_joint_density(kernel,roots,y,z,phase,mu_interval)
        densities.append({k:native_to_half(got[k],coordinates,logC,y) for k in ('C','C_Z')})
        rows.append(dict(original_native_source_log_bases=bases,actual_joint_source=got['record']))
    union=rc.density.local.same_source_union;c=coordinates.ctx
    den={k:union([row[k] for row in densities]) for k in ('C','C_Z')}
    mass=coordinates.scalar(iv(c,cell['positive_own_rate_final_endpoint_masses']['k']))
    return dict(C=den['C']*mass,C_Z=den['C_Z']*mass,record=dict(original_source_cell_index=index,
        exact_y_cell=cell['source']['exact_y_cell'],actual_original_phase_joint_sources=rows,
        original_rate_three_halves_mass_to_slope_exit=mass.record(),joint_source_union=terminal.native.records(den),
        actual_joint_C_contribution=(den['C']*mass).record(),actual_joint_C_Z_contribution=(den['C_Z']*mass).record()))


def pre_slope_attribution(coordinates,payload,pressure_tile):
    """Identify the inherited obstruction using the original 13-cell route."""
    c=coordinates.ctx;old=weighted.replay.inputs(c,payload)[0]
    lift=lambda rec:coordinates.rebase(source.restore_common_source(rec,old),payload['source_family'])
    flat,rows=weighted.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])
    local=[];final={k:[] for k in KEYS}
    decays={k:coordinates.scalar(1) for k in KEYS}
    for label,row in reversed(rows):
        gkey,skey,out,bkey,own,inherited=weighted.accepted.row_fields(label)
        geo=row[gkey];width=lift(geo['positive_true_log_radius_width'])
        g=dict(width=width,regular=iv(c,geo['regular_true_log_radius_width_cover']),scalar_cover=iv(c,geo['scalar_width_cover_used_only_for_directed_kernel_bounds']))
        records={}
        for k,rate in rc.RATES.items():
            v=lift(row['true_log_radius_signed_C0_contributions'][k]);z=lift(row['true_log_radius_signed_Z_contributions'][k])
            if label=='inner_reference' and k=='p':v=lift(pressure_tile['pressure_integral_proof']['selected_pressure_integral'])
            fv=v*decays[k];fz=z*decays[k];final[k].append((label,fv,fz))
            records[k]=dict(actual_C0_at_O2_inlet=fv.record(),actual_Z_at_O2_inlet=fz.record())
            decays[k]=decays[k]*rc.transfer.true_width_kernel(coordinates,g,rate)['decay']
        local.append(dict(label=label,source_contributions=records))
    dominant={}
    for k in KEYS:
        biggest=max(final[k],key=lambda row:ep(row[2].scale.evaluate()+c.ln(c.mpf(max(abs(x) for x in ep(row[2].coefficient)))))[1] if not row[2].zero else -c.inf)
        dominant[k]=dict(label=biggest[0],actual_inherited_Z_cap=terminal.magnitude(biggest[2]))
    first=payload['twelve_source_primitive_cover_records'][0]
    branch=first['genuine_same_source_conditional_paired_cover']
    return dict(actual_original_pre_slope_source_contributions=list(reversed(local)),dominant_inherited_Z_source=dominant,
        firstbridge_original_roots=branch['original_root_rows'],firstbridge_selected_primitive_C0_Z=first['selected_primitive_C0_Z_covers'],
        firstbridge_conditional_cutoff_branch_derivative_covers=[dict(name=b['name'],paired_derivative_covers=b['paired_derivative_covers']) for b in branch['conditional_cutoff_paired_branches']],
        source_defined_initial_flat_correction_zero=True,weighted_pressure_C0_refinement_retained=True,
        no_upstream_source_contribution_discarded_or_relabelled_as_slope_own=True)


def terminal_composition(coordinates,old,bd,ad,report,primary,C,CZ,amp,trace):
    """Separate pre-slope memory, slope own, axial, buffer and transition."""
    c=coordinates.ctx;zero=coordinates.scalar(0)
    # All primary upstream records are in the common Pstar-squared basis.
    common_bases=(c.mpf(0),coordinates.logP_squared,c.mpf(0),c.mpf(0),c.mpf(0))
    ledger=coordinates.ledger
    lift=lambda rec:coordinates.rebase(native_restore(rec,common_bases,ledger),coordinates.family)
    v={k:lift(primary['actual_refined_upstream_C0'][k]) for k in KEYS}
    z={k:lift(primary['actual_refined_upstream_Z'][k]) for k in KEYS}
    own={k:coordinates.scalar(iv(c,report['five_original_C0_integral_contributions'][k])) for k in KEYS}
    ownZ={k:coordinates.scalar(iv(c,report['five_genuine_ordinary_Z_integral_contributions'][k])) for k in KEYS}
    zbox=c.mpf(ad['exact_Z_range']);base=coordinates.scalar(c.exp(c.mpf('.3'))/(1+zbox*zbox))
    ae=base+base*amp['epsilon'];aeZ=ae*amp['ratio']
    incoming=v['k']-ae*v['m'];incomingZ=z['k']-aeZ*v['m']-ae*z['m']
    # Original axial width exp(Md)-1, then 11+1+2; the slope adds one.
    after=coordinates.logP_squared/2+2;all_width=after+1
    decay=coordinates.decay(after,'3/2');all_decay=coordinates.decay(all_width,'3/2')
    pairs=dict(pre_slope_inherited_component_cover=dict(C=incoming*all_decay,C_Z=incomingZ*all_decay),
        original_slope_joint_source=dict(C=C*decay,C_Z=CZ*decay))
    for name in ('original_axial_joint_source','original_buffer_joint_source','original_transition_joint_source'):
        row=old['source_contributions_to_divided_joint_defect'][name]
        pairs[name]=dict(C=terminal.complete.restore_half_source(row['joint_C'],coordinates),
            C_Z=terminal.complete.restore_half_source(row['joint_C_Z'],coordinates))
    total=sum((row['C'] for row in pairs.values()),zero);totalZ=sum((row['C_Z'] for row in pairs.values()),zero)
    history=bd['actual_updated_Rd_Rc_correction_and_own_history']
    values=terminal.restored(history['actual_updated_Rc_correction_C0'],coordinates)
    jets=terminal.restored(history['actual_updated_Rc_correction_Z'],coordinates)
    targets=terminal.normalize(values,jets,total,totalZ,amp);parts={}
    for name,rows,width in (('pre_slope_inherited_component_cover',(v,z),all_width),('original_slope_joint_source',(own,ownZ),after)):
        factors={k:coordinates.decay(width,rate) for k,rate in rc.RATES.items()}
        pv={k:rows[0][k]*factors[k] for k in KEYS};pz={k:rows[1][k]*factors[k] for k in KEYS}
        part=terminal.normalize(pv,pz,pairs[name]['C'],pairs[name]['C_Z'],amp)
        parts[name]=dict(joint_C=pairs[name]['C'].record(),joint_C_Z=pairs[name]['C_Z'].record(),five_target_diagnostics=terminal.summary(part))
    for name in ('original_axial_joint_source','original_buffer_joint_source','original_transition_joint_source'):
        parts[name]=old['source_contributions_to_divided_joint_defect'][name]
    chart_rows=[]
    for row in trace['actual_original_pre_slope_source_contributions']:
        chartv={k:terminal.complete.restore_half_source(encode(row['source_contributions'][k]['actual_C0_at_O2_inlet']),coordinates) for k in KEYS}
        chartz={k:terminal.complete.restore_half_source(encode(row['source_contributions'][k]['actual_Z_at_O2_inlet']),coordinates) for k in KEYS}
        chartC=chartv['k']-ae*chartv['m'];chartCZ=chartz['k']-aeZ*chartv['m']-ae*chartz['m']
        chart_rows.append(dict(label=row['label'],actual_joint_C_at_Rc=(chartC*all_decay).record(),actual_joint_C_Z_at_Rc=(chartCZ*all_decay).record()))
    return dict(original_Rc_amplitude_and_parameter=amp['record'],exact_post_slope_to_Rc_width='logPstar+2=exp(Md)+13',
        exact_pre_slope_to_Rc_width='logPstar+3=exp(Md)+14',original_post_slope_rate_three_halves_decay=decay.record(),
        actual_pre_slope_joint_component_at_O2_inlet=incoming.record(),actual_pre_slope_joint_component_Z_at_O2_inlet=incomingZ.record(),
        actual_joint_terminal_C=total.record(),actual_joint_terminal_C_Z=totalZ.record(),
        actual_signed_five_terminal_defect_C0=native.records(targets['values']),actual_signed_five_terminal_defect_Z=native.records(targets['Z_derivatives']),
        actual_N_scaled_targets_C0=native.records({k:value*N for k,value in targets['values'].items()}),
        actual_N_scaled_targets_Z=native.records({k:value*N for k,value in targets['Z_derivatives'].items()}),
        actual_five_target_diagnostics=terminal.summary(targets),five_separate_source_terminal_contributions=parts,
        pre_slope_chart_joint_component_contributions_to_Rc=chart_rows,
        preceding_composite_slope_inlet_diagnostics=old['actual_five_target_diagnostics'],
        original_five_control_linear_response_enclosures=terminal.linear_control_enclosures(coordinates,targets,amp),
        actual_pre_slope_root_correlation_not_recovered=True,actual_slope_own_joint_before_phase_hulls=True,
        unchanged_original_component_histories_used_for_other_four_rows=True,pressure_rate_zero_memory_and_separate_original_P0_retained=True,
        actual_nonlinear_controls_or_terminal_function_identity_admitted=False)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();hashes={};accepted=json.loads((HERE/terminal.NAME).read_bytes())
    common.attach_receipt(hashes,terminal,accepted['source_family']);wm,_=common.attach_receipt(hashes,weighted,accepted['source_family'])
    cm,_=common.attach_receipt(hashes,source,accepted['source_family']);c=MPIntervalContext();c.dps=240;archives=[]
    for oldarchive in accepted['actual_original_Rc_joint_terminal_defect_archives']:
        old=terminal.load_archive(oldarchive);bd=terminal.load_archive(old['accepted_original_buffer_archive']);ad=terminal.load_archive(bd['accepted_original_axial_archive'])
        coordinates=native.HalfPstarCoordinates(c,iv(c,old['common_directed_coordinate_theorem']['common_log_bases'][1]),accepted['source_family'])
        amp=terminal.amplitude(coordinates,ad);logC=iv(c,ad['original_logC']);dstar=iv(c,ad['original_dstar_log'])
        report=next(row for row in cm['original_complete_same_N_source_incoming_O2_integral_refinements'] if row['ordered_source_cells']==2048 and weighted.same_Z(row['exact_Z_range'],ad['exact_Z_range']))
        pressure_tile=next(row for row in wm['actual_weighted_pressure_tiles'] if weighted.same_Z(row['exact_Z_range'],ad['exact_Z_range']))
        primary=next(row for row in wm['genuine_original_O2_weighted_pressure_replays'] if row['ordered_source_cells']==2048 and weighted.same_Z(row['exact_Z_range'],ad['exact_Z_range']))
        payload=terminal.load_archive(pressure_tile['original_source_archive']);trace=pre_slope_attribution(coordinates,payload,pressure_tile)
        rows=[];C=coordinates.scalar(0);CZ=coordinates.scalar(0);query_count=0
        for archive in report['whole_original_density_and_density_Z_source_archives']:
            data=terminal.load_archive(archive);left,right=archive['source_cell_index_range']
            assert data['source_family']==accepted['source_family'] and data['candidate_N']==N and data['source_cell_count']==2048
            assert weighted.same_Z(data['exact_Z_range'],ad['exact_Z_range'])
            for index,cell in enumerate(data['records'],start=left):
                got=source_row(coordinates,cell,index,logC,dstar,amp['mu_interval']);C=C+got['C'];CZ=CZ+got['C_Z'];rows.append(got['record'])
                query_count+=len(cell['actual_source_phase_held_Z_records'])
            assert len(rows)==right
            print('Actual original slope joint source:',ad['exact_Z_range'],right,'/ 2048',flush=True)
        composition=terminal_composition(coordinates,old,bd,ad,report,primary,C,CZ,amp,trace)
        raw=dict(source_family=accepted['source_family'],candidate_N=N,exact_Z_range=ad['exact_Z_range'],
            common_directed_coordinate_theorem=coordinates.record(),original_effective_amplitude_identity='a_eff(y)=qi*exp(3/10-y/2-5mu/2)',
            exact_original_source_to_Rc_distance='logPstar+3-y=exp(Md)+14-y',
            original_complete_2048_cell_joint_source_rows=rows,actual_joint_source_phase_queries=query_count,
            slope_own_joint_C_at_slope_exit=C.record(),slope_own_joint_C_Z_at_slope_exit=CZ.record(),
            actual_updated_Rc_terminal_defects_and_source_ledger=composition,
            actual_pre_slope_C1_contribution_attribution=trace,accepted_original_slope_source_archives=report['whole_original_density_and_density_Z_source_archives'],
            accepted_original_terminal_archive=oldarchive,accepted_original_pre_slope_source_archive=pressure_tile['original_source_archive'],
            slope_own_is_separate_from_actual_upstream_incoming=True,actual_new_inlet_terminal_control_or_functional_closure_admitted=False,
            **dict.fromkeys(common.current.FLAGS,False))
        encoded=json.dumps(encode(raw),indent=2).encode()+b'\n';compressed=gzip.compress(encoded,compresslevel=9,mtime=0)
        tag='positive' if Fraction(ad['exact_Z_range'][0])>0 else 'negative';name=PREFIX+'current_slope_joint_source_'+tag+'.json.gz'
        (HERE/name).write_bytes(compressed);hashes[name]=sha(name)
        archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(encoded),lossless_original_json_sha256=hashlib.sha256(encoded).hexdigest(),exact_Z_range=ad['exact_Z_range']))
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=accepted['source_family'],candidate_N=N,actual_original_slope_joint_source_archives=archives,
        entire_original_slope_y_zero_to_one_joint_source_on_two_tiles_executed=True,
        first_inherited_C1_obstruction_attributed_to_original_pre_slope_charts=True,
        actual_inverse_solvers_and_upstream_producers_not_rerun=True,
        actual_original_new_terminal_function_oracle_controls_global_N_heat_stress_recursion_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
