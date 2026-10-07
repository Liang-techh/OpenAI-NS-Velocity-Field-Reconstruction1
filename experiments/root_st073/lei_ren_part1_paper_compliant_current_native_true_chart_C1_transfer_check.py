"""Independent radius/integral references and original true-chart C1 replay."""
from fractions import Fraction
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_true_chart_C1_transfer as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
require=checks.require;contains=checks.contains;history=current.history;density=current.density


def independent_geometry_kernel_checks(c):
    """Finite scalar fixtures exercise physical endpoint formulas independently."""
    p=mp.mp.clone();p.dps=c.dps+30;family={key:'independent-true-width-'+key for key in packets.FAMILY_KEYS}
    coordinates=history.CommonSourceCoordinates(c,40,family)
    # No native constructor or source graph: these are expressly scalar fixtures.
    binder=object.__new__(current.spatial.NativeSpatialPhase)
    binder.ctx=c;binder.family=family;binder.loghB=c.mpf(-2);binder.loghS=c.mpf(-4);binder.sc=c.mpf('.1')
    binder.bases=tuple(c.mpf(0) for _ in range(5));binder.ledger=coordinates.ledger
    binder.fixed=dict(logP=c.mpf(20),logC=c.mpf(50),T=c.mpf(10),Tw=c.mpf(4))
    binder.seed=SimpleNamespace(params=SimpleNamespace(Md=c.ln(9)))
    binder.identity=current.spatial.affine_identity_theorem()
    geometry=current.TrueChartGeometry(binder,coordinates)
    P,C,T,W,sc=p.mpf(20),p.mpf(50),p.mpf(10),p.mpf(4),p.mpf('.1');M=p.log(9)
    B,S=p.exp(-2),p.exp(-4);Ra=p.log(4)-4*P-1000;ref=p.log(110)+10*(C+P)
    # Original ABSOLUTE radius assignments, not the implementation's collected
    # width formulas. At these bounded fixtures endpoint subtraction is safe.
    physical={
        'bridge_first':lambda s:Ra+B*s,'bridge_second':lambda s:Ra+B*s,
        'bridge_macro':lambda f:Ra+(p.log(100)-Ra)*f+2*B*(1-f),
        'switch_first':lambda s:p.log(100)+S*s,'switch_second':lambda s:p.log(100)+S*s,
        'switch_power':lambda f:p.log(100)+(p.log(110)-p.log(100))*f+2*S*(1-f),
        'reshape':lambda f:p.log(110)+T*f,
        'inner_reference':lambda f:p.log(110)+T+(10*(C+P)-T-8)*f,
        'axial_restore':lambda s:ref-8+s,'restore_buffer':lambda s:ref+s,
        'actual_patch':lambda x:ref-6+p.log(x),'Rh_reference':lambda s:ref+s,
        'O2_slope':lambda s:ref+s,'O2_axial':lambda f:ref+p.exp(M*f),
        'O2_buffer':lambda s:ref+p.exp(M)+s,'O3_slope_mu':lambda s:ref+P+s,
        'O3_power':lambda f:ref+P+1+W*f}
    endpoints=dict(bridge_first=('.1','.3'),bridge_second=('1.1','1.3'),bridge_macro=('.2','.3'),
        switch_first=('.1','.3'),switch_second=('1.1','1.3'),switch_power=('.2','.3'),
        reshape=('.2','.3'),inner_reference=('.2','.3'),axial_restore=('.2','.3'),
        restore_buffer=('-6.8','-6.7'),actual_patch=('1.1','1.2'),Rh_reference=('-3','-2'),
        O2_slope=('.2','.3'),O2_axial=('.2','.3'),O2_buffer=('2','3'),O3_slope_mu=('.2','.3'),O3_power=('.2','.3'))
    radius_rows=0;kernel_rows=0
    cv=lambda value:c.mpf(p.nstr(value,c.dps+20))
    for chart,(a,b) in endpoints.items():
        got=geometry.cell(chart,a,b);w=physical[chart](p.mpf(b))-physical[chart](p.mpf(a))
        require(w>0 and contains(got['width'].finite_interval(),cv(w)),'Independent original physical width failed: '+chart)
        radius_rows+=1
        for rate in (Fraction(0),Fraction(1),Fraction(3,2)):
            rr=p.mpf(rate.numerator)/rate.denominator;factors=current.true_width_kernel(coordinates,got,rate)
            decay=p.exp(-rr*w);mass=w if not rr else -p.expm1(-rr*w)/rr
            require(contains(factors['decay'].finite_interval(),cv(decay)),'Independent true-width attenuation failed: '+chart)
            require(contains(factors['mass'].finite_interval(),cv(mass)),'Independent true-width Duhamel mass failed: '+chart)
            kernel_rows+=2
    collar=geometry.cell('bridge_first',{'selected_sc_multiple':'1/2'},{'selected_sc_multiple':'3/4'})
    require(contains(collar['width'].finite_interval(),cv(B*sc/4)),'Original fractional-sc collar width failed')
    power=geometry.cell('O3_power',{'original_power_offset':'.53'},{'original_power_offset':'.54'})
    require(contains(power['width'].finite_interval(),c.mpf('.01')),'Original power offsets must give exact difference width')
    # Retain unmaterializable microscopic widths, and do not declare positive
    # rate attenuation exactly one. Huge regular quiet widths retain decay.
    binder.loghB=c.mpf('-1e40');tiny=geometry.cell('bridge_first','.1','.3')
    require(not tiny['width'].zero and ep(tiny['width'].coefficient)[0]>0,'Original microscopic width became zero')
    for rate in (Fraction(1),Fraction(3,2)):
        factors=current.true_width_kernel(coordinates,tiny,rate)
        require(not factors['mass'].zero and ep(factors['mass'].coefficient)[0]>0,'Microscopic positive mass lost')
        require(ep(factors['decay'].finite_interval())[0]<1 and ep(factors['decay'].finite_interval())[1]==1,'Positive-rate tiny attenuation declared exact1')
    require(current.true_width_kernel(coordinates,tiny,0)['mass'].scale.offset._mpi_==tiny['width'].scale.offset._mpi_,'Rate0 lost formal width factor')
    huge=geometry.build('independent_huge_quiet',c.mpf('1e40'),{})
    for rate in (Fraction(1),Fraction(3,2)):
        factors=current.true_width_kernel(coordinates,huge,rate);mass=factors['mass'].finite_interval();limit=c.mpf(rate.denominator)/rate.numerator
        require(not factors['decay'].zero and ep(factors['decay'].coefficient)[0]>0,'Huge positive decay became zero')
        require(ep(mass)[0]<ep(limit)[0] and contains(mass,limit),'Huge mass directed exponential tail lost')
    # Nonzero C0/Z incoming memory across nonlinear-coordinate and quiet cells.
    z=c.mpf(('.49','.51'));operator=history.C1DuhamelOperator(coordinates);analytic={key:p.mpf(0) for key in current.RATES};total=p.mpf(0)
    for chart,a,b,quiet in (('actual_patch','1.1','1.2',False),('O2_axial','.2','.3',False),('O2_buffer','2','3',True)):
        got=geometry.cell(chart,a,b);w=physical[chart](p.mpf(b))-physical[chart](p.mpf(a));values={};jets={}
        for index,(key,rate) in enumerate(current.RATES.items()):
            rr=p.mpf(rate.numerator)/rate.denominator;sign=(-1 if index%2 else 1)*(index+1)
            mass=current.true_width_kernel(coordinates,got,rate)['mass']
            values[key]=coordinates.scalar(0) if quiet else mass*(sign*(1+z+z*z))
            jets[key]=coordinates.scalar(0) if quiet else mass*(sign*(1+2*z))
            analytic[key]=p.exp(-rr*w)*analytic[key]+(0 if quiet else sign*(w if not rr else -p.expm1(-rr*w)/rr))
        before=(operator.increments['p'].record(),operator.Z_increments['p'].record())
        current.append_true_cell(operator,got,values,jets,family);total+=w
        if quiet:require(packets.encode(before)==packets.encode((operator.increments['p'].record(),operator.Z_increments['p'].record())),'Quiet pressure C0/Z memory reset')
    incoming={key:coordinates.scalar((index+2)*(1+z)) for index,key in enumerate(current.RATES)}
    incoming_Z={key:coordinates.scalar(index+2) for index,key in enumerate(current.RATES)}
    got=operator.apply(incoming,incoming_Z,family);transfer_rows=0
    for point in ('.49','.5','.51'):
        zz=p.mpf(point)
        for index,(key,rate) in enumerate(current.RATES.items()):
            decay=p.exp(-p.mpf(rate.numerator)/rate.denominator*total)
            references=(decay*(index+2)*(1+zz)+analytic[key]*(1+zz+zz*zz),decay*(index+2)+analytic[key]*(1+2*zz))
            for value,reference in zip((got['values'][key],got['Z_derivatives'][key]),references):
                require(contains(value.finite_interval(),cv(reference)),'Independent nonlinear-coordinate C1 transport failed: '+key);transfer_rows+=1
    return dict(passed=True,independent_original_absolute_radius_difference_comparisons=radius_rows,
        independent_true_width_mass_and_decay_comparisons=kernel_rows,
        independent_signed_nonzero_incoming_nonlinear_coordinate_C0_Z_comparisons=transfer_rows,
        original_fractional_sc_and_power_offset_widths_checked=True,
        unmaterializable_positive_width_mass_and_decay_retained=True,positive_rate_tiny_decay_not_exact_unit=True,
        quiet_pressure_C0_Z_memory_retained=True,scalar_fixtures_do_not_replace_original_native_functions=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['positive_native_chart_length_count']==17,'All17 true chart lengths required')
    require(not any(saved.get(k) for k in packets.OPEN),'True-chart local transport cannot complete global reconstruction')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed true-chart prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        c1=density.NativeDensityC1LocalIntegrals(density.first.NativePhaseFirstJets(density.slow.NativeQSlowJets(density.current.NativeCorrelatedShearQ(history.prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))))
        owner=current.NativeTrueChartC1Transfer(history.NativeC1HistoryTransfer(c1))
        independent=independent_geometry_kernel_checks(owner.ctx)
        lengths=owner.geometry.full_lengths();require(packets.encode(lengths)==saved['original_full_chart_positive_length_records'],'Original true-chart length ledger changed')
        require(all(row['positive_true_log_radius_width']['sign']=='positive' for row in lengths.values()),'Original true lengths must all be positive')
        require(len(owner.owner.binder.identity['original_same_radius_periodic_phase_seam_identities'])==16,'All16 original radius seams required')
        collar_rows=0;contribution_rows=0;regions={}
        for name,old in saved['actual_initial_collar_C1_history_records'].items():
            Z=packets.interval(owner.ctx,old['Z_box']);got=owner.initial_collar(Z,N=saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual initial collar C1 history changed: '+name)
            require(all(row.zero for row in got['correction']['values'].values()) and all(row.zero for row in got['correction']['Z_derivatives'].values()),'Exact flat collar correction lost')
            require(old['original_background_histories_and_P0_not_reset'] and old['remaining_bridge_after_selected_collar_not_assumed_quiet'],'Collar cannot reset background or skip active bridge')
            require(not old['global_inlet_to_Rc_histories_admitted'] and not any(old.get(k) for k in packets.OPEN),'Initial collar cannot admit global gates')
            require(old['original_right_background_and_separate_P0_Z']['P0_not_merged_into_pressure_history'],'Separate original P0/P0_Z required')
            collar_rows+=20;regions[name]=dict(actual_initial_correction_C0_Z_rows=10,actual_own_history_C0_Z_rows=10)
            print('Original actual initial collar C1 histories checked:',name,flush=True)
        for chart,old in saved['actual_true_chart_C1_contribution_records'].items():
            spec=old['actual_true_chart_cell_geometry']['original_endpoint_specification']
            got=owner.contribution(chart=chart,Z=packets.interval(owner.ctx,old['Z_box']),left=spec[0],right=spec[1],N=saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual true-chart signed C1 contribution changed: '+chart)
            require(old['log_radius_Jacobian_integrated_in_true_width_and_not_multiplied_again'] and old['fixed_Z_independent_endpoints_and_phase_allow_differentiation_under_integral'],'True-width/Z integral contract changed')
            for value in [*got['contributions'].values(),*got['Z_derivatives'].values()]:
                require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,'True chart C0/Z common basis changed')
            contribution_rows+=10;regions[chart]=dict(signed_C0_Z_contribution_rows=10,global_incoming_not_assumed=True)
            print('Original true-chart C1 integral checked:',chart,flush=True)
        for left,right in (('.2','.1'),('.1','.1'),('-.1','.2'),(owner.ctx.mpf('.1'),'.2')):
            try:owner.geometry.cell('O2_slope',left,right)
            except ValueError:pass
            else:raise ArithmeticError('Invalid/non-exact true-chart endpoints accepted')
        for mult in ('1/2','1','2'):
            try:owner.initial_collar(right_multiple=mult)
            except ValueError:pass
            else:raise ArithmeticError('Quiet collar inferred beyond certified interval')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        original_positive_chart_lengths_checked=17,original_exact_radius_seams_reused=16,
        actual_initial_collar_C0_Z_history_rows_checked=collar_rows,actual_initial_collar_Z_queries_checked=2,
        actual_true_chart_signed_C0_Z_contribution_rows_checked=contribution_rows,actual_true_chart_C1_cells_checked=5,
        independent_true_width_and_C1_checks=independent,regions=regions,
        original_background_and_P0_Z_retained=True,microscopic_and_long_chart_lengths_not_replaced_by_unit_width=True,
        actual_inlet_to_O2_local_incoming_correction_installed=False,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('True-chart C1 PASS:17 lengths;',collar_rows,'initial history rows;',contribution_rows,'signed integral rows',flush=True)
    return result


if __name__=='__main__':run()
