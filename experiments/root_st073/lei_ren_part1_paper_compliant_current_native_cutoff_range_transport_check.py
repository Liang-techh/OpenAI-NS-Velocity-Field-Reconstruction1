"""Check continuous fixed-N route, quiet memory and independent Rc targets."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_cutoff_range_transport as current
import lei_ren_part1_paper_compliant_current_native_phase_first_jets_check as reference

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets=current.packets;ep=current.ep;require=reference.require;contains=reference.contains


def independent_target_checks(c):
    z=sy.symbols('z');a=2+z/3;mu=sy.Rational(1,10)
    f=dict(m=sy.Rational(1,5)+z/7,h=-sy.Rational(2,5)+z*z/11,
        e=sy.Rational(3,10)+z/13,p=-sy.Rational(1,10)+z*z/17)
    d=sy.Rational(3,10)+z/19
    f['k']=a*f['m']+mu*a*a*d
    expected={'M':f['m']/a,'I':f['h']/a,'S':f['e']/(a*a),'Cp':f['p']/(a*a),
        current.current.repair.ROWS[1]:d}
    bases=tuple(c.mpf(0) for _ in range(5));ledger=reference.qchecks.new_ledger()
    def scalar(expr,point):
        q=sy.simplify(expr.subs(z,point))
        return current.density.cutoff.prior.ScaledEnclosure(current.density.cutoff.prior.FormalScale(bases),
            c.mpf(int(q.p))/int(q.q),ledger)
    comparisons=0
    for point in (sy.Rational(-1,5),sy.Integer(0),sy.Rational(1,5)):
        history={key:scalar(expr,point) for key,expr in f.items()}
        jets={key:scalar(sy.diff(expr,z),point) for key,expr in f.items()}
        A,AZ=scalar(a,point),scalar(sy.diff(a,z),point)
        got=current.fixed_N_target_rows(history,jets,A,AZ,c.ln(c.mpf('1.9')),scalar(mu,point),c.ln(c.mpf('.1')))
        for group,order in (('values',0),('Z_derivatives',1)):
            for key,expr in expected.items():
                target=scalar(sy.diff(expr,z,order),point).finite_interval()
                require(contains(got[group][key].finite_interval()+c.mpf(('-1e-95','1e-95')),target),
                    'Independent normalized fixed-N target failed: '+group+' '+key)
                comparisons+=1
    return dict(passed=True,independent_exact_Rc_target_C0_Z_comparisons=comparisons,
        variable_nonzero_A_Z_and_joint_divided_row_used=True,
        reference='Independent rational source polynomials; exact symbolic differentiation of normalized targets',
        synthetic_reference_fixtures_not_original_field_values=True)


@current.native.inlet.source_precision
def run(role_owner,report=None,live=None):
    began=time.monotonic()
    owner=current.NativeCutoffRangeTransport(role_owner) if live is None else live[0]
    report=json.loads((HERE/current.NAME).read_bytes()) if report is None else report
    require(report[current.GATE] and report['source_family']==owner.family,'Same original fixed-N transport producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed fixed-N prerequisite: '+name)
    require(report['candidate_N']==2048 and report['enclosed_original_whole_cells']==24 and
        report['full24_original_C1_integral_range_transport_enclosed'] and not report['unresolved_original_whole_cell_labels'],
        'All24 original whole-cell fixed-N ranges required')
    cells=report['original_serial_cells']
    require([(q['label'],q['chart']) for q in cells]==[(q[0],q[1]) for q in current.current.ROUTE],
        'Original exact function route order differs')
    previous=None;rows=0;quiet_pressure=0;positive_masses=0
    for index,row in enumerate(cells):
        require(row['status']=='enclosed' and row['candidate_N']==2048,'Required continuous source cell unresolved')
        require(row['source_covers_entire_native_coordinate_and_Z_cell'] and
            row['cutoff_signed_u_and_phase_unions_before_one_positive_mass'],'Point provider or overlapping integral sum used')
        require(row['logarithmic_radius_measure_no_second_native_Jacobian'] and
            row['original_Z_independent_weights_and_endpoints_derivative_under_integral'],'Original integration units/Z rules differ')
        for group in ('C0_contributions','Z_contributions','incoming_C0','incoming_Z','outgoing_C0','outgoing_Z'):
            require(set(row[group])==set(current.RATES),'Missing original five-rate range rows');rows+=5
        for factor in row['true_log_radius_kernel_factors'].values():
            require(factor['mass']['sign']=='positive' and not factor['mass']['exact_zero'],'Original true-width mass is not positive')
            positive_masses+=1
        if previous is not None:
            require(row['incoming_C0']==previous['outgoing_C0'] and row['incoming_Z']==previous['outgoing_Z'],
                'Incoming correction history changed at route join')
        if index==0 or row['chart']=='O3_power':
            require(all(q['exact_zero'] for group in ('C0_contributions','Z_contributions') for q in row[group].values()),
                'Original quiet source modulation not exactly zero')
            require(row['incoming_C0']['p']==row['outgoing_C0']['p'] and row['incoming_Z']['p']==row['outgoing_Z']['p'],
                'Quiet zero-rate pressure memory was reset');quiet_pressure+=2
        if row['chart']!='bridge_first' or index!=0:
            source=row['original_whole_cell_source']
            require(source['status']=='enclosed' and 'separate_original_P0' in source and 'separate_original_P0_Z' in source,
                'Continuous original source or separate pressure datum missing')
        previous=row
    require(set(report['original_Rc_target_C0_ranges'])==set(current.current.repair.ROWS) and
        set(report['original_Rc_target_Z_ranges'])==set(current.current.repair.ROWS),'All five actual fixed-N target rows required')
    require(report['fixed_N_range_integrals_not_scalar_values_or_uniform_coefficients'] and
        report['independent_hulls_do_not_prove_joint_target_cancellation'],'Fixed-N range scope was widened')
    require(not any(report.get(k) for k in packets.OPEN) and not report['actual_five_controls_installed'] and
        not report['actual_terminal_Z_function_closure_installed'] and not report['current_whole_N_selected'],
        'Fixed-N range transport cannot complete controls/global N/recursion')
    live_rows=0
    if live is not None:
        got=live[1]
        require(got['record']['enclosed_original_whole_cells']==24 and got['target'] is not None,'Actual live route missing')
        for cell in got['cells']:
            for group in ('values','Z_derivatives'):
                for value in cell[group].values():
                    require(value.ctx is owner.ctx and value.ledger is owner.coordinates.ledger and value.scale.bases is owner.coordinates.bases,
                        'Live integral lost original common basis/context/ledger');live_rows+=1
        require(packets.encode(current.records(got['history']))==packets.encode(report['original_Rc_C0_history_ranges']),
            'Producer report differs from live fixed-N history')
    lower=owner.cell('additional_O2_at_lower_N','O2_slope','.12','.15',Z=(-1,1),N=160)
    require(lower['values'] is not None and lower['record']['status']=='enclosed','Actual lower-N continuous source unresolved')
    try:owner.route(N=159)
    except ValueError:pass
    else:raise AssertionError('Candidate N below original bound accepted')
    independent=independent_target_checks(owner.ctx)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=2048,
        continuous_original_route_cells_checked=24,exact_original_function_radius_joins=23,
        source_contribution_and_incoming_outgoing_C0_Z_rows_checked=rows,
        actual_live_common_basis_integral_rows_checked=live_rows,positive_original_kernel_mass_rows=positive_masses,
        exact_quiet_pressure_C0_Z_memory_rows=quiet_pressure,
        actual_original_lower_N_O2_source_and_integral=lower['record'],additional_original_candidate_N=160,
        independent_fixed_N_target_checks=independent,
        full24_original_C1_integral_range_transport_enclosed=True,
        fixed_N_original_Rc_target_C0_Z_ranges_installed=True,
        scalar_quadrature_or_quantitative_joint_cancellation_installed=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Continuous24-cell original fixed-N range transport, actual common-basis contributions, unchanged incoming/quiet pressure memory, lower-N O2 cell and independent normalized target derivatives. Broad ranges do not establish small defects, repair contraction, scalar values, terminal closure, global N or recursion.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original fixed-N range transport check PASS:24 cells,30 independent target comparisons',flush=True)
    return result
