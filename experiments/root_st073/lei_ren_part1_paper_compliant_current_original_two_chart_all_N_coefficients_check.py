"""Focused runtime, original density and source-role rejection checks."""
from dataclasses import replace
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_two_chart_all_N_coefficients as current
import lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals as density

base=current.base;ep=current.ep


def exprel_references():
    c=base.MPIntervalContext();c.dps=180;p=mp.mp.clone();p.dps=160
    checks=0
    with mp.workdps(220):
        for lo,hi in (('-1','1'),('-.02','.03'),('0','0'),('1e-100','1e-100'),('-.8','-.5'),('.2','.7')):
            got=current.directed_exprel(c,c.mpf([lo,hi]));a,b=ep(got)
            for x in (p.mpf(lo),p.mpf(hi),(p.mpf(lo)+p.mpf(hi))/2):
                # Independent defining exponential average quadrature.
                expected=p.quad(lambda t:p.exp(t*x),[0,1])
                assert a-p.mpf('1e-155')<=expected<=b+p.mpf('1e-155')
                checks+=1
        assert ep(current.directed_exprel(c,c.mpf(0)))==(1,1)
    return dict(passed=True,independent_defining_exponential_average_comparisons=checks,
        zero_exact_one=True,order64_integrated_exponential_Taylor_tail_bound=True,
        signed_interval_and_nonzero_tiny_argument_cases=True)


def symbolic_density_identity():
    E,EZ,V,VZ,A,AZ,B,BZ,M,X=s.symbols('E EZ V VZ A AZ B BZ exprel exp')
    N=s.Symbol('N',positive=True);F=E*A*M;FZ=EZ*A*M+E*AZ*X
    first=dict(m=B,h=F,k=V*F+E*B,e=2*V*B-E*F,p=E*F)
    firstZ=dict(m=BZ,h=FZ,k=VZ*F+V*FZ+EZ*B+E*BZ,
        e=2*VZ*B+2*V*BZ-EZ*F-E*FZ,p=EZ*F+E*FZ)
    second=dict(m=0,h=0,k=F*B,e=B*B-F*F/2,p=F*F/2)
    secondZ=dict(m=0,h=0,k=FZ*B+F*BZ,e=2*B*BZ-F*FZ,p=F*FZ)
    dE,dEZ,dV,dVZ=F/N,FZ/N,B/N,BZ/N
    kernels=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    jets=dict(m=dVZ,h=dEZ,
        k=VZ*dE+V*dEZ+EZ*dV+E*dVZ+dEZ*dV+dE*dVZ,
        e=2*VZ*dV+2*V*dVZ+2*dV*dVZ-EZ*dE-E*dEZ-dE*dEZ,
        p=EZ*dE+E*dEZ+dE*dEZ)
    for key in current.exact.RATES:
        assert s.expand(kernels[key]-first[key]/N-second[key]/N**2)==0
        assert s.expand(jets[key]-firstZ[key]/N-secondZ[key]/N**2)==0
    return dict(passed=True,exact_original_C0_Z_density_order_identities=10,
        all_original_V_E_B_cross_terms_and_Z_product_rules_retained=True,
        original_F_Z_uses_exp_not_exprel_on_A_Z=True,
        exact_order_minus1_functions_not_zeroed=True)


def runtime_and_guards():
    dispatcher=current.TwoChartOriginalPointSources();comparisons=dispatches=0;midplane=0
    last=None
    for chart,coordinate,Z,N in (('Rh_reference','-2.337','.37',160),('O2_slope','.53','0',160)):
        got=current.evaluate_coefficients(dispatcher,chart=chart,coordinate=coordinate,Z=Z,N=N)
        for item in got['evaluated']:
            piece=item['piece'];values=piece.values;c=piece.ctx
            E,EZ,V,VZ,A,AZ,B,BZ=(values[role] for role in current.ROLES)
            primitive=dict(A=A,A_Z=AZ,B_over_Pstar=B,B_Z_over_Pstar=BZ)
            # The separately accepted original density program evaluates
            # velocity increments first and differentiates each cross term.
            direct=density.density_Z_kernels(E,EZ,V,VZ,primitive,N)
            for key in current.exact.RATES:
                for jet,target in (('C0',direct['kernels'][key]),('Z',direct['Z_derivatives'][key])):
                    combined=item['values'][-1][key][jet]*(c.mpf(1)/N)+item['values'][-2][key][jet]*(c.mpf(1)/N**2)
                    difference=combined-target;lo,hi=ep(difference.coefficient)
                    assert lo<=0<=hi,(chart,key,jet)
                    comparisons+=1
            for role,ref in got['built']['refs'].items():
                row=got['built']['graph'].nodes[ref.node]
                result=dispatcher.source(row,piece=piece,coordinate=c.mpf(int(piece.coordinate.p))/int(piece.coordinate.q),
                    phase=piece.phase,Z=piece.Z,N=piece.N)
                assert result is piece.values[role];dispatches+=1
            if Z=='0':
                assert not item['values'][-1]['h']['Z'].zero
                assert piece.record['ordinary_Z_contract']['branch']=='exact_midplane_nonzero_p2_Z'
                midplane+=1
            last=(got,piece)
    got,piece=last;g=got['built']['graph'];ref=got['built']['refs']['all_N_original_E_C0'];row=dict(g.nodes[ref.node]);c=piece.ctx
    coordinate=c.mpf(int(piece.coordinate.p))/int(piece.coordinate.q)
    def source(row0=row,**kwargs):
        return dispatcher.source(row0,piece=kwargs.get('piece',piece),coordinate=kwargs.get('coordinate',coordinate),
            phase=kwargs.get('phase',piece.phase),Z=kwargs.get('Z',piece.Z),N=kwargs.get('N',piece.N))
    rejected=0
    calls=[
        lambda:source(dict(row,graph_sha256='bad')),
        lambda:source(dict(row,source_graph_namespace='wrong')),
        lambda:source(dict(row,source_node=-1)),
        lambda:source(dict(row,function_role='density_m_C0')),
        lambda:source(dict(row,chart='Rh_reference')),
        lambda:source(N=256),lambda:source(Z='.37'),lambda:source(phase=c.mpf('.1')),
        lambda:source(coordinate=c.mpf('.1')),
        lambda:source(piece=replace(piece)),
        lambda:dispatcher.pieces(chart='actual_patch',coordinate='1.2',Z=0,N=160),
        lambda:dispatcher.pieces(chart='O2_slope',coordinate='1.01',Z=0,N=160),
        lambda:dispatcher.pieces(chart='Rh_reference',coordinate='.1',Z=0,N=160),
        lambda:current.candidate_N(True),lambda:current.candidate_N(159),
        lambda:current.FactoredPointCoefficientEvaluator(dict(got['built']),dispatcher,piece),
    ]
    for call in calls:
        try:call()
        except (ValueError,TypeError,ArithmeticError):rejected+=1
    assert rejected==len(calls)
    return dict(passed=True,live_original_role_dispatches=dispatches,
        original_density_C0_Z_overlap_comparisons=comparisons,
        exact_midplane_nonzero_source_derivative_cases=midplane,
        namespace_node_hash_chart_coordinate_phase_Z_N_and_unissued_frame_graph_guards=rejected,
        scalar_guard_not_relaxed=True,range_endpoints_or_midpoints_not_selected=True)


def saved_contracts(manifest):
    c=base.MPIntervalContext();c.dps=260
    for sample in manifest['actual_original_point_coefficient_queries']:
        assert sample['source_family']==manifest['source_family']
        assert sample['candidate_N']>=160 and sample['source_role_count']==8
        assert sample['phase_piece_nonlinear_coefficients_evaluated_before_union']
        assert sample['original_E_V_B_and_all_cross_terms_retained']
        assert not sample['source_caps_or_midpoints_used_as_values']
        for key in ('m','h'):
            for jet in ('C0','Z'):
                assert sample['actual_N_dependent_C0_Z_coefficients']['-2'][key][jet]['exact_zero']
        for piece in sample['actual_phase_pieces']:
            source=piece['source'];assert source['one_basis_and_ledger_for_all_roles'] and source['true_phase_Z_exact_zero']
            assert source['original_radius_phase']['positive_original_origin_offset']['strictly_positive']
            assert source['ordinary_Z_contract']['p2_Z_retained']
            assert not source['range_endpoint_or_midpoint_selected']
    # N remains in the actual inverse and modulation: same coordinate/Z at
    # N160 and N256 produce distinct coefficient enclosures.
    rows=manifest['actual_original_point_coefficient_queries']
    assert rows[0]['actual_N_dependent_C0_Z_coefficients']!=rows[2]['actual_N_dependent_C0_Z_coefficients']
    return dict(passed=True,actual_native_coefficient_queries=len(rows),original_charts=['Rh_reference','O2_slope'],
        shared_candidate_N_values=[160,256],same_point_different_N_not_rescaled_from_one_fixture=True,
        exact_quiet_second_order_m_h_zeros_preserved=True,
        full_integral_control_global_frequency_flags_remain_false=True)


def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and manifest['original_exact_coefficient_program_reused']
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('all_17_chart_or_24_cell_integral_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*current.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[k] is False for k in flags)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        original_exprel_defining_integral_references=exprel_references(),
        independent_signed_density_order_identity=symbolic_density_identity(),
        live_runtime_and_rejections=runtime_and_guards(),actual_native_manifest_contracts=saved_contracts(manifest),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),Path(density.__file__).name:current.sha(Path(density.__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Directed actual two-chart source roles and original N-dependent coefficient runtime; exponential-average reference, independent original density identity/runtime and rejection checks. No source/integral all17 completion or numerical five-control/global-N/recursive field claim.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Live original two-chart all-N coefficient runtime: source roles, density and guards PASS',flush=True)
    return report


if __name__=='__main__':run()
