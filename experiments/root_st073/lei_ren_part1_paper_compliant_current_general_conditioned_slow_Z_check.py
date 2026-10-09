"""Independent scalar quadrature, finite differences and source-factor checks."""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_general_conditioned_slow_Z as current
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as scalar

base=current.base;ep=current.ep


def read_iv(c,row):
    return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                  mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))


def restore(row,bases,ledger):
    c=bases[0].ctx;s=row['formal_positive_scale']
    scale=current.prior.FormalScale(bases,tuple(s['source_exponents'])+(s['radius_power'],),
        read_iv(c,s['additional_log_interval']))
    return current.prior.ScaledEnclosure(scale,read_iv(c,row['coefficient_interval']),ledger)


def original_specialization(report):
    owner=current.old.OriginalO2ConditionedSlowZ();components=0;midplane=endpoint=False
    with mp.workdps(owner.owner.ctx.dps+40):
        for row in report['genuine_original_O2_general_backend_queries']:
            assert row['source_family']==owner.family and not row['native_caps_or_midpoints_selected_as_field_values']
            assert row['accepted_source_owner_and_auxiliary_logs_bound']
            assert row['original_source_basis_contract']['same_original_R_Pstar_delta_L_bound']
            query=owner.owner.query(y=row['y'],Z=row['Z']);k=query['kernel'];roots=query['roots'];c=k.c
            for record in row['general_derivatives_at_genuine_original_O2_points']:
                assert not record['free_fixed_phase_parameter']
                assert record['C0']['original_common_N_and_radius_phase_bound']
                assert not record['C0']['free_phase_parameter_not_spatial_phase']
                select=record['C0']['selected_inverse'];x=read_iv(c,select['coordinate_interval'])
                old,_=current.old.slow_values(k,roots,x,select['chart'])
                for key in old:
                    new=restore(record['slow_Z'][key],k.q.scale.bases,k.q.ledger)
                    assert ep((new-old[key]).coefficient)[0]<=0<=ep((new-old[key]).coefficient)[1],key
                    assert not record['slow_Z'][key]['point_value_selected'];components+=1
                proof=record['derivative_contract']
                assert proof['C0_q_and_differentiated_original_root_sigma_identical']
                assert proof['q_Z']['exact_zero'] and proof['t0_Z']['exact_zero']
                if row['Z']=='0':
                    assert not proof['u_Z']['exact_zero'] and not record['slow_Z']['A_Z_slow']['exact_zero'];midplane=True
                if row['y']=='1':
                    assert not record['C0']['geometry']['q']['exact_zero'];endpoint=True
    assert components==12 and midplane and endpoint
    return dict(passed=True,genuine_original_O2_specialization_component_overlaps=components,
        original_midplane_p2_Z_and_positive_tiny_endpoint_q_preserved=True,
        independent_accepted_reduced_derivative_backend_used=True)


def finite_references():
    scales=scalar.loop.GenericLoopScales(a_min='1',margin_min='2',boundary_kappa_excess_min='.02',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='1',dps=110)
    p=scales.ctx;eta=scales.eta
    cases=[('small_positive',dict(a='1.5',b='.2',p2='.002',E='2')),
           ('signed_positive',dict(a='1.5',b='.2',p2='.7',E='2')),
           ('signed_negative',dict(a='1.5',b='.2',p2='-.7',E='2')),
           ('midplane',dict(a='1.5',b='.2',p2='0',E='2')),
           ('cutoff',dict(a='1.5',b=p.sqrt(p.mpf('1.5')*(p.mpf('.5')+eta*p.mpf('.45'))),p2='.7',E='2')),
           ('kappa_two',dict(a='2',b='0',p2='.7',E='2')),
           ('flat',dict(a='2.1',b='.2',p2='.7',E='2'))]
    derivatives=dict(a='.03',b='-.04',p2='.02',E='.1')
    records=[];comparisons=finite_difference_checks=0;branches=set();cutoffs=set();maximum=p.mpf(0)
    with mp.workdps(220):
        for name,values in cases:
            values={key:p.mpf(v) for key,v in values.items()}
            ref=scalar.GenericLoopPointZ(scales,p1=6,**values,
                **{key+'_Z':v for key,v in derivatives.items()})
            c,roots=current.fixture(values,derivatives,dps=160)
            backend=current.GeneralConditionedSlowZ(roots,eta_log=c.ln(c.mpf(eta)),
                dstar_log=c.ln(c.mpf(scales.d_star)),log_a_lower=c.ln(c.mpf(scales.a_min)))
            for key,want in (('q_Z',ref.q_Z),('t0_Z',ref.t0_Z),('kappa_Z',ref.kappa_Z)):
                lo,hi=ep(backend.shear[key].finite_interval());allow=p.mpf('1e-100')*(1+abs(want))
                assert lo-allow<=want<=hi+allow,(name,key);comparisons+=1
            for phase in ('.137','.337','.663'):
                wanted=ref.evaluate(phase)
                charts=('psi','E') if backend.kernel.geometry=='signed_Mobius' else ('psi',)
                for chart in charts:
                    if backend.kernel.flat:x=c.mpf(phase)
                    else:x=backend.kernel.inverse_bracket(c.mpf(phase),chart,110)['coordinate_interval']
                    actual,proof=backend.values(x,chart);branches.add(proof['branch']);cutoffs.add(proof['original_cutoff_derivative'])
                    for key,want in (('psi_Z',wanted.angle_Z),('A_Z_slow',wanted.A_Z_slow),('B_Z_slow',wanted.B_Z_slow)):
                        lo,hi=ep(actual[key].finite_interval());allow=p.mpf('1e-92')*(1+abs(want))
                        assert lo-allow<=want<=hi+allow,(name,phase,chart,key,p.nstr(want,15),mp.nstr(lo,15),mp.nstr(hi,15))
                        maximum=max(maximum,lo-want,want-hi,p.mpf(0));comparisons+=1
                if name in ('small_positive','signed_positive','signed_negative','midplane','cutoff') and phase=='.337':
                    h=p.mpf('1e-24');points=[]
                    for sign in (-1,1):
                        moved={key:v+sign*h*p.mpf(derivatives[key]) for key,v in values.items()}
                        points.append(scalar.loop.GenericShearLoop(scales,p1=6,Utheta=moved.pop('E'),**moved).evaluate(phase))
                    for key,want in (('A',wanted.A_Z_slow),('B',wanted.B_Z_slow)):
                        estimate=(points[1][key]-points[0][key])/(2*h)
                        assert abs(estimate-want)<=p.mpf('1e-36')*(1+abs(want)),(name,key)
                        finite_difference_checks+=1
            for phase in ('0','.5','1'):
                actual,_=backend.values(c.mpf(phase));assert all(value.zero for value in actual.values())
            records.append(dict(case=name,manufactured_functions_only=True,geometry=backend.kernel.geometry,
                cutoff=backend.shear['cutoff_contract'],nonzero_b=not roots['b'][current.ZERO].zero))
    assert {'general_exact_midplane','general_small_r_Fourier','general_signed_psi','general_signed_E','exact_general_flat_or_symmetry'}<=branches
    assert 'original_directed_sigma_and_sigma_prime' in cutoffs
    return dict(passed=True,explicit_function_cases=records,independent_scalar_quadrature_component_comparisons=comparisons,
        independent_C0_finite_difference_checks=finite_difference_checks,branches=sorted(branches),
        cutoff_derivative_branches=sorted(cutoffs),maximum_reference_outside_discrepancy=maximum,
        finite_fixture_scales_not_used_as_native_parameters=True)


def uncertain_source_box():
    values=dict(a=('1.49','1.51'),b=('.19','.21'),p2=('.69','.71'),E=('1.99','2.01'))
    jets=dict(a=('.029','.031'),b=('-.041','-.039'),p2=('.019','.021'),E=('.099','.101'))
    scales=scalar.loop.GenericLoopScales(a_min='1',margin_min='2',boundary_kappa_excess_min='.02',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='1',dps=100)
    c,roots=current.fixture(values,jets,dps=150)
    backend=current.GeneralConditionedSlowZ(roots,eta_log=c.ln(c.mpf(scales.eta)),
        dstar_log=c.ln(c.mpf(scales.d_star)),log_a_lower=c.mpf(0))
    inverse=backend.kernel.inverse_bracket(c.mpf('.337'),'E',100)
    actual,_=backend.values(inverse['coordinate_interval'],'E');checks=0
    for pick in ((0,0,0,0),(1,1,1,1),(0,1,1,0),(1,0,0,1)):
        values1={key:values[key][i] for key,i in zip(('a','b','p2','E'),pick)}
        jets1={key:jets[key][i] for key,i in zip(('a','b','p2','E'),reversed(pick))}
        ref=scalar.GenericLoopPointZ(scales,p1=6,**values1,**{key+'_Z':v for key,v in jets1.items()}).evaluate('.337')
        for key,want in (('psi_Z',ref.angle_Z),('A_Z_slow',ref.A_Z_slow),('B_Z_slow',ref.B_Z_slow)):
            lo,hi=ep(actual[key].finite_interval());assert lo<=want<=hi,(pick,key);checks+=1
    return dict(passed=True,uncertain_source_and_jet_box_corner_comparisons=checks,
        directed_source_errors_propagated_without_point_selection=True)


def extreme_factors_and_guards():
    extremes=[]
    for sign in (1,-1):
        c,roots=current.fixture(dict(a=2,b=0,p2=str(sign*.4),E=1),dict(a='.03',b='.04',p2='.02',E='.1'),
            dps=200,bases_values=('1e20',0,0,0,0),source_powers={('a',current.ZROW):(-1,0,0,0,0),
                ('E',current.ZERO):(1,0,0,0,0),('E',current.ZROW):(1,0,0,0,0)})
        with mp.workdps(240):
            backend=current.GeneralConditionedSlowZ(roots,eta_log=c.mpf('-1e20'),dstar_log=c.mpf('-1e20'),log_a_lower=c.mpf(0))
            assert backend.kernel.geometry=='signed_Mobius' and not backend.kernel.q.zero and not backend.shear['q_Z'].zero
            for factor in (backend.kernel.rho,backend.kernel.s_source,backend.kernel.hinv):
                assert ep(factor.coefficient)[0]>0 and not factor.zero
            for chart in ('psi','E'):
                bracket=backend.kernel.inverse_bracket(c.mpf('.137'),chart,90)
                values,proof=backend.values(bracket['coordinate_interval'],chart)
                assert all(not value.zero for value in values.values())
                extremes.append(dict(sign=sign,chart=chart,proof=proof,derivatives={key:value.record() for key,value in values.items()},
                    source_geometry=backend.kernel.geometry_record(),no_huge_source_exponentials_materialized=True))
    c,roots=current.fixture(dict(a='1.5',b='.2',p2='.002',E=2),dict(a='.03',b='-.04',p2='.02',E='.1'))
    factory=lambda rows:current.GeneralConditionedSlowZ(rows,eta_log=c.ln(c.mpf('.01')),dstar_log=c.ln(c.mpf('.5')),log_a_lower=c.mpf(0))
    backend=factory(roots);guards=0
    def rejected(fn):
        nonlocal guards
        try:fn()
        except (ValueError,ArithmeticError,KeyError):guards+=1
        else:raise AssertionError('Invalid derivative request accepted')
    rejected(lambda:backend.values('-0.01'))
    rejected(lambda:backend.values('1.01'))
    rejected(lambda:backend.values('.2','E'))
    rejected(lambda:backend.evaluate('1.1'))
    rejected(lambda:backend.evaluate('.2',bits=3))
    _,alien=current.fixture(dict(a=1,b=0,p2=0,E=1),dict(a=0,b=0,p2=0,E=0))
    changed={key:dict(row) for key,row in roots.items()};changed['b'][current.ZROW]=alien['b'][current.ZROW]
    rejected(lambda:factory(changed))
    changed={key:dict(row) for key,row in roots.items()};changed['a'][current.ZERO]=roots['a'][current.ZERO].scalar(-1)
    rejected(lambda:factory(changed))
    changed={key:dict(row) for key,row in roots.items()};changed['a'][current.ZERO]=roots['a'][current.ZERO].scalar(('2.008','2.012'));changed['b'][current.ZERO]=roots['b'][current.ZERO].scalar(0)
    rejected(lambda:factory(changed))
    changed={key:dict(row) for key,row in roots.items()};changed['p2'][current.ZERO]=roots['p2'][current.ZERO].scalar(('-.7','.7'))
    rejected(lambda:factory(changed))
    backend.kernel.source=dict(backend.query);rejected(lambda:backend.values('.2'))
    assert guards==10
    return dict(passed=True,extreme_factored_derivative_records=extremes,invalid_requests_rejected=guards,
        positive_tiny_q_and_q_Z_rho_s_hinv_preserved=True,active_flat_mixed_boxes_require_refinement=True,
        exact_flat_branch_alone_sets_q_Z_zero=True)


def run():
    began=time.monotonic();report=json.loads((current.HERE/current.NAME).read_bytes())
    assert report[current.GATE] and report['no_ancestor_quadratures_or_producers_executed']
    assert report['generic_backend_requires_caller_owned_source_cone_and_scale_hypotheses']
    assert report['genuine_O2_path_binds_accepted_point_errors_family_logs_and_true_phase']
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    assert all(report[key] is False for key in current.OPEN)
    assert current.exact_identities()['passed'] and report['exact_general_derivative_identities']['passed']
    result=dict(all_passed=True,**{current.GATE:True},source_family=report['source_family'],
        genuine_original_O2_specialization=original_specialization(report),
        general_varying_source_scalar_references=finite_references(),
        uncertain_source_error_box_checks=uncertain_source_box(),
        extreme_source_factors_and_guards=extreme_factors_and_guards(),
        **dict.fromkeys(current.OPEN,False),input_hashes={**report['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),Path(scalar.__file__).name:current.sha(Path(scalar.__file__).name),
            Path(scalar.loop.__file__).name:current.sha(Path(scalar.loop.__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='General fixed-phase derivative backend and genuine O2 specialization. Independent finite scalar quadrature and C0 finite differences check variable a,b,q,t0; native active-patch defining histories, full integral closure and full NS reconstruction remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('General factored Z: independent varying-source references, O2 specialization and guards PASS',flush=True)
    return result


if __name__=='__main__':run()
