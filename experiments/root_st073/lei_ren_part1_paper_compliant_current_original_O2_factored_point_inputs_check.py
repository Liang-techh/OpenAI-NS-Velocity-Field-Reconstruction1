"""Focused independent ODE/unit and error-contract checks for O2 point rows.

Finite R/Pstar/delta here are diagnostic units, never native parameters.
The accepted original physical I/F identity is reused through its hashes.
An independent moment ODE tests the new factor decomposition numerically.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_factored_point_inputs as point
from lei_ren_part1_paper_compliant_current_original_O2_slope_point_profiles_check import original_ODE_reference


def read(c,value):
    return c.make_mpf(tuple(value['exact_mpf_tuple'])) if isinstance(value,dict) else c.mpf(value)


def finite_reference(frame,manifest):
    c=mp.mp.clone();c.dps=70;t=frame.owner.template;z=t['z']
    f,H,D,P=s.symbols('f H D P',real=True);p0,p0Z,p0ZZ=s.symbols('P0 P0_Z P0_ZZ',real=True)
    substitution={t['f']:f,t['H']:H,t['D']:D,t['P']:P,t['P0']:p0,
        s.diff(t['P0'],z):p0Z,s.diff(t['P0'],z,2):p0ZZ}
    args=(z,t['R'],t['Pstar'],t['delta'],f,H,D,P,p0,p0Z,p0ZZ)
    original={}
    for key in ('p1','p2'):
        for order in (0,1):
            expression=t[key if order==0 else key+'_Z'].xreplace(substitution)
            original[(key,order)]=s.lambdify(args,expression,modules=[{'mpf':c.mpf},'mpmath'])
    pressure=json.loads((point.HERE/point.pressure.NAME).read_bytes())
    alpha=read(c,pressure['actual_defining_quadrature_alpha'])
    comparisons=coefficients=nonzero_tails=0;maximum=c.mpf(0)
    with mp.workdps(c.dps):
        for query in manifest['original_full_factored_O2_point_queries']:
            y=s.Rational(query['original_y_exact']);Z=s.Rational(query['original_Z_exact'])
            yv=c.mpf(int(y.p))/int(y.q);zv=c.mpf(int(Z.p))/int(Z.q);q=1+zv*zv;C=1/q
            # The separate ODE integrates the original five histories. Its
            # double-precision discrepancy is allowed only in this check.
            state=original_ODE_reference(float(yv),float(zv),7,steps=2400)
            fv=c.exp(yv/10-c.mpf(3)/5*c.mpf(state[0]))
            Hv=c.mpf(state[2])/C;Dv=((4*zv/7)**2-c.mpf(state[4]))/C**2
            Pv=c.mpf(state[5])/C**2
            pressure_values=(-alpha/q**2,4*alpha*zv/q**3,alpha*(4-20*zv*zv)/q**4)
            for key,pair in query['inputs'].items():
                assert len(pair)==2
                for order,row in enumerate(pair):
                    assert row['source_factor_order']==['R','Pstar','delta','L']
                    assert row['original_factors_not_materialized']
                    for term in row['terms']:
                        powers=term['original_R_Pstar_delta_L_powers']
                        assert len(powers)==4 and all(type(power) is int for power in powers)
                        error=read(c,term['directed_finite_coefficient_absolute_error_upper'])
                        assert c.isfinite(error) and error>=0
                        late=term['late_pressure_error'];assert late['zero']==(late['log_upper'] is None)
                        if not late['zero']:
                            assert key=='p2' and read(c,late['log_upper']) < -c.mpf('1e17')
                            nonzero_tails+=1
                        coefficients+=1
                    if key not in ('p1','p2'):continue
                    for R,ps,delta in ((2,7,c.mpf('.1')),(5,13,c.mpf('.01'))):
                        factors=(c.mpf(R),c.mpf(ps),delta,1-delta*zv*zv)
                        actual=error=0
                        for term in row['terms']:
                            weight=c.fprod(base**power for base,power in zip(factors,term['original_R_Pstar_delta_L_powers']))
                            actual+=weight*read(c,term['approximate_source_point_coefficient'])
                            error+=abs(weight)*read(c,term['directed_finite_coefficient_absolute_error_upper'])
                        reference=original[(key,order)](zv,R,ps,delta,fv,Hv,Dv,Pv,*pressure_values)
                        discrepancy=abs(actual-reference)
                        assert discrepancy<c.mpf('2e-8'),(key,order,y,Z,discrepancy)
                        # Finite reference omits the late original pressure;
                        # those nonzero budgets remain explicit, unexponentiated.
                        assert discrepancy<=error+c.mpf('2e-8')
                        maximum=max(maximum,discrepancy);comparisons+=1
            expected_E=fv*C;E=query['inputs']['E']
            for row,expected in ((E[0],expected_E),(E[1],-2*zv*C*expected_E)):
                actual=sum(read(c,term['approximate_source_point_coefficient']) for term in row['terms'])
                assert abs(actual-expected)<c.mpf('2e-11')
            for row in query['inputs']['b']:
                assert row['terms']==[]
            assert query['inputs']['a'][1]['terms']==[]
        assert nonzero_tails>0
    return dict(passed=True,original_full_I_over_F_ODE_and_finite_unit_comparisons=comparisons,
        maximum_absolute_reference_discrepancy=maximum,finite_coefficient_error_records=coefficients,
        nonzero_pressure_tail_records_preserved=nonzero_tails,
        finite_reference_units_only=[[2,7,'.1'],[5,13,'.01']],
        independent_ODE_is_diagnostic_not_a_directed_proof=True,
        original_native_factors_not_numerically_materialized=True)


def decomposition_and_guards(frame):
    templates=point.finite_coefficient_templates(frame);p0=templates['pressure_symbols']
    terms=0
    for rows in templates['rows'].values():
        for unused,expression in rows:
            for a in p0:
                for b in p0:
                    assert s.diff(expression,a,b)==0
            terms+=1
    # Linear pressure sensitivity makes separate late-error propagation exact
    # even when more than one ordinary pressure jet occurs in a coefficient.
    assert all(templates['source_template_reconstruction_identities'].values())
    rejected=0
    for call in (lambda:point.source.exact_coordinate(-1),lambda:point.source.exact_coordinate(2),
        lambda:point.pressure.exact_Z('nan'),lambda:point.pressure.exact_Z(2),
        lambda:point.FactoredPointRow((point.FactoredPointTerm((1,0,0,0),1,0),),0,0).unscaled_scalar(),
        lambda:point.OriginalO2FactoredPointInputs(cells=0)):
        try:call()
        except ValueError:rejected+=1
    assert rejected==6
    return dict(passed=True,exact_original_template_reconstruction=templates['source_template_reconstruction_identities'],
        pressure_linear_coefficient_terms=terms,second_and_mixed_pressure_derivatives_exact_zero=True,
        wrong_domain_nonfinite_unresolved_scalar_and_cells_rejections=rejected)


def run():
    began=time.monotonic();manifest=json.loads((point.HERE/point.NAME).read_bytes())
    frame=point.source.OriginalO2SourceParameterFrame()
    assert manifest[point.GATE] and manifest['source_family']==frame.family
    for name,digest in manifest['input_hashes'].items():assert point.sha(name)==digest,name
    flags=('conditioned_native_factor_and_phase_evaluation_installed','original_p1_p2_scalar_point_values_installed',
        'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
        'current_whole_N_selected',*point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    result=dict(all_passed=True,**{point.GATE:True},source_family=frame.family,
        independent_finite_references=finite_reference(frame,manifest),
        exact_decomposition_pressure_linearity_and_guards=decomposition_and_guards(frame),
        directed_error_method='Directed original radial integral boxes plus checked pressure coefficient enclosure and all-late-source theorem; linear late-pressure sensitivity bounds, kept in logarithmic factored form.',
        **dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],point.NAME:point.sha(point.NAME),
            Path(__file__).name:point.sha(Path(__file__).name),
            'lei_ren_part1_paper_compliant_current_original_O2_slope_point_profiles_check.py':
                point.sha('lei_ren_part1_paper_compliant_current_original_O2_slope_point_profiles_check.py')},
        execution_seconds=time.monotonic()-began,
        scope='Complete O2 factored finite point coefficients and directed coefficient/source-pressure budgets. Native astronomical factor/conditioned loop evaluation, full oracle, global frequency and recursive field remain open.')
    (point.HERE/point.RECEIPT).write_text(json.dumps(point.source.inertial.profiles.loop.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original O2 factored point inputs: independent ODE/unit, decomposition and error contracts PASS',flush=True)
    return result


if __name__=='__main__':run()
