"""Independent full-density pressure reference and native error-budget checks."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_pressure_point_jets as point


def full_density_reference(owner):
    """Integrate the actual 14 defining densities in finite diagnostic units.

    This fixture checks the general all-late-stage lemma, not native parameter
    selection. It uses the original cutoffs/primitive and raw waiting root.
    """
    operator=owner.frame.owner.pressure;p=operator.partition;z,t=p['z'],p['t'];c=owner.ctx
    parameters={p['mu']:s.Rational(1,5),p['delta']:s.Rational(1,10),p['epsilon']:s.Rational(1,1000),
        p['yd']:s.Integer(3),p['Tw']:s.Rational(3,10),p['L']:s.Rational(2,5),p['Ts']:s.Integer(2)}
    sigma=lambda x:point.source.inertial.profiles.loop.flat_step(c,c.mpf(x))
    def J(x):
        x=c.mpf(x)
        if x<=0:return c.mpf(0)
        if x>=1:return x-c.mpf('.5')
        return owner.frame.owner.profiles.J(x)
    def phi(x):
        x=c.mpf(x)
        return c.mpf(0) if x>=3 else c.exp(-4/(3-x)**2)
    functions=[dict(J_sigma=J,sigma=sigma,phi=phi),'mpmath']
    def integrate(expression,variable,lo,hi):
        fn=s.lambdify(variable,expression,modules=functions)
        a=-mp.inf if lo==-s.oo else mp.mpf(str(lo))
        b=mp.inf if hi==s.oo else mp.mpf(str(hi))
        knots=[a]+[mp.mpf(v) for v in ('.5','1','2','3','10','30','50') if a<mp.mpf(v)<b]+[b]
        return mp.quad(fn,knots)
    with mp.workdps(25):
        root=operator.raw_waiting_root.subs(parameters)
        replacement={}
        for integral in root.atoms(s.Integral):
            variable,lo,hi=integral.limits[0]
            value=integrate(integral.function,variable,lo,hi)
            replacement[integral]=s.Float(str(value),55)
        W=mp.mpf(str(root.xreplace(replacement).evalf(50)))
        assert W>0
        parameters[p['W']]=s.Float(str(W),55)
        Z=s.Rational(37,100);coordinate=mp.mpf('.37');q=1+coordinate**2
        derivative_factors=(q**-2,-4*coordinate/q**3,(-4+20*coordinate**2)/q**4)
        errors=[];integrals=0
        for order,remainder_factor in enumerate((5,10,44)):
            actual=mp.mpf(0)
            for name,density in operator.original_densities.items():
                integrand=s.diff(density,z,order).subs(parameters).subs(z,Z)
                if integrand==0:continue
                lo,hi=(s.sympify(bound).subs(parameters).evalf(50) for bound in p['domains'][name])
                actual-=integrate(integrand,t,lo,hi);integrals+=1
            approximation=-owner.alpha*derivative_factors[order]
            error=abs(actual-approximation)
            bound=remainder_factor*mp.exp(mp.mpf('.6')-3)/(2*q*q)
            assert error<bound,(order,error,bound)
            errors.append(dict(order=order,full_original_normalized_pressure=actual,
                actual_error=error,all_late_source_bound=bound))
    return dict(passed=True,nonzero_original_density_integrals=integrals,
        all_fourteen_stages_in_value_reference=True,original_density_Z_integrals_through2=True,
        finite_diagnostic_parameters=dict(mu='.2',delta='.1',epsilon='.001',yd=3,Tw='.3',L='.4',Ts=2),
        approximate_original_raw_waiting_root=W,reference_Z='.37',comparisons=errors,
        finite_fixture_tests_general_source_lemma_not_native_parameter_selection=True,
        diagnostic_reference_decimal_precision=25)


def run():
    began=time.monotonic();owner=point.OriginalNormalizedPressurePointJets()
    report=json.loads((point.HERE/point.NAME).read_bytes())
    assert report[point.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():assert point.sha(name)==digest,name
    assert owner.proof==report['original_late_source_error_theorem']
    queries=[];maximum_error=mp.mpf(0)
    for Z in ('-.7','0','.37','1'):
        row=owner.evaluate(Z);queries.append(row)
        for jet in row['normalized_pressure_ordinary_Z_jets']:
            error=jet['coefficient_quadrature_and_arithmetic_absolute_error_upper']
            assert error>=0;maximum_error=max(maximum_error,error)
            tail=jet['finite_end_and_all_late_source_error']
            if tail['zero']:
                assert jet['ordinary_Z_order']==1 and Z=='0' and tail['representation']=='exact_zero'
                assert jet['approximate_normalized_value']==error==0
            else:assert tail['log_upper']<-mp.mpf('2e17')
        assert row['physical_Pstar_squared_error_amplification_not_discarded']
        assert not row['full_native_phase_and_global_scalar_oracle_installed']
    plus=owner.evaluate('.37')['normalized_pressure_ordinary_Z_jets']
    minus=owner.evaluate('-.37')['normalized_pressure_ordinary_Z_jets']
    for order in (0,1,2):assert plus[order]['approximate_normalized_value']==(-1)**order*minus[order]['approximate_normalized_value']
    rejected=0
    for fn in (lambda:owner.evaluate(2),lambda:owner.evaluate(-2),lambda:point.OriginalNormalizedPressurePointJets(cells=0),
        lambda:owner.evaluate(owner.ctx.nan),lambda:owner.evaluate(owner.ctx.inf)):
        try:fn()
        except ValueError:rejected+=1
    assert rejected==5
    independent=full_density_reference(owner)
    flags=('original_p1_p2_scalar_point_values_installed','numerical_original_source_point_or_integral_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*point.source.inertial.profiles.loop.OPEN)
    assert all(report[key] is False for key in flags)
    result=dict(all_passed=True,**{point.GATE:True},source_family=owner.family,
        original_late_source_error_theorem=owner.proof,independent_full_density_reference=independent,
        current_original_normalized_pressure_point_jet_service_installed=True,
        actual_native_normalized_pressure_queries=queries,
        maximum_coefficient_and_roundoff_absolute_error_upper=maximum_error,
        approximate_point_quadrature_distinct_from_directed_enclosure=True,
        exact_odd_symmetry_zero_uses_explicit_zero_budget=True,
        finite_axial_end_and_all_eleven_late_stages_not_discarded=True,
        physical_Pstar_squared_error_amplification_not_discarded=True,
        invalid_domain_or_quadrature_setting_rejections=rejected,**dict.fromkeys(flags,False),
        input_hashes={**owner.hashes,point.NAME:point.sha(point.NAME),Path(__file__).name:point.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original normalized pressure point approximations through ordinary Z order2 with directed coefficient/roundoff and source-linked all-late-stage errors. Finite full-density fixture tests the general bound only. Native physical amplification, conditioned p1/p2/phase, full numerical oracle, controls and global field remain open.')
    (point.HERE/point.RECEIPT).write_text(json.dumps(point.source.inertial.profiles.loop.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original pressure point jets: native budgets and full fourteen-density reference PASS',flush=True)
    return result


if __name__=='__main__':run()
