"""Independent differentiation of positive pressure factors and error budgets."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_preheat_interval_bounds import preheat_interval_bounds,q_power_derivative_bounds,positive_stage_mass_upper,positive_stage_pressure_bounds


def run():
    with mp.workdps(100):
        a=mp.mpf('.8');tiny=mp.exp(-1000)
        receipt=dict(Pstar_squared=mp.mpf(4),components={
            'prefix':dict(kind='q_power',beta=mp.mpf(2),value_at_Z0=mp.mpf('2.5')),
            'flatten':dict(kind='q_power_atoms',atoms=[dict(beta=mp.mpf('.7'),atom_at_Z0=mp.mpf('.2')),dict(beta=mp.mpf('.001'),atom_at_Z0=tiny)]),
            'exterior':dict(kind='exact_exterior_power_atom',Z_independent=True,value_at_Z0=mp.mpf('.1'))})
        bounds=preheat_interval_bounds(receipt,axial_radius=a)
        errors={label:row['mass']/100 for label,row in bounds['atoms'].items()}
        augmented=preheat_interval_bounds(receipt,axial_radius=a,mass_error_bounds=errors)
        checks=0
        for label,row in bounds['atoms'].items():
            for z in (-a,0,mp.mpf('.3'),a):
                for n in range(4):
                    derivative=mp.diff(lambda zz:row['mass']*(1+zz*zz)**(-row['beta']),z,n)
                    assert abs(derivative)<=row['derivative_bounds'][n]*(1+mp.mpf('1e-95'))
                    checks+=1
        assert bounds['atoms']['flatten:1']['mass']==tiny
        assert augmented['normalized_mass_error_derivative_bounds'][0]>0
        # g(y)=.2-.4y meets both slope bounds, on [0,3].
        mass=positive_stage_mass_upper(3,'.2','-1')
        direct=mp.quad(lambda y:mp.exp(2*(mp.mpf('.2')-mp.mpf('.4')*y))/2,[0,3])
        assert direct<=mass['mass_upper']
        variable=positive_stage_pressure_bounds(3,'.2','-1',axial_radius=a)
        for n in range(4):
            # beta(y)=2y/3 genuinely varies across the radial stage.
            val=mp.quad(lambda y:mp.exp(2*(mp.mpf('.2')-mp.mpf('.4')*y))/2*
                        mp.diff(lambda z:(1+z*z)**(-2*y/3),a,n),[0,3])
            assert abs(val)<=variable['derivative_bounds'][n]
        try:preheat_interval_bounds(receipt,axial_radius=a,mass_error_bounds={})
        except KeyError:pass
        else:raise AssertionError('Missing error budget was accepted')
        report=dict(independent_derivative_checks=checks,compact_radius='.8',tiny_atom_preserved=True,
            mass_error_budgets_separate=True,positive_stage_mass_envelope_passed=True,
            variable_beta_integral_derivatives_checked=True,
            normalized_C2_bound=mp.nstr(bounds['normalized_C2_bound'],60),
            actual_quadrature_enclosed=False,actual_source_norm_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('preheat analytic envelope checks passed',checks,flush=True)


if __name__=='__main__':run()
