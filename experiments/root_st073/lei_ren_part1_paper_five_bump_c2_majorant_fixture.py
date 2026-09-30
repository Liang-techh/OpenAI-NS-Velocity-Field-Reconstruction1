"""Resolved C2 analytic family versus a conditional interval tail bound."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse import iterate_inverse
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_five_bump_c2_majorant import compute_c2_majorant,amplitude_inverse_square_c2_bound


def run():
    with mp.workdps(120):
        m=FiveBumpMomentMap(precision=120,order=96);a=mp.mpf('.8');p=mp.mpf(2)
        seed=compute_c2_majorant(m,0,p,axial_radius=a)
        threshold=seed['raw']['condition_threshold'];e=threshold/100
        bound=compute_c2_majorant(m,e,p,axial_radius=a,inverse_steps=10)
        # d_i=c_i(1+Z^2): exact sum of derivative suprema on [-a,a].
        coeff=[e*i/(15*(2+a*a+2*a)) for i in range(1,6)]
        errors=[]
        for z in (-a,mp.mpf(0),mp.mpf('.3'),a):
            zd=AxialSecondJet(z,1,0,pressure_order=0,width_order=0)
            amplitude=p*mp.exp(mp.mpf('-.6'))/(1+zd*zd)
            defects=[c*(1+zd*zd) for c in coeff]
            finite=iterate_inverse(m,defects,amplitude,steps=10)
            reference=iterate_inverse(m,defects,amplitude,steps=40)
            error=mp.fsum(abs((x-y).value.evaluate())+abs((x-y).tangent.evaluate())+
                         abs((x-y).second.evaluate())/2 for x,y in zip(finite['h'],reference['h']))
            assert error<=bound['raw']['incremental_inverse_tail_bound']
            errors.append(mp.nstr(error,50))
        boundary=compute_c2_majorant(m,threshold,p,axial_radius=a)
        assert boundary['raw']['contraction_condition_passed']
        assert abs(boundary['raw']['lipschitz_bound']-mp.mpf('.5'))<mp.mpf('1e-110')
        assert not compute_c2_majorant(m,2*threshold,p,axial_radius=a)['raw']['contraction_condition_passed']
        # Derivative envelope endpoints agree with the elementary formula.
        expected=mp.exp(mp.mpf('1.2'))/p**2*((1+a*a)**2+4*a*(1+a*a)+2+6*a*a)
        assert abs(amplitude_inverse_square_c2_bound(p,a)-expected)<mp.mpf('1e-110')
        report=dict(resolved_family_interval_radius='.8',family_defect_norm=mp.nstr(e,50),
            finite_vs_40_update_pointwise_C2_errors=errors,
            conditional_C2_tail=mp.nstr(bound['raw']['incremental_inverse_tail_bound'],50),
            boundary_passed=True,above_threshold_rejected=True,
            independent_defect_interval_norm='analytic polynomial supremum sums',
            metadata=bound['metadata'],actual_source_norm_verified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('C2 family tail bound passed',report['conditional_C2_tail'],flush=True)


if __name__=='__main__':run()
