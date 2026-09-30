"""Full degree-six terminal moment tail of the nominal degree-three response.

Evaluate coefficients before summing, avoiding cancellation of large public
physical moments. This receipt does not bound the infinite response or source error.
"""
import json
import pickle
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_moment_reference_background import ReferenceDefectBackground
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_response import build_response, SparseDefectPolynomial
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'reference_defect_Z03.pkl'
    with cache.open('rb') as stream: snapshot=pickle.load(stream)
    if snapshot['version']!=1 or snapshot['Z']!='.3':raise ValueError('Unsupported source snapshot')
    with mp.workdps(snapshot['precision']):
        data=snapshot['data']; amplitude=data['Am'].evaluate()
        defects=[data['d'][i].evaluate() for i in range(1,6)]
        moment_map=FiveBumpMomentMap(precision=snapshot['precision'],order=96)
        response=build_response(moment_map,amplitude,degree=3)
        h=tuple(SparseDefectPolynomial({key:rows[i] for key,rows in response.coefficients.items()},
                                     max_degree=6) for i in range(5))
        actual=moment_map.apply(h,amplitude)
        residual=[]
        for i in range(5):
            unit=tuple(int(j==i) for j in range(5))
            residual.append(actual[i]+SparseDefectPolynomial({unit:mp.mpf(1)},max_degree=6))
        def evaluated(poly,degree):
            values=[]
            for key,coefficient in poly.terms.items():
                if sum(key)!=degree:continue
                term=coefficient
                for d,power in zip(defects,key):term*=d**power
                values.append((key,term))
            return values
        evaluated_h=[mp.fsum(rows[i]*mp.fprod(d**p for d,p in zip(defects,key))
                     for key,rows in response.coefficients.items()) for i in range(5)]
        direct=moment_map.apply(evaluated_h,amplitude)
        rows=[]; replay_errors=[]
        for i,poly in enumerate(residual):
            degrees={}
            tail=[]
            for degree in range(1,7):
                values=evaluated(poly,degree)
                total=mp.fsum(v for _,v in values)
                absolute=mp.fsum(abs(v) for _,v in values)
                degrees[str(degree)]=dict(term_count=len(values),sum=signed_log(total,70),
                                         absolute_term_sum=signed_log(absolute,70))
                if degree>=4:tail.extend(values)
            tail_sum=mp.fsum(v for _,v in tail);tail_l1=mp.fsum(abs(v) for _,v in tail)
            largest=sorted(tail,key=lambda kv:abs(kv[1]),reverse=True)[:5]
            full_terms=[v for degree in range(1,7) for _,v in evaluated(poly,degree)]
            replay_error=abs((direct[i]+defects[i])-mp.fsum(full_terms))
            replay_errors.append(replay_error)
            rows.append(dict(row=i+1,input_defect=signed_log(defects[i],70),degrees=degrees,
                aggregate_replay_absolute_error=signed_log(replay_error,70),
                unresolved_degree_4_to_6_tail=signed_log(tail_sum,70),
                absolute_term_sum=signed_log(tail_l1,70),
                largest_tail_terms=[dict(monomial=key,value=signed_log(v,70)) for key,v in largest]))
        assert max(replay_errors)<mp.mpf('1e-400'), 'Full quadratic replay mismatch'
        report=dict(Z='.3',nominal_evaluation='retained actual P9/W2 inputs at pressure=width=1',
            precision=snapshot['precision'],velocity_defect_degree=3,moment_defect_degree=6,
            full_finite_quadratic_products_retained=True,rows=rows,
            composition_mode='scalar response after nominal input evaluation; no additional P9/W2 product truncation',
            equivalent_to_nominal_truncated_field_jet_certified=False,
            aggregate_replay_max_absolute_error=signed_log(max(replay_errors),70),
            replay_scope='same quadrature map, scalar aggregate versus degree-six termwise evaluation; not independent quadrature',
            degrees_1_to_3='finite precision coefficient cancellation only',
            degrees_4_to_6='unrepaired finite-response tail, not discarded',
            infinite_response_tail_enclosed=False,source_error_enclosed=False,
            quadrature_enclosed=False,functional_closure=False,
            source_stage_products_expanded=False,first_Z_evaluated=False,
            temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        for row in rows:print('row',row['row'],'tail',row['unresolved_degree_4_to_6_tail'],flush=True)
        return report


if __name__=='__main__':run()
