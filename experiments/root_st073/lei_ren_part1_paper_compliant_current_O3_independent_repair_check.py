"""Independent acceptance of new source-bound implicit five-bump map scope."""
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_O3_independent_repair import (
    CurrentO3IndependentRepair,HERE,PREFIX,NAME,RECEIPT,GATES,OPEN,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import (
    exact_repair_theorem,scaled_quadratic,contraction_certificate)

@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3IndependentRepair(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    if encode(pack(field.manifest()))!=raw or any(raw[k] for k in GATES+OPEN):
        raise ValueError('Actual independent repair producer/scope differs')
    if encode(pack(exact_repair_theorem()))!=raw['exact_source_and_new_five_bump_map_theorem']:
        raise ValueError('Current source or exact physical/scaled nonlinear repair map differs')
    c=field.ctx;ci=field.weights['centers'];ell=field.weights['radius'];matrix=field.matrix
    if endpoints(ci[0]-ell)[0]<=0 or endpoints(ci[2]+ell)[1]>=1:
        raise ValueError('Independent bumps must lie strictly inside q[1,2]')
    for a,b in zip(ci[:-1],ci[1:]):
        if endpoints(b-a-2*ell)[0]<=0:raise ValueError('Distinct centers must have disjoint supports')
    B=matrix['linear'];inverse=matrix['inverse'];inverse_checks=0
    for i in range(5):
        for j in range(5):
            product=sum((B[i][k]*inverse[k][j] for k in range(5)),c.mpf(0))
            lo,hi=endpoints(product)
            if not lo<=int(i==j)<=hi:raise ValueError('New divided-difference inverse orientation fails')
            inverse_checks+=1
    cert=field.certificate;r=endpoints(cert['scaled_controls_ball_radius'])[1]
    if not cert['unique_exact_implicit_controls_exist'] or cert['common_cone_frequency_admitted']:
        raise ValueError('Independent repair-only contraction cannot admit complete cones')
    if endpoints(cert['strict_contraction_upper'])[1]>=1 or endpoints(cert['strict_image_radius_upper'])[1]>=r:
        raise ValueError('Strict scaled ball contraction/inclusion required')
    for control in field.controls:
        if max(abs(x) for x in endpoints(control))>=r:raise ValueError('Implicit controls do not lie strictly inside the ball')
    pure_axial=scaled_quadratic(c,field.mu,field.N,matrix,[c.mpf(1),c.mpf(0),c.mpf(0),c.mpf(0),c.mpf(0)])
    pure_swirl=scaled_quadratic(c,field.mu,field.N,matrix,[c.mpf(0),c.mpf(0),c.mpf(1),c.mpf(0),c.mpf(0)])
    overlap=scaled_quadratic(c,field.mu,field.N,matrix,[c.mpf(1),c.mpf(0),c.mpf(1),c.mpf(0),c.mpf(0)])
    if endpoints(pure_axial[3])[0]<=0 or endpoints(pure_axial[4])!=(0,0):
        raise ValueError('Actual axial-square kinetic term omitted')
    if endpoints(pure_swirl[3])[1]>=0 or endpoints(pure_swirl[4])[0]<=0:
        raise ValueError('Actual swirl-square energy/pressure terms omitted')
    if endpoints(overlap[1])[0]<=0:raise ValueError('Same-center angular/axial cross term omitted')
    try:contraction_certificate(c,field.mu,cert['repair_only_sufficient_integer_N']-1,matrix)
    except ValueError:pass
    else:raise ValueError('Below-bound frequency accepted as repair sufficient')
    point=field.coefficients()
    if not point['exact_controls_are_Z_independent'] or any(point[k] for k in OPEN):
        raise ValueError('Implicit constant controls exceed constructed field scope')
    result=dict(actual_five_defect_family_sha256=field.histories.family,
      implicit_source_sha256=field.histories.source,datum_enclosure_sha256=field.histories.datum_sha,
      parent_modified_source_definition_sha256=field.histories.modified_source_definition_sha256,
      repair_definition_sha256=field.repair_definition_sha256,finite_integer_N=field.N,
      new_exact_source_and_physical_scaled_map_identities=len(field.theorem['identities']),
      new_block_inverse_identity_enclosures=inverse_checks,
      fixed_independent_supports_before_common_N=True,
      actual_divided_axial_defect_formed_before_enclosure=True,
      overlapping_axial_swirl_cross_and_both_square_terms_retained=True,
      unique_constant_controls_define_exact_source_repair_not_midpoint_fields=True,
      original_incoming_moments_and_same_axis_datum_retained=True,
      strict_actual_N_repair_only_contraction_certified=True,
      repaired_partial_history_and_tensor_source_installation_remains_open=True,
      common_modified_cone_frequency_remains_open=True,
      input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      scope=raw['scope'],**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Independent source-correlated five-bump map, inverse and unique implicit constant controls PASS; installed repair/cone open',flush=True)
    return result

if __name__=='__main__':run()
