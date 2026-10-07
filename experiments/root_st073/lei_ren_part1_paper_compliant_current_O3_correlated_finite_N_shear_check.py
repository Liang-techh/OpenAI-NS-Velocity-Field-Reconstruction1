"""Independent joint-shear algebra, taper queries and scope acceptance."""
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_O3_correlated_finite_N_shear as source


def independent_correlated_samples():
    count=0
    # Arbitrary cutoff logarithmic derivatives test the feature discarded
    # by separate error caps. This is supplementary to the exact theorem.
    with mp.workdps(100):
        for muv in ('0.001','1e-20'):
            mu=mp.mpf(muv)
            for N in (1,37):
                for chi in (mp.mpf(1),mp.mpf('1e-30')):
                    for D in (mp.mpf('-1e60'),mp.mpf(-1),mp.mpf(0),mp.mpf('1e60')):
                        L=-mp.mpf('.5')-mu
                        for phase in (mp.mpf(0),mp.mpf('.137'),mp.mpf('.25'),mp.mpf('.499')):
                            sn,cs=mp.sin(2*mp.pi*phase),mp.cos(2*mp.pi*phase)
                            A=-mu*chi**2*mp.sin(4*mp.pi*phase)/(8*mp.pi)
                            # Evaluate the original two source expressions,
                            # rather than the completed-square expression.
                            dy=chi*D
                            slowA=-mu*chi*dy*mp.sin(4*mp.pi*phase)/(4*mp.pi)
                            slowB=-mp.sqrt(mu)*(chi*L+dy)*mp.cos(2*mp.pi*phase)/(2*mp.pi)
                            theta=mu*chi**2*mp.cos(4*mp.pi*phase)-2*slowA/N
                            axial=mp.exp(-A/N)*(2*mp.sqrt(mu)*chi*sn+2*slowB/N)
                            joint=theta+axial**2/(2+3*mu)
                            if joint < mu*chi**2/4:
                                raise ArithmeticError('Independent correlated primitive margin loses quarter reserve')
                            count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed joint primitive shear source: '+name)
    field=source.CurrentCorrelatedFiniteNShear(require_checked=False)
    if source.parameters.encoded(field.theorem)!=data['exact_correlated_finite_N_primitive_shear_theorem']:
        raise ValueError('Joint primitive source theorem differs')
    specifications=(('left_flat_edge','-2',1),('left_taper','-1.99',1),('negative_plateau','-1',1),
        ('old_seam','0',1),('right_taper','.49',1),('right_flat_edge','.5',1),('O3_terminal','1',1),
        ('high_frequency_taper','.337',10**12))
    for name,t,N in specifications:
        fresh=field.query(t,N)
        if source.parameters.encoded(fresh)!=data['examples'][name]:raise ValueError('Actual joint taper query differs')
    examples={name:field.query(t,N) for name,t,N in specifications}
    for name in ('left_flat_edge','right_flat_edge'):
        if not examples[name]['exact_flat_cutoff_value_and_first_derivative_zero']:
            raise ArithmeticError('Actual flat source endpoint changed')
        for label in ('actual_phase_theta_shear_error_absolute_cap','actual_phase_axial_shear_error_absolute_cap'):
            if examples[name][label]._mpi_!=field.ctx.mpf(0)._mpi_:
                raise ArithmeticError('Actual phase error cap does not vanish at the flat source endpoint')
    if examples['left_flat_edge']['direct_joint_margin_enclosure']._mpi_!=field.ctx.mpf(0)._mpi_:
        raise ArithmeticError('Original O2 exact zero endpoint margin was reset')
    if examples['right_flat_edge']['direct_joint_margin_enclosure']._mpi_!=field.data['mu']._mpi_:
        raise ArithmeticError('O3 exact right-edge margin differs from actual mu')
    if not source.parameters.numeric.transport.endpoints(examples['old_seam']['directed_scalar_margin_lower_bound'])[0]>0:
        raise ArithmeticError('Corrected old seam lacks a directed positive primitive floor')
    if any(data.get(key) is not False for key in source.OPEN):
        raise ValueError('Primitive inequality promoted completed repair/tensor/recursion')
    invalid=((0,0),(0,True),(-3,1),(2,1),(mp.inf,1))
    for t,N in invalid:
        try:field.query(t,N)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid actual primitive shear query accepted')
    tested=independent_correlated_samples()
    hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)}
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_correlated_primitive_source_and_square_identities=len(field.theorem['identities']),
        independent_original_shear_correlated_samples=tested,actual_taper_queries_recomputed=8,
        old_seam_positive_left_endpoint_zero_right_endpoint_actual_mu=True,
        actual_mu_interval_retained_without_dividing_by_taper=True,
        joint_primitive_shear_margin_certified_on_actual_O2_O3_profiles=True,
        all_real_phases_and_finite_integer_N_ge1_covered_by_exact_inequality=True,
        invalid_queries_rejected=len(invalid),**{key:False for key in source.OPEN},input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Joint primitive shear PASS:',len(field.theorem['identities']),'identities;',tested,'independent correlated samples;8 actual taper queries',flush=True)
    return result


if __name__=='__main__':run()
