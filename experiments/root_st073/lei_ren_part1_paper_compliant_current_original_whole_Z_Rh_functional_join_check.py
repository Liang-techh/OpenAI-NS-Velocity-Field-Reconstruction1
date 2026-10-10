"""Check the current functional join and genuine correction continuation."""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join as current


def run(field=None):
    began=time.monotonic();raw=json.loads((current.HERE/current.NAME).read_bytes())
    assert raw['all_passed'] and raw['candidate_current_Rh_functional_join_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.WholeZRhFunctionalJoin(require_checked=False)
        assert not field.acceptance_loaded and field.identity==raw['source_family']
        assert current.serialized(current.source_bindings())==raw['actual_source_bindings']
        assert current.exact_join_proof()==raw['exact_current_function_proof']
        # Independently derive the canonical five histories from physical
        # primitives. No interval subtraction is used for this equality.
        y,z=s.symbols('y Z',real=True);P=s.Symbol('Pstar',positive=True)
        E=s.exp((y-6)/10)/(1+z*z);V=4*z;x=s.exp(y)
        physical=dict(Mz=x*V,Mtheta=s.Rational(5,8)*s.exp(3*y/2)*E,
            Mtheta_z=s.Rational(5,8)*s.exp(3*y/2)*E*V,
            Mztheta=x*(V*V/P**2-s.Rational(5,12)*E*E),Mp=s.Rational(5,2)*E*E)
        histories=dict(m=physical['Mz']/x,h=physical['Mtheta']/s.exp(3*y/2),
            k=physical['Mtheta_z']/s.exp(3*y/2),e=physical['Mztheta']/x,p=physical['Mp'])
        expected=dict(m=V,h=s.Rational(5,8)*E,k=s.Rational(5,8)*E*V,
            e=V*V/P**2-s.Rational(5,12)*E*E,p=s.Rational(5,2)*E*E)
        assert all(s.simplify(histories[k]-expected[k])==0 for k in expected)
        # Duhamel inlet and differentiated equations, including rate zero.
        lam,Y,Z=s.symbols('lambda y Z',real=True);H0=s.Function('H0')(Z)
        D=s.Function('D');q=s.Symbol('s',real=True)
        H=s.exp(-lam*Y)*H0+s.Integral(s.exp(-lam*(Y-q))*D(q,Z),(q,0,Y))
        assert s.simplify(H.subs(Y,0)-H0)==0
        assert s.simplify(s.diff(H,Y)+lam*H-D(Y,Z))==0
        assert s.simplify(s.diff(H,Z).subs(Y,0)-s.diff(H0,Z))==0
        assert s.exp(-lam*Y).subs(lam,0)==1
        N,logRm,logRa,hb,sc,offset=s.symbols('N logRm logRa hb sc offset')
        phase_patch=N*(logRm+Y-logRa-hb*sc/2)
        phase_reference=N*(logRm+6+offset-logRa-hb*sc/2)
        assert s.expand(phase_patch-phase_reference.subs(offset,Y-6))==0
        assert phase_patch.subs(Y,1)==phase_reference.subs(offset,-5)
        counts=dict(current_whole_Z_cells=0,shared_implicit_owner_and_P0=0,
            actual_C1_inlet_handles=0,actual_reference_first_unit_histories=0,
            independent_physical_history_identities=5,exact_Duhamel_C1_identities=4,
            domain_rejections=0)
        for ends,stored in zip(current.CELLS,raw['source_cells']):
            join=field.seam(ends);assert current.serialized(join)==stored['join']
            op=field.outer.owner.owner(ends)
            assert field.leading.owner(ends).op is op
            assert join['exact_common_P0_axial5'] is op.P0
            assert any(not row.zero for row in op.P0)
            assert join['current_unique_implicit_owner_shared_by_both_branches']
            assert join['phase_same_at_Rh_and_not_restarted']
            assert join['current_global_phase']=='frac(N*(logRm+y-logRa-hb*s_c/2))'
            assert join['reference_global_phase']=='frac(N*(logRm+6+offset-logRa-hb*s_c/2))'
            assert join['leading_identity_does_not_zero_correction_incoming']
            assert not join['overlap_only_consistency_not_functional_identity_proof']
            assert current.outer.ep(op.c.exp(1))[0]>mp.mpf(71)/40
            p=field.leading.query(ends,'Rh')
            assert p['actual_partial_primitive_source_memory']['exact_terminal_identity_of_same_leading_map']
            assert all(current.outer.ep(v)==(0,0) for row in p['actual_gamma_ordinary_x_derivatives'] for v in row)
            inlet=field.continuation(ends);next_unit=field.continuation(ends,(-4,1))
            assert current.serialized(inlet)==stored['reference_inlet']
            assert current.serialized(next_unit)==stored['reference_first_unit']
            original=field.outer.normalized_Rh_incoming(ends)
            assert join['incoming_defining_source']['source_report']==current.outer.previous.NAME
            assert join['incoming_defining_source']['source_receipt']==current.outer.previous.RECEIPT
            assert join['incoming_defining_source']['method']=='WholeZAllNLongPatchFunctions.transport'
            assert join['incoming_defining_source']['original_core_zero_inlet_and_bridge_R110_predecessor_chain_retained']
            epsilon=current.outer.core.epsilon_cover(field.c,field.N0)
            # This is an all-N parameter cover, not epsilon fixed to 1/N0.
            assert current.outer.ep(epsilon)[0]==0
            assert current.outer.ep(epsilon)[1]==current.outer.ep(field.c.mpf(1)/field.N0)[1]
            assert inlet['epsilon']._mpi_==epsilon._mpi_ and next_unit['epsilon']._mpi_==epsilon._mpi_
            density=field.outer.query(ends,'Rh_reference',(-5,1),(-4,1))['N_scaled_density_C0_Z']
            for name,rate in current.RATES.items():
                graph=inlet['defining_correction_function_graph'][name]
                advanced=next_unit['defining_correction_function_graph'][name]
                assert graph['inlet_handle']==advanced['inlet_handle']
                assert graph['at_Rh_exactly_same_inlet_handle']
                assert not advanced['at_Rh_exactly_same_inlet_handle']
                assert graph['exact_log_distance']==[0,1] and advanced['exact_log_distance']==[1,1]
                assert graph['own_rate']==str(rate)
                assert graph['exact_Z_derivative_uses_same_inlet_Z_and_density_Z']
                assert graph['source_density']['method']=='normalized_drivers'
                assert graph['source_density']['actual_phase']=='frac(N*(logRh+s-logRa-hb*s_c/2))'
                assert graph['source_density']['module']==Path(current.outer.core.__file__).name
                assert current.serialized(current.outer.record(original[name]))==current.serialized(inlet['N_scaled_correction_range_C0_Z'][name])
                assert current.outer.previous.equivalent_rows([inlet['own_rate_incoming_memory'][name]],[op.flow.scalar(1)])
                lam0=field.c.mpf(rate.numerator)/rate.denominator
                want=op.flow.factor((0,0,0,0,0),-lam0)
                assert current.outer.previous.equivalent_rows([next_unit['own_rate_incoming_memory'][name]],[want])
                assert any(not row.zero for row in (original[name].value,original[name].Z))
                # Recompute both complete rows from the actual background,
                # actual local driver and one epsilon factor. This would
                # reject omission or double application of the correction.
                weight=field.outer.weights(ends,'Rh_reference',(-5,1),(-4,1),rate)
                advanced_correction=current.add(current.scale(original[name],weight['cell_decay']),
                    current.scale(density[name],weight['original_mass']))
                for offset0,packet,correction in (((-5,1),inlet,original[name]),((-4,1),next_unit,advanced_correction)):
                    background=field.outer.leading_packet(ends,'Rh_reference',offset0)['original_generic_source']['common_own_five_histories_axial5'][name]
                    complete=current.add(current.Pair(*background[:2]),current.scale(correction,epsilon))
                    assert current.serialized(current.outer.record(complete))==current.serialized(packet['complete_background_plus_epsilon_correction_range_C0_Z'][name])
                for row in inlet['N_scaled_local_integral_range_C0_Z'][name]:
                    assert row['exact_zero']
                counts['actual_C1_inlet_handles']+=2
                counts['actual_reference_first_unit_histories']+=2
            assert current.RATES['p']==0
            assert current.outer.previous.equivalent_rows([next_unit['own_rate_incoming_memory']['p']],[op.flow.scalar(1)])
            assert next_unit['exact_common_P0_axial5'] is op.P0
            assert next_unit['original_P0_independent_and_included_once']
            assert next_unit['ranges_are_enclosures_not_defining_function_values']
            assert not next_unit['full_finite_N_mixed4_or_global_field_installed']
            assert not any(next_unit[k] for k in current.OPEN)
            counts['current_whole_Z_cells']+=1;counts['shared_implicit_owner_and_P0']+=1
        for call in (lambda:field.seam(('0','0')),lambda:field.continuation(current.CELLS[0],(-6,1)),
                     lambda:field.continuation(current.CELLS[0],(1,1))):
            try:call()
            except ValueError:counts['domain_rejections']+=1
            else:raise AssertionError('Unadmitted cell/offset must reject')
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            fixed_leading_input_N0=field.N0,current_full_weight_map_and_source_AST_checked=True,
            canonical_leading_functions_and_implied_mixed4_join_checked=True,
            genuine_nonzero_correction_inlet_retained_and_C1_Duhamel_source_installed=True,
            original_independent_P0_and_actual_radius_phase_retained=True,
            interval_overlap_not_used_as_function_identity=True,replay_counts=counts,
            full_finite_N_mixed4_or_global_field_installed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.serialized(result),indent=2)+'\n',encoding='utf8')
    print('Current whole-Z functional Rh join/correction continuation checks passed',flush=True)
    return result


if __name__=='__main__':run()
