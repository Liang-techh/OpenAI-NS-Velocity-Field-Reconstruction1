"""Independent normalized density references and actual all-N bridge replay."""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

encoded=lambda value:current.current.encode(current.current.serialized(value))


def independent_density_fixture():
    c=MPIntervalContext();c.dps=150
    p=mp.mp.clone();p.dps=190
    f=MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    def sources(z):
        return [p.mpf('.9')+p.mpf('.2')*z+p.mpf('.03')*z*z,
            -p.mpf('.4')+p.mpf('.11')*z-p.mpf('.02')*z*z,
            p.mpf('.3')-p.mpf('.14')*z+p.mpf('.017')*z*z,
            -p.mpf('.25')+p.mpf('.13')*z-p.mpf('.019')*z*z]
    def numerator(z,epsilon):
        E,V,A,B=sources(z);newE=E*p.exp(epsilon*A);newV=V+epsilon*B
        return dict(m=newV-V,h=newE-E,k=newV*newE-V*E,
            e=newV*newV-V*V-(newE*newE-E*E)/2,p=(newE*newE-E*E)/2)
    def reference(z,epsilon,key):
        if not epsilon:return p.diff(lambda eps:numerator(z,eps)[key],0)
        return numerator(z,epsilon)[key]/epsilon
    count=0
    for ztext in ('-.3','0','.41'):
        z=p.mpf(ztext)
        pair=lambda i:current.Pair(f.scalar(c.mpf(p.nstr(sources(z)[i],185))),
            f.scalar(c.mpf(p.nstr(p.diff(lambda value:sources(value)[i],z),185))))
        E,V,A,B=[pair(i) for i in range(4)]
        primitives=dict(A=A.value,A_Z=A.Z,B_over_Pstar=B.value,B_Z_over_Pstar=B.Z)
        for denominator in (None,256,512):
            eps=p.mpf(0) if denominator is None else p.mpf(1)/denominator
            for mode in ('point','uniform'):
                argument=c.mpf(0) if denominator is None else c.mpf(1)/denominator
                got=current.normalized_drivers(E,V,primitives,256,argument if mode=='point' else None)
                for key,q in got['drivers'].items():
                    for order,row in enumerate((q.value,q.Z)):
                        want=p.diff(lambda value:reference(value,eps,key),z,order)
                        lo,hi=current.ep(row.finite_interval(max_log=1000))
                        assert lo<=want<=hi,(ztext,denominator,mode,key,order)
                        count+=1
    for bad in (c.mpf('-.01'),c.mpf('.1')):
        try:current.normalized_drivers(E,V,primitives,256,bad)
        except ValueError:pass
        else:raise AssertionError('Unadmitted epsilon accepted')
    # The extension retains the exact zero source and does not divide by eps.
    zero=f.scalar(0)
    got=current.normalized_drivers(current.Pair(zero,zero),V,
        dict(A=zero,A_Z=zero,B_over_Pstar=zero,B_Z_over_Pstar=zero),256,None)
    assert all(q.value.zero and q.Z.zero for q in got['drivers'].values())
    return dict(passed=True,independent_velocity_moment_difference_C0_Z_comparisons=count,
        epsilon_zero_limit_and_positive_samples_used=True,point_and_full_epsilon_range_compared=True,
        every_quadratic_cross_term_and_nonzero_Z_source_used=True,
        synthetic_reference_not_actual_source_claim=True)


def phase_cover(alternatives):
    pieces=sorted(set((Fraction(*row['exact_phase_left']),Fraction(*row['exact_phase_right']))
        for row in alternatives))
    assert pieces and pieces[0][0]==0 and pieces[-1][1]==1
    assert all(a<b for a,b in pieces)
    assert all(before[1]==after[0] for before,after in zip(pieces,pieces[1:]))
    return len(pieces)


def independent_union_fixture():
    c=MPIntervalContext();c.dps=150
    p=mp.mp.clone();p.dps=190
    f=MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    prior=current.backend.signed.density.local.prior
    values=[prior.ScaledEnclosure(prior.FormalScale(f.logs,offset=c.mpf(('-100000','100000'))),c.mpf(('-2','3')),f.ledger),
        prior.ScaledEnclosure(prior.FormalScale(f.logs,offset=c.mpf(('99999','100001'))),c.mpf(('-.7','-.2')),f.ledger),f.scalar(0)]
    union=current.uniform_source_union(values)
    lo,hi=current.ep(union.coefficient);count=0
    upper=current.ep(union.scale.offset)[0]
    for value in values[:-1]:
        logs=current.ep(value.scale.offset);coefficients=current.ep(value.coefficient)
        for log in (logs[0],sum(logs)/2,logs[1]):
            for coefficient in (coefficients[0],sum(coefficients)/2,coefficients[1]):
                want=coefficient*p.exp(log-upper)
                assert lo<=want<=hi
                count+=1
    assert lo<=0<=hi and current.uniform_source_union([f.scalar(0)]).zero
    assert union.scale.bases is f.logs and union.ledger is f.ledger
    return dict(passed=True,independent_signed_wide_log_union_endpoint_comparisons=count,
        independent_wide_offsets_and_tiny_directed_tails_retained=True,
        same_basis_ledger_and_exact_zero_retained=True)


def forbidden(*args,**kwargs):raise AssertionError('Unchanged local source must not be recomputed')


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed prerequisite: '+name)
    with mp.workdps(540):
        fixture=independent_density_fixture()
        union_fixture=independent_union_fixture()
        owner=current.WholeZAllNBridgeSourceFunctions()
        assert saved['source_family']==owner.identity and saved['fixed_leading_input_N0']==owner.N0
        assert saved['all_N_source_density_binding']==owner.binding
        cells=phase_pieces=0
        with patch.object(current.current.WholeZSharpBridgeFunctions,'query',forbidden), \
                patch.object(current.current.WholeZSharpBridgeFunctions,'local_integral',forbidden):
            for transport in saved['actual_four_Z_all_N_bridge_transports']:
                ends=tuple(transport['exact_Z_cell']);f=owner.owner.owner(ends).flow
                decode=lambda value:current.relative.restore_row(f,value)
                pairs=lambda rows:{name:current.Pair(*[decode(v) for v in pair]) for name,pair in rows.items()}
                assert all(current.relative.exact_row(a,decode(b)) for a,b in
                    zip(owner.owner.owner(ends).P0,transport['exact_common_P0_axial5']))
                incoming={name:current.Pair(f.scalar(0),f.scalar(0)) for name in current.RATES}
                assert all(row['exact_zero'] for pair in transport['actual_N_scaled_zero_inlet_C0_Z'].values() for row in pair)
                for window in transport['actual_all_N_bridge_windows']:
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_incoming_C0_Z']
                    local={name:current.Pair(f.scalar(0),f.scalar(0)) for name in current.RATES}
                    for cell in window['actual_all_N_source_cells']:
                        phase_pieces+=phase_cover(cell['full_phase_inverse_function_alternatives'])
                        assert cell['physical_Jacobian_applied_once'] and cell['all_N_phase_cover_not_fixed_N_phase_reuse']
                        assert current.ep(current.bridge.source_module.read(owner.c,cell['epsilon']))==current.ep(current.epsilon_cover(owner.c,owner.N0))
                        density=pairs(cell['N_scaled_density_C0_Z']);contributions={}
                        for name,rate in current.RATES.items():
                            mass,decay,suffix=owner.owner.weights(ends,window['actual_chart'],tuple(cell['exact_left']),tuple(cell['exact_right']),rate)
                            contribution=current.scale(density[name],mass*suffix)
                            contributions[name]=contribution;local[name]=current.add(local[name],contribution)
                            assert encoded(dict(original_mass=mass,cell_decay=decay,suffix_decay=suffix))==cell['original_own_rate_weights'][name]
                        assert encoded({k:current.record(v) for k,v in contributions.items()})==cell['actual_N_scaled_contribution_C0_Z']
                        cells+=1
                    assert encoded({k:current.record(v) for k,v in local.items()})==window['actual_N_scaled_local_C0_Z']
                    memory={name:decode(v) for name,v in window['original_incoming_own_rate_memory'].items()}
                    assert current.ep(memory['p'].finite_interval())==(1,1)
                    incoming={name:current.add(current.scale(incoming[name],memory[name]),local[name]) for name in current.RATES}
                    assert encoded({k:current.record(v) for k,v in incoming.items()})==window['actual_N_scaled_outgoing_C0_Z']
                assert encoded({k:current.record(v) for k,v in incoming.items()})==transport['actual_uniform_N_scaled_R100_correction_C0_Z']
        assert cells==40
        fresh=owner.query(current.CELLS[0],'first_micro',(3,4),(1,1))
        stored=saved['actual_four_Z_all_N_bridge_transports'][0]['actual_all_N_bridge_windows'][0]['actual_all_N_source_cells'][3]
        assert encoded({k:current.record(v) for k,v in fresh['N_scaled_density_C0_Z'].items()})==stored['N_scaled_density_C0_Z']
        assert phase_cover(fresh['source']['full_phase_inverse_function_alternatives'])==phase_cover(stored['full_phase_inverse_function_alternatives'])
        velocity=fresh['source']['original_E_V_C0_Z'];comparisons=0
        for alternative in fresh['source']['full_phase_inverse_function_alternatives']:
            primitive=alternative['primitive_C0_Z_phi']
            direct=current.bridge.SIGNED_DENSITY(velocity['E'].value,velocity['E'].Z,
                velocity['V'].value,velocity['V'].Z,primitive,owner.N0)
            normalized=current.normalized_drivers(velocity['E'],velocity['V'],primitive,owner.N0,owner.c.mpf(1)/owner.N0)
            for key,q in normalized['drivers'].items():
                for part,row in (('kernels',q.value),('Z_derivatives',q.Z)):
                    difference=direct[part][key]*owner.N0-row
                    lo,hi=current.ep(difference.coefficient)
                    assert lo<=0<=hi,(key,part)
                    comparisons+=1
        try:owner.at_N(current.CELLS[0],'first_micro',(3,4),(1,1),N=owner.N0-1)
        except ValueError:pass
        else:raise AssertionError('Unadmitted frequency accepted')
        for flag in (*current.OPEN,'actual_uniform_N_scaled_Rc_targets_installed',
                'actual_all_regions_all_N_adapter_installed','actual_N_squared_averaging_cancellation_proved',
                'actual_terminal_controls_installed','actual_global_frequency_admitted',
                'physical_original_exterior_five_targets_closed'):
            assert saved[flag] is False
        hashes=dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,
            fixed_leading_input_N0=owner.N0,independent_normalized_density_fixture=fixture,
            independent_uniform_source_union_fixture=union_fixture,
            actual_uniform_bridge_integration_cells_replayed=cells,
            exact_full_phase_pieces_checked=phase_pieces,
            actual_first_active_uniform_C1_source_cell_recomputed=True,
            actual_fixed_N0_same_primitive_density_consistency_comparisons=comparisons,
            all_four_Z_original_P0_and_zero_inlet_and_pressure_memory_retained=True,
            unchanged_local_source_integrations_not_repeated=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,actual_all_regions_all_N_adapter_installed=False,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: original whole-Z full-phase all-N bridge source functions and R100 transport',flush=True)
    return receipt


if __name__=='__main__':run()
