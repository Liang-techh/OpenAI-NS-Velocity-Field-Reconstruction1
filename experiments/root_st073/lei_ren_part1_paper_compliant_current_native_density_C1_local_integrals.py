"""Original signed density Z jets and actual local C1 Duhamel contributions.

All roots/velocity jets share the same native packet, formal basis and ledger.
The original radius phase and fixed O2 cell/weights are Z independent. Local
function contributions do not supply global incoming histories or Rc repair.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_phase_first_jets as first
import lei_ren_part1_paper_compliant_current_native_local_signed_integrals as local

slow=first.slow;current=first.current;prior=first.prior;native=first.native;packets=first.packets
spatial=first.spatial;density=spatial.density;phase=first.phase
HERE,PREFIX,sha=first.HERE,first.PREFIX,first.sha;ep=first.ep;ZERO=(0,0)
NAME=PREFIX+'current_native_density_C1_local_integrals.json'
RECEIPT=PREFIX+'current_native_density_C1_local_integrals_check.json'
GATE='current_original_native_signed_density_Z_and_local_C1_Duhamel_functions_executed'
RATES=local.RATES


def density_Z_kernels(E,E_Z,V,V_Z,primitives,N):
    """Differentiate the original velocity increments and every signed term."""
    N=density.candidate_integer(N);c=E.ctx;factor=c.mpf(1)/N
    A=primitives['A']*factor;A_Z=primitives['A_Z']*factor
    increment=density.factored_expm1(A)
    exponential=c.exp(phase.bounded_value(A))
    deltaE=E*increment;deltaE_Z=E_Z*increment+E*exponential*A_Z
    deltaV=primitives['B_over_Pstar']*factor;deltaV_Z=primitives['B_Z_over_Pstar']*factor
    kernels=density.signed_density_kernels(E,V,deltaE,deltaV)
    theta_cross_Z=E_Z*deltaE+E*deltaE_Z;theta_square_half_Z=deltaE*deltaE_Z
    derivatives=dict(m=deltaV_Z,h=deltaE_Z,
        k=V_Z*deltaE+V*deltaE_Z+E_Z*deltaV+E*deltaV_Z+deltaE_Z*deltaV+deltaE*deltaV_Z,
        e=V_Z*deltaV*2+V*deltaV_Z*2+deltaV*deltaV_Z*2-theta_cross_Z-theta_square_half_Z,
        p=theta_cross_Z+theta_square_half_Z)
    velocities=dict(original_E=E,original_E_Z=E_Z,original_V=V,original_V_Z=V_Z,
        deltaE=deltaE,deltaE_Z=deltaE_Z,deltaV=deltaV,deltaV_Z=deltaV_Z,
        candidate_E=E+deltaE,candidate_E_Z=E_Z+deltaE_Z,
        candidate_V=V+deltaV,candidate_V_Z=V_Z+deltaV_Z)
    return dict(kernels=kernels,Z_derivatives=derivatives,velocities=velocities)


class NativeDensityC1LocalIntegrals:
    def __init__(self,owner):
        if type(owner) is not first.NativePhaseFirstJets:raise ValueError('Same original first-jet owner required')
        receipt=json.loads((HERE/first.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[first.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Accepted original phase/A/B first jets required')
        if {key:Fraction(value) for key,value in packets.recovery.RATES.items()}!=RATES:
            raise ValueError('Unchanged original five Duhamel rates required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (first.RECEIPT,Path(__file__).name,
            PREFIX+'current_native_candidate_densities.py',PREFIX+'current_native_local_signed_integrals.py',
            PREFIX+'current_generic_shear_signed_jets.py',PREFIX+'current_generic_shear_moment_recovery.py')})

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        query=self.owner.spatial_query(chart,Z,coordinate,N);source=query['source']['source'];packet=source['packet']
        E=source['roots']['E'][ZERO];E_Z=source['roots']['E'][(0,1)]
        signed_owner=self.owner.owner.owner.owner
        # Packet axial velocity is already in the common normalized Pstar
        # unit. k! Taylor coefficient_k is the ordinary Z derivative.
        leaf=lambda k:prior.signed.expressions.RadiusPolynomial(packet.algebra,
            {0:prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)})
        V,V_Z=[signed_owner.leaf(leaf(k),E.scale.bases,E.ledger) for k in (0,1)]
        cells=[]
        for original in query['cells']:
            if original['values'] is None:
                cells.append(dict(record=dict(status=original['record']['status'],original_phase_first_jet_source=original['record']),values=None))
                continue
            got=density_Z_kernels(E,E_Z,V,V_Z,original['values'],N)
            record=dict(status='enclosed',original_phase_first_jet_source=original['record'],
                original_and_candidate_normalized_velocity_Z_enclosures={k:v.record() for k,v in got['velocities'].items()},
                five_original_signed_density_C0_enclosures={k:v.record() for k,v in got['kernels'].items()},
                five_original_signed_density_Z_derivative_enclosures={k:v.record() for k,v in got['Z_derivatives'].items()},
                tiny_original_expm1_increment_retained=True,original_axial_V_and_every_cross_term_retained=True,
                original_radius_phase_Z_derivative_exactly_zero=True,same_original_factor_basis_and_ledger=True,
                source_native_width_or_Pstar_conversion_not_reapplied=True,
                C1_integrals_or_global_cumulative_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
            cells.append(dict(record=record,values=got))
        record=dict(source_family=self.family,chart=chart,source_provenance=packet.provenance,candidate_N=N,
            actual_original_radius_phase=query['geometry']['record'],original_q_slow_jet_source=query['source']['record'],
            original_V_and_V_Z_source='same packet.velocity axial row0; ordinary_axial_coefficient(k=0,1); same E factor basis/ledger',
            actual_spatial_signed_density_Z_cells=[cell['record'] for cell in cells],
            five_signed_density_Z_functions_installed=all(cell['values'] is not None for cell in cells),
            whole_chart_or_global_C1_history_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,source=source,geometry=query['geometry'])

    @native.inlet.source_precision
    def contribution(self,*,Z,left,right,N,chart='O2_slope'):
        if chart!='O2_slope':raise ValueError('Local C1 integral requires the original installed dy=dcoordinate O2_slope chart')
        lo=spatial.exact_coordinate(left);hi=spatial.exact_coordinate(right)
        if lo is None or hi is None or not 0<=lo<hi<=1:
            raise ValueError('Explicit exact ordered O2_slope endpoints in[0,1] required')
        c=self.ctx;cv=lambda value:c.mpf(value.numerator)/value.denominator
        coordinate=c.mpf((ep(cv(lo))[0],ep(cv(hi))[1]));width=cv(hi-lo)
        query=self.spatial_query(chart,Z,coordinate,N)
        if any(cell['values'] is None for cell in query['cells']):
            raise ArithmeticError('Refine original source/phase before C1 integration')
        kernels={key:local.same_source_union([cell['values']['kernels'][key] for cell in query['cells']]) for key in RATES}
        jets={key:local.same_source_union([cell['values']['Z_derivatives'][key] for cell in query['cells']]) for key in RATES}
        masses={key:local.positive_kernel_mass(c,width,rate) for key,rate in RATES.items()}
        contributions={key:kernels[key]*masses[key] for key in RATES}
        derivatives={key:jets[key]*masses[key] for key in RATES}
        record=dict(source_family=self.family,chart=chart,Z_box=query['geometry']['raw']['Z'],candidate_N=N,
            exact_native_coordinate_endpoints=[spatial.fractional_record(lo),spatial.fractional_record(hi)],
            exact_log_radius_width_fraction=spatial.fractional_record(hi-lo),log_radius_width_enclosure=width,
            original_coordinate_identity='dy=dcoordinate on O2_slope; original dy/dcoordinate=1',
            original_whole_cell_signed_density_Z_source=query['record'],
            whole_cell_signed_density_C0_covers={k:v.record() for k,v in kernels.items()},
            whole_cell_signed_density_Z_derivative_covers={k:v.record() for k,v in jets.items()},
            original_positive_Duhamel_kernel_mass_covers=masses,
            five_actual_local_signed_Duhamel_contribution_C0_enclosures={k:v.record() for k,v in contributions.items()},
            five_actual_local_signed_Duhamel_contribution_Z_derivative_enclosures={k:v.record() for k,v in derivatives.items()},
            original_integral='I_j(Z)=integral_left^right exp(-lambda_j*(right-s))*f_j(s,Z) ds',
            original_integral_Z_derivative='I_j_Z(Z)=integral_left^right exp(-lambda_j*(right-s))*f_j_Z(s,Z) ds',
            derivative_under_integral_theorem='Original source and cutoff are smooth on the covered compact radial/Z cell; fixed endpoints/radius phase and Duhamel weights are Z independent. Whole-cell continuous derivative covers permit differentiation under the finite integral.',
            recovery_rates={key:str(value) for key,value in RATES.items()},
            entire_radial_and_Z_cell_source_covers_used=True,actual_phase_union_not_samples=True,
            positive_kernel_signed_range_integration=True,Z_independent_endpoints_and_weights_no_boundary_terms=True,
            original_P0_same_common_unit_unchanged=query['geometry']['raw']['P0'],
            no_extra_N_or_radius_or_Pstar_unit_factor=True,pressure_rate_zero_and_memory_not_reset=True,
            actual_local_C1_Z_integral_functions_installed=True,local_contributions_are_not_global_defect_histories=True,
            original_incoming_history_not_assumed_or_reset=True,
            global_C1_histories_or_Rc_targets_or_repair_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,contributions=contributions,Z_derivatives=derivatives,
            kernels=kernels,kernel_Z_derivatives=jets,width=width,source=query['source'],geometry=query['geometry'])


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeDensityC1LocalIntegrals(first.NativePhaseFirstJets(slow.NativeQSlowJets(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))))
        N=1024;points=dict(inner_reference='.1337',O2_slope='.1337',O2_buffer='5.337',O3_slope_mu='.537',O3_power={'original_power_offset':'.537'})
        records={}
        for chart,coordinate in points.items():
            records[chart]=owner.spatial_query(chart,('.5','.5'),coordinate,N)['record']
            print('Original actual signed density Z jets:',chart,flush=True)
        integrals={}
        for name,Z in (('Z_half',('.5','.5')),('Z_interval',('.49','.51'))):
            got=owner.contribution(Z=Z,left='.13369999',right='.13370001',N=N);integrals[name]=got['record']
            print('Actual local C1 signed integrals:',name,{k:v.record()['sign'] for k,v in got['Z_derivatives'].items()},flush=True)
        if not all(record['five_signed_density_Z_functions_installed'] for record in records.values()):
            raise ArithmeticError('Declared native signed-density Z source unresolved')
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=N,
        native_spatial_signed_density_Z_records=records,actual_local_C1_signed_integral_records=integrals,
        native_point_query_count=5,native_local_integral_Z_query_count=2,
        actual_signed_density_Z_functions_installed=True,actual_local_signed_integral_C1_Z_functions_installed=True,
        first_Z_derivative_order_only=True,global_cumulative_histories_or_Rc_repair_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Original five signed density C0/Z functions at actual candidate N*y phase on five native boxes and two whole O2 radial/Z cells, with ten local signed integral values and ten genuine Z derivative enclosures. First Z only; no whole-chart coverage, incoming histories, cumulative Rc functions, repair/common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
