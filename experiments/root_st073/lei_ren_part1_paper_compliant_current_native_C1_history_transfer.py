"""Common directed source coordinates, original inlet and serial C1 transfer.

Variable width/amplitude/radius log covers move to formal offsets; only the
same Pstar squared log remains a shared symbolic base. Adjacent-cell signed
increments compose as an affine Duhamel operator with an explicit incoming
correction argument. The original background memory and P0 never reset.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals as density

first=density.first;slow=density.slow;current=density.current;prior=density.prior
native=density.native;packets=density.packets;spatial=density.spatial
HERE,PREFIX,sha=density.HERE,density.PREFIX,density.sha;ep=density.ep;RATES=density.RATES
NAME=PREFIX+'current_native_C1_history_transfer.json'
RECEIPT=PREFIX+'current_native_C1_history_transfer_check.json'
GATE='current_original_native_inlet_common_basis_and_adjacent_C1_Duhamel_operator_executed'


class CommonSourceCoordinates:
    """An arithmetic coordinate choice, never a selected source field value."""
    def __init__(self,c,logP_squared,family):
        if set(family)!=set(packets.FAMILY_KEYS):raise ValueError('Explicit original source family/datum required')
        self.ctx=c;self.family=dict(family);self.logP_squared=c.mpf(logP_squared)
        self.bases=(c.mpf(0),self.logP_squared,c.mpf(0),c.mpf(0),c.mpf(0))
        self.ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
            positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    def scalar(self,value):return prior.ScaledEnclosure(prior.FormalScale(self.bases),value,self.ledger)
    def require_family(self,family):
        if family!=self.family:raise ValueError('Same original source family and analytic axis datum required')
    def rebase(self,value,family):
        self.require_family(family)
        source_bases=tuple(self.ctx.mpf(base) for base in value.scale.bases)
        if source_bases[1]._mpi_!=self.logP_squared._mpi_:
            raise ValueError('Same accepted Pstar squared log cover required')
        if value.zero:return self.scalar(0)
        # Every entire source log cover is retained. No exponential, radius
        # cap, endpoint, or midpoint is converted into a defining field value.
        p=value.scale.powers;offset=self.ctx.mpf(value.scale.offset)+sum((source_bases[j]*p[j]
            for j in (0,2,3,4) if p[j]),self.ctx.mpf(0))
        return prior.ScaledEnclosure(prior.FormalScale(self.bases,(0,p[1],0,0,0),offset),self.ctx.mpf(value.coefficient),self.ledger)
    def decay(self,width,rate):
        w=self.ctx.mpf(width);r=Fraction(rate)
        if ep(w)[0]<0 or r<0 or any(not mp.isfinite(x) for x in ep(w)):raise ValueError('Finite nonnegative true log-radius width/rate required')
        if not r:return self.scalar(1)
        return prior.ScaledEnclosure(prior.FormalScale(self.bases,offset=-w*self.ctx.mpf(r.numerator)/r.denominator),1,self.ledger)
    def record(self):
        return dict(source_family=self.family,common_log_bases=self.bases,
            shared_symbolic_factor='same original Pstar squared log only',
            variable_width_F0_swirl_and_radius_log_covers_retained_in_formal_offsets=True,
            no_logarithmic_source_factor_or_radius_exponential_materialized=True,
            entire_directed_log_and_coefficient_covers_copied_to_common_context=True,
            common_basis_is_arithmetic_not_a_new_source_function=True,
            rebase_can_widen_dependency_covers_not_establish_joint_correlations=True)


class C1DuhamelOperator:
    """L(H_in)=decay*H_in+signed_increment for all five C0/Z rows."""
    def __init__(self,coordinates):
        self.coordinates=coordinates;self.coefficients={key:coordinates.scalar(1) for key in RATES}
        # These zeros define the affine operator's additive identity, not an
        # assertion that the actual incoming correction/history is zero.
        self.increments={key:coordinates.scalar(0) for key in RATES}
        self.Z_increments={key:coordinates.scalar(0) for key in RATES};self.steps=0
    def append(self,width,values,Z_values,family):
        self.coordinates.require_family(family)
        if set(values)!=set(RATES) or set(Z_values)!=set(RATES):raise ValueError('All five original C0/Z increment functions required')
        for key,rate in RATES.items():
            decay=self.coordinates.decay(width,rate)
            value=self.coordinates.rebase(values[key],family);jet=self.coordinates.rebase(Z_values[key],family)
            self.coefficients[key]=decay*self.coefficients[key]
            self.increments[key]=decay*self.increments[key]+value
            self.Z_increments[key]=decay*self.Z_increments[key]+jet
        self.steps+=1
    def apply(self,incoming,incoming_Z,family):
        self.coordinates.require_family(family)
        if set(incoming)!=set(RATES) or set(incoming_Z)!=set(RATES):raise ValueError('All five actual incoming C0/Z functions required')
        values={};jets={}
        for key in RATES:
            value=self.coordinates.rebase(incoming[key],family);jet=self.coordinates.rebase(incoming_Z[key],family)
            values[key]=self.coefficients[key]*value+self.increments[key]
            jets[key]=self.coefficients[key]*jet+self.Z_increments[key]
        return dict(values=values,Z_derivatives=jets)
    def record(self):
        return dict(steps=self.steps,original_rates={key:str(rate) for key,rate in RATES.items()},
            incoming_C0_Z_decay_coefficients={k:v.record() for k,v in self.coefficients.items()},
            cumulative_signed_increment_C0_enclosures={k:v.record() for k,v in self.increments.items()},
            cumulative_signed_increment_Z_enclosures={k:v.record() for k,v in self.Z_increments.items()},
            affine_transfer='H_out=exp(-lambda*total_width)*H_in+sum_cell exp(-lambda*downstream_width)*I_cell; same for Z rows',
            incoming_argument_not_assumed_or_reset=True,quiet_cells_preserve_original_rate0_pressure_memory=True,
            actual_global_incoming_correction_not_supplied=True,global_cumulative_histories_admitted=False,
            **dict.fromkeys(packets.OPEN,False))


def packet_history_functions(packet,coordinates,signed_owner):
    """Original normalized histories/P0 and ordinary Z rows, without recovery."""
    coordinates.require_family(packet.source_family);c=coordinates.ctx
    bases=tuple(c.mpf(base) for base in packet.algebra.logs)+(c.mpf(0),);ledger=dict(directed_small_exponential_tails=0,
        positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    def convert(row,k):
        ordinary=prior.signed.ordinary_axial_coefficient(row,k)
        leaf=prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:ordinary})
        return coordinates.rebase(signed_owner.leaf(leaf,bases,ledger),packet.source_family)
    originals={key:convert(packet.histories[key][0],0) for key in RATES}
    jets={key:convert(packet.histories[key][0],1) for key in RATES}
    P0=convert(packet.P0,0);P0_Z=convert(packet.P0,1)
    pressure=P0+originals['p'];pressure_Z=P0_Z+jets['p']
    record=dict(source_family=packet.source_family,source_provenance=packet.provenance,
        original_normalized_history_C0_enclosures={k:v.record() for k,v in originals.items()},
        original_normalized_history_Z_enclosures={k:v.record() for k,v in jets.items()},
        original_separate_P0_over_Pstar_squared=P0.record(),original_separate_P0_Z_over_Pstar_squared=P0_Z.record(),
        original_absolute_pressure_over_Pstar_squared=pressure.record(),original_absolute_pressure_Z_over_Pstar_squared=pressure_Z.record(),
        original_units=packets.recovery.UNITS,ordinary_Z_factorial_conversion=True,
        original_histories_not_zeroed_or_reset=True,P0_not_merged_into_pressure_history=True,
        native_Pstar_width_or_radial_factor_not_reapplied=True,common_directed_coordinate_theorem=coordinates.record())
    return dict(record=record,originals=originals,Z_derivatives=jets,P0=P0,P0_Z=P0_Z,pressure=pressure,pressure_Z=pressure_Z)


class NativeC1HistoryTransfer:
    def __init__(self,owner):
        if type(owner) is not density.NativeDensityC1LocalIntegrals:raise ValueError('Same original signed-density/local C1 owner required')
        receipt=json.loads((HERE/density.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[density.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Accepted original density/local C1 functions required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.native=owner.owner.native;self.service=owner.service
        self.signed_owner=owner.owner.owner.owner.owner
        self.coordinates=CommonSourceCoordinates(self.ctx,self.native.seed.logP*2,self.family)
        self.binder=owner.owner.binder
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (density.RECEIPT,Path(__file__).name,
            PREFIX+'current_native_generic_left_inlet.py',PREFIX+'current_generic_shear_source_packets.py')})

    @native.inlet.source_precision
    def inlet(self,Z=(-1,1),N=1024):
        inlet=self.native.left_inlet(Z);packet=inlet['original_source_packet']
        geometry=self.binder.query('bridge_first',Z,{'selected_sc_multiple':'1/2'},N)
        if geometry['raw']['coordinate']._mpi_!=inlet['source_query_coordinate']._mpi_ or not geometry['offset'].zero:
            raise ArithmeticError('Same actual positive sc/2 inlet and exact zero radius offset required')
        if not inlet['left_loop_q_A_B_exact_zero_by_existing_checked_collar'] or not inlet['original_P0_and_all_five_incoming_histories_retained']:
            raise ArithmeticError('Original flat inlet and incoming memory theorem required')
        got=packet_history_functions(packet,self.coordinates,self.signed_owner)
        zeros={key:self.coordinates.scalar(0) for key in RATES}
        if any(row.terms for row in inlet['generic_inlet_defect_exact_zero'].values()):raise ArithmeticError('Original inlet correction initial condition changed')
        record=dict(got['record'],actual_original_left_inlet_geometry=geometry['record'],
            original_checked_flat_inlet_argument=inlet['flat_left_branch_source_argument'],
            original_inlet_correction_C0_Z_exact_zero={key:zero.record() for key,zero in zeros.items()},
            inlet_correction_zero_is_checked_initial_condition_not_downstream_history_reset=True,
            actual_native_incoming_functions_and_P0_Z_installed=True,
            global_inlet_to_Rc_histories_or_repair_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,functions=got,initial_defect=zeros,initial_defect_Z=dict(zeros),geometry=geometry)

    @native.inlet.source_precision
    def serial(self,*,Z,endpoints,N=1024):
        points=[spatial.exact_coordinate(point) for point in endpoints]
        if len(points)<2 or any(point is None for point in points) or not 0<=points[0]<points[-1]<=1 or any(a>=b for a,b in zip(points,points[1:])):
            raise ValueError('Exact strictly increasing adjacent O2_slope endpoints in[0,1] required')
        operator=C1DuhamelOperator(self.coordinates);cells=[]
        text=lambda point:str(point.numerator)+'/'+str(point.denominator)
        for left,right in zip(points,points[1:]):
            got=self.owner.contribution(Z=Z,left=text(left),right=text(right),N=N)
            operator.append(got['width'],got['contributions'],got['Z_derivatives'],self.family)
            cells.append(got['record'])
        # Original background history comes directly from the unchanged
        # endpoint packet. The affine correction at this local left endpoint
        # remains an explicit input; no missing inlet-to-local gap is skipped.
        packet=self.signed_owner.packet('O2_slope',Z,self.ctx.mpf(points[-1].numerator)/points[-1].denominator)
        background=packet_history_functions(packet,self.coordinates,self.signed_owner)
        record=dict(source_family=self.family,Z_box=self.ctx.mpf(Z),candidate_N=N,
            exact_adjacent_coordinate_partition=[spatial.fractional_record(point) for point in points],
            exact_total_log_radius_width_fraction=spatial.fractional_record(points[-1]-points[0]),
            actual_whole_cell_C1_signed_contributions=cells,common_directed_coordinate_theorem=self.coordinates.record(),
            serial_signed_C0_Z_affine_Duhamel_operator=operator.record(),
            original_right_endpoint_background_history_and_P0_Z=background['record'],
            completion_formula='own_history_right=original_background_history_right+operator(actual_incoming_correction_at_local_left)',
            actual_incoming_correction_at_local_left_not_computed_or_assumed_zero=True,
            missing_actual_inlet_to_local_gap_not_treated_as_quiet=True,
            adjacent_local_cumulative_C1_operator_installed=True,
            global_inlet_to_Rc_histories_or_repair_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,operator=operator,background=background)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeC1HistoryTransfer(density.NativeDensityC1LocalIntegrals(first.NativePhaseFirstJets(slow.NativeQSlowJets(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge)))))))
        inlets={};serial={};partition=['.12','.13','.14','.15'];N=1024
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            inlets[name]=owner.inlet(Z,N)['record'];print('Actual original inlet histories/P0_Z:',name,flush=True)
        for name,Z in (('Z_half',('.5','.5')),('Z_interval',('.49','.51'))):
            got=owner.serial(Z=Z,endpoints=partition,N=N);serial[name]=got['record']
            print('Actual adjacent C1 Duhamel operator:',name,len(partition)-1,'cells, true width.03',flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=N,
        actual_original_inlet_C1_history_records=inlets,actual_adjacent_C1_serial_operator_records=serial,
        original_inlet_Z_query_count=2,adjacent_radial_cell_count=3,serial_Z_query_count=2,
        actual_original_inlet_history_and_P0_Z_functions_installed=True,
        actual_common_basis_and_adjacent_signed_C1_operator_installed=True,
        covered_O2_radial_coordinate_interval=['.12','.15'],exact_total_width='3/100',
        whole_chart_cancellation_or_global_cumulative_history_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Original actual sc/2 inlet five histories/Z and separate P0/P0_Z, plus common-coordinate affine signed C1 Duhamel transfer over three adjacent O2 cells[.12,.15] at Z=.5/[.49,.51]. Incoming local correction remains explicit; no missing global gap/reset, whole chart cancellation, Rc/repair/common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
