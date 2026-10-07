"""True native chart lengths, factored C1 integration and inlet collar transport.

Widths are collected from the unchanged affine-radius identities before
endpoint subtraction. Microscopic bridge/switch lengths stay nonzero formal
factors. Whole-cell density covers integrate in log radius, so no second
native-coordinate Jacobian is applied. The missing active global route is
still explicit; only the accepted initial flat collar has known incoming data.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_C1_history_transfer as history

density=history.density;prior=history.prior;spatial=history.spatial;native=history.native
packets=history.packets;phase=density.phase;local=density.local
HERE,PREFIX,sha=history.HERE,history.PREFIX,history.sha;ep=history.ep;RATES=history.RATES
NAME=PREFIX+'current_native_true_chart_C1_transfer.json'
RECEIPT=PREFIX+'current_native_true_chart_C1_transfer_check.json'
GATE='current_original_true_chart_lengths_factored_C1_integrals_and_initial_collar_transport_executed'


class TrueChartGeometry:
    def __init__(self,binder,coordinates):
        if type(binder) is not spatial.NativeSpatialPhase:raise ValueError('Same original radius binder required')
        coordinates.require_family(binder.family)
        if not binder.identity['passed']:raise ValueError('Original radius/Jacobian/seam identities required')
        self.binder=binder;self.coordinates=coordinates;self.ctx=coordinates.ctx
    def cv(self,q):return self.ctx.mpf(q.numerator)/q.denominator
    def factored(self,value):
        # Radius binder bases are exactly zero; copy full cover to the common
        # arithmetic coordinates, without inventing a Pstar power.
        if any(ep(base)!=(0,0) for base in value.scale.bases):raise ValueError('Original geometry zero bases required')
        return prior.ScaledEnclosure(prior.FormalScale(self.coordinates.bases,offset=self.ctx.mpf(value.scale.evaluate())),
            self.ctx.mpf(value.coefficient),self.coordinates.ledger)
    def build(self,chart,regular,micro,coordinate=None,endpoint_record=None):
        c=self.ctx;regular=c.mpf(regular)
        parts={key:self.factored(self.binder.width(key,coefficient)) for key,coefficient in micro.items()}
        width=sum(parts.values(),self.coordinates.scalar(regular))
        if ep(width.coefficient)[0]<=0:raise ArithmeticError('Strict positive true chart length unresolved: '+chart)
        # This cover is used only in attenuation/mass functions. The true
        # microscopic positive width itself remains factored in every record.
        tiny=sum(parts.values(),self.coordinates.scalar(0))
        bounded_tiny=phase.bounded_value(tiny)
        scalar_cover=regular+bounded_tiny
        if ep(scalar_cover)[0]<0:raise ArithmeticError('Negative directed chart length cover')
        record=dict(chart=chart,source_family=self.coordinates.family,
            native_coordinate_box=coordinate,original_endpoint_specification=endpoint_record,
            regular_true_log_radius_width_cover=regular,
            microscopic_true_log_radius_width_components={k:v.record() for k,v in parts.items()},
            positive_true_log_radius_width=width.record(),
            width_collected_before_huge_absolute_radius_endpoint_subtraction=True,
            exact_original_radius_Jacobian_identities=self.binder.identity,
            scalar_width_cover_used_only_for_directed_kernel_bounds=scalar_cover,
            width_and_endpoints_independent_of_Z=True,
            no_radius_cap_or_selected_midpoint_used=True,native_velocity_width_conversion_not_reapplied=True)
        return dict(record=record,width=width,regular=regular,scalar_cover=scalar_cover,coordinate=coordinate)
    def regular_slope(self,chart):
        c=self.ctx;P=self.binder.fixed['logP'];C=self.binder.fixed['logC'];T=self.binder.fixed['T']
        return dict(bridge_macro=4*P+c.ln(100)-c.ln(4)+1000,
            switch_power=c.ln(110)-c.ln(100),reshape=T,inner_reference=10*(C+P)-T-8,
            axial_restore=c.mpf(1),restore_buffer=c.mpf(1),Rh_reference=c.mpf(1),O2_slope=c.mpf(1),
            O2_buffer=c.mpf(1),O3_slope_mu=c.mpf(1),O3_power=self.binder.fixed['Tw']).get(chart)
    def cell(self,chart,left,right):
        c=self.ctx
        special=isinstance(left,dict) or isinstance(right,dict)
        if special:
            if not isinstance(left,dict) or not isinstance(right,dict) or set(left)!=set(right):raise ValueError('Matching original endpoint expressions required')
            if chart=='bridge_first' and set(left)=={'selected_sc_multiple'}:
                lo=Fraction(left['selected_sc_multiple']);hi=Fraction(right['selected_sc_multiple'])
                if not 0<=lo<hi<=1:raise ValueError('Selected original collar fractions in[0,1] required')
                a=self.binder.sc*self.cv(lo);b=self.binder.sc*self.cv(hi)
                box=c.mpf((ep(a)[0],ep(b)[1]));return self.build(chart,0,{'bridge':self.binder.sc*self.cv(hi-lo)},box,[left,right])
            if chart=='O3_power' and set(left)=={'original_power_offset'}:
                lo=Fraction(left['original_power_offset']);hi=Fraction(right['original_power_offset'])
                if lo<0 or hi<=lo or ep(self.binder.fixed['Tw']-self.cv(hi))[0]<0:raise ValueError('Original ordered power offsets within Tw required')
                a=self.cv(lo)/self.binder.fixed['Tw'];b=self.cv(hi)/self.binder.fixed['Tw']
                box=c.mpf((ep(a)[0],ep(b)[1]));return self.build(chart,self.cv(hi-lo),{},box,[left,right])
            raise ValueError('Unsupported original endpoint expression')
        lo=spatial.exact_coordinate(left);hi=spatial.exact_coordinate(right)
        if chart not in native.DOMAINS or lo is None or hi is None or lo>=hi:raise ValueError('Exact ordered original native chart endpoints required')
        dl,dh=native.DOMAINS[chart];upper=ep(c.exp(1))[1] if dh=='e' else dh
        if lo<dl or ep(self.cv(hi))[1]>upper:raise ValueError('Original native chart domain required')
        a,b=self.cv(lo),self.cv(hi);box=c.mpf((ep(a)[0],ep(b)[1]));d=self.cv(hi-lo)
        regular=c.mpf(0);micro={}
        if chart in ('bridge_first','bridge_second'):micro={'bridge':d}
        elif chart in ('switch_first','switch_second'):micro={'switch':d}
        elif chart=='bridge_macro':regular=self.regular_slope(chart)*d;micro={'bridge':-2*d}
        elif chart=='switch_power':regular=self.regular_slope(chart)*d;micro={'switch':-2*d}
        elif chart=='actual_patch':regular=c.ln(b/a)
        elif chart=='O2_axial':
            M=c.mpf(self.binder.seed.params.Md)
            regular=c.exp(M*a)*c.expm1(M*d)
        else:regular=self.regular_slope(chart)*d
        return self.build(chart,regular,micro,box,[str(lo),str(hi)])
    def full_lengths(self):
        result={}
        for chart,(left,right) in native.DOMAINS.items():
            if chart=='actual_patch':got=self.build(chart,1,{},None,['1','original exp(1)'])
            else:got=self.cell(chart,left,right)
            result[chart]=got['record']
        return result


def true_width_kernel(coordinates,geometry,rate):
    """Directed mass/attenuation for tiny, moderate and huge true widths."""
    c=coordinates.ctx;r=Fraction(rate);width=geometry['width'];w=geometry['scalar_cover']
    if r<0:raise ValueError('Nonnegative original recovery rate required')
    if not r:return dict(decay=coordinates.scalar(1),mass=width,branch='exact_rate0')
    rr=c.mpf(r.numerator)/r.denominator;x=rr*w
    if ep(geometry['regular'])==(0,0):
        # Tiny source width stays in mass. The scalar factor is only a
        # directed analytic cover, never an assertion exp(-lambda*w)=1.
        factor=c.exp(c.mpf((-ep(x)[1],0)))
        return dict(decay=coordinates.scalar(factor),mass=width*factor,branch='formal_micro_width_directed_average')
    decay=coordinates.decay(w,r)
    if ep(x)[1]<ep(c.mpf('.5'))[0]:
        factor=c.exp(c.mpf((-ep(x)[1],0)))
        mass=width*factor;branch='small_width_directed_average'
    else:
        # exp(-x) may be unmaterializable; bounded_exp supplies an outward
        # nonzero upper tail rather than replacing the exponential by zero.
        tail=decay.bounded_exp(decay.scale.evaluate())
        mass=coordinates.scalar((1-tail)/rr);branch='directed_one_minus_exponential_mass'
    if ep(mass.coefficient)[0]<=0:raise ArithmeticError('Strict positive Duhamel mass cover failed')
    return dict(decay=decay,mass=mass,branch=branch)


def append_true_cell(operator,geometry,values,jets,family):
    operator.coordinates.require_family(family)
    if set(values)!=set(RATES) or set(jets)!=set(RATES):raise ValueError('All original signed C0/Z contribution rows required')
    for key,rate in RATES.items():
        decay=true_width_kernel(operator.coordinates,geometry,rate)['decay']
        operator.coefficients[key]=decay*operator.coefficients[key]
        operator.increments[key]=decay*operator.increments[key]+operator.coordinates.rebase(values[key],family)
        operator.Z_increments[key]=decay*operator.Z_increments[key]+operator.coordinates.rebase(jets[key],family)
    operator.steps+=1


class NativeTrueChartC1Transfer:
    def __init__(self,owner):
        if type(owner) is not history.NativeC1HistoryTransfer:raise ValueError('Accepted original inlet/C1 transfer owner required')
        receipt=json.loads((HERE/history.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[history.GATE] or receipt['source_family']!=owner.family:raise ValueError('Same original checked C1 transfer required')
        self.owner=owner;self.ctx=owner.ctx;self.coordinates=owner.coordinates;self.family=owner.family;self.service=owner.service
        self.geometry=TrueChartGeometry(owner.binder,self.coordinates)
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (history.RECEIPT,Path(__file__).name,
            spatial.RECEIPT,PREFIX+'current_generic_shear_loop_domain.json',PREFIX+'current_generic_shear_uniform_inputs.json')})
    @native.inlet.source_precision
    def contribution(self,*,chart,Z,left,right,N=1024):
        geometry=self.geometry.cell(chart,left,right)
        query=self.owner.owner.spatial_query(chart,Z,geometry['coordinate'],N)
        if any(cell['values'] is None for cell in query['cells']):raise ArithmeticError('Refine unresolved original source/phase before true-chart C1 integration: '+chart)
        kernels={key:local.same_source_union([cell['values']['kernels'][key] for cell in query['cells']]) for key in RATES}
        jets={key:local.same_source_union([cell['values']['Z_derivatives'][key] for cell in query['cells']]) for key in RATES}
        factors={key:true_width_kernel(self.coordinates,geometry,rate) for key,rate in RATES.items()}
        values={key:self.coordinates.rebase(kernels[key],self.family)*factors[key]['mass'] for key in RATES}
        derivatives={key:self.coordinates.rebase(jets[key],self.family)*factors[key]['mass'] for key in RATES}
        record=dict(source_family=self.family,chart=chart,Z_box=self.ctx.mpf(Z),candidate_N=N,
            actual_true_chart_cell_geometry=geometry['record'],original_whole_cell_signed_density_Z_source=query['record'],
            original_true_log_radius_kernel_factors={key:dict(branch=item['branch'],mass=item['mass'].record(),decay=item['decay'].record()) for key,item in factors.items()},
            signed_C0_contribution_enclosures={k:v.record() for k,v in values.items()},
            signed_Z_contribution_enclosures={k:v.record() for k,v in derivatives.items()},
            density_already_in_original_log_radius_coordinate=True,
            log_radius_Jacobian_integrated_in_true_width_and_not_multiplied_again=True,
            fixed_Z_independent_endpoints_and_phase_allow_differentiation_under_integral=True,
            positive_Duhamel_mass_times_whole_signed_function_covers=True,
            incoming_correction_not_assumed_or_reset=True,global_inlet_to_Rc_histories_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,geometry=geometry,contributions=values,Z_derivatives=derivatives)
    @native.inlet.source_precision
    def initial_collar(self,Z=(-1,1),right_multiple='3/4',N=1024):
        right=Fraction(right_multiple)
        if not Fraction(1,2)<right<=Fraction(3,4):raise ValueError('Transport stays inside accepted two-sided left collar[sc/4,3sc/4]')
        inlet=self.owner.inlet(Z,N)
        domain=json.loads((HERE/(PREFIX+'current_generic_shear_loop_domain.json')).read_bytes())['current_original_generic_loop_domain']
        uniform=json.loads((HERE/(PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
        c=self.ctx;read=lambda row:packets.interval(c,row)
        if domain['source_family']!=self.family or not domain['original_both_loop_edge_collars_strict'] or ep(read(domain['left_collar_fraction_box']))!=(ep(c.mpf('0.25'))[0],ep(c.mpf('0.75'))[1]):
            # Family and whole collar are required, not a neighboring point.
            raise ValueError('Same accepted complete left collar required')
        excess=read(domain['left_collar_kappa_minus2_lower']);eta_log=read(uniform['selected_positive_eta_log'])
        if ep(eta_log)[1]>ep(c.ln(excess)-c.ln(2))[0]:raise ValueError('Original eta/kappa flat collar inequality required')
        geometry=self.geometry.cell('bridge_first',{'selected_sc_multiple':'1/2'},{'selected_sc_multiple':str(right)})
        zeros={key:self.coordinates.scalar(0) for key in RATES};operator=history.C1DuhamelOperator(self.coordinates)
        append_true_cell(operator,geometry,zeros,zeros,self.family)
        correction=operator.apply(inlet['initial_defect'],inlet['initial_defect_Z'],self.family)
        coordinate=self.owner.binder.sc*self.geometry.cv(right)
        packet=self.owner.native.query('bridge_first',Z,coordinate)
        background=history.packet_history_functions(packet,self.coordinates,self.owner.signed_owner)
        own={key:background['originals'][key]+correction['values'][key] for key in RATES}
        own_Z={key:background['Z_derivatives'][key]+correction['Z_derivatives'][key] for key in RATES}
        record=dict(source_family=self.family,Z_box=c.mpf(Z),candidate_N=N,
            exact_collar_fraction_endpoints=['1/2',str(right)],actual_true_chart_cell_geometry=geometry['record'],
            original_whole_left_collar_fraction_box=read(domain['left_collar_fraction_box']),
            original_kappa_minus2_lower=excess,original_eta_log_cover=eta_log,
            exact_source_argument='On the entire accepted left collar, kappa-2>=excess>=2*eta and D>=0. Original cutoff q and all source/phase/Z derivatives are exactly zero, hence A=B=deltaE=deltaV=f_j=f_j_Z=0.',
            original_correction_initial_condition_is_the_checked_sc_half_inlet=True,
            actual_correction_C0_enclosures={k:v.record() for k,v in correction['values'].items()},
            actual_correction_Z_enclosures={k:v.record() for k,v in correction['Z_derivatives'].items()},
            original_right_background_and_separate_P0_Z=background['record'],
            actual_own_history_C0_enclosures={k:v.record() for k,v in own.items()},
            actual_own_history_Z_enclosures={k:v.record() for k,v in own_Z.items()},
            true_width_affine_operator=operator.record(),
            original_background_histories_and_P0_not_reset=True,
            actual_initial_inlet_to_collar_endpoint_C1_histories_installed=True,
            remaining_bridge_after_selected_collar_not_assumed_quiet=True,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,geometry=geometry,operator=operator,correction=correction,background=background,own=own,own_Z=own_Z)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        c1=density.NativeDensityC1LocalIntegrals(density.first.NativePhaseFirstJets(density.slow.NativeQSlowJets(density.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))))
        owner=NativeTrueChartC1Transfer(history.NativeC1HistoryTransfer(c1));N=1024
        lengths=owner.geometry.full_lengths();print('True original positive chart lengths:',len(lengths),flush=True)
        collars={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            collars[name]=owner.initial_collar(Z,N=N)['record'];print('Actual initial C1 collar transport:',name,flush=True)
        requests=dict(inner_reference=('.12','.15'),O2_slope=('.12','.15'),O2_buffer=('5.33','5.34'),
            O3_slope_mu=('.53','.54'),O3_power=({'original_power_offset':'.53'},{'original_power_offset':'.54'}))
        cells={}
        for chart,(left,right) in requests.items():
            cells[chart]=owner.contribution(chart=chart,Z=('.49','.51'),left=left,right=right,N=N)['record']
            print('Native true-chart signed C1 contribution:',chart,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=N,
        original_full_chart_positive_length_records=lengths,actual_initial_collar_C1_history_records=collars,
        actual_true_chart_C1_contribution_records=cells,positive_native_chart_length_count=len(lengths),
        native_true_chart_C1_cell_count=len(cells),initial_collar_Z_query_count=len(collars),
        actual_initial_C1_route_scope='true original sc/2 ->3sc/4 on bridge_first; whole Z[-1,1] and[.49,.51]',
        actual_inlet_to_O2_local_incoming_correction_installed=False,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='All17 original positive chart lengths and true-width C1 local contributions on five declared native cells, plus actual original inlet-to-three-quarter-sc C1 histories through the proven flat collar. No active bridge/global route, terminal repair, common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
