"""Directed C0/C1 enclosures of the true complete Rp-to-Rv pressure.

The selected native source is identified before bounds are transferred.
Quadratic derivative energy bounds already have physical pressure units;
the old extra A/reference multiplication is deliberately not inherited.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_pulse_native_numeric_identity as native
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C1 as caps
import lei_ren_part1_paper_compliant_outer_pulse_map as original_map
import sympy as sy
from types import SimpleNamespace
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

true=native.true
HERE,PREFIX,sha=native.HERE,native.PREFIX,native.sha
ep,record=caps.ep,caps.record
NAME=PREFIX+'current_original_n1_pulse_numeric_bridge.json.gz'
RECEIPT=PREFIX+'current_original_n1_pulse_numeric_bridge_check.json'
GATE='current_original_true_pulse_pressure_C1_numeric_bridge_installed'
OPEN=tuple(k for k in native.OPEN if k!=GATE)


def corrected_derivative_errors(numerical):
    """epsilon is already normalized by reference*A*Rv**p, once."""
    c=numerical.c
    q=numerical.quadratic
    errors={key:caps._scalar_bound(c,q[name]/2) for key,name in
        (('m13','reduced_S0_Z_absolute_bound'),('m15','reduced_S1_Z_absolute_bound'))}
    errors['m15']+=numerical.viscous_Z
    return errors


def normalization_functions(pulse):
    """True physical and normalized C1 function handles, no interval leaves."""
    g=pulse.g;a=true.parent.C1Algebra(g);own=pulse.own
    Gamma=own.jets.complete.Gamma;W=Gamma.pulse.functions
    mu=Gamma.schedule['mu'];logRp=Gamma.schedule['logRp']
    logP=true.source.FunctionRef(g,Gamma.pulse.preRp.operator.background['parameters']['logP'])
    u=own.jets.native['u'];U0=a.at(true.source.C1Function(u.value,u.Z),g.symbol('Z'),g.zero).value
    logU0=g.unary('log',U0)
    reference=g.mul(g.constant('1099/100'),g.unary('sqrt',g.mul(g.constant(2),g.node('mathematical_pi'))),
        g.sub(g.unary('exp',g.constant(5)),g.unary('exp',g.constant(3))))
    logref=g.unary('log',reference)
    divide=lambda x:g.quotient(x,mu,'positive actual original mu')
    logRv=g.add(logRp,divide(g.constant(13)))
    logA=g.add(logP,logU0,W['log_end_scale'],divide(g.constant('-13/2')),g.constant(-13))
    logs={key:g.add(logref,logP,logU0,W['log_end_scale'],g.mul(g.constant(p),logRp),
        divide(g.constant(f'{26*p-13}/2')),
        g.constant(-13)) for key,p in (('m13',1),('m15',2))}
    functions={key:a.div(pulse.functions[key],a.fixed(g.unary('exp',log)),
        'positive reference*A*Rv power; all factors Z independent') for key,log in logs.items()}
    lower=pulse.bindings['global_finite_endpoints'][0]
    scaled_lower=a.div(lower['R_q'],a.fixed(g.unary('exp',g.add(logA,logRv))),
        'positive A*Rv; incoming Rq uses no extra reference or radius ratio')
    return functions,dict(reference=reference,logA=logA,logRv=logRv,
        physical_unit_logs=logs,actual_native_U0=U0,log_native_U0=logU0,scaled_lower_R_q=scaled_lower)


def exact_scale_and_incoming_endpoint_identity(owner):
    """Bind actual physical scale recipes and lower endpoint function bodies."""
    g=owner.g;source=owner.numeric.pulse;frame=owner.native.selected.inlet.exact.frame
    Gamma=owner.pulse.own.jets.complete.Gamma;W=Gamma.pulse.functions
    importer=true.parent.GraphImport(g,frame.graph);ref=frame.functions['original_logRp']
    if source.units['source_Rp_node']!=ref.node or importer(ref)!=Gamma.schedule['logRp']:
        raise ValueError('Actual absolute Rp source nodes differ')
    mu,logC=sy.symbols('mu logC',positive=True);c=SimpleNamespace(ln=sy.log,exp=sy.exp,sqrt=sy.sqrt)
    logP=sy.exp(40)+11
    parameters=frame.bridge.leading.parameters
    parameter=lambda key:importer(true.source.FunctionRef(frame.graph,parameters[key].node))
    bindings={parameter('mu').node:mu,parameter('logC').node:logC}
    # The absolute radius is represented as phase-origin plus offset. Its
    # hbB/sc terms cancel exactly; keep them symbolic until that identity.
    for key in ('T','hbB','hbS','sc'):bindings[parameter(key).node]=sy.Symbol('original_'+key,positive=True)
    actual_logP=native.scalar_expression(g,parameter('logP'),{})
    if sy.simplify(actual_logP-logP)!=0:raise ValueError('Actual Pstar recipe differs')
    bindings[parameter('logP').node]=logP
    actual_logRp=native.scalar_expression(g,Gamma.schedule['logRp'],bindings)
    env=dict(c=c,logC=logC,self=SimpleNamespace(mu=mu))
    numeric=native.execute_assignments(caps.end.CurrentOriginalN1EndPressureSource.__init__,
        ('self.logRc','self.logRp','self.logP'),env)['self']
    if sy.simplify(actual_logRp-numeric.logRp)!=0 or sy.simplify(numeric.logP-logP)!=0:
        raise ValueError('True and numeric Rp/Pstar scale recipes differ')
    pref=native.execute_assignments(original_map.pulse_rows,('k0','u0','L','logpref'),
        dict(c=c,mu=mu))['logpref']
    old=native.execute_assignments(original_map.SharedOuterPulseMap.__init__,('self.logscale',),
        dict(c=c,self=SimpleNamespace(mu=mu,rows={'common_logpref':pref})))['self'].logscale
    actual=native.scalar_expression(g,W['log_end_scale'],{W['mu'].node:mu})
    if sy.simplify(actual-old)!=0:raise ValueError('True and numeric end amplitude scale recipes differ')
    # C3 native source identities were proved by the required native join.
    # Project the actual lower Rq graph using those exact functions.
    z,de,U0,Cnative,LP,LR,LE=sy.symbols('Z delta U0 Cnative logP logRp logE',positive=True)
    native_functions=owner.pulse.own.jets.native
    b={g.symbol('Z').node:z,Gamma.schedule['mu'].node:mu,Gamma.schedule['delta'].node:de,
        Gamma.schedule['logRp'].node:LR,W['log_end_scale'].node:LE,
        owner.normalization['actual_native_U0'].node:U0,
        true.source.FunctionRef(g,Gamma.pulse.preRp.operator.background['parameters']['logP']).node:LP}
    for key,body in (('u',U0/(1+z*z)),('m1',Cnative*z*(1+z*z))):
        jet=native_functions[key]
        for n,ref in enumerate((jet.value,jet.Z,jet.ZZ,jet.ZZZ)):b[ref.node]=sy.diff(body,z,n)
    pair=owner.normalization['scaled_lower_R_q']
    value=native.scalar_expression(g,pair.value,b);derivative=native.scalar_expression(g,pair.Z,b)
    factor=sy.exp(-13/(2*mu)+13-LE)
    if sy.simplify(value+Cnative*factor)!=0 or sy.simplify(derivative)!=0:
        raise ValueError('Actual normalized lower Rq/Z1 endpoint differs from native qinit')
    return dict(passed=True,actual_absolute_Rp_node_import_identity=True,
        actual_Pstar_Rp_and_end_amplitude_log_recipe_identities=3,
        native_U0_matches_required_exact_scalar_join=True,
        actual_lower_Rq_over_A_Rv_value_and_Z1_function_identities=2,
        native_qinit_is_minus_Cnative_times_mu_incoming_factor=True,
        native_qinit_Z1_exact_zero=True,finite_C0_weighted_viscous_endpoint_not_zeroed=True)


@dataclass(frozen=True,eq=False)
class OriginalTruePulseNumericPacket:
    scope:str='true complete original Rp-to-Rv pressure value/Z1 enclosures'


class CurrentOriginalN1PulseNumericBridge:
    @source_precision
    def __init__(self,native_identity,numerical,require_checked=True):
        if type(native_identity) is not native.CurrentOriginalN1PulseNativeNumericIdentity or not native_identity.acceptance_loaded:
            raise ValueError('Checked actual native-to-selected source identity required')
        if type(numerical) is not caps.CurrentOriginalN1PulsePressureC1Source or not numerical.acceptance_loaded:
            raise ValueError('Live checked same-source pulse C1 enclosure provider required')
        if numerical.pulse.coefficients.selected is not native_identity.selected:
            raise ValueError('Numerical pressure provider uses another selected source owner')
        self.native,self.numeric=native_identity,numerical
        self.pulse=native_identity.pulse;self.g=self.pulse.g;self.c=numerical.c
        self.family=copy.deepcopy(self.pulse.family)
        if numerical.family!=self.family:raise ValueError('True and numeric pulse families differ')
        self.functions,self.pieces=self.pulse.functions,self.pulse.pieces
        # The true pulse owner protects its entire original prefix. Snapshot
        # the remaining graph too, including normalization nodes reused by
        # an accepted second constructor rather than newly appended nodes.
        start=len(self.pulse.prefix)
        self.normalized_functions,self.normalization=normalization_functions(self.pulse)
        self._normalization_span=(start,len(self.g.nodes))
        self._normalization_nodes=copy.deepcopy(self.g.nodes[start:])
        self.bindings=true.own.encode(dict(physical=self.functions,normalized=self.normalized_functions,
            normalization=self.normalization,pieces=self.pieces))
        self._bindings=copy.deepcopy(self.bindings)
        self.derivative_errors=corrected_derivative_errors(numerical)
        self.hashes=dict(native_identity.hashes)
        for name,digest in numerical.hashes.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting pulse source '+name)
            self.hashes[name]=digest
        for module in (native,true,caps,caps.pulse):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.hashes[Path(original_map.__file__).name]=sha(Path(original_map.__file__).name)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed pulse enclosure source '+name)
        self.source_identity=self.assert_source_identity()
        self._identity=copy.deepcopy(self.source_identity)
        self._state=record((self.derivative_errors,numerical.energy,numerical.quadratic,numerical.source_cells,
            numerical.pulse._cells,numerical.pulse.Fprime,numerical.pulse.Wplus,numerical.pulse.rho,
            numerical.pulse.bounds,numerical.pulse.units))
        self._owners=(self.native,self.numeric,self.pulse,self.g,self.c)
        self._hashes=dict(self.hashes);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=true.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or archive!=sha(NAME) or receipt.get('source_family')!=self.family
                or receipt.get('source_identity')!=self.source_identity):
                raise ValueError('Checked true complete pulse numerical bridge required')
            self.acceptance_loaded=True

    @source_precision
    def assert_source_identity(self):
        self.native.assert_graph();self.numeric.assert_graph()
        own=self.pulse.own;complete=own.jets.complete;Gamma=complete.Gamma
        before=len(self.g.nodes)
        ends,certificate=true.actual_outer_endpoints(own)
        replay,pieces,bindings=true.pulse_pressure_functions(self.g,Gamma.schedule['delta'],Gamma.schedule['mu'],own.fields,ends)
        if (true.own.encode(replay)!=true.own.encode(self.functions) or true.own.encode(pieces)!=true.own.encode(self.pieces)
            or true.own.encode(bindings)!=true.own.encode(self.pulse.bindings)
            or certificate!=self.pulse.endpoint_certificate or len(self.g.nodes)!=before):
            raise ValueError('Actual pulse chart/endpoint/function replay differs')
        selected,equations=Gamma.pulse.original_coefficient_functions(complete.functions['weighted'])
        if (true.parent.encode(selected)!=true.parent.encode(complete.pulse_functions)
            or true.parent.encode(equations)!=true.parent.encode(complete.pulse_equations) or len(self.g.nodes)!=before):
            raise ValueError('Actual finite-mu positive selector or complete future changed')
        identity=complete.source_identity;selector=complete.pulse_branch_certificate
        if (not identity['passed'] or not identity['actual_C5_future_has_same_five_signed_complete_pieces']
            or not identity['graph_mu_exact_recipe_matches_original_selected_parameter']
            or not identity['complete_Gamma_value_and_first_Z_source_integrands_identified']
            or not selector['passed'] or not selector['actual_nu_strictly_positive_by_exact_product_mu_Kj_exp_2logE']):
            raise ValueError('Exact full future and nonzero-energy positive selector required')
        if not all(self.native.selected.current_selection_source_proof[k] for k in
            ('passed','checked_defining_parameter_identity_consumed','nonzero_formal_end_energy_and_incoming_terms_retained',
             'original_source_integrals_not_replaced_by_enclosure_values')):
            raise ValueError('Current numerical complete future/selector source proof required')
        source=self.numeric.pulse;c=self.c
        reference=c.mpf('10.99')*c.sqrt(2*c.pi)*(c.exp(5)-c.exp(3))
        if (source.reference._mpi_!=reference._mpi_ or ep(reference)[0]<=0
            or not source.units['all_original_pulse_charts_combined']
            or not source.units['terminal_M_at_Rv_exact_zero']
            or not source.units['integer_powers_combined_before_boxing']
            or len(self.numeric.source_cells)!=4 or len(source._cells)!=4):
            raise ValueError('Original pulse coverage or physical units differ')
        if tuple(tuple(row['ends']) for row in source._cells)!=tuple(caps.CELLS):
            raise ValueError('Original coefficient cell coverage differs')
        for row in self.numeric.source_cells:
            if row['source_order']<3 or row['incoming_mean_order']<3 or len(row['B'])!=4 or len(row['Tinitial'])!=4:
                raise ValueError('qZZ needs actual third-order source bounds')
        for energy in (source.energy,self.numeric.energy):
            if not energy['source_energy_bounds_nonzero'] or not energy['exact_units_combined_before_boxing']:
                raise ValueError('Original physical energy or unit source lost')
            if any(ep(value)[1]>=-1000 for row in energy['log_caps'] for value in row.values()):
                raise ValueError('Actual normalized pulse energy bound exceeds positive cover')
        if ep(source.Kmain)[0]<0 or ep(source.Kmain)[1]<=0 or ep(source.Rp_over_Rv)[1]<=0:
            raise ValueError('Positive main source or finite lower-radius ratio lost')
        if source.index.source_binding['selected_branch_unchanged'] is not True:
            raise ValueError('Actual original Fprime source changed')
        scale_proof=exact_scale_and_incoming_endpoint_identity(self)
        return dict(passed=True,exact_main_gap_end_pressure_C1_graph_and_endpoints_replayed=True,
            exact_native_numeric_source_identity=native.GATE,
            exact_finite_mu_selector_and_full_future_functions_replayed=True,
            complete_future_source_identity=identity,
            numerical_selection_source_proof=self.native.selected.current_selection_source_proof,
            all_three_charts_and_native_mean_histories_retained=True,
            physical_measure_main_gap_R_over_mu_end_R=True,
            auxiliary_q_is_Ur_over_sqrt_R_half_not_Ur_over_R=True,
            Pstar_applied_once_and_signed_m15_once=True,
            original_physical_units=copy.deepcopy(source.units),
            reference_is_positive_and_used_once=True,
            energy_epsilon_already_in_reference_A_Rv_pressure_units=True,
            corrected_C1_error_does_not_multiply_A_or_divide_reference_again=True,
            old_C1_local_grid_accuracy_not_used_as_transfer_evidence=True,
            all_nonzero_quadratic_finite_mu_and_weighted_viscous_errors_retained=True,
            exact_normalized_function_bindings=true.own.encode(self.normalized_functions),
            actual_physical_scale_and_incoming_endpoint_identities=scale_proof,
            inner_lower_preRp_and_reserved_repair_not_added=True)

    @source_precision
    def assert_graph(self):
        begin,end=self._normalization_span
        if (self._owners!=(self.native,self.numeric,self.pulse,self.g,self.c) or self.hashes!=self._hashes
            or self.source_identity!=self._identity or self.g.nodes[begin:end]!=self._normalization_nodes
            or true.own.encode(dict(physical=self.functions,normalized=self.normalized_functions,
                normalization=self.normalization,pieces=self.pieces))!=self._bindings
            or record((self.derivative_errors,self.numeric.energy,self.numeric.quadratic,self.numeric.source_cells,
                self.numeric.pulse._cells,self.numeric.pulse.Fprime,self.numeric.pulse.Wplus,self.numeric.pulse.rho,
                self.numeric.pulse.bounds,self.numeric.pulse.units))!=self._state):
            raise ValueError('Actual true pulse numerical attachment changed')
        self.native.assert_graph();self.numeric.assert_graph()

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();c=self.c;z=c.mpf(Z);lo,hi=ep(z)
        if lo<-1 or hi>1:raise ValueError('Original real Z domain [-1,1] required')
        sources=[];rows={key:[] for key in ('m13','m15')};active=[]
        for index,cell in enumerate(self.numeric.pulse._cells):
            left,right=map(c.mpf,cell['ends']);a=max(lo,ep(left)[0]);b=min(hi,ep(right)[1])
            if a>b:continue
            query=[a,b];active.append(index)
            packet=self.numeric.pulse.evaluate(query);value=self.numeric.pulse.report(packet)
            dp=self.numeric.evaluate(query);derivative=self.numeric.report(dp)
            for key in rows:
                linear=derivative['linear_pressure_Taylor'][key][1]
                error=ep(self.derivative_errors[key])[1]
                dz=linear+c.mpf([-error,error])
                if key=='m15':dz+=derivative['weighted_viscous_endpoint_derivative']
                rows[key].append([value['full_pressure_values'][key],dz])
            sources.append(dict(query=c.mpf(query),exact_source_coefficient_cell=value['source_coefficient_cell'],
                source_linear_inputs=value['source_linear_inputs'],
                C0_quadratic_errors=value['quadratic_source_errors'],
                corrected_C1_errors=self.derivative_errors,
                weighted_viscous_endpoint=value['weighted_viscous_endpoint'],
                weighted_viscous_endpoint_derivative=derivative['weighted_viscous_endpoint_derivative']))
        if not active:raise ValueError('Missing actual pulse coefficient cell coverage')
        hull=lambda values:c.mpf([min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)])
        view=dict(source_family=self.family,Z=z,radial_scope='original complete Rp to Rv',
            true_function_bindings=self.bindings,active_cells=active,
            normalized_true_pulse_C1_enclosures={k:[hull([p[j] for p in values]) for j in (0,1)] for k,values in rows.items()},
            original_source_cell_views=sources,
            normalization=dict(physical_unit_log_function_handles=true.own.encode(self.normalization['physical_unit_logs']),
                numerical_source_units=self.source_identity['original_physical_units'],
                actual_numeric_reference=self.numeric.pulse.reference,
                denominator='reference*A*Rv^p; p=1 for m13, p=2 for m15; Z-independent'),
            C1_energy_normalization_double_scaling_corrected=True,
            ordinary_Z0_Z1_enclosures_of_true_normalized_functions=True,
            complete_history_and_finite_mu_endpoints_retained=True,
            inner_lower_and_preRp_pressure_not_added=True,
            no_physical_hyper_exponential_radius_or_amplitude_materialized=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalTruePulseNumericPacket();self._packets[id(packet)]=(packet,view,record(view));return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalTruePulseNumericPacket or entry is None or entry[0] is not packet or record(entry[1])!=entry[2]:
            raise ValueError('Live true pulse C1 enclosure packet required')
        view=caps.end.copy_view(self.c,entry[1]);view[GATE]=self.acceptance_loaded;return view


@source_precision
def run(native_identity,numerical):
    began=time.monotonic();owner=CurrentOriginalN1PulseNumericBridge(native_identity,numerical,False)
    packets={name:owner.evaluate(Z) for name,Z in (('point','.371'),('boundary','0'),('whole_Z',['-1','1']))}
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,true_function_bindings=owner.bindings,
        normalization_graph_span=owner._normalization_span,normalization_graph_nodes=owner._normalization_nodes,
        samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_COMPLETE_PULSE_C1_NUMERIC_ENCLOSURES_DEFINED',len(packets),flush=True)
    return owner,raw,packets
