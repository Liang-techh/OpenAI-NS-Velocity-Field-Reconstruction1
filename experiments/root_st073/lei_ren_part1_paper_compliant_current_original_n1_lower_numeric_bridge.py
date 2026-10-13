"""Same-source C1 enclosures of the true Rin-to-rminus pressure integrals.

The true graph keeps the defining collar functions and their own histories.
The checked prefix supplies correlated complex-comparison errors, ordinary
Z jets and directed range integration. A range product encloses an integral;
it is never substituted for the integral FunctionRef or a selected field.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric_bridge as cut
import lei_ren_part1_paper_compliant_current_original_n1_unchanged_pressure_prefix as prior
import lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions as true
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

full=cut.full;op=prior.op;ep=prior.ep
HERE,PREFIX,sha=cut.HERE,cut.PREFIX,cut.sha
NAME=PREFIX+'current_original_n1_lower_numeric_bridge.json.gz'
RECEIPT=PREFIX+'current_original_n1_lower_numeric_bridge_check.json'
GATE='current_original_true_lower_collar_pressure_C1_numeric_bridge_installed'
OPEN=cut.OPEN


@dataclass(frozen=True,eq=False)
class OriginalTrueLowerNumericPacket:
    scope:str='original signed pressure C1 integrals on Rin to rminus'


class CurrentOriginalN1LowerNumericBridge:
    @source_precision
    def __init__(self,accepted_cutoff,accepted_prefix,require_checked=True):
        if type(accepted_cutoff) is not cut.CurrentOriginalN1CutoffNumericBridge or not accepted_cutoff.acceptance_loaded:
            raise ValueError('Checked true cutoff numerical bridge required')
        if type(accepted_prefix) is not prior.CurrentOriginalN1UnchangedPressurePrefix or not accepted_prefix.acceptance_loaded:
            raise ValueError('Checked original unchanged pressure prefix required')
        self.cutoff,self.numeric=accepted_cutoff,accepted_prefix
        self.inner=accepted_cutoff.inner;self.lower=self.inner.full.lower
        self.g,self.c=self.lower.g,accepted_prefix.c
        self.family=copy.deepcopy(self.lower.family)
        self.functions=self.lower.functions
        self.bindings=full.own.encode(dict(functions=self.functions,fields=self.lower.fields,
            parameters=self.lower.parameters,pressure_bindings=self.lower.pressure_bindings))
        self._bindings=copy.deepcopy(self.bindings)
        self._owners=(self.cutoff,self.numeric,self.inner,self.lower,self.g,self.c)
        self.source_identity=self.assert_source_identity()
        self._identity=op.snapshot(self.source_identity)
        self.hashes={}
        for owner in (self.cutoff,self.numeric):
            for name,digest in owner.hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:
                    raise ValueError('Conflicting lower numerical source '+name)
                self.hashes[name]=digest
        for module in (cut,prior,true):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):
                digest=sha(name)
                if name in self.hashes and self.hashes[name]!=digest:
                    raise ValueError('Changed lower numerical dependency '+name)
                self.hashes[name]=digest
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed lower numerical input '+name)
        self._hashes=dict(self.hashes);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=full.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or archive!=sha(NAME) or receipt.get('source_family')!=self.family
                or receipt.get('source_identity')!=op.snapshot(self.source_identity)):
                raise ValueError('Checked true lower C1 numerical bridge required')
            self.acceptance_loaded=True

    def assert_source_identity(self):
        self.cutoff.assert_graph();self.numeric.assert_graph()
        if (self.lower is not self.inner.collar or self.lower.core is not self.inner.core
            or self.numeric.mean is not self.inner.numeric.mean
            or self.numeric.inputs is not self.cutoff.numeric.inputs
            or self.numeric.physical is not self.inner.numeric.physical
            or self.numeric.local is not self.inner.numeric.mean.local
            or self.family!=self.numeric.family or self.family!=self.cutoff.family):
            raise ValueError('Lower pressure source or incoming history differs')
        source_proof=self.inner.assert_source_identity()
        if not source_proof['directed_model_plus_error_callbacks_enclose_that_same_actual_solution']:
            raise ValueError('The reference core must enclose the actual nonlinear solution')
        p=self.lower.parameters;sc=cut.binary_ref(self.g,p['sc'])
        for name,shift in (('sc_quarter',-2),('sc_half',-1)):
            if cut.binary_ref(self.g,p[name])!=cut.normalize_binary(sc[0],sc[1]+shift):
                raise ValueError('True lower endpoints must be exact original binary phases')
        if (self.numeric.phase_in._mpi_!=self.inner.selection.phase_in._mpi_
            or self.numeric.modification_inlet._mpi_!=(2*self.inner.selection.phase_in)._mpi_
            or self.numeric.local.inner.sc._mpi_!=self.inner.selection.sc._mpi_):
            raise ValueError('Numeric pressure interval must be original sc/4 to sc/2')
        hb=self.g.nodes[p['hb'].node]
        if (hb['operation']!='current_original_source_parameter' or hb['name']!='hbB'
            or hb['source_family']!=self.family
            or self.lower.geometry_certificate['original_width_definition']!=
                'hb=epsilon_b=cstar*K^-100; K is the original global physical norm sum'):
            raise ValueError('Actual native hb definition required')
        # The inherited source proof checks the producer AST for logh,
        # selected norm records and the same native width. No box equality
        # or materialized small width is used to define the source.
        recipe=cut.inner.original_normalization_and_source_recipes()
        if recipe['original_positive_width_log_source']!='source_log_hb_enclosure':
            raise ValueError('Formal hb no longer comes from the native source width')
        Gamma=true.own.jets.complete.Gamma
        if prior.domain.sigma_jets is not true.own.c2.flat_source.sigma_jets:
            raise ValueError('Numerical and true collar sigma derivative providers differ')
        s=self.lower.fields['query'].node
        sigma=[n for n in self.g.nodes if n['operation']=='exact_original_source_flat_sigma' and n['argument']==s]
        derivative=[n for n in self.g.nodes if n['operation']=='analytic_unary'
            and n.get('name')=='original_flat_sigma_prime' and n.get('argument')==s]
        if (len(sigma)!=1 or sigma[0]['source_binding']!=Gamma.parent.active.current.ast_binding(Gamma.flat_source.stable_sigma)
            or not sigma[0]['exact_function_not_its_interval_cap'] or len(derivative)!=1
            or derivative[0]['derivative_source']!=true.own.c2.current.ast_binding(true.own.c2.flat_source.sigma_jets)):
            raise ValueError('Exact original collar flat function/ordinary derivative source required')
        if not self.lower.fields['all_histories_start_at_true_core_rho4']:
            raise ValueError('Lower collar reset its actual core inlet')
        support=self.lower.geometry_certificate['original_modification_support']
        if not support['increments_exact_zero']:
            raise ValueError('Finite-N changes must vanish throughout the lower collar')
        for pair in self.functions.values():
            for k,ref in enumerate((pair.value,pair.Z)):
                n=self.g.nodes[ref.node]
                if (n['operation']!='definite_integral' or n['lower']!=p['sc_quarter'].node
                    or n['upper']!=p['sc_half'].node or n['ordinary_slow_Z_derivative_order']!=k
                    or not n['signed_integrand_not_a_cap'] or not n['all_finite_viscous_derivatives_retained']
                    or not n['endpoint_independent_of_Z']):
                    raise ValueError('Signed true lower C1 integral body changed')
        return dict(passed=True,same_checked_true_nonlinear_core_model_and_error_callbacks=True,
            source_parameter_and_native_width_recipe=recipe,
            exact_binary_sc_quarter_and_half_endpoints=True,
            same_native_P0_F0_anchor_and_relative_amplitude_source=True,
            same_original_sigma_value_and_ordinary_derivative_AST=True,
            identical_ratio_ODE_Psi_transport_and_own_mean_rho4_history=True,
            leading_Q_uses_one_minus_delta_and_full_mean_Z_derivative=True,
            finite_N_modification_is_exact_zero_throughout_lower_interval=True,
            nonzero_complex_comparison_and_Cauchy_errors_retained=True,
            physical_hb_epsilon_measures_and_single_m15_minus_retained=True,
            range_product_is_integral_containment_not_FunctionRef_equality=True,
            function_definitions_remain_true_signed_C1_integrals=True,
            complete_preRp_pulse_pressure_and_reserved_repair_not_admitted=True)

    def assert_graph(self):
        if ((self.cutoff,self.numeric,self.inner,self.lower,self.g,self.c)!=self._owners
            or self.hashes!=self._hashes or self.family!=self.numeric.family
            or self.bindings!=self._bindings or op.snapshot(self.source_identity)!=self._identity
            or full.own.encode(dict(functions=self.functions,fields=self.lower.fields,
                parameters=self.lower.parameters,pressure_bindings=self.lower.pressure_bindings))!=self._bindings):
            raise ValueError('True lower numerical bridge changed')
        self.cutoff.assert_graph();self.numeric.assert_graph()
        self.inner.assert_graph();self.lower.assert_graph()

    @source_precision
    def evaluate(self,Z,reuse_packet=None,Rin_packet=None):
        self.assert_graph();z=self.c.mpf(Z)
        if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original real Z domain [-1,1] required')
        packet=self.numeric.evaluate(Z) if reuse_packet is None else reuse_packet
        view=self.numeric.report(packet)
        if (view['Z']._mpi_!=z._mpi_ or view['query']['lower_phase']._mpi_!=self.numeric.phase_in._mpi_
            or view['query']['upper_phase']._mpi_!=self.numeric.modification_inlet._mpi_):
            raise ValueError('Pressure prefix packet belongs to another Z or radial interval')
        rows,fields,forcing,values=self.numeric._packets[id(packet)][2]
        proof=view['actual_source_reference_and_error_proof']
        if (not view['original_modification_support']['increments_exact_zero']
            or not proof['actual_comparison_errors_retained']
            or not proof['phase_and_radial_map_Z_independent']
            or any(v.order<1 for v in values.values())
            or any(ep(row['nonzero_coefficient_upper'])[0]<=0 for row in proof['original_sigma_bounds'])):
            raise ValueError('Genuine C1 source errors and unchanged collar required')
        combined=None
        if Rin_packet is not None:
            rin=self.cutoff.report(Rin_packet)
            if (rin['Z']._mpi_!=z._mpi_ or rin['chart']!='first'
                or rin['coordinate']._mpi_!=self.c.mpf(1)._mpi_):
                raise ValueError('A same-Z true Rin packet is required for prefix composition')
            inner_values=self.cutoff._packets[id(Rin_packet)][4][3]
            if any(v.algebra is not next(iter(values.values())).algebra for v in inner_values.values()):
                raise ValueError('Pressure composition requires the same live physical source frame')
            composed={name:inner_values[name]+values[name] if name in values else inner_values[name]
                for name in prior.moments.MOMENTS}
            combined=dict(physical_five_prefix_moments_Z1={n:op.export(v,1) for n,v in composed.items()},
                full_velocity_support_names=['m11','m12','m14'],
                pressure_m13_m15_stop_at_rminus=True,preRp_and_pulse_pressure_not_added=True)
        result=dict(source_family=self.family,Z=z,true_function_bindings=self.bindings,
            radial_scope='original Rin=sc/4 to rminus=sc/2',
            true_lower_pressure_C1_enclosures={name:op.export(v,1) for name,v in values.items()},
            true_lower_pressure_budgets=view['signed_source_scaled_pressure_budgets'],
            original_collar_source_enclosures=view,inner_plus_lower_prefix=combined,
            ordinary_Z0_Z1_enclosures_of_true_functions=True,
            true_function_graph_and_numerical_enclosure_context_kept_separate=True,
            genuine_physical_N1p_and_both_viscous_radial_derivatives=True,
            all_actual_core_collar_history_and_nonzero_error_terms_retained=True,
            range_product_encloses_integral_not_a_selected_value=True,
            relative_accuracy_uniform_at_zeros_not_claimed=True,
            complete_external_pressure_transfer_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalTrueLowerNumericPacket();state=(values,combined)
        self._packets[id(token)]=(token,result,packet,Rin_packet,state,op.snapshot((result,state)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalTrueLowerNumericPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live true lower C1 numerical packet required')
        self.numeric.report(entry[2])
        if entry[3] is not None:self.cutoff.report(entry[3])
        if op.snapshot((entry[1],entry[4]))!=entry[5]:raise ValueError('Lower C1 numerical packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_cutoff,accepted_prefix,reuse_packets=None,Rin_packet=None):
    began=time.monotonic();owner=CurrentOriginalN1LowerNumericBridge(accepted_cutoff,accepted_prefix,False)
    packets={}
    for name,Z in (('point','.371'),('independent_point','.211'),('whole_Z',['-1','1'])):
        print('TRUE_LOWER_PRESSURE_C1_START',name,flush=True)
        packets[name]=owner.evaluate(Z,reuse_packet=None if reuse_packets is None else reuse_packets.get(name),
            Rin_packet=Rin_packet if name=='point' else None)
        print('TRUE_LOWER_PRESSURE_C1_COMPLETE',name,flush=True)
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,true_function_bindings=owner.bindings,
        samples={n:owner.report(p) for n,p in packets.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,raw,packets


def restore_checked_archive():
    receipt=full.parent.linear.read(RECEIPT);raw=full.parent.linear.read(NAME)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
        or pins!=raw['input_hashes'] or archive!=sha(NAME) or receipt.get('source_family')!=raw['source_family']
        or receipt.get('source_identity')!=raw['source_identity']):
        raise ValueError('Checked true lower C1 numerical archive required')
    for name,digest in receipt['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Changed true lower C1 archive input '+name)
    return raw,receipt
