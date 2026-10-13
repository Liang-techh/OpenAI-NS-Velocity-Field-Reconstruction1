"""Same-source numerical enclosures of the true unrepaired common-inner n1.

The accepted true Volterra graph defines the functions. The separately
accepted regular-series owner only encloses them after source, physical
conjugation, complete-tail and exact first-chart guards. No center becomes
a FunctionRef and no n1 repair or global recursion gate is closed here.
"""
import ast
import copy
from dataclasses import dataclass
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as full
import lei_ren_part1_paper_compliant_current_original_n1_regular_series_function as regular
import lei_ren_part1_paper_compliant_current_original_core_hierarchy_source as hierarchy
import lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions as true_core
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE, PREFIX, sha = full.HERE, full.PREFIX, full.sha
NAME = PREFIX+'current_original_n1_inner_numeric_bridge.json.gz'
RECEIPT = PREFIX+'current_original_n1_inner_numeric_bridge_check.json'
GATE = 'current_original_true_unrepaired_n1_common_inner_numeric_enclosure_bridge_installed'
OPEN = tuple(dict.fromkeys(full.OPEN+(
    'n1_repair_velocity_C3_functions_installed', 'n1_repaired_own_M_V_P_installed',
    'positive_order_five_moment_repair_closed', 'uniform_full_pressure_derivative_accuracy_certified',
    'temporal_recursion', 'full_corrected_Navier_Stokes_solution')))


def callable_AST(fn):
    return hashlib.sha256(ast.dump(ast.parse(inspect.getsource(fn).strip())).encode()).hexdigest()


def original_normalization_and_source_recipes():
    """Pin exact producing recipes; differently rounded boxes are not identities."""
    fn = hierarchy.original.OriginalWholeZBridgeSource.__init__
    tree = ast.parse(inspect.getsource(fn).strip())
    assignments = {ast.unparse(target): ast.unparse(node.value)
        for node in ast.walk(tree) if isinstance(node, ast.Assign) for target in node.targets}
    expected = {'q': "scalar(major, 'scaled_map_Lipschitz_upper')",
        'core.sigma': 'core.j / 500',
        'self.logh': "scalar(original, 'source_log_hb_enclosure')"}
    if any(assignments.get(k) != v for k,v in expected.items()):
        raise ValueError('Original parameter/error/width recipes changed')
    core = next(node.value for node in ast.walk(tree) if isinstance(node, ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='core' for t in node.targets))
    fields = {kw.arg: ast.unparse(kw.value) for kw in core.keywords}
    if fields['correction'] != "scalar(major, 'scaled_map_size_upper') / (1 - q)":
        raise ValueError('Legacy enclosure no longer uses the actual-map correction norm')
    expected_fields = {'j': "scalar(major, 'required_j')", 'h': "scalar(tube, 'Xh_parameter')",
        'logC': "scalar(norm, 'selected_logCstar')", 'datum': 'admission.pressure'}
    if any(fields.get(k)!=v for k,v in expected_fields.items()):
        raise ValueError('Original core source leaves changed')
    return dict(passed=True,original_parameter_error_and_width_callable_AST_sha256=callable_AST(fn),
        original_correction_recipe='scaled_map_size_upper/(1-scaled_map_Lipschitz_upper)',
        original_sigma_recipe='j/500',original_positive_width_log_source='source_log_hb_enclosure',
        directed_rounding_of_the_same_recipe_not_required_to_be_identical=True)


@dataclass(frozen=True,eq=False)
class OriginalTrueN1NumericPacket:
    scope: str = 'true canonical unrepaired common-inner n1'


class CurrentOriginalN1InnerNumericBridge:
    @source_precision
    def __init__(self, function_owner, numerical_owner, hierarchy_owner, require_checked=True):
        if type(function_owner) is not full.CurrentOriginalN1FullDefectLinearFunctions or not function_owner.acceptance_loaded:
            raise ValueError('Checked true full defect/function graph required')
        if type(numerical_owner) is not regular.CurrentOriginalN1RegularSeriesFunction or not numerical_owner.acceptance_loaded:
            raise ValueError('Checked original complete regular-series numerical owner required')
        if type(hierarchy_owner) is not hierarchy.CurrentOriginalCoreHierarchySource or not hierarchy_owner.acceptance_loaded:
            raise ValueError('Checked actual numerical hierarchy/source owner required')
        self.full, self.numeric, self.hierarchy = function_owner, numerical_owner, hierarchy_owner
        self.inner = function_owner.cutoff.inner
        self.collar = self.inner.collar
        self.core = self.collar.core
        self.g, self.family = function_owner.g, copy.deepcopy(function_owner.family)
        self.functions = self.inner.functions
        self.recovery = self.inner.recovered
        self.c = numerical_owner.c
        self.common = self.core._oracle or self.core.attach_oracle()
        self.selection = self.numeric.mean.local.inner
        self.hashes = {}
        for owner in (self.full,self.numeric,self.hierarchy,self.common):
            for name,digest in owner.hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:
                    raise ValueError('Conflicting function/enclosure source pin '+name)
                self.hashes[name]=digest
        for module in (full,regular,hierarchy,true_core,true_core.common):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):
                digest=sha(name)
                if name in self.hashes and self.hashes[name]!=digest:
                    raise ValueError('Conflicting accepted module pin '+name)
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        checker=Path(__file__).stem+'_check.py'
        self.hashes[checker]=sha(checker)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed numerical bridge input '+name)
        self._owners=(self.full,self.numeric,self.hierarchy,self.inner,self.collar,self.core,self.common,self.selection,self.c)
        self._function_bindings=full.own.encode(dict(functions=self.functions,recovery=self.recovery))
        self._common_bindings=regular.op.snapshot((self.common.model_bindings,self.common.consumer_bindings,
            self.common.measure_bindings,self.common.pressure_bindings))
        self.source_identity=self.assert_source_identity()
        self._source_identity_snapshot=regular.op.snapshot(self.source_identity)
        self._hashes=copy.deepcopy(self.hashes)
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=full.parent.linear.read(RECEIPT)
            inputs=dict(receipt.get('input_hashes',{}))
            archive_sha=inputs.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE)
                or receipt.get('source_family')!=self.family or any(receipt.get(k) for k in OPEN)
                or inputs!=self.hashes or archive_sha!=sha(NAME)):
                raise ValueError('Checked true n1 numerical bridge receipt required')
            self.acceptance_loaded=True

    @source_precision
    def assert_source_identity(self):
        self.full.assert_graph();self.numeric.assert_graph();self.hierarchy.assert_graph()
        if (self.numeric.family!=self.family or self.hierarchy.family!=self.family
            or self.numeric.physical is not self.hierarchy.core
            or self.numeric.weighted.source is not self.hierarchy.source
            or self.numeric.mean.inputs.source is not self.hierarchy.source
            or self.core._oracle is not self.common or not self.common.acceptance_loaded):
            raise ValueError('True graph and numerical enclosure source chain differs')
        actual=self.core.parameter_identity(self.common)
        source_graph=self.common.source.current_source_graph()
        if not all(source_graph.values()) or not all(self.common.consumer_bindings.values()):
            raise ValueError('Actual nonlinear common source graph/consumers differ')
        if not self.common.model_bindings['passed'] or not self.common.measure_bindings['passed']:
            raise ValueError('Actual model/tube/pressure measure proof required')
        old=self.hierarchy.source;new=self.common.core
        # These exact saved producer records define the same scalar sources,
        # models, domains and norm premises. Their interval contexts can
        # recompute sigma/correction with different outward rounding.
        for key in ('core_transfer','core_uniform_bounds','physical_norm_family','shared_analytic_tube'):
            if old.records[key]!=new.records[key]:
                raise ValueError('Numerical and actual core producer records differ: '+key)
        if not all(self.hierarchy.parameter_checks.values()):
            raise ValueError('Numerical scalar recipes do not enclose the original exact parameters')
        recipe=original_normalization_and_source_recipes()
        read=lambda row,key:hierarchy.original.read(self.c,row[key])
        major=old.records['core_transfer']
        correction=read(major,'scaled_map_size_upper')/(1-read(major,'scaled_map_Lipschitz_upper'))
        if self.numeric.physical.correction._mpi_!=correction._mpi_ or old.core.sigma._mpi_!=(old.core.j/500)._mpi_:
            raise ValueError('Numeric source changed its directed actual-map correction/sigma recipes')
        amplitude=old.amplitude;roots=amplitude.roots
        restored_roots=hierarchy.original.restore(self.c,old.records['anchored_axis_amplitude']['root_certificates'])
        root_flags=('unique_real_anchor_certified','uniform_interval_parameter_root_certificates',
            'nominal_parameters_used_only_for_provisional_centers','all_six_poles_certified')
        anchor_lo,anchor_hi=regular.ep(roots['anchor_interval'])
        if (not all(roots[k] for k in root_flags)
            or not -regular.ep(old.core.j)[1]<=anchor_lo<=anchor_hi<=0
            or regular.ep(roots['anchor_Hprime_lower'])[0]<=0
            or amplitude.core is not old.core or amplitude.rebuild is not old.rebuild
            or regular.op.snapshot(roots)!=regular.op.snapshot(restored_roots)
            or amplitude.gradient_jets.__func__ is not self.common.source.amplitude.gradient_jets.__func__
            or any(amplitude.params[k]._mpi_!=getattr(old.core,k)._mpi_ for k in ('j','delta','sigma'))
            or amplitude.params['a']._mpi_!=((9-old.core.delta)/2)._mpi_):
            raise ValueError('Numeric F0 must use the same unique original H root and fixed scalar recipes')
        physical_AST=callable_AST(hierarchy.core_module.CompliantCorePhysicalField.normalized_jets)
        if self.numeric.physical.normalized_jets.__func__ is not hierarchy.core_module.CompliantCorePhysicalField.normalized_jets:
            raise ValueError('Numerical actual-model-plus-error core callback changed')
        binding=self.hierarchy.P0_binding
        flag_keys=('same_fourteen_continuous_atoms_and_analytic_flatten_function',
            'current_order6_first_six_and_native_order5_same_defining_rows',
            'common_P0_kept_independent_of_cumulative_pressure')
        if (not binding['passed'] or not all(binding[k] for k in flag_keys)
            or callable_AST(old.core.datum.normalized_jets.__func__)!=binding['inherited_normalized_jets_AST_sha256']
            or callable_AST(new.datum.normalized_jets.__func__)!=binding['inherited_normalized_jets_AST_sha256']):
            raise ValueError('Original independent fourteen-atom P0 function identity differs')
        if (list(self.selection.sc._mpi_[0])!=self.collar.sc_bits
            or self.selection.sc._mpi_[0]!=self.selection.sc._mpi_[1]
            or self.selection.phase_in._mpi_!=(self.selection.sc/4)._mpi_):
            raise ValueError('True physical first chart and numerical common-inner endpoint differ')
        hb=self.g.nodes[self.collar.parameters['hb'].node]
        if (hb['operation']!='current_original_source_parameter' or hb['name']!='hbB'
            or hb['source_family']!=self.family
            or not hb['quantity_is_exact_original_parameter_not_range_endpoint']
            or self.collar.geometry_certificate['original_width_definition']!=
                'hb=epsilon_b=cstar*K^-100; K is the original global physical norm sum'):
            raise ValueError('Numerical collar no longer encloses the original exact positive hb')
        regular_receipt=full.parent.linear.read(regular.RECEIPT)
        history_flags=('genuine_core_actual_first_inlet_function_bindings_preserved',
            'each_exact_V_word_uses_one_zero_to_x_path_and_inherits_source_seam')
        if not all(regular_receipt['checks'].get(k) for k in history_flags):
            raise ValueError('Actual first-collar own-history and continuous Volterra path proof required')
        if (not self.inner.tail_admission['all_k_derivative_words_allowed']
            or not self.numeric.contract()['all_remaining_words_enclosed_not_set_zero']):
            raise ValueError('Complete original feedback tail required')
        return dict(passed=True,true_graph_parameters_and_anchor=actual,
            exact_saved_core_parameter_model_and_norm_records_identical=True,
            numerical_actual_core_callback_AST_sha256=physical_AST,
            independent_P0_callable_AST_sha256=binding['inherited_normalized_jets_AST_sha256'],
            same_fourteen_atom_implicit_P0_function_not_hydrated_center=True,
            original_source_recipe_binding=recipe,
            common_unique_twenty_term_fixed_point_and_actual_tail_consumer_proof_bound=True,
            numerical_F0_uses_same_unique_H_root_not_provisional_center=True,
            legacy_full_root_certificate_matches_its_exact_saved_producer=True,
            equal_anchor_by_same_H_polynomial_and_certified_unique_root_not_box_overlap=True,
            directed_model_plus_error_callbacks_enclose_that_same_actual_solution=True,
            same_prescribed_collar_ratio_shear_and_own_integral_history=True,
            exact_first_fraction_to_phase='s=q*sc/4; core coordinate is rho',
            original_positive_hb_retained_in_log_frame=True,
            common_oracle_alone_not_used_as_first_collar_history=True,
            function_definitions_remain_true_graph_FunctionRefs=True,
            numerical_centers_and_boxes_not_imported_as_function_nodes=True,
            box_overlap_not_used_as_source_identity=True)

    def assert_graph(self):
        if ((self.full,self.numeric,self.hierarchy,self.inner,self.collar,self.core,self.common,self.selection,self.c)!=self._owners
            or self.hashes!=self._hashes
            or self.family!=self.full.family or self.family!=self.numeric.family
            or regular.op.snapshot(self.source_identity)!=self._source_identity_snapshot
            or self._function_bindings!=full.own.encode(dict(functions=self.functions,recovery=self.recovery))
            or self._common_bindings!=regular.op.snapshot((self.common.model_bindings,self.common.consumer_bindings,
                self.common.measure_bindings,self.common.pressure_bindings))):
            raise ValueError('True n1 numerical bridge owner/source state changed')
        self.full.assert_graph();self.numeric.assert_graph();self.hierarchy.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',reuse_packet=None):
        self.assert_graph();c=self.c;value=c.mpf(coordinate);lo,hi=regular.ep(value)
        if chart=='core':
            if lo<0 or hi>4:raise ValueError('Core rho must be a subset of [0,4]')
            phase=None
        elif chart=='first':
            if lo<0 or hi>1:raise ValueError('First coordinate is a fraction q in [0,1], not physical s')
            phase=value*self.selection.sc/4
        else:raise ValueError('Only the original common-inner core/first charts are available')
        packet=self.numeric.evaluate(Z,coordinate,chart) if reuse_packet is None else reuse_packet
        view=self.numeric.report(packet)
        if (view['source_family']!=self.family or not view['same_original_common_inner_function_for_every_query']
            or view['Z']._mpi_!=c.mpf(Z)._mpi_ or view['query']['chart']!=chart
            or view['query']['coordinate']._mpi_!=value._mpi_):
            raise ValueError('Numerical query escaped the admitted true common-inner family')
        result=dict(source_family=self.family,chart=chart,
            numerical_coordinate=value,original_true_phase_enclosure=phase,
            exact_coordinate_identity='rho=coordinate' if chart=='core' else 's=coordinate*sc/4; rho=4*exp(hb*s)',
            defining_true_function_nodes=self._function_bindings,
            numerical_enclosures=view,
            all_six_physical_function_ordinary_Z0_to_Z3_enclosed=True,
            supported_first_radial_mixed_Z0_to_Z2_enclosed=True,
            odd_last_two_rho_derivatives_at_axis_not_reported=view['odd_last_two_rho_derivatives_at_axis_not_reported'],
            actual_original_query_executed=True,
            retained_all_word_tail_and_nonzero_source_errors=True,
            physical_outputs_keep_unmaterializable_amplitudes_in_factored_log_frames=True,
            scoped_to_unrepaired_common_inner=True,
            repaired_five_moment_or_global_accuracy_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalTrueN1NumericPacket()
        self._packets[id(token)]=(token,result,packet,regular.op.snapshot(result))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalTrueN1NumericPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live true n1 numerical bridge query packet required')
        self.numeric.report(entry[2])
        if regular.op.snapshot(entry[1])!=entry[3]:raise ValueError('Changed numerical bridge packet')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(function_owner,numerical_owner,hierarchy_owner,reuse_packets=None):
    began=time.monotonic()
    owner=CurrentOriginalN1InnerNumericBridge(function_owner,numerical_owner,hierarchy_owner,require_checked=False)
    queries=[];packets=[]
    for i,(Z,coordinate,chart) in enumerate((('.125','0','core'),('.371','2','core'),('.371','.625','first'))):
        print('TRUE_N1_NUMERIC_QUERY_START',chart,Z,coordinate,flush=True)
        packet=owner.evaluate(Z,coordinate,chart,reuse_packet=None if reuse_packets is None else reuse_packets[i]);packets.append(packet)
        queries.append(owner.report(packet))
        print('TRUE_N1_NUMERIC_QUERY_COMPLETE',chart,Z,coordinate,flush=True)
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,
        true_function_bindings=owner._function_bindings,
        full_true_graph_node_count=len(owner.g.nodes),
        numerical_queries=queries,input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(regular.op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,raw,packets


def restore_checked_archive():
    receipt=full.parent.linear.read(RECEIPT);raw=full.parent.linear.read(NAME)
    if not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN):
        raise ValueError('Checked original true n1 numerical bridge archive required')
    if receipt.get('source_family')!=raw.get('source_family'):
        raise ValueError('Original numeric bridge archive source family differs')
    inputs=dict(receipt['input_hashes']);inputs.pop(NAME,None)
    if raw.get('input_hashes')!=inputs:
        raise ValueError('Original numeric archive and accepted input pins differ')
    for name,digest in receipt['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Changed numeric bridge archive input '+name)
    if not receipt['physical_conjugation_proof']['passed'] or not raw['source_identity']['passed']:
        raise ValueError('Exact same-source conjugation proof required')
    if raw['source_identity']!=receipt['source_identity']:
        raise ValueError('Saved archive and checked source identities differ')
    return raw,receipt
