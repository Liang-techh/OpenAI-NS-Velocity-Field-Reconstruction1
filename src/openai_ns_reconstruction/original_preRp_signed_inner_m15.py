"""Original preRp m15 with the computed signed inner/lower cumulative inlet.

Replace both legacy inner and lower magnitude components by the complete
direct source result at rminus. Its physical log sectors stay separate;
the remaining actual active/compact/quiet source terms are unchanged.
"""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX
from .original_preRp_m15 import OriginalPreRpM15Numeric


INNER_TERMS=('core_axial_coupling','core_pressure_integral_and_finite_boundary',
    'first_lower_axial_width','first_lower_regular_pressure_width',
    'first_lower_singular_pressure','first_lower_singular_width_error')
OUTER_TERMS=('active','compact_finite_N','compact_exact_leading','own_swirl')


class OriginalPreRpSignedInnerM15Numeric:
    def __init__(self, *, preRp=None, **kwargs):
        if preRp is not None and kwargs:
            raise ValueError('Supply an existing preRp owner or constructor options')
        self.pre=preRp if preRp is not None else OriginalPreRpM15Numeric(**kwargs)
        self.ctx,self.iv,self.family=self.pre.ctx,self.pre.iv,self.pre.family
        self.directory=self.pre.directory
        self.hashes,self.code_hashes=dict(self.pre.hashes),dict(self.pre.code_hashes)
        self._engine,self._packets=None,{}
        path=self.directory/(PREFIX+'compliant_current_original_n1_inner_m15_numeric.json')
        raw=path.read_bytes();data=json.loads(raw)
        if (data['source_family']!=self.family or not data['checks']['passed']
                or not data['limits']['signed_core_first_lower_m15_C0_C1']
                or data['limits']['independent_preRp_signed_numeric_seam']):
            raise ValueError('Accepted same-source signed core/first/lower m15 evidence required')
        required=data['checks']['independent_actual_source']
        if (required['total']!=6 or not required['entire_one_minus_h_in_singular_coefficient']
                or data['checks']['velocity_stop_and_nonzero_pressure_history_conditions']<12):
            raise ValueError('Original pressure source, physical scales and retained lower history required')
        self._merge_pins(data['input_hashes'],self.directory,self.hashes)
        self._merge_pins(data['implementation_hashes'],ROOT,self.code_hashes)
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self.saved={Fraction(p['Z']):p for p in data['samples']
            if Fraction(p['q'])==2 and p['cells_per_piece']==16}
        self._report_keys=set(self.saved)
        selection_name=PREFIX+'compliant_current_inner_exit_strict_collar.json'
        p=self.directory/selection_name;selection_raw=p.read_bytes()
        if self.hashes.get(selection_name)!=hashlib.sha256(selection_raw).hexdigest():
            raise ValueError('Original exact common-inner phase source must be pinned')
        selection=json.loads(selection_raw)
        endpoints=selection['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint']
        if endpoints['lower_exact_mpf_tuple']!=endpoints['upper_exact_mpf_tuple']:
            raise ValueError('Original selected sc singleton required')
        self.sc=endpoints['lower_exact_mpf_tuple']
        self.phase=list(self.sc);self.phase[2]-=2
        self.source=dict(direct_numeric_report=path.name,
            direct_numeric_method='OriginalN1InnerM15Numeric.evaluate(Z,q=2)',
            original_source='signed core Volterra plus actual cutoff/lower pressure derivative',
            old_inner_scope='axis toRin, phase sc/4',
            old_lower_scope='Rin to rminus, phase sc/4 tosc/2',
            replacement_scope='axis through rminus, actual q=2 and rho=4exp(hb*sc/2)',
            archived_FunctionRef_reuse_claimed=False,
            source_agreement_not_inferred_from_enclosure_overlap=True,
            preserved_core_and_collar_physical_units=True)
        name='src/openai_ns_reconstruction/original_preRp_signed_inner_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    @staticmethod
    def _merge_pins(incoming,root,target):
        for name,digest in incoming.items():
            if target.get(name,digest)!=digest:
                raise ValueError('Conflicting actual signed-inner/preRp source: '+name)
            if name not in target and hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed actual signed-inner/preRp source: '+name)
            target[name]=digest

    def inner_lower(self,Z):
        key=Fraction(Z)
        if key in self.saved:return self.saved[key],key in self._report_keys
        if self._engine is None:
            from .original_n1_inner_m15_numeric import OriginalN1InnerM15Numeric
            self._engine=OriginalN1InnerM15Numeric(cells=16,precision=80,directory=self.directory)
            if self._engine.family!=self.family:raise ValueError('Original signed inner source family changed')
            self._merge_pins(self._engine.hashes,self.directory,self.hashes)
            self._merge_pins(self._engine.code_hashes,ROOT,self.code_hashes)
        packet=self._engine.evaluate(Z,'2')
        self.saved[key]=packet
        return packet,False

    def evaluate(self,Z,chart,coordinate):
        # The original external owner enforces the full native query domain.
        legacy=self.pre.evaluate(Z,chart,coordinate)
        key=Z,chart,coordinate
        if key in self._packets:return self._packets[key]
        inner,saved=self.inner_lower(Z)
        if (inner['source_family']!=self.family or Fraction(inner['Z'])!=Fraction(Z)
                or Fraction(inner['q'])!=2 or inner['chart']!='selected_lower_collar'
                or inner['exact_sc_tuple']!=self.sc or inner['exact_phase_tuple']!=self.phase
                or tuple(t['name'] for t in inner['m15_C0_C1_terms'])!=INNER_TERMS):
            raise ValueError('Complete same-Z signed core/first/lower rminus inlet required')
        terms=legacy['m15_C0_C1_terms']
        expected=('inner','lower')+OUTER_TERMS+(('quiet_exact_leading',) if chart=='quiet' else ())
        if tuple(t['name'] for t in terms)!=expected:
            raise ValueError('Original source partition must contain inner/lower exactly once')
        out=copy.deepcopy(legacy)
        out['m15_C0_C1_terms']=copy.deepcopy(inner['m15_C0_C1_terms']+terms[2:])
        out['signed_inner_lower_source']=dict(self.source,used_saved_numeric_packet=saved,
            exact_sc_tuple=copy.deepcopy(self.sc),exact_phase_tuple=copy.deepcopy(self.phase),
            actual_Z=Z,native_q='2')
        out['source_partition_identity_reference']=out.pop('source_partition_identity')
        out['signed_inner_lower_canonical_function_reference']=dict(
            incoming_function=out['incoming_true_function'],compact_cumulative_function=out['compact_cumulative_true_function'],
            direct_numerical_source_registered_as_archived_graph_nodes=False,
            reference_does_not_establish_independent_signed_numeric_cancellation=True)
        out['contract'].update(signed_inner_lower_actual_cumulative_consumed_once=True,
            both_legacy_magnitude_components_removed=True,
            exact_selected_sc_and_rminus_coordinate_retained=True,
            original_inner_log_sectors_not_renormalized_again=True,
            active_compact_own_quiet_terms_preserved_unchanged=True,
            actual_direct_source_provenance_distinct_from_archived_graph_reference=True)
        out['limits'].update(signed_inner_lower_C0_C1_approximation=True,
            inherited_inner_lower_symmetric_proxies_replaced=True,
            signed_active_and_compact_finite_N_values=False,
            resolved_signed_point_centers=False,independent_signed_numeric_seam_residual=False,
            relative_cancellation_accuracy=False,global_five_moment_numeric_join=False)
        self._packets[key]=out
        return out
