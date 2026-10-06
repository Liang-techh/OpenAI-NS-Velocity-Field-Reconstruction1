"""Source-bound global physical T/E atlas cover with conservative unions."""
import ast
from dataclasses import dataclass,asdict
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_current_correlated_tensor_locator import (
    CurrentCorrelatedTensorLocator,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_correlated_radius_operator import SYMBOLS,LOG_HB
from lei_ren_part1_paper_compliant_current_global_tensor_cover_operator import (
    independent_radial_cover_theorem,independent_implicit_cover_theorem,
    directed_candidate_cover_argument,strict_source_margins)

NAME=PREFIX+'current_global_tensor_cover.json'
RECEIPT=PREFIX+'current_global_tensor_cover_check.json'
GATES=('current_original_33_chart_source_radius_ordering_certified','current_global_tensor_physical_cover_certified')
OPEN=tuple(key for key in EARLIER_OPEN if key not in GATES)


@dataclass(frozen=True)
class ExactSharedWidthWitness:
    family:str
    source:str
    datum:str
    uniform_Cstar_family:str
    analytic_core_family:str
    inner_parameter_family:str
    ledger_sha256:str
    exact_width_source:str='exp(log(cstar)-100*log(K))=cstar*K^-100=epsilon_b'


def class_assignment(stem,cls,target,wanted,hashes):
    name=PREFIX+stem+'.py';path=HERE/name;hashes[name]=sha(name)
    tree=ast.parse(path.read_text(encoding='utf8'))
    owner=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
    fn=next(n for n in owner.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    expected=ast.dump(ast.parse(wanted,mode='eval').body)
    rows=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
    if len(rows)!=1 or ast.dump(rows[0])!=expected:raise ValueError('Exact source assignment changed: '+stem+'.'+cls+'.'+target)
    return dict(target=target,exact_source_assignment=wanted,source_sha256=hashes[name])


def shared_width_source(field):
    p=field.locator.physical;b=p.dispatch.provider('bridge_first');w=p.dispatch.provider('switch_first')
    upstream=b.upstream;adapter=w.history.bridge;native=upstream.bridge
    if not (upstream is w.history.upstream and adapter.upstream is upstream and w.switch is w.history.switch
        and w.switch.bridge is adapter and b.bridge is adapter):
        raise ValueError('One original upstream and shared-width source path required')
    ledger=native.records['K1_ledger'];name=PREFIX+'K1_ledger.json'
    if ledger!=json.loads((HERE/name).read_bytes()):raise ValueError('Exact current shared-width ledger changed')
    _verify_hashes(ledger);definition=ledger['parameter_definition'];cstar=read_interval(field.ctx,ledger['cstar'])
    if not (definition['h_b']=='cstar*K^-100' and definition['epsilon_b']=='h_b, the identical scalar'
        and definition['K']=='actual physical norm sum in9.16, bounded by Cstar*Kbar'
        and ledger['h_b_equals_epsilon_b_by_definition'] and ledger['positive_width_not_materialized']
        and endpoints(cstar)[0]>0 and all(mp.isfinite(v) for v in endpoints(cstar))):
        raise ValueError('Original positive formal shared-width definition required')
    if (ledger['uniform_Cstar_family_sha256']!=p.core.records['physical_norm_family']['uniform_Cstar_family_sha256']
        or ledger['base_analytic_core_family_sha256']!=p.core.records['core_transfer']['analytic_core_family_sha256']):
        raise ValueError('Shared exact width does not belong to the checked current core')
    hashes={name:sha(name)};bindings={}
    for stem,cls,target,value in (
      ('inner_bridge_profiles','CompliantInnerBridgeProfiles','self.logh',"read_interval(c,ledger['shared_positive_width_log_enclosure'])"),
      ('inner_switch_profiles','CompliantInnerSwitchProfiles','self.logh','self.bridge.logh'),
      ('actual_bridge_integrals','CompliantActualBridgeIntegrals','self.logh','self.bridge.logh'),
      ('actual_bridge_switch','ActualR100BridgeAdapter','self.logh','c.mpf(upstream.logh)'),
      ('actual_bridge_switch','ActualInnerSwitch','self.logh','bridge.logh'),
      ('actual_switch_mixed_C4','CompliantActualSwitchMixedC4','self.logh','self.switch.logh'),
      ('current_actual_bridge_mixed_C4','CurrentActualBridgeMixedC4','self.logh','self.upstream.logh')):
        bindings[stem+'.'+cls]=class_assignment(stem,cls,target,value,hashes)
    # The formal recipe and original norm/family definition are part of the
    # witness. Equal endpoints merely confirm its numerical enclosure.
    for stem in ('K1_ledger','tolerance_exit_parameters'):
        name=PREFIX+stem+'.py';hashes[name]=sha(name)
        tree=ast.parse((HERE/name).read_text(encoding='utf8'))
        definitions=[kw.value for n in ast.walk(tree) if isinstance(n,ast.Call) for kw in n.keywords if kw.arg=='h_b']
        if not any(isinstance(v,ast.Constant) and v.value=='cstar*K^-100' for v in definitions):
            raise ValueError('Original formal h_b recipe changed: '+stem)
    witness=ExactSharedWidthWitness(field.family,field.source,field.datum_sha,
        ledger['uniform_Cstar_family_sha256'],ledger['base_analytic_core_family_sha256'],
        ledger['admitted_inner_parameter_family_sha256'],hashes[PREFIX+'K1_ledger.json'])
    witness_id=hashlib.sha256(json.dumps(asdict(witness),sort_keys=True).encode()).hexdigest()
    return dict(immutable_exact_shared_width_witness=asdict(witness),witness_sha256=witness_id,
        witness_transport_source_assignments=bindings,
        bridge_and_switch_same_upstream_source_object=True,
        source_path_witnesses=dict(bridge_provider=witness_id,actual_R100_adapter=witness_id,actual_inner_switch=witness_id,microswitch_provider=witness_id),
        original_K_norm_definition_and_same_core_family_bound=True,
        interval_object_identity_or_endpoint_overlap_not_the_exact_source_proof=True,
        exact_source_log_hb='log(cstar)-100*log(K); enclosed by the original shared positive width log ledger',
        input_hashes=hashes,passed=True)


def audit_source_hypotheses(ctx,values,delta,evaluate):
    if any(not mp.isfinite(v) for box in values.values() for v in endpoints(box)):
        raise ValueError('Finite original source parameters and exact log(hb) required')
    dl,dh=endpoints(delta);ml,mh=endpoints(values[SYMBOLS['mu']])
    if not (0<=dl<=dh<1 and 0<ml<=mh<mp.mpf('.25')):
        raise ValueError('Original global cover requires 0<=delta<1 and 0<mu<1/4')
    result={}
    for name,expr in strict_source_margins().items():
        box=evaluate(expr,values)
        if name=='hb':
            # Finite log(hb) proves exp(logh)>0 even when the cap touches 0.
            result[name]=dict(exact_source=str(expr),strictly_positive_by_finite_exact_log=True,enclosure=box)
        else:
            if endpoints(box)[0]<=0:raise ValueError('Strict source radius cover margin unproved: '+name)
            result[name]=dict(exact_source=str(expr),strictly_positive_directed_lower=True,enclosure=box)
    return dict(all_original_33_radial_monotonicity_premises_audited=True,
        strict_positive_source_margins=result,finite_exact_log_hb=values[LOG_HB],source_delta=delta,
        source_mu=values[SYMBOLS['mu']],parameters_are_enclosures_not_selected_midpoints=True,passed=True)


class CurrentGlobalTensorCover:
    @source_precision
    def __init__(self,field=None,require_checked=True):
        self.field=field if field is not None else CurrentCorrelatedTensorLocator()
        if type(self.field) is not CurrentCorrelatedTensorLocator or not self.field.acceptance_loaded:
            raise ValueError('Checked correlated source locator on the complete original graph required')
        self.locator=self.field.locator;self.registry=self.field.registry;self.ctx=self.field.ctx
        self.family=self.field.family;self.source=self.field.source;self.datum_sha=self.field.datum_sha
        self.radial=independent_radial_cover_theorem();self.implicit=independent_implicit_cover_theorem()
        self.candidate=directed_candidate_cover_argument();self.width=shared_width_source(self.field)
        self.hypotheses=self.audit_hypotheses()
        self.hashes={**self.field.hashes,**self.candidate['input_hashes'],**self.width['input_hashes']}
        for stem in ('current_global_tensor_cover_operator','current_global_tensor_cover'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[1]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Global tensor cover admission exceeds source-function scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def audit_hypotheses(self):
        values=self.field._values(self.field.radius('Ra'))
        return audit_source_hypotheses(self.ctx,values,self.locator.physical.delta,self.field._eval)

    def assert_graph(self):
        self.field.assert_graph()
        if not self.field.acceptance_loaded or self.locator is not self.field.locator or self.registry is not self.field.registry:
            raise ValueError('Foreign global physical source graph')
        if encode(pack(self.audit_hypotheses()))!=encode(pack(self.hypotheses)) or shared_width_source(self.field)!=self.width:
            raise ValueError('Current global cover source premises changed')

    def _attach(self,value,request=None):
        self.assert_graph();location=value.get('location',value)
        identity=(location.get('actual_five_defect_family_sha256',location.get('source_family')),
            location.get('implicit_source_sha256'),location.get('datum_enclosure_sha256'))
        if identity!=(self.family,self.source,self.datum_sha):raise ValueError('Foreign global cover family/source/datum')
        if request is None:
            expected=self.locator.locate_log_radius(location['physical_log_radius'],location['physical_z'],
                location['requested_log_tau'],location['theta'],location['physical_viscosity'])
            if encode(pack(location))!=encode(pack(expected)):raise ValueError('Foreign or changed global absolute source location')
        else:self.field._replay(request,location)
        ql,qh=endpoints(location['actual_log_lambda']);zl,zh=endpoints(location['Z'])
        if not (mp.isfinite(ql) and mp.isfinite(qh) and ql<=qh and -1<=zl<=zh<=1
            and location['solver_status'] in self.candidate['accepted_solver_statuses']
            and location['true_finite_point_has_abs_Z_strictly_below_one'] and location['candidates']):
            raise ValueError('Finite directed root bracket and nonempty complete candidate union required')
        return dict(original_physical_source_result=value,
            exact_cover_scope='Every finite physical point at tau>0, constant nu>0, on the fixed admitted original source family; candidate union retains parameter uncertainty.',
            independent_analytic_cover_applies=True,finite_directed_root_bracket_retained=True,
            source_dependent_or_absolute_input_uncertainty_retained=True,
            exact_shared_width_witness_sha256=self.width['witness_sha256'],
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def locate_log_radius(self,*args,**kwargs):return self._attach(self.locator.locate_log_radius(*args,**kwargs))

    @source_precision
    def tensor_log_radius(self,*args,**kwargs):
        location=self.locator.locate_log_radius(*args,**kwargs)
        return self._attach(self.locator.tensor_from_location(location))

    @source_precision
    def cartesian(self,*args,**kwargs):return self._attach(self.locator.cartesian(*args,**kwargs))

    @source_precision
    def correlated_tensor(self,*args,**kwargs):
        request=args[0] if args else kwargs['request']
        return self._attach(self.field.tensor(*args,**kwargs),request=request)

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            independent_original_radial_cover_theorem=self.radial,independent_global_implicit_coordinate_cover_theorem=self.implicit,
            source_bound_directed_candidate_union_cover_argument=self.candidate,exact_shared_width_source_witness=self.width,
            actual_current_source_cover_hypotheses=self.hypotheses,
            actual_current_tensor_regions_available=list(self.registry.registry),
            admitted_adjacent_common_source_traces=list(self.registry.adjacent),admitted_internal_common_source_traces=list(self.registry.internal),
            axis_is_separate_same_core_analytic_extension=True,unbounded_Gamma_is_not_a_finite_radial_truncation=True,
            scope='Independent source-bound global physical chart/T/E cover for every finite physical point at tau>0 and finite constant nu>0, for the admitted fixed construction family. Analytic radius ordering and unique axial inverse; directed candidate unions remain ambiguous. Unique fixed-coordinate seam selection, resolved u/v/w/p, cone/lift, NS, energy, completed flatness and coefficient recursion remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentGlobalTensorCover(require_checked=False)
    field.assert_graph();result=field.manifest()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Independent global physical tensor cover: unique axial inverse, ordered 33 original radial charts, axis/Gamma and conservative unions',flush=True)
    return result


if __name__=='__main__':run()
