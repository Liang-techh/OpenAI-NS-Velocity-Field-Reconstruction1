"""Physical finite-N viscous trace contraction on proved active seams.

The raw repaired state and its independent pulse target are immutable.
Only shared signed traces are eliminated; surviving outer traces keep
positive source budgets and are not assigned numerical centers.
"""
import ast
import copy
import hashlib
import json
import sys
from fractions import Fraction

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX,_bounds
from .original_micro_viscous_m15 import OriginalMicroViscousM15,EndpointFingerprint,contract_trace_chain


PARENT=PREFIX+'compliant_current_original_micro_viscous_m15_numeric.json'
PARENT_SHA256='0903f9100aa46abbfd9b69c667c91b147b364cb0a33cbcc526f4549a53d7a5b8'
CHARTS=('first_micro','second_micro','frozen_macro','first_switch','second_switch','post_power',
        'long_reshape','reference','restoration','postrestore')
VISC_INDEX=dict(zip(CHARTS,(12,16,27,50,68,87,105,123,139,155)))
SEAMS=(('restoration','postrestore'),)


class ActiveEndpointFingerprint(EndpointFingerprint):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.linear_forms={};self.product_forms={}

    def linear(self,value):
        digest,rational=value
        return {None:rational} if rational is not None else self.linear_forms.get(digest,{digest:Fraction(1)})

    def finish_linear(self,terms):
        terms={d:c for d,c in terms.items() if c}
        if not terms or set(terms)=={None}:return self.rational(terms.get(None,0))
        if len(terms)==1:
            digest,coefficient=next(iter(terms.items()))
            if coefficient==1:return digest,None
        rows=sorted(((d or '',c.numerator,c.denominator) for d,c in terms.items()))
        digest=self.digest(['exact_rational_linear_combination',rows])
        self.linear_forms[digest]=terms
        return digest,None

    def combine(self,op,values):
        # Exact additive group arithmetic, including a+((b-a)-8)=b-8.
        # Only rational scalar products distribute over sums. Products of
        # two genuine sums stay opaque; no unbounded polynomial expansion.
        if op=='sum':
            terms={}
            for value in values:
                for d,c in self.linear(value).items():terms[d]=terms.get(d,0)+c
            return self.finish_linear(terms)
        coefficient=Fraction(1);factors=[]
        for digest,rational in values:
            if rational is not None:coefficient*=rational
            else:
                form=self.linear_forms.get(digest)
                if form and len(form)==1 and None not in form:
                    digest,c=next(iter(form.items()));coefficient*=c
                factors.extend(self.product_forms.get(digest,[digest]))
        if not coefficient:return self.zero
        if not factors:return self.rational(coefficient)
        if len(factors)==1:
            return self.finish_linear({d:coefficient*c for d,c in self.linear((factors[0],None)).items()})
        digest=self.digest(['exact_commutative_product',sorted(factors)])
        self.product_forms[digest]=factors
        if coefficient==1:return digest,None
        return self.finish_linear({digest:coefficient})

    @staticmethod
    def name(name):
        if name in {'true_loop_C5_angle_'+n for n in CHARTS}:
            return 'actual_active_common_bound_loop_angle'
        return EndpointFingerprint.name(name)


class OriginalActiveViscousM15(OriginalMicroViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=120):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA256:raise ValueError('Accepted micro trace state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA256
        self.modified=self.archive('modified_Q4_pressure_C3_functions')
        self.loop=self.archive('prescribed_loop_C5_functions')
        self.bridge=self.archive('prescribed_bridge_C6_functions')
        self.outer=self.archive('prescribed_R100_Rm_C6_functions')
        self.restore=self.saved('restore_m15_numeric')
        self.nodes=self.modified['exact_graph_nodes'];self.m=self.modified['modified_Q4_pressure_C3_bindings']
        self.l=self.loop['prescribed_loop_C5_bindings']
        self.background=dict(self.outer['prescribed_R100_Rm_C6_bindings']['charts'])
        bridge=self.bridge['prescribed_bridge_C6_bindings']
        self.background.update(first_micro=bridge['micro'],second_micro=bridge['micro'],frozen_macro=bridge['macro'])
        self.background['post_power']=self.background['post_switch']
        self.background['postrestore']=self.background['post_restore']
        for source in (self.loop,self.bridge,self.outer):
            nodes=source['exact_graph_nodes']
            if self.nodes[:len(nodes)]!=nodes:raise ValueError('Identical original source prefixes required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        self.proofs=[self.prove_seam(*pair) for pair in SEAMS]
        self.inlet_probe=self.prove_seam('reference','restoration',same_native_measure=False,strict=False)
        self.budget_source=self.bind_budget_source()
        name='src/openai_ns_reconstruction/original_active_viscous_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def prove_seam(self,left,right,*,same_native_measure=True,strict=True):
        a,b=[self.l['charts'][n] for n in (left,right)]
        f=ActiveEndpointFingerprint(self.nodes,dict(coordinate=a['coordinate'],loop=self.l['charts']))
        scope={a['coordinate']:f.at(a['domain'][1],{}),b['coordinate']:f.at(b['domain'][0],{})}
        rows=[];unproved=[]
        def compare(group,key,order,ln,rn):
            lh,rh=f.at(ln,scope)[0],f.at(rn,scope)[0]
            if lh!=rh:
                if strict:raise ArithmeticError('Unproved actual seam '+str((left,right,group,key,order,ln,rn)))
                unproved.append(dict(group=group,quantity=key,ordinary_Z_order=order,left_node=ln,right_node=rn,
                    left_endpoint_sha256=lh,right_endpoint_sha256=rh))
                return
            rows.append(dict(group=group,quantity=key,ordinary_Z_order=order,left_node=ln,right_node=rn,canonical_endpoint_sha256=lh))
        for key in ('radius','global_offset','actual_phase'):compare('geometry',key,0,a[key],b[key])
        if same_native_measure:compare('geometry','physical_measure',0,a['physical_measure'],b['physical_measure'])
        radius_dependencies=[]
        for chart in (a,b):
            visited=set()
            def radius_Z_independent(i):
                if i in visited:return
                visited.add(i);row=self.nodes[i];op=row['operation']
                if op=='bound_variable':
                    if i!=chart['coordinate']:raise ValueError('Unexpected physical radius variable')
                elif op in ('original_source_parameter','current_original_source_parameter','exact_positive_integer_expression'):
                    if not row.get('Z_independent'):raise ValueError('Radius parameter must be Z-independent')
                elif op in ('sum','product'):
                    for j in row['arguments']:radius_Z_independent(j)
                elif op=='negative':radius_Z_independent(row['argument'])
                elif op=='positive_quotient':
                    radius_Z_independent(row['numerator']);radius_Z_independent(row['denominator'])
                elif op=='analytic_unary':radius_Z_independent(row['argument'])
                elif op not in ('exact_rational','mathematical_pi'):raise ValueError('Unissued physical radius operation '+op)
            radius_Z_independent(chart['radius'])
            radius_Z_independent(chart['physical_measure'])
            radius_dependencies.append(sorted(visited))
        if not a['total_Z_phase_derivative_exact_zero'] or not b['total_Z_phase_derivative_exact_zero']:
            raise ValueError('Original fixed slow-Z phase required')
        for key in ('U','U_y','m','m_y','Q','Q_y'):
            for k in (0,1):compare('modified',key,k,self.m['charts'][left]['physical'][key][k],self.m['charts'][right]['physical'][key][k])
        for k in (0,1):compare('modified','R_Q_y',k,self.m['charts'][left]['finite_boundary']['R_Q_y'][k],self.m['charts'][right]['finite_boundary']['R_Q_y'][k])
        for name,chart in ((left,a),(right,b)):
            for k in (0,1):
                node=self.nodes[self.m['charts'][name]['finite_boundary']['R_Q_y'][k]]
                if node['operation']!='product' or chart['radius'] not in node['arguments']:
                    raise ValueError('Actual physical R factor required in each RQy row')
        # These original leading rows suffice for the same physical Qy and
        # its ordinary Z1: Qy uses U_y, U-m, and (U-m)_Z.
        for group,key,orders in (('fields','V',3),('fields','m',3),('radial','V_y',2)):
            for k in range(orders):compare('original_leading',key,k,self.background[left][group][key][k],self.background[right][group][key][k])
        if a['physical_measure']!=self.background[left]['physical_measure'] or b['physical_measure']!=self.background[right]['physical_measure']:
            raise ValueError('Physical native Jacobian source required on each side')
        return dict(passed=not unproved,left=left,right=right,left_coordinate=a['coordinate'],right_coordinate=b['coordinate'],
            left_endpoint=a['domain'][1],right_endpoint=b['domain'][0],identity_rows=rows,
            unproved_rows=unproved,unproved_fingerprints_are_not_a_discontinuity_theorem=True,
            empty_integrals_eliminated=sorted(f.empty),global_phase_not_reset=True,
            actual_modified_and_original_leading_traces_equal=not unproved,
            native_Jacobians=[a['physical_measure'],b['physical_measure']],
            common_ordinary_Z_independent_radius=True,radius_and_measure_dependency_nodes=radius_dependencies,
            native_Jacobians_proved_equal=same_native_measure,physical_R_factor_in_both_ordinary_Z_rows=True,
            Qy_recovery='[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2)',
            independent_of_interval_caps=True,actual_signed_centers_evaluated=False)

    def bind_budget_source(self):
        name='src/openai_ns_reconstruction/original_restore_m15.py';raw=(ROOT/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=self.code_hashes[name]:raise ValueError('Pinned positive budget source required')
        tree=ast.parse(raw);cls=next(v for v in tree.body if isinstance(v,ast.ClassDef) and v.name=='OriginalRestoreM15')
        fn=next(v for v in cls.body if isinstance(v,ast.FunctionDef) and v.name=='finite_errors')
        statements=["upper=self.correction_densities(Z,packet,base)",
            "lower=self.correction_densities(Z,lower_packet,lower_base)",
            "result[name]=[self.rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'], None if lower[name][k] is None else lower[name][k]+power*lower_geo['logR'])) for k in (0,1)]"]
        dumps={ast.dump(v) for v in ast.walk(fn)}
        for text in statements:
            if ast.dump(ast.parse(text).body[0]) not in dumps:raise ValueError('Original outer endpoint budget changed')
        return dict(passed=True,source=name,source_sha256=self.code_hashes[name],method='OriginalRestoreM15.finite_errors',
            exact_source_statements=statements,physical_radius_power_for_viscous=1,
            positive_upper_lower_union_not_difference=True,source_cap_not_function_value=True)

    def decode(self,row):
        values=row.get('source_enclosure_exact_mpf_tuples')
        if values is None:values=[row[k+'_exact_mpf_tuple'] for k in ('lower','upper')]
        return self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in values])
    def record(self,value):
        lo,hi=_bounds(self.ctx,value)
        return dict(lower_approximation=self.ctx.nstr(lo,45),upper_approximation=self.ctx.nstr(hi,45),
            source_enclosure_exact_mpf_tuples=[list(v) for v in value._mpi_],
            positive_error_bound_atom_only=True,source_value_selected=False)
    def positive_union(self,values):
        values=[v for v in values if v is not None]
        if not values:raise ValueError('Nonzero remaining endpoint budget required')
        top=max(_bounds(self.ctx,v)[1] for v in values)
        return self.iv.mpf(top)+self.iv.log(len(values))

    def largest_caps(self,terms,target):
        # Normalize by the same independent physical m15 unit, never by P1.
        unit=sum((self.decode(v) for v in target['log_scale_terms']),self.iv.mpf(0))
        logs=[sum((self.decode(v) for v in t['log_scale_terms']),self.iv.mpf(0))-unit for t in terms]
        result=[]
        for k in (0,1):
            candidates=[]
            for t,log in zip(terms,logs):
                magnitude=_bounds(self.ctx,abs(self.decode(t['coefficient_C0_C1'][k])))[1]
                if magnitude:
                    candidates.append((_bounds(self.ctx,log+self.iv.log(self.iv.mpf(magnitude)))[1],t['name']))
            value,name=max(candidates)
            result.append(dict(ordinary_Z_order=k,name=name,
                normalized_log_cap_upper_approximation=self.ctx.nstr(value,40),
                normalized_log_cap_upper_exact_mpf_tuple=list(value._mpf_),
                cap_is_not_an_actual_residual_value=True))
        return result

    def evaluate(self,Z):
        parent=next((p for p in self.parent['samples'] if p['Z']==Z),None)
        if parent is None:raise ValueError('Saved axial query required; no ancestor replay')
        raw=parent['raw_m15_C0_C1_terms'];terms=parent['correlated_reduced_m15_C0_C1_terms']
        r=next(p for p in self.restore['samples'] if p['Z']==Z and p['chart']=='restoration' and p['coordinate']['value']=='1')
        p=next(p for p in self.restore['samples'] if p['Z']==Z and p['chart']=='postrestore' and p['coordinate']['value']=='1')
        if raw[:156]!=p['m15_C0_C1_terms'] or raw[:140]!=r['m15_C0_C1_terms']:
            raise ValueError('Exact original restoration/postrestore states required')
        if raw[139]['name']!='restoration_finite_N_viscous_error' or raw[155]['name']!='postrestore_finite_N_viscous_error':
            raise ValueError('Original signed viscous sectors required')
        # The saved lower packet is reference native1, not independently
        # certified as restoration native0. Do not reinterpret it as a new
        # surviving-endpoint bound. Keep BOTH entire accepted aggregates.
        old_lower,old_upper=raw[139],raw[155]
        for term in (old_lower,old_upper):
            if len(term['log_scale_terms'])!=1 or not term['log_scale_terms'][0]['positive_error_bound_atom_only']:
                raise ValueError('Original physical endpoint aggregates required')
        bound=self.positive_union([self.decode(t['log_scale_terms'][0]) for t in (old_lower,old_upper)])
        replacement=copy.deepcopy(old_upper);replacement['name']='restoration_postrestore_finite_N_viscous_external_endpoints'
        replacement['log_scale_terms']=[self.record(bound)]
        replacement['actual_function_binding']=dict(error_enclosure_only=True,
            function='delta(RQy)(restoration native0)-delta(RQy)(postrestore native1)',
            delta_definition='delta(RQy)=R*(Qy_modified-Qy_original_leading)',
            internal_seam='restoration native1 = postrestore native0',
            retained_external_endpoints=['restoration native0','postrestore native1 = Rm'],
            actual_internal_trace_proved_equal=True,both_external_traces_retained_without_zero_assumption=True,
            positive_budget_source=self.budget_source,ordinary_Z1=True,no_error_cap_subtraction=True,
            original_both_full_positive_aggregates_retained=True,
            separate_restoration_inlet_bound_claimed=False,
            inherited_chart_budget_validity_not_independently_recertified=True)
        reduced=[]
        for t in terms:
            if t['name']==raw[139]['name']:continue
            reduced.append(copy.deepcopy(replacement if t['name']==raw[155]['name'] else t))
        chain=contract_trace_chain([dict(name='restoration',lower='restoration_native0',upper='restoration_end'),
            dict(name='postrestore',lower='postrestore_start',upper='Rm')],{'postrestore_start':'restoration_end'})
        return dict(Z=Z,source_family=self.family,raw_m15_C0_C1_terms=copy.deepcopy(raw),
            correlated_reduced_m15_C0_C1_terms=reduced,original_raw_terms=213,previous_reduced_terms=200,
            correlated_reduced_terms=len(reduced),source_seam_proofs=self.proofs,signed_trace_chain=chain,
            external_inlet_trace_probe=self.inlet_probe,
            largest_remaining_individual_log_caps=self.largest_caps(reduced,parent['independent_target']),
            old_budget_terms_retained_as_provenance=[copy.deepcopy(raw[i]) for i in (139,155)],
            external_endpoint_budget_witness=dict(source=self.budget_source,
                restoration_positive_aggregate_retained=copy.deepcopy(old_lower['log_scale_terms'][0]),
                postrestore_positive_aggregate_retained=copy.deepcopy(old_upper['log_scale_terms'][0]),
                original_reference_lower_packet_provenance_only=copy.deepcopy(r['actual_lower_finite_N_endpoint_budget']),
                reference_lower_packet_not_reinterpreted_as_restoration_inlet=True,
                strategy='positive union of both full original sectors 139 and 155; inherited chart budgets only',
                no_new_separate_external_endpoint_bound=True,positive_union=self.record(bound)),
            replacement_outer_endpoints=replacement,independent_target=copy.deepcopy(parent['independent_target']),
            contract=dict(raw_213_verbatim=True,only_proved_viscous_seam_contracted=True,
                modified_and_leading_traces_both_verified=True,both_external_traces_retained=True,
                original_both_full_positive_aggregates_retained=True,
                reference_inlet_reinterpretation_not_used=True,
                ordinary_Z1_and_physical_R_once=True,all_other_source_volume_width_errors_retained=True,
                independent_target_unchanged=True,ancestors_not_reexecuted=True),
            limits=dict(restoration_postrestore_viscous_internal_trace_resolved=True,
                reference_restoration_finite_N_inlet_trace_resolved=False,
                inherited_external_budget_source_identity_independently_certified=False,
                remaining_external_trace_signed_values_resolved=False,independent_numeric_Rp_seam=False,
                global_five_moment_numeric_join=False,exact_heat_exterior=False,stress_cone=False,
                flat_remainder=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
