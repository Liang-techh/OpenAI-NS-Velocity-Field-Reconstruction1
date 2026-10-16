"""Source-bound finite-N viscous telescoping on the two micro windows.

The signed traces are functions, not their absolute error budgets. Only a
proved common physical trace can be contracted. Raw accepted terms remain
available; all other numerical, volume and width errors are retained.
"""
import copy
import ast
from fractions import Fraction
import gzip
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT, DIRECTORY, PREFIX, _bounds


REPAIRED = PREFIX+'compliant_current_original_repaired_Rp_m15_numeric.json'
REPAIRED_SHA256 = 'bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871'


def contract_trace_chain(pieces, aliases):
    """Contract lower-minus-upper source IDs; numerical bounds are ignored."""
    def resolve(name):
        seen=set()
        while name in aliases:
            if name in seen:raise ValueError('Cyclic trace alias')
            seen.add(name);name=aliases[name]
        return name
    coefficients={};visits={}
    for piece in pieces:
        for side,sign in (('lower',1),('upper',-1)):
            name=resolve(piece[side]);coefficients[name]=coefficients.get(name,0)+sign
            visits[name]=visits.get(name,0)+1
    return dict(signed_coefficients={k:v for k,v in coefficients.items() if v},
        cancelled_interfaces=[k for k,v in coefficients.items() if not v and visits[k]>1])


class EndpointFingerprint:
    """Exact endpoint substitution and alpha-equivalence on the seam collar.

    Rationals, empty integrals and exp(0) are simplified exactly. The only
    removed lazy predicate is the original first-only collar s<=3sc/4,
    proved false at s=1 by the pinned exact selector sc<1/4.
    """
    def __init__(self, nodes, bindings):
        self.nodes,self.b=nodes,bindings;self.cache={};self.empty=set();self.collar=set()
        self.zero=self.rational(0);self.one=self.rational(1)
        self.variables={row['name']:i for i,row in enumerate(nodes) if row['operation']=='bound_variable'}

    @staticmethod
    def digest(body):return hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    def rational(self,value):
        value=Fraction(value)
        return self.digest(['rational',value.numerator,value.denominator]),value
    @staticmethod
    def name(name):
        return 'actual_common_loop_angle' if name in ('true_loop_C5_angle_first_micro','true_loop_C5_angle_second_micro') else name
    def combine(self,op,values):
        rational=Fraction(0 if op=='sum' else 1);other=[]
        for digest,value in values:
            if value is None:other.append(digest)
            elif op=='sum':rational+=value
            else:rational*=value
        if op=='product' and not rational:return self.zero
        if not other:return self.rational(rational)
        if rational != (0 if op=='sum' else 1):other.append(self.rational(rational)[0])
        if len(other)==1:return other[0],None
        return self.digest([op,sorted(other)]),None
    def at(self,node,scope=None):
        scope=scope if scope is not None else {self.b['coordinate']:self.one}
        key=node,tuple(sorted((k,v[0]) for k,v in scope.items()))
        if node in scope:return scope[node]
        if key in self.cache:return self.cache[key]
        row=self.nodes[node];op=row['operation'];F=lambda i:self.at(i,scope)
        if op=='exact_rational':out=self.rational(Fraction(row['numerator'],row['denominator']))
        elif op=='bound_variable':out=self.digest([op,self.name(row['name'])]),None
        elif op in ('sum','product'):out=self.combine(op,[F(i) for i in row['arguments']])
        elif op=='negative':out=self.combine('product',[self.rational(-1),F(row['argument'])])
        elif op=='positive_quotient':
            a,b=F(row['numerator']),F(row['denominator'])
            if b[1]==0:raise ValueError('Undefined endpoint quotient outside a lazy branch')
            if a[1]==0:out=self.zero
            elif b[1]==1:out=a
            elif a[1] is not None and b[1] is not None:out=self.rational(a[1]/b[1])
            else:out=self.digest([op,a[0],b[0]]),None
        elif op in ('analytic_unary','exact_original_source_flat_sigma'):
            value=F(row['argument']);name=row.get('name')
            if (op=='exact_original_source_flat_sigma' or name=='original_flat_sigma') and value[1] is not None and (value[1]<=0 or value[1]>=1):
                out=self.zero if value[1]<=0 else self.one
            elif name in ('original_flat_sigma_prime','original_flat_sigma_second') and value[1] is not None and (value[1]<=0 or value[1]>=1):out=self.zero
            elif name=='exp' and value[1]==0:out=self.one
            elif name in ('expm1','sin') and value[1]==0:out=self.zero
            elif name=='cos' and value[1]==0:out=self.one
            elif name=='log' and value[1]==1:out=self.zero
            else:out=self.digest([op,name,value[0]]),None
        elif op=='function_substitution':
            if not row.get('Z_independent_substitution'):raise ValueError('Fixed-Z endpoint substitution required')
            changed=dict(scope);changed[row['variable']]=F(row['value'])
            out=self.at(row['expression'],changed)
        elif op=='definite_integral':
            lo,hi=F(row['lower']),F(row['upper'])
            if lo[0]==hi[0]:out=self.zero;self.empty.add(node)
            else:
                changed=dict(scope);changed.pop(self.variables[row['variable']],None)
                body=self.at(row['integrand'],changed)
                out=self.digest([op,body[0],self.name(row['variable']),lo[0],hi[0]]),None
        elif op=='original_monotone_phase_inverse':
            out=self.digest([op,F(row['phase_function'])[0],self.name(row['angle_variable']),
                F(row['target_phase'])[0],F(row['angle_lower'])[0],F(row['angle_upper'])[0]]),None
        elif op=='substitute_original_inverse_angle':
            out=self.digest([op,F(row['body'])[0],F(row['inverse_angle'])[0],self.name(row['angle_variable'])]),None
        elif op=='exact_real_comparison':
            left,right=F(row['left']),F(row['right'])
            if left[1] is not None and right[1] is not None:
                tests={'<=':lambda a,b:a<=b,'>=':lambda a,b:a>=b,'<':lambda a,b:a<b,'>':lambda a,b:a>b,'==':lambda a,b:a==b}
                out=self.rational(int(tests[row['operator']](left[1],right[1])))
            else:out=self.digest([op,row['operator'],left[0],right[0]]),None
        elif op=='logical_or':
            values=[F(i) for i in row['arguments']]
            if any(v[1]==1 for v in values):out=self.one
            elif all(v[1]==0 for v in values):out=self.zero
            else:out=self.digest([op,sorted(v[0] for v in values if v[1]!=0)]),None
        elif op=='original_lazy_flat_branch':
            predicate=F(row['flat_predicate']);pred=predicate[0]
            if predicate[1] in (0,1):
                out=F(row['flat_value'] if predicate[1] else row['active_body'])
                self.cache[key]=out;return out
            if row.get('original_collar_coordinate')==self.b['coordinate'] and scope.get(self.b['coordinate'])==self.one:
                predicates=list(self.nodes[row['flat_predicate']]['arguments'])
                if row['original_collar_upper']!=self.b['loop']['first_micro']['loop']['original_collar_upper']:
                    raise ValueError('Unexpected first-only collar predicate')
                target=dict(operation='exact_real_comparison',operator='<=',left=self.b['coordinate'],right=row['original_collar_upper'])
                remove=[i for i in predicates if self.nodes[i]==target]
                if len(remove)!=1:raise ValueError('Actual first-only collar guard required')
                predicates.remove(remove[0]);self.collar.add(node)
                pred=F(predicates[0])[0] if len(predicates)==1 else self.digest(['logical_or',sorted(F(i)[0] for i in predicates)])
            out=self.digest([op,pred,F(row['active_body'])[0],F(row['flat_value'])[0]]),None
        else:
            admitted=('current_original_unique_analytic_core_component','mathematical_pi','mathematical_positive_infinity',
                'original_source_parameter','current_original_source_parameter','exact_positive_integer_expression',
                'current_original_leading_function_recipe','current_original_axis_pressure_function_derivative')
            if op not in admitted:raise ValueError('Unissued source fingerprint operation '+op)
            dependencies={'current_original_unique_analytic_core_component':('rho','Z_variable'),
                'current_original_leading_function_recipe':('coordinate','Z_variable'),
                'current_original_axis_pressure_function_derivative':('original_function','Z_variable')}.get(op,())
            # Function leaves still depend on their actual argument nodes.
            # A substitution may not turn an opaque f(s=0) into f(s=1).
            out=self.digest([op,row,[(key,F(row[key])[0]) for key in dependencies]]),None
        self.cache[key]=out;return out


class OriginalMicroViscousM15:
    def __init__(self, *, directory=DIRECTORY, precision=120):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30
        self.hashes={};self.code_hashes={}
        raw=(directory/REPAIRED).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=REPAIRED_SHA256:raise ValueError('Accepted repaired Rp report required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked repaired Rp predecessor required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes)
        self.hashes[REPAIRED]=REPAIRED_SHA256
        self.active=self.archive('micro_active_radial4_pressure_C3_functions')
        self.common=self.archive('lower_active_common_source_join_functions')
        self.first=self.saved('first_micro_full_m15_numeric')
        self.second=self.saved('second_micro_m15_numeric')
        self.proof=self.prove_interface()
        self.positive_budget_source=self.bind_positive_budget_source()
        name='src/openai_ns_reconstruction/original_micro_viscous_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def merge(self,group,base,target):
        for name,digest in group.items():
            if target.get(name,digest)!=digest or hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed trace source '+name)
            target[name]=digest
    def saved(self,suffix):
        name=PREFIX+'compliant_current_original_'+suffix+'.json';raw=(self.directory/name).read_bytes()
        if self.hashes.get(name)!=hashlib.sha256(raw).hexdigest():raise ValueError('Pinned numerical ancestor required')
        data=json.loads(raw)
        if data['source_family']!=self.family or not data['checks']['passed']:raise ValueError('Checked same-source numerical ancestor required')
        return data
    def archive(self,suffix):
        stem=PREFIX+'compliant_current_original_'+suffix;name=stem+'.json.gz'
        raw=(self.directory/name).read_bytes();receipt_raw=(self.directory/(stem+'_check.json')).read_bytes()
        data=json.loads(gzip.decompress(raw));receipt=json.loads(receipt_raw)
        if not receipt['all_passed'] or data['source_family']!=self.family or receipt['source_family']!=self.family:
            raise ValueError('Checked original trace source required')
        self.merge(receipt['input_hashes'],self.directory,self.hashes)
        self.hashes[name]=hashlib.sha256(raw).hexdigest()
        self.hashes[stem+'_check.json']=hashlib.sha256(receipt_raw).hexdigest()
        data['_receipt']=receipt;return data
    def bind_positive_budget_source(self):
        """Bind the actual positive upper/lower aggregation before using it.

        This proves the second aggregate bounds its upper trace. It need
        not dominate a different, looser first-chart bound on the trace
        that was already eliminated by the source identity.
        """
        requests={
            ('original_first_micro_m15.py','OriginalFirstMicroM15','finite_N_errors'):[
                "viscous=[None if v is None else v+lr for v in errors['Q_y']]"],
            ('original_second_micro_m15.py','OriginalSecondMicroM15','evaluate'):[
                "upper=self.finite_N_errors(Z,s,packet)",
                "lower=self.finite_N_errors(Z,iv.mpf(1),lower_packet)",
                "error[name]=[self.rows.sum_logs((upper[name][k],lower[name][k])) for k in (0,1)]",
                "bound=self.rows.sum_logs(error[name])",
                "top=iv.mpf(_bounds(self.audit,bound)[1])"],
            ('original_first_micro_pressure_rows.py','OriginalFirstMicroPressureRows','sum_logs'):[
                "values=[v for v in values if v is not None]",
                "high=max(_bounds(self.ctx,v)[1] for v in values)",
                "return self.iv.mpf(high)+self.iv.log(len(values))"],
            ('original_first_micro_m15.py','OriginalFirstMicroM15','evaluate'):[
                "lower=self.rows.evaluate(Z,'.5',selected_sc_multiple=True)",
                "if not all(v['exact_zero'] for name in ('U','Q','Q_y') for v in lower['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds']): raise ValueError('Actual flat lower finite-N pressure endpoint must be zero')"]}
        evidence=[]
        for (file,cls,fn),statements in requests.items():
            name='src/openai_ns_reconstruction/'+file;raw=(ROOT/name).read_bytes()
            if self.code_hashes.get(name)!=hashlib.sha256(raw).hexdigest():raise ValueError('Pinned budget producer required')
            tree=ast.parse(raw);body=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
            function=next(n for n in body.body if isinstance(n,ast.FunctionDef) and n.name==fn)
            dumps={ast.dump(n) for n in ast.walk(function)}
            for text in statements:
                node=ast.parse(text).body[0]
                if ast.dump(node) not in dumps:raise ValueError('Original positive budget formula changed '+text)
            evidence.append(dict(source=name,source_sha256=self.code_hashes[name],method=cls+'.'+fn,
                exact_source_statements=statements))
        return dict(passed=True,source=evidence,
            theorem='exp(max(log_upper,log_lower)+log2) >= upper_bound; subsequent positive C0/C1 union preserves domination',
            only_surviving_upper_trace_requires_a_bound=True,
            removed_first_chart_cap_need_not_be_dominated=True,cap_subtraction_used=False)
    def prove_interface(self):
        b=self.active['micro_active_radial4_bindings'];nodes=self.active['exact_graph_nodes']
        rec=self.active['_receipt'];checks=rec['checks']
        if not rec['current_original_n1_first_second_micro_full_mixed_source_join_installed']:
            raise ValueError('Actual mixed micro join certificate required')
        sc=self.common['lower_active_common_source_bindings']['sc_binding']
        common_nodes=self.common['exact_graph_nodes']
        if nodes[:len(common_nodes)]!=common_nodes:raise ValueError('Identical shared source graph prefix required')
        sign,mantissa,exponent,bits=sc['exact_mpf_tuple']
        if (not sc['exact_singleton_recomputed_not_interval_overlap'] or sign or mantissa<=0
                or mantissa.bit_length()!=bits or exponent+bits>-2):raise ValueError('Original positive exact sc<1/4 required')
        for p in self.first['samples']:
            if p['original_sc_tuple']!=sc['exact_mpf_tuple']:raise ValueError('Actual first phase selector changed')
        for key in ('coordinate','radius','physical_measure','actual_phase','roots'):
            if b['loop']['first_micro'][key]!=b['loop']['second_micro'][key]:raise ValueError('Actual micro source mismatch '+key)
        for name in ('first_micro','second_micro'):
            if b['modified'][name]['radial_background']!=b['background_y']:raise ValueError('Shared leading background required')
        if not (b['loop']['first_micro']['domain'][1]==b['loop']['second_micro']['domain'][0]==b['seam']):
            raise ValueError('Actual native seam must be s=1')
        if nodes[b['seam']]!=dict(operation='exact_rational',numerator=1,denominator=1):raise ValueError('Fixed seam1 required')
        visited=set()
        def radius_Z_independent(i):
            if i in visited:return
            visited.add(i);row=nodes[i];op=row['operation']
            if op=='bound_variable':
                if i!=b['coordinate']:raise ValueError('Unexpected physical radius variable')
            elif op in ('original_source_parameter','current_original_source_parameter','exact_positive_integer_expression'):
                if not row.get('Z_independent'):raise ValueError('Radius parameter must be Z-independent')
            elif op in ('sum','product'):
                for j in row['arguments']:radius_Z_independent(j)
            elif op=='negative':radius_Z_independent(row['argument'])
            elif op=='positive_quotient':
                radius_Z_independent(row['numerator']);radius_Z_independent(row['denominator'])
            elif op=='analytic_unary':radius_Z_independent(row['argument'])
            elif op not in ('exact_rational','mathematical_pi'):raise ValueError('Unissued physical radius operation '+op)
        radius_Z_independent(b['radius'])
        for key in ('actual_leading_micro_sources','actual_leading_recovery_pressure_and_FTC','periodic_lift','source_contract'):
            if not checks[key]['passed']:raise ValueError('Required source proof failed '+key)
        if not checks['actual_leading_micro_sources']['actual_first_complete_Duhamel_history_is_second_incoming']:
            raise ValueError('Complete actual mean history must cross seam')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        f=EndpointFingerprint(nodes,b);rows=[]
        for field in ('U','m','Q'):
            for radial in (0,1):
                for k in (0,1):
                    left=b['charts']['first_micro']['physical'][field][radial][k]
                    right=b['charts']['second_micro']['physical'][field][radial][k]
                    lh,rh=f.at(left)[0],f.at(right)[0]
                    if lh!=rh:raise ArithmeticError('Actual endpoint source mismatch '+str((field,radial,k,left,right)))
                    rows.append(dict(field=field,physical_y_order=radial,ordinary_Z_order=k,
                        left_node=left,right_node=right,canonical_endpoint_sha256=lh))
        for k in (0,1):
            left=b['modified']['first_micro']['finite_boundary']['R_Q_y'][k]
            right=b['modified']['second_micro']['finite_boundary']['R_Q_y'][k]
            if nodes[left].get('operation')!='product' or b['radius'] not in nodes[left]['arguments']:
                raise ValueError('Actual left physical radius trace factor required')
            if nodes[right].get('operation')!='product' or b['radius'] not in nodes[right]['arguments']:
                raise ValueError('Actual right physical radius trace factor required')
            lh,rh=f.at(left)[0],f.at(right)[0]
            if lh!=rh:raise ArithmeticError('Actual physical RQy trace mismatch')
            rows.append(dict(field='R_Q_y',physical_y_order=0,ordinary_Z_order=k,
                left_node=left,right_node=right,canonical_endpoint_sha256=lh))
        for name in ('m','h','k','e','p'):
            if not all(v==0 for v in b['loop']['first_micro']['incoming'][name]):raise ValueError('Original relative history inlet zero required')
        lower_samples=[p for p in self.first['samples'] if p['coordinate'].get('kind')=='selected_sc_multiple'
            and p['coordinate'].get('multiple')=='.5']
        if {p['Z'] for p in lower_samples}!={p['Z'] for p in self.parent['samples']}:
            raise ValueError('Saved exact sc/2 flat source packets required at every axial query')
        for p in lower_samples:
            if not p['exact_empty_increment'] or p['coordinate']['exact_sc_mpf_tuple']!=sc['exact_mpf_tuple']:
                raise ValueError('Actual fixed sc/2 lower endpoint required')
            if not all(row['exact_zero'] for rows in p['finite_N_error_logs_before_physical_measures'].values() for row in rows):
                raise ValueError('Accepted actual lower flat endpoint errors must be exactly zero')
        integer=nodes[b['N']]
        for p in self.parent['samples']:
            if integer!=p['actual_selected_integer']:raise ValueError('Same actual selected N required')
        return dict(passed=True,coordinate_node=b['coordinate'],seam_node=b['seam'],radius_node=b['radius'],
            physical_measure_node=b['jacobian'],phase_node=b['phase'],selected_integer_node=b['N'],
            exact_sc_tuple=sc['exact_mpf_tuple'],source_identity_rows=rows,
            empty_integrals_eliminated=sorted(f.empty),first_collar_predicates_proved_false=sorted(f.collar),
            common_ordinary_Z_independent_radius=True,leading_background_and_modified_trace_both_shared=True,
            original_lower_flat_source_zero=dict(saved_Z=sorted(p['Z'] for p in lower_samples),
                exact_sc_half_packets_verified=True,all_five_relative_incoming_histories_zero=True,
                source_producer_checks_actual_U_Q_Qy_errors=True),
            phase_lift_from_actual_periodic_composites=True,
            definition='delta(RQy)=R*(Qy_modified-Qy_original_leading), at the same physical s=1',
            signed_trace_centers_evaluated=False,upper_caps_not_used_as_function_values=True)

    def evaluate(self,Z):
        p=next((v for v in self.parent['samples'] if v['Z']==Z),None)
        if p is None:raise ValueError('Saved axial query required; no ancestor replay')
        second=next(v for v in self.second['samples'] if v['Z']==Z and v['coordinate']['value']=='2')
        terms=p['m15_C0_C1_terms'];a,b=terms[12],terms[16]
        if (a['name']!='first_micro_finite_N_viscous_error' or b['name']!='second_micro_finite_N_viscous_error'
                or not a['actual_function_binding']['lower_endpoint_exact_zero_checked']
                or not b['actual_function_binding']['both_upper_and_lower_pressure_endpoints_retained']
                or terms[:17]!=second['m15_C0_C1_terms']):raise ValueError('Actual signed first/second viscous partition required')
        if a['log_scale_terms'][0]['source_enclosure_exact_mpf_tuples']!=b['log_scale_terms'][0]['source_enclosure_exact_mpf_tuples']:
            raise ValueError('Same physical epsilon factor required')
        for row in second['finite_N_error_logs_before_physical_measures']['viscous']:
            if row['exact_zero']:raise ValueError('Remaining upper budget may not be declared zero')
        if len(b['log_scale_terms'])!=2 or not b['log_scale_terms'][1]['positive_error_bound_atom_only']:
            raise ValueError('Original physical epsilon times positive common budget required')
        def restore(row):
            raw=row.get('source_enclosure_exact_mpf_tuples')
            if raw is None:raw=[row[k+'_exact_mpf_tuple'] for k in ('lower','upper')]
            return self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in raw])
        cap=restore(b['log_scale_terms'][1]);cap_upper=_bounds(self.ctx,cap)[1]
        witness=[]
        for k,row in enumerate(second['finite_N_error_logs_before_physical_measures']['viscous']):
            combined=restore(row['log_absolute_upper'])
            if cap_upper<_bounds(self.ctx,combined)[1]:raise ValueError('Retained cap must dominate each positive combined endpoint budget')
            witness.append(dict(ordinary_Z_order=k,combined_endpoint_log_bound=copy.deepcopy(row),
                retained_common_log_cap=copy.deepcopy(b['log_scale_terms'][1]),
                common_cap_dominates_combined_budget=True))
        chain=contract_trace_chain([dict(name='first_micro',lower='rminus_exact_zero',upper='first_s1'),
                                    dict(name='second_micro',lower='second_s1',upper='second_s2')],{'second_s1':'first_s1'})
        if chain['signed_coefficients']!={'rminus_exact_zero':1,'second_s2':-1}:raise ValueError('Signed physical trace orientations changed')
        replacement=copy.deepcopy(b);replacement['name']='micro_shared_finite_N_viscous_outer_endpoint_error'
        replacement['actual_function_binding']=dict(error_enclosure_only=True,
            function='-delta(RQy)(second_micro s=2)',ordinary_Z1=True,
            internal_trace_proof=self.proof['source_identity_rows'],
            remaining_upper_bound='original positive upper+lower second-micro budget; no cap subtraction',
            bound_is_conservative_not_actual_signed_value=True,no_other_error_caps_removed=True)
        reduced=[]
        for term in p['reduced_m15_C0_C1_terms']:
            if term['name']==a['name']:continue
            reduced.append(copy.deepcopy(replacement if term['name']==b['name'] else term))
        # The previous normalized ledger is retained. Only the cancelled
        # shared-trace pair is represented by its surviving outer endpoint.
        comparisons=[];allowed={term['name'] for term in p['reduced_m15_C0_C1_terms']} - {a['name']}
        for k in (0,1):
            candidates=[]
            for sector in p['comparison']['normalized_forward_sectors']:
                if sector['name'] not in allowed:continue
                row=sector['coefficient_C0_C1'][k]
                values=row['source_enclosure_exact_mpf_tuples'];coefficient=self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in values])
                cap=_bounds(self.ctx,abs(coefficient))[1]
                if not cap:continue
                values=sector['normalized_log_scale_enclosure']['source_enclosure_exact_mpf_tuples']
                log=self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in values])+self.iv.log(self.iv.mpf(cap))
                candidates.append((_bounds(self.ctx,log)[1],replacement['name'] if sector['name']==b['name'] else sector['name']))
            value,name=max(candidates)
            comparisons.append(dict(ordinary_Z_order=k,name=name,
                normalized_log_cap_upper_approximation=self.ctx.nstr(value,40),
                normalized_log_cap_upper_exact_mpf_tuple=list(value._mpf_),
                cap_is_not_an_actual_residual_value=True))
        return dict(Z=Z,source_family=self.family,raw_m15_C0_C1_terms=copy.deepcopy(terms),
            correlated_reduced_m15_C0_C1_terms=reduced,original_raw_terms=213,previous_reduced_terms=201,
            correlated_reduced_terms=len(reduced),signed_trace_chain=chain,trace_identity=self.proof,
            old_budget_terms_retained_as_provenance=[copy.deepcopy(a),copy.deepcopy(b)],
            replacement_outer_endpoint=replacement,largest_remaining_individual_log_caps=comparisons,
            positive_upper_budget_witness=dict(source=self.positive_budget_source,serialized_C0_C1=witness,
                interpretation='the combined second bound covers D2; no domination of cancelled D1 cap is required'),
            independent_target=copy.deepcopy(p['independent_target']),
            contract=dict(raw_213_verbatim=True,proved_physical_trace_contracted=True,
                ordinary_Z1_and_radius_preserved=True,only_first_second_viscous_representation_changed=True,
                second_positive_aggregate_budget_retained_for_upper=True,all_volume_width_and_source_errors_retained=True,
                no_interval_cap_subtraction=True,ancestors_not_reexecuted=True),
            limits=dict(first_second_finite_N_viscous_internal_trace_resolved=True,
                remaining_outer_trace_signed_value_resolved=False,other_internal_traces_resolved=False,
                independent_numeric_Rp_seam=False,global_five_moment_numeric_join=False,
                exact_heat_exterior=False,stress_cone=False,flat_remainder=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
