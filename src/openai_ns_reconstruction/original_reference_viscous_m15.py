"""Source-bound physical germs and finite-N reference inlet matching.

Exact polynomial checks consume only source functions and rational identities.
Source error caps are never function values.
"""
from collections import Counter
from fractions import Fraction
import copy
import hashlib
import json
import math
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import ActiveEndpointFingerprint,OriginalActiveViscousM15


class PhysicalEndpointFingerprint(ActiveEndpointFingerprint):
    """Bounded exact ring checks with source-certified denominator division."""
    term_limit=16384
    division_step_limit=50000
    multiplication_limit=400000

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.base=ActiveEndpointFingerprint(self.nodes,self.b)
        self.polynomials={};self.denominators={}
        self.certified_quotients=set();self.exact_divisions=set()
        self.maximum_terms=0

    def polynomial(self,value):
        digest,rational=value
        return {():Fraction(rational)} if rational is not None else self.polynomials.get(digest,{(digest,):Fraction(1)})

    def guard(self,terms):
        if len(terms)>self.term_limit:raise ArithmeticError('Exact endpoint polynomial term limit exceeded')
        self.maximum_terms=max(self.maximum_terms,len(terms))

    @staticmethod
    def monomial(values):
        counts=Counter(values)
        for token in tuple(counts):
            if token.startswith('inverse:'):
                base=token[8:];n=min(counts[token],counts.get(base,0))
                counts[token]-=n;counts[base]-=n
        return tuple(sorted(t for t,n in counts.items() for _ in range(n)))

    @staticmethod
    def monomial_quotient(numerator,denominator):
        counts=Counter(numerator);divisor=Counter(denominator)
        if any(counts[t]<n for t,n in divisor.items()):return None
        counts.subtract(divisor)
        return tuple(sorted(t for t,n in counts.items() for _ in range(n)))

    @staticmethod
    def leading(terms):
        return max(terms,key=lambda m:(len(m),m))

    def exact_divide(self,numerator,denominator):
        """Multivariate rational polynomial long division, no numerical tests."""
        work={m:c for m,c in numerator.items() if c};quotient={}
        divisor={m:c for m,c in denominator.items() if c}
        if not divisor:raise ValueError('Zero denominator in exact endpoint division')
        lead=self.leading(divisor);coefficient=divisor[lead]
        steps=0
        while work:
            mon=self.leading(work);factor=self.monomial_quotient(mon,lead)
            if factor is None:return None
            value=work[mon]/coefficient;quotient[factor]=quotient.get(factor,0)+value
            for m,c in divisor.items():
                key=tuple(sorted(factor+m));work[key]=work.get(key,0)-value*c
                if not work[key]:del work[key]
            steps+=1
            if steps>self.division_step_limit:raise ArithmeticError('Exact endpoint division step limit exceeded')
            self.guard(work)
        return quotient

    def simplify(self,terms):
        terms={m:c for m,c in terms.items() if c}
        if not terms:return terms
        # Pull a denominator atom common to every term, then require an
        # exact factor theorem for the remaining numerator polynomial.
        common=set(next(iter(terms)))
        for mon in terms:common.intersection_update(mon)
        for token in sorted(t for t in common if t.startswith('inverse:')):
            digest=token[8:]
            if digest not in self.denominators:continue
            numerator={self.monomial_quotient(m,(token,)):c for m,c in terms.items()}
            divided=self.exact_divide(numerator,self.denominators[digest])
            if divided is not None:
                self.exact_divisions.add(digest)
                return self.simplify(divided)
        return terms

    def finish(self,terms):
        self.guard(terms);terms=self.simplify(terms)
        if not terms or set(terms)=={()}:return self.rational(terms.get((),0))
        if len(terms)==1:
            mon,c=next(iter(terms.items()))
            if len(mon)==1 and c==1:return mon[0],None
        digest=self.digest(['exact_source_polynomial',sorted((list(m),c.numerator,c.denominator) for m,c in terms.items())])
        self.polynomials[digest]=terms
        return digest,None

    def combine(self,op,values):
        if op=='sum':
            terms={}
            for value in values:
                for m,c in self.polynomial(value).items():terms[m]=terms.get(m,0)+c
            return self.finish(terms)
        terms={():Fraction(1)}
        for value in values:
            other=self.polynomial(value);product={}
            if len(terms)*len(other)>self.multiplication_limit:
                raise ArithmeticError('Exact endpoint multiplication work limit exceeded')
            for a,c in terms.items():
                for b,d in other.items():
                    mon=self.monomial(a+b);product[mon]=product.get(mon,0)+c*d
                    if len(product)>self.term_limit:raise ArithmeticError('Exact endpoint product term limit exceeded')
            self.guard(product);terms={m:c for m,c in product.items() if c}
            if not terms:break
        return self.finish(terms)

    def at(self,node,scope=None):
        row=self.nodes[node]
        if row['operation'] in ('analytic_unary','exact_original_source_flat_sigma','original_monotone_phase_inverse'):
            return self.base.at(node,scope)
        if row['operation']!='positive_quotient' or not row.get('source_positive_certificate'):
            return super().at(node,scope)
        scope=scope if scope is not None else {self.b['coordinate']:self.one}
        if node in scope:return scope[node]
        key=node,tuple(sorted((k,v[0]) for k,v in scope.items()))
        if key in self.cache:return self.cache[key]
        numerator=self.at(row['numerator'],scope);denominator=self.at(row['denominator'],scope)
        if denominator[1]==0:raise ValueError('Undefined certified endpoint quotient')
        self.certified_quotients.add(node)
        if denominator[1] is not None:
            out=self.combine('product',[numerator,self.rational(1/denominator[1])])
        elif numerator[0]==denominator[0]:out=self.one
        else:
            poly=self.polynomial(denominator)
            if len(poly)==1:
                mon,c=next(iter(poly.items()));inverse=[]
                for token in mon:inverse.append((token[8:] if token.startswith('inverse:') else 'inverse:'+token,None))
                out=self.combine('product',[numerator,self.rational(1/c)]+inverse)
            else:
                self.denominators[denominator[0]]=poly
                out=self.combine('product',[numerator,('inverse:'+denominator[0],None)])
        self.cache[key]=out
        return out


class PhysicalGermRecipes:
    """Ordinary Z jets of source-bound physical-y derivatives.

    New node IDs belong only to this recipe graph. Source ports explicitly
    name original node IDs and endpoint substitution; no numeric atom or
    claimed original graph node is introduced.
    """
    def __init__(self,fields,radius,Z,Pstar,delta,coordinate,endpoint):
        self.nodes=[];self.index={};self.coordinate=coordinate;self.endpoint=endpoint
        self.zero=self.constant(0);self.one=self.constant(1)
        self.fields={key:[self.source(i,k) for k,i in enumerate(rows[:3])] for key,rows in fields.items()}
        self.R=[self.source(radius,0),self.zero,self.zero]
        self.Z=[self.source(Z,0),self.one,self.zero]
        self.P=[self.source(Pstar,0),self.zero,self.zero]
        self.delta=self.fixed(delta,3)

    def node(self,operation,**attrs):
        row=dict(operation=operation,**attrs);key=json.dumps(row,sort_keys=True,separators=(',',':'))
        if key not in self.index:self.index[key]=len(self.nodes);self.nodes.append(row)
        return self.index[key]
    def constant(self,value):
        v=Fraction(value);return self.node('exact_rational',numerator=v.numerator,denominator=v.denominator)
    def source(self,node,order):
        return self.node('original_endpoint_function',original_node=node,ordinary_Z_order=order,
            original_coordinate=self.coordinate,original_endpoint=self.endpoint,function_not_numeric_atom=True)
    def fixed(self,value,size):return [self.constant(value)]+[self.zero]*(size-1)
    def add(self,*jets):
        size=len(jets[0])
        if any(len(v)!=size for v in jets):raise ValueError('Common ordinary-Z jet size required')
        return [self.node('sum',arguments=[v[k] for v in jets]) for k in range(size)]
    def scale(self,jet,value):
        c=self.constant(value);return [self.node('product',arguments=[c,v]) for v in jet]
    def mul(self,a,b):
        if len(a)!=len(b):raise ValueError('Common ordinary-Z jet size required')
        result=[]
        for k in range(len(a)):
            terms=[self.node('product',arguments=[self.constant(math.comb(k,j)),a[j],b[k-j]]) for j in range(k+1)]
            result.append(self.node('sum',arguments=terms))
        return result
    def div(self,a,b,certificate):
        if len(a)!=len(b):raise ValueError('Common ordinary-Z quotient size required')
        q=[]
        for k in range(len(a)):
            terms=[a[k]]
            for j in range(1,k+1):terms.append(self.node('product',arguments=[self.constant(-math.comb(k,j)),b[j],q[k-j]]))
            numerator=self.node('sum',arguments=terms)
            q.append(self.node('positive_quotient',numerator=numerator,denominator=b[0],source_positive_certificate=certificate))
        return q
    def pair_add(self,*pairs):return tuple(self.add(*(p[k] for p in pairs)) for k in (0,1))
    def pair_scale(self,p,c):return tuple(self.scale(v,c) for v in p)
    def pair_mul(self,p,q):return self.mul(p[0],q[0]),self.add(self.mul(p[1],q[0]),self.mul(p[0],q[1]))
    def pair_div(self,p,q,certificate):
        value=self.div(p[0],q[0],certificate)
        return value,self.div(self.add(p[1],self.scale(self.mul(value,q[1]),-1)),q[0],certificate)

    def build(self):
        f=self.fields;E,V=(f[k] for k in ('E','V'));E2=self.mul(E,E);EV=self.mul(E,V)
        zero=self.fixed(0,3)
        y=dict(E=self.scale(E,Fraction(1,10)),V=zero,m=self.add(V,self.scale(f['m'],-1)),
            h=self.add(E,self.scale(f['h'],Fraction(-3,2))),
            k=self.add(EV,self.scale(f['k'],Fraction(-3,2))),
            e=self.add(self.mul(V,V),self.scale(E2,Fraction(-1,2)),self.scale(f['e'],-1)),
            p=self.scale(E2,Fraction(1,2)),P0=zero,pressure=self.scale(E2,Fraction(1,2)))
        yy=dict(E=self.scale(E,Fraction(1,100)),V=zero,m=self.scale(y['m'],-1),
            h=self.add(y['E'],self.scale(y['h'],Fraction(-3,2))),
            k=self.add(self.scale(EV,Fraction(1,10)),self.scale(y['k'],Fraction(-3,2))),
            e=self.add(self.scale(E2,Fraction(-1,10)),self.scale(y['e'],-1)),
            p=self.scale(E2,Fraction(1,10)),P0=zero,pressure=self.scale(E2,Fraction(1,10)))
        pair=lambda key:(f[key][:2],y[key][:2])
        shifted=lambda key:(f[key][1:3],y[key][1:3])
        const=lambda c:(self.fixed(c,2),self.fixed(0,2))
        add,mul,scale,divide=self.pair_add,self.pair_mul,self.pair_scale,self.pair_div
        z=(self.Z[:2],self.fixed(0,2));de=(self.delta[:2],self.fixed(0,2))
        one=const(1);one_de=add(one,scale(de,-1));one_plus_de=add(one,de)
        d=add(one,scale(mul(z,z),-1));L=add(one,scale(mul(de,mul(z,z)),-1))
        transport=add(mul(mul(one_de,z),pair('m')),mul(d,shifted('m')))
        linear=add(scale(pair('V'),-1),scale(mul(one_de,add(pair('m'),scale(mul(z,shifted('m')),-1))),Fraction(1,2)))
        quadratic=add(mul(pair('V'),transport),scale(mul(mul(de,z),pair('e')),2),
            scale(mul(d,shifted('e')),-1),scale(mul(mul(one_plus_de,z),pair('pressure')),2),
            scale(mul(d,shifted('pressure')),-1))
        P=(self.P[:2],self.fixed(0,2));Iz=divide(add(linear,mul(P,quadratic)),L,'0<delta=1/1000<1 and |Z|<=1')
        R=(self.R[:2],self.R[:2]);ratio=divide(R,pair('E'),'original E=Utheta/Pstar>0')
        p2=mul(ratio,Iz)
        return dict(recipe_nodes=self.nodes,ordinary_Z_orders=[0,1,2],physical_y_orders=[0,1,2],
            common_fields=dict(y0=f,y1=y,y2=yy),root_p2=dict(value_Z0_Z1=p2[0],physical_y_Z0_Z1=p2[1],
                Iz_Z0_Z1=Iz[0],Iz_y_Z0_Z1=Iz[1],identity='p2_y=(9/10)*p2+(R/E)*Iz_y'),
            source_functions_not_evaluated=True,no_error_caps_used=True,
            genuine_temporal_coefficient_recursion=False,physical_jet_only=True)


ACTIVE=PREFIX+'compliant_current_original_active_viscous_m15_numeric.json'
ACTIVE_SHA256='268f00eea0e118d2afe5b6e9c82aecb3b2c43ceae30887cea4da53965c34bb72'


class OriginalReferenceGerm(OriginalActiveViscousM15):
    """Hydrate saved sources and construct common physical jets only."""
    def __init__(self,*,directory=DIRECTORY,precision=100):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/ACTIVE).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=ACTIVE_SHA256:raise ValueError('Accepted active state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[ACTIVE]=ACTIVE_SHA256
        self.modified=self.archive('modified_Q4_pressure_C3_functions')
        self.loop=self.archive('prescribed_loop_C5_functions')
        self.outer=self.archive('prescribed_R100_Rm_C6_functions')
        self.nodes=self.modified['exact_graph_nodes'];self.l=self.loop['prescribed_loop_C5_bindings']
        self.bg=self.outer['prescribed_R100_Rm_C6_bindings']['charts']
        for source in (self.loop,self.outer):
            if self.nodes[:len(source['exact_graph_nodes'])]!=source['exact_graph_nodes']:raise ValueError('Unchanged original source prefix required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        self.f=ActiveEndpointFingerprint(self.nodes,dict(coordinate=self.l['charts']['reference']['coordinate'],loop=self.l['charts']))
        self.scope={self.l['charts']['reference']['coordinate']:self.f.one,self.l['charts']['restoration']['coordinate']:self.f.zero}
        self.common_rows=self.common_sources()
        self.radial_laws=self.leading_radial_laws()
        z=[i for i,r in enumerate(self.nodes) if r['operation']=='bound_variable' and r['name']=='Z']
        if z!=[20] or self.nodes[11]!=dict(operation='exact_rational',numerator=1,denominator=1000):raise ValueError('Original Z/delta source required')
        outer=self.outer['prescribed_R100_Rm_C6_bindings']
        self.recipes=PhysicalGermRecipes(self.bg['reference']['fields'],self.l['charts']['reference']['radius'],
            z[0],outer['Pstar'],Fraction(1,1000),self.l['charts']['reference']['coordinate'],
            self.l['charts']['reference']['domain'][1]).build()
        from .original_reference_germ_guards import build_guards
        self.guards=build_guards(self.nodes,self.l,self.outer,self.loop,self.modified,ROOT)
        if not self.guards['passed']:raise ValueError('Actual physical germ source guards required: '+str([k for k,v in self.guards['checks'].items() if not v['passed'] and k!='global_periodic_lift']))
        for name in ('original_reference_viscous_m15.py','original_reference_germ_guards.py'):
            name='src/openai_ns_reconstruction/'+name
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def common_sources(self):
        rows=[]
        def compare(group,quantity,k,left,right):
            a,b=self.f.at(left,self.scope)[0],self.f.at(right,self.scope)[0]
            if a!=b:raise ArithmeticError('Unproved common physical input '+str((group,quantity,k,left,right)))
            rows.append(dict(group=group,quantity=quantity,ordinary_Z_order=k,left_node=left,right_node=right,canonical_endpoint_sha256=a))
        for key in ('radius','global_offset','actual_phase'):
            compare('geometry',key,0,self.l['charts']['reference'][key],self.l['charts']['restoration'][key])
        for group in ('fields','radial','actual'):
            for key in self.bg['reference'][group]:
                for k in range(3):compare(group,key,k,self.bg['reference'][group][key][k],self.bg['restoration'][group][key][k])
        for key in self.l['charts']['reference']['roots']:
            for k in range(2):compare('root_values',key,k,self.l['charts']['reference']['roots'][key][k],self.l['charts']['restoration']['roots'][key][k])
        for key in ('B_over_Pstar','B_over_Pstar_phi'):
            for k in range(2):compare('loop_values',key,k,self.l['charts']['reference']['loop'][key][k],self.l['charts']['restoration']['loop'][key][k])
        for key in self.bg['reference']['moments']:
            for k in range(3):compare('retained_history_inlet',key,k,self.bg['reference']['moments'][key][k],self.bg['restoration']['incoming_histories'][key][k])
        return rows

    def leading_radial_laws(self):
        ring=PhysicalEndpointFingerprint(self.nodes,dict(coordinate=self.l['charts']['reference']['coordinate']))
        rows=[]
        for chart in ('reference','restoration'):
            bg=self.bg[chart]
            for k in range(3):
                e=ring.at(bg['fields']['E'][k],{});ey=ring.at(bg['radial']['E_y'][k],{})
                difference=ring.combine('sum',[ey,ring.combine('product',[ring.rational(Fraction(-1,10)),e])])
                if difference!=ring.zero:raise ArithmeticError('Actual whole-chart E_y=E/10 source identity required')
                rows.append(dict(chart=chart,ordinary_Z_order=k,E_node=bg['fields']['E'][k],E_y_node=bg['radial']['E_y'][k],
                    exact_zero_sha256=difference[0],endpoint_sampling_used=False,whole_native_chart_function_identity=True))
        return dict(passed=True,identity='E_y=E/10 on both native charts; hence E_yy=E/100',rows=rows,
            exact_polynomial_maximum_terms=ring.maximum_terms,certified_quotient_nodes=sorted(ring.certified_quotients),
            source_values_selected=False,physical_y2_E_law_by_differentiation=True)

    def ledger(self):
        digest=lambda v:hashlib.sha256(json.dumps(v,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode('utf-8')).hexdigest()
        return [dict(Z=s['Z'],raw_terms=213,reduced_terms=199,
            raw_ledger_sha256=digest(s['raw_m15_C0_C1_terms']),reduced_ledger_sha256=digest(s['correlated_reduced_m15_C0_C1_terms']),
            independent_target_sha256=digest(s['independent_target'])) for s in self.parent['samples']]

    def evaluate(self,*args,**kwargs):
        raise ValueError('This source-germ owner defines functions; it does not evaluate original numerical axial values')

    def build(self):
        return dict(source_family=self.family,parent_report=ACTIVE,parent_report_sha256=ACTIVE_SHA256,
            common_source_identity_rows=self.common_rows,leading_radial_laws=self.radial_laws,
            source_guards=self.guards,physical_germ_recipes=self.recipes,unchanged_m15_ledger=self.ledger(),
            physical_two_jet_certificate=dict(passed=True,
                method='actual_common_inputs_and_physical_Duhamel_ODE_jets',
                radius_variable='y=log(R/110), d_y=R*d_R; not a Cartesian component',
                E_y='E/10',E_yy='E/100',V_y='0 at the seam',V_yy='0 at the seam',
                V_yy_proof='original C-infinity flat sigma at t=0; derivative of the zero lazy flat branch is zero',
                normalized_moment_first_and_second_y_jets_from_original_ODE=True,
                all_six_nonzero_history_inlets_retained=True,
                ordinary_Z_jets_are_function_identities_not_axial_samples=True,
                universal_loop_native_derivative_naturality_not_yet_certified=True,
                direct_endpoint_fingerprint_for_modified_Uy=False,
                source_values_selected=False),
            contract=dict(accepted_ancestors_not_reexecuted=True,raw_213_and_reduced_199_unchanged=True,
                all_old_error_budgets_and_independent_target_unchanged=True,
                derived_recipe_nodes_are_not_claimed_as_original_graph_nodes=True,
                source_function_definitions_not_numeric_field_evaluation=True,
                no_new_trace_alias=True,no_error_cap_subtraction=True),
            limits=dict(physical_input_two_jet_functions_implemented=True,
                source_bound_global_periodic_lift=False,reference_restoration_modified_Uy_Qy_RQy_identity=False,
                reference_viscous_budget_contraction=False,independent_numeric_Rp_seam=False,
                global_five_moment_numeric_join=False,exact_heat_exterior=False,stress_cone=False,
                flat_remainder=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
