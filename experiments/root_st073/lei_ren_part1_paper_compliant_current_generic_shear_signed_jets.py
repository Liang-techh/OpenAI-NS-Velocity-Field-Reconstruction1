"""Signed ordinary source-jet expressions with deferred positive quotients.

Leaves carry signed original function derivative covers, not chosen field
values. Saved covers and explicitly injected live queries stay distinct.
The derivative DAG retains real denominators, fixed log factors and R.
"""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_loop_jet_bounds as current
import lei_ren_part1_paper_compliant_current_generic_shear_inputs as expressions

packets,bounds=current.packets,current.bounds
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
NAME=PREFIX+'current_generic_shear_signed_jets.json'
RECEIPT=PREFIX+'current_generic_shear_signed_jets_check.json'
VIEWS=PREFIX+'current_generic_shear_signed_jets_views.json.gz'
GATE='current_original_signed_generic_input_jet_expressions_y2_Z1_available'
OPEN=current.OPEN
ORDERS=bounds.ORDERS
ZERO=(0,0)
O3_OWNER_CLASSES={
    'O3_slope_mu':('current_O3_transition_background_tensor','CurrentO3TransitionBackgroundTensor'),
    'O3_power':('current_pulse_entrance_incoming_background_tensor','CurrentPulseEntranceIncomingBackgroundTensor')}


def validate_source_receipt(checked,family):
    if (not checked.get('all_passed') or not checked.get(current.GATE)
        or checked.get('source_family')!=family or any(checked.get(k) for k in OPEN)):
        raise ValueError('Checked actual loop/source derivative bounds with original scope required')


def ordinary_axial_coefficient(row,k):
    """k! coefficient_k encloses D_Z^k; fixed base factors are unchanged."""
    if type(k) is not int or not 0<=k<=row.order:raise ValueError('Available original axial derivative order required')
    c=row.ctx
    return packets.FactoredJet(row.algebra,{mode:packets.IntervalTaylor.constant(c,jet[k]*math.factorial(k),0)
        for mode,jet in row.terms.items()},0)


def source_leaves(packet,delta):
    """Same original full rows and remaining inertial half-shift as checked norms."""
    c=packet.algebra.ctx;algebra=packet.algebra;field=packet.recover_original(delta)
    E=packet.velocity['theta'];V=packet.velocity['axial']
    rows=dict(E=E,C=tuple(E[j]-2*E[j+1] for j in range(3)),
        B=tuple(2*V[j+1] for j in range(3)))
    full=field['full_signed_stress_ordinary_y_rows']
    for name,component in (('nt','theta'),('nz','axial')):
        physical=[left+algebra.shift(right,(0,.5,0,0)) for left,right in
            zip(full['inertial_'+component+'_linear'],full['inertial_'+component+'_quadratic'])]
        # physical rows already differentiate sqrt(R). Divide by the same
        # shear F: the numerator is R*shape, so apply only the remaining .5.
        rows[name]=tuple(sum((physical[i]*(math.comb(j,i)*c.mpf('.5')**(j-i))
            for i in range(j+1)),algebra.lift(0)) for j in range(3))
    leaves={}
    for name,values in rows.items():
        for j,k in ORDERS:
            coefficient=ordinary_axial_coefficient(values[j],k)
            leaves[name+'_'+'y%d_Z%d'%(j,k)]=expressions.RadiusPolynomial(algebra,
                {1 if name in ('nt','nz') else 0:coefficient})
    return leaves


class SignedJetDAG:
    """Signed arithmetic instructions; no inversion of coefficient enclosures."""
    def __init__(self,leaves,certificates):
        self.leaves=dict(leaves);self.certificates=dict(certificates);self.nodes=[];self.keys={}
        self.algebra=next(iter(leaves.values())).algebra
        if any(v.algebra is not self.algebra for v in leaves.values()):raise ValueError('Same original factored source algebra required')
        self.zero=self.constant(0)
    def node(self,operation,**attributes):
        value=dict(operation=operation,**attributes);key=json.dumps(value,sort_keys=True)
        if key not in self.keys:self.keys[key]=len(self.nodes);self.nodes.append(value)
        return self.keys[key]
    def constant(self,value):
        if type(value) is not int:raise ValueError('Exact integer jet weights required')
        return self.node('constant',value=value)
    def leaf(self,name):
        if name not in self.leaves:raise ValueError('Original derivative source leaf required')
        if not self.leaves[name].terms:return self.zero
        return self.node('source_derivative',name=name)
    def add(self,*args):
        args=tuple(q for q in args if q!=self.zero)
        if not args:return self.zero
        return args[0] if len(args)==1 else self.node('sum',arguments=list(args))
    def negate(self,arg):return self.zero if arg==self.zero else self.node('negative',argument=arg)
    def multiply(self,*args):
        if self.zero in args:return self.zero
        args=tuple(q for q in args if q!=self.constant(1))
        if not args:return self.constant(1)
        return args[0] if len(args)==1 else self.node('product',arguments=list(args))
    def divide(self,numerator,denominator,certificate):
        if certificate not in self.certificates:raise ValueError('Actual positive-function denominator certificate required')
        if denominator==self.zero:raise ValueError('Identically zero defining source denominator')
        return self.zero if numerator==self.zero else self.node('positive_function_quotient',
            numerator=numerator,denominator=denominator,positive_certificate=certificate)
    def quotient_jet(self,numerator,denominator,certificate):
        """Differentiated n=d*g with ordinary binomial weights and minus signs."""
        result={};d0=denominator[ZERO]
        for j,k in sorted(ORDERS,key=lambda key:(sum(key),key)):
            correction=[]
            for i in range(j+1):
                for ell in range(k+1):
                    if i==ell==0:continue
                    correction.append(self.multiply(self.constant(math.comb(j,i)*math.comb(k,ell)),
                        denominator[(i,ell)],result[(j-i,k-ell)]))
            result[(j,k)]=self.divide(self.add(numerator[(j,k)],self.negate(self.add(*correction))),d0,certificate)
        return result
    def product_jet(self,left,right):
        return {(j,k):self.add(*[self.multiply(self.constant(math.comb(j,i)*math.comb(k,ell)),
            left[(i,ell)],right[(j-i,k-ell)]) for i in range(j+1) for ell in range(k+1)]) for j,k in ORDERS}
    def record(self,roots):
        return dict(nodes=self.nodes,roots={key:{'y%d_Z%d'%order:node for order,node in value.items()} for key,value in roots.items()},
            source_derivative_leaves={key:value.record() for key,value in self.leaves.items()},
            positive_denominator_certificates=self.certificates,
            signed_ordinary_derivative_orders=[list(k) for k in ORDERS],
            source_caps_not_substituted_for_function_values=True,
            actual_source_factors_radius_or_reciprocals_not_evaluated=True)


def from_packet(packet,delta,positive):
    """Construct exact signed jet recipes over original enclosing leaves."""
    leaves=source_leaves(packet,delta)
    c=packet.algebra.ctx
    certificates=dict(E=dict(function='E=Utheta/Pstar',log_positive_lower=positive['log_E_positive_lower']),
        C=dict(function='C=a*E',log_positive_lower=positive['log_C_positive_lower']),
        a=dict(function='a=C/E',log_positive_lower=positive['log_actual_a_positive_lower']))
    dag=SignedJetDAG(leaves,certificates)
    table=lambda name:{key:dag.leaf(name+'_'+'y%d_Z%d'%key) for key in ORDERS}
    E,C,B,nt,nz=(table(key) for key in ('E','C','B','nt','nz'))
    roots=dict(E=E,a=dag.quotient_jet(C,E,'E'),b=dag.quotient_jet(B,E,'E'),
        p1=dag.quotient_jet(nt,E,'E'),p2=dag.quotient_jet(nz,E,'E'),
        t0=dag.quotient_jet({k:dag.negate(v) for k,v in B.items()},C,'C'))
    bb=dag.product_jet(roots['b'],roots['b']);b2_over_a=dag.quotient_jet(bb,roots['a'],'a')
    kappa={k:dag.add(roots['a'][k],b2_over_a[k]) for k in ORDERS}
    roots['kappa_minus2']={k:dag.add(v,dag.negate(dag.constant(2))) if k==ZERO else v for k,v in kappa.items()}
    result=dict(chart=packet.chart,source_family=packet.source_family,
        source_provenance=packet.provenance,original_radius_source=packet.provenance['original_radius_source'],
        fixed_original_log_basis_names=packets.LOG_NAMES,fixed_original_log_bases=packet.algebra.logs,
        jet_expression_dag=dag.record(roots),
        defining_source_kind='original analytic source derivatives with signed enclosing coefficients',
        original_five_histories_and_separate_P0_retained=True,
        pressure_energy_meridional_terms_in_p1_p2_retained=True,
        inertial_sqrt_R_to_R_remaining_half_shift_applied_once=True,
        R_not_materialized_or_applied_twice=True,
        quotient_expression_not_materialized_point_value=True,
        actual_branch_Delta_stays_signed_source_expression=True,
        original_O3_small_excess_definition=positive.get('exact_a_minus2_source'),
        original_O3_positive_mu=positive.get('actual_positive_mu'),
        normalized_units=dict(E='Utheta/Pstar',p1='Itheta/F',p2='Iz/F',F='Utheta/sqrt(2R)',
            derivatives='ordinary y=logR, ordinary Z; coefficient factorials included'))
    if packet.algebra.proofs or packet.algebra.final_rows:raise ValueError('Production signed source factors were resolved')
    return result


def exact_theorem():
    y,Z=s.symbols('y Z');n=s.Function('n')(y,Z);d=s.Function('d')(y,Z);g=n/d
    checks={}
    for j,k in ORDERS:
        correction=sum(math.comb(j,i)*math.comb(k,ell)*s.diff(d,y,i,Z,ell)*s.diff(g,y,j-i,Z,k-ell)
            for i in range(j+1) for ell in range(k+1) if (i,ell)!=(0,0))
        rhs=(s.diff(n,y,j,Z,k)-correction)/d
        if s.cancel(s.diff(g,y,j,Z,k)-rhs)!=0:raise ArithmeticError('Signed quotient derivative identity failed')
        checks['y%d_Z%d'%(j,k)]=True
    R=s.exp(y);f=s.Function('inertial_shape')(y,Z)
    physical=[s.diff(s.exp(y/2)*f,y,j)/s.exp(y/2) for j in range(3)]
    for j in range(3):
        converted=sum(math.comb(j,i)*s.Rational(1,2)**(j-i)*physical[i] for i in range(j+1))
        if s.simplify(converted-s.diff(R*f,y,j)/R)!=0:raise ArithmeticError('Signed numerator radius conversion failed')
        checks['remaining_radius_half_shift_y%d'%j]=True
    return dict(passed=True,exact_signed_ordinary_quotient_and_radius_identities=checks,
        original_axial_Taylor_derivative_rule='D_Z^k f=k!*coefficient_k at the same basepoint',
        fixed_source_log_bases_not_differentiated_again=True,
        no_absolute_value_or_positive_part_in_signed_expression_recurrence=True,
        source_function_positivity_not_coefficient_box_positivity=True)


class CurrentSignedInputJets:
    def __init__(self,*,service=None,O3_owners=None):
        self.provider=current.source.domain.source.CurrentO3Sources(service)
        self.service=self.provider.service;self.ctx=self.service.ctx;self.family=self.service.family
        self.O3_owners=dict(O3_owners or {})
        if set(self.O3_owners)-set(O3_OWNER_CLASSES):raise ValueError('Known original O3 owner chart required')
        checked=json.loads((HERE/current.RECEIPT).read_bytes())
        validate_source_receipt(checked,self.family)
        self.service.bind_hashes(checked['input_hashes']);self.service.bind_hashes({current.RECEIPT:sha(current.RECEIPT)})
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name),Path(expressions.__file__).name:sha(Path(expressions.__file__).name)})
        self.inventory=json.loads((HERE/current.NAME).read_bytes())['current_actual_loop_jet_log_bounds_by_chart']
        self.theorem=exact_theorem()
    def positive(self,chart,c=None):
        if chart not in self.inventory:raise ValueError('Admitted positive original quotient chart required; axis/core separate')
        target=self.ctx if c is None else c
        def decode(v):
            if isinstance(v,dict) and 'lower' in v and 'upper' in v:return packets.interval(target,v)
            if isinstance(v,dict):return {k:decode(q) for k,q in v.items()}
            if isinstance(v,list):return [decode(q) for q in v]
            return v
        return decode(self.inventory[chart]['actual_positive_denominator_theorem'])
    def saved(self,chart):
        positive=self.positive(chart)
        packet=self.service.saved(chart) if chart in packets.CHARTS else self.provider.saved(chart)
        return from_packet(packet,self.service.data['delta'],positive)
    def query(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.positive(chart)
        if chart in packets.CHARTS:
            packet=self.service.query(chart,Z,coordinate,log_tau,theta,viscosity)
        else:
            if chart not in self.O3_owners:raise ValueError('Arbitrary O3 coordinates require an injected checked owner; saved fallback forbidden')
            stem,cls=O3_OWNER_CLASSES[chart];owner=self.O3_owners[chart]
            if type(owner).__name__!=cls or type(owner).__module__!=PREFIX+stem or not owner.acceptance_loaded:
                raise ValueError('Exact existing checked original O3 owner required')
            if dict(zip(packets.FAMILY_KEYS,(owner.family,owner.source,owner.datum_sha)))!=self.family:
                raise ValueError('Injected O3 owner source/axis datum differs')
            owner.assert_graph();self.service.bind_hashes(owner.hashes)
            # Validate/bind the checked original regional source, without creating owners.
            self.provider.saved_view(chart)
            view=owner.chart(chart,Z,coordinate,log_tau,theta,viscosity)
            proxy=SimpleNamespace(ctx=owner.ctx,family=self.family,service=self.service)
            packet=type(self.provider).adapt(proxy,chart,view,dict(mode='injected_current_O3_source_cover',
                cache_cover=False,arbitrary_coordinates_evaluated=True,owner_class=cls,
                receipt=PREFIX+stem+'_check.json',source_graph_asserted=True,
                adaptation_context='original checked owner; global constants and positive proof copied exactly'))
        return from_packet(packet,packet.algebra.ctx.mpf(self.service.data['delta']),self.positive(chart,packet.algebra.ctx))
    def run(self):
        records={chart:self.saved(chart) for chart in self.inventory}
        (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(packets.encode(records),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
        self.service.bind_hashes({VIEWS:sha(VIEWS)})
        result=dict(source_family=self.family,source_cover_chart_count=len(records),
            signed_source_jet_root_count=sum(len(v['jet_expression_dag']['roots'])*len(ORDERS) for v in records.values()),
            signed_source_derivative_expression_views=VIEWS,
            exact_signed_source_jet_theorem=self.theorem,
            signed_ordinary_input_jet_inventory={chart:dict(input_roots=list(v['jet_expression_dag']['roots']),
                node_count=len(v['jet_expression_dag']['nodes']),source_leaf_count=len(v['jet_expression_dag']['source_derivative_leaves']),
                mode=v['source_provenance']['mode']) for chart,v in records.items()},
            **{GATE:True},**dict.fromkeys(OPEN,False),
            signed_current_point_loop_or_inverse_jets_installed=False,
            signed_source_query_requires_injected_exact_checked_owner=True,
            successful_arbitrary_coordinate_live_owner_query_tested=False,
            source_graph_ancestor_constructors_called=False,
            scope='Exact signed derivative expression DAG for a,b,p1,p2,E,t0,Delta through y2/Z1 over17 original coefficient covers, with actual positive denominator proofs and explicit existing-owner query routing. Saved covers are not point functions; no physical loop/inverse, changed histories/repair/N or true recursion admission.',
            input_hashes=self.service.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        print('Signed original input derivative expressions PASS:17 charts,714 roots; source factors deferred',flush=True)
        return result


def run():return CurrentSignedInputJets().run()


if __name__=='__main__':run()
