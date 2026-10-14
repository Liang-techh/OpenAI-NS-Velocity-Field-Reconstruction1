"""Directed physical C0/C1 ranges for the actual reserved n1 coefficients.

The finite matrix bounds enclose the selected operator; they do not replace
its functions. Absolute source scales stay logarithmic throughout.
"""
import copy
import ast
import inspect
import math
import textwrap
from types import SimpleNamespace
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge as full
import lei_ren_part1_paper_compliant_current_transition_factored_native_source as transition
import lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N as original_O3
import lei_ren_part1_paper_compliant_current_original_C3_source_targets as arc_targets

ep,interval,record,source_precision=full.ep,full.interval,full.record,full.source_precision
COEFFICIENTS=('alpha1','alpha2','beta1','beta2','beta3')
HERE,PREFIX,sha=full.HERE,full.PREFIX,full.sha
NAME=PREFIX+'current_original_n1_reserved_coefficient_enclosures.json.gz'
RECEIPT=PREFIX+'current_original_n1_reserved_coefficient_enclosures_check.json'
GATE='current_original_n1_actual_reserved_coefficient_C1_numeric_enclosures_installed'
OPEN=full.OPEN


def exact_original_arc_source_identity(operator):
    """Original Rd -> Rw -> Rc body; never infer it from range overlap."""
    proof=transition.canonical_amplitude_proof()
    if not proof['passed'] or not proof['exact_original_O2_exit_J_one_half']:
        raise ValueError('Actual canonical Rd amplitude and J(1) identity required')
    fn=ast.parse(textwrap.dedent(inspect.getsource(original_O3.background_cell))).body[0]
    branch=next(n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.unparse(n.test)=="chart == 'transition'")
    value=lambda rows:next(n.value for n in rows if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='factor' for t in n.targets))
    mu,LP,z=sy.symbols('original_mu logP Z',real=True)
    ctx=SimpleNamespace(mpf=lambda value:sy.Rational(str(value)),exp=sy.exp)
    env=dict(c=ctx,mu_cover=mu,K={'J':sy.Rational(1,2)})
    first=eval(compile(ast.Expression(value(branch.body)),'<actual-O2-factor>','eval'),{},dict(env,t=sy.Integer(1)))
    second=eval(compile(ast.Expression(value(branch.orelse)),'<actual-O3-factor>','eval'),{},dict(env,t=sy.Integer(2)))
    body=sy.exp(-LP/2+sy.Rational(3,10))*first*second/(1+z*z)
    scale=-LP/2+sy.Rational(1,2)
    finite=sy.exp(-sy.Rational(17,10)-sy.Rational(5,2)*mu)/(1+z*z)
    if sy.simplify(body-sy.exp(scale)*finite)!=0:
        raise ValueError('Actual Rc angular amplitude source factorization differs')
    E=next(n.value for n in fn.body if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='E' for t in n.targets))
    if ast.unparse(E)!='f.scale(u1, factor)':raise ValueError('Actual amplitude Z jets use another multiplier')
    g=operator.graph;canonical=full.true.parent.linear.read(arc_targets.NAME)
    canonical_nodes=canonical['actual_C3_density_transport_targets']['actual_terminal_amplitude']
    identity_keys=('operation','source_family','native_chart','recipe','quantity','Z_order','coordinate','Z_variable')
    for order,node in enumerate(operator.background['actual_A_rc_ordinary_Z0_Z3']):
        leaf=g.nodes[node]
        if (leaf['operation']!='current_original_leading_function_recipe' or leaf['native_chart']!='O3_power'
            or leaf['quantity']!='E' or leaf['Z_order']!=order
            or leaf['recipe']['function']!='WholeZAllNOuterRcFunctions.leading_packet'
            or leaf['recipe']['quantity_paths']['E']!='original_generic_source.common_velocity_E_axial5[Z_order]'
            or leaf.get('Taylor_coefficient_factorial',1)!=(1,1,2,6)[order]
            or g.nodes[leaf['coordinate']]!=dict(operation='exact_rational',numerator=2,denominator=1)):
            raise ValueError('Actual A_rc C1 source leaf differs')
        saved=canonical['exact_graph_nodes'][canonical_nodes[order]]
        if any(leaf.get(k)!=saved.get(k) for k in identity_keys):
            raise ValueError('Actual angular source family/provider/coordinate/variable differs')
    return dict(passed=True,actual_canonical_Rd_amplitude_body_bound=True,actual_J_one_half=True,
        actual_O2_O3_finite_mu_multiplier_ASTs_executed=True,
        actual_Rc_ordinary_C0_to_C3_source_leaves_bound=True,
        actual_O3_background_cell_AST_sha256=hashlib.sha256(ast.dump(fn,include_attributes=False).encode()).hexdigest(),
        exact_fixed_source_log_scale='-logP/2+1/2',
        exact_finite_coefficient='exp(-17/10-5*mu/2)/(1+Z^2)',
        actual_nonzero_finite_mu_factor_retained=True,source_range_overlap_not_used_as_identity=True)


def actual_arc_coefficient_rows(c,Z,mu):
    """Ordinary physical A_rc derivatives divided by its one fixed source scale."""
    z=c.mpf(Z);den=1+z*z
    if ep(den)[0]<=0:raise ValueError('Split at the original four Z cells before rational evaluation')
    constant=c.exp(-c.mpf(17)/10-c.mpf(5)*mu/2)
    return [dict(coefficient_interval=constant/den,encloses_original_source_function=True),
            dict(coefficient_interval=-2*z*constant/(den*den),encloses_original_source_function=True)]


def divided_matrix_factors(c,values):
    mu=interval(c,values['mu'])
    a,b=[[interval(c,v) for v in row] for row in values['axial_divided_rows']]
    axial=interval(c,values['positive_negative_axial_divided_determinant'])
    angular=interval(c,values['positive_angular_determinant'])
    if ep(mu)[0]<=0 or ep(axial)[0]<=0 or ep(angular)[0]<=0:
        raise ValueError('Actual positive finite-mu divided determinants required')
    T=[[interval(c,v) for v in row] for row in values['angular_rows']]
    inverse=[]
    for i in range(3):
        row=[]
        for j in range(3):
            ii=[k for k in range(3) if k!=j];jj=[k for k in range(3) if k!=i]
            row.append((-1)**(i+j)*(T[ii[0]][jj[0]]*T[ii[1]][jj[1]]
                -T[ii[0]][jj[1]]*T[ii[1]][jj[0]])/angular)
        inverse.append(row)
    # y0 is the first inverse input. Combine the exact numerator BEFORE
    # the physical 1/mu factor is placed in the fixed alpha log unit.
    alpha=[[(a[1]+mu*b[1])/axial,-a[1]/axial],
           [-(a[0]+mu*b[0])/axial,a[0]/axial]]
    return dict(mu=mu,alpha_signed_factors=alpha,angular_inverse=inverse,
        actual_axial_divided_determinant=axial,actual_angular_determinant=angular)


def divide_arc_pair(pair,arc):
    value,Z=pair;a,az=arc
    return [value/a,(Z*a-value*az)/(a*a)]


def signed_cover(c,coefficient,log):
    """A nonzero finite cover, with a sharper absolute log bound retained."""
    lo,hi=ep(coefficient);magnitude=ep(abs(coefficient))[1]
    if magnitude==0:return dict(value=c.mpf(0),absolute_log_upper=None,exact_zero=True)
    bound=log+c.ln(c.mpf(magnitude))
    if ep(bound)[1]>=-1000:raise ValueError('Actual coefficient sector is outside the fixed small cover')
    b=ep(c.exp(-1000))[1]
    value=c.mpf([0,b]) if lo>=0 else c.mpf([-b,0]) if hi<=0 else c.mpf([-b,b])
    return dict(value=value,absolute_log_upper=bound,exact_zero=False)


def physical_log_bound(c,terms,unit_log,order):
    logs=[row['rows'][order]['absolute_log_upper'] for row in terms
          if not row['rows'][order]['exact_zero']]
    if not logs:return None
    upper=max(ep(v)[1] for v in logs)
    return unit_log+c.mpf(upper)+c.ln(len(logs))


@source_precision
def coefficient_cell_enclosures(c,view,matrix,arc_rows,logP,L,logSigma,logMu):
    """Transform same-source signed sectors; no boxed y0/ym subtraction."""
    a,az=[interval(c,r['coefficient_interval']) for r in arc_rows[:2]]
    if ep(a)[0]<=0:raise ValueError('Actual positive A_rc source coefficient required')
    unit_alpha=-logMu-L/2
    unit_beta=view['actual_physical_pressure_unit_logs']['m15']-logP-L/2-logSigma
    unit_beta=interval(c,unit_beta)
    arc=[a,az];alpha=[]
    for factors in matrix['alpha_signed_factors']:
        terms=[];covers=[c.mpf(0),c.mpf(0)]
        for key,factor,shift,divide in (
            ('m11',factors[0],c.mpf(0),False),
            ('m14',factors[1],-logP-L/2-logSigma,True)):
            for term in view['full_physical_velocity_moment_C1_sectors'][key]['terms']:
                pair=[interval(c,v) for v in term['ordinary_Z_derivatives'][:2]]
                if divide:pair=divide_arc_pair(pair,arc)
                pair=[factor*v for v in pair]
                log=interval(c,term['combined_source_log'])+shift
                rows=[signed_cover(c,value,log) for value in pair]
                for order,row in enumerate(rows):covers[order]+=row['value']
                terms.append(dict(moment=key,source_exponents=term['source_exponents'],
                    original_absolute_source_log=interval(c,term['combined_source_log']),
                    exact_physical_log_shift=shift,transformed_ordinary_coefficients=pair,rows=rows))
        alpha.append(dict(normalized_C0_C1_enclosures=covers,physical_log_unit=unit_alpha,
            physical_C0_C1_absolute_log_upper=[physical_log_bound(c,terms,unit_alpha,k) for k in (0,1)],
            signed_original_source_terms=terms,correlated_source_numerator_before_nonzero_mu_division=True))
    pressure=view['normalized_full_pressure_moment_C1_enclosures']
    p13=divide_arc_pair([interval(c,v) for v in pressure['m13']],arc)
    p15=divide_arc_pair([interval(c,v) for v in pressure['m15']],arc)
    # theta1/unit_beta=-exp(L-logRv)*p13/(2*a). The exact Rp recipe
    # cancels logRc before boxing. Do not materialize exp(-13/mu).
    ratio_log=c.ln(2)+2+60*logMu-13/matrix['mu']
    theta1=[signed_cover(c,-v/2,ratio_log) for v in p13]
    theta0=[c.mpf(0),c.mpf(0)];theta0_sources=[]
    for term in view['full_physical_velocity_moment_C1_sectors']['m12']['terms']:
        log=interval(c,term['combined_source_log'])-L-unit_beta
        rows=[signed_cover(c,-interval(c,v),log) for v in term['ordinary_Z_derivatives'][:2]]
        for k,row in enumerate(rows):theta0[k]+=row['value']
        theta0_sources.append(dict(source_exponents=term['source_exponents'],rows=rows,
            original_absolute_source_log=interval(c,term['combined_source_log']),normalized_log_shift=-L-unit_beta))
    theta=[theta0,[r['value'] for r in theta1],p15]
    beta=[]
    for row in matrix['angular_inverse']:
        values=[sum((B*pair[k] for B,pair in zip(row,theta)),c.mpf(0)) for k in (0,1)]
        logs=[]
        for value in values:
            magnitude=ep(abs(value))[1]
            logs.append(None if magnitude==0 else unit_beta+c.ln(c.mpf(magnitude)))
        beta.append(dict(normalized_C0_C1_enclosures=values,physical_log_unit=unit_beta,
            physical_C0_C1_absolute_log_upper=logs,all_signed_theta_inputs_retained=True))
    return dict(coefficients=dict(zip(COEFFICIENTS,alpha+beta)),actual_A_rc_C0_C1_coefficient_enclosures=arc,
        actual_A_rc_shared_fixed_log_scale=logSigma,normalized_theta_inputs=theta,
        nonzero_theta0_original_sectors=theta0_sources,nonzero_theta1_ratio_cover=theta1,
        actual_Astar_Z_in_all_quotient_derivatives=True,
        coefficient_intervals_and_source_units_not_selected_as_function_values=True,
        no_physical_hyper_exponential_materialized=True)


def from_archive(c,value):
    if isinstance(value,dict):
        if 'lower_exact_mpf_tuple' in value and 'upper_exact_mpf_tuple' in value:return interval(c,value)
        return {key:from_archive(c,item) for key,item in value.items()}
    if isinstance(value,list):return [from_archive(c,item) for item in value]
    return value


def exact_reserved_radius_identity(owner):
    native=full.pulse.native;g=owner.g
    frame=owner.full.pulse.native.selected.inlet.exact.frame
    importer=full.true.parent.GraphImport(g,frame.graph)
    parameters=frame.bridge.leading.parameters
    parameter=lambda key:importer(full.true.source.FunctionRef(frame.graph,parameters[key].node))
    mu,logC=sy.symbols('original_mu original_logCstar',positive=True)
    LP=sy.exp(40)+11
    bindings={parameter('mu').node:mu,parameter('logC').node:logC,parameter('logP').node:LP}
    for key in ('T','hbB','hbS','sc'):bindings[parameter(key).node]=sy.Symbol('original_'+key,positive=True)
    L=native.scalar_expression(g,owner.operator.functions['log2Rc'],bindings)
    expected=sy.log(2)+sy.log(110)+10*logC+11*sy.exp(40)+124
    if sy.simplify(L-expected)!=0:raise ValueError('Actual reserved log(2Rc) source recipe differs')
    actual_LP=native.scalar_expression(g,full.true.source.FunctionRef(g,owner.operator.background['parameters']['logP']),{})
    actual_mu=native.scalar_expression(g,owner.operator.functions['mu'],{})
    if sy.simplify(actual_LP-LP)!=0 or sy.simplify(actual_mu-sy.exp(-4*LP-sy.log(1000)))!=0:
        raise ValueError('Actual reserved mu/Pstar scalar source recipes differ')
    desired=owner.scales['logMu']
    for log in (interval(owner.c,owner.operator.background['actual_mu_log']),owner.c.ln(owner.matrix['mu'])):
        if ep(log)[0]>ep(desired)[0] or ep(log)[1]<ep(desired)[1]:
            raise ValueError('Reserved numeric mu/log-mu does not enclose the original scalar recipe')
    actual_LP_box=owner.c.exp(40)+11
    if actual_LP_box._mpi_!=owner.scales['logP']._mpi_:
        raise ValueError('Reserved and pulse numeric Pstar source differs')
    return dict(passed=True,actual_reserved_Rc_and_pulse_Rc_recipe_identical=True,
        actual_graph_mu_and_logP_recipes_and_numeric_enclosures_bound=True,
        phase_origin_hb_sc_terms_cancel_before_absolute_radius_interpretation=True)


@dataclass(frozen=True,eq=False)
class OriginalReservedCoefficientPacket:
    scope:str='actual selected two axial and three angular physical C0/C1 coefficient ranges'


class CurrentOriginalN1ReservedCoefficientEnclosures:
    @source_precision
    def __init__(self,accepted_full,require_checked=True):
        if type(accepted_full) is not full.CurrentOriginalN1FullFiveNumericBridge or not accepted_full.acceptance_loaded:
            raise ValueError('Live checked complete true five numerical defects required')
        self.full=accepted_full;self.g=accepted_full.g;self.c=c=accepted_full.c
        accepted_full.assert_graph();self.operator=accepted_full.full.linear
        self.functions=accepted_full.full.coefficients;self.family=copy.deepcopy(accepted_full.family)
        self.matrix=divided_matrix_factors(c,self.operator.values)
        end=accepted_full.pulse.numeric.pulse.end
        self.scales=dict(logP=end.logP,L=c.ln(2)+end.logRc,logMu=-4*end.logP-c.ln(1000),
            logSigma=-end.logP/2+c.mpf(1)/2)
        a=full.true.parent.C1Algebra(self.g);W=self.operator.functions
        LP=full.true.source.FunctionRef(self.g,self.operator.background['parameters']['logP'])
        half=self.g.constant('1/2');L=W['log2Rc'];lm=self.g.unary('log',W['mu'])
        sigma=self.g.add(self.g.mul(self.g.constant('-1/2'),LP),half)
        logs=dict(alpha=self.g.add(self.g.neg(lm),self.g.mul(self.g.constant('-1/2'),L)),
            beta=self.g.add(self.full.pulse.normalization['physical_unit_logs']['m15'],self.g.neg(LP),
                self.g.mul(self.g.constant('-1/2'),L),self.g.neg(sigma)))
        self.unit_functions=logs
        self.normalized_functions={key:a.div(value,a.fixed(self.g.unary('exp',logs['alpha' if key.startswith('alpha') else 'beta'])),
            'actual positive fixed original coefficient physical unit') for key,value in self.functions.items()}
        start=self.full._normalization_span[1]
        self._span=(start,len(self.g.nodes));self._nodes=copy.deepcopy(self.g.nodes[start:])
        self.bindings=full.true.own.encode(dict(physical=self.functions,normalized=self.normalized_functions,
            actual_fixed_log_unit_functions=self.unit_functions))
        self.archive=full.true.parent.linear.read(full.NAME)
        if self.archive['input_hashes']!=accepted_full.hashes or self.archive['true_function_bindings']!=accepted_full.bindings:
            raise ValueError('Complete numerical archive belongs to another true function source')
        self.hashes=dict(accepted_full.hashes)
        for module in (full,transition,original_O3):
            names=(Path(module.__file__).name,)
            if module is full:names+=(module.NAME,module.RECEIPT)
            for name in names:self.hashes[name]=sha(name)
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.source_identity=self.assert_source_identity()
        self._identity=copy.deepcopy(self.source_identity);self._bindings=copy.deepcopy(self.bindings)
        self._hashes=dict(self.hashes);self._state=record((self.matrix,self.scales,self.archive))
        self._owners=(self.full,self.operator,self.g,self.c);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=full.true.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or archive!=sha(NAME) or receipt.get('source_family')!=self.family
                or receipt.get('source_identity')!=self.source_identity):
                raise ValueError('Checked actual reserved coefficient C1 numerical bounds required')
            self.acceptance_loaded=True

    @source_precision
    def assert_source_identity(self):
        self.full.assert_graph();self.operator.assert_graph();before=len(self.g.nodes)
        coefficients=self.operator.coefficient_functions(self.full.functions)
        if full.true.own.encode(dict(zip(COEFFICIENTS,coefficients)))!=full.true.own.encode(self.functions) or len(self.g.nodes)!=before:
            raise ValueError('Actual selected coefficient function action differs')
        amplitude=exact_original_arc_source_identity(self.operator)
        radius=exact_reserved_radius_identity(self)
        for cell in self.operator.background['amplitude_positive_source_cells']:
            expected=actual_arc_coefficient_rows(self.c,cell['exact_Z_cell'],self.matrix['mu'])
            for row,new in zip(cell['actual_A_rc_ordinary_Z0_Z3'][:2],expected):
                scale=row['formal_positive_scale']
                if scale['source_exponents']!=[0,0,0,0] or scale['radius_power']!=0:
                    raise ValueError('Actual Rc amplitude uses another physical formal basis')
                old=interval(self.c,row['coefficient_interval']);value=interval(self.c,new['coefficient_interval'])
                if ep(old)[0]>ep(value)[0] or ep(old)[1]<ep(value)[1]:
                    raise ValueError('Exact same-source amplitude range exceeds accepted C3 source enclosure')
                sigma=interval(self.c,scale['additional_log_interval'])
                if ep(sigma)[0]>ep(self.scales['logSigma'])[0] or ep(sigma)[1]<ep(self.scales['logSigma'])[1]:
                    raise ValueError('Actual fixed amplitude scale differs from its source enclosure')
        return dict(passed=True,actual_selected_five_coefficient_graph_action_replayed=True,
            actual_arc_source_identity=amplitude,actual_reserved_radius_identity=radius,
            exact_arc_model_contained_in_accepted_ordinary_C1_source_rows=True,
            actual_fixed_arc_scale_contained_in_original_source_bounds=True,
            alpha_fixed_unit='1/(mu*sqrt(2Rc))',
            beta_fixed_unit='reference*A*Rv^2/(Pstar*sqrt(2Rc)*exp(-logP/2+1/2))',
            all_units_Z_independent=True,actual_Astar_Z_retained=True,
            complete_signed_source_ledger_used_before_mu_division=True,
            no_global_inner_power_cancellation_assumed=True,
            coefficient_and_higher_order_repair_closure_not_claimed=True)

    @source_precision
    def assert_graph(self):
        if (self._owners!=(self.full,self.operator,self.g,self.c) or self.hashes!=self._hashes
            or self.source_identity!=self._identity or record((self.matrix,self.scales,self.archive))!=self._state
            or self.g.nodes[self._span[0]:self._span[1]]!=self._nodes
            or full.true.own.encode(dict(physical=self.functions,normalized=self.normalized_functions,
                actual_fixed_log_unit_functions=self.unit_functions))!=self._bindings):
            raise ValueError('Actual coefficient numerical source attachment changed')
        self.full.assert_graph();self.operator.assert_graph()

    @source_precision
    def evaluate(self,Z):
        c=self.c;z=c.mpf(Z)
        if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original real Z domain [-1,1] required')
        self.assert_graph();sample=None
        for name,saved in self.archive['samples'].items():
            if interval(c,saved['Z'])._mpi_==z._mpi_:sample=from_archive(c,saved);source=name;break
        if sample is None:
            sample=self.full.report(self.full.evaluate(Z));source='new checked complete moment query'
        cells=[]
        for row in self.operator.background['amplitude_positive_source_cells']:
            ends=c.mpf(row['exact_Z_cell']);left=max(ep(z)[0],ep(ends)[0]);right=min(ep(z)[1],ep(ends)[1])
            if left>right:continue
            clip=c.mpf([left,right]);rows=actual_arc_coefficient_rows(c,clip,self.matrix['mu'])
            view=coefficient_cell_enclosures(c,sample,self.matrix,rows,**self.scales)
            cells.append(dict(original_source_Z_cell=row['exact_Z_cell'],clipped_Z=clip,**view))
        if not cells:raise ValueError('No actual amplitude source cell covers this query')
        result={}
        for key in COEFFICIENTS:
            pieces=[row['coefficients'][key] for row in cells]
            pair=[c.mpf([min(ep(p['normalized_C0_C1_enclosures'][k])[0] for p in pieces),
                         max(ep(p['normalized_C0_C1_enclosures'][k])[1] for p in pieces)]) for k in (0,1)]
            logs=[]
            for k in (0,1):
                candidates=[p['physical_C0_C1_absolute_log_upper'][k] for p in pieces
                    if p['physical_C0_C1_absolute_log_upper'][k] is not None]
                logs.append(None if not candidates else c.mpf(max(ep(v)[1] for v in candidates)))
            result[key]=dict(normalized_C0_C1_enclosures=pair,physical_log_unit=pieces[0]['physical_log_unit'],
                physical_C0_C1_absolute_log_upper=logs)
        view=dict(Z=z,source_family=self.family,true_function_bindings=self.bindings,
            coefficients=result,cells=cells,complete_moment_numeric_source=source,
            actual_five_coefficient_C0_C1_enclosures_present=True,actual_Astar_Z_retained=True,
            all_nonzero_source_terms_and_finite_mu_preserved=True,
            physical_coefficients_not_replaced_by_interval_centers=True,
            repaired_own_terminal_zero_identities_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalReservedCoefficientPacket();self._packets[id(packet)]=(packet,view,record(view));return packet

    @source_precision
    def report(self,packet):
        entry=self._packets.get(id(packet))
        if type(packet) is not OriginalReservedCoefficientPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live actual reserved coefficient packet required')
        self.assert_graph()
        if record(entry[1])!=entry[2]:raise ValueError('Actual coefficient packet changed')
        view=full.pulse.caps.end.copy_view(self.c,entry[1]);view[GATE]=self.acceptance_loaded;return view


@source_precision
def run(accepted_full):
    began=time.monotonic();owner=CurrentOriginalN1ReservedCoefficientEnclosures(accepted_full,False)
    packets={}
    for name,Z in (('point','.371'),('whole_Z',['-1','1'])):
        print('ACTUAL_RESERVED_COEFFICIENT_C1_START',name,flush=True)
        packets[name]=owner.evaluate(Z)
        print('ACTUAL_RESERVED_COEFFICIENT_C1_QUERY_DONE',name,flush=True)
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,
        true_function_bindings=owner.bindings,samples={name:owner.report(p) for name,p in packets.items()},
        coefficient_normalization_graph_span=owner._span,coefficient_normalization_graph_nodes=owner._nodes,
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,raw,packets
