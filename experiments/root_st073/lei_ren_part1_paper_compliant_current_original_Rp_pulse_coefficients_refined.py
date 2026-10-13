"""Refine the actual selected pulse coefficients without changing the branch.

The original current selected owner supplies genuine C5 incoming and full
future functions. Refined enclosures of the same fixed pulse/bump integrals
feed the original divided-difference inverse and positive quadratic equation.
No ancestor cache, source constant, field definition or artifact is replaced.
"""
import copy
import ast
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_basis_energy as basis
import lei_ren_part1_paper_compliant_current_original_Rp_selected_pulse as selected
import lei_ren_part1_paper_compliant_axial_high_jets as implicit
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE,PREFIX,sha,ep = basis.HERE,basis.PREFIX,basis.sha,basis.ep
source_precision=basis._source_precision
NAME=PREFIX+'current_original_Rp_pulse_coefficients_refined.json.gz'
RECEIPT=PREFIX+'current_original_Rp_pulse_coefficients_refined_check.json'
GATE='current_original_Rp_selected_coefficients_refined_installed'
OPEN=tuple(k for k in basis.OPEN if k not in ('current_original_Rp_selected_amplitude_installed',
    'current_original_Rp_end_coefficients_installed'))+('uniform_coefficient_derivative_accuracy_certified',)
FIELDS=('ap','C1','C2')


def rebox(c,value):
    if isinstance(value,IntervalTaylor):
        return IntervalTaylor(c,[rebox(c,v) for v in value.coefficients])
    if hasattr(value,'_mpi_'):return c.make_mpf(value._mpi_)
    return c.mpf(value)


def record(value):
    if isinstance(value,IntervalTaylor):
        return dict(Taylor_order=value.order,coefficients=[basis.record(v) for v in value.coefficients])
    if isinstance(value,dict):return {k:record(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [record(v) for v in value]
    return basis.record(value)


def function_binding(fn):
    """Pin the actual callable and its defining code, including decorators."""
    function = getattr(fn, '__func__', fn)
    source = textwrap.dedent(inspect.getsource(function))
    tree = ast.dump(ast.parse(source), include_attributes=False)
    path = Path(inspect.getsourcefile(function)).resolve()
    if path.parent != HERE.resolve():
        raise ValueError('Pulse selection function outside accepted source directory')
    return dict(module=function.__module__, qualname=function.__qualname__,
                source_file=path.name, source_sha256=sha(path.name),
                source_ast_sha256=hashlib.sha256(tree.encode()).hexdigest())


def intersect(c,left,right):
    lo=max(ep(left)[0],ep(right)[0]);hi=min(ep(left)[1],ep(right)[1])
    if lo>hi:raise ValueError('Refined coefficients conflict with original source enclosure')
    return c.mpf([lo,hi])


def intersect_jet(c,left,right):
    if left.order!=right.order:raise ValueError('Same original Taylor order required')
    return IntervalTaylor(c,[intersect(c,a,b) for a,b in zip(left.coefficients,right.coefficients)])


def linear_inverse(A,D,det,r1,r2):
    return [(r1*D[1]-r2*A[1])/det,(r2*A[0]-r1*D[0])/det]


def positive_selection(c,K,nu,u,v,target):
    """Original completed-square C0 and implicit derivative recurrence."""
    if any(isinstance(item,IntervalTaylor) for item in [K]+list(nu)+list(v)):
        raise ValueError('Original K, nu and v must be Z-independent scalar intervals')
    if len(nu)!=2 or len(u)!=2 or len(v)!=2 or not isinstance(target,IntervalTaylor):
        raise ValueError('Original two-coefficient quadratic source required')
    if any(not isinstance(j,IntervalTaylor) or j.order!=target.order for j in u):
        raise ValueError('Matching genuine incoming and target Taylor orders required')
    A2=K;A1=target*0;A0=-target
    for weight,uj,vj in zip(nu,u,v):
        A2+=weight*vj*vj;A1+=uj*(2*weight*vj);A0+=uj*uj*weight
    if ep(A2)[0]<=0 or ep(A0[0])[1]>=0:raise ValueError('Original positive root conditions lost')
    half=A1[0]/(2*A2);disc=half*half-A0[0]/A2
    if ep(disc)[0]<=0:raise ValueError('Positive completed-square discriminant required')
    a0=-half+c.sqrt(disc)
    if ep(a0)[0]<=ep(c.mpf('.9'))[1] or ep(a0)[1]>=ep(c.mpf('1.2'))[0]:
        raise ValueError('Original positive branch in (.9,1.2) required')
    ap,denominator=implicit.positive_quadratic_jets(c,a0,A2,A1,A0)
    C=[uj+ap*vj for uj,vj in zip(u,v)]
    return ap,C,dict(A2=A2,A1=A1,A0=A0,positive_derivative_denominator=denominator)


def value_budget(c,jet):
    value=jet[0];lo,hi=ep(value)
    if lo<=0<=hi:return dict(sign='zero_or_unresolved',relative_width=None,value_accuracy_certified=False)
    norm=c.mpf([min(abs(lo),abs(hi)),min(abs(lo),abs(hi))])
    width=(value.b-value.a)/norm
    return dict(sign='negative' if hi<0 else 'positive',relative_width=width,
        value_accuracy_target='.001',value_accuracy_certified=ep(width)[1]<ep(c.mpf('.001'))[0])


class CurrentOriginalRpPulseCoefficientsRefined:
    @source_precision
    def __init__(self,selected_owner=None,basis_owner=None,require_checked=True):
        self.selected=selected_owner if selected_owner is not None else selected.CurrentOriginalRpSelectedPulse()
        self.basis=basis_owner if basis_owner is not None else basis.CurrentOriginalRpPulseBasisEnergy()
        if type(self.selected) is not selected.CurrentOriginalRpSelectedPulse or not self.selected.acceptance_loaded:
            raise ValueError('Accepted current typed selected pulse required')
        if type(self.basis) is not basis.CurrentOriginalRpPulseBasisEnergy or not self.basis.acceptance_loaded:
            raise ValueError('Accepted refined original basis/energy required')
        self.selected.assert_graph();self.basis.assert_graph()
        self.family=copy.deepcopy(self.selected.family_record)
        if self.family!=self.basis.family:raise ValueError('Same actual source family required')
        self.c=self.basis.ctx
        self.mu=intersect(self.c,rebox(self.c,self.selected.pulse.mu),rebox(self.c,self.basis.mu))
        self.hashes=dict(self.selected.hashes)
        for source in (self.basis,self.basis.saddle_owner):
            for name,digest in source.hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting source '+name)
                self.hashes[name]=digest
        for module in (selected,basis,basis.saddle):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.source_bindings={name:function_binding(fn) for name,fn in (
            ('original_linear_inverse',self.selected.amplitude.linear_inverse),
            ('original_positive_amplitude',self.selected.amplitude.select),
            ('original_current_C5_selection',self.selected.fifth.select),
            ('original_current_C5_future',self.selected.fifth.future),
            ('implicit_quadratic_jets',implicit.positive_quadratic_jets),
            ('current_selected_owner',selected.CurrentOriginalRpSelectedPulse.__init__),
            ('refined_linear_inverse',linear_inverse),
            ('refined_positive_selection',positive_selection))}
        for binding in self.source_bindings.values():
            self.hashes[binding['source_file']]=binding['source_sha256']
        self.acceptance_loaded=False;self._packets={}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Refined coefficient receipt/scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed refined coefficient source '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._owners=(self.selected,self.basis,self.c);self._family=copy.deepcopy(self.family)
        self._hashes=copy.deepcopy(self.hashes)
        self._bindings=copy.deepcopy(self.source_bindings)

    @source_precision
    def assert_graph(self):
        if (self._owners!=(self.selected,self.basis,self.c) or self.family!=self._family or
                self.hashes!=self._hashes or self.source_bindings!=self._bindings):
            raise ValueError('Refined selected coefficient graph changed')
        self.selected.assert_graph();self.basis.assert_graph()

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();began=time.monotonic();c=self.c;z=c.mpf(Z)
        if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original real Z domain required')
        # These are the current copied owners and genuine original functions,
        # including all five derivatives and the complete future source.
        old=self.selected.fifth.select(rebox(self.selected.ctx,z))
        amplitude=self.selected.amplitude
        incoming=[rebox(c,j) for j in old['incoming']['moment_Taylor']]
        energy=rebox(c,old['incoming']['energy_Taylor'])
        future=rebox(c,old['actual_weighted_future_energy_Taylor'])
        if any(j.order!=5 for j in incoming+[energy,future]):raise ValueError('Genuine source C5 required')
        factors=[rebox(c,v) for v in amplitude.incoming_factor_caps]
        rows=[j*f for j,f in zip(incoming,factors)]
        A=[rebox(c,v) for v in self.basis.basis['first_row']]
        D=[rebox(c,v) for v in self.basis.basis['exact_divided_difference_row']]
        det=rebox(c,self.basis.basis['divided_difference_determinant'])
        Pi=[rebox(c,v['refined']) for v in self.basis.saddle_owner.refined_rows]
        u=linear_inverse(A,D,det,-rows[0]*self.mu,-(rows[1]-rows[0]))
        v=linear_inverse(A,D,det,-self.mu*Pi[0],-(Pi[1]-Pi[0]))
        # Old positive energy caps remain enclosures of the same exact nu_j.
        # They are not replaced by zero when their physical source is tiny.
        nu=[rebox(c,w) for w in amplitude.nu]
        if any(ep(w)[1]<=0 for w in nu):raise ValueError('Original nonzero end energy must be retained')
        target=future-energy*self.mu+rebox(c,amplitude.base)
        ap,controls,quadratic=positive_selection(c,rebox(c,self.basis.Kpulse),nu,u,v,target)
        prior_ap=rebox(c,old['selected_ap_Taylor'])
        prior_C=[rebox(c,j) for j in old['selected_scaled_end_coefficient_Taylor']]
        ap=intersect_jet(c,ap,prior_ap)
        controls=[intersect_jet(c,j,p) for j,p in zip(controls,prior_C)]
        jets=dict(zip(FIELDS,[ap]+controls));budgets={n:value_budget(c,j) for n,j in jets.items()}
        if budgets['ap']['sign']!='positive' or budgets['C1']['sign']!='negative' or budgets['C2']['sign']!='positive':
            raise ValueError('Original amplitude and end-coefficient signs lost')
        constraints=dict(first_linear=sum((j*a for j,a in zip(controls,A)),rows[0]*0)+(rows[0]+ap*Pi[0])*self.mu,
            divided_linear=sum((j*d for j,d in zip(controls,D)),rows[0]*0)+(rows[1]-rows[0])+ap*(Pi[1]-Pi[0]),
            energy=ap*ap*rebox(c,self.basis.Kpulse)+sum((j*j*w for j,w in zip(controls,nu)),target*0)-target)
        view=dict(source_family=self.family,Z=z,refined_coefficients=record(jets),
            value_budgets=record(budgets),original_coefficient_intersections_used=True,
            original_genuine_C5_incoming_and_future_retained=True,
            implicit_C5_recovered_from_source_equation_not_root_iteration=True,
            original_positive_branch_preserved=True,original_source_end_scale=record(rebox(c,amplitude.log_end_scale)),
            exact_energy_equation=old['exact_selected_energy_equation'],
            source_linear_constraints='A*C=-mu*(incoming1+ap*F1); D*C=-(incoming2-incoming1+ap*(F2-F1))',
            source_bindings=copy.deepcopy(self.source_bindings),
            source_inputs=record(dict(incoming_moment_Taylor=incoming,
                incoming_energy_Taylor=energy,weighted_complete_future_Taylor=future,
                fixed_incoming_row_factor_caps=factors,
                fixed_incoming_row_factor_log_definitions=[rebox(c,v) for v in amplitude.incoming_factor_log_definitions],
                positive_incoming_row_factor_log_margins=[rebox(c,v) for v in amplitude.incoming_factor_margins],
                positive_end_energy_log_margins=[rebox(c,v) for v in amplitude.energy_log_margins],
                end_energy_definition='nu_j=mu*Kj*exp(2*log_end_scale)>0; finite caps are enclosures',
                original_incoming_shapes=old['incoming']['exact_shapes'])),
            nonzero_energy_caps_retained=record(nu),quadratic=record(quadratic),
            source_constraint_enclosures=record(constraints),
            prior_value_budgets=record({n:value_budget(c,j) for n,j in zip(FIELDS,[prior_ap]+prior_C)}),
            source_cache_and_ancestor_artifacts_not_replaced=True,
            derivative_ranges_available_through_order=5,derivative_accuracy_not_claimed=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=object();state=record((view,jets,constraints))
        self._packets[id(packet)]=(packet,view,jets,constraints,state)
        self.assert_graph()
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();item=self._packets.get(id(packet))
        if item is None or item[0] is not packet:raise ValueError('Live refined coefficient packet required')
        if record((item[1],item[2],item[3]))!=item[4]:
            raise ValueError('Refined coefficient packet changed')
        return copy.deepcopy(item[1])


@source_precision
def run(selected_owner=None,basis_owner=None):
    owner=CurrentOriginalRpPulseCoefficientsRefined(selected_owner,basis_owner,False)
    packets={}
    for name,z in (('negative_outer',['-1','-.5']),('negative_inner',['-.5','0']),
            ('positive_inner',['0','.5']),('positive_outer',['.5','1']),('named_real_point','.371')):
        packets[name]=owner.evaluate(z)
        print('REFINED_ORIGINAL_PULSE_COEFFICIENT_QUERY',name,flush=True)
    raw=dict(source_family=owner.family,samples={n:owner.report(p) for n,p in packets.items()},
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_REFINED_PULSE_COEFFICIENTS_COMPUTED',len(packets),flush=True)
    return owner,packets
