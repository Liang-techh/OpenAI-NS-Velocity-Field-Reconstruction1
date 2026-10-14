"""Complete true signed n1 defects with directed ordinary C0/C1 enclosures.

Pressure components share the true pulse physical unit. Velocity moments
retain their physical log sectors. No original physical exponential is
materialized and no coefficient or terminal-zero claim is made here.
"""
import copy
import ast
import inspect
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
import textwrap
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_pulse_numeric_bridge as pulse
import lei_ren_part1_paper_compliant_current_original_n1_preRp_numeric_bridge as preRp
import lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric_bridge as cutoff
import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as old
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

lower=preRp.lower
true=pulse.true
HERE,PREFIX,sha=pulse.HERE,pulse.PREFIX,pulse.sha
ep=pulse.ep
NAME=PREFIX+'current_original_n1_full_five_numeric_bridge.json.gz'
RECEIPT=PREFIX+'current_original_n1_full_five_numeric_bridge_check.json'
GATE='current_original_true_complete_five_defect_C1_numeric_bridge_installed'
OPEN=pulse.OPEN
VELOCITY=('m11','m12','m14')
PRESSURE=('m13','m15')


def record(value):
    """Preserve scalar endpoints as exact tuples as well as interval records."""
    if hasattr(value,'_mpf_'):
        return dict(lower_exact_mpf_tuple=list(value._mpf_),upper_exact_mpf_tuple=list(value._mpf_),
            exact_scalar_endpoint=True)
    if isinstance(value,dict):return {key:record(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [record(item) for item in value]
    return pulse.record(value)


def interval(c,value):
    if hasattr(value,'_mpi_'):return c.make_mpf(value._mpi_)
    return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),tuple(value['upper_exact_mpf_tuple'])))


def exact_pressure_unit_source_identity(owner):
    """Lift actual numeric recipes and compare the true graph scale handles."""
    c=SimpleNamespace(mpf=lambda v:sy.Rational(v) if isinstance(v,str) else sy.sympify(v),
        exp=sy.exp,ln=lambda v:sy.expand_log(sy.log(v),force=True),sqrt=sy.sqrt,pi=sy.pi)
    logC,lu=sy.symbols('original_logCstar log_native_U0',real=True)
    lp=sy.exp(40)+11;lm=-4*lp-sy.log(1000);mu=sy.exp(lm)
    native=pulse.native
    pref=native.execute_assignments(pulse.original_map.pulse_rows,('k0','u0','L','logpref'),
        dict(c=c,mu=mu))['logpref']
    logE=native.execute_assignments(pulse.original_map.SharedOuterPulseMap.__init__,('self.logscale',),
        dict(c=c,self=SimpleNamespace(mu=mu,rows={'common_logpref':pref})))['self'].logscale
    end=pulse.caps.end
    formulas=native.execute_assignments(end.CurrentOriginalN1EndPressureSource.__init__,
        ('self.logRc','self.logRp','self.logP','self.logUp','logE','self.logA','self.logRv'),
        dict(c=c,logC=logC,self=SimpleNamespace(mu=mu,coefficients=SimpleNamespace(selected=SimpleNamespace(
            constants={'U':sy.exp(lu)},amplitude=SimpleNamespace(log_end_scale=logE))))),
        {'rebox':lambda ctx,value:value})['self']
    tree=ast.parse(textwrap.dedent(inspect.getsource(old.pressure_unit_log)))
    namespace=dict(old.pressure_unit_log.__globals__);namespace['interval']=lambda ctx,value:value
    exec(compile(tree,'<actual-combined-pressure-unit-recipe>','exec'),namespace)
    lift=namespace['pressure_unit_log'];g=owner.g
    Gamma=owner.pulse.pulse.own.jets.complete.Gamma;W=Gamma.pulse.functions
    bindings={true.source.FunctionRef(g,Gamma.pulse.preRp.operator.background['parameters']['logP']).node:formulas.logP,
        owner.pulse.normalization['log_native_U0'].node:formulas.logUp,
        Gamma.schedule['logRp'].node:formulas.logRp,W['log_end_scale'].node:logE,
        Gamma.schedule['mu'].node:mu}
    cache=dict(bindings)
    def scalar(node):
        if node in cache:return cache[node]
        row=g.nodes[node];op=row['operation']
        if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op=='mathematical_pi':value=sy.pi
        elif op=='sum':value=sum((scalar(i) for i in row['arguments']),sy.Integer(0))
        elif op=='product':value=sy.prod(scalar(i) for i in row['arguments'])
        elif op=='negative':value=-scalar(row['argument'])
        elif op=='positive_quotient':value=scalar(row['numerator'])/scalar(row['denominator'])
        elif op=='analytic_unary' and row['name'] in ('exp','log','sqrt'):
            value=getattr(sy,row['name'])(scalar(row['argument']))
        else:raise ValueError('Unsupported physical unit scalar operation '+op)
        cache[node]=value;return value
    for key,power in (('m13',1),('m15',2)):
        actual=sy.expand_log(scalar(owner.pulse.normalization['physical_unit_logs'][key].node),force=True)
        combined=lift(c,dict(logP=formulas.logP,logUp=formulas.logUp,logRc=formulas.logRc),power)
        direct=c.ln(c.mpf('10.99')*c.sqrt(2*c.pi)*(c.exp(5)-c.exp(3)))+formulas.logA+power*formulas.logRv
        if sy.simplify(sy.expand(actual-combined))!=0 or sy.simplify(sy.expand(direct-combined))!=0:
            raise ValueError('Actual true and numerical physical pressure unit recipes differ')
    physical=owner.pulse.numeric.pulse.end;ctx=owner.c
    if physical.coefficients.selected is not owner.pulse.native.selected:
        raise ValueError('Pressure scales use another actual native scalar owner')
    lc=end.coefficients.basis.read_interval(ctx,physical.radius_binding['selected_logCstar'])
    expected=dict(logRc=ctx.ln(110)+10*lc+11*ctx.exp(40)+124,logP=ctx.exp(40)+11,
        logUp=ctx.ln(interval(ctx,physical.coefficients.selected.constants['U'])))
    for key,value in expected.items():
        if value._mpi_!=getattr(physical,key)._mpi_:raise ValueError('Physical scalar scale recipe changed '+key)
    inner_receipt=true.parent.linear.read(cutoff.inner.RECEIPT)
    amplitude=inner_receipt['physical_conjugation_proof']
    if (not inner_receipt['all_passed'] or not inner_receipt[cutoff.inner.GATE] or not amplitude['passed']
        or amplitude['original_relative_amplitude_ordinary_Z0_to_Z3_rows']!=16
        or not amplitude['single_original_F0_Z_gauge_verified']):
        raise ValueError('Accepted actual amplitude Taylor derivative proof required')
    return dict(passed=True,actual_physical_pressure_unit_graph_identities=2,
        actual_direct_and_combined_log_recipe_identities=2,
        actual_native_U0_and_selected_logCstar_scalar_constructor_recipes_bound=True,
        ordinary_pressure_unit_Z1_exact_zero=True,
        accepted_actual_relative_amplitude_ordinary_Z0_to_Z3_rows=16,
        source_bases_fixed_and_amplitude_Z_derivatives_in_coefficients=True)


def sector_C1_in_fixed_pressure_unit(c,sectors,unit_log):
    """Enclose already differentiated physical sectors in a Z-independent unit."""
    if sectors['ordinary_Z_derivative_order']<1:raise ValueError('Actual ordinary C1 sector rows required')
    values=[c.mpf(0),c.mpf(0)];sources=[];beta=c.exp(-1000)
    for term in sectors['terms']:
        log=interval(c,term['combined_source_log']);rows=[]
        for order in (0,1):
            coefficient=interval(c,term['ordinary_Z_derivatives'][order]);lo,hi=ep(coefficient)
            magnitude=ep(abs(coefficient))[1]
            if magnitude==0:
                value=c.mpf(0);bound=None
            else:
                bound=log-unit_log+c.ln(c.mpf(magnitude))
                if ep(bound)[1]>=-1000:raise ValueError('Actual sector exceeds the original pressure-unit cover')
                b=ep(beta)[1]
                value=c.mpf([0,b]) if lo>=0 else c.mpf([-b,0]) if hi<=0 else c.mpf([-b,b])
            values[order]+=value
            rows.append(dict(ordinary_Z_order=order,source_coefficient=coefficient,
                coefficient_exact_zero_in_this_query=magnitude==0,
                normalized_log_absolute_upper=bound,normalized_enclosure=value))
        sources.append(dict(source_exponents=term['source_exponents'],combined_source_log=log,rows=rows))
    return dict(C1_enclosures=values,physical_sector_sources=sources,
        original_source_bases_are_fixed_within_Taylor_jets=True,
        source_ordinary_Z1_already_contains_F0_amplitude_derivatives=True,
        positive_unit_is_Z_independent=True,nonzero_sector_not_selected_zero=True)


@dataclass(frozen=True,eq=False)
class OriginalTrueFullFiveNumericPacket:
    scope:str='complete unrepaired true signed five n1 moment C0/C1 defects'


class CurrentOriginalN1FullFiveNumericBridge:
    @source_precision
    def __init__(self,accepted_cutoff,accepted_lower,accepted_preRp,accepted_pulse,require_checked=True):
        specs=((accepted_cutoff,cutoff.CurrentOriginalN1CutoffNumericBridge),
            (accepted_lower,lower.CurrentOriginalN1LowerNumericBridge),
            (accepted_preRp,preRp.CurrentOriginalN1PreRpNumericBridge),
            (accepted_pulse,pulse.CurrentOriginalN1PulseNumericBridge))
        if any(type(o) is not cls or not o.acceptance_loaded for o,cls in specs):
            raise ValueError('Four live checked true numerical component owners required')
        self.cut,self.lower,self.preRp,self.pulse=(o for o,cls in specs)
        self.full=self.cut.inner.full;self.g=self.full.g;self.c=self.pulse.c
        self.family=copy.deepcopy(self.full.family)
        if (self.lower.cutoff is not self.cut or self.preRp.lower is not self.lower
            or self.lower.lower is not self.full.lower or self.preRp.preRp is not self.full.preRp
            or self.pulse.pulse is not self.full.pulse or any(o.family!=self.family for o,cls in specs)
            or any(o.g is not self.g for o,cls in specs)):
            raise ValueError('Components must belong to one actual full original function graph')
        self.functions=self.full.moments
        a=true.parent.C1Algebra(self.g)
        self.normalized_pressure_functions={key:a.div(self.functions[key],
            a.fixed(self.g.unary('exp',self.pulse.normalization['physical_unit_logs'][key])),
            'same positive reference*A*Rv power, Z independent') for key in PRESSURE}
        start=self.pulse._normalization_span[1]
        self._normalization_span=(start,len(self.g.nodes))
        self._normalization_nodes=copy.deepcopy(self.g.nodes[start:])
        self.bindings=true.own.encode(dict(physical_five=self.functions,
            normalized_pressure=self.normalized_pressure_functions))
        end=self.pulse.numeric.pulse.end
        self.scales=dict(logP=end.logP,logUp=end.logUp,logRc=end.logRc)
        self.unit_logs={key:old.pressure_unit_log(self.c,self.scales,p) for key,p in (('m13',1),('m15',2))}
        self.hashes={}
        for owner,cls in specs:
            for name,digest in owner.hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting complete moment source '+name)
                self.hashes[name]=digest
        for module in (cutoff,lower,preRp,pulse,old):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed complete moment source '+name)
        self.source_identity=self.assert_source_identity()
        self._identity=copy.deepcopy(self.source_identity);self._bindings=copy.deepcopy(self.bindings)
        self._hashes=dict(self.hashes);self._state=record((self.scales,self.unit_logs))
        self._owners=(self.cut,self.lower,self.preRp,self.pulse,self.full,self.g,self.c)
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=true.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or archive!=sha(NAME) or receipt.get('source_family')!=self.family
                or receipt.get('source_identity')!=self.source_identity):
                raise ValueError('Checked complete true five C1 numerical defects required')
            self.acceptance_loaded=True

    @source_precision
    def assert_source_identity(self):
        for owner in (self.cut,self.lower,self.preRp,self.pulse):owner.assert_graph()
        before=len(self.g.nodes)
        replay=cutoff.full.compose_full_moments(self.g,self.full.cutoff.endpoint,self.full.lower.functions,
            self.full.preRp.functions,self.full.pulse.functions)
        if true.own.encode(replay)!=true.own.encode(self.functions) or len(self.g.nodes)!=before:
            raise ValueError('Actual complete five moment function partition differs')
        units=self.preRp.source_identity['original_physical_units']
        pulse_units=self.pulse.source_identity['original_physical_units']
        for key in ('source_Rp_node','source_native_frame','source_A','source_Rv_formula','source_logRp_formula',
            'source_logE_formula','normalization','m13','m15'):
            if units[key]!=pulse_units[key]:raise ValueError('Pre-Rp and pulse physical normalization differs '+key)
        if not self.pulse.source_identity['actual_physical_scale_and_incoming_endpoint_identities']['passed']:
            raise ValueError('Actual physical scale identities required')
        zero=self.full.beyond_Rv_support
        support=cutoff.full.support
        if (zero['endpoint']['definition']!='Rv=exp(original_logRp+13/original_mu)'
            or not zero['endpoint']['angular_and_pressure_not_zeroed']
            or set(zero['future_charts'])!=set(support.FUTURE)
            or any(row['exact_zero_unweighted_and_R_weighted_forcing_integrals']!=[0,0]
                   for row in zero['future_charts'].values())):
            raise ValueError('Actual beyond-Rv pressure support theorem required')
        return dict(passed=True,one_same_true_graph_and_live_source_owner_chain=True,
            actual_pressure_unit_source_identity=exact_pressure_unit_source_identity(self),
            exact_five_signed_function_composition_replayed=True,
            inner_lower_preRp_pulse_pressure_partitions_complete=True,
            beyond_Rv_true_pressure_forcing_zero=self.full.beyond_Rv_support,
            velocity_m11_m12_m14_cover_full_support_at_Rin=True,
            preRp_pulse_same_original_positive_reference_A_Rv_units=True,
            physical_pressure_unit_logs_use_combined_inverse_mu_powers=True,
            factored_source_Z1_already_includes_F0_relative_amplitude_derivatives=True,
            pressure_unit_Z1_exact_zero=True,
            m15_minus_already_in_original_signed_components=True,
            actual_linear_coefficients_and_repaired_own_fields_not_yet_enclosed=True)

    @source_precision
    def assert_graph(self):
        if (self._owners!=(self.cut,self.lower,self.preRp,self.pulse,self.full,self.g,self.c)
            or self.hashes!=self._hashes or self.source_identity!=self._identity
            or self.g.nodes[self._normalization_span[0]:self._normalization_span[1]]!=self._normalization_nodes
            or true.own.encode(dict(physical_five=self.functions,normalized_pressure=self.normalized_pressure_functions))!=self._bindings
            or record((self.scales,self.unit_logs))!=self._state):
            raise ValueError('Actual complete numerical moment attachment changed')
        end=self.pulse.numeric.pulse.end
        if record(self.scales)!=record(dict(logP=end.logP,logUp=end.logUp,logRc=end.logRc)):
            raise ValueError('Actual physical scale source changed')
        for owner in (self.cut,self.lower,self.preRp,self.pulse):owner.assert_graph()
        self.full.assert_graph()

    @source_precision
    def evaluate(self,Z,reuse=None):
        self.assert_graph();reuse=reuse or {};c=self.c;z=c.mpf(Z)
        if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original real Z domain [-1,1] required')
        packets={}
        packets['inner']=reuse.get('inner') or self.cut.evaluate(Z,'1','first',cells=4)
        packets['lower']=reuse.get('lower') or self.lower.evaluate(Z,Rin_packet=packets['inner'])
        packets['preRp']=reuse.get('preRp') or self.preRp.evaluate(Z)
        packets['pulse']=reuse.get('pulse') or self.pulse.evaluate(Z)
        owners=dict(inner=self.cut,lower=self.lower,preRp=self.preRp,pulse=self.pulse)
        views={name:owners[name].report(packet) for name,packet in packets.items()}
        for name,view in views.items():
            if view['Z']._mpi_!=owners[name].c.mpf(Z)._mpi_:
                raise ValueError('A component packet belongs to another original Z query')
        inner=views['inner']
        if inner['chart']!='first' or inner['coordinate']._mpi_!=self.cut.c.mpf(1)._mpi_:
            raise ValueError('Actual full-support true Rin moment packet required')
        sectors=inner['physical_five_inner_cumulative_Z1']
        pressure={};parts={}
        for key in PRESSURE:
            inner_part=sector_C1_in_fixed_pressure_unit(c,sectors[key],self.unit_logs[key])
            lower_part=sector_C1_in_fixed_pressure_unit(c,views['lower']['true_lower_pressure_C1_enclosures'][key],self.unit_logs[key])
            other=[interval(c,v) for v in views['preRp']['normalized_true_preRp_C1_enclosures'][key]]
            last=[interval(c,v) for v in views['pulse']['normalized_true_pulse_C1_enclosures'][key]]
            pressure[key]=[inner_part['C1_enclosures'][k]+lower_part['C1_enclosures'][k]+other[k]+last[k] for k in (0,1)]
            parts[key]=dict(inner=inner_part,lower=lower_part,preRp=other,pulse=last)
        view=dict(source_family=self.family,Z=z,moment_order=['m11','m12','m13','m14','m15'],
            true_function_bindings=self.bindings,
            full_physical_velocity_moment_C1_sectors={key:sectors[key] for key in VELOCITY},
            full_velocity_moment_source_budgets={key:inner['five_inner_cumulative_budgets'][key] for key in VELOCITY},
            normalized_full_pressure_moment_C1_enclosures=pressure,pressure_component_enclosures=parts,
            actual_physical_pressure_unit_logs=self.unit_logs,pressure_normalization=self.pulse.source_identity['original_physical_units'],
            all_five_true_unrepaired_moment_C1_enclosures_present=True,all_pressure_radial_support_included=True,
            velocity_moments_remain_factored_physical_not_pressure_units=True,
            signed_m15_not_negated_again=True,nonzero_source_errors_preserved=True,
            actual_coefficients_and_terminal_zero_identities_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalTrueFullFiveNumericPacket()
        self._packets[id(packet)]=(packet,view,packets,record(view));return packet

    @source_precision
    def report(self,packet):
        entry=self._packets.get(id(packet))
        if type(packet) is not OriginalTrueFullFiveNumericPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live true complete five C1 numerical packet required')
        self.assert_graph()
        if record(entry[1])!=entry[3]:raise ValueError('Changed true complete five C1 numerical packet')
        view=pulse.caps.end.copy_view(self.c,entry[1]);view[GATE]=self.acceptance_loaded;return view


@source_precision
def run(accepted_cutoff,accepted_lower,accepted_preRp,accepted_pulse,reuse=None):
    began=time.monotonic();owner=CurrentOriginalN1FullFiveNumericBridge(accepted_cutoff,accepted_lower,accepted_preRp,accepted_pulse,False)
    packets={}
    for name,Z in (('point','.371'),('whole_Z',['-1','1'])):
        print('TRUE_COMPLETE_FIVE_C1_START',name,flush=True)
        packets[name]=owner.evaluate(Z,None if reuse is None else reuse.get(name))
        print('TRUE_COMPLETE_FIVE_C1_QUERY_DONE',name,flush=True)
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,true_function_bindings=owner.bindings,
        normalization_graph_span=owner._normalization_span,normalization_graph_nodes=owner._normalization_nodes,
        samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,raw,packets
