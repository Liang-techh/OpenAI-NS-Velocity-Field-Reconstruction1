"""Original correlated shear jets and live Section11 q enclosures.

Logarithmic velocity rows are exposed from original source locals, before
amplitude/width caps. Their analytic identity a=1-2*logU_y sharpens the
previous independent C/E interval quotient. The same signed inertial roots
and positive source theorems remain. q is an interval function backend,
including unresolved branch boxes, before phase inversion and integrals.
"""
import ast
import copy
import dataclasses
import json
import math
from pathlib import Path
from types import FunctionType
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_signed_input_enclosures as prior
import lei_ren_part1_paper_compliant_current_bridge_stress_operator as bridge_ops
import lei_ren_part1_paper_compliant_current_microswitch_stress_operator as micro_ops
import lei_ren_part1_paper_compliant_microswitch_mixed_C4 as micro_original
import lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 as bridge_original
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import copy_jet,SourceAST

packets=prior.packets;native=prior.native;HERE,PREFIX,sha=prior.HERE,prior.PREFIX,prior.sha
NAME=PREFIX+'current_native_correlated_shear_q.json'
RECEIPT=PREFIX+'current_native_correlated_shear_q_check.json'
GATE='current_live_correlated_generic_shear_and_q_enclosures_executed'
KEY='original_unresolved_logU_coordinate_rows'
ORDERS=prior.signed.ORDERS;ZERO=(0,0);ep=prior.ep


def original_log_row_compilers():
    """Change only returned metadata, never a defining source assignment."""
    asts=SourceAST();fn=micro_ops.expose_phase_rows(asts);before=copy.deepcopy(fn)
    returned=next(n for n in fn.body if isinstance(n,ast.Return))
    returned.value.keywords.append(ast.keyword(arg=KEY,value=ast.parse('dict(logU=logU,coordinate_scale=h)',mode='eval').body))
    compare=copy.deepcopy(fn);next(n for n in compare.body if isinstance(n,ast.Return)).value.keywords.pop()
    if ast.dump(compare)!=ast.dump(before):raise ValueError('Original phase math changed while exposing log rows')
    env=dict(vars(micro_original));phase=micro_ops.compile_function(fn,env,'<original phase; unresolved logU output only>')
    env=dict(vars(micro_original));env['phase_physical']=phase
    micro=micro_ops.compile_function(asts.method('microswitch_mixed_C4','evaluate'),env,
        '<original microswitch with output-only logU callback>')
    replay=bridge_original.REPLAY;bridge_env=dict(replay.__globals__);bridge_env['phase_physical']=phase
    exposed=FunctionType(replay.__code__,bridge_env,replay.__name__,replay.__defaults__,replay.__closure__)
    exposed.__kwdefaults__=replay.__kwdefaults__
    env=dict(vars(bridge_original));env['REPLAY']=exposed
    bridge=micro_ops.compile_function(asts.method('current_actual_bridge_mixed_C4','evaluate'),env,
        '<unchanged current bridge; logU output callback only>')
    asts.method('bridge_mixed_C4','evaluate');asts.method('current_actual_bridge_mixed_C4','current_coordinate_replay')
    return bridge,micro,dict(passed=True,original_phase_assignments_unchanged=True,
        only_existing_logU_and_h_return_metadata_added=True,
        original_current_bridge_replay_code_identical=exposed.__code__ is replay.__code__,
        original_micro_and_current_bridge_wrappers_unchanged=True,input_hashes=asts.hashes)


def analytic_shear_identity():
    y,Z=sy.symbols('y Z');L=sy.Function('actual_logU')(y,Z);E=sy.exp(L)
    a=(E-2*sy.diff(E,y))/E;checks={}
    for j,k in ORDERS:
        expected=(1 if (j,k)==ZERO else 0)-2*sy.diff(L,y,j+1,Z,k)
        if sy.simplify(sy.diff(a,y,j,Z,k)-expected)!=0:raise ArithmeticError('Correlated ordinary log-shear identity failed')
        checks['y%d_Z%d'%(j,k)]=True
    return dict(passed=True,exact_ordinary_shear_identities=checks,
        frozen_amplitude_and_width_not_materialized=True)


def square(value):
    """The same source squared; do not multiply independent signed boxes."""
    return prior.ScaledEnclosure(value.scale+value.scale,value.coefficient**2,value.ledger)


def multiply_jet(left,right):
    result={}
    for j,k in ORDERS:
        parts=[]
        for i in range(j+1):
            for ell in range(k+1):
                a,b=left[(i,ell)],right[(j-i,k-ell)]
                part=square(a) if a is b else a*b
                parts.append(part*(math.comb(j,i)*math.comb(k,ell)))
        result[(j,k)]=sum(parts,left[ZERO].scalar(0))
    return result


def quotient_jet(numerator,denominator,lower):
    result={}
    for j,k in sorted(ORDERS,key=lambda o:(sum(o),o)):
        correction=[]
        for i in range(j+1):
            for ell in range(k+1):
                if i==ell==0:continue
                correction.append(denominator[(i,ell)]*result[(j-i,k-ell)]*(math.comb(j,i)*math.comb(k,ell)))
        value=numerator[(j,k)]-sum(correction,numerator[ZERO].scalar(0))
        result[(j,k)]=value.positive_divide(denominator[ZERO],lower)
    return result


def nonnegative_sqrt(value):
    lo,hi=ep(value.coefficient)
    if lo<0:raise ValueError('Nonnegative enclosing source required for square root')
    # Collect the log before halving so a quarter power need not be materialized.
    scale=prior.FormalScale(value.scale.bases,offset=value.scale.evaluate()/2)
    return prior.ScaledEnclosure(scale,value.ctx.sqrt(value.coefficient),value.ledger)


def q_enclosure(a,Delta,eta_log,log_a_lower):
    """Same original lazy cutoff on full active/flat/mixed source boxes."""
    c=a.ctx;scalar=a.scalar
    eta=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=eta_log),1,a.ledger)
    diff=Delta-eta;lo,hi=ep(diff.coefficient)
    if lo>=0:
        return dict(q=scalar(0),branch='flat',sigma_interval=c.mpf(0),
            active_body_evaluated=False,branch_difference=diff,
            cutoff_scope='original Delta>=eta flat branch',positive_q_lower_retained=False)
    branch='active' if hi<0 else 'requires_source_box_refinement'
    gamma=eta*2-Delta
    if ep(gamma.coefficient)[1]<=0:raise ArithmeticError('Mixed branch cover must contain positive active gamma')
    if ep(gamma.coefficient)[0]<=0:
        # On the active part Delta<eta implies gamma>eta. The flat part is
        # evaluated separately as exact zero and included below in the hull.
        gamma=gamma.positive_intersection(eta_log)
    ratio=gamma.positive_divide(a*2,log_a_lower+c.ln(2))
    root=nonnegative_sqrt(ratio)
    if Delta.zero or ep(Delta.coefficient)[1]<=0:
        sigma=c.mpf(1);sigma_scope='Delta<=0 implies original sigma argument>=1 exactly'
    elif branch=='active':
        try:
            d=Delta.positive_divide(eta,eta_log)
            try:normalized=d.finite_interval()
            except ArithmeticError:
                if ep(d.scale.evaluate())[1]<-1000:normalized=d.coefficient*d.bounded_exp(d.scale.evaluate())
                else:raise
            sigma=prior.sigma_jets(c,1-normalized)[0]
            sigma_scope='original monotone flat sigma interval at1-Delta/eta'
        except ArithmeticError:
            sigma=c.mpf([0,1]);sigma_scope='conservative original sigma range; active argument needs refinement'
    else:
        sigma=c.mpf([0,1]);sigma_scope='conservative original sigma range on mixed branch box'
    q=root*sigma
    if branch!='active':
        q=prior.ScaledEnclosure(q.scale,c.mpf([0,ep(q.coefficient)[1]]),q.ledger)
    return dict(q=q,branch=branch,sigma_interval=sigma,active_body_evaluated=True,
        branch_difference=diff,positive_active_gamma_enclosure=gamma,
        cutoff_scope=sigma_scope,positive_q_lower_retained=ep(q.coefficient)[0]>0,
        exact_zero_not_substituted_for_small_positive_q=True)


class NativeCorrelatedShearQ:
    def __init__(self,owner):
        if type(owner) is not prior.NativeSignedInputEnclosures:raise ValueError('Existing native signed-root backend required')
        self.owner=owner;self.native=owner.native;self.service=owner.service;self.ctx=owner.ctx;self.family=owner.family
        admitted=json.loads((HERE/prior.RECEIPT).read_bytes())
        if not admitted['all_passed'] or not admitted[prior.GATE] or admitted['source_family']!=self.family:
            raise ValueError('Same accepted native signed-root stage required')
        self.service.bind_hashes(admitted['input_hashes']);self.service.bind_hashes({prior.RECEIPT:sha(prior.RECEIPT)})
        self.bridge,self.micro,self.compiler=original_log_row_compilers()
        self.service.bind_hashes(self.compiler['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.identity=analytic_shear_identity();self.queries=[]
    def correlated_a_rows(self,packet,coordinate):
        chart=packet.chart;c=self.ctx;algebra=packet.algebra
        adapter=lambda value:packets.decode_row(algebra,packets.factored_rows_record(value))
        rule=None
        if chart in ('bridge_first','bridge_second','bridge_macro','switch_first','switch_second'):
            original=(self.bridge(self.native.bridge,packet.provenance['Z_box'],coordinate,chart.replace('bridge_',''))
                if chart.startswith('bridge_') else self.micro(self.native.switch,packet.provenance['Z_box'],coordinate,chart.replace('switch_','')))
            logrows=original[KEY]['logU']
            source=original[micro_ops.SOURCE_KEY];source_algebra=source['algebra']
            if len(logrows)!=4 or any(a._mpi_!=b._mpi_ for a,b in zip(source_algebra.logs,algebra.logs)):
                raise ValueError('Original four logarithmic rows and same source bases required')
            ordinary=[]
            for j,row in enumerate(logrows[:3]):
                converted=source_algebra.lift(row) if chart=='bridge_macro' else source_algebra.width(row,-j-1)
                ordinary.append(adapter(converted))
            rule='unchanged original phase logU rows; hb^-(j+1) exactly once, macro ordinary rows unchanged'
        elif chart in ('switch_power','reshape'):
            source=self.native.original_raw(chart,packet.provenance['Z_box'],coordinate)['source_packet']
            ordinary=[adapter(row) for row in source['actual_log_Utheta_y_derivative_axial5'][:3]]
            rule='original variable ordinary logU_y axial5 rows, amplitude independent'
        elif chart in ('inner_reference','axial_restore','restore_buffer'):
            ordinary=[algebra.lift(c.mpf('.1')),algebra.lift(0),algebra.lift(0)]
            rule='accepted original raw_restore_rows: theta_y^j=u*(1/10)^j, constant logU_y=1/10'
        elif chart in ('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_slope_mu','O3_power'):
            source=self.native.original_raw(chart,packet.provenance['Z_box'],coordinate)['source_packet']
            ordinary=[adapter(copy_jet(c,row)) for row in source['log_Utheta_ordinary_y_derivatives'][:3]]
            rule='original pre ordinary logU_y source rows; original context copied exactly'
        else:
            return None,dict(correlated_log_source_available=False,source='original signed quotient retained for implicit patch')
        a=[algebra.lift(1)-ordinary[0]*2]+[-row*2 for row in ordinary[1:]]
        return a,dict(correlated_log_source_available=True,source=rule,
            original_logU_ordinary_y_rows=ordinary,
            shear_ordinary_y_rows=a,amplitude_quotient_not_used_to_define_shear=True)
    @native.inlet.source_precision
    def query(self,chart,Z,coordinate):
        incoming=self.owner.query(chart,Z,coordinate);packet=incoming['packet'];old=incoming['roots']
        arows,evidence=self.correlated_a_rows(packet,coordinate)
        roots={key:{order:row['y%d_Z%d'%order] for order in ORDERS} for key,row in old.items()}
        bases=roots['a'][ZERO].scale.bases;ledger=roots['a'][ZERO].ledger
        if arows is not None:
            a={}
            for j,k in ORDERS:
                row=prior.signed.ordinary_axial_coefficient(arows[j],k)
                a[(j,k)]=self.owner.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),bases,ledger)
            roots['a']=a
        positive=self.owner.decode(self.owner.inventory[chart]['actual_positive_denominator_theorem'])
        roots['a'][ZERO]=roots['a'][ZERO].positive_intersection(positive['log_actual_a_positive_lower'])
        roots['t0']=quotient_jet({k:-v for k,v in roots['b'].items()},roots['a'],positive['log_actual_a_positive_lower'])
        b2=multiply_jet(roots['b'],roots['b']);quotient=quotient_jet(b2,roots['a'],positive['log_actual_a_positive_lower'])
        if not chart.startswith('O3_'):
            roots['kappa_minus2']={k:roots['a'][k]+quotient[k]-(2 if k==ZERO else 0) for k in ORDERS}
        else:
            # Exact Delta is the original source, even when adding mu to2 rounds.
            roots['a']={k:roots['kappa_minus2'][k]+(2 if k==ZERO else 0) for k in ORDERS}
        eta=packets.interval(self.ctx,self.owner.scales['selected_positive_eta_log'])
        loop=q_enclosure(roots['a'][ZERO],roots['kappa_minus2'][ZERO],eta,positive['log_actual_a_positive_lower'])
        record=dict(chart=chart,source_family=self.family,source_provenance=packet.provenance,
            original_correlated_shear_source=evidence,
            previous_independent_root_branch=incoming['record']['branch_against_actual_eta'],
            correlated_shear_and_signed_root_enclosures={key:{'y%d_Z%d'%k:v.record() for k,v in row.items()} for key,row in roots.items()},
            original_q_enclosure={key:(v.record() if isinstance(v,prior.ScaledEnclosure) else v) for key,v in loop.items()},
            numerical_arithmetic_ledger=dict(ledger),
            exact_O3_excess_derivatives_retained=True,correlated_b_squared_nonnegative=True,
            radius_phase_correlations_evaluated=False,source_caps_or_midpoints_used_as_field_values=False,
            actual_phase_inverse_installed=False,changed_density_integrals_installed=False,
            **dict.fromkeys(packets.OPEN,False))
        self.queries.append(dict(chart=chart,branch=loop['branch'],old_branch=incoming['record']['branch_against_actual_eta'],
            correlated_shear=evidence['correlated_log_source_available'],positive_q_lower=loop['positive_q_lower_retained']))
        return dict(record=record,packet=packet,roots=roots,q=loop['q'],loop=loop,evidence=evidence)


def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge)))
        old=json.loads((HERE/prior.NAME).read_bytes())['current_native_signed_source_root_records']
        records={}
        for chart,previous in old.items():
            p=previous['source_provenance'];Z=packets.interval(owner.ctx,p['Z_box']);v=packets.interval(owner.ctx,p['coordinate_box'])
            records[chart]=owner.query(chart,Z,v)['record']
            print('Live correlated shear/q:',chart,records[chart]['original_q_enclosure']['branch'],flush=True)
    result=dict(source_family=owner.family,**{GATE:True},live_native_chart_count=len(records),
        actual_correlated_shear_and_q_records=records,query_inventory=owner.queries,
        unchanged_original_log_row_output_compilers=owner.compiler,analytic_shear_identity=owner.identity,
        actual_live_q_enclosure_functions_available=True,
        actual_phase_inverse_installed=False,actual_changed_defect_integral_functions_installed=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Live correlated original shear jets and Section11 q interval functions on17 original query boxes; unresolved branch boxes retain conservative q covers. No radius/phase correlation, inverse, changed density integrals, common N, repair or recursive/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
