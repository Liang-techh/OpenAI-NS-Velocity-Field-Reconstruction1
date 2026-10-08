"""Source-selected scales and true-radius O2 conditioned C0 primitives.

Genuine defining point coefficients and their errors replace broad chart
covers. The original selected auxiliary scales remain logarithmically
factored. The accepted two-angle kernel encloses the original inverse and
A/B at the actual radius phase, without an ancestor construction.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_radius_phase_points as radius
import lei_ren_part1_paper_compliant_current_native_conditioned_phase as conditioned

point=radius.point;prior=conditioned.prior;current=conditioned.current
HERE,PREFIX,sha=radius.HERE,radius.PREFIX,radius.sha
NAME=PREFIX+'current_original_O2_conditioned_primitives.json'
RECEIPT=PREFIX+'current_original_O2_conditioned_primitives_check.json'
GATE='original_O2_source_point_and_true_radius_conditioned_C0_inverse_A_B_enclosures_connected'
UNIFORM=PREFIX+'current_generic_shear_uniform_inputs.json'
UNIFORM_CHECK=PREFIX+'current_generic_shear_uniform_inputs_check.json'
FUNCTION=PREFIX+'current_generic_loop_function_sources.json'
FUNCTION_CHECK=PREFIX+'current_generic_loop_function_sources_check.json'


def encoded(value):
    """Serialize cloned scalar and interval contexts without losing exact bits."""
    if hasattr(value,'_mpi_'):
        lo,hi=point.endpoints(value)
        return dict(lower=mp.nstr(lo,80),upper=mp.nstr(hi,80),
            lower_exact_mpf_tuple=list(lo._mpf_),upper_exact_mpf_tuple=list(hi._mpf_))
    if hasattr(value,'_mpf_'):
        return dict(decimal=mp.nstr(value,80),exact_mpf_tuple=list(value._mpf_))
    if isinstance(value,dict):return {key:encoded(row) for key,row in value.items()}
    if isinstance(value,(list,tuple)):return [encoded(row) for row in value]
    return value


def bits(record):
    if record['lower_exact_mpf_tuple']!=record['upper_exact_mpf_tuple']:
        raise ValueError('Explicit original selected auxiliary log constant required')
    return tuple(record['lower_exact_mpf_tuple'])


def factored_row_enclosure(row,bases,ledger):
    """Restore a genuine source row/error in one shared positive-factor basis."""
    c=bases[0].ctx;scalar=lambda value:prior.ScaledEnclosure(prior.FormalScale(bases),value,ledger)
    result=scalar(0)
    for term in row.terms:
        r,p,d,ell=term.factor_powers
        center=c.mpf(term.coefficient);error=c.mpf(term.coefficient_error_upper)
        coefficient=center+c.mpf([-point.endpoints(error)[1],point.endpoints(error)[1]])
        if term.late_pressure_error_log_upper is not None:
            bound=scalar(1).bounded_exp(c.mpf(term.late_pressure_error_log_upper))
            upper=point.endpoints(bound)[1];coefficient+=c.mpf([-upper,upper])
        result+=prior.ScaledEnclosure(prior.FormalScale(bases,(p,d,ell,0,r)),coefficient,ledger)
    return result


class FactoredGenericLoopScaleFrame:
    """Accepted whole-input auxiliary choices, never source field samples."""
    def __init__(self,family,ctx):
        self.family=family;self.ctx=ctx;self.hashes={}
        for name,gate in ((UNIFORM_CHECK,'current_original_whole_generic_input_margins_and_log_scales_certified'),
            (FUNCTION_CHECK,'current_original_generic_loop_velocity_and_five_increment_function_graphs_defined'),
            (conditioned.RECEIPT,conditioned.GATE)):
            receipt=json.loads((HERE/name).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate) or receipt['source_family']!=family:
                raise ValueError('Same accepted original scale/function/inverse family required: '+name)
            for filename,digest in {**receipt['input_hashes'],name:sha(name)}.items():
                if filename in self.hashes and self.hashes[filename]!=digest:
                    raise ValueError('Selected-scale dependencies disagree: '+filename)
                if sha(filename)!=digest:raise ValueError('Selected-scale source changed: '+filename)
                self.hashes[filename]=digest
        uniform=json.loads((HERE/UNIFORM).read_bytes())
        self.scales=uniform['current_actual_logarithmic_loop_scales']
        if not uniform['current_whole_source_loop_scales_in_log_form_certified']:
            raise ValueError('Whole-source selected logarithmic scale contract required')
        if self.scales['constant_interpretation']!='Each selected positive constant is exp(its saved log); None upper is exact zero. No source exponentiation is performed.':
            raise ValueError('Original auxiliary constant interpretation changed')
        lower=self.scales['logarithmic_selected_positive_lower_constants']
        self.logs={name:ctx.mpf(mp.mp.make_mpf(bits(record))) for name,record in lower.items()}
        self.logs['eta']=ctx.mpf(mp.mp.make_mpf(bits(self.scales['selected_positive_eta_log'])))
        self.upper_logs={name:None if record is None else ctx.mpf(mp.mp.make_mpf(bits(record)))
            for name,record in self.scales['logarithmic_conservative_upper_constants'].items()}
        graph=json.loads((HERE/FUNCTION).read_bytes());view=graph['actual_generic_loop_function_graph_views']
        self.graph=json.loads(gzip.decompress((HERE/view).read_bytes()))['O2_slope']
        if self.graph['source_family']!=family:raise ValueError('Original O2 function graph family differs')
        nodes={node['name']:node for node in self.graph['function_graph_nodes']
            if node['operation']=='original_positive_log_scale'}
        if set(nodes)!={'eta','d_star'}:raise ValueError('Original selected loop scale graph roots required')
        for name,record in (('eta',self.scales['selected_positive_eta_log']),('d_star',lower['d_star'])):
            if bits(nodes[name]['selected_source_log'])!=bits(record) or not nodes[name]['source_scale_never_materialized']:
                raise ValueError('Exact function graph and auxiliary choice differ: '+name)
        binding=self.graph['common_phase_binding']
        if binding['phi']!='fractional_part(N*log(R/r_minus))' or not binding['total_Z_phase_derivative_exact_zero']:
            raise ValueError('Original shared radius phase convention required')
        self.graph_hash=sha(view)
        for name in (UNIFORM,FUNCTION,view):self.hashes[name]=sha(name)

    def record(self):
        return dict(source_family=self.family,selected_positive_logs=self.logs,
            selected_auxiliary_upper_logs=self.upper_logs,
            accepted_scale_recipe=self.scales['original_Section11_scale_recipe'],
            O2_function_graph_sha256=self.graph_hash,
            whole_original_input_margins_and_scale_choices_attached=True,
            source_field_caps_not_selected_as_values=True,
            selected_eta_dstar_not_materialized=True,global_N_not_selected=True)


class OriginalO2ConditionedPrimitives:
    mode='original_O2_true_radius_conditioned_C0_primitive_enclosures'
    def __init__(self,dps=260,cells=4096,coefficient_dps=50,phase_dps=80):
        if type(dps) is not int or dps<260:raise ValueError('At least 260 interval digits preserve selected log bits')
        if type(coefficient_dps) is not int or coefficient_dps<40:raise ValueError('At least 40 defining-coefficient digits required')
        self.inputs=point.OriginalO2FactoredPointInputs(dps=coefficient_dps,cells=cells)
        self.radius=radius.OriginalO2RadiusPhasePoints(dps=phase_dps);self.family=self.inputs.family
        definitions=self.inputs.frame.definitions;logP=definitions['logPstar']
        if definitions['Md']!=40 or s.simplify(logP-s.exp(40)-11)!=0:
            raise ValueError('Original O2 logarithmic pressure parameter recipe differs')
        if s.simplify(definitions['log_delta']+4*logP+30)!=0:
            raise ValueError('Original O2 delta factor recipe differs')
        if s.expand(self.inputs.frame.radius_log(0)-s.log(110)-10*(definitions['logCstar']+logP))!=0:
            raise ValueError('Original O2 radius factor recipe differs')
        if self.inputs.frame.selected_logCstar_mpf_tuple!=self.radius.frame.selected_logCstar_mpf_tuple:
            raise ValueError('Original source point and radius phase selected parameters differ')
        self.precisions=dict(interval_digits=dps,defining_point_coefficient_digits=coefficient_dps,
            directed_radial_coefficient_digits=coefficient_dps+40,radius_phase_digits=phase_dps,
            directed_radial_cells=cells,separate_precisions_with_all_errors_propagated=True)
        self.ctx=MPIntervalContext();self.ctx.dps=dps
        self.scales=FactoredGenericLoopScaleFrame(self.family,self.ctx)
        self.hashes={}
        for closure in (self.inputs.hashes,self.radius.hashes,self.scales.hashes):
            for name,digest in closure.items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original loop dependencies disagree: '+name)
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={}

    def query(self,*,y,Z):
        incoming=self.inputs.evaluate(y=y,Z=Z);c=self.ctx
        y=point.source.exact_coordinate(y);Z=point.pressure.exact_Z(Z)
        ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
            positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        with mp.workdps(c.dps+40):
            iz=c.mpf(int(Z.p))/int(Z.q);iy=c.mpf(int(y.p))/int(y.q)
            logP=c.exp(40)+11;logdelta=-4*logP-30
            delta=c.exp(logdelta);L=1-delta*iz*iz
            if point.endpoints(L)[0]<=0:raise ArithmeticError('Original positive L source lost')
            logC=c.mpf(mp.mp.make_mpf(self.inputs.frame.selected_logCstar_mpf_tuple))
            logR=c.ln(110)+10*(logC+logP)+iy
            bases=(logP,logdelta,c.ln(L),c.mpf(0),logR)
            scalar=lambda value:prior.ScaledEnclosure(prior.FormalScale(bases),value,ledger)
            roots={key:{(0,order):factored_row_enclosure(row,bases,ledger) for order,row in enumerate(pair)}
                for key,pair in incoming['inputs'].items()}
            # The original endpoint identity is exactly a(1)=2. A forward
            # rational-interval rounding box must not invent a sign of a-2.
            if y==1:roots['a'][(0,0)]=scalar(2)
            if not roots['b'][(0,0)].zero or not roots['b'][(0,1)].zero:
                raise ValueError('Exact original O2 b=0 required for this source adapter')
            roots['t0']={(0,0):scalar(0),(0,1):scalar(0)}
            Delta=roots['a'][(0,0)]-2
            qsource=current.q_enclosure(roots['a'][(0,0)],Delta,self.scales.logs['eta'],self.scales.logs['a_min'])
            kernel=conditioned.ConditionedPhase(dict(q=qsource['q'],roots=roots),self.scales.logs['d_star'])
        return dict(point=incoming,roots=roots,q=qsource,kernel=kernel,ledger=ledger,
            basis_contract=dict(order=['logPstar','logdelta','logL','zero','logR'],
                same_original_R_Pstar_delta_L_bound=True,
                source_L='1-original_delta*Z²',source_R='Rref*exp(y)',
                true_point_radius_log_enclosure_not_selected_as_field_value=True,
                original_positive_eta_delta_and_Pstar_inverse_sectors_retained=True,
                exact_a1_equals2_intersection_used=y==1))

    def evaluate(self,*,y,Z,N,bits=80):
        y=point.source.exact_coordinate(y);Z=point.pressure.exact_Z(Z)
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Explicit inverse bisection bits in[4,256] required')
        phase=self.radius.evaluate(y=y,N=N);key=(y,Z,N,bits)
        if key in self.cache:return self.cache[key]
        query=self.query(y=y,Z=Z);kernel=query['kernel'];c=self.ctx
        rows=[]
        with mp.workdps(c.dps+40):
            for box in phase['true_original_phase_directed_boxes']:
                result=kernel.evaluate(c.mpf([box['lower'],box['upper']]),bits=bits)
                if result['status']!='enclosed':raise ArithmeticError('Refine actual O2 point inputs for signed inverse')
                result.update(free_phase_parameter_not_spatial_phase=False,
                    original_common_N_and_radius_phase_bound=True,
                    actual_O2_defining_point_inputs_and_errors_consumed=True,
                    original_positive_auxiliary_scales_attached=True,
                    source_caps_or_midpoints_used_as_field_values=False)
                rows.append(result)
        result=dict(source_family=self.family,mode=self.mode,original_y_exact=str(y),original_Z_exact=str(Z),
            explicit_candidate_N=N,actual_original_radius_phase=phase,
            actual_original_O2_factored_inputs=point.record_query(query['point']),
            selected_scale_contract=self.scales.record(),source_factor_basis=query['basis_contract'],
            effective_numerical_precisions=self.precisions,
            conditioned_source_geometry=kernel.geometry_record(),actual_radius_phase_C0_inverse_and_A_B_enclosures=rows,
            numerical_arithmetic_ledger=query['ledger'],
            original_C0_inverse_and_primitives_enclose_true_O2_point_source=True,
            scalar_A_B_point_values_selected=False,slow_Z_primitives_installed=False,
            actual_changed_five_moment_integral_evaluated=False,
            numerical_original_source_point_or_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False)
        self.cache[key]=result;return result


def run():
    began=time.monotonic();owner=OriginalO2ConditionedPrimitives()
    records=[]
    for y,Z,N in (('.53','-.37',7),('.53','0',7),('.53','.37',7),('1','.37',257)):
        records.append(owner.evaluate(y=y,Z=Z,N=N))
        print('Live true-radius O2 conditioned primitives:',y,Z,N,records[-1]['conditioned_source_geometry']['branch'],flush=True)
    report=dict(**{GATE:True},source_family=owner.family,mode=owner.mode,
        selected_original_whole_input_scale_contract=owner.scales.record(),
        actual_original_O2_conditioned_point_queries=records,
        original_C0_inverse_and_primitives_enclose_true_O2_point_source=True,
        original_auxiliary_scale_and_function_graph_binding_installed=True,
        native_source_parameter_recipes_checked_against_common_frame=True,
        effective_numerical_precisions=owner.precisions,
        slow_Z_primitives_installed=False,actual_changed_five_moment_integral_evaluated=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='True original O2 source-point/error and actual-radius phase C0 inverse/A-B enclosures using the same accepted original selected logarithmic auxiliary scales. No scalar field point selected from an interval, slow-Z primitives, signed integral, all-chart oracle, common N, controls or recursion is installed.')
    (HERE/NAME).write_text(json.dumps(encoded(report),indent=2)+'\n',encoding='utf8')
    print('True-radius source-selected-scale O2 C0 inverse/A-B enclosure adapter connected',flush=True)
    return report


if __name__=='__main__':run()
