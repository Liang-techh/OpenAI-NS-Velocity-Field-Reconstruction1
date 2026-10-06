"""Original cone attached to signed sectors of the checked current graph.

Current pulse main/exit/end have source-correlated shear/log reductions.
Other charts return complete signed records with an explicit missing shear
adapter. No source envelope is promoted to a point or a global lift.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_global_tensor_cover import (
    CurrentGlobalTensorCover,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN as EARLIER_OPEN)
from lei_ren_part1_paper_compliant_current_original_cone_operator import (
    original_cone_theorem,cone_margins,reference_covariance)
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import nonpositive_exp
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST

NAME=PREFIX+'current_original_cone.json'
RECEIPT=PREFIX+'current_original_cone_check.json'
VIEWS_NAME=PREFIX+'current_original_cone_views.json.gz'
GATES=('current_original_paper_cone_variable_and_scale_map_available',
       'current_signed_pulse_cone_margin_adapter_available',
       'current_whole_pulse_end_signed_two_vector_cone_certified',
       'current_original_covariance_and_signed_linear_algebra_available')
OPEN=EARLIER_OPEN
REGIONS=('pulse_main','pulse_exit','pulse_end')
TGROUP='physical_cylindrical_stress_mixed3/'


def current_source_bindings():
    asts=SourceAST();bound={}
    for stem,method in (
        ('current_pulse_end_background_tensor','end'),
        ('current_pulse_main_exit_background_tensor','chart'),
        ('pulse_end_physical_C2','pulse_velocity_rows'),
        ('pulse_end_physical_C2','lift_physical_packet'),
        ('pulse_main_exit_physical_C2','main_exit_to_physical_packet'),
        ('pulse_main_exit_physical_C2','compiled_main_exit_lift'),
        ('current_tensor_registry','native')):
        fn=asts.method(stem,method)
        import ast
        bound[stem+'.'+method]=hashlib.sha256(ast.dump(fn,include_attributes=False).encode()).hexdigest()
    asts.expression('current_pulse_end_background_tensor','end','velocity',wanted='pulse_velocity_rows(c,delta,mu,z,C,Bh,m)')
    asts.expression('current_pulse_end_background_tensor','end','physical',wanted='lift_physical_packet(c,packet,delta,velocity,lt,theta,nu)')
    asts.expression('current_pulse_main_exit_background_tensor','chart','velocity',wanted="dict(local=rows['velocity_local'],incoming=rows['velocity_incoming'])")
    asts.expression('current_pulse_main_exit_background_tensor','chart','point',wanted='self.lift(c,mapped,delta,velocity,lt,theta,nu)')
    asts.expression('pulse_end_physical_C2','lift_physical_packet','beta',wanted='-2-delta')
    import sympy as s
    zz,delta,mu=s.symbols('Z delta mu',real=True);CC=1/(1+zz**2);rate=1-mu
    actual=asts.evaluate(asts.expression('pulse_end_stress_C3','pulse_coefficients','equilibrium',
        wanted='(C*((mu-delta/2)/r)-z*axial_derivative(C)*(b/r))/L'),
        dict(C=CC,mu=mu,delta=delta,r=rate,z=zz,b=(1-delta)/2,L=1-delta*zz**2,axial_derivative=lambda x:s.diff(x,zz)))
    positive=(CC*(mu-delta/2)/rate+(1-delta)*zz**2*CC**2/rate)/(1-delta*zz**2)
    if s.cancel(actual-positive)!=0:raise ValueError('Original equilibrium positive rewrite differs')
    # This is the actual original y-derivative source, not a generic shear fitted to T.
    import ast
    fn=asts.method('pulse_end_physical_C2','pulse_velocity_rows')
    expected=ast.dump(ast.parse("[C*(-(c.mpf('.5')+mu))**j for j in range(5)]",mode='eval').body)
    theta=[kw.value for node in ast.walk(fn) if isinstance(node,ast.Call) and ast.unparse(node.func)=='dict'
           for kw in node.keywords if kw.arg=='theta']
    if len(theta)!=1 or ast.dump(theta[0])!=expected:raise ValueError('Original pulse pure-power theta jets changed')
    expected_axial=ast.dump(ast.parse("shifted_rows([C*v for v in Bh],-(c.mpf('.5')+mu),4)",mode='eval').body)
    axial=[kw.value for node in ast.walk(fn) if isinstance(node,ast.Call) and ast.unparse(node.func)=='dict'
           for kw in node.keywords if kw.arg=='axial']
    if len(axial)!=1 or ast.dump(axial[0])!=expected_axial:raise ValueError('Original normalized axial y jets changed')
    native_end=asts.method('axial_pulse_field','end');native_main=asts.method('axial_pulse_field','main')
    for node,key,wanted in ((native_end,'Uz_y_over_Utheta',"(Byhat-Bhat*(c.mpf('.5')+self.mu))*self.Ecap"),
                            (native_main,'Uz_y_over_Utheta',"By-B*(c.mpf('.5')+self.mu)")):
        rows=[kw.value for call in ast.walk(node) if isinstance(call,ast.Call) for kw in call.keywords if kw.arg==key]
        if len(rows)!=1 or ast.dump(rows[0])!=ast.dump(ast.parse(wanted,mode='eval').body):
            raise ValueError('Native main/end axial shear factor source changed')
    asts.expression('axial_pulse_field','end','B',wanted='Bhat*self.Ecap')
    asts.expression('current_pulse_end_background_tensor','end','logD',wanted='c.mpf(self.pulse.logE)')
    asts.bindings['normalized_axial_rows_and_first_y_shear_AST']=ast.unparse(axial[0])
    asts.bindings['one_exact_end_D_factor_from_same_native_logE']=True
    return dict(defining_method_AST_sha256=bound,exact_source_assignment_bindings=asts.bindings,
        theta_source='Utheta=B*C; Utheta_y=-(1/2+mu)*B*C',
        axial_source='End: Uz=B*D*C*Bhat. Main/exit: Uz=B*C*Bhat; incoming velocity only contributes to radial recovery.',
        positive_F='B>0 as exp(finite exact logs); C=1/(1+Z^2)>0; R>0. F=B*C/sqrt(2R)>0.',
        a_minus2_source='2*mu, retained before adding it to 2; numerical subtraction a-2 is not used',
        exact_equilibrium_positive_rewrite='(C*(mu-delta/2)/(1-mu)+(1-delta)*Z^2*C^2/(1-mu))/(1-delta*Z^2); C=1/(1+Z^2)',
        same_velocity_and_stress_source_packet_bound=True,input_hashes=asts.hashes,passed=True)


def coefficient_scaled_exp(c,value,coefficient,target):
    """Enclose exp(value) relative to a positive requested error scale.

    A fixed exp(-1024) cap can dominate an astronomical small mu after
    multiplication by a large coefficient. This cap is chosen from the
    coefficient and target; the original exponent remains the source.
    """
    magnitude=max(abs(x) for x in endpoints(coefficient))
    if not magnitude:return nonpositive_exp(c,value),None
    if endpoints(value)[1]>0 or endpoints(target)[0]<=0:raise ValueError('Nonpositive source log and positive target required')
    cap=c.ln(target)-c.ln(c.mpf(magnitude))-1000
    cutoff=endpoints(cap)[0];lo,hi=endpoints(value)
    if hi<cutoff:return c.mpf([0,endpoints(c.exp(cap))[1]]),cap
    upper=endpoints(c.exp(c.mpf(hi)))[1]
    lower=0 if lo<cutoff else endpoints(c.exp(c.mpf(lo)))[0]
    return c.mpf([lower,upper]),cap


def normalized_signed_pulse(c,view,mu,delta):
    """Cancel the common positive sqrt(R/2)*B source before enclosure."""
    raw=view['original_complete_view'];packet=raw['current_actual_source_stress_packet']
    sectors=packet['full_meridional_stress_log_sectors'];result={};retained={}
    logs=(packet['exact_logR'],sum(packet['exact_pulse_reference_logB_parts'].values(),c.mpf(0)),
          packet['exact_logD'],packet['exact_logH'])
    for label in ('theta','axial'):
        physical=raw['physical_cylindrical_stress_mixed3'][label]
        if set(physical)!=set(sectors[label]):raise ValueError('Current original signed stress sector lost')
        total=c.mpf(0);retained[label]={}
        canonical=view['canonical_signed_component_groups'][TGROUP+label+'/r0_z0']
        if len(canonical)!=len(physical):raise ValueError('Canonical signed sector count differs')
        for name,part in sectors[label].items():
            row=physical[name]['r0_z0'];coefficient=part['full_stress_mixed3_coefficient_enclosures']['s0_Z0']
            if encode(pack(row['signed_coefficient']))!=encode(pack(coefficient)):
                raise ValueError('Physical order-zero signed coefficient differs from actual source')
            if row['actual_source_log_parts']!=part['exact_source_log_parts']:
                raise ValueError('Physical source factors differ from the original pulse packet')
            original_coefficient=coefficient
            if label=='theta' and name=='equilibrium':
                Z=c.mpf(raw['Z']);square=Z**2;C=1/(1+square)
                positive=(C*(mu-delta/2)/(1-mu)+(1-delta)*square*C**2/(1-mu))/(1-delta*square)
                lo=max(endpoints(coefficient)[0],endpoints(positive)[0]);hi=min(endpoints(coefficient)[1],endpoints(positive)[1])
                if lo>hi:raise ValueError('Source-identical equilibrium enclosures are inconsistent')
                coefficient=c.mpf([lo,hi])
            # Subtract powers exactly, never enormous interval log boxes.
            relative=tuple(c.mpf(v)-base for v,base in zip(part['mode'],(c.mpf('.5'),1,0,0)))
            relative_log=sum((power*log for power,log in zip(relative,logs) if endpoints(power)!=(0,0)),c.mpf(0))
            if 'extra_source' in part:relative_log+=packet['exact_source_logs']['extra'][part['extra_source']]
            if endpoints(coefficient)==(0,0):factor=None;signed=c.mpf(0)
            else:
                if endpoints(relative_log)[1]>0:raise ValueError('Positive relative factor requires an additional correlated scale reduction: '+name)
                factor,cap=coefficient_scaled_exp(c,relative_log,coefficient,mu);signed=coefficient*factor
            total+=signed
            retained[label][name]=dict(original_signed_physical_row=row,
                source_signed_coefficient=coefficient,exact_relative_log_mode=relative,
                original_unreduced_signed_coefficient=original_coefficient,
                source_positive_equilibrium_rewrite_used=(label=='theta' and name=='equilibrium'),
                exact_relative_log_enclosure=relative_log,positive_factor_enclosure=factor,
                coefficient_relative_cap_log=None if endpoints(coefficient)==(0,0) else cap,
                normalized_signed_contribution=signed)
        result[label]=total
    return dict(common_positive_source='nu*lambda^(-2-delta)*sqrt(R/2)*B',
        common_scale_canceled_algebraically_before_interval_enclosure=True,
        signed_totals=result,all_signed_sectors=retained,
        cap_is_an_enclosure_of_exp_not_a_replacement_source=True,
        interval_midpoints_or_absolute_envelopes_not_used_as_stress=True)


class CurrentOriginalCone:
    @source_precision
    def __init__(self,cover=None,require_checked=True):
        self.cover=cover if cover is not None else CurrentGlobalTensorCover()
        if type(self.cover) is not CurrentGlobalTensorCover or not self.cover.acceptance_loaded:
            raise ValueError('Checked current global T/E cover required')
        self.registry=self.cover.registry;self.ctx=self.cover.ctx
        self.family=self.cover.family;self.source=self.cover.source;self.datum_sha=self.cover.datum_sha
        self.theorem=original_cone_theorem();self.bindings=current_source_bindings()
        self.hashes={**self.cover.hashes,**self.theorem['input_hashes'],**self.bindings['input_hashes']}
        for stem in ('current_original_cone','current_original_cone_operator'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Original cone operator receipt exceeds its regional/algebraic scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.cover.assert_graph()
        if self.registry is not self.cover.registry or current_source_bindings()!=self.bindings:
            raise ValueError('Current original cone graph/source changed')

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();view=self.registry.native(region,Z,coordinate,log_tau,theta,viscosity)
        if (view['actual_five_defect_family_sha256'],view['implicit_source_sha256'],view['datum_enclosure_sha256'])!=(self.family,self.source,self.datum_sha):
            raise ValueError('Foreign signed tensor source')
        if region not in REGIONS:
            return dict(region=region,original_complete_signed_tensor_view=view,
                status='same_source_shear_adapter_pending',scope='Full source records retained; no inferred cone status',
                module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
                **dict.fromkeys(GATES+OPEN,False))
        c=self.ctx;raw=view['original_complete_view'];packet=raw['current_actual_source_stress_packet']
        owner=self.registry.owners['end' if region=='pulse_end' else 'main'];mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta)
        velocity=raw['current_source_three_component_velocity_rows']
        if region!='pulse_end':velocity=velocity['local']
        C=velocity['theta'][0][0]
        if endpoints(C)[0]<=0:raise ValueError('Original positive theta coefficient required')
        b_coefficient=2*velocity['axial'][1][0]/C
        if region=='pulse_end':Dfactor,Dcap=coefficient_scaled_exp(c,c.mpf(packet['exact_logD']),b_coefficient,mu)
        else:Dfactor,Dcap=c.mpf(1),None
        bs=Dfactor*b_coefficient;a=2+2*mu
        km=2*mu+bs**2/a
        reduced=normalized_signed_pulse(c,view,mu,delta);tt=reduced['signed_totals']['theta'];tz=reduced['signed_totals']['axial']
        tested=cone_margins(c,a,bs,tt,tz,km)
        reference=reference_covariance(c,a,bs,km,tt,tz)
        return dict(region=region,source_family=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            original_complete_signed_tensor_view=view,normalized_current_signed_stress=reduced,
            current_source_shear=dict(a=a,bs=bs,vs_minus2=km,positive_F_by_original_recipe=True,
                exact_a_minus2_source=2*mu,nonnegative_axial_factor_enclosure=Dfactor,
                exact_axial_factor_positive_by_finite_source_log=True,
                exact_axial_factor_log=packet['exact_logD'],coefficient_relative_axial_cap_log=Dcap,
                full_axial_shear_retained=True,source_velocity_jets=velocity),
            original_leading_two_vector_cone=tested,reference_covariance_diagnostic=reference,
            original_positive_lambda_and_viscosity_scale_cancels=True,
            current_completed_diagonal_divergence_and_all_remainder_groups_retained=True,
            higher_order_signed_corrections_require_Proposition_7_6=True,
            scope='Directed current source cone margins on this requested box; not whole annulus, smooth edge amplitudes, or actual wave covariance realization',
            module_level_construction_gates={k:self.acceptance_loaded for k in GATES},
            **{k:(self.acceptance_loaded and (k!='current_whole_pulse_end_signed_two_vector_cone_certified' or region=='pulse_end')) for k in GATES},
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            original_paper_variable_cone_and_covariance_theorem=self.theorem,current_actual_source_bindings=self.bindings,
            current_source_correlated_shear_and_signed_margin_regions=list(REGIONS),
            full_signed_tensor_record_routes=list(self.registry.registry),
            actual_current_tensor_region_count=33,actual_current_completed_tensor_adjacent_interface_count=32,
            actual_current_completed_tensor_internal_interface_count=14,
            scope='Original source/equation map and current signed pulse cone adapter; reference and general matrix algebra only. All-region cone repair, actual homogeneous pulses/H, uniform edge weights, true coefficient recovery, corrected NS, energy and flatness remain open.',
            complete_fresh_signed_view_storage=VIEWS_NAME,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentOriginalCone(require_checked=False)
    result=field.manifest();(HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original cone source map and signed current pulse margin adapter built; actual covariance waves remain open',flush=True)
    return result
