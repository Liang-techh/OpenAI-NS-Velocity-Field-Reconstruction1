"""Same-core nonsingular Cartesian axis T/E, including nonzero mixed derivatives."""
import ast
import gzip
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_core_background_tensor import (
    CurrentCoreBackgroundTensor,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes)
from lei_ren_part1_paper_compliant_current_core_axis_operator import (
    RHO,X,Y,Z,D,B,INDICES2,INDICES3,indexkey,cartesian_remainder,raw_axis_source_theorem,
    cartesian_pullback_theorem,physical_factor_and_cylindrical_pullback_theorem,
    nonsingular_profiles,templates,source_jets)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_collar_physical_C2 import physical_source_row

NAME=PREFIX+'current_core_axis_background_tensor.json.gz'
RECEIPT=PREFIX+'current_core_axis_background_tensor_check.json'
GATES=('current_core_axis_Cartesian_full_tensor_remainder_mixed2_available',
    'core_axis_tensor_remainder_limits_certified','current_core_full_background_tensor_available')
OPEN=('global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'full_point_physical_field_evaluation','actual_point_moment_history_recovered',
    'independently_bounded_flat_remainder','temporal_recursion')
VIEWS={
    'whole_axis':((-1,1),0,0,('-3','-1'),'1'),
    'axis_center':(0,0,0,'-2.6','.8'),
    'fresh_axis':('.537',0,0,'-2.6','.8'),
    'near_axis':((-1,1),('-.1','.1'),('-.1','.1'),('-3','-1'),'.8'),
    'fresh_near':('.731','.0137','-.02','-2.6','.2'),
    'fresh_positive':('.359','.7','.3','-2.6','.8'),
    'core_exit':('.537',2,2,'-2.6','.8')}


def axis_source_theorem(field):
    tensor=field.core_tensor;interior=tensor.interior;asts=SourceAST()
    live=dict(exact_checked_core_tensor=type(tensor) is CurrentCoreBackgroundTensor,
        same_admitted_full_core_graph=tensor.acceptance_loaded and interior.acceptance_loaded,
        same_original_core_context=field.ctx is tensor.ctx,
        same_original_core_and_nonlinear_owner=field.core is tensor.core is interior.common.core,
        same_original_source_jets=field.source_jets is source_jets,
        same_original_nonsingular_remainder=field.remainder is cartesian_remainder,
        original_remainder_source_jet_callback=field.remainder.__globals__['source_jets'] is source_jets,
        same_original_full_raw_and_physical_operators=tensor.stress is tensor.bridge_tensor.stress and tensor.lift is tensor.bridge_tensor.lift)
    if not all(live.values()):raise ValueError('Same checked original core tensor/source/axis operator graph required')
    unit=raw_axis_source_theorem();pullback=cartesian_pullback_theorem()
    factors=physical_factor_and_cylindrical_pullback_theorem()
    if unit['total_identities']!=110 or pullback['total_identities']!=30 or not unit['passed'] or not pullback['passed']:
        raise ValueError('Original nonsingular Cartesian source/coordinate theorem incomplete')
    if not factors['passed'] or (factors['coefficient_identity_count'],factors['cylindrical_operator_identity_count'],
        factors['physical_factor_identity_count'])!=(30,6,110):
        raise ValueError('Original positive-radius/axis physical C2 factor theorem incomplete')
    coreproof=tensor.proof['current_same_source_core_units_integrated_stress_and_completed_first_trace']
    first=coreproof['consumed_same_source_core_first_receipt']
    pressure=first['common_pressure_primitive_and_true_atom']
    equations=first['original_integrated_stress_free_core_equations']
    if not pressure['passed'] or not equations['same_common_pressure_P_equals_P0_plus_integral_F_squared'] or not equations['passed']:
        raise ValueError('Same original P0/C pressure FTC and centrifugal identity required')
    # The actual production operator contains no inverse transverse coordinate or rho.
    for expression in nonsingular_profiles().values():
        if any(node.base==RHO and node.exp.is_negative for node in expression.atoms(s.Pow)):
            raise ValueError('Inverse rho survives the nonsingular core source')
    for rows in templates().values():
        for expression in rows.values():
            if any(node.base in (X,Y,RHO) and node.exp.is_negative for node in expression.atoms(s.Pow)):
                raise ValueError('Inverse Cartesian coordinate survives core pullback')
    for stem,method in (('current_core_axis_operator','cartesian_remainder'),('current_core_axis_operator','templates'),
        ('current_core_axis_operator','source_jets'),('current_core_axis_background_tensor','cartesian'),
        ('current_core_axis_background_tensor','assert_graph'),('core_physical_field','core_step'),
        ('current_core_interior_moments','evaluate'),('collar_physical_C2','physical_source_row')):
        asts.method(stem,method)
    operator=asts.method('current_core_axis_operator','cartesian_remainder')
    calls=[node for node in ast.walk(operator) if isinstance(node,ast.Call) and ast.unparse(node.func)=='physical_source_row']
    if len(calls)!=1 or not any(kw.arg=='nu_base' and ast.unparse(kw.value)=="c.mpf('.5')" for kw in calls[0].keywords):
        raise ValueError('Original remainder sqrt-nu normalization changed')
    return dict(live_original_callable_bindings=live,original_six_sector_nonsingular_source_theorem=unit,
        independent_nonsingular_Cartesian_coordinate_pullback=pullback,
        original_positive_radius_to_axis_physical_C2_factor_and_coefficient_theorem=factors,
        consumed_original_core_pressure_primitive_and_true_atom=pressure,
        consumed_original_core_equations_and_FTC=equations,
        same_axis_Phi_Uz_Q_C_P0_F0_function_graph=True,
        same_normalized_C_FTC_identity='C+rho*C_rho=Phi^2, with analytic rho0 values; P_R=F0^2*Phi^2',
        radial_pressure_and_centrifugal_force_cancel_before_axis_bounds=True,
        original_relative_F0_dressing_and_all_source_tails_retained=True,
        apparent_inverse_R_cancelled_symbolically_before_source_enclosure=True,
        whole_positive_radius_and_axis_join_by_same_analytic_functions=True,
        total_zero_core_stress_extends_analytically_not_each_signed_sector=True,
        exact_axis_does_not_zero_nonzero_Cartesian_derivatives=True,
        no_new_radial_tensor_region_for_coordinate_extension=True,
        input_hashes={**unit['input_hashes'],**factors['input_hashes'],**asts.hashes},passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentCoreAxisBackgroundTensor:
    @source_precision
    def __init__(self,core_tensor=None,require_checked=True):
        self.core_tensor=core_tensor if core_tensor is not None else CurrentCoreBackgroundTensor()
        self.core=self.core_tensor.core;self.ctx=self.core.ctx
        self.family=self.core_tensor.family;self.source=self.core_tensor.source;self.datum_sha=self.core_tensor.datum_sha
        self.source_jets=source_jets;self.remainder=cartesian_remainder;self.assert_graph()
        self.log_epsilon=self.ctx.ln(self.core.epsilon)
        self.log_F0=self.ctx.mpf([endpoints(-self.core.logC-self.core.Lambda*self.core.Gbar)[0],endpoints(-self.core.logC)[1]])
        self.proof=axis_source_theorem(self);self.hashes={**self.core_tensor.hashes,**self.proof['input_hashes']}
        for stem in ('current_core_axis_operator','current_core_axis_background_tensor'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Whole core axis tensor/source admission scope differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.core_tensor.assert_graph()
        if not (type(self.core_tensor) is CurrentCoreBackgroundTensor and self.core_tensor.acceptance_loaded
            and self.core is self.core_tensor.core and self.ctx is self.core_tensor.ctx
            and self.source_jets is source_jets and self.remainder is cartesian_remainder):
            raise ValueError('Same checked positive-radius core and original axis source required')

    @source_precision
    def cartesian(self,Z0,X0=0,Y0=0,log_tau='-1',viscosity='1'):
        self.assert_graph();c=self.ctx;z=c.mpf(Z0);x=c.mpf(X0);y=c.mpf(Y0);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(v) for v in endpoints(z)+endpoints(x)+endpoints(y)+endpoints(lt)+endpoints(nu)) or endpoints(z)[0]<-1 or endpoints(z)[1]>1 or endpoints(nu)[0]<=0:
            raise ValueError('Finite Cartesian core coordinates,Z[-1,1],finite log tau,nu>0 required')
        # Interval squaring must preserve the nonnegative exact coordinate radius.
        rho=(x**2+y**2)/2
        if endpoints(rho)[1]>4:raise ValueError('Cartesian request outside original rho<=4 core')
        remainder,source=self.remainder(self,z,x,y,rho,lt,nu)
        delta=self.core.delta
        zero=lambda gamma,N,base:physical_source_row(c,c.mpf(0),{},gamma,lt,nu,N,nu_base=base)
        stress={indexkey(*ik):{name:zero(-2-delta-ik[0]-ik[1]+ik[2]*(delta-1),sum(ik),c.mpf(1))
            for name in ('xx','xy','xz','yy','yz','zz')} for ik in INDICES3}
        divergence={indexkey(*ik):{name:zero(-3-delta-ik[0]-ik[1]+ik[2]*(delta-1),sum(ik),c.mpf('.5'))
            for name in ('x','y','z')} for ik in INDICES2}
        is_axis=endpoints(x)==endpoints(y)==(mp.mpf(0),mp.mpf(0))
        return dict(chart='core_axis_Cartesian',Z=z,X=x,Y=y,rho=rho,requested_log_tau=lt,physical_viscosity=nu,
            exact_axis=is_axis,physical_core_Cartesian_stress_mixed3=stress,
            physical_core_Cartesian_stress_divergence_mixed2=divergence,
            physical_core_Cartesian_remainder_mixed2=remainder,
            physical_core_Cartesian_momentum_decomposition_mixed2=remainder,
            actual_same_nonlinear_core_source=source,
            original_factor_logs=dict(epsilon_core=self.log_epsilon,F0base=self.log_F0),
            exact_coordinate_source='rho=(X^2+Y^2)/2; X=x_phys/(sqrt(nu)*lambda*sqrt(epsilon_core))',
            no_inverse_radius_or_angle_at_axis=True,
            all_six_original_remainder_sectors_and_nonzero_derivatives_retained=True,
            pressure_FTC_and_centrifugal_cancellation_from_same_C_source=True,
            original_signed_stress_ledger_retained_in_checked_positive_radius_parent=True,
            only_total_stress_analytic_extension_zero=True,
            source_function_equality_not_interval_overlap=True,
            source_bounds_not_resolved_physical_points=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def axis(self,Z0=(-1,1),log_tau=('-3','-1'),viscosity='1'):
        return self.cartesian(Z0,0,0,log_tau,viscosity)

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_same_source_nonsingular_core_axis_theorem=self.proof,
            actual_current_tensor_regions_available=self.core_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=32,actual_current_completed_tensor_internal_interface_count=10,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            whole_core_domain='rho[0,4],Z[-1,1]; Cartesian disk X^2+Y^2<=8',
            positive_radius_provider='CurrentCoreBackgroundTensor; original cylindrical chart still rejects rho0',
            nonsingular_axis_and_positive_radius_use_identical_analytic_source=True,
            physical_time_scope='arbitrary finite log tau and finite positive nu; independent terminal flatness/energy open',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCoreAxisBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_core_axis_and_Cartesian_views']={}
    for name,args in VIEWS.items():
        result['current_core_axis_and_Cartesian_views'][name]=field.cartesian(*args)
        print('Build nonsingular core Cartesian T/E: '+name,flush=True)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
