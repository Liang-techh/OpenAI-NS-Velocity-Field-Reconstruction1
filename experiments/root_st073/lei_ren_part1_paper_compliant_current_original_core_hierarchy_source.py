"""Original real core functions and the genuine order-one known forcing.

The radial Taylor coefficients of the leading fixed point are not temporal
hierarchy coefficients. No positive-order solution or complex-domain theorem
is manufactured by this source adapter.
"""
import copy
from dataclasses import dataclass
import gzip
import ast
import hashlib
import inspect
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_native_constants as native
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source as original
import lei_ren_part1_paper_compliant_core_physical_field as core_module
import lei_ren_part1_paper_compliant_current_core_first_interface as first
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame import source_assignment

HERE,PREFIX,sha=original.HERE,original.PREFIX,original.sha
NAME=PREFIX+'current_original_core_hierarchy_source.json.gz'
RECEIPT=PREFIX+'current_original_core_hierarchy_source_check.json'
GATE='current_original_core_real_functions_and_genuine_n1_known_forcing_installed'
OPEN=('coefficient_F1_Uz1_P1_solved','positive_order_five_moment_repair_closed',
      'common_holomorphic_extension_assumption14_1_certified','temporal_recursion',
      'regional_admissible_stress_certified','flat_remainder_certified',
      'full_corrected_Navier_Stokes_solution')
FIELDS=('delta','logLambda','Lambda','epsilon','logP','logC','j','h','sigma',
        'correction','phi_floor','phi_ceiling','Gbar')
SAMPLES={'whole_real_core':(['0','4.1'],['-1','1']),
         'whole_axis':('0',['-1','1']),'axis_offcenter':('0','.371'),
         'inner_offcenter':('2','.371'),'beyond_Ra':('4.1','.371')}


class OriginalCoreContext:
    """Independent cache; unchanged callbacks on the canonical source tuple."""
    axis_inputs=core_module.CompliantCorePhysicalField.axis_inputs
    normalized_jets=core_module.CompliantCorePhysicalField.normalized_jets
    profiles=core_module.CompliantCorePhysicalField.profiles

    def __init__(self,source):
        for key in FIELDS:setattr(self,key,getattr(source.core,key))
        self.ctx=source.c;self.datum=source.core.datum;self.cache={}


@dataclass(frozen=True,eq=False)
class OriginalCoreHierarchySource:
    coordinate:str='rho=Lambda*R'


def ends(value):return original.ep(value)


def contains(outer,inner):
    a,b=ends(outer);x,y=ends(inner)
    return a<=x and y<=b


def axial_second(c,grid,rho,z,delta,weight):
    """Z_(a-1+delta) Z_a, with R d_R = rho d_rho.

    Grid entries are ordinary mixed derivatives of the true profile, with
    only its common F0 amplitude factored out when appropriate.
    """
    at=lambda i,k:grid[core_module.gridkey(i,k)]
    f,fz,fr=at(0,0),at(0,1),at(1,0)
    d=1-z*z;L=1-delta*z*z
    numerator=weight*z*f+d*fz-2*z*rho*fr
    value=numerator/L
    derivative_z=((weight*f+(weight-2)*z*fz+d*at(0,2)-2*rho*fr
                   -2*z*rho*at(1,1))*L+2*delta*z*numerator)/(L*L)
    derivative_rho=((weight-2)*z*fr+d*at(1,1)-2*z*rho*at(2,0))/L
    return ((weight-1+delta)*z*value+d*derivative_z-2*z*rho*derivative_rho)/L


def regular_radial_forcing(core,profile):
    """-Omega0/(2R) using V0=R Q, including R=0 directly."""
    g=profile['ordinary_mixed_profile_grids'];rho,z=profile['rho'],profile['Z']
    at=lambda i,k:g[core_module.Q][core_module.gridkey(i,k)]
    q,qz=at(0,0),at(0,1);rho_qr=rho*at(1,0)
    uz=g[core_module.V][core_module.gridkey(0,0)]
    L=1-core.delta*z*z;d=1-z*z
    time=(q+(1-core.delta)*z*qz/2+rho_qr)/L
    radial_transport=q*(q/2+rho_qr)
    axial_transport=uz*(-2*z*q+d*qz-2*z*rho_qr)/L
    radial_viscosity=-core.Lambda*(4*at(1,0)+2*rho*at(2,0))
    return -(time+radial_transport+axial_transport+radial_viscosity)/2


class CurrentOriginalCoreHierarchySource:
    @source_precision
    def __init__(self,constants,source=None,require_checked=True):
        if type(constants) is not native.CurrentOriginalRpNativeConstants or not constants.acceptance_loaded:
            raise ValueError('Accepted live original Rp constant-function owner required')
        self.constants=constants
        self.source=source if source is not None else original.OriginalWholeZBridgeSource(500)
        if type(self.source) is not original.OriginalWholeZBridgeSource:
            raise ValueError('Canonical live original whole-Z core source required')
        self.family=copy.deepcopy(constants.identity);self.ctx=c=self.source.c
        got=dict(implicit_source_sha256=self.source.source,
                 actual_five_defect_family_sha256=self.source.family,datum_enclosure_sha256=self.source.datum)
        if got!=self.family:raise ValueError('Core and current original Rp source families differ')
        ps=self.source.records['pressure_source']['compliant_source']
        pressure=constants.frame.reports[native.frame.native.NAME]['current_P0_source_binding']
        if pressure['analytic_definition']!=ps['implicit_source_definition']:
            raise ValueError('Actual inner and outer independent P0 functions differ')
        if not pressure['passed'] or not all(pressure[key] for key in (
                'same_fourteen_continuous_atoms_and_analytic_flatten_function',
                'current_order6_first_six_and_native_order5_same_defining_rows',
                'common_P0_kept_independent_of_cumulative_pressure')):
            raise ValueError('Accepted exact P0 function/callable identity required')
        callable_hash=hashlib.sha256(ast.dump(ast.parse(inspect.getsource(
            self.source.core.datum.normalized_jets.__func__).strip())).encode()).hexdigest()
        if callable_hash!=pressure['inherited_normalized_jets_AST_sha256']:
            raise ValueError('Live original normalized pressure callable changed')
        self.P0_binding=copy.deepcopy(pressure)
        self.hashes=dict(constants.hashes)
        def bind(rows):
            for name,digest in rows.items():original.bind(self.hashes,name,digest)
        native_receipt=json.loads((HERE/native.RECEIPT).read_bytes())
        if not native_receipt['all_passed'] or not all(native_receipt[key] for key in native.GATES):
            raise ValueError('Accepted original Rp defining-function receipt required')
        bind(native_receipt['input_hashes'])
        for name in (native.NAME,native.RECEIPT,original.NAME,first.NAME):
            original.bind(self.hashes,name,sha(name))
        bind(self.source.hashes)
        self.parents={}
        for name,gate in ((original.RECEIPT,original.GATE),
                          (first.RECEIPT,first.GATE)):
            receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Accepted core source/interface required '+name)
            for key,wanted in self.family.items():
                if receipt[key]!=wanted:raise ValueError('Core receipt belongs to another source '+name)
            bind(receipt['input_hashes']);original.bind(self.hashes,name,sha(name));self.parents[name]=receipt
        self.core=OriginalCoreContext(self.source)
        scale_receipt=PREFIX+'current_core_recurrence_source_check.json'
        scale=json.loads((HERE/scale_receipt).read_bytes())
        if not scale['all_passed'] or not scale['pressure_and_seed_bindings']['Lambda_times_core_epsilon_exactly_one']:
            raise ValueError('Accepted exact Lambda epsilon identity required')
        if any(scale[key]!=value for key,value in self.family.items()):
            raise ValueError('Core radial scale proof belongs to another source')
        bind(scale['input_hashes']);original.bind(self.hashes,scale_receipt,sha(scale_receipt))
        self.parameter_bindings=[source_assignment('lei_ren_part1_paper_compliant_core_transfer.py',
            'run',target,expression) for target,expression in
            (('logLambda','4*logP+1000'),('lam','ctx.exp(logLambda)'),('eps','ctx.exp(-logLambda)'))]
        exact_logP=c.exp(40)+11
        checks=dict(original_logP_contains_defining_function=contains(self.core.logP,exact_logP),
            original_logLambda_contains_defining_function=contains(self.core.logLambda,4*exact_logP+1000),
            original_delta_contains_defining_function=contains(self.core.delta,c.exp(-4*exact_logP-30)),
            original_Lambda_contains_defining_function=contains(self.core.Lambda,c.exp(4*exact_logP+1000)),
            original_core_epsilon_contains_defining_function=contains(self.core.epsilon,c.exp(-4*exact_logP-1000)),
            original_Lambda_epsilon_interval_contains_one=contains(self.core.Lambda*self.core.epsilon,c.mpf(1)),
            positive_original_delta=ends(self.core.delta)[0]>0,
            positive_original_core_epsilon=ends(self.core.epsilon)[0]>0,
            independent_P0_definition_identified=True)
        checks['exact_Lambda_times_core_epsilon_identity_inherited']=True
        if not all(checks.values()):raise ValueError('Original core parameter binding failed '+str(checks))
        self.parameter_checks=checks;self._snapshot=self._parameters()
        self._datum_snapshot=self._datum_parameters()
        self._amplitude_parameters={key:value._mpi_ for key,value in self.source.amplitude.params.items()}
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self._fields={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Original core hierarchy source receipt/scope differs')
            for name in (Path(__file__).name,Path(__file__).stem+'_check.py',NAME):
                if receipt['input_hashes'].get(name)!=sha(name):raise ValueError('Unbound original core hierarchy source '+name)
            bind(receipt['input_hashes']);self.acceptance_loaded=True

    def _parameters(self):
        return {key:(getattr(self.core,key)._mpi_,getattr(self.source.core,key)._mpi_) for key in FIELDS}

    def _datum_parameters(self):
        datum=self.source.core.datum
        return dict(values={key:getattr(datum,key)._mpi_ for key in ('m0','m2','rho','flatten_complex_upper')},
            flatten_mass=datum.stages['z_flatten']['mass']._mpi_,logP=datum.parameters.logPstar._mpi_,
            source=datum.source_sha,datum=datum.datum_sha)

    def assert_graph(self):
        checks=dict(unchanged_core_parameter_tuple=self._snapshot==self._parameters(),
            same_actual_pressure=self.core.datum is self.source.core.datum is self.source.admission.pressure,
            unchanged_original_pressure_atoms=self._datum_snapshot==self._datum_parameters(),
            original_profiles_callback=self.core.profiles.__func__ is core_module.CompliantCorePhysicalField.profiles,
            original_axis_callback=self.core.axis_inputs.__func__ is core_module.CompliantCorePhysicalField.axis_inputs,
            original_nonlinear_bounds_callback=self.core.normalized_jets.__func__ is core_module.CompliantCorePhysicalField.normalized_jets,
            independent_core_cache=self.core.cache is not self.source.cache,
            same_interval_context=self.core.ctx is self.ctx is self.source.c,
            unchanged_anchored_amplitude_parameters=self._amplitude_parameters=={
                key:value._mpi_ for key,value in self.source.amplitude.params.items()},
            actual_amplitude_and_rebuild_core=self.source.amplitude.core is self.source.rebuild.core is self.source.core,
            original_source_family=(self.source.source,self.source.family,self.source.datum)==tuple(self.family[key]
                for key in ('implicit_source_sha256','actual_five_defect_family_sha256','datum_enclosure_sha256')),
            accepted_Rp_function_source=self.constants.acceptance_loaded,
            identical_source_family=self.constants.identity==self.family)
        if not all(checks.values()):raise ValueError('Original core function graph differs '+str(checks))
        return checks

    @source_precision
    def evaluate(self,rho,Z):
        self.assert_graph();c=self.ctx;rho,z=c.mpf(rho),c.mpf(Z)
        if ends(rho)[0]<0 or ends(rho)[1]>ends(c.mpf('4.1'))[1] or ends(z)[0]<-1 or ends(z)[1]>1:
            raise ValueError('Original core domain rho[0,4.1], Z[-1,1] required')
        profile=self.core.profiles(rho,z)
        amplitude=self.source.amplitude.evaluate(z)
        g=profile['ordinary_mixed_profile_grids'];logF=amplitude['logF0']
        fields={}
        for name,key,base_log in (('F_0',core_module.F,logF),('Uz_0',core_module.V,c.mpf(0)),
                                  ('V_0_over_R',core_module.Q,c.mpf(0))):
            fields[name]={}
            for i in range(5):
                for k in range(5-i):
                    index=core_module.gridkey(i,k)
                    fields[name][index]=dict(coefficient=g[key][index],log_scale=base_log+i*self.core.logLambda,
                        exact_scale_formula=('logF0+' if name=='F_0' else '')+str(i)+'*logLambda',
                        radial_derivative=i,ordinary_Z_derivative=k)
        fields['P_0']={}
        for i in range(5):
            for k in range(5-i):
                index=core_module.gridkey(i,k)
                fields['P_0'][index]=dict(axis_datum=dict(coefficient=g[core_module.PD][index],
                    log_scale=2*self.core.logP+i*self.core.logLambda,exact_scale_formula='2logP+'+str(i)+'*logLambda'),
                    centrifugal_increment=dict(coefficient=g[core_module.PI][index],
                    log_scale=2*logF+(i-1)*self.core.logLambda,exact_scale_formula='2logF0+'+str(i-1)+'*logLambda'))
        forcing=dict(N1theta=dict(coefficient=-axial_second(c,g[core_module.F],rho,z,self.core.delta,-2-self.core.delta),
                                 log_scale=logF,exact_scale_formula='logF0'),
                     N1z=dict(coefficient=-axial_second(c,g[core_module.V],rho,z,self.core.delta,-1-self.core.delta),
                              log_scale=c.mpf(0),exact_scale_formula='0'),
                     N1p=dict(coefficient=regular_radial_forcing(self.core,profile),
                              log_scale=c.mpf(0),exact_scale_formula='0'))
        result=dict(rho=rho,Z=z,actual_leading_function_derivatives=fields,genuine_order_one_known_forcing=forcing,
            exact_weights=dict(e1='2delta',b0='-2-delta',c0='-1-delta',p0='-2-2delta',
                b1='-2+delta',c1='-1+delta',p1='-2'),
            regular_recovery='Q=[2Z Uz0-(1-delta)Z(Mz0/R)-(1-Z^2)d_Z(Mz0/R)]/(1-delta Z^2)',
            radial_forcing_has_no_division_by_R=True,pressure_axis_datum_retained_separately=True,
            radial_derivatives_converted_by_Lambda_power=True,all_output_jets_are_ordinary_derivatives=True,
            ordinary_radial_Taylor_order_not_used_as_hierarchy_order=True,
            real_common_domain=dict(rho=['0','4.1'],Z=['-1','1'],R_in='4.1/Lambda',R_a='4/Lambda'),
            exact_original_parameter_recipes=dict(logP='exp(40)+11',logdelta='-4logP-30',
                logLambda='4logP+1000',epsilon_core='exp(-logLambda)',epsilon_pressure='delta/1000',
                logF0='-selected_logCstar-Lambda*G(Z)'),
            original_P0_analytic_definition=copy.deepcopy(self.source.records['pressure_source']['compliant_source']['implicit_source_definition']),
            original_P0_exact_function=self.P0_binding['exact_analytic_function'],
            exact_P0_callable_and_hydration_AST_identity_inherited=True,
            outputs_are_enclosures_of_original_functions_not_selected_values=True,
            **dict.fromkeys(OPEN,False))
        value=OriginalCoreHierarchySource()
        self._fields[id(value)]=(value,profile,result,copy.deepcopy(first.encode(original.serialized(result))))
        return value

    @source_precision
    def report(self,value):
        entry=self._fields.get(id(value))
        if type(value) is not OriginalCoreHierarchySource or entry is None or entry[0] is not value:
            raise ValueError('Live core hierarchy source packet issued by this owner required')
        self.assert_graph()
        if first.encode(original.serialized(entry[2]))!=entry[3]:raise ValueError('Original source values changed')
        return dict(source_family=copy.deepcopy(self.family),**copy.deepcopy(entry[2]),**{GATE:self.acceptance_loaded})


@source_precision
def run(constants,source=None):
    began=time.monotonic();owner=CurrentOriginalCoreHierarchySource(constants,source,require_checked=False)
    values={name:owner.evaluate(rho,z) for name,(rho,z) in SAMPLES.items()}
    result=dict(source_family=owner.family,actual_core_function_and_forcing_samples={name:owner.report(value) for name,value in values.items()},
        original_parameter_checks=owner.parameter_checks,original_parameter_AST_bindings=owner.parameter_bindings,
        original_core_function_graph=owner.assert_graph(),input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    data=json.dumps(first.encode(original.serialized(result)),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('ORIGINAL_CORE_N1_FORCING_SOURCE_READY',len(values),flush=True)
    return owner,values
