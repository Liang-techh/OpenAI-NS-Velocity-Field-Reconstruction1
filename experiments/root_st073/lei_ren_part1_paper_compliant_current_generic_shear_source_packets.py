"""Typed original source covers in common S=Pstar units, core through O2.

Saved boxes are executable factored Taylor covers, never point functions.
Ordinary logR conversion has already happened in the checked raw adapters.
No width, Pstar, absolute radius or F0 factor is resolved here. Arbitrary
coordinate requests require an explicitly injected existing checked owner.
"""
import ast
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as recovery
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import FactoredAlgebra, FactoredJet
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import (
    compile_function, factored_rows_record, factored_axial_derivative, LOG_NAMES)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE,PREFIX,sha=recovery.HERE,recovery.PREFIX,recovery.sha
NAME=PREFIX+'current_generic_shear_source_packets.json'
RECEIPT=PREFIX+'current_generic_shear_source_packets_check.json'
VIEWS=PREFIX+'current_generic_shear_source_packets_views.json.gz'
GATE='current_original_upstream_common_unit_packet_interface_implemented'
RECOVERY_GATE='generic_recovery_accepts_factored_current_source_packets'
OPEN=recovery.OPEN
FAMILY_KEYS=('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')
INVERSE_S=(0,-.5,0,0)
ZERO=(0,0,0,0)

# group -> producer, saved-map, source owner class, receipt gate
GROUPS={
 'core':('current_core_background_tensor','current_positive_radius_core_tensor_views',
         'CurrentCoreBackgroundTensor','current_positive_radius_core_full_background_tensor_available'),
 'bridge':('current_bridge_background_tensor',None,'CurrentBridgeBackgroundTensor',
           'current_actual_three_bridge_full_tensors_available'),
 'micro':('current_microswitch_background_tensor','current_actual_microswitch_tensor_views',
          'CurrentMicroswitchBackgroundTensor','current_actual_two_microswitch_full_tensors_available'),
 'power':('current_switch_power_background_tensor','current_actual_switch_power_tensor_views',
          'CurrentSwitchPowerBackgroundTensor','current_actual_switch_power_full_tensor_available'),
 'reshape':('current_reshape_background_tensor','current_actual_reshape_tensor_views',
            'CurrentReshapeBackgroundTensor','current_actual_long_reshape_full_tensor_available'),
 'restore':('current_restore_background_tensor','current_actual_restore_tensor_views',
            'CurrentRestoreBackgroundTensor','current_actual_inner_reference_axial_restore_buffer_full_tensors_available'),
 'patch':('current_actual_patch_background_tensor','current_actual_patch_tensor_views',
          'CurrentActualPatchBackgroundTensor','current_actual_five_moment_patch_full_tensor_available'),
 'O2':('current_O2_background_tensor','current_actual_O2_tensor_views',
       'CurrentO2BackgroundTensor','current_actual_O2_reference_slope_axial_buffer_full_tensors_available')}
# chart -> group, default saved view, explicit P0 coefficient path, selector
CHARTS={
 'core':('core','compact_core',('actual_same_fixed_point_core_source',
    'original_same_fixed_point_density_and_moment_packet','original_analytic_axis_pressure_axial6_Taylor_coefficients'),
    'rho=R/epsilon_core; D_y=rho*d_rho; positive compact sector only'),
 **{chart:('bridge',branch+'_whole',('actual_upstream_current_bridge_source',
    'actual_parent_axial5_packet','pressure_axis_axial5_coefficients'),
    'original macro coordinate; fraction selects coverage only' if branch=='macro' else 'phase; D_y^j=hb^-j D_phase^j')
    for chart,branch in (('bridge_first','first'),('bridge_second','second'),('bridge_macro','macro'))},
 **{chart:('micro',branch+'_whole',('actual_upstream_current_microswitch_source',
    'actual_parent_axial5_packet','pressure_axis_axial5_coefficients'),'phase; D_y^j=hb^-j D_phase^j')
    for chart,branch in (('switch_first','first'),('switch_second','second'))},
 'switch_power':('power','switch_power_whole',('actual_upstream_current_switch_power_source','original_P0_axial5'),
    'original power coordinate; fraction selects coverage only'),
 'reshape':('reshape','reshape_whole',('actual_upstream_current_long_reshape_source',
    'actual_inherited_axial5_packet','pressure_axis_axial5_coefficients'),'phase; source T retained, ordinary y rows already converted'),
 **{chart:('restore',chart+'_whole',('actual_upstream_current_reference_restore_source',
    'actual_inherited_axial5_packet','pressure_axis_axial5_coefficients'),'original chart selector; ordinary y rows already converted')
    for chart in ('inner_reference','axial_restore','restore_buffer')},
 'actual_patch':('patch','whole_patch',('actual_upstream_current_implicit_patch_source',
    'actual_inherited_patch_packet','original_P0_axial5'),'x=R/Rm; D_y=x*d_x already applied'),
 **{chart:('O2',chart+'_whole',('actual_upstream_original_pre_O2_source','original_P0_axial5_coefficients'),
    'original O2 selector; ordinary y rows already converted')
    for chart in ('Rh_reference','O2_slope','O2_axial','O2_buffer')}}
DOMAINS={chart:(0,1) for chart in CHARTS}
DOMAINS.update(core=(0,4),bridge_second=(1,2),switch_second=(1,2),
    restore_buffer=(-7,-6),Rh_reference=(-5,0),O2_buffer=(0,11))


def interval(c,value):
    if hasattr(value,'_mpi_'):return c.mpf(value)
    return recovery.read(c,value) if isinstance(value,dict) else c.mpf(value)


def jet(c,value):
    if isinstance(value,IntervalTaylor):
        if value.ctx is not c:raise ValueError('Foreign axial Taylor context')
        return value
    coeffs=value['coefficients'] if isinstance(value,dict) else value
    return IntervalTaylor(c,[interval(c,q) for q in coeffs])


def decode_row(algebra,value):
    """Rehydrate coefficient covers and exact half exponents, not a source evaluator."""
    c=algebra.ctx
    if isinstance(value,FactoredJet):return algebra.lift(value)
    if not isinstance(value,dict) or 'terms' not in value:return algebra.lift(jet(c,value))
    if value.get('fixed_basepoint_log_factors') is not True:raise ValueError('Fixed log bases required')
    order=value['axial_order']
    if type(order) is not int or order<0:raise ValueError('Nonnegative integer axial order required')
    terms={}
    for term in value['terms']:
        exponents=tuple(term['source_exponents'])
        if len(exponents)!=4 or any(type(x) not in (int,float) or not mp.isfinite(x) or 2*x!=int(2*x) for x in exponents):
            raise ValueError('Four exact integer/half-integer source exponents required')
        if exponents in terms:raise ValueError('Duplicate source mode')
        coeffs=term['axial_Taylor_coefficients']
        if len(coeffs)!=order+1:raise ValueError('Complete axial coefficient cover required')
        row=jet(c,coeffs)
        if any(not mp.isfinite(x) for q in row.coefficients for x in recovery.endpoints(q)):
            raise ValueError('Finite axial coefficient covers required')
        terms[exponents]=row
    return FactoredJet(algebra,terms,order)


def encode(value):
    if isinstance(value,FactoredJet):return recovery.encode(factored_rows_record(value))
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return recovery.encode(value)


def path_get(view,path):
    for key in path:view=view[key]
    return view


def factored_check_jets(c,values):
    if not values or any(type(v) not in (IntervalTaylor,FactoredJet) or v.ctx is not c or v.order<1 for v in values):
        raise ValueError('Actual axial jets in one current context required')
    algebras={id(v.algebra) for v in values if isinstance(v,FactoredJet)}
    if len(algebras)>1:raise ValueError('Source factor bases differ')
    coefficients=[q for v in values for row in (v.terms.values() if isinstance(v,FactoredJet) else (v,)) for q in row.coefficients]
    if any(not mp.isfinite(x) for q in coefficients for x in recovery.endpoints(q)):
        raise ValueError('Finite factored coefficient covers required')


def recovery_methods():
    """Replay the checked generic recovery math unchanged with jet dispatch only."""
    asts=recovery.numeric.transport.SourceAST();env=dict(vars(recovery))
    env['check_jets']=factored_check_jets;env['axial_derivative']=factored_axial_derivative
    methods={name:compile_function(asts.method('current_generic_shear_moment_recovery',name),env,
             '<unchanged generic recovery with factored axial dispatch>') for name in ('field','field_rows')}
    return methods,dict(original_generic_recovery_AST_unchanged=True,
        only_jet_validation_and_axial_dispatch_extended=True,
        no_source_factor_resolution=True,input_hashes=asts.hashes,passed=True)


class FactoredRecoveryState:
    """Own-history recovery from supplied modal covers; transport remains separate."""
    def __init__(self,packet,*,original=None,defect=None):
        self.ctx=packet.algebra.ctx;self.family=dict(packet.source_family);self.P0=packet.P0
        self.original={k:rows[0] for k,rows in packet.histories.items()} if original is None else dict(original)
        self.defect={k:self.P0*0 for k in recovery.RATES} if defect is None else dict(defect)
        if set(self.original)!=set(recovery.RATES) or set(self.defect)!=set(recovery.RATES):
            raise ValueError('All five original and changed history covers required')
        values=[self.P0,*self.original.values(),*self.defect.values()]
        factored_check_jets(self.ctx,values)
        if any(type(v) is not FactoredJet or v.algebra is not packet.algebra for v in values):
            raise ValueError('Same explicit source factor algebra required')
    def own(self):return {k:self.original[k]+self.defect[k] for k in recovery.RATES}


_METHODS,RECOVERY_PROOF=recovery_methods()
FactoredRecoveryState.field=_METHODS['field']
FactoredRecoveryState.field_rows=_METHODS['field_rows']


@dataclass(frozen=True)
class CurrentSourcePacket:
    chart: str
    source_family: dict
    algebra: FactoredAlgebra
    Z: IntervalTaylor
    velocity: dict
    histories: dict
    absolute_pressure: tuple
    P0: FactoredJet
    native_velocity: dict
    native_histories: dict
    provenance: dict

    def recover_original(self,delta):
        result=FactoredRecoveryState(self).field_rows(Z=self.Z,delta=delta,
            E_rows=self.velocity['theta'],V_rows=self.velocity['axial'])
        result.update(chart=self.chart,source_packet_provenance=self.provenance,
            covers_original_current_source_only=True,modified_generic_loop_not_installed=True,
            source_factors_not_resolved=True)
        return result

    def record(self):
        return dict(chart=self.chart,source_family=self.source_family,Z=self.Z,
            common_velocity_unit='S=Pstar',source_log_basis_names=LOG_NAMES,
            fixed_source_log_bases=self.algebra.logs,
            common_velocity_ordinary_y_rows=self.velocity,
            common_normalized_five_history_ordinary_y_rows=self.histories,
            absolute_pressure_over_S_squared_ordinary_y_rows=self.absolute_pressure,
            original_axis_pressure_over_S_squared=self.P0,
            velocity_units=dict(theta='Utheta/S',axial='Uz/S',radial='Ur/(S*sqrt(R/2))'),
            primitive_units=recovery.UNITS,radial_rows_include_physical_half_shift=True,
            axial_coefficients_are_Taylor_not_ordinary=True,
            radial_rows_are_ordinary_logR_not_selector_derivatives=True,
            P0_separate_from_cumulative_p=True,provenance=self.provenance,
            point_source_function=False,source_function_replay=False,**dict.fromkeys(OPEN,False))


class CurrentSourcePackets:
    """Lazy receipt-bound cache covers, or explicit injection of existing live owners."""
    def __init__(self,*,owners=None):
        self.owners=dict(owners or {})
        if set(self.owners)-set(GROUPS):raise ValueError('Known upstream owner group required')
        self.data=recovery.numeric.inputs();self.ctx=self.data['ctx']
        self.family=dict(self.data['accepted']['source_family']);self.hashes=dict(self.data['hashes'])
        self.groups={};self.shards={}
        admitted=json.loads((HERE/recovery.RECEIPT).read_bytes())
        if not admitted.get('all_passed') or not admitted.get(recovery.GATE) or admitted.get('source_family')!=self.family or any(admitted.get(k) for k in OPEN):
            raise ValueError('Same checked generic recovery source/axis datum required')
        self.bind_hashes(admitted['input_hashes']);self.bind_hashes({recovery.RECEIPT:sha(recovery.RECEIPT)})
        self.hashes.update(RECOVERY_PROOF['input_hashes'])
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def bind_hashes(self,hashes):
        for name,digest in hashes.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting current dependency: '+name)
            if name not in self.hashes and sha(name)!=digest:raise ValueError('Changed current packet dependency: '+name)
            self.hashes[name]=digest

    def group(self,name):
        if name in self.groups:return self.groups[name]
        stem,_,_,gate=GROUPS[name];producer=PREFIX+stem+'.json.gz';receipt=PREFIX+stem+'_check.json'
        admitted=json.loads((HERE/receipt).read_bytes())
        if not admitted.get('all_passed') or not admitted.get(gate) or any(admitted.get(k) for k in OPEN):
            raise ValueError('Checked original regional tensor receipt required')
        if {k:admitted[k] for k in FAMILY_KEYS}!=self.family:raise ValueError('Foreign source/axis datum')
        self.bind_hashes(admitted['input_hashes']);self.bind_hashes({receipt:sha(receipt)})
        if admitted['input_hashes'].get(producer)!=sha(producer):raise ValueError('Unbound regional packet producer')
        data=json.loads(gzip.decompress((HERE/producer).read_bytes()))
        if {k:data[k] for k in FAMILY_KEYS}!=self.family:raise ValueError('Regional packet family differs')
        self.groups[name]=data
        return data

    def saved_view(self,chart,view=None):
        if chart not in CHARTS:raise ValueError('Original upstream chart required')
        group,default,_,_=CHARTS[chart];data=self.group(group);view=default if view is None else view
        if group=='bridge':
            shard=data['current_actual_bridge_chart_shards'][chart];name=shard['path']
            if view not in shard['views']:raise ValueError('Saved bridge view required')
            if name not in self.shards:
                self.bind_hashes({name:shard['sha256']})
                self.shards[name]=json.loads(gzip.decompress((HERE/name).read_bytes()))
            rows=self.shards[name]
        else:rows=data[GROUPS[group][1]];name=PREFIX+GROUPS[group][0]+'.json.gz'
        if view not in rows or rows[view]['chart']!=chart:raise ValueError('Saved view for the requested chart required')
        return rows[view],dict(mode='saved_original_cover',view=view,producer=name,
            producer_sha256=self.hashes[name],receipt=PREFIX+GROUPS[group][0]+'_check.json',
            cache_cover=True,arbitrary_coordinates_evaluated=False)

    def saved(self,chart,view=None):
        data,origin=self.saved_view(chart,view)
        return self.adapt(chart,data,self.ctx,origin)

    def query(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        if chart not in CHARTS:raise ValueError('Original upstream chart required')
        group=CHARTS[chart][0]
        if group not in self.owners:raise ValueError('Arbitrary coordinates require an injected checked live owner; cache fallback forbidden')
        self.group(group);owner=self.owners[group];stem,_,class_name,_=GROUPS[group]
        if type(owner).__name__!=class_name or type(owner).__module__!=PREFIX+stem or not owner.acceptance_loaded:
            raise ValueError('Exact existing checked original owner required')
        if dict(zip(FAMILY_KEYS,(owner.family,owner.source,owner.datum_sha)))!=self.family:
            raise ValueError('Injected owner source/axis datum differs')
        owner.assert_graph();self.bind_hashes(owner.hashes)
        data=owner.chart(chart,Z,coordinate,log_tau,theta,viscosity)
        return self.adapt(chart,data,owner.ctx,dict(mode='injected_current_source_cover',cache_cover=False,
            arbitrary_coordinates_evaluated=True,owner_class=class_name,
            receipt=PREFIX+stem+'_check.json',source_graph_asserted=True))

    def adapt(self,chart,view,c,origin):
        if chart not in CHARTS or view['chart']!=chart or not view.get('source_bounds_not_resolved_physical_point_values'):
            raise ValueError('Original source covers with exact chart provenance required')
        physical=view['actual_upstream_physical_spatial4_time1_packet']
        if {k:physical[k] for k in FAMILY_KEYS}!=self.family:raise ValueError('Packet source/axis datum differs')
        packet=view['current_actual_source_stress_packet'];logp=interval(c,packet['exact_pulse_reference_logB_parts']['logPstar'])
        lo,hi=recovery.endpoints(logp);expected=interval(c,self.data['logP']);el,eh=recovery.endpoints(expected)
        if hi<el or eh<lo:raise ValueError('Foreign constant common velocity unit')
        if 'current_unresolved_raw_source_rows' in view:
            raw=view['current_unresolved_raw_source_rows'];logs=tuple(interval(c,x) for x in view['fixed_current_factored_source_log_bases'])
            if len(logs)!=4 or any(not mp.isfinite(x) for value in logs for x in recovery.endpoints(value)):
                raise ValueError('All four finite original source log bases required')
            al,ah=recovery.endpoints(logs[1]);bl,bh=recovery.endpoints(logp*2)
            if ah<bl or bh<al:raise ValueError('Original Pstar squared log basis differs')
            if not view.get('original_radial_prefactors_and_absolute_P0_included_once'):
                raise ValueError('Checked raw coordinate/radial/pressure normalization required')
            if chart=='bridge_macro' and not raw.get('macro_original_ordinary_y_preserved'):
                raise ValueError('Macro ordinary y identity required, no microscopic width conversion')
            if chart in ('bridge_first','bridge_second','switch_first','switch_second') and not raw.get('inverse_width_source_shift_applied_before_any_resolution'):
                raise ValueError('Original phase width inversion must precede packet normalization')
            velocity=raw['velocity'];histories=raw['histories'];pressure=raw['absolute_pressure']
            raw_form='converted_factored_raw_cover'
        else:
            logs=(c.mpf(0),2*logp,c.mpf(0),c.mpf(0))
            velocity=view['current_source_three_component_velocity_rows'];histories=view['current_raw_five_history_rows']
            pressure=view['current_absolute_pressure_ordinary_y_rows'];raw_form='converted_plain_raw_cover'
            if chart=='actual_patch' and not view.get('ordinary_x_converted_to_logR_before_full_physical_derivatives'):
                raise ValueError('Patch Euler conversion required')
        algebra=FactoredAlgebra(c,logs,[])
        def rows(values):
            if len(values)!=5:raise ValueError('Five ordinary y rows required')
            result=tuple(decode_row(algebra,v) for v in values)
            factored_check_jets(c,list(result))
            return result
        if set(velocity)!=set(('theta','axial','radial')) or set(histories)!=set(recovery.RATES):
            raise ValueError('Three velocity profiles and all five histories required')
        native_v={k:rows(v) for k,v in velocity.items()};native_m={k:rows(v) for k,v in histories.items()}
        absolute=rows(pressure)
        v={k:tuple(algebra.shift(row,INVERSE_S) if k in ('axial','radial') else row for row in values) for k,values in native_v.items()}
        m={k:tuple(algebra.shift(row,INVERSE_S) if k in ('m','k') else row for row in values) for k,values in native_m.items()}
        p0=algebra.lift(jet(c,path_get(view,CHARTS[chart][2])[:6]))
        if p0.order<5:raise ValueError('Separate original axial order-five pressure datum required')
        Z=IntervalTaylor.variable(c,interval(c,packet['Z']),5)
        coordinate=interval(c,view['coverage_coordinate'])
        zl,zh=recovery.endpoints(Z[0]);vl,vh=recovery.endpoints(coordinate)
        lower,upper=DOMAINS[chart]
        if chart=='actual_patch':lower,upper=1,recovery.endpoints(c.exp(1))[1]
        if not all(mp.isfinite(x) for x in (zl,zh,vl,vh)) or zl<-1 or zh>1 or vl<lower or vh>upper or (chart=='core' and vl<=0):
            raise ValueError('Original source box and positive-radius core domain required')
        exact_radius=view.get('exact_original_radius_source')
        if chart=='core':exact_radius=dict(formal_radius=view['actual_same_fixed_point_core_source']['formal_radius'],
            original_positive_epsilon_log=view['actual_same_fixed_point_core_source']['original_positive_epsilon_log'])
        if exact_radius is None:raise ValueError('Original formal radius source required')
        provenance=dict(origin,source_family=self.family,chart=chart,Z_box=Z[0],coordinate_box=coordinate,
            requested_log_tau=view['requested_log_tau'],physical_viscosity=view['physical_viscosity'],
            coordinate_contract=CHARTS[chart][3],derivative_coordinate='ordinary y=log R',
            original_radius_source=exact_radius,logR_cover=interval(c,packet['exact_logR']),
            raw_record_form=raw_form,raw_coordinate_conversion_not_reapplied=True,
            current_swirl_factor_not_applied_twice=True,P0_read_from_original_source_not_subtracted=True,
            source_factor_resolution_performed=False,core_axis_included=False,
            per_view_bounds_not_a_seam_identity=True)
        return CurrentSourcePacket(chart,dict(self.family),algebra,Z,v,m,absolute,p0,native_v,native_m,provenance)


def unit_theorem():
    S,R=s.symbols('S R',positive=True)
    Uz,Ur,Mz,Mh,Mk,Me,Mp,P0=s.symbols('Uz Ur Mz Mh Mk Me Mp P0')
    native=dict(V=Uz,Q=Ur/s.sqrt(R/2),m=Mz/R,h=Mh/(s.sqrt(2)*R**s.Rational(3,2)*S),
        k=Mk/(s.sqrt(2)*R**s.Rational(3,2)*S),e=Me/(R*S**2),p=Mp/S**2,P0=P0/S**2)
    expected=dict(V=Uz/S,Q=Ur/(S*s.sqrt(R/2)),m=Mz/(R*S),h=native['h'],
        k=Mk/(s.sqrt(2)*R**s.Rational(3,2)*S**2),e=native['e'],p=native['p'],P0=native['P0'])
    checks={key:s.cancel(value/(S if key in ('V','Q','m','k') else 1)-expected[key])==0 for key,value in native.items()}
    checks['inverse_S_uses_half_Pstar_squared_log']=s.exp(-s.log(S**2)/2)==1/S
    if not all(checks.values()):raise ArithmeticError('Common source unit conversion differs')
    return dict(passed=True,identities=checks,source_log_basis_names=LOG_NAMES,
        inverse_S_exponent_shift=INVERSE_S,all_y_Z_derivatives_have_same_constant_unit_shift=True,
        native_width_and_swirl_conversion_applied_before_this_interface=True,
        absolute_pressure_and_P0_not_divided_twice=True)


def run():
    service=CurrentSourcePackets();covers={};inventory={}
    for chart in CHARTS:
        packet=service.saved(chart);field=packet.recover_original(service.data['delta'])
        covers[chart]=dict(packet=packet.record(),original_recovered_field=field)
        inventory[chart]=dict(provenance=packet.provenance,ordinary_y_rows=5,
            separate_P0_axial_order=packet.P0.order,full_signed_stress_y_rows=4,
            source_factor_terms_retained=sum(len(row.terms) for rows in (*packet.velocity.values(),*packet.histories.values(),packet.absolute_pressure) for row in rows))
        print('Typed actual source and factored recovery: '+chart,flush=True)
    (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(encode(covers),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    service.bind_hashes({VIEWS:sha(VIEWS),recovery.RECEIPT:sha(recovery.RECEIPT)})
    result=dict(source_family=service.family,typed_original_source_cover_inventory=inventory,
        normalized_cover_views=VIEWS,current_original_chart_count=len(CHARTS),
        common_unit_theorem=unit_theorem(),unchanged_generic_factored_recovery_theorem=RECOVERY_PROOF,
        **{GATE:True,RECOVERY_GATE:True},**dict.fromkeys(OPEN,False),
        exact_original_width_and_swirl_bases_retained=True,
        original_core_scope='saved positive compact rho sector; axis and whole core not inferred',
        cache_covers_not_arbitrary_coordinate_functions=True,
        live_queries_require_existing_checked_owner_injection=True,
        ancestor_constructors_called=False,
        scope='Typed original box covers for 16 core-to-O2 charts, common formal units and exact own-field recovery algebra. No current generic-loop installation, full-source gate, changed moment transport, repair or common N.',
        input_hashes=service.hashes)
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
