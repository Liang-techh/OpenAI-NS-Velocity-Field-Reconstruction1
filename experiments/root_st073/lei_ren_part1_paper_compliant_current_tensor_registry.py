"""One checked graph, 33 full tensor native routes and 46 source trace routes.

This registry preserves the original full views and signed sectors. Regional
coordinates are explicit; inverse physical-coordinate location remains open.
"""
import ast
import hashlib
import importlib
import inspect
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_angular_internal_background_tensor import (
    CurrentAngularInternalBackgroundTensor,canonical_angular_groups,
    HERE,PREFIX,OPEN as BACKGROUND_OPEN,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,SEAMS as ANGULAR_SEAMS)
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import canonical_tensor_groups
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor import EDGES as PULSE_EDGES

NAME=PREFIX+'current_tensor_registry.json'
RECEIPT=PREFIX+'current_tensor_registry_check.json'
GATES=('current_33_region_full_tensor_native_registry_available',
       'current_32_adjacent_and_14_internal_tensor_trace_routes_bound')
OPEN=BACKGROUND_OPEN+('current_global_physical_tensor_coordinate_locator_certified',
    'current_global_tensor_physical_cover_certified')

# Paths are relative to the previous milestone's checked root. No new owners.
CHAIN=('core_tensor','bridge_tensor','micro_tensor','power_tensor','reshape_tensor',
    'restore_tensor','patch_tensor','o2_tensor','o3_tensor','entrance_incoming_tensor',
    'main_tensor','gap_tensor','end_tensor','flatten_power','heat_tensor','o7','entry_source','angular')
ALIASES=('core','bridge','micro','power','reshape','restore','patch','o2','o3',
    'incoming','main','gap','end','flatten','heat','o7','entry','angular')
ROUTES=(
    ('core_positive_radius','core','chart','core','rho in (0,4]; exact rho=0 uses the same Cartesian axis source'),
    ('bridge_first','bridge','chart','bridge_first','phase [0,1]; ordinary logR derivatives include hb^-j'),
    ('bridge_second','bridge','chart','bridge_second','phase [1,2]; ordinary logR derivatives include hb^-j'),
    ('bridge_macro','bridge','chart','bridge_macro','fraction [0,1] labels coverage; original ordinary logR rows unchanged'),
    ('switch_first','micro','chart','switch_first','phase [0,1]'),
    ('switch_second','micro','chart','switch_second','phase [1,2]'),
    ('switch_power','power','chart','switch_power','fraction [0,1]'),
    ('reshape','reshape','chart','reshape','phase [0,1]'),
    ('inner_reference','restore','chart','inner_reference','phase [0,1]'),
    ('axial_restore','restore','chart','axial_restore','phase [0,1]'),
    ('restore_buffer','restore','chart','restore_buffer','offset [-7,-6]'),
    ('actual_patch','patch','chart','actual_patch','x=R/Rm in [1,e]; named Rh/support source limits remain available'),
    ('Rh_reference','o2','chart','Rh_reference','offset [-5,0]'),
    ('O2_slope','o2','chart','O2_slope','phase [0,1]'),
    ('O2_axial','o2','chart','O2_axial','phase [0,1]; original axial phase keyword'),
    ('O2_buffer','o2','chart','O2_buffer','buffer offset [0,11]; original buffer_offset keyword'),
    ('O3_slope_mu','o3','chart','O3_slope_mu','offset [0,1]; variable log-amplitude slope'),
    ('O3_power','incoming','chart','O3_power','phase [0,1]'),
    ('pulse_entrance','incoming','chart','pulse_entrance','xi [0,.02]'),
    ('pulse_main','main','chart','pulse_main','xi [.02,10]'),
    ('pulse_exit','main','chart','pulse_exit','xi [10,11]'),
    ('pulse_gap','gap','chart','pulse_gap','xi [11,12]'),
    ('pulse_gap_end','gap','chart','pulse_gap_end','phase [0,1]; d in [4mu,1], exact original s=-d/mu'),
    ('pulse_end','end','end',None,'s [-4,0]'),
    ('flatten','flatten','chart','flatten','offset [0,100]'),
    ('outer_power','flatten','chart','outer_power','phase [0,1]; y=(Lrel-4)*phase'),
    ('outer_angular','angular','angular',None,'s [-4,0]; exact support limits use the internal trace routes'),
    ('steep_entry','entry','entry',None,'offset [0,1]'),
    ('steep_power','o7','chart','steep_power','phase [0,1]'),
    ('steep_exit','o7','chart','steep_exit','phase [0,1]'),
    ('waiting','o7','chart','waiting','offset [0,1]'),
    ('heat_collar','heat','chart','heat_collar','offset [0,3]'),
    ('heat_exterior','heat','chart','heat_exterior','offset >=3; separate same-source unbounded Gamma identity'))

SOURCE_COORDINATES=(
    'rho','bridge_phase','bridge_phase','macro_coverage_fraction',
    'microscopic_switch_phase','microscopic_switch_phase','power_coverage_fraction',
    'reshape_phase','reference_phase','restoration_phase','postrestore_offset',
    'R_over_Rm','reference_offset','slope_phase','axial_phase','axial_buffer_offset',
    'slope_mu_offset','incoming_power_phase','pulse_xi','pulse_xi','pulse_xi','pulse_xi',
    'reciprocal_gap_phase','end_s','flatten_offset','outer_power_phase','angular_s',
    'steep_entry_offset','steep_power_phase','steep_exit_phase','waiting_offset','heat_offset','heat_offset')

def derivative_pullback(region,alias,coordinate):
    if region=='core_positive_radius':return 'D_y=rho*d_rho via the checked Euler/Stirling order4 source conversion; rho0 uses the separate Cartesian core pullback'
    if region in ('bridge_first','bridge_second'):return 'D_y^j=hb^-j*d_bridge_phase^j in the original bridge adapter; registry applies no second width factor'
    if region=='bridge_macro':return 'Original ordinary-logR rows retained; macro fraction selects coverage only, with no inverse macro-length factor'
    if region=='pulse_gap_end':return 'Exact s=-d/mu with d in [4mu,1]; original ordinary-logR mixed rows recovered before enclosure, not phase derivatives'
    if region=='outer_power':return 'Original y=(Lrel-4)*phase; full mixed rows already differentiate ordinary y before the physical pullback'
    if region=='actual_patch':return 'Original x=R/Rm source rows converted to ordinary logR by the checked full patch adapter'
    return 'The original '+alias+' operator maps '+coordinate+' to ordinary logR and physical mixed rows; source method AST fixes its width/offset conversion, registry applies no derivative rescaling'

# Each entry points at actual admitted rows, not a cumulative count alone.
ADJACENT=(
    ('core','core_bridge','current_core_bridge_completed_tensor_interface_rows',None),
    ('bridge','bridge_first_second','current_actual_bridge_tensor_interface_rows','bridge_first_second'),
    ('bridge','bridge_second_macro','current_actual_bridge_tensor_interface_rows','bridge_second_macro'),
    ('bridge','bridge_switch_R100','current_actual_bridge_tensor_interface_rows','bridge_switch_R100'),
    ('micro','first_second','current_actual_phase1_R2_tensor_interface_rows','first_second'),
    ('micro','second_power','current_actual_phase1_R2_tensor_interface_rows','second_power'),
    ('power','power_reshape','current_actual_R110_tensor_interface_rows','power_reshape'),
    ('reshape','reshape_reference','current_actual_Rsh_tensor_interface_rows','reshape_reference'),
    ('restore','reference_restore','current_actual_three_restore_tensor_interface_rows','reference_restore'),
    ('restore','restore_buffer','current_actual_three_restore_tensor_interface_rows','restore_buffer'),
    ('restore','buffer_patch','current_actual_three_restore_tensor_interface_rows','buffer_patch'),
    ('patch','patch_Rh','current_actual_patch_Rh_and_six_support_tensor_interface_rows','patch_Rh'),
    ('o2','reference_slope','current_actual_four_O2_tensor_interface_rows','reference_slope'),
    ('o2','slope_axial','current_actual_four_O2_tensor_interface_rows','slope_axial'),
    ('o2','axial_buffer','current_actual_four_O2_tensor_interface_rows','axial_buffer'),
    ('o2','buffer_transition','current_actual_four_O2_tensor_interface_rows','buffer_transition'),
    ('o3','transition_power','current_actual_O3_transition_power_interface_rows','transition_power'),
    ('incoming','incoming_entrance','current_actual_two_entrance_incoming_tensor_interface_rows','incoming_entrance'),
    ('incoming','entrance_main','current_actual_two_entrance_incoming_tensor_interface_rows','entrance_main'),
    ('main','main_exit','current_actual_two_main_exit_tensor_interface_rows','main_exit'),
    ('main','exit_gap','current_actual_two_main_exit_tensor_interface_rows','exit_gap'),
    ('gap','gap_coordinate','current_actual_two_gap_tensor_interface_rows','gap_coordinate'),
    ('gap','gap_end','current_actual_two_gap_tensor_interface_rows','gap_end'),
    ('end','end_flatten','current_actual_end_flatten_tensor_interface_rows',None),
    ('flatten','flatten_power','current_actual_two_flatten_power_tensor_interface_rows','flatten_power'),
    ('flatten','power_angular','current_actual_two_flatten_power_tensor_interface_rows','power_angular'),
    ('entry','angular_entry','actual_angular_entry_common_tensor_rows_checked',None),
    ('o7','steep_entry_power','current_actual_O7_three_common_tensor_interface_rows','steep_entry_power'),
    ('o7','steep_power_exit','current_actual_O7_three_common_tensor_interface_rows','steep_power_exit'),
    ('o7','steep_exit_waiting','current_actual_O7_three_common_tensor_interface_rows','steep_exit_waiting'),
    ('heat','waiting_collar','current_actual_two_heat_tensor_interface_rows','waiting_collar'),
    ('heat','collar_exterior','current_actual_two_heat_tensor_interface_rows','collar_exterior'))

def defining_method(owner,method):
    fn=getattr(owner,method).__func__;asts=SourceAST()
    source=asts.method(owner.__class__.__module__.removeprefix(PREFIX),fn.__name__)
    return dict(provider=owner.__class__.__module__+'.'+owner.__class__.__name__,method=method,
        original_signature=str(inspect.signature(getattr(owner,method))),
        defining_method_AST_sha256=hashlib.sha256(ast.dump(source,include_attributes=False).encode()).hexdigest(),
        defining_method_source=ast.unparse(source),input_hashes=asts.hashes)

def cartesian_core_groups(view):
    groups={}
    for name in ('physical_core_Cartesian_stress_mixed3','physical_core_Cartesian_stress_divergence_mixed2',
                 'physical_core_Cartesian_remainder_mixed2','physical_core_Cartesian_momentum_decomposition_mixed2'):
        for index,components in view[name].items():
            for component,row in components.items():
                groups[name+'/'+index+'/'+component]=list(row.values()) if 'signed_coefficient' not in row else [row]
    return groups

class CurrentTensorRegistry:
    @source_precision
    def __init__(self,background=None,require_checked=True):
        self.background=background if background is not None else CurrentAngularInternalBackgroundTensor()
        self.physical=self.background.physical;self.ctx=self.background.ctx
        self.family=self.background.family;self.source=self.background.source;self.datum_sha=self.background.datum_sha
        self.owners={};obj=self.background.core_axis
        for alias,attribute in zip(ALIASES,CHAIN):
            obj=getattr(obj,attribute);self.owners[alias]=obj
        self.owners['axis']=self.background.core_axis;self.owners['angular_internal']=self.background
        self.owner_ids={name:id(owner) for name,owner in self.owners.items()}
        self.routes={row[0]:row for row in ROUTES};self.receipts={};self.bindings={};self.methods={}
        self.hashes=dict(self.background.hashes);self.assert_graph()
        for alias,owner in self.owners.items():
            receipt=owner.__class__.__module__+'_check.json'
            if receipt not in self.hashes or self.hashes[receipt]!=sha(receipt):raise ValueError('Unbound regional receipt: '+receipt)
            record=json.loads((HERE/receipt).read_bytes())
            if not record['all_passed'] or (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'],record['datum_enclosure_sha256'])!=(self.family,self.source,self.datum_sha):
                raise ValueError('Regional source receipt differs: '+receipt)
            self.receipts[alias]=(receipt,record)
        for region,alias,method,chart,domain in ROUTES:
            self._bind(alias,method)
        for alias,method in (('axis','axis'),('axis','cartesian'),('heat','unbounded_exterior'),('end','support_interface'),('angular_internal','interface')):
            self._bind(alias,method)
        for alias,_,_,_ in ADJACENT:self._bind(alias,'interface')
        self.registry=self._region_records();self.adjacent=self._adjacent_records();self.internal=self._internal_records()
        regions=self.background.manifest()['actual_current_tensor_regions_available']
        if list(self.registry)!=regions or len(self.adjacent)!=32 or len(self.internal)!=14:
            raise ValueError('Exact same-source 33/32/14 inventory required')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Native tensor registry receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def _bind(self,alias,method):
        key=alias+'.'+method
        if key not in self.bindings:
            owner=self.owners[alias];self.bindings[key]=defining_method(owner,method)
            self.methods[key]=getattr(owner,method).__func__;self.hashes.update(self.bindings[key]['input_hashes'])

    def assert_graph(self):
        self.background.assert_graph();self.background.assert_programs()
        if self.routes!={row[0]:row for row in ROUTES}:raise ValueError('Original native route table changed')
        if type(self.background) is not CurrentAngularInternalBackgroundTensor or not self.background.acceptance_loaded:
            raise ValueError('Checked full core-axis/angular tensor root required')
        obj=self.background.core_axis
        for alias,attribute in zip(ALIASES,CHAIN):
            obj=getattr(obj,attribute)
            if self.owners[alias] is not obj:raise ValueError('Foreign tensor owner path: '+alias)
        if self.owners['axis'] is not self.background.core_axis or self.owners['angular_internal'] is not self.background:
            raise ValueError('Foreign axis/internal source owner')
        for alias,owner in self.owners.items():
            physical=owner.core_tensor.physical if alias=='axis' else owner.physical
            if id(owner)!=self.owner_ids[alias] or not owner.acceptance_loaded or physical is not self.physical or owner.ctx is not self.ctx or (owner.family,owner.source,owner.datum_sha)!=(self.family,self.source,self.datum_sha):
                raise ValueError('One admitted physical/core/history/pressure graph required')
        for key,fn in self.methods.items():
            alias,method=key.split('.')
            if getattr(self.owners[alias],method).__func__ is not fn:raise ValueError('Foreign native operator: '+key)

    def _evidence(self,alias,key,child=None):
        name,record=self.receipts[alias];value=record[key] if child is None else record[key][child]
        if not value:raise ValueError('Missing actual tensor trace admission')
        return dict(receipt=name,receipt_sha256=self.hashes[name],receipt_field=key,
            receipt_child=child,actual_admitted_trace_evidence=value)

    def _region_records(self):
        result={}
        for (region,alias,method,chart,domain),coordinate in zip(ROUTES,SOURCE_COORDINATES):
            result[region]=dict(owner_alias=alias,owner_path='core_axis.'+'.'.join(CHAIN[:ALIASES.index(alias)+1]),
                method=method,native_chart_argument=chart,native_domain=domain,source_Z_domain=[-1,1],
                source_coordinate=coordinate,derivative_pullback=derivative_pullback(region,alias,coordinate),
                defining_method_binding=alias+'.'+method,acceptance_receipt=self.receipts[alias][0],
                receipt_sha256=self.hashes[self.receipts[alias][0]])
        return result

    def _adjacent_records(self):
        regions=list(self.routes);result={}
        for j,(alias,name,key,child) in enumerate(ADJACENT):
            result[name]=dict(owner_alias=alias,method='interface',named_argument=alias not in ('end','entry'),
                left_region=regions[j],right_region=regions[j+1],defining_endpoint_recipe=alias+'.interface',
                **self._evidence(alias,key,child))
        return result

    def _internal_records(self):
        result={}
        for edge in (49,51,59,61,69,71):
            name='patch_support_'+str(edge);result[name]=dict(owner_alias='patch',method='interface',region='actual_patch',
                source_argument=name,defining_endpoint_recipe='patch.interface',
                **self._evidence('patch','current_actual_patch_Rh_and_six_support_tensor_interface_rows',name))
        for edge in PULSE_EDGES:
            name='pulse_end_'+edge['exact_edge'];result[name]=dict(owner_alias='end',method='support_interface',region='pulse_end',
                exact_source_edge=edge['exact_edge'],source_argument=edge,defining_endpoint_recipe='end.support_interface',
                **self._evidence('end','current_actual_four_pulse_end_support_tensor_interface_rows',edge['exact_edge']))
        for name,q in ANGULAR_SEAMS.items():
            result[name]=dict(owner_alias='angular_internal',method='interface',region='outer_angular',
                exact_source_edge=q,source_argument=name,defining_endpoint_recipe='angular_internal.interface',
                **self._evidence('angular_internal','current_four_angular_internal_full_tensor_view_counts',name+'_whole'))
        return result

    def _view(self,region,view,layout):
        groups=cartesian_core_groups(view) if layout=='Cartesian_core' else (
            canonical_angular_groups(self.ctx,view) if layout=='postpulse' else canonical_tensor_groups(view))
        return dict(region=region,layout=layout,original_complete_view=view,canonical_signed_component_groups=groups,
            original_full_source_sectors_and_factors_retained=True,return_kind='Directed source-function enclosure, not a resolved physical point value',
            actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            **dict.fromkeys(OPEN,False))

    @source_precision
    def native(self,region,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph()
        if region not in self.routes:raise ValueError('Unknown full tensor region')
        _,alias,method,chart,_=self.routes[region];owner=self.owners[alias]
        if region=='core_positive_radius' and endpoints(self.ctx.mpf(coordinate))==endpoints(self.ctx.mpf(0)):
            if theta is not None and not all(mp.isfinite(v) for v in endpoints(self.ctx.mpf(theta))):raise ValueError('Finite angle required')
            return self._view(region,self.owners['axis'].axis(Z,log_tau,viscosity),'Cartesian_core')
        fn=getattr(owner,method)
        view=fn(chart,Z,coordinate,log_tau,theta,viscosity) if method=='chart' else fn(Z,coordinate,log_tau,theta,viscosity)
        return self._view(region,view,'postpulse' if alias in ('flatten','angular','entry','o7','heat') else 'cylindrical')

    @source_precision
    def core_cartesian(self,Z,X=0,Y=0,log_tau='-1',viscosity='1'):
        self.assert_graph();return self._view('core_positive_radius',self.owners['axis'].cartesian(Z,X,Y,log_tau,viscosity),'Cartesian_core')

    @source_precision
    def trace(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph();args=(Z,log_tau,theta,viscosity)
        if name in self.adjacent:
            route=self.adjacent[name];fn=self.owners[route['owner_alias']].interface
            value=fn(name,*args) if route['named_argument'] else fn(*args)
        elif name in self.internal:
            route=self.internal[name];value=getattr(self.owners[route['owner_alias']],route['method'])(route['source_argument'],*args)
        else:raise ValueError('Unknown admitted full tensor trace')
        key='common_actual_angular_entry_tensor_rows' if route['owner_alias']=='entry' else 'common_actual_tensor_rows'
        return dict(seam=name,common_actual_tensor_rows=value[key],original_complete_trace=value,
            defining_endpoint_recipe=route['defining_endpoint_recipe'],receipt=route['receipt'],receipt_sha256=route['receipt_sha256'],
            original_source_identity_proof_and_complete_common_bounds_retained=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def exterior(self,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph();return self._view('heat_exterior',self.owners['heat'].unbounded_exterior(Z,log_tau,theta,viscosity),'postpulse')

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_full_tensor_native_region_registry=self.registry,current_full_tensor_adjacent_trace_registry=self.adjacent,
            current_full_tensor_internal_trace_registry=self.internal,current_original_callable_source_bindings=self.bindings,
            actual_current_tensor_regions_available=list(self.registry),actual_current_completed_tensor_adjacent_interface_count=32,
            actual_current_completed_tensor_internal_interface_count=14,internal_count_composition=dict(actual_patch=6,pulse_end=4,outer_angular=4),
            core_axis_is_separate_analytic_extension_not_a_counted_internal_support=True,
            current_core_axis_dispatch=dict(region='core_positive_radius',provider='axis',method='cartesian/axis',
                source_coordinates='(Z,X,Y), rho=(X^2+Y^2)/2; whole rho[0,4] disk',
                derivative_pullback='Original physical Cartesian mixed3 stress and mixed2 divergence/remainder; exact epsilon/F0/nu/lambda factors, no inverse radius or angle at axis',
                defining_method_binding='axis.cartesian'),
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            current_leading_origin_remainder_time_obstruction=self.background.leading_obstruction,
            scope='One checked graph; original 33 native full tensor APIs and 32 adjacent/14 support trace APIs with actual receipt rows. Explicit regional coordinates; compact requested positive times. Physical-coordinate locator/global cover, cone/lift, resolved points, completed flatness and coefficient recovery remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentTensorRegistry(require_checked=False)
    field.assert_graph();result=field.manifest()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current full tensor registry built: 33 native regions, 32 adjacent and 14 internal support routes; same checked graph',flush=True)
    return result

if __name__=='__main__':run()
