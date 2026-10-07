"""Six registry-named physical patch interfaces without a graph rebuild.

Queries return bounds on the admitted source functions. A physical query
uses actual lambda and constant viscosity; it does not resolve primitive
coefficients. Other global velocity interfaces still require providers.
"""
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_physical_numeric_bounds as numeric
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    CompliantGlobalPhysicalAssembly as BASE,cartesian_source_row,time_source_row,log_row,BASES)
from lei_ren_part1_paper_compliant_current_modified_physical_velocity_operator import (
    actual_physical_velocity_bounds,viscosity_power)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces_operator import physical_contribution_groups

t = numeric.transport
HERE,PREFIX,sha = numeric.HERE,numeric.PREFIX,numeric.sha
NAME = PREFIX+'current_patch_velocity_interface_atlas.json'
RECEIPT = PREFIX+'current_patch_velocity_interface_atlas_check.json'
VIEWS_NAME = PREFIX+'current_patch_velocity_interface_atlas_views.json.gz'
SEAMS = {'patch_support_'+str(edge):edge for edge in t.EDGES}


def defining_germs(source,edge):
    record = source['accepted']['records'][t.SOURCE_NAME]['exact_actual_patch_support_velocity_pressure_theorem']['exact_support_source_theorems'][str(edge)]
    return {side:dict(named_source_side=side,exact_x=record['exact_x'],
        active_bump=record['active_bump'],entering=record['entering'],
        inside_strip_source_symbols=record['inside_strip_source_symbols_by_side'][i],
        last_edge_uses_current_five_functional_implicit_equations=record['last_edge_uses_current_five_functional_implicit_equations'],
        source_receipt=t.SOURCE_RECEIPT,physical_source_receipt=t.RECEIPT)
        for i,side in enumerate(('left','right'))}


def request_arguments(c,Z,log_tau,theta,viscosity,physical=False):
    z,lt,nu = c.mpf(Z),c.mpf(log_tau),c.mpf(viscosity)
    for value in (z,lt,nu):
        if any(not mp.isfinite(v) for v in t.endpoints(value)): raise ValueError('Finite Z, log_tau and constant viscosity required')
    lo,hi = t.endpoints(z)
    if lo < -1 or hi > 1 or (physical and (lo<=-1 or hi>=1)):
        raise ValueError('Physical points require |Z|<1; source closure requires |Z|<=1')
    if t.endpoints(nu)[0] <= 0: raise ValueError('One positive constant viscosity required')
    if theta is None: cosine=sine=c.mpf([-1,1])
    else:
        angle=c.mpf(theta)
        if any(not mp.isfinite(v) for v in t.endpoints(angle)): raise ValueError('Finite angle required')
        cosine,sine=c.cos(angle),c.sin(angle)
    return z,lt,nu,cosine,sine


def requested_native_view(source,edge,Z,lt,cosine,sine):
    c=source['ctx'];packet=source['primitive']['source_bump_edge_packets'][t.EDGES.index(edge)]
    x=c.mpf(edge)/40;logr=source['logRref']-6+c.ln(x)
    raw={label:{key:t.read_interval(c,value) for key,value in rows.items()}
        for label,rows in packet['physical_velocity_pressure_y_Z_mixed4'].items()}
    grids,logs,amplitudes=BASE.normalized_sources(source['field'],'actual_patch',Z,x,
        {'physical_velocity_pressure_y_Z_mixed4':raw},None,logr)
    spatial={};time={}
    for i,j,b in t.INDICES:
        key='x%d_y%d_z%d'%(i,j,b);spatial[key]={}
        for component in t.COMPONENTS:
            parts=cartesian_source_row(c,grids,component,i,j,b,Z,source['delta'],cosine,sine,amplitudes)
            spatial[key][component]={label:log_row(c,rows,logs,gamma,lt/2) for label,(rows,gamma) in parts.items()}
    basis={'ux':{t.UR:cosine,t.UT:-sine},'uy':{t.UR:sine,t.UT:cosine},'uz':{t.UZ:c.mpf(1)},'p':{t.P:c.mpf(1)}}
    for component,labels in basis.items():
        time[component]={}
        for label,angular in labels.items():
            rows,gamma=time_source_row(c,grids,label,Z,source['delta'],amplitudes)
            time[component][label]=log_row(c,[(powers,value*angular) for powers,value in rows],logs,gamma,lt/2)
    return dict(chart='actual_patch',exact_support_x=str(s.Rational(edge,40)),Z=Z,requested_log_tau=lt,
        cos_theta=cosine,sin_theta=sine,physical_layout='four_label_Cartesian',
        physical_spatial_cartesian_mixed4=spatial,first_fixed_x_physical_time_derivative=time,
        source_logR_enclosure=logr,positive_source_log_bases=dict(zip(BASES,logs)),
        primitive_coefficients_enclose_whole_Z_not_resolved_at_requested_Z=True,
        normalized_source_rows_not_refitted_or_interpolated=True,
        physical_log_bounds_use_lambda_ge_sqrt_tau=True)


def constant_viscosity_view(c,view,nu):
    result=dict(view);lnnu=c.ln(nu);spatial={};time={}
    def scale(row,component,order):
        power=c.mpf(str(viscosity_power(component,order)));out=dict(row)
        out.update(constant_viscosity_power=power,
            constant_viscosity_factor_held_separately_from_native_signed_terms=True)
        if not row['exact_zero']: out['log_absolute_upper']=row['log_absolute_upper']+power*lnnu
        return out
    for i,j,b in t.INDICES:
        key='x%d_y%d_z%d'%(i,j,b)
        spatial[key]={component:{label:scale(row,component,i+j+b) for label,row in rows.items()}
            for component,rows in view['physical_spatial_cartesian_mixed4'][key].items()}
    for component,rows in view['first_fixed_x_physical_time_derivative'].items():
        time[component]={label:scale(row,component,0) for label,row in rows.items()}
    result.update(physical_spatial_cartesian_mixed4=spatial,first_fixed_x_physical_time_derivative=time,
        physical_constant_viscosity=nu,viscosity_constant_in_space_and_time=True)
    return result


class CurrentPatchVelocityInterfaceAtlas:
    def __init__(self,require_checked=True):
        self.source=numeric.inputs();self.ctx=self.source['ctx']
        record=json.loads((HERE/numeric.RECEIPT).read_bytes())
        if not record['all_passed'] or not record['current_six_patch_support_physical_numeric_trace_bounds_available']:
            raise ValueError('Fresh checked primitive-source physical bounds required')
        for path,digest in record['input_hashes'].items():
            if sha(path)!=digest: raise ValueError('Changed fresh physical input: '+path)
        if record['fresh_patch_physical_bound_definition_sha256']!=numeric.definition(self.source):
            raise ValueError('Foreign primitive/scalar source bound definition')
        self.hashes={**record['input_hashes'],numeric.RECEIPT:sha(numeric.RECEIPT),Path(__file__).name:sha(Path(__file__).name)}
        self.theorem=exact_atlas_theorem(self)
        self.hashes.update(self.theorem['input_hashes'])
        identity=dict(parent_fresh_physical_bound=numeric.definition(self.source),
            adapter_sha256=sha(Path(__file__).name),physical_operator_sha256=sha(PREFIX+'current_modified_physical_velocity_operator.py'))
        self.definition=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            for path,digest in checked['input_hashes'].items():
                if sha(path)!=digest: raise ValueError('Changed accepted patch interface atlas: '+path)
            if not checked['all_passed'] or checked['patch_velocity_interface_atlas_definition_sha256']!=self.definition:
                raise ValueError('Checked current patch interface atlas required')
            self.acceptance_loaded=True

    def _edge(self,name):
        if name not in SEAMS: raise ValueError('One of six registry patch_support names required')
        return SEAMS[name]

    def interface(self,name,Z=(-1,1),log_tau=(-3,-1),theta=None,viscosity=1):
        edge=self._edge(name);z,lt,nu,cosine,sine=request_arguments(self.ctx,Z,log_tau,theta,viscosity)
        native=requested_native_view(self.source,edge,z,lt,cosine,sine)
        view=constant_viscosity_view(self.ctx,native,nu)
        germs=defining_germs(self.source,edge);common=physical_contribution_groups(view)
        return dict(seam=name,Z=z,requested_log_tau=lt,physical_constant_viscosity=nu,
            left_defining_source_germ=germs['left'],right_defining_source_germ=germs['right'],
            left_complete_velocity_pressure_source=dict(view,named_source_side='left'),
            right_complete_velocity_pressure_source=dict(view,named_source_side='right'),
            common_physical_velocity_pressure_rows=common,current_common_physical_contribution_count=len(common),
            patch_velocity_interface_atlas_definition_sha256=self.definition,
            source_function_and_FTC_identities_precede_physical_bounds=True,
            interval_overlap_or_midpoint_not_used_as_function_identity=True,
            all_finite_log_tau_sectors_and_positive_constant_viscosities_supported=True,
            resolved_physical_point_values_available=False,all_global_velocity_interfaces_certified=False)

    def physical_interface(self,name,Z,log_tau,theta=None,viscosity=1):
        edge=self._edge(name);z,lt,nu,cosine,sine=request_arguments(self.ctx,Z,log_tau,theta,viscosity,physical=True)
        native=requested_native_view(self.source,edge,z,lt,cosine,sine)
        q=(lt-self.ctx.ln(1-z*z))/2
        location=dict(actual_log_lambda=q,requested_log_tau=lt,physical_viscosity=nu)
        bounds=actual_physical_velocity_bounds(SimpleNamespace(ctx=self.ctx),native,location)
        coordinates=dict(actual_log_lambda=q,
            physical_log_r=self.ctx.ln(nu)/2+q+(native['source_logR_enclosure']+self.ctx.ln(2))/2,
            physical_z_source_coefficient=z,
            physical_z_log_scale=self.ctx.ln(nu)/2+(1-self.source['delta'])*q,
            implicit_relation='lambda^2-(z_phys^2/nu)*lambda^(2delta)=tau',
            radius_relation='r_phys=sqrt(nu)*lambda*sqrt(2*Rm*x)')
        return dict(seam=name,left_and_right_defining_source_germs=defining_germs(self.source,edge),
            canonical_native_velocity_pressure_trace=native,physical_coordinates=coordinates,
            actual_physical_velocity_pressure_bounds=bounds,
            patch_velocity_interface_atlas_definition_sha256=self.definition,
            same_source_function_bounds_for_both_physical_traces=True,
            resolved_physical_point_values_available=False,all_global_velocity_interfaces_certified=False)

    def manifest(self):
        return dict(source_family=self.source['accepted']['source_family'],
            patch_velocity_interface_atlas_definition_sha256=self.definition,
            parent_fresh_physical_bound_definition_sha256=numeric.definition(self.source),
            registry_named_internal_patch_interfaces=list(SEAMS),exact_source_transport_theorem=self.theorem,
            complete_views=VIEWS_NAME,finite_time_sector_restriction_removed_by_fresh_source_mapping=True,
            positive_constant_viscosity_and_actual_lambda_query_bridge_available=True,
            current_patch_physical_velocity_interface_provider_available=True,
            resolved_point_values_global_interfaces_energy_common_N_cones_recursion_NS_remain_open=True,
            input_hashes=self.hashes)


def exact_atlas_theorem(atlas):
    asts=t.SourceAST();checks={}
    asts.expression('current_tensor_registry','_internal_records','name',wanted="'patch_support_'+str(edge)")
    for method in ('cartesian_source_row','time_source_row','normalized_sources'): asts.method('global_physical_assembly',method)
    for method in ('actual_physical_velocity_bounds','viscosity_power'): asts.method('current_modified_physical_velocity_operator',method)
    asts.expression('current_patch_velocity_interface_atlas','physical_interface','q',
        wanted='(lt-self.ctx.ln(1-z*z))/2')
    asts.expression('current_patch_velocity_interface_atlas','physical_interface','bounds',
        wanted='actual_physical_velocity_bounds(SimpleNamespace(ctx=self.ctx),native,location)')
    asts.expression('current_patch_velocity_interface_atlas','requested_native_view','logr',
        wanted="source['logRref']-6+c.ln(x)")
    z,delta,lt,lnnu,lr=s.symbols('Z delta log_tau log_nu log_R',real=True)
    q=(lt-s.log(1-z*z))/2
    lam,nu=s.symbols('lambda nu',positive=True)
    physical_z=s.sqrt(nu)*z*lam**(1-delta)
    identities=dict(actual_implicit_lambda_relation=s.exp(2*q)*(1-z*z)-s.exp(lt),
        actual_log_radius_relation=(lnnu/2+q+(lr+s.log(2))/2)-(lnnu+s.log(2)+lr+2*q)/2,
        actual_signed_axial_implicit_relation=lam**2-(physical_z**2/nu)*lam**(2*delta)-lam**2*(1-z*z))
    for key,value in identities.items():
        if s.simplify(value)!=0: raise ArithmeticError('Physical patch coordinate identity: '+key)
        checks[key]=True
    for component in t.COMPONENTS:
        for order in range(5):
            if viscosity_power(component,order)!=s.Rational((2 if component=='p' else 1)-order,2):
                raise ArithmeticError('Constant viscosity source power differs')
            checks['same_actual_nu_power_'+component+'_'+str(order)]=True
    hashes={**atlas.source['hashes'],**asts.hashes}
    return dict(passed=True,identities=checks,source_bindings=asts.bindings,input_hashes=hashes,
        consumed_six_independent_source_and_physical_join_receipts=[t.SOURCE_RECEIPT,t.RECEIPT,numeric.RECEIPT],
        only_six_registry_named_internal_velocity_interfaces_installed=True,
        underlying_primitive_grid_uncertainty_not_resolved_or_interpolated=True)


def run():
    atlas=CurrentPatchVelocityInterfaceAtlas(require_checked=False);views={}
    for name in SEAMS:
        views[name]=atlas.interface(name)
        print('Registry-named actual patch velocity interface:',name,flush=True)
    views['fresh_time_nu_sector']=atlas.interface('patch_support_59',Z=('-0.7','0.8'),log_tau=('-12','-11'),viscosity='.8')
    views['actual_physical_point']=atlas.physical_interface('patch_support_71','.237','-2.337','.337','1.3')
    views['actual_physical_box']=atlas.physical_interface('patch_support_49',('-.999','-.997'),('-9','-8'),viscosity='.8')
    payload=(json.dumps(numeric.encode(numeric.pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result=atlas.manifest();result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_bytes((json.dumps(numeric.encode(numeric.pack(result)),indent=2)+'\n').encode())
    return result


if __name__=='__main__': run()
