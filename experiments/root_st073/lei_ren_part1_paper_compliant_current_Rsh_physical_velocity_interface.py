"""Reshape/reference physical trace adapter from both accepted source germs.

The original flat boundary theorem is consumed before merging bounds.
Frozen basepoint normalization is transported without re-differentiation.
Neither source coefficients nor the finite left neighborhood are replaced.
"""
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_velocity_interface_atlas as patch
import lei_ren_part1_paper_compliant_actual_Rsh_source_join as join
from lei_ren_part1_paper_compliant_current_interface_atlas import merge_trace_rows
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import original_native_radius_recipes,SYMBOLS,V
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value

numeric,t = patch.numeric,patch.t
HERE,PREFIX,sha = patch.HERE,patch.PREFIX,patch.sha
NAME=PREFIX+'current_Rsh_physical_velocity_interface.json'
RECEIPT=PREFIX+'current_Rsh_physical_velocity_interface_check.json'
VIEWS_NAME=PREFIX+'current_Rsh_physical_velocity_interface_views.json.gz'
JOIN_NAME=PREFIX+'actual_Rsh_source_join.json'
JOIN_RECEIPT=PREFIX+'actual_Rsh_source_join_check.json'
PACKETS={'left':PREFIX+'actual_long_reshape_mixed_C4.json',
         'right':PREFIX+'actual_reference_restore_mixed_C4.json'}
LABELS=('Utheta_over_current_Utheta','Uz','Ur_over_current_sqrt_R_over_2','P_over_Pstar2')
KEYS={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}


def inputs():
    source=numeric.inputs();c=source['ctx'];hashes=dict(source['hashes'])
    records={name:json.loads((HERE/name).read_bytes()) for name in (JOIN_NAME,JOIN_RECEIPT,*PACKETS.values())}
    receipt=records[JOIN_RECEIPT];proof=records[JOIN_NAME]['boundary_source_proof']
    if not receipt['all_passed'] or not receipt['current_Rsh_source_functional_join_certified'] or not proof['passed']:
        raise ValueError('Checked current two-sided Rsh source theorem required')
    if proof['finite_left_neighborhood_replaced'] or not proof['normalization_is_frozen_after_physical_differentiation']:
        raise ValueError('Original finite source and frozen post-derivative units required')
    for name,record in records.items():
        for key,value in source['accepted']['source_family'].items():
            if record.get(key)!=value: raise ValueError('Foreign Rsh source or analytic P0: '+name+' '+key)
        for path,digest in record['input_hashes'].items():
            if path in hashes and hashes[path]!=digest:raise ValueError('Conflicting Rsh/current source: '+path)
            if sha(path)!=digest:raise ValueError('Changed Rsh source dependency: '+path)
            hashes[path]=digest
        hashes[name]=sha(name)
    if join.source_bindings()!=records[JOIN_NAME]['current_Rsh_source_bindings']:
        raise ValueError('Original Rsh defining source call sites changed')
    packets={side:restore_value(c,records[name]['actual_Rsh_exit']) for side,name in PACKETS.items()}
    for side,packet in packets.items():
        if t.endpoints(packet['Z'])!=(-1,1) or set(packet['physical_velocity_pressure_y_Z_mixed4'])!=set(LABELS):
            raise ValueError('Whole-Z four-label primitive source rows required')
        if any(set(row)!=KEYS for row in packet['physical_velocity_pressure_y_Z_mixed4'].values()):
            raise ValueError('Complete ordinary mixed4 source derivative grid required')
    if packets['left']['phase']._mpi_!=c.mpf(1)._mpi_ or packets['right']['chart']!='Rsh_to_Rz':
        raise ValueError('Actual source boundary endpoints required')
    T=packets['left']['original_T']
    if T._mpi_!=packets['left']['actual_inherited_axial5_packet']['log_R_over_110_enclosure']._mpi_:
        raise ValueError('Original selected T differs from the defining Rsh parent')
    logR=c.ln(110)+T
    return dict(source=source,ctx=c,packets=packets,records=records,logR=logR,original_T=T,hashes=hashes)


def native_trace(data,side,Z,lt,cosine,sine):
    c=data['ctx'];packet=data['packets'][side];chart='reshape' if side=='left' else 'inner_reference'
    grids,logs,amplitudes=patch.BASE.normalized_sources(data['source']['field'],chart,Z,
        c.mpf(1 if side=='left' else 0),packet,None,data['logR'])
    spatial={};time={}
    for i,j,b in t.INDICES:
        key='x%d_y%d_z%d'%(i,j,b);spatial[key]={}
        for component in t.COMPONENTS:
            parts=patch.cartesian_source_row(c,grids,component,i,j,b,Z,data['source']['delta'],cosine,sine,amplitudes)
            spatial[key][component]={label:patch.log_row(c,rows,logs,gamma,lt/2,False) for label,(rows,gamma) in parts.items()}
    basis={'ux':{t.UR:cosine,t.UT:-sine},'uy':{t.UR:sine,t.UT:cosine},'uz':{t.UZ:c.mpf(1)},'p':{t.P:c.mpf(1)}}
    for component,labels in basis.items():
        time[component]={}
        for label,angle in labels.items():
            rows,gamma=patch.time_source_row(c,grids,label,Z,data['source']['delta'],amplitudes)
            time[component][label]=patch.log_row(c,[(powers,value*angle) for powers,value in rows],logs,gamma,lt/2,False)
    return dict(chart=chart,named_source_side=side,Z=Z,requested_log_tau=lt,cos_theta=cosine,sin_theta=sine,
        source_logR_enclosure=data['logR'],positive_source_log_bases=dict(zip(patch.BASES,logs)),
        physical_spatial_cartesian_mixed4=spatial,first_fixed_x_physical_time_derivative=time,
        physical_layout='four_label_Cartesian',primitive_source=PACKETS[side]+'#/actual_Rsh_exit',
        whole_Z_frozen_primitive_normalization_retained=True,
        already_differentiated_physical_rows_not_differentiated_again=True)


def exact_theorem(data):
    asts=t.SourceAST();checks={}
    for method in ('normalized_sources','cartesian_source_row','time_source_row'):
        asts.method('global_physical_assembly',method)
    asts.expression('current_Rsh_physical_velocity_interface','native_trace','(grids, logs, amplitudes)',
        wanted="patch.BASE.normalized_sources(data['source']['field'],chart,Z,c.mpf(1 if side=='left' else 0),packet,None,data['logR'])")
    recipes,radius_proof=original_native_radius_recipes();q=SYMBOLS
    left=recipes['reshape'].subs(V,1);right=recipes['inner_reference'].subs(V,0)
    if s.expand(left-right)!=0 or s.expand(left-q['T']-s.log(110))!=0:
        raise ArithmeticError('Original two-sided Rsh radius source differs')
    checks['two_original_native_radius_recipes_at_Rsh_identical']=True
    # Canonical Cartesian/time maps are linear in the already physical
    # source rows. This audits their source order and exact fixed units;
    # the admitted Rsh function theorem supplies two-sided equality.
    c=t.exact_context();z,delta,logP,logU,logR=s.symbols('Z delta logP logU logR',real=True)
    atoms={label:{key:s.Symbol(label+'_'+key) for key in KEYS} for label in LABELS}
    packet=dict(physical_velocity_pressure_y_Z_mixed4=atoms,
        actual_inherited_axial5_packet=dict(log_Utheta_over_Pstar_axial5_coefficients=[logU-logP]))
    field=SimpleNamespace(ctx=c,logP=logP)
    grids,logs,amp=patch.BASE.normalized_sources(field,'reshape',z,1,packet,None,logR)
    if logs[4]!=logU or amp[t.UT]!={4:1} or amp[t.UR]!={5:s.Rational(1,2),6:-s.Rational(1,2)} or amp[t.P]!={1:2}:
        raise ArithmeticError('Frozen Rsh physical normalization changed')
    checks['original_frozen_swirl_radial_pressure_units_retained']=True
    source_atoms={v for row in atoms.values() for v in row.values()}
    def linear_source_value(value):
        active=value.free_symbols & source_atoms
        if active and s.Poly(value,*active).total_degree()>1:
            raise ArithmeticError('Physical Rsh operator is not linear in source derivatives')
    count=0
    for i,j,b in t.INDICES:
        for component in t.COMPONENTS:
            parts=patch.cartesian_source_row(c,grids,component,i,j,b,z,delta,t.CS,t.SN,amp)
            for label,(rows,gamma) in parts.items():
                for powers,value in rows:
                    if not value.free_symbols.issubset({z,delta,t.CS,t.SN,*[v for row in atoms.values() for v in row.values()]}):
                        raise ArithmeticError('Physical Rsh map requests unknown source derivatives')
                    linear_source_value(value)
                checks['linear_spatial_map_'+component+'_'+label+'_%d_%d_%d'%(i,j,b)]=True;count+=1
    for label in (t.UT,t.UZ,t.UR,t.P):
        rows,gamma=patch.time_source_row(c,grids,label,z,delta,amp)
        if any(not value.free_symbols.issubset({z,delta,*[v for row in atoms.values() for v in row.values()]}) for powers,value in rows):
            raise ArithmeticError('Fixed-position Rsh time source derivatives changed')
        for powers,value in rows:linear_source_value(value)
        checks['fixed_position_time_map_'+label]=True
    if count!=210:raise ArithmeticError('Complete Rsh spatial source contributions required')
    checks['accepted_source_function_identity_precedes_common_enclosure']=True
    return dict(passed=True,identities=checks,source_bindings=asts.bindings,
        consumed_Rsh_function_theorem=JOIN_RECEIPT,original_source_function_mixed4_rows=135,
        physical_spatial_contribution_count=count,physical_time_contribution_count=6,
        same_frozen_normalization_not_reapplied_as_product_derivative=True,
        input_hashes={**asts.hashes,**radius_proof['input_hashes'],Path(__file__).name:sha(Path(__file__).name)})


class CurrentRshPhysicalVelocityInterface:
    def __init__(self,require_checked=True):
        self.data=inputs();self.ctx=self.data['ctx'];self.theorem=exact_theorem(self.data)
        self.hashes={**self.data['hashes'],**self.theorem['input_hashes']}
        self.definition=hashlib.sha256(json.dumps(self.hashes,sort_keys=True).encode()).hexdigest()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or receipt['Rsh_physical_velocity_interface_definition_sha256']!=self.definition:
                raise ValueError('Checked current Rsh physical interface required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed accepted Rsh interface: '+name)

    def interface(self,name='reshape_reference',Z=(-1,1),log_tau=(-3,-1),theta=None,viscosity=1):
        if name!='reshape_reference':raise ValueError('Only reshape_reference is provided')
        z,lt,nu,cs,sn=patch.request_arguments(self.ctx,Z,log_tau,theta,viscosity)
        traces={side:patch.constant_viscosity_view(self.ctx,native_trace(self.data,side,z,lt,cs,sn),nu) for side in PACKETS}
        groups={side:patch.physical_contribution_groups(view) for side,view in traces.items()}
        common=merge_trace_rows(self.ctx,list(groups.values()))
        return dict(seam=name,source_family=self.data['source']['accepted']['source_family'],
            Z=z,requested_log_tau=lt,physical_constant_viscosity=nu,
            left_complete_velocity_pressure_source=traces['left'],right_complete_velocity_pressure_source=traces['right'],
            left_physical_contribution_rows=groups['left'],right_physical_contribution_rows=groups['right'],
            common_physical_velocity_pressure_rows=common,current_common_physical_contribution_count=len(common),
            Rsh_physical_velocity_interface_definition_sha256=self.definition,
            independent_source_germs_identified_by_current_function_theorem=True,
            matching_enclosures_not_used_to_prove_source_equality=True,
            all_global_velocity_interfaces_certified=False,resolved_physical_point_values_available=False)

    def physical_interface(self,Z,log_tau,theta=None,viscosity=1):
        z,lt,nu,cs,sn=patch.request_arguments(self.ctx,Z,log_tau,theta,viscosity,physical=True)
        q=(lt-self.ctx.ln(1-z*z))/2;location=dict(actual_log_lambda=q,requested_log_tau=lt,physical_viscosity=nu)
        traces={side:native_trace(self.data,side,z,lt,cs,sn) for side in PACKETS}
        return dict(seam='reshape_reference',actual_log_lambda=q,
            physical_log_r=self.ctx.ln(nu)/2+q+(self.data['logR']+self.ctx.ln(2))/2,
            physical_z_source_coefficient=z,physical_z_log_scale=self.ctx.ln(nu)/2+(1-self.data['source']['delta'])*q,
            independent_left_right_actual_physical_bounds={side:patch.actual_physical_velocity_bounds(SimpleNamespace(ctx=self.ctx),trace,location) for side,trace in traces.items()},
            Rsh_physical_velocity_interface_definition_sha256=self.definition,
            same_actual_source_function_from_checked_Rsh_theorem=True,
            all_global_velocity_interfaces_certified=False,resolved_physical_point_values_available=False)


def run():
    field=CurrentRshPhysicalVelocityInterface(require_checked=False)
    views=dict(whole_Z=field.interface(),fresh_time_nu=field.interface(Z=('.1','.4'),log_tau=(-12,-11),viscosity='.8'),
        actual_physical_point=field.physical_interface('.237','-2.337','.337','1.3'))
    payload=(json.dumps(numeric.encode(numeric.pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        Rsh_physical_velocity_interface_definition_sha256=field.definition,exact_physical_adapter_theorem=field.theorem,
        complete_views=VIEWS_NAME,current_Rsh_physical_spatial4_time1_interface_available=True,
        arbitrary_finite_time_positive_constant_nu_actual_lambda_queries_available=True,
        all_global_velocity_interfaces_energy_cones_recursion_NS_certified=False,
        input_hashes={**field.hashes,VIEWS_NAME:sha(VIEWS_NAME)})
    (HERE/NAME).write_bytes((json.dumps(numeric.encode(numeric.pack(result)),indent=2)+'\n').encode())
    print('Rsh two-sided physical interface generated: source theorem;216 common groups; actual lambda/nu',flush=True)
    return result


if __name__=='__main__':run()
