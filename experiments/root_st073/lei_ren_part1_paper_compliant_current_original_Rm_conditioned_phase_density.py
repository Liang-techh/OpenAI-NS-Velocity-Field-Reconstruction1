"""Same actual Rm owner: conditioned Z-only phase and five signed densities.

The original first-dual formulas are restricted to Z/angle directions,
without inventing y derivatives. Exact Rm radius phase and candidate N
remain separate from global frequency, stress-cone and integral admission.
"""
import ast
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rm_positive_quotients_q as previous
import lei_ren_part1_paper_compliant_current_native_phase_first_jets as first
import lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals as densities
import lei_ren_part1_paper_compliant_current_original_all_chart_point_phase as radius_phase

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_conditioned_phase_density.json.gz'
RECEIPT=PREFIX+'current_original_Rm_conditioned_phase_density_check.json'
GATE='original_actual_Rm_conditioned_Z_phase_and_five_signed_density_interface_installed'
OUTPUTS=('A','B_over_Pstar','A_Z','B_Z_over_Pstar','A_phi','B_phi_over_Pstar')


def compile_Z_only():
    """Restrict derivative axes; retain every original defining assignment."""
    tree=ast.parse(Path(first.__file__).read_text(encoding='utf8'))
    wanted=('First','source_dual','small_series','conditioned_first_jets')
    nodes=[copy.deepcopy(next(node for node in tree.body if getattr(node,'name',None)==name)) for name in wanted]
    edits=[]
    class Restrict(ast.NodeTransformer):
        def visit_Tuple(self,node):
            self.generic_visit(node)
            if all(isinstance(v,ast.Constant) for v in node.elts) and [v.value for v in node.elts]==['y','Z']:
                node.elts=node.elts[1:];edits.append('Z-only slow direction tuple')
            elif all(isinstance(v,ast.Constant) for v in node.elts) and [v.value for v in node.elts]==['A','B_over_Pstar','A_y','A_Z','B_y_over_Pstar','B_Z_over_Pstar']:
                node.elts=[v for v in node.elts if v.value not in ('A_y','B_y_over_Pstar')];edits.append('Z-only symmetry outputs')
            return node
        def visit_Call(self,node):
            self.generic_visit(node)
            if isinstance(node.func,ast.Attribute) and node.func.attr=='update' and [kw.arg for kw in node.keywords]==['y','Z']:
                node.keywords=node.keywords[1:];edits.append('only genuine Z source dual')
            elif isinstance(node.func,ast.Name) and node.func.id=='dict' and [kw.arg for kw in node.keywords]==['y','Z','x']:
                node.keywords=node.keywords[1:];edits.append('Z-only independent angle dual')
            return node
    nodes=[Restrict().visit(node) for node in nodes]
    expected=['only genuine Z source dual','Z-only independent angle dual','Z-only slow direction tuple',
        'Z-only slow direction tuple','Z-only symmetry outputs','Z-only slow direction tuple']
    if edits!=expected:raise ValueError('Original first-dual derivative sites changed: '+repr(edits))
    env=dict(vars(first));env.update(DIRECTIONS=('Z','x'),OUTPUTS=OUTPUTS)
    code=compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),'<original first dual; genuine Z and angle only>','exec')
    exec(code,env)
    return env['conditioned_first_jets'],dict(original_defining_phase_primitive_and_density_formulas_unchanged=True,
        only_derivative_axes_restricted=True,no_fabricated_y_source_derivatives=True,exact_axis_restrictions=edits)


Z_FIRST,BINDINGS=compile_Z_only()


def parameter_bindings():
    assignment=radius_phase.parameters.source_assignment
    checks=[assignment('lei_ren_part1_paper_logarithmic_outer_parameters.py','__init__','self.md','c.mpf(Md)','LogarithmicOuterParameters'),
        assignment('lei_ren_part1_paper_logarithmic_outer_parameters.py','__init__','self.logPstar','c.exp(self.md)+11','LogarithmicOuterParameters'),
        assignment(PREFIX+'core_transfer.py','run','datum',"CompliantPressureDatum('40')"),
        assignment(PREFIX+'core_transfer.py','run','logLambda','4*logP+1000'),
        assignment(PREFIX+'physical_norm_family.py','run','selected_logC','b.hi(2*max(endpoints(v)[1] for v in restrictions.values()))'),
        assignment(PREFIX+'global_physical_assembly.py','__init__','self.logRref','c.ln(110)+10*(self.logC+self.logP)','CompliantGlobalPhysicalAssembly'),
        assignment(PREFIX+'current_original_Rm_defect_patch_inverse.py','__init__','self.Rm_factor',
            'f.factor((0,5,0,0,0),10*reference.logC+c.ln(110)-6)','_ActualRmLeadingPatch')]
    return dict(original_parameter_defining_AST_assignments=checks,
        original_selected_parameter_recipe=first.spatial.source_assignment_theorem(),
        original_radius_and_inlet_affine_identity=first.spatial.affine_identity_theorem())


def candidate_N(N):
    if type(N) is not int or N<160 or N.bit_length()>4096:
        raise ValueError('Explicit original candidate integer N>=160 up to4096 bits required')
    return N


def source_frame(op,generic,quotients):
    f=op.flow
    if generic['common_original_P0_axial5'] is not op.P0 or quotients['common_original_P0_axial5'] is not op.P0:
        raise ValueError('Same actual Rm P0 owner required for phase/density')
    if generic['source_geometry']!=quotients['source_geometry']:
        raise ValueError('Same actual source geometry required')
    if not quotients['source_context_basis_ledger_geometry_and_P0_bound_to_same_actual_owner']:
        raise ValueError('Accepted exact owner/source binding required')
    qjet=quotients['original_shear_q']['q_axial_coefficients']
    if qjet is None:raise ValueError('Source q cutoff derivatives need refinement')
    shear=quotients['shear_quotients_axial5']
    rows=dict(a=shear['a'],t0=shear['t0'],E=generic['common_velocity_E_axial5'],
        p2=quotients['actual_factored_p1_p2_axial4']['p2'])
    for row in (*rows.values(),qjet):previous.same_source(f,row)
    roots={key:{(0,0):row[0],(0,1):row[1]} for key,row in rows.items()}
    qrows={(0,0):qjet[0],(0,1):qjet[1]}
    return dict(q=qjet[0],roots=roots),qrows


class RmRadiusPhase:
    """Direct specialisation of the accepted original actual_patch map."""
    def __init__(self,op,family,parameter_family):
        self.op=op;self.family=family;self.c=op.c;f=op.flow;c=op.c
        if f is not op.reference.flow or op.P0 is not op.reference.P0:
            raise ValueError('Same actual reference/Rm owner required')
        logP=c.exp(40)+11
        lo,hi=ep(f.logs[1]);rlo,rhi=ep(logP*2)
        if not lo<=rlo<=rhi<=hi:
            raise ValueError('Original analytic logPstar recipe must match the actual flow')
        if f.logs[3]._mpi_!=(c.ln(4)-2*f.logs[1]-1000)._mpi_:
            raise ValueError('Original logRa must match the actual flow')
        core=json.loads((HERE/(PREFIX+'core_transfer.json')).read_bytes())
        bridge=json.loads((HERE/(PREFIX+'actual_bridge_integrals.json')).read_bytes())
        read=fields.previous.read_interval
        if f.logs[0]._mpi_!=read(c,bridge['source_log_hb_enclosure'])._mpi_ or f.logs[1]._mpi_!=(2*read(c,core['logPstar']))._mpi_:
            raise ValueError('Same defining bridge width and core logPstar source tuples required')
        if any(bridge[key]!=parameter_family[key] for key in first.packets.FAMILY_KEYS) or parameter_family['actual_five_defect_family_sha256']!=family:
            raise ValueError('Same actual bridge and Rm family/source/datum required')
        self.logC=op.reference.logC
        if ep(self.logC)[0]!=ep(self.logC)[1]:raise ValueError('Same selected dyadic logC required')
        parameters=json.loads((HERE/(PREFIX+'physical_norm_family.json')).read_bytes())
        if self.logC._mpi_!=read(c,parameters['selected_logCstar'])._mpi_:
            raise ValueError('Same defining selected physical logCstar source tuple required')
        expected=f.factor((0,5,0,0,0),10*self.logC+c.ln(110)-6)
        if op.Rm_factor.scale.powers!=expected.scale.powers or op.Rm_factor.scale.offset._mpi_!=expected.scale.offset._mpi_:
            raise ValueError('Same original factored Rm reference radius required')
        collar=json.loads((HERE/first.spatial.COLLAR).read_bytes())
        collar_family=collar.get('source_family') or {key:collar[key] for key in first.packets.FAMILY_KEYS}
        if any(collar_family[key]!=parameter_family[key] for key in first.packets.FAMILY_KEYS):
            raise ValueError('Same selected inner inlet source family required')
        self.sc=first.packets.interval(c,collar['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint'])
        if ep(self.sc)[0]!=ep(self.sc)[1] or ep(self.sc)[0]<=0:
            raise ValueError('Same positive selected singleton s_c required')
        collar_width=first.packets.interval(c,collar['explicit_current_inner_exit_strict_collar']['exact_source_width_log_enclosure'])
        if collar_width._mpi_!=f.logs[0]._mpi_:
            raise ValueError('Exact same actual bridge/collar width source tuple required')
        ledger=json.loads((HERE/(PREFIX+'K1_ledger.json')).read_bytes())
        if not ledger.get('h_b_equals_epsilon_b_by_definition') or read(c,ledger['shared_positive_width_log_enclosure'])._mpi_!=f.logs[0]._mpi_:
            raise ValueError('Same defining K1 positive bridge width required')
        self.parameter_binding=dict(same_actual_bridge_width_and_core_logP_source_tuples=True,
            same_selected_logCstar_source_tuple=True,same_selected_collar_source_family=True,
            actual_logP_analytic_refinement_contained_in_same_source=True)

    def point(self,coordinate,N,decimal_digits=80):
        N=candidate_N(N)
        if type(decimal_digits) is not int or not 40<=decimal_digits<=1000:raise ValueError('Original phase precision in[40,1000] required')
        c=MPIntervalContext();c.dps=max(260,decimal_digits+(N.bit_length()*30103+99999)//100000+60)
        with mp.workdps(c.dps+40):
            if coordinate=='Rh':x=c.exp(1);logx=c.mpf(1);geometry=dict(point=True,exact_x_source='exp(1)')
            else:
                if not isinstance(coordinate,tuple) or len(coordinate)!=2:raise ValueError('Exact rational patch coordinate required')
                q=Fraction(*coordinate);x=c.mpf(q.numerator)/q.denominator;logx=c.ln(x)
                if q<1 or ep(x)[1]>ep(c.exp(1))[0]:raise ValueError('Exact actual patch x in[1,e] required')
                geometry=dict(point=True,exact_x=[q.numerator,q.denominator])
            projection=dict(full_period=False,boxes=[c.mpf(0)]);proofs=[]
            for name,value in (('constant',c.ln(c.mpf(110)/4)+994+logx),('logP',14*(c.exp(40)+11))):
                wrapped=first.spatial.ordinary_mod_one(c,value*N)
                if wrapped['full_period']:raise ArithmeticError('Actual analytic phase needs higher precision')
                pieces=[]
                for part in wrapped['boxes']:
                    added=radius_phase.periodic_add(c,projection['boxes'],part)
                    if added['full_period']:raise ArithmeticError('Actual analytic phase sum needs higher precision')
                    pieces.extend(added['boxes'])
                projection=radius_phase.periodic_add(c,pieces,c.mpf(0))
                proofs.append(dict(component=name,method='directed_analytic_modulus',projection=wrapped))
            value,proof=first.spatial.binary_mod_one(c,c.mpf(self.logC),Fraction(10*N))
            projection=radius_phase.periodic_add(c,projection['boxes'],value)
            epsilon=c.mpf(10)**-(decimal_digits+20)
            microscopic_log=c.mpf(self.op.flow.logs[0])+c.ln(N)+c.ln(c.mpf(self.sc)/2)
            if ep(microscopic_log)[1]>=ep(c.ln(epsilon))[0]:raise ArithmeticError('Actual inlet width needs phase refinement')
            error=c.mpf([-ep(epsilon)[1],0]);projection=radius_phase.periodic_add(c,projection['boxes'],error)
            if projection['full_period']:raise ArithmeticError('Actual radius phase needs refinement')
            return dict(source_family=self.family,source_geometry=geometry,candidate_N=N,
                exact_original_radius_phase='frac(N*(log(110/4)+14logPstar+10logCstar+994+logx-hb*s_c/2))',
                phase_boxes=[self.c.mpf(part) for part in projection['boxes']],
                original_regular_component_proofs=proofs,selected_dyadic_logC_modulus_proof=proof,
                signed_actual_inlet_width_phase_error=self.c.mpf(error),actual_positive_width_not_zeroed=True,
                exact_source_Rm_factor=self.op.Rm_factor,exact_source_logRa=self.op.flow.logs[3],
                actual_parameter_source_binding=self.parameter_binding,
                actual_Rm_radius_phase_Z_independent=True,phase_value_not_selected=True,
                candidate_N_not_global_frequency_admission=True)

    def cell(self,left,right,N):
        left_value=previous.source_coordinate(self.c,dict(point=True,exact_x=list(left)))
        right_value=self.c.exp(1) if right=='Rh' else previous.source_coordinate(self.c,dict(point=True,exact_x=list(right)))
        if Fraction(*left)<1 or (right!='Rh' and Fraction(*left)>=Fraction(*right)) or ep(right_value)[1]>ep(self.c.exp(1))[1]:
            raise ValueError('Strictly ordered actual patch source cell required')
        result=self.point(left,N);width=self.c.ln(right_value/left_value)
        if ep(width)[0]<=0:raise ValueError('Positive actual logarithmic cell width required')
        variation=self.c.mpf([0,ep(width*N)[1]])
        projection=radius_phase.periodic_add(self.c,result['phase_boxes'],variation)
        geometry=dict(point=False,exact_left=list(left))
        if right=='Rh':geometry['exact_right_source']='exp(1)'
        else:geometry['exact_right']=list(right)
        result.update(source_geometry=geometry,phase_boxes=projection['boxes'],whole_period_cover=projection['full_period'],
            source_radial_logarithmic_cell_width=width,actual_phase_cell_from_monotone_radius_image=True,
            source_radial_measure='dy=dx/x',phase_midpoint_not_selected=True)
        return result


class OriginalRmConditionedPhaseDensity:
    mode='same_actual_Rm_conditioned_Z_only_phase_and_original_five_signed_density_interface'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmPositiveQuotientsQ(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.bindings=BINDINGS;self.mappers={}
        self.parameter_family=self.upstream.upstream.upstream.upstream.repair
        self.parameter_bindings=parameter_bindings()
        for module in (first,first.phase,densities,densities.density,radius_phase,first.spatial):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE):
                raise ValueError('Accepted original stateless phase/density source required')
            if any(receipt['source_family'][key]!=self.parameter_family[key] for key in first.packets.FAMILY_KEYS):
                raise ValueError('Same original phase/density parameter family required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,module.RECEIPT,sha(module.RECEIPT))
        collar_receipt=first.spatial.COLLAR.replace('.json','_check.json')
        checked=json.loads((HERE/collar_receipt).read_bytes())
        if not checked.get('all_passed') or not checked.get('current_inner_exit_strict_collar_attached_to_current_source_graph_certified'):
            raise ValueError('Accepted exact selected collar parameter source required')
        for key in first.packets.FAMILY_KEYS:
            if checked[key]!=self.parameter_family[key]:raise ValueError('Same exact selected collar source family required')
        for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
        fields.previous.bind(self.hashes,collar_receipt,sha(collar_receipt))
        for name in (Path(__file__).name,Path(first.__file__).name,Path(first.phase.__file__).name,
            Path(densities.__file__).name,Path(densities.density.__file__).name,Path(radius_phase.__file__).name,
            Path(first.spatial.__file__).name,first.spatial.COLLAR,
            PREFIX+'core_transfer.json',PREFIX+'actual_bridge_integrals.json',PREFIX+'physical_norm_family.json',
            'lei_ren_part1_paper_logarithmic_outer_parameters.py',PREFIX+'core_transfer.py',
            PREFIX+'physical_norm_family.py',PREFIX+'global_physical_assembly.py',
            PREFIX+'current_original_Rm_defect_patch_inverse.py',PREFIX+'K1_ledger.json'):
            fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm conditioned phase/density receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def fixed_phase(self,label,coordinate,phi,N):
        N=candidate_N(N)
        return self.source_phase(label,self.upstream.upstream.evaluate(label,coordinate),self.upstream.evaluate(label,coordinate),phi,N)

    def source_phase(self,label,generic,quotients,phi,N):
        N=candidate_N(N);op=self.upstream.upstream.upstream.owner(label).op
        with mp.workdps(self.c.dps+40):
            source,qrows=source_frame(op,generic,quotients)
            got=Z_FIRST(source,qrows,self.upstream.dstar_log,self.c.mpf(phi));values=got['values'];density=None
            if values is not None:
                E,V=generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']
                density=densities.density_Z_kernels(E[0],E[1],V[0],V[1],values,N)
            return dict(source_geometry=generic['source_geometry'],source_family=self.family,candidate_N=N,
                original_phase_Z_only_result=got['record'],actual_original_primitive_Z_values=values,
                actual_original_five_signed_density_C0_Z=density,
                original_common_P0_axial5=op.P0,exact_same_shared_radius_factor=quotients['exact_same_shared_positive_radius_factor'],
                source_q_Z_includes_all_b_and_b_Z_terms=True,no_fabricated_y_derivative_exported=True,
                original_source_V_already_in_common_Pstar_units=True,
                source_radius_Z_independent=True,phase_is_independent_candidate_parameter=True,
                full_stress_cone_or_global_N_or_density_integrals_admitted=False,
                **dict.fromkeys(fields.previous.OPEN,False))

    def spatial_point(self,label,coordinate,N):
        N=candidate_N(N);op=self.upstream.upstream.upstream.owner(label).op
        if label not in self.mappers:self.mappers[label]=RmRadiusPhase(op,self.family,self.parameter_family)
        geometry=self.mappers[label].point(coordinate,N);cells=[]
        for phi in geometry['phase_boxes']:
            cell=self.fixed_phase(label,coordinate,phi,N)
            if cell['source_geometry']!=geometry['source_geometry']:raise ValueError('Same actual Rm radius/density geometry required')
            cell['phase_is_independent_candidate_parameter']=False;cells.append(cell)
        return dict(source_family=self.family,actual_original_Rm_radius_phase=geometry,
            actual_source_bound_phase_density_cells=cells,
            actual_original_spatial_Z_density_interface_installed=bool(cells) and all(cell['actual_original_five_signed_density_C0_Z'] is not None for cell in cells),
            actual_spatial_y_derivative_or_radial_integral_installed=False,
            **dict.fromkeys(fields.previous.OPEN,False))

    def spatial_cell(self,label,left,right,N):
        N=candidate_N(N);op=self.upstream.upstream.upstream.owner(label).op
        if label not in self.mappers:self.mappers[label]=RmRadiusPhase(op,self.family,self.parameter_family)
        geometry=self.mappers[label].cell(left,right,N)
        generic=self.upstream.upstream.cell(label,left,right);quotients=self.upstream.cell(label,left,right);cells=[]
        for phi in geometry['phase_boxes']:
            cell=self.source_phase(label,generic,quotients,phi,N)
            if cell['source_geometry']!=geometry['source_geometry']:raise ValueError('Same exact actual radial source cell required')
            cell['phase_is_independent_candidate_parameter']=False;cells.append(cell)
        return dict(source_family=self.family,actual_original_Rm_radius_phase=geometry,
            actual_source_bound_phase_density_cells=cells,actual_closed_radial_source_cell=True,
            actual_original_spatial_Z_density_interface_installed=bool(cells) and all(cell['actual_original_five_signed_density_C0_Z'] is not None for cell in cells),
            actual_spatial_y_derivative_or_radial_integral_installed=False,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmConditionedPhaseDensity(require_checked=False);frames={};N=257
    for label in ('0','.5'):
        frames[label]=dict(active=owner.spatial_point(label,(5,4),N),
            terminal=owner.spatial_point(label,(2,1),N),Rh=owner.spatial_point(label,'Rh',N),
            active_cell=owner.spatial_cell(label,(124999999,100000000),(125000001,100000000),N),
            terminal_cell=owner.spatial_cell(label,(71,40),'Rh',N))
        if not all(row['actual_original_spatial_Z_density_interface_installed'] for row in frames[label].values()):
            raise ArithmeticError('Requested actual Rm phase/density source cells require subdivision')
        print('Actual Rm source-bound phase and five signed Z densities',label,flush=True)
    result=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=N,
        original_Z_only_compiler_binding=owner.bindings,original_radius_parameter_source_bindings=owner.parameter_bindings,frames=fields.serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(result),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
