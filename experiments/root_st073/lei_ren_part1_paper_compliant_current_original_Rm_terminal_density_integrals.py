"""Actual Rm signed C0/Z Duhamel integrals with the original dx/x measure.

Every integration cell encloses the actual source and radius phase. Local
contributions transport to Rh with their original own rates. Incoming
finite-N histories stay an explicit unsupplied affine argument; P0 stays
separate from the pressure row.
"""
import ast
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_conditioned_phase_density as previous
import lei_ren_part1_paper_compliant_current_native_local_signed_integrals as local
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as recovery

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_terminal_density_integrals.json.gz'
RECEIPT=PREFIX+'current_original_Rm_terminal_density_integrals_check.json'
GATE='original_actual_Rm_terminal_five_density_C0_Z_Duhamel_contributions_installed'
RATES={key:Fraction(value) for key,value in recovery.RATES.items()}
PARTITION=((71,40),(2,1),(9,4),(5,2),'Rh')


def source_bindings():
    if RATES!={key:Fraction(value) for key,value in previous.densities.density.RATES.items()} or RATES!=local.RATES:
        raise ValueError('Original common-unit own five rates required')
    wanted={
        recovery:('increment_densities','dict(m=u,h=e,k=V*e+E*u+e*u,e=2*V*u+u*u-E*e-e*e/2,p=E*e+e*e/2)'),
        previous.densities.density:('signed_density_kernels','dict(m=deltaV,h=deltaE,k=V*deltaE+E*deltaV+deltaE*deltaV,e=V*deltaV*2+axial_square-theta_cross-theta_square*E.ctx.mpf(".5"),p=theta_cross+theta_square*E.ctx.mpf(".5"))')}
    checks=[]
    for module,(name,expression) in wanted.items():
        tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
        fn=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name==name)
        returned=next(node for node in ast.walk(fn) if isinstance(node,ast.Return))
        if ast.dump(returned.value)!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original signed density defining expression changed')
        checks.append(dict(module=Path(module.__file__).name,function=name,expression=expression))
    return dict(original_five_increment_density_defining_returns_bound=True,defining_returns=checks,
        original_own_rates={key:str(value) for key,value in RATES.items()},
        common_units=recovery.UNITS,no_patch_repair_tail_rate_substitution=True,
        radial_measure='dy=dx/x',incoming_decay='exp(-lambda*log(right/left))',
        cell_to_terminal_suffix='exp(-lambda*log(target/right))')


def coordinate(c,value):
    if value=='Rh':return c.exp(1)
    if not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value):
        raise ValueError('Exact rational source coordinate or Rh required')
    q=Fraction(*value)
    if q<Fraction(71,40):raise ValueError('Actual terminal support x>=71/40 required')
    x=c.mpf(q.numerator)/q.denominator
    if ep(x)[1]>ep(c.exp(1))[0]:raise ValueError('Actual terminal support x<=e required')
    return x


def exact_geometry(left,right):
    out=dict(point=False,exact_left=list(left))
    if right=='Rh':out['exact_right_source']='exp(1)'
    else:out['exact_right']=list(right)
    return out


def validate_partition(c,partition):
    if not isinstance(partition,(tuple,list)) or len(partition)<2 or partition[-1]!='Rh':
        raise ValueError('Nonempty terminal partition ending at exact Rh required')
    xs=[coordinate(c,v) for v in partition]
    if partition[0]=='Rh' or any(ep(xs[j+1])[0]<=ep(xs[j])[1] for j in range(len(xs)-1)):
        raise ValueError('Strictly ordered exact terminal partition required')
    return xs


def integrate_source_cell(f,source,left,right,target,*,source_family,P0,Rm_factor):
    """Positive-weight range integral of actual whole-cell signed functions."""
    if source['source_family']!=source_family or source.get('actual_closed_radial_source_cell') is not True or source['actual_original_spatial_Z_density_interface_installed'] is not True:
        raise ValueError('Enclosed actual spatial source density required')
    geometry=source['actual_original_Rm_radius_phase']
    previous.previous.same_source(f,[Rm_factor,*P0])
    if geometry['exact_source_Rm_factor'] is not Rm_factor:
        raise ValueError('Same live actual Rm radius owner required')
    if geometry['source_geometry']!=exact_geometry(left,right):
        raise ValueError('Exact same integral/source cell required')
    if not geometry['actual_Rm_radius_phase_Z_independent']:
        raise ValueError('Z-independent original radius phase required')
    c=f.c;xl,xr,xt=(coordinate(c,v) for v in (left,right,target))
    if ep(xt)[0]<ep(xr)[0] or ep(xr)[0]<=ep(xl)[1]:raise ValueError('Ordered cell and downstream target required')
    width=c.ln(xr/xl);suffix=c.mpf(0) if right==target else c.ln(xt/xr)
    if width._mpi_!=geometry['source_radial_logarithmic_cell_width']._mpi_:
        raise ValueError('Same original logarithmic cell-width source required')
    cells=source['actual_source_bound_phase_density_cells']
    if not cells:raise ValueError('Nonempty actual phase-cell union required')
    expected_R=Rm_factor*previous.previous.source_coordinate(c,geometry['source_geometry'])
    for cell in cells:
        if cell['source_family']!=source_family or cell['source_geometry']!=geometry['source_geometry'] or cell['original_common_P0_axial5'] is not P0:
            raise ValueError('Same actual phase-cell family, geometry and P0 owner required')
        R=cell['exact_same_shared_radius_factor'];previous.previous.same_source(f,[R])
        if R.scale.powers!=expected_R.scale.powers or R.scale.offset._mpi_!=expected_R.scale.offset._mpi_ or R.coefficient._mpi_!=expected_R.coefficient._mpi_:
            raise ValueError('Same exact factored actual Rm*x source radius required')
    kernels={};jets={}
    for key in RATES:
        values=[cell['actual_original_five_signed_density_C0_Z']['kernels'][key] for cell in cells]
        derivatives=[cell['actual_original_five_signed_density_C0_Z']['Z_derivatives'][key] for cell in cells]
        previous.previous.same_source(f,values+derivatives)
        kernels[key]=local.same_source_union(values);jets[key]=local.same_source_union(derivatives)
    masses={key:local.positive_kernel_mass(c,width,rate) for key,rate in RATES.items()}
    decays={key:c.mpf(1) if not rate else c.exp(-suffix*c.mpf(rate.numerator)/rate.denominator) for key,rate in RATES.items()}
    local_rows={key:[kernels[key]*masses[key],jets[key]*masses[key]] for key in RATES}
    target_rows={key:[row*decays[key] for row in rows] for key,rows in local_rows.items()}
    return dict(source_geometry=geometry['source_geometry'],source_family=source['source_family'],
        actual_phase_source=geometry,candidate_N=geometry['candidate_N'],target_coordinate=target,
        source_radial_logarithmic_cell_width=width,source_logarithmic_suffix_to_target=suffix,
        original_positive_own_rate_masses=masses,original_own_rate_suffix_decays=decays,
        whole_cell_signed_density_C0=kernels,whole_cell_signed_density_Z=jets,
        actual_local_cell_integral_C0_Z=local_rows,actual_cell_to_target_integral_C0_Z=target_rows,
        actual_integral='integral_left^right exp(-lambda*log(target/x))*f(x,Z)*dx/x',
        source_cells_and_phase_unions_not_samples=True,no_extra_R_Pstar_or_N_factor=True,
        Z_independent_endpoints_weights_and_phase_no_boundary_terms=True,
        conditional_two_axial_source_frames_only=True,
        **dict.fromkeys(fields.previous.OPEN,False))


def affine_transport(f,incoming,contribution,width):
    """Conditional transport of a supplied five-vector; no default zero inlet."""
    if set(incoming)!=set(RATES) or set(contribution)!=set(RATES):raise ValueError('All five incoming/contribution histories required')
    c=f.c;width=c.mpf(width);result={}
    if ep(width)[0]<0 or any(not mp.isfinite(v) for v in ep(width)):
        raise ValueError('Finite nonnegative original logarithmic transport width required')
    for key,rate in RATES.items():
        if len(incoming[key])!=2 or len(contribution[key])!=2:raise ValueError('Genuine C0 and ordinary Z history rows required')
        previous.previous.same_source(f,[*incoming[key],*contribution[key]])
        decay=c.mpf(1) if not rate else c.exp(-width*c.mpf(rate.numerator)/rate.denominator)
        result[key]=[incoming[key][n]*decay+contribution[key][n] for n in range(2)]
    return result


class OriginalRmTerminalDensityIntegrals:
    mode='same_actual_Rm_terminal_C0_Z_signed_density_Duhamel_contributions_and_conditional_inlet_transport'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmConditionedPhaseDensity(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.bindings=source_bindings();self.cache={}
        for module in (local,recovery):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE):raise ValueError('Accepted original signed integration contract required')
            family=receipt.get('source_family') or {key:receipt[key] for key in previous.first.packets.FAMILY_KEYS}
            if any(family[key]!=self.upstream.parameter_family[key] for key in previous.first.packets.FAMILY_KEYS):
                raise ValueError('Same original normalized density parameter family required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,module.RECEIPT,sha(module.RECEIPT))
        for name in (Path(__file__).name,Path(local.__file__).name,Path(recovery.__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm terminal integral receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def contribution(self,label,partition=PARTITION,N=257):
        N=previous.candidate_N(N);xs=validate_partition(self.c,partition);key=(label,tuple(partition),N)
        if key in self.cache:return self.cache[key]
        op=self.upstream.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        with mp.workdps(c.dps+40):
            sources=[]
            for left,right in zip(partition,partition[1:]):
                source=self.upstream.spatial_cell(label,left,right,N)
                sources.append(integrate_source_cell(f,source,left,right,'Rh',source_family=self.family,P0=op.P0,Rm_factor=op.Rm_factor))
            rows={key:[sum((source['actual_cell_to_target_integral_C0_Z'][key][n] for source in sources),f.scalar(0)) for n in range(2)] for key in RATES}
            generic=self.upstream.upstream.upstream
            incoming=generic.evaluate(label,partition[0]);outgoing=generic.evaluate(label,'Rh')
            if incoming['common_original_P0_axial5'] is not op.P0 or outgoing['common_original_P0_axial5'] is not op.P0:
                raise ValueError('Same actual leading inlet/endpoint P0 required')
            width=c.mpf(1)-c.ln(xs[0])
            decays={key:c.mpf(1) if not rate else c.exp(-width*c.mpf(rate.numerator)/rate.denominator) for key,rate in RATES.items()}
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_radial_partition=list(partition),
                actual_cells_to_Rh=sources,actual_terminal_local_defect_integral_C0_Z=rows,
                full_original_logarithmic_interval_width=width,original_incoming_to_Rh_decays=decays,
                exact_common_P0_axial5=op.P0,
                actual_leading_inlet_memory=incoming['common_own_five_histories_axial5'],
                actual_leading_Rh_memory=outgoing['common_own_five_histories_axial5'],
                finite_N_incoming_defect_history_is_unsupplied_affine_argument=True,
                original_pressure_incoming_decay_exactly_one=True,
                exact_incoming_recipe='deltaH_j(Rh)=decay_j*deltaH_j(left)+local_integral_j',
                complete_prefix_or_corrected_Rh_history_not_claimed=True,
                actual_terminal_density_local_C0_Z_integral_enclosures_installed=True,
                signed_cell_contributions_have_downstream_suffixes=True,
                integral_enclosures_are_bounds_not_selected_values=True,
                full_axis_or_finite_N_repair_or_global_N_admitted=False,
                **dict.fromkeys(fields.previous.OPEN,False))
            self.cache[(label,tuple(partition),N)]=result;return result

    def transport_supplied_incoming(self,label,incoming,partition=PARTITION,N=257):
        packet=self.contribution(label,partition,N);f=self.upstream.upstream.upstream.upstream.owner(label).op.flow
        if incoming is None:raise ValueError('Actual or explicitly supplied incoming history cannot be silently reset')
        with mp.workdps(self.c.dps+40):
            rows=affine_transport(f,incoming,packet['actual_terminal_local_defect_integral_C0_Z'],packet['full_original_logarithmic_interval_width'])
        return dict(source_frame=label,source_family=self.family,conditional_transported_defect_C0_Z=rows,
            exact_common_P0_axial5=packet['exact_common_P0_axial5'],
            supplied_incoming_enclosures_not_a_complete_prefix_source_proof=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmTerminalDensityIntegrals(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Actual Rm terminal five-density C0/Z integrals',label,flush=True)
    result=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_density_and_rate_bindings=owner.bindings,frames=fields.serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(result),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
