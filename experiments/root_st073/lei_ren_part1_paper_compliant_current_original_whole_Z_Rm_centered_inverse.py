"""Source-centered whole-Z Rm defects and the original leading implicit map.

Cancel the prescribed affine source BEFORE interval evaluation. Retain all
coupled core rows/tails, micro/macro/first-switch increments and finite long
incoming memories. Finite-N correction histories are not leading defects.
"""
import ast
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_finite_N as upstream
import lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse as leading

current=upstream.current;source=current.source;macro=current.moments
HERE,PREFIX,sha,read,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.read,upstream.bind,upstream.ep
NAME=PREFIX+'current_original_whole_Z_Rm_centered_inverse.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_Rm_centered_inverse_check.json'
GATE='current_original_whole_Z_source_centered_Rm_defects_and_leading_axial5_inverse_enclosed'
OPEN=upstream.OPEN


def assignment(module,method,target,expression):
    tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    functions=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method]
    want=ast.dump(ast.parse(expression,mode='eval').body)
    matches=[n for fn in functions for n in ast.walk(fn) if isinstance(n,ast.Assign)
             and any(ast.unparse(t)==target for t in n.targets) and ast.dump(n.value)==want]
    if len(matches)!=1:raise ValueError('Original centered-source assignment changed: '+method+'/'+target)
    return hashlib.sha256(ast.dump(matches[0]).encode()).hexdigest()


def source_bindings():
    return dict(
        core_affine_seed=assignment(source.rebuild_module,'seed','u','4*z+self.core.j'),
        macro_V_seed=assignment(macro,'evaluate','Vin','f.add(f.V0,f.deltaV_in)'),
        macro_field=assignment(macro,'evaluate','field_V','f.add(V.evaluate(S,R0,R1),f.error_rows(V_error))'),
        first_field_increment=assignment(current.first,'evaluate','deltaV','f.add(*Vparts.values())'),
        centered_mean_identity='V=4Z+E; M_Rsh=V+mean_memory, hence M_Rsh-4Z=E+mean_memory',
        centered_axial_identity='A_Rsh=V^2+axial_memory; A_Rsh-8Z*M_Rsh+16Z^2=E^2+axial_memory-8Z*mean_memory',
        centered_mixed_identity='K_Rsh=K_inherited+V*Ktheta; H_Rsh=H_inherited+Ktheta; K_Rsh-4Z*H_Rsh=E*H_Rsh+K_inherited-V*H_inherited',
        same_first_switch_increment_and_second_post_V_constant=True,
        exact_affine_cancellation_before_interval_extension=True)


def centered_core_exit(original,Z):
    c=original.c;packet=original.rebuild.rebuild(Z,24,6)
    # The same admitted source has row0=4Z+j. It is cancelled symbolically,
    # rather than by subtracting two independent intervals for that row.
    z=source.IntervalTaylor.variable(c,Z,30);seed=4*z+original.core.j
    assert all(a._mpi_==b._mpi_ for a,b in zip(packet['fixed']['U0_Z_taylor'],seed.coefficients))
    assert all(a._mpi_==b._mpi_ for a,b in zip(packet['rows']['Uz'][0],seed.coefficients))
    core=original.core;N=packet['radial_degree'];r=c.mpf(4);values=[];tails=[]
    for k in range(6):
        finite=(core.j if k==0 else c.mpf(0))+sum(
            (packet['rows']['Uz'][n][k]*r**n for n in range(1,N+1)),c.mpf(0))
        factor=source.rebuild_module.tail_factor(c,degree=N,radial_order=0,axial_order=k,
            radius=r,h=core.h)['tail_per_Xh_norm']
        tail=core.epsilon*core.correction*factor/math.factorial(k)
        values.append(finite+source.bridge.symmetric(c,tail));tails.append(tail)
    return source.IntervalTaylor(c,values),dict(radial_degree=N,axial_order=5,
        original_j=core.j,actual_coupled_nonaffine_rows=[row[:6] for row in packet['rows']['Uz'][1:]],
        original_differentiated_Uz_tail_bounds=tails,affine_4Z_cancelled_by_source_identity=True,
        source_coupled_rows_not_replaced_by_linear_model=True)


def macro_centered_exit(f,E0):
    """The unchanged actual macro V polynomial, with only its affine seed cancelled."""
    c=f.c;S,R0,R1,empty=f.geometry((1,1));const=lambda row:macro.ExponentialPolynomial.constant(f,row)
    one=const([f.scalar(1)]+[f.scalar(0)]*5);ell=const(f.ell_in)
    for l,row in enumerate(f.d):
        primitive=macro.ExponentialPolynomial(f,{(1-l,0):f.scale(row,R0)}).primitive()
        ell=ell-primitive.scale(f.h*c.mpf('.5'))
    exp_poly=one+ell+(ell*ell).scale(c.mpf('.5'))
    E=const(f.add(f.jet(E0),f.deltaV_in));source_mass=f.scalar(0)
    for part,p in (('hydro',1),('pressure',1),('swirl',2)):
        scale=f.factor((1,int(part=='pressure'),int(part=='swirl'),0,0))
        for j,row in enumerate(f.drive[part]):
            product=f.multiply(f.q,row)
            drive=macro.ExponentialPolynomial(f,{(p-j,0):f.scale(product,f.radial_power(R0,p))})
            E=E-(drive*exp_poly).primitive().scale(scale)
            I=f.I(p,j,S,R0,R1,empty)
            mass=macro.prior.ScaledEnclosure(I.scale,macro.fields.magnitude(c,I.coefficient),f.ledger)
            source_mass+=f.norm(product)*mass*scale
    # Exact same third-order exponential remainder used by ActualMacroMoments.
    full=f.evaluate((1,1));B=full['full_prefix_log_jet_norm'];upper=f.small_upper(B)
    if ep(upper)[1]>.5:raise ValueError('Same macro full log-jet remainder admission required')
    error=source_mass*B*B*B*(c.exp(upper)/6)
    return f.add(E.evaluate(S,R0,R1),f.error_rows(error)),dict(
        actual_centered_source_polynomial_terms={str(key):row for key,row in E.terms.items()},
        original_complete_field_error_norm=error,actual_micro_deltaV_incoming=f.deltaV_in,
        exact_original_macro_geometry=dict(S=S,R0=R0,R1=R1),
        all_hydro_pressure_swirl_and_nonlinear_remainders_retained=True)


def first_increment(first_op):
    f,c=first_op.flow,first_op.c;delta={part:[f.scalar(0)]*6 for part in current.fields.PARTS};cells=[]
    for n in range(first_op.cells):
        left,right=Fraction(n,first_op.cells),Fraction(n+1,first_op.cells)
        t=first_op.box(left,right);at=first_op.angular(t)
        dforce,_,_=first_op.force_delta(t,at);mass=first_op.mass(left,right)
        for part in delta:delta[part]=f.add(delta[part],f.scale(dforce[part],mass))
        cells.append(dict(left=[left.numerator,left.denominator],right=[right.numerator,right.denominator],
            original_weighted_cutoff_mass=mass,signed_original_force_difference=dforce))
    pieces={part:f.scale(f.add(f.scale(first_op.force0[part],c.mpf('.5')),delta[part]),
                         -f.h*f.h*first_op.scales[part]) for part in delta}
    return f.add(*pieces.values()),dict(signed_original_first_switch_V_parts=pieces,source_cells=cells,
        exact_total_sigma_mass=c.mpf('.5'),original_source_measure_applied_once=True)


class WholeZRmCenteredInverse:
    def __init__(self,dps=500):
        self.upstream=upstream.WholeZLongRmFiniteN(dps);self.c=c=self.upstream.c
        self.identity=self.upstream.identity;self.hashes=dict(self.upstream.hashes);self.N=self.upstream.N
        receipt=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(upstream.GATE) or receipt['source_family']!=self.identity:
            raise ValueError('Checked actual whole-Z Rm background source required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,upstream.RECEIPT,sha(upstream.RECEIPT))
        micro_receipt=json.loads((HERE/current.RECEIPT).read_bytes())
        if not micro_receipt.get('all_passed') or not micro_receipt.get(current.GATE):
            raise ValueError('Checked same-source current R100 functions required')
        self.micro=json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
        if self.micro['source_family']!=self.identity:raise ValueError('Same whole-Z R100 family required')
        self.micro_cells={tuple(row['exact_Z_cell']):row for row in self.micro['source_cells']}
        self.repair=json.loads((HERE/(PREFIX+'five_moment_repair.json')).read_bytes())
        checked=json.loads((HERE/(PREFIX+'five_moment_repair_check.json')).read_bytes())
        admitted=json.loads((HERE/upstream.background.ADMISSION).read_bytes())
        for key,want in self.identity.items():
            if self.repair[key]!=want or admitted[key]!=want:raise ValueError('Same original fixed leading map source required')
        for gate in ('actual_implicit_functional_five_moment_closure_independently_checked',
            'corrected_partial_moments_and_same_pressure_independently_checked','residual_zero_containment_not_used_as_closure_proof'):
            if not checked.get(gate):raise ValueError('Checked original leading scalar/matrix map required')
        self.t=read(c,self.repair['axial_box_scale']);self.admitted_logP=read(c,admitted['logPstar'])
        wn='lei_ren_part1_paper_bump_integral_enclosures_check.json'
        self.wraw=json.loads((HERE/wn).read_bytes());self.W=leading.weights(c,self.wraw)
        self.fixed=json.loads((HERE/'lei_ren_part1_paper_shared_bump_constants.json').read_bytes())
        if not self.wraw.get('directed_integrals_certified') or any(ep(v)!=ep(read(c,self.fixed['fixed_matrix'][i][j]))
            for i,row in enumerate(self.W['L']) for j,v in enumerate(row)):
            raise ValueError('Original complete directed beta weights and fixed matrix required')
        self.bump=leading.SharedFiveMomentRepair.__new__(leading.SharedFiveMomentRepair);self.bump.ctx=c
        self.bump.normalization=leading.restore_value(c,self.wraw['normalization'])
        self.bump.beta_sup=read(c,self.fixed['bump_sup_bound']);self.bump.beta_deriv_sup=read(c,self.fixed['bump_first_derivative_sup_bound'])
        self.bump.full_weights={(Fraction(r['center']),Fraction(r['power']),r['multiplicity']):leading.restore_value(c,r['weight_interval']) for r in self.wraw['weight_records'].values()}
        self.bindings=source_bindings();self.owners={}
        for record in (micro_receipt,self.repair,checked):
            for name,digest in record['input_hashes'].items():bind(self.hashes,name,digest)
        for name in (current.NAME,current.RECEIPT,PREFIX+'five_moment_repair.json',PREFIX+'five_moment_repair_check.json',
            wn,'lei_ren_part1_paper_shared_bump_constants.json',Path(leading.__file__).name,Path(__file__).name):bind(self.hashes,name,sha(name))

    def owner(self,ends):
        key=tuple(ends)
        if key in self.owners:return self.owners[key]
        if key not in self.micro_cells:raise ValueError('Current whole-Z source cell required')
        op=self.upstream.owner(ends);f,c=op.flow,op.c;original=self.upstream.source.background.source
        if not ep(self.admitted_logP)[0]<=ep(op.reference.logP)[0]<=ep(op.reference.logP)[1]<=ep(self.admitted_logP)[1]:
            raise ValueError('Same admitted canonical Pstar source required')
        Ecore,core_proof=centered_core_exit(original,op.Z)
        prepared=original.actual.prepare(op.Z)
        series=current.first.FirstSwitchFunctions.__new__(current.first.FirstSwitchFunctions);series.flow,series.c=f,c
        decoder=current.micro.HydratedComparison.__new__(current.micro.HydratedComparison)
        decoder.flow,decoder.ctx,decoder.weight=f,c,f.weight;decoder.cap,decoder.Z,decoder.delta=original.cap,op.Z,original.core.delta
        decoder.data,decoder.inputs=prepared['data'],original.inputs(op.Z);decoder.comparison,decoder.cache=decoder,{}
        decoder.proof=dict(live_original_interval_inlet=True,source_frame_values_not_used=True,
            source_L3_norm_upper=decoder.data['L3_norm'],source_V3_norm_upper=decoder.data['V3_norm'])
        mic=current.micro.MicroFunctions(f,series,decoder,source.bridge._encode(op.source_proof['fresh_actual_core']))
        micro_exit=mic.prefix('second_micro',c.mpf(2));f.ell_in=micro_exit['ell'];f.deltaV_in=f.add(*micro_exit['deltaV'].values())
        E100,macro_proof=macro_centered_exit(f,Ecore)
        axis=upstream.upstream.decode(f,self.micro_cells[key]['actual_physical_R100'])
        first=current.first.FirstSwitchFunctions(f,axis['normalized_actual_fields'],axis['normalized_actual_own_six_moments'])
        increment,first_proof=first_increment(first);E110=f.add(E100,increment)
        # Enclose the SAME actual centered long endpoint; every finite
        # incoming memory remains. No zero target replaces a defect.
        term=op.terminal;shapes=term['actual_terminal_normalized_six_history_shapes'];memory=term['actual_inherited_terminal_history_contributions']
        zrows=op.reference.zrows
        mean=f.add(E110,term['actual_mean_incoming_memory'])
        axial=f.add(f.multiply(E110,E110),term['actual_axial_incoming_memory'],
                    f.scale(f.multiply(zrows,term['actual_mean_incoming_memory']),-8))
        mixed=f.add(f.multiply(E110,shapes['theta']),memory['theta_z'],
                    f.scale(f.multiply(op.R110_fields['V'],memory['theta']),-1))
        ref=upstream.background.ActualReferenceRestoreFunctions(f,op.Z,original.core.delta,shapes,
            op.R110_fields['V'],op.reference.P0,self.upstream.T,self.upstream.logC,self.upstream.provider)
        ref.E=E110;ref.E2=f.multiply(E110,E110)
        ref.centered.update(mean_error=mean,mixed_error=mixed,axial_square=axial)
        patch=leading._ActualRmLeadingPatch(ref,self.t,self.W,self.bump)
        result=SimpleNamespace(upstream=op,reference=ref,patch=patch,actual_E110=E110,
            source_centering_proof=dict(core=core_proof,macro=macro_proof,first_switch=first_proof,
                actual_E100=E100,actual_E110=E110,actual_centered_long_endpoint=ref.centered,
                original_finite_long_memory_retained=True,leading_background_not_finite_N_correction=True))
        self.owners[key]=result;print('Whole-Z actual centered Rm implicit axial5 inverse: '+str(ends),flush=True)
        return result


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZRmCenteredInverse();cells=[]
        for ends in source.CELLS:
            op=owner.owner(ends);p=op.patch
            cells.append(upstream.reference.serialized(leading.pack(dict(exact_Z_cell=list(ends),source_identity=owner.identity,
                source_centering_proof=op.source_centering_proof,coefficient_solution=p.coefficients(),
                original_partial_patch_endpoint_queries=[p.evaluate(q) for q in ((1,1),(2,1))],
                exact_common_P0_axial5=p.P0,same_common_normalized_initial_box=True))))
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N_not_used_for_leading_map=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(ends) for ends in source.CELLS],source_cells=cells,
            original_source_bindings=owner.bindings,
            whole_Z_leading_Rm_implicit_controls_solved=True,common_box_unique_smooth_leading_family=True,
            full_radial_mixed4_leading_patch_installed=False,current_same_N_finite_correction_after_Rm_installed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Original whole-Z source-centered leading Rm defects, same-common-box implicit axial5 '
                  'controls and original patch endpoint functions. Finite-N corrections stay separate. '
                  'Whole radial mixed4, same-N post-Rm/Rc transport, terminal repair, cone/heat and recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(source.bridge._encode(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z source-centered leading Rm implicit controls generated',flush=True);return report


if __name__=='__main__':run()
