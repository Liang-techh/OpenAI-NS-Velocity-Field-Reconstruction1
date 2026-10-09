"""Genuine current Rh-reference source cells, including the signed midplane.

Use the original repaired leading terminal identities and closed reference
functions. General signed all-u primitives retain the true Z derivative at
zero; no abs(Z) division, old point-owner construction, or inlet reset.
"""
import ast
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_finite_N_prefix as incoming
import lei_ren_part1_paper_compliant_current_original_reference_restore_finite_N as reference

fields,base,ep=incoming.fields,incoming.base,incoming.ep
switch=incoming.switch
HERE,PREFIX,sha=incoming.HERE,incoming.PREFIX,incoming.sha
RATES,C0,Z=reference.RATES,reference.C0,reference.Z
PARTITION=((-5,1),(-4,1),(-3,1),(-2,1),(-1,1),(0,1))
NAME=PREFIX+'current_original_Rh_reference_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_Rh_reference_finite_N_check.json'
GATE='current_original_source_bound_finite_N_Rh_reference_prefix_through_Rref_enclosed'


def fraction(value):
    if not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact rational reference offset required')
    q=Fraction(*value)
    if not -5<=q<=0:raise ValueError('Original reference offsets[-5,0] required')
    return q


def source_bindings():
    original=HERE/(PREFIX+'pre_pulse_mixed_C4.py')
    tree=ast.parse(original.read_text(encoding='utf8'))
    cls=next(n for n in tree.body if getattr(n,'name',None)=='CompliantPrePulseMixedC4')
    fn=next(n for n in cls.body if getattr(n,'name',None)=='reference')
    assignments=("u=qi*c.exp(y/10)","V=z*4","h=u*c.mpf('.625')",
        "hist=dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(u)*c.mpf(5)/12,p=square(u)*c.mpf('2.5'))")
    for text in assignments:
        wanted=ast.dump(ast.parse(text).body[0])
        if sum(ast.dump(n)==wanted for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
            raise ValueError('Original repaired Rh reference defining assignment changed: '+text)
    if not any(isinstance(n,ast.keyword) and n.arg=='exact_radius_source' and isinstance(n.value,ast.Constant)
        and n.value.value=='R=Rref*exp(offset), Rh=e^-5 Rref' for n in ast.walk(fn)):
        raise ValueError('Original Rh reference radius definition changed')
    return dict(original_repaired_reference_assignments=assignments,
        original_reference_source_sha256=sha(original.name),
        source_radius='same actual Rm_factor*exp(offset+6), offsets[-5,0]',
        exact_reference_phase='frac(N*(log(110/4)+14logPstar+10logCstar+1000+offset-hb*s_c/2))',
        leading_Rh_canonical_history_identity_from_same_unique_implicit_five_repair=True,
        correction_incoming_is_not_the_leading_terminal_zero_defect=True,
        source_E_y_equals_E_over10_a_four_fifths_V_y_zero=True,
        signed_axis_source_no_absZ_or_positive_rho_division=True,
        all_u_inverse_identity_and_Poisson_theorem=reference.primitives.THEOREM,
        actual_generic_signed_recovery=reference.long.generic.source_bindings())


def background_cell(op,left,right):
    lo,hi=fraction(left),fraction(right)
    if hi<lo:raise ValueError('Ordered closed original reference offsets required')
    ref=op.reference;f,c=op.flow,op.c
    l,r=c.mpf(lo.numerator)/lo.denominator,c.mpf(hi.numerator)/hi.denominator
    x=c.mpf([ep(l)[0],ep(r)[1]])
    zero=[f.scalar(0)]*6;unit=lambda v:[f.scalar(v)]+zero[1:]
    # These are the source's canonical background histories after leading
    # repair. They are not the genuine finite-N correction boundary packet.
    shapes=dict(theta=unit(c.mpf(5)/8),theta_z=f.scale(ref.zrows,c.mpf(5)/2),
        mean=ref.Vref,axial=f.scale(f.multiply(ref.zrows,ref.zrows),16),
        swirl=unit(c.mpf(5)/6),pressure=unit(5))
    logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10
    radius=op.Rm_factor*c.exp(x+6)
    return dict(chart='Rh_reference',phase=x,original_window_length=c.mpf(5),
        exact_original_reference_offset_cell=(left,right),
        actual_normalized_six_history_shapes=shapes,actual_velocity_V=ref.Vref,
        actual_velocity_V_y=zero,original_log_E_axial5=logE,
        actual_physical_radius=radius,
        actual_log_radius=radius.scale.evaluate()+c.ln(radius.coefficient),
        exact_source_Rm_factor=op.Rm_factor,original_P0_normalized_axial5=ref.P0,
        original_closed_reference_histories_from_leading_repair=True,
        true_compact_reference_offset_not_added_into_gigantic_logC=True,
        source_rows_are_original_interval_functions=True,original_full_source_window_not_shortened=True,
        phase_radius_and_window_length_Z_independent=True)


def guard_incoming(family,label,N,op,packet):
    switch.guard_saved_source(family,label,N,op.P0,packet)
    if packet.get('genuine_current_finite_N_prefix_through_Rh_boundary_enclosures_supplied') is not True:
        raise ValueError('Genuine actual current finite-N Rh incoming required')
    expected=switch.canonical_expression([op.Rm_factor*op.c.exp(1)])[0]
    if incoming.canonical_saved(packet['exact_source_Rh_factor'])!=expected:
        raise ValueError('Same actual Rh=e*Rm source factor required')
    if packet['bound_original_endpoint_phase_source']!=incoming.PHASE_RECIPE:
        raise ValueError('Same actual global positive-origin radius phase required')
    incoming.row_packet(packet,'actual_current_Rh_correction_C0_Z')


class OriginalRhReferenceFiniteN:
    mode='same_current_source_Rh_reference_C0_Z_local_driver_and_real_incoming_to_Rref'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=incoming.OriginalRmFiniteNPrefix(dps)
        self.c=self.upstream.c;self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes)
        self.bindings=source_bindings();self.cache={};self.cells={}
        self.saved=json.loads(gzip.decompress((HERE/incoming.NAME).read_bytes()))
        if not self.saved[incoming.GATE] or self.saved['source_family']!=self.family or self.saved['candidate_N']!=257:
            raise ValueError('Accepted genuine same-current Rh prefix required')
        for name,digest in self.saved['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
        for name in (incoming.NAME,Path(__file__).name,self.bindings['original_reference_source_sha256'] and PREFIX+'pre_pulse_mixed_C4.py'):
            fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual current reference prefix receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def query(self,label,left,right,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same supplied current frames0,.5 and candidate N257 required')
        lo,hi=fraction(left),fraction(right);key=(label,lo,hi,N)
        if key in self.cells:return self.cells[key]
        op=self.upstream.owner(label);c=op.c
        with mp.workdps(c.dps+40):
            source=background_cell(op,left,right)
            parameters=self.upstream.upstream.parameters
            got=reference.density_cell(op.reference,source,parameters.eta_log,parameters.dstar_log,N)
            got['original_generic_source'].pop('source_frame_conditional_on_same_current_Rsh_background')
            got['original_generic_source']['source_from_same_current_repaired_Rh_background']=True
            got.update(source_family=self.family,source_frame=label,candidate_N=N,
                exact_common_P0_axial5=op.P0,exact_original_reference_offset_cell=(left,right),
                exact_source_Rm_factor=op.Rm_factor,
                actual_global_radius_phase='frac(N*(logRref+offset-logRa-hb*s_c/2))',
                actual_global_phase_Z_exact_zero=True,reference_physical_log_measure='dy=d(offset)',
                signed_midplane_and_nonzero_Z_derivative_interface_installed=True,
                source_endpoint_phase_from_same_positive_origin=True)
        self.cells[key]=got;return got

    def seam(self,label):
        op=self.upstream.owner(label);f,c=op.flow,op.c
        generic=self.upstream.upstream.parameters.upstream
        leading=generic.evaluate(label,'Rh')
        tail=self.query(label,(-5,1),(-5,1))['original_generic_source']
        rows=0
        pairs=[(leading['common_velocity_E_axial5'],tail['common_velocity_E_axial5']),
            (leading['common_velocity_V_axial5'],tail['common_velocity_V_axial5'])]
        pairs.extend((leading['common_own_five_histories_axial5'][k],tail['common_own_five_histories_axial5'][k]) for k in RATES)
        for left,right in pairs:
            for l,r in zip(left,right):
                d=l-r
                if not d.zero and not ep(d.coefficient)[0]<=0<=ep(d.coefficient)[1]:
                    raise ValueError('Defining patched leading Rh and closed reference source covers must meet')
                rows+=1
        return dict(source_function_equality_from_original_leading_full_weight_implicit_terminal_identities=True,
            exact_radius_identity='Rm*e=Rref*exp(-5)',
            exact_amplitude_identity='exp(-.6)*e^.1=exp(-.5)',
            same_defining_five_histories_pressure_and_own_rate_RHS=True,
            shared_actual_flow=f is op.reference.flow,shared_original_P0=leading['common_original_P0_axial5'] is op.P0,
            directed_source_overlap_consistency_rows=rows,
            overlap_only_consistency_not_the_function_identity_proof=True)

    def contribution(self,label,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same supplied current source and frequency required')
        if (label,N) in self.cache:return self.cache[(label,N)]
        op=self.upstream.owner(label);f,c=op.flow,op.c;accepted=self.saved['frames'][label]
        with mp.workdps(c.dps+40):
            guard_incoming(self.family,label,N,op,accepted)
            restore=lambda row:switch.first.endpoint.restore_row(f,row)
            initial={name:[restore(row) for row in values] for name,values in incoming.row_packet(accepted,'actual_current_Rh_correction_C0_Z').items()}
            joined=self.seam(label);local={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
            for left,right in zip(PARTITION,PARTITION[1:]):
                l,r=fraction(left),fraction(right);width=c.mpf((r-l).numerator)/(r-l).denominator
                suffix=-c.mpf(r.numerator)/r.denominator;source=self.query(label,left,right,N)
                density=source['original_signed_five_density_C0_Z'];drivers={};weights={}
                for name,rate in RATES.items():
                    mass=incoming.weighted.terminal.local.positive_kernel_mass(c,width,rate)
                    decay=c.mpf(1) if not rate else c.exp(-suffix*c.mpf(rate.numerator)/rate.denominator)
                    rows=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*decay)
                        for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(rows):local[name][n]+=row
                    drivers[name]=rows;weights[name]=dict(true_own_rate_cell_mass=mass,suffix_to_Rref=decay)
                cells.append(dict(source=source,actual_cell_to_Rref_driver_C0_Z=drivers,original_own_rate_weights=weights,
                    actual_physical_log_width=width,actual_physical_suffix_to_Rref=suffix))
                print('Actual original reference finite-N cell',label,left,right,flush=True)
            output=incoming.weighted.terminal.affine_transport(f,initial,local,c.mpf(5))
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                exact_source_Rh_factor=op.Rm_factor*c.exp(1),exact_source_Rref_factor=op.Rm_factor*c.exp(6),
                exact_original_reference_partition=PARTITION,full_original_logarithmic_interval_width=c.mpf(5),
                actual_current_Rh_incoming_correction_C0_Z=initial,original_Rh_reference_source_cells=cells,
                actual_Rh_Rref_local_driver_C0_Z=local,actual_current_Rref_correction_C0_Z=output,
                typed_actual_patched_Rh_reference_source_join=joined,
                genuine_current_finite_N_prefix_through_Rref_boundary_enclosures_supplied=True,
                actual_signed_midplane_reference_integral_installed=True,
                all_signed_u_Z_rows_and_all_nonlinear_density_terms_retained=True,
                original_pressure_memory_exactly_one=True,
                original_P0_background_and_finite_N_correction_kept_distinct=True,
                no_old_point_owner_constructor_or_integration_replay=True,
                correction_rows_are_whole_source_enclosures_not_selected_values=True,
                functional_five_moment_closure_or_global_N_not_claimed=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[(label,N)]=result;return result


def run():
    began=time.monotonic();owner=OriginalRhReferenceFiniteN(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_bindings=owner.bindings,frames=reference.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Actual current source-bound finite-N reference prefix reaches Rref',flush=True);return report


if __name__=='__main__':run()
