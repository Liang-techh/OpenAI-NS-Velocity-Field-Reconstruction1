"""Original O2 slope source and real current finite-N Rref input transport.

Only pure directed scalar J/mass functions are reused. Field and history
rows live in the same current source algebra; old O2 owners are not run.
"""
import ast
from fractions import Fraction
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N as previous
import lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 as original

fields,base,ep=previous.fields,previous.base,previous.ep
incoming,reference,switch=previous.incoming,previous.reference,previous.switch
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
RATES,C0,Z=previous.RATES,previous.C0,previous.Z
PARTITION=((0,1),(1,4),(1,2),(3,4),(1,1))
SCALAR_CELLS=128
NAME=PREFIX+'current_original_O2_slope_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_O2_slope_finite_N_check.json'
GATE='current_original_source_bound_finite_N_O2_slope_prefix_through_slope_exit_enclosed'


def fraction(value):
    if not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact rational original O2 slope coordinate required')
    q=Fraction(*value)
    if not 0<=q<=1:raise ValueError('Original O2 slope y in[0,1] required')
    return q


def source_bindings():
    tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if getattr(n,'name',None)=='CompliantPrePulseMixedC4')
    fn=next(n for n in cls.body if getattr(n,'name',None)=='slope')
    assignments=("J,mass=slope_masses(c,y,self.cells)","factor=c.exp(y/10-c.mpf('.6')*J)",
        "u=qi*factor","V=z*4","h=qi*(c.mpf('.625')+mass[0])*c.exp(-c.mpf('1.5')*y)",
        "hist=dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),p=square(qi)*(c.mpf('2.5')+mass[1]/2))")
    for text in assignments:
        wanted=ast.dump(ast.parse(text).body[0])
        if sum(ast.dump(n)==wanted for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
            raise ValueError('Original O2 slope field/history source changed: '+text)
    return dict(original_O2_slope_assignments=assignments,
        exact_source_E='exp(y/10-.6*J(y))/(1+Z^2)',exact_source_V='4Z',
        actual_E_y='(.1-.6*sigma(y))*E',actual_a='.8+1.2*sigma(y)',actual_b='0',
        original_dimensionless_mass_order=('theta_1.6_power1','pressure_.2_power2','energy_1.2_power2'),
        directed_scalar_helpers_only_no_old_source_owner=True,
        source_own_histories_start_from_original_Rref_background_not_finite_N_zero=True,
        actual_finite_N_input_from_same_current_Rref_prefix=True,
        exact_source_radius='same actual Rm_factor*exp(6+y)',
        exact_global_phase='frac(N*(logRref+y-logRa-hb*s_c/2))',
        source_full_signed_generic_recovery=reference.long.generic.source_bindings())


def background_cell(op,left,right,scalar_cells=SCALAR_CELLS):
    lo,hi=fraction(left),fraction(right)
    if hi<lo:raise ValueError('Ordered original O2 slope cell required')
    if type(scalar_cells) is not int or not 16<=scalar_cells<=4096:
        raise ValueError('Bounded directed scalar partition16..4096 required')
    f,c=op.flow,op.c;ref=op.reference
    y=c.mpf([ep(c.mpf(lo.numerator)/lo.denominator)[0],ep(c.mpf(hi.numerator)/hi.denominator)[1]])
    J,mass=original.slope_masses(c,y,scalar_cells)
    sig=original.sigma_jets(c,y)[0]
    factor=c.exp(y/10-c.mpf('.6')*J)
    qi=f.jet((1+ref.z*ref.z).reciprocal());qi2=f.multiply(qi,qi)
    E=f.scale(qi,factor);V=ref.Vref;zero=[f.scalar(0)]*6
    h=f.scale(qi,(c.mpf(5)/8+mass[0])*c.exp(-c.mpf('1.5')*y))
    hist=dict(m=V,h=h,k=f.multiply(h,V),
        e=f.add(f.scale(f.multiply(V,V),f.factor((0,-1,0,0,0))),
            f.scale(qi2,-(c.mpf(5)/12+mass[2]/2)*c.exp(-y))),
        p=f.scale(qi2,c.mpf(5)/2+mass[1]/2))
    Ey=f.scale(E,c.mpf('.1')-c.mpf('.6')*sig)
    dy=dict(m=f.add(V,f.scale(hist['m'],-1)),h=f.add(E,f.scale(h,-c.mpf('1.5'))),
        k=f.add(f.multiply(E,V),f.scale(hist['k'],-c.mpf('1.5'))),
        e=f.add(f.scale(f.multiply(V,V),f.factor((0,-1,0,0,0))),f.scale(hist['e'],-1),f.scale(f.multiply(E,E),-c.mpf('.5'))),
        p=f.scale(f.multiply(E,E),c.mpf('.5')))
    a=[f.scalar(c.mpf('.8')+c.mpf('1.2')*sig)]+zero[1:]
    if lo==hi==1:
        # Preserve the defining critical endpoint exactly. A tiny arithmetic
        # interval around2 would overwhelm the genuine positive eta term.
        a=[f.scalar(2)]+zero[1:]
        Ey=f.scale(E,-c.mpf('.5'))
    radius=op.Rm_factor*c.exp(6+y)
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=dict(
        histories={name:[row,dy[name]] for name,row in hist.items()},
        velocity=dict(theta=[E,Ey],axial=[V,zero])),original_P0_normalized_axial5=op.P0,
        geometry=dict(chart='O2_slope',exact_source_coordinate_cell=(left,right),actual_physical_radius=radius,
            actual_compact_source_coordinate=y,original_window_length=c.mpf(1),phase_and_radius_Z_independent=True))
    return dict(raw=packet,actual_a_axial5=a,original_J=J,original_dimensionless_masses=mass,
        original_sigma=sig,original_scalar_partition=scalar_cells,
        exact_source_Rm_factor=op.Rm_factor,actual_source_radius=radius,
        true_compact_tail_coordinate_not_added_to_gigantic_logC=True,
        original_field_and_histories_are_closed_source_functions=True)


def guard_incoming(family,label,N,op,packet):
    switch.guard_saved_source(family,label,N,op.P0,packet)
    if packet.get('genuine_current_finite_N_prefix_through_Rref_boundary_enclosures_supplied') is not True:
        raise ValueError('Genuine current finite-N Rref input required')
    radius=switch.canonical_expression([op.Rm_factor*op.c.exp(6)])[0]
    if incoming.canonical_saved(packet['exact_source_Rref_factor'])!=radius:
        raise ValueError('Same actual canonical Rref radius required')
    incoming.row_packet(packet,'actual_current_Rref_correction_C0_Z')


class OriginalO2SlopeFiniteN:
    mode='current_original_O2_slope_whole_source_and_actual_finite_N_Rref_incoming'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRhReferenceFiniteN(dps)
        self.c=self.upstream.c;self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes)
        self.bindings=source_bindings();self.cache={};self.cells={}
        self.saved=json.loads(gzip.decompress((HERE/previous.NAME).read_bytes()))
        if not self.saved[previous.GATE] or self.saved['source_family']!=self.family or self.saved['candidate_N']!=257:
            raise ValueError('Accepted same-current actual Rref prefix required')
        for name,digest in self.saved['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
        for module in (original,reference,switch):fields.previous.bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        for name in (previous.NAME,Path(__file__).name,'lei_ren_part1_paper_interval_outer_slope_field.py'):
            fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual O2 slope finite-N prefix receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):return self.upstream.upstream.owner(label)

    def query(self,label,left,right,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same current frames0,.5,N257 required')
        lo,hi=fraction(left),fraction(right);key=(label,lo,hi,N)
        if key in self.cells:return self.cells[key]
        op=self.owner(label);f,c=op.flow,op.c
        with mp.workdps(c.dps+40):
            source=background_cell(op,left,right);a=source['actual_a_axial5']
            E=source['raw']['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][0]
            proxy=SimpleNamespace(flow=f,c=c,reference=op.reference,zrows=op.zrows,P0=op.P0,
                Pstar=op.Pstar,source_radius=source['actual_source_radius'],correlated_C=f.multiply(a,E))
            generic=reference.long.RECOVER(proxy,source['raw'])
            generic.pop('source_frame_conditional_on_same_actual_Rm_inlet')
            generic['source_from_same_current_original_O2_slope_background']=True
            parameters=self.upstream.upstream.upstream.parameters
            roots,qr,quotients=switch.general_quotients(proxy,generic,a,parameters.eta_log)
            got=reference.primitives.all_u_primitive_bounds(f,roots,qr,parameters.dstar_log,c.mpf([0,1]))
            if qr[C0].zero and qr[Z].zero:
                got['values']={name:f.scalar(0) for name in got['values']}
                got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
            V=generic['common_velocity_V_axial5']
            density=reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                original_closed_O2_slope_background=source,original_generic_source=generic,
                original_full_signed_quotients=quotients,original_roots=roots['roots'],original_q_rows=qr,
                original_primitive_values=got['values'],original_primitive_proof=got['record'],
                original_signed_five_density_C0_Z=density,
                exact_source_radius=source['actual_source_radius'],actual_compact_source_coordinate_cell=(left,right),
                actual_global_radius_phase='frac(N*(logRref+y-logRa-hb*s_c/2))',
                actual_phase_full_period_cover=c.mpf([0,1]),actual_phase_Z_exact_zero=True,
                actual_physical_log_measure='dy=d(y)',all_signed_u_midplane_Z_derivative_retained=True,
                original_variable_a_and_nonzero_inertia_retained=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cells[key]=result;return result

    def seam(self,label):
        op=self.owner(label);ref=self.upstream.query(label,(0,1),(0,1))['original_generic_source']
        slope=self.query(label,(0,1),(0,1))['original_generic_source'];rows=0
        pairs=[(ref['common_velocity_E_axial5'],slope['common_velocity_E_axial5']),
            (ref['common_velocity_V_axial5'],slope['common_velocity_V_axial5'])]
        pairs.extend((ref['common_own_five_histories_axial5'][k],slope['common_own_five_histories_axial5'][k]) for k in RATES)
        for left,right in pairs:
            for l,r in zip(left,right):
                d=l-r
                if not d.zero and not ep(d.coefficient)[0]<=0<=ep(d.coefficient)[1]:
                    raise ValueError('Original shared Rref background function covers must meet')
                rows+=1
        return dict(exact_original_function_join_by_J_and_masses_zero_at_y0=True,
            same_original_Rref_P0_flow_and_five_history_seed=True,
            source_overlap_consistency_rows=rows,overlap_not_function_identity_proof=True,
            exact_radius_and_phase_join='Rref*exp(0)=Rref; y=0 is reference offset0')

    def contribution(self,label,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same supplied current frame/frequency required')
        if (label,N) in self.cache:return self.cache[(label,N)]
        op=self.owner(label);f,c=op.flow,op.c;accepted=self.saved['frames'][label]
        with mp.workdps(c.dps+40):
            guard_incoming(self.family,label,N,op,accepted)
            restore=lambda row:switch.first.endpoint.restore_row(f,row)
            initial={name:[restore(row) for row in values] for name,values in incoming.row_packet(accepted,'actual_current_Rref_correction_C0_Z').items()}
            joined=self.seam(label);local={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
            for left,right in zip(PARTITION,PARTITION[1:]):
                l,r=fraction(left),fraction(right);width=c.mpf((r-l).numerator)/(r-l).denominator
                suffix=1-c.mpf(r.numerator)/r.denominator;source=self.query(label,left,right,N)
                density=source['original_signed_five_density_C0_Z'];drivers={};weights={}
                for name,rate in RATES.items():
                    mass=incoming.weighted.terminal.local.positive_kernel_mass(c,width,rate)
                    decay=c.mpf(1) if not rate else c.exp(-suffix*c.mpf(rate.numerator)/rate.denominator)
                    pair=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*decay) for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):local[name][n]+=row
                    drivers[name]=pair;weights[name]=dict(actual_own_rate_mass=mass,actual_suffix_to_slope_exit=decay)
                cells.append(dict(source=source,actual_cell_to_slope_exit_driver_C0_Z=drivers,
                    actual_original_log_width=width,actual_original_suffix=suffix,original_own_rate_weights=weights))
                print('Actual current O2 slope finite-N source cell',label,left,right,flush=True)
            output=incoming.weighted.terminal.affine_transport(f,initial,local,c.mpf(1))
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                exact_source_Rref_factor=op.Rm_factor*c.exp(6),exact_source_O2_slope_exit_factor=op.Rm_factor*c.exp(7),
                exact_original_slope_partition=PARTITION,full_original_logarithmic_interval_width=c.mpf(1),
                actual_current_Rref_incoming_correction_C0_Z=initial,actual_O2_slope_source_cells=cells,
                actual_O2_slope_local_driver_C0_Z=local,actual_current_O2_slope_exit_correction_C0_Z=output,
                actual_O2_slope_exit_background_source=self.query(label,(1,1),(1,1)),
                typed_actual_Rref_O2_source_join=joined,
                genuine_current_finite_N_prefix_through_O2_slope_exit_supplied=True,
                original_scalar_J_and_all_three_history_integrals_retained=True,
                variable_a_critical_end_and_nonzero_midplane_Z_inertia_retained=True,
                original_P0_background_and_actual_correction_kept_separate=True,
                original_pressure_memory_exactly_one=True,no_old_O2_owner_or_source_replay=True,
                source_enclosures_not_chosen_field_values_or_functional_closure=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[(label,N)]=result;return result


def run():
    began=time.monotonic();owner=OriginalO2SlopeFiniteN(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_bindings=owner.bindings,frames=reference.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Actual current finite-N prefix reaches original O2 slope exit',flush=True);return report


if __name__=='__main__':run()
