"""Actual whole-Z Rm patch generic source and same-N correction to Rh.

The leading patch supplies the background. Its actual angular perturbation
is differentiated before common amplitude cancellation. Genuine incoming
finite-N corrections are carried separately by the original five kernels.
"""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_mixed4 as upstream

prefix = upstream.upstream.upstream
long, switch = prefix.long, prefix.upstream.switch
parameters, primitives, phase, bounds = prefix.parameters, prefix.primitives, prefix.phase, prefix.bounds
HERE, PREFIX, sha, read, bind, ep = upstream.HERE, upstream.PREFIX, upstream.sha, upstream.read, upstream.bind, upstream.ep
RATES, C0, Z, OPEN = prefix.RATES, prefix.C0, prefix.Z, upstream.OPEN
PARTITION = upstream.PARTITION
NAME = PREFIX+'current_original_whole_Z_Rm_patch_finite_N.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_Rm_patch_finite_N_check.json'
GATE = 'current_original_whole_Z_same_N_actual_Rm_patch_source_and_correction_through_Rh_enclosed'


def serialized(value):
    return upstream.serialized(value)


def raw_packet(packet):
    """Keep genuine source axial5 and only the consumed first two y rows."""
    raw = packet['raw_current_radius_y_derivative_axial_coefficients']
    return dict(raw_current_radius_y_derivative_axial_coefficients=dict(
        histories={name:rows[:2] for name, rows in raw['histories'].items()},
        velocity={name:raw['velocity'][name][:2] for name in ('theta', 'axial')}),
        original_P0_normalized_axial5=packet['original_P0_normalized_axial5'],
        geometry=packet['geometry'])


def recover_patch(op, packet):
    """Original recovery with explicit Rm*x and source-correlated C."""
    f, c = op.flow, op.c
    if packet['original_P0_normalized_axial5'] is not op.P0:
        raise ValueError('Same independent live patch P0 required')
    a, shear = parameters.correlated_shear_a(op, packet)
    raw = raw_packet(packet)
    E = raw['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][0]
    C = parameters.multiply(f, a, E)
    x = parameters.source_coordinate(c, packet['geometry'])
    radius = op.Rm_factor*x
    proxy = SimpleNamespace(flow=f, c=c, reference=op.reference, zrows=op.zrows,
        P0=op.P0, Pstar=op.Pstar, source_radius=radius, correlated_C=C)
    recovered = long.RECOVER(proxy, raw)
    recovered.pop('source_frame_conditional_on_same_actual_Rm_inlet')
    recovered.update(actual_patch_source_owner=type(op).__name__,
        source_frame_conditional_on_checked_whole_Z_Rm_leading_family=True,
        actual_source_C_is_same_function_E_minus2Ey=True,
        correlated_C_identity='E=am*(x^.1+f); C=am*(.8*x^.1+f-2*x*f_x)=a*E',
        leading_background_not_finite_N_correction=True)
    return proxy, raw, recovered, a, shear


def patch_weights(f, left, right, rate):
    """Exact log widths and suffixes for the original x=R/Rm partition."""
    c=f.c
    l=c.mpf(left[0])/left[1]
    r=None if right == 'Rh' else c.mpf(right[0])/right[1]
    yl=c.ln(l);yr=c.mpf(1) if r is None else c.ln(r)
    width=yr-yl;suffix=c.mpf(1)-yr
    mass, decay, tail=long.kernel_weight(f,width,suffix,rate)
    return width,suffix,mass,decay,tail


class WholeZRmPatchFiniteN:
    def __init__(self,dps=500):
        self.source=upstream.WholeZRmPatchMixed4(dps)
        self.c=self.source.c;self.identity=self.source.identity;self.N=self.source.N
        self.hashes=dict(self.source.hashes)
        for module, gate in ((upstream,upstream.GATE),(prefix,prefix.GATE)):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate) or receipt['source_family'] != self.identity:
                raise ValueError('Checked current whole-Z source receipt required: '+module.RECEIPT)
            for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
            bind(self.hashes,module.RECEIPT,sha(module.RECEIPT))
        if not receipt['actual_R110_incoming_all_five_C0_Z_drivers_true_weights_and_Rm_boundary_checked'] \
                or receipt['candidate_N'] != self.N:
            raise ValueError('Genuine same-N Rm correction prefix required')
        self.saved=json.loads(gzip.decompress((HERE/prefix.NAME).read_bytes()))
        if self.saved['candidate_N'] != self.N or self.saved['source_family'] != self.identity:
            raise ValueError('Same current source and unchanged candidate N required')
        self.saved_cells={tuple(row['exact_Z_cell']):row for row in self.saved['source_cells']}
        if set(self.saved_cells) != set(upstream.upstream.source.CELLS):
            raise ValueError('Complete original whole-Z Rm correction cover required')
        current_prefix=self.source.source.upstream
        self.eta_log=current_prefix.eta_log;self.dstar_log=current_prefix.dstar_log
        self.sc=read(self.c,current_prefix.source.saved['source_owned_left_inlet_proof']['selected_positive_s_c'])
        self.cache={};self.inlets={}
        self.bindings=dict(original_generic_recovery=long.RECOVERY_BINDING,
            original_generic_equations=long.generic.source_bindings(),
            original_mixed_source=upstream.mixed.source_bindings(),
            exact_correlated_shear_identity='a=.8-2*(x*f_x-f/10)/(x^.1+f)',
            original_general_quotients=Path(switch.__file__).name,
            original_five_own_rates={name:str(rate) for name,rate in RATES.items()})
        for module in (upstream,prefix,long,switch,parameters,primitives,phase,bounds,long.generic):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def owner(self,ends):
        return self.source.owner(ends).op

    def inlet(self,ends):
        key=tuple(ends)
        if key in self.inlets:return self.inlets[key]
        op=self.owner(ends);f=op.flow;saved=self.saved_cells[key]
        if saved['source_identity'] != self.identity or saved['candidate_N'] != self.N:
            raise ValueError('Same current Rm incoming source and N required')
        P0=prefix.upstream.decode(f,saved['exact_common_P0_axial5'])
        radius=prefix.upstream.decode(f,saved['exact_Rm_radius'])
        for live,old in zip(op.P0,P0):
            if live.coefficient._mpi_ != old.coefficient._mpi_ or live.scale.powers != old.scale.powers \
                    or live.scale.offset._mpi_ != old.scale.offset._mpi_:
                raise ValueError('Exact common incoming P0 source tuples required')
        if len(P0) != 6 or radius.coefficient._mpi_ != op.Rm_factor.coefficient._mpi_ \
                or radius.scale.powers != op.Rm_factor.scale.powers or radius.scale.offset._mpi_ != op.Rm_factor.scale.offset._mpi_:
            raise ValueError('Same exact actual Rm radius source required')
        incoming=prefix.upstream.decode(f,saved['actual_Rm_correction_C0_Z'])
        if set(incoming) != set(RATES) or any(len(row)!=2 for row in incoming.values()):
            raise ValueError('Genuine correction-only five C0/Z rows required')
        parameters.same_source(f,[v for rows in incoming.values() for v in rows])
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,
            actual_Rm_correction_C0_Z=incoming,
            incoming_source_report=prefix.NAME,
            common_original_function_definitions_with_source_centered_leading_owner=True,
            correction_only_not_complete_own_history=True,
            no_saved_label_or_N257_density_transplant=True)
        self.inlets[key]=result;return result

    def query(self,ends,left,right=None):
        key=(tuple(ends),left,right)
        if key in self.cache:return self.cache[key]
        op=self.owner(ends);f,c=op.flow,op.c
        packet=self.source.query(ends,left,right)
        proxy,raw,recovered,a,shear=recover_patch(op,packet)
        roots,qr,proof=switch.general_quotients(proxy,recovered,a,self.eta_log)
        got=primitives.all_u_primitive_bounds(f,roots,qr,self.dstar_log,c.mpf([0,1]))
        if qr[C0].zero and qr[Z].zero:
            got['values']={name:f.scalar(0) for name in got['values']}
            got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        A=got['values']['A'];logA=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
        if logA is not None and ep(c.ln(self.N))[0] <= logA:
            raise ValueError('Same N fails actual patch A/N budget at '+str(key))
        got=prefix.upstream.upstream.bounded_exponent_cover(f,got,self.N)
        E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
        density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            source_geometry=packet['geometry'],exact_common_P0_axial5=op.P0,
            original_raw_source=raw,original_generic_source=recovered,
            actual_correlated_shear_source=shear,original_full_source_quotients=proof,
            original_roots=roots['roots'],original_q_C0_Z=qr,
            original_primitive_values=got['values'],original_primitive_proof=got['record'],
            actual_original_A_log_absolute_upper=None if logA is None else c.mpf(logA),
            actual_A_below_same_candidate_N=True,
            original_signed_five_density_C0_Z=density,
            actual_radius_phase_definition='frac(N*(logR-logRa-hb*s_c/2))',
            actual_phase_origin_s_c=self.sc,actual_phase_Z_exact_zero=True,
            actual_radius_phase_full_period_cover=c.mpf([0,1]),
            phase_average_cancellation_not_claimed=True,
            whole_Z_actual_patch_finite_N_C0_Z_density_oracle_installed=True,
            actual_high_order_finite_N_correction_jets_installed=False)
        self.cache[key]=result;return result

    def contribution(self,ends):
        op=self.owner(ends);f,c=op.flow,op.c;inlet=self.inlet(ends)
        incoming=inlet['actual_Rm_correction_C0_Z']
        total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
        for left,right in zip(PARTITION,PARTITION[1:]):
            source=self.query(ends,left,right);density=source['original_signed_five_density_C0_Z']
            rows={};weights={}
            for name,rate in RATES.items():
                width,suffix,mass,decay,tail=patch_weights(f,left,right,rate)
                parameters.positive_source(f,mass,'actual_Rm_patch_Duhamel_mass')
                pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                    for part in ('kernels','Z_derivatives')]
                for n,row in enumerate(pair):total[name][n]+=row
                rows[name]=pair;weights[name]=dict(true_full_mass=mass,true_cell_decay=decay,
                    true_suffix_decay=tail,own_rate=str(rate),true_radial_measure_applied_once=True)
            cells.append(dict(source=source,exact_left=list(left),
                exact_right=right if right=='Rh' else list(right),
                actual_log_width=width,actual_suffix_to_Rh=suffix,
                signed_cell_driver_C0_Z=rows,own_rate_weights=weights))
            print('Whole-Z actual same-N Rm patch driver: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
        memory={name:long.kernel_weight(f,c.mpf(1),c.mpf(0),rate)[1] for name,rate in RATES.items()}
        inherited={name:[row*memory[name] for row in incoming[name]] for name in RATES}
        outgoing={name:[inherited[name][n]+total[name][n] for n in range(2)] for name in RATES}
        terminal=self.query(ends,'Rh')
        original={name:list(rows[:2]) for name,rows in terminal['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete={name:[original[name][n]+outgoing[name][n] for n in range(2)] for name in RATES}
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_Rh_radius=op.Rm_factor*c.exp(1),
            actual_Rm_incoming_binding=inlet,actual_source_cells=cells,
            exact_total_log_length=1,exact_length_identity='sum log(x_right/x_left)=log(e/1)=1',
            actual_Rm_correction_incoming_C0_Z=incoming,incoming_own_rate_memory=memory,
            actual_local_driver_C0_Z=total,actual_retained_incoming_C0_Z=inherited,
            actual_Rh_correction_C0_Z=outgoing,actual_original_Rh_background_C0_Z=original,
            actual_Rh_complete_own_history_C0_Z=complete,actual_original_Rh_source=terminal,
            actual_nonzero_incoming_memory_retained=True,
            leading_controls_not_replaced_by_finite_N_correction=True,
            full_same_N_correction_prefix_through_Rh_installed=True,
            **dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZRmPatchFiniteN()
        cells=[serialized(owner.contribution(ends)) for ends in upstream.upstream.source.CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(ends) for ends in upstream.upstream.source.CELLS],
            exact_x_domain=['1','exp(1)'],exact_x_partition=serialized(PARTITION),source_cells=cells,
            original_source_bindings=owner.bindings,original_eta_log=owner.eta_log,original_dstar_log=owner.dstar_log,
            full_same_N_correction_prefix_through_Rh_installed=True,
            actual_patch_finite_N_C0_Z_density_oracle_installed=True,
            actual_high_order_finite_N_correction_jets_installed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original Rm patch generic source, all-u signed C0/Z densities '
                  'and genuine unchanged-N correction incoming through Rh. Leading source and '
                  'correction stay distinct. Later Rref/Rc repair, high correction jets, global '
                  'N/cone/heat and genuine n-dependent recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(upstream.upstream.source.bridge._encode(report),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine same-N correction prefix reaches Rh',flush=True)
    return report


if __name__=='__main__':run()
