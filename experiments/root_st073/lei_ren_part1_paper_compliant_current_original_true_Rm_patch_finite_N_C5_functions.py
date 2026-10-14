"""Actual finite-N Rm patch C5, total radial functions and signed pressure C3.

The original nonzero Rm correction and pressure-integral inlets are carried
through the exact dx/x kernels. Leading controls and corrections stay separate.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions as previous
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_long_patch_functions as original_transport
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_finite_N as original_patch

radial=previous.previous
engine=radial.previous
source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=engine.LoopAlgebra
RATES=engine.RATES
NAME=PREFIX+'current_original_true_Rm_patch_finite_N_C5_functions.json.gz'
RECEIPT=PREFIX+'current_original_true_Rm_patch_finite_N_C5_functions_check.json'
GATES=('current_original_true_Rm_patch_finite_N_loop_density_history_C5_functions_installed',
       'current_original_true_Rm_patch_finite_N_C5_transport_installed',
       'current_original_true_Rm_patch_modified_Q4_total_radial_C5_functions_installed',
       'current_original_true_rminus_to_Rh_signed_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+('current_original_true_Rm_patch_finite_N_C5_numeric_containment_installed',
    'current_original_true_six_outer_charts_finite_N_C5_functions_installed')))


def transport_log(g,densities,incoming,coordinate,lower,label):
    """Exact physical log-radius Duhamel transport, dy=dx/x."""
    a=Algebra(g,5);t=g.symbol('true_patch_C5_transport_'+label)
    width=g.sub(g.unary('log',coordinate),g.unary('log',lower));out={};parts={};memory={}
    suffix=g.sub(g.unary('log',coordinate),g.unary('log',t))
    measure=g.quotient(g.one,t,'positive original patch integration coordinate')
    for name,rate in RATES.items():
        decay=g.unary('exp',g.mul(g.constant(-rate),width)) if rate else g.one
        kernel=g.unary('exp',g.mul(g.constant(-rate),suffix)) if rate else g.one
        body=a.scale(a.at(densities[name],coordinate,t),g.mul(measure,kernel))
        contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True,
            original_rate=str(rate),original_chart=label,physical_Jacobian_applied_once=True,
            physical_measure='dt/t',exact_original_log_suffix=True,
            nonzero_original_Rm_correction_inlet_retained=True,predecessor_memory_never_reset=True)
        out[name]=a.add(a.scale(incoming[name],decay),contribution)
        parts[name]=contribution;memory[name]=decay
    return out,parts,memory


def native_derivative(leading,phase):
    g=leading.g;chart=leading.chart
    D=radial.NativeDerivative(g,leading.x,phase);D.shared['family']=leading.family
    D.shared['free'].update({i:set(v) for i,v in leading.free_source_seeds.items()})
    background_radial=radial.background_radial_rows(g,chart)
    for key,q in chart['fields'].items():
        for ref,dr in zip(rows(q),rows(background_radial[key])):
            D.seed(ref,g.mul(chart['physical_measure'],dr))
    for gamma in leading.gamma:
        D.seed(gamma[0],gamma[1]);D.seed(gamma[1],gamma[2])
    return D,background_radial


def pressure_continuation(g,S,physical,Rm,x,incoming,initial):
    a=Algebra(g,3);Q,Qy=a.trim(physical['Q']),a.trim(physical['Q_y']);R=g.mul(Rm,x)
    boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,R))
    local={};cumulative={};t=g.symbol('true_patch_pressure_C3_t')
    for key,power,weight in (('S0',1,Rm),('S1',2,g.mul(Rm,Rm,x))):
        body=a.at(a.scale(S,weight),x,t)
        local[key]=a.integral(body,t,g.one,x,original_signed_reduced_pressure=True,original_chart='Rm_patch',
            physical_radius_power=power,physical_Jacobian_applied_once=True,
            exact_patch_measure='R*dx/x=Rm*dx' if power==1 else 'R^2*dx/x=Rm^2*x*dx',
            nonzero_preceding_S0_S1_retained=True,independent_P0_retained_in_original_source=True)
        cumulative[key]=a.add(incoming[key],local[key])
    pressure=dict(m13=a.add(a.scale(cumulative['S0'],'-1/2'),a.sub(boundary['Q_plus_Q_y'],initial['Q_plus_Q_y'])),
        m15=a.sub(a.scale(cumulative['S1'],'1/2'),a.sub(boundary['R_Q_y'],initial['R_Q_y'])))
    return boundary,local,cumulative,pressure


class CurrentOriginalTrueRmPatchFiniteNC5Functions:
    def __init__(self,leading,require_checked=True):
        if type(leading) is not previous.CurrentOriginalTrueRmPatchC6Functions or not leading.acceptance_loaded:
            raise TypeError('Checked true original leading Rm patch C6/generic C5 source required')
        leading.assert_graph();self.leading=leading;self.modified=leading.modified;self.loop=self.modified.loop
        self.g=g=leading.g;self.family=copy.deepcopy(leading.family);self.hashes=dict(leading.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (previous.NAME,previous.RECEIPT,Path(original_transport.__file__).name,Path(original_patch.__file__).name):self._pin(name,sha(name))
        self.bindings=dict(original_all_N_source=original_transport.source_binding(),
            original_log_kernel=parent.active.current.ast_binding(original_patch.patch_weights),
            original_loop=self.loop.bindings['original_loop'],original_density=self.loop.bindings['original_density_graph'],
            original_pressure=self.modified.bindings['original_pressure_reduction'])
        self.N,self.eta,self.dstar,self.sc=self.loop.N,self.loop.eta,self.loop.dstar,self.loop.sc
        self.Pstar,self.P0,self.delta=leading.Pstar,leading.P0,leading.background.delta
        self.coordinate=x=leading.x;self.source_chart=bg=leading.chart;self.incoming=self.loop.Rm_correction
        self.phase_origin=self.loop.phase_origin
        Ra=leading.background.bridge.Ra
        offset=g.sub(g.unary('log',g.quotient(bg['radius'],Ra,'same original positive actual R/Ra')),self.phase_origin)
        phase=g.unary('fractional_part',g.mul(self.N,offset));self.phase=phase
        roots=engine.generic_roots(g,bg,self.Pstar)
        loop=engine.loop_functions(g,roots,phase,self.eta,self.dstar,'Rm_patch',x,self.sc)
        finite=engine.finite_N_functions(g,roots,loop,self.N)
        history,parts,memory=transport_log(g,finite['five_signed_density_increments'],self.incoming,x,g.one,'Rm_patch')
        a=Algebra(g,5);fields={key:a.trim(q) for key,q in bg['fields'].items()}
        fields.update(E=finite['E_N'],V=finite['V_N'])
        for key in RATES:fields[key]=a.add(a.trim(bg['fields'][key]),history[key])
        fields['pressure']=a.add(fields['P0'],fields['p'])
        self.chart=dict(source_chart=bg,coordinate=x,domain=bg['domain'],radius=bg['radius'],
            physical_measure=bg['physical_measure'],global_offset=offset,actual_phase=phase,
            roots=roots,loop=loop,finite_N=finite,incoming=self.incoming,history=history,
            contributions=parts,memory=memory,modified_fields=fields)
        self.derivative,self.radial_background=native_derivative(leading,phase)
        self.radial=radial.total_radial_functions(g,self.chart,self.derivative,self.N)
        self.generic=radial.generic_functions(g,fields,self.delta,self.Pstar,bg['radius'])
        self.physical=radial.physical_pressure_fields(g,fields,self.radial,self.delta,self.Pstar)
        self.S=radial.reduced_pressure_S(g,self.physical,self.delta)
        self.initial=self.modified.initial;self.incoming_integrals=self.modified.Rm['integrals']
        self.boundary,self.local_integrals,self.cumulative_integrals,self.pressure=pressure_continuation(
            g,self.S,self.physical,leading.inlet['radius'],x,self.incoming_integrals,self.initial)
        self.endpoints={key:self.at(value) for key,value in (('Rm',g.one),('Rh',bg['domain'][1]))}
        self.Rh=self.endpoints['Rh'];self.Rh_correction=self.Rh['correction']
        self.definition=dict(actual_same_N_Section11_loop_inverse_density_and_history_C5=True,
            nonzero_true_Rm_correction_and_signed_pressure_inlets_retained=True,
            exact_log_coordinate_transport_and_rate_zero_pressure_memory=True,
            original_global_phase_Ra_hbB_sc_half_not_new_patch_phase_origin=True,
            original_collar_not_applied_to_Rm_patch=True,
            separate_leading_control_and_finite_N_correction_layers=True,
            full_fast_and_slow_total_radial_C5_and_modified_generic_C4=True,
            original_gamma_first_second_and_same_implicit_control_native_sources=True,
            physical_U_Uy_mean_C5_Q_Qy_C4_and_signed_rminus_to_Rh_pressure_C3=True,
            exact_physical_Rm_and_Rm_squared_x_pressure_measures_once=True,
            independent_P0_and_original_signed_density_sectors_retained=True,
            six_outer_charts_full_preRp_second_radial_and_numeric_norms_open=True,
            all17_controls_global_n1_heat_stress_recursion_and_corrected_NS_open=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.leading,self.modified,self.loop,self.g,self.incoming,self.incoming_integrals);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or receipt['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['Rm_patch_finite_N_C5_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual finite-N patch C5/radial/pressure C3 source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed finite-N patch source '+name)
            self.acceptance_loaded=True
    def at(self,value):
        g=self.g;x=self.coordinate
        group=lambda data,order:{key:Algebra(g,order).at(q,x,value) for key,q in data.items()}
        return dict(correction=group(self.chart['history'],5),modified_fields=group(self.chart['modified_fields'],5),
            generic=group(self.generic,4),radial={key:Algebra(g,5).at(self.radial[key],x,value)
                for key in ('E_y','V_y','m_y','slow_A_y','slow_B_y')},
            physical={key:Algebra(g,4 if key in ('Q','Q_y') else 5).at(q,x,value) for key,q in self.physical.items()},
            integrals=group(self.cumulative_integrals,3),pressure=group(self.pressure,3),boundary=group(self.boundary,3),
            radius=g.node('function_substitution',expression=self.source_chart['radius'].node,
                variable=x.node,value=value.node,Z_independent_substitution=True))
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed patch input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        chart={k:v for k,v in self.chart.items() if k!='source_chart'};loop=dict(chart['loop'])
        loop['phase_partial_coefficients']=[dict(Z=k,angle=j,node=v.node) for (k,j),v in loop['phase_partial_coefficients'].items()]
        chart=dict(chart,loop=loop)
        return encode(dict(chart=chart,N=self.N,eta=self.eta,dstar=self.dstar,sc=self.sc,phase_origin=self.phase_origin,
            incoming=self.incoming,radial_background=self.radial_background,radial=self.radial,generic=self.generic,
            physical=self.physical,S=self.S,boundary=self.boundary,initial=self.initial,
            incoming_integrals=self.incoming_integrals,local_integrals=self.local_integrals,
            cumulative_integrals=self.cumulative_integrals,pressure=self.pressure,endpoints=self.endpoints,
            definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.leading.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self.snapshot!=self._snapshot() or self.hashes!=self._hashes
            or self._owners!=(self.leading,self.modified,self.loop,self.g,self.incoming,self.incoming_integrals)
            or self.source_chart is not self.leading.chart or self.chart['source_chart'] is not self.leading.chart
            or self.incoming is not self.loop.Rm_correction or self.incoming_integrals is not self.modified.Rm['integrals']
            or self.initial is not self.modified.initial):
            raise ValueError('Actual patch or nonzero original predecessor changed')


def run(leading):
    began=time.monotonic();owner=CurrentOriginalTrueRmPatchFiniteNC5Functions(leading,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,Rm_patch_finite_N_C5_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_RM_PATCH_FINITE_N_C5_RADIAL_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(leading) reuses checked true leading Rm C6 owner')
