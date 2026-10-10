"""Complete current Rc chain to the original native Rp inlet functions.

This downstream adapter requires the full Rc bridge and actual compact/quiet
chain. True Z4/Z5 inlet rows come from the identified analytic functions,
never from padding the previous C3 packet. Pulse owner injection stays open.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge as bridge
import lei_ren_part1_paper_compliant_current_limit_Rp_native_identity as native
import lei_ren_part1_paper_compliant_current_native_spatial_phase as spatial

leading,quiet=bridge.leading,bridge.leading.quiet
current,source,target=bridge.current,bridge.source,bridge.target
HERE,PREFIX,sha=bridge.HERE,bridge.PREFIX,bridge.sha
NAME=PREFIX+'current_original_C3_Rp_native_frame.json.gz'
RECEIPT=PREFIX+'current_original_C3_Rp_native_frame_check.json'
GATES=('current_original_complete_C3_chain_to_Rp_native_frame_identity_installed',
    'current_original_Rp_native_true_Z0_to_Z5_inlet_functions_installed',
    'current_original_Rp_native_exact_radius_identity_installed')
OPEN=bridge.OPEN
FRAME_KEYS=frozenset(('u','m1','m2','X','energy','Mp','P0','pressure'))


def inputs():
    hashes={};reports={}
    for module in (bridge,quiet,native,spatial):
        receipt=bridge.outer.rh.read(module.RECEIPT)
        gates=module.GATES if hasattr(module,'GATES') else (module.GATE,)
        if not receipt['all_passed'] or not all(receipt[k] for k in gates):
            raise ValueError('Accepted actual Rp prerequisite required '+module.RECEIPT)
        for name,value in receipt['input_hashes'].items():
            if name in hashes and hashes[name]!=value:raise ValueError('Conflicting Rp input '+name)
            hashes[name]=value
        hashes[module.NAME]=sha(module.NAME);hashes[module.RECEIPT]=sha(module.RECEIPT)
        reports[module.NAME]=bridge.outer.rh.read(module.NAME)
    identity=reports[bridge.NAME]['source_family']
    if any(row['source_family']!=identity for row in reports.values()):raise ValueError('Same original Rp source family required')
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,value in hashes.items():
        if sha(name)!=value:raise ValueError('Changed actual Rp prerequisite '+name)
    return reports,hashes


def build(field):
    g=field.graph;b=field.bridge;lf=b.leading.functions;params=b.leading.parameters
    A=b.functions['actual_outer_Rc_amplitude_C3'];H=b.functions['actual_outer_Rc_leading_C3']
    mu=params['mu'];two=g.constant(2);log2=g.unary('log',two)
    Tw=g.mul(g.constant(-60),g.unary('log',mu));length=g.sub(g.sub(Tw,two),log2)
    first=quiet.power_flow(g,A,H,mu,log2)
    A2=first['profiles_y0_y1_y2_y3_y4_C3'][0]['E'];H2=first['complete_histories_y0_y1_y2_y3_y4_C3'][0]
    flow=quiet.power_flow(g,A2,H2,mu,length)
    Ep=flow['profiles_y0_y1_y2_y3_y4_C3'][0]['E'];Hp=flow['complete_histories_y0_y1_y2_y3_y4_C3'][0]
    P0=b.functions['actual_independent_P0_C3']
    c=lf['charts']['O3_power'];native_power=quiet.power_flow(g,c['seeds']['A'],c['seeds']['H'],mu,Tw)
    phase_rc=source.FunctionRef(g,b.leading.bridge.source_integer_and_Rc_binding['outer_Rc_offset_node'])
    # The source transport offset is log(R/r_minus), not log(R).
    # Restore the same original microscopic origin, without materializing h_B.
    logminus=g.add(g.unary('log',g.constant(4)),g.neg(g.mul(g.constant(4),params['logP'])),
        g.constant(-1000),g.mul(params['hbB'],params['sc'],g.constant('1/2')))
    rc=g.add(phase_rc,logminus)
    rw=g.sub(rc,two);rp=g.add(rw,Tw)
    return dict(actual_leading_2Rc_E_C3=A2,actual_leading_2Rc_histories_C3=H2,
        actual_identified_Rp_E_C3=Ep,actual_identified_Rp_histories_C3=Hp,
        actual_identified_Rp_native_frame_C3=quiet.pulse_input(g,Ep,Hp,P0),
        original_power_at_Tw_C3=native_power,Tw=Tw,quiet_length=length,
        original_Rc_phase_offset=phase_rc,original_log_r_minus=logminus,
        original_logRc=rc,original_logRw=rw,original_logRp=rp,
        source_identity_requires_actual_compact_exit_not_an_outer_shortcut=True)


class CurrentOriginalC3RpNativeFrame:
    def __init__(self,require_checked=True):
        self.reports,self.hashes=inputs();self.bridge=bridge.CurrentOuterLeadingC3Bridge()
        self.quiet=quiet.CurrentC3PowerToRp();self.identity=self.bridge.identity
        self.graph=self.bridge.graph;self.prefix=copy.deepcopy(self.graph.nodes)
        self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=bridge.outer.rh.read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual native Rp frame required')
            for name,value in receipt['input_hashes'].items():
                if sha(name)!=value:raise ValueError('Changed native Rp frame '+name)
            self.acceptance_loaded=True

    def symbolic_frame(self):
        """Exact underlying inlet functions; scalar mass covers stay bounds."""
        q=bridge.leading_check.Interpreter(self.bridge.leading)
        # Keep the actual Tw definition separate while proving the arbitrary
        # positive-duration semigroup. The checker binds it to -60 log(mu).
        Tw=s.Symbol('Tw',positive=True);q.bindings[self.functions['Tw'].node]=Tw
        z=q.z;P0=s.Function('original_P0')(z)
        rawP=self.bridge.functions['actual_independent_P0_C3']
        pressure=self.reports[native.NAME]['current_P0_source_binding']
        analytic=s.sympify(pressure['exact_analytic_function'],locals={'Z':z,'F_flat':s.Function('F_flat')})
        for j,row in enumerate(target.rows(rawP)):q.bindings[row.node]=s.diff(analytic,z,j)
        packet=self.functions['actual_identified_Rp_native_frame_C3']
        frame={key:q.at(target.rows(value)[0]) for key,value in packet.items()
            if isinstance(value,target.C3Function)}
        if set(frame)!=FRAME_KEYS:raise ValueError('Complete eight-function Rp frame required')
        allowed={q.z,q.mu,q.logP,Tw,s.Symbol('pressure_M0'),s.Symbol('pressure_M2')}
        if any(value.free_symbols-allowed for value in frame.values()):
            raise ValueError('Unbound graph coordinate or parameter in native inlet')
        return q,frame

    def incoming_exact_taylor(self,Z):
        """Six real source Taylor rows for a mathematical native inlet callback.

        This returns exact SymPy expressions, not MPIntervalTaylor point values.
        A downstream numerical owner must evaluate these same defining kernels
        and analytic pressure datum before injection into its interval runtime.
        """
        q,frame=self.symbolic_frame();value=s.sympify(Z)
        if value.is_number and (value<-1 or value>1):raise ValueError('Native Rp inlet requires Z in [-1,1]')
        return {key:[s.diff(expr,q.z,j).subs(q.z,value)/math.factorial(j) for j in range(6)]
            for key,expr in frame.items()}


def run():
    began=time.monotonic();field=CurrentOriginalC3RpNativeFrame(require_checked=False)
    q,frame=field.symbolic_frame()
    report=dict(candidate_current_original_C3_Rp_native_frame_constructed=True,source_family=field.identity,
        exact_graph_nodes=field.graph.nodes,accepted_complete_Rc_graph_prefix_length=len(field.prefix),
        actual_Rp_native_frame_functions=target.encoded(field.functions),
        actual_analytic_native_inlet_functions={key:s.sstr(value) for key,value in frame.items()},
        actual_native_Z0_to_Z5_Taylor_rows={key:[s.sstr(s.diff(value,q.z,j)/math.factorial(j)) for j in range(6)]
            for key,value in frame.items()},
        actual_common_P0_definition=field.reports[native.NAME]['current_P0_source_binding'],
        native_inlet_callback_output='exact symbolic Taylor coefficients; numerical owner injection pending',
        selected_native_pulse_constructor_consumes_current_C3_frame=False,
        native_interval_inlet_callback_installed=False,current_numeric_point_field_oracle_installed=False,
        global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(report,separators=(',',':'))+'\n').encode(),mtime=0))
    print('Current complete C3 chain and true six-row native Rp inlet functions constructed',flush=True)
    return field


if __name__=='__main__':run()
