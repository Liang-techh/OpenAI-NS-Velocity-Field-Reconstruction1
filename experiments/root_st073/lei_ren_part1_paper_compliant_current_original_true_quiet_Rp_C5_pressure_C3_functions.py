"""Actual compact exit, quiet C5 power flow and rminus-to-Rp pressure C3.

The original full composer keeps inner, lower collar, preRp and pulse
components separate. A partitioned Rin-to-Rp sum is exposed for inspection;
it is not a proof of the unresolved lower/active field seam or inner C3.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions as previous
import lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp as original_flow
import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives as original_quiet
import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_lower_pressure_C3_functions as lower
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as original_full

source,parent,controls=previous.source,previous.parent,previous.controls
engine,radial=previous.engine,previous.radial
outer=previous.previous
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=previous.Algebra
RATES=previous.RATES
NAME=PREFIX+'current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_true_quiet_Rp_C5_pressure_C3_functions_check.json'
GATES=('current_original_true_2Rc_to_Rp_quiet_C5_pressure_C3_functions_installed',
       'current_original_true_rminus_to_Rp_signed_pressure_C3_functions_installed',
       'current_original_true_quiet_K_matches_original_native_source_installed',
       'current_original_true_quiet_y4_C5_and_direct_Qyy_C4_functions_installed',
       'current_original_true_Rin_to_Rp_partitioned_pressure_C3_sum_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_true_quiet_Rp_C5_pressure_C3_numeric_containment_installed',
    'current_original_true_lower_active_rminus_high_order_field_join_installed',
    'current_original_n1_complete_inner_to_Rp_pressure_C3_functions_installed',
    'current_original_true_Rp_C5_selected_native_pulse_inlet_bridge_installed')))


def power_flow(g,A,own,mu,s,order):
    """Original finite-mu flow, with genuine ordinary rows through order."""
    a=Algebra(g,order);A=a.trim(A);alpha=g.add(g.constant('1/2'),mu)
    decay=lambda rate:g.unary('exp',g.neg(g.mul(g.constant(rate),s)))
    mass=lambda rate:g.mul(s,g.unary('exprel',g.neg(g.mul(rate,s))))
    E=a.scale(A,g.unary('exp',g.neg(g.mul(alpha,s))));AA=a.mul(A,A);zero=a.fixed(0)
    theta=g.mul(decay('3/2'),s,g.unary('exprel',g.mul(g.sub(g.one,mu),s)))
    H=dict(m=a.scale(a.trim(own['m']),decay(1)),
        h=a.add(a.scale(a.trim(own['h']),decay('3/2')),a.scale(A,theta)),
        k=a.scale(a.trim(own['k']),decay('3/2')),
        e=a.add(a.scale(a.trim(own['e']),decay(1)),a.scale(AA,g.mul(g.constant('-1/2'),decay(1),mass(g.mul(g.constant(2),mu))))),
        p=a.add(a.trim(own['p']),a.scale(AA,g.mul(g.constant('1/2'),mass(g.add(g.one,g.mul(g.constant(2),mu)))))))
    profiles=[dict(E=a.scale(E,g.mul(g.constant((-1)**j),*([alpha]*j))),V=zero) for j in range(5)]
    histories=[H];EE=a.mul(E,E)
    for j in range(4):
        square=g.mul(g.constant((-2)**j),*([alpha]*j))
        densities=dict(m=zero,h=profiles[j]['E'],k=zero,e=a.scale(EE,g.mul(g.constant('-1/2'),square)),
            p=a.scale(EE,g.mul(g.constant('1/2'),square)))
        histories.append({k:a.sub(densities[k],a.scale(histories[-1][k],rate)) for k,rate in RATES.items()})
    return dict(profiles_y0_to_y4=profiles,histories_y0_to_y4=histories,alpha=alpha,order=order)


def pulse_frame(g,E,H,P0):
    a=Algebra(g,5);EE=a.mul(E,E)
    div=lambda v,d:a.div(v,d,'same actual strictly positive terminal swirl')
    return dict(u=E,m1=div(H['m'],E),m2=div(H['k'],EE),X=div(H['h'],E),
        energy=div(H['e'],EE),Mp=H['p'],P0=a.trim(P0),pressure=a.add(a.trim(P0),H['p']))


class CurrentOriginalTrueQuietRpC5PressureC3Functions:
    def __init__(self,compact,lower_owner,require_checked=True):
        if type(compact) is not previous.CurrentOriginalTrueCompactC5Functions or not compact.acceptance_loaded:
            raise TypeError('Checked actual compact C5 exit required')
        if type(lower_owner) is not lower.CurrentOriginalN1AnalyticCoreLowerPressureC3Functions or not lower_owner.acceptance_loaded:
            raise TypeError('Checked actual lower pressure C3 component required')
        compact.assert_graph();lower_owner.assert_graph()
        if compact.g is not lower_owner.g or compact.family!=lower_owner.family:raise ValueError('Same original source graph required')
        self.compact,self.lower=compact,lower_owner;self.g=g=compact.g;self.family=copy.deepcopy(compact.family);self.hashes=dict(compact.hashes)
        for module in (previous,lower):
            receipt=parent.linear.read(module.RECEIPT)
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT):self._pin(name,sha(name))
        for module in (original_flow,original_quiet,original_full):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        for module in (original_flow,original_quiet):
            receipt=parent.linear.read(module.RECEIPT)
            if not receipt['all_passed'] or receipt['source_family']!=self.family:raise ValueError('Same original quiet source required')
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT):self._pin(name,sha(name))
        self.bindings=dict(original_power_flow=parent.active.current.ast_binding(original_flow.power_flow),
            original_pulse_frame=parent.active.current.ast_binding(original_flow.pulse_input),
            original_quiet_primitive=parent.active.current.ast_binding(original_quiet.CurrentOriginalN1QuietPressurePrimitives.primitive),
            original_quiet_mass_recovery=parent.active.current.ast_binding(original_quiet.CurrentOriginalN1QuietPressurePrimitives.recover_physical_q),
            original_full_composer=parent.active.current.ast_binding(original_full.compose_full_moments),
            true_compact_seed=compact.definition,original_Rh_leading_identity=compact.outer.seed_identity,
            original_six_outer_profiles=compact.outer.bindings,
            true_lower_source=lower_owner.definition)
        self.N,self.mu,self.Pstar,self.P0,self.delta=compact.N,compact.mu,compact.Pstar,compact.P0,compact.delta
        self.initial=compact.initial;self.incoming=compact.terminal
        self.leading_incoming=compact.terminal['background_fields'];self.relative_incoming=compact.terminal['correction']
        self.incoming_integrals=compact.terminal['integrals'];self.incoming_pressure=compact.terminal['pressure']
        self.s=s=g.symbol('true_post_2Rc_quiet_s');self.Tw=compact.outer.Tw
        self.length=g.sub(g.sub(self.Tw,g.constant(2)),g.unary('log',g.constant(2)))
        self.R2=self.incoming['radius'];self.R=g.mul(self.R2,g.unary('exp',s))
        self.Rp=outer.substitute(g,self.R,s,self.length)
        self.flow=power_flow(g,self.incoming['modified_fields']['E'],self.incoming['modified_fields'],self.mu,s,5)
        self.leading_flow=power_flow(g,self.leading_incoming['E'],self.leading_incoming,self.mu,s,6)
        a=Algebra(g,5);H=self.flow['histories_y0_to_y4'][0]
        self.relative={k:a.scale(v,g.unary('exp',g.mul(g.constant(-rate),s))) for k,rate in RATES.items() for v in (self.relative_incoming[k],)}
        self.fields=dict(E=self.flow['profiles_y0_to_y4'][0]['E'],V=self.flow['profiles_y0_to_y4'][0]['V'],**H,P0=a.trim(self.P0))
        self.fields['pressure']=a.add(self.fields['P0'],self.fields['p'])
        self.radial=dict(E_y=self.flow['profiles_y0_to_y4'][1]['E'],V_y=self.flow['profiles_y0_to_y4'][1]['V'],
            E_yy=self.flow['profiles_y0_to_y4'][2]['E'],V_yy=self.flow['profiles_y0_to_y4'][2]['V'],
            m_y=self.flow['histories_y0_to_y4'][1]['m'],m_yy=self.flow['histories_y0_to_y4'][2]['m'])
        self.generic=radial.generic_functions(g,self.fields,self.delta,self.Pstar,self.R)
        self.physical=radial.physical_pressure_fields(g,self.fields,self.radial,self.delta,self.Pstar)
        second=radial.physical_pressure_fields(g,dict(V=self.radial['V_y'],m=self.radial['m_y']),
            dict(V_y=self.radial['V_yy']),self.delta,self.Pstar)
        self.physical['Q_yy']=second['Q_y']
        self.S=radial.reduced_pressure_S(g,self.physical,self.delta)
        b=Algebra(g,3);self.N1p=b.add(b.scale(self.S,'-1/2'),
            b.scale(b.add(b.trim(self.physical['Q_y']),b.trim(self.physical['Q_yy'])),g.quotient(g.one,self.R,'positive actual quiet radius')))
        self.boundary,self.pressure_local,self.integrals,self.pressure=outer.pressure_geometry(g,self.S,self.physical,
            self.R,s,g.zero,g.one,self.incoming_integrals,self.initial,'true_quiet_2Rc_Rp')
        # The leading Rh mass is exactly4Z/Pstar, subsequent leading V is
        # linear in Z, and its own linear Duhamel equation preserves that
        # property. The compact relative correction is exactly zero at2Rc.
        # Evaluate the actual leading mass atZ=1 to retain its scalar source;
        # no tiny kappa or complete absolute history is replaced by zero.
        self.kappa2=g.mul(self.Pstar,outer.substitute(g,self.leading_incoming['m'].value,g.symbol('Z'),g.one))
        self.K=g.mul(self.R2,self.kappa2)
        self.linear_mass_theorem=g.node('proved_original_quiet_linear_mass_source_identity',
            actual_compact_complete_mass_rows=encode(self.incoming['modified_fields']['m']),
            actual_leading_mass_rows=encode(self.leading_incoming['m']),physical_scale=self.Pstar.node,
            scalar_kappa2=self.kappa2.node,constant_K=self.K.node,
            original_Rh_identity=encode(compact.outer.seed_identity),
            original_outer_profile_binding=encode(compact.outer.bindings),
            compact_relative_zero=compact.terminal_certificates['m'].node,
            theorem='m_Rh=4Z/Pstar; linear-Z V and own linear Duhamel preserve m=coefficient*Z; compact relative terminal mass is zero',
            current_quiet_identity='physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q',
            original_quiet_identity_receipt=original_quiet.RECEIPT,
            source_constant_not_a_saved_range_endpoint=True)
        K2=g.mul(self.K,self.K);invR2=g.quotient(g.one,self.R2,'positive physical 2Rc');invR=g.quotient(g.one,self.R,'positive quiet radius')
        self.quiet_closed=dict(m13=b.fixed(g.mul(g.constant('1/4'),K2,g.sub(invR2,invR))),
            m15=b.fixed(g.mul(g.constant('-1/4'),K2,s)))
        self.closed_pressure={k:b.add(self.incoming_pressure[k],q) for k,q in self.quiet_closed.items()}
        self.Rp_endpoint=self.at(self.length);self.inlet=self.at(g.zero)
        self.Rp_frame=pulse_frame(g,self.Rp_endpoint['fields']['E'],self.Rp_endpoint['fields'],self.P0)
        self.lower_pressure={k:b.trim(q) for k,q in lower_owner.pressure.items()}
        self.partitioned_Rin_Rp={k:b.add(self.lower_pressure[k],self.Rp_endpoint['pressure'][k]) for k in self.pressure}
        self.definition=dict(actual_compact_exit_E_complete_histories_pressure_used=True,
            quiet_original_finite_mu_C5_fields_and_y0_to_y4_histories=True,
            original_R2_and_Tw_minus2_minus_log2_geometry=True,
            actual_relative_correction_endpoint_expressions_not_reset=True,
            quiet_direct_Qyy_C4_and_N1p_C3_functions=True,
            original_linear_Z_mass_and_nonzero_scalar_K_source_retained=True,
            original_native_K_source_equality_required_by_exact_kernel_checker=True,
            complete_original_rminus_to_Rp_pressure_component_C3=True,
            lower_Rin_rminus_pressure_C3_plus_preRp_partition_sum=True,
            original_full_composer_lower_preRp_pulse_entries_not_changed=True,
            inner_0_Rin_C3_and_lower_active_field_seam_not_proved=True,
            selected_native_pulse_C5_consumes_this_Rp_frame=False,
            quantitative_high_order_norms_global_frequency_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.compact,self.lower,self.g,self.incoming,self.incoming_integrals,self.initial);self.acceptance_loaded=False
        if require_checked:
            raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not rec.get('all_passed') or not all(rec.get(k) for k in GATES) or any(rec.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['true_quiet_Rp_C5_C3_bindings']!=self.snapshot or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual quiet Rp C5/pressure C3 required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed quiet source '+name)
            self.acceptance_loaded=True
    def at(self,value):
        group=lambda data,order:{k:Algebra(self.g,order).at(q,self.s,value) for k,q in data.items()}
        return dict(fields=group(self.fields,5),leading_fields=group(self.leading_flow['histories_y0_to_y4'][0],6),
            relative=group(self.relative,5),radial=group(self.radial,5),generic=group(self.generic,4),
            physical={k:Algebra(self.g,5 if k in ('U','U_y','m','m_y') else 4).at(q,self.s,value) for k,q in self.physical.items()},
            integrals=group(self.integrals,3),pressure=group(self.pressure,3),closed_pressure=group(self.closed_pressure,3),
            boundary=group(self.boundary,3),radius=outer.substitute(self.g,self.R,self.s,value))
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed quiet input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(N=self.N,mu=self.mu,Pstar=self.Pstar,P0=self.P0,delta=self.delta,initial=self.initial,
            incoming=self.incoming,leading_incoming=self.leading_incoming,relative_incoming=self.relative_incoming,
            incoming_integrals=self.incoming_integrals,incoming_pressure=self.incoming_pressure,s=self.s,Tw=self.Tw,
            length=self.length,R2=self.R2,R=self.R,Rp=self.Rp,flow=self.flow,leading_flow=self.leading_flow,
            relative=self.relative,fields=self.fields,radial=self.radial,generic=self.generic,physical=self.physical,S=self.S,
            N1p=self.N1p,boundary=self.boundary,pressure_local=self.pressure_local,integrals=self.integrals,pressure=self.pressure,
            kappa2=self.kappa2,K=self.K,linear_mass_theorem=self.linear_mass_theorem,quiet_closed=self.quiet_closed,
            closed_pressure=self.closed_pressure,Rp_endpoint=self.Rp_endpoint,inlet=self.inlet,Rp_frame=self.Rp_frame,
            lower_pressure=self.lower_pressure,partitioned_Rin_Rp=self.partitioned_Rin_Rp,definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.compact.assert_graph();self.lower.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.compact,self.lower,self.g,self.incoming,self.incoming_integrals,self.initial)
            or self.incoming is not self.compact.terminal or self.incoming_integrals is not self.compact.terminal['integrals']
            or self.initial is not self.compact.initial):raise ValueError('Actual quiet source or compact exit changed')


def run(compact,lower_owner):
    began=time.monotonic();owner=CurrentOriginalTrueQuietRpC5PressureC3Functions(compact,lower_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,true_quiet_Rp_C5_C3_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_QUIET_RP_C5_AND_PRE_RP_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(compact,lower_owner) reuses accepted source owners')
