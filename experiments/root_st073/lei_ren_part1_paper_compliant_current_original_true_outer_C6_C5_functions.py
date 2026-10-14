"""Six true outer C6 sources and actual finite-N C5 continuation to Rc.

The actual Rh leading and correction inlets remain separate. Original
physical log-radius maps govern all kernels, including exponential O2.
The compact Rc-to-2Rc repair and quiet extension to Rp remain subsequent.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions as previous
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity as original_leading
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions as original_outer
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection as original_geometry
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join as original_Rh

engine,radial=previous.engine,previous.radial
source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=engine.LoopAlgebra
RATES=engine.RATES
CHARTS=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')
NAME=PREFIX+'current_original_true_outer_C6_C5_functions.json.gz'
RECEIPT=PREFIX+'current_original_true_outer_C6_C5_functions_check.json'
GATES=('current_original_true_six_outer_prescribed_C6_generic_C5_functions_installed',
       'current_original_all_17_chart_axial6_generic_source_installed',
       'current_original_true_six_outer_charts_finite_N_C5_functions_installed',
       'current_original_all_17_chart_finite_N_C5_functions_installed',
       'current_original_true_all17_modified_Q4_total_radial_C5_functions_installed',
       'current_original_true_rminus_to_Rc_signed_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_true_outer_C6_C5_numeric_containment_installed',
    'current_original_true_Rc_to_Rp_compact_C5_pressure_C3_functions_installed')))


def substitute(g,ref,variable,value):
    return g.node('function_substitution',expression=ref.node,variable=variable.node,value=value.node,Z_independent_substitution=True)


def geometry(g,x,chart):
    c=g.constant;Y=g.unary('exp',c(40))
    offset={'Rh_reference':g.add(c(6),x),'O2_slope':g.add(c(6),x),
        'O2_axial':g.add(c(6),g.unary('exp',g.mul(c(40),x))),
        'O2_buffer':g.add(c(6),Y,x),'O3_transition':g.add(c(17),Y,x),
        'O3_power':g.add(c(18),Y,x)}[chart]
    jac=g.mul(c(40),g.unary('exp',g.mul(c(40),x))) if chart=='O2_axial' else g.one
    domains={'Rh_reference':(c(-5),g.zero),'O2_buffer':(g.zero,c(11)),'O3_power':(g.zero,c(2))}
    return offset,jac,domains.get(chart,(g.zero,g.one))


def transport_geometry(g,density,incoming,x,lower,offset,jacobian,label,order):
    """Duhamel with actual physical log coordinate, once in each measure."""
    a=Algebra(g,order);t=g.symbol('true_outer_transport_'+label+'_'+str(order))
    origin=substitute(g,offset,x,lower);at_t=substitute(g,offset,x,t)
    width=g.sub(offset,origin);suffix=g.sub(offset,at_t);measure=substitute(g,jacobian,x,t)
    out={};local={};memory={}
    for key,rate in RATES.items():
        decay=g.unary('exp',g.mul(g.constant(-rate),width)) if rate else g.one
        kernel=g.unary('exp',g.mul(g.constant(-rate),suffix)) if rate else g.one
        body=a.scale(a.at(density[key],x,t),g.mul(measure,kernel))
        local[key]=a.integral(body,t,lower,x,original_outer_Duhamel=True,original_chart=label,
            original_rate=str(rate),physical_Jacobian_applied_once=True,
            physical_log_offset=offset.node,physical_Jacobian=jacobian.node,
            predecessor_memory_never_reset=True,leading_background=order==6,finite_N_correction=order==5)
        memory[key]=decay;out[key]=a.add(a.scale(incoming[key],decay),local[key])
    return out,local,memory


def background_density(g,E,V):
    a=Algebra(g,6);EE=a.mul(E,E)
    return dict(m=V,h=E,k=a.mul(E,V),e=a.sub(a.mul(V,V),a.scale(EE,'1/2')),p=a.scale(EE,'1/2'))


def pressure_geometry(g,S,physical,R,x,lower,jac,incoming,initial,label):
    a=Algebra(g,3);Q,Qy=a.trim(physical['Q']),a.trim(physical['Q_y'])
    boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,R));local={};cumulative={}
    t=g.symbol('true_outer_pressure_'+label)
    for key,power in (('S0',1),('S1',2)):
        weight=g.mul(jac,*([R]*power));body=a.at(a.scale(S,weight),x,t)
        local[key]=a.integral(body,t,lower,x,original_signed_reduced_pressure=True,original_chart=label,
            physical_radius_power=power,physical_Jacobian_applied_once=True,
            nonzero_preceding_S0_S1_retained=True,independent_P0_retained_in_original_source=True)
        cumulative[key]=a.add(incoming[key],local[key])
    pressure=dict(m13=a.add(a.scale(cumulative['S0'],'-1/2'),a.sub(boundary['Q_plus_Q_y'],initial['Q_plus_Q_y'])),
        m15=a.sub(a.scale(cumulative['S1'],'1/2'),a.sub(boundary['R_Q_y'],initial['R_Q_y'])))
    return boundary,local,cumulative,pressure


class CurrentOriginalTrueOuterC6C5Functions:
    def __init__(self,patch,require_checked=True):
        if type(patch) is not previous.CurrentOriginalTrueRmPatchFiniteNC5Functions or not patch.acceptance_loaded:
            raise TypeError('Checked actual Rh leading/correction/pressure source required')
        patch.assert_graph();self.patch=patch;self.leading=patch.leading;self.g=g=patch.g
        self.family=copy.deepcopy(patch.family);self.hashes=dict(patch.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for module in (original_leading,original_outer,original_geometry,original_Rh):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        self.bindings=dict(original_leading=original_leading.source_contracts(),
            original_outer_frontends=original_outer.FRONTEND_BINDINGS,
            original_geometry=parent.active.current.ast_binding(original_geometry.exact_geometry),
            original_Rh_leading_identity=original_Rh.source_bindings(),
            original_normalized_history_ODE=parent.active.current.ast_binding(original_leading.own_first_y),
            original_mu_admission=parent.active.current.ast_binding(original_outer.o3.source_mu_admission),
            original_power=parent.active.current.ast_binding(original_leading.quiet.power_flow),
            original_loop=patch.bindings['original_loop'],original_density=patch.bindings['original_density'],
            original_pressure=patch.bindings['original_pressure'])
        self.N,self.eta,self.dstar,self.sc=patch.N,patch.eta,patch.dstar,patch.sc
        self.Pstar,self.P0,self.delta=patch.Pstar,patch.P0,patch.delta
        self.phase_origin=patch.phase_origin;self.Rm=self.leading.inlet['radius'];self.Ra=self.leading.background.bridge.Ra
        self.logmu=g.sub(g.unary('log',g.constant('1/1000')),g.mul(g.constant(4),g.add(g.unary('exp',g.constant(40)),g.constant(11))))
        self.mu=g.unary('exp',self.logmu);self.Tw=g.mul(g.constant(-60),self.logmu)
        c=self.leading.background.bridge.inlet.common.ctx
        proof=g.nodes[self.eta.node]['canonical_parameter_proof']
        eta_log=previous.previous.leading.restore_value(c,proof['original_canonical_eta_log'])
        logmu,Tw,margin=original_outer.o3.source_mu_admission(c,eta_log)
        self.mu_proof=dict(source_family=self.family,eta_proof=copy.deepcopy(proof),
            exact_original_mu_recipe='mu=1/1000*exp(-4*(exp(40)+11))',
            positive_mu_and_power_2mu_ge_eta_and_Rp_reservation=True,
            original_admission_binding=self.bindings['original_mu_admission'],bounds_not_field_values=True)
        self.leading_Rh=self.leading.endpoints['Rh']
        # An identity of smooth functions supplies all ordinary derivatives.
        # The old C1 correction is deliberately not an input to this adapter.
        Rh_receipt=parent.linear.read(original_Rh.RECEIPT)
        Rh_source=parent.linear.read(original_Rh.NAME)
        leading_receipt=parent.linear.read(previous.previous.RECEIPT)
        if (not Rh_receipt.get('all_passed') or not Rh_receipt.get(original_Rh.GATES[0])
            or Rh_receipt['source_family']!=self.family or Rh_receipt['fixed_leading_input_N0']!=self.leading.fixed_N0
            or Rh_source['actual_source_bindings']!=original_Rh.serialized(self.bindings['original_Rh_leading_identity'])
            or not leading_receipt.get('same_exact_original_unique_branch_theorem_bound')
            or leading_receipt['source_family']!=self.family):
            raise ValueError('Same original Rh leading function identity and true C6 branch required')
        for name in (original_Rh.RECEIPT,original_Rh.NAME):self._pin(name,sha(name))
        for name,digest in Rh_receipt['input_hashes'].items():self._pin(name,digest)
        self.seed_identity=dict(actual_Rh_fields=self.leading_Rh['fields'],
            true_C6_endpoint_receipt=previous.previous.RECEIPT,
            same_exact_unique_branch_and_C6_implicit_derivatives=True,
            original_whole_Z_leading_identity_receipt=original_Rh.RECEIPT,
            original_full_weight_function_identity=original_Rh.exact_join_proof(),
            canonical_normalized_Rh='E=exp(-1/2)/(1+Z^2); V=m=4Z/Pstar; h=5E/8; k=h*V; e=V^2-5E^2/12; p=5E^2/2',
            original_normalized_history_ODE_binding=self.bindings['original_normalized_history_ODE'],
            propagated_original_histories_identified_by_linear_ODE_uniqueness=True,
            legacy_C1_or_all_N_correction_never_used_as_selected_N_inlet=True,
            function_identity_does_not_admit_numeric_C6_norms=True)
        self.initial_background={key:self.leading_Rh['fields'][key] for key in RATES}
        self.initial_correction=patch.Rh_correction;self.initial_integrals=patch.Rh['integrals'];self.initial=patch.initial
        self.background={};self.background_endpoints={};self.charts={};self.endpoints={};self.derivatives={}
        a=Algebra(g,6);incoming=self.initial_background;correction=self.initial_correction;integrals=self.initial_integrals
        sigma=engine.previous.sigma;prime=engine.previous.sigma_prime
        for name in CHARTS:
            x=g.symbol('true_outer_'+name);offset,jac,domain=geometry(g,x,name)
            if name=='Rh_reference':
                E=a.scale(self.leading_Rh['fields']['E'],g.unary('exp',g.mul(g.constant('1/10'),g.add(x,g.constant(5)))))
                V=a.trim(self.leading_Rh['fields']['V']);Ey=a.scale(E,'1/10');Vy=a.fixed(0)
            elif name=='O2_slope':
                t=g.symbol('true_outer_slope_sigma_t')
                J=controls.integral(g,sigma(g,t),g.nodes[t.node]['name'],g.zero,x,exact_original_sigma_integral=True)
                E=a.scale(self.background_endpoints['Rh_reference']['fields']['E'],
                    g.unary('exp',g.sub(g.mul(g.constant('1/10'),x),g.mul(g.constant('3/5'),J))))
                V=a.trim(self.leading_Rh['fields']['V'])
                Ey=a.scale(E,g.sub(g.constant('1/10'),g.mul(g.constant('3/5'),sigma(g,x))));Vy=a.fixed(0)
            elif name=='O2_axial':
                distance=g.sub(g.unary('exp',g.mul(g.constant(40),x)),g.one)
                E=a.scale(self.background_endpoints['O2_slope']['fields']['E'],g.unary('exp',g.mul(g.constant('-1/2'),distance)))
                beta=sigma(g,g.sub(g.one,x));V=a.scale(self.leading_Rh['fields']['V'],beta)
                Ey=a.scale(E,'-1/2')
                Vy=a.scale(self.leading_Rh['fields']['V'],g.quotient(g.neg(prime(g,g.sub(g.one,x))),jac,'positive original O2 exponential Jacobian'))
            elif name=='O2_buffer':
                E=a.scale(self.background_endpoints['O2_axial']['fields']['E'],g.unary('exp',g.mul(g.constant('-1/2'),x)))
                V=a.fixed(0);Ey=a.scale(E,'-1/2');Vy=a.fixed(0)
            elif name=='O3_transition':
                t=g.symbol('true_outer_transition_sigma_t')
                J=controls.integral(g,sigma(g,t),g.nodes[t.node]['name'],g.zero,x,exact_original_sigma_integral=True)
                E=a.scale(self.background_endpoints['O2_buffer']['fields']['E'],
                    g.unary('exp',g.sub(g.mul(g.constant('-1/2'),x),g.mul(self.mu,J))))
                V=a.fixed(0);Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),g.mul(self.mu,sigma(g,x)))));Vy=a.fixed(0)
            else:
                E=a.scale(self.background_endpoints['O3_transition']['fields']['E'],
                    g.unary('exp',g.neg(g.mul(g.add(g.constant('1/2'),self.mu),x))))
                V=a.fixed(0);Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),self.mu)));Vy=a.fixed(0)
            H,parts,memory=transport_geometry(g,background_density(g,E,V),incoming,x,domain[0],offset,jac,name,6)
            fields=dict(E=E,V=V,**H,P0=a.trim(self.P0));fields['pressure']=a.add(fields['P0'],fields['p'])
            R=g.mul(self.Rm,g.unary('exp',offset));generic=radial.generic_functions(g,fields,self.delta,self.Pstar,R,order=5)
            bg=dict(phase=x,domain=list(domain),radius=R,physical_measure=jac,log_offset=offset,
                fields=fields,radial=dict(E_y=Ey,V_y=Vy),generic=generic,incoming=incoming,
                history_contributions=parts,memory=memory)
            self.background[name]=bg;self.background_endpoints[name]=self.background_at(name,domain[1]);incoming=self.background_endpoints[name]['history']
            roots=engine.generic_roots(g,bg,self.Pstar)
            phase=g.unary('fractional_part',g.mul(self.N,g.sub(g.unary('log',g.quotient(R,self.Ra,'same original positive actual R/Ra')),self.phase_origin)))
            loop=engine.loop_functions(g,roots,phase,self.eta,self.dstar,name,x,self.sc)
            finite=engine.finite_N_functions(g,roots,loop,self.N)
            HN,local,mem=transport_geometry(g,finite['five_signed_density_increments'],correction,x,domain[0],offset,jac,name,5)
            b=Algebra(g,5);modified={key:b.trim(q) for key,q in fields.items()}
            modified.update(E=finite['E_N'],V=finite['V_N'])
            for key in RATES:modified[key]=b.add(b.trim(fields[key]),HN[key])
            modified['pressure']=b.add(modified['P0'],modified['p'])
            chart=dict(source_chart=bg,coordinate=x,domain=list(domain),radius=R,physical_measure=jac,
                actual_phase=phase,roots=roots,loop=loop,finite_N=finite,incoming=correction,
                history=HN,contributions=local,memory=mem,modified_fields=modified)
            D=radial.NativeDerivative(g,x,phase);D.shared['family']=self.family
            D.shared['free'].update({i:set(v) for i,v in self.leading.free_source_seeds.items()})
            bg_radial=radial.background_radial_rows(g,bg)
            for key,q in fields.items():
                for ref,dr in zip(rows(q),rows(bg_radial[key])):D.seed(ref,g.mul(jac,dr))
            self.derivatives[name]=D
            total=radial.total_radial_functions(g,chart,D,self.N)
            modified_generic=radial.generic_functions(g,modified,self.delta,self.Pstar,R)
            physical=radial.physical_pressure_fields(g,modified,total,self.delta,self.Pstar)
            S=radial.reduced_pressure_S(g,physical,self.delta)
            boundary,pressure_local,cumulative,pressure=pressure_geometry(g,S,physical,R,x,domain[0],jac,integrals,self.initial,name)
            chart.update(radial_background=bg_radial,radial=total,generic=modified_generic,physical=physical,S=S,
                boundary=boundary,incoming_integrals=integrals,local_integrals=pressure_local,cumulative_integrals=cumulative,pressure=pressure)
            self.charts[name]=chart;self.endpoints[name]=self.at(name,domain[1])
            correction=self.endpoints[name]['correction'];integrals=self.endpoints[name]['integrals']
        self.Rc=self.endpoints['O3_power'];self.Rc_background=self.background_endpoints['O3_power']
        self.definition=dict(all_six_true_prescribed_background_C6_and_full_generic_C5=True,
            leading_Rh_and_actual_selected_N_correction_Rh_sources_separate=True,
            all_background_and_finite_N_predecessors_retained=True,
            exact_original_mu_and_exponential_O2_log_radius_Jacobian=True,
            all17_actual_loop_density_history_modified_C5_sources=True,
            total_radial_C5_modified_signed_generic_C4_and_physical_Q_Qy_C4=True,
            signed_rminus_to_Rc_pressure_C3_and_original_global_boundary=True,
            O3_power_quiet_local_density_does_not_reset_correction_memory=True,
            original_six_chart_terminal_is_Rc_not_Rp=True,
            compact_Rc_to_2Rc_controls_and_quiet_2Rc_to_Rp_remain_open=True,
            original_scale_C6_C5_C3_mixed_bounds_and_global_frequency_not_admitted=True,
            full_preRp_global_n1_heat_stress_genuine_recursion_and_NS_open=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.patch,self.leading,self.g,self.initial_correction,self.initial_integrals);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or receipt['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['true_outer_C6_C5_bindings']!=self.snapshot or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked six true outer C6/C5 and Rc pressure C3 required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed outer source '+name)
            self.acceptance_loaded=True
    def background_at(self,name,value):
        bg=self.background[name];x=bg['phase'];a=Algebra(self.g,6)
        fields={key:a.at(q,x,value) for key,q in bg['fields'].items()}
        return dict(fields=fields,history={key:fields[key] for key in RATES},
            radius=substitute(self.g,bg['radius'],x,value),radial={key:a.at(q,x,value) for key,q in bg['radial'].items()},
            generic={key:Algebra(self.g,5).at(q,x,value) for key,q in bg['generic'].items()})
    def at(self,name,value):
        chart=self.charts[name];x=chart['coordinate'];group=lambda data,order:{key:Algebra(self.g,order).at(q,x,value) for key,q in data.items()}
        return dict(correction=group(chart['history'],5),modified_fields=group(chart['modified_fields'],5),
            generic=group(chart['generic'],4),radial={key:Algebra(self.g,5).at(chart['radial'][key],x,value)
                for key in ('E_y','V_y','m_y','slow_A_y','slow_B_y')},
            physical={key:Algebra(self.g,4 if key in ('Q','Q_y') else 5).at(q,x,value) for key,q in chart['physical'].items()},
            integrals=group(chart['cumulative_integrals'],3),pressure=group(chart['pressure'],3),
            boundary=group(chart['boundary'],3),radius=substitute(self.g,chart['radius'],x,value))
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed outer input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        charts={}
        for name,before in self.charts.items():
            chart={k:v for k,v in before.items() if k!='source_chart'};loop=dict(chart['loop'])
            loop['phase_partial_coefficients']=[dict(Z=k,angle=j,node=v.node) for (k,j),v in loop['phase_partial_coefficients'].items()]
            charts[name]=dict(chart,loop=loop)
        return encode(dict(background=self.background,background_endpoints=self.background_endpoints,charts=charts,endpoints=self.endpoints,
            initial_background=self.initial_background,initial_correction=self.initial_correction,initial_integrals=self.initial_integrals,
            initial=self.initial,N=self.N,eta=self.eta,dstar=self.dstar,sc=self.sc,phase_origin=self.phase_origin,
            logmu=self.logmu,mu=self.mu,Tw=self.Tw,mu_proof=self.mu_proof,seed_identity=self.seed_identity,
            definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.patch.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.patch,self.leading,self.g,self.initial_correction,self.initial_integrals)
            or self.initial_correction is not self.patch.Rh_correction or self.initial_integrals is not self.patch.Rh['integrals']
            or self.initial is not self.patch.initial):raise ValueError('Actual outer source or Rh predecessor changed')
        bg,corr,pressure=self.initial_background,self.initial_correction,self.initial_integrals
        for name in CHARTS:
            if (self.background[name]['incoming'] is not bg or self.charts[name]['incoming'] is not corr
                or self.charts[name]['incoming_integrals'] is not pressure or self.charts[name]['source_chart'] is not self.background[name]):
                raise ValueError('An original outer predecessor was reset')
            bg=self.background_endpoints[name]['history'];corr=self.endpoints[name]['correction'];pressure=self.endpoints[name]['integrals']


def run(patch):
    began=time.monotonic();owner=CurrentOriginalTrueOuterC6C5Functions(patch,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,true_outer_C6_C5_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_SIX_OUTER_C6_C5_AND_RC_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(patch) reuses the accepted actual Rh source')
