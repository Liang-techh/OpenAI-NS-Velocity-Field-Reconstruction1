"""Actual Rc targets, same compact implicit C5 branch, and pressure to 2Rc.

The selected-N incoming correction is already physical and is not divided
by N again. Exact original bump integrals define every weight. Function
identities and original C1 invertibility do not admit new C5 numeric norms.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as previous
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as compact
import lei_ren_part1_paper_compliant_current_native_repair_band_cone_budget as original_band
import lei_ren_part1_paper_compliant_current_original_whole_Z_limit_repair_band as original_limit
import lei_ren_part1_paper_compliant_current_original_C3_source_targets as original_targets
import lei_ren_part1_paper_compliant_current_original_C3_limit_controls as original_controls

source,parent,controls=previous.source,previous.parent,previous.controls
engine,radial=previous.engine,previous.radial
leading=previous.previous.previous
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=engine.LoopAlgebra
RATES=engine.RATES
NAME=PREFIX+'current_original_true_compact_C5_functions.json.gz'
RECEIPT=PREFIX+'current_original_true_compact_C5_functions_check.json'
GATES=('current_original_true_all17_normalized_compact_targets_C5_functions_installed',
       'current_original_true_same_branch_compact_controls_C5_functions_installed',
       'current_original_true_Rc_to_2Rc_compact_modified_history_C5_functions_installed',
       'current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed',
       'current_original_true_rminus_to_2Rc_signed_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_true_compact_C5_numeric_containment_installed',
    'current_original_true_2Rc_to_Rp_quiet_C5_pressure_C3_functions_installed')))


def normalized_targets(g,correction,A,mu,N):
    """d=N*r, from actual selected-N corrections; retain the joint numerator."""
    a=Algebra(g,5);A=a.trim(A);AA=a.mul(A,A)
    scaled={k:a.scale(v,N) for k,v in correction.items()}
    d={row:a.div(scaled[k],A if degree==1 else AA,'original positive Rc amplitude')
        for row,k,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2))}
    joint=a.sub(scaled['k'],a.mul(A,scaled['m']))
    d[compact.ROWS[1]]=a.div(joint,a.scale(AA,mu),'original positive mu times Rc amplitude squared')
    return {key:d[key] for key in compact.ROWS},scaled,joint


def nonlinear(a,W,h):
    zero=a.fixed(0);cross=a.add(a.scale(a.mul(h[0],h[2]),W['cross'][0]),a.scale(a.mul(h[1],h[4]),W['cross'][1]))
    axial=a.add(a.scale(a.mul(h[0],h[0]),W['energy'][0]),a.scale(a.mul(h[1],h[1]),W['energy'][2]))
    swirl=a.add(*(a.scale(a.mul(h[i+2],h[i+2]),W['energy'][i]) for i in range(3)))
    pressure=a.add(*(a.scale(a.mul(h[i+2],h[i+2]),W['pressure'][i]) for i in range(3)))
    return [zero,a.scale(cross,a.g.quotient(a.g.one,W['mu'],'same original positive mu')),zero,
        a.sub(axial,a.scale(swirl,'1/2')),a.scale(pressure,'1/2')]


def implicit_rows(g,W,N,d,values,order=5):
    """Ordinary derivatives of B*h+d+Q(h)/N, with one full Jacobian."""
    invN=g.quotient(g.one,N,'same original exact positive integer N')
    columns=[compact.quadratic_action(g,W,values,[g.one if i==j else g.zero for i in range(5)]) for j in range(5)]
    J=[[g.add(W['B'][i][j],g.mul(invN,columns[j][i])) for j in range(5)] for i in range(5)]
    inverse,det=leading.inverse_matrix(g,J);out=[[v] for v in values];recurrence=[]
    for n in range(1,order+1):
        a=Algebra(g,n);known=[FunctionJet(tuple(q)+(g.zero,)) for q in out]
        Q=nonlinear(a,W,known)
        rhs=[g.add(rows(d[i])[n],g.mul(invN,rows(Q[i])[n])) for i in range(5)]
        solved=[g.neg(g.add(*(g.mul(inverse[i][j],rhs[j]) for j in range(5)))) for i in range(5)]
        for i in range(5):out[i].append(solved[i])
        recurrence.append(dict(order=n,known_rhs=rhs,solved_ordinary_row=solved,
            only_unknown_highest_row_temporarily_zero=True,same_full_Jacobian=True))
    return [FunctionJet(tuple(q)) for q in out],dict(Jacobian=J,inverse=inverse,determinant=det,recurrence=recurrence)


def bump_rows(g,W,x):
    """Actual compact beta, beta_y and beta_yy, y=log(x)."""
    out=[];logx=g.unary('log',x)
    invell=g.quotient(g.one,W['ell'],'original positive log bump width')
    scale=g.quotient(g.one,g.mul(W['ell'],W['normal'],x),'original positive ell J0 x')
    for center in W['centers']:
        argument=g.mul(g.sub(logx,center),invell);raw=leading.raw_beta(g,argument)
        first,second=leading.beta_first_two(g,argument,raw)
        b=g.mul(scale,raw)
        by=g.mul(scale,g.sub(g.mul(first,invell),raw))
        byy=g.mul(scale,g.add(g.mul(second,invell,invell),g.mul(g.constant(-2),first,invell),raw))
        out.append((b,by,byy))
    return out


def profiles(g,W,N,A,h,x):
    a=Algebra(g,5);bumps=bump_rows(g,W,x);invN=g.quotient(g.one,N,'same original exact positive integer N')
    alpha=g.add(g.constant('1/2'),W['mu']);power=g.unary('exp',g.neg(g.mul(alpha,g.unary('log',x))))
    radial=[]
    for j in range(3):
        F=a.scale(a.add(*(a.scale(h[i+2],bumps[i][j]) for i in range(3))),invN)
        G=a.scale(a.add(a.scale(h[0],bumps[0][j]),a.scale(h[1],bumps[2][j])),invN)
        fixed=g.mul(power,g.constant((-1)**j),*([alpha]*j))
        E=a.scale(a.trim(A),fixed);dE=a.mul(a.trim(A),F);dV=a.mul(a.trim(A),G)
        radial.append(dict(F=F,G=G,original_E=E,delta_E=dE,delta_V=dV,E=a.add(E,dE),V=dV))
    return dict(bumps=bumps,radial_y_rows=radial,alpha=alpha,power=power)


def density(g,E,dE,dV):
    a=Algebra(g,5);EF=a.mul(E,dE);FF=a.mul(dE,dE);VV=a.mul(dV,dV)
    return dict(m=dV,h=dE,k=a.mul(a.add(E,dE),dV),
        e=a.sub(a.sub(VV,EF),a.scale(FF,'1/2')),p=a.add(EF,a.scale(FF,'1/2')))


class CurrentOriginalTrueCompactC5Functions:
    def __init__(self,outer,require_checked=True):
        if type(outer) is not previous.CurrentOriginalTrueOuterC6C5Functions or not outer.acceptance_loaded:
            raise TypeError('Checked actual all17 Rc source required')
        outer.assert_graph();self.outer=outer;self.g=g=outer.g;self.family=copy.deepcopy(outer.family);self.hashes=dict(outer.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for module in (compact,original_band,original_limit,original_targets,original_controls):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        self.N,self.mu,self.Pstar,self.P0,self.delta=outer.N,outer.mu,outer.Pstar,outer.P0,outer.delta
        self.Rc=outer.Rc['radius'];self.A=outer.Rc_background['fields']['E']
        self.incoming_background=outer.Rc_background['history'];self.incoming_correction=outer.Rc['correction']
        self.incoming_integrals=outer.Rc['integrals'];self.initial=outer.initial
        self.bindings=dict(original_weights=parent.active.current.ast_binding(compact.exact_weights),
            original_quadratic=parent.active.current.ast_binding(compact.quadratic),
            original_quadratic_action=parent.active.current.ast_binding(compact.quadratic_action),
            original_inverse_action=parent.active.current.ast_binding(compact.inverse_action),
            original_target_transport=parent.active.current.ast_binding(original_targets.build_transport),
            original_N_scaled_density=parent.active.current.ast_binding(original_targets.coefficient_C3),
            original_profile=parent.active.current.ast_binding(original_band.profiles_at),
            original_signed_density=parent.active.current.ast_binding(original_band.density_pairs),
            original_partial_history_terminal_identity=parent.active.current.ast_binding(original_band.exact_partial_band_graph),
            original_limit=parent.active.current.ast_binding(original_limit.exact_limit_adapter),
            original_C3_branch=original_controls.source_bindings(),
            original_all17_profiles=outer.bindings,original_Rh_seed=outer.seed_identity)
        admission=parent.linear.read(previous.original_geometry.NAME)
        receipt=parent.linear.read(previous.original_geometry.RECEIPT)
        proof=admission['actual_whole_Z_frequency_connection']
        c3_receipt=parent.linear.read(original_controls.RECEIPT)
        old_targets=parent.linear.read(original_targets.NAME)
        if (not receipt['all_passed'] or not receipt[previous.original_geometry.GATE]
            or admission['source_family']!=self.family or c3_receipt['source_family']!=self.family
            or not c3_receipt['all_passed'] or not all(c3_receipt[k] for k in original_controls.GATES)
            or proof['selected_integer']['J']!=g.nodes[self.N.node]['J']
            or not proof['C1_limit']['exact_C1_limit_exists_unique_for_this_repair_only_integer']
            or proof['C1_limit']['contraction_at_most']!='1/2'):
            raise ValueError('Same original selected-N unique compact branch required')
        if old_targets['source_family']!=self.family:raise ValueError('Same original target-function source required')
        for rec in (receipt,c3_receipt):
            for name,digest in rec['input_hashes'].items():self._pin(name,digest)
        for name in (previous.original_geometry.NAME,previous.original_geometry.RECEIPT,original_controls.RECEIPT):self._pin(name,sha(name))
        for name in (original_targets.NAME,original_targets.RECEIPT):self._pin(name,sha(name))
        old=admission['exact_current_integral_and_control_graph']
        self.original_prefix=dict(target_archive=original_targets.NAME,
            actual_C3_target_functions=old_targets['actual_C3_density_transport_targets']['actual_N_scaled_five_targets'],
            actual_C3_N_scaled_incoming=old_targets['actual_C3_density_transport_targets']['terminal_N_scaled_histories'],
            actual_C3_amplitude=old_targets['actual_C3_density_transport_targets']['actual_terminal_amplitude'],
            original_mu_node=old['parameters']['mu'],original_N_node=old['N'],
            original_parameter_archive=previous.original_geometry.NAME,
            prefix_orders=[0,1,2,3],normalization_expressions_not_numeric_ranges=True)
        self.W=compact.exact_weights(g,self.mu)
        self.targets,self.N_scaled_incoming,self.joint=normalized_targets(g,self.incoming_correction,self.A,self.mu,self.N)
        self.branch_transfer=dict(source_family=self.family,selected_integer=proof['selected_integer'],
            original_C1_limit_proof=proof['C1_limit'],original_admission_receipt=previous.original_geometry.RECEIPT,
            exact_source_correspondence=dict(same_all17_prescribed_source_and_normalized_profile_functions=True,
                original_N_global_phase_eta_dstar_and_initial_zero_correction_retained=True,
                source_density_identity='N*delta=(B,F_N,V*F_N+E*B,2V*B-E*F_N,E*F_N)+(0,0,F_N*B,B^2-F_N^2/2,F_N^2/2)/N; F_N=N*E*expm1(A_loop/N)',
                physical_selected_N_incoming_not_divided_again=True,
                own_rate_original_Duhamel_maps_and_inlet_bound_by_predecessor=True,
                original_positive_Rc_amplitude_and_joint_mu_target_normalization=True),
            same_original_B_Q_N_and_exact_bump_integrals=True,
            invertibility='J=B+DQ(h*)/N; norm(B^-1*DQ/N)<=1/2 on the original branch',
            C5_extension='Actual targets are C5; the same invertible J differentiates the existing unique branch by IFT.',
            original_C1_contraction_not_new_C5_numeric_norm_admission=True)
        generic=[g.symbol('true_compact_branch_h%d'%i) for i in range(5)]
        Q=compact.quadratic(g,self.W,generic)
        vec=[g.add(self.targets[key].value,g.quotient(q,self.N,'same original exact positive integer N')) for key,q in zip(compact.ROWS,Q)]
        control_map=[g.neg(v) for v in compact.inverse_action(g,self.W,vec)]
        self.system=g.node('original_compact_unique_Banach_limit_system',variables=[q.node for q in generic],
            exact_map=[q.node for q in control_map],actual_target_values=[q.value.node for q in self.targets.values()],
            actual_selected_N=self.N.node,original_weights=encode(self.W),initial_vector=[g.zero.node]*5,
            source_family=self.family,original_branch_transfer=self.branch_transfer,
            definition='same original unique C1 Banach limit; no finite iterate or enclosure chosen as its value')
        values=[g.node('original_compact_unique_Banach_limit_component',system=self.system.node,component=i,
            Z_variable=g.symbol('Z').node,source_family=self.family,mathematical_function_not_numeric_interval=True) for i in range(5)]
        self.h,self.implicit=implicit_rows(g,self.W,self.N,list(self.targets.values()),values)
        a=Algebra(g,5);Q=nonlinear(a,self.W,self.h);invN=g.quotient(g.one,self.N,'same original exact positive integer N')
        self.residual={key:a.add(self.targets[key],a.add(*(a.scale(q,b) for q,b in zip(self.h,row))),a.scale(qq,invN))
            for key,row,qq in zip(compact.ROWS,self.W['B'],Q)}
        self.free_source_seeds={i:set(v) for i,v in outer.leading.free_source_seeds.items()}
        for q in values:self.free_source_seeds[q.node]={'Z'}
        for row in g.nodes:
            if row['operation']=='compact_raw_beta':
                # All exact-weight raw bumps have this normalized scalar argument.
                argument=g.nodes[row['argument']]
                if argument['operation']=='bound_variable':
                    ref=g.node('compact_raw_beta',**{k:v for k,v in row.items() if k!='operation'})
                    self.free_source_seeds[ref.node]={argument['name']}
        self.x=x=g.symbol('true_compact_repair_x');offset=g.unary('log',x);jac=g.quotient(g.one,x,'original positive repair radius x')
        self.profiles=profiles(g,self.W,self.N,self.A,self.h,x)
        lead_alg=Algebra(g,6);E=lead_alg.scale(self.A,self.profiles['power']);V=lead_alg.fixed(0)
        histories,parts,memory=previous.transport_geometry(g,previous.background_density(g,E,V),self.incoming_background,
            x,g.one,offset,jac,'true_compact_background',6)
        fields=dict(E=E,V=V,**histories,P0=lead_alg.trim(self.P0));fields['pressure']=lead_alg.add(fields['P0'],fields['p'])
        R=g.mul(self.Rc,x);Ey=lead_alg.scale(E,g.neg(self.profiles['alpha']))
        self.background=dict(phase=x,domain=[g.one,g.constant(2)],radius=R,fields=fields,radial=dict(E_y=Ey,V_y=V),
            history=histories,incoming=self.incoming_background,history_contributions=parts,memory=memory,
            generic=radial.generic_functions(g,fields,self.delta,self.Pstar,R,5))
        f=self.profiles['radial_y_rows'];self.density=density(g,f[0]['original_E'],f[0]['delta_E'],f[0]['delta_V'])
        self.history,self.contributions,self.memory=previous.transport_geometry(g,self.density,self.incoming_correction,
            x,g.one,offset,jac,'true_compact_correction',5)
        self.modified_fields={k:a.trim(v) for k,v in fields.items()};self.modified_fields.update(E=f[0]['E'],V=f[0]['V'])
        for k in RATES:self.modified_fields[k]=a.add(a.trim(fields[k]),self.history[k])
        self.modified_fields['pressure']=a.add(self.modified_fields['P0'],self.modified_fields['p'])
        self.radial=dict(E_y=f[1]['E'],V_y=f[1]['V'],m_y=a.sub(self.modified_fields['V'],self.modified_fields['m']),
            E_yy=f[2]['E'],V_yy=f[2]['V'])
        self.radial['m_yy']=a.sub(self.radial['V_y'],self.radial['m_y'])
        self.generic=radial.generic_functions(g,self.modified_fields,self.delta,self.Pstar,R)
        self.physical=radial.physical_pressure_fields(g,self.modified_fields,self.radial,self.delta,self.Pstar)
        self.S=radial.reduced_pressure_S(g,self.physical,self.delta)
        self.boundary,self.pressure_local,self.integrals,self.pressure=previous.pressure_geometry(g,self.S,self.physical,R,
            x,g.one,jac,self.incoming_integrals,self.initial,'true_compact')
        self.terminal=self.at(g.constant(2));self.inlet=self.at(g.one)
        AA=a.mul(a.trim(self.A),a.trim(self.A));self.terminal_rhs={}
        for key,row,amplitude in (('m','M',a.trim(self.A)),('h','I',a.trim(self.A)),('e','S',AA),('p','Cp',AA)):
            decay=g.unary('exp',g.mul(g.constant(-RATES[key]),g.unary('log',g.constant(2))))
            self.terminal_rhs[key]=a.scale(a.mul(amplitude,self.residual[row]),g.mul(decay,invN))
        joint=a.add(self.residual['M'],a.scale(self.residual[compact.ROWS[1]],self.mu))
        decay=g.unary('exp',g.mul(g.constant('-3/2'),g.unary('log',g.constant(2))))
        self.terminal_rhs['k']=a.scale(a.mul(AA,joint),g.mul(decay,invN))
        self.terminal_certificates={key:g.node('proved_original_compact_relative_C5_terminal_identity',
            actual_endpoint_rows=encode(self.terminal['correction'][key]),equivalent_residual_rows=encode(self.terminal_rhs[key]),
            source_system=self.system.node,ordinary_implicit_recurrence=encode(self.implicit['recurrence']),
            original_signed_endpoint_binding=self.bindings['original_partial_history_terminal_identity'],
            zero_rows=[g.zero.node]*6,relative_correction_only=True,absolute_leading_histories_and_pressure_not_zeroed=True,
            theorem='original normalized bump integrals give endpoint=residual; same exact branch and ordinary derivatives make residual zero') for key in RATES}
        self.definition=dict(genuine_Rc_normalized_targets_and_same_branch_controls_C5=True,
            actual_selected_N_Rc_correction_used_without_additional_N_division=True,
            exact_original_compact_profiles_and_signed_history_C5=True,
            compact_modified_total_radial_and_second_velocity_rows_C5=True,
            compact_modified_generic_C4_and_signed_pressure_C3=True,
            terminal_relative_five_C5_identity_leaves_absolute_background_and_P0=True,
            original_C1_invertibility_transferred_not_C5_numeric_norms=True,
            quiet_2Rc_to_Rp_full_preRp_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.outer,self.g,self.A,self.incoming_background,self.incoming_correction,self.incoming_integrals);self.acceptance_loaded=False
        if require_checked:
            raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not rec.get('all_passed') or not all(rec.get(k) for k in GATES) or any(rec.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['true_compact_C5_bindings']!=self.snapshot or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual compact C5 and pressure C3 required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed compact source '+name)
            self.acceptance_loaded=True
    def at(self,value):
        group=lambda data,order:{k:Algebra(self.g,order).at(q,self.x,value) for k,q in data.items()}
        return dict(correction=group(self.history,5),modified_fields=group(self.modified_fields,5),
            background_fields=group(self.background['fields'],6),radial=group(self.radial,5),generic=group(self.generic,4),
            integrals=group(self.integrals,3),pressure=group(self.pressure,3),boundary=group(self.boundary,3),
            radius=previous.substitute(self.g,self.background['radius'],self.x,value))
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed compact input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(N=self.N,mu=self.mu,Pstar=self.Pstar,P0=self.P0,delta=self.delta,Rc=self.Rc,A=self.A,
            incoming_background=self.incoming_background,incoming_correction=self.incoming_correction,
            incoming_integrals=self.incoming_integrals,initial=self.initial,targets=self.targets,
            N_scaled_incoming=self.N_scaled_incoming,joint=self.joint,W=self.W,h=self.h,system=self.system,
            implicit=self.implicit,residual=self.residual,branch_transfer=self.branch_transfer,bindings=self.bindings,
            original_prefix=self.original_prefix,
            profiles=self.profiles,background=self.background,density=self.density,history=self.history,
            contributions=self.contributions,memory=self.memory,modified_fields=self.modified_fields,
            radial=self.radial,generic=self.generic,physical=self.physical,S=self.S,boundary=self.boundary,
            pressure_local=self.pressure_local,integrals=self.integrals,pressure=self.pressure,terminal=self.terminal,
            inlet=self.inlet,terminal_rhs=self.terminal_rhs,terminal_certificates=self.terminal_certificates,definition=self.definition))
    def assert_graph(self):
        self.outer.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.outer,self.g,self.A,self.incoming_background,self.incoming_correction,self.incoming_integrals)
            or self.A is not self.outer.Rc_background['fields']['E'] or self.incoming_correction is not self.outer.Rc['correction']
            or self.incoming_background is not self.outer.Rc_background['history'] or self.incoming_integrals is not self.outer.Rc['integrals']
            or self.initial is not self.outer.initial):raise ValueError('Actual compact source or Rc inlet changed')


def run(outer):
    began=time.monotonic();owner=CurrentOriginalTrueCompactC5Functions(outer,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,true_compact_C5_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_COMPACT_C5_AND_2RC_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(outer) reuses the checked actual all17 Rc source')
