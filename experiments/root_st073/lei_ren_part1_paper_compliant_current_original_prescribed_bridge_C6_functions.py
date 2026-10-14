"""True prescribed background bridge C6, before finite-N loop corrections.

The known comparison direction remains distinct from actual field histories.
Both micro controls, frozen macro comparison and complete finite endpoints
are original functions. No comparison enclosure defines a function value.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source as inlet
import lei_ren_part1_paper_compliant_actual_bridge_integrals as bridge
import lei_ren_part1_paper_compliant_bridge_mixed_C4 as comparison
import lei_ren_part1_paper_compliant_inner_bridge_profiles as profiles

source,parent,controls=inlet.source,inlet.parent,inlet.controls
HERE,PREFIX,sha=inlet.HERE,inlet.PREFIX,inlet.sha
NAME=PREFIX+'current_original_prescribed_bridge_C6_functions.json.gz'
RECEIPT=PREFIX+'current_original_prescribed_bridge_C6_functions_check.json'
GATES=('current_original_true_comparison_bridge_C7_functions_installed',
    'current_original_true_prescribed_micro_macro_C6_functions_installed',
    'current_original_true_prescribed_bridge_Q_inertial_p2_C5_functions_installed')
OPEN=tuple(dict.fromkeys(inlet.OPEN+('current_original_actual_modified_bridge_C6_functions_installed',
    'current_original_prescribed_bridge_C6_numeric_containment_installed')))
FunctionJet,rows,encode=inlet.FunctionJet,inlet.rows,inlet.encode
RATES=bridge.RATES


class BridgeAlgebra(inlet.JetAlgebra):
    def exp(self,jet):
        q=rows(self.trim(jet));g=self.g;out=[g.unary('exp',q[0])]
        for n in range(self.order):
            out.append(g.add(*(g.mul(g.constant(math.comb(n,j)),q[j+1],out[n-j]) for j in range(n+1))))
        return FunctionJet(tuple(out))
    def log(self,jet,certificate):
        q=rows(self.trim(jet));a=BridgeAlgebra(self.g,self.order-1)
        quotient=a.div(a.shift(jet),a.trim(jet),certificate)
        return FunctionJet((self.g.unary('log',q[0]),)+rows(quotient))


def extend(derivative,jet,order):
    q=list(rows(jet))
    for left,right in zip(q,q[1:]):derivative.seed(left,right)
    while len(q)<=order:q.append(derivative.derivative(q[-1]))
    if len(q)!=order+1:raise ValueError('Exact ordinary source order required')
    return FunctionJet(tuple(q))


def densities(a,phi,U):
    return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),
        A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))


def histories(a,phi,U,initial,phase,jacobian,label):
    """Six complete Duhamel functions, with nonzero original incoming memory."""
    g=a.g;t=g.symbol(label+'_history_t');bodies=densities(a,phi,U);out={}
    for name,rate in RATES.items():
        decay=g.unary('exp',g.mul(g.constant(-rate),jacobian,phase))
        kernel=g.unary('exp',g.mul(g.constant(-rate),jacobian,g.sub(phase,t)))
        body=a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))
        out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,
            actual_six_history_rate=rate,actual_history_name=name,
            original_incoming_memory_retained=True,source_measure_applied_once=True,
            prescribed_bridge_chart=label))
    return out


def comparison_directions(g,phi,U,J,F0,P0,delta):
    """Exact (9.13) functions; sixth direction consumes seventh mean rows."""
    a=BridgeAlgebra(g,6);b=BridgeAlgebra(g,7);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5)
    one=a.fixed(1);de=a.fixed(delta);d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    M,A=(a.trim(J[name]) for name in ('M','A'));Mz=a.shift(J['M']);Az=a.shift(J['A'])
    f,h,k=(b.mul(F0,jet) for jet in (phi,J['H'],J['K']))
    F2=b.mul(F0,F0);B,C=(b.mul(F2,J[name]) for name in ('B','C'))
    angular=a.add(a.mul(a.sub(one,a.scale(de,'1/2')),a.trim(h)),
        a.neg(a.scale(a.mul(a.sub(one,de),a.mul(z,a.shift(h))),'1/2')),
        a.neg(a.mul(d,a.shift(k))),a.mul(a.sub(a.scale(de,2),one),a.mul(z,a.trim(k))))
    W=a.sub(one,a.add(a.mul(a.sub(one,de),a.mul(z,M)),a.mul(d,Mz)))
    halfL=a.scale(L,2);divide=lambda jet:a.div(jet,halfL,'same positive 2*(1-delta*Z^2)')
    D=a.div(a.add(a.neg(W),a.div(angular,a.scale(a.trim(f),2),
        'same strictly positive full comparison F0*phi')),L,'same positive 1-delta*Z^2')
    hydro=divide(a.add(a.neg(a.mul(W,a.trim(U))),
        a.scale(a.mul(a.sub(one,de),a.sub(M,a.mul(z,Mz))),'1/2'),
        a.scale(a.mul(de,a.mul(z,A)),2),a.neg(a.mul(d,Az))))
    pressure=divide(a.add(a.scale(a.mul(a.add(one,de),a.mul(z,a.trim(P0))),2),a.neg(a.mul(d,a.shift(P0)))))
    dressed_swirl=divide(a.add(a.neg(a.scale(a.mul(de,a.mul(z,a.trim(B))),2)),
        a.mul(d,a.shift(B)),a.scale(a.mul(a.add(one,de),a.mul(z,a.trim(C))),2),a.neg(a.mul(d,a.shift(C)))))
    # The external physical drive multiplies F0^2 once. These derivatives
    # are dressed first, then normalized back by the SAME full F0^2.
    swirl=a.div(dressed_swirl,a.trim(F2),'same strictly positive full F0 squared')
    return dict(D_over_R=D,drive_hydro=hydro,drive_pressure=pressure,drive_swirl=swirl)


def actual_fields(a,known,initial,phase,jacobian,chi,radius,F0,Pstar,barphi,label):
    g=a.g;t=g.symbol(label+'_actual_log_t')
    angular=a.scale(known['D_over_R'],g.mul(chi,radius,jacobian,g.constant('-1/2')))
    logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase,
        actual_prescribed_angular_equation=True,both_original_micro_controls_retained=True,
        prescribed_bridge_chart=label)
    phi=a.mul(initial['phi'],a.exp(logratio))
    full_drive=a.add(a.scale(known['drive_hydro'],radius),
        a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar)),
        a.mul(a.scale(a.mul(F0,F0),g.mul(radius,radius)),known['drive_swirl']))
    ratio=a.div(phi,a.trim(barphi),'same positive actual comparison phi')
    Vy=a.scale(a.mul(ratio,full_drive),g.neg(chi))
    v=g.symbol(label+'_actual_U_t')
    U=a.sub(initial['U'],a.integral(a.scale(a.at(Vy,phase,v),g.neg(jacobian)),v,g.zero,phase,
        actual_prescribed_axial_equation=True,all_hydro_pressure_swirl_sectors_retained=True,
        source_measure_applied_once=True,prescribed_bridge_chart=label))
    phiy=a.scale(a.mul(known['D_over_R'],phi),g.mul(chi,radius,g.constant('-1/2')))
    return dict(phi=phi,U=U,phi_y=phiy,U_y=Vy,logratio=logratio,ratio=ratio,
        angular_density=angular,full_signed_axial_drive=full_drive)


def generic_fields(a,actual,J,F0,P0,Pstar,radius,delta):
    g=a.g;c=g.constant;inverse=g.quotient(g.one,Pstar,'same strictly positive Pstar')
    W=a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))
    FF=a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))
    fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse),
        h=a.scale(a.mul(W,J['H']),'1/2'),k=a.scale(a.mul(W,J['K']),g.mul(c('1/2'),inverse)),
        e=a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B'])),p=a.mul(FF,J['C']),P0=a.trim(P0))
    fields['pressure']=a.add(fields['P0'],fields['p'])
    radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))
    return fields,radial,inlet.generic_functions(g,fields,delta,Pstar,radius)


class CurrentOriginalPrescribedBridgeC6Functions:
    def __init__(self,inlet_owner,require_checked=True):
        if type(inlet_owner) is not inlet.CurrentOriginalCoreInletAxial6GenericSource or not inlet_owner.acceptance_loaded:
            raise TypeError('Checked true core inlet C6 source required')
        inlet_owner.assert_graph();self.inlet=inlet_owner;self.core=inlet_owner.core;self.g=g=inlet_owner.g
        self.family=copy.deepcopy(inlet_owner.family);self.hashes=dict(inlet_owner.hashes)
        for name,digest in parent.linear.read(inlet.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (inlet.NAME,inlet.RECEIPT):self._pin(name,sha(name))
        for module in (bridge,comparison,profiles):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        self.bindings=bridge.bridge_control_bindings()
        self.derivative=d=inlet.previous.HigherCoreZDerivative(inlet_owner.previous.lower)
        self.core7={name:extend(d,jet,7) for name,jet in inlet_owner.core_functions.items()}
        self.radial7={name:[extend(d,jet,7) for jet in jets]
            for name,jets in inlet_owner.previous.core_radial_rows.items()}
        self.initial7={name:extend(d,jet,7) for name,jet in inlet_owner.inlet_moments.items()}
        self.P0=self.pressure7();a=BridgeAlgebra(g,7);b=BridgeAlgebra(g,6);c=g.constant
        self.hb=hb=inlet_owner.previous.lower.parameters['hb'];self.Ra=Ra=inlet_owner.radius
        self.Y=g.sub(g.unary('log',c(100)),g.unary('log',Ra));self.length=length=g.sub(self.Y,g.mul(c(2),hb))
        s=g.symbol('original_prescribed_C6_micro_s');u=g.symbol('original_prescribed_C6_comparison_t')
        sigma=lambda q:inlet.previous.lower.own.jets.complete.Gamma.sigma(g,q)
        rho=lambda q:g.mul(c(4),g.unary('exp',g.mul(hb,q)))
        atcore=lambda jet,q:a.at(jet,self.core.rho,rho(q))
        phi0=a.at(self.core7['Phi'],self.core.rho,c(4))
        U0=a.add(self.core7['U0'],a.scale(a.at(self.core7['Psi'],self.core.rho,c(4)),self.core.parameters['epsilon']))
        alpha=g.sub(g.one,sigma(g.sub(u,g.one)))
        ratio=a.div(atcore(self.radial7['Phi'][1],u),atcore(self.core7['Phi'],u),
            'same positive unique analytic Phi on rho4..4exp(2hb)')
        logbar=a.add(a.log(phi0,'same positive core inlet Phi'),a.integral(a.scale(ratio,g.mul(alpha,hb,rho(u))),
            u,g.zero,s,exact_original_comparison_smoothing=True,comparison_not_actual_field=True))
        barphi=a.exp(logbar)
        barU=a.add(U0,a.integral(a.scale(atcore(self.radial7['Psi'][1],u),
            g.mul(alpha,hb,rho(u),self.core.parameters['epsilon'])),u,g.zero,s,
            exact_original_comparison_smoothing=True,comparison_not_actual_field=True))
        barJ=histories(a,barphi,barU,self.initial7,s,hb,'comparison_micro')
        known=comparison_directions(g,barphi,barU,barJ,self.core7['F0'],self.P0,self.core.parameters['delta'])
        self.comparison_micro=dict(phase=s,phi=barphi,U=barU,logphi=logbar,moments=barJ,directions=known)
        chi=g.add(g.one,g.mul(g.sub(hb,g.one),sigma(s)))
        R=g.mul(Ra,g.unary('exp',g.mul(hb,s)));F0=b.trim(self.core7['F0']);Pstar=inlet_owner.Pstar
        actual=actual_fields(b,known,dict(phi=b.trim(phi0),U=b.trim(U0)),s,hb,chi,R,F0,Pstar,barphi,'micro')
        actualJ=histories(b,actual['phi'],actual['U'],{k:b.trim(v) for k,v in self.initial7.items()},s,hb,'actual_micro')
        fields,radial,generic=generic_fields(b,actual,actualJ,F0,self.P0,Pstar,R,self.core.parameters['delta'])
        self.micro=dict(phase=s,domain=[0,2],subcharts=dict(first=[0,1],second=[1,2]),radius=R,
            physical_measure=hb,chi=chi,actual=actual,moments=actualJ,fields=fields,radial=radial,generic=generic)
        self.actual_micro_end={name:b.at(actual[name],s,c(2)) for name in ('phi','U')}
        self.actual_moment_end={name:b.at(jet,s,c(2)) for name,jet in actualJ.items()}
        self.comparison_micro_end=dict(phi=a.at(barphi,s,c(2)),U=a.at(barU,s,c(2)),
            moments={name:a.at(jet,s,c(2)) for name,jet in barJ.items()})
        m=g.symbol('original_prescribed_C6_macro_fraction');theta=g.unary('exp',g.neg(g.mul(length,m)))
        barphi2,barU2=(self.comparison_micro_end[name] for name in ('phi','U'))
        equilibrium=dict(H=barphi2,M=barU2,K=a.mul(barphi2,barU2),A=a.mul(barU2,barU2),
            B=a.scale(a.mul(barphi2,barphi2),'1/2'),C=a.mul(barphi2,barphi2))
        macroJ={}
        for name,rate in RATES.items():
            decay=g.unary('exp',g.mul(c(-rate),length,m))
            macroJ[name]=a.add(a.scale(self.comparison_micro_end['moments'][name],decay),
                a.scale(equilibrium[name],g.sub(g.one,decay)))
        macroknown=comparison_directions(g,barphi2,barU2,macroJ,self.core7['F0'],self.P0,self.core.parameters['delta'])
        self.comparison_macro=dict(phase=m,phi=barphi2,U=barU2,moments=macroJ,directions=macroknown,
            frozen_only_the_comparison_profile=True,theta=theta)
        Rmacro=g.mul(Ra,g.unary('exp',g.add(g.mul(c(2),hb),g.mul(length,m))))
        actualmacro=actual_fields(b,macroknown,self.actual_micro_end,m,length,hb,Rmacro,F0,Pstar,barphi2,'macro')
        actualmacroJ=histories(b,actualmacro['phi'],actualmacro['U'],self.actual_moment_end,m,length,'actual_macro')
        macrofields,macroradial,macrogeneric=generic_fields(b,actualmacro,actualmacroJ,F0,self.P0,Pstar,Rmacro,self.core.parameters['delta'])
        self.macro=dict(phase=m,domain=[0,1],radius=Rmacro,physical_measure=length,chi=hb,
            actual=actualmacro,moments=actualmacroJ,fields=macrofields,radial=macroradial,generic=macrogeneric)
        self.R100=dict(radius=c(100),actual={k:b.at(v,m,g.one) for k,v in actualmacro.items()},
            moments={k:b.at(v,m,g.one) for k,v in actualmacroJ.items()},
            fields={k:b.at(v,m,g.one) for k,v in macrofields.items()},
            generic={k:BridgeAlgebra(g,5).at(v,m,g.one) for k,v in macrogeneric.items()})
        self.definition=dict(true_prescribed_background_before_finite_N_modification=True,
            actual_twenty_term_core_and_own_six_histories_retained=True,
            comparison_C7_before_direction_C6_before_actual_C6=True,
            known_direction_not_replaced_by_actual_moments=True,
            both_original_micro_chi_controls_and_macro_hb_retained=True,
            exact_micro_dy_hb_ds_and_macro_dy_length_dm_once=True,
            F0_derivatives_dressed_before_swirl_normalization_and_scale_once=True,
            nonzero_all_six_incoming_histories_and_micro_macro_memory_retained=True,
            finite_N_modified_velocity_and_correction_histories_not_installed=True,
            actual_C6_numeric_enclosures_and_global_mixed_norms_not_certified=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.inlet,self.core,self.g);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES)
                or any(receipt.get(k) for k in OPEN) or receipt['source_family']!=self.family
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or archive!=sha(NAME)
                or checker!=sha(Path(__file__).stem+'_check.py') or raw['prescribed_bridge_C6_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual prescribed bridge C6 functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed prescribed bridge source '+name)
            self.acceptance_loaded=True

    def pressure7(self):
        q=list(rows(self.inlet.P0));attrs=copy.deepcopy(self.g.nodes[q[6].node]);operation=attrs.pop('operation')
        attrs.update(Z_order=7,ordinary_slow_Z_derivative_order=7,Taylor_coefficient_factorial=5040,
            source_projection_binding=parent.active.current.ast_binding(self.project_P0_seventh),
            pressure_provider_binding=parent.active.current.ast_binding(self.P0_rows),
            source_projection='factorial(7)*P0_rows(Z).coefficients[7]',true_order7_pressure_provider=True)
        q.append(self.g.node(operation,**attrs));return FunctionJet(tuple(q))
    def P0_rows(self,Z):
        self.inlet.hierarchy.assert_graph();provider=self.inlet.hierarchy.source.admission.pressure
        packet=provider.normalized_jets(Z,7)
        if not packet['all_14_pressure_atoms_included'] or len(packet['normalized_pressure_coefficients'])!=8:
            raise ValueError('Same full eight-row independent P0 source required')
        return packet['normalized_pressure_coefficients']
    def project_P0_seventh(self,ref,Z):
        self.assert_graph()
        if ref!=self.P0.rows[7]:raise ValueError('Exact registered seventh pressure derivative required')
        return self.P0_rows(Z)[7]*5040
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed prescribed bridge input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(core7=self.core7,radial7=self.radial7,initial7=self.initial7,P0=self.P0,
            comparison_micro=self.comparison_micro,comparison_macro=self.comparison_macro,
            micro=self.micro,macro=self.macro,R100=self.R100,hb=self.hb,Ra=self.Ra,Y=self.Y,length=self.length,
            actual_micro_end=self.actual_micro_end,actual_moment_end=self.actual_moment_end,
            comparison_micro_end=self.comparison_micro_end,definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.inlet.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.inlet,self.core,self.g)):
            raise ValueError('True prescribed bridge source or histories changed')


def run(inlet_owner):
    began=time.monotonic();owner=CurrentOriginalPrescribedBridgeC6Functions(inlet_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,prescribed_bridge_C6_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_PRESCRIBED_MICRO_MACRO_C6_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':raise SystemExit('run(inlet_owner) reuses the admitted actual core')
