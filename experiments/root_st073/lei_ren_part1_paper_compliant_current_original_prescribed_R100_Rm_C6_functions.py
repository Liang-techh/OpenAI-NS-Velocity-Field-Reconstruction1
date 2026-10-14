"""Original prescribed background C6 from R100 through switches to Rm.

The full finite long reshape and all nonzero moment memories remain. This
constructs the background before finite-N loop/correction-history composition.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions as previous
import lei_ren_part1_paper_compliant_current_original_first_switch_functions as first_source
import lei_ren_part1_paper_compliant_current_original_second_switch_R110 as second_source
import lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint as long_source
import lei_ren_part1_paper_compliant_current_original_reference_restore_functions as restore_source
import lei_ren_part1_paper_compliant_inner_switch_profiles as switch_definition
import lei_ren_part1_paper_compliant_long_reshape_profiles as long_definition
import lei_ren_part1_paper_compliant_reference_restore_profiles as restore_definition
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection as connection
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_long_patch_functions as native_long

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_prescribed_R100_Rm_C6_functions.json.gz'
RECEIPT=PREFIX+'current_original_prescribed_R100_Rm_C6_functions_check.json'
GATES=('current_original_true_prescribed_R100_R110_C6_functions_installed',
    'current_original_true_prescribed_finite_long_reference_restore_Rm_C6_functions_installed',
    'current_original_true_prescribed_R100_Rm_Q_inertial_p2_C5_functions_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+('current_original_actual_modified_R100_Rm_C6_functions_installed',
    'current_original_prescribed_R100_Rm_C6_numeric_containment_installed')))
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=previous.BridgeAlgebra
CHARTS=('first_switch','second_switch','post_switch','long_reshape','reference','restoration','post_restore')


def source_bindings():
    return dict(first_switch=first_source.fixed_comparison_denominator_binding(),second_switch=second_source.source_bindings(),
        finite_long=long_source.source_bindings(),reference_restoration=restore_source.source_bindings())


def sigma(g,x):return previous.inlet.previous.lower.own.jets.complete.Gamma.sigma(g,x)


def sigma_prime(g,x):
    own=previous.inlet.previous.lower.own
    return g.node('analytic_unary',name='original_flat_sigma_prime',argument=x.node,
        derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets),
        ordinary_derivative='same original sigma Taylor coefficient[1]',
        exact_flat_source_derivative_not_a_bound=True)


def full_drive(a,known,radius,F0,Pstar):
    g=a.g
    return a.add(a.scale(known['drive_hydro'],radius),
        a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar)),
        a.mul(a.scale(a.mul(F0,F0),g.mul(radius,radius)),known['drive_swirl']))


def first_switch_fields(a,known,initial,phase,hb,radius,F0,Pstar,barphi):
    g=a.g;t=g.symbol('true_first_switch_angular_t');v=g.symbol('true_first_switch_axial_t')
    angular=a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius))
    logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase,
        original_first_switch_full_current_radius_direction=True)
    phi=a.mul(initial['phi'],a.exp(logratio));ratio=a.div(phi,a.trim(barphi),'same positive frozen comparison phi')
    drive=full_drive(a,known,radius,F0,Pstar);chi_axial=g.mul(hb,g.sub(g.one,sigma(g,phase)))
    Uy=a.scale(a.mul(ratio,drive),g.neg(chi_axial))
    U=a.add(initial['U'],a.integral(a.scale(a.at(Uy,phase,v),hb),v,g.zero,phase,
        original_first_switch_axial_control=True,source_measure_applied_once=True))
    phiy=a.scale(a.mul(known['D_over_R'],phi),g.mul(g.constant('-1/2'),hb,radius))
    return dict(phi=phi,U=U,phi_y=phiy,U_y=Uy,logratio=logratio,
        angular_phase_density=angular,full_signed_axial_drive=drive)


def second_switch_fields(a,known,initial,phase,hb,radius):
    g=a.g;t=g.symbol('true_second_switch_angular_t');sig=sigma(g,phase)
    angular=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))),
        a.fixed(g.mul(g.constant('-2/5'),hb,sig)))
    logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase,
        original_second_switch_full_current_radius_direction_and_power=True)
    phi=a.mul(initial['phi'],a.exp(logratio))
    slope=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,radius,g.sub(g.one,sig))),
        a.fixed(g.mul(g.constant('-2/5'),sig)))
    return dict(phi=phi,U=initial['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0),
        logratio=logratio,angular_phase_density=angular)


class CurrentOriginalPrescribedR100RmC6Functions:
    def __init__(self,bridge_owner,require_checked=True):
        if type(bridge_owner) is not previous.CurrentOriginalPrescribedBridgeC6Functions or not bridge_owner.acceptance_loaded:
            raise TypeError('Checked true prescribed bridge through R100 required')
        bridge_owner.assert_graph();self.bridge=bridge_owner;self.core=bridge_owner.core;self.g=g=bridge_owner.g
        self.family=copy.deepcopy(bridge_owner.family);self.hashes=dict(bridge_owner.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (previous.NAME,previous.RECEIPT):self._pin(name,sha(name))
        for module in (first_source,second_source,long_source,restore_source,switch_definition,long_definition,restore_definition,connection,native_long):
            self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        self.bindings=source_bindings();self.a=a=Algebra(g,6);self.b=b=Algebra(g,7);c=g.constant
        self.bridge_hb=bridge_owner.hb
        self.parameter_binding=parent.active.current.ast_binding(connection.build_graph)
        self.hb=hb=self.native_parameter('hbS','original switch_owner.flow.h; not merged with hbB by name')
        self.F0=F0=a.trim(bridge_owner.core7['F0']);self.P0=bridge_owner.P0
        self.Pstar=Pstar=bridge_owner.inlet.Pstar;self.delta=delta=self.core.parameters['delta']
        self.T=T=self.native_parameter('T','same original T=400*A singleton; all-N long window_length exact tuple guards')
        self.geometry_certificate=self.parameter_geometry_guard()
        self.charts={};self.endpoints={}
        s=g.symbol('true_prescribed_first_switch_s');R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))
        known,barphi=self.comparison_at_radius(R)
        initial={key:bridge_owner.R100['actual'][key] for key in ('phi','U')}
        actual=first_switch_fields(a,known,initial,s,hb,R,F0,Pstar,barphi)
        self.add_chart('first_switch',s,R,hb,actual,bridge_owner.R100['moments'],known=known)
        self.endpoints['R1']=self.end('first_switch',g.one)
        t=g.symbol('true_prescribed_second_switch_t');R=g.mul(c(100),g.unary('exp',g.mul(hb,g.add(g.one,t))))
        known,_=self.comparison_at_radius(R);incoming=self.endpoints['R1']
        actual=second_switch_fields(a,known,incoming['actual'],t,hb,R)
        self.add_chart('second_switch',t,R,hb,actual,incoming['moments'],known=known)
        self.endpoints['R2']=self.end('second_switch',g.one)
        p=g.symbol('true_prescribed_post_switch_fraction');self.post_length=g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))
        y=g.mul(self.post_length,p);R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))
        incoming=self.endpoints['R2'];phi=a.scale(incoming['actual']['phi'],g.unary('exp',g.mul(c('-2/5'),y)))
        actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0))
        self.add_chart('post_switch',p,R,self.post_length,actual,incoming['moments'])
        self.endpoints['R110']=self.end('post_switch',g.one)
        z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5);self.z=z
        self.logq=logq=a.log(a.add(a.fixed(1),a.mul(z,z)),'positive original 1+Z^2')
        incoming=self.endpoints['R110'];self.B=B=a.add(a.scale(a.trim(bridge_owner.core7['G']),g.neg(self.core.parameters['Lambda'])),
            a.log(incoming['actual']['phi'],'positive actual R110 phi'),a.fixed(g.mul(c('1/2'),g.unary('log',c(220)))),logq)
        q=g.symbol('true_prescribed_long_reshape_fraction');y=g.mul(T,q);R=g.mul(c(110),g.unary('exp',y))
        logUtheta=a.add(a.scale(B,g.sub(g.one,sigma(g,q))),a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))
        Utheta=a.exp(logUtheta);phi=self.phi_from_Utheta(Utheta,R)
        slope=a.add(a.fixed('-2/5'),a.scale(B,g.neg(g.quotient(sigma_prime(g,q),T,'same positive original T'))))
        actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0),
            physical_Utheta=Utheta,log_Utheta=logUtheta)
        self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])
        self.endpoints['Rsh']=self.end('long_reshape',g.one)
        self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))
        self.reference_length=length=g.sub(g.sub(logref,T),c(8))
        q=g.symbol('true_prescribed_reference_fraction');y=g.add(T,g.mul(length,q));R=g.mul(c(110),g.unary('exp',y))
        logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))
        Utheta=a.exp(logUtheta);phi=self.phi_from_Utheta(Utheta,R);incoming=self.endpoints['Rsh']
        actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),
            physical_Utheta=Utheta,log_Utheta=logUtheta)
        self.add_chart('reference',q,R,length,actual,incoming['moments'])
        self.endpoints['Rz']=self.end('reference',g.one)
        self.E=E=a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))
        t=g.symbol('true_prescribed_restoration_t');offset=g.add(c(-8),t);R=g.mul(c(110),g.unary('exp',g.add(logref,offset)))
        logUtheta=a.add(a.neg(logq),a.fixed(g.add(self.core.parameters['logP'],g.mul(c('1/10'),offset))))
        Utheta=a.exp(logUtheta);phi=self.phi_from_Utheta(Utheta,R)
        U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))
        actual=dict(phi=phi,U=U,phi_y=a.scale(phi,'-2/5'),U_y=a.scale(E,g.neg(sigma_prime(g,t))),
            physical_Utheta=Utheta,log_Utheta=logUtheta)
        self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])
        self.endpoints['Rrestore']=self.end('restoration',g.one)
        p=g.symbol('true_prescribed_post_restore_fraction');offset=g.add(c(-7),p);R=g.mul(c(110),g.unary('exp',g.add(logref,offset)))
        logUtheta=a.add(a.neg(logq),a.fixed(g.add(self.core.parameters['logP'],g.mul(c('1/10'),offset))))
        Utheta=a.exp(logUtheta);phi=self.phi_from_Utheta(Utheta,R)
        actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),
            physical_Utheta=Utheta,log_Utheta=logUtheta)
        self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])
        self.endpoints['Rm']=self.end('post_restore',g.one)
        self.definition=dict(true_prescribed_background_before_finite_N_loop=True,
            exact_original_both_R100_switch_controls_and_current_radius_retained=True,
            switch_hbS_and_preceding_bridge_hbB_source_nodes_not_aliased=True,
            exact_shared_upstream_K1_width_identity_and_positive_geometry_guarded=True,
            frozen_comparison_direction_not_actual_moment_replacement=True,
            full_original_finite_T_400A_not_shortened_or_infinite=True,
            original_B_from_actual_R110_and_anchored_G=True,
            all_six_nonzero_incoming_and_every_chart_memory_retained=True,
            exact_Rref_Rz_Rrestore_Rm_offsets=[-8,-7,-6],
            original_axial_restoration_4Z_plus_actual_E_times_one_minus_sigma=True,
            original_independent_P0_and_physical_Pstar_F0_factors_once=True,
            all_seven_source_rows_and_six_generic_output_rows_genuine=True,
            finite_N_modified_fields_correction_histories_and_all17_not_installed=True,
            new_original_scale_C6_numeric_containment_not_certified=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.bridge,self.core,self.g);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or pins!=self.hashes or raw['input_hashes']!=self.hashes
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['prescribed_R100_Rm_C6_bindings']!=self.snapshot or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked true prescribed R100-to-Rm functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed downstream prescribed source '+name)
            self.acceptance_loaded=True

    def comparison_at_radius(self,R):
        b=self.b;g=self.g;c=g.constant;old=self.bridge.comparison_macro;phase=old['phase']
        initial={name:b.at(jet,phase,g.one) for name,jet in old['moments'].items()}
        phi,U=old['phi'],old['U'];distance=g.unary('log',g.quotient(R,c(100),'same positive current radius over100'))
        eq=dict(H=phi,M=U,K=b.mul(phi,U),A=b.mul(U,U),B=b.scale(b.mul(phi,phi),'1/2'),C=b.mul(phi,phi))
        J={}
        for name,rate in previous.RATES.items():
            decay=g.unary('exp',g.mul(c(-rate),distance))
            J[name]=b.add(b.scale(initial[name],decay),b.scale(eq[name],g.sub(g.one,decay)))
        return previous.comparison_directions(g,phi,U,J,self.bridge.core7['F0'],self.P0,self.delta),phi
    def native_parameter(self,name,definition):
        found=[source.FunctionRef(self.g,i) for i,row in enumerate(self.g.nodes)
            if row['operation']=='current_original_source_parameter' and row.get('name')==name]
        if len(found)!=1:raise ValueError('Unique original source parameter required: '+name)
        row=self.g.nodes[found[0].node]
        original=self.g.nodes[self.bridge_hb.node]
        if (not row['Z_independent'] or row['definition']!=definition or row['source_family']!=self.family
            or row['fixed_leading_input_N0']!=original['fixed_leading_input_N0']
            or row['current_source_receipt']!=connection.current.RECEIPT
            or row['current_source_receipt_sha256']!=sha(connection.current.RECEIPT)):
            raise ValueError('Exact native parameter recipe/family/receipt required: '+name)
        return found[0]
    def parameter_geometry_guard(self):
        """Source identities and directed bounds; none defines a field value."""
        common=self.bridge.inlet.common;wb=common.source.bridge;sw=wb.switch;shared=wb.upstream
        if not (wb.acceptance_loaded and shared is sw.history.upstream and wb.logh is shared.logh
            and shared.logh._mpi_==sw.logh._mpi_ and wb.family==sw.family==common.family
            and wb.source==sw.source==common.source_sha):
            raise ValueError('Same admitted bridge/switch upstream and width source required')
        ledger=PREFIX+'K1_ledger.json';digest=sha(ledger)
        if not all(owner.hashes.get(ledger)==digest for owner in (wb,sw,shared)):
            raise ValueError('Identical actual K1 width definition and ledger required')
        self._pin(ledger,digest)
        for name in (native_long.NAME,native_long.RECEIPT):self._pin(name,sha(name))
        checked,raw=parent.linear.read(native_long.RECEIPT),parent.linear.read(native_long.NAME)
        if not checked['all_passed'] or raw['source_family']!=self.family:
            raise ValueError('Accepted canonical exact long-window tuple source required')
        lengths=[]
        for view in raw['actual_four_Z_all_N_long_patch_transports']:
            for window in view['actual_all_N_long_patch_windows']:
                if window['actual_chart']=='long_reshape':lengths.append(window['actual_original_window_length'])
        if len(lengths)!=4 or any(row!=lengths[0] for row in lengths):
            raise ValueError('Same frozen T across all original source cells required')
        exact=lengths[0]
        if exact['lower_exact_mpf_tuple']!=exact['upper_exact_mpf_tuple']:
            raise ValueError('Actual T singleton source, not a chosen interval endpoint')
        c=common.ctx;Tvalue=c.mpf(mp.make_mpf(tuple(exact['lower_exact_mpf_tuple'])))
        native=wb.anchor.patch.reference_mixed.long_mixed.reshape
        if native.T._mpi_!=Tvalue._mpi_ or native.core.logC._mpi_!=common.core.logC._mpi_:
            raise ValueError('Actual T/logC source tuples must equal the admitted native source')
        ep=previous.inlet.macro.ep;gap=10*(common.core.logC+common.core.logP)-Tvalue-8
        hcap=c.mpf('1/100');post_lower=c.ln(c.mpf('11/10'))-2*hcap
        if not (ep(Tvalue)[0]>0 and ep(gap)[0]>0 and ep(shared.logh)[1]<ep(c.ln(hcap))[0]
            and ep(post_lower)[0]>0 and ep(common.core.phi_floor)[0]>0):
            raise ValueError('Positive exact T/gap/post length and initial Phi required')
        self._width_sources=(wb,sw,shared,native)
        self._width_signature=(shared.logh._mpi_,native.T._mpi_,native.core.logC._mpi_)
        return dict(shared_upstream_is_same_object=True,bridge_logh_is_actual_upstream_object=True,
            switch_width_is_exact_source_reboxing=True,shared_K1_ledger_sha256=digest,
            actual_T_singleton_tuple=exact,all_four_native_T_source_tuples_equal=True,
            actual_native_T_and_logC_match_common_source=True,bridge_hbB_switch_hbS_kept_distinct=True,
            directed_T_positive=True,directed_reference_gap_positive=True,directed_post_length_positive=True,
            width_upper_bound='hb<1/100 is a bound only, never a source width',
            positive_R110_phi='positive true core inlet Phi times finite real exponentials',
            complete_radius_order='100<R1<R2<110<Rsh<Rz<Rrestore<Rm',
            bounds_are_not_defining_field_parameters=True)
    def phi_from_Utheta(self,Utheta,R):
        a=self.a;g=self.g
        return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R))),
            'same positive sqrt(2R)*F0; F0 is Z-only')
    def add_chart(self,name,phase,R,jacobian,actual,incoming,known=None):
        a=self.a
        J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)
        fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)
        self.charts[name]=dict(phase=phase,domain=[0,1],radius=R,physical_measure=jacobian,actual=actual,
            incoming_histories=incoming,moments=J,fields=fields,radial=radial,generic=generic,known_comparison_direction=known)
    def end(self,name,value):
        chart=self.charts[name];a=self.a;at=lambda jet:a.at(jet,chart['phase'],value)
        return dict(radius=self.g.node('function_substitution',expression=chart['radius'].node,
            variable=chart['phase'].node,value=value.node,Z_independent_substitution=True),
            actual={key:at(chart['actual'][key]) for key in ('phi','U','phi_y','U_y')},moments={k:at(jet) for k,jet in chart['moments'].items()},
            fields={k:at(jet) for k,jet in chart['fields'].items()},
            generic={k:Algebra(self.g,5).at(jet,chart['phase'],value) for k,jet in chart['generic'].items()})
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed prescribed R100-Rm input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(charts=self.charts,endpoints=self.endpoints,hb=self.hb,bridge_hb=self.bridge_hb,
            parameter_binding=self.parameter_binding,T=self.T,F0=self.F0,P0=self.P0,Pstar=self.Pstar,
            post_length=self.post_length,reference_length=self.reference_length,logref=self.logref,B=self.B,E=self.E,
            logq=self.logq,geometry_certificate=self.geometry_certificate,definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.bridge.assert_graph()
        wb,sw,shared,native=self._width_sources
        if not (wb.upstream is shared is sw.history.upstream and wb.logh is shared.logh
            and shared.logh._mpi_==sw.logh._mpi_
            and self._width_signature==(shared.logh._mpi_,native.T._mpi_,native.core.logC._mpi_)):
            raise ValueError('Native width/T/logC source identity changed')
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.bridge,self.core,self.g)):
            raise ValueError('True prescribed downstream field or incoming history changed')


def run(bridge_owner):
    began=time.monotonic();owner=CurrentOriginalPrescribedR100RmC6Functions(bridge_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,prescribed_R100_Rm_C6_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_PRESCRIBED_R100_RM_C6_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':raise SystemExit('run(bridge_owner) reuses the actual checked R100 functions')
