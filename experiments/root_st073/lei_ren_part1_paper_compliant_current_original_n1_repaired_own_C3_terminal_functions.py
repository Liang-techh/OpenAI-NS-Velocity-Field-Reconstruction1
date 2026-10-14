"""Actual repaired reserved C3 histories, Q/V C2 and terminal integral zeros.

The full current compact pressure supplies the unpatched prefix. Its fixed
leading N1p is retained. Terminal zeros require actual integral reduction;
global joins, quantitative norms and higher-order recursion remain open.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_complete_C3_defect_repair_functions as previous
import lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities as original_terminal
import lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields as original_own

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode,Algebra=previous.FunctionJet,previous.rows,previous.encode,previous.Algebra
MOMENTS=previous.MOMENTS
NAME=PREFIX+'current_original_n1_repaired_own_C3_terminal_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_repaired_own_C3_terminal_functions_check.json'
GATES=('current_original_n1_reserved_own_C3_histories_QV_C2_functions_installed',
       'current_original_n1_C3_repaired_own_terminal_five_function_zeros',
       'current_original_n1_mass_pressure_C3_QV_C2_tail_functions_installed',
       'current_original_n1_reserved_band_total_pressure_y_C3_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_repaired_own_C3_numeric_containment_installed',
    'current_original_n1_repaired_own_global_mixed_norms_installed')))


def patch_histories(g,complete):
    """Actual Section15 physical densities, ordinary Z0..3 and dR=R ds."""
    a=Algebra(g,3);patch=complete.repair_velocity;W=complete.linear.functions
    bg=complete.linear.background;R=patch['physical_radius'];s=patch['log_R_over_Rc']
    A=FunctionJet(tuple(source.FunctionRef(g,i) for i in bg['Astar_ordinary_Z0_Z3']))
    power=g.unary('exp',g.neg(g.mul(g.add(g.constant('1/2'),W['mu']),
        g.add(g.unary('log',g.constant(2)),source.FunctionRef(g,bg['logRc']),s))))
    theta0=a.scale(A,power);U,theta,F=(patch[k] for k in ('Uz1','Utheta1','F1'))
    sqrt2R=g.unary('sqrt',g.mul(g.constant(2),R))
    per_dR=dict(m11=U,m12=a.scale(F,g.mul(g.constant(2),R)),
        m13=a.scale(a.mul(theta0,theta),g.quotient(g.one,R,'positive reserved physical radius')),
        m14=a.scale(a.mul(theta0,U),sqrt2R),m15=a.neg(a.mul(theta0,theta)))
    per_ds={k:a.scale(q,R) for k,q in per_dR.items()}
    t=g.symbol('true_repaired_n1_C3_reserved_history_t')
    bodies={k:a.at(q,s,t) for k,q in per_ds.items()}
    histories={k:a.integral(q,t,g.zero,s,physical_measure='dR=R ds',paper_moment=k,
        original_selected_alpha_beta_functions=True,history_not_reset_between_bumps=True,
        fixed_leading_N1p_not_subtracted=True,actual_signed_repaired_n1_history=True) for k,q in bodies.items()}
    return histories,dict(per_dR=per_dR,per_ds=per_ds,bodies=bodies,variable=t,
                         leading_Utheta0=theta0,physical_F1=F)


def own_Q(g,U,mean,delta):
    """Positive-order recovery: genuine mean C3 supplies Q C2."""
    a=Algebra(g,2);z=FunctionJet((g.symbol('Z'),g.one,g.zero));one=a.fixed(1)
    de=a.fixed(delta);d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    numerator=a.sub(a.scale(a.mul(z,a.trim(U)),2),
        a.add(a.mul(a.add(one,de),a.mul(z,a.trim(mean))),a.mul(d,a.shift(mean))))
    return a.div(numerator,L,'original positive L=1-delta*Z^2 on common axial domain')


def recovery_fields(g,patch,mass,pressure,delta,patch_pressure_y,total_pressure_y):
    a=Algebra(g,3);b=Algebra(g,2);R=patch['physical_radius']
    mean=a.div(mass,a.fixed(R),'positive nonaxis physical radius R');Q=own_Q(g,patch['Uz1'],mean,delta)
    return dict(**{k:patch[k] for k in ('Uz1','Utheta1','F1')},
        mean1=mean,K1=a.sub(mean,patch['Uz1']),M1=mass,P1=pressure,
        Q1=Q,V1=b.scale(Q,R),mean1_y=a.sub(patch['Uz1'],mean),
        M1_y=a.scale(patch['Uz1'],R),P1_y_patch=patch_pressure_y,P1_y_total=total_pressure_y)


class CurrentOriginalN1RepairedOwnC3TerminalFunctions:
    def __init__(self,complete,legacy_terminal,require_checked=True):
        if type(complete) is not previous.CurrentOriginalN1CompleteC3DefectRepairFunctions or not complete.acceptance_loaded:
            raise TypeError('Checked actual complete C3 repair source required')
        if type(legacy_terminal) is not original_terminal.CurrentOriginalN1TerminalIntegralIdentities or not legacy_terminal.acceptance_loaded:
            raise TypeError('Checked actual physical bump integral theorem required')
        complete.assert_graph();legacy_terminal.assert_graph()
        self.complete,self.legacy_terminal=complete,legacy_terminal;self.g=g=complete.g
        self.family=copy.deepcopy(complete.family);self.hashes=dict(complete.hashes)
        if legacy_terminal.g is not g or legacy_terminal.family!=self.family or legacy_terminal.full is not complete.full:
            raise ValueError('Same actual graph, source and reserved operator required')
        for module in (previous,original_terminal,original_own):
            rec=parent.linear.read(module.RECEIPT)
            if not rec['all_passed'] or rec['source_family']!=self.family:raise ValueError('Admitted same-source theorem required')
            for name,digest in rec['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        integral_check=parent.linear.read(original_terminal.RECEIPT)['checks']['actual_integrals']
        if (not integral_check['passed'] or integral_check['basis_rows']!=25
            or not integral_check['actual_beta_normalizer_and_support_domain_bound']):
            raise ValueError('Actual translated physical integral theorem required')
        self.quiet=complete.quiet;self.compact=self.quiet.compact;self.patch=complete.repair_velocity
        self.s=self.patch['log_R_over_Rc'];self.R=self.patch['physical_radius'];self.delta=complete.cutoff.inner.parameters['delta']
        a=Algebra(g,3);self.patch_histories,self.patch_bodies=patch_histories(g,complete)
        self.band_pressure={k:a.at(q,self.compact.x,g.unary('exp',self.s)) for k,q in self.compact.pressure.items()}
        self.active_pressure={k:a.trim(q) for k,q in self.compact.inlet['pressure'].items()}
        self.compact_local_pressure={k:a.sub(q,self.active_pressure[k]) for k,q in self.band_pressure.items()}
        self.incoming={k:complete.moments[k] for k in ('m11','m12','m14')}
        for k in ('m13','m15'):self.incoming[k]=a.add(complete.segments['inner'][k],complete.segments['lower'][k],self.active_pressure[k])
        self.cumulative={k:a.add(self.incoming[k],self.patch_histories[k],self.compact_local_pressure[k])
            if k in ('m13','m15') else a.add(self.incoming[k],self.patch_histories[k]) for k in MOMENTS}
        # Direct second radial recovery uses genuine compact V_y,V_yy,m_y
        # rows; the complete pressure primitive itself already retains its
        # original S0/S1 integral and finite boundaries.
        radial=previous.previous.radial
        second=radial.physical_pressure_fields(g,
            dict(V=self.compact.radial['V_y'],m=self.compact.radial['m_y']),
            dict(V_y=self.compact.radial['V_yy']),self.compact.delta,self.compact.Pstar)
        self.compact_Qyy=second['Q_y']
        compact_R=self.compact.background['radius']
        self.fixed_N1p=a.add(a.scale(self.compact.S,'-1/2'),
            a.scale(a.add(a.trim(self.compact.physical['Q_y']),a.trim(self.compact_Qyy)),
                g.quotient(g.one,compact_R,'positive actual compact radius')))
        self.fixed_pressure_y=a.at(a.scale(self.fixed_N1p,compact_R),self.compact.x,g.unary('exp',self.s))
        self.total_pressure_y=a.add(self.fixed_pressure_y,self.patch_bodies['per_ds']['m13'])
        self.fields=recovery_fields(g,self.patch,self.cumulative['m11'],self.cumulative['m13'],self.delta,
                                    self.patch_bodies['per_ds']['m13'],self.total_pressure_y)
        end=g.unary('log',g.constant(2))
        self.band_endpoint={k:a.at(q,self.s,end) for k,q in self.cumulative.items()}
        self.patch_endpoint={k:a.at(q,self.s,end) for k,q in self.patch_histories.items()}
        self.quiet_increment={k:a.sub(a.trim(complete.preRp_pressure[k]),a.trim(self.compact.terminal['pressure'][k]))
                              for k in ('m13','m15')}
        self.terminal={k:a.add(self.band_endpoint[k],self.quiet_increment[k],complete.segments['pulse'][k])
            if k in ('m13','m15') else self.band_endpoint[k] for k in MOMENTS}
        self.reduced_patch={k:a.add(*(a.mul(B,C) for B,C in zip(complete.response_matrix[k],complete.coefficients.values())))
                            for k in MOMENTS}
        self.reduced_terminal={k:a.add(complete.moments[k],self.reduced_patch[k]) for k in MOMENTS}
        self.integral_identity={k:g.node('proved_original_n1_C3_transported_patch_integral_identity',
            actual_endpoint_rows=encode(self.patch_endpoint[k]),equivalent_physical_response_rows=encode(self.reduced_patch[k]),
            actual_signed_density_rows=encode(self.patch_bodies['per_ds'][k]),
            original_basis_integral_receipt=original_terminal.RECEIPT,
            original_basis_integral_receipt_sha256=sha(original_terminal.RECEIPT),
            physical_measure='dR=R ds; s=center+ell*t',ordinary_Z_orders=[0,1,2,3],
            supports_and_endpoints_independent_of_Z=True,
            theorem='same exact 25 physical basis integrals; linearity and ordinary differentiation of smooth compact integrands')
            for k in MOMENTS}
        self.tail_origins=dict(legacy_terminal.tail_origins)
        self.tail_origins['pressure_source_support']=copy.deepcopy(legacy_terminal.tail_origins['pressure_source_support'])
        self.mass_coordinate=g.symbol('true_repaired_n1_C3_log_R_over_2Rc')
        self.pressure_coordinate=g.symbol('true_repaired_n1_C3_log_R_over_Rv')
        Rt=g.unary('exp',g.add(self.tail_origins['mass_log_radius'],self.mass_coordinate));mean=a.scale(self.terminal['m11'],
            g.quotient(g.one,Rt,'positive actual mass-tail radius'))
        zero=a.fixed(0);b=Algebra(g,2);Q=own_Q(g,zero,mean,self.delta)
        self.mass_tail=dict(M1=self.terminal['m11'],mean1=mean,K1=mean,Uz1=zero,Utheta1=zero,F1=zero,
                            Q1=Q,V1=b.scale(Q,Rt))
        self.pressure_tail=dict(P1=self.terminal['m13'])
        self.terminal_zero_functions={k:a.fixed(0) for k in MOMENTS}
        self.mass_tail_zero_functions={k:Algebra(g,len(rows(q))-1).fixed(0) for k,q in self.mass_tail.items()}
        self.pressure_tail_zero_functions=dict(P1=a.fixed(0))
        self.bindings=dict(original_signed_density=parent.active.current.ast_binding(original_own.patch_histories),
            original_positive_order_recovery=parent.active.current.ast_binding(original_own.recovery_fields),
            original_tail_recovery=parent.active.current.ast_binding(original_terminal.mass_tail_fields),
            original_physical_basis_receipt=original_terminal.RECEIPT,
            actual_compact_pressure=previous.previous.previous.NAME,
            actual_compact_second_radial_recovery=parent.active.current.ast_binding(radial.physical_pressure_fields))
        self.definition=dict(actual_signed_integral_histories_C3_not_caps=True,
            actual_current_compact_pressure_prefix_C3=True,
            active_rminus_Rc_pressure_explicitly_retained_in_incoming=True,
            local_compact_increment_subtracts_same_active_prefix_once=True,
            inner_lower_pressure_datum_retained_once=True,fixed_leading_N1p_retained=True,
            mean_M_K_pressure_C3_and_QV_C2_using_true_mean_Z3=True,
            own_recovery_weight_one_plus_delta=True,leading_pressure_weight_one_minus_delta=True,
            compact_direct_Qyy_C4_and_total_P1_y_C3=True,
            actual_transported_terminal_integrals_C3_equal_selected_response=True,
            mass_mean_K_C3_QV_C2_zero_after_2Rc=True,n1_pressure_C3_zero_after_Rv=True,
            leading_angular_and_heat_pressure_not_zeroed=True,
            terminal_identity_not_interval_overlap_or_algebraic_response_alone=True,
            queries_restricted_to_chart_and_checked_zero_certificates=True,
            lower_active_field_join_and_global_n1_mixed_bounds_open=True,
            quantitative_original_scale_C3_norms_and_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(complete,legacy_terminal,self.quiet,self.compact,g,self.patch);self.acceptance_loaded=False
        if require_checked:
            raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not rec.get('all_passed') or not all(rec.get(k) for k in GATES) or any(rec.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['repaired_own_C3_terminal_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual C3 own fields and terminal integral proof required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed repaired own C3 dependency '+name)
            self.acceptance_loaded=True
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed repaired own C3 source '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(patch_histories=self.patch_histories,patch_bodies=self.patch_bodies,
            band_pressure=self.band_pressure,active_pressure=self.active_pressure,
            compact_local_pressure=self.compact_local_pressure,incoming=self.incoming,cumulative=self.cumulative,
            compact_Qyy=self.compact_Qyy,fixed_N1p=self.fixed_N1p,fixed_pressure_y=self.fixed_pressure_y,
            total_pressure_y=self.total_pressure_y,fields=self.fields,band_endpoint=self.band_endpoint,
            patch_endpoint=self.patch_endpoint,quiet_increment=self.quiet_increment,terminal=self.terminal,
            reduced_patch=self.reduced_patch,reduced_terminal=self.reduced_terminal,integral_identity=self.integral_identity,
            tail_origins=self.tail_origins,mass_coordinate=self.mass_coordinate,pressure_coordinate=self.pressure_coordinate,
            mass_tail=self.mass_tail,pressure_tail=self.pressure_tail,terminal_zero_functions=self.terminal_zero_functions,
            mass_tail_zero_functions=self.mass_tail_zero_functions,pressure_tail_zero_functions=self.pressure_tail_zero_functions,
            bindings=self.bindings,definition=self.definition))
    def assert_graph(self):
        self.complete.assert_graph();self.legacy_terminal.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.complete,self.legacy_terminal,self.complete.quiet,self.complete.quiet.compact,self.g,self.complete.repair_velocity)):
            raise ValueError('Actual repaired own C3 source or transported history changed')
    def at_band_fraction(self,fraction):
        self.assert_graph()
        if not isinstance(fraction,(int,str,Fraction)) or isinstance(fraction,bool):raise TypeError('Exact rational fraction required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Own-field chart requires 0<=s/log(2)<=1')
        value=self.g.mul(self.g.constant(q),self.g.unary('log',self.g.constant(2)))
        fields={k:Algebra(self.g,len(rows(v))-1).at(v,self.s,value) for k,v in self.fields.items()}
        return dict(exact_s_over_log2=str(q),fields=fields,mean_M_K_P_C3=True,QV_C2=True,total_pressure_y_C3=True)
    def tail_at_ratio(self,kind,ratio):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Tail zeros require completed actual integral proof')
        if kind not in ('mass','pressure'):raise ValueError('Choose mass or pressure tail')
        if not isinstance(ratio,(int,str,Fraction)) or isinstance(ratio,bool):raise TypeError('Exact rational ratio required')
        q=Fraction(ratio)
        if q<1:raise ValueError('Tail ratio must be at least one')
        fields=self.mass_tail if kind=='mass' else self.pressure_tail
        zeros=self.mass_tail_zero_functions if kind=='mass' else self.pressure_tail_zero_functions
        x=self.mass_coordinate if kind=='mass' else self.pressure_coordinate;value=self.g.unary('log',self.g.constant(q))
        return dict(kind=kind,exact_radius_ratio=str(q),origin='2Rc' if kind=='mass' else 'Rv',
            actual_source_functions={k:Algebra(self.g,len(rows(v))-1).at(v,x,value) for k,v in fields.items()},
            certified_equal_zero_functions=zeros,actual_terminal_integral_C3_proof_consumed=True,
            leading_angular_and_heat_pressure_unchanged=True)


def run(complete,legacy_terminal):
    began=time.monotonic();owner=CurrentOriginalN1RepairedOwnC3TerminalFunctions(complete,legacy_terminal,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,repaired_own_C3_terminal_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_REPAIRED_OWN_C3_TERMINAL_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(complete,legacy_terminal) reuses accepted source owners')
