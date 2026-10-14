"""Complete signed n1 C3 partition, selected inverse and repair velocities.

Reuse the already admitted inner cutoff C3 endpoints. The numerical inner
bridge's C1 scope does not reduce those function rows. Algebraic response
cancellation does not install the global joins, repaired own fields or norms.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions as previous
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as original
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions as pulse

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode,Algebra=previous.FunctionJet,previous.rows,previous.encode,previous.Algebra
MOMENTS,COEFFICIENTS=original.MOMENTS,original.COEFFICIENTS
NAME=PREFIX+'current_original_n1_complete_C3_defect_repair_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_complete_C3_defect_repair_functions_check.json'
GATES=('current_original_n1_complete_signed_partition_C3_functions_installed',
       'current_original_n1_selected_reserved_inverse_C3_functions_installed',
       'n1_repair_velocity_C3_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_complete_C3_defect_numeric_containment_installed',
    'current_original_n1_selected_reserved_C3_coefficient_norms_installed',
    'n1_repaired_own_M_V_P_installed','positive_order_five_moment_repair_closed',
    'current_original_n1_C3_repaired_own_terminal_five_function_zeros')))


def compose_full_moments(g,endpoints,lower,preRp,pulse_pressure):
    """0..Rin endpoint plus separately signed exterior pressure pieces."""
    a=Algebra(g,3);out={k:a.trim(v) for k,v in endpoints.items()}
    for k in ('m13','m15'):
        out[k]=a.add(out[k],*(a.trim(segment[k]) for segment in (lower,preRp,pulse_pressure)))
    return out


def coefficient_functions(g,W,background,moments):
    """The original finite-mu selected inverse, with all ordinary Z0..3."""
    if set(moments)!=set(MOMENTS):raise ValueError('All five signed C3 moments required')
    a=Algebra(g,3)
    A=FunctionJet(tuple(source.FunctionRef(g,i) for i in background['Astar_ordinary_Z0_Z3']))
    if len(rows(A))!=4:raise ValueError('Actual amplitude C3 source required')
    divA=lambda q:a.div(a.trim(q),A,'same actual strictly positive original Astar')
    factor=lambda k:g.unary('exp',g.mul(k,W['log2Rc']))
    mu=W['mu'];invsr=factor(g.constant('-1/2'))
    y0=a.scale(moments['m11'],g.neg(invsr))
    ym=a.scale(divA(moments['m14']),g.neg(g.mul(invsr,factor(mu))))
    yd=a.scale(a.sub(ym,y0),g.quotient(g.one,mu,'same original strictly positive mu'))
    theta=(a.scale(moments['m12'],g.neg(factor(g.constant(-1)))),
           a.scale(divA(moments['m13']),g.mul(g.constant('-1/2'),factor(g.add(g.one,mu)))),
           a.scale(divA(moments['m15']),factor(mu)))
    rhs=(y0,yd,*theta)
    columns=[original.linear.inverse_normalized_action(g,W,[rows(q)[k] for q in rhs]) for k in range(4)]
    return [FunctionJet(tuple(column[j] for column in columns)) for j in range(5)],dict(Astar=A,normalized_rhs=rhs)


def physical_response_matrix(g,W,background):
    """Lift the unchanged original physical response through ordinary Z3."""
    a=Algebra(g,3);fixed=a.fixed;c=g.constant
    A=FunctionJet(tuple(source.FunctionRef(g,i) for i in background['Astar_ordinary_Z0_Z3']))
    exp=lambda q:g.unary('exp',g.mul(q,W['log2Rc']))
    mu=W['mu'];sr=exp(c('1/2'));inv=exp(g.neg(mu));zero=fixed(0)
    z0=[fixed(g.mul(sr,q)) for q in W['axial_rows'][0]]
    zm=[fixed(g.mul(sr,inv,g.add(q,g.mul(mu,d)))) for q,d in zip(*W['axial_rows'])]
    t0=[fixed(g.mul(exp(g.one),q)) for q in W['angular_rows'][0]]
    t1=[a.scale(A,g.mul(c(2),exp(g.neg(g.add(g.one,mu))),q)) for q in W['angular_rows'][1]]
    t2=[a.scale(A,g.neg(g.mul(inv,q))) for q in W['angular_rows'][2]]
    return dict(m11=z0+[zero]*3,m12=[zero]*2+t0,m13=[zero]*2+t1,
                m14=[a.mul(A,q) for q in zm]+[zero]*3,m15=[zero]*2+t2)


def response_residuals(g,moments,coefficients,matrix):
    a=Algebra(g,3)
    return {k:a.add(moments[k],*(a.mul(B,C) for B,C in zip(matrix[k],coefficients))) for k in MOMENTS}


def repair_velocity_functions(g,W,coefficients,background):
    a=Algebra(g,3);zero=a.fixed(0)
    U=a.add(zero,*(a.scale(q,W['physical_bumps'][j]) for q,j in zip(coefficients[:2],original.linear.AXIAL)))
    theta=a.add(zero,*(a.scale(q,W['physical_bumps'][j]) for q,j in zip(coefficients[2:],original.linear.ANGULAR)))
    s=W['bump_log_radius_variable'];logRc=source.FunctionRef(g,background['logRc'])
    R=g.unary('exp',g.add(logRc,s));sqrt2R=g.unary('sqrt',g.mul(g.constant(2),R))
    F=a.scale(theta,g.quotient(g.one,sqrt2R,'positive reserved physical radius'))
    return dict(Uz1=U,Utheta1=theta,F1=F,physical_radius=R,log_R_over_Rc=s)


class CurrentOriginalN1CompleteC3DefectRepairFunctions:
    def __init__(self,quiet,full,pulse_owner,require_checked=True):
        expected=(previous.CurrentOriginalTrueQuietRpC5PressureC3Functions,
                  original.CurrentOriginalN1FullDefectLinearFunctions,
                  pulse.CurrentOriginalN1PulsePressureC3Functions)
        for owner,cls in zip((quiet,full,pulse_owner),expected):
            if type(owner) is not cls or not owner.acceptance_loaded:raise TypeError('Checked original source owners required')
            owner.assert_graph()
        self.quiet,self.full,self.pulse=quiet,full,pulse_owner
        self.cutoff=full.cutoff;self.lower=quiet.lower;self.linear=full.linear
        self.g=g=quiet.g;self.family=copy.deepcopy(quiet.family);self.hashes=dict(quiet.hashes)
        self.sources=(quiet,full,pulse_owner,self.cutoff,self.lower,self.linear)
        for owner in self.sources:
            og=owner.graph if owner is self.linear else owner.g
            if og is not g or owner.family!=self.family or not owner.acceptance_loaded:
                raise ValueError('Same checked source graph, family and selected operator required')
        for module in (previous,original,pulse):
            rec=parent.linear.read(module.RECEIPT)
            if not rec['all_passed'] or rec['source_family']!=self.family:raise ValueError('Checked component receipt required')
            for name,digest in rec['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        a=Algebra(g,3)
        self.inner_endpoint=self.cutoff.endpoint
        self.lower_pressure=self.lower.pressure
        self.preRp_pressure=quiet.Rp_endpoint['pressure']
        self.pulse_pressure=pulse_owner.pressure
        self.segments=dict(inner={k:a.trim(v) for k,v in self.inner_endpoint.items()},
            lower={k:a.trim(v) for k,v in self.lower_pressure.items()},
            preRp={k:a.trim(v) for k,v in self.preRp_pressure.items()},
            pulse={k:a.trim(v) for k,v in self.pulse_pressure.items()},
            beyond_Rv={k:a.fixed(0) for k in ('m13','m15')})
        self.moments=compose_full_moments(g,self.inner_endpoint,self.lower_pressure,self.preRp_pressure,self.pulse_pressure)
        values,self.inverse_bindings=coefficient_functions(g,self.linear.functions,self.linear.background,self.moments)
        self.coefficients=dict(zip(COEFFICIENTS,values))
        self.response_matrix=physical_response_matrix(g,self.linear.functions,self.linear.background)
        self.linear_response_residuals=response_residuals(g,self.moments,values,self.response_matrix)
        self.repair_velocity=repair_velocity_functions(g,self.linear.functions,values,self.linear.background)
        self.reserved_background=copy.deepcopy(full.reserved_background)
        self.beyond_Rv_support=copy.deepcopy(full.beyond_Rv_support)
        self.partition=dict(inner='0 to Rin: all five admitted cutoff C3 endpoints',
            lower='Rin to rminus',preRp='rminus to Rp: actual active/compact/quiet',
            pulse='Rp to Rv',beyond_Rv='exact zero meridional N1p, angular/heat retained',
            velocity_moments='m11/m12/m14 constant after inner cutoff support',
            lower_component_added_once=True,partitioned_Rin_Rp_not_used_as_preRp=True,
            m15_signed_components_added_without_second_negation=True)
        self.bindings=dict(original_C1_coefficient_action=parent.active.current.ast_binding(original.linear.coefficient_action),
            original_selected_inverse=parent.active.current.ast_binding(original.linear.inverse_normalized_action),
            original_physical_response=parent.active.current.ast_binding(original.physical_response_matrix),
            original_repair_velocity=parent.active.current.ast_binding(original.repair_velocity_functions),
            inner_endpoint_source=original.cutoff.RECEIPT,lower_source=previous.lower.RECEIPT,
            preRp_source=previous.RECEIPT,pulse_source=pulse.RECEIPT,
            inner_C3_was_already_admitted_before_this_milestone=True)
        self.definition=dict(complete_signed_partition_C3_functions=True,
            same_selected_finite_mu_inverse_C3_functions=True,actual_Astar_Z0_to_Z3_used=True,
            repair_Uz_Utheta_F_C3_functions=True,ordinary_derivatives_not_Taylor_coefficients=True,
            physical_F_equals_Utheta_divided_by_sqrt_2R_once=True,
            response_zero_is_algebraic_not_terminal_repaired_field_proof=True,
            lower_active_field_join_installed=False,selected_native_full_C5_pulse_bridge_installed=False,
            quantitative_C3_bounds_transferred=False,repaired_own_M_V_P_installed=False,
            global_n1_recursion_and_corrected_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(quiet,full,pulse_owner,self.cutoff,self.lower,self.linear,g,
                      self.inner_endpoint,self.lower_pressure,self.preRp_pressure,self.pulse_pressure)
        self.acceptance_loaded=False
        if require_checked:
            raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not rec.get('all_passed') or not all(rec.get(k) for k in GATES) or any(rec.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['complete_C3_defect_repair_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked complete signed C3 partition and inverse required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed complete C3 source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed complete C3 input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(segments=self.segments,moments=self.moments,coefficients=self.coefficients,
            inverse_bindings=self.inverse_bindings,response_matrix=self.response_matrix,
            linear_response_residuals=self.linear_response_residuals,repair_velocity=self.repair_velocity,
            reserved_background=self.reserved_background,beyond_Rv_support=self.beyond_Rv_support,
            partition=self.partition,bindings=self.bindings,definition=self.definition))
    def assert_graph(self):
        for owner in self.sources:owner.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.quiet,self.full,self.pulse,self.cutoff,self.lower,self.linear,self.g,
                self.cutoff.endpoint,self.lower.pressure,self.quiet.Rp_endpoint['pressure'],self.pulse.pressure)):
            raise ValueError('Actual complete C3 source ownership or graph changed')


def run(quiet,full,pulse_owner):
    began=time.monotonic();owner=CurrentOriginalN1CompleteC3DefectRepairFunctions(quiet,full,pulse_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,complete_C3_defect_repair_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('COMPLETE_SIGNED_C3_DEFECT_AND_RESERVED_REPAIR_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(quiet,full,pulse_owner) reuses checked sources')
