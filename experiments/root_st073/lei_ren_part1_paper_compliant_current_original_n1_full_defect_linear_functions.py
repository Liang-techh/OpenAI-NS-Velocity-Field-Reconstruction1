"""Actual full n1 signed C1 defects and five selected linear coefficients.

All pressure intervals use the retained same-graph source owners. The
repair velocities are C1 functions on the reserved compact supports.
Higher derivative bounds and repaired own M/V/P remain separate work.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions as cutoff
import lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator as linear
import lei_ren_part1_paper_compliant_current_original_n1_pressure_forcing_support as support

parent,source,controls,own=cutoff.parent,cutoff.source,cutoff.controls,cutoff.own
HERE,PREFIX,sha=cutoff.HERE,cutoff.PREFIX,cutoff.sha
NAME=PREFIX+'current_original_n1_full_defect_linear_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_full_defect_linear_functions_check.json'
GATES=('n1_complete_signed_defect_C1_function_adapter_installed',
       'n1_selected_linear_coefficient_functions_installed',
       'current_original_n1_linear_repair_velocity_C1_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in cutoff.OPEN if k not in GATES)+(
    'n1_repair_velocity_C3_functions_installed',
    'n1_repaired_own_M_V_P_installed','positive_order_five_moment_repair_closed',
    'uniform_full_pressure_derivative_accuracy_certified')))
MOMENTS=cutoff.MOMENTS
COEFFICIENTS=('alpha1','alpha2','beta1','beta2','beta3')


def compose_full_moments(g,endpoints,lower,preRp,pulse):
    """Rin endpoints plus the three exact exterior pressure intervals."""
    moments={key:source.C1Function(q.value,q.Z) for key,q in endpoints.items()}
    for key in ('m13','m15'):
        for segment in (lower,preRp,pulse):moments[key]=g.c1add(moments[key],segment[key])
    return moments


def physical_response_matrix(g,W,background):
    """Paper physical five-moment response, coefficient order alpha,beta."""
    fixed=lambda q:source.C1Function(q,g.zero)
    scale=g.c1scale;mul=g.c1mul;c=g.constant
    Av,Az=(source.FunctionRef(g,n) for n in background['Astar_ordinary_Z0_Z3'][:2])
    A=source.C1Function(Av,Az);L=W['log2Rc'];mu=W['mu']
    exp=lambda q:g.unary('exp',g.mul(q,L))
    sr=exp(c('1/2'));inverse_power=exp(g.neg(mu))
    z0=[fixed(g.mul(sr,a)) for a in W['axial_rows'][0]]
    zm=[fixed(g.mul(sr,inverse_power,g.add(a,g.mul(mu,d))))
        for a,d in zip(*W['axial_rows'])]
    t0=[fixed(g.mul(exp(g.one),q)) for q in W['angular_rows'][0]]
    t1=[mul(scale(g.mul(c(2),exp(g.neg(g.add(g.one,mu)))),A),fixed(q)) for q in W['angular_rows'][1]]
    t2=[mul(scale(g.neg(inverse_power),A),fixed(q)) for q in W['angular_rows'][2]]
    zero=fixed(g.zero)
    return dict(m11=z0+[zero]*3,m12=[zero]*2+t0,m13=[zero]*2+t1,
        m14=[mul(A,q) for q in zm]+[zero]*3,m15=[zero]*2+t2)


def response_residuals(g,moments,coefficients,matrix):
    out={}
    for key in MOMENTS:
        value=moments[key]
        for B,C in zip(matrix[key],coefficients):value=g.c1add(value,g.c1mul(B,C))
        out[key]=value
    return out


def repair_velocity_functions(g,W,coefficients,background):
    zero=source.C1Function(g.zero,g.zero)
    U=zero;theta=zero
    for q,j in zip(coefficients[:2],linear.AXIAL):U=g.c1add(U,g.c1scale(W['physical_bumps'][j],q))
    for q,j in zip(coefficients[2:],linear.ANGULAR):theta=g.c1add(theta,g.c1scale(W['physical_bumps'][j],q))
    s=W['bump_log_radius_variable'];logRc=source.FunctionRef(g,background['logRc'])
    R=g.unary('exp',g.add(logRc,s));sqrt2R=g.unary('sqrt',g.mul(g.constant(2),R))
    F=g.c1scale(g.quotient(g.one,sqrt2R,'positive reserved physical radius'),theta)
    return dict(Uz1=U,Utheta1=theta,F1=F,physical_radius=R,log_R_over_Rc=s)


def zero_from_old_support(g,node,zero_nodes,cache):
    """Exact zero propagation, only after disjoint flat support is proved."""
    if node in zero_nodes:return True
    if node in cache:return cache[node]
    row=g.nodes[node];op=row['operation'];at=lambda i:zero_from_old_support(g,i,zero_nodes,cache)
    if op=='exact_rational':value=row['numerator']==0
    elif op=='sum':value=all(at(i) for i in row['arguments'])
    elif op=='product':value=any(at(i) for i in row['arguments'])
    elif op=='negative':value=at(row['argument'])
    elif op=='positive_quotient':value=at(row['numerator'])
    else:value=False
    cache[node]=value;return value


def reserved_background_guard(operator):
    """Actual leading Uz=0 and pure swirl power on every new support.

    Nonlocal mean and fixed leading Omega/N1p are not set to zero.
    Their unchanged terms affect the defect, not its linear variation.
    """
    g=operator.graph;p=operator.background['old_profile_functions']
    if operator.support!=linear.exact_supports():raise ValueError('Original disjoint support changed')
    if not p['compact_support_all_endpoint_jets_zero']:raise ValueError('Actual old flat support required')
    zero_nodes=set(p['bumps']+p['bump_y']+p['bump_yy']);cache={};checked=[]
    for i,row in enumerate(p['radial_y_rows']):
        # Preserved C2 rows can use different arithmetic nodes for V and
        # delta_V. Prove zero of the actual V expression on this support.
        for key in ('V','delta_E'):
            for k,node in enumerate(row[key]):
                if not zero_from_old_support(g,node,zero_nodes,cache):
                    raise ValueError('Nonzero leading local field on reserved support '+str((key,i,k)))
                checked.append(dict(quantity=key,ordinary_y_order=i,ordinary_Z_order=k,node=node))
    return dict(exact_new_and_old_supports=operator.support,
        actual_local_Uz_and_angular_increment_zero_rows=checked,
        leading_swirl_equals_original_power_on_new_support=True,
        nonlocal_mean_pressure_and_leading_N1p_not_zeroed=True,
        fixed_Omega_cancels_from_linear_pressure_variation=True)


class CurrentOriginalN1FullDefectLinearFunctions:
    def __init__(self,cutoff_owner=None,require_checked=True):
        self.cutoff=cutoff_owner or cutoff.CurrentOriginalN1InnerCutoffFunctions()
        if not self.cutoff.acceptance_loaded:raise ValueError('Checked actual inner cutoff required')
        self.cutoff.assert_graph();self.g=g=self.cutoff.g;self.family=copy.deepcopy(self.cutoff.family)
        self.hashes=dict(self.cutoff.hashes)
        for name,digest in parent.linear.read(cutoff.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (cutoff.NAME,cutoff.RECEIPT):self._pin(name,sha(name))
        self.lower=self.cutoff.inner.collar;self.pulse=self.lower.core.pulse
        self.Gamma=self.pulse.own.jets.complete.Gamma;self.preRp=self.Gamma.pulse.preRp
        self.linear=self.preRp.operator
        self.sources=(self.cutoff,self.lower,self.preRp,self.pulse,self.linear)
        for owner in self.sources:
            owner.assert_graph()
            og=owner.graph if owner is self.linear else owner.g
            if og is not g or owner.family!=self.family or not owner.acceptance_loaded:
                raise ValueError('Exact same checked source graph and family required')
        for module in (linear,support):
            receipt=parent.linear.read(module.RECEIPT)
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Checked actual linear/support source required')
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        zero=parent.linear.read(support.NAME)
        if (zero['endpoint']['definition']!='Rv=exp(original_logRp+13/original_mu)'
            or not zero['endpoint']['angular_and_pressure_not_zeroed']
            or set(zero['future_charts'])!=set(support.FUTURE)
            or any(row['exact_zero_unweighted_and_R_weighted_forcing_integrals']!=[0,0]
                   for row in zero['future_charts'].values())):
            raise ValueError('Actual full beyond-Rv meridional forcing support required')
        self.beyond_Rv_support=copy.deepcopy(zero)
        self.reserved_background=reserved_background_guard(self.linear)
        self.segments=dict(inner={key:source.C1Function(q.value,q.Z) for key,q in self.cutoff.endpoint.items()},
            lower=self.lower.functions,preRp=self.preRp.functions,pulse=self.pulse.functions,
            beyond_Rv={key:source.C1Function(g.zero,g.zero) for key in ('m13','m15')})
        self.moments=compose_full_moments(g,self.cutoff.endpoint,self.lower.functions,self.preRp.functions,self.pulse.functions)
        values=self.linear.coefficient_functions(self.moments)
        self.coefficients=dict(zip(COEFFICIENTS,values))
        self.response_matrix=physical_response_matrix(g,self.linear.functions,self.linear.background)
        self.linear_response_residuals=response_residuals(g,self.moments,values,self.response_matrix)
        self.repair_velocity=repair_velocity_functions(g,self.linear.functions,values,self.linear.background)
        self.partition=dict(velocity_m11_m12_m14='full cutoff support, constant after q=3/4',
            inner='0 to Rin',lower='Rin to r_minus',preRp='r_minus to Rp: active/compact/quiet',
            pulse='Rp to Rv: main/gap/end',beyond_Rv='exact zero N1p; angular/pressure heat field retained',
            Gamma_future_is_upstream_pulse_input_not_an_added_pressure_segment=True,
            m15_signed_contributions_added_without_second_negation=True,
            no_graph_import_or_old_numeric_sector_selection=True)
        self.definition=dict(full_five_signed_C1_functions=True,selected_five_coefficients_are_actual_C1_functions=True,
            repair_velocity_C1_functions=True,physical_swirl_repair_divided_by_sqrt_2R_once=True,
            original_nonzero_mu_retained_in_divided_axial_inverse=True,
            coefficient_norms_and_higher_derivative_bounds_transferred=False,
            repaired_own_M_V_P_installed=False,positive_order_five_moment_repair_closed=False,
            original_numeric_n1_solution_oracle_executed=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or raw['exact_graph_nodes']!=g.nodes
                or raw['full_defect_linear_bindings']!=self.snapshot):raise ValueError('Checked actual full linear functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual full defect source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed full defect dependency '+name)
        self.hashes[name]=digest

    def _snapshot(self):return own.encode(dict(segments=self.segments,moments=self.moments,coefficients=self.coefficients,
        response_matrix=self.response_matrix,linear_response_residuals=self.linear_response_residuals,
        repair_velocity=self.repair_velocity,partition=self.partition,definition=self.definition,
        reserved_background=self.reserved_background,beyond_Rv_support=self.beyond_Rv_support))

    def assert_graph(self):
        self.cutoff.assert_graph()
        for owner in self.sources:owner.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('Full original n1 defect/coefficient graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1FullDefectLinearFunctions(require_checked=False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,full_defect_linear_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_FULL_DEFECT_AND_LINEAR_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':run()
