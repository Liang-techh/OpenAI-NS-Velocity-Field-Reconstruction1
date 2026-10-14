"""Actual pulse third/fourth radial sources and own pressure yy/yyy.

Scalar flat branches are retained before evaluating interior quotients.
Ordinary Z rows and physical-y derivatives use independent conventions.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_pulse_internal_radial_source_joins as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,radial,rp=previous.outer,previous.radial,previous.rp
flat=rp.pulse5.true.own.c2.flat_source
NAME=PREFIX+'current_original_pulse_higher_radial_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_pulse_higher_radial_pressure_C3_functions_check.json'
GATES=('current_original_pulse_Uyyy_Uyyyy_C5_Qyyy_Qyyyy_C4_source_functions_installed',
    'current_original_pulse_own_pressure_second_radial_C3_source_functions_installed',
    'current_original_pulse_own_pressure_third_radial_C3_source_functions_installed',
    'current_original_pulse_internal_higher_radial_pressure_C3_source_joins_installed',
    'current_original_n1_Rp_Rv_higher_radial_pressure_C3_source_joins_installed')
OPEN=tuple(dict.fromkeys(tuple(key for key in previous.OPEN if key not in GATES)+(
    'current_original_n1_2Rc_higher_radial_pressure_C3_source_join_installed',)))


def power(g,x,n):
    return g.mul(*([x]*n)) if n else g.one


def sigma_high(g,x,n,pin):
    if n not in (3,4):raise ValueError('Actual sigma third/fourth row required')
    c=g.constant;s=rp.pulse5.true.own.jets.complete.Gamma.sigma(g,x)
    one_minus=g.sub(g.one,x);q=g.mul(s,g.sub(g.one,s));a=g.sub(g.one,g.mul(c(2),s))
    inv=lambda k,sign:g.add(g.quotient(c(sign),power(g,x,k),'interior 0<x<1 only'),
        g.quotient(c(abs(sign)),power(g,one_minus,k),'interior 0<x<1 only'))
    L1,L2,L3,L4=inv(3,2),inv(4,-6),inv(5,24),inv(6,-120)
    if n==3:
        body=g.mul(q,g.add(L3,g.mul(c(3),a,L1,L2),g.mul(g.sub(g.one,g.mul(c(6),q)),power(g,L1,3))))
    else:
        body=g.mul(q,g.add(L4,g.mul(a,g.add(g.mul(c(4),L1,L3),g.mul(c(3),L2,L2))),
            g.mul(c(6),g.sub(g.one,g.mul(c(6),q)),L1,L1,L2),
            g.mul(a,g.sub(g.one,g.mul(c(12),q)),power(g,L1,4))))
    return g.node('analytic_unary',name='original_flat_sigma_order'+str(n),argument=x.node,
        interior_body=body.node,active_interval=[0,1],flat_value=0,
        active_body_not_evaluated_on_flat_branch=True,derivative_order=n,
        ordinary_from_Taylor_factor=math.factorial(n),derivative_source=pin,
        exact_function_not_a_numerical_bound=True)


def beta_high(g,x,n,pin):
    if n not in (3,4):raise ValueError('Actual beta third/fourth row required')
    coefficients=list(flat.BETA_POLYNOMIALS[n]);c=g.constant;poly=g.zero
    for coefficient in reversed(coefficients):poly=g.add(g.mul(poly,x),c(coefficient))
    w=g.sub(g.one,g.mul(x,x));beta=g.node('compact_raw_beta',argument=x.node,
        definition='exp(-1/(1-r^2)) on abs(r)<1; zero otherwise')
    body=g.quotient(g.mul(beta,poly),power(g,w,2*n),'interior abs(r)<1 only')
    return g.node('analytic_unary',name='original_compact_beta_order'+str(n),argument=x.node,
        interior_body=body.node,active_interval=[-1,1],flat_value=0,polynomial_coefficients=coefficients,
        active_body_not_evaluated_on_flat_branch=True,derivative_order=n,
        ordinary_from_Taylor_factor=math.factorial(n),derivative_source=pin,
        exact_function_not_a_numerical_bound=True)


def higher_pressure_density(g,f,delta,R,old_density):
    a=outer.Algebra(g,3);z=outer.FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero));one=a.fixed(1)
    de=a.fixed(delta);d=a.sub(one,a.mul(z,z));h=a.sub(one,de);L=a.sub(one,a.mul(de,a.mul(z,z)))
    div=lambda value:a.div(value,L,'same admitted positive leading L')
    U,Uy,Uyy=[a.trim(f[key]) for key in ('U','U_y','U_yy')]
    Q,Qy,Qyy,Qyyy,Qyyyy=[a.trim(f[key]) for key in ('Q','Q_y','Q_yy','Q_yyy','Q_yyyy')]
    Qz,Qyz,Qyyz=[a.shift(f[key]) for key in ('Q','Q_y','Q_yy')]
    A1=a.add(Qy,a.scale(a.mul(h,a.mul(z,Qyz)),'1/2'),Qyy)
    A2=a.add(Qyy,a.scale(a.mul(h,a.mul(z,Qyyz)),'1/2'),Qyyy)
    C=a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2))
    C1=a.add(a.scale(a.mul(z,Qy),-2),a.mul(d,Qyz),a.scale(a.mul(z,Qyy),-2))
    C2=a.add(a.scale(a.mul(z,Qyy),-2),a.mul(d,Qyyz),a.scale(a.mul(z,Qyyy),-2))
    B1=a.add(a.mul(Q,Qy),a.mul(Qy,Qy),a.mul(Q,Qyy))
    B2=a.add(a.mul(Qy,Qy),a.mul(Q,Qyy),a.scale(a.mul(Qy,Qyy),3),a.mul(Q,Qyyy))
    Sy=a.add(div(a.add(A1,a.mul(Uy,C),a.mul(U,C1))),B1)
    Syy=a.add(div(a.add(A2,a.mul(Uyy,C),a.scale(a.mul(Uy,C1),2),a.mul(U,C2))),B2)
    S,Py=old_density['S'],old_density['P1_y']
    Pyy=a.add(a.scale(a.add(S,Sy),g.mul(g.constant('-1/2'),R)),Qyy,Qyyy)
    Pyyy=a.add(a.scale(a.add(S,a.scale(Sy,2),Syy),g.mul(g.constant('-1/2'),R)),Qyyy,Qyyyy)
    invR=g.quotient(g.one,R,'positive actual pulse radius')
    return dict(S_y=Sy,S_yy=Syy,P1_yy=Pyy,P1_yyy=Pyyy,
        N1p_y=a.scale(a.sub(Pyy,Py),invR),N1p_yy=a.scale(a.add(Pyyy,a.scale(Pyy,-2),Py),invR),
        m15_yy=a.scale(a.add(Py,Pyy),g.neg(R)),
        m15_yyy=a.scale(a.add(Py,a.scale(Pyy,2),Pyyy),g.neg(R)))


class CurrentOriginalPulseHigherRadialPressureC3Functions:
    def __init__(self,internal=None,require_checked=True):
        self.internal=internal if internal is not None else previous.CurrentOriginalPulseInternalRadialSourceJoins()
        if not self.internal.acceptance_loaded:raise ValueError('Checked actual internal pulse joins required')
        self.internal.assert_graph();self.pulse=self.internal.pulse;self.inlet=self.internal.inlet;self.joins=self.internal.joins
        self.g=g=self.internal.g;self.family=copy.deepcopy(self.internal.family);self.hashes=dict(self.internal.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        self.sigma_source=parent.active.current.ast_binding(flat._sigma_left)
        self.beta_source=parent.active.current.ast_binding(flat.beta_jets)
        self.beta_polynomial_source=parent.active.current.ast_binding(flat.beta_polynomials)
        for pin in (self.sigma_source,self.beta_source,self.beta_polynomial_source):
            if self.hashes[pin['module']]!=pin['module_sha256']:raise ValueError('Pinned original flat source required')
        self.charts={};a=outer.Algebra(g,5);bp=g.add(g.constant('1/2'),self.pulse.mu)
        for name,data in self.pulse.charts.items():
            x=data['coordinate'];B=[data[key] for key in ('normalized_axial','normalized_axial_y','normalized_axial_yy')]
            for n in (3,4):
                value=a.fixed(0)
                if name=='main':
                    I=self.pulse.at_scalar(self.pulse.primitive,x)
                    left=g.mul(g.constant(50),x);right=g.sub(g.constant(11),x)
                    sl=list(self.pulse.sigma_rows(left))+[sigma_high(g,left,3,self.sigma_source),sigma_high(g,left,4,self.sigma_source)]
                    sr=list(self.pulse.sigma_rows(right))+[sigma_high(g,right,3,self.sigma_source),sigma_high(g,right,4,self.sigma_source)]
                    terms=[]
                    for k in range(n+1):
                        primitive_derivative=I if k==0 else g.mul(g.constant(50**(k-1)),sl[k-1])
                        terms.append(g.mul(g.constant(math.comb(n,k)*(-1)**(n-k)),primitive_derivative,sr[n-k]))
                    value=a.scale(self.pulse.selected['ap'],g.mul(power(g,self.pulse.mu,n),g.add(*terms)))
                elif name=='end':
                    ell,norm=self.pulse.W['ell'],self.pulse.W['normalization']
                    for key,center in (('C1',-3),('C2',-1)):
                        argument=g.quotient(g.sub(x,g.constant(center)),ell,'positive original ell')
                        derivative=beta_high(g,argument,n,self.beta_source)
                        value=a.add(value,a.scale(self.pulse.selected[key],g.quotient(derivative,
                            g.mul(power(g,ell,n+1),norm),'positive original ell and beta normalization')))
                    value=a.scale(value,g.unary('exp',self.pulse.W['log_end_scale']))
                elif name!='gap':raise ValueError('Actual original pulse chart required')
                B.append(value)
            high_U=[]
            for n in (3,4):
                normalized=a.add(*(a.scale(B[k],g.mul(g.constant(math.comb(n,k)*(-1)**(n-k)),power(g,bp,n-k)))
                    for k in range(n+1)))
                high_U.append(a.mul(data['fields']['Utheta'],normalized))
            fields=dict(data['fields'],U_yyy=high_U[0],U_yyyy=high_U[1])
            fields['m_yyy']=a.sub(fields['U_yy'],fields['m_yy'])
            fields['m_yyyy']=a.sub(fields['U_yyy'],fields['m_yyy'])
            fields['Q_yyy']=radial.leading_Q(g,fields['U_yyy'],fields['m_yyy'],self.pulse.delta)
            fields['Q_yyyy']=radial.leading_Q(g,fields['U_yyyy'],fields['m_yyyy'],self.pulse.delta)
            density=higher_pressure_density(g,fields,self.pulse.delta,data['radius'],data['density'])
            self.charts[name]=dict(coordinate=x,domain=data['domain'],radius=data['radius'],
                physical_y_Jacobian=data['physical_y_Jacobian'],fields=fields,higher_density=density,
                normalized_axial_yyy=B[3],normalized_axial_yyyy=B[4])
        self.endpoints={name:{edge:self.at(name,value) for edge,value in zip(('left','right'),data['domain'])}
            for name,data in self.charts.items()}
        # The admitted quiet source already saves real y0..4 power-flow
        # rows. Recover physical higher fields from those rows, rather
        # than filling the quiet side from a shared pulse endpoint.
        q=copy.deepcopy(self.joins.quiet);P=self.ref(q['Pstar']);R=self.ref(q['R'])
        for key,order,count in (('flow',5,6),('leading_flow',6,7)):
            if q[key]['order']!=order:raise ValueError('Exact admitted quiet ordinary order required')
            for packet,field in (('profiles_y0_to_y4','V'),('histories_y0_to_y4','m')):
                if len(q[key][packet])!=5 or any(len(row[field])!=count for row in q[key][packet]):
                    raise ValueError('Exact quiet y0..4 ordinary source row lengths required')
        physical={key:self.jet(ids) for key,ids in q['physical'].items()}
        for suffix,k in (('yy',2),('yyy',3),('yyyy',4)):
            V=q['flow']['profiles_y0_to_y4'][k]['V'];m=q['flow']['histories_y0_to_y4'][k]['m']
            physical['U_'+suffix]=a.scale(self.jet(V),P)
            physical['m_'+suffix]=a.scale(self.jet(m),P)
        for suffix in ('yyy','yyyy'):
            physical['Q_'+suffix]=radial.leading_Q(g,physical['U_'+suffix],physical['m_'+suffix],self.pulse.delta)
        a3=outer.Algebra(g,3);quiet_old=dict(S=self.jet(q['S']),P1_y=a3.scale(self.jet(q['N1p']),R))
        quiet_high=higher_pressure_density(g,physical,self.pulse.delta,R,quiet_old)
        self.quiet=dict(source=q,coordinate=self.ref(q['s']),radius=R,fields=physical,
            old_density=quiet_old,higher_density=quiet_high)
        keys=('P1_yy','P1_yyy','N1p_y','N1p_yy','m15_yy','m15_yyy')
        self.outer_pressure_traces={
            'Rp':dict(left={key:a3.at(quiet_high[key],self.ref(q['s']),self.ref(q['length'])) for key in keys},
                right={key:self.endpoints['main']['left']['higher_density'][key] for key in keys}),
            'Rv':dict(left={key:self.endpoints['end']['right']['higher_density'][key] for key in keys},
                right={key:a3.fixed(0) for key in keys})}
        self.retained_tail_pressure=self.joins.tail_pressure
        self.definition=dict(actual_sigma_beta_third_fourth_interior_bodies_and_flat_branches=True,
            ordinary_source_rows_not_interval_bounds_or_Taylor_coefficients=True,
            actual_B_U_mass_Q_third_fourth_radial_source_functions=True,
            actual_P1_yy_P1_yyy_derivatives_of_same_retained_pressure_FTC=True,
            actual_quiet_y3_y4_saved_flow_rows_and_once_Pstar_used=True,
            actual_both_Rp_pressure_traces_and_constant_Rv_derivative_tail_preserved=True,
            actual_R_y_equals_R_no_extra_mu=True,
            leading_one_minus_delta_recovery_kept_separate_from_own_one_plus_delta=True,
            original_nonzero_means_normalizations_and_histories_retained=True,
            outer_higher_joins_full_mixed_numeric_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.internal,self.pulse,self.inlet,self.joins,g)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return outer.FunctionJet(tuple(self.ref(i) for i in ids))
    def at(self,name,value):
        data=self.charts[name];x=data['coordinate']
        group=lambda packet:{key:outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,x,value) for key,jet in packet.items()}
        return dict(fields=group(data['fields']),higher_density=group(data['higher_density']),
            radius=outer.Algebra(self.g,3).at(outer.Algebra(self.g,3).fixed(data['radius']),x,value).value)
    def _snapshot(self):
        return outer.encode(dict(sigma_source=self.sigma_source,beta_source=self.beta_source,
            beta_polynomial_source=self.beta_polynomial_source,charts=self.charts,endpoints=self.endpoints,
            quiet=self.quiet,outer_pressure_traces=self.outer_pressure_traces,
            retained_tail_pressure=self.retained_tail_pressure,definition=self.definition))
    def assert_graph(self):
        self.internal.assert_graph()
        if (self._owners!=(self.internal,self.pulse,self.inlet,self.joins,self.g)
            or self.hashes!=self._hashes or self._snapshot()!=self.snapshot
            or self.g.nodes[:len(self.prefix)]!=self.prefix or self.numerical_oracle_available):
            raise ValueError('Actual higher radial pressure source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed higher radial source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py'
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or archive!=sha(NAME) or pins.pop(checker)!=sha(checker) or pins!=self.hashes
            or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['pulse_higher_radial_pressure_C3_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked higher pulse radial pressure required')
        self.acceptance_loaded=True


def run(internal=None):
    began=time.monotonic();owner=CurrentOriginalPulseHigherRadialPressureC3Functions(internal,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        pulse_higher_radial_pressure_C3_bindings=owner.snapshot,
        **dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_PULSE_HIGHER_RADIAL_PRESSURE_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
