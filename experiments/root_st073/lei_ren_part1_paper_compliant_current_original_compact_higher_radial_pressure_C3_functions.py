"""Actual compact/own y4 sources and total pressure y2/y3 at 2Rc.

Modified background and own positive-order recovery keep separate signs.
All original branch, controls, nonzero histories and pressure anchors remain.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_pulse_higher_radial_pressure_C3_functions as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,radial,rp=previous.outer,previous.radial,previous.rp
compact=rp.quiet.previous
NAME=PREFIX+'current_original_compact_higher_radial_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_compact_higher_radial_pressure_C3_functions_check.json'
GATES=('current_original_compact_modified_profile_history_y4_C5_source_functions_installed',
    'current_original_compact_Qyyy_Qyyyy_C4_pressure_yy_yyy_C3_source_functions_installed',
    'current_original_n1_reserved_own_y4_C3_QV_C2_source_functions_installed',
    'current_original_n1_flat_outer_own_y4_C3_QV_C2_source_functions_installed',
    'current_original_n1_2Rc_higher_radial_pressure_C3_source_join_installed',
    'current_original_n1_2Rc_own_y4_C3_QV_C2_source_joins_installed')
OPEN=tuple(dict.fromkeys(k for k in previous.OPEN if k not in GATES))


def product_row(a,left,right,j):
    return a.add(*(a.scale(a.mul(left[k],right[j-k]),math.comb(j,k)) for k in range(j+1)))


class ReservedHighRadialDerivative(rp.pulse3.OrdinaryZDerivative):
    """The actual ordinary radial chain rule, with exact guarded beta rows."""
    def __init__(self,g,phase,coefficients,beta_source):
        super().__init__(g);self.derivatives[self.z.node]=g.zero;self.z=phase
        self.derivatives[phase.node]=g.one;self.beta_source=beta_source;self.beta_chains={}
        for jet in coefficients.values():
            for ref in outer.rows(jet):self.seed(ref,g.zero)
    def rule(self,i):
        row=self.g.nodes[i]
        if row['operation']=='compact_raw_beta':
            argument=row['argument'];dx=self.derivatives[argument]
            if dx==self.g.zero:return self.g.zero
            arg=source.FunctionRef(self.g,argument);raw=source.FunctionRef(self.g,i)
            first,second=compact.leading.beta_first_two(self.g,arg,raw)
            third=previous.beta_high(self.g,arg,3,self.beta_source)
            fourth=previous.beta_high(self.g,arg,4,self.beta_source)
            self.seed(first,self.g.mul(second,dx));self.seed(second,self.g.mul(third,dx))
            self.seed(third,self.g.mul(fourth,dx))
            self.beta_chains[i]=(raw,first,second,third,fourth)
            return self.g.mul(first,dx)
        return super().rule(i)


class CurrentOriginalCompactHigherRadialPressureC3Functions:
    def __init__(self,high=None,require_checked=True):
        self.high=high if high is not None else previous.CurrentOriginalPulseHigherRadialPressureC3Functions()
        if not self.high.acceptance_loaded:raise ValueError('Checked higher pulse/quiet source required')
        self.high.assert_graph();self.g=g=self.high.g;self.family=copy.deepcopy(self.high.family)
        self.hashes=dict(self.high.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.c=copy.deepcopy(self.high.inlet.inputs[compact.NAME]);c=self.c
        self.own_outer=self.high.joins.outer;self.x=self.ref(c['background']['phase'])
        self.W=self.decode(c['W']);self.mu=self.ref(c['mu']);self.delta=self.ref(c['delta'])
        self.Pstar=self.ref(c['Pstar']);self.R=self.ref(c['background']['radius'])
        self.beta_source=self.high.beta_source;self.beta_polynomial_source=self.high.beta_polynomial_source
        a=outer.Algebra(g,5);A=a.trim(self.jet(c['A']));h=[self.jet(ids) for ids in c['h']]
        if len(h)!=5 or any(len(outer.rows(jet))!=6 for jet in h):raise ValueError('Actual compact C5 controls required')
        self.profiles=[self.jets(row) for row in c['profiles']['radial_y_rows']]
        if len(self.profiles)!=3:raise ValueError('Exact admitted compact radial0..2 prefix required')
        bumps=[list(row) for row in compact.bump_rows(g,self.W,self.x)]
        invell=g.quotient(g.one,self.W['ell'],'original positive log bump width')
        scale=g.quotient(g.one,g.mul(self.W['ell'],self.W['normal'],self.x),'original ell J0 x')
        self.compact_beta_chains=[]
        for i,center in enumerate(self.W['centers']):
            arg=g.mul(g.sub(g.unary('log',self.x),center),invell)
            beta=compact.leading.raw_beta(g,arg);first,second=compact.leading.beta_first_two(g,arg,beta)
            chain=[beta,first,second]+[previous.beta_high(g,arg,n,self.beta_source) for n in (3,4)]
            self.compact_beta_chains.append(chain)
            for n in (3,4):
                bumps[i].append(g.mul(scale,g.add(*(g.mul(g.constant(math.comb(n,k)*(-1)**(n-k)),
                    chain[k],previous.power(g,invell,k)) for k in range(n+1)))))
        self.bumps=bumps;alpha=self.ref(c['profiles']['alpha']);fixed_power=self.ref(c['profiles']['power'])
        invN=g.quotient(g.one,self.ref(c['N']),'same original exact positive integer N')
        for n in (3,4):
            F=a.scale(a.add(*(a.scale(h[i+2],bumps[i][n]) for i in range(3))),invN)
            G=a.scale(a.add(a.scale(h[0],bumps[0][n]),a.scale(h[1],bumps[2][n])),invN)
            E0=a.scale(A,g.mul(fixed_power,g.constant((-1)**n),previous.power(g,alpha,n)))
            dE,dV=a.mul(A,F),a.mul(A,G)
            self.profiles.append(dict(F=F,G=G,original_E=E0,delta_E=dE,delta_V=dV,E=a.add(E0,dE),V=dV))
        E,V=([row[key] for row in self.profiles] for key in ('E','V'))
        self.histories=[{key:self.jet(c['modified_fields'][key]) for key in compact.RATES}]
        for n in range(4):
            EE=product_row(a,E,E,n)
            density=dict(m=V[n],h=E[n],k=product_row(a,E,V,n),
                e=a.sub(product_row(a,V,V,n),a.scale(EE,'1/2')),p=a.scale(EE,'1/2'))
            self.histories.append({key:a.sub(density[key],a.scale(self.histories[-1][key],rate))
                for key,rate in compact.RATES.items()})
        self.physical={}
        for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
            self.physical['U'+suffix]=a.scale(V[j],self.Pstar)
            self.physical['m'+suffix]=a.scale(self.histories[j]['m'],self.Pstar)
            self.physical['Q'+suffix]=radial.leading_Q(g,self.physical['U'+suffix],self.physical['m'+suffix],self.delta)
        self.density=radial.pressure_density(g,self.physical,self.delta,self.R)
        self.higher_density=previous.higher_pressure_density(g,self.physical,self.delta,self.R,self.density)
        self.retained_compact_pressure=self.jets(c['pressure'])
        self.compact_endpoint=self.at_compact(g.constant(2))
        # The own fields are a separate source, with genuine C3 input rows.
        self.phase=self.own_outer.reserved_phase;self.own_R=self.own_outer.reserved_radius
        fields=self.own_outer.reserved_fields;coefficients=self.own_outer.coefficients
        if any(len(outer.rows(jet))!=4 for jet in coefficients.values()):raise ValueError('Actual own C3 coefficients required')
        self.coefficient_phase_independence=copy.deepcopy(self.own_outer.coefficient_phase_independence)
        D=ReservedHighRadialDerivative(g,self.phase,coefficients,self.beta_source)
        own={key:[fields[key]] for key in ('Uz1','Utheta1','F1')}
        for j in range(4):
            for key in own:
                own[key].append(outer.FunctionJet(tuple(D.derivative(ref) for ref in outer.rows(own[key][-1]))))
        self.own_beta_chains={str(i):chain for i,chain in D.beta_chains.items()}
        b=outer.Algebra(g,3);v=outer.Algebra(g,2)
        own['mean1']=[fields['mean1']]
        for j in range(4):own['mean1'].append(b.sub(own['Uz1'][j],own['mean1'][-1]))
        own['M1']=[fields['M1']]
        for j in range(1,5):
            own['M1'].append(b.scale(b.add(*(b.scale(own['Uz1'][k],math.comb(j-1,k)) for k in range(j))),self.own_R))
        own['K1']=[b.sub(mean,U) for mean,U in zip(own['mean1'],own['Uz1'])]
        own['Q1']=[outer.own_Q(g,U,mean,self.delta) for U,mean in zip(own['Uz1'],own['mean1'])]
        own['V1']=[v.scale(v.add(*(v.scale(own['Q1'][k],math.comb(j,k)) for k in range(j+1))),self.own_R)
            for j in range(5)]
        self.own_reserved=own
        # Preserve the actual own pressure value and accepted first FTC.
        repaired=self.own_outer.inputs[outer.intermediate.previous.NAME]
        self.repaired=copy.deepcopy(repaired)
        theta0=self.jet(repaired['patch_bodies']['leading_Utheta0']);bp=g.add(g.constant('1/2'),self.mu)
        self.theta0_radial=[b.scale(theta0,g.mul(g.constant((-1)**j),previous.power(g,bp,j))) for j in range(3)]
        self.patch_pressure_density=[self.jet(repaired['patch_bodies']['per_ds']['m13'])]
        self.patch_pressure_density.extend(product_row(b,self.theta0_radial,own['Utheta1'],j) for j in (1,2))
        at_own=lambda jet:b.at(jet,self.x,g.unary('exp',self.phase))
        self.own_pressure={
            'P1':self.jet(repaired['fields']['P1']),
            'P1_y':self.jet(repaired['total_pressure_y']),
            'P1_yy':b.add(at_own(self.higher_density['P1_yy']),self.patch_pressure_density[1]),
            'P1_yyy':b.add(at_own(self.higher_density['P1_yyy']),self.patch_pressure_density[2])}
        self.own_m15_higher={
            'm15_yy':b.scale(b.add(self.own_pressure['P1_y'],self.own_pressure['P1_yy']),g.neg(self.own_R)),
            'm15_yyy':b.scale(b.add(self.own_pressure['P1_y'],b.scale(self.own_pressure['P1_yy'],2),self.own_pressure['P1_yyy']),g.neg(self.own_R))}
        self.own_outside={}
        packets={'quiet':self.own_outer.charts['quiet']}
        packets.update({name:self.high.joins.current_own[name] for name in ('main','gap','end')})
        for name,packet in packets.items():
            f=packet['fields'];rows={key:[f[key]]+[b.fixed(0) for _ in range(4)] for key in ('Uz1','Utheta1','F1','M1')}
            for key in ('mean1','K1'):
                rows[key]=[b.scale(f[key],(-1)**j) for j in range(5)]
            rows['Q1']=[v.scale(f['Q1'],(-1)**j) for j in range(5)]
            rows['V1']=[f['V1']]+[v.fixed(0) for _ in range(4)]
            self.own_outside[name]=rows
        self.own_2Rc_traces=dict(left=self.at_own(g.unary('log',g.constant(2))),right=self.at_quiet_own(g.zero))
        self.pressure_2Rc_traces=dict(
            left={key:b.at(jet,self.phase,g.unary('log',g.constant(2))) for key,jet in self.own_pressure.items() if key in ('P1_yy','P1_yyy')},
            right={key:b.at(self.high.quiet['higher_density'][key],self.high.quiet['coordinate'],g.zero) for key in ('P1_yy','P1_yyy')})
        self.m15_2Rc_traces=dict(left={key:b.at(jet,self.phase,g.unary('log',g.constant(2))) for key,jet in self.own_m15_higher.items()},
            right={key:b.at(self.high.quiet['higher_density'][key],self.high.quiet['coordinate'],g.zero) for key in self.own_m15_higher})
        self.definition=dict(actual_compact_original_radial_prefix_preserved=True,
            actual_guarded_beta3_beta4_log_bump_functions=True,
            actual_modified_complete_nonzero_history_y4_C5_functions=True,
            actual_leading_one_minus_delta_Q4_pressure_y2_y3_C3=True,
            actual_own_positive_order_one_plus_delta_y4_C3_QV_C2=True,
            actual_own_pressure_patch_theta0_theta1_product_derivatives_kept=True,
            original_nonzero_mass_and_pressure_values_not_reset=True,
            full_outer_mixed_numeric_global_n1_recursion_corrected_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.high,self.own_outer,g);self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return outer.FunctionJet(tuple(self.ref(i) for i in ids))
    def jets(self,packet):return {key:self.jet(ids) for key,ids in packet.items()}
    def decode(self,value):
        if isinstance(value,dict):return {key:self.decode(v) for key,v in value.items()}
        if isinstance(value,list):return [self.decode(v) for v in value]
        return self.ref(value)
    def at_compact(self,x):
        at=lambda packet:{key:outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,self.x,x) for key,jet in packet.items()}
        return dict(profiles=[at(row) for row in self.profiles],histories=[at(row) for row in self.histories],
            physical=at(self.physical),density=at(self.density),higher_density=at(self.higher_density))
    def at_own(self,x):
        return {key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,self.phase,x) for jet in row]
            for key,row in self.own_reserved.items()}
    def at_quiet_own(self,s):
        phase=self.high.quiet['coordinate']
        return {key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,phase,s) for jet in row]
            for key,row in self.own_outside['quiet'].items()}
    def _snapshot(self):
        return outer.encode(dict(c=self.c,W=self.W,x=self.x,R=self.R,mu=self.mu,delta=self.delta,Pstar=self.Pstar,
            profiles=self.profiles,bumps=self.bumps,compact_beta_chains=self.compact_beta_chains,
            histories=self.histories,physical=self.physical,density=self.density,higher_density=self.higher_density,
            retained_compact_pressure=self.retained_compact_pressure,compact_endpoint=self.compact_endpoint,
            phase=self.phase,own_R=self.own_R,own_beta_chains=self.own_beta_chains,
            coefficient_phase_independence=self.coefficient_phase_independence,
            own_reserved=self.own_reserved,repaired=self.repaired,theta0_radial=self.theta0_radial,
            patch_pressure_density=self.patch_pressure_density,own_pressure=self.own_pressure,
            own_m15_higher=self.own_m15_higher,own_outside=self.own_outside,
            own_2Rc_traces=self.own_2Rc_traces,pressure_2Rc_traces=self.pressure_2Rc_traces,
            m15_2Rc_traces=self.m15_2Rc_traces,definition=self.definition))
    def assert_graph(self):
        self.high.assert_graph()
        if (self._owners!=(self.high,self.own_outer,self.g) or self.hashes!=self._hashes
            or self._snapshot()!=self.snapshot or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.numerical_oracle_available):raise ValueError('Actual compact higher source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed compact higher source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py'
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or archive!=sha(NAME) or pins.pop(checker)!=sha(checker) or pins!=self.hashes
            or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['compact_higher_radial_pressure_C3_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked compact higher radial source required')
        self.acceptance_loaded=True


def run(high=None):
    began=time.monotonic();owner=CurrentOriginalCompactHigherRadialPressureC3Functions(high,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        compact_higher_radial_pressure_C3_bindings=owner.snapshot,
        **dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_COMPACT_HIGHER_RADIAL_PRESSURE_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
