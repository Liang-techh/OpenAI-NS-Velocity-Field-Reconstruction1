"""Actual leading collar radial4, pressure rho2 and full finite own Rin join.

The original analytic core and nonzero collar histories define every row.
The high own cutoff extension is restricted to its genuine flat collar.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_n1_Rin_source_join as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.outer
rc_source=previous.previous
high=rc_source.high
NAME=PREFIX+'current_original_lower_higher_radial_Rin_join.json.gz'
RECEIPT=PREFIX+'current_original_lower_higher_radial_Rin_join_check.json'
GATES=('current_original_actual_lower_radial4_C5_Q_radial4_C4_source_functions_installed',
    'current_original_actual_lower_pressure_rho2_C3_source_functions_installed',
    'current_original_n1_Rin_full_y4_pressure_y3_mixed_source_join_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)
CONTRACT=copy.deepcopy(rc_source.previous.CONTRACT)


def product_row(a,left,right,n):
    return a.add(*(a.scale(a.mul(left[k],right[n-k]),math.comb(n,k)) for k in range(n+1)))


def core_definition(row):
    result=copy.deepcopy(row)
    for key in ('radial_order','ordinary_slow_Z_derivative_order'):result.pop(key)
    return result


class CurrentOriginalLowerHigherRadialRinJoin:
    def __init__(self,rin=None,require_checked=True):
        self.rin=rin if rin is not None else previous.CurrentOriginalN1RinSourceJoin()
        if not self.rin.acceptance_loaded:raise ValueError('Accepted actual finite Rin source required')
        self.rin.assert_graph();self.g=g=self.rin.g;self.family=copy.deepcopy(self.rin.family)
        self.hashes=dict(self.rin.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.old=copy.deepcopy(self.rin.leading);self.parameters=dict(self.rin.parameters)
        self.coordinate=self.rin.coordinate;self.rho=self.rin.rho;self.radius=self.rin.radius
        self.core_rho=self.ref(self.rin.inner['coordinates']['core'])
        p=self.parameters;eps,hb,delta=p['epsilon'],p['hb'],p['delta']
        a=outer.Algebra(g,5);b=outer.Algebra(g,3);qalg=outer.Algebra(g,4)
        self.core_templates={key:copy.deepcopy(g.nodes[ids[0][0]]) for key,ids in self.old['core_radial_rows'].items()}
        self.core_radial={key:[self.jet(ids) for ids in data]+[self.core_jet(key,n) for n in (3,4)]
            for key,data in self.old['core_radial_rows'].items()}
        self.collar_core={key:[a.at(jet,self.core_rho,self.rho) for jet in jets] for key,jets in self.core_radial.items()}
        flat=previous.cutoff_source.own.c2.flat_source
        self.sigma_pin=parent.active.current.ast_binding(flat.sigma_jets)
        self.sigma_interior_pin=parent.active.current.ast_binding(flat._sigma_left)
        for pin in (self.sigma_pin,self.sigma_interior_pin):
            if self.hashes[pin['module']]!=pin['module_sha256']:raise ValueError('Same original sigma source required')
        sigma0=previous.cutoff_source.own.jets.complete.Gamma.sigma(g,self.coordinate)
        sigma1=g.node('analytic_unary',name='original_flat_sigma_prime',argument=self.coordinate.node,
            derivative_source=self.sigma_pin,ordinary_derivative='same original sigma Taylor coefficient[1]',
            exact_flat_source_derivative_not_a_bound=True)
        sigma2=g.node('analytic_unary',name='original_flat_sigma_second',argument=self.coordinate.node,
            derivative_source=self.sigma_pin,ordinary_derivative='2*same original sigma Taylor coefficient[2]',
            Taylor_to_ordinary_factor=2,exact_flat_source_derivative_not_a_bound=True)
        sigma3=high.sigma_high(g,self.coordinate,3,self.sigma_interior_pin)
        self.sigma=[sigma0,sigma1,sigma2,sigma3]
        div=lambda x,y:g.quotient(x,y,'same positive original collar hb and rho factors')
        hm1=g.sub(hb,g.one);r2=g.mul(self.rho,self.rho);h2=g.mul(hb,hb)
        self.chi=[a.fixed(self.ref(self.rin.lower['fields']['chi'])),a.fixed(self.ref(self.rin.lower['fields']['chi_rho'])),
            a.fixed(div(g.mul(hm1,g.sub(sigma2,g.mul(hb,sigma1))),g.mul(h2,r2))),
            a.fixed(div(g.mul(hm1,g.add(sigma3,g.mul(g.constant(-3),hb,sigma2),g.mul(g.constant(2),h2,sigma1))),g.mul(h2,hb,r2,self.rho)))]
        cp=self.collar_core['Phi']
        self.A=[a.div(cp[1],cp[0],'same nonzero analytic core Phi')]
        self.A.append(a.div(a.sub(cp[2],a.mul(self.A[0],cp[1])),cp[0],'same nonzero analytic core Phi'))
        self.A.append(a.div(a.sub(cp[3],a.add(a.mul(self.A[0],cp[2]),a.scale(a.mul(self.A[1],cp[1]),2))),cp[0],'same nonzero analytic core Phi'))
        cm1=[a.sub(self.chi[0],a.fixed(1)),*self.chi[1:]]
        self.g_rows=[product_row(a,cm1,self.A,n) for n in range(3)]
        self.ratio=[self.jet(self.old['fields']['ratio'])]
        for n in range(1,4):self.ratio.append(product_row(a,self.ratio,self.g_rows,n-1))
        self.B=[product_row(a,self.ratio,self.chi,n) for n in range(4)]
        self.reconstructed={};self.radial={}
        for name in ('Phi','Psi'):
            old=[self.jet(ids) for ids in self.old['radial_fields'][name+'_rho_rows']]
            calculated=[old[0]]+[a.add(*(a.scale(a.mul(self.B[k],self.collar_core[name][n-k]),math.comb(n-1,k)) for k in range(n))) for n in range(1,5)]
            self.reconstructed[name]=calculated;self.radial[name]=old+calculated[3:]
        old_mean=[self.jet(ids) for ids in self.old['radial_fields']['mean_Psi_rho_rows']]
        means=[old_mean[0]]
        for n in range(1,5):means.append(a.scale(a.sub(self.radial['Psi'][n-1],a.scale(means[-1],n)),div(g.one,self.rho)))
        self.reconstructed['mean_Psi']=means;self.radial['mean_Psi']=old_mean+means[3:]
        for key,physical in (('Phi','F'),('Psi','U'),('mean_Psi','m')):
            old=[self.jet(ids) for ids in self.old['radial_fields'][physical+'_rho_rows']]
            if physical=='F':extra=[a.mul(self.jet(self.old['core_functions']['F0']),jet) for jet in self.radial[key][3:]]
            else:extra=[a.scale(jet,eps) for jet in self.radial[key][3:]]
            self.radial[physical]=old+extra
        self.recovered_Q=[rc_source.radial.leading_Q(g,U,m,delta) for U,m in zip(self.radial['U'],self.radial['m'])]
        self.radial['Q']=[self.jet(ids) for ids in self.old['radial_fields']['Q_rho_rows']]+self.recovered_Q[3:]
        z=outer.FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero));one=b.fixed(1);de=b.fixed(delta)
        d=b.sub(one,b.mul(z,z));h=b.sub(one,de);L=b.sub(one,b.mul(de,b.mul(z,z)))
        Q=[b.trim(jet) for jet in self.radial['Q']];Qz=[b.shift(jet) for jet in self.radial['Q']]
        U=[b.trim(jet) for jet in self.radial['U']]
        y=[b.add(b.scale(Q[n+1],self.rho),b.scale(Q[n],n)) for n in range(3)]
        C=[b.add(b.scale(b.mul(z,Q[n]),-2),b.mul(d,Qz[n]),b.scale(b.mul(z,y[n]),-2)) for n in range(3)]
        self.S=[];self.N_reconstructed=[]
        for n in range(3):
            linear=b.div(b.add(Q[n],b.scale(b.mul(h,b.mul(z,Qz[n])),'1/2'),y[n]),L,'same positive leading L')
            quadratic=b.add(b.scale(product_row(b,Q,Q,n),'1/2'),product_row(b,Q,y,n))
            axial=b.div(product_row(b,U,C,n),L,'same positive leading L')
            self.S.append(b.add(linear,quadratic,axial))
            self.N_reconstructed.append(b.add(b.scale(self.S[-1],'-1/2'),b.scale(b.add(b.scale(Q[n+1],n+2),b.scale(Q[n+2],self.rho)),p['Lambda'])))
        self.N=[self.jet(self.old['N1p']),*self.N_reconstructed[1:]]
        N,Nr,Nrr=self.N
        Py=b.scale(N,self.radius)
        Pyy=b.scale(b.add(N,b.scale(Nr,self.rho)),self.radius)
        Pyyy=b.scale(b.add(N,b.scale(Nr,g.mul(g.constant(3),self.rho)),b.scale(Nrr,r2)),self.radius)
        self.pressure_y=[Py,Pyy,Pyyy]
        negative_R2=g.neg(g.mul(self.radius,self.radius))
        self.m15_y=[b.scale(N,negative_R2),b.scale(b.add(b.scale(N,2),b.scale(Nr,self.rho)),negative_R2),
            b.scale(b.add(b.scale(N,4),b.scale(Nr,g.mul(g.constant(5),self.rho)),b.scale(Nrr,r2)),negative_R2)]
        self.signed_bodies=dict(m13=[b.scale(N,g.mul(eps,hb,self.rho)),
            b.scale(b.add(N,b.scale(Nr,self.rho)),g.mul(eps,hb)),
            b.scale(b.add(b.scale(Nr,2),b.scale(Nrr,self.rho)),g.mul(eps,hb))],
            m15=[b.scale(N,g.neg(g.mul(eps,eps,hb,r2))),
                b.scale(b.add(b.scale(N,g.mul(g.constant(2),self.rho)),b.scale(Nr,r2)),g.neg(g.mul(eps,eps,hb))),
                b.scale(b.add(b.scale(N,2),b.scale(Nr,g.mul(g.constant(4),self.rho)),b.scale(Nrr,r2)),g.neg(g.mul(eps,eps,hb)))])
        self.flat_domain=[g.mul(g.constant('3/4'),self.rin.endpoint),self.rin.endpoint]
        self.left=self.extend_own(self.rin.left);self.right=self.extend_own(self.rin.right)
        self.join={key:self.rin.at(packet,self.rin.endpoint) for key,packet in (('left',self.left),('right',self.right))}
        self.contract=copy.deepcopy(CONTRACT)
        self.definition=dict(actual_analytic_core_radial_orders_0_to_4_ordinary_Z0_to_Z5=True,
            actual_lower_Phi_Psi_mean_F_U_m_radial4_C5_and_Q_radial4_C4=True,
            old_radial0_1_2_and_all_nonzero_primitives_preserved=True,
            ratio_chi_and_Psi_transport_ODE_not_reset=True,
            pressure_N_rho0_rho1_rho2_C3_from_same_leading_one_minus_delta_source=True,
            high_own_left_scope_is_original_flat_cutoff_neighborhood=True,
            physical_y_derivatives_and_m15_primitive_order_distinct=True,
            full_finite_Rin_246_row_contract=True,numerical_global_n1_other_joins_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def guard_core(self,row):
        key=row.get('quantity')
        if (key not in self.core_templates or core_definition(row)!=core_definition(self.core_templates[key])
            or row.get('source_family')!=self.family or not 0<=row.get('radial_order',-1)<=4
            or not 0<=row.get('ordinary_slow_Z_derivative_order',-1)<=5):raise ValueError('Unadmitted analytic radial core component')
    def core_jet(self,key,n):
        jets=[]
        for k in range(6):
            row=copy.deepcopy(self.core_templates[key]);row.update(radial_order=n,ordinary_slow_Z_derivative_order=k)
            self.guard_core(row);op=row.pop('operation');jets.append(self.g.node(op,**row))
        return outer.FunctionJet(tuple(jets))
    def extend_own(self,packet):
        g=self.g;a=outer.Algebra(g,3);v=outer.Algebra(g,2)
        fields={key:list(values) for key,values in packet['fields'].items()}
        for key in rc_source.previous.STATE:
            alg=v if key in ('Q1','V1') else a
            fields[key]+=[alg.scale(fields[key][0],(-1)**n) if key in ('mean1','K1','Q1') else alg.fixed(0) for n in (3,4)]
        fields['P1']+=self.pressure_y[1:]
        moments={key:list(values) for key,values in packet['cumulative'].items()}
        for key in moments:moments[key]+=[a.fixed(0),a.fixed(0)]
        moments['m13']=list(fields['P1']);moments['m15'][2:]=self.m15_y[1:]
        return dict(fields=fields,cumulative=moments)
    def _snapshot(self):return outer.encode(dict(old=self.old,core_templates=self.core_templates,
        core_radial=self.core_radial,collar_core=self.collar_core,sigma=self.sigma,chi=self.chi,A=self.A,g_rows=self.g_rows,
        ratio=self.ratio,B=self.B,reconstructed=self.reconstructed,radial=self.radial,recovered_Q=self.recovered_Q,
        S=self.S,N_reconstructed=self.N_reconstructed,N=self.N,pressure_y=self.pressure_y,m15_y=self.m15_y,
        signed_bodies=self.signed_bodies,flat_domain=self.flat_domain,left=self.left,right=self.right,join=self.join,
        contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.rin.assert_graph()
        if (self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix or self.hashes!=self._hashes
            or self.numerical_oracle_available):raise ValueError('Actual higher lower/Rin source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed higher lower/Rin source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['lower_higher_radial_Rin_bindings']!=self.snapshot or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual higher lower/Rin source required')
        self.acceptance_loaded=True


def run(rin=None):
    began=time.monotonic();owner=CurrentOriginalLowerHigherRadialRinJoin(rin,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        lower_higher_radial_Rin_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_LOWER_RADIAL4_PRESSURE_RHO2_RIN_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
