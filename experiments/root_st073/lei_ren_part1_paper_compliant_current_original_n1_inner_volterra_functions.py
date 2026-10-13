"""Actual raw n1 six-row system and its complete regular Volterra functions.

The Green bodies, all B0/B1 entries and finite words are executable DAGs.
Solution nodes denote the unique convergent complete series, never a center
of an old prefix/tail enclosure. Raw units and the zero-axis extension stay
separate from the normalized potential kernel.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions as collar
import lei_ren_part1_paper_compliant_current_original_n1_inner_operator as raw
import lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step as coupled
import lei_ren_part1_paper_compliant_current_original_n1_regular_series_function as regular

parent,source,controls,c3,own=collar.parent,collar.source,collar.controls,collar.c3,collar.own
HERE,PREFIX,sha=collar.HERE,collar.PREFIX,collar.sha
NAME=PREFIX+'current_original_n1_inner_volterra_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_inner_volterra_functions_check.json'
GATE='current_original_n1_true_inner_Volterra_function_adapter_installed'
OPEN=collar.OPEN
FIELDS=regular.FIELDS
D=(0,0,2,0,3,1)
POWERS=(2,2,2,2,1,1)


class AnalyticZDerivative(own.jets.OrdinaryZDerivative):
    """Differentiate the actual analytic limits, not numerical enclosures."""
    COMPONENTS=('current_original_unique_analytic_core_component','current_original_unique_regular_n1_component')
    def operands(self,i):
        if self.g.nodes[i]['operation'] in self.COMPONENTS:return []
        return super().operands(i)
    def rule(self,i):
        row=self.g.nodes[i]
        if row['operation'] in self.COMPONENTS:
            definition=copy.deepcopy(row);operation=definition.pop('operation')
            definition['ordinary_slow_Z_derivative_order']+=1
            return self.g.node(operation,**definition)
        return super().rule(i)
    def c3(self,value):
        out=[value]
        for _ in range(3):out.append(self.derivative(out[-1]))
        return c3.C3Function(*out)


def raw_matrix_expressions():
    names='rho x Z delta epsilon F Fr FZ U Ur UZ Q'.split()
    v={name:sy.Symbol(name) for name in names}
    r,x,z,de,eps,f,fr,fz,w,wr,wz,q=[v[name] for name in names]
    d=1-z*z;L=1-de*z*z
    H=z*(1-de)/2+d*w;B=q+(1-2*z*w)/L;E=z*w-sy.Rational(1,2)
    Cf=f+r*fr;Zf=((-2-de)*z*f+d*fz-2*z*r*fr)/L
    Zw=((-1-de)*z*w+d*wz-2*z*r*wr)/L
    B0={'1,5':sy.Integer(1),'2,6':sy.Integer(1),'3,6':sy.Integer(-1),
        '4,1':eps*4*x*f,
        '5,1':eps*2*(q+E*(-2+de)/L),
        '5,2':eps*(2*Zf-2*(-1+de)*z*Cf/L),
        '5,3':-eps*2*(1+de)*z*Cf/L,
        '5,5':eps*x*B,
        '6,1':-eps**2*8*z*r*f/L,
        '6,2':eps*2*(E*(-1+de)/L+Zw-wr*z*(-1+de)*r/L),
        '6,3':-eps*2*(1+de)*z*r*wr/L,'6,4':-eps*4*z/L,'6,6':eps*x*B}
    B1={'5,1':eps*2*H/L,'5,2':-eps*2*d*Cf/L,'5,3':-eps*2*d*Cf/L,
        '6,2':eps*2*(H-wr*d*r)/L,'6,3':-eps*2*d*r*wr/L,'6,4':eps*2*d/L}
    return dict(B0=B0,B1=B1)


def axial_second(g,rows,rho,delta,weight):
    """The actual Z_(weight-1+delta) Z_weight differential composition."""
    z=g.symbol('Z');c=g.constant;de=delta;d=g.sub(g.one,g.mul(z,z));L=g.sub(g.one,g.mul(de,z,z))
    f,fr,frr=rows;div=lambda n:g.quotient(n,L,'positive original L')
    A=g.add(g.mul(weight,z,f.value),g.mul(d,f.Z),g.neg(g.mul(c(2),z,rho,fr.value)))
    AZ=g.add(g.mul(weight,f.value),g.mul(g.sub(weight,c(2)),z,f.Z),g.mul(d,f.ZZ),
        g.neg(g.mul(c(2),rho,fr.value)),g.neg(g.mul(c(2),z,rho,fr.Z)))
    AR=g.add(g.mul(g.sub(weight,c(2)),z,fr.value),g.mul(d,fr.Z),g.neg(g.mul(c(2),z,rho,frr.value)))
    return div(g.add(div(g.mul(g.add(g.sub(weight,g.one),de),z,A)),
        g.mul(d,g.add(div(AZ),g.quotient(g.mul(c(2),de,z,A),g.mul(L,L),'positive L squared'))),
        g.neg(div(g.mul(c(2),z,rho,AR)))))


def raw_system(g,fields,parameters,derivative):
    a=parent.C1Algebra(g);c=g.constant;r=fields['rho'];de=parameters['delta'];eps=parameters['epsilon']
    F,U,Q=fields['F_rho_rows'],fields['U_rho_rows'],fields['Q_rho_rows']
    c1=lambda q:source.C1Function(q.value,q.Z)
    z=source.C1Function(g.symbol('Z'),g.one);dd=a.sub(a.fixed(1),a.mul(z,z))
    delta=a.fixed(de);L=a.sub(a.fixed(1),a.mul(delta,z,z));h=a.sub(a.fixed(1),delta)
    q,qz,u=c1(Q[0]),source.C1Function(Q[0].Z,Q[0].ZZ),c1(U[0])
    qy=a.scale(r,c1(Q[1]))
    S=a.add(a.div(a.add(q,a.scale('1/2',a.mul(h,z,qz)),qy),L),
        a.mul(q,a.add(a.scale('1/2',q),qy)),
        a.div(a.mul(u,a.add(a.scale(-2,a.mul(z,q)),a.mul(dd,qz),a.scale(-2,a.mul(z,qy)))),L))
    # This regular rho form is valid at the axis; no division by R=0.
    Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))
    Nt=g.neg(axial_second(g,F,r,de,g.sub(c(-2),de)))
    Nz=g.neg(axial_second(g,U,r,de,g.sub(c(-1),de)))
    x=g.unary('sqrt',r);zv=z.value;Lv=L.value
    values=dict(rho=r,x=x,Z=zv,delta=de,epsilon=eps,
        F=F[0].value,Fr=F[1].value,FZ=F[0].Z,U=U[0].value,Ur=U[1].value,UZ=U[0].Z,Q=Q[0].value)
    expressions=raw_matrix_expressions()
    matrices={key:{name:derivative.c3(source.expression(g,body,values,x)) for name,body in rows.items()}
        for key,rows in expressions.items()}
    g4=g.mul(c(2),eps,x,Np.value);g5=g.mul(c(2),eps,Nt)
    g6=g.sub(g.mul(c(2),eps,Nz),g.quotient(g.mul(c(4),zv,r,eps,eps,Np.value),Lv,'positive original L'))
    return dict(**matrices,C0={key:q for key,q in matrices['B0'].items() if key not in ('1,5','2,6','3,6')},
        g={key:derivative.c3(value) for key,value in (('4',g4),('5',g5),('6',g6))},
        N1theta=derivative.c3(Nt),N1z=derivative.c3(Nz),N1p=derivative.c3(Np.value),x=x,rho=r,D=D)


def core_fields(owner):
    g=owner.g;core=owner.core;a=c3.C3Algebra(g);eps=core.parameters['epsilon']
    rows=owner.core_radial_rows
    U=[a.add(core.functions['U0'],a.scale(rows['Psi'][0],eps)),*[a.scale(q,eps) for q in rows['Psi'][1:]]]
    mean=[a.add(core.functions['U0'],a.scale(rows['mean_Psi'][0],eps)),*[a.scale(q,eps) for q in rows['mean_Psi'][1:]]]
    return dict(rho=core.rho,query=core.rho,F_rho_rows=[a.mul(core.functions['F0'],q) for q in rows['Phi']],
        U_rho_rows=U,m_rho_rows=mean,Q_rho_rows=[own.radial_Q(g,core.parameters['delta'],u,m) for u,m in zip(U,mean)])


def green_apply(g,vectors,chart,coordinates,parameters):
    """(I+J)G on both charts, with one unrestarted core-to-collar path."""
    c=g.constant;a=c3.C3Algebra(g);hb=parameters['hb'];r=coordinates['core'];s=coordinates['first']
    t=g.symbol('original_n1_Green_t');v=g.symbol('original_n1_Green_phase_v')
    integral=lambda q,var,lo,hi: c3.C3Function(*(controls.integral(g,row,g.nodes[var.node]['name'],lo,hi,
        actual_raw_regular_n1_Green_body=True,ordinary_slow_Z_derivative_order=k,
        endpoint_independent_of_Z=True,zero_endpoint_integrand_continuation=0,
        complete_history_not_reset_at_core_seam=True,signed_body_not_magnitude_cap=True)
        for k,row in enumerate(c3.rows(q))))
    if chart=='core':
        X=g.unary('sqrt',r);sample=g.mul(r,t,t);measure=X
        K=g.mul(c('1/2'),X,t,g.sub(g.one,g.mul(t,t)))
        log=g.mul(X,t,g.unary('log',g.quotient(g.one,t,'positive t inside unit integral')))
        kernels=(K,log,g.neg(K),g.one,g.mul(t,t,t),t)
    elif chart=='first':
        X=g.mul(c(2),g.unary('exp',g.mul(c('1/2'),hb,s)));sample=g.mul(c(4),t,t);measure=c(2)
        distance=g.mul(hb,s);tail=g.mul(distance,g.unary('exprel',g.neg(distance)))
        K=g.mul(t,g.add(g.sub(g.one,g.mul(t,t)),g.mul(t,t,tail)))
        log=g.mul(c(2),t,g.add(g.mul(c('1/2'),hb,s),g.unary('log',g.quotient(g.one,t,'positive t inside unit integral'))))
        ratio=g.mul(t,g.unary('exp',g.neg(g.mul(c('1/2'),hb,s))))
        kernels=(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)
    else:raise ValueError('Regular n1 owns only core/first charts')
    out=[]
    for row,kernel in zip(('5','6','6','4','5','6'),kernels):
        body=collar.at(g,vectors['core'][row],r,sample)
        out.append(integral(a.scale(body,g.mul(kernel,measure)),t,g.zero,g.one))
    if chart=='first':
        x=g.mul(c(2),g.unary('exp',g.mul(c('1/2'),hb,v)))
        distance=g.mul(hb,g.sub(s,v));ratio=g.unary('exp',g.neg(g.mul(c('1/2'),distance)))
        K=g.mul(c('1/2'),x,distance,g.unary('exprel',g.neg(distance)))
        log=g.mul(c('1/2'),x,distance);measure=g.mul(c('1/2'),hb,x)
        kernels=(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)
        for i,(row,kernel) in enumerate(zip(('5','6','6','4','5','6'),kernels)):
            body=collar.at(g,vectors['first'][row],s,v)
            out[i]=a.add(out[i],integral(a.scale(body,g.mul(kernel,measure)),v,g.zero,s))
    return dict(zip(FIELDS,out))


def feedback_vectors(g,systems,state,derivative):
    a=c3.C3Algebra(g);out={}
    for chart,system in systems.items():
        vector={key:a.fixed(0) for key in ('4','5','6')}
        for label in ('C0','B1'):
            for key,B in system[label].items():
                row,column=key.split(',');q=state[chart][FIELDS[int(column)-1]]
                if label=='B1':q=derivative.c3(q.Z)
                vector[row]=a.add(vector[row],a.mul(B,q))
        out[chart]=vector
    return out


def recovery(g,state,rho,delta):
    a=c3.C3Algebra(g);b=own.c2.C2Algebra(g);z=own.c2.C2Function(g.symbol('Z'),g.one,g.zero)
    U,K=state['Uz1'],state['K1'];mean=a.add(U,K);c2=lambda q:own.c2.C2Function(q.value,q.Z,q.ZZ)
    d=b.add(b.fixed(1),b.neg(b.mul(z,z)));de=b.fixed(delta)
    L=b.add(b.fixed(1),b.neg(b.mul(de,b.mul(z,z))))
    mz=own.c2.C2Function(mean.Z,mean.ZZ,mean.ZZZ)
    Q=b.div(b.add(b.mul(b.add(b.fixed(1),b.neg(de)),b.mul(z,c2(U))),
        b.neg(b.mul(b.add(b.fixed(1),de),b.mul(z,c2(K)))),b.neg(b.mul(d,mz))),L)
    R=g.mul(rho['epsilon'],rho['rho'])
    return dict(mean1=mean,radial_primitive_M1=a.scale(mean,R),Q1=Q,V1=b.scale(Q,R),
        temporal_n1_mean_weight_is_1_plus_delta=True,raw_physical_units_retained=True)


class CurrentOriginalN1InnerVolterraFunctions:
    def __init__(self,collar_owner=None,require_checked=True):
        self.collar=collar_owner or collar.CurrentOriginalN1LowerCollarFunctions()
        if not self.collar.acceptance_loaded:raise ValueError('Checked true core/collar functions required')
        self.collar.assert_graph();self.g=g=self.collar.g;self.family=copy.deepcopy(self.collar.family)
        self.hashes=dict(self.collar.hashes)
        for name,digest in parent.linear.read(collar.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (collar.NAME,collar.RECEIPT):self._pin(name,sha(name))
        for module in (raw,coupled,regular):
            receipt=parent.linear.read(module.RECEIPT)
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual raw n1 source/series prerequisite required')
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        packet=next(iter(parent.linear.read(coupled.NAME)['samples'].values()))
        self.tail_admission=copy.deepcopy(packet['analytic_preconditioned_tail'])
        tail=self.tail_admission;ctx=collar.core.common.MPIntervalContext();ctx.dps=500
        read=lambda row:collar.selection.read(ctx,row)
        if (not tail['all_k_derivative_words_allowed'] or tail['half_derivative_block_reduction_used']
            or tail['old_unpreconditioned_tail_reused'] or tail['retained_preconditioned_indices']!=[0,1]
            or not 0<=collar.core.common.endpoints(read(tail['geometric_tail_ratio_upper']))[0]
                <=collar.core.common.endpoints(read(tail['geometric_tail_ratio_upper']))[1]<1):
            raise ValueError('Complete raw preconditioned all-word convergence required')
        self.parameters=self.collar.parameters;self.coordinates=dict(core=self.collar.core.rho,first=self.collar.fields['query'])
        self.derivative=derivative=AnalyticZDerivative(g)
        # The anchored G primitive has a moving Z endpoint. Its exact FTC
        # jet is already accepted; do not re-differentiate it as a fixed-end
        # integral. Higher derivatives use the rational gradient body.
        for name in ('G','logF0','F0'):derivative.seed_jet(self.collar.core.functions[name])
        self.leading=dict(core=core_fields(self.collar),first=self.collar.fields)
        self.systems={chart:raw_system(g,fields,self.parameters,derivative) for chart,fields in self.leading.items()}
        operator_bindings=own.encode(self.systems);self.functions={}
        for chart,coordinate in self.coordinates.items():
            self.functions[chart]={name:c3.C3Function(*(g.node('current_original_unique_regular_n1_component',
                quantity=name,chart=chart,coordinate=coordinate.node,Z_variable=g.symbol('Z').node,
                coordinate_domain_lower=g.zero.node,
                coordinate_domain_upper=(g.constant(4).node if chart=='core' else self.parameters['sc_quarter'].node),
                ordinary_slow_Z_derivative_order=k,source_family=self.family,
                actual_raw_operator_bindings=operator_bindings,zero_axis_regular_power=POWERS[i],
                exact_function_definition='W=sum_k>=0 H^k W0; W0=(I+J)Gg; H=(I+J)G(C0+B1*dZ)',
                all_words_retained_in_analytic_limit=True,unique_regular_analytic_solution=True,
                raw_units_not_normalized_potential_kernel=True,source_parameter_bounds_not_selected_as_values=True,
                defining_Green_provider=parent.active.current.ast_binding(green_apply),
                defining_feedback_provider=parent.active.current.ast_binding(feedback_vectors),
                convergence_provider=parent.active.current.ast_binding(coupled.preconditioned_tail))
                for k in range(4))) for i,name in enumerate(FIELDS)}
        self.words={0:{chart:green_apply(g,{key:system['g'] for key,system in self.systems.items()},chart,self.coordinates,self.parameters)
            for chart in self.coordinates}}
        self.words[1]=self.apply(self.words[0])
        correction=self.apply(self.functions);a=c3.C3Algebra(g)
        self.fixed_point_rhs={chart:{name:a.add(self.words[0][chart][name],correction[chart][name])
            for name in FIELDS} for chart in self.coordinates}
        self.recovered={chart:recovery(g,state,dict(rho=self.systems[chart]['rho'],epsilon=self.parameters['epsilon']),self.parameters['delta'])
            for chart,state in self.functions.items()}
        self.definition=dict(D=D,regular_axis_powers=POWERS,kinematic_entries={'1,5':1,'2,6':1,'3,6':-1},
            raw_C0_is_B0_minus_Akin=True,one_core_to_collar_path_no_seam_reset=True,
            common_first_phase_upper=self.parameters['sc_quarter'],
            common_inner_endpoint_is_Rin_not_rminus=True,
            leading_source_pressure_domain_is_distinct_sc_quarter_to_sc_half=True,
            full_regular_solution_distinct_from_finite_words_and_tail_enclosures=True,
            normalized_Kchi_not_used_with_raw_matrix=True,numerical_original_n1_point_accuracy_claimed=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw_data,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or raw_data['exact_graph_nodes']!=g.nodes
                or raw_data['inner_Volterra_bindings']!=self.snapshot):raise ValueError('Checked true n1 Volterra functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed inner n1 function source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed n1 dependency '+name)
        self.hashes[name]=digest

    def apply(self,state):
        vectors=feedback_vectors(self.g,self.systems,state,self.derivative)
        return {chart:green_apply(self.g,vectors,chart,self.coordinates,self.parameters) for chart in self.coordinates}

    def finite_word(self,index):
        if isinstance(index,bool) or not isinstance(index,int) or not 0<=index<=8:raise ValueError('Finite executable word index 0..8 required')
        for i in range(2,index+1):
            if i not in self.words:self.words[i]=self.apply(self.words[i-1])
        return self.words[index]

    def _snapshot(self):return own.encode(dict(functions=self.functions,systems=self.systems,leading=self.leading,
        words={str(i):self.words[i] for i in (0,1)},fixed_point_rhs=self.fixed_point_rhs,recovered=self.recovered,
        parameters=self.parameters,coordinates=self.coordinates,definition=self.definition,tail_admission=self.tail_admission))

    def assert_graph(self):
        self.collar.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('True n1 Volterra graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1InnerVolterraFunctions(require_checked=False)
    raw_data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,inner_Volterra_bindings=owner.snapshot,
        actual_raw_six_row_operator_and_complete_regular_limit=True,
        **dict.fromkeys((GATE,)+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw_data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_TRUE_INNER_VOLTERRA_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw_data


if __name__=='__main__':run()
