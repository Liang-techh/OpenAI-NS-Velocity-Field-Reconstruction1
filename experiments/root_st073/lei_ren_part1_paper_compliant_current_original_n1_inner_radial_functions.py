"""True n1 even radial jets from the complete raw equation and axis limits.

F1/Uz1 have ordinary rho derivatives through order two. Their own mean,
K1 and radial velocity share the same regular continuation. Pressure has
its exact first derivative only. No odd rho jet at the axis is invented.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions as inner

parent,source,controls,c3,own=inner.parent,inner.source,inner.controls,inner.c3,inner.own
HERE,PREFIX,sha=inner.HERE,inner.PREFIX,inner.sha
NAME=PREFIX+'current_original_n1_inner_radial_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_inner_radial_functions_check.json'
GATE='current_original_n1_true_inner_even_radial_functions_installed'
OPEN=inner.OPEN


class CoreRhoDerivative(inner.AnalyticZDerivative):
    """Differentiate actual core components and arithmetic at fixed Z."""
    def __init__(self,g,rho,core_functions=None):
        self.g=g;self.z=rho
        self.derivatives={g.zero.node:g.zero,g.one.node:g.zero,rho.node:g.one}
        for jet in (core_functions or {}).values():
            for ref in c3.rows(jet) if type(jet) is c3.C3Function else own.c2.rows(jet):
                # These axis/gauge functions are independent of rho.
                self.derivatives[ref.node]=g.zero

    def operands(self,i):
        if self.g.nodes[i]['operation']=='current_original_unique_analytic_core_component':return []
        return super().operands(i)

    def rule(self,i):
        row=self.g.nodes[i]
        if row['operation']=='current_original_unique_regular_n1_component':
            raise ValueError('n1 radial derivatives must come from the full equation')
        if row['operation']=='current_original_unique_analytic_core_component':
            if row['rho']!=self.z.node:raise ValueError('Different core radial coordinate')
            definition=copy.deepcopy(row);operation=definition.pop('operation')
            definition['radial_order']+=1
            return self.g.node(operation,**definition)
        return super().rule(i)


def equation_velocity_rows(g,system,state,derivative):
    """Even-component rho1/rho2 on rho>0, from the full six-row system."""
    a=c3.C3Algebra(g);r,x=system['rho'],system['x']
    H={row:system['g'][row] for row in ('5','6')}
    for label in ('B0','B1'):
        for key,B in system[label].items():
            row,column=key.split(',')
            if row not in H:continue
            q=state[inner.FIELDS[int(column)-1]]
            if label=='B1':q=derivative.c3(q.Z)
            H[row]=a.add(H[row],a.mul(B,q))
    f,u=state['d_x_F1'],state['d_x_Uz1']
    div=lambda q,d:a.div(q,a.fixed(d))
    f1=div(f,g.mul(g.constant(2),x));u1=div(u,g.mul(g.constant(2),x))
    f2=a.add(div(H['5'],g.mul(g.constant(4),r)),a.neg(div(f,g.mul(r,x))))
    u2=a.add(div(H['6'],g.mul(g.constant(4),r)),a.neg(div(u,g.mul(g.constant(2),r,x))))
    return dict(F1=[state['F1'],f1,f2],Uz1=[state['Uz1'],u1,u2]),H


def regular_axis_rows(g,system,leading,parameters,rho_derivative,Z_derivative):
    """Exact rho1/rho2 axis coefficients of the same regular solution.

    All even values vanish at rho=0. The row-5 and row-6 rho coefficients
    yield 24*f2 and 16*u2; ordinary second derivatives are 2*f2/2*u2.
    """
    a=c3.C3Algebra(g);r=system['rho'];eps=parameters['epsilon'];de=parameters['delta']
    at=lambda q:inner.collar.at(g,q,r,g.zero)
    scale=lambda q,s:a.scale(q,g.constant(s))
    f1=scale(at(system['g']['5']),'1/8');u1=scale(at(system['g']['6']),'1/4')
    p1=a.scale(at(system['N1p']),eps);k1=scale(u1,'-1/2')
    slopes=dict(F1=f1,Uz1=u1,K1=k1,P1=p1)
    z=c3.C3Function(g.symbol('Z'),g.one,g.zero,g.zero);L=a.add(a.fixed(1),a.neg(a.scale(a.mul(z,z),de)))
    B=a.add(leading['Q_rho_rows'][0] if type(leading['Q_rho_rows'][0]) is c3.C3Function else
        Z_derivative.c3(leading['Q_rho_rows'][0].value),
        a.div(a.add(a.fixed(1),a.neg(scale(a.mul(z,leading['U_rho_rows'][0]),2))),L))
    beta=at(a.scale(B,eps))
    coefficient={row:at(Z_derivative.c3(rho_derivative.derivative(system['g'][row].value))) for row in ('5','6')}
    for label in ('B0','B1'):
        for key,Bij in system[label].items():
            row,column=key.split(',');j=int(column)-1
            if row not in coefficient or j>=4:continue
            value=slopes[inner.FIELDS[j]]
            if label=='B1':value=Z_derivative.c3(value.Z)
            coefficient[row]=a.add(coefficient[row],a.mul(at(Bij),value))
    coefficient['5']=a.add(coefficient['5'],scale(a.mul(beta,f1),2))
    coefficient['6']=a.add(coefficient['6'],scale(a.mul(beta,u1),2))
    f2=scale(coefficient['5'],'1/12');u2=scale(coefficient['6'],'1/8')
    zero=a.fixed(0)
    return dict(F1=[zero,f1,f2],Uz1=[zero,u1,u2],
        mean1=[zero,scale(u1,'1/2'),scale(u2,'1/3')],
        K1=[zero,k1,scale(u2,'-2/3')],P1_rho=p1),dict(
        row5_rho_coefficient=coefficient['5'],row6_rho_coefficient=coefficient['6'],
        diagonal_B0_over_x_axis=beta,source_radial_derivatives_are_actual_functions=True,
        ordinary_second_derivative_factors=dict(F1='1/12',Uz1='1/8'))


def own_mean_rows(g,state,velocity_rows,rho):
    a=c3.C3Algebra(g);U=velocity_rows['Uz1'];mean=a.add(state['Uz1'],state['K1'])
    m1=a.div(a.add(U[0],a.neg(mean)),a.fixed(rho))
    m2=a.div(a.add(U[1],a.neg(a.scale(m1,2))),a.fixed(rho))
    return dict(mean1=[mean,m1,m2],K1=[state['K1'],a.add(m1,a.neg(U[1])),a.add(m2,a.neg(U[2]))])


def axis_extend(g,interior,axis,rho,quantity,order):
    """Explicit regular continuation: the interior formula and exact limit."""
    return c3.C3Function(*(g.node('regular_even_radial_axis_continuation',
        interior_expression=q.node,axis_expression=p.node,rho=rho.node,
        quantity=quantity,ordinary_rho_derivative_order=order,
        ordinary_slow_Z_derivative_order=k,interior_domain='rho>0',axis_domain='rho=0',
        axis_value_is_exact_equation_limit=True,numerical_enclosure_not_selected=True)
        for k,(q,p) in enumerate(zip(c3.rows(interior),c3.rows(axis)))))


def recovered_radial_rows(g,rows,rho,epsilon,delta):
    """Q/V jets from the same own U and mean, with n1 temporal weight."""
    out=[]
    for U,M in zip(rows['Uz1'],rows['mean1']):
        state=dict(Uz1=U,K1=c3.C3Algebra(g).add(M,c3.C3Algebra(g).neg(U)))
        out.append(inner.recovery(g,state,dict(rho=rho,epsilon=epsilon),delta)['Q1'])
    b=own.c2.C2Algebra(g);R=g.mul(epsilon,rho)
    V=[b.scale(out[0],R),b.add(b.scale(out[0],epsilon),b.scale(out[1],R)),
        b.add(b.scale(out[1],g.mul(g.constant(2),epsilon)),b.scale(out[2],R))]
    return dict(Q1=out,V1=V)


def pressure_radial_function(g,leadingF,F1,N1p,epsilon):
    a=c3.C3Algebra(g)
    return a.scale(a.add(a.scale(a.mul(leadingF,F1),2),N1p),epsilon)


class CurrentOriginalN1InnerRadialFunctions:
    def __init__(self,inner_owner=None,require_checked=True):
        self.inner=inner_owner or inner.CurrentOriginalN1InnerVolterraFunctions()
        if not self.inner.acceptance_loaded:raise ValueError('Checked full original n1 functions required')
        self.inner.assert_graph();self.g=g=self.inner.g;self.family=copy.deepcopy(self.inner.family)
        self.hashes=dict(self.inner.hashes)
        for name,digest in parent.linear.read(inner.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (inner.NAME,inner.RECEIPT):self._pin(name,sha(name))
        self.parameters=self.inner.parameters;self.derivative=self.inner.derivative
        core=self.inner.collar.core
        fixed={key:core.functions[key] for key in ('U0','gradient','chi','G','logF0','F0','P0')}
        rd=CoreRhoDerivative(g,core.rho,fixed)
        self.axis,self.axis_binding=regular_axis_rows(g,self.inner.systems['core'],self.inner.leading['core'],
            self.parameters,rd,self.derivative)
        self.rows={};self.interior_rows={};self.equation_bindings={};self.recovered={};self.pressure_rho={}
        for chart,system in self.inner.systems.items():
            state=self.inner.functions[chart];r=system['rho'];a=c3.C3Algebra(g)
            velocity,H=equation_velocity_rows(g,system,state,self.derivative)
            rows={**velocity,**own_mean_rows(g,state,velocity,r)}
            self.interior_rows[chart]=rows;self.equation_bindings[chart]=H
            if chart=='core':
                rows={key:[jet if i==0 else axis_extend(g,jet,self.axis[key][i],r,key,i)
                    for i,jet in enumerate(value)] for key,value in rows.items()}
            self.rows[chart]=rows
            self.recovered[chart]=recovered_radial_rows(g,rows,r,self.parameters['epsilon'],self.parameters['delta'])
            self.pressure_rho[chart]=pressure_radial_function(g,self.inner.leading[chart]['F_rho_rows'][0],
                state['F1'],system['N1p'],self.parameters['epsilon'])
        self.definition=dict(full_raw_B0_and_B1_dZ_retained=True,
            ordinary_even_velocity_and_mean_K_rho_orders=[0,1,2],ordinary_Q_V_rho_orders=[0,1,2],
            ordinary_Q_V_axial_orders=[0,1,2],ordinary_even_axial_orders=[0,1,2,3],
            pressure_rho_order=1,pressure_second_radial_function_installed=False,
            odd_d_x_component_rho_derivatives_at_axis_not_reported=True,
            axis_limits_derived_from_regular_raw_equation=True,
            core_first_seam_inherits_one_regular_solution=True,
            analytic_higher_mixed_source_derivatives_are_definitions_not_numeric_oracles=True,
            original_numeric_n1_solution_oracle_executed=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw_data,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or raw_data['exact_graph_nodes']!=g.nodes
                or raw_data['inner_radial_bindings']!=self.snapshot):raise ValueError('Checked true n1 radial functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed true n1 radial source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed radial dependency '+name)
        self.hashes[name]=digest

    def _snapshot(self):return own.encode(dict(rows=self.rows,interior_rows=self.interior_rows,
        axis=self.axis,axis_binding=self.axis_binding,recovered=self.recovered,
        pressure_rho=self.pressure_rho,equation_bindings=self.equation_bindings,definition=self.definition))

    def assert_graph(self):
        self.inner.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('True n1 radial graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1InnerRadialFunctions(require_checked=False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,inner_radial_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys((GATE,)+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_TRUE_INNER_RADIAL_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':run()
