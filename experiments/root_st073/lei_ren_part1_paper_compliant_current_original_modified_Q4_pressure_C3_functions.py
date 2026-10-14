"""True ten-chart modified generic C4, total radial C5 and pressure C3.

The reduced pressure integral is the original exact integration-by-parts
identity. It retains both physical finite endpoints and every chart history.
The actual Rm patch/six outer charts and quantitative mixed bounds stay open.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions as previous
import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check as lexical
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as original_pressure

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_modified_Q4_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_modified_Q4_pressure_C3_functions_check.json'
GATES=('current_original_true_ten_chart_modified_full_generic_C4_functions_installed',
       'current_original_true_ten_chart_modified_total_Uy_mean_C5_Qy_C4_functions_installed',
       'current_original_true_rminus_to_Rm_signed_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+('current_original_modified_Q4_pressure_C3_numeric_containment_installed',
    'current_original_true_ten_chart_modified_Qyy_C4_functions_installed',
    'current_original_complete_preRp_pressure_C3_functions_installed')))
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=previous.LoopAlgebra
CHARTS=previous.CHARTS


def generic_functions(g,fields,delta,Pstar,radius,order=4):
    """Unchanged full signed leading recovery with order-aware true inputs."""
    a=Algebra(g,order);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*(order-1))
    one=a.fixed(1);de=a.fixed(delta);d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    E,V,m,h,k,e=(a.trim(fields[name]) for name in ('E','V','m','h','k','e'))
    pressure=a.trim(fields['pressure'])
    mz,hz,kz,ez,pz=(a.shift(fields[name]) for name in ('m','h','k','e','pressure'))
    div=lambda q:a.div(q,L,'same original positive L=1-delta*Z^2')
    transport=a.add(a.mul(a.sub(one,de),a.mul(z,m)),a.mul(d,mz))
    Q=div(a.sub(a.scale(a.mul(z,V),2),transport))
    itl=div(a.sub(a.add(a.neg(E),a.mul(a.sub(one,a.scale(de,'1/2')),h)),
        a.scale(a.mul(a.sub(one,de),a.mul(z,hz)),'1/2')))
    itq=div(a.add(a.mul(a.sub(a.scale(de,2),one),a.mul(z,k)),a.neg(a.mul(d,kz)),a.mul(E,transport)))
    izl=div(a.add(a.neg(V),a.scale(a.mul(a.sub(one,de),a.sub(m,a.mul(z,mz))),'1/2')))
    izq=div(a.add(a.mul(V,transport),a.scale(a.mul(de,a.mul(z,e)),2),a.neg(a.mul(d,ez)),
        a.scale(a.mul(a.add(one,de),a.mul(z,pressure)),2),a.neg(a.mul(d,pz))))
    It=a.scale(a.add(itl,a.scale(itq,Pstar)),radius)
    Iz=a.scale(a.add(izl,a.scale(izq,Pstar)),radius)
    p2=a.div(Iz,E,'same actual positive modified swirl E_N')
    return dict(Q=Q,transport=transport,theta_linear=itl,theta_quadratic=itq,
        axial_linear=izl,axial_quadratic=izq,inertial_theta=It,inertial_axial=Iz,p2=p2)


def background_radial_rows(g,chart):
    """Original normalized leading moment ODEs, all seven genuine rows."""
    a=Algebra(g,6);f=chart['fields'];E,V=f['E'],f['V'];E2=a.mul(E,E)
    result=dict(E=chart['radial']['E_y'],V=chart['radial']['V_y'],
        m=a.sub(V,f['m']),h=a.sub(E,a.scale(f['h'],'3/2')),
        k=a.sub(a.mul(E,V),a.scale(f['k'],'3/2')),
        e=a.sub(a.sub(a.mul(V,V),a.scale(E2,'1/2')),f['e']),
        p=a.scale(E2,'1/2'),P0=a.fixed(0))
    result['pressure']=result['p'];return result


class NativeDerivative:
    """Exact partial differentiation, including moving integral endpoints.

    The issued fast phase is held fixed here. Its total-y contribution is
    added explicitly from the admitted A_phi/B_phi functions afterwards.
    """
    def __init__(self,g,variable,phase,shared=None):
        self.g,self.variable,self.phase=g,variable,phase
        self.name=g.nodes[variable.node]['name']
        self.shared=shared if shared is not None else dict(free={},factories={})
        self.shared['factories'][variable.node]=self
        self.cache={phase.node:g.zero};self.core_extensions=[]
    def other(self,variable):
        if variable.node not in self.shared['factories']:
            NativeDerivative(self.g,variable,self.phase,self.shared)
        return self.shared['factories'][variable.node]
    def seed(self,left,right):
        self.g.ids([left,right])
        if left.node in self.cache and self.cache[left.node]!=right:
            raise ValueError('Conflicting native source derivative')
        if left in (self.g.zero,self.g.one) and right!=self.g.zero:
            raise ValueError('A constant cannot have a nonzero radial derivative')
        self.cache[left.node]=right
    def sub(self,body,variable,value):
        return self.g.node('function_substitution',expression=body.node,variable=variable.node,value=value.node,
            Z_independent_substitution=True,exact_native_chain_rule_substitution=True)
    def sigma_prime(self,x):
        g=self.g;c=g.constant;sig=previous.previous.sigma(g,x);comp=g.sub(g.one,x)
        first=g.add(g.quotient(c(2),g.mul(x,x,x),'interior original sigma x>0'),
            g.quotient(c(2),g.mul(comp,comp,comp),'interior original sigma x<1'))
        predicate=g.node('logical_or',arguments=[
            g.node('exact_real_comparison',operator='<=',left=x.node,right=g.zero.node).node,
            g.node('exact_real_comparison',operator='>=',left=x.node,right=g.one.node).node])
        return g.node('original_lazy_flat_branch',active_body=g.mul(sig,g.sub(g.one,sig),first).node,
            Delta=x.node,eta=g.one.node,flat_value=g.zero.node,flat_predicate=predicate.node,
            active_body_not_evaluated_on_flat_branch=True,original_sigma_exact_logistic_first_derivative=True)
    def derivative(self,ref):
        g=self.g;g.ids([ref]);i=ref.node
        if i in self.cache:return self.cache[i]
        if self.name not in lexical.free_variables(g,ref,self.shared['free']):
            self.cache[i]=g.zero;return g.zero
        row=g.nodes[i];op=row['operation'];D=lambda j:self.derivative(source.FunctionRef(g,j))
        R=lambda j:source.FunctionRef(g,j)
        if op=='bound_variable':value=g.one if i==self.variable.node else g.zero
        elif op=='sum':value=g.add(*(D(j) for j in row['arguments']))
        elif op=='product':
            ids=row['arguments'];value=g.add(*(g.mul(*(D(v) if k==j else R(v) for k,v in enumerate(ids))) for j in range(len(ids))))
        elif op=='negative':value=g.neg(D(row['argument']))
        elif op=='positive_quotient':
            u,v=R(row['numerator']),R(row['denominator'])
            value=g.quotient(g.sub(g.mul(D(u.node),v),g.mul(u,D(v.node))),g.mul(v,v),row['source_positive_certificate'])
        elif op in ('analytic_unary','exact_original_source_flat_sigma'):
            x=R(row['argument']);dx=D(x.node);name='original_flat_sigma' if op=='exact_original_source_flat_sigma' else row['name']
            if name in ('exp','expm1'):outer=g.unary('exp',x)
            elif name=='log':outer=g.quotient(g.one,x,'same original positive logarithm input')
            elif name=='log1p':outer=g.quotient(g.one,g.add(g.one,x),'same original positive 1+x')
            elif name in ('sqrt','positive_sqrt'):outer=g.quotient(g.one,g.mul(g.constant(2),g.unary('positive_sqrt',x)),
                'same original active positive square root')
            elif name=='cos':outer=g.neg(g.unary('sin',x))
            elif name=='sin':outer=g.unary('cos',x)
            elif name=='original_flat_sigma':outer=self.sigma_prime(x)
            elif name=='original_flat_sigma_prime':
                value=self.derivative(self.sigma_prime(x));self.cache[i]=value;return value
            else:raise ValueError('Unissued native unary derivative '+name)
            value=g.mul(outer,dx)
        elif op=='function_substitution':
            body,var,arg=(R(row[key]) for key in ('expression','variable','value'))
            direct=g.zero if var==self.variable else self.sub(self.derivative(body),var,arg)
            chain=g.mul(self.sub(self.other(var).derivative(body),var,arg),self.derivative(arg))
            value=g.add(direct,chain)
        elif op=='definite_integral':
            body,lower,upper=(R(row[key]) for key in ('integrand','lower','upper'))
            var=g.symbol(row['variable']);value=g.zero
            if row['variable']!=self.name:
                partial=self.derivative(body)
                if partial!=g.zero:
                    attrs=copy.deepcopy(row);operation=attrs.pop('operation');attrs['integrand']=partial.node
                    attrs['exact_native_partial_of_same_integrand']=True
                    value=g.node(operation,**attrs)
            value=g.add(value,g.mul(self.sub(body,var,upper),self.derivative(upper)),
                g.neg(g.mul(self.sub(body,var,lower),self.derivative(lower))))
        elif op=='substitute_original_inverse_angle':
            body,inverse=R(row['body']),R(row['inverse_angle']);angle=g.symbol(row['angle_variable'])
            at=lambda q:previous.at_inverse(g,q,inverse,angle)
            direct=g.zero if angle==self.variable else at(self.derivative(body))
            value=g.add(direct,g.mul(at(self.other(angle).derivative(body)),self.derivative(inverse)))
        elif op=='original_monotone_phase_inverse':
            Phi=R(row['phase_function']);angle=g.symbol(row['angle_variable']);target=R(row['target_phase'])
            at=lambda q:previous.at_inverse(g,q,ref,angle)
            held=g.zero if angle==self.variable else self.derivative(Phi)
            lam=at(self.other(angle).derivative(Phi))
            value=g.quotient(g.sub(self.derivative(target),at(held)),lam,'original positive phase lambda; exact implicit native derivative')
        elif op=='original_lazy_flat_branch':
            attrs=copy.deepcopy(row);operation=attrs.pop('operation')
            attrs['active_body']=D(row['active_body']).node;attrs['flat_value']=D(row['flat_value']).node
            attrs['same_original_lazy_native_derivative']=True;value=g.node(operation,**attrs)
        elif op=='current_original_unique_analytic_core_component':
            if (row['source_family'] is None or (self.shared.get('family') is not None and row['source_family']!=self.shared['family'])
                or not row['unique_analytic_Banach_solution']
                or not row['actual_twenty_term_operator'] or not row['admitted_solution_not_a_model_or_enclosure_center']):
                raise ValueError('Native derivative requires the admitted true analytic core')
            dr,dz=D(row['rho']),D(row['Z_variable']);terms=[]
            for key,factor in (('radial_order',dr),('ordinary_slow_Z_derivative_order',dz)):
                if factor!=g.zero:
                    attrs=copy.deepcopy(row);operation=attrs.pop('operation');attrs[key]+=1
                    if attrs['radial_order']>2:raise ValueError('An additional true radial core source is required')
                    new=g.node(operation,**attrs);terms.append(g.mul(new,factor));self.core_extensions.append(new.node)
            value=g.add(*terms)
        else:raise ValueError('Unissued actual native derivative operation '+op)
        self.cache[i]=value;return value


def total_radial_functions(g,chart,derivative,N):
    a=Algebra(g,5);bg=chart['source_chart'];jac=chart['physical_measure'];invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')
    slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref)) for ref in rows(chart['loop'][key])))
          for key in ('A','B_over_Pstar')}
    invN=g.quotient(g.one,N,'same original exact positive common N')
    x=a.scale(chart['loop']['A'],invN);ex=a.exp(x)
    xy=a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])
    Ey=a.mul(ex,a.add(a.trim(bg['radial']['E_y']),a.mul(chart['roots']['E'],xy)))
    Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])
    mean_y=a.sub(chart['modified_fields']['V'],chart['modified_fields']['m'])
    return dict(E_y=Ey,V_y=Vy,m_y=mean_y,slow_A_y=slow['A'],slow_B_y=slow['B_over_Pstar'],
        total_fast_phase_terms_included=True,physical_y_equals_native_partial_divided_by_Jacobian=True)


def physical_pressure_fields(g,fields,radial,delta,Pstar):
    a=Algebra(g,5);b=Algebra(g,4);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*3)
    U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)
    m=a.scale(fields['m'],Pstar);my=a.sub(U,m)
    one=b.fixed(1);de=b.fixed(delta);d=b.sub(one,b.mul(z,z));L=b.sub(one,b.mul(de,b.mul(z,z)))
    div=lambda q:b.div(q,L,'same original L=1-delta*Z^2')
    numerator=lambda u,mean: b.sub(b.scale(b.mul(z,b.trim(u)),2),
        b.add(b.mul(b.sub(one,de),b.mul(z,b.trim(mean))),b.mul(d,b.shift(mean))))
    Q=div(numerator(U,m));Qy=div(numerator(Uy,my))
    return dict(U=U,U_y=Uy,m=m,m_y=my,Q=Q,Q_y=Qy)


def reduced_pressure_S(g,physical,delta):
    a=Algebra(g,3);z=FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero));one=a.fixed(1);de=a.fixed(delta)
    d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    U,Q,Qy=(a.trim(physical[key]) for key in ('U','Q','Q_y'));Qz=a.shift(physical['Q'])
    div=lambda q:a.div(q,L,'same original positive meridional L')
    linear=div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy))
    quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))
    axial=div(a.mul(U,a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2))))
    return a.add(linear,quadratic,axial)


class CurrentOriginalModifiedQ4PressureC3Functions:
    def __init__(self,loop_owner,require_checked=True):
        if type(loop_owner) is not previous.CurrentOriginalPrescribedLoopC5Functions or not loop_owner.acceptance_loaded:
            raise TypeError('Checked true ten-chart loop/density/history C5 source required')
        loop_owner.assert_graph();self.loop=loop_owner;self.background=loop_owner.background;self.g=g=loop_owner.g
        self.family=copy.deepcopy(loop_owner.family);self.hashes=dict(loop_owner.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (previous.NAME,previous.RECEIPT,Path(lexical.__file__).name,Path(original_pressure.__file__).name):self._pin(name,sha(name))
        self.bindings=dict(original_generic=previous.previous.previous.inlet.generic.source_bindings(),
            original_pressure_reduction=parent.active.current.ast_binding(original_pressure.pressure_functions),
            original_endpoint_reduction=parent.active.current.ast_binding(original_pressure.q_plus_q_y))
        self.charts={};self.derivatives={};self.endpoints={};a=Algebra(g,3)
        self.initial=None;incoming=dict(S0=a.fixed(0),S1=a.fixed(0));self.zero_pressure_integrals=incoming
        for name in CHARTS:
            chart=loop_owner.charts[name];bg=chart['source_chart'];coordinate=chart['coordinate'];jac=chart['physical_measure']
            derivative=NativeDerivative(g,coordinate,chart['actual_phase']);self.derivatives[name]=derivative
            derivative.shared['family']=self.family
            radial_bg=background_radial_rows(g,bg)
            for key,jet in bg['fields'].items():
                for ref,dr in zip(rows(jet),rows(radial_bg[key])):derivative.seed(ref,g.mul(jac,dr))
            radial=total_radial_functions(g,chart,derivative,loop_owner.N)
            generic=generic_functions(g,chart['modified_fields'],self.background.delta,self.background.Pstar,chart['radius'])
            physical=physical_pressure_fields(g,chart['modified_fields'],radial,self.background.delta,self.background.Pstar)
            S=reduced_pressure_S(g,physical,self.background.delta)
            Q,Qy=a.trim(physical['Q']),a.trim(physical['Q_y'])
            boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))
            at=lambda jets,value:{key:a.at(jet,coordinate,value) for key,jet in jets.items()}
            if self.initial is None:self.initial=at(boundary,chart['domain'][0])
            local={};cumulative={}
            t=g.symbol('true_modified_pressure_C3_'+name)
            for key,power in (('S0',1),('S1',2)):
                weight=g.mul(jac,*([chart['radius']]*power))
                body=a.at(a.scale(S,weight),coordinate,t)
                local[key]=a.integral(body,t,chart['domain'][0],coordinate,
                    original_signed_reduced_pressure=True,original_chart=name,
                    physical_radius_power=power,physical_Jacobian_applied_once=True,
                    independent_P0_retained_in_original_source=True,finite_endpoints_independent_of_Z=True)
                cumulative[key]=a.add(incoming[key],local[key])
            pressure=dict(m13=a.add(a.scale(cumulative['S0'],'-1/2'),a.sub(boundary['Q_plus_Q_y'],self.initial['Q_plus_Q_y'])),
                m15=a.sub(a.scale(cumulative['S1'],'1/2'),a.sub(boundary['R_Q_y'],self.initial['R_Q_y'])))
            self.charts[name]=dict(source_chart=chart,radial_background=radial_bg,radial=radial,
                generic=generic,physical=physical,S=S,finite_boundary=boundary,
                incoming_integrals=incoming,local_integrals=local,cumulative_integrals=cumulative,pressure=pressure)
            endpoint=dict(integrals=at(cumulative,chart['domain'][1]),pressure=at(pressure,chart['domain'][1]),
                boundary=at(boundary,chart['domain'][1]))
            self.endpoints[name]=endpoint;incoming=endpoint['integrals']
        self.Rm=self.endpoints['postrestore']
        self.definition=dict(full_modified_Q_and_signed_generic_inertial_pressure_C4=True,
            actual_slow_native_loop_derivatives_from_true_function_bodies=True,
            moving_integral_limits_and_implicit_inverse_chain_rule_retained=True,
            total_modified_Ey_Vy_mean_y_C5_with_fast_phase_terms=True,
            physical_U_mean_Uy_C5_and_Q_Qy_C4_Pstar_once=True,
            leading_Q_weight='1-delta',positive_order_n1_own_weight_unchanged='1+delta',
            reduced_original_S_and_signed_m13_m15_C3_through_Rm=True,
            complete_ten_finite_intervals_and_global_original_rminus_boundary_retained=True,
            physical_measures='S0:R*Jacobian*dcoordinate; S1:R^2*Jacobian*dcoordinate',
            m15_minus_endpoint_sign_applied_once=True,
            pressure_components_not_replacing_existing_full_preRp_C1=True,
            direct_Qyy_or_N1p_C3_not_installed=True,
            remaining_Rm_patch_and_six_outer_charts_not_installed=True,
            numeric_C3_C5_mixed_bounds_and_global_frequency_not_admitted=True,
            full_n1_matching_heat_stress_recursion_and_corrected_NS_open=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.loop,self.background,self.g);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or pins!=self.hashes or raw['input_hashes']!=self.hashes
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['modified_Q4_pressure_C3_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked true modified Q4/pressure C3 source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed modified pressure source '+name)
            self.acceptance_loaded=True
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed source '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        charts={name:{key:value for key,value in chart.items() if key!='source_chart'} for name,chart in self.charts.items()}
        return encode(dict(charts=charts,endpoints=self.endpoints,initial=self.initial,Rm=self.Rm,
            zero_pressure_integrals=self.zero_pressure_integrals,definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.loop.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.hashes!=self._hashes or self._owners!=(self.loop,self.background,self.g)):
            raise ValueError('True modified radial/pressure source changed')
        previous=self.zero_pressure_integrals
        for name in CHARTS:
            if self.charts[name]['source_chart'] is not self.loop.charts[name] or self.charts[name]['incoming_integrals'] is not previous:
                raise ValueError('Original pressure source or incoming signed integrals changed')
            previous=self.endpoints[name]['integrals']


def run(loop_owner):
    began=time.monotonic();owner=CurrentOriginalModifiedQ4PressureC3Functions(loop_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,modified_Q4_pressure_C3_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_MODIFIED_Q4_TOTAL_RADIAL_AND_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':raise SystemExit('run(loop_owner) reuses checked true C5 loop source')
