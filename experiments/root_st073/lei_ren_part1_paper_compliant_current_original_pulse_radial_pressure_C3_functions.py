"""Actual current Rp source, pulse second radial recovery and pressure FTC.

Ordinary Z rows are not Taylor coefficients. Main/gap physical-y derivative
is mu*d/dxi; end derivative is d/ds. Numerical oracles remain separate.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from types import FunctionType

import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins as previous
import lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions as kernels

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_pulse_radial_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_pulse_radial_pressure_C3_functions_check.json'
GATES=('current_original_pulse_Uyy_C5_Qyy_C4_source_functions_installed',
    'current_original_pulse_radial_N1p_pressure_C3_source_functions_installed',
    'current_original_n1_partial_pulse_pressure_physical_y_FTC_C3_functions_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+(
    'current_original_n1_Rp_Rv_total_radial_pressure_C3_source_joins_installed',
    'current_original_pulse_radial_pressure_numeric_bounds_installed')))


def leading_Q(g,U,mean,delta):
    """Leading 1-delta recovery; mean C5 supplies Q C4."""
    a=previous.Algebra(g,4);z=previous.FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero,g.zero))
    one=a.fixed(1);de=a.fixed(delta);d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    numerator=a.sub(a.scale(a.mul(z,a.trim(U)),2),
        a.add(a.mul(a.sub(one,de),a.mul(z,a.trim(mean))),a.mul(d,a.shift(mean))))
    return a.div(numerator,L,'same original positive leading L=1-delta*Z^2')


def pressure_density(g,fields,delta,R):
    a=previous.Algebra(g,3);z=previous.FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero))
    one=a.fixed(1);de=a.fixed(delta);d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    U,Q,Qy=[a.trim(fields[k]) for k in ('U','Q','Q_y')];Qz=a.shift(fields['Q'])
    div=lambda q:a.div(q,L,'same original positive leading L')
    linear=div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy))
    quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))
    axial=div(a.mul(U,a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2))))
    S=a.add(linear,quadratic,axial)
    Py=a.add(a.scale(S,g.mul(g.constant('-1/2'),R)),a.trim(fields['Q_y']),a.trim(fields['Q_yy']))
    N1p=a.scale(Py,g.quotient(g.one,R,'positive actual pulse radius'))
    return dict(S=S,N1p=N1p,P1_y=Py,m15_y=a.scale(Py,g.neg(R)))


def aligned_pressure_call(original_map,*args):
    """Keep a pinned formula while aligning restored ordinary row classes."""
    aligned_globals=dict(original_map.__globals__,FunctionJet=previous.FunctionJet,
        Algebra=previous.Algebra)
    aligned_map=FunctionType(original_map.__code__,aligned_globals,original_map.__name__,
        original_map.__defaults__,original_map.__closure__)
    if (original_map.__closure__ is not None or aligned_map.__code__ is not original_map.__code__
        or aligned_globals['FunctionJet'] is not previous.FunctionJet
        or aligned_globals['Algebra'] is not previous.Algebra):
        raise ValueError('Unchanged pressure body with the actual restored source algebra required')
    return aligned_map(*args)


def pressure_FTC(g,fields,delta,R):
    """Physical-y Leibniz/FTC inputs to the actual signed pressure map."""
    a=previous.Algebra(g,3);U,Uy,Q,Qy,Qyy=[a.trim(fields[k]) for k in ('U','U_y','Q','Q_y','Q_yy')]
    Uz,Qz=a.shift(fields['U']),a.shift(fields['Q']);UU=a.mul(U,U);QQ=a.mul(Q,Q)
    products=dict(I0_q_Z=Qz,I0_U_q_Z=a.mul(U,Qz),I0_U_squared=UU,I0_U_U_Z=a.mul(U,Uz),
        I1_q=Q,I1_q_Z=Qz,I1_q_squared=QQ,I1_U_q_Z=a.mul(U,Qz),I1_U_squared=UU,I1_U_U_Z=a.mul(U,Uz))
    integral_y={k:a.scale(q,R if k.startswith('I0_') else g.mul(R,R)) for k,q in products.items()}
    qsquare_y=a.scale(a.mul(Q,Qy),2);usquare_y=a.scale(a.mul(U,Uy),2)
    upper_y=dict(R_q=a.scale(a.add(Q,Qy),R),R_squared_q=a.scale(a.add(a.scale(Q,2),Qy),g.mul(R,R)),
        R_q_squared=a.scale(a.add(QQ,qsquare_y),R),
        R_squared_q_squared=a.scale(a.add(a.scale(QQ,2),qsquare_y),g.mul(R,R)),
        R_U_squared=a.scale(a.add(UU,usquare_y),R),
        R_squared_U_squared=a.scale(a.add(a.scale(UU,2),usquare_y),g.mul(R,R)),
        q_plus_q_y=a.add(Qy,Qyy),R_q_y=a.scale(a.add(Qy,Qyy),R))
    zero={k:a.fixed(0) for k in upper_y}
    # Only the linear pressure map is differentiated. Nonlinear products
    # have their actual endpoint derivatives above, and complete preceding
    # charts/global lower tensors have no current physical query.
    # Restore can replace the core FunctionJet class while a previously
    # imported pressure module still holds its old class alias. Execute the
    # unchanged, pinned pressure body with this owner's actual jet algebra;
    # do not mutate the accepted module or convert derivatives numerically.
    signed=aligned_pressure_call(previous.intermediate.pressure_functions,g,delta,integral_y,zero,upper_y)
    return dict(actual_integral_physical_y_derivatives=integral_y,
        actual_upper_tensor_physical_y_derivatives=upper_y,actual_lower_tensor_physical_y_derivatives=zero,
        signed_pressure_physical_y_derivatives=signed,actual_products=products)


class CurrentOriginalPulseRadialPressureC3Functions:
    def __init__(self,outer=None,require_checked=True):
        self.outer=outer or previous.CurrentOriginalN1OuterOwnSourceJoins(require_checked=False)
        self.outer.assert_graph();self.inlet=self.outer.inlet;self.g=g=self.outer.g
        self.family=copy.deepcopy(self.outer.family);self.hashes=dict(self.outer.hashes)
        rec=parent.linear.read(previous.CONSTRUCTION_RECEIPT);raw=parent.linear.read(previous.NAME)
        if (rec['source_family']!=self.family or not rec['all_three_ordinary_axial_trace_identities_passed']
            or rec['ordinary_own_field_and_cumulative_trace_rows']!=204
            or raw['outer_own_source_join_bindings']!=self.outer.snapshot
            or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual three ordinary outer source traces required')
        for name,digest in rec['input_hashes'].items():
            if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed outer trace source '+name)
            self.hashes[name]=digest
        self.hashes[previous.CONSTRUCTION_RECEIPT]=sha(previous.CONSTRUCTION_RECEIPT)
        kernel=parent.linear.read(kernels.NAME);kr=parent.linear.read(kernels.RECEIPT)
        if (not kr['all_passed'] or not kr[kernels.GATE] or kr['source_family']!=self.family
            or kernel['source_family']!=self.family
            or kernel['exact_graph_nodes']!=g.nodes[:len(kernel['exact_graph_nodes'])]):
            raise ValueError('Same checked original pulse kernels required')
        for name,digest in kr['input_hashes'].items():
            if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed actual pulse kernel '+name)
            self.hashes[name]=digest
        for name in (kernels.NAME,kernels.RECEIPT,Path(kernels.__file__).name):self.hashes[name]=sha(name)
        self.W={k:self.ref(v) for k,v in kernel['exact_pulse_kernel_functions'].items()
            if k in ('original_gp','original_gp_variable','ell','normalization','log_end_scale','mu')}
        self.selected=self.inlet.selected_pulse;self.mu=self.inlet.geometry['mu'];self.delta=self.inlet.geometry['delta']
        if self.W['mu']!=self.mu:raise ValueError('Same exact original pulse mu required')
        self.fixed_geometry_source_independence=[]
        for key,ref in self.W.items():
            if key in ('original_gp','original_gp_variable'):continue
            closure=self.outer.operand_closure(ref.node)
            if g.symbol('Z').node in closure:raise ValueError('Pulse scalar geometry depends on Z '+key)
            self.fixed_geometry_source_independence.append(dict(key=key,node=ref.node,closed_source_nodes=len(closure)))
        legacy=previous.previous.pulse5.true.own
        self.sigma_source=parent.active.current.ast_binding(legacy.c2.flat_source.sigma_jets)
        self.sigma_second_source=parent.active.current.ast_binding(legacy.c2.flat_source._sigma_left)
        self.pressure_map_source=parent.active.current.ast_binding(previous.intermediate.pressure_functions)
        self.pressure_endpoint_source=parent.active.current.ast_binding(previous.intermediate.endpoint_functions)
        if self.hashes[self.pressure_map_source['module']]!=self.pressure_map_source['module_sha256']:
            raise ValueError('Pinned original signed pressure map required')
        for proof in (self.sigma_source,self.sigma_second_source):
            if self.hashes[proof['module']]!=proof['module_sha256']:raise ValueError('Pinned actual sigma derivative source required')
        self.primitive=self.ref(next(j for j in g.nodes[self.W['original_gp'].node]['arguments']
            if g.nodes[j]['operation']=='definite_integral'))
        self.charts={};a=previous.Algebra(g,5);bp=g.add(g.constant('1/2'),self.mu)
        for name,original in self.inlet.leading_fields.items():
            coordinate=self.inlet.pulse_coordinates[name];x=coordinate['query']
            jac=g.one if name=='end' else self.mu
            for key,jet in self.selected.items():
                for ref in previous.rows(jet):
                    if previous_check_free(self.outer,ref.node,x.node):raise ValueError('Selected coefficient depends on radial query')
            B=a.fixed(0);By=a.fixed(0);Byy=a.fixed(0)
            if name=='main':
                gp=self.at_scalar(self.W['original_gp'],x)
                I=self.at_scalar(self.primitive,x)
                left=g.mul(g.constant(50),x);right=g.sub(g.constant(11),x)
                sl,sl1,_=self.sigma_rows(left);sr,sr1,sr2=self.sigma_rows(right)
                first=g.sub(g.mul(sl,sr),g.mul(I,sr1))
                second=g.add(g.mul(g.constant(50),sl1,sr),g.mul(g.constant(-2),sl,sr1),g.mul(I,sr2))
                B=a.scale(self.selected['ap'],gp);By=a.scale(self.selected['ap'],g.mul(self.mu,first))
                Byy=a.scale(self.selected['ap'],g.mul(self.mu,self.mu,second))
            elif name=='end':
                ell=self.W['ell'];normal=self.W['normalization']
                for key,center in (('C1',-3),('C2',-1)):
                    local=g.sub(x,g.constant(center));argument=g.quotient(local,ell,'positive original ell')
                    raw_beta=g.node('compact_raw_beta',argument=argument.node,
                        definition='exp(-1/(1-r^2)) on abs(r)<1; zero otherwise')
                    first,second=previous.previous.compact.leading.beta_first_two(g,argument,raw_beta)
                    B=a.add(B,a.scale(self.selected[key],g.quotient(raw_beta,g.mul(ell,normal),'positive original ell and beta normalization')))
                    By=a.add(By,a.scale(self.selected[key],g.quotient(first,g.mul(ell,ell,normal),'positive ell squared and beta normalization')))
                    Byy=a.add(Byy,a.scale(self.selected[key],g.quotient(second,g.mul(ell,ell,ell,normal),'positive ell cubed and beta normalization')))
                factor=g.unary('exp',self.W['log_end_scale']);B,By,Byy=[a.scale(q,factor) for q in (B,By,Byy)]
            elif name!='gap':raise ValueError('Original main/gap/end chart required')
            theta=original['Utheta'];Uyy=a.mul(theta,a.add(Byy,a.scale(By,g.mul(g.constant(-2),bp)),a.scale(B,g.mul(bp,bp))))
            myy=a.sub(original['U_y'],original['m_y']);Qyy=leading_Q(g,Uyy,myy,self.delta)
            fields=dict(original,U_yy=Uyy,m_yy=myy,Q_yy=Qyy)
            density=pressure_density(g,fields,self.delta,coordinate['radius'])
            FTC=pressure_FTC(g,fields,self.delta,coordinate['radius'])
            self.charts[name]=dict(coordinate=x,radius=coordinate['radius'],physical_y_Jacobian=jac,
                domain=self.outer.charts['pulse_'+name]['domain'],fields=fields,
                normalized_axial=B,normalized_axial_y=By,normalized_axial_yy=Byy,
                density=density,FTC=FTC,original_gp_integral_retained=True)
        # Rebuild actual partial primitives on the current Rp source,
        # with the field query name as the integral binder. Retain every
        # complete preceding chart and the finite global lower tensors.
        a3=previous.Algebra(g,3);self.complete_integrals={}
        for name,data in self.charts.items():
            x=data['coordinate'];lo,hi=data['domain'];R=data['radius'];jac=data['physical_y_Jacobian']
            dx=g.quotient(R,jac,'positive actual physical-y chart Jacobian')
            bodies={k:a3.scale(v,dx if k.startswith('I0_') else g.mul(R,dx))
                for k,v in data['FTC']['actual_products'].items()}
            def integrate(body,upper,key):
                return previous.FunctionJet(tuple(previous.intermediate.controls.integral(g,ref,
                    g.nodes[x.node]['name'],lo,upper,endpoint_independent_of_Z=True,
                    ordinary_slow_Z_derivative_order=n,physical_radius_power=1 if key.startswith('I0_') else 2,
                    physical_Jacobian_applied_once=True,actual_current_Rp_pulse_body=True,
                    complete_original_history_no_cutoff=True) for n,ref in enumerate(previous.rows(body))))
            data['actual_integrands']=bodies
            data['partial_integrals']={k:integrate(body,x,k) for k,body in bodies.items()}
            self.complete_integrals[name]={k:integrate(body,hi,k) for k,body in bodies.items()}
        first=self.charts['main'];x=first['coordinate'];lo=first['domain'][0]
        initial={k:a3.at(a3.trim(first['fields'][k]),x,lo) for k in ('U','Q','Q_y')}
        r0=a3.at(a3.fixed(first['radius']),x,lo).value
        self.global_lower_tensors=aligned_pressure_call(previous.intermediate.endpoint_functions,g,
            initial['U'],initial['Q'],initial['Q_y'],r0)
        names=('main','gap','end')
        for name,data in self.charts.items():
            data['preceding_complete_charts']=list(names[:names.index(name)])
            data['total_integrals']={k:a3.add(*(self.complete_integrals[earlier][k]
                for earlier in data['preceding_complete_charts']),v) for k,v in data['partial_integrals'].items()}
            f=data['fields'];data['upper_tensors']=aligned_pressure_call(previous.intermediate.endpoint_functions,
                g,f['U'],f['Q'],f['Q_y'],data['radius'])
            data['signed_pressure_increment']=aligned_pressure_call(previous.intermediate.pressure_functions,
                g,self.delta,data['total_integrals'],self.global_lower_tensors,data['upper_tensors'])
        self.endpoints={name:{side:self.at(name,end) for side,end in zip(('left','right'),data['domain'])}
            for name,data in self.charts.items()}
        self.definition=dict(actual_current_Rp_C5_source_inputs_consumed=True,
            original_finite_mu_gp_integral_and_beta_normalization_retained=True,
            sigma_beta_second_derivative_source_not_numerical_caps=True,
            field_Z_rows_already_ordinary_not_multiplied_by_factorial_again=True,
            sigma_second_is_two_times_source_Taylor_row2=True,
            leading_one_minus_delta_not_own_one_plus_delta=True,
            physical_y_Jacobian_main_gap_mu_end_one=True,
            actual_partial_pressure_Leibniz_and_FTC_inputs_constructed=True,
            actual_partial_primitive_binder_is_the_field_query_name=True,
            preceding_complete_pulse_charts_and_finite_global_lower_tensors_retained=True,
            outer_radial_pressure_joins_higher_mixed_numeric_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.outer,self.inlet,g);self._family=copy.deepcopy(self.family)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def sigma_rows(self,x):
        g=self.g;c=g.constant
        sigma=previous.previous.pulse5.true.own.jets.complete.Gamma.sigma(g,x)
        first=g.node('analytic_unary',name='original_flat_sigma_prime',argument=x.node,
            definition='sigma(x)*(1-sigma(x))*(2/x^3+2/(1-x)^3) for 0<x<1; zero otherwise',
            derivative_source=self.sigma_source,ordinary_derivative='same original sigma Taylor coefficient[1]',
            exact_flat_source_derivative_not_a_bound=True)
        second=g.node('analytic_unary',name='original_flat_sigma_second',argument=x.node,
            definition='sigma*(1-sigma)*(6/(1-x)^4-6/x^4+(1-2*sigma)*(2/x^3+2/(1-x)^3)^2) inside (0,1); zero outside',
            derivative_source=self.sigma_second_source,ordinary_derivative='2*same original sigma Taylor coefficient[2]',
            exact_flat_source_derivative_not_a_bound=True)
        return sigma,first,second
    def at_scalar(self,body,value):
        return self.g.node('function_substitution',expression=body.node,variable=self.W['original_gp_variable'].node,
            value=value.node,Z_independent_substitution=True)
    def at(self,name,value):
        data=self.charts[name];x=data['coordinate']
        group=lambda rows:{k:previous.Algebra(self.g,len(previous.rows(q))-1).at(q,x,value) for k,q in rows.items()}
        return dict(fields=group(data['fields']),density=group(data['density']),
            pressure_FTC=group(data['FTC']['signed_pressure_physical_y_derivatives']),
            signed_pressure_increment=group(data['signed_pressure_increment']),
            radius=previous.Algebra(self.g,3).at(previous.Algebra(self.g,3).fixed(data['radius']),x,value).value)
    def _snapshot(self):
        return previous.encode(dict(W=self.W,selected=self.selected,mu=self.mu,delta=self.delta,
            fixed_geometry_source_independence=self.fixed_geometry_source_independence,
            sigma_source=self.sigma_source,sigma_second_source=self.sigma_second_source,
            pressure_map_source=self.pressure_map_source,
            pressure_endpoint_source=self.pressure_endpoint_source,
            complete_integrals=self.complete_integrals,global_lower_tensors=self.global_lower_tensors,
            primitive=self.primitive,charts=self.charts,endpoints=self.endpoints,definition=self.definition))
    def assert_graph(self):
        self.outer.assert_graph()
        if (self.hashes!=self._hashes or self.family!=self._family or self._owners!=(self.outer,self.inlet,self.g)
            or self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.numerical_oracle_available):
            raise ValueError('Actual pulse radial source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed pulse radial source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT)
        pins=dict(rec['input_hashes']);archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py'
        if (not rec['all_passed'] or not all(rec[k] for k in GATES) or any(rec[k] for k in OPEN)
            or pins.pop(checker)!=sha(checker) or archive!=sha(NAME) or pins!=self.hashes
            or raw['input_hashes']!=self.hashes or raw['source_family']!=self.family or rec['source_family']!=self.family
            or raw['pulse_radial_pressure_C3_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual pulse radial/FTC source required')
        self.acceptance_loaded=True
    def at_fraction(self,name,fraction):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Checked pulse radial source required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Exact chart fraction must lie in [0,1]')
        lo,hi=self.charts[name]['domain'];value=self.g.add(lo,self.g.mul(self.g.constant(q),self.g.sub(hi,lo)))
        return dict(chart=name,exact_fraction=str(q),**self.at(name,value),numerical_oracle_available=False,global_n1=False)


def previous_check_free(owner,node,query):
    # Import the already pinned binder-aware source dependency algorithm.
    import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_check as check
    return check.free_query_dependency(owner,node,query)


def run(outer=None):
    began=time.monotonic();owner=CurrentOriginalPulseRadialPressureC3Functions(outer,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,pulse_radial_pressure_C3_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_PULSE_SECOND_RADIAL_AND_PRESSURE_FTC_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
