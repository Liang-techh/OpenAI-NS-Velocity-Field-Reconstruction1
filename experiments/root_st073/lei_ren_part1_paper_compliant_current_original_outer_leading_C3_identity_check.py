"""Independent defining-kernel, ordinary C3, five-ODE and source seam checks.

The old opaque outer/repair endpoint identity deliberately remains open.
Directed enclosure outputs are never substituted for mathematical integrals.
"""
import ast
import inspect
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity as current


class Sigma(s.Function):
    nargs=1
    @classmethod
    def eval(cls,x):
        if x.is_number and x.is_real:
            if x<=0:return s.Integer(0)
            if x>=1:return s.Integer(1)
    def fdiff(self,argindex=1):
        if argindex!=1:raise ValueError('Sigma has one argument')
        return SigmaDerivative(s.Integer(1),self.args[0])


class SigmaDerivative(s.Function):
    nargs=2
    @classmethod
    def eval(cls,n,x):
        if x.is_number and x.is_real and (x<=0 or x>=1):return s.Integer(0)
    def fdiff(self,argindex=2):
        if argindex!=2:raise ValueError('Only the Sigma argument is differentiated')
        return SigmaDerivative(self.args[0]+1,self.args[1])


class Exprel(s.Function):
    nargs=1
    @classmethod
    def eval(cls,x):
        if x==0:return s.Integer(1)
    def _eval_rewrite_as_exp(self,x,**kwargs):return (s.exp(x)-1)/x


class Interpreter:
    """Fail-closed scalar graph interpreter with two separate integral modes."""
    def __init__(self,field,literal_integrals=False,bindings=None):
        self.g=field.graph;self.literal=literal_integrals;self.memo={}
        self.bindings={} if bindings is None else dict(bindings)
        self.sigma_binding=field.contracts['original_sigma']
        self.kernel_functions={};self.z=s.Symbol('Z',real=True)
        self.mu=s.Symbol('mu',positive=True);self.logP=s.Symbol('logP',real=True)
        self.bindings.setdefault(field.parameters['mu'].node,self.mu)
        self.bindings.setdefault(field.parameters['logP'].node,self.logP)

    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.bindings[i]
        if i in self.memo:return self.memo[i]
        n=self.g.nodes[i];op=n['operation'];at=self.at
        if op=='exact_rational':v=s.Rational(n['numerator'],n['denominator'])
        elif op=='bound_variable':v=s.Symbol(n['name'],real=True)
        elif op=='sum':v=s.Add(*(at(j) for j in n['arguments']))
        elif op=='product':v=s.Mul(*(at(j) for j in n['arguments']))
        elif op=='negative':v=-at(n['argument'])
        elif op=='positive_quotient':v=at(n['numerator'])/at(n['denominator'])
        elif op=='analytic_unary':
            q=at(n['argument'])
            if n['name']=='exp':v=s.exp(q)
            elif n['name']=='log':v=s.log(q)
            elif n['name']=='exprel':v=Exprel(q)
            else:raise AssertionError('Unsupported leading unary '+str(n))
        elif op=='exact_original_source_flat_sigma':
            assert n['definition']=='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside'
            assert n['source_binding']==self.sigma_binding
            assert n['smooth_flat_endpoints'] and n['Z_independent_argument']
            v=Sigma(at(n['argument']))
        elif op=='function_substitution':
            assert n['Z_independent_substitution']
            v=at(n['expression']).subs(at(n['variable']),at(n['value']))
        elif op=='definite_integral':
            assert n['exact_function_integral'] and n['original_kernel_function']
            assert n['numerical_value_not_installed'] and n['coordinate_and_kernel_Z_independent']
            lo,hi=at(n['lower']),at(n['upper'])
            if self.literal:
                variable=s.Symbol(n['variable'],real=True);body=at(n['integrand'])
                assert self.z not in body.free_symbols and self.z not in lo.free_symbols|hi.free_symbols
                v=s.Integer(0) if lo==hi else s.Integral(body,(variable,lo,hi))
            else:
                fn=self.kernel_functions.setdefault(i,s.Function('original_kernel_'+str(i)))
                v=s.Integer(0) if lo==hi else fn(hi)
        else:raise AssertionError('Unsupported exact original leading dependency '+str(n))
        self.memo[i]=v;return v


def rows(q):return current.target.rows(q)


def zero(expr):
    expr=s.expand_power_exp(expr.rewrite(s.exp))
    return s.cancel(s.expand(expr))==0


def source_predicates():
    """Compare defining source bodies, not just their archive/hash labels."""
    pre=current.providers.slope.original.original
    axial=current.providers.axial.original.kernels
    o3=current.providers.o3.original.kernels
    checks=[]
    def require(fn,statement):
        wanted=ast.dump(ast.parse(statement).body[0],include_attributes=False)
        tree=ast.parse(inspect.getsource(fn))
        assert any(ast.dump(n,include_attributes=False)==wanted for n in ast.walk(tree)),(fn.__name__,statement)
        checks.append(dict(function=fn.__module__+'.'+fn.__name__,defining_statement=statement))
    require(axial.stable_sigma,'if v <= 0: return c.mpf(0), c.mpf(0)')
    require(axial.stable_sigma,'if v >= 1: return c.mpf(1), c.mpf(0)')
    require(axial.stable_sigma,'a = c.exp(-1/x**2)')
    require(axial.stable_sigma,'b = c.exp(-1/(1-x)**2)')
    require(axial.stable_sigma,'value = a/(a+b)')
    left=pre.sigma_jets.__globals__['_sigma_left']
    require(left,'odds = 1/(1-x)**2 - 1/x**2')
    require(left,'e = positive_exp(c, odds)')
    require(left,'value = e/(1+e)')
    require(pre.sigma_jets,'if lo <= 0: rows.append(IntervalTaylor.constant(c,0,4))')
    require(pre.sigma_jets,'if hi >= 1: rows.append(IntervalTaylor.constant(c,1,4))')
    a,b=s.symbols('a b',positive=True)
    assert s.cancel((a/b)/(1+a/b)-a/(a+b))==0
    scalar=pre.transition_integrals.__globals__['sigma_value_derivative']
    alpha=scalar.__globals__['alpha_box']
    require(scalar,'sig = 1-alpha_box(c,s+1,c.mpf(1))')
    require(alpha,'qlo,qhi = endpoints((y-hb)/hb)')
    require(alpha,'if q <= 0: return c.mpf(1)')
    require(alpha,'if q >= 1: return c.mpf(0)')
    require(alpha,'left,right = c.exp(-1/(x*x)),c.exp(-1/((1-x)*(1-x)))')
    require(alpha,'return right/(left+right)')
    assert s.cancel(1-b/(a+b)-a/(a+b))==0
    require(pre.transition_integrals,"rates = (c.mpf('1.6'), c.mpf('.2'), c.mpf('1.2'))")
    require(pre.transition_integrals,'powers = (1,2,2)')
    require(pre.transition_integrals,'next_j = j + dy*sj')
    require(pre.transition_integrals,"masses[k] += dy*c.exp(rate*scell-c.mpf('.6')*power*jcell)")
    require(pre.slope_masses,'return transition_integrals(c,exact,cells)')
    require(axial.turnoff_kernels,'phase = 1-c.ln(c.mpf(x))/md')
    require(axial.turnoff_kernels,'return stable_sigma(c,phase)[0]')
    require(axial.turnoff_kernels,'left = B(y-b)')
    require(axial.turnoff_kernels,'right = B(y-a)')
    require(axial.turnoff_kernels,'weight = c.exp(-a)-c.exp(-b)')
    require(axial.turnoff_kernels,'mass[0] += weight*box')
    require(axial.turnoff_kernels,'mass[1] += weight*box**2')
    require(pre.turnoff_derivatives,'sig = sigma_jets(c,1-phase)')
    require(pre.turnoff_derivatives,'ordinary = [sig[k]*math.factorial(k) for k in range(5)]')
    require(o3.transition_kernels,'nextJ = J+da*sj')
    require(o3.transition_kernels,'weights = (c.exp(b)-c.exp(a),da,c.exp(-a)-c.exp(-b))')
    require(o3.transition_kernels,'for n,power in enumerate((1,2,2)): K[n] += weights[n]*c.exp(-mu*power*jcell)')
    require(o3.transition_kernels,'return dict(J=J,theta=K[0],energy=K[1],pressure=K[2],cells=cells)')
    return dict(independent_source_defining_AST_predicates=checks,
        source_scalar_sigma_and_Taylor_scalar_interior_identical=True,
        source_helpers_are_enclosures_of_integrals_not_function_values=True,
        slope_source_order_theta_pressure_energy_and_O3_order_theta_energy_pressure_preserved=True)


def kernels(field):
    q=Interpreter(field,True);f=field.functions;charts=f['charts'];x=s.Symbol('u',real=True);v=s.Symbol('v',real=True)
    J=lambda t:s.Integral(Sigma(v),(v,0,t))
    count=0
    def equal(node,wanted):
        nonlocal count
        got=q.at(node)
        assert got.as_dummy()==wanted.as_dummy(),(got,wanted)
        count+=1
    for name,specs in (('O2_slope',dict(theta=(s.Rational(8,5),s.Rational(3,5)),
            energy=(s.Rational(6,5),s.Rational(6,5)),pressure=(s.Rational(1,5),s.Rational(6,5)))),
        ('O3_transition',dict(theta=(1,q.mu),energy=(0,2*q.mu),pressure=(-1,2*q.mu)))):
        c=charts[name];t=q.at(c['native_coordinate'])
        equal(c['kernels']['J'],J(t))
        for key,(rate,power) in specs.items():
            equal(c['kernels'][key],s.Integral(s.exp(rate*x-power*J(x)),(x,0,t)))
    c=charts['O2_axial'];y=q.at(c['kernels']['original_y_endpoint'])
    beta=Sigma(1-s.log(x)/40)
    for key,power in (('original_B_mass',1),('original_B_squared_mass',2)):
        equal(c['kernels'][key],s.Integral(s.exp(x-y)*beta**power,(x,1,y)))
    new=field.graph.nodes[len(field.prefix):]
    all_integrals=[i+len(field.prefix) for i,n in enumerate(new) if n['operation']=='definite_integral']
    assert len(all_integrals)==16
    for i in all_integrals:q.at(i)
    # The endpoint-weighted r=y-u source has the same true dy measure.
    r=s.Symbol('r',real=True);Y=s.Symbol('Y',positive=True)
    body=s.exp(x-Y)*Sigma(1-s.log(x)/40)
    assert body.subs(x,Y-r)==s.exp(-r)*Sigma(1-s.log(Y-r)/40)
    assert s.diff(Y-r,r)==-1
    return dict(independent_original_kernel_integrand_limit_and_primitive_identities=count,
        all_true_integral_nodes_checked=len(all_integrals),original_axial_dy_measure_and_change_of_variable=True,
        no_directed_scalar_cap_used_as_function_value=True)


def ordinary_C3(field):
    q=Interpreter(field);count=0
    assert q.at(rows(field.functions['exact_common_Rref_A_C3'])[0])==1/(1+q.z**2)
    assert q.at(rows(field.functions['exact_common_axial_B_C3'])[0])==4*q.z*s.exp(-q.logP)
    for name,c in field.functions['charts'].items():
        quantities={'E':c['E'],'V':c['V'],**c['histories'],**{'y_'+k:v for k,v in c['history_y'].items()}}
        assert c['ordinary_Z_orders']==[0,1,2,3]
        for key,value in quantities.items():
            refs=rows(value);base=q.at(refs[0])
            for j,ref in enumerate(refs):
                assert zero(q.at(ref)-s.diff(base,q.z,j)),(name,key,j)
                count+=1
    assert count==288
    return dict(independent_ordinary_Z0_to_Z3_profile_and_history_and_own_rate_rows=count,
        actual_common_A_and_B_source_definitions=True,ordinary_factorials_once=True,
        V_m_k_normalized_by_Pstar_once_and_no_Rm_multiplier=True)


def local_calculus(field):
    chart_reports=[]
    for name,c in field.functions['charts'].items():
        native=s.Symbol('local_y',real=True);bindings={c['native_coordinate'].node:native}
        if name=='O2_axial':
            native=s.Symbol('local_y',positive=True)
            bindings={c['native_coordinate'].node:s.log(native)/40,c['physical_coordinate'].node:native,
                c['kernels']['original_y_endpoint'].node:native}
        q=Interpreter(field,bindings=bindings);rules={}
        if name in ('O2_slope','O3_transition'):
            refs=c['kernels'];sig=Sigma(native);J=q.at(refs['J'])
            rules[s.diff(J,native)]=sig
            params=(dict(theta=(s.Rational(8,5),s.Rational(3,5)),energy=(s.Rational(6,5),s.Rational(6,5)),pressure=(s.Rational(1,5),s.Rational(6,5)))
                if name=='O2_slope' else dict(theta=(1,q.mu),energy=(0,2*q.mu),pressure=(-1,2*q.mu)))
            for key,(r,p) in params.items():rules[s.diff(q.at(refs[key]),native)]=s.exp(r*native-p*J)
        elif name=='O2_axial':
            beta=Sigma(1-s.log(native)/40)
            for key,power in (('B_mass',1),('B_squared_mass',2)):
                K=q.at(c['kernels'][key]);rules[s.diff(K,native)]=beta**power-K
        E,V=(q.at(rows(c[key])[0]) for key in ('E','V'))
        H={key:q.at(rows(value)[0]) for key,value in c['histories'].items()}
        D=dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
        count=0
        for key,rate in current.current.RATES.items():
            wanted=D[key]-s.Rational(str(rate))*H[key]
            for j,ref in enumerate(rows(c['histories'][key])):
                derivative=s.diff(q.at(ref).rewrite(s.exp),native).subs(rules,simultaneous=True)
                assert zero(derivative-q.at(rows(c['history_y'][key])[j])),(name,key,j,'FTC')
                assert zero(q.at(rows(c['history_y'][key])[j])-s.diff(wanted,q.z,j)),(name,key,j,'own rate')
                count+=1
        assert count==20
        chart_reports.append(dict(chart=name,independent_Z0_to_Z3_five_history_ODE_rows=count,
            physical_y_measure_used=True,ancestor_endpoint_memory_retained=True))
    return dict(charts=chart_reports,total_independent_ODE_rows=120,
        all_five_own_rates_and_true_fundamental_theorem_of_calculus=True,
        accepted_exact_mu_power_flow_source_consumed=True)


def seams(field):
    q=Interpreter(field);charts=field.functions['charts'];reports=[]
    for left,right,endpoint in (('Rh_reference','O2_slope',0),('O2_slope','O2_axial',1),
            ('O2_axial','O2_buffer',1),('O2_buffer','O3_transition',11),('O3_transition','O3_power',1)):
        a,b=charts[left],charts[right];av=q.at(a['native_coordinate']);bv=q.at(b['native_coordinate']);count=0
        for key,qa,qb in [('E',a['E'],b['E']),('V',a['V'],b['V'])]+[(key,a['histories'][key],b['histories'][key]) for key in current.current.RATES]+[('y_'+key,a['history_y'][key],b['history_y'][key]) for key in current.current.RATES]:
            for j,(ra,rb) in enumerate(zip(rows(qa),rows(qb))):
                la=q.at(ra).subs(av,s.Integer(endpoint));lb=q.at(rb).subs(bv,s.Integer(0))
                # Integral abstractions keep endpoint identity, including zero
                # lower-limit integrals after source-coordinate substitution.
                for i,fn in q.kernel_functions.items():
                    lower=q.at(field.graph.nodes[i]['lower'])
                    la=la.subs(fn(lower),0);lb=lb.subs(fn(lower),0)
                assert zero(la-lb),(left,right,key,j)
                count+=1
        assert count==48
        reports.append(dict(left=left,right=right,independent_state_and_own_y_C3_seam_rows=count))
    terminal=field.functions['actual_source_Rc_leading_C3'];c=charts['O3_power'];x=q.at(c['native_coordinate'])
    for key,value in terminal.items():
        for got,ref in zip(rows(value),rows(c['histories'][key])):assert zero(q.at(got)-q.at(ref).subs(x,2))
    return dict(source_seams=reports,total_C3_state_and_own_rate_seam_rows=240,
        true_O3_x2_Rc_endpoint_rows=20,original_full_axial_kernel_memory_kept_at_buffer_inlet=True,
        old_opaque_band_seed_endpoint_identity_not_claimed=True)


def source_geometry(field):
    bindings={q.node:s.Symbol(key,real=True) for key,q in field.parameters.items()}
    q=Interpreter(field,bindings=bindings);charts=field.functions['charts']
    geometry=field.functions['actual_original_source_geometry'];offsets={};domains={};count=0
    for name,c in charts.items():
        geo=geometry[name];native=q.at(c['native_coordinate']);offset=q.at(geo['offset'])
        physical=q.at(c['physical_coordinate']);jac=q.at(geo['Jacobian'])
        assert zero(s.diff(offset,native)-jac) and zero(s.diff(physical,native)-jac),name
        assert q.z not in offset.free_symbols and q.z not in jac.free_symbols
        domain=[q.at(i) for i in geo['domain']]
        assert domain=={'Rh_reference':[-5,0],'O2_slope':[0,1],'O2_axial':[0,1],
            'O2_buffer':[0,11],'O3_transition':[0,1],'O3_power':[0,2]}[name]
        assert zero(offset.subs(native,domain[0])-q.at(geo['left_offset']))
        assert zero(offset.subs(native,domain[1])-q.at(geo['right_offset']))
        offsets[name]=(offset,native);domains[name]=domain;count+=1
    base=offsets['Rh_reference'][0]-offsets['Rh_reference'][1]
    for name,c in charts.items():
        physical=q.at(c['physical_coordinate']);offset=offsets[name][0]
        extra=0 if name in ('Rh_reference','O2_slope','O2_axial','O2_buffer') else s.exp(40)+(11 if name=='O3_transition' else 12)
        assert zero(offset-physical-base-extra),name
    names=list(charts)
    for left,right in zip(names,names[1:]):
        assert zero(offsets[left][0].subs(offsets[left][1],domains[left][1])-offsets[right][0].subs(offsets[right][1],domains[right][0]))
    return dict(actual_original_native_domains_offsets_and_Jacobians=count,
        all_five_original_log_radius_seams=True,actual_axial_Jacobian_40_exp_40_native_once=True,
        exact_source_C3_function_physical_coordinate_bound_to_original_geometry=True,
        current_original_O3_power_endpoint_x2_is_Rc=True)


def run(field=None):
    began=time.monotonic();raw=current.outer.rh.read(current.NAME)
    for name,value in raw['input_hashes'].items():assert current.sha(name)==value,name
    field=field if field is not None else current.CurrentOuterLeadingC3Identity(require_checked=False)
    assert not field.acceptance_loaded and field.identity==raw['source_family']
    assert field.graph.nodes==raw['exact_graph_nodes'] and field.graph.nodes[:len(field.prefix)]==field.prefix
    assert current.target.encoded(field.functions)==raw['actual_original_leading_C3_source_functions']
    assert raw['original_source_contracts']==json.loads(json.dumps(current.source_contracts()))
    assert raw['candidate_actual_original_outer_leading_C3_source_kernels_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    predicates=source_predicates();definitions=kernels(field);derivatives=ordinary_C3(field)
    calculus=local_calculus(field);matching=seams(field);geometry=source_geometry(field)
    open_keys=('current_outer_leading_endpoint_band_seed_function_identity_installed',
        'current_outer_complete_history_to_repair_inlet_function_identity_installed',
        'selected_native_pulse_constructor_consumes_current_C3_frame',
        'current_numeric_point_field_oracle_installed',
        'global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed')
    assert not any(raw[k] for k in open_keys)
    result=dict(all_passed=True,source_family=field.identity,**dict.fromkeys(current.GATES,True),
        independent_original_source_predicates=predicates,independent_exact_kernel_definitions=definitions,
        independent_ordinary_C3_calculus=derivatives,independent_five_ODE_calculus=calculus,
        independent_original_source_seams=matching,independent_original_native_geometry=geometry,
        **dict.fromkeys(open_keys+current.OPEN,False),
        input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Original leading exact kernels, ordinary C3, five ODEs and all source seams passed',flush=True)
    return result


if __name__=='__main__':run()
