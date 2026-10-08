"""Safe same-original O3 microscopic kernels, source amplitude and history rows.

Copied original source/physical operators run in directed factored arithmetic.
No accepted provider is monkeypatched; inherited histories and P0 are kept.
"""
import ast
import copy
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source as previous
import lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 as original
import lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator as stress
import lei_ren_part1_paper_compliant_actual_Rp_source_join as binding

rc,common,prior,slow=previous.rc,previous.common,previous.prior,previous.slow
first,current=previous.first,previous.current
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
ep,iv,encode=previous.ep,previous.iv,previous.encode
ZERO=(0,0);N=previous.N
NAME=PREFIX+'current_transition_factored_native_source.json'
RECEIPT=PREFIX+'current_transition_factored_native_source_check.json'
GATE='original_microscopic_transition_safe_kernels_amplitude_and_native_history_rows_executed'


class HalfPstarCoordinates(rc.history.CommonSourceCoordinates):
    """Same original Pstar, with sqrt(Pstar) as its shared arithmetic unit."""
    def __init__(self,c,logP_squared,family):
        super().__init__(c,logP_squared,family)
        self.bases=(c.mpf(0),self.logP_squared,c.mpf(0),self.logP_squared/4,c.mpf(0))
    def rebase(self,value,family):
        self.require_family(family);c=self.ctx;p=value.scale.powers
        if c.mpf(value.scale.bases[1])._mpi_!=self.logP_squared._mpi_:raise ValueError('Same original Pstar required')
        known=c.mpf(value.scale.bases[3])._mpi_==self.bases[3]._mpi_
        offset=c.mpf(value.scale.offset)+sum((c.mpf(value.scale.bases[j])*p[j] for j in (0,2,4) if p[j]),c.mpf(0))
        if not known and p[3]:offset+=c.mpf(value.scale.bases[3])*p[3]
        return prior.ScaledEnclosure(prior.FormalScale(self.bases,(0,0,0,4*p[1]+(p[3] if known else 0),0),offset),c.mpf(value.coefficient),self.ledger)
    def record(self):
        return dict(source_family=self.family,common_log_bases=self.bases,
            shared_symbolic_factor='same original sqrt(Pstar); base3=base1/4 exactly',
            exact_original_Pstar_unit_relation='Pstar**2=(sqrt(Pstar))**4',
            common_basis_is_arithmetic_not_a_new_source_function=True,
            entire_original_parameter_covers_retained=True)


def records(value):
    if isinstance(value,prior.ScaledEnclosure):return value.record()
    if isinstance(value,ScaledTaylor):return [v.record() for v in value.coefficients]
    if isinstance(value,dict):return {k:records(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [records(v) for v in value]
    return value


class ScaledTaylor:
    """Ordinary axial Taylor coefficients with original formal source factors."""
    def __init__(self,template,coefficients):
        self.template=template;self.ctx=template.ctx
        self.coefficients=tuple(v if isinstance(v,prior.ScaledEnclosure) else template.scalar(v) for v in coefficients)
        self.order=len(self.coefficients)-1
    def __getitem__(self,n):return self.coefficients[n]
    def truncate(self,n):return ScaledTaylor(self.template,self.coefficients[:min(n,self.order)+1])
    def pair(self,other):
        if not isinstance(other,ScaledTaylor):
            other=ScaledTaylor(self.template,[other]+[0]*self.order)
        n=min(self.order,other.order)
        return self.truncate(n),other.truncate(n)
    def __neg__(self):return ScaledTaylor(self.template,[-v for v in self.coefficients])
    def __add__(self,other):
        a,b=self.pair(other);return ScaledTaylor(self.template,[a[k]+b[k] for k in range(a.order+1)])
    __radd__=__add__
    def __sub__(self,other):return self+-self.pair(other)[1]
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):
        a,b=self.pair(other)
        rows=[]
        for n in range(a.order+1):
            terms=[current.square(a[k]) if a[k] is b[n-k] else a[k]*b[n-k] for k in range(n+1)]
            rows.append(sum(terms,self.template.scalar(0)))
        return ScaledTaylor(self.template,rows)
    __rmul__=__mul__
    def reciprocal(self):
        v=self[0];lo,hi=ep(v.coefficient)
        if lo<=0:raise ArithmeticError('Positive original axial denominator cover required')
        lower=v.scale.evaluate()+self.ctx.ln(self.ctx.mpf(lo))
        rows=[v.scalar(1).positive_divide(v,lower)]
        for n in range(1,self.order+1):
            terms=sum((self[k]*rows[n-k] for k in range(1,n+1)),v.scalar(0))
            rows.append((-terms).positive_divide(v,lower))
        return ScaledTaylor(self.template,rows)
    def __truediv__(self,other):return self*self.pair(other)[1].reciprocal()
    def __rtruediv__(self,other):return self.reciprocal()*other
    def __pow__(self,n):
        if n<0:return self.reciprocal()**(-n)
        out=self*0+1
        for _ in range(n):out=out*self
        return out
    def derivative(self):return ScaledTaylor(self.template,[(n+1)*self[n+1] for n in range(self.order)])


class ProtocolConstant(prior.ScaledEnclosure):
    """Original scalar literals in a copied source's factored context."""
    def __truediv__(self,other):
        if type(other) in (int,float):return self*(self.ctx.mpf(1)/other)
        lo=ep(other.coefficient)[0]
        if lo<=0:raise ArithmeticError('Positive source literal denominator required')
        return self.positive_divide(other,other.scale.evaluate()+self.ctx.ln(self.ctx.mpf(lo)))
    def __pow__(self,n):
        if type(n) is not int or n<0:raise ValueError('Exact nonnegative scalar power')
        out=self.scalar(1)
        for _ in range(n):out=out*self
        return out


def bounded(value):
    if value.zero:return value.ctx.mpf(0)
    return value.coefficient*value.bounded_exp(value.scale.evaluate())


def exponential_factor(exponent,ledger):
    """Exact exp(x)=1+expm1(x); retain the nonzero correction separately."""
    c=exponent.ctx;x=bounded(exponent)
    if max(abs(v) for v in ep(x))>1:raise ArithmeticError('Original microscopic exponent cover required')
    increment=rc.density.density.factored_expm1(exponent)
    value=exponent.scalar(1)+increment
    ledger.append(dict(original_exponent=exponent.record(),original_expm1_increment=increment.record(),
        enclosing_exponential=value.record(),exact_identity='exp(x)=1+x*integral_0^1 exp(t*x)dt',
        nonzero_source_increment_not_replaced_by_zero=not exponent.zero))
    return value


def positive_hull(left,right):
    if left.zero and right.zero:return left
    c=left.ctx
    # Whole positive source range [left.lower,right.upper]; arithmetic upper
    # coordinates never select the source. Endpoint source expressions remain.
    ref=max(ep(x.scale.evaluate())[1] for x in (left,right) if not x.zero)
    coeffs=[]
    for x in (left,right):
        coeffs.append(c.mpf(0) if x.zero else x.coefficient*x.bounded_exp(x.scale.evaluate()-c.mpf(ref)))
    lo=max(0,ep(coeffs[0])[0]);hi=max(0,ep(coeffs[1])[1])
    if lo>hi:raise ArithmeticError('Original monotone positive source bracket')
    return prior.ScaledEnclosure(prior.FormalScale(left.scale.bases,offset=c.mpf(ref)),c.mpf((lo,hi)),left.ledger)


def original_sigma_rows(template,s_cover):
    """Same logistic, retaining positive log factors through ordinary y4."""
    c=template.ctx;s=c.mpf(s_cover)
    if ep(s)==(0,0):return [template.scalar(0)]*5
    if ep(s)[0]<=0:raise ArithmeticError('Positive source endpoint or exact zero required')
    got=previous.previous.positive_log_sigma(template,s)
    sigma=got['values'][ZERO];inv=iv(c,encode(got['record']['positive_denominator_inverse']))
    L1=2/(1-s)**3+2/s**3;L2=6/(1-s)**4-6/s**4
    L3=24/(1-s)**5+24/s**5;L4=120/(1-s)**6-120/s**6
    prod=sigma*inv;a=2*inv-1
    out=[got['values'][(j,0)] for j in range(3)]
    out.append(prod*((1-prod*6)*L1**3+(L3+3*a*L1*L2)))
    out.append(prod*((1-prod*6)*(L1**2*L2*6)+(1-prod*12)*(a*L1**4)+(L4+a*(4*L1*L3+3*L2**2))))
    return out


def safe_transition_kernels(template,s,s_cover,mu,cells=16):
    """Unchanged monotone integral rectangles, cancellation-free true masses."""
    c=template.ctx
    if type(cells) is not int or cells<1:raise ValueError('Positive original cell count required')
    if s.zero:return dict(values={k:s for k in ('J','theta','energy','pressure')},record=dict(
        exact_original_s_zero_all_transition_kernels_zero=True,cells=cells,cell_records=[]))
    if ep(s_cover)[0]<0 or ep(s_cover)[1]>ep(c.mpf('.5'))[0]:raise ValueError('Original microscopic left-half source required')
    J=template.scalar(0);K=[template.scalar(0) for _ in range(3)];rows=[];exp_records=[]
    width=s*(c.mpf(1)/cells)
    for index in range(cells):
        a=s*(c.mpf(index)/cells);b=s*(c.mpf(index+1)/cells)
        ac=s_cover*(c.mpf(index)/cells);bc=s_cover*(c.mpf(index+1)/cells)
        sa=original_sigma_rows(template,ac)[0];sb=original_sigma_rows(template,bc)[0]
        sj=positive_hull(sa,sb)
        nextJ=J+width*sj;jcell=positive_hull(J,nextJ)
        # exp(b)-exp(a) and exp(-a)-exp(-b), both exact positive masses.
        positive_increment=rc.density.density.factored_expm1(width)
        theta_weight=exponential_factor(a,exp_records)*positive_increment
        pressure_weight=exponential_factor(-b,exp_records)*positive_increment
        weights=(theta_weight,width,pressure_weight)
        for n,power in enumerate((1,2,2)):
            K[n]=K[n]+weights[n]*exponential_factor(-jcell*(mu*power),exp_records)
        rows.append(dict(exact_normalized_cell=[str(Fraction(index,cells)),str(Fraction(index+1,cells))],
            original_positive_cell_width=width.record(),original_left_s=a.record(),original_right_s=b.record(),
            original_left_sigma=sa.record(),original_right_sigma=sb.record(),
            same_original_sigma_monotone_cell_cover=sj.record(),original_J_cell_cover=jcell.record(),
            original_next_J=nextJ.record(),original_positive_kernel_weights=[v.record() for v in weights]))
        J=nextJ
    values=dict(J=J,theta=K[0],energy=K[1],pressure=K[2])
    return dict(values=values,record=dict(cells=cells,cell_records=rows,original_exponential_corrections=exp_records,
        original_kernel_values=records(values),positive_width_collected_before_endpoint_subtraction=True,
        exact_weight_identities=['exp(b)-exp(a)=exp(a)*expm1(b-a)','exp(-a)-exp(-b)=exp(-b)*expm1(b-a)'],
        original_sigma_and_J_integral_formulas_unchanged=True,
        all_positive_kernel_masses_retained=True,original_kernel_enclosures_not_selected_as_source_values=True))


def compile_original(stem,name,callbacks,method=False):
    module=__import__(PREFIX+stem)
    tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    candidates=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name]
    if len(candidates)!=1:raise ValueError('One original source function required')
    fn=copy.deepcopy(candidates[0]);before=ast.dump(fn);sites=[];divisions=[]
    # Only commutative orientation; mpf intervals must not try to coerce
    # a factored source jet before its own arithmetic method can run.
    class Orient(ast.NodeTransformer):
        def visit_BinOp(self,node):
            self.generic_visit(node)
            if isinstance(node.op,ast.Mult) and isinstance(node.left,ast.Name) and node.left.id=='mu':
                sites.append(ast.unparse(node));node.left,node.right=node.right,node.left
            if isinstance(node.op,ast.Div) and isinstance(node.right,ast.Constant) and type(node.right.value) is int and node.right.value>0:
                divisions.append(ast.unparse(node))
                inverse=ast.parse('c.mpf(1)/'+str(node.right.value),mode='eval').body
                node=ast.BinOp(left=node.left,op=ast.Mult(),right=inverse)
            return node
    Orient().visit(fn)
    environment=dict(vars(module));environment.update(callbacks)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<same original source, directed arithmetic callbacks>','exec'),environment)
    return environment[name],dict(source_module=Path(module.__file__).name,source_function=name,
        source_sha256=sha(Path(module.__file__).name),original_AST_sha256=hashlib.sha256(before.encode()).hexdigest(),
        commutative_mu_multiplications_reoriented=sites,
        exact_positive_integer_divisions_reoriented_as_right_scalar_factors=divisions,
        defining_source_assignments_not_replaced=True,original_provider_not_monkeypatched=True)


def canonical_amplitude_proof():
    parameter=binding.actual_parameter_formula_bindings()
    wanted={
        ('slope','factor'):"c.exp(y/10-c.mpf('.6')*J)",
        ('slope','u'):'qi*factor',('inlet','t_dummy'):None,
        ('axial','t'):'y-1',('axial','u'):'u1*root',
        ('axial','root'):'c.exp(-t/2)'}
    tree=ast.parse((HERE/(PREFIX+'pre_pulse_mixed_C4.py')).read_text(encoding='utf8'))
    for (name,target),text in wanted.items():
        if text is None:continue
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(text,mode='eval').body):raise ValueError('Original amplitude definition changed')
    slope_tree=ast.parse((HERE/'lei_ren_part1_paper_interval_outer_slope_field.py').read_text(encoding='utf8'))
    fn=next(n for n in slope_tree.body if isinstance(n,ast.FunctionDef) and n.name=='transition_integrals')
    endpoint=next(n for n in fn.body if isinstance(n,ast.If) and ast.unparse(n.test)=='y == 1')
    if ast.unparse(endpoint.body[0])!="j = c.mpf('.5')":raise ValueError('Original J(1)=1/2 identity changed')
    M,z=sy.symbols('Md Z')
    original_log=sy.Rational(1,10)-sy.Rational(3,5)/2-(sy.exp(M)+10)/2-sy.log(1+z*z)
    canonical=-(sy.exp(M)+11)/2+sy.Rational(3,10)-sy.log(1+z*z)
    if sy.simplify(original_log-canonical)!=0:raise ArithmeticError('Original Rd amplitude common factor identity')
    return dict(passed=True,original_parameter_binding=parameter,
        exact_original_Ud_identity='Utheta(Rd,Z)/Pstar=Pstar**(-1/2)*exp(3/10)/(1+Z**2)',
        exact_source_logPstar_equals_exp_Md_plus_11=True,
        exact_original_O2_exit_J_one_half=True,
        factor_correlations_collected_before_evaluation=True,
        source_functions_not_fitted_or_capped=True,
        input_hashes={**parameter['input_hashes'],PREFIX+'pre_pulse_mixed_C4.py':sha(PREFIX+'pre_pulse_mixed_C4.py'),
            'lei_ren_part1_paper_interval_outer_slope_field.py':sha('lei_ren_part1_paper_interval_outer_slope_field.py')})


class FactoredOriginalTransition:
    def __init__(self,owner,coordinates,mu,eta):
        self.owner=owner;self.c=coordinates.ctx;self.coordinates=coordinates;self.template=coordinates.scalar(1)
        self.pre=owner.transfer.owner.owner.owner.native.pre
        self.mu=mu;self.eta=eta;self.proof=canonical_amplitude_proof();self.parents={}
    def lift(self,row):return ScaledTaylor(self.template,[self.c.mpf(ep(v)) for v in row.coefficients])
    def source(self,Z,geometry,cells=16):
        c=self.c;t=self.template
        if geometry['kind']=='xi':
            x=c.mpf((ep(c.mpf(geometry['left'].numerator)/geometry['left'].denominator)[0],
                     ep(c.mpf(geometry['right'].numerator)/geometry['right'].denominator)[1]))
            s=previous.OriginalTransitionCoordinates(self.coordinates,self.mu,self.eta).W*x
        else:
            coord=previous.OriginalTransitionCoordinates(self.coordinates,self.mu,self.eta)
            k=c.mpf((ep(coord.cv(geometry['right']))[0],ep(coord.cv(geometry['left']))[1]))
            s=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=-c.ln(coord.D+k)/2),1,t.ledger)
        key=tuple(Z)
        if key not in self.parents:self.parents[key]=self.pre.axial(Z,buffer_offset=11)
        rawparent=self.parents[key]
        parent=dict(rawparent)
        parent['log_Utheta_over_Pstar_base_source']=t.scalar(c.mpf(rawparent['log_Utheta_over_Pstar_base_source']))
        parent['actual_normalized_primitive_y_derivative_axial5']={k:[self.lift(v) for v in rows] for k,rows in rawparent['actual_normalized_primitive_y_derivative_axial5'].items()}
        z=ScaledTaylor(t,[c.mpf(Z),1,0,0,0,0]);qi=(1+z*z).reciprocal()
        amp=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,-1,0)),c.exp(c.mpf('.3')),t.ledger)
        ud=qi*amp;parent['Utheta_over_Pstar_axial5_coefficients']=ud.coefficients
        kernels=safe_transition_kernels(t,s,geometry['coordinate'],self.mu,cells)
        sigma=original_sigma_rows(t,geometry['coordinate'])
        exponent_records=[]
        class Context:
            def mpf(_,value):return s if value is marker else ProtocolConstant(prior.FormalScale(t.scale.bases),c.mpf(value),t.ledger)
            def exp(_,value):return exponential_factor(value if isinstance(value,prior.ScaledEnclosure) else t.scalar(value),exponent_records)
        ctx=Context();marker=object()
        def factory(_,values):return ScaledTaylor(t,values)
        factory.variable=lambda _,value,order:ScaledTaylor(t,[value,1]+[0]*(order-1))
        phys,physproof=compile_original('pre_pulse_mixed_C4','physical_mixed',dict(
            IntervalTaylor=factory,square=lambda v:v*v,derivative=lambda v:v.derivative()))
        emitted={}
        def packet(_,ZZ,chart,offset,u,logu,logU,V,history,extra):
            datum=ScaledTaylor(t,[c.mpf(ep(v)) for v in rawparent['original_P0_axial5_coefficients']])
            invP2=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,-4,0)),1,t.ledger)
            emitted.update(physical=phys(ctx,ZZ,c.mpf(self.pre.delta),u,logU,V,history,datum,invP2),
                u=u,logU=logU,V=V,history=history,P0=datum,original_logu_cover=logu,extra=extra)
            return emitted
        proxy=SimpleNamespace(ctx=ctx,params=self.pre.params,delta=c.mpf(self.pre.delta),cells=cells,
            axial=lambda ZZ,buffer_offset:parent,coordinates=lambda ZZ:(z,qi),packet=None)
        proxy.packet=lambda *args:packet(proxy,*args)
        fn,proof=compile_original('pre_pulse_mixed_C4','slope_mu',dict(
            IntervalTaylor=factory,
            square=lambda value:value*value,
            transition_kernels=lambda ctx2,offset,mu,n:kernels['values'],
            sigma_jets=lambda ctx2,offset:[value*(c.mpf(1)/math.factorial(j)) for j,value in enumerate(sigma)],
            endpoints=lambda value:ep(bounded(value)) if isinstance(value,prior.ScaledEnclosure) else ep(value)))
        fn(proxy,Z,marker)
        emitted['original_s_source']=s
        return dict(values=emitted,record=dict(source_family=self.owner.family,exact_Z_range=list(Z),
            actual_original_typed_geometry=geometry['record'],same_original_source_program=proof,
            same_original_physical_mixed_program=physproof,canonical_original_amplitude_source=self.proof,
            lossless_actual_native_parent_inputs=dict(
                exact_Z_range=list(Z),original_mu=self.mu,original_eta_log=self.eta,original_delta=self.pre.delta,
                original_log_Utheta_over_Pstar_base_source=rawparent['log_Utheta_over_Pstar_base_source'],
                original_history_axial_Taylor_coefficients={key:[[c.mpf(ep(v)) for v in row.coefficients] for row in rows]
                    for key,rows in rawparent['actual_normalized_primitive_y_derivative_axial5'].items()},
                original_P0_axial_Taylor_coefficients=[c.mpf(ep(v)) for v in rawparent['original_P0_axial5_coefficients']]),
            original_inherited_Rd_history_rows=records(parent['actual_normalized_primitive_y_derivative_axial5']),
            original_separate_P0=records(emitted['P0']),original_safe_transition_kernels=kernels['record'],
            original_profile_exponential_corrections=exponent_records,
            original_current_amplitude=records(emitted['u']),original_logU_ordinary_y1_through_y4=records(emitted['logU']),
            actual_original_ordinary_velocity_history_and_pressure_rows=records(emitted['physical']),
            all_five_inherited_histories_and_analytic_P0_retained=True,
            original_Ud_common_Pstar_factor_collected=True,original_provider_graph_unmodified=True,
            ordinary_derivative_coordinate='y=log R and axial Z; no typed Jacobian reapplied'))


def signed_root_source(adapter,got,normalized,eta,log_a_lower,seed):
    c=adapter.c;t=adapter.template;v=got['values'];phys=v['physical']
    u=phys['physical_velocity_pressure_y_derivative_Taylor']['Utheta_over_Pstar']
    V=phys['physical_velocity_pressure_y_derivative_Taylor']['Uz']
    history=phys['actual_normalized_primitive_y_derivative_axial5']
    P=[v['P0']+history['p'][0]]+history['p'][1:]
    fn,proof=compile_original('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(axial_derivative=lambda x:x.derivative()))
    parts=fn(c,c.mpf(adapter.pre.delta),ScaledTaylor(t,[c.mpf(got['record']['exact_Z_range']),1,0,0,0,0]),u,V,history,P)
    Ps=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,2,0)),1,t.ledger)
    invPs=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,-2,0)),1,t.ledger)
    R0=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,22,0),
        10*c.mpf(seed.core.logC)+c.ln(110)),1,t.ledger)
    geometry=got['record']['actual_original_typed_geometry']
    scover=iv(c,encode(geometry['native_coordinate_box']))
    logRdelta=v['original_s_source'];radexp=[]
    radius=R0*exponential_factor(logRdelta,radexp)
    ns={}
    for label,omit in (('theta','variable_radial_shear'),('axial','axial_radial_shear')):
        raw=[u[0]*0 for _ in range(3)]
        for name,part in parts[label].items():
            if name==omit:continue
            factor=1 if label=='theta' else invPs if part['mode'][1]==0 else Ps
            for j in range(3):raw[j]=raw[j]+part['shape'][j]*factor
        ns[label]=[sum((raw[i]*(math.comb(j,i)) for i in range(j+1)),raw[0]*0)*radius for j in range(3)]
    roots={}
    for name,rows in (('E',u),('p1_numerator',ns['theta']),('p2_numerator',ns['axial'])):
        roots[name]={(j,k):rows[j][k]*math.factorial(k) for j,k in slow.ORDERS}
    lower=roots['E'][ZERO].scale.evaluate()+c.ln(c.mpf(ep(roots['E'][ZERO].coefficient)[0]))
    roots['p1']=current.quotient_jet(roots.pop('p1_numerator'),roots['E'],lower)
    roots['p2']=current.quotient_jet(roots.pop('p2_numerator'),roots['E'],lower)
    q=previous.correlated_original_q(normalized,eta,log_a_lower)
    roots.update(q['roots']);roots['b']={k:t.scalar(0) for k in slow.ORDERS};roots['t0']=dict(roots['b'])
    record=dict(source_family=adapter.owner.family,source_provenance=got['record']['same_original_source_program'],
        original_signed_stress_source_program=proof,
        original_signed_roots={name:{'y%d_Z%d'%k:value.record() for k,value in rows.items()} for name,rows in roots.items()},
        original_formal_radius=radius.record(),original_nonzero_radius_exponential_correction=radexp,
        full_pressure_P0_plus_p_added_exactly_once=True,
        original_p2_complete_energy_pressure_and_retained_m_terms_not_assumed_positive=True)
    return dict(roots=roots,q=q['rows'][ZERO],record=record,qrows=q['rows'],source_record=got['record'])


def original_local_density_integrals(adapter,got,roots,geometry,phase,dstar):
    """Original C0/Z density covers integrated over the true typed width."""
    t=adapter.template;c=adapter.c
    invPs=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,-2,0)),1,t.ledger)
    axial=got['values']['physical']['physical_velocity_pressure_y_derivative_Taylor']['Uz'][0]
    V,VZ=axial[0]*invPs,axial[1]*invPs
    E,EZ=roots['roots']['E'][ZERO],roots['roots']['E'][(0,1)]
    inverses=[];densities=[];density_records=[]
    for phi in phase['actual_phase_boxes']:
        inv=first.conditioned_first_jets(roots,roots['qrows'],dstar,phi)
        inverses.append(inv['record'])
        if inv['values'] is None:continue
        density=rc.density.density_Z_kernels(E,EZ,V,VZ,inv['values'],N)
        densities.append(density)
        density_records.append(dict(actual_original_fraction_phase_box=phi,
            original_candidate_common_Pstar_velocity_C0_Z=records(density['velocities']),
            original_five_signed_density_C0=records(density['kernels']),
            original_five_signed_density_Z=records(density['Z_derivatives'])))
    if len(densities)!=len(inverses):
        return dict(inverses=inverses,record=dict(status='requires_original_phase_source_refinement',
            exact_phase_union_not_partially_integrated=True,original_density_records=density_records))
    kernels={k:rc.density.local.same_source_union([d['kernels'][k] for d in densities]) for k in rc.RATES}
    jets={k:rc.density.local.same_source_union([d['Z_derivatives'][k] for d in densities]) for k in rc.RATES}
    factors={k:rc.transfer.true_width_kernel(adapter.coordinates,geometry,rate) for k,rate in rc.RATES.items()}
    values={k:kernels[k]*factors[k]['mass'] for k in rc.RATES}
    derivatives={k:jets[k]*factors[k]['mass'] for k in rc.RATES}
    return dict(inverses=inverses,record=dict(status='enclosed',
        original_density_records=density_records,original_complete_phase_union_density_C0=records(kernels),
        original_complete_phase_union_density_Z=records(jets),
        original_true_width_kernel_factors={k:dict(branch=x['branch'],mass=x['mass'].record(),decay=x['decay'].record()) for k,x in factors.items()},
        original_local_five_signed_C0_Duhamel_contributions=records(values),
        original_local_five_signed_Z_Duhamel_contributions=records(derivatives),
        original_common_velocity_unit='S=Pstar',original_five_common_velocity_units=first.packets.recovery.UNITS,
        original_axial_velocity_divided_by_Pstar_once=True,
        original_source_width_Jacobian_installed_once=True,ordinary_Z_endpoint_and_phase_derivatives_zero=True,
        full_original_prefix_incoming_or_Rc_targets_admitted=False))


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/previous.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,previous,accepted['source_family'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(previous.previous.generic.tail.baseline.build_owner(bridge)))
        native=previous.previous.PositiveLogTransition(owner)
        coordinates=HalfPstarCoordinates(owner.ctx,owner.coordinates.logP_squared,owner.family)
        coord=previous.OriginalTransitionCoordinates(coordinates,native.mu,native.eta)
        adapter=FactoredOriginalTransition(owner,coordinates,native.mu,native.eta)
        archives=[];summaries=[]
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            rows=[]
            for label,kind,left,right in (('bulk','xi','1/4','3/4'),('seam','k','7/4','5/3')):
                geometry=coord.geometry(kind,left,right);got=adapter.source(Z,geometry)
                norm=coord.normalized_excess(adapter.template,geometry)
                roots=signed_root_source(adapter,got,norm['rows'],native.eta,native.positive['log_actual_a_positive_lower'],adapter.owner.transfer.geometry.binder.seed)
                phase=previous.typed_original_phase(owner.transfer.geometry.binder,coordinates,geometry,Z,N)
                local=original_local_density_integrals(adapter,got,roots,geometry,phase,native.dstar)
                inverses=local['inverses']
                rows.append(dict(label=label,original_typed_native_source=got['record'],
                    original_correlated_q_source=previous.correlated_original_q(norm['rows'],native.eta,native.positive['log_actual_a_positive_lower'])['record'],
                    original_signed_root_source=roots['record'],actual_original_typed_phase=phase,
                    actual_original_conditioned_first_jet_queries=inverses,
                    actual_original_local_five_C0_Z_density_and_integral=local['record'],
                    original_full_prefix_integral_not_claimed=True))
                print('Safe original typed native source:',tag,label,[r['status'] for r in inverses],flush=True)
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                common_directed_coordinate_theorem=coordinates.record(),original_source_rows=rows,
                original_logC=adapter.owner.transfer.geometry.binder.seed.core.logC,
                original_dstar_log=native.dstar,original_positive_a_lower_log=native.positive['log_actual_a_positive_lower'])
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_transition_factored_native_source_'+tag+'.json.gz'
            (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
            archives.append(dict(filename=filename,compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
                lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),exact_Z_range=list(Z)))
            summaries.append(dict(exact_Z_range=list(Z),actual_source_queries=[dict(label=r['label'],
                p2_sign=r['original_signed_root_source']['original_signed_roots']['p2']['y0_Z0']['sign'],
                inverse_statuses=[x['status'] for x in r['actual_original_conditioned_first_jet_queries']],
                local_five_C0_Z_integral_status=r['actual_original_local_five_C0_Z_density_and_integral']['status']) for r in rows]))
        hashes.update(owner.service.hashes);hashes.update(adapter.proof['input_hashes'])
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        for name in (PREFIX+'outer_buffer.py',PREFIX+'outer_initial.py',PREFIX+'current_pre_pulse_stress_operator.py'):
            hashes[name]=sha(name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            original_factored_native_source_archives=archives,actual_original_query_summaries=summaries,
            same_original_source_amplitude_kernel_and_history_programs_executed=True,
            original_typed_local_five_signed_C0_and_Z_integral_contributions_executed=True,
            full_prefix_density_integrals_or_Rc_targets_admitted=False,
            **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,
            execution_seconds=time.monotonic()-began,source_seed_evidence=seed,
            scope='Same-original microscopic kernel and amplitude/history/ordinary physical rows; original complete signed p1/p2 roots, actual phase/conditioned first jets and local five signed C0/Z density/integral covers on two symmetric interior Z tiles. Undetermined p2 sign is retained; small-r inversion needs no sign assumption. Full prefix/axial/targets/controls/global N/stress/recursion/full NS remain open.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
