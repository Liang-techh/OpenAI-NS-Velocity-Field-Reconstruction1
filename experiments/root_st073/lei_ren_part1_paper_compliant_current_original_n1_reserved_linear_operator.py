"""Original positive-order linear repair on a fresh, disjoint Rc support.

This constructs normalized bump functions, physical moment matrices and a
linear C1 inverse action. It does not select the enclosed defect values as
functions, install repaired fields, or close a pressure derivative gate.
"""
import copy
from fractions import Fraction
from functools import wraps
import gzip
import json
import math
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as defects
import lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp as quiet
import lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame as frame
import lei_ren_part1_paper_compliant_current_native_Rc_C1_integral_realization as realization
import lei_ren_part1_paper_compliant_current_generic_shear_source_packets as packets

HERE,PREFIX,sha,ep=defects.HERE,defects.PREFIX,defects.sha,defects.ep
record,interval=defects.record,defects.interval
source,controls=realization.source,realization.controls
NAME=PREFIX+'current_original_n1_reserved_linear_operator.json.gz'
RECEIPT=PREFIX+'current_original_n1_reserved_linear_operator_check.json'
GATE='current_original_n1_reserved_normalized_linear_repair_operator_installed'
OPEN=tuple(dict.fromkeys(defects.OPEN+frame.OPEN+(
    'n1_complete_signed_defect_C1_function_adapter_installed',
    'n1_selected_linear_coefficient_functions_installed',
    'n1_repaired_own_M_V_P_installed',
    'positive_order_five_moment_repair_closed',
    'uniform_full_pressure_derivative_accuracy_certified')))
CENTERS=(12,13,14,15,16)
AXIAL=(0,4)
ANGULAR=(1,2,3)
DPS=100


def precision(fn):
    @wraps(fn)
    def call(*a,**k):
        with mp.workdps(DPS):return fn(*a,**k)
    return call


def read(name):
    data=(HERE/name).read_bytes()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)


def inputs():
    reports,hashes={},{}
    for module in (defects,defects.external,quiet.band,quiet.band.control,quiet,frame):
        raw,receipt=read(module.NAME),read(module.RECEIPT)
        gates=module.GATES if hasattr(module,'GATES') else (module.GATE,)
        if not receipt['all_passed'] or not all(receipt.get(k) for k in gates):
            raise ValueError('Checked same-family input required '+module.NAME)
        for name,digest in receipt['input_hashes'].items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Conflicting source '+name)
            hashes[name]=digest
        for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):hashes[name]=sha(name)
        reports[module.NAME]=raw
    family=reports[defects.NAME]['source_family']
    if any(r['source_family']!=family for r in reports.values()):raise ValueError('Same original source family required')
    for module in (realization,controls,source,packets):hashes[Path(module.__file__).name]=sha(Path(module.__file__).name)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Changed accepted source '+name)
    return reports,hashes,family


def exact_supports():
    """All endpoints are multiples of log(2), so rational tests are exact."""
    new=[(Fraction(j,40)-Fraction(1,400),Fraction(j,40)+Fraction(1,400)) for j in CENTERS]
    old=[(Fraction(j-1,40),Fraction(j+1,40)) for j in (8,20,32)]
    gap=(Fraction(11,40),Fraction(17,40))
    if not all(gap[0]<a<b<gap[1] for a,b in new):raise ValueError('Support left exact reserved gap')
    if not all(new[i][1]<new[i+1][0] for i in range(4)):raise ValueError('New supports overlap')
    if not all(b<c or d<a for a,b in new for c,d in old):raise ValueError('Old correction intersects n1 support')
    return dict(new_log_supports_over_log2=[[str(a),str(b)] for a,b in new],
        old_log_supports_over_log2=[[str(a),str(b)] for a,b in old],
        exact_reserved_gap_over_log2=list(map(str,gap)),autonomous_support_choice=True,
        axial_indices=list(AXIAL),angular_indices=list(ANGULAR),
        old_corrections_and_all_local_derivatives_zero_on_new_support=True,
        nonlocal_histories_pressure_and_density_not_zeroed=True)


def arithmetic(nodes,node,bindings):
    if node in bindings:return bindings[node]
    row=nodes[node];at=lambda i:arithmetic(nodes,i,bindings);op=row['operation']
    if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
    if op=='sum':return sy.Add(*map(at,row['arguments']))
    if op=='product':return sy.Mul(*map(at,row['arguments']))
    if op=='negative':return -at(row['argument'])
    if op=='positive_quotient':return at(row['numerator'])/at(row['denominator'])
    if op=='analytic_unary' and row['name'] in ('log','exp'):return getattr(sy,row['name'])(at(row['argument']))
    raise ValueError('Unsupported support/constant expression '+op)


def bind_background(reports,g):
    b=reports[quiet.band.NAME];q=reports[quiet.NAME]
    nodes=b['exact_graph_nodes'];f=b['actual_C3_repair_band_functions'];p=f['profiles_at_band_x']
    x=sy.Symbol('x',positive=True)
    raw_beta=[]
    for j,bump in zip((8,20,32),p['bumps']):
        row=nodes[bump]
        if row['operation']!='positive_quotient':raise ValueError('Actual compact bump recipe changed')
        raw=nodes[row['numerator']]
        if raw['operation']!='compact_raw_beta':raise ValueError('Actual lazy compact source required')
        arg=arithmetic(nodes,raw['argument'],{f['original_band_variable']:x})
        expected=(sy.log(x)-sy.Rational(j,40)*sy.log(2))/(sy.log(2)/40)
        if sy.simplify(arg-expected)!=0:raise ValueError('Actual old support does not match reserved gap')
        raw_beta.append(row['numerator'])
    binding=q['actual_range_source_node_binding'];amplitudes=binding['actual_amplitude_source_nodes']
    for order,node in enumerate(amplitudes):
        leaf=nodes[node]
        if (leaf['operation']!='current_original_leading_function_recipe' or leaf['source_family']!=b['source_family']
            or leaf['native_chart']!='O3_power' or leaf['quantity']!='E' or leaf['Z_order']!=order
            or leaf['recipe']['quantity_paths']['E']!='original_generic_source.common_velocity_E_axial5[Z_order]'
            or leaf.get('Taylor_coefficient_factorial',1)!=(1,1,2,6)[order]
            or nodes[leaf['coordinate']]!=dict(operation='exact_rational',numerator=2,denominator=1)):
            raise ValueError('Genuine original ordinary amplitude rows required')
    base=reports[quiet.band.control.NAME]['actual_C3_limit_controls_and_Picard_functions']['exact_C2_limit_functions']['exact_C1_limit_adapter']
    ref=lambda i:source.FunctionRef(g,i);params={k:ref(i) for k,i in base['parameters'].items()}
    logminus=g.add(g.unary('log',g.constant(4)),g.neg(g.mul(g.constant(4),params['logP'])),
        g.constant(-1000),g.mul(params['hbB'],params['sc'],g.constant('1/2')))
    logRc=g.add(ref(base['original_Rc_offset']),logminus)
    log2Rc=g.add(g.unary('log',g.constant(2)),logRc)
    alpha=g.add(g.constant('1/2'),params['mu'])
    factor=g.unary('exp',g.add(params['logP'],g.mul(alpha,log2Rc)))
    A=[g.mul(factor,ref(i)) for i in amplitudes]
    # E is already Utheta/Pstar in the actual common-source frontend.
    # Thus Utheta=Pstar*A_rc*x^-alpha=Astar*(2R)^-alpha.
    return dict(parameters=base['parameters'],logRc=logRc.node,log2Rc=log2Rc.node,
        Astar_ordinary_Z0_Z3=[r.node for r in A],actual_A_rc_ordinary_Z0_Z3=amplitudes,
        amplitude_positive_source_cells=[dict(exact_Z_cell=r['exact_Z_cell'],
            actual_A_rc_ordinary_Z0_Z3=r['actual_endpoint_source']['ordinary_terminal_amplitude_C3'])
            for r in b['actual_four_Z_C3_band_ranges']],
        physical_common_unit='Pstar',common_scale_applied_once=True,
        original_power_definition='Utheta0=Pstar*A_rc*(R/Rc)^(-1/2-mu)=Astar*(2R)^(-1/2-mu)',
        old_compact_raw_beta_nodes=raw_beta,old_profile_functions=p,
        actual_reserved_band_variable=f['original_band_variable'],
        actual_mu_log=binding['actual_mu_log'],same_source_actual_angular_amplitude=True,
        nonlocal_histories_retained=True)


def graph_weights(g,background):
    ref=lambda i:source.FunctionRef(g,i)
    mu=ref(background['parameters']['mu']);L=g.unary('log',g.constant(2));ell=g.mul(L,g.constant('1/400'))
    centers=[g.mul(L,g.constant(Fraction(j,40))) for j in CENTERS]
    t=g.symbol('n1_reserved_normalized_bump_t');s=g.symbol('n1_reserved_log_R_over_Rc')
    beta=g.node('compact_raw_beta',argument=t.node,definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
        defining_module=PREFIX+'outer_pulse_map.py',defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))
    integ=lambda a:controls.integral(g,a,'n1_reserved_normalized_bump_t',g.constant(-1),g.one,
        measure='dt; the translated physical normalizer retains ell*exp(center/2)')
    J=lambda k:integ(g.mul(beta,g.unary('exp',g.mul(k,ell,t))))
    normal=J(g.constant('1/2'));ratio=lambda a:g.quotient(a,normal,'positive integral beta(t)*exp(ell*t/2) dt')
    H=lambda k:ratio(J(g.add(k,g.constant('1/2')))) # M_j(k)=exp(k*center_j)*H(k).
    half=g.constant('1/2');H0=H(half);Hm=H(g.sub(half,mu))
    rawD=ratio(integ(g.mul(beta,g.unary('exp',g.mul(ell,t)),g.neg(g.mul(ell,t)),
        g.unary('exprel',g.neg(g.mul(mu,ell,t))))))
    axial=[g.mul(g.unary('exp',g.mul(half,centers[j])),H0) for j in AXIAL]
    D=[g.mul(g.unary('exp',g.mul(half,centers[j])),g.add(rawD,
        g.mul(Hm,g.neg(centers[j]),g.unary('exprel',g.neg(g.mul(mu,centers[j])))))) for j in AXIAL]
    gap=g.sub(centers[AXIAL[1]],centers[AXIAL[0]])
    positive_axial=g.mul(H0,Hm,g.unary('exp',g.sub(g.mul(half,g.add(centers[0],centers[4])),
        g.mul(mu,centers[0]))),gap,g.unary('exprel',g.neg(g.mul(mu,gap))))
    powers=[g.one,g.neg(g.add(g.one,mu)),g.neg(mu)]
    angular=[[g.mul(H(k),g.unary('exp',g.mul(k,centers[j]))) for j in ANGULAR] for k in powers]
    step=g.sub(centers[2],centers[1])
    angular_det=g.mul(*(H(k) for k in powers),g.unary('exp',g.mul(g.add(*powers),centers[1])))
    for i,j in ((0,1),(0,2),(1,2)):
        delta=g.sub(powers[j],powers[i])
        angular_det=g.mul(angular_det,g.unary('exp',g.mul(powers[i],step)),delta,step,
            g.unary('exprel',g.mul(delta,step)))
    inverse=[]
    for i in range(3):
        row=[]
        for j in range(3):
            ii=[k for k in range(3) if k!=j];jj=[k for k in range(3) if k!=i]
            cofactor=g.mul(g.constant((-1)**(i+j)),g.sub(g.mul(angular[ii[0]][jj[0]],angular[ii[1]][jj[1]]),
                g.mul(angular[ii[0]][jj[1]],angular[ii[1]][jj[0]])))
            row.append(g.quotient(cofactor,angular_det,'positive new translated angular Vandermonde determinant'))
        inverse.append(row)
    logRc=ref(background['logRc']);log2Rc=ref(background['log2Rc'])
    physical_bumps=[]
    for center in centers:
        arg=g.quotient(g.sub(s,center),ell,'positive new support half-width')
        raw=g.node('compact_raw_beta',argument=arg.node,definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
            defining_module=PREFIX+'outer_pulse_map.py',defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))
        pref=g.unary('exp',g.mul(half,g.sub(L,logRc)))
        physical_bumps.append(g.quotient(g.mul(pref,raw),g.mul(ell,g.unary('exp',g.mul(half,center)),normal),
            'actual positive physical dR/sqrt(2R) normalization'))
    return dict(mu=mu,ell=ell,centers=centers,normal=normal,physical_bumps=physical_bumps,
        bump_log_radius_variable=s,axial_rows=[axial,D],angular_rows=angular,
        positive_axial_divided_determinant=positive_axial,positive_angular_determinant=angular_det,
        inverse_angular=inverse,log2Rc=log2Rc,
        axial_physical_row_log_factors=[g.mul(half,log2Rc),g.mul(g.sub(half,mu),log2Rc)],
        angular_physical_row_log_factors=[log2Rc,g.neg(g.mul(g.add(g.one,mu),log2Rc)),g.neg(g.mul(mu,log2Rc))])


def inverse_normalized_action(g,W,y):
    g.ids(y)
    if len(y)!=5:raise ValueError('Five signed function handles required')
    a,b=W['axial_rows'];den=W['positive_axial_divided_determinant'];ym,yd,*theta=y
    result=[g.quotient(g.sub(g.mul(a[1],yd),g.mul(b[1],ym)),den,'positive negative-of-axial-divided determinant'),
        g.quotient(g.sub(g.mul(b[0],ym),g.mul(a[0],yd)),den,'positive negative-of-axial-divided determinant')]
    return result+[g.add(*(g.mul(v,q) for v,q in zip(row,theta))) for row in W['inverse_angular']]


def coefficient_action(g,W,background,moments):
    """Physical Section 15 signs, for supplied genuine signed C1 handles.

    The cancellation before division by mu stays an exact graph expression.
    Calling this with caps, affine midpoint values or another graph fails.
    """
    if set(moments)!=set(('m11','m12','m13','m14','m15')):raise ValueError('All five signed moments required')
    if any(type(p) is not source.C1Function for p in moments.values()):raise TypeError('C1 function handles required')
    for p in moments.values():g.ids([p.value,p.Z])
    ref=lambda i:source.FunctionRef(g,i);Av,Az=map(ref,background['Astar_ordinary_Z0_Z3'][:2])
    divideA=lambda p:source.C1Function(g.quotient(p.value,Av,'actual positive Astar'),
        g.quotient(g.sub(g.mul(p.Z,Av),g.mul(p.value,Az)),g.mul(Av,Av),'actual positive Astar squared'))
    scaled=lambda p,k:g.c1scale(k,p);log2Rc=W['log2Rc'];mu=W['mu']
    factor=lambda k:g.unary('exp',g.mul(k,log2Rc))
    invsr=factor(g.constant('-1/2'))
    y0=scaled(moments['m11'],g.neg(invsr))
    ym=scaled(divideA(moments['m14']),g.neg(g.mul(invsr,factor(mu))))
    yd=scaled(g.c1add(ym,scaled(y0,g.constant(-1))),g.quotient(g.one,mu,'same original strictly positive mu'))
    theta=[scaled(moments['m12'],g.neg(factor(g.constant(-1)))),
        scaled(divideA(moments['m13']),g.mul(g.constant('-1/2'),factor(g.add(g.one,mu)))),
        scaled(divideA(moments['m15']),factor(mu))]
    pairs=[y0,yd,*theta]
    values=inverse_normalized_action(g,W,[p.value for p in pairs]);derivatives=inverse_normalized_action(g,W,[p.Z for p in pairs])
    return [source.C1Function(v,d) for v,d in zip(values,derivatives)]


def exprel(c,x):
    """Directed entire-function evaluation without tiny exp(x)-1 subtraction."""
    lo,hi=ep(x);bound=c.mpf(max(abs(lo),abs(hi)))
    if ep(bound)[1]>.1:return (c.exp(x)-1)/x
    value=sum((x**j/math.factorial(j+1) for j in range(9)),c.mpf(0))
    err=c.exp(bound)*bound**9/math.factorial(10)
    return value+c.mpf([-ep(err)[1],ep(err)[1]])


def integral_enclosure(c,k,ell,first_moment=False):
    """Finite Taylor integration with whole-cell omitted derivatives and tails."""
    cut=Fraction(127,128);order=12;initial=48;total=c.mpf(0);errors=c.mpf(0);leaves=0
    rat=lambda q:c.mpf(q.numerator)/q.denominator
    def integrand(t):
        result=(t*(k*ell)-1/(1-t*t)).exp()
        return result*t if first_moment else result
    def cell(a,b,depth):
        nonlocal total,errors,leaves
        mid=(a+b)/2;radius=rat((b-a)/2)
        point=integrand(IntervalTaylor.variable(c,rat(mid),order))
        box=c.mpf([ep(rat(a))[0],ep(rat(b))[1]])
        remainder=integrand(IntervalTaylor.variable(c,box,order))[order]
        absolute=c.mpf(max(map(abs,ep(remainder))))
        err=2*absolute*radius**(order+1)/(order+1)
        if ep(err)[1]>mp.mpf('1e-18')*ep(rat(b-a))[0] and depth<12:
            cell(a,mid,depth+1);cell(mid,b,depth+1);return
        if ep(err)[1]>mp.mpf('1e-16')*ep(rat(b-a))[0]:raise ArithmeticError('New bump Taylor budget exceeded')
        total+=sum((2*point[j]*radius**(j+1)/(j+1) for j in range(0,order,2)),c.mpf(0))
        errors+=err;leaves+=1
    for i in range(initial):cell(-cut+2*cut*i/initial,-cut+2*cut*(i+1)/initial,0)
    tail=2*rat(1-cut)*c.exp(-1/(1-rat(cut)**2)+abs(k*ell))
    error=errors+tail
    value=total+c.mpf([-ep(error)[1],ep(error)[1]])
    return value,dict(Taylor_order=order,leaves=leaves,interior_remainder=errors,
        endpoint_strip_bound=tail,whole_cell_omitted_derivative=True,
        nonzero_quadrature_error=error,first_moment=first_moment)


def numeric_weights(c,mu):
    if ep(mu)[0]<0 or ep(mu)[1]>.5:raise ValueError('mu in [0,1/2] required')
    L=c.ln(2);ell=L/400;centers=[j*L/40 for j in CENTERS];cache={};diags={}
    def J(k,first=False):
        key=(k._mpi_,first)
        if key not in cache:cache[key],diags[str(len(diags))]=integral_enclosure(c,k,ell,first)
        return cache[key]
    normal=J(c.mpf('.5'));H=lambda k:J(k+c.mpf('.5'))/normal
    if ep(normal)[0]<=0:raise ArithmeticError('Positive new normalizer lost')
    H0=H(c.mpf('.5'));Hm=H(c.mpf('.5')-mu)
    axial=[c.exp(centers[j]/2)*H0 for j in AXIAL]
    # -D is a positive integral of s*exp(s/2)*exprel(-mu*s).
    D=[]
    for j in AXIAL:
        positive=c.exp(centers[j]/2)*(centers[j]*J(c.mpf(1))+ell*J(c.mpf(1),True))/normal
        factor=c.mpf([ep(c.exp(-mu*(centers[j]+ell)))[0],1])
        D.append(-positive*factor)
    gap=centers[4]-centers[0]
    det_ax=H0*Hm*c.exp((centers[0]+centers[4])/2-mu*centers[0])*gap*exprel(c,-mu*gap)
    powers=[c.mpf(1),-1-mu,-mu]
    angular=[[H(k)*c.exp(k*centers[j]) for j in ANGULAR] for k in powers]
    step=centers[2]-centers[1]
    det_ang=math.prod(H(k) for k in powers)*c.exp(sum(powers)*centers[1])
    for i,j in ((0,1),(0,2),(1,2)):
        d=powers[j]-powers[i];det_ang*=c.exp(powers[i]*step)*d*step*exprel(c,d*step)
    if ep(det_ax)[0]<=0 or ep(det_ang)[0]<=0:raise ArithmeticError('Fresh matrix determinant not positive')
    return dict(mu=mu,ell=ell,centers=centers,base_normalizer=normal,
        physical_normalizer_definition='J_j=ell*exp(center_j/2)*J(1/2); b_j=sqrt(2/Rc)*beta/J_j',
        axial_divided_rows=[axial,D],angular_rows=angular,
        positive_negative_axial_divided_determinant=det_ax,positive_angular_determinant=det_ang,
        actual_axial_determinant_definition='-mu*(2Rc)^(1-mu)*positive_negative_axial_divided_determinant',
        actual_angular_determinant_definition='(2Rc)^(-2mu)*positive_angular_determinant',
        Taylor_integral_diagnostics=diags,quadrature_error_nonzero=True,
        finite_mu_preserved=True,small_mu_divided_difference_evaluated_without_cancellation=True)


def encode(value):
    if type(value) is source.FunctionRef:return value.node
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value


class CurrentOriginalN1ReservedLinearOperator:
    @precision
    def __init__(self,require_checked=True):
        self.reports,self.hashes,self.family=inputs();self.c=MPIntervalContext();self.c.dps=DPS
        band=self.reports[quiet.band.NAME];self.graph=realization.restore_graph(band['exact_graph_nodes'])
        self.prefix=copy.deepcopy(self.graph.nodes);self.support=exact_supports()
        self.background=bind_background(self.reports,self.graph)
        self.functions=graph_weights(self.graph,self.background)
        for row in self.background['amplitude_positive_source_cells']:
            value=row['actual_A_rc_ordinary_Z0_Z3'][0]
            if value['sign']!='positive' or ep(interval(self.c,value['coefficient_interval']))[0]<=0:
                raise ArithmeticError('Actual Astar positivity missing on axial cell')
        mu=self.c.exp(interval(self.c,self.background['actual_mu_log']))
        self.values=numeric_weights(self.c,mu)
        self._immutable_nodes=copy.deepcopy(self.graph.nodes)
        self._immutable_bindings=json.dumps((self.background,encode(self.functions),self.support,record(self.values)),sort_keys=True)
        self.acceptance_loaded=False
        if require_checked:
            receipt=read(RECEIPT);raw=read(NAME)
            if not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt[k] for k in OPEN):
                raise ValueError('Checked original n1 linear basis required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed linear repair source '+name)
            if raw['source_family']!=self.family or raw['exact_graph_nodes']!=self.graph.nodes:
                raise ValueError('Accepted original linear source graph differs')
            self.acceptance_loaded=True

    @precision
    def assert_graph(self):
        if (self.graph.nodes[:len(self._immutable_nodes)]!=self._immutable_nodes
            or json.dumps((self.background,encode(self.functions),self.support,record(self.values)),sort_keys=True)!=self._immutable_bindings
            or sha(Path(__file__).name)!=self.hashes[Path(__file__).name]):
            raise ValueError('Original n1 linear operator or source bindings changed')

    @precision
    def coefficient_functions(self,moments):
        self.assert_graph()
        return coefficient_action(self.graph,self.functions,self.background,moments)


@precision
def run():
    began=time.monotonic();owner=CurrentOriginalN1ReservedLinearOperator(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.graph.nodes,
        accepted_background_graph_prefix_length=len(owner.prefix),
        actual_background_binding=owner.background,support_certificate=owner.support,
        exact_normalized_basis_and_matrix_functions=encode(owner.functions),
        actual_normalized_matrix_enclosures=record(owner.values),
        physical_moment_equations=['Bz*alpha=-[m11,m14/Astar]',
            'Btheta*beta=-[m12,m13/(2*Astar),-m15/Astar]'],
        full_signed_defect_values_archive=defects.NAME,
        coefficient_API='coefficient_functions requires five same-graph genuine signed C1 functions',
        enclosed_values_not_selected_as_functions=True,selected_coefficient_functions_installed=False,
        complete_pressure_C1_accuracy_installed=False,ancestor_producers_and_checkers_not_rerun=True,
        **dict.fromkeys((GATE,)+OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('Original n1 physical bump matrices and stable linear C1 inverse action constructed',flush=True)
    return owner,raw


if __name__=='__main__':run()
