"""Four exact current angular support traces with full weighted histories.

The source limit is reduced before interval evaluation. Flat beta input does
not erase the exponential background or accumulated angular/pressure/energy.
This stage identifies traces; it does not certify global stress or recursion.
"""
import ast
import copy
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import (
    CurrentPulseSupportInterfaces,weighted_support_FTC_theorem,physical_summary,
    fixed_normalization_and_angular_weight_source_proof,FlatPulseDerivatives,read_interval,
    EDGES,OPEN,HERE,PREFIX,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import (
    statement,current_postpulse_trace_proof,physical_operator_trace_proof)
from lei_ren_part1_paper_compliant_current_pulse_interfaces import EndpointDispatchView,BASE
from lei_ren_part1_paper_compliant_current_selected_energy_source import function,binding
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_current_physical_interfaces import supported_beta_edges
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_angular_support_interfaces.json'
RECEIPT=PREFIX+'current_angular_support_interfaces_check.json'
GATES=('current_angular_four_support_functional_mixed4_interfaces_identified',
       'current_angular_support_physical_spatial4_time1_traces_identified',
       'all_eight_current_internal_support_source_traces_identified')
WEIGHTS={'A':('1-mu',1),'B':('-1-2mu',1),'D':('-1-2mu',2),
         'E':('-2mu',1),'F':('-2mu',2)}


def angular_endpoint_theorem():
    """Actual source expressions and mixed recurrences for arbitrary C5 data."""
    z,mu,x,L,Ev2=s.symbols('Z mu source_s Lrel Ev2',real=True)
    bp=s.Rational(1,2)+mu;rate=1-mu;prate=1+2*mu
    arbitrary=lambda label:s.Function('current_'+label)(z)
    Xf,Pf,Post=map(arbitrary,('flatten_X','absolute_flatten_P','complete_post_future'))
    A,Q,H=map(arbitrary,('retained_past_A','retained_past_B_D','retained_future_E_F'))
    F=s.Symbol('F',positive=True);D=lambda a,y:(1-s.exp(-a*y))/a
    Xbase=1/rate+(Xf-1/rate)*s.exp(-rate*(L+x))
    env={'F[0]':F,'Xbase':Xbase,'pastA':A,'pastP':Q,'futureE':H,
         'self.rate':rate,'self.mu':mu,'self.bp':bp,'self.prate':prate,
         'self.Lrel':L,'self.flatten.Ev2':Ev2,'s':x,'y':L+x,
         "data['post']":Post,"data['flatten_exit_pressure']['P_over_Pstar_squared']":Pf,
         'decay_integral(c, 2 * self.mu, -s)':D(2*mu,-x),
         'baseline':D(prate,L+x)*Ev2*s.exp(-100*prate)/8}
    actual={label:assignment('compliant_power_angular_C4','angular',target,env)
            for label,target in (('theta','theta'),('X','X'),('energy','energy'),('pressure','pressure'))}
    actual['pressure']+=assignment('compliant_power_angular_C4','angular','pressure',env,augmented=True)
    outside={'theta':s.exp(-bp*(100+L+x))/2,
             'X':Xbase+A*s.exp(-rate*x),
             'energy':((Post+H)*s.exp(2*mu*x)+D(2*mu,-x))/2,
             'pressure':Pf+env['baseline']+Q*Ev2*s.exp(-prate*(100+L))/4}
    # Quotient recurrence and the actual mixed transport statements are
    # bound separately from symbolic boundary substitutions.
    bindings={}
    specs=(('compliant_power_angular_C4','angular','pastA+=dj*(c.exp(self.rate*center)*past[\'A\'])'),
        ('compliant_power_angular_C4','angular','pastP+=(dj*past[\'B\']+dj*dj*past[\'D\']/2)*c.exp(-self.prate*center)'),
        ('compliant_power_angular_C4','angular','futureE+=(dj*(2*future[\'E\'])+dj*dj*future[\'F\'])*c.exp(-2*self.mu*center)'),
        ('compliant_power_angular_C4','quotient_log_rates','value=rows[n+1]'),
        ('compliant_power_angular_C4','quotient_log_rates','for j in range(1,n+1):value-=rows[j]*rates[n-j]*math.comb(n,j)'),
        ('compliant_power_angular_C4','quotient_log_rates','rates.append(value/rows[0])'),
        ('compliant_flatten_mixed_C4','flatten_mixed','lograte[0]=lograte[0]-bp'),
        ('compliant_flatten_mixed_C4','flatten_mixed','angularrate[0]=angularrate[0]+1-mu'),
        ('compliant_flatten_mixed_C4','flatten_mixed','energyrate[0]=energyrate[0]+2*mu'),
        ('compliant_flatten_mixed_C4','flatten_mixed','th+=lograte[j]*theta_rows[k-j]*math.comb(k,j)'),
        ('compliant_flatten_mixed_C4','flatten_mixed','xx-=angularrate[j]*Xrows[k-j]*math.comb(k,j)'),
        ('compliant_flatten_mixed_C4','flatten_mixed','ee+=energyrate[j]*erows[k-j]*math.comb(k,j)'),
        ('compliant_flatten_mixed_C4','flatten_mixed','square+=theta_rows[j]*theta_rows[k-j]*math.comb(k,j)'),
        ('compliant_flatten_mixed_C4','flatten_mixed','theta_rows.append(th); Xrows.append(xx); erows.append(ee); prows.append(square*(Ev2/2))'))
    for stem,method,text in specs:
        for part in text.split(';'):bindings[stem+'.'+method+':'+part.strip()]=statement(stem,method,part.strip())
    binding('compliant_power_angular_C4','_packet','mixed',
        'flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    # F=1, derivatives1..4=0 implies the four logarithmic rate rows vanish.
    rates=[]
    for n in range(4):
        value=s.Integer(0)
        for j in range(1,n+1):value-=s.Integer(0)*rates[n-j]*math.comb(n,j)
        rates.append(value)
    if any(rates):raise ArithmeticError('Flat angular source did not eliminate rate jets')
    primitive={};velocity={}
    for edge in EDGES:
        v=s.Rational(edge['exact_edge'])
        a={label:[value.subs(F,1).subs(x,v)] for label,value in actual.items()}
        b={label:[value.subs(x,v)] for label,value in outside.items()}
        for order in range(4):
            for rows in (a,b):
                rows['theta'].append(-bp*rows['theta'][order])
                rows['X'].append((1 if order==0 else 0)-rate*rows['X'][order])
                rows['energy'].append((-s.Rational(1,2) if order==0 else 0)+2*mu*rows['energy'][order])
                rows['pressure'].append(Ev2/2*sum(math.comb(order,j)*rows['theta'][j]*rows['theta'][order-j] for j in range(order+1)))
        for label in a:
            for order in range(5):
                diff=s.simplify(s.expand_power_exp(a[label][order]-b[label][order]))
                if diff!=0:raise ArithmeticError('Actual angular endpoint history differs: '+label)
                for n in range(5-order):
                    if s.diff(diff,z,n)!=0:raise ArithmeticError('Axial angular trace differs')
                    primitive[edge['exact_edge']+'_'+label+'_y%d_Z%d'%(order,n)]=True
        for label in ('UZ','Utheta','UR','P'):
            source='theta' if label=='Utheta' else 'pressure' if label=='P' else None
            for order in range(5):
                for n in range(5-order):
                    diff=s.Integer(0) if source is None else a[source][order]-b[source][order]
                    if s.simplify(s.diff(s.expand_power_exp(diff),z,n))!=0:raise ArithmeticError('Angular velocity/absolute pressure trace differs')
                    velocity[edge['exact_edge']+'_'+label+'_y%d_Z%d'%(order,n)]=True
    return dict(actual_angular_weighted_translation_and_mixed_recurrence_AST=bindings,
        current_primitive_mixed4_identities=primitive,current_velocity_pressure_mixed4_identities=velocity,
        arbitrary_current_C5_coefficient_and_nonzero_history_functions_retained=True,
        five_weighted_integral_endpoint_values_and_FTC_derivatives_substituted_before_axial_differentiation=True,
        rates_vanish_but_exponential_background_derivatives_do_not=True,
        quadratic_D_F_pressure_energy_terms_preserved=True,absolute_pressure_uses_Pstar_squared=True,passed=True)


def source_edge_integrals(owner,center,edge):
    """Exact support ordering; no rounded clipping selects a source value."""
    ell=s.Rational(3,20);position=s.Rational(edge)
    zero={key:owner.ctx.mpf(0) for key in owner.weights}
    if position<=s.Integer(center)-ell:return zero,owner.weights
    if position>=s.Integer(center)+ell:return owner.weights,zero
    raise ValueError('Source edge unexpectedly lies inside another support')


def exact_endpoint_program():
    """Private replay of the unchanged recipe after proved source limits.

    Only beta jets and the numerical past/future selection are substituted.
    All formulas, owners, fixed units and retained histories stay original.
    """
    fn=copy.deepcopy(function('compliant_power_angular_C4','angular'))
    fn.name='_source_edge_angular';fn.args.args.append(ast.arg(arg='source_edge'))
    beta_replacements=integral_replacements=0
    class Reduce(ast.NodeTransformer):
        def visit_Assign(self,node):
            nonlocal beta_replacements
            if len(node.targets)==1 and ast.unparse(node.targets[0])=='beta':
                if ast.unparse(node.value)!='self.flat.beta(local)':raise ValueError('Angular beta acquisition changed')
                beta_replacements+=1
                return ast.parse('beta=[c.mpf(0) for _ in range(5)]').body[0]
            return self.generic_visit(node)
        def visit_If(self,node):
            nonlocal integral_replacements
            if ast.unparse(node.test)=='hi <= -endpoints(ell)[1]':
                integral_replacements+=1
                return ast.parse('past,future=source_edge_integrals(self,center,source_edge)').body[0]
            return self.generic_visit(node)
    fn=Reduce().visit(fn)
    if (beta_replacements,integral_replacements)!=(1,1):raise ValueError('Source endpoint reduction differs')
    tree=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))
    from lei_ren_part1_paper_compliant_power_angular_C4 import CompliantPowerAngularC4
    env=dict(CompliantPowerAngularC4.angular.__globals__);env['source_edge_integrals']=source_edge_integrals
    exec(compile(tree,'<source-proved-current-angular-support>','exec'),env)
    return env[fn.name]


@source_precision
def current_angular_support_proof(physical):
    physical.assert_graph();outer=physical.history.outer;c=outer.ctx
    normalization=fixed_normalization_and_angular_weight_source_proof(physical)
    if set(outer.weights)!=set(WEIGHTS):raise ValueError('All five original angular weights required')
    if outer.fifth is not physical.pulse.fifth or outer.future is not outer.fifth.fourth.energy.base:
        raise ValueError('Current angular coefficients/future owner changed')
    if not 0<endpoints(outer.mu)[0]<=endpoints(outer.mu)[1]<mp.mpf('.25'):
        raise ValueError('Current positive angular decay parameter required')
    flat=physical.pulse.flat
    if outer.normalization is not outer.flat.normalization or outer.flat.beta.__func__ is not FlatPulseDerivatives.beta:
        raise ValueError('Actual current angular beta does not use its declared flat reader')
    saved=json.loads((HERE/(PREFIX+'axial_pulse_field.json')).read_bytes())
    if endpoints(outer.normalization)!=endpoints(read_interval(c,saved['bump_normalization'])):
        raise ValueError('Current angular normalized beta reader differs from its pinned publication')
    if any(endpoints(value)!=endpoints(c.mpf(endpoints(outer.future.repair.weights[key]))) for key,value in outer.weights.items()):
        raise ValueError('Current angular full weight copy differs from its defining repair source')
    if endpoints(outer.flat.normalization)!=endpoints(c.mpf(endpoints(flat.normalization))):
        raise ValueError('Angular flat normalization differs from the same pinned source')
    for edge in EDGES:
        for center in (-3,-1):source_edge_integrals(outer,center,edge['exact_edge'])
    exact_endpoint_program()
    current=current_postpulse_trace_proof(physical)
    return dict(current_postpulse_defining_function_transfer=current,
        exact_normalization_and_five_weight_defining_source=normalization,
        weighted_integral_FTC_and_flat_beta_source_theorem=weighted_support_FTC_theorem(),
        five_current_weight_definitions=WEIGHTS,
        current_angular_endpoint_source_theorem=angular_endpoint_theorem(),
        exact_private_endpoint_reduction=dict(original_source_sha256=sha(PREFIX+'power_angular_C4.py'),
            beta='All exact source support endpoints have beta derivatives0..4 zero',
            integrals='Past/future exact empty/full support integrals by rational center +/-3/20 ordering',
            preserved='Current data(Z), translations A/B/D/E/F, nonlinear D/F terms, Post, flatten X/P, exponentials, _packet and original physical maps',
            public_guards_and_live_owners_unchanged=True),
        original_linear_physical_trace_transfer=physical_operator_trace_proof(physical),
        full_current_angular_local_stress_error_difference_bounds_remain_open=True,passed=True)


class CurrentAngularSupportInterfaces:
    @source_precision
    def __init__(self,pulse_support=None,require_checked=True):
        self.pulse_support=pulse_support if pulse_support is not None else CurrentPulseSupportInterfaces()
        if not self.pulse_support.acceptance_loaded:raise ValueError('Checked current four pulse support transfers required')
        self.physical=self.pulse_support.physical;self.outer=self.physical.history.outer;self.ctx=self.outer.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.proof=current_angular_support_proof(self.physical);self.hashes=dict(self.pulse_support.hashes)
        for name,digest in self.proof['exact_normalization_and_five_weight_defining_source']['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current angular support source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.endpoint_recipe=exact_endpoint_program();self.endpoint_cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current angular support admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def endpoint(self,edge,Z=(-1,1),log_tau='-1',theta=None):
        if edge not in EDGES:raise ValueError('Unknown angular source edge')
        self.physical.assert_graph();c=self.ctx;Z=c.mpf(Z);q=s.Rational(edge['exact_edge'])
        coordinate=c.mpf(str(q.p))/q.q;key=(edge['exact_edge'],tuple(endpoints(Z)))
        if key not in self.endpoint_cache:
            exact=self.endpoint_recipe(self.outer,Z,coordinate,edge['exact_edge'])
            diagnostic=self.outer.angular(Z,coordinate)
            self.endpoint_cache[key]=(exact,diagnostic)
        exact,diagnostic=self.endpoint_cache[key]
        view=copy.copy(self.physical);view.dispatch=EndpointDispatchView(self.physical,'outer_angular',exact)
        physical=BASE.evaluate(view,'outer_angular',Z,coordinate,log_tau=log_tau,theta=theta)
        return dict(edge=edge,Z=Z,exact_source_edge=edge['exact_edge'],source_coordinate_enclosure=coordinate,
            exact_source_endpoint_packet=exact,unchanged_native_coordinate_consistency_packet=diagnostic,
            original_physical_spatial4_time1_endpoint_log_bounds=physical_summary(physical),
            current_complete_future_and_absolute_pressure_retained=True,
            rounded_coordinate_consistency_is_not_a_function_identity=True,
            all_current_angular_local_stress_error_difference_bounds_certified=False,**dict.fromkeys(OPEN,False))

    def manifest(self):
        ledger=supported_beta_edges(self.physical)
        for row in ledger.values():row['current_full_field_interface_trace_admitted']=True
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_angular_support_source_proof=self.proof,
            current_internal_support_ledger=ledger,current_source_identified_adjacent_interface_count=14,
            current_source_identified_internal_support_count=8,
            all_current_angular_local_stress_error_difference_bounds_certified=False,
            quantitative_all_interface_admission=False,input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentAngularSupportInterfaces(require_checked=False)
    result=field.manifest();result['current_whole_Z_angular_support_endpoint_views']=[]
    for edge in EDGES:
        result['current_whole_Z_angular_support_endpoint_views'].append(field.endpoint(edge))
        print('Current angular exact support source generated: '+edge['exact_edge'],flush=True)
    result['input_hashes']=field.hashes
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
