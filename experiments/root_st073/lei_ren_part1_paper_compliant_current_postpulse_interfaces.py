"""Current complete-energy instantiation of eight postpulse mixed4 seams.

Function identities precede interval consistency. This does not admit all
interfaces, a full tensor/NS solution, resolved points or time recursion.
"""
import ast
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_physical_interfaces import (
    CurrentPhysicalInterfaces,SEAMS,OPEN,HERE,PREFIX,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_postpulse_energy_history import energy_transport_source_proof
from lei_ren_part1_paper_compliant_current_full_exterior_stress import current_full_history_transfer
from lei_ren_part1_paper_compliant_power_angular_C4_check import functional_source_identities as power_theorem
from lei_ren_part1_paper_compliant_steep_waiting_C4_check import functional_source_identities as steep_theorem
from lei_ren_part1_paper_compliant_collar_Gamma_C4_check import functional_source_identities as heat_theorem
from lei_ren_part1_paper_compliant_current_selected_energy_source import function,binding
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import phi_jets
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_postpulse_interfaces.json';RECEIPT=PREFIX+'current_postpulse_interfaces_check.json'
GATES=('current_complete_postpulse_functional_mixed4_interfaces_certified',
       'current_complete_postpulse_physical_spatial4_time1_traces_identified')
POSTSEAMS=SEAMS[6:]


def statement(stem,method,text):
    wanted=ast.dump(ast.parse(text).body[0])
    if sum(ast.dump(n)==wanted for n in ast.walk(function(stem,method)))!=1:
        raise ValueError('Postpulse production statement changed: '+stem+'.'+method+' '+text)
    return text


def ordinary_coordinate_proof(physical):
    """Bind actual phase recipes; derivatives are already ordinary logR jets."""
    statements={}
    specs=(('compliant_power_angular_C4','power','y=(self.Lrel-4)*phase'),
        ('compliant_steep_waiting_C4','steep_power','t=self.Ts*phase'),
        ('compliant_steep_waiting_C4','waiting','t=self.wait*phase'),
        ('compliant_global_physical_assembly','radius',
         "outer=provider;offset=100+(outer.Lrel-4)*v if chart=='outer_power' else 100+outer.Lrel+v"),
        ('compliant_global_physical_assembly','radius',
         "offset=(origin+v if chart=='steep_entry' else origin+1+steep.Ts*v if chart=='steep_power' else origin+1+steep.Ts+v if chart=='steep_exit' else origin+2+steep.Ts+steep.wait*v)"),
        ('compliant_global_physical_assembly','radius','offset=100+steep.outer.Lrel+2+steep.Ts+steep.wait+v'))
    for stem,method,text in specs:
        # A semicolon is two separate production statements.
        for item in text.split(';'):
            statements[stem+'.'+method+':'+item]=statement(stem,method,item)
    for target,value in (('self.wait','box(f.angular.waiting)'),('self.Ts','box(f.params.Ts)')):
        class_assignment('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',target,value)
        statements['current_steep_constructor:'+target]=value
    c=physical.ctx;h=physical.history
    copies={'Lrel':(h.outer.Lrel,h.selected.future.Lrel),
            'Ts':(h.steep.Ts,h.selected.future.params.Ts),
            'wait':(h.steep.wait,h.selected.future.angular.waiting)}
    for name,(value,source) in copies.items():
        if endpoints(value)!=endpoints(c.mpf(endpoints(source))):
            raise ValueError('Current phase length differs from its defining source: '+name)
    if any(endpoints(v)[0]<=0 for v in (h.outer.Lrel-4,h.steep.Ts,h.steep.wait)):
        raise ValueError('Ordinary phase conversion requires positive source lengths')
    p,o=s.symbols('phase origin',real=True);J=s.symbols('J',positive=True)
    dy,dz=s.symbols('ordinary_y_derivative ordinary_Z_derivative',commutative=True)
    if s.diff(o+J*p,p)!=J:raise ArithmeticError('Source phase Jacobian differs')
    chains={}
    for order in range(5):
        for n in range(5-order):
            # Constant source lengths give D_phase=J D_y. These smooth
            # derivative operators commute because J is independent of Z.
            if s.simplify((J*dy)**order*dz**n/J**order-dy**order*dz**n)!=0:
                raise ArithmeticError('Phase to ordinary logR chain differs')
            chains['phase_y%d_Z%d'%(order,n)]=True
    rv,L,T,W=s.symbols('logRv L Ts wait',real=True)
    radius=( (rv+100,rv+100),(rv+100+L-4,rv+100+L-4),
        (rv+100+L,rv+100+L),(rv+100+L+1,rv+100+L+1),
        (rv+100+L+1+T,rv+100+L+1+T),
        (rv+100+L+2+T,rv+100+L+2+T),
        (rv+100+L+2+T+W,rv+100+L+2+T+W),
        (rv+100+L+2+T+W+3,rv+100+L+2+T+W+3))
    equal={row[0]:s.expand(a-b)==0 for row,(a,b) in zip(POSTSEAMS,radius)}
    if not all(equal.values()):raise ArithmeticError('Postpulse exact radius differs')
    return dict(actual_coordinate_statements=statements,same_current_defining_lengths=True,
        ordinary_chain_identities=chains,phase_derivatives_to_y_factor='(dy/dphase)^(-order)',
        production_native_mixed_grids_already_use_y_not_phase=True,
        exact_common_radius_identities=equal,source_lengths_independent_of_Z=True,passed=True)


def flatten_power_current_theorem():
    """Actual common flatten exit is substituted before axial differentiation."""
    for text in ('theta=F/q*c.exp(-self.bp*t)',
        'theta=IntervalTaylor.constant(c,c.exp(-100*self.bp)/2,5)',
        'energy=(ev-Eint/2)*c.exp(2*self.mu*t)/(F*F)'):
        statement('compliant_flatten_mixed_C4','flatten',text)
    binding('compliant_current_power_angular_source','data','f','self.flatten.flatten(Z,100)')
    fn=function('compliant_current_power_angular_source','data')
    for key,text in (('flatten_exit_X',"f['angular_Taylor']"),
                     ('flatten_exit_pressure',"f['pressure']"),('flatten_exit_energy',"f['energy_Taylor']")):
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.keyword) and n.arg==key]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(text,mode='eval').body):
            raise ValueError('Actual current flatten boundary publication changed: '+key)
    binding('compliant_power_angular_C4','_packet','mixed',
        'flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    z,mu,L,Ev2=s.symbols('Z mu L Ev2',real=True);q=1+z*z;bp=s.Rational(1,2)+mu
    EF=s.Function('current_exact_flatten_piece')(z);post=s.Function('current_complete_post_future')(z)
    change=s.Function('current_signed_angular_change')(z)
    Xf=s.Function('actual_flatten_exit_X')(z);Pf=s.Function('actual_flatten_exit_absolute_P')(z)
    D=lambda r,t:(1-s.exp(-r*t))/r
    total=EF+q*q*s.exp(-200*mu)/4*(D(2*mu,L)+s.exp(-2*mu*L)*(post+change))
    flat=[q/2/q*s.exp(-100*bp),Xf,(total/2-EF/2)*s.exp(200*mu)/(q*q/4),Pf]
    env={'one':s.Integer(1),'self.bp':bp,'self.rate':1-mu,'self.mu':mu,'y':s.Integer(0),
        "data['flatten_exit_X']":Xf,"data['flatten_exit_pressure']['P_over_Pstar_squared']":Pf,
        'deltaP':s.Integer(0),"data['post']":post,"data['full_angular_change']":change,
        'D':L,'decay_integral(c, 2 * self.mu, D)':D(2*mu,L)}
    power=[assignment('compliant_power_angular_C4','power',v,env) for v in ('theta','X','energy','pressure')]
    power[2]+=assignment('compliant_power_angular_C4','power','energy',env,augmented=True)
    proofs={};a=flat;b=power
    for order in range(5):
        for n in range(5-order):
            for label,left,right in zip(('theta','X','energy','pressure'),a,b):
                if s.simplify(s.diff(left-right,z,n))!=0:
                    raise ArithmeticError('Current flatten/power functional mixed4 trace differs')
                proofs['flatten_power_'+label+'_y%d_Z%d'%(order,n)]=True
        a=[-bp*a[0],(1 if order==0 else 0)-(1-mu)*a[1],
           (-s.Rational(1,2) if order==0 else 0)+2*mu*a[2],Ev2*flat[0]**2/2*(-2*bp)**order]
        b=[-bp*b[0],(1 if order==0 else 0)-(1-mu)*b[1],
           (-s.Rational(1,2) if order==0 else 0)+2*mu*b[2],Ev2*power[0]**2/2*(-2*bp)**order]
    return proofs


def physical_operator_trace_proof(physical):
    """Source-equal rows pass through the same exact linear physical map.

    This is a function/operator proof before log_row's triangle enclosures.
    It never treats two overlapping physical upper bounds as equal fields.
    """
    bindings={}
    value='BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)'
    class_assignment('current_complete_physical_assembly','CurrentCompletePhysicalAssembly','evaluate','packet',value)
    bindings['current_complete_evaluate_call']=value
    specs=(('compliant_current_pulse_physical_assembly','normalized_sources',
            'grids,logs,amplitudes=BASE.normalized_sources(self,chart,Z,value,packet,provider,logR)'),
        ('compliant_current_pulse_physical_assembly','normalized_sources',
            'if chart not in PULSE_CHARTS:return grids,logs,amplitudes'),
        ('compliant_global_physical_assembly','normalized_sources',
            "grids={label:ordinary_terms(raw[label]) for label in (UZ,UT,UR,P)}"),
        ('compliant_global_physical_assembly','evaluate',
            'grids,logs,amplitudes=self.normalized_sources(chart,Z,coordinate,packet,provider,logR)'),
        ('compliant_global_physical_assembly','evaluate',
            'parts=cartesian_source_row(c,grids,component,i,j,b,Z,self.delta,cosine,sine,amplitudes)'),
        ('compliant_global_physical_assembly','evaluate',
            'terms,gamma=time_source_row(c,grids,label,Z,self.delta,amplitudes)'),
        ('compliant_global_physical_assembly','cartesian_source_row',
            'coefficient=angular_factor*interval_expression(c,expression,Z,delta,beta[label])'),
        ('compliant_global_physical_assembly','cartesian_source_row',
            'by_label.setdefault(label,[]).append((shift(powers,powers_to_add),source_coefficient*coefficient))'),
        ('compliant_global_physical_assembly','time_source_row',
            'for index,coefficient in (((0,0),-beta/(2*L)),((0,1),(1-delta)*Z/(2*L)),((1,0),1/L)):\n    for powers,value in grids[label][index]:terms.append((shift(powers,amplitudes[label]),value*coefficient))'))
    for stem,method,text in specs:bindings[stem+'.'+method+':'+text]=statement(stem,method,text)
    if not all(physical.operator_bindings.values()):raise ValueError('Original physical operators changed')
    a,b,C=s.symbols('left_native_row right_native_row common_operator_coefficient',real=True)
    if s.expand(a*C-b*C-C*(a-b))!=0:raise ArithmeticError('Exact linear physical difference transfer failed')
    return dict(actual_original_physical_operator_AST=bindings,
        all_nine_postpulse_charts_bypass_pulse_rebase=True,
        same_positive_Ev0_and_absolute_Pstar_squared_units=physical.proof['original_post_fixed_unit_binding'],
        actual_spatial_and_fixed_x_time_rows_linear_in_native_mixed_rows=True,
        coefficients_depend_on_shared_Z_delta_angle_radius_units_not_chart_histories=True,
        source_equal_mixed_rows_give_zero_exact_physical_operator_difference=True,
        physical_tau_is_a_shared_symbol_not_the_selected_diagnostic_time=True,
        domain='finite physical traces: R>0, abs(Z)<1, tau>0; Z=+/-1 covered as limiting source bounds',
        triangle_log_bounds_are_not_defining_field_values=True,passed=True)


@source_precision
def current_postpulse_trace_proof(physical):
    physical.assert_graph();h=physical.history;c=physical.ctx
    if not physical.physical_acceptance_loaded:raise ValueError('Checked current complete source required')
    coordinate=ordinary_coordinate_proof(physical)
    if h.heat.Ev2 is not h.flatten.Ev2:
        raise ValueError('Heat and flatten absolute pressure use different Ev0 sources')
    energy=energy_transport_source_proof(h);full=current_full_history_transfer(h)
    canonical={'power':power_theorem(),'steep':steep_theorem(),'heat':heat_theorem()}
    if not energy['passed'] or not full['passed'] or not all(all(rows.values()) for rows in canonical.values()):
        raise ValueError('Current source / arbitrary-function theorem unavailable')
    flats={}
    for x in (0,1):
        jets=sigma_jets(c,x)
        if endpoints(jets[0])!=(mp.mpf(x),mp.mpf(x)) or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in jets[1:]):
            raise ValueError('Original sigma endpoint mixed4 flatness lost')
        flats['sigma_'+str(x)]=True
    if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in phi_jets(c,3)):
        raise ValueError('Original collar phi mixed4 terminal flatness lost')
    flats['phi_3']=True
    for endpoint in (-4,0):
        for center in (-3,-1):
            if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in h.outer.flat.beta(endpoint-center)):
                raise ValueError('Angular endpoint original beta support is nonempty')
    flats['angular_endpoint_beta_jets_zero']=True
    identities=flatten_power_current_theorem()
    mappings=(('power_angular','power','new_provider_endpoint'),
        ('angular_steep','steep','angular_entry'),('steep_entry_power','steep','entry_power'),
        ('steep_power_exit','steep','power_exit'),('steep_exit_waiting','steep','exit_waiting'),
        ('waiting_collar','heat','waiting_collar'),('collar_exterior','heat','collar_Gamma'))
    for seam,group,prefix in mappings:
        for label in ('theta','X','energy','pressure'):
            for order in range(5):
                for n in range(5-order):
                    key=prefix+'_'+label+'_y%d_Z%d'%(order,n)
                    if not canonical[group].get(key):raise ValueError('Arbitrary-function seam identity missing: '+key)
                    identities[seam+'_'+label+'_y%d_Z%d'%(order,n)]=True
    if len(identities)!=480:raise ValueError('Eight current primitive seam grids incomplete')
    return dict(current_owner_graph=h.assert_graph(),current_complete_energy_substitution=energy,
        current_original_angular_pressure_function_transfer=full,
        recomputed_arbitrary_function_theorems=canonical,current_instantiated_seam_mixed4_identities=identities,
        original_flat_endpoint_jet_evidence=flats,ordinary_coordinates=coordinate,
        original_fixed_Ev0_Pstar_squared_units=physical.proof['original_post_fixed_unit_binding'],
        original_physical_operators=physical.operator_bindings,
        source_bound_exact_linear_physical_trace_transfer=physical_operator_trace_proof(physical),
        exact_source_substitutions=dict(Post='current complete post-angular future, with all transition and full Gamma pieces',
            H='current complete E0(Z)/(1-epsilon)^2',EF='same exact full flatten energy piece',
            X='actual forward angular history, unchanged by energy replacement',
            P='actual original P0 plus forward pressure history in Pstar^2',
            meridional='selected terminal zero histories propagated by their homogeneous equations and FTC'),
        arbitrary_functions_instantiated_before_axial_differentiation=True,
        all_current_postpulse_primitive_mixed4_traces_identified=True,
        physical_trace_argument='Same source-identified mixed4 velocity/absolute pressure at the same exact radii; source-bound exact linear Cartesian/moving-basis/fixed-x time maps preserve the traces for every positive tau and abs(Z)<1. Limiting bounds cover Z=+/-1. This is structural function identification, not global quantitative or resolved-point certification.',
        interval_overlap_is_not_the_function_identity_proof=True,passed=True)


class CurrentPostpulseInterfaces:
    @source_precision
    def __init__(self,interfaces=None,require_checked=True):
        self.interfaces=interfaces if interfaces is not None else CurrentPhysicalInterfaces()
        if not self.interfaces.acceptance_loaded:raise ValueError('Checked current interface prerequisite required')
        self.physical=self.interfaces.physical;self.ctx=self.physical.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.proof=current_postpulse_trace_proof(self.physical);self.hashes=dict(self.interfaces.hashes)
        for stem in ('power_angular_C4_check','steep_waiting_C4_check','collar_Gamma_C4_check',
                     'current_postpulse_energy_history','current_full_exterior_stress'):
            name=PREFIX+stem+'.py'
            if name in self.hashes and self.hashes[name]!=sha(name):raise ValueError('Current postpulse dependency conflict: '+name)
            self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current eight-interface admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def endpoint(self,name,Z,log_tau='-1',theta=None):
        rows={row[0]:row for row in POSTSEAMS}
        if name not in rows:raise ValueError('Unknown current postpulse seam: '+name)
        _,l,lc,r,rc,*_=rows[name];p=self.physical;p.assert_graph()
        a=p.dispatch.evaluate(l,Z,lc)['source_packet'];b=p.dispatch.evaluate(r,Z,rc)['source_packet']
        differences={label:{key:value-b['physical_mixed_derivatives_total_order_le4'][label][key] for key,value in grid.items()}
            for label,grid in a['physical_mixed_derivatives_total_order_le4'].items()}
        for label in ('angular','energy'):
            left=a[label+'_y_derivative_Taylor'];right=b[label+'_y_derivative_Taylor']
            differences[label]={'y%d_Z%d'%(k,n):(left[k][n]-right[k][n])*math.factorial(n)
                for k in range(5) for n in range(5-k)}
        return dict(seam=name,Z=self.ctx.mpf(Z),same_current_native_mixed4_difference_enclosures=differences,
            left_physical_trace=p.evaluate(l,Z,lc,log_tau=log_tau,theta=theta),
            right_physical_trace=p.evaluate(r,Z,rc,log_tau=log_tau,theta=theta),
            source_function_identity_precedes_consistency=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        ledger={name:dict(row) for name,row in self.interfaces.ledger.items()}
        for name in [row[0] for row in POSTSEAMS]:ledger[name]['current_functional_mixed4_trace_admitted']=True
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_postpulse_source_trace_proof=self.proof,
            current_source_identified_interface_ledger=ledger,
            source_identified_adjacent_interface_count=9,
            remaining_current_adjacent_trace_interfaces=[name for name,row in ledger.items() if not row['current_functional_mixed4_trace_admitted']],
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPostpulseInterfaces(require_checked=False)
    result=field.manifest();result['current_whole_Z_postpulse_trace_views']={name:field.endpoint(name,[-1,1]) for name,*_ in POSTSEAMS}
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current eight postpulse source-identified mixed4 and physical traces generated',flush=True)
    return result


if __name__=='__main__':run()
