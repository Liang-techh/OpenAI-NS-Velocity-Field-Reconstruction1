"""Current interface inventory and source-identified physical end/flatten trace."""
import ast
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_complete_physical_assembly import (
    CurrentCompletePhysicalAssembly,HERE,PREFIX,sha,pack,encode,endpoints,binding)
from lei_ren_part1_paper_compliant_current_selected_energy_source import function
from lei_ren_part1_paper_compliant_current_pulse_flatten_source import current_terminal_source_binding
from lei_ren_part1_paper_compliant_pulse_end_flatten_join import functional_join_identities
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UT
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_physical_interfaces.json';RECEIPT=PREFIX+'current_physical_interfaces_check.json'
GATES=('current_complete_affected_interface_ledger_available',
    'current_selected_end_flatten_functional_mixed4_join_certified',
    'current_selected_end_flatten_physical_spatial4_time1_trace_identified')
OPEN=('current_complete_quantitative_physical_interfaces_certified',
    'uniform_pulse_C4_chart_interface_certificate_available',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'independently_bounded_flat_remainder','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','temporal_recursion')
# Each coordinate is a source expression, not a selected floating endpoint.
SEAMS=(
    ('entrance_main','pulse_entrance','.02/mu','pulse_main','.02','mu*t=xi','1','1/mu'),
    ('main_exit','pulse_main','10','pulse_exit','10','same main callable','1/mu','1/mu'),
    ('exit_gap','pulse_exit','11','pulse_gap','11','same logR','1/mu','1/mu'),
    ('gap_coordinate','pulse_gap','12','pulse_gap_end','-1/mu','xi=13+mu*s','1/mu','1'),
    ('gap_end','pulse_gap_end','-4','pulse_end','-4','same end-scale s','1','1'),
    ('end_flatten','pulse_end','0','flatten','0','s=t=log(R/Rv)','1','1'),
    ('flatten_power','flatten','100','outer_power','0','y=(Lrel-4)*phase','1','Lrel-4'),
    ('power_angular','outer_power','1','outer_angular','-4','y=Lrel+s','Lrel-4','1'),
    ('angular_steep','outer_angular','0','steep_entry','0','same Rv*exp(100+Lrel)','1','1'),
    ('steep_entry_power','steep_entry','1','steep_power','0','t=Ts*phase','1','Ts'),
    ('steep_power_exit','steep_power','1','steep_exit','0','same logR at origin+1+Ts','Ts','1'),
    ('steep_exit_waiting','steep_exit','1','waiting','0','t=wait*phase','1','wait'),
    ('waiting_collar','waiting','1','heat_collar','0','same Rtail','wait','1'),
    ('collar_exterior','heat_collar','3','heat_exterior','3','same log(R/Rtail)','1','1'))


def interface_ledger(field):
    field.assert_graph();rows={}
    for name,left,lc,right,rc,relation,lj,rj in SEAMS:
        rows[name]=dict(left_chart=left,left_source_coordinate=lc,right_chart=right,
            right_source_coordinate=rc,exact_coordinate_relation=relation,
            dy_d_left_coordinate=lj,dy_d_right_coordinate=rj,
            left_owner=field.source_owners[left],right_owner=field.source_owners[right],
            native_grid_coordinate='ordinary y=logR,Z derivatives, not phase derivatives',
            required_velocity_pressure_mixed_order=4,required_primitive_histories=5,
            absolute_pressure_unit='Pstar^2',velocity_unit='pulse Pstar*exp(-bp*log(R/Rp)); postpulse Ev0',
            current_functional_mixed4_trace_admitted=name=='end_flatten',
            quantitative_all_interface_admission=False)
    return rows


def supported_beta_edges(field):
    """Interior compact-support boundaries are separate from adjacent charts."""
    return {chart+'_'+str(center)+'_'+str(side):dict(chart=chart,
        current_source_owner=field.source_owners[chart],ordinary_coordinate='s=logR offset',
        exact_source_edge=str(s.Rational(center)+side*s.Rational(3,20)),
        normalized_beta_source='flat.beta(s-center), ell=3/20',
        required_velocity_pressure_mixed_order=4,
        current_full_field_interface_trace_admitted=False,
        next_source_obligation='Same boundary histories and current selected/angular amplitudes; beta forcing jets vanish, but full current physical stress/remainder bounds must be transferred separately')
        for chart in ('pulse_end','outer_angular') for center in (-3,-1) for side in (-1,1)}


def native_end_flatten_trace_proof(field):
    """Bind native source recipes before deriving all60 endpoint identities."""
    terminal=current_terminal_source_binding(field);c=field.ctx
    if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in sigma_jets(c,0)):
        raise ArithmeticError('Original flatten inlet sigma jets are not flat')
    bindings={}
    specs={('compliant_pulse_radial_C4','pressure_moment'):{
        'log_decay':'sum(log_parts.values(),c.mpf(0))','decay':'self.factor(log_decay)',
        'kernel':"(decay_integral(c,self.prate,t) if endpoints(self.prate*t)[1]<mp.mpf('1e-20') else (1-decay)/self.prate)",
        'p':"IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)"},
        ('compliant_flatten_mixed_C4','flatten_mixed'):{
        'bp':"c.mpf('.5')+mu",'theta_rows':'[theta]','lograte[0]':'lograte[0]-bp',
        'rows':'{UZ:[zero for _ in range(5)],UT:theta_rows,UR:[zero for _ in range(5)],P:prows}'}}
    for (stem,method),values in specs.items():
        for target,value in values.items():
            binding(stem,method,target,value);bindings[stem+'.'+method+':'+target]=True
    # Raw pulse UT has ALREADY divided out the common theta radial factor.
    # Only U converts it to the flatten Ev0 unit; S must not be divided out
    # again. Pressure keeps its original absolute Pstar^2 unit throughout.
    binding('compliant_pulse_mixed_C4','transport_mixed','physical',
        "{'Uz_over_Pstar_without_common_theta_radial_factor':[u*binomial_product(Brows,-bp,k) for k in range(5)],'Utheta_over_Pstar_without_common_theta_radial_factor':[u*(-bp)**k for k in range(5)],'Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor':[u.truncate(4)*binomial_product(A,-mu,k) for k in range(5)],'P_over_Pstar_squared':pressure}")
    bindings['actual_raw_velocity_grid_excludes_common_radial_factor']=True
    pressure=function('compliant_pulse_radial_C4','pressure_moment')
    value=[n.value for n in ast.walk(pressure) if isinstance(n,ast.keyword) and n.arg=='P_y_over_Pstar_squared']
    if len(value)!=1 or ast.dump(value[0])!=ast.dump(ast.parse('u*u*(decay/2)',mode='eval').body):
        raise ValueError('Actual native pressure radial derivative source changed')
    end=function('compliant_axial_pulse_field','end')
    expression='self.pressure_moment(Z,t,inlet,u,dict(inverse_mu_term=-13*self.prate/self.mu,finite_offset=-self.prate*s))'
    if sum(ast.dump(n)==ast.dump(ast.parse(expression,mode='eval').body) for n in ast.walk(end) if isinstance(n,ast.Call))!=1:
        raise ValueError('Actual native terminal pressure exact log source changed')
    bindings['actual_native_pressure_y_and_end_log_parts_bound']=True
    binding('compliant_axial_pulse_field','__init__','self.prate','1+2*self.mu')
    binding('compliant_flatten_mixed_C4','flatten','mixed',
        'flatten_mixed(c,self.mu,rho,sj,theta,X,energy,pressure,self.Ev2)')
    binding('compliant_power_inlet_C4','incoming','raw',
        "self.datum.normalized_jets(Z,5)['normalized_pressure_coefficients']")
    binding('compliant_pulse_radial_C4','pressure_moment','rows',
        "self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),5)['normalized_pressure_coefficients']")
    for stem,method,key,expression in (
        ('compliant_power_inlet_C4','incoming','Mp','invq*invq*self.inlet_P'),
        ('compliant_pulse_radial_C4','data','Mp_over_Pstar_squared','[self.inlet_P*v for v in q_jets(c,Z,5)]')):
        values=[n.value for n in ast.walk(function(stem,method)) if isinstance(n,ast.keyword) and n.arg==key]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Actual common Pin publication changed: '+stem)
    bindings['actual_shared_Pin_P0_and_flat_Ev2_argument_bound']=True
    # Bind the actual production recurrence, including its binomial terms.
    fn=function('compliant_flatten_mixed_C4','flatten_mixed')
    for text in ('square+=theta_rows[j]*theta_rows[k-j]*math.comb(k,j)',
            'prows.append(square*(Ev2/2))',
            'th+=lograte[j]*theta_rows[k-j]*math.comb(k,j)'):
        target=ast.dump(ast.parse(text).body[0])
        if sum(ast.dump(n)==target for n in ast.walk(fn))!=1:
            raise ValueError('Native endpoint recurrence changed: '+text)
        bindings[text]=True
    z,U,mu,S,Pin=s.symbols('Z U mu exact_decay Pin',real=True)
    q=1+z*z;bp=s.Rational(1,2)+mu;prate=1+2*mu
    datum=s.Function('same_original_P0')(z)
    identities={}
    for k in range(5):
        left_theta=(U/q)*(-bp)**k/U
        right_theta=(-bp)**k/q
        if k==0:
            left_p=right_p=Pin/q**2+U**2/q**2*(1-S)/(2*prate)+datum
        else:
            left_p=U**2*S/(2*q*q)*(-prate)**(k-1)
            square=sum(s.binomial(k-1,j)*(-bp)**j/q*(-bp)**(k-1-j)/q for j in range(k))
            right_p=U**2*S*square/2
        for n in range(5-k):
            for label,left,right in (('Utheta',left_theta,right_theta),('P',left_p,right_p),('Ur',0,0),('Uz',0,0)):
                if s.simplify(s.diff(left-right,z,n))!=0:
                    raise ArithmeticError('Current native endpoint differs: '+label)
                identities['%s_y%d_Z%d'%(label,k,n)]=True
    lp,mu0,uu=s.symbols('logPstar mu U',positive=True)
    end_unit=lp-13/(2*mu0)-13
    flat_unit=lp+s.log(uu)-13/(2*mu0)-13
    unit_identities=dict(raw_theta_factor_excluded_and_U_rebase=s.expand(end_unit+s.log(uu)-flat_unit)==0,
        native_terminal_pressure_decay_is_same_S=s.expand(-(1+2*mu0)*13/mu0+13/mu0+26)==0,
        exact_flat_Ev0_squared_is_U_squared_times_same_S=s.expand(2*(s.log(uu)-13/(2*mu0)-13)-(2*s.log(uu)-13/mu0-26))==0)
    unit_identities['native_prate_equals_twice_flatten_bp']=s.expand(1+2*mu0-2*(s.Rational(1,2)+mu0))==0
    if not all(unit_identities.values()):
        raise ArithmeticError('Original endpoint positive unit differs')
    return dict(current_terminal_source_binding=terminal,actual_native_recurrence_bindings=bindings,
        mixed4_identities=identities,identity_count=len(identities),
        exact_positive_velocity_and_pressure_log_source_identities=unit_identities,
        raw_native_UT_grid_excludes_the_common_radial_exponential=True,
        original_sigma_inlet_jets_exactly_zero=True,
        original_positive_Ev0_equals_U_times_pulse_endpoint_unit=True,
        original_pressure_decay_source='S=exp(-13/mu-26), Ev0^2/Pstar^2=U^2*S',
        exact_pressure_source_bound_before_its_numeric_cap_enclosure=True,
        exact_current_energy_half_and_empty_support_histories=field.history.proof,
        shared_current_cartesian_map=field.operator_bindings,
        physical_trace_argument='Equal ordinary y,Z mixed4 functions at the same Rv, same positive units and datum; identical original Cartesian/moving-basis and fixed-x time operators preserve these traces for every positive tau',
        source_caps_are_enclosures_not_defining_endpoint_values=True,passed=True)


class CurrentPhysicalInterfaces:
    @source_precision
    def __init__(self,physical=None,require_checked=True):
        self.physical=physical if physical is not None else CurrentCompletePhysicalAssembly()
        if not self.physical.physical_acceptance_loaded:raise ValueError('Checked current33 physical source required')
        self.physical.assert_graph();self.ctx=self.physical.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.ledger=interface_ledger(self.physical);self.proof=native_end_flatten_trace_proof(self.physical)
        self.canonical=functional_join_identities();self.hashes=dict(self.physical.hashes)
        for name,digest in self.canonical['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:
                raise ValueError('Current physical/canonical interface source conflict: '+name)
            self.hashes[name]=digest
        for stem in ('current_pulse_flatten_source','pulse_end_flatten_join','flat_pulse_derivatives'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current endpoint interface admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    @source_precision
    def endpoint(self,Z,log_tau='-1',theta=None):
        p=self.physical;p.assert_graph()
        left=p.dispatch.evaluate('pulse_end',Z,0)['source_packet'];right=p.history.flatten.flatten(Z,0)
        c=self.ctx;U=p.flatten.U
        lg=left['physical_mixed_derivatives_total_order_le4'];rg=right['physical_mixed_derivatives_total_order_le4']
        differences={label:{key:lv/(U if label==UT else 1)-rg[label][key] for key,lv in rows.items()}
            for label,rows in lg.items()}
        return dict(Z=c.mpf(Z),native_fixed_Ev0_velocity_and_Pstar_squared_pressure_difference_enclosures=differences,
            left_physical_trace=p.evaluate('pulse_end',Z,0,log_tau=log_tau,theta=theta),
            right_physical_trace=p.evaluate('flatten',Z,0,log_tau=log_tau,theta=theta),
            exact_function_identity_precedes_interval_consistency=True,
            source_identity_scope='one end/flatten physical spatial4/fixed-x time1 trace, not all interfaces',
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_affected_interface_ledger=self.ledger,
            current_supported_beta_edge_inventory=supported_beta_edges(self.physical),
            current_end_flatten_native_source_proof=self.proof,
            recomputed_arbitrary_shared_function_endpoint_theorem=self.canonical,
            exact_end_flatten_source_boundary='R=Rv, s=t=0; d_s=d_t=d_logR',
            other_current_mixed4_interfaces_remaining=[name for name in self.ledger if name!='end_flatten'],
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPhysicalInterfaces(require_checked=False)
    result=field.manifest();result['current_end_flatten_views']={
        'whole_Z':field.endpoint([-1,1]),'fresh':field.endpoint('.419','-2.31','.41')}
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current14 interface ledger and current end/flatten physical traces generated',flush=True)
    return result


if __name__=='__main__':run()
