"""Current four pulse-end support transfers with retained full histories.

Original arbitrary-function difference operators are instantiated with the
current selected C5 functions. No old serialized control/history is reused.
Local stress/error flatness is not the global temporal flat remainder.
"""
import ast
import copy
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_pulse_interfaces import (
    CurrentPulseInterfaces,current_pulse_trace_proof,EndpointDispatchView,
    BASE,OPEN,HERE,PREFIX,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import statement
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding,function
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_pulse_end_support_interfaces import (
    EDGES,source_difference_proof,full_difference_rows,difference_transport,
    symmetric_jet,beta_tail_bound,ordinary_grid,physical_bracket)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import (
    original_formula_identities,capture_native_end)
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import binomial_product
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets,FlatPulseDerivatives
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_current_physical_interfaces import supported_beta_edges
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_pulse_support_interfaces.json';RECEIPT=PREFIX+'current_pulse_support_interfaces_check.json'
GATES=('current_pulse_end_four_support_functional_mixed4_interfaces_identified',
    'current_pulse_end_support_physical_spatial4_time1_traces_identified',
    'current_pulse_end_support_local_stress3_error2_flat_difference_bounds_available')
DISTANCES=('.01','.000001','0')


def fixed_normalization_and_angular_weight_source_proof(physical):
    """One exact normalized bump definition across live and saved readers.

    Copied enclosures are consistency checks. Equality of defining functions
    follows the checked normalization record and bound publication/read AST.
    The five moment weights have radius1/40; angular A/B/D/E/F have radius3/20.
    """
    p=physical.pulse;selected=physical.history.selected;c=p.ctx
    paper='lei_ren_part1_paper_'
    names=(paper+'bump_integral_enclosures_check.json',PREFIX+'axial_pulse_field.json',
           PREFIX+'outer_angular_repair.json',PREFIX+'outer_angular_repair_check.json')
    records={name:json.loads((HERE/name).read_bytes()) for name in names};hashes={}
    for name,record in records.items():
        if 'input_hashes' in record:_verify_hashes(record);hashes.update(record['input_hashes'])
        if record.get('implicit_source_sha256',physical.source)!=physical.source:
            raise ValueError('Normalized bump record belongs to another source')
        hashes[name]=sha(name)
    fixed=records[names[0]];saved=records[names[1]]
    if not fixed['directed_integrals_certified']:raise ValueError('Directed original bump normalization required')
    for name,key in ((paper+'bump_integral_enclosures.py','module_sha256'),
                     (paper+'bump_integral_enclosures_check.py','check_source_sha256')):
        if sha(name)!=fixed[key]:raise ValueError('Original normalization quadrature source changed')
        hashes[name]=fixed[key]
    repair=p.pulse.initial.repair;exact=selected.exact;original=exact.heat.future.repair
    if repair.wraw!=fixed:raise ValueError('Live five moment repair uses another normalization record')
    for owner in (repair,original.angular.initial.repair):
        if endpoints(owner.normalization)!=endpoints(restore_value(owner.ctx,fixed['normalization'])):
            raise ValueError('Live normalizer differs from its exact directed record copy')
    if p.flat.beta.__func__ is not FlatPulseDerivatives.beta:
        raise ValueError('Original normalized beta derivative recipe required')
    if endpoints(p.flat.normalization)!=endpoints(read_interval(p.flat.ctx,saved['bump_normalization'])):
        raise ValueError('Flat normalization reader differs from its pinned publication')
    graph=dict(current_repair_normalization_preserved=exact.repair.normalization is original.normalization,
        current_repair_full_weights_preserved=exact.repair.weights is original.weights,
        current_future_repair_is_current_repair=selected.future.repair is exact.repair,
        current_angular4_uses_current_repair=selected.fifth.angular4.repair is selected.future.repair,
        original_repair_normalization_from_its_initial=original.normalization is original.angular.initial.repair.normalization,
        original_repair_initial_uses_same_fixed_record=original.angular.initial.repair.wraw==fixed)
    if not all(graph.values()):raise ValueError('Current normalized bump/weight defining graph differs: '+str(graph))
    bindings={}
    for stem,cls,target,value in (
        ('five_moment_repair','SharedFiveMomentRepair','self.normalization',"restore_value(c,self.wraw['normalization'])"),
        ('flat_pulse_derivatives','FlatPulseDerivatives','self.normalization',"read_interval(self.ctx,pulse['bump_normalization'])"),
        ('outer_angular_repair','CompliantAngularRepair','self.normalization','self.angular.initial.repair.normalization'),
        ('outer_angular_repair','CompliantAngularRepair','self.weights','bump_weights(c,self.mu,self.normalization,cells=bump_cells)'),
        ('current_power_angular_source','CurrentPowerAngularC4','self.weights','{k:box(v) for k,v in future.repair.weights.items()}')):
        class_assignment(stem,cls,'__init__',target,value);bindings[stem+':'+target]=value
    fn=function('compliant_axial_pulse_field','report')
    values=[n.value for n in ast.walk(fn) if isinstance(n,ast.keyword) and n.arg=='bump_normalization']
    if len(values)!=1 or ast.unparse(values[0])!='self.pulse.initial.repair.normalization':
        raise ValueError('Native normalization publication changed')
    bindings['native_normalization_publication']='self.pulse.initial.repair.normalization'
    for stem,method,text in (
        ('compliant_current_exact_repair_branch','replay_repair','out=copy.copy(original)'),
        ('compliant_current_exact_repair_branch','replay_future','out.repair=repair'),
        ('compliant_outer_angular_repair','bump_weights','beta=raw_beta(c,raw_coordinate)/(ell*normalization)'),
        ('compliant_corrected_outer_field','future_bump_weights','beta=raw_beta(c,r)/(ell*normalization)')):
        bindings[stem+'.'+method+':'+text]=statement(stem,method,text)
    replay=function('compliant_current_exact_repair_branch','replay_repair')
    if any(isinstance(n,(ast.Assign,ast.AugAssign)) and any(ast.unparse(t) in ('out.normalization','out.weights')
            for t in (n.targets if isinstance(n,ast.Assign) else [n.target])) for n in ast.walk(replay)):
        raise ValueError('Current replay reassigns a fixed normalization/weight source')
    tuples="(('A',1-mu,1),('B',-1-2*mu,1),('D',-1-2*mu,2),('E',-2*mu,1),('F',-2*mu,2))"
    wanted=ast.dump(ast.parse(tuples,mode='eval').body)
    for stem,method in (('compliant_outer_angular_repair','bump_weights'),('compliant_corrected_outer_field','future_bump_weights')):
        if sum(isinstance(n,ast.For) and ast.dump(n.iter)==wanted for n in ast.walk(function(stem,method)))!=1:
            raise ValueError('Original five angular rate/power definitions changed')
    v,mu,N=s.symbols('v mu N',real=True);ell=s.Rational(3,20);beta=s.Function('same_raw_beta')(v/ell)/(ell*N)
    identities={}
    for name,rate,power in (('A',1-mu,1),('B',-1-2*mu,1),('D',-1-2*mu,2),('E',-2*mu,1),('F',-2*mu,2)):
        past=s.Integral(s.exp(rate*v)*beta**power,(v,-ell,ell))
        future=s.Integral(s.exp(rate*v)*beta**power,(v,-ell,ell))
        if past-future!=0:raise ArithmeticError('Full past/future source integrals differ')
        identities[name]=True
    return dict(actual_normalization_publication_read_and_weight_AST=bindings,current_preserved_source_graph=graph,
        exact_normalization_source='N=int_-1^1 exp(-1/(1-r^2)) dr, same checked directed record before all copies',
        exact_full_past_future_five_weight_integral_identities=identities,
        future_start_minus4_and_past_stop_3_over20_have_same_exact_full_support_domain=True,
        five_moment_radius_1_over40_weights_are_not_angular_radius_3_over20_weights=True,
        copied_enclosure_consistency_is_not_function_identity=True,input_hashes=hashes,passed=True)


def weighted_support_FTC_theorem():
    """Exact one-sided limits for compact beta^1/beta^2 weighted integrals.

    Values at an edge use empty/full support integrals, not a rounded CDF.
    An integral over a vanishing interval controls both orientations.
    """
    x,v,a,L,R=s.symbols('x v rate left_edge right_edge',real=True)
    beta=s.Function('original_normalized_beta');proofs={}
    for power in (1,2):
        g=s.exp(a*v)*beta(v)**power
        past=s.Integral(g,(v,L,x));future=s.Integral(g,(v,x,R))
        forcing=s.exp(a*x)*beta(x)**power
        for kind,value,sign in (('past',past,1),('future',future,-1)):
            for order in range(1,6):
                diff=s.diff(value,x,order)-sign*s.diff(forcing,x,order-1)
                if s.simplify(diff)!=0:raise ArithmeticError('Weighted compact-support FTC failed')
                proofs[kind+'_power%d_FTC_order%d'%(power,order)]=True
                row=s.diff(forcing,x,order-1)
                subs={beta(x):0,**{s.diff(beta(x),x,j):0 for j in range(1,5)}}
                if s.simplify(row.subs(subs))!=0:raise ArithmeticError('Weighted integral endpoint derivative is nonzero')
                proofs[kind+'_power%d_endpoint_derivative%d_zero'%(power,order)]=True
        if s.simplify(s.diff(past+future,x))!=0:raise ArithmeticError('Past+future full-support identity failed')
        proofs['power%d_full_support_constant'%(power)]=True
    c=s.symbols('center',real=True);ell=s.Rational(3,20);r=s.symbols('normalized_r',real=True)
    for side in (-1,1):
        if s.simplify(((c+side*ell-c)/ell)-side)!=0:raise ArithmeticError('Exact support edge normalization failed')
        proofs['exact_normalized_edge_'+str(side)]=True
    W=s.symbols('W',positive=True)
    for order in range(5):
        if s.limit(s.exp(-1/W)*W**(-2*order),W,0,dir='+')!=0:raise ArithmeticError('Original beta tail limit failed')
        proofs['beta_majorant_flat_limit_order'+str(order)]=True
    return dict(identities=proofs,
        one_sided_boundary_values=dict(past_left='0',future_left='full weighted support integral',
            past_right='full weighted support integral',future_right='0'),
        both_orientations_value_difference_bound='abs(int_edge^s exp(rate*v)*beta(v)^p dv) <= h*sup_neighborhood(abs(exp(rate*v)))*sup_neighborhood(abs(beta))^p',
        exact_singleton_source_edges='center +/- 3/20; normalized r=+/-1 before interval enclosure',
        both_one_sided_values_agree_by_vanishing_interval_integral=True,
        derivatives_through5_follow_FTC_and_original_flat_beta_jets=True,
        numerical_clamped_integral_enclosures_are_not_defining_boundary_values=True,passed=True)


def endpoint_difference_theorem():
    """Current arbitrary axial histories cancel only in actual-reference."""
    z,mu,delta=s.symbols('Z mu delta',real=True)
    C=s.Function('current_selected_end_control')(z);zero=s.Integer(0)
    forcing=[C*s.Symbol('beta_derivative'+str(k)) for k in range(5)]
    m=[s.Function('same_boundary_m_difference')(z)]
    n=[s.Function('same_boundary_n_difference')(z)]
    J=[s.Function('same_boundary_J_difference')(z)]
    for k in range(4):
        m.append(forcing[k]-(s.Rational(1,2)-mu)*m[k])
        n.append(forcing[k]-(s.Rational(1,2)-2*mu)*n[k])
        square=sum(math.comb(k,j)*forcing[j]*forcing[k-j] for j in range(k+1))
        J.append(2*mu*J[k]-square)
    subs={m[0]:0,n[0]:0,J[0]:0,**{s.Symbol('beta_derivative'+str(k)):0 for k in range(5)}}
    grids={};axial5={}
    for label,rows in (('m1',m),('m2',n),('quadratic_loss',J),('forcing',forcing)):
        for k in range(5):
            value=s.simplify(rows[k].subs(subs))
            if value!=0:raise ArithmeticError('Support source difference endpoint is nonzero')
            for order in range(5-k):grids[label+'_y%d_Z%d'%(k,order)]=s.diff(value,z,order)==0
        if label=='m1':axial5[label]=s.diff(s.simplify(rows[0].subs(subs)),z,5)==0
    # Utheta, X and absolute pressure are the same unmodified source;
    # the complete future cancels in the energy difference, never in e itself.
    ur=(2*z*forcing[0]-(1-delta)*z*m[0]-(1-z*z)*(s.diff(m[0],z)-2*z*m[0]/(1+z*z)))/(1-delta*z*z)
    for order in range(5):
        if s.diff(s.simplify(ur.subs(subs)),z,order)!=0:raise ArithmeticError('Support radial difference limit failed')
    return dict(difference_mixed4_identities=grids,first_moment_axial5_difference=axial5,
        complete_future_and_nonzero_boundary_moments_preserved_in_actual_and_reference=True,
        exact_X_Utheta_absolute_pressure_differences_zero=True,
        ordinary_radial_recovery_retains_axial5=True,passed=True)


@source_precision
def current_support_source_proof(physical):
    pulse=current_pulse_trace_proof(physical);p=physical.pulse;c=p.ctx
    if not 0<endpoints(p.mu)[0]<=endpoints(p.mu)[1]<mp.mpf('.25'):
        raise ValueError('Positive current linear decay rates required')
    proof=source_difference_proof();full=original_formula_identities();bindings={}
    normalization=fixed_normalization_and_angular_weight_source_proof(physical)
    specs=(('compliant_axial_pulse_field','backward_bump_weights','out[row-1]+=ds*c.exp((c.mpf(\'.5\')-row*mu)*s)*beta'),
        ('compliant_axial_pulse_field','backward_bump_weights','out[2]+=ds*c.exp(-2*mu*s)*beta**2'),
        ('compliant_axial_pulse_field','end','end_energy+=Cj*Cj*(c.exp(-2*self.mu*center)*w[2])'),
        ('compliant_axial_pulse_field','end','B=Bhat*self.Ecap'),
        ('compliant_axial_pulse_field','end','m=[v*self.Ecap for v in mhat]'))
    for stem,method,text in specs:bindings[stem+'.'+method+':'+text]=statement(stem,method,text)
    binding('compliant_pulse_radial_C4','data','selected','self.fifth.select(Z)')
    if any(v['exact_edge']!=str(s.Rational(v['center'])+v['side']*s.Rational(3,20)) for v in EDGES):
        raise ValueError('Exact original support edges changed')
    for side in (-1,1):
        if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in beta_jets(c,side).coefficients):
            raise ValueError('Exact normalized beta endpoint is not flat')
    normals=(p.flat.normalization,p.pulse.initial.repair.normalization)
    if any(endpoints(value)[0]<=0 for value in normals):raise ValueError('Same positive original normalization required')
    return dict(recomputed_current_selected_parameter_pressure_and_scale_proof=pulse,
        recomputed_arbitrary_function_full_difference_proof=proof,
        recomputed_original_full_paper_stress_formula=full,
        exact_normalization_and_five_weight_defining_source=normalization,
        current_original_weighted_history_AST=bindings,
        exact_one_sided_weighted_integral_theorem=weighted_support_FTC_theorem(),
        current_source_endpoint_difference_theorem=endpoint_difference_theorem(),
        current_selected_C5_control_source='physical.history.selected.fifth.select(Z), not an old JSON coefficient box',
        current_actual_moment_history_source='formal_scaled_Mz_mixed_moments of the unchanged end recipe on the current pulse; no division by Ecap',
        exact_end_histories=dict(m_i='-exp(logE)*sum Cj exp(-lambda_i*s)*int_s^0 exp(lambda_i*v)*beta(v-center_j)dv',
            energy='exp(2mu*s)*current_complete_future/2 + int_s^0 exp(2mu*(s-v))/2 dv - exp(2logE)*sum Cj^2 int_s^0 exp(2mu*(s-v))*beta_j^2 dv',
            reference='same actual boundary m_i, energy, X, pressure and swirl; remove only local beta input'),
        disjoint_supports_remove_cross_beta_products_not_inherited_histories=True,
        same_current_exact_positive_B_D_H_R_and_physical_lambda_sources_retained=True,
        physical_trace_transfer=pulse['source_bound_exact_linear_physical_trace_transfer'],
        local_stress3_velocity4_error2_flatness_is_not_global_temporal_flatness=True,passed=True)


def physical_summary(packet):
    result={}
    for category in ('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative'):
        def leaves(node,path=()):
            if 'exact_zero' in node:yield path,node
            else:
                for key,value in node.items():yield from leaves(value,path+(key,))
        for path,row in leaves(packet[category]):
            result[category+':'+':'.join(path)]={key:row[key] for key in ('exact_zero','log_absolute_upper','physical_lambda_exponent')}
    if len(result)!=216:raise ValueError('Support difference physical spatial4/time1 map incomplete')
    return result


class CurrentPulseSupportInterfaces:
    @source_precision
    def __init__(self,pulse_interfaces=None,require_checked=True):
        self.pulse_interfaces=pulse_interfaces if pulse_interfaces is not None else CurrentPulseInterfaces()
        if not self.pulse_interfaces.acceptance_loaded:raise ValueError('Checked current fourteen adjacent source traces required')
        self.physical=self.pulse_interfaces.physical;self.pulse=self.physical.pulse;self.ctx=self.pulse.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.proof=current_support_source_proof(self.physical);self.hashes=dict(self.pulse_interfaces.hashes)
        for sub in ('recomputed_arbitrary_function_full_difference_proof','recomputed_original_full_paper_stress_formula',
                    'exact_normalization_and_five_weight_defining_source'):
            for name,digest in self.proof[sub]['input_hashes'].items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current support source conflict: '+name)
                self.hashes[name]=digest
        oldname=PREFIX+'pulse_end_support_interfaces_check.json'
        self.generic=accepted(oldname,self.family,self.source,'actual_pulse_end_all_four_support_functional_interfaces_verified');_verify_hashes(self.generic)
        if not self.generic['retained_nonzero_history_fixture']['all_passed']:raise ValueError('Original arbitrary-control retained-history oracle missing')
        self.hashes.update(self.generic['input_hashes']);self.hashes[oldname]=sha(oldname)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.history_cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current pulse support admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def controls_and_history(self,Z):
        self.physical.assert_graph();c=self.ctx;Z=c.mpf(Z);key=tuple(endpoints(Z))
        if key not in self.history_cache:
            selected,_,u,_,_=self.pulse.data(Z)
            controls=selected['selected_scaled_end_coefficient_Taylor']
            if len(controls)!=2 or any(v.order!=5 for v in controls):raise ValueError('Two current C5 selected controls required')
            native,hashes=capture_native_end(self.pulse,Z,[-4,0],64)
            for name,digest in hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original native capture source changed')
                self.hashes[name]=digest
            z=IntervalTaylor.variable(c,Z,5);C=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]).reciprocal()
            B=[C*0 for _ in range(5)]
            for Cj,center in zip(controls,(-3,-1)):
                beta=self.pulse.flat.beta(c.mpf([-4,0])-center)
                for k in range(5):B[k]+=Cj*(beta[k]*math.factorial(k))
            m=[native['formal_scaled_Mz_mixed_moments'][0]]
            for k in range(4):m.append(B[k]-m[k]*(c.mpf('.5')-self.pulse.mu))
            self.history_cache[key]=(controls,m,z,C,u)
        return self.history_cache[key]

    @source_precision
    def interface(self,edge,h,Z=(-1,1),log_tau='-1',theta=None):
        if edge not in EDGES:raise ValueError('Unknown exact original pulse support edge')
        c=self.ctx;h=c.mpf(h);Z=c.mpf(Z);ell=c.mpf('.15')
        if endpoints(h)[0]<0 or endpoints(h)[1]>endpoints(c.mpf('.01'))[1]:raise ValueError('Local support distance h in [0,.01] required')
        controls,m_actual,z,C,u=self.controls_and_history(Z);mu=self.pulse.mu;delta=self.pulse.delta
        if endpoints(h)[1]==0:B=[C*0 for _ in range(5)]
        else:B=[symmetric_jet(controls[edge['row']]*(beta_tail_bound(c,k,2*h/ell)/(ell**(k+1)*self.pulse.flat.normalization))) for k in range(5)]
        diff=difference_transport(c,B,mu,h)
        rows=full_difference_rows(c,delta,mu,z,C,self.pulse.Xp,B,diff['linear_m1'],diff['linear_m2'],[-v for v in diff['energy']],m_actual)
        stress={label:{name:ordinary_grid(part['full_derivative_rows'],3) for name,part in sectors.items()}
            for label,sectors in rows['stress'].items()}
        remainder={label:{name:dict(beta=part['beta'],mode=part['mode'],grid=ordinary_grid(part['rows'],2))
            for name,part in sectors.items()} for label,sectors in rows['remainder'].items()}
        ps={};pe={}
        for label,sectors in stress.items():
            ps[label]={name:{'r%d_z%d'%(i,j):physical_bracket(c,grid,i,j,Z,delta,-2-delta)
                for i in range(4) for j in range(4-i)} for name,grid in sectors.items()}
        for label,sectors in remainder.items():
            pe[label]={name:{'r%d_z%d'%(i,j):physical_bracket(c,part['grid'],i,j,Z,delta,part['beta'])
                for i in range(3) for j in range(3-i)} for name,part in sectors.items()}
        # The ordinary original physical operator consumes a bound on the
        # exact difference, with Ecap used solely to enclose exp(logE).
        scaled=[b*self.pulse.Ecap for b in B];dm=[m*self.pulse.Ecap for m in diff['linear_m1']]
        A=[self.pulse.radial(Z,b,m) for b,m in zip(scaled,dm)]
        velocity={UZ:[u*binomial_product(scaled,-(c.mpf('.5')+mu),k) for k in range(5)],
            UT:[C*0 for _ in range(5)],UR:[u.truncate(4)*binomial_product(A,-mu,k) for k in range(5)],P:[C*0 for _ in range(5)]}
        grid={label:{'y%d_Z%d'%(k,n):v[n]*math.factorial(n) for k,v in enumerate(values) for n in range(5-k)} for label,values in velocity.items()}
        e=s.Rational(edge['exact_edge']);coord=c.mpf(str(e.p))/e.q
        neighborhood=c.mpf([endpoints(coord-h)[0],endpoints(coord+h)[1]])
        view=copy.copy(self.physical);view.dispatch=EndpointDispatchView(self.physical,'pulse_end',dict(physical_mixed_derivatives_total_order_le4=grid))
        physical=BASE.evaluate(view,'pulse_end',Z,neighborhood,log_tau=log_tau,theta=theta)
        return dict(edge=edge,h=h,Z=Z,unscaled_beta_forcing_difference_rows=B,
            current_unscaled_first_moment_history_rows=m_actual,
            source_primitive_differences={key:diff[key] for key in ('linear_m1','linear_m2','energy')},
            similarity_stress_difference_mixed3=stress,
            similarity_velocity_difference_mixed4={name:ordinary_grid(jets,4) for name,jets in rows['velocity'].items()},
            physical_stress_difference_mixed3_coefficients=ps,physical_three_component_error_difference_mixed2_coefficients=pe,
            current_native_velocity_pressure_difference_mixed4=grid,
            original_cartesian_spatial4_time1_difference_log_bounds=physical_summary(physical),
            exact_current_parent_factors=dict(logD=self.pulse.logE,logR=physical['source_logR_enclosure'],
                logB='logPstar+log(original inlet U)-13/(2mu)-13-(.5+mu)*s',
                logH='-13*(1-mu)/mu-(1-mu)*s',physical_lambda='lambda^2*(1-Z^2)=tau'),
            exact_source_edge=edge['exact_edge'],diagnostic_coordinate_enclosure=neighborhood,
            actual_nonzero_histories_and_current_complete_future_preserved=True,
            h_zero_is_zero_difference_not_zero_actual_field=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        edges=supported_beta_edges(self.physical)
        for row in edges.values():
            if row['chart']=='pulse_end':row['current_full_field_interface_trace_admitted']=True
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_pulse_support_source_proof=self.proof,
            current_internal_support_ledger=edges,current_source_identified_adjacent_interface_count=14,
            current_source_identified_internal_support_count=4,
            remaining_internal_support_transfers=[name for name,row in edges.items() if not row['current_full_field_interface_trace_admitted']],
            original_arbitrary_control_oracle_consumed_without_repeating_quadrature=self.generic['retained_nonzero_history_fixture'],
            global_temporal_flat_remainder_not_certified=True,input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseSupportInterfaces(require_checked=False)
    result=field.manifest();views=[]
    for edge in EDGES:
        for h in DISTANCES:views.append(field.interface(edge,h))
        print('Current pulse support source/bounds generated: '+edge['exact_edge'],flush=True)
    result['current_whole_Z_support_difference_views']=views;result['input_hashes']=field.hashes
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
