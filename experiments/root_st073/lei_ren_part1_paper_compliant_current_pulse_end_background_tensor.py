"""Actual current meridional pulse-end tensor and its five tensor traces.

Current selected controls, full future and closed absolute pressure supply
all original sectors. Small scale ratios are cancelled before enclosure.
"""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_flatten_power_background_tensor import (
    CurrentFlattenPowerBackgroundTensor,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,
    copy_jet,source_precision,accepted,_verify_hashes,SourceAST,flatten_remaining_kernels,decay_integral)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import pulse_coefficients,capture_native_end
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import (
    full_physical_identities,pulse_velocity_rows,lift_physical_packet)
from lei_ren_part1_paper_compliant_pulse_end_flatten_join import functional_join_identities
from lei_ren_part1_paper_compliant_pulse_end_support_interfaces import EDGES
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding,function
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_axial_pulse_field import backward_bump_weights
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_pulse_end_background_tensor.json'
RECEIPT=PREFIX+'current_pulse_end_background_tensor_check.json'
GATES=('current_actual_pulse_end_full_meridional_tensor_available',
    'current_actual_pulse_end_three_component_physical_decomposition_available',
    'current_actual_end_flatten_completed_tensor_join_certified',
    'current_actual_four_pulse_end_support_tensor_traces_certified')
TENSOR_KEYS=('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_divergence_mixed2',
    'completed_theta_theta_stress_mixed2','physical_three_component_remainder_mixed2',
    'physical_completed_stress_tensor_cartesian','physical_completed_stress_divergence_cartesian',
    'physical_remainder_cartesian','physical_momentum_residual_decomposition_cartesian')
VIEWS={'whole_end':((-1,1),(-4,0),('-3','-1'),None,'1'),
    'end_inlet':((-1,1),'-4','-1',None,'1'),
    'first_center':(('-.8','.8'),'-3',('-5','-2'),None,'.2'),
    'gap_between_supports':(('-.8','.8'),'-2',('-5','-2'),None,'.2'),
    'second_center':(('-.8','.8'),'-1',('-5','-2'),None,'.2'),
    'flatten_attachment':((-1,1),'0','-1',None,'1'),
    'fresh_first_support':('.537','-2.937','-2.6','.41','.8'),
    'fresh_second_support':('.731','-.947','-2.6','.41','.8')}


def reduced_pressure_at_Rv(field,Z):
    """Same remaining pressure / Ev0^2, without division by a tiny C0."""
    c=field.ctx;h=field.flatten_power.heat;mu=c.mpf(field.pulse.mu)
    L=c.mpf(field.flatten_power.outer.Lrel);pp=1+2*mu;d=c.mpf(h.a)-mu
    kernels=flatten_remaining_kernels(c,mu,Z,c.mpf(0),field.flatten_power.flatten.cells)
    right=copy_jet(c,field.flatten_power.right_data(Z)['P'])
    tail=right*c.exp(8*d-pp*(L+96))/4
    power=decay_integral(c,pp,L-4)*c.exp(-100*pp)/8
    return -(kernels['pressure']/4+tail+power)


def current_pressure_and_amplitude_theorem(field):
    """Exact reduction of the checked flatten/power remaining pressure."""
    a,mu,L,z=s.symbols('a mu Lrel Z',real=True);d=a-mu;p=1+2*a;pp=1+2*mu
    C0=2*s.exp(-d*(L+100));right=s.Function('same_current_angular_minus4_P')(z)
    JP=s.Function('same_flatten_suffix_P')(z)
    I=(1-s.exp(-pp*(L-4)))/pp
    power=s.exp(-p*(L-4))*right+s.exp(-2*d*L)*I/2
    full=s.exp(-100*p)*power+s.exp(-2*d*(L+100))*JP
    reduced=JP/4+s.exp(-100*pp)*I/8+s.exp(8*d-pp*(L+96))*right/4
    checks={}
    if s.simplify(s.expand_power_exp(full-C0**2*reduced))!=0:raise ArithmeticError('Current absolute pulse pressure reduction differs')
    checks['same_full_remaining_P_at_flatten0_equals_C0_squared_times_reduced_future']=True
    v=-L-100;bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a
    thetaR=s.exp(-bp*(100+L))/2
    if s.simplify(s.expand_power_exp(C0*thetaR*s.exp(-bh*v)))!=1:raise ArithmeticError('Current pulse/flatten reference amplitude differs')
    checks['same_B_KR_times_C0_equals_actual_Ev0_at_Rv']=True
    q=1+z*z
    if s.simplify(C0/q-s.exp(d*v)*2/q)!=0:raise ArithmeticError('Current flatten K0 differs from same C0/q')
    checks['same_current_flatten_K_at_zero_is_C0_over_q']=True
    delta=s.symbols('same_current_delta',real=True)
    for i in range(3):
        for j in range(3-i):
            pulse=-3+delta-i+j*(delta-1)
            flatten=-1-delta-i+(j+2)*(delta-1)
            if s.expand(pulse-flatten)!=0:raise ArithmeticError('Current physical remainder exponent differs')
            checks['same_theta_remainder_lambda_power_r%d_z%d'%(i,j)]=True
    checks['same_stress_divergence_diagonal_physical_powers_from_identical_generic_operator']=True
    lp,U,mu0=s.symbols('logP U mu',positive=True)
    for label,left,right in (('Ev0',lp+s.log(U)-13/(2*mu0)-13,lp+s.log(U)-13/(2*mu0)-13),
            ('same_radius',s.Symbol('logRp')+13/mu0,s.Symbol('logRp')+13/mu0)):
        if s.expand(left-right)!=0:raise ArithmeticError('Current pulse source log unit differs')
        checks[label]=True
    asts=SourceAST()
    asts.expression('current_pulse_end_background_tensor','reduced_pressure_at_Rv','tail',
        wanted='right*c.exp(8*d-pp*(L+96))/4')
    asts.expression('current_pulse_end_background_tensor','reduced_pressure_at_Rv','power',
        wanted='decay_integral(c,pp,L-4)*c.exp(-100*pp)/8')
    class_assignment('current_pulse_flatten_source','CurrentFlattenMixedC4','__init__','self.U',"self.inlet.constants['U']")
    terminal=field.atlas.interfaces.proof
    if not terminal['passed'] or not all(terminal['mixed4_identities'].values()):raise ValueError('Current native end-flatten source attachment missing')
    support=field.atlas.pulse_support.proof
    if not support['current_source_endpoint_difference_theorem']['complete_future_and_nonzero_boundary_moments_preserved_in_actual_and_reference']:
        raise ValueError('Current support history source missing')
    return dict(identities=checks,exact_reduced_signed_P_at_Rv='-JP_flatten/4-exp(-100pp)*Ipp(Lrel-4)/8-exp(8d-pp*(Lrel+96))*Pangular_minus4/4',
        current_complete_pressure_source=field.flatten_power.units,
        current_native_terminal_function_and_units=terminal,current_selected_support_source=support,
        same_nonzero_Pin_P0_forward_history_retained=True,
        absolute_pressure_equal_by_same_density_FTC_and_current_Cp_zero=True,
        reduction_before_enclosure_not_tiny_amplitude_cap_division=True,input_hashes=asts.hashes,passed=True)


def exact_edge_capture(field,Z,coordinate,edge,cells):
    """Evaluate source flat beta and full/empty suffixes at an exact edge."""
    asts=SourceAST();fn=copy.deepcopy(asts.method('axial_pulse_field','end'))
    class EdgeEvaluation(ast.NodeTransformer):
        def visit_Assign(self,node):
            targets=[ast.unparse(v) for v in node.targets]
            if 'beta' in targets or 'beta_y' in targets:node.value=ast.parse('c.mpf(0)',mode='eval').body
            if 'w' in targets:node.value=ast.parse('edge_weights(c,self.mu,normal,center,cells)',mode='eval').body
            return self.generic_visit(node)
    fn=EdgeEvaluation().visit(fn);ret=next(v for v in ast.walk(fn) if isinstance(v,ast.Return))
    for key,expression in (('formal_backward_energy_loss_Taylor','end_energy*c.exp(2*self.mu*s)'),
            ('formal_unperturbed_energy_Taylor','future*c.exp(2*self.mu*s)+baseline')):
        ret.value.keywords.append(ast.keyword(arg=key,value=ast.parse(expression,mode='eval').body))
    rational=s.Rational(edge['exact_edge'])
    def weights(c,mu,normal,center,n):
        left=s.Rational(center)-s.Rational(3,20);right=s.Rational(center)+s.Rational(3,20)
        if rational<=left:return backward_bump_weights(c,mu,normal,-4,n)
        if rational>=right:return [c.mpf(0)]*3
        raise ValueError('Declared exact edge lies inside another support')
    env=dict(mp=mp,endpoints=endpoints,edge_weights=weights)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<current exact beta edge, unchanged full histories>','exec'),env)
    return env['end'](field,Z,coordinate,cells),asts.hashes


class CurrentPulseEndBackgroundTensor:
    @source_precision
    def __init__(self,flatten_power=None,require_checked=True,cells=64):
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive directed cell count required')
        self.flatten_power=flatten_power if flatten_power is not None else CurrentFlattenPowerBackgroundTensor()
        if not self.flatten_power.acceptance_loaded:raise ValueError('Checked actual current flatten/power tensor required')
        self.physical=self.flatten_power.physical;self.history=self.physical.history;self.pulse=self.history.selected.pulse
        self.atlas=self.flatten_power.angular.atlas;self.ctx=self.physical.ctx;self.cells=cells
        self.family=self.flatten_power.family;self.source=self.flatten_power.source;self.datum_sha=self.flatten_power.datum_sha
        self.units=current_pressure_and_amplitude_theorem(self)
        self.formulas=self.atlas.pulse_support.proof['recomputed_original_full_paper_stress_formula']
        self.physical_proof=full_physical_identities();self.join= functional_join_identities()
        if not self.join['arbitrary_shared_axial_future_and_pressure_functions'] or not all(self.join['identities'].values()):raise ValueError('Original full end-flatten tensor function join missing')
        self.hashes=dict(self.flatten_power.hashes)
        for proof in (self.units,self.formulas,self.physical_proof,self.join):self.hashes.update(proof['input_hashes'])
        for stem in ('pulse_end_stress_C3','pulse_end_physical_C2','pulse_end_flatten_join','current_pulse_end_background_tensor'):
            name=PREFIX+stem+'.py';digest=sha(name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Checked pulse-end source dependency changed')
            self.hashes[name]=digest
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current pulse-end receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.flatten_power.assert_graph()
        if not (self.flatten_power.acceptance_loaded and self.physical is self.flatten_power.physical and self.history is self.physical.history and self.pulse is self.history.selected.pulse and self.pulse is self.physical.pulse and self.pulse is self.flatten_power.flatten.pulse and self.atlas is self.flatten_power.angular.atlas and self.ctx is self.physical.ctx):
            raise ValueError('Pulse-end tensor must use same checked current selected pulse, not old native dispatcher')

    @source_precision
    def end(self,Z,offset,log_tau='-1',theta=None,viscosity='1',edge=None):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(offset);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<-4 or endpoints(v)[1]>0 or endpoints(nu)[0]<=0:raise ValueError('Current end s[-4,0],Z[-1,1],finite logtau,nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        if edge is not None:
            if edge not in EDGES:raise ValueError('Declared exact current pulse edge required')
            rational=s.Rational(edge['exact_edge']);exact=c.mpf(str(rational.p))/rational.q
            if endpoints(v)!=endpoints(exact):raise ValueError('Exact support source coordinate required')
            native,hashes=exact_edge_capture(self.pulse,Z,v,edge,self.cells)
        else:native,hashes=capture_native_end(self.pulse,Z,v,self.cells)
        self.hashes.update(hashes);selected=self.pulse.data(Z)[0]
        z=IntervalTaylor.variable(c,Z,5);C=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]).reciprocal()
        Bh=[C*0 for _ in range(5)]
        if edge is None:
            for Cj,center in zip(selected['selected_scaled_end_coefficient_Taylor'],(-3,-1)):
                control=copy_jet(c,Cj);beta=self.pulse.flat.beta(v-center)
                for j in range(5):Bh[j]+=control*(c.mpf(beta[j])*math.factorial(j))
        m=[copy_jet(c,native['formal_scaled_Mz_mixed_moments'][0])]
        n=[copy_jet(c,native['formal_scaled_Mz_mixed_moments'][1])]
        e0=[copy_jet(c,native['formal_unperturbed_energy_Taylor'])]
        J=[copy_jet(c,native['formal_backward_energy_loss_Taylor'])]
        key=endpoints(Z)
        if key not in self.cache:self.cache[key]=reduced_pressure_at_Rv(self,Z)
        mu=c.mpf(self.pulse.mu);delta=c.mpf(self.pulse.delta);pp=1+2*mu
        P=[self.cache[key]*c.exp(pp*v)+C*C*(c.expm1(pp*v)/(2*pp))]
        for j in range(4):
            m.append(Bh[j]-m[j]*(c.mpf('.5')-mu));n.append(Bh[j]-n[j]*(c.mpf('.5')-2*mu))
            square=sum((Bh[l]*Bh[j-l]*math.comb(j,l) for l in range(j+1)),C*0)
            e0.append(e0[j]*(2*mu)-(C*0+1/2 if j==0 else C*0));J.append(J[j]*(2*mu)-square)
            P.append(P[j]*pp+(C*C/2 if j==0 else C*0))
        parts=pulse_coefficients(delta,mu,z,C,c.mpf(self.pulse.Xp),Bh,m,n,e0,J,P)
        logB=dict(logPstar=self.physical.logP,actual_log_inlet_U=c.ln(c.mpf(self.flatten_power.flatten.U)),inverse_mu=-13/(2*mu),finite=-13-(c.mpf('.5')+mu)*v)
        logR,radius_source=self.physical.radius('pulse_end',v,{},self.pulse)
        logD=c.mpf(self.pulse.logE);logH=-13*(1-mu)/mu-(1-mu)*v;sectors={}
        for label,rows in parts.items():
            sectors[label]={}
            for name,part in rows.items():
                rp,bp,dp,hp=part['mode'];logs=dict(source_logR=rp*logR,**{key:bp*value for key,value in logB.items()},
                    selected_log_end_scale=dp*logD,signed_original_memory_log=hp*logH,normalization=-c.ln(2)/2)
                grid={'s%d_Z%d'%(j,n):jet[n]*math.factorial(n) for j,jet in enumerate(part['full_derivative_rows']) for n in range(4-j)}
                sectors[label][name]=dict(mode=part['mode'],exact_source_log_parts=logs,source_logR_rate=part['source_logR_rate'],full_stress_mixed3_coefficient_enclosures=grid)
        packet=dict(Z=Z,s=v,full_meridional_stress_log_sectors=sectors,exact_pulse_reference_logB_parts=logB,exact_logR=logR,exact_logD=logD,exact_logH=logH)
        velocity=pulse_velocity_rows(c,delta,mu,z,C,Bh,m)
        physical=lift_physical_packet(c,packet,delta,velocity,lt,theta,nu)
        return dict(physical,current_actual_source_stress_packet=packet,
            source_formal_Bhat_rows=Bh,source_formal_Mz_rows=m,source_formal_Mtheta_z_rows=n,
            source_unperturbed_complete_energy_rows=e0,source_selected_backward_energy_loss_rows=J,
            signed_pressure_relative_pure_swirl_reference_rows=P,original_native_forward_absolute_pressure=native['pressure'],
            original_native_nonzero_complete_future=native['positive_terminal_future_energy_Taylor'],source_radius=radius_source,
            current_source_three_component_velocity_rows=velocity,actual_full_stress_not_local_difference=True,
            current_selected_complete_history_graph_retained=True,pressure_remaining_reduction_before_enclosure=True,
            exact_support_boundary_source=edge,source_bounds_not_resolved_physical_point_values=True,
            actual_completed_tensor_regional_decomposition_not_global_temporal_flatness=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph();left=self.end(Z,0,log_tau,theta,viscosity);right=self.flatten_power.chart('flatten',Z,0,log_tau,theta,viscosity)
        c=self.ctx
        def bound(rows):
            values=[endpoints(row['log_absolute_upper'])[1] for row in rows if not row['exact_zero']]
            if not values:return dict(exact_zero=True,log_absolute_upper=None)
            # Triangle enclosure after exact source function equality. A
            # directed log(n)+max avoids materializing extreme source logs.
            return dict(exact_zero=False,log_absolute_upper=c.mpf(max(values))+c.ln(len(values)))
        for label in ('radial','axial'):
            if any(not row['exact_zero'] for sector in left['physical_three_component_remainder_mixed2'][label].values() for row in sector.values()):raise ValueError('Original empty terminal supports must give Er=Ez=0')
        rows={}
        for key,endkey in (('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_mixed3'),('physical_stress_divergence_mixed2','physical_cylindrical_stress_divergence_mixed2')):
            for label,grid in right[key].items():
                for name,row in grid.items():rows[key+'/'+label+'/'+name]=bound([row]+[v[name] for v in left[endkey][label].values()])
        for name,row in right['completed_theta_theta_stress_mixed2'].items():rows['completed_theta_theta_stress_mixed2/'+name]=bound([row]+[v[name] for v in left['completed_theta_theta_stress_mixed2'].values()])
        for name,row in right['physical_axial_viscosity_remainder_mixed2'].items():rows['physical_axial_viscosity_remainder_mixed2/'+name]=bound([row,left['physical_three_component_remainder_mixed2']['theta']['axial_viscosity'][name]])
        maps=(('completed_background_tensor_cartesian_components','physical_completed_stress_tensor_cartesian'),
            ('physical_completed_stress_divergence_cartesian','physical_completed_stress_divergence_cartesian'),
            ('physical_remainder_cartesian','physical_remainder_cartesian'),
            ('physical_momentum_residual_decomposition_cartesian','physical_momentum_residual_decomposition_cartesian'))
        nt=len(left['physical_cylindrical_stress_mixed3']['theta'])
        nz=len(left['physical_cylindrical_stress_mixed3']['axial'])
        # Preserve the canonical theta/diagonal or divergence/remainder
        # contributions individually. Original sector sums have already
        # been source-identified with each canonical common function.
        for key,endkey in maps:
            for component,values in right[key].items():
                other=left[endkey][component]
                if key=='completed_background_tensor_cartesian_components' and component in ('xx','xy','yy'):
                    groups=(other[:nt],other[nt:])
                elif key=='physical_momentum_residual_decomposition_cartesian':
                    n=nt if component in ('x','y') else nz;groups=(other[:n],other[n:])
                else:groups=(other,)
                if len(groups)!=len(values):raise ValueError('Original source sector/canonical contribution routing differs')
                for i,(row,group) in enumerate(zip(values,groups)):
                    for source in group:
                        if not row['exact_zero'] and not source['exact_zero']:
                            # The equivalent source lambda expressions are
                            # proved above before enclosure. Their directed
                            # copies have different precision/dependency;
                            # byte equality of rounded powers is not proof.
                            power='physical_viscosity_exponent'
                            if encode(pack(row[power]))!=encode(pack(source[power])):raise ValueError('Common source physical viscosity exponent differs')
                    rows[key+'/'+component+'/'+str(i)]=bound([row]+group)
        if len(rows)!=65:raise ValueError('Current end-flatten common physical layout incomplete')
        return dict(seam='end_flatten',common_actual_tensor_rows=rows,current_common_tensor_contribution_count=len(rows),
            same_actual_current_complete_moments_and_pressure=True,
            original_arbitrary_function_full_tensor_join=self.join,
            current_native_terminal_and_amplitude_pressure_reduction=self.units,
            original_full_meridional_physical_operator_identity=self.physical_proof,
            source_function_equality_precedes_common_triangle_bounds=True,interval_overlap_not_used_as_function_identity=True,
            **dict.fromkeys(OPEN,False))

    @source_precision
    def support_interface(self,edge,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if edge not in EDGES:raise ValueError('Exact original pulse support edge required')
        rational=s.Rational(edge['exact_edge']);v=self.ctx.mpf(str(rational.p))/rational.q
        endpoint=self.end(Z,v,log_tau,theta,viscosity,edge=edge)
        def collect(value,path):
            if isinstance(value,dict) and 'signed_coefficient' in value:return {path:dict(exact_zero=value['exact_zero'],log_absolute_upper=value['log_absolute_upper'])}
            result={}
            for key,row in (value.items() if isinstance(value,dict) else enumerate(value)):result.update(collect(row,path+'/'+str(key)))
            return result
        rows={}
        for key in TENSOR_KEYS:rows.update(collect(endpoint[key],key))
        proof=self.atlas.pulse_support.proof
        return dict(edge=edge,common_actual_tensor_rows=rows,current_common_tensor_contribution_count=len(rows),
            exact_empty_full_support_evaluation_before_enclosure=True,actual_nonzero_shared_boundary_histories_preserved=True,
            current_weighted_support_FTC_function_theorem=proof['exact_one_sided_weighted_integral_theorem'],
            current_mixed4_and_axial5_zero_difference_theorem=proof['current_source_endpoint_difference_theorem'],
            current_actual_full_tensor_and_meridional_remainder_endpoint=endpoint,
            both_one_sided_actual_tensor_traces_identified=True,local_difference_not_substituted_for_actual_tensor=True,
            global_temporal_flatness_not_inferred=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_pulse_complete_pressure_units_and_five_history_theorem=self.units,
            original_full_meridional_paper_stress_AST_theorem=self.formulas,original_full_meridional_physical_decomposition_theorem=self.physical_proof,
            original_arbitrary_complete_future_pressure_end_flatten_join_theorem=self.join,
            actual_current_tensor_regions_available=['pulse_end']+self.flatten_power.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=9,actual_current_completed_tensor_internal_interface_count=4,
            current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
            source_domain='end s[-4,0],exact edges(-63/20,-57/20,-23/20,-17/20),Z[-1,1],R>0,|Z|<1,tau>0,finite compact logtau,constant nu>0',
            directed_end_suffix_cells=self.cells,full_actual_selected_meridional_and_radial_remainder_sectors_retained=True,
            remaining_pulse_core_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseEndBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_pulse_end_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_pulse_end_tensor_views'][name]=field.end(*args)
        print('Current actual pulse-end tensor: '+name,flush=True)
    result['current_actual_end_flatten_tensor_interface']=field.interface()
    result['current_actual_four_pulse_end_support_tensor_interfaces']={edge['exact_edge']:field.support_interface(edge) for edge in EDGES}
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
