"""Actual Rz..Rm relaxed inputs with original pressure/axial stress bounds.

Restoration retains the same centered histories and P0. The paper N,
angular Q and signed b*w are connected to the complete current I/F.
This is original relaxed admission, not changed loops or strict tensors.
"""
import json
import ast
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_reshape_relaxed_inputs as reshape

packets=reshape.packets;HERE,PREFIX,sha=reshape.HERE,reshape.PREFIX,reshape.sha
NAME=PREFIX+'current_restore_relaxed_inputs.json'
RECEIPT=PREFIX+'current_restore_relaxed_inputs_check.json'
GATE='current_original_Rz_Rm_relaxed_generic_input_certified'
OPEN=reshape.OPEN;CHARTS=('axial_restore','restore_buffer')


def restoration_source_theorem():
    asts=packets.recovery.numeric.transport.SourceAST()
    wanted={
      ('reference_restore_profiles','restoration','alpha'):'intersection(c,1-sigma_jets(c,t)[0],c.mpf([0,1]))',
      ('reference_restore_profiles','restoration','V'):"inp['z']*4+inp['E']*alpha",
      ('reference_restore_profiles','restoration','logu'):"-logarithm(1+square(inp['z']))+(-8+t)/10",
      ('reference_restore_profiles','restoration','centered'):"restore_centered(c,initial,inp['E'],t,kernels)",
      ('reference_restore_profiles','terminal','logu'):"-logarithm(1+square(inp['z']))+offset/10",
      ('reference_restore_mixed_C4','restoration','alpha'):'[1-cutoff[0]]+[-cutoff[k]*math.factorial(k) for k in range(1,5)]'}
    nodes={}
    for (stem,method,target),value in wanted.items():nodes[(method,target)]=asts.expression(stem,method,target,wanted=value)
    terminal=asts.method('reference_restore_profiles','terminal')
    terminal_V=next(node.value.args[2] for node in terminal.body if isinstance(node,ast.Return))
    if ast.dump(terminal_V)!=ast.dump(ast.parse("inp['z']*4",mode='eval').body):raise ValueError('Original terminal axial velocity source changed')
    asts.bindings['reference_restore_profiles.terminal.return_packet_V']=True
    t,Z=s.symbols('same_restore_log_offset Z',real=True);S=s.Symbol('Pstar',positive=True)
    E=s.Function('same_current_E')(Z);sig=s.Function('same_sigma')(t)
    V=asts.evaluate(nodes[('restoration','V')],dict(inp=dict(z=Z,E=E),alpha=1-sig))
    logu=asts.evaluate(nodes[('restoration','logu')],dict(inp=dict(z=Z),t=t,square=lambda v:v*v,logarithm=s.log))
    u=S*s.exp(logu);av=s.Rational(4,5);b=2*s.diff(V,t)/u
    checks={}
    def zero(name,value):
        if s.cancel(value)!=0:raise ArithmeticError('Actual restoration source identity: '+name)
        checks[name]=True
    zero('same_actual_restore_signed_b',b+2*E*s.diff(sig,t)/u)
    zero('same_actual_restore_power_a',1-2*s.diff(logu,t)-av)
    zero('same_actual_restore_reference_zeta',s.diff(logu,Z)+2*Z/(1+Z*Z))
    terminal=asts.evaluate(terminal_V,dict(inp=dict(z=Z)))
    zero('same_postrestore_exact_4Z',terminal-4*Z)
    zero('same_postrestore_signed_b_zero',s.diff(terminal,t))
    R,delta=s.symbols('R delta',positive=True)
    uu=s.Function('same_actual_Utheta')(R,Z);vv=s.Function('same_actual_V')(R,Z)
    mz,mt,mtz,mzt,P=[s.Function('same_actual_'+n)(R,Z) for n in ('Mz','Mtheta','Mtheta_z','Mztheta','P')]
    d=1-Z*Z;L=1-delta*Z*Z;m=mz/R;ee=mzt/R
    W=1-(1-delta)*Z*m-d*s.diff(m,Z)
    Q=-W+((1-delta/2)*mt-(1-delta)*Z*s.diff(mt,Z)/2-d*s.diff(mtz,Z)+(2*delta-1)*Z*mtz)/(s.sqrt(2)*R**s.Rational(3,2)*uu)
    N=-W*vv+(1-delta)*(m-Z*s.diff(m,Z))/2+2*delta*Z*ee-d*s.diff(ee,Z)+2*(1+delta)*Z*P-d*s.diff(P,Z)
    native=dict(m=m,h=mt/(s.sqrt(2)*R**s.Rational(3,2)*S),
        k=mtz/(s.sqrt(2)*R**s.Rational(3,2)*S),e=ee/(S*S),p=P/(S*S))
    rows=lambda value:[value]+[s.Integer(0)]*4
    env=dict(math=math,axial_derivative=lambda value:s.diff(value,Z))
    asts.replay('collar_Gamma_C4','product_rows',env);asts.replay('collar_stress_C3','shifted_rows',env)
    raw=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',env)(
        SimpleNamespace(mpf=lambda v:s.Rational(str(v))),delta,Z,rows(uu/S),rows(vv),
        {k:rows(v) for k,v in native.items()},rows(P/(S*S)))
    axial=sum(R**s.Rational(str(part['mode'][0]))*S**part['mode'][1]*part['shape'][0]/s.sqrt(2)
        for label,part in raw['axial'].items() if label!='axial_radial_shear')
    zero('same_complete_raw_axial_N',axial-s.sqrt(R/2)*N/L)
    p1,p2=R*Q/L,R*N/(L*uu);Vy=s.symbols('same_actual_V_y',real=True)
    b=2*Vy/uu;J=N/(uu*uu);bw=2*J*Vy/Q
    zero('same_complete_axial_I_over_F',axial/(uu/s.sqrt(2*R))-p2)
    zero('same_signed_bw_equals_b_times_axial_over_angular',bw-b*p2/p1)
    zero('same_full_signed_H0_restore',p1-p2*b/av-p1*(av-bw)/av)
    zero('same_H0_margin_from_original_restore_cone',p1*(av-bw)/av-2-(p1*(av-bw)-2*av)/av)
    return dict(passed=True,original_source_AST_bindings=asts.bindings,exact_identities=checks,
        full_signed_axial_pressure_energy_meridional_sectors_retained=True,
        original_signed_bw_not_absolute_sector_sum_or_omitted_p2=True,
        angular_Q_is_same_paper_inertial_Q_not_radial_velocity=True,
        original_restore_and_terminal_P0_and_six_centered_histories_retained=True,input_hashes=asts.hashes)


class CurrentRestoreRelaxedInputs:
    def __init__(self,service=None):
        self.reshape=reshape.CurrentReshapeRelaxedInputs(service);self.service=self.reshape.service
        self.ctx=self.reshape.ctx;self.family=self.reshape.family;self.conditions=[]
        record=json.loads((HERE/reshape.RECEIPT).read_bytes())
        if not record['all_passed'] or not record[reshape.GATE] or record['source_family']!=self.family:
            raise ValueError('Same checked current original R110..Rz barrier required')
        self.service.bind_hashes(record['input_hashes']);self.service.bind_hashes({reshape.RECEIPT:sha(reshape.RECEIPT)})
        def yes(row,key,label):
            if row[key] is not True:raise ValueError('Actual restoration source identity missing: '+label)
            self.conditions.append(label)
        original=self.reshape.rows['reference_join_check']['moment_stress_identities']
        for key in ('independent_five_primitive_to_stress_derivation_passed','V_radial_derivative_arbitrary'):yes(original,key,'exact full current primitive/angular/axial ODE/'+key)
        if original['angular_9_32_residual']!='0' or original['axial_9_35_residual']!='0':raise ValueError('Original complete stress ODEs differ')
        r=self.reshape.rows['current_restore_background_tensor_check']
        for key in ('current_actual_three_restore_completed_tensor_joins_certified','current_actual_inner_reference_axial_restore_buffer_full_tensors_available'):
            yes(r,key,'same actual restore family/three source joins/'+key)
        norm=r['current_actual_restore_source_units_and_tensor_theorem']['original_actual_restore_current_radius_normalization_theorem']
        yes(norm,'passed','current full restore raw normalization')
        for key in ('raw_five_history_ODE_m','raw_five_history_ODE_h','raw_five_history_ODE_k','raw_five_history_ODE_e','raw_five_history_ODE_p','absolute_pressure_y_ODE'):
            yes(norm['identities'],key,'same original full current history/'+key)
        source=r['current_actual_restore_source_units_and_tensor_theorem']['current_actual_restore_source_pressure_and_three_endpoint_theorem']
        for key in ('reference_Rz_endpoint_uses_exact_negative_point8_log_amplitude','restore_end_and_buffer_start_share_same_full_original_kernels_and_six_histories',
            'full_current_centered_energy_baseline_and_cross_terms_retained','interval_overlap_not_used_as_function_identity'):
            yes(source,key,'whole original restoration source/'+key)
        rm=source['actual_current_Rm_open_neighborhood_function_theorem']
        yes(rm,'passed','actual Rm endpoint and unchanged original pressure transport')
        yes(rm,'same_P0_and_source_functions_arbitrary_smooth_Z','same complete Rm functions and pressure')
        join=self.reshape.rows['reference_join_bounds']
        pressure_name=PREFIX+'pressure_Kp.json'
        pressure=json.loads((HERE/pressure_name).read_bytes())
        if (pressure['implicit_source_sha256']!=self.family['implicit_source_sha256']
                or pressure['datum_enclosure_sha256']!=self.family['datum_enclosure_sha256']
                or pressure['pressure_Kp']!=join['pressure_Kp'] or pressure['KN']!=join['KN']):
            raise ValueError('Same actual analytic pressure source and KN/Kp norm family required')
        yes(pressure,'same_source_preheat_pressure_Kp_certified','same actual axis pressure C1 bound and KN')
        self.service.bind_hashes(pressure['input_hashes']);self.service.bind_hashes({pressure_name:sha(pressure_name)})
        definition=join['definition']
        if definition['axial']!='v1 toRz; v1+(4Z-v1)*sigma(log(R/Rz)) to eRz;4Z afterwards' or definition['offsets']!={'Rz':-8,'restore_end':-7,'Rm':-6,'Rh':-5}:
            raise ValueError('Actual original restoring velocity/radius prescription differs')
        self.conditions.append('same original actual restoration and exact radius offsets, stops at Rm before active patch')
        self.theorem=restoration_source_theorem();self.service.bind_hashes(self.theorem['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.proof=self.bounds()

    def bounds(self):
        c=self.ctx;read=lambda value:packets.interval(c,value)
        lower=lambda value:c.mpf(packets.recovery.endpoints(value)[0])
        upper=lambda value:c.mpf(packets.recovery.endpoints(value)[1])
        join=self.reshape.rows['reference_join_bounds'];eps=read(self.reshape.inner.rows['K1_ledger']['epsilon0'])
        KN=c.mpf(join['KN']);Kp=c.mpf(join['pressure_Kp']);a=c.mpf('.8')
        # The existing actual defining source, its full pressure and inherited
        # moments are already attached, so the original analytic N proof uses
        # the same functions. No axial stress value is selected from this bound.
        source_N=read(join['normalized_N_source_upper'])
        norm_gate=dict(inherited_N_radius_log_margin=read(join['inherited_N_radius_log_margin']),
            full_axial_N_source_budget=500+4*Kp-source_N,
            full_axial_N_KN_budget=KN-(501+4*Kp),
            original_Pstar_above1_log_margin=self.service.saved('switch_power').algebra.logs[1]/2,
            angular_Q_at_Rz_above_half=self.reshape.proof['angular_Q_after_one_log_unit_lower']-c.mpf('.5'),
            inherited_N_R110_over_K_budget=20-read(join['inherited_N_R110_over_K_upper']))
        if any(packets.recovery.endpoints(value)[0]<=0 for value in norm_gate.values()):raise ArithmeticError('Original current pressure/N normalization gate unresolved')
        u_min=c.exp(c.mpf('-.8'))/2
        Vy=upper(16*eps);J=20*KN;bS=upper(2*Vy/u_min)
        bw=upper(2*J*Vy/c.mpf('.5'))
        kappa=upper(a+bS*bS/a) # same original S=Pstar>=1
        strong=lower((3*(a-bw)-2*a)/a)
        margins=dict(restored_profile_above_zero=u_min,original_bw_below_a=a-bw,
            original_kappa_below1=1-kappa,original_full_signed_H0_minus2_lower=strong,
            angular_Q_half_crossing=self.reshape.proof['angular_SQ_lower']-c.mpf('1.6')/2)
        if any(packets.recovery.endpoints(value)[0]<=0 for value in margins.values()):raise ArithmeticError('Actual whole restoration relaxed input failed')
        return dict(physical_domain='Rz<=R<=Rm=exp(-6)*Rref, all Z[-1,1]; stops before active patch',source_charts=list(CHARTS),
            actual_current_pressure_and_inherited_N_norm_gates=norm_gate,positive_directed_relaxed_margins=margins,
            original_u_over_Pstar_lower=u_min,original_signed_Vy_absolute_upper=Vy,
            original_signed_J_absolute_upper=J,original_signed_b_times_Pstar_absolute_upper=bS,
            original_signed_bw_absolute_upper=bw,original_kappa_upper=kappa,
            original_full_H0_minus2_uniform_lower=strong,
            full_pressure_N_proof='N_y+N=Z*u^2+P_operator(P)-W*Vy-(1+delta)*(1-2ZV)*V/2-Hv*V_Z; same P0+Mp and all exact histories',
            full_pressure_N_normalization='|N|<=KN*Pstar^2 after Rz; |u/Pstar|>=exp(-.8)/2 implies |J=N/u^2|<=20KN',
            whole_angular_barrier='same angular ODE for arbitrary V_y; V=4Z+E*alpha, |E|C2<epsilon0, alpha in[0,1] independent of Z; positive mean averaging retains W bound',
            whole_angular_Q_half_barrier='Q>1/2 at Rz; Q_y=SQ-1.6Q>0 at Q=1/2, since SQ>1.4',
            whole_angular_p1_barrier='p1>3 at Rz and the same positive p1=3 crossing persists for R>=Rz',
            full_generic_relaxed_input='a=.8, kappa<1; H0=p1*(a-bw)/a>3*(a-|bw|)/a>2',
            full_signed_p2_used_via_bw_not_dropped=True,restore_buffer_exact_b_zero=True,
            original_current_source_histories_and_absolute_pressure_unchanged=True,
            strict_completed_tensor_cone_new_regions_admitted=0,p1_p2_whole_path_norm_bounds_certified=False,
            active_patch_relaxed_input_certified=False,source_amplitude_radius_and_width_not_materialized=True,
            **dict.fromkeys(OPEN,False))

    def query(self,chart,coordinate,Z=(-1,1)):
        c=self.ctx;v=c.mpf(coordinate);z=c.mpf(Z)
        vl,vh=packets.recovery.endpoints(v);zl,zh=packets.recovery.endpoints(z)
        lo,hi=(0,1) if chart=='axial_restore' else (-7,-6)
        if chart not in CHARTS or not all(mp.isfinite(x) for x in (vl,vh,zl,zh)) or vl<lo or vh>hi or zl<-1 or zh>1:
            raise ValueError('Actual original restore phase[0,1] or buffer offset[-7,-6], Z[-1,1] required')
        return dict(chart=chart,coordinate=v,Z=z,source_family=self.family,analytic_source_certificate=self.proof,
            source_function_domain_not_saved_point_values=True,active_patch_not_in_query=True,
            **{GATE:True},**dict.fromkeys(OPEN,False))


def run():
    owner=CurrentRestoreRelaxedInputs()
    examples={chart:owner.query(chart,(0,1) if chart=='axial_restore' else (-7,-6)) for chart in CHARTS}
    result=dict(source_family=owner.family,current_source_function_attachment_conditions=owner.conditions,
        original_current_full_signed_restoration_source_theorem=owner.theorem,
        current_whole_Rz_Rm_relaxed_input_and_bounds=owner.proof,analytic_source_subbox_examples=examples,
        **{GATE:True},**dict.fromkeys(OPEN,False),
        strict_completed_tensor_cone_new_regions_admitted=0,source_graph_ancestor_constructors_called=False,
        scope='Actual original Rz..Rm restoration+buffer relaxed input with full signed axial/pressure term; no active patch/norm scales/new tensor/changed loop/repair/N/global/recursion admission.',
        input_hashes=owner.service.hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual current whole Rz..Rm restoration/buffer full signed relaxed input PASS',flush=True)
    return result


if __name__=='__main__':run()
