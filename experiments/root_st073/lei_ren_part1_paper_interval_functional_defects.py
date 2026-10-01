"""Directed physical five-row defects on the fresh axial center family.

Fixed accepted construction parameters, not formal pressure/width atoms.
Full flat-shape kernel bounded analytically for negative B; no saddle sampling.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box,tail_jet,pack
from lei_ren_part1_paper_interval_exit_continuation_enclosure import (
    restore_value,restore_jet,load_receipts,continue_exit,_frozen_state,_driver)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,constant
from lei_ren_part1_paper_candidate_general_center_factory import PARAMETERS
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
HERE=Path(__file__).parent


def abs_upper(c,v):
    lo,hi=endpoints(v)
    return c.mpf(max(abs(lo),abs(hi)))


def logjet(v):
    if v.order<1 or endpoints(v[0])[0]<=0:
        raise ValueError('positive C1 log input required')
    return IntervalTaylor(v.ctx,[v.ctx.ln(v[0]),v[1]/v[0]])


def flat_kernel(c,k,m,B,T):
    """Bound I=int exp(-k*l)expm1(m B sigma(l/T)) dl and I_Z.

    B<=0: abs(expm1(m B sigma))<=m abs(B) sigma;
    dI/dB<=m int exp(-k*l)sigma dl. Split at L<=T/2.
    """
    if endpoints(B[0])[1]>0 or endpoints(T)[0]<=0 or endpoints(k)[0]<=0:
        raise ValueError('negative B and positive k,T required')
    tl,kl=c.mpf(endpoints(T)[0]),c.mpf(endpoints(k)[0])
    candidate=(tl*tl/kl)**(c.mpf(1)/3)
    L=c.mpf(min(endpoints(candidate)[0],endpoints(tl/2)[0]))
    q=L/tl
    J=(c.exp(-(tl/L)**2+1/(1-q)**2)+c.exp(-kl*L))/kl
    E=abs_upper(c,B[0])*m*J
    EZ=abs_upper(c,B[1])*m*J
    jet=IntervalTaylor(c,[c.mpf([-endpoints(E)[1],0]),
        c.mpf([-endpoints(EZ)[1],endpoints(EZ)[1]])])
    return dict(jet=jet,weighted_cutoff_integral_upper=J,split=L,
        entire_kernel_and_first_axial_derivative_enclosed=True)


def restore_constants(c,cells=512):
    values={name:c.mpf(0) for name in ('K1','K16','K2')}
    step=c.mpf(1)/cells
    for i in range(cells):
        t=c.mpf([endpoints(step*i)[0],endpoints(step*(i+1))[1]])
        cutoff=alpha_box(c,t+1,c.mpf(1))
        values['K1']+=c.exp(t)*cutoff*step
        values['K16']+=c.exp(t*c.mpf('1.6'))*cutoff*step
        values['K2']+=c.exp(t)*cutoff*cutoff*step
    return values


def centered_axial_exit(calc,bridge,continuation):
    c=calc.ctx
    zero=calc.u[0].truncate(1)*0
    core=zero+c.mpf(PARAMETERS['j'])
    for n,row in enumerate(calc.u[1:],1):
        core+=row.truncate(1)*c.mpf(4)**n
    core+=tail_jet(calc,'Psi_tail',scale=calc.eps).truncate(1)
    bridge_increment=restore_jet(c,bridge['U_increment'],order=1)
    correction=core+bridge_increment+continuation['delta_u_target']
    # The only additional U change is in the first source switch. Enclose
    # its complete RHS and integrate over hb; the second switch/tail keep U.
    bar=continuation['endpoint_comparison_state']
    hb=c.mpf('.005');se=continuation['endpoint_scaled_radius']
    s=c.mpf([endpoints(c.mpf(100)/calc.eps)[0],endpoints(c.exp(hb)*100/calc.eps)[1]])
    comparison=_frozen_state(bar,bar['phi'],bar['U'],s-se,s*s-se*se)
    driver=_driver(calc,bar['phi'],bar['U'],comparison,s)
    epsilon=continuation['exit_shear_epsilon']
    gr=-(driver['D']*epsilon*c.mpf([0,endpoints(hb)[1]])/2)
    ratio=continuation['normalized_exit_state']['phi']/bar['phi'].truncate(1)
    rhs=-(ratio*gr.exp()*driver['I_z']*epsilon*c.sqrt(s*calc.eps/2)*c.mpf([0,1]))
    switch_increment=rhs*hb
    return correction+switch_increment,dict(core_centered=core,
        bridge_increment=bridge_increment,continuation_increment=continuation['delta_u_target'],
        switch_increment=switch_increment,baseline_4Z_cancelled_algebraically=True)


def assemble(calc,switch,bridge,continuation,restore_cells=512):
    c=calc.ctx
    with mp.workdps(calc.precision+60):
        z=calc.z.truncate(1)
        logC=c.mpf('5e151');logP=c.mpf(14);T=c.mpf('4e152')
        rm_log=10*(logC+logP)-6
        Rm=110*c.exp(rm_log);rho=c.exp(-rm_log);alpha=c.exp(T-rm_log)
        Am=(1+z*z).reciprocal()*c.exp(logP-c.mpf('.6'))
        uref=(1+z*z).reciprocal()*c.exp(-logC)
        scale=Am*(c.sqrt(2)*Rm**c.mpf('1.5'))
        phi=restore_jet(c,switch['normalized_exit_state']['phi'],order=1)
        B0=c.ln(calc.F0)+c.ln(phi[0])+logC+c.ln((1+z*z)[0])+c.ln(c.sqrt(220))
        B1=calc.ell[0]+phi[1]/phi[0]+(1+z*z)[1]/(1+z*z)[0]
        B=IntervalTaylor(c,[B0,B1])
        kernels={name:flat_kernel(c,c.mpf(k),c.mpf(m),B,T) for name,k,m in
            (('theta','1.6','1'),('energy','1.2','2'),('pressure','.2','2'))}
        K=restore_constants(c,restore_cells)
        g,gparts=centered_axial_exit(calc,bridge,continuation)
        moments={key:restore_jet(c,value,order=1) for key,value in switch['physical_moments'].items()}
        Am2=Am*Am
        seed={1:(moments['z']-z*440)/Rm,
            2:(moments['theta_z']-z*moments['theta']*4)/scale,
            3:(moments['theta']-uref*(c.mpf(5)/8*c.sqrt(2)*c.mpf(110)**c.mpf('1.5')))/scale,
            4:(moments['z_theta']-z*moments['z']*8+z*z*1760+uref*uref*(c.mpf(5)/12*110))/(Am2*Rm),
            5:(moments['p']-uref*uref*c.mpf('2.5'))/Am2}
        direct={1:g*(c.exp(-2)-rho),
            2:g*((c.exp(c.mpf('-3.2'))-rho**c.mpf('1.6'))/c.mpf('1.6')),
            3:g*0,4:g*g/Am2*(c.exp(-2)-rho),5:g*0}
        flat={1:g*0,2:g*kernels['theta']['jet']*alpha**c.mpf('1.6'),
            3:kernels['theta']['jet']*alpha**c.mpf('1.6'),
            4:kernels['energy']['jet']*(-c.mpf('.5')*alpha**c.mpf('1.2')),
            5:kernels['pressure']['jet']*(c.mpf('.5')*alpha**c.mpf('.2'))}
        restoration={1:g*(c.exp(-2)*K['K1']),2:g*(c.exp(c.mpf('-3.2'))*K['K16']),
            3:g*0,4:g*g/Am2*(c.exp(-2)*K['K2']),5:g*0}
        rows={str(n):seed[n]+direct[n]+flat[n]+restoration[n] for n in range(1,6)}
        return dict(center_family=calc.center_family,state_sha256=calc.state_hash,
            source_profile_R=110,rows=rows,B=B,g=g,centered_axial_parts=gparts,
            parts={name:{str(n):v for n,v in data.items()} for name,data in
                (('seed',seed),('direct',direct),('flat',flat),('restoration',restoration))},
            restoration_constants=K,restoration_cells=restore_cells,flat_kernel_bounds=kernels,
            P0=calc.p0.truncate(1),accepted_axis_pressure_preserved=True,
            fixed_parameter_functional_family_defects=True,retained_axial_order=1,
            formal_pressure_width_atoms=False,pressure_parameter_remainders_enclosed=False,
            whole_axis_defects_enclosed=False,terminal_five_moment_closure=False,
            five_bump_controls_solved=False,temporal_recursion=False)


def run():
    calc=IntervalComparisonJets(4)
    path=HERE/'lei_ren_part1_paper_interval_exit_switch_enclosure.json'
    raw=json.loads(path.read_text())
    if raw['state_sha256']!=calc.state_hash:raise ValueError('switch core mismatch')
    for name,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('switch dependency changed: '+name)
    bridge,cells=load_receipts(calc)
    continuation=continue_exit(calc,bridge,cells['endpoint'],cells['initial_phi'])
    result=assemble(calc,raw,bridge,continuation)
    result['input_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
        (Path(__file__).name,path.name,'lei_ren_part1_paper_interval_comparison_enclosure.py',
        'lei_ren_part1_paper_interval_exit_continuation_enclosure.py',
        'lei_ren_part1_paper_centered_component_defects.py')}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Five functional defect rows and first axial derivatives enclosed on',calc.center_family,flush=True)
    print('Centered axial correction:',[mp.nstr(v,18) for v in endpoints(result['g'][0])],flush=True)
    return result

if __name__=='__main__':run()
