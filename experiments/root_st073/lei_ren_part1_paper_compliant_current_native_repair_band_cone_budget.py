"""Source-bound partial repair histories and quiet/power cone budgets.

Exact finite Picard functions are kept separate from conditional control-ball
bounds. The latter apply to the genuine fixed point if it is instantiated at
a compatible N; finite Picard terminal residuals are never declared zero.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_generic_cone_state_errors as state

allN, controls = state.allN, state.allN.controls
source, repair, packets = allN.source, allN.repair, allN.packets
HERE, PREFIX, sha, ep = allN.HERE, allN.PREFIX, allN.sha, allN.ep
LogUpper, Poly = state.LogUpper, state.Poly
NAME = PREFIX+'current_native_repair_band_cone_budget.json'
RECEIPT = PREFIX+'current_native_repair_band_cone_budget_check.json'
GATE = 'current_original_repair_band_partial_C1_functions_and_quiet_power_cone_budget_bound'


def profiles_at(built, x):
    """Same original compact bumps, same finite controls and Rc amplitude."""
    g, W, N = built['graph'], built['repair_weights'], built['N']
    h = built['finite_picard_sequence'][-1]
    logx = g.unary('log', x)
    bumps = []
    for center in W['centers']:
        t = g.quotient(g.sub(logx, center), W['ell'], 'original positive log bump width')
        raw = g.node('compact_raw_beta', argument=t.node,
            definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
            defining_module=PREFIX+'outer_pulse_map.py',
            defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))
        bumps.append(g.quotient(raw, g.mul(W['ell'], W['normal'], x), 'original positive ell J0 x'))
    overN = lambda v: g.quotient(v, N, 'one common original positive integer N')
    F = source.C1Function(*[overN(g.add(*(g.mul(b, getattr(h[i+2], row)) for i,b in enumerate(bumps))))
                           for row in ('value','Z')])
    G = source.C1Function(*[overN(g.add(g.mul(bumps[0], getattr(h[0], row)),
                                     g.mul(bumps[2], getattr(h[1], row)))) for row in ('value','Z')])
    power = g.unary('exp', g.neg(g.mul(g.add(g.constant('1/2'), W['mu']), logx)))
    E = g.c1mul(built['amplitude'], source.C1Function(power, g.zero))
    dE, dV = g.c1mul(built['amplitude'], F), g.c1mul(built['amplitude'], G)
    return dict(original_E=E, delta_E=dE, delta_V=dV, F=F, G=G,
                E=g.c1add(E,dE), V=dV, bumps=bumps)


def density_pairs(g, E, dE, dV):
    """Full signed ORIGINAL density differences on the reserved V=0 power."""
    half = lambda q: g.c1scale(g.constant('1/2'), q)
    EF, FF, VV = g.c1mul(E,dE), g.c1mul(dE,dE), g.c1mul(dV,dV)
    return dict(m=dV, h=dE, k=g.c1mul(g.c1add(E,dE),dV),
                e=g.c1add(VV,g.c1scale(g.constant(-1),g.c1add(EF,half(FF)))),
                p=g.c1add(EF,half(FF)))


def exact_partial_band_graph(built):
    """D(x)=x^-rate*(D(Rc)+int_1^x t^(rate-1)*density(t)dt).

    Endpoint and integration variables are distinct. All rates, incoming
    pressure memory and ordinary-Z product rules remain exact functions.
    """
    g = built['graph']; x = g.symbol(built['band_variable'])
    variable = 'repair_partial_integration_t'; t = g.symbol(variable)
    profiles = profiles_at(built,t)
    densities = density_pairs(g,profiles['original_E'],profiles['delta_E'],profiles['delta_V'])
    partial, contributions = {}, {}
    for key, rate in allN.RATES.items():
        kernel = g.unary('exp',g.mul(g.constant(rate-1),g.unary('log',t)))
        decay = g.unary('exp',g.neg(g.mul(g.constant(rate),g.unary('log',x))))
        integral_rows = [controls.integral(g,g.mul(kernel,getattr(densities[key],row)),
            variable,g.one,x,measure='dt; exact t^(rate-1) original recovery weight',
            original_recovery_rate=str(rate),partial_endpoint='repair_x in[1,2]',
            endpoint_independent_of_Z=True) for row in ('value','Z')]
        contribution = source.C1Function(*integral_rows)
        contributions[key] = g.c1scale(decay,contribution)
        partial[key] = g.c1scale(decay,g.c1add(built['history'][key],contribution))
    residual = dict(zip(controls.ROWS,built['control_residual']))
    invN = g.quotient(g.one,built['N'],'one common original positive integer N')
    A = built['amplitude']; AA = g.c1mul(A,A)
    endpoints = {}
    for key,row,amplitude in (('m','M',A),('h','I',A),('e','S',AA),('p','Cp',AA)):
        decay = g.unary('exp',g.neg(g.mul(g.constant(allN.RATES[key]),g.unary('log',g.constant(2)))))
        endpoints[key] = g.c1scale(g.mul(decay,invN),g.c1mul(amplitude,residual[row]))
    joint = g.c1add(residual['M'],g.c1scale(built['parameters']['mu'],residual[controls.ROWS[1]]))
    decay = g.unary('exp',g.neg(g.mul(g.constant(allN.RATES['k']),g.unary('log',g.constant(2)))))
    endpoints['k'] = g.c1scale(g.mul(decay,invN),g.c1mul(AA,joint))
    return dict(**built,partial_band_histories=partial,partial_band_contributions=contributions,
                band_density_pairs=densities,finite_terminal_residual_identities=endpoints,
                partial_integration_variable=variable)


def bump_caps(c, weights):
    """|beta|<=e^-1, |beta'|<=8/e^2, g=beta/(ell J0 x)."""
    ell, J0 = weights['radius'], weights['raw_normalization']
    if ep(ell)[0]<=0 or ep(J0)[0]<=0:
        raise ValueError('Original positive bump width and normalization required')
    value = c.exp(-1)/(ell*J0)
    derivative = 8*c.exp(-2)/(ell*ell*J0)+value
    return dict(value=value,log_radius_derivative=derivative,
        raw_beta_upper=c.exp(-1),raw_beta_derivative_upper=8*c.exp(-2),
        within_swirl_and_within_axial_supports_pairwise_disjoint=True,
        source_normalization_is_enclosure_not_defining_value=True)


def density_envelopes(c, E, EZ, dE, dEZ, dV, dVZ):
    """Uniform C0/Z density caps; source amplitude and its Z row retained."""
    C = lambda x: LogUpper.constant(c,x)
    linear = dE.scale(E); linearZ = dE.scale(EZ)+dEZ.scale(E)
    square = (dE*dE).scale(C('.5')); squareZ = dE*dEZ
    return (dict(m=dV,h=dE,k=dV.scale(E)+dE*dV,
                 e=dV*dV+linear+square,p=linear+square),
            dict(m=dVZ,h=dEZ,k=dV.scale(EZ)+dVZ.scale(E)+dEZ*dV+dE*dVZ,
                 e=(dV*dVZ).scale(C(2))+linearZ+squareZ,p=linearZ+squareZ))


def partial_history_envelopes(c, incoming, incomingZ, density, densityZ):
    """Full band width log2, each OWN positive Duhamel mass, no resets."""
    masses = {}; D = {}; DZ = {}
    for key,rate in allN.RATES.items():
        rr = c.mpf(rate.numerator)/rate.denominator
        mass = c.ln(2) if not rate else -c.expm1(-rr*c.ln(2))/rr
        if ep(mass)[0]<=0:raise ValueError('Own positive band mass required')
        masses[key]=mass
        cap = lambda rows: Poly(c,{p:allN.target_module.magnitude(v) for p,v in rows[key].items()})
        D[key]=cap(incoming)+density[key].scale(LogUpper.constant(c,mass))
        DZ[key]=cap(incomingZ)+densityZ[key].scale(LogUpper.constant(c,mass))
    return D,DZ,masses


def power_margin(c,logRw,offset,mu,p1,p2,proof):
    """Whole original O3 power D>=R*theta_floor, Q/D^2>=2-direction.

    b=0, a=2+2mu; G=(a,a*2mu,a*D,a^3*Q). Never replace
    the signed directional theorem by a small-axial-only assertion.
    """
    if ep(mu)[0]<=0 or ep(mu)[1]>=ep(c.mpf(1)/6)[0]:
        raise ValueError('Same original 0<mu<1/6 required for power/shear caps')
    theta = packets.interval(c,proof['theta_lower_constants']['theta_floor'])
    qfactor = 2-packets.interval(c,proof['full_directional_expression_upper'])
    if ep(theta)[0]<=0 or ep(qfactor)[0]<=0:raise ValueError('Original full signed D,Q reserve required')
    logD = c.mpf(ep(logRw+offset+c.ln(theta))[0])
    logs = dict(G1=c.ln(2),G2=c.ln(4)+c.ln(mu),
                G3=c.ln(2)+logD,G4=c.ln(8)+2*logD+c.ln(qfactor))
    logg = c.mpf(min(ep(v)[0] for v in logs.values()))
    maxstate = c.mpf(max(ep(c.ln(3))[1],*[ep(v.log)[1] for v in (p1,p2) if v.log is not None]))
    logH = LogUpper.add(c,[LogUpper.constant(c,2),LogUpper(c,maxstate)]).log
    logrho = c.mpf(min(0,ep(logg-c.ln(512)-5*logH)[0]))
    return dict(log_original_D_positive_lower=logD,original_Q_over_D_squared_lower=qfactor,
        G_log_positive_lowers=logs,log_g_positive_lower=logg,log_H_state_upper=logH,
        log_rho_stability_positive=logrho,original_power_a='2+2mu',original_power_b_exact_zero=True,
        full_original_energy_pressure_and_directional_correlation_retained=True,
        tolerance='max state error <= min(1,g/(512H^5)) gives every G_j >= g/2')


def band_normalized_state(c,errors,RoverE,S,p1,p2,relative,shape_derivative,shape_value):
    """|F/power|<=1/2; actual positivity supplies reciprocal factor<=2."""
    C = lambda x: LogUpper.constant(c,x)
    p = {}
    for key,axis,old in (('p1','theta',p1),('p2','axial',p2)):
        p[key]=((errors[axis+'_linear']+errors[axis+'_quadratic'].scale(S)).scale(RoverE)
                +relative.scale(old)).scale(C(2))
    p['a']=(shape_derivative+shape_value.scale(C(c.mpf(2)/3))).scale(C(4)*C(c.power(2,c.mpf(2)/3)))
    p['b']=shape_derivative.scale(C(4)*C(c.power(2,c.mpf(2)/3)))
    return p


class NativeRepairBandConeBudget:
    def __init__(self,owner,cone_owner):
        if type(owner) is not allN.NativeRcAllNFunctionControls or type(cone_owner) is not state.NativeGenericConeStateErrors:
            raise TypeError('Accepted live all-N and cone-state owners required')
        if cone_owner.owner is not owner:raise ValueError('Same original all-N owner required')
        checked=json.loads((HERE/state.RECEIPT).read_bytes())
        if not checked['all_passed'] or not checked[state.GATE] or checked['source_family']!=owner.family:
            raise ValueError('Checked same-original cone-state route required')
        powername=PREFIX+'current_O3_power_cone.json'; powercheck=PREFIX+'current_O3_power_cone_check.json'
        powerreceipt=json.loads((HERE/powercheck).read_bytes())
        if not powerreceipt['all_passed']:raise ValueError('Checked original full power cone required')
        self.proof=json.loads((HERE/powername).read_bytes())['current_whole_O3_power_correlated_bounds']
        self.owner,self.ctx,self.family=owner,owner.ctx,owner.family
        self.hashes={**cone_owner.hashes,**checked['input_hashes'],state.NAME:sha(state.NAME),
            state.RECEIPT:sha(state.RECEIPT),powername:sha(powername),powercheck:sha(powercheck),
            Path(__file__).name:sha(Path(__file__).name)}
        owner.service.bind_hashes(self.hashes)

    @allN.paired.native.inlet.source_precision
    def compute(self,live,cone_live):
        c=self.ctx;owner=self.owner;root=owner.target.q_owner.owner.owner
        domain,reservation=owner.target.owner.domain,owner.target.owner.reservation
        if not reservation['original_whole_power_cone_consumed_on_r_plus_through_twice_Rc']:
            raise ValueError('Exact original r_plus,Rc,2Rc reservation required')
        Tw=owner.target.transfer.geometry.binder.fixed['Tw']
        if ep(Tw-2-c.ln(2))[0]<=0:raise ValueError('Whole repair band must be inside original O3 power chart')
        packet_query=owner.target.q_owner.query('O3_power',(-1,1),c.mpf((2,ep(2+c.ln(2))[1]))/Tw)
        roots=packet_query['source']['roots'];packet=packet_query['source']['packet']
        original=state.history_transfer.packet_history_functions(packet,owner.target.coordinates,
            owner.target.transfer.owner.signed_owner)
        mag=allN.target_module.magnitude;C=lambda x:LogUpper.constant(c,x)
        poly=lambda v:Poly(c,{-1:v})
        bump=bump_caps(c,live['weights'])
        radius_record=live['conditions']['formal_control_C1_ball_radius_log']
        if radius_record['exact_zero']:raise ValueError('Positive accepted C1 control ball required')
        rho=LogUpper(c,packets.interval(c,radius_record['log_absolute_upper']))
        shape=poly(rho*C(bump['value']));shapeY=poly(rho*C(bump['log_radius_derivative']))
        A,AZ=mag(live['amplitude']),mag(live['amplitude_Z'])
        dE=dV=shape.scale(A);dEZ=dVZ=shape.scale(LogUpper.add(c,[A,AZ]))
        density,densityZ=density_envelopes(c,A,AZ,dE,dEZ,dV,dVZ)
        D,DZ,masses=partial_history_envelopes(c,live['history'],live['Z_derivatives'],density,densityZ)
        T=LogUpper.add(c,[mag(original['originals']['m']),mag(original['Z_derivatives']['m'])])
        errors=state.inertial_error_envelopes(c,A,C(0),T,dE,dV,D,DZ)
        logRw=packets.interval(c,domain['exact_log_Rw']);mu=live['mu']
        logEfloor=live['logA']-c.mpf(2)/3*c.ln(2)
        logRupper=logRw+2+c.ln(2)
        relative=shape.scale(C(c.power(2,c.mpf(2)/3)))
        p1,p2=mag(roots['p1'][allN.ZERO]),mag(roots['p2'][allN.ZERO])
        bandstate=band_normalized_state(c,errors,LogUpper(c,logRupper-logEfloor),
            LogUpper(c,c.mpf(owner.service.data['logP'])),p1,p2,relative,shapeY,shape)
        margin=power_margin(c,logRw,c.mpf(2),mu,p1,p2,self.proof)
        requirements=state.negative_power_threshold(c,bandstate,margin['log_rho_stability_positive'])
        bandN=c.mpf(max(ep(v['sufficient_log_N_lower'])[1] for v in requirements.values()))
        positivity=c.mpf(max(0,ep(c.ln(2)+rho.log+c.ln(bump['value'])+c.mpf(2)/3*c.ln(2))[1]))
        quiet=[]
        for old,row in zip(live['cells'],cone_live['cells']):
            if row['record']['label']!='power_to_Rc':continue
            if len(row['branches'])!=1 or row['record']['chart']!='O3_power':
                raise ValueError('Whole original quiet power branch required')
            branch=row['branches'][0];source_roots=old['branches'][0]['conditional']['source']['roots']
            qm=power_margin(c,logRw,c.mpf(1),mu,mag(source_roots['p1'][allN.ZERO]),
                mag(source_roots['p2'][allN.ZERO]),self.proof)
            req=state.negative_power_threshold(c,branch['state'],qm['log_rho_stability_positive'])
            qN=c.mpf(max(ep(v['sufficient_log_N_lower'])[1] for v in req.values()))
            quiet.append(dict(record=dict(label='r_plus_to_Rc',original_cell_label='power_to_Rc',
                original_geometry=old['geometry']['record'],original_q_flat_exact_profile_unchanged=True,
                incoming_five_history_and_pressure_memory_not_zeroed=True,
                margin=qm,state_error_polynomials={k:v.record() for k,v in branch['state'].items()},
                requirements=req,sufficient_log_N_lower=qN),margin=qm,state=branch['state'],requirements=req,logN=qN))
        if len(quiet)!=1:raise ValueError('Exact whole quiet interval must be retained')
        combined=c.mpf(max(ep(cone_live['combined_logN'])[1],ep(bandN)[1],ep(positivity)[1],ep(quiet[0]['logN'])[1]))
        return dict(bump=bump,control_radius=rho,density=density,densityZ=densityZ,D=D,DZ=DZ,masses=masses,
            errors=errors,state=bandstate,requirements=requirements,margin=margin,bandN=bandN,
            positivityN=positivity,quiet=quiet,combined_logN=combined,original_query=packet_query,
            original_log_E_positive_lower=logEfloor,original_log_R_cover=logRupper,
            relative=relative,shape=shape,shapeY=shapeY,amplitude_cap=A,amplitude_Z_cap=AZ)


@allN.paired.native.inlet.source_precision
def run(owner,cone_owner,live,cone_live,*,return_live=False):
    began=time.monotonic();adapter=NativeRepairBandConeBudget(owner,cone_owner)
    built=exact_partial_band_graph(owner.build());got=adapter.compute(live,cone_live)
    record=dict(source_family=owner.family,**{GATE:True},exact_function_graph_nodes=built['graph'].nodes,
        exact_partial_band_history_roots=controls.pair_roots(built['partial_band_histories'].values(),built['partial_band_histories'].keys()),
        exact_partial_band_contribution_roots=controls.pair_roots(built['partial_band_contributions'].values(),built['partial_band_contributions'].keys()),
        exact_band_density_roots=controls.pair_roots(built['band_density_pairs'].values(),built['band_density_pairs'].keys()),
        exact_finite_terminal_residual_identity_roots=controls.pair_roots(built['finite_terminal_residual_identities'].values(),built['finite_terminal_residual_identities'].keys()),
        finite_Picard_depth=len(built['finite_picard_sequence'])-1,
        exact_partial_history_formula='x^-rate*(D(Rc)+integral_1^x t^(rate-1)*signed_density(t)dt)',
        partial_endpoint='x=R/Rc in[1,2], Z in[-1,1]',integration_variable=built['partial_integration_variable'],
        finite_terminal_identity='Dm=2^-1*A*res_M/N; Dh=2^-3/2*A*res_I/N; De=2^-1*A^2*res_S/N; Dp=A^2*res_Cp/N; Dk=2^-3/2*A^2*(res_M+mu*res_D)/N',
        bump_C0_and_log_radius_derivative_bounds=got['bump'],control_C1_ball_radius=got['control_radius'].record(),
        band_density_C0_error_polynomials={k:v.record() for k,v in got['density'].items()},
        band_density_Z_error_polynomials={k:v.record() for k,v in got['densityZ'].items()},
        original_own_rate_band_masses=got['masses'],
        band_partial_history_C0_error_polynomials={k:v.record() for k,v in got['D'].items()},
        band_partial_history_Z_error_polynomials={k:v.record() for k,v in got['DZ'].items()},
        band_full_inertial_and_radial_error_polynomials={k:v.record() for k,v in got['errors'].items()},
        band_normalized_state_error_polynomials={k:v.record() for k,v in got['state'].items()},
        original_band_log_E_positive_lower=got['original_log_E_positive_lower'],
        original_band_log_R_cover=got['original_log_R_cover'],original_power_query=got['original_query']['record'],
        original_band_power_margin=got['margin'],band_state_frequency_requirements=got['requirements'],
        band_sufficient_log_N_lower=got['bandN'],band_positivity_sufficient_log_N_lower=got['positivityN'],
        quiet_power_records=[row['record'] for row in got['quiet']],
        source_repair_active_quiet_band_sufficient_log_N_lower=got['combined_logN'],
        band_budget_conditional_on_same_source_control_C1_ball=True,
        existing_original_source_live_owners_reused=True,original_source_ancestor_constructors_called=False,
        full_original_P0_cancels_only_in_difference=True,
        local_quiet_and_band_cone_persistence_conditional_on_recorded_N_and_control_ball=True,
        all_chart_q_flat_input_cone_margins_certified=False,
        finite_Picard_functions_equal_fixed_point=False,finite_terminal_residual_declared_zero=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=adapter.hashes,execution_seconds=time.monotonic()-began,
        scope='Exact source-bound finite repair-band partial C1 functions and terminal residual identities; conditional whole-band control-ball histories/full state error bounds; original signed whole-power margins and sufficient quiet+band frequency budgets. Global q-flat input, chosen N, actual solved controls/tail/terminal closure, outer joins and recursion remain open.')
    (HERE/NAME).write_text(json.dumps(packets.encode(record),indent=2)+'\n',encoding='utf8')
    print('Original repair-band partial functions and quiet/power cone budget produced; terminal closure remains open',flush=True)
    return (record,adapter,built,got) if return_live else record
