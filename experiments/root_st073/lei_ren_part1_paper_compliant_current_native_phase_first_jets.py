"""Conditioned implicit inverse and original A/B first slow/phase derivatives.

First signed duals retain the same native factor basis. Both signed Mobius
angle charts and the small-r Fourier series are differentiated. Slow y/Z
derivatives hold phi fixed; the fast phi derivative is reported separately.
"""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_q_slow_jets as slow
import lei_ren_part1_paper_compliant_current_native_conditioned_phase as phase
import lei_ren_part1_paper_compliant_current_native_spatial_phase as spatial

current=slow.current;prior=slow.prior;native=slow.native;packets=slow.packets
HERE,PREFIX,sha=slow.HERE,slow.PREFIX,slow.sha;ep=slow.ep;ZERO=slow.ZERO
NAME=PREFIX+'current_native_phase_first_jets.json'
RECEIPT=PREFIX+'current_native_phase_first_jets_check.json'
GATE='current_original_native_conditioned_phase_and_A_B_first_slow_jets_executed'
DIRECTIONS=('y','Z','x')
OUTPUTS=('A','B_over_Pstar','A_y','A_Z','B_y_over_Pstar','B_Z_over_Pstar','A_phi','B_phi_over_Pstar')


class First:
    def __init__(self,value,derivatives=None):
        self.v=value;self.d={key:value.scalar(0) for key in DIRECTIONS} if derivatives is None else derivatives
    def coerce(self,other):return other if isinstance(other,First) else First(self.v.scalar(other))
    def __neg__(self):return First(-self.v,{key:-value for key,value in self.d.items()})
    def __add__(self,other):
        other=self.coerce(other);return First(self.v+other.v,{key:self.d[key]+other.d[key] for key in DIRECTIONS})
    __radd__=__add__
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        other=self.coerce(other)
        value=current.square(self.v) if self.v is other.v else self.v*other.v
        return First(value,{key:self.d[key]*other.v+self.v*other.d[key] for key in DIRECTIONS})
    __rmul__=__mul__
    def square(self):return First(current.square(self.v),{key:self.v*value*2 for key,value in self.d.items()})
    def divide(self,other,lower):
        other=self.coerce(other);value=self.v.positive_divide(other.v,lower)
        return First(value,{key:(self.d[key]-value*other.d[key]).positive_divide(other.v,lower) for key in DIRECTIONS})
    def sine(self):
        c=self.v.ctx;argument=phase.bounded_value(self.v)
        return First(self.v.scalar(c.sin(argument)),{key:value*c.cos(argument) for key,value in self.d.items()})


def lower_log(value):
    lo,_=ep(value.coefficient)
    if lo<=0:raise ArithmeticError('A retained positive original source factor is required')
    return value.ctx.mpf(ep(value.scale.evaluate()+value.ctx.ln(value.ctx.mpf(lo)))[0])


def source_dual(value,rows=None):
    derivatives={key:value.scalar(0) for key in DIRECTIONS}
    if rows is not None:derivatives.update(y=rows[(1,0)],Z=rows[(0,1)])
    return First(value,derivatives)


def periodic_trig(c,fraction,half=False,cosine=False):
    lo,hi=ep(fraction)
    if lo==hi and lo in (0,mp.mpf('.5'),1):
        if half:
            return c.mpf((1,0,-1)[(0,mp.mpf('.5'),1).index(lo)] if cosine else (0,1,0)[(0,mp.mpf('.5'),1).index(lo)])
        return c.mpf(1 if lo in (0,1) else -1) if cosine else c.mpf(0)
    argument=c.pi*fraction*(1 if half else 2)
    return c.cos(argument) if cosine else c.sin(argument)


def original_angle_denominator(loop,fraction,inverse):
    c=loop.c;rho2=current.square(loop.rho);ar=loop.scalar(abs(loop.r))
    # Forward: r>0 has its small denominator at psi=0. Inverse: r<0
    # has its small denominator at E=0. Keep rho^2, never rounded1-|r|.
    use_cosine=(loop.sign>0) if inverse else (loop.sign<0)
    trig=loop.scalar(periodic_trig(c,fraction,half=True,cosine=use_cosine))
    result=rho2+ar*current.square(trig)*4
    return result.positive_intersection(lower_log(rho2)),lower_log(rho2)


def small_series(loop,r,s,psi):
    """Original W1/W2 plus rigorous first partial derivative tail covers."""
    c=loop.c;M=48;R=c.mpf(max(abs(value) for value in ep(loop.r)));smin=c.mpf('0.9')
    one=First(loop.scalar(1));powers=[one]
    for k in range(M):powers.append(powers[-1]*r)
    W1=First(loop.scalar(0));W2=psi.divide(s*2,c.ln(smin*2))
    for k in range(1,M+1):
        sine=(psi*k).sine()*(c.mpf(1)/k)
        W1+=powers[k-1]*sine
        coefficient=powers[k].divide(s,c.ln(smin))
        if k>=2:coefficient+=powers[k-2]*(c.mpf(k-1)/2)
        W2+=coefficient*sine
    tails1=dict(value=R**M/((M+1)*(1-R)),r=R**(M-1)/(1-R),psi=R**M/(1-R),s=c.mpf(0))
    tails2=dict(value=R**(M+1)/(smin*(M+1)*(1-R))+R**(M-1)/(2*(1-R)),
        r=R**M/(smin*(1-R))+R**(M-2)/2*((M-1)/(1-R)+R/(1-R)**2),
        psi=R**(M+1)/(smin*(1-R))+R**(M-1)/2*(M/(1-R)+R/(1-R)**2),
        s=R**(M+1)/(smin*smin*(M+1)*(1-R)))
    for W,tails in ((W1,tails1),(W2,tails2)):
        symmetric=lambda bound:c.mpf((-ep(bound)[1],ep(bound)[1]))
        W.v=W.v+loop.scalar(symmetric(tails['value']))
        for key in DIRECTIONS:W.d[key]=W.d[key]+r.d[key]*symmetric(tails['r'])+psi.d[key]*symmetric(tails['psi'])+s.d[key]*symmetric(tails['s'])
    return W1,W2,dict(order=M,r_absolute_cover=R,one_minus_r_squared_lower=smin,
        original_series_value_and_first_partial_tail_bounds={'W1':tails1,'W2':tails2})


def conditioned_first_jets(source,qrows,dstar_log,phi):
    loop=phase.ConditionedPhase(source,dstar_log);c=loop.c;value=loop.evaluate(phi,bits=60)
    if value['status']!='enclosed':return dict(record=value,values=None)
    if qrows is None:return dict(record=dict(status='requires_original_q_slow_source_subdivision'),values=None)
    if loop.flat:
        if any(not row.zero for row in qrows.values()):
            return dict(record=dict(status='requires_original_flat_q_jet_consistency'),values=None)
        zero=loop.scalar(0)
        outputs={name:zero for name in OUTPUTS}
        return dict(record=dict(status='enclosed',geometry='flat',original_C0_inverse=value,
            exact_flat_A_B_and_all_first_derivatives_zero=True,
            original_A_B_first_derivative_enclosures={key:v.record() for key,v in outputs.items()},
            slow_derivatives_hold_actual_phi_fixed=True,fast_phi_derivative_reported_separately=True,
            actual_spatial_fast_N_chain_installed=False,density_C1_integrals_installed=False,**dict.fromkeys(packets.OPEN,False)),
            values=outputs,loop=loop)
    selected=value['selected_inverse'];x=selected['coordinate_interval'];chart=selected['chart']
    variables={name:source_dual(rows[ZERO],rows) for name,rows in source['roots'].items()}
    a,t0,E=variables['a'],variables['t0'],variables['E'];q=source_dual(loop.q,qrows)
    zero=loop.scalar(0);one=loop.scalar(1)
    u=variables['p2']*q
    u=u.divide(First(loop.dstar),dstar_log)
    hinv=First(loop.hinv,{key:-loop.u*loop.hinv*loop.hinv*loop.hinv*u.d[key] for key in DIRECTIONS})
    s=hinv.square();r_value=loop.scalar(loop.r) if loop.geometry=='signed_Mobius' else loop.u*loop.hinv
    r=First(r_value,{key:loop.hinv*loop.hinv*loop.hinv*u.d[key] for key in DIRECTIONS})
    nu=1+t0.square()+q.square()*2
    xf=First(loop.scalar(x),dict(y=zero,Z=zero,x=one));tail=None
    if loop.geometry=='signed_Mobius':
        psif,Ef=loop.angles(x,chart)
        D,logD=original_angle_denominator(loop,x,inverse=chart=='E')
        ds=s.v.positive_divide(D,logD)
        dr=loop.scalar(periodic_trig(c,x)*( -1 if chart=='E' else 1)/c.pi).positive_divide(D,logD)
        mapped=First(loop.scalar(psif if chart=='E' else Ef),{key:dr*r.d[key]+(ds if key=='x' else zero) for key in DIRECTIONS})
        psi,angle=(mapped,xf) if chart=='E' else (xf,mapped)
        divr=lambda numerator:numerator.divide(r if loop.sign>0 else -r,c.ln(c.mpf('.1')))*(1 if loop.sign>0 else -1)
        difference=angle-psi;sinE=(angle*(2*c.pi)).sine()
        phase_numerator=(1+t0.square())*psi+divr(t0*q*hinv*difference*2)
        phase_numerator+=divr(divr(q.square()*((2-s*3)*angle+s*psi+r*sinE*(1/c.pi))))
        F=phase_numerator.divide(nu,0)
        numerator=divr(t0*q*hinv*difference*2)+divr(divr(q.square()*((2-s*3)*difference+r*sinE*(1/c.pi))))
        A=(a*numerator).divide(nu,0)*c.mpf('.5')
        B=E*(t0*A-divr(a*q*hinv*difference)*c.mpf('.5'))
        if chart=='E':
            w=t0.v*loop.hinv+loop.q*(periodic_trig(c,x,cosine=True)+loop.r)*2
            numerator=s.v+current.square(w)
            numerator=numerator.positive_intersection(lower_log(s.v))
            Fx=numerator.positive_divide(nu.v*D,logD)
            logFx=lower_log(s.v)-c.ln(4)-c.mpf(ep(nu.v.record()['log_absolute_upper'])[1])
        else:
            cosine=periodic_trig(c,x,cosine=True)
            w=loop.scalar(cosine)-r.v
            direction=loop.t0+(loop.q*loop.hinv*w*2).positive_divide(D,logD)
            Fx=(1+current.square(direction)).positive_divide(nu.v,0)
            logFx=-c.mpf(ep(nu.v.record()['log_absolute_upper'])[1])
    else:
        psi=xf;rad=psi*(2*c.pi);W1,W2,tail=small_series(loop,r,s,rad)
        F=((1+t0.square())*rad+t0*q*hinv*W1*4+q.square()*s*W2*4).divide(nu,0)*(1/(2*c.pi))
        numerator=t0*q*hinv*W1*4+q.square()*(s*W2*4-psi*(4*c.pi))
        A=(a*numerator).divide(nu,0)*(1/(4*c.pi))
        B=E*(t0*A-a*q*hinv*W1*(1/(2*c.pi)))
        cosine=periodic_trig(c,x,cosine=True)
        D=1-r.v*cosine*2+current.square(r.v);logD=c.ln(c.mpf('.5'))
        direction=loop.t0+(loop.q*loop.hinv*(loop.scalar(cosine)-r.v)*2).positive_divide(D,logD)
        Fx=(1+current.square(direction)).positive_divide(nu.v,0)
        logFx=-c.mpf(ep(nu.v.record()['log_absolute_upper'])[1])
    implicit={key:(-F.d[key]).positive_divide(Fx,logFx) for key in ('y','Z')}
    outputs=dict(A=A.v,B_over_Pstar=B.v)
    for primitive,name in ((A,'A'),(B,'B')):
        for key in ('y','Z'):outputs[name+'_'+key+('_over_Pstar' if name=='B' else '')]=primitive.d[key]+primitive.d['x']*implicit[key]
        outputs[name+'_phi'+('_over_Pstar' if name=='B' else '')]=primitive.d['x'].positive_divide(Fx,logFx)
    philo,phihi=ep(c.mpf(phi))
    if philo==phihi and philo in (0,mp.mpf('.5'),1):
        for name in ('A','B_over_Pstar','A_y','A_Z','B_y_over_Pstar','B_Z_over_Pstar'):outputs[name]=zero
        implicit={key:zero for key in ('y','Z')}
    record=dict(status='enclosed',geometry=loop.geometry,original_C0_inverse=value,
        selected_chart=chart,original_positive_fraction_phase_derivative=Fx.record(),
        positive_fraction_phase_derivative_lower_log=logFx,
        implicit_inverse_first_slow_derivative_enclosures={key:v.record() for key,v in implicit.items()},
        original_A_B_first_derivative_enclosures={key:v.record() for key,v in outputs.items()},
        small_r_value_and_first_derivative_tail_theorem=tail,
        slow_derivatives_hold_actual_phi_fixed=True,fast_phi_derivative_reported_separately=True,
        original_direction_is_cosine_Poisson_shape=True,positive_rho_and_s_not_subtracted_from_one=True,
        actual_spatial_fast_N_chain_installed=False,density_C1_integrals_installed=False,**dict.fromkeys(packets.OPEN,False))
    return dict(record=record,values=outputs,loop=loop)


class NativePhaseFirstJets:
    def __init__(self,owner):
        if type(owner) is not slow.NativeQSlowJets:raise ValueError('Same original q slow-jet source owner required')
        receipt=json.loads((HERE/slow.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[slow.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Accepted same-source q slow jets required')
        self.owner=owner;self.native=owner.native;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (slow.RECEIPT,Path(__file__).name,
            PREFIX+'current_native_conditioned_phase.py')})
        self.dstar=packets.interval(self.ctx,owner.owner.owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        self.binder=spatial.NativeSpatialPhase(self.native)

    @native.inlet.source_precision
    def query(self,chart,Z,coordinate,phi):
        qsource=self.owner.query(chart,Z,coordinate)
        got=conditioned_first_jets(qsource['source'],qsource['rows'],self.dstar,phi)
        record=dict(got['record'],source_family=self.family,chart=chart,
            source_provenance=qsource['source']['packet'].provenance,original_q_slow_jet_source=qsource['record'],
            declared_fractional_phase=self.ctx.mpf(phi),phase_is_independent_input_on_this_query=True,
            original_ordinary_y_Z_derivatives=True,source_native_width_or_Pstar_conversion_not_reapplied=True,
            global_C1_histories_or_Rc_repair_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,values=got['values'],source=qsource,loop=got.get('loop'))

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        geometry=self.binder.query(chart,Z,coordinate,N)
        # Query the original q/source once at the exact coordinate decoded by
        # the radius binder. All derivative rows share its factor basis/ledger.
        qsource=self.owner.query(chart,Z,geometry['raw']['coordinate']);cells=[]
        for phi in geometry['phase_boxes']:
            got=conditioned_first_jets(qsource['source'],qsource['rows'],self.dstar,phi)
            values=got['values'];chain=None
            if values is not None:
                chain=dict(A_y_total=values['A_y']+values['A_phi']*N,
                    B_y_total_over_Pstar=values['B_y_over_Pstar']+values['B_phi_over_Pstar']*N,
                    A_Z_total=values['A_Z'],B_Z_total_over_Pstar=values['B_Z_over_Pstar'])
            record=dict(got['record'],derived_spatial_fractional_phase_box=phi,
                phase_is_independent_input_on_this_query=False,actual_spatial_fast_N_chain_installed=chain is not None,
                original_total_spatial_first_derivative_enclosures=None if chain is None else {k:v.record() for k,v in chain.items()},
                chain_rule='D_y P(y,Z,N*y)=P_y|phi+N*P_phi; D_Z P=P_Z|phi since radius phase is Z independent',
                global_C1_histories_or_Rc_repair_admitted=False,**dict.fromkeys(packets.OPEN,False))
            cells.append(dict(record=record,values=values,chain=chain))
        record=dict(source_family=self.family,chart=chart,source_provenance=qsource['source']['packet'].provenance,
            original_q_slow_jet_source=qsource['record'],actual_original_radius_phase=geometry['record'],
            actual_spatial_phase_and_first_derivative_cells=[cell['record'] for cell in cells],
            candidate_N=N,original_radius_phase_is_Z_independent=True,
            explicit_candidate_N_not_global_admission=True,global_C1_histories_or_Rc_repair_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,source=qsource,geometry=geometry)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativePhaseFirstJets(slow.NativeQSlowJets(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge)))))
        points=dict(inner_reference='.1337',O2_slope='.1337',O2_buffer='5.337',O3_slope_mu='.537',O3_power={'original_power_offset':'.537'})
        fixed={};actual={};N=1024
        for chart,coordinate in points.items():
            geometry=owner.binder.query(chart,('.5','.5'),coordinate,N)
            fixed[chart]=owner.query(chart,('.5','.5'),geometry['raw']['coordinate'],'.137')['record']
            actual[chart]=owner.spatial_query(chart,('.5','.5'),coordinate,N)['record']
            print('Original A/B first jets and actual spatial chain:',chart,fixed[chart]['status'],flush=True)
        endpoints={phi:owner.query('O2_slope',('.5','.5'),'.1337',phi)['record'] for phi in ('0','.5','1')}
        local=owner.spatial_query('O2_slope',('.49','.51'),owner.ctx.mpf(('.13369999','.13370001')),N)['record']
        print('Whole local original radial/Z cell A/B first jets and spatial chain:',local['actual_spatial_phase_and_first_derivative_cells'][0]['status'],flush=True)
        if not all(record['status']=='enclosed' for record in fixed.values()):raise ArithmeticError('Declared native first jet source unresolved')
        if not all(cell['status']=='enclosed' for record in [*actual.values(),local] for cell in record['actual_spatial_phase_and_first_derivative_cells']):
            raise ArithmeticError('Declared native spatial first jet source unresolved')
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=N,
        native_fixed_phase_first_jet_records=fixed,native_actual_spatial_first_jet_records=actual,
        native_symmetry_phase_first_jet_records=endpoints,actual_whole_integral_cell_spatial_first_jet_record=local,
        original_A_B_first_slow_and_fast_derivative_functions_installed=True,
        actual_candidate_spatial_fast_N_chain_installed=True,first_derivative_order_only=True,
        density_C1_integrals_or_global_C1_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Same-source q slow jets feed signed Mobius/small-r conditioned implicit inverse and A/B first y/Z/phi derivatives, with actual N*y spatial chain at candidate N=1024 on five declared boxes and one whole local O2 radial/Z cell. Three exact symmetry phases retained. No higher phase jets, density C1 integration, cumulative inlet-to-Rc histories, repair/common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
