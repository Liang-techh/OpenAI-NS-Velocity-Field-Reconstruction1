"""Current full gap/gap-end tensors with exact reciprocal source endpoints.

The zero axial input retains the selected nonzero incoming histories.
Ordinary derivatives and reduced source logs precede physical enclosure.
"""
import ast
import json
import math
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor import (
    CurrentPulseEndBackgroundTensor,HERE,PREFIX,OPEN,TENSOR_KEYS,sha,pack,encode,
    endpoints,copy_jet,source_precision,accepted,_verify_hashes,SourceAST,
    IntervalTaylor,backward_bump_weights,reduced_pressure_at_Rv,
    pulse_velocity_rows,lift_physical_packet)
from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import gap_shapes,source_proof
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import pulse_remainder_sectors
from lei_ren_part1_paper_compliant_flatten_stress_C3 import axial_derivative

NAME=PREFIX+'current_pulse_gap_background_tensor.json'
RECEIPT=PREFIX+'current_pulse_gap_background_tensor_check.json'
GATES=('current_actual_gap_gapend_full_tensors_available',
    'current_actual_gap_gapend_three_component_decomposition_available',
    'current_actual_gapend_end_completed_tensor_join_certified',
    'current_actual_reciprocal_gap_completed_tensor_join_certified')
VIEWS={'whole_gap':('pulse_gap',(-1,1),(11,12),('-3','-1'),None,'1'),
    'whole_gap_end':('pulse_gap_end',(-1,1),(0,1),('-3','-1'),None,'1'),
    'gap_left':('pulse_gap',(-1,1),11,'-1',None,'1'),
    'gap_coordinate':('pulse_gap',(-1,1),12,'-1',None,'1'),
    'reciprocal_gap_end':('pulse_gap_end',(-1,1),0,'-1',None,'1'),
    'end_attachment':('pulse_gap_end',(-1,1),1,'-1',None,'1'),
    'fresh_gap':('pulse_gap','.537','11.337','-2.6','.41','.8'),
    'fresh_gap_end':('pulse_gap_end','.731','.947','-2.6','.41','.8')}
SEAMS=('gap_coordinate','gap_end')


def current_gap_source_theorem(field):
    """Instantiate the arbitrary-history splitter on the admitted current graph."""
    generic=source_proof({})
    if not all(generic['identities'].values()):raise ValueError('Original arbitrary gap history splitter failed')
    current=field.atlas.pulse_interfaces.proof
    mixed=current['source_bound_five_current_pulse_mixed4_theorem']
    if not current['passed'] or not mixed['passed'] or not all(current['exact_coordinate_reductions'].values()):
        raise ValueError('Checked current reciprocal and end primitive function traces required')
    asts=SourceAST()
    asts.expression('axial_pulse_field','_gap','future',wanted="self.selection.future.future(Z)['complete_future_energy_Taylor']/2")
    asts.expression('axial_pulse_field','_gap','full_weights',wanted='backward_bump_weights(c,self.mu,normal,-4)')
    asts.expression('axial_pulse_field','_gap','end_energy',wanted='Cj*Cj*Kj',augmented=True)
    asts.expression('axial_pulse_field','gap','D',wanted='13-xi')
    asts.expression('axial_pulse_field','gap_from_end','D',wanted='-self.mu*s')
    asts.expression('current_pulse_gap_background_tensor','chart','future',wanted="copy_jet(c,self.pulse.selection.future.future(Z)['complete_future_energy_Taylor'])/2")
    asts.expression('current_pulse_gap_background_tensor','chart','P0',wanted='reduced_pressure_at_Rv(self.end_tensor,Z)')
    fn=asts.method('current_pulse_gap_background_tensor','chart')
    expected=ast.dump(ast.parse("Cj*(c.exp((c.mpf('.5')-i*mu)*center)*self.weights[i-1])",mode='eval').body)
    terms=[node for node in ast.walk(fn) if isinstance(node,ast.AugAssign) and isinstance(node.op,ast.Sub)
        and ast.unparse(node.target)=='M[i - 1]' and ast.dump(node.value)==expected]
    if len(terms)!=1:raise ValueError('Current selected negative full beta moment source changed')
    asts.bindings['current_selected_negative_full_linear_beta_weight']=True
    asts.expression('current_pulse_gap_background_tensor','chart','J0',wanted='Cj*Cj*(c.exp(-2*mu*center)*self.weights[2])',augmented=True)
    mu,z,delta,phase=s.symbols('mu Z delta phase',real=True)
    d=4*mu+(1-4*mu)*(1-phase);lam=s.Rational(1,2)-mu
    checks={}
    selected=field.history.selected
    if not (field.pulse.pulse is selected.amplitude.pulse and selected.amplitude.log_end_scale is field.pulse.pulse.logscale):
        raise ValueError('Current selected logE must retain the actual shared outer pulse map logscale')
    asts.expression('axial_amplitude_selection','__init__','self.log_end_scale',wanted='self.pulse.logscale')
    asts.expression('current_selected_energy_source','__init__','self.pulse.pulse',wanted='self.amplitude.pulse')
    L0,u0=s.symbols('same_current_saddle_L same_current_saddle_u0',positive=True)
    cs=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),ln=s.log,exp=s.exp)
    logpref=asts.evaluate(asts.expression('outer_pulse_map','pulse_rows','logpref'),dict(k0=1/(2*mu),L=L0,u0=u0,mu=mu,c=cs))
    logscale=asts.evaluate(asts.expression('outer_pulse_map','__init__','self.logscale'),
        dict(self=SimpleNamespace(rows={'common_logpref':logpref},mu=mu),c=cs))
    if s.simplify(logscale-(-1/mu-3*L0+2*s.log(u0)-s.log(6)/2-2*s.log(mu)))!=0:
        raise ArithmeticError('Actual current pulse.logE does not equal the reduced D0 source')
    asts.expression('current_pulse_gap_background_tensor','chart','finite',wanted='-3*L+2*c.ln(u0)-c.ln(6)/2-2*c.ln(mu)')
    checks['current_live_logE_to_same_outer_map_to_reduced_D0_AST_chain']=True
    # Every cell count encloses the same defining integral. In particular
    # 64 and 256 are partitions, not different exact beta functions. The
    # current end and gap also share the actual bound evaluation cells.
    asts.expression('axial_pulse_field','backward_bump_weights','beta',wanted='raw_beta(c,r)/(ell*normalization)')
    asts.expression('axial_pulse_field','backward_bump_weights','out[2]',wanted='ds*c.exp(-2*mu*s)*beta**2',augmented=True)
    asts.expression('outer_pulse_map','correction_basis','gram',wanted='(c.mpf(2)/cells)*beta**2*c.exp(-2*mu*s)/(ell*normalization**2)',augmented=True)
    rr,N=s.symbols('normalized_support_coordinate same_fixed_normalization',real=True)
    ell=s.Rational(3,20);raw_beta=s.Function('same_raw_beta')(rr)
    if s.simplify(ell*s.exp(-2*mu*ell*rr)*(raw_beta/(ell*N))**2-s.exp(-2*mu*ell*rr)*raw_beta**2/(ell*N**2))!=0:
        raise ArithmeticError('Full backward square weight and selected energy Gram defining integral differ')
    k,ncells=s.symbols('cell cells',integer=True,positive=True)
    mass=s.summation(ell*((-1+2*(k+1)/ncells)-(-1+2*k/ncells)),(k,0,ncells-1))
    if s.simplify(mass-2*ell)!=0:raise ArithmeticError('Weighted support partition does not cover the whole defining integral')
    checks['same_exact_beta_square_integral_by_change_of_variable']=True
    checks['positive_cell_partition_additivity_independent_of_64_or_256']=True
    for label,value in (('phase0_exact_reciprocal_distance',d.subs(phase,0)-1),
            ('phase1_exact_end_distance',d.subs(phase,1)-4*mu),
            ('same_original_ordinary_coordinate',-mu*s.diff(-s.Symbol('distance')/mu,s.Symbol('distance'))-1)):
        if s.expand(value)!=0:raise ArithmeticError('Current reduced gap coordinate differs')
        checks[label]=True
    # Instantiate the original full histories with one actual current set
    # of selected controls and weighted beta integrals. The end and gap use
    # the same quadrature source/cell count; these remain enclosures of the
    # original exact defining integrals, not coefficient point values.
    controls=[s.Function('same_current_selected_C'+str(j))(z) for j in (1,2)]
    weights=s.symbols('same_W1 same_W2 same_Wsquare',real=True)
    future=s.Function('same_current_complete_C5_future_over2')(z)
    P0=s.Function('same_current_reduced_absolute_Rv_P')(z)
    Xp=s.symbols('same_current_native_forward_Xp',real=True);pp=1+2*mu
    def equal(label,left,right):
        if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError('Current actual gap endpoint history differs: '+label)
        checks[label]=True
    M=[]
    for i in (1,2):
        rate=s.Rational(1,2)-i*mu
        base=-sum(Cj*s.exp(rate*center)*weights[i-1] for Cj,center in zip(controls,(-3,-1)))
        end=-sum(Cj*s.exp(rate*(center+4))*weights[i-1] for Cj,center in zip(controls,(-3,-1)))
        M.append(base)
        for j in range(5):equal('same_current_m'+str(i)+'_ordinary_row'+str(j),base*s.exp(4*rate)*(-rate)**j,end*(-rate)**j)
    J0=sum(Cj**2*s.exp(-2*mu*center)*weights[2] for Cj,center in zip(controls,(-3,-1)))
    end_J=sum(Cj**2*s.exp(-2*mu*(center+4))*weights[2] for Cj,center in zip(controls,(-3,-1)))
    distance=s.symbols('distance',real=True);C=1/(1+z*z)
    e=future*s.exp(-2*distance)+(1-s.exp(-2*distance))/(4*mu)
    J=J0*s.exp(-2*distance)
    P=-C*C/(2*pp)+(P0+C*C/(2*pp))*s.exp(-pp*distance/mu)
    X=1/(1-mu)+(Xp-1/(1-mu))*s.exp(-(1-mu)*(13-distance)/mu)
    for j in range(5):
        equal('same_current_unperturbed_energy_ordinary_row'+str(j),
            ((-mu)**j*s.diff(e,distance,j)).subs(distance,4*mu),
            future*s.exp(-8*mu)+(1-s.exp(-8*mu))/(4*mu) if j==0 else (future-1/(4*mu))*s.exp(-8*mu)*(2*mu)**j)
        equal('same_current_quadratic_loss_ordinary_row'+str(j),
            ((-mu)**j*s.diff(J,distance,j)).subs(distance,4*mu),end_J*(2*mu)**j)
        equal('same_current_pressure_ordinary_row'+str(j),
            ((-mu)**j*s.diff(P,distance,j)).subs(distance,4*mu),
            P0*s.exp(-4*pp)+C*C*(s.exp(-4*pp)-1)/(2*pp) if j==0 else (P0+C*C/(2*pp))*s.exp(-4*pp)*pp**j)
        equal('same_current_X_ordinary_row'+str(j),
            ((-mu)**j*s.diff(X,distance,j)).subs(distance,4*mu),
            1/(1-mu)+(Xp-1/(1-mu))*s.exp(-13*(1-mu)/mu+4*(1-mu)) if j==0 else
            (Xp-1/(1-mu))*s.exp(-13*(1-mu)/mu+4*(1-mu))*(-(1-mu))**j)
    finite,logRp,logP,logU=s.symbols('same_finite same_logRp same_logP same_logU',real=True)
    logD0=-1/mu+finite
    for i in (1,2):
        logDi=(distance/2-1)/mu+finite-i*distance
        equal('same_current_logD'+str(i)+'_end_scale_mode_bridge',logDi.subs(distance,4*mu),logD0+4*(s.Rational(1,2)-i*mu))
    equal('same_current_gap_end_logR',logRp+(13-distance)/mu,logRp+13/mu-distance/mu)
    equal('same_current_gap_end_logB',logP+logU+(-13+distance)/(2*mu)-13+distance,
        logP+logU-13/(2*mu)-13+(s.Rational(1,2)+mu)*distance/mu)
    equal('same_current_gap_end_logH',-(1-mu)*(13-distance)/mu,-13*(1-mu)/mu+(1-mu)*distance/mu)
    equal('same_current_Q_at_end',(-pp*distance/mu).subs(distance,4*mu),-4*pp)
    # The full radial velocity and every remainder sector are homogeneous
    # in the original first moment. D1 is a frozen basepoint factor; rows
    # were already differentiated in ordinary logR, not in phase.
    M=s.Function('same_current_full_linear_beta_weight')(z);C=1/(1+z*z)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)));zero=[s.Integer(0)]*5
    asts.replay('collar_Gamma_C4','product_rows',env:={})
    asts.replay('collar_stress_C3','shifted_rows',env)
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_operators
    env.update(mp=SimpleNamespace(mpf=lambda v:s.Rational(str(v))),math=math,
        axial_derivative=lambda v:s.diff(v,z),physical_operators=physical_operators)
    for fn in ('pulse_velocity_rows','axial_n','axial_operator_rows','pulse_remainder_sectors'):
        asts.replay('pulse_end_physical_C2',fn,env)
    velocity=env['pulse_velocity_rows'](c,delta,mu,z,C,zero,[M*(-lam)**j for j in range(5)])
    endpoint=env['pulse_velocity_rows'](c,delta,mu,z,C,zero,[M*s.exp(4*lam)*(-lam)**j for j in range(5)])
    for label,rows in velocity.items():
        factor=s.exp(4*lam) if label in ('radial','axial') else 1
        for j,value in enumerate(rows):
            if s.simplify(s.expand_power_exp(endpoint[label][j]-value*factor))!=0:raise ArithmeticError('Actual gap end velocity row differs')
            checks['same_actual_'+label+'_ordinary_row'+str(j)]=True
    errors=env['pulse_remainder_sectors'](c,delta,mu,z,velocity)
    endpoint_errors=env['pulse_remainder_sectors'](c,delta,mu,z,endpoint)
    for label,sectors in errors.items():
        for name,sector in sectors.items():
            power=s.Rational(sector['mode'][2])
            for j,value in enumerate(sector['rows']):
                diff=endpoint_errors[label][name]['rows'][j]-value*s.exp(4*lam*power)
                if s.simplify(s.expand_power_exp(diff))!=0:raise ArithmeticError('Actual gap end remainder source differs')
                checks['same_'+label+'_'+name+'_ordinary_row'+str(j)]=True
    return dict(original_arbitrary_full_gap_moment_stress_splitter=generic,
        current_selected_five_primitive_velocity_pressure_trace_theorem=current,
        current_complete_pressure_and_units=field.end_tensor.units,
        source_normalization_and_full_weight_theorem=field.atlas.pulse_support.proof['exact_normalization_and_five_weight_defining_source'],
        original_full_paper_stress_AST_theorem=field.end_tensor.formulas,
        original_full_physical_decomposition_theorem=field.end_tensor.physical_proof,
        current_exact_velocity_and_remainder_scale_identities=checks,
        current_selected_parameter_algorithms_bound_to_both_sides=True,
        current_full_stress_rows_identified_by_generic_splitter_with_bound_current_histories=True,
        same_current_end_weight_cells=field.cells,
        weighted_quadrature_scope='Directed partitions enclose the same exact defining normalized beta integrals; cell coefficients are not field point values.',
        current_shared_outer_map_logscale_to_reduced_D0_AST_bridge=True,
        same_actual_end_logD_source=field.atlas.pulse_interfaces.proof['current_logE_and_transport_source_copies'],
        selected_source_functions_substituted_before_axial_differentiation=True,
        same_native_complete_future_and_full_squared_beta_source=True,
        ordinary_logR_rows_not_phase_derivatives=True,
        reciprocal_coordinate_reduced_before_enclosure=True,
        current_absolute_pressure_source_is_checked_reduced_Rv_datum=True,
        physical_tensor_trace_from_same_original_linear_pullback=True,
        input_hashes={**generic['input_hashes'],**asts.hashes},passed=True)


def canonical_tensor_groups(view):
    """Full component functions, independent of each chart's sector split."""
    result={}
    for key in TENSOR_KEYS:
        value=view[key]
        if key in TENSOR_KEYS[:2]:
            for label,sectors in value.items():
                for row in next(iter(sectors.values())):
                    result[key+'/'+label+'/'+row]=[part[row] for part in sectors.values()]
        elif key=='completed_theta_theta_stress_mixed2':
            for row in next(iter(value.values())):result[key+'/'+row]=[part[row] for part in value.values()]
        elif key=='physical_three_component_remainder_mixed2':
            for label,sectors in value.items():
                for row in next(iter(sectors.values())):result[key+'/'+label+'/'+row]=[part[row] for part in sectors.values()]
        else:
            for label,rows in value.items():result[key+'/'+label]=rows
    if len(result)!=71:raise ValueError('Complete common tensor component layout required')
    return result


class CurrentPulseGapBackgroundTensor:
    @source_precision
    def __init__(self,end_tensor=None,require_checked=True,cells=None):
        self.end_tensor=end_tensor if end_tensor is not None else CurrentPulseEndBackgroundTensor()
        if not self.end_tensor.acceptance_loaded:raise ValueError('Checked actual current pulse-end required')
        self.physical=self.end_tensor.physical;self.history=self.physical.history
        self.pulse=self.history.selected.pulse;self.ctx=self.physical.ctx;self.atlas=self.end_tensor.atlas
        self.family=self.end_tensor.family;self.source=self.end_tensor.source;self.datum_sha=self.end_tensor.datum_sha
        self.cells=self.end_tensor.cells if cells is None else cells
        if self.cells!=self.end_tensor.cells:raise ValueError('Both sides must share the current checked end defining weight cells')
        self.assert_graph();self.proof=current_gap_source_theorem(self)
        c=self.ctx;mu=c.mpf(self.pulse.mu)
        normal=c.mpf(self.pulse.pulse.initial.repair.normalization)
        self.weights=backward_bump_weights(c,mu,normal,-4,self.cells)
        self.cache={};self.hashes=dict(self.end_tensor.hashes);self.hashes.update(self.proof['input_hashes'])
        for stem in ('pulse_gap_similarity_C4','pulse_gap_physical_C2','current_pulse_gap_background_tensor'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current gap tensor receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.end_tensor.assert_graph()
        if not (self.end_tensor.acceptance_loaded and self.physical is self.end_tensor.physical and self.history is self.physical.history
                and self.pulse is self.history.selected.pulse and self.pulse is self.physical.pulse and self.pulse is self.end_tensor.pulse
                and self.atlas is self.end_tensor.atlas and self.ctx is self.physical.ctx and self.cells==self.end_tensor.cells):
            raise ValueError('Actual gap must share checked current selected pulse/complete pressure graph')

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in ('pulse_gap','pulse_gap_end'):raise ValueError('Actual gap/gap-end chart required')
        if not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(nu)[0]<=0:raise ValueError('Z[-1,1], finite coordinates/logtau and nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        mu=c.mpf(self.pulse.mu);delta=c.mpf(self.pulse.delta);pp=1+2*mu
        if chart=='pulse_gap':
            if endpoints(v)[0]<11 or endpoints(v)[1]>12:raise ValueError('Original main gap xi[11,12] required')
            d=13-v;ss=-d/mu;xi=v
        else:
            if endpoints(v)[0]<0 or endpoints(v)[1]>1:raise ValueError('Gap-end exact reduced phase[0,1] required')
            if endpoints(v)==(mp.mpf(0),mp.mpf(0)):d=c.mpf(1);ss=-1/mu
            elif endpoints(v)==(mp.mpf(1),mp.mpf(1)):d=4*mu;ss=c.mpf(-4)
            else:d=4*mu+(1-4*mu)*(1-v);ss=-d/mu
            xi=13-d
        key=endpoints(Z)
        if key not in self.cache:
            selected,inlet,u,_,_=self.pulse.data(Z);controls=[copy_jet(c,value) for value in selected['selected_scaled_end_coefficient_Taylor']]
            zero=controls[0]*0;M=[zero,zero];J0=zero
            for Cj,center in zip(controls,(-3,-1)):
                for i in (1,2):M[i-1]-=Cj*(c.exp((c.mpf('.5')-i*mu)*center)*self.weights[i-1])
                J0+=Cj*Cj*(c.exp(-2*mu*center)*self.weights[2])
            future=copy_jet(c,self.pulse.selection.future.future(Z)['complete_future_energy_Taylor'])/2
            if future.order!=5 or any(value.order!=5 for value in controls+M+[J0]):raise ValueError('Current selected full axial5 history required')
            P0=reduced_pressure_at_Rv(self.end_tensor,Z)
            self.cache[key]=dict(controls=controls,M=M,J0=J0,future=future,P0=P0,inlet=inlet,u=u)
        data=self.cache[key];z=IntervalTaylor.variable(c,Z,5)
        C=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]).reciprocal()
        rows=gap_shapes(c,mu,delta,z,C,c.mpf(self.pulse.Xp),*data['M'],data['future'],data['J0'],data['P0'],d)
        L=c.mpf(self.pulse.pulse.rows['saddle_L']);u0=c.mpf(self.pulse.pulse.rows['saddle_u0'])
        finite=-3*L+2*c.ln(u0)-c.ln(6)/2-2*c.ln(mu)
        logB=dict(logPstar=self.physical.logP,actual_log_inlet_U=c.ln(c.mpf(self.end_tensor.flatten_power.flatten.U)),
            inverse_mu=(-13+d)/(2*mu),finite=-13+d)
        # The end chart preserves its finite original offset at the exact
        # attachment, instead of recovering -4 by rounded division.
        if chart=='pulse_gap_end':logB=dict(logPstar=logB['logPstar'],actual_log_inlet_U=logB['actual_log_inlet_U'],inverse_mu=-13/(2*mu),finite=-13-(c.mpf('.5')+mu)*ss)
        logR=self.physical.logRp+13/mu+ss
        D={0:-1/mu+finite,1:(d/2-1)/mu+finite-d,2:(d/2-1)/mu+finite-2*d}
        H=-(1-mu)*(13-d)/mu;Q=-pp*d/mu;sectors={}
        for label,parts in rows['stress'].items():
            sectors[label]={}
            for name,part in parts.items():
                rp,bp,dp,hp=part['mode'];which=int(part['selected_D_recipe'][-1])
                logs=dict(source_logR=rp*logR,**{k:bp*value for k,value in logB.items()},
                    selected_log_end_scale=dp*D[which],signed_original_memory_log=hp*H,normalization=-c.ln(2)/2)
                if part['pressure_memory']:logs['same_absolute_pressure_memory']=Q
                grid={'s%d_Z%d'%(j,n):jet[n]*math.factorial(n) for j,jet in enumerate(part['full_derivative_rows']) for n in range(4-j)}
                sectors[label][name]=dict(mode=part['mode'],selected_D_recipe=part['selected_D_recipe'],
                    exact_source_log_parts=logs,full_stress_mixed3_coefficient_enclosures=grid)
        packet=dict(Z=Z,s=ss,full_meridional_stress_log_sectors=sectors,exact_pulse_reference_logB_parts=logB,
            exact_logR=logR,exact_logD=D[1],exact_logH=H)
        point=lift_physical_packet(c,packet,delta,rows['velocity'],lt,theta,nu)
        native_pressure=self.pulse.pressure_moment(Z,xi/self.pulse.mu,data['inlet'],data['u'],
            dict(inverse_mu_term=-self.pulse.prate*xi/self.pulse.mu))
        return dict(point,chart=chart,coordinate=v,distance=d,ordinary_radial_derivative='d_logR=d_s=-mu*d_distance',
            source_reduced_endpoint=('D=1,xi=12,s=-1/mu' if chart=='pulse_gap_end' and endpoints(v)==(mp.mpf(0),mp.mpf(0)) else
                'D=4mu,s=-4' if chart=='pulse_gap_end' and endpoints(v)==(mp.mpf(1),mp.mpf(1)) else None),
            current_actual_source_stress_packet=packet,source_five_full_history_rows=rows,
            current_selected_controls=data['controls'],complete_current_future=data['future'],full_selected_end_energy_loss=data['J0'],
            full_signed_absolute_Rv_pressure=data['P0'],original_native_forward_absolute_pressure=native_pressure,
            exact_source_D0_D1_D2_Q_logs={**{'D'+str(k):value for k,value in D.items()},'Q':Q},current_source_three_component_velocity_rows=rows['velocity'],
            actual_full_stress_not_local_difference=True,current_selected_complete_history_graph_retained=True,
            pressure_remaining_reduction_before_enclosure=True,effective_radial_remainder_scale_is_D1=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name not in SEAMS:raise ValueError('Current gap_coordinate or gap_end tensor seam required')
        if name=='gap_coordinate':left=self.chart('pulse_gap',Z,12,log_tau,theta,viscosity);right=self.chart('pulse_gap_end',Z,0,log_tau,theta,viscosity)
        else:left=self.chart('pulse_gap_end',Z,1,log_tau,theta,viscosity);right=self.end_tensor.end(Z,-4,log_tau,theta,viscosity)
        a=canonical_tensor_groups(left);b=canonical_tensor_groups(right);common={}
        if set(a)!=set(b):raise ValueError('Actual common tensor functions differ')
        for key in a:
            values=[endpoints(row['log_absolute_upper'])[1] for row in a[key]+b[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_source_function_tensor_trace_theorem=self.proof,source_function_equality_precedes_common_triangle_bounds=True,
            source_reduced_reciprocal_or_finite_endpoint=True,actual_nonzero_shared_boundary_histories_preserved=True,
            interval_overlap_not_used_as_function_identity=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_gap_five_history_full_tensor_and_source_join_theorem=self.proof,
            actual_current_tensor_regions_available=['pulse_gap','pulse_gap_end']+self.end_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=11,actual_current_completed_tensor_internal_interface_count=4,
            current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),source_domain='gap xi[11,12];gap_end phase[0,1] -> d[4mu,1], original s=-d/mu;R>0,|Z|<1,tau>0,constant nu>0,compact logtau',
            current_shared_end_gap_weight_cells=self.cells,original_gap_selected_full_beta_weights=self.weights,
            current_gap_actual_three_component_remainder_retained=True,remaining_pulse_core_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseGapBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_gap_gapend_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_gap_gapend_tensor_views'][name]=field.chart(*args)
        print('Current actual gap tensor: '+name,flush=True)
    result['current_actual_two_gap_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
