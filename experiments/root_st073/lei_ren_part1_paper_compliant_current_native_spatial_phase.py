"""Original affine log-radius offsets and candidate spatial periodic phase.

Selected binary logCstar/T scales admit exact modular arithmetic. Analytic
logPstar/Tw remain intervals. Microscopic offsets never enter a rounded huge
absolute logR. This binds N*log(R/r_minus), before global frequency admission.
"""
import ast
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_candidate_densities as density

phase=density.phase;current=phase.current;prior=phase.prior;packets=phase.packets;native=phase.native
HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha;ep=phase.ep
NAME=PREFIX+'current_native_spatial_phase.json'
RECEIPT=PREFIX+'current_native_spatial_phase_check.json'
GATE='current_original_native_radius_offsets_and_candidate_spatial_phase_executed'
DOMAIN=PREFIX+'current_generic_shear_loop_domain.json'
COLLAR=PREFIX+'current_inner_exit_strict_collar.json'


def fractional_record(q):return dict(numerator=q.numerator,denominator=q.denominator)


def exact_coordinate(value):
    """Only explicitly supplied exact literals; never an interval midpoint."""
    if isinstance(value,Fraction):return value
    if type(value) in (str,int):return Fraction(value)
    return None


def binary_mod_one(c,selected,coefficient):
    """Fraction of a selected binary scalar without allocating its integer."""
    lo,hi=ep(selected)
    if lo!=hi:raise ValueError('Only a source-selected singleton binary constant can use exact modulus')
    if not isinstance(coefficient,Fraction):raise ValueError('Explicit exact rational multiplier required')
    sign,man,exponent,bc=lo._mpf_;num=coefficient.numerator*man*(-1 if sign else 1);den=coefficient.denominator
    if exponent>=0:
        remainder=(num*pow(2,exponent,den))%den
    elif -exponent<=8192:
        den <<= -exponent;remainder=num%den
    else:raise ArithmeticError('Tiny binary constant needs a factored fractional representation')
    q=Fraction(remainder,den)
    return c.mpf(q.numerator)/q.denominator,dict(source_exact_mpf_tuple=[sign,man,exponent,bc],
        exact_multiplier=fractional_record(coefficient),exact_fraction=fractional_record(q),
        huge_integer_materialized=False,integer_period_retained_as_source_expression=True)


def ordinary_mod_one(c,value):
    lo,hi=ep(value);width=c.mpf(hi)-c.mpf(lo)
    if ep(width)[1]>=1:return dict(full_period=True,boxes=[c.mpf((0,1))])
    p=mp.mp.clone();p.dps=c.dps+20
    low=p.floor(p.mpf(lo));high=p.floor(p.mpf(hi))
    if low==high:
        boxes=[phase.clipped(c,value-c.mpf(low),0,1)]
    else:
        boxes=[c.mpf((ep(value-c.mpf(low))[0],1)),c.mpf((0,ep(value-c.mpf(high))[1]))]
    return dict(full_period=False,boxes=boxes,integer_period_cover=c.mpf((low,high)),
        floor_used_only_for_directed_endpoint_bracketing=True)


def source_assignment_theorem():
    """Bind selected binary parameters to original assignments, not overlaps."""
    specs={
        'physical_norm_family.py':('selected_logC','b.hi(2 * max(endpoints(v)[1] for v in restrictions.values()))'),
        'actual_long_reshape_mixed_C4.py':('self.T','400 * self.A')}
    hashes={};checks={}
    for stem,(target,expected) in specs.items():
        name=PREFIX+stem;tree=ast.parse((HERE/name).read_text(encoding='utf8'))
        found=[n.value for n in ast.walk(tree) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        if len(found)!=1 or ast.dump(found[0])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Original selected radius parameter assignment changed: '+target)
        hashes[name]=sha(name);checks[target]=True
    tree=ast.parse((HERE/(PREFIX+'physical_norm_family.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='LogBounds')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='hi')
    returned=next(n for n in fn.body if isinstance(n,ast.Return))
    if ast.dump(returned.value)!=ast.dump(ast.parse('self.c.mpf(endpoints(self.c.mpf(value))[1])',mode='eval').body):
        raise ValueError('Original upper scalar parameter selector changed')
    name=PREFIX+'pre_pulse_mixed_C4.py';tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantPrePulseMixedC4')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='axial')
    found=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='y' for t in n.targets)]
    expected=[ast.parse(v,mode='eval').body for v in ('c.exp(md * phase)','c.exp(md) + selector')]
    if sorted(ast.dump(v) for v in found)!=sorted(ast.dump(v) for v in expected):
        raise ValueError('Original axial/buffer radius coordinate assignments changed')
    hashes[name]=sha(name)
    return dict(passed=True,original_selected_scalar_assignments=checks,original_endpoint_selection_bound=True,
        original_axial_and_buffer_radius_coordinates_bound=True,
        analytic_logPstar_and_Tw_not_selected=True,input_hashes=hashes)


def affine_identity_theorem():
    P,C,T,W,f,s,sc,B,S,x,M=sy.symbols('P C T W f s sc hbB hbS x Md',real=True)
    a,b,d=sy.symbols('log4 log100 log110');Ra=a-4*P-1000;ref=d+10*(C+P);minus=Ra+B*sc/2
    actual={
        'bridge_first':Ra+B*s,'bridge_second':Ra+B*s,'bridge_macro':Ra+(b-Ra)*f+2*B*(1-f),
        'switch_first':b+S*s,'switch_second':b+S*s,'switch_power':b+(d-b)*f+2*S*(1-f),
        'reshape':d+T*f,'inner_reference':d+T+(10*(C+P)-T-8)*f,
        'axial_restore':ref-8+s,'restore_buffer':ref+s,'actual_patch':ref-6+sy.log(x),
        'Rh_reference':ref+s,'O2_slope':ref+s,'O2_axial':ref+sy.exp(M*f),
        'O2_buffer':ref+sy.exp(M)+s,'O3_slope_mu':ref+P+s,'O3_power':ref+P+1+W*f}
    rewritten={
        'bridge_first':B*(s-sc/2),'bridge_second':B*(s-sc/2),
        'bridge_macro':4*P*f+(b-a+1000)*f+B*(2*(1-f)-sc/2),
        'switch_first':4*P+b-a+1000+S*s-B*sc/2,'switch_second':4*P+b-a+1000+S*s-B*sc/2,
        'switch_power':4*P+b-a+1000+(d-b)*f+2*S*(1-f)-B*sc/2,
        'reshape':4*P+d-a+1000+T*f-B*sc/2,
        'inner_reference':(4+10*f)*P+10*f*C+(1-f)*T+d-a+1000-8*f-B*sc/2,
        'axial_restore':14*P+10*C+d-a+992+s-B*sc/2,
        'restore_buffer':14*P+10*C+d-a+1000+s-B*sc/2,
        'actual_patch':14*P+10*C+d-a+994+sy.log(x)-B*sc/2,
        'Rh_reference':14*P+10*C+d-a+1000+s-B*sc/2,
        'O2_slope':14*P+10*C+d-a+1000+s-B*sc/2,
        'O2_axial':14*P+10*C+d-a+1000+sy.exp(M*f)-B*sc/2,
        'O2_buffer':14*P+10*C+d-a+1000+sy.exp(M)+s-B*sc/2,
        'O3_slope_mu':15*P+10*C+d-a+1000+s-B*sc/2,
        'O3_power':15*P+10*C+d-a+1001+W*f-B*sc/2}
    jacobian=dict(bridge_first=B,bridge_second=B,bridge_macro=4*P+b-a+1000-2*B,
        switch_first=S,switch_second=S,switch_power=d-b-2*S,reshape=T,
        inner_reference=10*(C+P)-T-8,axial_restore=1,restore_buffer=1,actual_patch=1/x,
        Rh_reference=1,O2_slope=1,O2_axial=M*sy.exp(M*f),O2_buffer=1,O3_slope_mu=1,O3_power=W)
    coordinate={key:f if key in ('bridge_macro','switch_power','reshape','inner_reference','O2_axial','O3_power')
        else x if key=='actual_patch' else s for key in actual}
    for key,value in actual.items():
        if sy.expand(value-minus-rewritten[key])!=0:raise ArithmeticError('Exact source radius cancellation failed: '+key)
        if sy.simplify(sy.diff(value,coordinate[key])-jacobian[key])!=0:
            raise ArithmeticError('Original coordinate Jacobian identity failed: '+key)
    seams=[('bridge_first',s,1,'bridge_second',s,1),('bridge_second',s,2,'bridge_macro',f,0),
        ('bridge_macro',f,1,'switch_first',s,0),('switch_first',s,1,'switch_second',s,1),
        ('switch_second',s,2,'switch_power',f,0),('switch_power',f,1,'reshape',f,0),
        ('reshape',f,1,'inner_reference',f,0),('inner_reference',f,1,'axial_restore',s,0),
        ('axial_restore',s,1,'restore_buffer',s,-7),('restore_buffer',s,-6,'actual_patch',x,1),
        ('actual_patch',x,sy.E,'Rh_reference',s,-5),('Rh_reference',s,0,'O2_slope',s,0),
        ('O2_slope',s,1,'O2_axial',f,0),('O2_axial',f,1,'O2_buffer',s,0),
        ('O2_buffer',s,11,'O3_slope_mu',s,0),('O3_slope_mu',s,1,'O3_power',f,0)]
    checked_seams=[]
    for left,lc,lv,right,rc,rv in seams:
        delta=rewritten[left].subs(lc,lv)-rewritten[right].subs(rc,rv)
        if sy.simplify(delta.subs(P,sy.exp(M)+11))!=0:
            raise ArithmeticError('Original radius/periodic phase seam failed: '+left+' -> '+right)
        checked_seams.append(left+' -> '+right)
    return dict(passed=True,exact_source_radius_minus_same_left_radius_identities=list(actual),
        original_native_coordinate_Jacobian_identities=list(jacobian),
        original_same_radius_periodic_phase_seam_identities=checked_seams,
        seams_do_not_admit_velocity_or_stress_derivative_matching=True,
        same_logRa_cancelled_before_numeric_arithmetic=True,
        microscopic_width_derivative_not_reapplied_to_native_rows=True)


def spatial_candidate_record(candidate,geometry):
    """Annotate a kernel supplied with the derived phase, without changing it."""
    result=dict(candidate['record']);result['derived_spatial_phase_box']=result.pop('explicit_free_phase_parameter')
    loop=dict(result['original_phase_inverse_and_primitives'])
    loop['free_phase_parameter_not_spatial_phase']=False
    loop['original_radius_phase_bound_at_candidate_N']=True
    result.update(original_phase_inverse_and_primitives=loop,spatial_phase_binding_installed=True,
        original_radius_phase_source=geometry,phase_box_is_derived_not_independently_selected=True)
    return result


def spatial_candidate_functions(owner,binder,chart,Z,coordinate,N):
    if type(owner) is not density.NativeCandidateDensities:raise ValueError('Existing native candidate density owner required')
    if owner.family!=binder.family or owner.owner.owner.owner.native.seed is not binder.seed:
        raise ValueError('Same original native density/radius seed and family required')
    geometry=binder.query(chart,Z,coordinate,N)
    source,loop,V=owner.source(chart,Z,geometry['raw']['coordinate']);records=[]
    for box in geometry['phase_boxes']:
        candidate=density.candidate_at_phase(loop,V,box,N)
        if candidate['values'] is None:
            records.append(dict(status=candidate['record']['status'],spatial_phase_binding_installed=True,
                derived_spatial_phase_box=box,original_radius_phase_source=geometry['record']))
        else:records.append(spatial_candidate_record(candidate,geometry['record']))
    return dict(original_native_source=source['record'],original_normalized_V=V.record(),
        actual_source_radius_and_phase=geometry['record'],actual_spatial_candidate_cells=records,
        explicit_candidate_N_not_global_admission=True)


class NativeSpatialPhase:
    def __init__(self,backend):
        if type(backend) is not native.NativeGenericSourcePackets:raise ValueError('Same original17 native source backend required')
        self.native=backend;self.ctx=c=backend.ctx;self.seed=backend.seed;self.family=backend.family;self.service=backend.service
        if self.seed.radius.__func__ is not prior.CompliantGlobalPhysicalAssembly.radius:raise ValueError('Unchanged original radius operator required')
        self.theorem=source_assignment_theorem();self.service.bind_hashes(self.theorem['input_hashes'])
        for stem,gate in (('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'),
                          ('current_inner_exit_strict_collar','current_inner_exit_strict_collar_attached_to_current_source_graph_certified')):
            name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
            family=receipt.get('source_family') or {k:receipt[k] for k in packets.FAMILY_KEYS}
            if not receipt.get('all_passed') or not receipt.get(gate) or family!=self.family:raise ValueError('Same original radius/domain/collar source required')
            self.service.bind_hashes(receipt['input_hashes']);self.service.bind_hashes({name:sha(name)})
        left=json.loads((HERE/COLLAR).read_bytes())['explicit_current_inner_exit_strict_collar']
        self.sc=packets.interval(c,left['selected_first_phase_endpoint']);self.loghB=c.mpf(backend.bridge.logh);self.loghS=c.mpf(backend.switch.logh)
        self.fixed=dict(logP=self.seed.logP,logC=self.seed.core.logC,T=backend.reshape.reshape.T,Tw=self.seed.params.Tw)
        for name in ('logC','T'):
            if ep(self.fixed[name])[0]!=ep(self.fixed[name])[1]:raise ValueError('Original selected singleton '+name+' required for modular arithmetic')
        selected=self.seed.core.records['physical_norm_family']
        if ep(packets.interval(c,selected['selected_logCstar']))!=ep(self.fixed['logC']):raise ValueError('Same original selected logCstar required')
        if ep(packets.interval(c,selected['A_upper'])*400)!=ep(self.fixed['T']):raise ValueError('Same original selected400*Abar scale required')
        self.bases=tuple(c.mpf(0) for _ in range(5));self.ledger=dict(directed_small_exponential_tails=0,
            positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        self.service.bind_hashes({name:sha(name) for name in (DOMAIN,COLLAR,PREFIX+'global_physical_assembly.py',Path(__file__).name)})
        self.identity=affine_identity_theorem()

    def scalar(self,value):return prior.ScaledEnclosure(prior.FormalScale(self.bases),value,self.ledger)
    def width(self,kind,coefficient):
        return prior.ScaledEnclosure(prior.FormalScale(self.bases,offset=self.loghB if kind=='bridge' else self.loghS),coefficient,self.ledger)

    def query(self,chart,Z,coordinate,N):
        N=density.candidate_integer(N);c=self.ctx;tag=exact_coordinate(coordinate);source_expression=None
        if isinstance(coordinate,dict):
            if set(coordinate)=={'original_power_offset'} and chart=='O3_power':
                offset=Fraction(coordinate['original_power_offset']);v=(c.mpf(offset.numerator)/offset.denominator)/self.fixed['Tw'];source_expression='phase=explicit original offset / same Tw'
            elif set(coordinate)=={'selected_sc_multiple'} and chart=='bridge_first':
                mult=Fraction(coordinate['selected_sc_multiple']);v=self.sc*(c.mpf(mult.numerator)/mult.denominator);source_expression='phase=explicit multiple of same selected s_c'
            else:raise ValueError('Unsupported original source coordinate expression')
        elif tag is not None:v=c.mpf(tag.numerator)/tag.denominator
        else:v=c.mpf(coordinate)
        raw=self.native.original_raw(chart,Z,v);packet=raw['source_packet']
        # The original radius descriptor is retained as source evidence. Its
        # cap-based numerical logR is deliberately never used in arithmetic.
        _,descriptor=self.seed.radius(chart,v,packet,raw['native_provider'])
        if chart.startswith('bridge_') and ep(packet['source_width_log'])!=ep(self.loghB):raise ValueError('Same actual bridge width required')
        if chart.startswith('switch_') and ep(packet['exact_positive_width_log'])!=ep(self.loghS):raise ValueError('Same actual switch width required')
        cv=lambda q:c.mpf(q.numerator)/q.denominator if isinstance(q,Fraction) else c.mpf(q)
        f=tag if tag is not None else v
        coeff={};constant=c.mpf(0);micro={'bridge':-self.sc/2};jac_regular=c.mpf(0);jac_micro={}
        L100=c.ln(100)-c.ln(4)+1000;L110=c.ln(110)-c.ln(4)+1000
        if chart in ('bridge_first','bridge_second'):
            micro['bridge']=v-self.sc/2;jac_micro['bridge']=c.mpf(1)
            if source_expression and 'selected_sc_multiple' in coordinate:micro['bridge']=self.sc*cv(mult-Fraction(1,2))
        elif chart=='bridge_macro':
            coeff['logP']=4*f;constant=L100*cv(f);micro['bridge']=2*(1-v)-self.sc/2
            jac_regular=4*self.fixed['logP']+L100;jac_micro['bridge']=c.mpf(-2)
        elif chart in ('switch_first','switch_second'):
            coeff['logP']=Fraction(4);constant=L100;micro['switch']=v;jac_micro['switch']=c.mpf(1)
        elif chart=='switch_power':
            coeff['logP']=Fraction(4);constant=L100+(c.ln(110)-c.ln(100))*v
            micro['switch']=2*(1-v);jac_regular=c.ln(110)-c.ln(100);jac_micro['switch']=c.mpf(-2)
        elif chart=='reshape':
            coeff={'logP':Fraction(4),'T':f};constant=L110;jac_regular=self.fixed['T']
        elif chart=='inner_reference':
            coeff={'logP':4+10*f,'logC':10*f,'T':1-f};constant=L110-8*v
            jac_regular=10*(self.fixed['logC']+self.fixed['logP'])-self.fixed['T']-8
        else:
            coeff={'logP':Fraction(14),'logC':Fraction(10)};constant=L110;jac_regular=c.mpf(1)
            if chart=='axial_restore':constant+=v-8
            elif chart in ('restore_buffer','Rh_reference','O2_slope'):constant+=v
            elif chart=='actual_patch':constant+=c.ln(v)-6;jac_regular=1/v
            elif chart in ('O2_axial','O2_buffer'):
                constant+=packet['actual_y'];jac_regular=packet['actual_y']*c.mpf(self.seed.params.Md) if chart=='O2_axial' else c.mpf(1)
            elif chart=='O3_slope_mu':coeff['logP']=Fraction(15);constant+=v
            elif chart=='O3_power':
                coeff['logP']=Fraction(15);constant+=1;jac_regular=self.fixed['Tw']
                if source_expression:constant+=cv(offset)
                else:coeff['Tw']=f
            else:raise ValueError('Original modification chart required')
        regular=constant+sum((self.fixed[k]*cv(vv) for k,vv in coeff.items()),c.mpf(0))
        pieces={k:self.width(k,x) for k,x in micro.items()};offset=sum(pieces.values(),self.scalar(regular))
        jac=sum((self.width(k,x) for k,x in jac_micro.items()),self.scalar(jac_regular))
        residual=c.mpf(0);component_proofs=[];full=False
        for key,value in [('constant',constant),*[(k,self.fixed[k]) for k in coeff]]:
            multiplier=Fraction(N) if key=='constant' else coeff[key]*N
            if key in ('logC','T') and isinstance(multiplier,Fraction):
                box,proof=binary_mod_one(c,value,multiplier);residual+=box;component_proofs.append(dict(component=key,method='exact_selected_binary_modulus',**proof))
            else:
                term=value*cv(multiplier);wrapped=ordinary_mod_one(c,term)
                component_proofs.append(dict(component=key,method='directed_analytic_interval_modulus',**wrapped))
                if wrapped['full_period']:full=True;break
                # Multiple cells are a genuine union; a hull is sufficient for
                # the next modulo projection and does not select a cell.
                residual+=c.mpf((min(ep(q)[0] for q in wrapped['boxes']),max(ep(q)[1] for q in wrapped['boxes'])))
        tiny_parts={k:x*N for k,x in pieces.items()};tiny=sum(tiny_parts.values(),self.scalar(0))
        if full:projected=dict(full_period=True,boxes=[c.mpf((0,1))])
        else:
            try:projected=ordinary_mod_one(c,residual+phase.bounded_value(tiny))
            except ArithmeticError:projected=dict(full_period=True,boxes=[c.mpf((0,1))],requires_source_phase_refinement=True)
        record=dict(chart=chart,source_family=self.family,Z_box=raw['Z'],coordinate_box=raw['coordinate'],
            declared_exact_coordinate=None if tag is None else fractional_record(tag),original_coordinate_expression=source_expression,
            original_radius_source_descriptor=descriptor,source_radius_caps_not_consumed=True,
            original_same_Ra_and_reference_logs_collected_before_evaluation=True,
            exact_regular_source_coefficients={k:fractional_record(x) if isinstance(x,Fraction) else x for k,x in coeff.items()},
            regular_constant=constant,regular_log_radius_offset_cover=regular,
            microscopic_log_radius_offset_components={k:x.record() for k,x in pieces.items()},
            actual_log_R_over_same_r_minus_enclosure=offset.record(),
            original_dy_dnative_coordinate_enclosure=jac.record(),native_width_conversion_not_reapplied=True,
            candidate_N=N,phase_component_proofs=component_proofs,
            microscopic_N_y_components={k:x.record() for k,x in tiny_parts.items()},
            microscopic_N_y_source_sum=tiny.record(),periodic_projection=projected,
            original_spatial_phase_bound_at_candidate_N=True,phase_independent_of_Z=True,
            selected_binary_constants_not_midpoints=True,analytic_logPstar_and_Tw_not_selected=True,
            integer_cycles_not_allocated_as_huge_integers=True,
            global_common_N_or_actual_integrals_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,phase_boxes=projected['boxes'],offset=offset,jacobian=jac,raw=raw)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeSpatialPhase(native.NativeGenericSourcePackets(bridge));N=1024
        points=dict(bridge_first='.1337',bridge_second='1.831',bridge_macro='.537',switch_first='.537',
            switch_second='1.337',switch_power='.537',reshape='.537',inner_reference='.1337',axial_restore='.537',
            restore_buffer='-6.337',actual_patch='1.337',Rh_reference='-1.337',O2_slope='.1337',O2_axial='.1337',
            O2_buffer='5.337',O3_slope_mu='.537',O3_power={'original_power_offset':'.537'})
        records={}
        for chart,coordinate in points.items():
            result=owner.query(chart,('.5','.5'),coordinate,N);records[chart]=dict(query_coordinate=coordinate,result=result['record'])
            print('Original spatial phase:',chart,'full period' if result['record']['periodic_projection']['full_period'] else 'bounded cells',flush=True)
        inlet=owner.query('bridge_first',(-1,1),{'selected_sc_multiple':'1/2'},N)
        candidate_owner=density.NativeCandidateDensities(phase.NativeConditionedPhase(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(owner.native))))
        spatial={}
        for chart in ('inner_reference','O2_slope','O2_buffer','O3_slope_mu','O3_power'):
            spatial[chart]=spatial_candidate_functions(candidate_owner,owner,chart,('.5','.5'),points[chart],N)
            print('Actual spatial candidate velocity/densities:',chart,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_chart_count=len(records),candidate_N=N,
        current_original_radius_and_spatial_phase_records=records,actual_same_left_inlet=inlet['record'],
        actual_spatial_candidate_velocity_and_density_records=spatial,
        original_selected_parameter_source_theorem=owner.theorem,original_affine_radius_identity_theorem=owner.identity,
        spatial_phase_is_actual_N_log_R_over_r_minus=True,actual_spatial_velocity_density_functions_installed=True,
        actual_spatial_velocity_density_scope='five explicitly declared source-coordinate boxes at candidate N=1024',
        actual_changed_defect_integral_functions_installed=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual original radius offsets/Jacobians and candidate N*y periodic enclosures on17 native charts plus real inlet, feeding spatial E_N/V_N and five signed density functions on five exact coordinate requests. Selected binary constants reduced exactly; analytic constants remain intervals. No global frequency/whole inverse coverage/integrated histories/repair/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
