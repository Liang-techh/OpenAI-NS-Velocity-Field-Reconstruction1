"""Current waiting -> original collar -> complete Gamma exterior sources.

Only acquisition changes. Exact Gamma, infinite tails and original forward
angular/absolute-pressure histories are preserved; caps remain enclosures.
"""
import ast
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_steep_waiting_source import (
    CurrentSteepWaitingSourceAssembly, CHARTS as PRIOR_CHARTS,
    UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4 as BASE
from lei_ren_part1_paper_compliant_collar_Gamma_C4_check import functional_source_identities
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

METHODS={'heat_collar':'collar','heat_exterior':'exterior'}
NEW_CHARTS=tuple(METHODS)
CHARTS=PRIOR_CHARTS+NEW_CHARTS
RECEIPT=PREFIX+'current_heat_source_check.json'
GATE='current_heat_source_ownership_certified'
PROVED='current_waiting_collar_and_exact_Gamma_functional_joins_proved'
THEOREM=PREFIX+'collar_Gamma_C4_check.json'
GAMMA_THEOREM=PREFIX+'heat_pressure_C4_check.json'
VIEWS={'heat_collar':{'inlet':0,'whole_domain':[0,3],'exit':3,'phi_crossing':['2.99','3']},
       'heat_exterior':{'inlet':3,'whole_domain':[3,mp.inf]}}


class CurrentCollarGammaC4(BASE):
    @source_precision
    def __init__(self,source,cells=64):
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive directed cell count required')
        self.steep=source.steep;self.outer=self.steep.outer;self.future=self.steep.future
        # Do not name this alias `heat`: the original physical radius mapper
        # uses that name for a companion's collar provider, not the raw datum.
        self.exact_heat=self.future.heat;self.repair=self.future.repair
        self.ctx=c=self.steep.ctx;self.mu=self.steep.mu;self.delta=self.steep.delta;self.a=self.delta/2
        self.eps=self.steep.epsilon;self.k=self.steep.k;self.bh=self.steep.bh;self.prate=1+self.delta
        self.family=source.family;self.source=source.source;self.Scap=self.steep.S
        self.S=c.mpf([0,endpoints(self.Scap)[1]]);self.cells=cells;self.hashes=dict(source.hashes)
        self.cache={};self.shape_cache={};self.tail_cache={};self.gamma_cache={}
        self.theta_base=self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)
        self.Ev2=self.outer.flatten.Ev2;self.pressure_scale=self.Ev2*self.theta_base**2

    def packet(self,*args,**kwargs):
        packet=BASE.packet(self,*args,**kwargs)
        packet.update(full_pulse_C4_installed=False,current_live_waiting_heat_source_used=True,
            source_caps_used_as_defining_field_values=False,
            current_nonzero_angular_and_absolute_pressure_histories_retained=True)
        return packet


def same_exact_Gamma_future_bindings():
    """Bind production integrands before different directed cell enclosures."""
    bindings={
        'current_prefix_gamma_call':assignment_source_bindings('fifth_axial_jets','angular',{
            'rawheat':'heat_future_jets(c,Z,box(r.delta)/2,box(r.heat.epsilon),box(r.strong_S_cap),5,self.angular4.cells)'}),
        'prefix_collar_source':assignment_source_bindings('angular_high_jets','heat_future_jets',{
            'D':'gamma_deficit_jets(c,Z,a,S_cap,order,t)*(sig*C)',
            'squaredifference':'D*(D*(-a*Sbox)+2*pre)',
            'tails':'integrated_gamma_tails(c,Z,a,S_cap,order)'}),
        'current_collar_source':assignment_source_bindings('collar_Gamma_C4','collar_tails',{
            'shape':'self.shape(Z,v,False)',"D":"shape['D_rows'][0]",'square':'D*(pre*2-D*(self.a*self.S))',
            'tails':'integrated_gamma_tails(c,Z,self.a,self.Scap,5,3)'}),
        'current_terminal_source':assignment_source_bindings('collar_Gamma_C4','data',{
            'terminal':'self.steep.waiting(Z,1)','heat0':'self.collar_tails(Z,0)',
            'Xtail':"terminal['angular_Taylor']",'Ptail':"terminal['pressure_over_Pstar_squared_Taylor']",
            'defect':"Xtail*(1-self.eps)-heat0['angular_numerator']"})}
    # Both helpers import the identical infinite-tail callable. Their Gamma
    # point/derivative enclosures are retained in GAMMA_THEOREM, not assumed
    # equal as boxes and not replaced by an S polynomial or cap endpoint.
    parsed=ast.parse((HERE/(PREFIX+'collar_Gamma_C4.py')).read_text(encoding='utf8'))
    if not any(isinstance(n,ast.ImportFrom) and n.module=='lei_ren_part1_paper_compliant_angular_high_jets'
            and any(v.name=='integrated_gamma_tails' for v in n.names) for n in parsed.body):
        raise ValueError('Actual common infinite Gamma tail import changed')
    syntax={
        ('angular_high_jets','heat_future_jets','pre'):'1-eps*(1-sig+sig*phi_box)',
        ('angular_high_jets','heat_future_jets','C'):'1-eps*phi_box',
        ('angular_high_jets','heat_future_jets','sig'):'stable_sigma(c,t)[0]'}
    for (stem,method,target),expression in syntax.items():
        assignment_source_bindings(stem,method,{target:expression})
    # Bind actual flat phi C0: sigma source equivalence is already admitted
    # by current O7. The two written phi formulas are exactly the same.
    parsed=ast.parse((HERE/(PREFIX+'angular_high_jets.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(parsed) if isinstance(n,ast.FunctionDef) and n.name=='phi')
    expected='c.exp(-4/(3-c.mpf(v))**2) if v<3 else c.mpf(0)'
    if len(fn.body)!=1 or ast.dump(fn.body[0].value)!=ast.dump(ast.parse(expected,mode='eval').body):
        raise ValueError('Actual future phi definition changed')
    assignment_source_bindings('collar_Gamma_C4','phi_jets',{
        'log_upper':'-4/D**2','dist':'IntervalTaylor(c,[3-t,-1,0,0,0])','jets':'(dist**(-2)*(-4)).exp()'})
    proofs={}
    def zero(name,value):
        if s.simplify(value)!=0:raise ValueError('Shared Gamma future identity differs: '+name)
        proofs[name]=True
    a,S,eps,sig,phi,dt,v=s.symbols('a S eps sig phi dt v',real=True)
    H=s.Function('H')(v);D=sig*(1-eps*phi)*(1-H)/(a*S);W=1-sig+sig*phi;pre=1-eps*W
    square=D*(2*pre-a*S*D)
    env={'D':D,'pre':pre,'a':a,'Sbox':S}
    zero('actual_prefix_squared_difference',assignment('compliant_angular_high_jets','heat_future_jets','squaredifference',env)-square)
    env={'D':D,'pre':pre,'self.a':a,'self.S':S}
    zero('actual_collar_squared_difference',assignment('compliant_collar_Gamma_C4','collar_tails','square',env)-square)
    for label,target in (('theta',a*D*s.exp((1-a)*v)),('pressure',square*s.exp(-(1+2*a)*v)/2),
            ('energy',square*s.exp(-2*a*v))):
        prefix=assignment('compliant_angular_high_jets','heat_future_jets',label,{
            'D':D,'squaredifference':square,'dt':dt,'a':a,'t':v},augmented=True)
        collar=assignment('compliant_collar_Gamma_C4','collar_tails','angular' if label=='theta' else label,{
            'D':D,'square':square,'ds':dt,'self.a':a,'self.k':1-a,'self.prate':1+2*a,'self.delta':2*a,'v':v},augmented=True)
        zero('actual_prefix_'+label+'_exact_integrand',prefix-dt*target)
        zero('actual_collar_'+label+'_exact_integrand',collar-dt*target)
    # Exact integral additivity, for the admitted full Gamma (including its
    # error envelopes): int_0^infty = int_0^3 + int_3^infty, never a cutoff.
    F=s.Function('F');integral=F(s.oo)-F(0)
    zero('complete_tail_additivity',(F(3)-F(0))+(F(s.oo)-F(3))-integral)
    EW,EW2,EP,EP2,J,hatE,hatP,hatA,k,prate=s.symbols('EW EW2 EP EP2 J hatE hatP hatA k prate',real=True)
    collarenv={'one':s.Integer(1),'self.k':k,'self.delta':2*a,'self.prate':prate,'t':s.Integer(0),
        'self.eps':eps,'self.S':S,'self.a':a,'angular':hatA,'energy':hatE,'pressure':hatP,
        "atoms['JW']":J,"atoms['EW']":EW,"atoms['EW2']":EW2,"atoms['PW']":EP,"atoms['PW2']":EP2}
    expected={'A':1/k+S*hatA+eps*J,'E':1/(2*a)-2*eps*EW+eps**2*EW2-a*S*hatE,
        'P':1/(2*prate)-eps*EP+eps**2*EP2/2-a*S*hatP}
    for name,value in expected.items():
        zero('actual_collar0_'+name+'_full_future_units',assignment('compliant_collar_Gamma_C4','collar_tails',name,collarenv)-value)
    e0=expected['E'];lone=s.symbols('logone',real=True)
    preheat=assignment('compliant_current_steep_waiting_source','data','preheat_tail',{
        'IntervalTaylor.constant(c, self.preheat_scalar, 5)':1/(2*a)-2*eps*EW+eps**2*EW2,'heat':hatE,
        'c.mpf([0, endpoints(self.delta * self.S / 2)[1]])':a*S})
    zero('current_waiting_and_collar0_same_complete_energy',preheat-e0)
    actualH=assignment('compliant_current_steep_waiting_source','data','H',{
        'preheat_tail':preheat,'self.tail_normalization':s.exp(-2*lone)})
    zero('current_waiting_inlet_energy_epsilon_units',actualH.subs(lone,s.log(1-eps))-e0/(1-eps)**2)
    return dict(actual_source_assignments=bindings,exact_full_Gamma_future_identities=proofs,
        same_original_sigma_source_from_current_O7=True,same_original_flat_phi_definition=True,
        same_canonical_positive_Gamma_expectation_and_derivative_enclosures=True,
        same_infinite_Gamma_tail_callable=True,different_cells_not_assumed_equal_boxes=True,
        all_three_scaled_Gamma_future_integrands_bound=True,
        both_energy_epsilon_atoms_retained=True,interval_overlap_used_as_join_proof=False,passed=True)


def current_source_bindings(owner):
    p=owner.heat;steep=owner.before.steep;f=p.future
    graph=dict(checked_current_waiting=owner.before.acceptance_loaded,
        same_current_waiting_object=p.steep is steep,same_current_outer=p.outer is steep.outer,
        same_actual_C4_C5_prefix_future=f is p.outer.fifth.fourth.energy.base,
        same_future_heat_and_repair=p.exact_heat is f.repair.heat,
        same_native_context=p.ctx is steep.ctx,same_native_mu=p.mu is steep.mu,same_native_delta=p.delta is steep.delta,
        original_shape_callable=p.shape.__func__ is BASE.shape,original_local_Gamma_callable=p.local_Gamma.__func__ is BASE.local_Gamma,
        original_full_collar_tails_callable=p.collar_tails.__func__ is BASE.collar_tails,
        original_live_terminal_callable=p.data.__func__ is BASE.data,
        original_forward_pressure_callable=p.forward_pressure.__func__ is BASE.forward_pressure,
        original_collar_callable=p.collar.__func__ is BASE.collar,original_exterior_callable=p.exterior.__func__ is BASE.exterior,
        current_prefix_parameter_bridge=owner.before.before.bindings['current_native_parameter_source_bridge']['passed'],
        current_C4_C5_prefix_bridge=owner.before.before.bindings['current_prefix_and_future_decomposition']['passed'])
    if not all(graph.values()):raise ValueError('Current heat graph or original algorithms differ')
    ctor={target:class_assignment('current_heat_source','CurrentCollarGammaC4','__init__',target,expression)
        for target,expression in {'self.steep':'source.steep','self.outer':'self.steep.outer','self.future':'self.steep.future',
            'self.exact_heat':'self.future.heat','self.repair':'self.future.repair','self.mu':'self.steep.mu',
            'self.delta':'self.steep.delta','self.a':'self.delta/2','self.eps':'self.steep.epsilon',
            'self.k':'self.steep.k','self.bh':'self.steep.bh','self.prate':'1+self.delta',
            'self.Scap':'self.steep.S','self.S':'c.mpf([0,endpoints(self.Scap)[1]])',
            'self.theta_base':'self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)',
            'self.Ev2':'self.outer.flatten.Ev2','self.pressure_scale':'self.Ev2*self.theta_base**2'}.items()}
    return dict(current_defining_object_graph=graph,current_constructor_source_bindings=ctor,
        shared_exact_Gamma_future_binding=same_exact_Gamma_future_bindings(),
        exact_inverse_radius_source='S=1/Rtail; current raw future/repaired datum, never Scap',
        retained_canonical_Gamma_enclosure_theorem=GAMMA_THEOREM,
        canonical_Gamma_enclosure_evidence_reused=owner.gamma_evidence,
        current_waiting_energy_and_heat0_same_full_function=True,
        current_angular_tail_constant_defect_retained=True,current_absolute_pressure_inlet_not_reset=True,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,passed=True)


class CurrentHeatSourceAssembly:
    @source_precision
    def __init__(self,require_checked=True):
        self.before=CurrentSteepWaitingSourceAssembly()
        self.family=self.before.family;self.source=self.before.source;self.datum_sha=self.before.datum_sha
        self.heat=CurrentCollarGammaC4(self.before);self.hashes=dict(self.heat.hashes);self.registry=dict(self.before.registry)
        theorem=accepted(THEOREM,self.family,self.source,'waiting_collar_and_collar_Gamma_joins_certified')
        self.functional_proof=functional_source_identities()
        if self.functional_proof!=theorem['functional_production_source_identities'] or not all(self.functional_proof.values()):
            raise ValueError('Canonical collar/Gamma theorem differs')
        gamma=accepted(GAMMA_THEOREM,self.family,self.source,'absolute_pressure_same_source_mixed4_available')
        self.gamma_evidence=gamma['defining_source_bridge']['defining_function_bridge']['gamma_derivative_enclosure']
        if not self.gamma_evidence['verified']:raise ValueError('Canonical full Gamma enclosure theorem omitted')
        for record in (theorem,gamma):
            for name,digest in record['input_hashes'].items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current heat source conflict: '+name)
                self.hashes[name]=digest
        self.hashes[THEOREM]=sha(THEOREM);self.hashes[GAMMA_THEOREM]=sha(GAMMA_THEOREM)
        self.bindings=current_source_bindings(self);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.registry.update({chart:dict(provider=PREFIX+'current_heat_source.CurrentCollarGammaC4',
            method=METHODS[chart],coverage_coordinate='t=log(R/Rtail)',domain='[0,3]' if chart=='heat_collar' else '[3,infinity)',
            acceptance_receipt=RECEIPT) for chart in NEW_CHARTS})
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current heat datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def provider(self,chart):return self.heat if chart in NEW_CHARTS else self.before.provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart not in NEW_CHARTS:return self.before.evaluate(chart,Z,coordinate)
        packet=getattr(self.heat,METHODS[chart])(Z,coordinate)
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,source_provider=self.registry[chart]['provider'],
            acceptance_receipt=RECEIPT,source_packet=packet,
            physical_mixed_grids={'physical_mixed_derivatives_total_order_le4':packet['physical_mixed_derivatives_total_order_le4']},
            **{GATE:self.acceptance_loaded,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
            current_heat_physical_owner_installed=False,
            output_kind='current collar/full Gamma derivative source enclosures; no production point selection',
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_downstream_chart_owner_count=29,
            ordered_current_chart_registry=self.registry,current_heat_source_bindings=self.bindings,
            retained_actual_canonical_functional_identities=self.functional_proof,
            retained_current_steep_waiting_source_certified=True,
            **{GATE:self.acceptance_loaded,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
            current_heat_physical_owner_installed=False,all_profile_source_charts_callable=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();packets={}
        for chart,views in VIEWS.items():
            packets[chart]={}
            for name,value in views.items():
                packets[chart][name]=self.evaluate(chart,[-1,1],value)
                print('Current heat source: '+chart+' '+name,flush=True)
        result.update(whole_current_heat_source_maps=packets,
            current_waiting_terminal=self.before.evaluate('waiting',[-1,1],1),input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentHeatSourceAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current heat generated: twenty-nine profile source owners, entire Gamma exterior',flush=True)
    return result


if __name__=='__main__':run()
