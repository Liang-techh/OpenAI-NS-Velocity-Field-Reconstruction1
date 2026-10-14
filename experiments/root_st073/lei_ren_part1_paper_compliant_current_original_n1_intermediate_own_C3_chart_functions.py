"""Actual intermediate n1 own C3 fields and partial pulse pressure functions.

Reuse regional leading pressure primitives, full signed histories and the
checked reserved C3 repair. This installs regional functions, not numerical
norms or the unresolved global leading-field joins/native C5 inlet bridge.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_repaired_own_C3_terminal_functions as previous
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as pressure_source
import lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions as modified_source
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions as patch_source
import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as outer_source
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions as pulse_source

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode,Algebra=previous.FunctionJet,previous.rows,previous.encode,previous.Algebra
MOMENTS=previous.MOMENTS
NAME=PREFIX+'current_original_n1_intermediate_own_C3_chart_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_intermediate_own_C3_chart_functions_check.json'
GATES=('current_original_n1_lower_and_17_preRc_own_C3_QV_C2_chart_functions_installed',
       'current_original_n1_quiet_and_pulse_own_pressure_C3_chart_functions_installed',
       'current_original_n1_actual_partial_pulse_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_intermediate_chart_numeric_containment_installed',
    'current_original_n1_global_field_join_C3_installed')))


def endpoint_functions(g,U,Q,Qy,R):
    a=Algebra(g,3);U,Q,Qy=(a.trim(q) for q in (U,Q,Qy));r=a.fixed(R);r2=a.mul(r,r)
    return dict(R_q=a.mul(r,Q),R_squared_q=a.mul(r2,Q),R_q_squared=a.mul(r,a.mul(Q,Q)),
        R_squared_q_squared=a.mul(r2,a.mul(Q,Q)),R_U_squared=a.mul(r,a.mul(U,U)),
        R_squared_U_squared=a.mul(r2,a.mul(U,U)),q_plus_q_y=a.add(Q,Qy),R_q_y=a.mul(r,Qy))


def pressure_functions(g,delta,integrals,lower,upper):
    """Unchanged original signed reduced-S formula, ordinary Z0..3."""
    a=Algebra(g,3);z=FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero));one=a.fixed(1)
    c=a.fixed(delta);d=a.sub(one,a.mul(z,z));h=a.sub(one,c);L=a.sub(one,a.mul(c,a.mul(z,z)))
    L2=a.mul(L,L);half=lambda q:a.scale(q,'1/2');bd=lambda k:a.sub(upper[k],lower[k])
    mul=lambda *v:_product(a,v);div=lambda q:a.div(q,L,'same original positive leading L')
    div2=lambda q:a.div(q,L2,'same original positive leading L squared');i=integrals
    S0=a.add(div(bd('R_q')),half(div(mul(h,z,i['I0_q_Z']))),half(bd('R_q_squared')),
        div(mul(d,i['I0_U_q_Z'])),a.scale(div2(mul(z,z,bd('R_U_squared'))),-2),
        a.scale(div2(mul(c,z,z,i['I0_U_squared'])),-2),a.scale(div2(mul(z,d,i['I0_U_U_Z'])),2))
    S1=a.add(div(bd('R_squared_q')),a.neg(div(i['I1_q'])),half(div(mul(h,z,i['I1_q_Z']))),
        half(bd('R_squared_q_squared')),a.neg(half(i['I1_q_squared'])),div(mul(d,i['I1_U_q_Z'])),
        a.scale(div2(mul(z,z,bd('R_squared_U_squared'))),-2),
        a.scale(div2(mul(h,z,z,i['I1_U_squared'])),2),a.scale(div2(mul(z,d,i['I1_U_U_Z'])),2))
    return dict(m13=a.add(a.neg(half(S0)),bd('q_plus_q_y')),m15=a.sub(half(S1),bd('R_q_y')))


def _product(a,values):
    out=a.fixed(1)
    for value in values:out=a.mul(out,value)
    return out


def partial_integrals(g,integrals,upper,bodies):
    """Actual body derivatives, with the original measure/variable/domain."""
    out={}
    for key,q in integrals.items():
        template=copy.deepcopy(g.nodes[q.value.node])
        if template.pop('operation')!='definite_integral' or not template['endpoint_independent_of_Z']:
            raise ValueError('Actual fixed-Z integral source required')
        refs=[]
        for k,ref in enumerate(rows(bodies[key])):
            attrs=dict(template,integrand=ref.node,upper=upper.node,ordinary_slow_Z_derivative_order=k)
            refs.append(g.node('definite_integral',**attrs))
        if len(refs)!=4:raise ValueError('Genuine C3 integral rows required')
        out[key]=FunctionJet(tuple(refs))
    return out


def pulse_integrands(g,pulse,chart):
    """Genuine U C5/Q C4 supply every nonradial product through Z3."""
    a=Algebra(g,3);f=pulse.leading_fields[chart];U,Q=a.trim(f['U']),a.trim(f['Q'])
    Uz,Qz=a.shift(f['U']),a.shift(f['Q']);UU=a.mul(U,U);Uqz=a.mul(U,Qz);Uuz=a.mul(U,Uz)
    products=dict(I0_q_Z=Qz,I0_U_q_Z=Uqz,I0_U_squared=UU,I0_U_U_Z=Uuz,
        I1_q=Q,I1_q_Z=Qz,I1_q_squared=a.mul(Q,Q),I1_U_q_Z=Uqz,I1_U_squared=UU,I1_U_U_Z=Uuz)
    meta=pulse.pulse.own.fields[chart];R=meta['radius'];mu=pulse.pulse.own.jets.complete.Gamma.schedule['mu']
    jac=R if chart=='end' else g.quotient(R,mu,'same original positive pulse mu')
    bodies={k:a.scale(q,jac if k.startswith('I0_') else g.mul(R,jac)) for k,q in products.items()}
    if set(products)!=set(pulse.integrals[chart]):raise ValueError('Original ten signed nonradial products required')
    return bodies,products


def pulse_pressure_at(g,pulse,chart,value):
    a=Algebra(g,3);names=('main','gap','end')
    if chart not in names:raise ValueError('Original pulse chart required')
    meta=pulse.pulse.own.fields[chart];x=meta['query'];fields=pulse.leading_fields[chart]
    bodies,products=pulse_integrands(g,pulse,chart)
    current=partial_integrals(g,pulse.integrals[chart],value,bodies)
    total={key:a.add(*(a.trim(pulse.integrals[name][key]) for name in names[:names.index(chart)]),q)
           for key,q in current.items()}
    at=lambda q:a.at(q,x,value)
    R=a.at(a.fixed(meta['radius']),x,value).value
    upper=endpoint_functions(g,at(fields['U']),at(fields['Q']),at(fields['Q_y']),R)
    lower={k:a.trim(q) for k,q in pulse.global_endpoints[0].items()}
    delta=pulse.pulse.own.jets.complete.Gamma.schedule['delta']
    pressure=pressure_functions(g,delta,total,lower,upper)
    return pressure,dict(current_integrals=current,current_integrands=bodies,products=products,
        total_integrals=total,lower=lower,upper=upper,
        coordinate=x,value=value,radius=R,preceding_complete_charts=list(names[:names.index(chart)]))


def zero_velocity_fields(g,cumulative,R,delta):
    """Outside the inner/repair supports, keep the complete own mass datum."""
    a=Algebra(g,3);b=Algebra(g,2);zero=a.fixed(0)
    mass=a.trim(cumulative['m11']);mean=a.scale(mass,g.quotient(g.one,R,'positive actual nonaxis regional radius'))
    Q=previous.own_Q(g,zero,mean,delta)
    return dict(M1=mass,mean1=mean,K1=mean,P1=a.trim(cumulative['m13']),
        Uz1=zero,Utheta1=zero,F1=zero,Q1=Q,V1=b.scale(Q,R),
        M1_y=zero,mean1_y=a.neg(mean),Q1_y=b.neg(Q),V1_y=b.fixed(0))


class CurrentOriginalN1IntermediateOwnC3ChartFunctions:
    def __init__(self,own_owner,require_checked=True):
        if type(own_owner) is not previous.CurrentOriginalN1RepairedOwnC3TerminalFunctions or not own_owner.acceptance_loaded:
            raise TypeError('Checked repaired own C3 source required')
        own_owner.assert_graph();self.own=own_owner;self.complete=own_owner.complete
        self.g=g=own_owner.g;self.family=copy.deepcopy(own_owner.family);self.hashes=dict(own_owner.hashes)
        self.quiet=own_owner.quiet;self.compact=own_owner.compact;self.outer=self.compact.outer
        self.patch=self.outer.patch;self.modified=self.patch.modified
        self.lower=self.complete.lower;self.pulse=self.complete.pulse
        self.sources=(own_owner,self.quiet,self.compact,self.outer,self.patch,self.modified,self.lower,self.pulse)
        for owner in self.sources:
            if owner.g is not g or owner.family!=self.family or not owner.acceptance_loaded:
                raise ValueError('Same actual checked function graph and family required')
        for module in (previous,modified_source,patch_source,outer_source,pulse_source):
            rec=parent.linear.read(module.RECEIPT)
            if not rec['all_passed'] or rec['source_family']!=self.family:raise ValueError('Checked regional pressure source required')
            for name,digest in rec['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        self._pin(Path(pressure_source.__file__).name,sha(Path(pressure_source.__file__).name))
        a=Algebra(g,3);self.delta=own_owner.delta;self.inner=self.complete.segments['inner'];self.charts={}
        lower_meta=self.lower.lower.fields;x=lower_meta['query'];domain=lower_meta['domain']
        self.lower_partial=partial_integrals(g,self.lower.pressure,x,self.lower.signed_bodies)
        self._before('lower',x,[domain['lower'],domain['upper']],lower_meta['radius'],self.lower_partial,False)
        self.preRc_names=[]
        for name,data in self.modified.charts.items():
            meta=data['source_chart'];self._before(name,meta['coordinate'],meta['domain'],meta['radius'],data['pressure'],True)
            self.preRc_names.append(name)
        self._before('Rm_patch',self.patch.coordinate,self.patch.source_chart['domain'],
                     self.patch.source_chart['radius'],self.patch.pressure,True);self.preRc_names.append('Rm_patch')
        for name,data in self.outer.charts.items():
            self._before(name,data['coordinate'],data['domain'],data['radius'],data['pressure'],True);self.preRc_names.append(name)
        if len(self.preRc_names)!=17:raise ValueError('All original17 pre-Rc charts required')
        cumulative=dict(own_owner.band_endpoint)
        for key in ('m13','m15'):
            cumulative[key]=a.add(cumulative[key],a.sub(a.trim(self.quiet.pressure[key]),a.trim(self.quiet.incoming_pressure[key])))
        self.charts['quiet']=dict(coordinate=self.quiet.s,domain=[g.zero,self.quiet.length],radius=self.quiet.R,
            fixed_pressure=self.quiet.pressure,cumulative=cumulative,
            fields=zero_velocity_fields(g,cumulative,self.quiet.R,self.delta),region='after2Rc',
            complete_incoming_mass_not_reset=True)
        self.Rp_cumulative={k:a.at(q,self.quiet.s,self.quiet.length) for k,q in cumulative.items()}
        self.pulse_bindings={}
        for name in ('main','gap','end'):
            meta=self.pulse.pulse.own.fields[name];x=meta['query'];fixed,binding=pulse_pressure_at(g,self.pulse,name,x)
            cumulative=dict(self.Rp_cumulative)
            for key in ('m13','m15'):cumulative[key]=a.add(cumulative[key],fixed[key])
            self.charts['pulse_'+name]=dict(coordinate=x,domain=[meta['domain']['lower'],meta['domain']['upper']],
                radius=meta['radius'],fixed_pressure=fixed,cumulative=cumulative,
                fields=zero_velocity_fields(g,cumulative,meta['radius'],self.delta),region='after2Rc',
                complete_incoming_mass_not_reset=True)
            self.pulse_bindings[name]=binding
        self.endpoints={name:{edge:self.at(name,value) for edge,value in zip(('left','right'),data['domain'])}
                        for name,data in self.charts.items()}
        self.bindings=dict(original_pressure_formula=parent.active.current.ast_binding(pressure_source.pressure_functions),
            original_finite_endpoint_formula=parent.active.current.ast_binding(pressure_source.endpoint_functions),
            original_positive_order_recovery=parent.active.current.ast_binding(previous.own_Q),
            actual_pulse_C3_integral_source=pulse_source.RECEIPT,
            original_pulse_global_endpoint_certificate=copy.deepcopy(self.pulse.pulse.endpoint_certificate))
        self.definition=dict(lower_plus_17_preRc_actual_own_C3_QV_C2_functions=True,
            quiet_and_main_gap_end_actual_own_pressure_C3_functions=True,
            original_positive_order_one_plus_delta_recovery=True,
            all_complete_nonzero_mass_pressure_histories_retained=True,
            pulse_partial_integrals_keep_original_bodies_measures_and_all_preceding_charts=True,
            pulse_global_startup_endpoint_keeps_native_mean=True,
            leading_pulse_velocity_not_misidentified_as_positive_order_velocity=True,
            own_velocities_zero_only_outside_inner_and_reserved_compact_supports=True,
            after2Rc_actual_mass_source_uses_checked_terminal_integral_not_forced_zero=True,
            reserved_band_queries_delegated_to_checked_previous_owner=True,
            numerical_C3_bounds_and_global_leading_joins_not_claimed=True,
            current_native_full_C5_inlet_bridge_and_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(own_owner,self.complete,self.quiet,self.compact,self.outer,self.patch,self.modified,self.lower,self.pulse,g)
        self.acceptance_loaded=False
        if require_checked:
            raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec.get('input_hashes',{}))
            archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not rec.get('all_passed') or not all(rec.get(k) for k in GATES) or any(rec.get(k) for k in OPEN)
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['intermediate_own_C3_chart_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual intermediate own C3 chart functions required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed intermediate own source '+name)
            self.acceptance_loaded=True
    def _before(self,name,x,domain,R,pressure,include_lower):
        a=Algebra(self.g,3);cumulative={k:a.trim(q) for k,q in self.inner.items()}
        for key in ('m13','m15'):
            cumulative[key]=a.add(cumulative[key],a.trim(pressure[key]),
                a.trim(self.complete.segments['lower'][key]) if include_lower else a.fixed(0))
        self.charts[name]=dict(coordinate=x,domain=list(domain),radius=R,fixed_pressure=pressure,cumulative=cumulative,
            fields=zero_velocity_fields(self.g,cumulative,R,self.delta),region='beforeRc',
            complete_incoming_mass_not_reset=True)
    def at(self,name,value):
        data=self.charts[name];x=data['coordinate'];g=self.g
        group=lambda fields:{k:Algebra(g,len(rows(q))-1).at(q,x,value) for k,q in fields.items()}
        return dict(fields=group(data['fields']),cumulative=group(data['cumulative']),
            radius=Algebra(g,3).at(Algebra(g,3).fixed(data['radius']),x,value).value)
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed intermediate input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(charts=self.charts,preRc_names=self.preRc_names,lower_partial=self.lower_partial,
            Rp_cumulative=self.Rp_cumulative,pulse_bindings=self.pulse_bindings,endpoints=self.endpoints,
            bindings=self.bindings,definition=self.definition))
    def assert_graph(self):
        self.own.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or self._owners!=(self.own,self.own.complete,self.own.quiet,self.own.compact,self.own.compact.outer,
                self.own.compact.outer.patch,self.own.compact.outer.patch.modified,self.own.complete.lower,self.own.complete.pulse,self.g)):
            raise ValueError('Actual intermediate own chart or source graph changed')
    def at_fraction(self,name,fraction):
        self.assert_graph()
        if name=='reserved':return self.own.at_band_fraction(fraction)
        if name not in self.charts:raise ValueError('Unknown original own-field chart')
        if not isinstance(fraction,(int,str,Fraction)) or isinstance(fraction,bool):raise TypeError('Exact rational coordinate fraction required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Native coordinate fraction must lie in[0,1]')
        lo,hi=self.charts[name]['domain'];value=self.g.add(lo,self.g.mul(self.g.constant(q),self.g.sub(hi,lo)))
        return dict(chart=name,exact_native_coordinate_fraction=str(q),**self.at(name,value),
            actual_own_functions_not_numeric_ranges=True,mean_mass_K_pressure_C3=True,QV_C2=True,
            global_join_or_numeric_accuracy_claimed=False)


def run(own_owner):
    began=time.monotonic();owner=CurrentOriginalN1IntermediateOwnC3ChartFunctions(own_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,intermediate_own_C3_chart_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('INTERMEDIATE_N1_OWN_C3_CHART_FUNCTIONS_DEFINED',len(owner.charts),len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(own_owner) reuses accepted source owners')
