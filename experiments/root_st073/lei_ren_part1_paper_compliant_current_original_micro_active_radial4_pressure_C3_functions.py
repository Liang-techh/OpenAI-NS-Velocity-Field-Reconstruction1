"""Actual active micro radial4 sources and full own micro seam.

Total derivatives act on periodic lifted loop composites, never on a raw
fractional-part function. The original inverse/integral bodies are retained.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_rminus_mixed_source_join as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,high,rc_source=previous.outer,previous.high,previous.rc_source
modified_source=previous.modified_source
NAME=PREFIX+'current_original_micro_active_radial4_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_micro_active_radial4_pressure_C3_functions_check.json'
GATES=('current_original_actual_active_micro_U_mean_y4_C5_Q_y4_C4_sources_installed',
    'current_original_actual_active_micro_pressure_N_y2_C3_sources_installed',
    'current_original_n1_first_second_micro_full_mixed_source_join_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)
CHARTS=('first_micro','second_micro')


class PeriodicTotalNativeDerivative(modified_source.NativeDerivative):
    """Local lifted derivative; periodic composite identities glue wraps."""
    def __init__(self,g,variable,phase,main,jacobian,N,guard,shared=None):
        self.main,self.jacobian,self.N,self.guard=main,jacobian,N,guard
        super().__init__(g,variable,phase,shared)
        if variable==main:self.cache[phase.node]=g.mul(N,jacobian)
        self.sigma_extensions=[]
    def other(self,variable):
        if variable.node not in self.shared['factories']:
            PeriodicTotalNativeDerivative(self.g,variable,self.phase,self.main,
                self.jacobian,self.N,self.guard,self.shared)
        return self.shared['factories'][variable.node]
    def sigma_order(self,x,order):
        if not 1<=order<=4:raise ValueError('Actual sigma source through order four required')
        t=self.g.symbol('actual_active_micro_sigma_argument')
        local=self.other(t);body=local.sigma_prime(t)
        for _ in range(1,order):body=local.derivative(body)
        result=self.sub(body,t,x);self.sigma_extensions.append(dict(argument=x.node,order=order,source=result.node))
        return result
    def derivative(self,ref):
        g=self.g;i=ref.node
        if i in self.cache:return self.cache[i]
        row=g.nodes[i];op=row['operation']
        if op=='current_original_unique_analytic_core_component':
            self.guard(row);terms=[]
            for coordinate,order in (('rho','radial_order'),('Z_variable','ordinary_slow_Z_derivative_order')):
                factor=self.derivative(source.FunctionRef(g,row[coordinate]))
                if factor!=g.zero:
                    attrs=copy.deepcopy(row);attrs[order]+=1;self.guard(attrs)
                    operation=attrs.pop('operation');new=g.node(operation,**attrs)
                    terms.append(g.mul(new,factor));self.core_extensions.append(new.node)
            value=g.add(*terms)
        elif op=='analytic_unary' and row['name']=='original_flat_sigma_second':
            x=source.FunctionRef(g,row['argument'])
            value=g.mul(self.sigma_order(x,3),self.derivative(x))
        else:return super().derivative(ref)
        self.cache[i]=value;return value


class CurrentOriginalMicroActiveRadial4PressureC3Functions:
    def __init__(self,base=None,require_checked=True):
        self.base=base if base is not None else previous.CurrentOriginalN1RminusMixedSourceJoin()
        if not self.base.acceptance_loaded:raise ValueError('Checked actual own rminus source required')
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.parameters=dict(self.base.parameters);self.lower=self.base.lower
        self.loop_input=copy.deepcopy(self.base.inputs['prescribed_loop_C5_bindings'])
        self.modified_input=copy.deepcopy(self.base.inputs['modified_Q4_pressure_C3_bindings'])
        self.loop={name:self.loop_input['charts'][name] for name in CHARTS}
        self.modified={name:self.modified_input['charts'][name] for name in CHARTS}
        self.background=copy.deepcopy(self.base.bridge['micro'])
        self.coordinate=self.ref(self.background['phase']);self.radius=self.ref(self.background['radius'])
        self.jacobian=self.ref(self.background['physical_measure']);self.N=self.ref(self.loop_input['N'])
        self.phase=self.ref(self.loop['first_micro']['actual_phase']);self.seam=g.one
        self.collar=[g.constant('7/8'),g.constant('9/8')]
        self.core_templates=copy.deepcopy(self.lower.core_templates)
        self.factories={};self.requests={};self.charts={};a=outer.Algebra(g,5);b=outer.Algebra(g,3)
        self.background_y={key:self.jet(ids) for key,ids in self.modified['first_micro']['radial_background'].items()}
        for name in CHARTS:
            loop=self.loop[name];old=self.modified[name];J=self.ref(loop['physical_measure'])
            D=PeriodicTotalNativeDerivative(g,self.ref(loop['coordinate']),self.ref(loop['actual_phase']),
                self.coordinate,J,self.N,self.guard_core)
            D.shared['family']=self.family;self.factories[name]=D;self.requests[name]=[]
            for key,ids in self.background['fields'].items():
                for node,dy in zip(ids,outer.rows(self.background_y[key])):D.seed(self.ref(node),g.mul(J,dy))
            invJ=g.quotient(g.one,J,'same positive original micro physical Jacobian')
            physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])] for key in ('U','m','Q')}
            for n in range(2,5):
                refs=[]
                for node in outer.rows(physical['U'][-1]):
                    native=D.derivative(node);refs.append(g.mul(invJ,native))
                    self.requests[name].append(dict(input=node.node,native_output=native.node,physical_output=refs[-1].node,physical_y_order=n))
                physical['U'].append(outer.FunctionJet(tuple(refs)))
                physical['m'].append(a.sub(physical['U'][n-1],physical['m'][n-1]))
                physical['Q'].append(rc_source.radial.leading_Q(g,physical['U'][n],physical['m'][n],self.parameters['delta']))
            f={key+suffix:physical[key][n] for key in ('U','m','Q') for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
            density=rc_source.radial.pressure_density(g,f,self.parameters['delta'],self.radius)
            higher=high.higher_pressure_density(g,f,self.parameters['delta'],self.radius,density)
            Ny=[density['N1p'],higher['N1p_y'],higher['N1p_yy']]
            raw_pressure_y=[density['P1_y'],higher['P1_yy'],higher['P1_yyy']]
            raw_m15_y=[b.scale(raw_pressure_y[0],g.neg(self.radius)),higher['m15_yy'],higher['m15_yyy']]
            n0,n1,n2=Ny;negative_R2=g.neg(g.mul(self.radius,self.radius))
            pressure_y=[b.scale(n0,self.radius),b.scale(b.add(n0,n1),self.radius),b.scale(b.add(n0,b.scale(n1,2),n2),self.radius)]
            m15_y=[b.scale(n0,negative_R2),b.scale(b.add(b.scale(n0,2),n1),negative_R2),
                b.scale(b.add(b.scale(n0,4),b.scale(n1,4),n2),negative_R2)]
            intermediate=copy.deepcopy(self.base.lower.rin.rc.own.inputs[previous.common_source.previous.NAME]['charts'][name])
            values={key:self.jet(ids) for key,ids in intermediate['fields'].items()};fields={}
            for key in rc_source.previous.STATE:
                alg=outer.Algebra(g,2 if key in ('Q1','V1') else 3)
                fields[key]=[values[key]]+[alg.scale(values[key],(-1)**n) if key in ('mean1','K1','Q1') else alg.fixed(0) for n in range(1,5)]
            fields['P1']=[values['P1'],*pressure_y]
            cumulative={key:[self.jet(ids)]+[b.fixed(0) for _ in range(3)] for key,ids in intermediate['cumulative'].items()}
            cumulative['m13']=list(fields['P1']);cumulative['m15'][1:]=m15_y
            domain=[self.collar[0],self.seam] if name=='first_micro' else [self.seam,self.collar[1]]
            self.charts[name]=dict(domain=domain,physical=physical,density=density,higher_density=higher,N_y=Ny,pressure_y=pressure_y,
                raw_pressure_y=raw_pressure_y,raw_m15_y=raw_m15_y,
                m15_y=m15_y,intermediate=intermediate,own=dict(fields=fields,cumulative=cumulative))
        self.contract=copy.deepcopy(previous.CONTRACT)
        self.join={side:self.at(self.charts[name]['own'],self.seam) for side,name in zip(('left','right'),CHARTS)}
        self.definition=dict(actual_original_U_y_sources_differentiated_three_more_times=True,
            actual_native_core_radial4_Z6_from_same_unique_analytic_solution=True,
            actual_sigma_lazy_source_derivatives_not_flat_collar_zero_extension=True,
            total_native_phase_lift_N_times_J_only_in_main_native_factory=True,
            fractional_part_itself_has_no_global_ordinary_derivative_rule=True,
            periodic_composite_A_B_and_all_required_phase_mixed_rows_glue_wraps=True,
            same_original_mean_ODE_and_distinct_leading_one_minus_delta_Q_recovery=True,
            actual_signed_pressure_primitives_and_nonzero_own_histories_preserved=True,
            full_finite_micro_own_mixed_seam_not_global_n1_or_quantitative_bounds=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,node):return source.FunctionRef(self.g,node)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def decode_background(self):
        return dict(fields={k:self.jet(ids) for k,ids in self.background['fields'].items()},
            radial={k:self.jet(ids) for k,ids in self.background['radial'].items()})
    def guard_core(self,row):
        key=row.get('quantity');core_definition=previous.previous.core_definition
        if (key not in self.core_templates or core_definition(row)!=core_definition(self.core_templates[key])
            or row.get('source_family')!=self.family or not 0<=row.get('radial_order',-1)<=4
            or not 0<=row.get('ordinary_slow_Z_derivative_order',-1)<=6):
            raise ValueError('Same actual analytic core radial4/Z6 derivative required')
    def at(self,packet,value):return self.base.at(packet,value)
    def _snapshot(self):
        return outer.encode(dict(coordinate=self.coordinate,radius=self.radius,jacobian=self.jacobian,N=self.N,phase=self.phase,
            seam=self.seam,collar=self.collar,loop=self.loop,modified=self.modified,background=self.background,
            background_y=self.background_y,charts=self.charts,requests=self.requests,contract=self.contract,join=self.join,
            definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):raise ValueError('Actual active micro source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed active micro source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['micro_active_radial4_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):raise ValueError('Checked active micro source required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalMicroActiveRadial4PressureC3Functions(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        micro_active_radial4_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_ACTIVE_MICRO_RADIAL4_SOURCE_DEFINED',len(owner.g.nodes),flush=True);return owner,raw
