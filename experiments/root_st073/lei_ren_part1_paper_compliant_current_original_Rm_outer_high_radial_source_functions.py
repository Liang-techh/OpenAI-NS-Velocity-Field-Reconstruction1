"""Actual high physical-log-radius source functions on Rm and five outer charts.

Native Jacobians need not be constant. Every derivative applies J^-1 D_x
again, with the lifted phase derivative N*J and exact compact-beta calculus.
The pre-existing O3-power/Rc and compact/quiet sources are reused separately.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_annular_mixed_source_joins as previous
import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions as high_source
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions as patch_leading
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions as patch_source
import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as outer_source

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,high,radial=high_source.outer,high_source.high,high_source.rc_source.radial
PeriodicDerivative=high_source.previous.previous.PeriodicTotalNativeDerivative
NAME=PREFIX+'current_original_Rm_outer_high_radial_source_functions.json.gz'
RECEIPT=PREFIX+'current_original_Rm_outer_high_radial_source_functions_check.json'
CHARTS=('Rm_patch','Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition')
GATES=('current_original_Rm_and_five_outer_U_mean_y4_C5_Q_y4_C4_sources_installed',
       'current_original_Rm_and_five_outer_N_y2_pressure_y3_C3_sources_installed',
       'current_original_Rm_and_five_outer_own_high_radial_sources_installed',
       'current_original_nonconstant_J_total_physical_y4_source_calculus_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)


class RmOuterTotalNativeDerivative(PeriodicDerivative):
    """The admitted total-native calculus plus the original raw beta rule."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.beta_extensions=[]
    def other(self,variable):
        if variable.node not in self.shared['factories']:
            RmOuterTotalNativeDerivative(self.g,variable,self.phase,self.main,self.jacobian,self.N,self.guard,self.shared)
        return self.shared['factories'][variable.node]
    def derivative(self,ref):
        g=self.g;i=ref.node
        if i in self.cache:return self.cache[i]
        row=g.nodes[i]
        if row['operation']=='compact_raw_beta':
            x=source.FunctionRef(g,row['argument'])
            if g.nodes[patch_leading.raw_beta(g,x).node]!=row:raise ValueError('Same exact compact beta source required')
            first,_=patch_leading.beta_first_two(g,x,ref)
            result=g.mul(first,self.derivative(x));self.cache[i]=result
            self.beta_extensions.append(dict(input=i,argument=x.node,first_argument_derivative=first.node,native_output=result.node))
            return result
        return super().derivative(ref)


def fraction_packet(g,chart,name,fraction):
    q=Fraction(fraction)
    if not 0<=q<=1:raise ValueError('Fraction outside complete original native chart')
    lo,hi=chart['domain'];value=g.add(lo,g.mul(g.constant(q),g.sub(hi,lo)));x=chart['coordinate']
    at=lambda jets:[outer.Algebra(g,len(outer.rows(jet))-1).at(jet,x,value) for jet in jets]
    own={group:{field:at(jets) for field,jets in fields.items()} for group,fields in chart['own'].items()}
    return dict(chart=name,native_query=value,
        radius=outer.Algebra(g,0).at(outer.FunctionJet((chart['radius'],)),x,value),
        physical={field:at(jets) for field,jets in chart['physical'].items()},N_y=at(chart['N_y']),
        pressure_y=at(chart['pressure_y']),m15_y=at(chart['m15_y']),own=own,numerical_oracle_available=False)


class CurrentOriginalRmOuterHighRadialSourceFunctions:
    def __init__(self,base=None,require_checked=True):
        fresh=base is None
        self.base=base if base is not None else previous.CurrentOriginalAnnularMixedSourceJoins()
        if not self.base.acceptance_loaded:raise ValueError('Accepted annular source interfaces required')
        if fresh and require_checked:previous._restore_checked_layout(self.base,NAME,RECEIPT)
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        inputs={}
        for module,key in ((patch_leading,'true_Rm_patch_C6_bindings'),(patch_source,'Rm_patch_finite_N_C5_bindings'),
                           (outer_source,'true_outer_C6_C5_bindings')):
            raw,rec=parent.linear.read(module.NAME),parent.linear.read(module.RECEIPT)
            if (not rec['all_passed'] or raw['source_family']!=self.family
                or not all(rec[k] for k in module.GATES)
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Same accepted original source required: '+module.NAME)
            for filename,digest in rec['input_hashes'].items():
                if self.hashes.get(filename)!=digest or sha(filename)!=digest:raise ValueError('Changed Rm/outer input '+filename)
            inputs[key]=copy.deepcopy(raw[key])
        self.leading_input=inputs['true_Rm_patch_C6_bindings'];self.patch_input=inputs['Rm_patch_finite_N_C5_bindings']
        self.outer_input=inputs['true_outer_C6_C5_bindings'];self.parameters=dict(self.base.parameters);self.N=self.base.N
        if self.patch_input['N']!=self.outer_input['N'] or self.patch_input['N']!=self.N.node:
            raise ValueError('One common actual finite N required')
        self.backgrounds=dict(Rm_patch=self.leading_input['chart'],
            **{name:self.outer_input['background'][name] for name in CHARTS[1:]})
        self.modified=dict(Rm_patch=self.patch_input,**{name:self.outer_input['charts'][name] for name in CHARTS[1:]})
        self.loop=dict(Rm_patch=self.patch_input['chart'],**{name:self.outer_input['charts'][name] for name in CHARTS[1:]})
        self.free_seeds={row['node']:set(row['variables']) for row in self.leading_input['free_source_seeds']}
        bump=patch_leading.ExactBumpWeights(g);x=self.ref(self.loop['Rm_patch']['coordinate'])
        self.gamma=[bump.gamma(x,i) for i in range(3)]
        self.factories={};self.background_y={};self.requests={};self.charts={}
        for name in CHARTS:
            self.charts[name]=self.build_chart(name)
            print('ACTUAL_RM_OUTER_HIGH_SOURCE_DEFINED',name,len(g.nodes),flush=True)
        self.definition=dict(actual_six_complete_Rm_and_outer_charts_not_ordinary_Z_order_relabelling=True,
            every_high_U_row_is_repeated_J_inverse_times_total_native_derivative=True,
            variable_J_derivatives_and_main_NJ_phase_retained=True,
            Z_only_exact_Rm_implicit_controls_and_admitted_common_core_guard=True,
            original_flat_compact_beta_derivatives_not_numeric_or_zero_extension=True,
            complete_nonzero_background_correction_pressure_and_own_histories_retained=True,
            O3_power_Rc_compact_quiet_and_combined_outer_sources_reused_not_rebuilt=True,
            high_source_coverage_not_interface_closure_numeric_bounds_or_global_n1=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()

    def __getattr__(self,name):return getattr(self.base,name)
    def ref(self,node):return source.FunctionRef(self.g,node)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def build_chart(self,name):
        g=self.g;old=self.modified[name];loop=self.loop[name];bg=self.backgrounds[name]
        x,R,J,phase=(self.ref(loop[key]) for key in ('coordinate','radius','physical_measure','actual_phase'))
        self.background_y[name]={key:self.jet(ids) for key,ids in old['radial_background'].items()}
        D=RmOuterTotalNativeDerivative(g,x,phase,x,J,self.N,self.base.active.guard_core)
        D.shared['family']=self.family;D.shared['free'].update({i:set(v) for i,v in self.free_seeds.items()})
        self.factories[name]=D;self.requests[name]=[]
        for key,ids in bg['fields'].items():
            for node,dy in zip(ids,outer.rows(self.background_y[name][key])):D.seed(self.ref(node),g.mul(J,dy))
        if name=='Rm_patch':
            for gamma in self.gamma:D.seed(gamma[0],gamma[1]);D.seed(gamma[1],gamma[2])
        invJ=g.quotient(g.one,J,'positive actual original variable log-radius Jacobian')
        a=outer.Algebra(g,5);b=outer.Algebra(g,3)
        physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])] for key in ('U','m','Q')}
        for n in range(2,5):
            refs=[]
            for ref in outer.rows(physical['U'][-1]):
                native=D.derivative(ref);physical_ref=g.mul(invJ,native);refs.append(physical_ref)
                self.requests[name].append(dict(input=ref.node,native_output=native.node,physical_output=physical_ref.node,
                    physical_y_order=n,Jacobian_inverse=invJ.node))
            physical['U'].append(outer.FunctionJet(tuple(refs)))
            physical['m'].append(a.sub(physical['U'][n-1],physical['m'][n-1]))
            physical['Q'].append(radial.leading_Q(g,physical['U'][n],physical['m'][n],self.parameters['delta']))
        f={key+suffix:physical[key][n] for key in ('U','m','Q') for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        density=radial.pressure_density(g,f,self.parameters['delta'],R)
        higher=high.higher_pressure_density(g,f,self.parameters['delta'],R,density)
        Ny=[density['N1p'],higher['N1p_y'],higher['N1p_yy']];n0,n1,n2=Ny
        P=[b.scale(n0,R),b.scale(b.add(n0,n1),R),b.scale(b.add(n0,b.scale(n1,2),n2),R)]
        minusR2=g.neg(g.mul(R,R));m15=[b.scale(n0,minusR2),b.scale(b.add(b.scale(n0,2),n1),minusR2),
            b.scale(b.add(b.scale(n0,4),b.scale(n1,4),n2),minusR2)]
        packet=copy.deepcopy(self.base.intermediate['charts'][name]);values={key:self.jet(ids) for key,ids in packet['fields'].items()}
        fields={}
        for key in high_source.rc_source.previous.STATE:
            alg=outer.Algebra(g,2 if key in ('Q1','V1') else 3)
            fields[key]=[values[key]]+[alg.scale(values[key],(-1)**n) if key in ('mean1','K1','Q1') else alg.fixed(0) for n in range(1,5)]
        fields['P1']=[values['P1'],*P]
        cumulative={key:[self.jet(ids)]+[b.fixed(0) for _ in range(3)] for key,ids in packet['cumulative'].items()}
        cumulative['m13']=list(fields['P1']);cumulative['m15'][1:]=m15
        return dict(domain=[self.ref(node) for node in loop['domain']],coordinate=x,radius=R,jacobian=J,phase=phase,
            physical=physical,density=density,higher_density=higher,N_y=Ny,pressure_y=P,m15_y=m15,
            raw_pressure_y=[density['P1_y'],higher['P1_yy'],higher['P1_yyy']],intermediate=packet,
            own=dict(fields=fields,cumulative=cumulative))
    def at_fraction(self,name,fraction):
        if name not in CHARTS:raise ValueError('Actual Rm/outer original chart required')
        return fraction_packet(self.g,self.charts[name],name,fraction)
    def _snapshot(self):
        return outer.encode(dict(leading_input=self.leading_input,patch_input=self.patch_input,outer_input=self.outer_input,
            backgrounds=self.backgrounds,modified=self.modified,loop=self.loop,free_seeds=[dict(node=i,variables=sorted(v)) for i,v in self.free_seeds.items()],
            gamma=self.gamma,background_y=self.background_y,charts=self.charts,requests=self.requests,definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self._snapshot()!=self.snapshot or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):raise ValueError('Actual Rm/outer high sources changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed Rm/outer high input '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or len(rec.get('checks',{}))!=6 or not all(rec[key] for key in GATES)
            or any(rec[key] for key in OPEN) or pins.pop(NAME)!=sha(NAME)
            or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['Rm_outer_high_radial_source_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual Rm/outer high sources required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalRmOuterHighRadialSourceFunctions(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        Rm_outer_high_radial_source_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_SIX_RM_OUTER_HIGH_RADIAL_SOURCE_CHARTS_DEFINED',len(owner.prefix),flush=True)
    return owner,raw
