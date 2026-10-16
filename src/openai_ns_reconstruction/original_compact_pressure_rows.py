"""Source-bound C3 compact pressure and exact quiet continuation to Rp."""
import copy
from fractions import Fraction
import gzip
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_post_axial_pressure_rows import OriginalPostAxialPressureRows
from .original_compact_numeric import supported_mass
from .original_Rm_finite_N_numeric import RATES


class OriginalCompactPressureRows(OriginalPostAxialPressureRows):
    DOMAINS={'compact':(Fraction(1),Fraction(2))}
    LOOP_CHARTS={'compact':'actual_selected_N_compact_profile_physical_y'}
    ORIGINS={'compact':14}
    START='1'

    def _archives(self,specifications):
        out={}
        for stem,gates in specifications:
            path=self.directory/(PREFIX+'compliant_'+stem+'_check.json');raw=path.read_bytes();receipt=json.loads(raw)
            if receipt['source_family']!=self.family or not receipt['all_passed'] or not all(receipt[k] for k in gates):
                raise ValueError('Accepted same original compact/quiet source required: '+stem)
            for name,digest in receipt['input_hashes'].items():
                if self.hashes.get(name,digest)!=digest or hashlib.sha256((self.directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed compact/quiet source: '+name)
                self.hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
            path=self.directory/(PREFIX+'compliant_'+stem+'.json.gz');raw=path.read_bytes()
            self.hashes[path.name]=hashlib.sha256(raw).hexdigest();out[stem]=json.loads(gzip.decompress(raw))
            if out[stem]['source_family']!=self.family:raise ValueError('Conflicting source family')
        return out

    def __init__(self,*,previous=None,**kwargs):
        if previous is not None and kwargs:raise ValueError('Supply Rc predecessor or constructor options')
        p=previous if previous is not None else OriginalPostAxialPressureRows(**kwargs).next_field().next_field()
        if p.chart!='O3_power':raise ValueError('Complete original Rc pressure predecessor required')
        self.previous=p;self.chart='compact';self.incoming_coordinate='2'
        for name in ('ctx','iv','J','audit','directory','family','o','t','micro','owner','reference','cover_owner','patch_rows',
            'axial','mean_source_summary','mean_source_proof','mean_point','mean_truth'):
            setattr(self,name,getattr(p,name))
        self.hashes,self.code_hashes=dict(p.hashes),dict(p.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_post_axial_m15_numeric.json');raw=path.read_bytes();data=json.loads(raw)
        if not data['checks']['passed'] or data['source_family']!=self.family:raise ValueError('Accepted complete Rc pressure state required')
        for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,digest in group.items():
                if target.get(name,digest)!=digest or hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed original Rc pressure input: '+name)
                target[name]=digest
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self.rc_inlets={Fraction(p['Z']):p for p in data['samples'] if p['chart']=='O3_power' and p['native_coordinate']=='2'}
        sources=self._archives((
            ('current_original_C3_limit_controls',('current_original_same_Banach_limit_C3_controls_installed','current_original_same_repair_integer_quantitative_C3_Picard_tail_installed')),
            ('current_original_C3_repair_band',('current_original_same_limit_repair_band_C3_velocity_histories_installed','current_original_same_limit_repair_band_C3_relative_terminal_identities_installed','current_original_same_limit_repair_band_C3_directed_magnitude_ranges_installed')),
            ('current_original_true_compact_C5_functions',('current_original_true_same_branch_compact_controls_C5_functions_installed','current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed'))))
        control=sources['current_original_C3_limit_controls'];band=sources['current_original_C3_repair_band'];functions=sources['current_original_true_compact_C5_functions']
        self.ranges=band['actual_four_Z_C3_band_ranges'];self.binding=functions['true_compact_C5_bindings'];nodes=functions['exact_graph_nodes']
        selected=control['actual_quantitative_C3_bounds']['actual_same_repair_integer']
        if not selected['J']==band['actual_selected_repair_integer']['J']==nodes[self.binding['N']]['J']==self.owner.integer['J']:
            raise ValueError('Original exact selected N changed')
        self.certificates={}
        for name,node in self.binding['terminal_certificates'].items():
            v=nodes[node]
            if v['operation']!='proved_original_compact_relative_C5_terminal_identity' or not v['relative_correction_only'] or not v['absolute_leading_histories_and_pressure_not_zeroed'] or v['source_system']!=self.binding['system']:
                raise ValueError('Actual original relative terminal source identity required')
            self.certificates[name]=dict(node=node,source_system=v['source_system'],relative_only=True,ordinary_Z_orders=list(range(6)))
        self.L=self.iv.log(2);self.supports=[(self.L*(self.iv.mpf(v)-self.iv.mpf(1)/40),self.L*(self.iv.mpf(v)+self.iv.mpf(1)/40)) for v in ('.2','.5','.8')]
        self._reset()
        for name in ('original_compact_pressure_rows.py','original_compact_numeric.py'):
            name='src/openai_ns_reconstruction/'+name;self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _reset(self):
        self._histories={};self._incoming={};self._backgrounds={};self._packets={};self._errors={};self._scalars={};self._uniforms={};self._next=None

    def next_field(self):
        if self._next is None:self._next=OriginalQuietPressureRows(self)
        return self._next

    def _matched(self,Z):
        lo,hi=_bounds(self.audit,self.iv.mpf(Z))
        frames=[p for p in self.ranges if not (self.audit.mpf(p['exact_Z_cell'][1])<lo or self.audit.mpf(p['exact_Z_cell'][0])>hi)]
        if not frames or any(p['source_family']!=self.family or not p['source_derivative_units_and_independent_P0_retained'] or not p['selected_repair_integer_unchanged'] for p in frames):
            raise ValueError('Actual compact C3 source cover required')
        return frames

    def geometry(self,value):
        q,p,t=self.coordinate(value);c,iv=self.ctx,self.iv;sp,st=c.log(p),iv.log(t)
        terms=list(self.axial.geometry('1')['atomic_log_radius_terms'])+[(c.mpf(14),iv.mpf(14)),(sp,st)]
        return dict(native=q,local_point=sp,local=st,atomic_log_radius_terms=terms,
            logR_point=sum((a for a,b in terms),c.mpf(0)),logR=sum((b for a,b in terms),iv.mpf(0)),
            physical_jacobian_point=1/p,physical_jacobian=1/t,exact_start=q==1)

    def scalars(self,value):
        if value not in self._scalars:
            _,p,t=self.coordinate(value)
            self._scalars[value]=dict(point=dict(beta=self.ctx.mpf(0),beta_y=self.ctx.mpf(0),mean=self.mean_point*self.ctx.exp(-14)/p),
                truth=dict(beta=self.iv.mpf(0),beta_y=self.iv.mpf(0),mean=self.mean_truth*self.iv.exp(-14)/t))
        return self._scalars[value]

    def uniform_background(self,Z,value):
        key=Z,value
        if key not in self._uniforms:
            end=self.scalars(value)['truth']['mean'];start=self.scalars(self.START)['truth']['mean']
            mean=self.iv.mpf([max(0,_bounds(self.audit,end)[0]),_bounds(self.audit,start)[1]])
            self._uniforms[key]=self.fields(Z,self.iv,self.iv.mpf(0),mean,self.iv.mpf(0))
        return self._uniforms[key]

    def minimum_log(self,*values):
        if any(v is None for v in values):return None
        return self.iv.mpf(min(_bounds(self.audit,v)[1] for v in values))

    def cap(self,frames,group,name,k):return self.cover_owner._log_cover([p[group][name][k] for p in frames])

    def history_logs(self,Z,value):
        key=Z,value
        if key in self._histories:return self._histories[key]
        q,_,_=self.coordinate(value);iv=self.iv;geo=self.geometry(value);s=geo['local'];frames=self._matched(Z)
        incoming=self.incoming_history_logs(Z);hist={};prefix={};proof={}
        for name,rate in RATES.items():
            supports=[self.supports[i] for i in ((0,2) if name in ('m','k') else (0,1,2))]
            fm,fp=supported_mass(iv,s,self.L,rate,'forward',supports);bm,bp=supported_mass(iv,s,self.L,rate,'backward',supports)
            pm,pp=supported_mass(iv,s,self.L,0,'forward',supports)
            masslog=lambda v:None if _bounds(self.audit,v)[1]==0 else iv.log(v)
            fl,bl,pl=masslog(fm),masslog(bm),masslog(pm);r=iv.mpf(rate);values=[];sup=[]
            for k in range(4):
                density=self.cap(frames,'actual_signed_density_magnitude_C3_bounds',name,k)
                initial=self.minimum_log(incoming[name][k],self.cap(frames,'actual_incoming_correction_C3_bounds',name,k))
                forward=self.sum_logs((None if initial is None else initial-r*s,None if density is None or fl is None else density+fl))
                backward=None if density is None or bl is None else density+bl
                values.append(self.minimum_log(forward,backward))
                sup.append(self.sum_logs((initial,None if density is None or pl is None else density+pl)))
            hist[name],prefix[name]=values,sup
            proof[name]=dict(rate=rate,forward_mass=fm,backward_mass=bm,prefix_source_width=pm,
                forward_parts=fp,backward_parts=bp,terminal_source_certificate=self.certificates[name],
                incoming_and_density_already_contain_one_over_N=True)
        sl,sh=_bounds(self.audit,s);local=[];visited=[]
        for i,(a,b) in enumerate(self.supports):
            if not (sh<_bounds(self.audit,a)[0] or sl>_bounds(self.audit,b)[1]):local.append(i)
            if sh>=_bounds(self.audit,a)[0]:visited.append(i)
        out=dict(history=hist,prefix_history=prefix,frames=[p['exact_Z_cell'] for p in frames],source_frames=frames,
            local_bumps=local,visited_bumps=visited,groups=[],local=frames,visited_cells=frames,steps=[],supported_transport=proof,
            exact_relative_terminal=q==2)
        self._histories[key]=out;return out

    def error_data(self,Z,value,prefix=False):
        key=Z,value,prefix
        if key in self._errors:return self._errors[key]
        q,_,_=self.coordinate(value);h=self.history_logs(Z,value);frames=h['source_frames'];lp=self.o.core.lp_iv
        active=bool(set(h['visited_bumps' if prefix else 'local_bumps'])&{0,2})
        shift=lambda v:None if v is None else v+lp
        profile=lambda y,k:self.cover_owner._log_cover([p['actual_selected_N_profile_y0_y1_y2_C3_bounds'][y]['delta_V'][k] for p in frames]) if active else None
        errors=dict(U=[shift(profile(0,k)) for k in range(4)],m=[shift(v) for v in h['prefix_history' if prefix else 'history']['m']],
            U_y=[shift(profile(1,k)) for k in range(2)])
        components=dict(actual_physical_deltaV_C3=self.log_rows(errors['U']),actual_physical_deltaV_y_C1=self.log_rows(errors['U_y']),
            archived_profiles_already_include_one_over_N=True,extra_N_factor=False,extra_native_x_or_hb_factor=False)
        if prefix or q==1:
            old,_=self.previous.error_data(Z,self.incoming_coordinate)
            for name in ('U','U_y'):errors[name]=[self.sum_logs((a,b)) for a,b in zip(errors[name],old[name])]
            components['predecessor_and_own_chart_budgets_unioned']=True
        out=self.finish_errors(Z,errors),components;self._errors[key]=out;return out

    def evaluate(self,Z,value):
        out=super().evaluate(Z,value);h=self.history_logs(Z,value);geo=self.geometry(value)
        from .original_Rm_patch_numeric import encode
        out['geometry'].update(exact_start=geo['exact_start'],native_coordinate_is_physical_log_offset=False,
            physical_Jacobian=self.o._summary(geo['physical_jacobian_point'],geo['physical_jacobian']),native_definition='x=R/Rc; dy=dx/x')
        out['supported_C3_transport']=encode(h['supported_transport'])
        out['source_cover'].update(active_bumps=h['local_bumps'],visited_bumps=h['visited_bumps'],
            same_original_unique_compact_branch=True,relative_terminal_source_certificates=copy.deepcopy(self.certificates))
        out['contract'].update(actual_selected_compact_profiles_already_N_scaled=True,forward_backward_C3_source_bounds=True,
            only_relative_histories_zero_at_2Rc=True,full_prefix_not_terminal_zero_for_volume=True)
        return out


class OriginalQuietPressureRows(OriginalCompactPressureRows):
    DOMAINS={'quiet':(Fraction(0),Fraction(1))}
    LOOP_CHARTS={'quiet':'exact_original_quiet_source_identity'}
    ORIGINS={'quiet':'14+log2'}
    START='0'

    def __init__(self,previous):
        if previous.chart!='compact':raise ValueError('Original compact terminal required')
        self.previous=previous;self.chart='quiet';self.incoming_coordinate='2'
        for name in ('ctx','iv','J','audit','directory','family','o','t','micro','owner','reference','cover_owner','patch_rows',
            'axial','mean_source_summary','mean_source_proof','mean_point','mean_truth','certificates','ranges','rc_inlets'):
            setattr(self,name,getattr(previous,name))
        self.hashes,self.code_hashes=dict(previous.hashes),dict(previous.code_hashes)
        raw=self._archives((
            ('current_original_true_quiet_Rp_C5_pressure_C3_functions',('current_original_true_2Rc_to_Rp_quiet_C5_pressure_C3_functions_installed','current_original_true_quiet_K_matches_original_native_source_installed')),
            ('current_original_n1_repaired_own_C3_terminal_functions',('current_original_n1_mass_pressure_C3_QV_C2_tail_functions_installed',))))
        self.quiet_binding=raw['current_original_true_quiet_Rp_C5_pressure_C3_functions']['true_quiet_Rp_C5_C3_bindings']
        self.quiet_mass_theorem=raw['current_original_true_quiet_Rp_C5_pressure_C3_functions']['exact_graph_nodes'][self.quiet_binding['linear_mass_theorem']]
        self.tail=raw['current_original_n1_repaired_own_C3_terminal_functions']['repaired_own_C3_terminal_bindings']
        if not self.tail['definition']['mass_mean_K_C3_QV_C2_zero_after_2Rc'] or self.quiet_mass_theorem['operation']!='proved_original_quiet_linear_mass_source_identity' or self.quiet_mass_theorem['current_quiet_identity']!='physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q':
            raise ValueError('Full original quiet relative-zero and pressure identities required')
        self.logmu=self.ctx.log(self.ctx.mpf('.001'))-4*self.o.core.logP
        self.logmui=self.iv.log(self.iv.mpf('.001'))-4*self.o.core.lp_iv
        self.length=-60*self.logmu-2-self.ctx.log(2);self.lengthi=-60*self.logmui-2-self.iv.log(2)
        if _bounds(self.audit,self.lengthi)[0]<=0:raise ValueError('Original positive quiet length required')
        self._reset()
        name='src/openai_ns_reconstruction/original_quiet_numeric.py';self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):raise ValueError('Rp requires an independent reserved pressure/source seam')

    def geometry(self,value):
        q,p,t=self.coordinate(value);c,iv=self.ctx,self.iv;sp,st=self.length*p,self.lengthi*t
        terms=list(self.axial.geometry('1')['atomic_log_radius_terms'])+[(c.mpf(14),iv.mpf(14)),(c.log(2),iv.log(2)),(sp,st)]
        return dict(native=q,local_point=sp,local=st,atomic_log_radius_terms=terms,
            logR_point=sum((a for a,b in terms),c.mpf(0)),logR=sum((b for a,b in terms),iv.mpf(0)),
            physical_jacobian_point=self.length,physical_jacobian=self.lengthi,exact_start=q==0)

    def scalars(self,value):
        if value not in self._scalars:
            g=self.geometry(value)
            self._scalars[value]=dict(point=dict(beta=self.ctx.mpf(0),beta_y=self.ctx.mpf(0),mean=self.mean_point*(self.ctx.exp(-14)/2)*self.ctx.exp(-g['local_point'])),
                truth=dict(beta=self.iv.mpf(0),beta_y=self.iv.mpf(0),mean=self.mean_truth*(self.iv.exp(-14)/2)*self.iv.exp(-g['local'])))
        return self._scalars[value]

    def history_logs(self,Z,value):
        self.coordinate(value);incoming=self.incoming_history_logs(Z)
        if any(v is not None for rows in incoming.values() for v in rows):raise ValueError('Relative compact terminal must be exactly zero by source certificate')
        zero={name:[None]*4 for name in RATES}
        return dict(history=zero,prefix_history=zero,frames=[p['exact_Z_cell'] for p in self._matched(Z)],
            groups=[],local=[],visited_cells=[],steps=[],supported_transport={},local_bumps=[],visited_bumps=[])

    def error_data(self,Z,value,prefix=False):
        self.history_logs(Z,value)
        errors=dict(U=[None]*4,m=[None]*4,U_y=[None]*2)
        return self.finish_errors(Z,errors),dict(full_original_quiet_source_identity=True,
            relative_finite_N_histories_and_local_velocity_exact_zero=True,absolute_leading_mean_and_pressure_not_reset=True)

    def evaluate(self,Z,value):
        out=super().evaluate(Z,value)
        out['geometry']['native_definition']='q in[0,1], s=Lq; dy=Ldq once'
        out['quiet_source_identity']=copy.deepcopy(self.quiet_binding)
        out['contract'].update(full_quiet_not_just_leading_Qy_minus_Q=True,finite_N_zero_from_source_identity=True,
            incoming_cumulative_pressure_not_zeroed=True)
        return out
