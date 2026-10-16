"""Current Rp inlet and original selected main/gap/end pulse values.

Archive-backed coefficient enclosures feed the unchanged original segmented
integrators. No selected/future ancestor is reconstructed or mutated.
Interior P1 transport is a separate dependency; missing P1 is never zero.
"""
import copy
from fractions import Fraction
import gzip
import hashlib
import json
from types import SimpleNamespace, FunctionType, MethodType

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_quiet_numeric import OriginalQuietVelocity


class OriginalNativePulseNumeric:
    def __init__(self, *, precision=85, amplitude_precision=110, directory=None, cells=128):
        if type(cells) is not int or not 32<=cells<=1024: raise ValueError('Pulse integration cells32..1024 required')
        self.quiet = OriginalQuietVelocity(precision=precision, amplitude_precision=amplitude_precision, directory=directory)
        f = self.quiet
        # The numerical Taylor engine uses its own directed context; the
        # physical amplitude/time context is deliberately separate.
        backend=f.numeric.compact.numeric
        self.ctx,self.iv,self.J,self.amplitude = backend.ctx,backend.iv,backend.J,f.amplitude
        self.family,self.directory = f.family,f.numeric.directory
        self.hashes,self.code_hashes = dict(f.hashes),dict(f.code_hashes)
        self.cells = cells; self._data_cache = {}; self._frame_cache = {}; self._packets = {}
        path = self.directory/(PREFIX+'compliant_current_original_quiet_numeric.json')
        accepted = json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family: raise ValueError('Accepted actual Rp quiet numerical inlet required')
        for n,h in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h: raise ValueError('Changed quiet implementation: '+n)
            self.code_hashes[n] = h
        self.hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.inlets = {p['Z']:p for p in accepted['samples'] if p['native_quiet_fraction']=='1'}
        specs=(('current_original_Rp_interval_inlet','.json',('current_original_Rp_source_bound_interval_constant_view_installed',)),
            ('current_original_Rp_selected_pulse','.json',('current_original_Rp_interval_constants_consumed_by_selected_pulse',)),
            ('current_original_Rp_native_C5_pulse_source','.json.gz',('current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions',)),
            ('current_original_Rp_pulse_coefficients_refined','.json.gz',('current_original_Rp_selected_coefficients_refined_installed',)),
            ('current_original_Rp_pulse_basis_energy','.json.gz',('current_original_Rp_pulse_basis_energy_refined_installed',)))
        raw = {}
        for stem,ext,gates in specs:
            receipt_path=self.directory/(PREFIX+'compliant_'+stem+'_check.json')
            receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not all(receipt[g] for g in gates) or receipt['source_family']!=self.family: raise ValueError('Checked original pulse input required: '+stem)
            for n,h in receipt['input_hashes'].items():
                if n in self.hashes:
                    if self.hashes[n]!=h: raise ValueError('Conflicting pulse input: '+n)
                elif hashlib.sha256((self.directory/n).read_bytes()).hexdigest()!=h: raise ValueError('Changed pulse input: '+n)
                self.hashes[n]=h
            self.hashes[receipt_path.name]=hashlib.sha256(receipt_path.read_bytes()).hexdigest()
            p=self.directory/(PREFIX+'compliant_'+stem+ext);data=p.read_bytes()
            self.hashes[p.name]=hashlib.sha256(data).hexdigest();raw[stem]=json.loads(gzip.decompress(data) if ext.endswith('.gz') else data)
            if raw[stem]['source_family']!=self.family: raise ValueError('Same original pulse family required')
        self.constants={n:self._decode(v) for n,v in raw['current_original_Rp_interval_inlet']['current_interval_constant_functions'].items()}
        self.selection=raw['current_original_Rp_pulse_coefficients_refined']['samples']
        self.binding=raw['current_original_Rp_native_C5_pulse_source']['Rp_native_C5_pulse_source_bindings']
        basis=raw['current_original_Rp_pulse_basis_energy']['view'];self.normal=self._decode(basis['original_normalization'])
        self.gram=self._decode(basis['bump_basis']['energy_gram'])
        self.mu,self.delta = f.numeric.outer.mui,self.iv.mpf(self.amplitude.delta_iv)
        from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField
        from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
        self.decay_integral=decay_integral
        # Only pure original segmented methods are reused. Their data callback
        # is this checked actual Rp frame and the same original selected jets.
        field=object.__new__(CompliantAxialPulseField);field.ctx=self.iv
        field.mu,field.delta=self.mu,self.delta
        field.rate,field.prate=1-self.mu,1+2*self.mu
        self.field=field
        field.data=self._data
        field.selection=SimpleNamespace(future=SimpleNamespace(future=self._future),
                                        K=self._decode(basis['fixed_pulse_energy_Kpulse']))
        field.pressure_moment=self._pressure
        # |d beta_j/ds| <= 2/(ell^2*N)*sup_(x>=1) x^2*exp(-x).
        field.beta_derivative_sup=8*self.iv.exp(-2)/(self.iv.mpf('.15')**2*self.normal)
        field.pulse=SimpleNamespace(initial=SimpleNamespace(repair=SimpleNamespace(normalization=self.normal)),rows={})
        logmu=f.numeric.outer.logmui;u0=self.iv.exp((self.iv.log(4)+logmu)/3)
        field.pulse.rows.update(saddle_L=self.iv.exp(-self.iv.mpf(2)/3*(self.iv.log(4)+logmu)),saddle_u0=u0)
        field.logE=self._decode(self.selection['positive_inner']['original_source_end_scale'])
        field.logEcap=logmu-1000;field.logE2cap=2*logmu-1000
        if _bounds(self.ctx,field.logEcap-field.logE)[0]<=0 or _bounds(self.ctx,field.logE2cap-2*field.logE)[0]<=0: raise ValueError('Actual end-scale cap required')
        field.Ecap=self.iv.mpf([0,_bounds(self.ctx,self.iv.exp(field.logEcap))[1]])
        field.E2cap=self.iv.mpf([0,_bounds(self.ctx,self.iv.exp(field.logE2cap))[1]])
        field.Xp=self.constants['Xp']
        field.full_end_energy_weights=[self.iv.exp(-2*self.mu*j)*self.gram for j in (-3,-1)]
        # Cache only pure Z-independent scalar kernels in a private copy of
        # the original method namespaces. No ancestor/module cache is edited.
        # The defining functions, arguments and directed results are unchanged.
        import lei_ren_part1_paper_compliant_axial_pulse_field as original
        namespace=dict(original.__dict__);self.kernel_caches={}
        def key(value):
            if hasattr(value,'_mpi_'): return ('interval',value._mpi_)
            if isinstance(value,(str,int,float,bool,type(None))): return value
            return ('context',id(value))
        def memo(name,fn):
            cache={};self.kernel_caches[name]=cache
            def call(*args,**kwargs):
                token=tuple(key(v) for v in args),tuple((n,key(v)) for n,v in sorted(kwargs.items()))
                if token not in cache:cache[token]=fn(*args,**kwargs)
                return cache[token]
            return call
        for name in ('gp','gp_energy','gp_future_energy','backward_bump_weights'):
            fn=getattr(original,name)
            clone=FunctionType(fn.__code__,namespace,fn.__name__,fn.__defaults__,fn.__closure__)
            namespace[name]=memo(name,clone)
        for name in ('main','_gap','end'):
            fn=getattr(CompliantAxialPulseField,name)
            setattr(field,name,MethodType(FunctionType(fn.__code__,namespace,fn.__name__,fn.__defaults__,fn.__closure__),field))
        self.function_source=PREFIX+'compliant_axial_pulse_field.py'
        self.hashes[self.function_source]=hashlib.sha256((self.directory/self.function_source).read_bytes()).hexdigest()
        name='src/openai_ns_reconstruction/original_native_pulse_numeric.py';self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _decode(self,row):
        return self.iv.mpf([mp.make_mpf(tuple(row[k+'_exact_mpf_tuple'])) for k in ('lower','upper')])

    def _jet(self,row): return self.J(self.iv,[self._decode(v) for v in row['coefficients']])

    def _rows(self,jet):
        ref=self.quiet.numeric.compact.numeric.patch.leading.reference
        return [ref._summary(sum(_bounds(self.ctx,v),self.ctx.mpf(0))/2,v) for v in jet.coefficients]

    def _row(self,value): return self._rows(self.J.constant(self.iv,value,0))[0]

    def _intersect(self,a,b):
        out=[]
        for x,y in zip(a.coefficients,b.coefficients):
            lo,hi=_bounds(self.ctx,x);ll,hh=_bounds(self.ctx,y)
            if max(lo,ll)>min(hi,hh): raise ArithmeticError('Actual quiet and selected native Rp frames are disjoint')
            out.append(self.iv.mpf([max(lo,ll),min(hi,hh)]))
        # Higher rows use the proved canonical native function identity.
        return self.J(self.iv,out+list(a.coefficients[len(out):]))

    def _data(self,Z):
        iv,J=self.iv,self.J;z=iv.mpf(Z);key=z._mpi_
        if key in self._data_cache: return self._data_cache[key]
        zl,zh=_bounds(self.ctx,z)
        label=next(n for n,ends in (('negative_outer',(-1,-.5)),('negative_inner',(-.5,0)),('positive_inner',(0,.5)),('positive_outer',(.5,1))) if ends[0]<=zl and zh<=ends[1])
        source=self.selection[label];ap=self._jet(source['refined_coefficients']['ap'])
        C=[self._jet(source['refined_coefficients'][n]) for n in ('C1','C2')]
        var=J.variable(iv,z,5);q=1+var*var;inverse=q.reciprocal();k=self.constants
        frame=dict(u=inverse*k['U'],m1=var*q*k['C1'],m2=var*q*k['C2'],
            X=J.constant(iv,k['Xp'],5),energy=var*var*q*q*k['C_E']+k['C0'],Mp=inverse*inverse*k['Pin'])
        if type(Z) is not str: raise ValueError('Original exact decimal Z must reach the pulse callback')
        Ztext=Z
        if Ztext not in self.inlets: self.inlets[Ztext]=self.quiet.numeric.evaluate(Ztext,'1')
        inlet=self.inlets[Ztext];actual=inlet['native_pulse_input_frame']
        checked={n:self._intersect(v,self.quiet.numeric._jet(actual[n],True)) for n,v in frame.items()}
        original_P0=self.quiet.numeric._jet(inlet['independent_P0_over_Pstar_squared'],True)
        if original_P0.order!=5:raise ValueError('All six original independent datum rows required, without zero padding')
        checked['P0']=self._intersect(original_P0,self.quiet.numeric._jet(actual['P0'],True))
        checked['pressure']=self._intersect(checked['P0']+checked['Mp'],self.quiet.numeric._jet(actual['pressure'],True))
        if set(checked)!={'u','m1','m2','X','energy','Mp','P0','pressure'} or any(v.order!=5 for v in checked.values()):
            raise ValueError('All eight true C5 current inlet functions required')
        self._frame_cache[key]=checked
        future=self._jet(source['source_inputs']['weighted_complete_future_Taylor'])*(2*iv.exp(26)/self.mu)
        selected=dict(selected_ap_Taylor=ap,selected_scaled_end_coefficient_Taylor=C,source_cell=label,complete_future=future)
        pin=dict(Mp_over_Pstar_squared=list(checked['Mp'].coefficients),P0_over_Pstar_squared=list(checked['P0'].coefficients))
        out=selected,pin,checked['u'],checked['energy'],[checked['m1'],checked['m2']]
        self._data_cache[key]=out
        return out

    def _future(self,Z): return dict(complete_future_energy_Taylor=self._data(Z)[0]['complete_future'])

    def _pressure(self,Z,t,inlet,u,log_parts=None):
        iv,J=self.iv,self.J
        if log_parts is None: log_parts=dict(time_term=-self.field.prate*t)
        decay=self.field.factor(sum(log_parts.values(),iv.mpf(0)))
        kernel=self.decay_integral(iv,self.field.prate,t) if _bounds(self.ctx,self.field.prate*t)[1]<iv.mpf('1e-20') else (1-decay)/self.field.prate
        Mp=J(iv,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)
        P0=J(iv,inlet['P0_over_Pstar_squared'])
        return dict(Mp_over_Pstar_squared=Mp,P0_over_Pstar_squared=P0,P_over_Pstar_squared=Mp+P0,
                    retained_swirl_pressure_decay_log_parts=log_parts,exact_pressure_source_unchanged=True)

    def evaluate(self,Z,chart,coordinate):
        c,iv,J=self.ctx,self.iv,self.J
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1: raise ValueError('Interior decimal Z required')
        if chart not in ('entrance','main','gap','gap_end','end') or type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)): raise ValueError('Original pulse chart and decimal coordinate required')
        key=Z,chart,coordinate
        if key in self._packets: return self._packets[key]
        selected,_,u,_,_=self._data(Z)
        with mp.workdps(iv.dps+30):
            if chart=='entrance': raw=self.field.entrance(Z,coordinate)
            elif chart=='main': raw=self.field.main(Z,coordinate,cells=self.cells)
            elif chart=='gap': raw=self.field.gap(Z,coordinate)
            elif chart=='gap_end': raw=self.field.gap_from_end(Z,coordinate)
            else: raw=self.field.end(Z,coordinate,cells=self.cells)
        z=J.variable(iv,iv.mpf(Z),4);d=1-z*z;L=1-z*z*self.delta
        q=1+z*z;thetaZ=-2*z/q
        B=raw.get('formal_scaled_Uz_over_Utheta',raw['Uz_over_Utheta'])
        m1,m2=raw.get('formal_scaled_Mz_mixed_moments',[raw['Mz_over_R_Utheta'],raw['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']])
        derivative=J(iv,[(i+1)*m1[i+1] for i in range(m1.order)])
        Q=(z*B.truncate(4)*2-z*m1.truncate(4)*(1-self.delta)-d*(derivative+thetaZ*m1.truncate(4)))/L
        source=self.inlets[Z]
        fields=dict(Uz_over_Utheta=B,Mz_over_R_Utheta=m1,Mtheta_z_over_sqrt2_R_3half_Utheta_squared=m2,
            Q_over_Utheta=Q,X=J.constant(iv,raw['Mtheta_over_sqrt2_R_3half_Utheta'],5),
            energy=raw['Mztheta_over_R_Utheta_squared'],pressure=raw['pressure']['P_over_Pstar_squared'],Mp=raw['pressure']['Mp_over_Pstar_squared'])
        if chart=='entrance': t=iv.mpf(coordinate);xi=self.mu*t;end_offset=iv.mpf(0)
        elif chart in ('main','gap'):xi=iv.mpf(coordinate);t=xi/self.mu;end_offset=iv.mpf(0)
        else:end_offset=iv.mpf(coordinate);xi=13+self.mu*end_offset;t=13/self.mu+end_offset
        # In end coordinates xi and t retain finite offsets as separate atoms.
        if chart in ('end','gap_end'):
            radial_log=[self._row(13/self.mu),self._row(end_offset)]
            angular_log=[self._row(-iv.mpf('6.5')/self.mu),self._row(-13),self._row(-end_offset/2),self._row(-self.mu*end_offset)]
            radial_mu_log=[self._row(-13),self._row(-self.mu*end_offset)]
            mass_log=[self._row(iv.mpf('6.5')/self.mu),self._row(-13),self._row(end_offset/2),self._row(-self.mu*end_offset)]
        else:
            radial_log=[self._row(t)];angular_log=[self._row(-t/2),self._row(-xi)];radial_mu_log=[self._row(-xi)]
            mass_log=[self._row(t/2),self._row(-xi)]
        lm=self.quiet.numeric.outer.logmui
        end_atoms=self._rows(J(iv,[-1/self.mu,-3*iv.exp(-iv.mpf(2)/3*(iv.log(4)+lm)),
            -iv.mpf(4)/3*lm,iv.mpf(2)/3*iv.log(4)-iv.log(6)/2])) if chart=='end' else []
        out=dict(Z=Z,native_pulse_chart=chart,native_pulse_coordinate=coordinate,source_family=self.family,
            selected_integer=source['selected_integer'],inlet_frame=source['native_pulse_input_frame'],
            native_selected_inlet_C5_source_callback=True,selected_source_cell=selected['source_cell'],
            native_callback_frame_C5={n:self._rows(v) for n,v in self._frame_cache[iv.mpf(Z)._mpi_].items()},
            native_u=self._rows(u),normalized_fields={n:self._rows(v) for n,v in fields.items()},
            selected_coefficients=dict(ap=self._rows(selected['selected_ap_Taylor']),end=[self._rows(v) for v in selected['selected_scaled_end_coefficient_Taylor']]),
            coordinate=dict(Rp_global_y_log_terms=copy.deepcopy(source['coordinate']['global_y_log_terms']),log_R_over_Rp_terms=radial_log),
            theta_radial_log_terms=angular_log,radial_mu_log_terms=radial_mu_log,end_source_log_terms=end_atoms,
            mass_radial_log_terms=mass_log,
            incoming_P1_at_Rp=copy.deepcopy(source['reserved_n1_pressure']),
            terminal_linear_zero_by_source_identity=raw.get('terminal_linear_zero_by_empty_future_supports',False),
            leading_pressure_includes_same_P0_once=True,original_source_method=self.function_source+'::CompliantAxialPulseField.'+chart,
            limits=dict(numerical_velocity_and_leading_pressure=True,interior_P1=False,
                signed_coefficients_are_function_enclosures=True,physical_gradients=False,relative_pressure_cancellation_accuracy=False,
                exact_heat_exterior=False,stress_cone=False,common_global_N=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out


class OriginalNativePulseVelocity:
    def __init__(self, **kwargs):
        self.numeric=OriginalNativePulseNumeric(**kwargs);f=self.numeric.quiet
        self.ctx,self.iv,self.J,self.amplitude=f.ctx,f.iv,f.J,f.amplitude
        self.family=f.family;self.hashes,self.code_hashes=dict(self.numeric.hashes),dict(self.numeric.code_hashes)
        self._summary,self._pair,self._restore,self._rotate=f._summary,f._pair,f._restore,f._rotate

    def evaluate(self,Z,chart,coordinate,*,log_tau='-1',scaled_time=None,azimuth='0'):
        p=self.numeric.evaluate(Z,chart,coordinate);c,iv,J=self.ctx,self.iv,self.J
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)): raise ValueError('Finite decimal azimuth required')
        ell,de,time=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        def scale(v,k):return v[0]*k,v[1]*k
        def logs(a,b,P=0,Rp=0,constant=None):
            out=[self._summary(scale(v,a)) for v in ell if a]+[self._summary(scale(v,b)) for v in de if b]
            if P:out.append(self._summary(scale((self.amplitude.lp,self.amplitude.lp_iv),P)))
            if Rp:
                out.extend(self._summary(scale(self._pair(v),c.mpf(Rp))) for v in p['coordinate']['Rp_global_y_log_terms'])
                out.append(self._summary((c.mpf(Rp)*c.log(110),iv.mpf(Rp)*iv.log(110))))
            if constant is not None:out.append(self._summary(constant))
            return out
        def product(a,b):
            aa,bb=self._pair(a),self._pair(b)
            return self._summary((aa[0]*bb[0],aa[1]*bb[1]))
        def term(name,row,weights,order=0):return dict(name=name,temporal_order=order,coefficient=copy.deepcopy(row),log_scale_terms=copy.deepcopy(weights))
        fields=p['normalized_fields'];u=p['native_u'][0];end=p['end_source_log_terms']
        zero=self._summary((c.mpf(0),iv.mpf(0)))
        swirl=[term('pulse_leading_swirl',u,logs(-1,-1,P=1)+p['theta_radial_log_terms']),term('n1_pulse_swirl_zero',zero,[],1)]
        axial=[term('pulse_leading_axial',product(u,fields['Uz_over_Utheta'][0]),logs(-1,-1,P=1)+p['theta_radial_log_terms']+end),term('n1_pulse_axial_zero',zero,[],1)]
        radial=[term('pulse_leading_radial',product(u,fields['Q_over_Utheta'][0]),logs(-1,0,P=1,Rp='.5',constant=(-c.log(2)/2,-iv.log(2)/2))+p['radial_mu_log_terms']+end),term('n1_pulse_radial_zero',zero,[],1)]
        stream=[term('pulse_leading_streamfunction',product(u,fields['Mz_over_R_Utheta'][0]),logs(1,-1,P=1,Rp=1)+p['mass_radial_log_terms']+end),term('n1_pulse_mass_zero',zero,[],1)]
        radius=term('physical_radius',self._summary((c.mpf(1),iv.mpf(1))),logs(1,0,Rp='.5',constant=(c.log(2)/2,iv.log(2)/2))+
                    [self._summary(scale(self._pair(v),c.mpf('.5'))) for v in p['coordinate']['log_R_over_Rp_terms']],None)
        pressure=[term('pulse_leading_pressure',fields['pressure'][0],logs(-2,-2,P=2))]
        angle,anglei=c.mpf(azimuth),iv.mpf(azimuth);cs,sn=(c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))
        physical=dict(source_family=self.family,time=time,native=dict(chart=chart,Z=Z,coordinate=coordinate,
            physical_viscosity=1,R_definition='Rp*exp(t); t=xi/mu or13/mu+s'),
            coordinates=dict(r=radius,x=self._rotate(radius,cs),y=self._rotate(radius,sn),
                z=term('physical_z',self._summary((c.mpf(Z),iv.mpf(Z))),logs(1,-1),None),
                t=dict(kind='one_minus_exp_of_log_tau',definition=time)),
            cylindrical=dict(ur=radial,utheta=swirl,uz=axial),Cartesian=dict(
                u=[self._rotate(v,cs) for v in radial]+[self._rotate(v,(-sn[0],-sn[1])) for v in swirl],
                v=[self._rotate(v,sn) for v in radial]+[self._rotate(v,cs) for v in swirl],w=copy.deepcopy(axial)),
            streamfunction=stream,pressure=pressure,coefficient_orders_included=[0,1],pressure_orders_included=[0],
            full_n1_pressure_included=False,interior_P1_not_substituted_by_zero=True,
            leading_pressure_includes_same_P0_once=True,n1_velocity_zero_from_source_identity=True,
            end_scale_kept_as_original_atomic_log_recipe=chart=='end',radial_half_power_cancelled_before_boxing=True,
            shared_streamfunction_meridional_recovery=True,physical_hyper_exponentials_materialized=False,
            global_field_available=False,genuine_n_ge_2_temporal_recursion=False,full_corrected_NS_solution=False)
        out=copy.deepcopy(p);out['physical']=physical
        return out
