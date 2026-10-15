"""Original selected-N second switch and post-power velocities to exact R110.

Six unmodified own histories are evaluated separately from the five
finite-N correction bounds. Native chart views share only source caches.
"""
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_collar_velocity import _InletGeometry
from .original_first_switch_numeric import OriginalFirstSwitchNumeric
from .original_macro_numeric import OriginalMacroNumeric
from .original_second_micro_numeric import OriginalSecondMicroNumeric


RATES=dict(H=2,Mpsi=1,Kpsi=2,Apsi=1,B=2,C=1)


def post_power_histories(ctx,initial,phi,psi,y):
    """Exact positive-kernel solution in original normalized moment units."""
    e1=ctx.exp(-y);e2=ctx.exp(-2*y)
    ephi=ctx.exp(-ctx.mpf(2)*y/5);ephi2=ctx.exp(-ctx.mpf(4)*y/5)
    whk=ctx.mpf(5)/4*(ephi-e2);wb=ctx.mpf(5)/6*(ephi2-e2)
    wc=5*(ephi2-e1)
    return dict(H=initial['H']*e2+phi*whk,
        Mpsi=initial['Mpsi']*e1+psi*(1-e1),
        Kpsi=initial['Kpsi']*e2+phi*psi*whk,
        Apsi=initial['Apsi']*e1+psi*psi*(1-e1),
        B=initial['B']*e2+phi*phi*wb,C=initial['C']*e1+phi*phi*wc)


class OriginalSecondSwitchNumeric(OriginalFirstSwitchNumeric):
    chart='second_switch'
    log_radius_definition='Y+hbS*(1+s)'
    radius_definition='100*exp(hbS*(1+s))'

    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        path=self.leading.core.directory/(PREFIX+'compliant_current_original_first_switch_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted actual original first-switch numerical inlet required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed accepted first-switch numerical source: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        iv,c=self.transport.iv,self.transport.ctx
        self.post_log=iv.ln(iv.mpf(11)/10);self.post_log_point=c.ln(c.mpf(11)/10)
        self.post_length=self.post_log-iv.mpf([0,_bounds(c,2*self.leading.hcap)[1]])
        if _bounds(c,self.post_length)[0]<=0:raise ValueError('Positive original post-power length required')
        self._second_switch_background={};self._second_switch_increments={};self._post_power_background={}
        self._six_cache={};self._post_power_view=None
        relative='src/openai_ns_reconstruction/original_r110_numeric.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()

    def post_power_field(self):
        """Independent chart cache, same immutable source and unmodified inlets."""
        if self._post_power_view is None:
            view=object.__new__(OriginalPostPowerNumeric)
            view.__dict__=dict(self.__dict__)
            view._micro_cache={};view._upgraded=set();view._post_power_view=None
            self._post_power_view=view
        return self._post_power_view

    def post_power(self,Z,q,**kwargs):
        return self.post_power_field().evaluate(Z,q,**kwargs)

    def _coordinate(self,coordinate,selected):
        c,iv=self.transport.ctx,self.transport.iv
        if selected or type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Finite decimal '+self.chart+' coordinate in[0,1] required')
        s=iv.mpf(coordinate);lo,hi=_bounds(c,s)
        if lo<0 or hi>1:raise ValueError(self.chart+' coordinate in[0,1] required')
        return s,False,dict(kind='exact_decimal_'+self.chart+'_coordinate',value=coordinate)

    def _covers(self,Z,s):
        c,iv=self.transport.ctx,self.transport.iv
        zl,zh=_bounds(c,iv.mpf(Z));sl,sh=_bounds(c,s)
        local=[];prefix=[];frames=[]
        order=('first_micro','second_micro','frozen_macro','first_switch','second_switch','post_power')
        for frame in self.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            frames.append(frame['exact_Z_cell'])
            for chart in order[:order.index(self.chart)+1]:
                current=next(w for w in frame['actual_17_chart_C3_transports'] if w['chart']==chart)
                for cell in current['cells']:
                    left,right=self._endpoint(cell['exact_left']),self._endpoint(cell['exact_right'])
                    ll,lh=_bounds(c,left);rl,rh=_bounds(c,right)
                    if chart==self.chart and not rh<sl and not ll>sh:local.append(cell)
                    end=min(sh,rh) if chart==self.chart else rh
                    if end>ll:
                        width=iv.mpf(_bounds(c,iv.mpf(end)-iv.mpf(ll))[1])
                        logmeasure=(iv.ln(self.Y) if chart=='frozen_macro' else
                                    iv.ln(self.post_log) if chart=='post_power' else self.leading.logh)
                        prefix.append((cell,width,logmeasure,chart))
        if not local:raise ValueError('Actual complete '+self.chart+' source cover required')
        return local,prefix,frames

    def _second_increment(self,Z,s):
        key=Z,s._mpi_
        if key in self._second_switch_increments:return self._second_switch_increments[key]
        o,t=self.leading,self.transport;b=o.bridge;c=b.ctx
        p=o._prepared[Z]
        module=__import__(PREFIX+'compliant_actual_bridge_integrals',fromlist=['jet_norm'])
        with mp.workdps(280):
            phase=c.mpf(_bounds(t.ctx,s)[1]);complement,sigma=b.pulse_masses(phase)
            norm=sum((module.jet_norm(row,b.weight)
                *module.upper(c,abs(p['R0']**j*c.mpf(100)**(1-j)))
                *c.exp(abs(1-j)*b.cap*(1+phase))
                for j,row in enumerate(p['modes']['D_over_R'])),c.mpf(0))
            shear=o._width_error(o._width_error(o._transfer(norm*complement/2)))
            power=o._width_error(o._transfer(c.mpf(2)*sigma/5))
        out=shear+power,dict(angular_shear_weight='1-sigma(s)',angular_power_weight='sigma(s)',
            native_shear_width_power=2,native_power_width_power=1,
            unmodified_U_constant=True,local_selected_N_B_still_evolves=True,
            complement_mass_exact_mpf_tuples=[list(v) for v in complement._mpi_],
            sigma_mass_exact_mpf_tuples=[list(v) for v in sigma._mpi_],
            physical_measure='hbS ds once',current_comparison_modes_retained=True)
        self._second_switch_increments[key]=out
        return out

    def _second_geometry(self,s):
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
        tiny=o._width_error(1+s)
        rho=iv.exp(iv.ln(4)+self.Y+iv.mpf([0,_bounds(c,tiny)[1]]))
        return rho,self.Y,tiny

    def _downstream_geometry(self,s):return self._second_geometry(s)

    def _point_log_radius(self,s):return self.Ypoint

    def _background(self,Z,s):
        key=Z,s._mpi_
        if key in self._second_switch_background:return self._second_switch_background[key]
        o,iv=self.leading,self.transport.iv
        point,inlet,_=OriginalFirstSwitchNumeric._background(self,Z,iv.mpf(1))
        ell,_=self._second_increment(Z,s)
        phi_error=o._norm(inlet['Phi'])*ell*iv.exp(ell)
        truth=dict(Phi=inlet['Phi']+o._norm_error_jet(phi_error,3),Psi=inlet['Psi'],
            Mpsi=inlet['Mpsi']+o._norm_error_jet(o._width_error(s*o._norm(inlet['Psi']-inlet['Mpsi'])),3))
        rho,_,_=self._second_geometry(s)
        self._second_switch_background[key]=point,truth,rho
        return point,truth,rho

    def _targets(self,values):
        phi,psi=values['Phi'],values['Psi']
        return dict(H=phi,Mpsi=psi,Kpsi=phi*psi,Apsi=psi*psi,B=phi*phi/2,C=phi*phi)

    def _target_errors(self,values,ephi,epsi):
        o=self.leading;np,nv=o._norm(values['Phi']),o._norm(values['Psi'])
        square=2*np*ephi+ephi*ephi
        return dict(H=ephi,Mpsi=epsi,Kpsi=np*epsi+nv*ephi+ephi*epsi,
                    Apsi=2*nv*epsi+epsi*epsi,B=square/2,C=square)

    def _r100_six(self,Z):
        key=Z,'R100'
        if key in self._six_cache:return self._six_cache[key]
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
        core=o.inlet_jets(Z,degree=24,order=3)
        point,truth={},{};a=core['truth']['values'];b=core['point']['values']
        p=o._prepared[Z]
        ell=o._width_error(o._transfer(p['angular_prefix_norm_per_hb']))
        ephi=o._norm(a['Phi'])*ell*iv.exp(ell)
        epsi=OriginalMacroNumeric._uniform_Psi_increment(self,Z,p)
        errors=self._target_errors(a,ephi,epsi);tp,tt=self._targets(b),self._targets(a)
        for name,rate in RATES.items():
            decayp=c.exp(-rate*self.Ypoint);decay=iv.exp(-rate*self.Y)
            point[name]=b[name]*decayp+tp[name]*(1-decayp)
            truth[name]=a[name]*decay+tt[name]*(1-decay)+o._norm_error_jet(errors[name]*(1-decay),3)
        self._six_cache[key]=point,truth
        return point,truth

    def _short_six(self,initial,values,ephi,epsi,s):
        o=self.leading;point,truth=initial;targets=self._targets(values)
        errors=self._target_errors(values,ephi,epsi)
        # For I'=rate*(target-I), subtract the actual incoming I before
        # bounding the positive kernel by1. No comparison moments enter.
        out={name:truth[name]+o._norm_error_jet(o._width_error(
            rate*s*(o._norm(targets[name]-truth[name])+errors[name])),3)
            for name,rate in RATES.items()}
        return point,out

    def _first_six(self,Z):
        key=Z,'R1'
        if key not in self._six_cache:
            o,iv=self.leading,self.transport.iv
            _,incoming,_=OriginalMacroNumeric._background(self,Z,iv.mpf(1))
            ell,epsi,_=self._switch_increments(Z,iv.mpf(1))
            ephi=o._norm(incoming['Phi'])*ell*iv.exp(ell)
            self._six_cache[key]=self._short_six(self._r100_six(Z),incoming,ephi,epsi,iv.mpf(1))
        return self._six_cache[key]

    def _second_six(self,Z,s):
        key=Z,'second',s._mpi_
        if key not in self._six_cache:
            o,iv=self.leading,self.transport.iv
            _,incoming,_=OriginalFirstSwitchNumeric._background(self,Z,iv.mpf(1))
            ell,_=self._second_increment(Z,s)
            ephi=o._norm(incoming['Phi'])*ell*iv.exp(ell)
            self._six_cache[key]=self._short_six(self._first_six(Z),incoming,ephi,iv.mpf(0),s)
        return self._six_cache[key]

    def _own_six(self,Z,s):return self._second_six(Z,s)

    def _numeric(self,Z,coordinate,selected):
        numeric=OriginalSecondMicroNumeric._numeric(self,Z,coordinate,selected)
        packet=numeric[0];s=numeric[4];post=self.chart=='post_power'
        charts=['first_micro','second_micro','frozen_macro','first_switch','second_switch']
        if post:charts.append('post_power')
        packet['source_cover'].update(native_source_domain=['0','1'],local_source_chart=self.chart,
            cumulative_history_charts=charts,true_phase='fractional_part(original_N*('+self.log_radius_definition+'-hbB*original_sc/2))')
        packet['background_control_contract']=(dict(unmodified_U_constant=True,
            angular_log_slope='-2/5',exact_length='log(11/10)-2hbS',physical_measure='exact_length dq once')
            if post else self._second_increment(Z,s)[1])
        packet['history_contract'].update(measure='micro/first/second hb ds; macro(Y-2hbB)dq; post(log(11/10)-2hbS)dq; each once',
            partial_prefix='All preceding complete charts plus current native prefix; positive kernels bounded by1 and physical lengths bounded above only inside error caps',
            all_nonzero_incoming_five_corrections_retained=True)
        packet['own_mean_contract']=dict(equation='dMpsi/dy=Psi-Mpsi',
            unmodified_Psi_fixed_at_actual_first_switch_endpoint=True,
            exact_source='exp(-y)*Mpsi_in+(1-exp(-y))*Psi_in',
            y='(log(11/10)-2hbS)*q' if post else 'hbS*s',
            local_selected_N_B_and_cumulative_m_correction_not_frozen=True)
        packet['limits'].pop('actual_modified_second_micro_velocity_Z3',None)
        packet['limits'].update(actual_modified_second_switch_velocity_Z3=not post,
            actual_modified_post_power_velocity_Z3=post,R110_endpoint_reached=post and _bounds(self.transport.ctx,s)==(1,1),
            all_later_charts=False,unmodified_own_six_Z2_available=True)
        pp,tt=self._own_six(Z,s);o=self.leading
        packet['unmodified_own_six']={name:o._jet_summary(pp[name].truncate(2),tt[name].truncate(2)) for name in RATES}
        packet['unmodified_six_contract']=dict(rates=RATES,source='actual core atoms, full finite bridge source, actual short-switch and exact post-power histories',
            raw_recovery=dict(H='H',M='U0+epsilon*Mpsi',K='U0*H+epsilon*Kpsi',
                A='U0^2+2*U0*epsilon*Mpsi+epsilon^2*Apsi',B='B',C='C'),
            selected_N_five_histories_remain_separate=True,comparison_moments_not_substituted=True,
            terminal_five_matching_not_claimed=True,pressure_not_exported=True)
        return numeric

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        numeric=self._numeric(Z,s,False)
        return dict(**numeric[0],physical=self._physical(numeric,Z,log_tau,scaled_time,azimuth))


    def _physical(self,numeric,Z,log_tau,scaled_time,azimuth):
        packet,point,truth,rho,s=numeric;c,iv=self.ctx,self.iv;o,t=self.leading,self.transport
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):raise ValueError('Finite decimal azimuth required')
        rows=packet['normalized_modified_leading'];z,zi=c.mpf(Z),iv.mpf(Z)
        de,di=self.amplitude.delta,self.amplitude.delta_iv
        u0,u0i=4*z+self.amplitude.j,4*zi+self.amplitude.j_iv
        profiles={name:self._pair(rows[source][0]) for name,source in
                  (('Phi','Phi'),('Psi','Psi'),('mean_Psi','Mpsi'),('Qpsi','Qpsi'))}
        profiles.update(U0=(u0,u0i),Qaxis=(((1+de)*z*u0-4*(1-z*z))/(1-de*z*z),
            ((1+di)*zi*u0i-4*(1-zi*zi))/(1-di*zi*zi)))
        if Z not in self._terminal_cache:self._terminal_cache[Z]=t.lower_collar(Z,'2')
        terminal=self._terminal_cache[Z]['normalized_own_primitive']['actual_primitive'];J=t.Jet
        pp=J(t.ctx,[t.ctx.mpf(row['approximation']) for row in terminal])
        pt=J(t.iv,[t.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]) for row in terminal])
        qp=t.ctx.mpf(packet['coordinate']['value'])
        rhop=4*t.ctx.exp(self._point_log_radius(s))
        meanp,meant=pp/rhop,pt/rho
        def recover(mean,ctx,delta):
            zz=J.variable(ctx,ctx.mpf(Z),2)
            return (-zz*mean*(1+delta)-(1-zz*zz)*o._derivative(mean))/(1-zz*zz*delta)
        q1p,q1t=recover(meanp,t.ctx,o.core.delta),recover(meant,t.iv,o.core.delta_iv)
        transfer=lambda p,b:(c.mpf(p),iv.mpf([mp.make_mpf(tuple(v)) for v in b._mpi_]))
        n1=dict(f=(c.mpf(0),iv.mpf(0)),u=(c.mpf(0),iv.mpf(0)),
            mean_u=transfer(meanp[0],meant[0]),q=transfer(q1p[0],q1t[0]))
        angle=c.mpf(azimuth);anglei=iv.mpf(azimuth)
        adapter=_InletGeometry(self,profiles,n1)
        out=adapter._assemble(Z,'4',log_tau,scaled_time,((c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))))
        _,base,tiny=self._downstream_geometry(s)
        def geometry(multiplier):
            value=base*t.iv.mpf(multiplier)
            bound=tiny*t.iv.mpf(multiplier)
            bounds=iv.mpf([mp.make_mpf(v) for v in (value+t.iv.mpf([0,_bounds(t.ctx,bound)[1]]))._mpi_])
            nominal=c.mpf(self._point_log_radius(s))*c.mpf(multiplier)
            row=self._summary((nominal,bounds))
            row.update(kind='original_'+self.chart+'_radius_log',definition='multiplier*('+self.log_radius_definition+')',
                multiplier=multiplier,coordinate=packet['coordinate'],original_log_hb=self._log_width(),
                original_hb_not_replaced_by_cap=True)
            return row
        half,whole=geometry('.5'),geometry('1')
        for key in ('ur','utheta'):
            for term in out['cylindrical'][key]:term['log_scale_terms'].append(half)
        for key in ('u','v'):
            for term in out['Cartesian'][key]:term['log_scale_terms'].append(half)
        for key in ('r','x','y'):out['coordinates'][key]['log_scale_terms'].append(half)
        for term in out['streamfunction']:term['log_scale_terms'].append(whole)
        out['native']=dict(chart='actual_modified_'+self.chart,Z=Z,phase_recipe=packet['coordinate'],
            rho_definition='4exp('+self.log_radius_definition+')',R_definition=self.radius_definition,
            left_endpoint_exact=_bounds(t.ctx,s)==(0,0),
            right_endpoint_exact=_bounds(t.ctx,s)==(1,1),
            R110_endpoint_exact=self.chart=='post_power' and _bounds(t.ctx,s)==(1,1),
            coefficient_geometry='rho4 physical weights plus full current radius log; directed n1 mean uses actual current rho, point omits only the bounded original width contribution')
        out['original_selected_integer']=packet['selected_integer']
        out['actual_modified_leading_errors_included']=True
        out['n1_own_history']=dict(actual_terminal_primitive=self._terminal_cache[Z]['normalized_own_primitive'],
            mean_u=self._summary(n1['mean_u']),q=self._summary(n1['q']),
            primitive_retained_through_current_radius=True,f_and_u_exact_zero=True,
            current_rho_not_rho4_used=True,directed_mean_uses_actual_current_rho=True,
            point_mean_radius_proxy='4exp(Y) in second_switch;4exp(Y+log(11/10)*q) in post_power; omitted hbS contribution fully enclosed',
            point_proxy_does_not_select_original_hb=True,
            no_later_reserved_repair_support_on_current_chart=True)
        out['physical_gradient_or_vorticity_available']=False
        return out


class OriginalPostPowerNumeric(OriginalSecondSwitchNumeric):
    chart='post_power'
    log_radius_definition='Y+log(11/10)*q+2hbS*(1-q)'
    radius_definition='100*exp(2hbS+(log(11/10)-2hbS)*q)'

    def _downstream_geometry(self,q):
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
        base=self.Y+self.post_log*q;tiny=o._width_error(2*(1-q))
        rho=iv.exp(iv.ln(4)+base+iv.mpf([0,_bounds(c,tiny)[1]]))
        return rho,base,tiny

    def _point_log_radius(self,q):
        return self.Ypoint+self.post_log_point*self.transport.ctx.make_mpf(q._mpi_[0])

    def _background(self,Z,q):
        key=Z,q._mpi_
        if key in self._post_power_background:return self._post_power_background[key]
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
        pp,tt,_=OriginalSecondSwitchNumeric._background(self,Z,iv.mpf(1))
        y=self.post_length*q;yp=self.post_log_point*c.make_mpf(q._mpi_[0])
        point=dict(Phi=pp['Phi']*c.exp(-c.mpf(2)*yp/5),Psi=pp['Psi'],
            Mpsi=pp['Mpsi']*c.exp(-yp)+pp['Psi']*(1-c.exp(-yp)))
        decay=iv.exp(-y)
        truth=dict(Phi=tt['Phi']*iv.exp(-iv.mpf(2)*y/5),Psi=tt['Psi'],
            Mpsi=tt['Mpsi']*decay+tt['Psi']*(1-decay))
        rho,_,_=self._downstream_geometry(q)
        self._post_power_background[key]=point,truth,rho
        return point,truth,rho

    def _own_six(self,Z,q):
        key=Z,'post',q._mpi_
        if key not in self._six_cache:
            o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
            hp,ht=self._second_six(Z,iv.mpf(1))
            pp,tt,_=OriginalSecondSwitchNumeric._background(self,Z,iv.mpf(1))
            yp=self.post_log_point*c.make_mpf(q._mpi_[0]);y=self.post_length*q
            self._six_cache[key]=(post_power_histories(c,hp,pp['Phi'],pp['Psi'],yp),
                post_power_histories(iv,ht,tt['Phi'],tt['Psi'],y))
        return self._six_cache[key]
