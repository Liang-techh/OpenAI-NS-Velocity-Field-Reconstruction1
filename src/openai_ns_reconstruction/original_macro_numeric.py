"""Actual selected-N macro values with true radius and own-mean relaxation.

Only the comparison profile is frozen. Numerical core proxies carry the
full actual finite-width source error and the nonzero modified histories.
"""
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_collar_velocity import _InletGeometry
from .original_second_micro_numeric import OriginalSecondMicroNumeric


class OriginalMacroNumeric(OriginalSecondMicroNumeric):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        path=self.leading.core.directory/(PREFIX+'compliant_current_original_second_micro_numeric.json')
        report=json.loads(path.read_text(encoding='utf-8'))
        if not report['checks']['passed'] or report['source_family']!=self.family:
            raise ValueError('Accepted actual second-micro numerical inlet required')
        for name,digest in report['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed accepted second-micro numerical source: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self.Y=self.leading.core.ll_iv+self.transport.iv.ln(25)
        self.Ypoint=self.leading.core.logLambda+self.transport.ctx.ln(25)
        self._macro_background={};self._macro_uniform={}
        relative='src/openai_ns_reconstruction/original_macro_numeric.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()

    def _coordinate(self,coordinate,selected):
        c,iv=self.transport.ctx,self.transport.iv
        if selected or type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Finite decimal macro fraction in[0,1] required')
        q=iv.mpf(coordinate);lo,hi=_bounds(c,q)
        if lo<0 or hi>1:raise ValueError('Macro fraction in[0,1] required')
        return q,False,dict(kind='exact_decimal_macro_fraction',value=coordinate)

    def _covers(self,Z,q):
        c,iv=self.transport.ctx,self.transport.iv
        zl,zh=_bounds(c,iv.mpf(Z));ql,qh=_bounds(c,q)
        local=[];prefix=[];frames=[]
        for frame in self.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            frames.append(frame['exact_Z_cell'])
            for chart in ('first_micro','second_micro','frozen_macro'):
                current=next(w for w in frame['actual_17_chart_C3_transports'] if w['chart']==chart)
                for cell in current['cells']:
                    left,right=self._endpoint(cell['exact_left']),self._endpoint(cell['exact_right'])
                    ll,lh=_bounds(c,left);rl,rh=_bounds(c,right)
                    if chart=='frozen_macro' and not rh<ql and not ll>qh:local.append(cell)
                    end=min(qh,rh) if chart=='frozen_macro' else rh
                    if end>ll:
                        width=iv.mpf(_bounds(c,iv.mpf(end)-iv.mpf(ll))[1])
                        # Lmacro=Y-2hb <=Y; this is an error-bound measure,
                        # not a replacement of the defining physical length.
                        logmeasure=iv.ln(self.Y) if chart=='frozen_macro' else self.leading.logh
                        prefix.append((cell,width,logmeasure,chart))
        if not local:raise ValueError('Actual frozen-macro complete source cover required')
        return local,prefix,frames

    def _correction_logs(self,local,prefix,flat):
        if not any(chart=='frozen_macro' for cell,width,logmeasure,chart in prefix):
            return super()._correction_logs(local,[(cell,width) for cell,width,logmeasure,chart in prefix],flat)
        iv,o=self.transport.iv,self.leading
        A=[];B=[];hist={name:[] for name in ('m','h','k','e','p')}
        for k in range(4):
            aa=self._sum_logs([self._source_log(cell['actual_source_and_primitive_range']['A'][k]) for cell in local])
            bb=self._sum_logs([self._source_log(cell['actual_source_and_primitive_range']['B'][k]) for cell in local])
            A.append(None if aa is None else aa-self.logN)
            B.append(None if bb is None else bb-self.logN+o.core.lp_iv+o.core.ll_iv)
            for name in hist:
                logs=[]
                for cell,width,logmeasure,chart in prefix:
                    v=self._source_log(cell['actual_density_ranges']['N_scaled'][name][k])
                    if v is not None:logs.append(v+iv.ln(width)+logmeasure)
                bound=self._sum_logs(logs)
                hist[name].append(None if bound is None else bound-self.logN)
        return A,B,hist

    def _macro_increments(self,prepared,q,order):
        o=self.leading;iv,J=o.iv,o.Jet
        with mp.workdps(280):angular,axial=o.bridge.contributions(prepared,'macro',q)
        Bell=iv.mpf(0);Bpsi=iv.mpf(0)
        for term in angular:
            row=J(iv,[o._scaled(v,term) for v in term['normalized_axial_coefficients'][:order+1]])
            Bell+=o._norm(row)
        for terms in axial.values():
            for term in terms:
                row=J(iv,[o._scaled(v,term,divide_eps=True) for v in term['normalized_axial_coefficients'][:order+1]])
                error=o._scaled(term['nonlinear_prefix_error_weighted_norm_per_next_hb'],term,
                    extra_power=1,divide_eps=True,symmetric=True)
                Bpsi+=o._norm(row)+o._abs(error)
        return Bell,Bpsi

    def _uniform_Psi_increment(self,Z,prepared):
        """hb times an absolute-modal norm, uniform over every macro prefix.

        The signed endpoint integral alone cannot bound earlier source
        values. Reproduce the original positive kernel majorant before
        applying epsilon normalization and the original width log.
        """
        if Z in self._macro_uniform:return self._macro_uniform[Z]
        o=self.leading;b=o.bridge;c=b.ctx
        module=__import__(PREFIX+'compliant_actual_bridge_integrals',fromlist=['jet_norm'])
        p=prepared;total=o.iv.mpf(0)
        with mp.workdps(280):
            for name,power,scale in (('hydro',1,c.mpf(0)),('pressure',1,b.pressure_log),('swirl',2,b.swirl_log)):
                n=(module.jet_norm(p['micros'][0]['quotient']*p['micros'][0]['direction']['drive_'+name],b.weight)
                    *module.upper(c,p['R0']**power)*(1+b.cap)
                    +module.jet_norm(p['micros'][1]['quotient']*p['micros'][1]['direction']['drive_'+name],b.weight)
                    *module.upper(c,p['R0']**power)*b.cap
                    +module.jet_norm(p['quotient'],b.weight)*module.weighted_mode_norm(p['modes']['drive_'+name],
                        module.radius_kernels(c,p['R0'],c.mpf(100),p['L'],power),b.weight))*c.exp(6*b.cap)
                term=b.source_term(p['data']['phi0']*0,1,scale)
                total+=o._abs(o._scaled(n,term,divide_eps=True))
        result=o._width_error(total)
        self._macro_uniform[Z]=result
        return result

    def _geometry(self,q):
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
        # Only hb is bounded by a small positive cover. Y is a directed
        # original analytic parameter; the macroscopic geometry is retained.
        length=self.Y-iv.mpf([0,_bounds(c,2*o.hcap)[1]])
        base=self.Y*q
        tiny=o._width_error(2*(1-q))
        logrho=iv.ln(4)+base+iv.mpf([0,_bounds(c,tiny)[1]])
        rho=iv.exp(logrho)
        theta=iv.exp(-length*q)
        if not q:theta=iv.mpf(1)
        mass=1-theta
        ml,mh=_bounds(c,mass);mass=iv.mpf([max(c.mpf(0),ml),min(c.mpf(1),mh)])
        return rho,theta,mass,base,tiny

    def _background(self,Z,q):
        key=Z,q._mpi_
        if key in self._macro_background:return self._macro_background[key]
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv;order=3
        bp,micro,_=OriginalSecondMicroNumeric._background(self,Z,iv.mpf(2))
        a=o.inlet_jets(Z,degree=24,order=order)['truth']['values']
        prepared=o._prepared[Z]
        Bell,Bpsi=self._macro_increments(prepared,_bounds(c,q)[1],order)
        Bphi=o._norm(a['Phi'])*Bell*iv.exp(o._width_error(Bell))
        truth={name:a[name]+o._norm_error_jet(o._width_error(bound),order)
               for name,bound in (('Phi',Bphi),('Psi',Bpsi))}
        rho,theta,mass,_,_=self._geometry(q)
        point=dict(Phi=bp['Phi'],Psi=bp['Psi'])
        tp=c.exp(-self.Ypoint*c.make_mpf(q._mpi_[0]))
        point['Mpsi']=bp['Mpsi']*tp+bp['Psi']*(1-tp)
        # M'=Psi-M. The nonzero true micro inlet is multiplied by its
        # actual decay; the full signed velocity source remains enclosed.
        truth['Mpsi']=micro['Mpsi']*theta+a['Psi']*mass\
            +o._norm_error_jet(self._uniform_Psi_increment(Z,prepared)*mass,order)
        self._macro_background[key]=point,truth,rho
        return point,truth,rho

    def _numeric(self,Z,coordinate,selected):
        numeric=super()._numeric(Z,coordinate,selected)
        packet=numeric[0]
        packet['source_cover'].update(native_source_domain=['0','1'],local_source_chart='frozen_macro',
            cumulative_history_charts=['first_micro','second_micro','frozen_macro'],
            true_phase='fractional_part(original_N*(2hb+(Y-2hb)*q-hb*original_sc/2))')
        packet['background_control_contract'].update(
            macro_control='original_hb, physical dy=(Y-2hb)dq; actual fields continue evolving',
            physical_measure='micro hb ds; macro (Y-2hb)dq, each once',
            only_comparison_profiles_frozen=True,actual_mean_incoming_retained=True)
        packet['history_contract'].update(measure='micro hb ds plus macro (Y-2hb)dq once',
            partial_prefix='Both complete micro windows plus macro0..q; positive source kernels bounded by1, macro length bounded above byY only in error caps')
        packet['limits'].pop('actual_modified_second_micro_velocity_Z3',None)
        endpoint=_bounds(self.transport.ctx,numeric[4])==(1,1)
        packet['limits'].update(actual_modified_macro_velocity_Z3=True,actual_radius_R100_reached=endpoint,
            all_later_charts=False)
        packet['own_mean_contract']=dict(equation='dMpsi/dy=Psi-Mpsi',
            exact_length='Y-2hb',Y='log(25)+logLambda',
            baseline='Mpsi_micro2*exp(-(Y-2hb)q)+Psi_core*(1-exp(-(Y-2hb)q))',
            actual_source_feedback='Integral of the true Psi-Psi_core, bounded by uniform absolute-modal norm times1-decay',
            comparison_mean_not_substituted=True,macro_density_not_multiplied_by_hb=True)
        return numeric

    def _physical(self,numeric,Z,log_tau,scaled_time,azimuth):
        packet,point,truth,rho,q=numeric;c,iv=self.ctx,self.iv;o,t=self.leading,self.transport
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
        rhop=4*t.ctx.exp(self.Ypoint*qp)
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
        _,_,_,base,tiny=self._geometry(q)
        def geometry(multiplier):
            value=base*t.iv.mpf(multiplier)
            bound=tiny*t.iv.mpf(multiplier)
            bounds=iv.mpf([mp.make_mpf(v) for v in (value+t.iv.mpf([0,_bounds(t.ctx,bound)[1]]))._mpi_])
            nominal=c.mpf(self.Ypoint*qp)*c.mpf(multiplier)
            row=self._summary((nominal,bounds))
            row.update(kind='original_macro_radius_log',definition='multiplier*(Y*q+2hb*(1-q))',
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
        out['native']=dict(chart='actual_modified_frozen_macro',Z=Z,fraction_recipe=packet['coordinate'],
            rho_definition='4exp(Y*q+2hb*(1-q))',R_definition='Ra*exp(2hb+(Y-2hb)*q)',
            R100_endpoint_exact=_bounds(t.ctx,q)==(1,1),
            coefficient_geometry='rho4 weights plus full original log radius; n1 mean uses actual current rho')
        out['original_selected_integer']=packet['selected_integer']
        out['actual_modified_leading_errors_included']=True
        out['n1_own_history']=dict(actual_terminal_primitive=self._terminal_cache[Z]['normalized_own_primitive'],
            mean_u=self._summary(n1['mean_u']),q=self._summary(n1['q']),
            primitive_retained_through_current_radius=True,f_and_u_exact_zero=True,
            current_rho_not_rho4_used=True,no_later_reserved_repair_support_on_macro=True)
        out['physical_gradient_or_vorticity_available']=False
        return out

    def evaluate(self,Z,q,*,log_tau='-1',scaled_time=None,azimuth='0'):
        numeric=self._numeric(Z,q,False)
        return dict(**numeric[0],physical=self._physical(numeric,Z,log_tau,scaled_time,azimuth))

    def selected(self,*args,**kwargs):
        raise ValueError('Use OriginalFirstMicroNumeric for selected-sc coordinates')
