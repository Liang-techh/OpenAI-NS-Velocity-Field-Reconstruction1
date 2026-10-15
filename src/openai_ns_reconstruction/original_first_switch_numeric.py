"""Original selected-N first switch from the actual R100 numerical inlet.

Angular control is unweighted; the axial source uses the true flat
1-sigma prefix mass. Both retain two original width factors.
"""
import gzip
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_collar_velocity import _InletGeometry
from .original_macro_numeric import OriginalMacroNumeric
from .original_second_micro_numeric import OriginalSecondMicroNumeric


class OriginalFirstSwitchNumeric(OriginalMacroNumeric):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        path=self.leading.core.directory/(PREFIX+'compliant_current_original_macro_numeric.json')
        report=json.loads(path.read_text(encoding='utf-8'))
        if not report['checks']['passed'] or report['source_family']!=self.family:
            raise ValueError('Accepted actual selected-N R100 numerical inlet required')
        for name,digest in report['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed accepted macro numerical source: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        stem=PREFIX+'compliant_current_original_prescribed_R100_Rm_C6_functions'
        path=self.leading.core.directory/(stem+'_check.json')
        checked=json.loads(path.read_text(encoding='utf-8'))
        if (not checked['all_passed'] or checked['source_family']!=self.family
                or not checked['current_original_true_prescribed_R100_R110_C6_functions_installed']):
            raise ValueError('Accepted actual prescribed R100 switch source required')
        for name,digest in checked['input_hashes'].items():
            if hashlib.sha256((self.leading.core.directory/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed prescribed switch source: '+name)
            if self.hashes.get(name,digest)!=digest:raise ValueError('Conflicting prescribed switch source')
            self.hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        path=self.leading.core.directory/(stem+'.json.gz')
        source=json.loads(gzip.decompress(path.read_bytes()))
        self.width_identity=source['prescribed_R100_Rm_C6_bindings']['geometry_certificate']
        for gate in ('shared_upstream_is_same_object','bridge_logh_is_actual_upstream_object',
                     'switch_width_is_exact_source_reboxing','bridge_hbB_switch_hbS_kept_distinct'):
            if not self.width_identity[gate]:raise ValueError('Original bridge/switch width identity required')
        ledger=PREFIX+'compliant_K1_ledger.json'
        if self.width_identity['shared_K1_ledger_sha256']!=self.hashes[ledger]:
            raise ValueError('Identical actual K1 width ledger required')
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        module=__import__(PREFIX+'compliant_first_switch_leading',fromlist=['switch_control_source_bridge'])
        self.control_binding=module.switch_control_source_bridge()
        if (self.control_binding['first_angular_weight']!='1'
                or self.control_binding['first_axial_weight']!='1-sigma(s)'):
            raise ValueError('Original unequal angular/axial switch controls required')
        for name,digest in self.control_binding['input_hashes'].items():
            if self.hashes.get(name,digest)!=digest:raise ValueError('Changed switch control source')
            self.hashes[name]=digest
        self._first_switch_background={};self._switch_increment_cache={}
        relative='src/openai_ns_reconstruction/original_first_switch_numeric.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()

    def _coordinate(self,coordinate,selected):
        c,iv=self.transport.ctx,self.transport.iv
        if selected or type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Finite decimal first-switch phase in[0,1] required')
        s=iv.mpf(coordinate);lo,hi=_bounds(c,s)
        if lo<0 or hi>1:raise ValueError('First-switch phase in[0,1] required')
        return s,False,dict(kind='exact_decimal_first_switch_phase',value=coordinate)

    def _covers(self,Z,s):
        c,iv=self.transport.ctx,self.transport.iv
        zl,zh=_bounds(c,iv.mpf(Z));sl,sh=_bounds(c,s)
        local=[];prefix=[];frames=[]
        for frame in self.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            frames.append(frame['exact_Z_cell'])
            for chart in ('first_micro','second_micro','frozen_macro','first_switch'):
                current=next(w for w in frame['actual_17_chart_C3_transports'] if w['chart']==chart)
                for cell in current['cells']:
                    left,right=self._endpoint(cell['exact_left']),self._endpoint(cell['exact_right'])
                    ll,lh=_bounds(c,left);rl,rh=_bounds(c,right)
                    if chart=='first_switch' and not rh<sl and not ll>sh:local.append(cell)
                    end=min(sh,rh) if chart=='first_switch' else rh
                    if end>ll:
                        width=iv.mpf(_bounds(c,iv.mpf(end)-iv.mpf(ll))[1])
                        logmeasure=iv.ln(self.Y) if chart=='frozen_macro' else self.leading.logh
                        prefix.append((cell,width,logmeasure,chart))
        if not local:raise ValueError('Actual complete first-switch source cover required')
        return local,prefix,frames

    def _switch_increments(self,Z,s):
        """Uniform prefix source errors, with the original angular/axial weights."""
        key=Z,s._mpi_
        if key in self._switch_increment_cache:return self._switch_increment_cache[key]
        o,t=self.leading,self.transport;iv=o.iv;b=o.bridge;c=b.ctx
        p=o._prepared[Z]
        module=__import__(PREFIX+'compliant_actual_bridge_integrals',fromlist=['jet_norm'])
        upper=_bounds(t.ctx,s)[1]
        with mp.workdps(280):
            phase=c.mpf(upper);cutmass,_=b.pulse_masses(phase)
            # Convert R0 modes ONCE to R100; comparison moments continue
            # at the current radius instead of freezing their R100 values.
            def modal_norm(name,power):
                return sum((module.jet_norm(row,b.weight)
                    *module.upper(c,abs(p['R0']**j*c.mpf(100)**(power-j)))
                    *c.exp(abs(power-j)*b.cap*phase)
                    for j,row in enumerate(p['modes'][name])),c.mpf(0))
            angular=o._transfer(modal_norm('D_over_R',1)*phase/2)
            ell=o._width_error(o._width_error(angular))
            macro_ell=o._width_error(o._transfer(p['angular_prefix_norm_per_hb']))
            ratio=o._transfer(module.jet_norm(p['quotient'],b.weight))*iv.exp(macro_ell+ell)
            psi=iv.mpf(0)
            for name,power,scale in (('hydro',1,c.mpf(0)),('pressure',1,b.pressure_log),('swirl',2,b.swirl_log)):
                norm=o._transfer(modal_norm('drive_'+name,power)*cutmass)*ratio
                term=b.source_term(p['data']['phi0']*0,2,scale)
                # _scaled removes one formal hb; _width_error applies
                # the other. Pressure, F0 and epsilon are normalized once.
                psi+=o._abs(o._scaled(norm,term,divide_eps=True))
            psi=o._width_error(psi)
        evidence=dict(original_hb_power_angular=2,original_hb_power_axial=2,
            angular_weight='1',axial_weight='1-sigma(s)',
            axial_prefix_mass_exact_mpf_tuples=[list(v) for v in cutmass._mpi_],
            angular_prefix_mass_exact_mpf_tuples=[list(v) for v in s._mpi_],
            original_source_width_identity=self.width_identity,
            current_comparison_modes_not_R100_constant=True,
            uniform_every_prefix_by_positive_modal_norms=True,
            all_hydro_pressure_swirl_parts=True,physical_hb_measure_applied_once=True)
        self._switch_increment_cache[key]=ell,psi,evidence
        return self._switch_increment_cache[key]

    def _switch_geometry(self,s):
        o,t=self.leading,self.transport;c,iv=t.ctx,t.iv
        tiny=o._width_error(s)
        rho=iv.exp(iv.ln(4)+self.Y+iv.mpf([0,_bounds(c,tiny)[1]]))
        return rho,self.Y,tiny

    def _background(self,Z,s):
        key=Z,s._mpi_
        if key in self._first_switch_background:return self._first_switch_background[key]
        o,t=self.leading,self.transport;iv=o.iv;order=3
        # Keep _geometry inherited as the MACRO geometry, because this
        # explicit parent evaluation must reach its true q=1 endpoint.
        point,inlet,_=OriginalMacroNumeric._background(self,Z,iv.mpf(1))
        ell,psi,_=self._switch_increments(Z,s)
        phi_error=o._norm(inlet['Phi'])*ell*iv.exp(ell)
        truth=dict(Phi=inlet['Phi']+o._norm_error_jet(phi_error,order),
            Psi=inlet['Psi']+o._norm_error_jet(psi,order))
        # Exact source: M(s)=exp(-hb*s)M100 + hb*int exp(-hb*(s-t))*Psi(t)dt.
        # Subtract the shared incoming mean before bounding; retain its
        # nonzero value. The positive kernel <=1 bounds this small change.
        mean_error=o._width_error(s*(o._norm(inlet['Psi']-inlet['Mpsi'])+psi))
        truth['Mpsi']=inlet['Mpsi']+o._norm_error_jet(mean_error,order)
        rho,_,_=self._switch_geometry(s)
        self._first_switch_background[key]=point,truth,rho
        return point,truth,rho

    def _numeric(self,Z,coordinate,selected):
        numeric=OriginalSecondMicroNumeric._numeric(self,Z,coordinate,selected)
        packet=numeric[0];s=numeric[4]
        packet['source_cover'].update(native_source_domain=['0','1'],local_source_chart='first_switch',
            cumulative_history_charts=['first_micro','second_micro','frozen_macro','first_switch'],
            true_phase='fractional_part(original_N*(Y+hbS*s-hbB*original_sc/2))')
        packet['background_control_contract']=self._switch_increments(Z,s)[2]
        packet['history_contract'].update(measure='micro hbB ds, macro (Y-2hbB)dq, first-switch hbS ds; each once',
            partial_prefix='Both micro windows, complete macro and current first-switch prefix; positive own-rate kernels bounded by1',
            actual_R100_incoming_all_five_retained=True)
        packet['own_mean_contract']=dict(equation='dMpsi/dy=Psi-Mpsi',
            exact_source='exp(-hbS*s)*Mpsi100 + hbS*integral_0^s exp(-hbS*(s-t))*Psi(t)dt',
            numerical_change_bound='hbS*s*(norm(Psi100-Mpsi100)+uniform norm(Psi-Psi100))',
            actual_nonzero_inlet_retained=True,comparison_mean_not_substituted=True)
        packet['limits'].pop('actual_modified_second_micro_velocity_Z3',None)
        packet['limits'].update(actual_modified_first_switch_velocity_Z3=True,
            first_switch_endpoint_R1_reached=_bounds(self.transport.ctx,s)==(1,1),all_later_charts=False)
        return numeric

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        numeric=self._numeric(Z,s,False)
        return dict(**numeric[0],physical=self._physical(numeric,Z,log_tau,scaled_time,azimuth))

    def selected(self,*args,**kwargs):
        raise ValueError('Use OriginalFirstMicroNumeric for selected-sc coordinates')

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
        rhop=4*t.ctx.exp(self.Ypoint)
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
        _,base,tiny=self._switch_geometry(s)
        def geometry(multiplier):
            value=base*t.iv.mpf(multiplier)
            bound=tiny*t.iv.mpf(multiplier)
            bounds=iv.mpf([mp.make_mpf(v) for v in (value+t.iv.mpf([0,_bounds(t.ctx,bound)[1]]))._mpi_])
            nominal=c.mpf(self.Ypoint)*c.mpf(multiplier)
            row=self._summary((nominal,bounds))
            row.update(kind='original_first_switch_radius_log',definition='multiplier*(Y+hbS*s)',
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
        out['native']=dict(chart='actual_modified_first_switch',Z=Z,phase_recipe=packet['coordinate'],
            rho_definition='100Lambda*exp(hbS*s)',R_definition='100*exp(hbS*s)',
            R100_inlet_exact=_bounds(t.ctx,s)==(0,0),
            R1_endpoint_exact=_bounds(t.ctx,s)==(1,1),
            coefficient_geometry='rho4 physical weights plus full original log radius; n1 directed mean uses current rho, with R100 point proxy and complete width error')
        out['original_selected_integer']=packet['selected_integer']
        out['actual_modified_leading_errors_included']=True
        out['n1_own_history']=dict(actual_terminal_primitive=self._terminal_cache[Z]['normalized_own_primitive'],
            mean_u=self._summary(n1['mean_u']),q=self._summary(n1['q']),
            primitive_retained_through_current_radius=True,f_and_u_exact_zero=True,
            current_rho_not_rho4_used=True,directed_mean_uses_actual_current_rho=True,
            point_mean_radius_proxy='100Lambda; omitted hbS*s factor fully enclosed by actual primitive/current_rho',
            point_proxy_does_not_select_original_hb=True,
            no_later_reserved_repair_support_on_first_switch=True)
        out['physical_gradient_or_vorticity_available']=False
        return out
