"""Selected-original-N second-micro values with cumulative source histories.

The first coupled window is retained completely. Both velocities admit
ordinary Z3 here, while the inherited mean history only admits useful Z2.
"""
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_micro_numeric import OriginalFirstMicroNumeric


class OriginalSecondMicroNumeric(OriginalFirstMicroNumeric):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        path=self.leading.core.directory/(PREFIX+'compliant_current_original_first_micro_numeric.json')
        report=json.loads(path.read_text(encoding='utf-8'))
        if not report['checks']['passed'] or report['source_family']!=self.family:
            raise ValueError('Accepted original selected-N first-micro value enclosure required')
        for name,digest in report['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed first-micro numerical prerequisite: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_second_micro_numeric.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._second_background={};self._upgraded=set()

    def _coordinate(self,coordinate,selected):
        c,iv=self.transport.ctx,self.transport.iv
        if selected or type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Second-micro decimal native phase1..2 required')
        s=iv.mpf(coordinate);lo,hi=_bounds(c,s)
        if lo<1 or hi>2:raise ValueError('Second-micro native phase in[1,2] required')
        return s,False,dict(kind='exact_decimal_phase',value=coordinate)

    def _covers(self,Z,s):
        c,iv=self.transport.ctx,self.transport.iv
        zl,zh=_bounds(c,iv.mpf(Z));sl,sh=_bounds(c,s)
        local=[];prefix=[];frames=[]
        for frame in self.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            frames.append(frame['exact_Z_cell'])
            for chart in ('first_micro','second_micro'):
                current=next(w for w in frame['actual_17_chart_C3_transports'] if w['chart']==chart)
                for cell in current['cells']:
                    left,right=self._endpoint(cell['exact_left']),self._endpoint(cell['exact_right'])
                    ll,lh=_bounds(c,left);rl,rh=_bounds(c,right)
                    if chart=='second_micro' and not rh<sl and not ll>sh:local.append(cell)
                    if sh>ll:
                        width=iv.mpf(max(c.mpf(0),_bounds(c,iv.mpf(min(sh,rh))-iv.mpf(ll))[1]))
                        if width:prefix.append((cell,width))
        if not local:raise ValueError('Actual second-micro complete source cover required')
        return local,prefix,frames

    def _second_increments(self,prepared,s,order):
        o=self.leading;iv,J=o.iv,o.Jet
        with mp.workdps(280):angular,axial=o.bridge.contributions(prepared,'second',s)
        ell=J.constant(iv,0,order);psi=J.constant(iv,0,order)
        Bell=iv.mpf(0);Bpsi=iv.mpf(0)
        for term in angular:
            row=J(iv,[o._scaled(v,term) for v in term['normalized_axial_coefficients'][:order+1]])
            ell+=row;Bell+=o._norm(row)
        for terms in axial.values():
            for term in terms:
                row=J(iv,[o._scaled(v,term,divide_eps=True) for v in term['normalized_axial_coefficients'][:order+1]])
                error=o._scaled(term['nonlinear_prefix_error_weighted_norm_per_next_hb'],term,
                    extra_power=1,divide_eps=True,symmetric=True)
                psi+=row+o._norm_error_jet(error,order);Bpsi+=o._norm(row)+o._abs(error)
        return ell,Bell,psi,Bpsi

    def _background(self,Z,s):
        key=Z,s._mpi_
        if key in self._second_background:return self._second_background[key]
        o,t=self.leading,self.transport;iv=o.iv;order=3
        inlet=o.inlet_jets(Z,degree=24,order=order)
        if Z not in o._prepared:
            with mp.workdps(280):o._prepared[Z]=o.bridge.prepare(Z)
        upper=_bounds(t.ctx,s)[1]
        _,Bell,_,Bpsi=self._second_increments(o._prepared[Z],upper,order)
        a=inlet['truth']['values']
        Bphi=o._norm(a['Phi'])*Bell*iv.exp(o._width_error(Bell))
        # The mean is transported cumulatively from the actual core. This
        # retains the first-window memory instead of inventing a new inlet.
        Bmean=iv.mpf(upper)*o._norm(a['Psi']-a['Mpsi'])+o._width_error(iv.mpf(upper)*Bpsi)
        truth={name:a[name]+o._norm_error_jet(o._width_error(bound),order)
               for name,bound in (('Phi',Bphi),('Psi',Bpsi),('Mpsi',Bmean))}
        rho=iv.mpf(4)+iv.mpf([0,_bounds(t.ctx,o._width_error(4*iv.mpf(upper)*iv.exp(o.hcap*iv.mpf(upper))))[1]])
        self._second_background[key]=inlet['point']['values'],truth,rho
        return self._second_background[key]

    @staticmethod
    def _third_exponential_increment_bound(x):
        return 3*(x[3]+x[1]*x[2]+x[1]**3/6)

    def _numeric(self,Z,coordinate,selected):
        numeric=super()._numeric(Z,coordinate,selected)
        key=Z,coordinate,selected
        if key in self._upgraded:return numeric
        packet,point,truth,rho,s=numeric
        c,iv,o=self.transport.ctx,self.transport.iv,self.leading
        local,prefix,frames=self._covers(Z,s);A,B,hist=self._correction_logs(local,prefix,False)
        bp,bt,_=self._background(Z,s);proof=packet['small_error_proofs']
        x=[self._cap(A[k],proof,'second_A_over_N_Z'+str(k))/math.factorial(k) for k in range(4)]
        exp_error=self._exponential_increment_bounds(x)+[self._third_exponential_increment_bound(x)]
        delta_phi=bt['Phi']*o.Jet(iv,[o._symmetric(v) for v in exp_error])
        delta_psi=o.Jet(iv,[o._symmetric(self._cap(B[k],proof,'second_delta_Psi_Z'+str(k))/math.factorial(k)) for k in range(4)])
        for name,change in (('Phi',delta_phi),('Psi',delta_psi)):
            point[name]=bp[name];truth[name]=bt[name]+change
            packet['normalized_modified_leading'][name]=o._jet_summary(point[name],truth[name])
            packet['modified_error_contributions'][name]=o._jet_summary(o.Jet.constant(c,0,3),change)
        packet['source_cover'].update(native_source_domain=['1','2'],local_source_chart='second_micro',
            cumulative_history_charts=['first_micro','second_micro'],
            whole_original_fast_phase_cover=packet['source_cover'].pop('whole_original_phase_cover'))
        packet['background_control_contract']=dict(
            first_complete_window='hb*wc(1)+hb^2*ws(1), wc(1)=ws(1)=1/2',
            second_additional_window='hb^2*(s-1) with second comparison direction and quotient',
            both_phi_and_U_controls_retained=True,all_three_axial_drive_parts_retained=True,
            first_values_and_all_own_histories_not_reset=True,physical_measure='original_hb*ds once')
        packet['history_contract']['partial_prefix']='first_micro sc/2..1 plus second_micro1..s, each original hb measure, positive own kernels retained in the source and bounded by1'
        packet['normalized_third_derivative_log_error_bounds']['delta_Psi_numeric_small_error_proved']=True
        packet['normalized_third_derivative_log_error_bounds']['delta_Mpsi_incoming_first_chart_still_unresolved']=True
        packet['coefficient_convention']='Phi/Psi Taylor Z0..3; Mpsi Z0..2; Qpsi Z0..1. k=D_Z^k/k!, no zero padding.'
        packet['limits'].pop('actual_modified_first_micro_values_and_Z2')
        packet['limits'].update(actual_modified_second_micro_velocity_Z3=True,actual_own_mean_Z2=True,
            actual_velocity_Z3_numeric=True,actual_mean_Z3_numeric=False,Qpsi_ZZ_numeric=False)
        self._upgraded.add(key)
        return numeric

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        numeric=self._numeric(Z,s,False)
        physical=self._physical(numeric,Z,log_tau,scaled_time,azimuth)
        physical['native']['chart']='actual_modified_second_micro'
        return dict(**numeric[0],physical=physical)

    def selected(self,*args,**kwargs):
        raise ValueError('Use OriginalFirstMicroNumeric for selected-sc first-chart coordinates')
