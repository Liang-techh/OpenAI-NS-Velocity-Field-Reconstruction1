"""Actual reference/restoration to Rm, with separated j/epsilon histories.

All incoming long histories survive. Unit offsets stay separate from the
huge log reference radius. No reference equilibrium replaces an inlet.
"""
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_long_numeric import OriginalLongNumeric


class OriginalReferenceNumeric:
    charts=('reference','restoration','postrestore')
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.previous=OriginalLongNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        p=self.previous;self.owner=p.owner;self.leading=p.leading;self.transport=p.transport
        self.ctx,self.iv,self.J=p.ctx,p.iv,p.J;self.family=p.family
        self.hashes=dict(p.hashes);self.code_hashes=dict(p.code_hashes)
        path=self.leading.core.directory/(PREFIX+'compliant_current_original_long_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted actual long terminal numerical inlet required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed long prerequisite: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self.T,self.Tp=p.T,p.Tp;c,iv=self.ctx,self.iv;a=self.owner.amplitude
        self.logC,self.logCp=iv.mpf(a.logC_iv),c.mpf(a.logC)
        self.logP,self.logPp=iv.mpf(a.lp_iv),c.mpf(a.lp)
        self.length=10*(self.logC+self.logP)-self.T-8
        self.lengthp=10*(self.logCp+self.logPp)-self.Tp-8
        if _bounds(c,self.length)[0]<=0:raise ValueError('Original positive finite reference gap required')
        import lei_ren_part1_paper_compliant_reference_restore_profiles as source
        from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
        self.source=source;self.sigma_jets=sigma_jets
        for name in (PREFIX+'compliant_reference_restore_profiles.py',PREFIX+'compliant_current_original_reference_restore_functions.py',
                     PREFIX+'compliant_current_original_prescribed_R100_Rm_C6_functions.py'):
            path=self.leading.core.directory/name;digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if self.hashes.get(name,digest)!=digest:raise ValueError('Changed reference source: '+name)
            self.hashes[name]=digest
        self._initial={};self._refs={};self._restores={};self._kernels={}
        name='src/openai_ns_reconstruction/original_reference_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _pair(self,row):return self.ctx.mpf(row['approximation']),self.previous._restore(row)
    def _summary(self,p,t):return self.leading._summary(p,t)
    def _rows(self,p,t):return self.leading._jet_summary(p,t)

    def inputs(self,Z):
        if Z not in self._initial:
            packet=self.previous.evaluate(Z,'1');pairs={}
            for name,rows in packet['unmodified_six_current_amplitude_units'].items():
                pairs[name]=(self.J(self.ctx,[self._pair(r)[0] for r in rows]),
                             self.J(self.iv,[self._pair(r)[1] for r in rows]))
            inp=self.previous.inputs(Z)
            self._initial[Z]=dict(packet=packet,six=pairs,Psi=(inp[0]['values']['Psi'],inp[1]['values']['Psi']),
                G=(inp[0]['G'],inp[1]['G']),logshape=(inp[0]['logshape'],inp[1]['logshape']))
        return self._initial[Z]

    def _decay(self,p,t,rate,gp,g,proof,label):
        c,iv,J=self.ctx,self.iv,self.J
        if _bounds(c,g)==(0,0):return p,t
        r=c.mpf(rate);ri=iv.mpf(rate)
        # exp(-g) is still representable for ordinary g such as10000 or
        # 1e18. A normalized incoming coefficient may have log magnitude
        # ~1e18, so a large-looking g alone does not justify a tiny tail.
        if gp<c.mpf('1e30'):return p*c.exp(-r*gp),t*iv.exp(-ri*g)
        errors=[]
        for k,v in enumerate(t.coefficients):
            mag=self.leading._abs(v);log=None if _bounds(c,mag)[1]==0 else iv.ln(mag)-ri*g
            errors.append(self.leading._symmetric(self.owner._cap(log,proof,label+'_Z'+str(k))))
        return J.constant(c,0,p.order),J(iv,errors)

    def _relax(self,initial,target,rate,gp,g,proof,label):
        dp,dt=self._decay(initial[0]-target[0],initial[1]-target[1],rate,gp,g,proof,label)
        return target[0]+dp,target[1]+dt

    def _reference_six(self,Z,gp,g,proof):
        inp=self.inputs(Z);c,iv,J=self.ctx,self.iv,self.J;v=inp['Psi'];six={}
        const=lambda x:(J.constant(c,c.mpf(x),2),J.constant(iv,iv.mpf(x),2))
        for name,target,rate in (('H_over_Phi',const('1.25'),'1.6'),
            ('Kpsi_over_Phi',(v[0].truncate(2)*c.mpf('1.25'),v[1].truncate(2)*iv.mpf('1.25')),'1.6'),
            ('B_over_Phi2',(J.constant(c,c.mpf(5)/6,2),J.constant(iv,iv.mpf(5)/6,2)),'1.2'),
            ('C_over_Phi2',const(5),'.2'),('Mpsi',(v[0].truncate(2),v[1].truncate(2)),1),
            ('Apsi',((v[0]*v[0]).truncate(2),(v[1]*v[1]).truncate(2)),1)):
            six[name]=self._relax(inp['six'][name],target,rate,gp,g,proof,name+'_reference_memory')
        one=c.mpf(1),iv.mpf(1)
        return six,dict(alpha=one,gamma1=one,gamma2=one,Kalpha_over_Phi=six['H_over_Phi'],Dpsi=six['Mpsi']),v

    def _alpha(self,tp,t):
        c,iv=self.ctx,self.iv
        if _bounds(c,t)==(0,0):return c.mpf(1),iv.mpf(1)
        if _bounds(c,t)==(1,1):return c.mpf(0),iv.mpf(0)
        if tp<=0:point=c.mpf(1)
        elif tp>=1:point=c.mpf(0)
        else:
            odds=1/(1-tp)**2-1/tp**2
            point=c.mpf(0) if odds>10000 else c.mpf(1) if odds<-10000 else 1/(1+c.exp(odds))
        truth=1-self.sigma_jets(iv,t)[0];lo,hi=_bounds(c,truth)
        return point,iv.mpf([max(c.mpf(0),lo),min(c.mpf(1),hi)])

    def kernels(self,coordinate):
        if coordinate not in self._kernels:
            c,iv=self.ctx,self.iv;t,ti=c.mpf(coordinate),iv.mpf(coordinate)
            if t==0:self._kernels[coordinate]=dict.fromkeys(('mean','mixed','square'),(c.mpf(0),iv.mpf(0)))
            else:
                truth=self.source.restoration_kernels(iv,ti,cells=128)
                def alpha(s):
                    if s<=0:return c.mpf(1)
                    if s>=1:return c.mpf(0)
                    odds=1/(1-s)**2-1/s**2
                    return c.mpf(0) if odds>10000 else c.mpf(1) if odds<-10000 else 1/(1+c.exp(odds))
                point={name:c.quad(lambda s:c.exp(-c.mpf(rate)*(t-s))*alpha(s)**power,[0,t/2,t])
                       for name,rate,power in (('mean',1,1),('mixed','1.6',1),('square',1,2))}
                self._kernels[coordinate]={n:(point[n],truth[n]) for n in point}
        return self._kernels[coordinate]

    def _restoration_six(self,Z,coordinate,proof):
        c,iv,J=self.ctx,self.iv,self.J;tp,t=c.mpf(coordinate),iv.mpf(coordinate)
        initial,_,v=self._reference_six(Z,self.lengthp,self.length,proof)
        kernels=self.kernels(coordinate);six={}
        const=lambda x:(J.constant(c,c.mpf(x),2),J.constant(iv,iv.mpf(x),2))
        for name,target,rate in (('H_over_Phi',const('1.25'),'1.6'),
            ('B_over_Phi2',(J.constant(c,c.mpf(5)/6,2),J.constant(iv,iv.mpf(5)/6,2)),'1.2'),
            ('C_over_Phi2',const(5),'.2')):
            six[name]=self._relax(initial[name],target,rate,tp,t,proof,name+'_restore_memory')
        kernel_for=dict(Mpsi=('mean',v),Kpsi_over_Phi=('mixed',(v[0]*2,v[1]*2)),
                        Apsi=('square',(v[0]*v[0],v[1]*v[1])))
        for name,(kind,source) in kernel_for.items():
            rate='1.6' if name=='Kpsi_over_Phi' else 1
            decayp,decayt=self._decay(*initial[name],rate,tp,t,proof,name+'_restore')
            kp,kt=kernels[kind]
            six[name]=decayp+source[0].truncate(2)*kp,decayt+source[1].truncate(2)*kt
        e,ei=c.exp(-tp),iv.exp(-t);alpha=self._alpha(tp,t)
        kp,kt=kernels['mixed'];hp,ht=initial['H_over_Phi']
        dp,dt=initial['Mpsi'];sp,st=kernels['square']
        baseline=dict(alpha=alpha,
            gamma1=(e+kernels['mean'][0],ei+kernels['mean'][1]),
            gamma2=(e+sp,ei+st),
            Kalpha_over_Phi=(hp*c.exp(-c.mpf('1.6')*tp)+2*kp,ht*iv.exp(-iv.mpf('1.6')*t)+2*kt),
            Dpsi=(dp*e+v[0].truncate(2)*sp,dt*ei+v[1].truncate(2)*st))
        return six,baseline,(v[0]*alpha[0],v[1]*alpha[1])

    def _post_six(self,Z,tp,t,proof):
        c,iv,J=self.ctx,self.iv,self.J;six,baseline,_=self._restoration_six(Z,'1',proof)
        const=lambda x:(J.constant(c,c.mpf(x),2),J.constant(iv,iv.mpf(x),2))
        for name,target,rate in (('H_over_Phi',const('1.25'),'1.6'),
            ('B_over_Phi2',(J.constant(c,c.mpf(5)/6,2),J.constant(iv,iv.mpf(5)/6,2)),'1.2'),('C_over_Phi2',const(5),'.2')):
            six[name]=self._relax(six[name],target,rate,tp,t,proof,name+'_post_memory')
        for name,rate in (('Kpsi_over_Phi','1.6'),('Mpsi',1),('Apsi',1)):
            six[name]=self._decay(*six[name],rate,tp,t,proof,name+'_post_memory')
        baseline['gamma1']=(baseline['gamma1'][0]*c.exp(-tp),baseline['gamma1'][1]*iv.exp(-t))
        baseline['gamma2']=(baseline['gamma2'][0]*c.exp(-tp),baseline['gamma2'][1]*iv.exp(-t))
        for name,rate in (('Kalpha_over_Phi','1.6'),('Dpsi',1)):
            baseline[name]=self._decay(*baseline[name],rate,tp,t,proof,name+'_post_memory')
        return six,baseline,(J.constant(c,0,3),J.constant(iv,0,3))

    def _covers(self,Z,chart,q):
        c,iv=self.ctx,self.iv;o=self.owner
        _,prefix,frames=self.previous._covers(Z,iv.mpf(1));local=[]
        zl,zh=_bounds(c,iv.mpf(Z));ql,qh=_bounds(c,q)
        for frame in o.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            for name in self.charts[:self.charts.index(chart)+1]:
                source=next(w for w in frame['actual_17_chart_C3_transports'] if w['chart']==name)
                for cell in source['cells']:
                    ll,lh=_bounds(c,o._endpoint(cell['exact_left']));rl,rh=_bounds(c,o._endpoint(cell['exact_right']))
                    if name==chart and not rh<ql and not ll>qh:local.append(cell)
                    end=min(qh,rh) if name==chart else rh
                    if end>ll:
                        width=iv.mpf(_bounds(c,iv.mpf(end)-iv.mpf(ll))[1])
                        prefix.append((cell,width,iv.ln(self.length) if name=='reference' else iv.mpf(0),name))
        if not local:raise ValueError('Actual reference/restoration source cover required')
        return local,prefix,frames

    def reference(self,Z,q):return self._evaluate(Z,'reference',q)
    def at_reference_gap(self,Z,gap):return self._evaluate(Z,'reference',gap,True)
    def restoration(self,Z,t):return self._evaluate(Z,'restoration',t)
    def postrestore(self,Z,p):return self._evaluate(Z,'postrestore',p)

    def _evaluate(self,Z,chart,coordinate,physical_gap=False):
        c,iv,J=self.ctx,self.iv,self.J;lead,o=self.leading,self.owner
        if type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):raise ValueError('Finite decimal native coordinate required')
        tp,t=c.mpf(coordinate),iv.mpf(coordinate)
        if physical_gap:gp,g=tp,t;tp,t=gp/self.lengthp,g/self.length
        else:gp,g=self.lengthp*tp,self.length*t
        lo,hi=_bounds(c,t)
        if lo<0 or hi>1:raise ValueError('Original native coordinate in[0,1] required')
        inp=self.inputs(Z);proof=[]
        if chart=='reference':six,baseline,v=self._reference_six(Z,gp,g,proof)
        elif chart=='restoration':six,baseline,v=self._restoration_six(Z,coordinate,proof)
        else:six,baseline,v=self._post_six(Z,tp,t,proof)
        # Atomic geometry keeps -8,-7,-6 and unit restoration distances
        # outside the enclosure of the enormous log reference radius.
        if chart=='reference':
            terms=[(self.Tp*(1-tp),self.T*(1-t)),(10*self.logCp*tp,10*self.logC*t),
                   (10*self.logPp*tp,10*self.logP*t),(-8*tp,-8*t)]
            logup=-inp['logshape'][0]+(self.Tp/10-self.logCp)*(1-tp)+(self.logPp-c.mpf('.8'))*tp
            logut=-inp['logshape'][1]+(self.T/10-self.logC)*(1-t)+(self.logP-iv.mpf('.8'))*t
        else:
            op,offset=(-8+tp,-8+t) if chart=='restoration' else (-7+tp,-7+t)
            terms=[(10*self.logCp,10*self.logC),(10*self.logPp,10*self.logP),(op,offset)]
            logup=-inp['logshape'][0]+self.logPp+op/10
            logut=-inp['logshape'][1]+self.logP+offset/10
        yp=sum((a for a,b in terms),c.mpf(0));y=sum((b for a,b in terms),iv.mpf(0))
        if chart=='reference' and not physical_gap and tp==0:
            yp,y=self.Tp,self.T
        if chart=='reference' and physical_gap:
            yp,y=self.Tp+gp,self.T+g
            terms=[(self.Tp,self.T),(gp,g)]
            logup=-inp['logshape'][0]+self.Tp/10-self.logCp+gp/10
            logut=-inp['logshape'][1]+self.T/10-self.logC+g/10
        # Once reference starts, the G/inlet terms in physical swirl have
        # already canceled. Its remaining Z-dependent angular factor is fixed.
        lp=inp['G'][0]*c.mpf(o.amplitude.Lambda)-inp['logshape'][0]-c.ln(220)/2
        lt=inp['G'][1]*iv.mpf(o.amplitude.Lambda_iv)-inp['logshape'][1]-iv.ln(220)/2
        for a,b in terms:lp-=a*c.mpf('.4');lt-=b*iv.mpf('.4')
        local,prefix,frames=self._covers(Z,chart,t);A,B,hist=o._correction_logs(local,prefix,False)
        def correction(logs,order,label):
            return J(iv,[lead._symmetric(o._cap(logs[k],proof,label+'_Z'+str(k))/math.factorial(k)) for k in range(order+1)])
        aj=correction(A,3,'reference_A_over_N');bj=correction(B,3,'reference_delta_Psi')
        meanlogs=[None if x is None else x+lead.core.lp_iv+lead.core.ll_iv for x in hist['m']]
        mj=correction(meanlogs,2,'reference_delta_Mpsi')
        point=dict(logPhi=lp,Psi=v[0],Mpsi=six['Mpsi'][0]);truth=dict(logPhi=lt+aj,Psi=v[1]+bj,Mpsi=six['Mpsi'][1]+mj)
        def recover(values,ctx,de):
            z=J.variable(ctx,ctx.mpf(Z),2);m=values['Mpsi'];psi=values['Psi'].truncate(2)
            return ((2*z*psi-z*m*(1-de)-(1-z*z)*lead._derivative(m))/(1-z*z*de)).truncate(1)
        point['Qpsi']=recover(point,c,lead.core.delta);truth['Qpsi']=recover(truth,iv,lead.core.delta_iv)
        base={name:(self._rows(*pair) if name in ('Kalpha_over_Phi','Dpsi') else self._summary(*pair)) for name,pair in baseline.items()}
        packet=dict(Z=Z,chart=chart,source_family=self.family,
            coordinate=dict(kind='exact_reference_gap' if physical_gap else 'exact_decimal_'+chart,value=coordinate,
                native=self._summary(tp,t),y=self._summary(yp,y),y_log_terms=[self._summary(a,b) for a,b in terms],
                physical_measure='original_reference_length*dq' if chart=='reference' else 'dt once'),
            original_T_singleton_exact_mpf_tuple=self.previous.T_tuple,
            normalized_modified_leading={name:self._rows(point[name],truth[name]) for name in point},
            physical_similarity_log_Utheta=self._rows(logup,logut+aj),
            baseline_axial_histories=base,
            unmodified_six_current_amplitude_units={name:self._rows(*pair) for name,pair in six.items()},
            baseline_contract=dict(V='4Z+j*alpha+epsilon*Psi_sector',M='4Z+j*gamma1+epsilon*Mpsi',
                K='4Z*H+j*Phi*Kalpha_over_Phi+epsilon*Phi*Kpsi_over_Phi',
                A='16Z^2+8Z*j*gamma1+j^2*gamma2+8Z*epsilon*Mpsi+2j*epsilon*Dpsi+epsilon^2*Apsi',
                B='Phi^2*B_over_Phi2',C='Phi^2*C_over_Phi2; undressed',
                auxiliary_inlet_identity='accepted raw K=U0*H+epsilon*Kpsi and A=U0^2+2U0*epsilon*Mpsi+epsilon^2*Apsi; alpha=1 before restore imply Kalpha=H and Dpsi=Mpsi',
                normalized_mixed_ODE='(Kalpha/Phi)_y+1.6*(Kalpha/Phi)=2*alpha; source2, equilibrium1.25',
                j_and_epsilon_sectors_never_merged=True,all_nonzero_histories_retained=True),
            selected_integer=inp['packet']['selected_integer'],
            source_cover=dict(exact_Z_cells=frames,local_source_cells=len(local),prefix_source_cells=len(prefix),
                original_fast_phase='frac(N*(log(rho/4)-hbB*sc/2))',fast_phase_not_chosen=True,
                reference_length_definition='10*(original_logCstar+original_logPstar)-original_T-8',
                physical_current_prefix_measure='reference L*dq; restoration dt; postrestore dp; each once',
                all_incoming_five_corrections_retained=True),
            original_five_correction_history_bounds={name:[o._source_summary(x) for x in rows] for name,rows in hist.items()},
            positive_order_n1_own_history=inp['packet']['positive_order_n1_own_history'],
            small_error_proofs=proof,
            limits=dict(logPhi_Psi_Z_order=3,own_mean_Z_order=2,Qpsi_Z_order=1,unmodified_six_Z_order=2,
                full_physical_vector_assembled=False,merged_modified_six=False,pressure_exported=False,
                mean_Z3_and_QZZ=False,resolved_signed_correction_integrals=False,resolved_fast_phase=False,
                total_radial_derivatives=False,actual_patch_connected=False,global_N_admission=False,
                global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        if chart!='reference':packet['restoration_kernels']={n:self._summary(*p) for n,p in self.kernels(coordinate if chart=='restoration' else '1').items()}
        return packet
