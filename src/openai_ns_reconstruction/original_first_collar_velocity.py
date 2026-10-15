"""Physical leading+n1 velocity on the actual admitted first-half collar.

Native phase is required. Original hb remains a defined positive source,
with its log enclosure: no width bound or midpoint becomes the field width.
Inlet, leading hb increments and full n1 increments stay separate.
"""
import copy
import hashlib
import json

from .original_core_log_amplitude import OriginalCoreLogAmplitude
from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_first_collar_transport import OriginalFirstCollarTransport


class _InletGeometry(OriginalCoreVelocity):
    """Reuse accepted physical weights at rho=4, before collar expansion."""
    def __init__(self, owner, profiles, n1):
        self.ctx,self.iv,self.amplitude=owner.ctx,owner.iv,owner.amplitude
        self.family,self.profiles,self.values=owner.family,profiles,n1
    def _leading(self,Z,rho):
        if rho!='4':raise ValueError('Inlet geometry requires rho=4')
        return self.profiles
    def _n1(self,Z,rho):
        if rho!='4':raise ValueError('Inlet geometry requires rho=4')
        return self.values


class OriginalFirstCollarVelocity:
    """Cartesian u,v,w at native (Z,s,azimuth,time), with directed errors.

This is a finite two-order coefficient field on 0<=s<=1/2. It is neither
an ordinary floating-point vector nor an admitted infinite temporal sum.
"""
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.transport=OriginalFirstCollarTransport(precision=precision,directory=directory)
        self.leading=self.transport.leading
        self.amplitude=OriginalCoreLogAmplitude(precision=amplitude_precision,directory=directory)
        self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.family=self.transport.family
        if self.amplitude.family!=self.family:raise ValueError('One original collar/amplitude family required')
        self.hashes=dict(self.transport.hashes)
        for name,digest in self.amplitude.hashes.items():
            if self.hashes.get(name,digest)!=digest:raise ValueError('Changed original amplitude input: '+name)
            self.hashes[name]=digest
        self.code_hashes=dict(self.transport.code_hashes)
        for name in ('original_core_velocity.py','original_core_log_amplitude.py','original_first_collar_velocity.py'):
            relative='src/openai_ns_reconstruction/'+name
            self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        name=PREFIX+'compliant_current_original_first_collar_transport.json'
        path=self.leading.core.directory/name
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted same-family finite-width n1 collar transport required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Accepted collar transport implementation changed: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._cache={}

    def _restore(self,row,key='source_enclosure_exact_mpf_tuples'):
        return self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in row[key]])

    def _pair(self,row,key='source_enclosure_exact_mpf_tuples'):
        return self.ctx.mpf(row['approximation']),self._restore(row,key)

    def _summary(self,pair):return self.amplitude.summarize(*pair)

    def _sources(self,Z,s):
        key=Z,s
        if key not in self._cache:
            leading=self.leading.evaluate(Z,s,degree=24,order=1)
            n1=self.transport.evaluate(Z,s)
            if leading['source_family']!=n1['source_family']:raise ValueError('Mixed collar source family')
            self._cache[key]=leading,n1
        return self._cache[key]

    def _log_width(self):
        return dict(kind='original_log_hb',definition='log(cstar)-100*log(actual_K_9_16)',
            source_definition=PREFIX+'compliant_K1_ledger.json:parameter_definition',
            source_enclosure_exact_mpf_tuples=[list(t) for t in self.leading.logh._mpi_],
            source_value_selected=False,exponential_materialized=False)

    def _geometry_log(self,s,multiplier):
        c,iv=self.ctx,self.iv
        size=self.leading._width_error(self.leading.iv.mpf(s)*multiplier)
        high=_bounds(c,size)[1]
        row=self._summary((c.mpf(0),iv.mpf([0,high])))
        row.update(kind='original_hb_linear_log',phase=s,multiplier=str(multiplier),
            definition='multiplier*hb*phase',original_log_hb=self._log_width(),
            zero_is_numeric_proxy_only=True,original_positive_width_preserved=True)
        return row

    def _rotate(self,row,pair):
        out=copy.deepcopy(row)
        p,t=self._pair(row['coefficient'])
        out['coefficient']=self._summary((p*pair[0],t*pair[1]))
        return out

    def _n1_increments(self,packet,Z):
        c,iv=self.ctx,self.iv
        v={}
        for name in ('f','u','k'):
            row=packet['factored_first_collar_increments'][name]
            if row['epsilon_power'] or row['Pstar_squared_power']:
                raise ValueError('Natural value increments must have scale one')
            v[name]=[self._pair(q) for q in row['coefficient']]
        def recover(i,z,de):
            du,dk=v['u'][0][i],v['k'][0][i]
            mean=du+dk
            meanZ=v['u'][1][i]+v['k'][1][i]
            return mean,(2*z*du-(1+de)*z*mean-(1-z*z)*meanZ)/(1-de*z*z)
        mp,qp=recover(0,c.mpf(Z),self.amplitude.delta)
        mt,qt=recover(1,iv.mpf(Z),self.amplitude.delta_iv)
        return dict(f=v['f'][0],u=v['u'][0],mean_u=(mp,mt),q=(qp,qt))

    def _comparison(self,packet,leading,n1,Z,log_tau,scaled_time):
        c,iv=self.ctx,self.iv
        row=leading['normalized_fields_and_histories']['Phi']
        phi,phii=self._pair(row['inlet'][0])
        delta=self._restore(row['increment_per_original_hb'][0])
        size=self.leading._width_error(self.leading._abs(self.leading._transfer(delta)))
        bound=_bounds(c,size)[1]
        phii+=iv.mpf([-bound,bound])
        f,fi=self._pair(n1['normalized_original_n1']['f'][0])
        lo,hi=_bounds(c,fi)
        if not f or not (lo>0 or hi<0) or _bounds(c,phii)[0]<=0:
            return dict(available=False,reason='Natural swirl sign unresolved')
        adapter=_InletGeometry(self,{}, {})
        _,de,_=adapter._time(Z,log_tau,scaled_time)
        point=self.amplitude.ll+2*sum((p[0] for p in de),c.mpf(0))+c.log(abs(f)/phi)
        truth=self.amplitude.ll_iv+2*sum((p[1] for p in de),iv.mpf(0))+iv.ln(abs(fi)/phii)
        lo,hi=_bounds(c,truth)
        return dict(available=True,log_absolute_ratio=self._summary((point,truth)),
            classification='n1_smaller' if hi<0 else 'n1_larger' if lo>0 else 'unresolved',
            formula='logLambda+2*delta*log(lambda)+log(abs(f_s/Phi_s))',
            F0_and_collar_radius_cancelled_exactly=True,full_finite_width_profiles_enclosed=True,
            proves_higher_temporal_recursion=False)

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        c,iv=self.ctx,self.iv
        if (any(type(v) is not str for v in (Z,s,azimuth))
                or not all(c.isfinite(c.mpf(v)) for v in (Z,s,azimuth))
                or not -1<c.mpf(Z)<1 or not 0<=c.mpf(s)<=c.mpf('.5')):
            raise ValueError('Finite decimal |Z|<1, first phase0..1/2 and azimuth required')
        leading,n1=self._sources(Z,s)
        fields=leading['normalized_fields_and_histories']
        profiles={name:self._pair(fields[source]['inlet'][0]) for name,source in
            (('Phi','Phi'),('Psi','Psi'),('mean_Psi','Mpsi'))}
        profiles['U0']=4*c.mpf(Z)+self.amplitude.j,4*iv.mpf(Z)+self.amplitude.j_iv
        for name,source in (('Qaxis','Qaxis'),('Qpsi','Qpsi_inlet')):
            profiles[name]=self._pair(leading['radial_recovery'][source][0])
        inlet={}
        for name in ('f','u','mean_u','q'):
            row=n1['incoming_core'][name];row=row[0] if isinstance(row,list) else row
            inlet[name]=self._pair(row,'original_source_enclosure_exact_mpf_tuples')
        angle,anglei=c.mpf(azimuth),iv.mpf(azimuth)
        cos,sin=(c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))
        adapter=_InletGeometry(self,profiles,inlet)
        out=adapter._assemble(Z,'4',log_tau,scaled_time,(cos,sin))
        def add(group,source_name,name,pair,width=False):
            old=next(row for row in group if row['name']==source_name)
            row=copy.deepcopy(old);row.update(name=name,coefficient=self._summary(pair))
            if width:row['log_scale_terms'].append(self._log_width())
            row['increment_semantics']='per original hb' if width else 'complete natural finite-width increment'
            group.append(row)
        uz,ur,ut,psi=(out['cylindrical']['uz'],out['cylindrical']['ur'],
                      out['cylindrical']['utheta'],out['streamfunction'])
        for group,old,name,source in (
            (uz,'leading_axial_correction','leading_axial_collar_increment','Psi'),
            (ut,'leading_swirl','leading_swirl_collar_increment','Phi'),
            (psi,'leading_mass_correction','leading_mass_collar_increment','Mpsi')):
            add(group,old,name,self._pair(fields[source]['increment_per_original_hb'][0]),True)
        add(ur,'leading_radial_correction','leading_radial_collar_increment',
            self._pair(leading['radial_recovery']['Qpsi_increment_per_hb'][0]),True)
        delta=self._n1_increments(n1,Z)
        for group,name,source in ((uz,'n1_axial','u'),(ut,'n1_swirl','f'),
                                  (ur,'n1_radial','q'),(psi,'n1_mass','mean_u')):
            add(group,name,name+'_collar_increment',delta[source])
        half=self._geometry_log(s,c.mpf('.5'));whole=self._geometry_log(s,c.mpf(1))
        for group in (ur,ut):
            for row in group:row['log_scale_terms'].append(copy.deepcopy(half))
        for row in psi:row['log_scale_terms'].append(copy.deepcopy(whole))
        out['coordinates']['r']['log_scale_terms'].append(copy.deepcopy(half))
        out['coordinates']['x']=self._rotate(out['coordinates']['r'],cos)
        out['coordinates']['y']=self._rotate(out['coordinates']['r'],sin)
        negsin=-sin[0],-sin[1]
        out['Cartesian']=dict(u=[self._rotate(q,cos) for q in ur]+[self._rotate(q,negsin) for q in ut],
            v=[self._rotate(q,sin) for q in ur]+[self._rotate(q,cos) for q in ut],w=copy.deepcopy(uz))
        out['native']=dict(chart='original_first_collar',Z=Z,phase=s,azimuth=azimuth,
            rho_definition='4*exp(hb*phase)',physical_R_definition='epsilon*rho',
            base_rho='4',original_log_hb=self._log_width())
        out['swirl_order_comparison']=self._comparison(out,leading,n1,Z,log_tau,scaled_time)
        out.pop('original_core_coefficient_and_log_errors_available',None)
        out.update(original_first_collar_coefficient_and_log_errors_available=True,
            original_positive_hb_not_selected_from_enclosure=True,
            leading_hb_increments_and_n1_natural_increments_kept_separate=True,
            original_collar_geometry_kept_in_separate_log_terms=True,
            native_phase_required=True,direct_Cartesian_inverse_available=False,
            actual_native_domain=dict(first_phase_min='0',first_phase_max='.5'),
            normalized_source_packets=dict(leading_axial_order=1,n1_axial_order=1,
                common_original_source_family=self.family,
                n1_finite_width_remainder_budgets=n1['finite_width_remainder_budgets']),
            value_semantics='sum coefficient*exp(sum separate log_scale_terms); original_log_hb and original_hb_linear_log are defined source scalars with directed enclosures, never selected parameter values')
        return out
