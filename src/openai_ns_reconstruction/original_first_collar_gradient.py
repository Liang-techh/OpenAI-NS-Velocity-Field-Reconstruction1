"""Original leading+n1 physical spatial gradient and vorticity on first collar.

Native first derivatives, original coordinate chains and F0 gauge derivative
are used before rotation. Physical exponentials stay factored in shared log
atoms. A core-boundary path supplies the independent native seam endpoint.
"""
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_log_amplitude import OriginalCoreLogAmplitude
from .original_core_velocity import OriginalCoreVelocity
from .original_first_collar_velocity import OriginalFirstCollarVelocity
from .original_first_collar_axial2 import OriginalFirstCollarAxial2


class OriginalFirstCollarGradient(OriginalFirstCollarVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.transport=OriginalFirstCollarAxial2(precision=precision,directory=directory)
        self.leading=self.transport.leading
        self.amplitude=OriginalCoreLogAmplitude(precision=amplitude_precision,directory=directory)
        self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.family=self.transport.family
        if self.amplitude.family!=self.family:
            raise ValueError('Same original gradient/amplitude family required')
        self.hashes=dict(self.transport.hashes);self.code_hashes=dict(self.transport.code_hashes)
        for name,digest in self.amplitude.hashes.items():
            if self.hashes.get(name,digest)!=digest:
                raise ValueError('Conflicting original physical gradient input: '+name)
            self.hashes[name]=digest
        for stem in ('compliant_current_original_first_collar_axial2',
                     'compliant_current_original_first_collar_velocity'):
            path=self.leading.core.directory/(PREFIX+stem+'.json')
            report=json.loads(path.read_text(encoding='utf-8'))
            if not report['checks']['passed'] or report['source_family']!=self.family:
                raise ValueError('Accepted original Z2 transport and physical weights required')
            for name,digest in report['implementation_hashes'].items():
                if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed physical gradient prerequisite: '+name)
                self.code_hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_first_collar_gradient.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._native_cache={}

    def _transfer(self,value,directed):
        if directed:
            return self.iv.mpf([mp.make_mpf(tuple(t)) for t in value._mpi_])
        return self.ctx.mpf(value)

    def _native(self,Z,s,core_boundary=False):
        key=Z,s,core_boundary
        if key in self._native_cache:return self._native_cache[key]
        o=self.leading
        if core_boundary:
            inlet=o.inlet_jets(Z,degree=24,order=5)
            packet=self.transport.n1.evaluate(Z,'4',degree=20,kernel_degree=28)
            native=None;row_key='original_source_enclosure_exact_mpf_tuples'
        else:
            packet=self.transport.evaluate(Z,s)
            native=self.transport._inputs(Z,s,order=5);inlet=native['inlet']
            row_key='source_enclosure_exact_mpf_tuples'
        result={}
        for directed,label in ((False,'point'),(True,'truth')):
            c=self.iv if directed else self.ctx;z=c.mpf(Z)
            de=self.amplitude.delta_iv if directed else self.amplitude.delta
            d,L=1-z*z,1-de*z*z
            a=native['actual'] if directed and native is not None else inlet[label]['values']
            values={name:[self._transfer(a[name][k],directed) for k in range(3)]
                    for name in ('Phi','Psi','Mpsi')}
            comparison=native['comparison'] if directed and native is not None else inlet[label]['radial']
            radial={name:[self._transfer(q[0],directed) for q in comparison[name][:2]]
                    for name in ('Phi','Psi')}
            rho=self._transfer(native['rho'][0],True) if directed and native is not None else c.mpf(4)
            x=c.sqrt(rho)
            if core_boundary:sig=c.mpf(0);h=c.mpf(0)
            elif directed:
                with mp.workdps(280):sig=self.transport.sigma_jets(c,c.mpf(s))[0]
                h=c.mpf([0,_bounds(self.ctx,o._width_error(o.iv.mpf(1)))[1]])
            else:
                t=c.mpf(s);sig=c.mpf(0) if not t else 1/(1+c.exp(1/t**2-1/(1-t)**2))
                h=c.mpf(0) # Numerical center; the directed branch keeps actual hb.
            chi=1-sig+h*sig
            phi,pz,_=values['Phi'];psi,sz,_=values['Psi'];mean,mz,mt=values['Mpsi']
            pr=phi*chi*radial['Phi'][1]/radial['Phi'][0]
            sr=phi/radial['Phi'][0]*chi*radial['Psi'][1]
            mr=(psi-mean)/rho;mrz=(sz-mz)/rho
            U0=4*z+(self.amplitude.j_iv if directed else self.amplitude.j)
            qa=((1+de)*z*U0-4*d)/L
            qaZ=((1+de)*(U0+4*z)+8*z)/L+qa*(2*de*z)/L
            N=2*z*psi-(1-de)*z*mean-d*mz
            NZ=2*psi+2*z*sz-(1-de)*(mean+z*mz)+2*z*mz-d*(2*mt)
            qp=N/L;qpZ=NZ/L+qp*(2*de*z)/L
            qpr=(2*z*sr-(1-de)*z*mr-d*mrz)/L
            def read(row):
                if directed:return self._restore(row,row_key)
                return c.mpf(row['approximation'])
            state={name:[read(q) for q in packet['normalized_original_n1'][name]]
                   for name in ('f','u','k','fx','ux','q')}
            fr=state['fx'][0]/(2*x);ur=state['ux'][0]/(2*x)
            qr=(z*state['ux'][0]/x+(1+de)*z*state['k'][0]/rho+d*state['k'][1]/rho)/L
            fixed,_,_,_=o.core._seed(Z,24,2,directed=directed)
            g=self._transfer(-fixed['ell_Z_taylor'][0],directed)
            profiles=dict(U0=(U0,c.mpf(4),c.mpf(0)),Phi=(phi,pz,pr),Psi=(psi,sz,sr),
                Qaxis=(qa,qaZ,c.mpf(0)),Qpsi=(qp,qpZ,qpr),
                f=(state['f'][0],state['f'][1],fr),u=(state['u'][0],state['u'][1],ur),
                q=(state['q'][0],state['q'][1],qr))
            result[label]=dict(profiles=profiles,rho=rho,g=g,delta=de,d=d,L=L)
        self._native_cache[key]=result
        return result

    @staticmethod
    def _specifications():
        # component, source, temporal order, lambda integer/delta powers,
        # logLambda power, rho power, radius geometry, F0 power.
        return [('ur','Qaxis',0,-1,0,'-.5','.5','radial_geometry',0),
                ('ur','Qpsi',0,-1,0,'-1.5','.5','radial_geometry',0),
                ('ur','q',1,-1,2,'-1.5','.5','radial_geometry',0),
                ('utheta','Phi',0,-1,-1,'-.5','.5','swirl_geometry',1),
                ('utheta','f',1,-1,1,'.5','.5','swirl_geometry',1),
                ('uz','U0',0,-1,-1,'0','0',None,0),
                ('uz','Psi',0,-1,-1,'-1','0',None,0),
                ('uz','u',1,-1,1,'-1','0',None,0)]

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        return self._evaluate(Z,s,log_tau,scaled_time,azimuth,False)

    def core_endpoint(self,Z,*,log_tau='-1',scaled_time=None,azimuth='0'):
        """Original rho=4 core endpoint, using inlet radial rows directly."""
        return self._evaluate(Z,'0',log_tau,scaled_time,azimuth,True)

    def _evaluate(self,Z,s,log_tau,scaled_time,azimuth,core_boundary):
        c,iv=self.ctx,self.iv
        if (any(type(v) is not str for v in (Z,s,azimuth))
                or not all(c.isfinite(c.mpf(v)) for v in (Z,s,azimuth))
                or not -1<c.mpf(Z)<1 or not 0<=c.mpf(s)<=c.mpf('.5')):
            raise ValueError('Finite native |Z|<1, first phase0..1/2 and azimuth required')
        ell,dell,time=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        data=self._native(Z,s,core_boundary)
        point,truth=data['point'],data['truth']
        atoms={}
        for name,pair in zip(('lambda_time','lambda_geometry','delta_lambda_time','delta_lambda_geometry'),ell+dell):
            atoms[name]=self._summary(pair)
        atoms['log_Lambda']=self._summary((self.amplitude.ll,self.amplitude.ll_iv))
        atoms['log_delta']=self._summary((-4*self.amplitude.lp-30,-4*self.amplitude.lp_iv-30))
        atoms['radial_geometry']=self._summary((c.log(2)/2,iv.ln(2)/2))
        atoms['swirl_geometry']=self._summary((c.log(8)/2,iv.ln(8)/2))
        atoms['half_hbs']=self._geometry_log('0' if core_boundary else s,c.mpf('.5'))
        z,zi=c.mpf(Z),iv.mpf(Z)
        # ln(1-a) lies in [-2a,-a] for 0<=a<1/2. This retains tiny delta.
        a=self.amplitude.delta_iv*zi*zi
        if _bounds(c,a)[1]>=c.mpf('.5'):raise ValueError('Original L positivity lost')
        atoms['log_L']=self._summary((c.log1p(-self.amplitude.delta*z*z),-a*iv.mpf([1,2])))
        amplitude=self.amplitude.global_amplitude(Z)
        atoms['F0_anchor']=dict(negative_of_exact_mpf_tuple=self.amplitude.selected_logC_tuple)
        atoms['F0_relative']=amplitude['log_F0_relative_to_anchor']
        pool={};values={n:[] for n in ('ur','utheta','uz')};dr={n:[] for n in values};dz={n:[] for n in values};over_r={n:[] for n in values}
        def logs(a,b,m,geometry,swirl):
            rows={name:c.mpf(value) for name,value in
                  (('lambda_time',a),('lambda_geometry',a),('delta_lambda_time',b),
                   ('delta_lambda_geometry',b),('log_Lambda',m)) if c.mpf(value)}
            if geometry:rows.update({geometry:c.mpf(1),'half_hbs':c.mpf(1)})
            if swirl:rows.update(F0_anchor=c.mpf(1),F0_relative=c.mpf(1))
            return rows
        def adjust(base,changes):
            row=dict(base)
            for atom,factor in changes.items():row[atom]=row.get(atom,c.mpf(0))+factor
            return {atom:factor for atom,factor in row.items() if factor}
        radius_logs=logs(1,0,'-.5','swirl_geometry',0)
        minus_radius={atom:-v for atom,v in radius_logs.items()}
        inverse_z={'lambda_time':-1,'lambda_geometry':-1,'delta_lambda_time':1,'delta_lambda_geometry':1,'log_L':-1}
        def emit(name,order,pair,scales,kind,source):
            if name in pool:raise ArithmeticError('Duplicate physical derivative term')
            pool[name]=dict(temporal_order=order,coefficient=self._summary(pair),
                log_scale_terms=[dict(atom=atom,multiplier=str(value)) for atom,value in scales.items()],
                kind=kind,source=source)
            return name
        def pair(fn):return fn(point,c,z),fn(truth,iv,zi)
        for component,source,order,A,B,M,p,geometry,swirl in self._specifications():
            base=logs(A,B,M,geometry,swirl);power=c.mpf(p)
            val=pair(lambda v,ctx,zz:v['profiles'][source][0])
            values[component].append(emit(component+':'+source,order,val,base,'value',source))
            radial=pair(lambda v,ctx,zz:2*v['rho']*v['profiles'][source][2]+2*ctx.mpf(p)*v['profiles'][source][0])
            dr[component].append(emit(component+'_r:'+source,order,radial,adjust(base,minus_radius),'radial_derivative',source))
            zr=pair(lambda v,ctx,zz:v['d']*v['profiles'][source][1]-2*v['rho']*zz*v['profiles'][source][2]+(A-2*ctx.mpf(p))*zz*v['profiles'][source][0])
            zlogs=adjust(base,inverse_z)
            dz[component].append(emit(component+'_z:'+source,order,zr,zlogs,'axial_regular_derivative',source))
            if B:
                delta_pair=pair(lambda v,ctx,zz:B*zz*v['profiles'][source][0])
                dz[component].append(emit(component+'_z_delta:'+source,order,delta_pair,adjust(zlogs,{'log_delta':1}),'axial_delta_weight_derivative',source))
            if swirl:
                gauge=pair(lambda v,ctx,zz:-v['d']*v['g']*v['profiles'][source][0])
                dz[component].append(emit(component+'_z_F0:'+source,order,gauge,adjust(zlogs,{'log_Lambda':1}),'original_F0_gauge_derivative',source))
            if component!='uz':
                over_r[component].append(emit(component+'_over_r:'+source,order,val,adjust(base,minus_radius),'cylindrical_connection',source))
        def refs(names,factor=1):return [dict(term=name,factor=self._summary((c.mpf(factor),iv.mpf(factor)))) for name in names]
        J=[[refs(dr['ur']),refs(over_r['utheta'],-1),refs(dz['ur'])],
           [refs(dr['utheta']),refs(over_r['ur']),refs(dz['utheta'])],
           [refs(dr['uz']),[],refs(dz['uz'])]]
        cp,sp=c.cos(c.mpf(azimuth)),c.sin(c.mpf(azimuth));ct,st=iv.cos(iv.mpf(azimuth)),iv.sin(iv.mpf(azimuth))
        R=[[(cp,ct),(-sp,-st),(c.mpf(0),iv.mpf(0))],
           [(sp,st),(cp,ct),(c.mpf(0),iv.mpf(0))],
           [(c.mpf(0),iv.mpf(0)),(c.mpf(0),iv.mpf(0)),(c.mpf(1),iv.mpf(1))]]
        def rotate_matrix():
            out=[]
            for i in range(3):
                row=[]
                for j in range(3):
                    entries=[]
                    for a in range(3):
                        for b in range(3):
                            f=R[i][a][0]*R[j][b][0],R[i][a][1]*R[j][b][1]
                            if f[0]==0 and _bounds(c,f[1])==(0,0):continue
                            for ref in J[a][b]:
                                q=self._pair(ref['factor'])
                                entries.append(dict(term=ref['term'],factor=self._summary((q[0]*f[0],q[1]*f[1]))))
                    row.append(entries)
                out.append(row)
            return out
        omega=[refs(dz['utheta'],-1),refs(dz['ur'])+refs(dr['uz'],-1),refs(dr['utheta'])+refs(over_r['utheta'])]
        omega_cart=[]
        for i in range(3):
            row=[]
            for j in range(3):
                f=R[i][j]
                if f[0]==0 and _bounds(c,f[1])==(0,0):continue
                for ref in omega[j]:
                    q=self._pair(ref['factor']);row.append(dict(term=ref['term'],factor=self._summary((q[0]*f[0],q[1]*f[1]))))
            omega_cart.append(row)
        divergence={}
        for name,U,Q,beta in (('leading_axis','U0','Qaxis',-1),('leading_correction','Psi','Qpsi',-1),('n1','u','q',1)):
            def defect(v,ctx,zz):
                u,uz,ur=v['profiles'][U];q,_,qr=v['profiles'][Q]
                return v['L']*(q+v['rho']*qr)+v['d']*uz-2*v['rho']*zz*ur+(-1+beta*v['delta'])*zz*u
            divergence[name]=self._summary(pair(defect))
        def axial_vorticity_coefficient(v,source):
            f,_,fr=v['profiles'][source]
            return 2*v['rho']*fr+2*f
        fp,ft=axial_vorticity_coefficient(point,'f'),axial_vorticity_coefficient(truth,'f')
        pp,pt=axial_vorticity_coefficient(point,'Phi'),axial_vorticity_coefficient(truth,'Phi')
        resolved=lambda q:(_bounds(c,q)[0]>0 or _bounds(c,q)[1]<0)
        comparison=dict(available=False,reason='axial vorticity coefficient sign unresolved')
        if fp and pp and resolved(ft) and resolved(pt):
            value=self.amplitude.ll+2*sum((q[0] for q in dell),c.mpf(0))+c.log(abs(fp/pp))
            box=self.amplitude.ll_iv+2*sum((q[1] for q in dell),iv.mpf(0))+iv.ln(abs(ft/pt))
            lo,hi=_bounds(c,box)
            comparison=dict(available=True,log_absolute_n1_over_leading=self._summary((value,box)),
                classification='n1_smaller' if hi<0 else 'n1_larger' if lo>0 else 'unresolved',
                formula='logLambda+2delta*log(lambda)+log(abs((rho*f_rho+f)/(rho*Phi_rho+Phi)))',
                original_F0_and_radius_cancelled_exactly=True,proves_temporal_n_ge_2_recursion=False)
        native_rows={name:dict(value=self._summary((point['profiles'][name][0],truth['profiles'][name][0])),
            Z=self._summary((point['profiles'][name][1],truth['profiles'][name][1])),
            rho=self._summary((point['profiles'][name][2],truth['profiles'][name][2]))) for name in point['profiles']}
        return dict(native=dict(chart='original_core_rho4_boundary' if core_boundary else 'original_first_collar',
                Z=Z,phase=s,azimuth=azimuth,first_phase_max='.5'),time=time,source_family=self.family,
            log_atoms=atoms,term_pool=pool,
            cylindrical_velocity={name:refs(rows) for name,rows in values.items()},
            gradient_cylindrical=J,gradient_Cartesian=rotate_matrix(),
            vorticity_cylindrical=dict(zip(('r','theta','z'),omega)),
            vorticity_Cartesian=dict(zip(('x','y','z'),omega_cart)),
            normalized_meridional_divergence_enclosures=divergence,native_source_rows=native_rows,
            axial_vorticity_order_comparison=comparison,
            physical_chain=dict(partial_r='(2rho/r)*partial_rho',
                partial_z='(d*partial_Z-2rho*Z*partial_rho)/(lambda^(1-delta)*L)',
                F0_log_derivative='-g/epsilon; its Lambda factor is kept in a separate log atom'),
            semantics='Every reference contributes factor*term_pool.coefficient*exp(sum multiplier*log_atom). Retain all source intervals and exact-source atoms; do not substitute missing approximations with zero.',
            original_F0_gauge_and_delta_derivatives_kept_separate=True,
            native_radial_rows_used_instead_of_differentiating_hb_increment_logs=True,
            limits=dict(native_phase_required=True,first_phase_max='.5',nonzero_radius=True,
                physical_first_spatial_derivatives=True,physical_second_spatial_derivatives=False,
                ordinary_float_velocity_gradient=False,useful_combined_physical_relative_error=False,
                complete_Cartesian_C2_join=False,global_field=False,
                temporal_n_ge_2_recursion=False,full_corrected_NS=False))
