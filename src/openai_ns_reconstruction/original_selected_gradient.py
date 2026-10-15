"""Actual selected inner/lower leading+n1 physical velocity and first gradient.

Only F1 and Uz1 are cut. Q is recovered from the retained own primitive.
The inverse hb*(sc/4) derivative scale remains an independent log sector.
This module stops at r_minus, where later leading modifications may start.
"""
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_log_amplitude import OriginalCoreLogAmplitude
from .original_core_velocity import OriginalCoreVelocity
from .original_first_collar_gradient import OriginalFirstCollarGradient
from .original_n1_cutoff_numeric import OriginalN1CutoffNumeric, stable_sigma_point


class OriginalSelectedGradient(OriginalFirstCollarGradient):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,control_cells=128):
        self.transport=OriginalN1CutoffNumeric(precision=precision,directory=directory,control_cells=control_cells)
        self.leading=self.transport.leading
        self.amplitude=OriginalCoreLogAmplitude(precision=amplitude_precision,directory=directory)
        self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.family=self.transport.family
        if self.amplitude.family!=self.family:raise ValueError('Same original selected/amplitude family required')
        self.hashes=dict(self.transport.hashes);self.code_hashes=dict(self.transport.code_hashes)
        for name,digest in self.amplitude.hashes.items():
            if self.hashes.get(name,digest)!=digest:raise ValueError('Conflicting selected physical source: '+name)
            self.hashes[name]=digest
        for stem in ('compliant_current_original_n1_cutoff_numeric',
                     'compliant_current_original_first_collar_gradient'):
            path=self.leading.core.directory/(PREFIX+stem+'.json')
            report=json.loads(path.read_text(encoding='utf-8'))
            if not report['checks']['passed'] or report['source_family']!=self.family:
                raise ValueError('Accepted same-family cutoff and physical chain reports required')
            for name,digest in report['implementation_hashes'].items():
                if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed selected physical prerequisite: '+name)
                self.code_hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_selected_gradient.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._native_cache={}

    def _selected_geometry(self,fraction):
        c,iv=self.ctx,self.iv
        phase=self.transport.phase_iv*self.transport.iv.mpf(fraction)
        size=self.leading._width_error(phase/2)
        high=_bounds(c,size)[1]
        row=self._summary((c.mpf(0),iv.mpf([0,high])))
        row.update(kind='original_hb_linear_log',fraction=fraction,
            exact_phase_in_mpf_tuple=self.transport.phase_tuple,
            definition='original_hb*original_phase_in*fraction/2',
            original_log_hb=self._log_width(),zero_is_numeric_proxy_only=True,
            original_positive_width_preserved=True)
        return row

    def _phase_log(self):
        return self._summary((self.ctx.log(self.transport.phase_exact),
            self.iv.ln(self._transfer(self.transport.phase_iv,True))))

    def _native(self,Z,fraction,lower=False):
        key=Z,fraction,lower
        if key in self._native_cache:return self._native_cache[key]
        o=self.leading;c=self.ctx;iv=self.iv;t=self.transport
        packet=t.lower_collar(Z,fraction) if lower else t.evaluate(Z,fraction)
        phase=t.phase_iv*t.iv.mpf(fraction);bound_phase=_bounds(t.ctx,phase)[1]
        native=t._inputs(Z,bound_phase,order=5);inlet=native['inlet']
        result={'cutoff_packet':packet}
        for directed,label in ((False,'point'),(True,'truth')):
            ctx=iv if directed else c;z=ctx.mpf(Z)
            de=self.amplitude.delta_iv if directed else self.amplitude.delta
            d,L=1-z*z,1-de*z*z
            a=native['actual'] if directed else inlet[label]['values']
            values={name:[self._transfer(a[name][k],directed) for k in range(3)]
                    for name in ('Phi','Psi','Mpsi')}
            comparison=native['comparison'] if directed else inlet[label]['radial']
            radial={name:[self._transfer(q[0],directed) for q in comparison[name][:2]]
                    for name in ('Phi','Psi')}
            rho=self._transfer(native['rho'][0],True) if directed else ctx.mpf(4)
            if directed:
                with mp.workdps(280):sig=t.sigma_jets(ctx,ctx.mpf([0,bound_phase]))[0]
                h=ctx.mpf([0,_bounds(c,o._width_error(o.iv.mpf(1)))[1]])
            else:
                sig=stable_sigma_point(ctx,ctx.mpf(bound_phase));h=ctx.mpf(0)
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
                return self._restore(row) if directed else ctx.mpf(row['approximation'])
            state={name:[read(q) for q in packet['normalized_original_n1'][name]]
                   for name in ('f','u','q')}
            regular={};singular={}
            for name in ('f','u','q'):
                rows=packet['native_first_radial_width_sectors'][name+'_rho']
                regular[name]=read(rows['0,0'][0]);singular[name]=read(rows['-1,-1'][0])
            fixed,_,_,_=o.core._seed(Z,24,2,directed=directed)
            g=self._transfer(-fixed['ell_Z_taylor'][0],directed)
            profiles=dict(U0=(U0,ctx.mpf(4),ctx.mpf(0)),Phi=(phi,pz,pr),Psi=(psi,sz,sr),
                Qaxis=(qa,qaZ,ctx.mpf(0)),Qpsi=(qp,qpZ,qpr),
                **{name:(state[name][0],state[name][1],regular[name]) for name in ('f','u','q')})
            result[label]=dict(profiles=profiles,radial_singular=singular,
                rho=rho,g=g,delta=de,d=d,L=L)
        self._native_cache[key]=result
        return result

    def evaluate(self,Z,q,*,log_tau='-1',scaled_time=None,azimuth='0'):
        """Actual selected inner fraction q in[0,1]."""
        return self._assemble(Z,q,log_tau,scaled_time,azimuth,False)

    def lower_collar(self,Z,fraction,*,log_tau='-1',scaled_time=None,azimuth='0'):
        """Selected lower fraction1..2, ending at r_minus."""
        return self._assemble(Z,fraction,log_tau,scaled_time,azimuth,True)

    def core_endpoint(self,*args,**kwargs):
        raise ValueError('Use the accepted core API for the independent rho=4 endpoint')

    def _assemble(self,Z,s,log_tau,scaled_time,azimuth,lower):
        c,iv=self.ctx,self.iv
        if (any(type(v) is not str for v in (Z,s,azimuth))
                or not all(c.isfinite(c.mpf(v)) for v in (Z,s,azimuth))
                or not -1<c.mpf(Z)<1 or not (1 if lower else 0)<=c.mpf(s)<=(2 if lower else 1)):
            raise ValueError('Finite native |Z|<1, selected inner fraction0..1 or lower1..2 and azimuth required')
        ell,dell,time=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        data=self._native(Z,s,lower)
        point,truth=data['point'],data['truth']
        atoms={}
        for name,pair in zip(('lambda_time','lambda_geometry','delta_lambda_time','delta_lambda_geometry'),ell+dell):
            atoms[name]=self._summary(pair)
        atoms['log_Lambda']=self._summary((self.amplitude.ll,self.amplitude.ll_iv))
        atoms['log_delta']=self._summary((-4*self.amplitude.lp-30,-4*self.amplitude.lp_iv-30))
        atoms['radial_geometry']=self._summary((c.log(2)/2,iv.ln(2)/2))
        atoms['swirl_geometry']=self._summary((c.log(8)/2,iv.ln(8)/2))
        atoms['half_hbs']=self._selected_geometry(s)
        atoms['original_log_hb']=self._log_width()
        atoms['log_phase_in']=self._phase_log()
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
            if source in ('f','u','q'):
                width_logs={'original_log_hb':-1,'log_phase_in':-1}
                radial_singular=pair(lambda v,ctx,zz:2*v['rho']*v['radial_singular'][source])
                axial_singular=pair(lambda v,ctx,zz:-2*v['rho']*zz*v['radial_singular'][source])
                dr[component].append(emit(component+'_r_width:'+source,order,radial_singular,
                    adjust(adjust(base,minus_radius),width_logs),'radial_cutoff_width_derivative',source))
                dz[component].append(emit(component+'_z_width:'+source,order,axial_singular,
                    adjust(zlogs,width_logs),'axial_cutoff_width_derivative',source))
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
        velocity_cart=[]
        for i in range(3):
            row=[]
            for j,name in enumerate(('ur','utheta','uz')):
                f=R[i][j]
                if f[0]==0 and _bounds(c,f[1])==(0,0):continue
                for ref in refs(values[name]):
                    p,t=self._pair(ref['factor'])
                    row.append(dict(term=ref['term'],factor=self._summary((p*f[0],t*f[1]))))
            velocity_cart.append(row)
        divergence={}
        for name,U,Q,beta in (('leading_axis','U0','Qaxis',-1),('leading_correction','Psi','Qpsi',-1),('n1','u','q',1)):
            def defect(v,ctx,zz):
                u,uz,ur=v['profiles'][U];q,_,qr=v['profiles'][Q]
                return v['L']*(q+v['rho']*qr)+v['d']*uz-2*v['rho']*zz*ur+(-1+beta*v['delta'])*zz*u
            divergence[name]={'0,0':self._summary(pair(defect))}
            if name=='n1':
                divergence[name]['-1,-1']=self._summary(pair(lambda v,ctx,zz:
                    v['L']*v['rho']*v['radial_singular']['q']
                    -2*v['rho']*zz*v['radial_singular']['u']))
        native_rows={}
        for name in point['profiles']:
            row=dict(value=self._summary((point['profiles'][name][0],truth['profiles'][name][0])),
                Z=self._summary((point['profiles'][name][1],truth['profiles'][name][1])),
                rho_width_sectors={'0,0':self._summary((point['profiles'][name][2],truth['profiles'][name][2]))})
            if name in ('f','u','q'):
                row['rho_width_sectors']['-1,-1']=self._summary((point['radial_singular'][name],truth['radial_singular'][name]))
            native_rows[name]=row
        packet=data['cutoff_packet']
        return dict(native=dict(chart=packet['chart'],Z=Z,fraction=s,azimuth=azimuth,
                fraction_definition='phase/(original_sc/4)',
                exact_sc_mpf_tuple=packet['geometry']['exact_sc_mpf_tuple'],
                exact_phase_in_mpf_tuple=packet['geometry']['exact_phase_in_mpf_tuple']),
            time=time,source_family=self.family,log_atoms=atoms,term_pool=pool,
            cylindrical_velocity={name:refs(rows) for name,rows in values.items()},
            Cartesian_velocity=dict(zip(('u','v','w'),velocity_cart)),
            gradient_cylindrical=J,gradient_Cartesian=rotate_matrix(),
            vorticity_cylindrical=dict(zip(('r','theta','z'),omega)),
            vorticity_Cartesian=dict(zip(('x','y','z'),omega_cart)),
            normalized_meridional_divergence_width_sectors=divergence,native_source_rows=native_rows,
            own_history=dict(primitive=packet['normalized_own_primitive'],
                cutoff=packet['cutoff'],only_f_and_u_cut=True,q_recovered_from_same_retained_mean=True),
            physical_chain=dict(partial_r='(2rho/r)*partial_rho',
                partial_z='(d*partial_Z-2rho*Z*partial_rho)/(lambda^(1-delta)*L)',
                F0_log_derivative='-g/epsilon; its Lambda factor is a separate log atom',
                radial_sector='a,b denotes original_hb^a*original_phase_in^b times coefficient'),
            semantics='Each reference is factor*coefficient*exp(sum multiplier*log_atom). Retain all source intervals and exact-source atoms; missing approximations are not zero. Point rows are numerical proxies, with omitted finite-width contributions charged to absolute errors.',
            original_F0_gauge_delta_and_cutoff_width_scales_kept_separate=True,
            limits=dict(actual_selected_cut_and_own_physical_field=True,nonzero_radius=True,
                inner_fraction_range=['0','1'],lower_fraction_range=['1','2'],
                stops_at_rminus=True,raw_uncut_domain_extended=False,
                physical_first_spatial_derivatives=True,physical_second_spatial_derivatives=False,
                ordinary_float_velocity_gradient=False,useful_combined_physical_relative_error=False,
                complete_Cartesian_C2_join=False,global_field=False,complete_five_moment_repair=False,
                temporal_n_ge_2_recursion=False,full_corrected_NS=False))
