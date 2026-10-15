"""Actual selected-N first-micro modification with directed numerical errors.

The signed periodic source remains the accepted Section11 function. Full
phase bounds enclose its numerical proxy, without choosing a fast angle.
The exact repair-only DoubleDyadicInteger is retained, never materialized.
"""
import gzip
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_collar_velocity import _InletGeometry
from .original_selected_gradient import OriginalSelectedGradient


class OriginalFirstMicroNumeric(OriginalSelectedGradient):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.records={}
        prerequisites=(
            ('compliant_current_original_C3_target_ranges',
             ('current_original_actual_17_chart_C3_primitive_density_ranges_installed',
              'current_original_actual_four_Z_C3_five_target_ranges_installed')),
            ('compliant_current_original_C3_repair_band',
             ('current_original_same_limit_repair_band_C3_velocity_histories_installed',)),
            ('compliant_current_original_prescribed_loop_C5_functions',
             ('current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed',)))
        for stem,gates in prerequisites:
            receipt=self.leading.core.directory/(PREFIX+stem+'_check.json')
            checked=json.loads(receipt.read_text(encoding='utf-8'))
            if not checked['all_passed'] or checked['source_family']!=self.family or not all(checked[k] for k in gates):
                raise ValueError('Accepted same-family modified source/range/integer required: '+stem)
            for name,digest in checked['input_hashes'].items():
                if hashlib.sha256((self.leading.core.directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed modified numerical dependency: '+name)
                if self.hashes.get(name,digest)!=digest:raise ValueError('Conflicting modified source: '+name)
                self.hashes[name]=digest
            self.hashes[receipt.name]=hashlib.sha256(receipt.read_bytes()).hexdigest()
            # The loop archive is source-pinned by its receipt. Only the two
            # quantitative archives are hydrated; no ancestor is rebuilt.
            if 'prescribed_loop' not in stem:
                path=self.leading.core.directory/(PREFIX+stem+'.json.gz')
                self.records[stem]=json.loads(gzip.decompress(path.read_bytes()))
                if self.records[stem]['source_family']!=self.family:raise ValueError('Range archive family changed')
                self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self.ranges=self.records['compliant_current_original_C3_target_ranges']['actual_four_Z_C3_target_transports']
        repair=self.records['compliant_current_original_C3_repair_band']
        self.selection=repair['actual_selected_repair_integer']
        integers=[row for row in repair['exact_graph_nodes'] if row.get('operation')=='exact_positive_integer_expression']
        if (len(integers)!=1 or integers[0]['type']!='DoubleDyadicInteger'
                or integers[0]['J']!=self.selection['J']
                or not integers[0]['one_common_N_for_all_charts_and_control_denominators']
                or self.selection['J']<12 or not self.selection['N_ge_N0']
                or self.selection['N0_bit_length']>4096):
            raise ValueError('One exact original repair-only integer and N>=N0 proof required')
        self.integer=integers[0]
        c,iv=self.transport.ctx,self.transport.iv
        logtwo=iv.ln(2);J=self.integer['J']
        self.logN=iv.mpf([mp.make_mpf((sgn,man,exp+J,bits)) for sgn,man,exp,bits in logtwo._mpi_])
        stored=self._decode(self.selection['log_N'])
        if max(_bounds(c,stored)[0],_bounds(c,self.logN)[0])>min(_bounds(c,stored)[1],_bounds(c,self.logN)[1]):
            raise ValueError('Selected integer recipe and stored logN disagree')
        self._micro_cache={};self._terminal_cache={}
        relative='src/openai_ns_reconstruction/original_first_micro_numeric.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()

    def _decode(self,row):
        return self.transport.iv.mpf([mp.make_mpf(tuple(row[k+'_exact_mpf_tuple'])) for k in ('lower','upper')])

    def _source_log(self,row):
        return None if row['exact_zero'] else self._decode(row['log_absolute_upper'])

    def _source_summary(self,value):
        return dict(exact_zero=value is None,log_absolute_upper=None if value is None else dict(
            lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1])))

    def _sum_logs(self,logs):
        logs=[v for v in logs if v is not None]
        if not logs:return None
        c,iv=self.transport.ctx,self.transport.iv
        high=max(_bounds(c,v)[1] for v in logs)
        # A positive sum bound. It does not imply signed cancellation.
        return iv.mpf(high)+iv.ln(len(logs))

    def _cap(self,log,proof,label):
        iv,c=self.transport.iv,self.transport.ctx
        if log is None:return iv.mpf(0)
        bound=iv.mpf('1e-1000')
        if _bounds(c,log)[1]>_bounds(c,iv.ln(bound))[0]:
            raise ArithmeticError('Modified normalized error is not bounded by 1e-1000: '+label)
        proof.append(dict(quantity=label,original_log_upper=self._source_summary(log),
            numerical_absolute_error_cap='1e-1000',cap_is_bound_not_source_value=True))
        return bound

    def _coordinate(self,coordinate,selected):
        c,iv,t=self.transport.ctx,self.transport.iv,self.transport
        if type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Finite decimal source coordinate required')
        if selected:
            if not c.mpf('.5')<=c.mpf(coordinate)<=1:raise ValueError('Selected sc multiple in[1/2,1] required')
            sc=t.phase_iv*4;s=sc*iv.mpf(coordinate)
            flat=_bounds(c,iv.mpf(coordinate))[1]<=c.mpf('.75')
            recipe=dict(kind='selected_sc_multiple',multiple=coordinate,exact_sc_mpf_tuple=t.sc_tuple)
        else:
            s=iv.mpf(coordinate);lo,hi=_bounds(c,s)
            if lo<_bounds(c,t.phase_iv*2)[1] or hi>1:raise ValueError('Actual first_micro phase in[original_sc/2,1] required')
            flat=hi<=_bounds(c,t.phase_iv*3)[0]
            recipe=dict(kind='exact_decimal_phase',value=coordinate)
        return s,flat,recipe

    def _endpoint(self,raw):
        iv=self.transport.iv
        return self.transport.phase_iv*2 if raw==['sc_half'] else iv.mpf(raw[0])/raw[1]

    def _covers(self,Z,s):
        iv,c=self.transport.iv,self.transport.ctx;z=iv.mpf(Z)
        zl,zh=_bounds(c,z);sl,sh=_bounds(c,s);local=[];prefix=[];frames=[]
        for frame in self.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            frames.append(frame['exact_Z_cell'])
            chart=next(q for q in frame['actual_17_chart_C3_transports'] if q['chart']=='first_micro')
            for cell in chart['cells']:
                left,right=self._endpoint(cell['exact_left']),self._endpoint(cell['exact_right'])
                ll,lh=_bounds(c,left);rl,rh=_bounds(c,right)
                if not rh<sl and not ll>sh:local.append(cell)
                if sh>ll:
                    width=iv.mpf(max(c.mpf(0),_bounds(c,iv.mpf(min(sh,rh))-iv.mpf(ll))[1]))
                    if width:prefix.append((cell,width))
        if not local:raise ValueError('No complete admitted first-micro source cover')
        return local,prefix,frames

    def _background(self,Z,s):
        o,t=self.leading,self.transport;iv=o.iv;order=3
        inlet=o.inlet_jets(Z,degree=24,order=order)
        if Z not in o._prepared:
            with mp.workdps(280):o._prepared[Z]=o.bridge.prepare(Z)
        upper=_bounds(t.ctx,s)[1]
        _,Bell,_,Bpsi=o._source_jet_increments(o._prepared[Z],upper,order)
        a=inlet['truth']['values']
        Bphi=o._norm(a['Phi'])*Bell*iv.exp(o._width_error(Bell))
        Bmean=iv.mpf(upper)*o._norm(a['Psi']-a['Mpsi'])+o._width_error(iv.mpf(upper)*Bpsi)
        truth={name:a[name]+o._norm_error_jet(o._width_error(bound),order)
               for name,bound in (('Phi',Bphi),('Psi',Bpsi),('Mpsi',Bmean))}
        rho=iv.mpf(4)+iv.mpf([0,_bounds(t.ctx,o._width_error(4*iv.mpf(upper)*iv.exp(o.hcap*iv.mpf(upper))))[1]])
        return inlet['point']['values'],truth,rho

    def _correction_logs(self,local,prefix,flat):
        iv,o=self.transport.iv,self.leading
        A=[];B=[];hist={name:[] for name in ('m','h','k','e','p')}
        for k in range(4):
            alogs=[self._source_log(cell['actual_source_and_primitive_range']['A'][k]) for cell in local]
            blogs=[self._source_log(cell['actual_source_and_primitive_range']['B'][k]) for cell in local]
            aa=self._sum_logs(alogs);bb=self._sum_logs(blogs)
            A.append(None if flat or aa is None else aa-self.logN)
            B.append(None if flat or bb is None else bb-self.logN+o.core.lp_iv+o.core.ll_iv)
            for name in hist:
                logs=[]
                for cell,width in prefix:
                    v=self._source_log(cell['actual_density_ranges']['N_scaled'][name][k])
                    if v is not None:logs.append(v+iv.ln(width))
                bound=self._sum_logs(logs)
                hist[name].append(None if flat or bound is None else bound+o.logh-self.logN)
        return A,B,hist

    @staticmethod
    def _exponential_increment_bounds(x):
        return [3*x[0],3*x[1],3*(x[2]+x[1]*x[1]/2)]

    def _numeric(self,Z,coordinate,selected):
        c,iv,o,t=self.transport.ctx,self.transport.iv,self.leading,self.transport
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:
            raise ValueError('Interior finite decimal |Z|<1 required')
        key=Z,coordinate,selected
        if key in self._micro_cache:return self._micro_cache[key]
        s,flat,recipe=self._coordinate(coordinate,selected)
        local,prefix,frames=self._covers(Z,s);bp,bt,rho=self._background(Z,s)
        proof=[];zero=o.Jet.constant(iv,0,2)
        A,B,hist=self._correction_logs(local,prefix,flat)
        # First two ordinary axial derivatives have a finite small error.
        # Third rows remain logarithmic; the selected repair-only N does not
        # make all conservative B_ZZZ/history caps numerically small.
        x=[self._cap(A[k],proof,'A_over_N_Z'+str(k))/math.factorial(k) for k in range(3)]
        exp_error=self._exponential_increment_bounds(x)
        delta_phi=bt['Phi'].truncate(2)*o.Jet(iv,[o._symmetric(v) for v in exp_error])
        delta_psi=o.Jet(iv,[o._symmetric(self._cap(B[k],proof,'delta_Psi_Z'+str(k))/math.factorial(k)) for k in range(3)])
        meanlogs=[None if v is None else v+o.core.lp_iv+o.core.ll_iv for v in hist['m']]
        delta_mean=o.Jet(iv,[o._symmetric(self._cap(meanlogs[k],proof,'delta_Mpsi_Z'+str(k))/math.factorial(k)) for k in range(3)])
        point={name:bp[name].truncate(2) for name in ('Phi','Psi','Mpsi')}
        truth={name:bt[name].truncate(2)+inc for name,inc in
               (('Phi',delta_phi),('Psi',delta_psi),('Mpsi',delta_mean))}
        def recover(values,ctx,de):
            z=o.Jet.variable(ctx,ctx.mpf(Z),2);d=1-z*z;L=1-z*z*de
            return ((2*z*values['Psi']-z*values['Mpsi']*(1-de)-d*o._derivative(values['Mpsi']))/L).truncate(1)
        qp=recover(point,c,o.core.delta);qt=recover(truth,iv,o.core.delta_iv)
        out={name:o._jet_summary(point[name],truth[name]) for name in point}
        out['Qpsi']=o._jet_summary(qp,qt)
        result=dict(Z=Z,coordinate=recipe,source_family=self.family,
            normalized_modified_leading=out,
            modified_error_contributions=dict(Phi=o._jet_summary(o.Jet.constant(c,0,2),delta_phi),
                Psi=o._jet_summary(o.Jet.constant(c,0,2),delta_psi),
                Mpsi=o._jet_summary(o.Jet.constant(c,0,2),delta_mean)),
            selected_integer=dict(recipe=self.integer,original_N_not_reselected=True,
                log_N_exact_mpf_tuples=[list(v) for v in self.logN._mpi_],
                N_ge_original_N0=True,global_high_derivative_stress_heat_admission=False),
            source_cover=dict(exact_Z_cells=frames,local_source_cells=len(local),prefix_source_cells=len(prefix),
                whole_original_phase_cover=[0,1],fast_phase_not_chosen=True,
                signed_source='Phi_N=Phi*exp(A/N), U_N=U+Pstar*B_over_Pstar/N',
                true_phase='fractional_part(original_N*original_hb*(s-original_sc/2))',
                defining_signed_C3_receipt=PREFIX+'compliant_current_original_C3_source_targets_check.json',
                true_prescribed_C5_receipt=PREFIX+'compliant_current_original_prescribed_loop_C5_functions_check.json',
                background_incoming_history_retained=True,exact_flat_collar=flat),
            original_five_correction_history_bounds={name:[self._source_summary(v) for v in rows] for name,rows in hist.items()},
            history_contract=dict(rate_m=1,measure='original_hb*ds once',
                partial_prefix='integral_sc/2^s exp(-hb*(s-t))*B(t,Z,phase)/N * hb dt',
                bound='positive rate kernels<=1; sum partial-cell lengths times complete original density caps',
                all_five_signed_source_rates_retained=True,only_original_rminus_correction_inlet_zero=True,
                caps_are_not_signed_point_values_or_exact_integrals=True),
            normalized_third_derivative_log_error_bounds=dict(A_over_N=self._source_summary(A[3]),
                delta_Psi=self._source_summary(B[3]),delta_Mpsi=self._source_summary(meanlogs[3])),
            small_error_proofs=proof,coefficient_convention='Taylor Zk=D_Z^k/k!; cap inputs ordinary Z derivatives',
            semantics='Original signed functions plus bounded numerical proxies. Full-phase caps are errors, never field values or freely chosen phases.',
            limits=dict(actual_modified_first_micro_values_and_Z2=True,Qpsi_Z1=True,
                resolved_fast_phase_or_loop_point=False,actual_velocity_Z3_numeric=False,Qpsi_ZZ_numeric=False,
                physical_radial_derivatives=False,physical_C2_join=False,global_field=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._micro_cache[key]=(result,point,truth,rho,s)
        return self._micro_cache[key]

    def _physical(self,numeric,Z,log_tau,scaled_time,azimuth):
        packet,point,truth,rho,s=numeric;c,iv=self.ctx,self.iv;o,t=self.leading,self.transport
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):raise ValueError('Finite decimal azimuth required')
        read=lambda row:self._pair(row)
        rows=packet['normalized_modified_leading']
        z,zi=c.mpf(Z),iv.mpf(Z);de,di=self.amplitude.delta,self.amplitude.delta_iv
        u0,u0i=4*z+self.amplitude.j,4*zi+self.amplitude.j_iv
        profiles={name:read(rows[source][0]) for name,source in
                  (('Phi','Phi'),('Psi','Psi'),('mean_Psi','Mpsi'),('Qpsi','Qpsi'))}
        profiles.update(U0=(u0,u0i),Qaxis=(((1+de)*z*u0-4*(1-z*z))/(1-de*z*z),
            ((1+di)*zi*u0i-4*(1-zi*zi))/(1-di*zi*zi)))
        if Z not in self._terminal_cache:self._terminal_cache[Z]=t.lower_collar(Z,'2')
        terminal=self._terminal_cache[Z]['normalized_own_primitive']['actual_primitive']
        J=t.Jet
        pp=J(t.ctx,[t.ctx.mpf(row['approximation']) for row in terminal])
        pt=J(t.iv,[t.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]) for row in terminal])
        mp0,mt=pp/4,pt/rho
        def q(mean,ctx,delta):
            zz=J.variable(ctx,ctx.mpf(Z),2)
            return (-zz*mean*(1+delta)-(1-zz*zz)*o._derivative(mean))/(1-zz*zz*delta)
        qp,qt=q(mp0,t.ctx,o.core.delta),q(mt,t.iv,o.core.delta_iv)
        transfer=lambda p,b:(c.mpf(p),iv.mpf([mp.make_mpf(tuple(v)) for v in b._mpi_]))
        n1=dict(f=(c.mpf(0),iv.mpf(0)),u=(c.mpf(0),iv.mpf(0)),
                mean_u=transfer(mp0[0],mt[0]),q=transfer(qp[0],qt[0]))
        adapter=_InletGeometry(self,profiles,n1)
        angle=c.mpf(azimuth);anglei=iv.mpf(azimuth)
        out=adapter._assemble(Z,'4',log_tau,scaled_time,((c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))))
        def geometry(multiplier):
            high=_bounds(t.ctx,o._width_error(s*t.iv.mpf(multiplier)))[1]
            row=self._summary((c.mpf(0),iv.mpf([0,high])))
            row.update(kind='original_hb_linear_log',definition='multiplier*original_hb*actual_phase',
                multiplier=multiplier,phase_recipe=packet['coordinate'],original_log_hb=self._log_width(),
                zero_is_numeric_proxy_only=True)
            return row
        half,whole=geometry('.5'),geometry('1')
        for key in ('ur','utheta'):
            for term in out['cylindrical'][key]:term['log_scale_terms'].append(half)
        for key in ('u','v'):
            for term in out['Cartesian'][key]:term['log_scale_terms'].append(half)
        for key in ('r','x','y'):out['coordinates'][key]['log_scale_terms'].append(half)
        for term in out['streamfunction']:term['log_scale_terms'].append(whole)
        out['native']=dict(chart='actual_modified_first_micro',Z=Z,phase_recipe=packet['coordinate'],
            rho_definition='4exp(original_hb*actual_phase)',rho4_in_coefficients_is_numeric_proxy=True)
        out['original_selected_integer']=packet['selected_integer']
        out['actual_modified_leading_errors_included']=True
        out['n1_own_history']=dict(actual_terminal_primitive=self._terminal_cache[Z]['normalized_own_primitive'],
            mean_u=self._summary(n1['mean_u']),q=self._summary(n1['q']),
            primitive_retained_through_current_radius=True,f_and_u_exact_zero=True,
            no_later_reserved_repair_support_on_first_micro=True)
        out['physical_gradient_or_vorticity_available']=False
        return out

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        """Actual first-micro native phase sc/2<=s<=1; original N fixed."""
        numeric=self._numeric(Z,s,False)
        return dict(**numeric[0],physical=self._physical(numeric,Z,log_tau,scaled_time,azimuth))

    def selected(self,Z,multiple,*,log_tau='-1',scaled_time=None,azimuth='0'):
        """Exact selected phase=multiple*sc, multiple in[1/2,1]."""
        numeric=self._numeric(Z,multiple,True)
        return dict(**numeric[0],physical=self._physical(numeric,Z,log_tau,scaled_time,azimuth))

    def lower_collar(self,*args,**kwargs):
        raise ValueError('Use OriginalSelectedGradient for the preceding selected lower collar')
