"""Original end-pulse n1 pressure values with correlated moment rows.

Actual coefficient functions and their C5 ranges are reused on four axial
cells. Linear pressure is evaluated from exact finite-mu row identities.
Quadratic source terms and the weighted viscous endpoint remain as genuine,
nonzero bounded contributions in fixed physical source units.
"""
import copy
import gzip
import json
import time
from pathlib import Path

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_coefficients_refined as coefficients
import lei_ren_part1_paper_compliant_current_original_n1_repair_band_pressure_source as radius
import lei_ren_part1_paper_compliant_axial_pulse_field as original_end
import lei_ren_part1_paper_compliant_pulse_high_jets as original_radial
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE,PREFIX,sha,ep=coefficients.HERE,coefficients.PREFIX,coefficients.sha,coefficients.ep
source_precision=coefficients.source_precision
record=coefficients.record
rebox=coefficients.rebox
NAME=PREFIX+'current_original_n1_end_pressure_source.json.gz'
RECEIPT=PREFIX+'current_original_n1_end_pressure_source_check.json'
GATE='current_original_n1_end_pressure_source_installed'
OPEN=('complete_signed_pressure_tails_certified','n1_second_and_third_repair_closed',
      'global_physical_pressure_defect_installed','uniform_full_pressure_derivative_accuracy_certified',
      'global_stress_certified','global_NS_solution_certified')
CELLS=(('-1','-.5'),('-.5','0'),('0','.5'),('.5','1'))
MODELS=256


def derivative(jet):
    return IntervalTaylor(jet.ctx,[(i+1)*jet[i+1] for i in range(jet.order)])


def read_jet(c,item):
    return IntervalTaylor(c,[coefficients.basis.read_interval(c,v) for v in item['coefficients']])


def copy_view(c,value):
    if isinstance(value,IntervalTaylor):return rebox(c,value)
    if hasattr(value,'_mpi_'):return rebox(c,value)
    if isinstance(value,dict):return {key:copy_view(c,item) for key,item in value.items()}
    if isinstance(value,(tuple,list)):return [copy_view(c,item) for item in value]
    return value


def abs_upper(value):
    lo,hi=ep(value)
    return max(abs(lo),abs(hi))


def scalar_bound(c,value):
    return c.mpf([0,ep(value)[1]])


def polynomials(c,Z,order,delta):
    z=c.mpf(Z)
    variable=IntervalTaylor.variable(c,z,order)
    def quadratic(value,first,second):
        return IntervalTaylor(c,[value,first,second]+[0]*(order-2))
    Q=quadratic(1+z**2,2*z,c.mpf(1))
    d=quadratic(1-z**2,-2*z,c.mpf(-1))
    L=quadratic(1-delta*z**2,-2*delta*z,-delta)
    return variable,Q.reciprocal(),d,L


def radial_q(c,Z,delta,U,m):
    order=min(U.order,m.order-1)
    z,_,d,L=polynomials(c,Z,order,delta)
    return (U.truncate(order)*z*2-m.truncate(order)*z*(1-delta)-derivative(m)*d)/L


def end_linear_jets(c,Z,delta,mu,ap,C,rows,Pi,Wplus,rho):
    """Exact linear source plus directed nonzero finite-mu X remainder."""
    z,ubar,_,L=polynomials(c,Z,5,delta)
    G1=rows[0]+ap*Pi[0]
    G2=rows[1]-rows[0]+ap*(Pi[1]-Pi[0])
    Wsum=-G1*mu
    rho_sum=sum((j*w for j,w in zip(C,rho)),ap*0)
    Xsum=Wsum*4+G2+rho_sum
    Wsum_plus=sum((j*w for j,w in zip(C,Wplus)),ap*0)
    I0U=ubar*Wsum;I0m=-ubar*Xsum
    I1U=ubar*Wsum_plus
    I1m=-ubar*(Wsum_plus-Wsum*c.exp(-4))
    q0=radial_q(c,Z,delta,I0U,I0m)
    q1=radial_q(c,Z,delta,I1U,I1m)
    qbase=radial_q(c,Z,delta,ap*0,-ubar*Wsum)
    L3=L.truncate(3);z3=z.truncate(3)
    S0=-qbase.truncate(3)/L3+derivative(q0)*z3*((1-delta)/2)/L3
    S1=(-qbase.truncate(3)*c.exp(-4)-q1.truncate(3))/L3
    S1+=derivative(q1)*z3*((1-delta)/2)/L3
    return dict(m13_linear=-S0/2,m15_linear=S1/2,qbase=qbase,
        I0U=I0U,I0m=I0m,I1U=I1U,I1m=I1m,I0q=q0,I1q=q1,
        first_row_sum=Wsum,Xsum=Xsum,X_finite_mu_remainder=rho_sum,
        lower_q_plus_R_qR_exact_zero=True,
        m15_viscous_coefficient=-qbase)


def end_field_bounds(c,delta,mu,C,normalization):
    """Genuine whole-band U,q,q_Z bounds in the common unit A."""
    zbox=['-1','1'];_,ubar,_,_=polynomials(c,zbox,5,delta)
    beta=scalar_bound(c,c.exp(-1)/(c.mpf('.15')*normalization))
    h=c.mpf([-ep(c.exp(4*(c.mpf('.5')-mu)))[1],0])
    decay=c.mpf([1,ep(c.exp(4*(c.mpf('.5')+mu)))[1]])
    U=ubar*sum((j*beta for j in C),C[0]*0)*decay
    m=ubar*sum((j*h for j in C),C[0]*0)*decay
    q=radial_q(c,zbox,delta,U,m)
    return dict(U=c.mpf(abs_upper(U[0])),UZ=c.mpf(abs_upper(derivative(U)[0])),
                q=c.mpf(abs_upper(q[0])),qZ=c.mpf(abs_upper(derivative(q)[0])),
                beta_upper=beta,h_range=h,radial_decay_range=decay)


def quadratic_error_bounds(c,bounds,qbase_bound,A_cap,normalization):
    """All nonzero quadratic reduced-S sources, including boundaries."""
    U,UZ,q,qZ=(bounds[key] for key in ('U','UZ','q','qZ'))
    # delta<1/2, L^-1<=2, |Z|<=1, 0<=d<=1; integrate e^s<=1,
    # and integrate e^(2s)<=1/2. U is zero at both finite endpoints.
    S0=c.exp(4)*qbase_bound**2/2+2*U*qZ+4*U**2+8*U*UZ
    S1=qbase_bound**2/2+q**2/4+U*qZ+4*U**2+4*U*UZ
    return dict(m13=scalar_bound(c,A_cap*S0/(2*normalization)),
                m15=scalar_bound(c,A_cap*S1/(2*normalization)),
                reduced_S0_absolute_bound=S0,reduced_S1_absolute_bound=S1,
                quadratic_sources_retained=True)


class CurrentOriginalN1EndPressureSource:
    @source_precision
    def __init__(self,coefficient_owner=None,require_checked=True):
        self.coefficients=coefficient_owner or coefficients.CurrentOriginalRpPulseCoefficientsRefined()
        if type(self.coefficients) is not coefficients.CurrentOriginalRpPulseCoefficientsRefined or not self.coefficients.acceptance_loaded:
            raise ValueError('Checked original refined coefficient owner required')
        self.coefficients.assert_graph();self.c=c=self.coefficients.c
        self.family=copy.deepcopy(self.coefficients.family);self.mu=self.coefficients.mu
        self.delta=rebox(c,self.coefficients.selected.pulse.delta)
        if not (0<ep(self.mu)[0]<=ep(self.mu)[1]<c.mpf('.01') and
                0<ep(self.delta)[0]<=ep(self.delta)[1]<c.mpf('.5')):
            raise ValueError('Actual positive mu/delta source regime required')
        self.radius_binding=radius._accepted_native_radius_binding(self.family)
        logC=coefficients.basis.read_interval(c,self.radius_binding['selected_logCstar'])
        self.logRc=c.ln(110)+10*logC+11*c.exp(40)+124
        self.logRp=self.logRc-2-60*c.ln(self.mu)
        self.logP=c.exp(40)+11
        self.logUp=c.ln(rebox(c,self.coefficients.selected.constants['U']))
        logE=rebox(c,self.coefficients.selected.amplitude.log_end_scale)
        self.logA=self.logP+self.logUp+logE-13/(2*self.mu)-13
        self.logRv=self.logRp+13/self.mu
        if ep(self.logA)[1]>=-1000 or ep(self.logRv)[0]<=1000:
            raise ValueError('Original quadratic/viscous source cap not proved')
        self.A_cap=c.exp(-1000);self.invRv=scalar_bound(c,c.exp(-1000))
        self.reference=c.mpf('10.99')*c.sqrt(2*c.pi)*(c.exp(5)-c.exp(3))
        B=self.coefficients.basis
        self.Pi=[rebox(c,v['refined']) for v in B.saddle_owner.refined_rows]
        A=[rebox(c,v) for v in B.basis['first_row']]
        J=B.basis['moments']['1.5'];ell=c.mpf('.15')
        remainder=scalar_bound(c,self.mu**2*ell**2*c.exp(self.mu*ell)*J['J']/2)
        Jplus=J['J']-J['J_r']*self.mu+c.mpf([-ep(remainder)[1],ep(remainder)[1]])
        self.Wplus=[c.exp((c.mpf('1.5')-self.mu)*(-d))*Jplus for d in (3,1)]
        self.rho=[scalar_bound(c,a*self.mu*b**2*c.exp(self.mu*b)/2)
                  for a,b in zip(A,(c.mpf('3.15'),c.mpf('1.15')))]
        self._cells=[];self._packets={}
        for ends in CELLS:
            packet=self.coefficients.evaluate(list(ends));view=self.coefficients.report(packet)
            jets=self.coefficients._packets[id(packet)][2]
            inputs=view['source_inputs'];factors=inputs['fixed_incoming_row_factor_caps']
            rows=[read_jet(c,j)*coefficients.basis.read_interval(c,f) for j,f in
                  zip(inputs['incoming_moment_Taylor'],factors)]
            cell=dict(ends=ends,ap=jets['ap'],C=[jets['C1'],jets['C2']],rows=rows)
            self._cells.append(cell)
        qbounds=[];fieldbounds=[]
        for cell in self._cells:
            terms=end_linear_jets(c,list(cell['ends']),self.delta,self.mu,cell['ap'],cell['C'],cell['rows'],self.Pi,self.Wplus,self.rho)
            qbounds.append(abs_upper(terms['qbase'][0]))
            fieldbounds.append(end_field_bounds(c,self.delta,self.mu,cell['C'],B.normalization))
        global_bounds={key:c.mpf(max(ep(item[key])[1] for item in fieldbounds))
                       for key in ('U','UZ','q','qZ')}
        self.qbase_bound=c.mpf(max(qbounds))
        self.errors=quadratic_error_bounds(c,global_bounds,self.qbase_bound,self.A_cap,self.reference)
        self.viscous_error=scalar_bound(c,self.invRv*self.qbase_bound/self.reference)
        self.full_errors={key:self.errors[key]+(self.viscous_error if key=='m15' else 0)
                          for key in ('m13','m15')}
        node=self.coefficients.selected.current_Rp_physical_log_radius.node
        if node!=self.radius_binding['functions']['original_logRp']:
            raise ValueError('Original absolute Rp source node differs')
        self.units=dict(source_Rp_node=node,source_native_frame=radius.NATIVE_FRAME_NAME,
            source_logRc_formula='log(110)+10*logCstar+11*exp(40)+124',
            source_logRp_formula='logRc-2-60*log(mu)',
            source_Rv_formula='Rp*exp(13/mu)',source_A='Pstar*U0*exp(logE-13/(2*mu)-13)',
            source_logE_formula='-1/mu-3/u0^2+2*log(u0)-log(6)/2-2*log(mu); u0=(4mu)^(1/3)',
            inverse_mu_integer_powers_include_logE_leading_term=True,
            normalization='10.99*sqrt(2*pi)*(exp(5)-exp(3))',
            m13=dict(A_power=1,Rv_power=1,inverse_mu_numerator=11,inverse_mu_denominator=2),
            m15=dict(A_power=1,Rv_power=2,inverse_mu_numerator=37,inverse_mu_denominator=2),
            quadratic_m13=dict(A_power=2,Rv_power=1,inverse_mu_numerator=-4,inverse_mu_denominator=2),
            quadratic_m15=dict(A_power=2,Rv_power=2,inverse_mu_numerator=22,inverse_mu_denominator=2),
            integer_powers_combined_before_boxing=True,no_physical_scale_materialization=True)
        self.source_bindings={name:coefficients.function_binding(fn) for name,fn in (
            ('original_end_profile',original_end.CompliantAxialPulseField.end),
            ('original_radial_recovery',original_radial.CompliantPulseHighJets.radial),
            ('end_linear_jets',end_linear_jets),('radial_q',radial_q),
            ('quadratic_error_bounds',quadratic_error_bounds))}
        self.bounds=dict(global_field_bounds=global_bounds,qbase_bound=self.qbase_bound,
            log_A_upper_source=self.logA,log_Rv_lower_source=self.logRv,
            positive_A_cap=self.A_cap,positive_inverse_Rv_enclosure=self.invRv,
            quadratic=self.errors,viscous_error=self.viscous_error,reference=self.reference)
        self.hashes=dict(self.coefficients.hashes)
        for name,digest in self.radius_binding['receipt']['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting radius dependency')
            self.hashes[name]=digest
        for name in (coefficients.NAME,coefficients.RECEIPT,radius.NATIVE_FRAME_NAME,
                     radius.NATIVE_FRAME_RECEIPT,radius.NORM_NAME,Path(radius.__file__).name,Path(__file__).name):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Checked end-pressure source receipt required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed end-pressure input '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._state=record((self.family,self.units,self.bounds,self.Wplus,self.rho,
            self._cells,self.mu,self.delta,self.Pi,self.source_bindings))
        self._hashes=copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        self.coefficients.assert_graph()
        if self._state!=record((self.family,self.units,self.bounds,self.Wplus,self.rho,
                self._cells,self.mu,self.delta,self.Pi,self.source_bindings)) or self.hashes!=self._hashes:
            raise ValueError('Original end-pressure source state changed')

    def _cell(self,Z):
        lo,hi=ep(self.c.mpf(Z))
        for cell in self._cells:
            left,right=map(mp.mpf,cell['ends'])
            if left<=lo<=hi<=right:return cell
        raise ValueError('Axial query must lie in an accepted original coefficient cell')

    def _make_view(self,Z):
        c=self.c;cell=self._cell(Z)
        terms=end_linear_jets(c,Z,self.delta,self.mu,cell['ap'],cell['C'],cell['rows'],self.Pi,self.Wplus,self.rho)
        linear={key:terms[key+'_linear']/self.reference for key in ('m13','m15')}
        full={key:linear[key][0]+c.mpf([-ep(self.errors[key])[1],ep(self.errors[key])[1]]) for key in linear}
        viscous=-terms['qbase'][0]*self.invRv/self.reference
        full['m15']+=viscous
        view=dict(Z=c.mpf(Z),source_family=self.family,source_coefficient_cell=list(cell['ends']),
            full_pressure_values=full,linear_pressure_Taylor=linear,source_linear_inputs=terms,
            original_physical_units=self.units,quadratic_source_errors=self.errors,
            weighted_viscous_endpoint=viscous,unweighted_viscous_boundary_exact_zero=True,
            exact_source_row_correlation_retained=True,finite_mu_X_remainder_retained=True,
            full_pressure_derivative_accuracy_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        return view

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();view=self._make_view(Z)
        packet=object();self._packets[id(packet)]=(packet,view,record(view))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();item=self._packets.get(id(packet))
        if item is None or item[0] is not packet or item[2]!=record(item[1]):
            raise ValueError('Live unchanged original end-pressure packet required')
        return copy_view(self.c,item[1])

    @source_precision
    def affine_models(self,count=MODELS):
        self.assert_graph()
        if count<4 or count%4:raise ValueError('Models must align with all four source cells')
        c=self.c;models=[]
        for i in range(count):
            left=c.mpf(-1)+c.mpf(2)*i/count;right=c.mpf(-1)+c.mpf(2)*(i+1)/count
            center=(left+right)/2;half=(right-left)/2
            point=self._make_view(center);box=self._make_view([ep(left)[0],ep(right)[1]])
            fields={}
            for key in ('m13','m15'):
                p=point['linear_pressure_Taylor'][key];whole=box['linear_pressure_Taylor'][key]
                remainder=c.mpf(abs_upper(whole[2]))*half**2
                ambiguity=(p[0].b-p[0].a)/2+(p[1].b-p[1].a)*half/2
                error=remainder+ambiguity+self.full_errors[key]
                fields[key]=dict(center_intercept=p[0],center_slope=p[1],linear_Z2_half_range=whole[2],
                    linear_Taylor_remainder=remainder,coefficient_ambiguity=ambiguity,
                    full_nonzero_source_error=self.full_errors[key],absolute_error=error,
                    absolute_accuracy_certified=ep(error)[1]<ep(c.mpf('.001'))[0])
            models.append(dict(left=left,right=right,center=center,half_width=half,
                approximation='exact midpoint(intercept endpoints) + exact midpoint(slope endpoints)*(Z-center)',
                source_coefficient_cell=point['source_coefficient_cell'],fields=fields))
        self.assert_graph()
        return models


@source_precision
def run(coefficient_owner=None):
    began=time.monotonic();owner=CurrentOriginalN1EndPressureSource(coefficient_owner,False)
    samples={z:owner.report(owner.evaluate(z)) for z in ('-1','-.639','0','.211','.371','.639','1')}
    models=owner.affine_models()
    raw=dict(source_family=owner.family,samples=record(samples),models=record(models),
        original_physical_units=owner.units,source_bounds=record(owner.bounds),
        source_bindings=owner.source_bindings,
        beta_Wplus=record(owner.Wplus),finite_mu_rho=record(owner.rho),input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_END_PRESSURE_COMPUTED',len(models),flush=True)
    return owner,samples,models
