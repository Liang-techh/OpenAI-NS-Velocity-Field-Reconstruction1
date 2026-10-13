"""Combined original Rp-to-Rv n1 pressure, preserving terminal moments.

Whole-pulse own means use exact terminal integration identities. Nonzero
main weighted source, finite lower endpoints and nonlinear energy terms are
retained; pointwise end-unit bounds are never extended to the main pulse.
"""
import copy
import gzip
import json
import math
import time
from pathlib import Path

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_n1_end_pressure_source as end
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_index_derivative as index

HERE,PREFIX,sha,ep=end.HERE,end.PREFIX,end.sha,end.ep
record,rebox,source_precision=end.record,end.rebox,end.source_precision
NAME=PREFIX+'current_original_n1_pulse_pressure_source.json.gz'
RECEIPT=PREFIX+'current_original_n1_pulse_pressure_source_check.json'
GATE='current_original_n1_combined_pulse_pressure_source_installed'
OPEN=end.OPEN
CELLS,MODELS=end.CELLS,end.MODELS


def pulse_linear_jets(c,Z,delta,mu,ap,C,rows,Pi,Fprime,Wplus,rho,Kmain,Rp_over_Rv):
    z,ubar,_,L=end.polynomials(c,Z,5,delta)
    G2=rows[1]-rows[0]+ap*(Pi[1]-Pi[0])
    rho_sum=sum((j*w for j,w in zip(C,rho)),ap*0)
    I0U=-ubar*rows[0]*mu
    I0m=ubar*(ap*Fprime+rows[0]*13-G2-rho_sum)
    I1U=ubar*(sum((j*w for j,w in zip(C,Wplus)),ap*0)+ap*Kmain)
    I1m=-I1U-ubar*rows[0]*(mu*Rp_over_Rv)
    q0=end.radial_q(c,Z,delta,I0U,I0m)
    q1=end.radial_q(c,Z,delta,I1U,I1m)
    qinit=end.radial_q(c,Z,delta,ap*0,ubar*rows[0]*mu)
    L3=L.truncate(3);z3=z.truncate(3)
    S0=-qinit.truncate(3)/L3+end.derivative(q0)*z3*((1-delta)/2)/L3
    S1=(-qinit.truncate(3)*Rp_over_Rv-q1.truncate(3))/L3
    S1+=end.derivative(q1)*z3*((1-delta)/2)/L3
    return dict(m13_linear=-S0/2,m15_linear=S1/2,qinit=qinit,
        I0U=I0U,I0m=I0m,I1U=I1U,I1m=I1m,I0q=q0,I1q=q1,
        index_derivative_source=Fprime,finite_mu_X_remainder=rho_sum,
        main_weighted_source=ubar*ap*Kmain,
        lower_own_mean=ubar*rows[0]*mu,
        lower_q_plus_R_qR_exact_zero=True,m15_viscous_coefficient=-qinit)


def pulse_energy_caps(c,mu,logA0,logE,reference,B,Tinitial):
    """Normalize actual whole-pulse energy estimates in their exact units."""
    length=13/mu;logs=[]
    for b,t in zip(B,Tinitial):
        lb=c.ln(b);lt=c.ln(t)
        common=logA0-logE-13/(2*mu)+13-c.ln(reference)
        row=dict(unweighted_U=common+2*lb-c.ln(2*mu),
            unweighted_initial=common+2*lt,
            weighted_U_with_length_squared=common-26+c.ln(2)+2*lb+c.ln(1+length**2),
            weighted_initial_with_length=logA0-logE-39/(2*mu)+13+2*lt+c.ln(1+length)-c.ln(reference))
        if any(ep(value)[1]>=-1000 for value in row.values()):
            raise ValueError('Original whole-pulse physical energy cap not proved')
        logs.append(row)
    epsilon=c.exp(-1000)
    # E(m_j)<=10 epsilon in either fixed unit, j=0,1,2.
    # q<=4U+2m+2mZ; qZ<=16U+4UZ+8m+8mZ+2mZZ.
    q_energy=288*epsilon;qZ_energy=7960*epsilon
    cross=c.sqrt(7960)*epsilon
    S0=(8+4+8)*epsilon+2*cross
    S1=(8+144+8+8)*epsilon+2*cross
    return dict(log_caps=logs,length=length,epsilon=epsilon,
        U_derivative_bounds=B,initial_mean_derivative_bounds=Tinitial,
        source_energy_bounds_nonzero=True,
        unweighted_Hardy_mean_factor=10,weighted_Cauchy_mean_factor=3,
        conservative_common_mean_factor=10,q_energy=q_energy,qZ_energy=qZ_energy,
        boundary_q_squared_bound=16*epsilon,
        m13=end.scalar_bound(c,S0/2),m15=end.scalar_bound(c,S1/2),
        exact_units_combined_before_boxing=True)


class CurrentOriginalN1PulsePressureSource(end.CurrentOriginalN1EndPressureSource):
    @source_precision
    def __init__(self,end_owner=None,index_owner=None,require_checked=True):
        self.end=end_owner or end.CurrentOriginalN1EndPressureSource()
        self.index=index_owner or index.CurrentOriginalRpPulseSaddleIndexDerivative()
        if not (self.end.acceptance_loaded and self.index.acceptance_loaded):
            raise ValueError('Checked original end pressure and index derivative required')
        self.end.assert_graph()
        self.coefficients=self.end.coefficients;self.c=c=self.end.c
        self.family=copy.deepcopy(self.end.family)
        if self.index.family!=self.family:raise ValueError('Same original pulse family required')
        for name in ('mu','delta','Pi','Wplus','rho','reference','invRv'):
            setattr(self,name,getattr(self.end,name))
        self.Fprime=rebox(c,self.index.index_derivative)
        logE=rebox(c,self.coefficients.selected.amplitude.log_end_scale)
        if ep(logE)[1]>=0:raise ValueError('Original positive end scale must be below one')
        self.Kmain_log=c.ln(11/(c.mpf('1.5')-self.mu))-logE-3/self.mu+2
        self.Rp_over_Rv_log=-13/self.mu
        if ep(self.Kmain_log)[1]>=-1000 or ep(self.Rp_over_Rv_log)[1]>=-1000:
            raise ValueError('Actual weighted main/lower source caps not proved')
        self.Kmain=end.scalar_bound(c,c.exp(-1000))
        self.Rp_over_Rv=end.scalar_bound(c,c.exp(-1000))
        self._cells=copy.copy(self.end._cells);self._packets={}
        # Actual source is U=A0*exp(-(.5+mu)t)*ubar*B(t,Z).
        # gp lies in[0,11]; end beta is bounded by its original normalization.
        _,ubar,_,_=end.polynomials(c,['-1','1'],5,self.delta)
        beta=c.exp(-1)/(c.mpf('.15')*self.coefficients.basis.normalization)
        B=[c.mpf(0)]*3;Tinitial=[c.mpf(0)]*3
        for cell in self._cells:
            p=self.coefficients.evaluate(list(cell['ends']));view=self.coefficients.report(p)
            first=end.read_jet(c,view['source_inputs']['incoming_moment_Taylor'][0])
            initial=ubar*first;main=ubar*cell['ap']*11
            ends=[ubar*j*beta for j in cell['C']]
            for j in range(3):
                b=(end.abs_upper(main[j])+sum(end.abs_upper(v[j]) for v in ends))*math.factorial(j)
                t=end.abs_upper(initial[j])*math.factorial(j)
                B[j]=c.mpf(max(ep(B[j])[1],b))
                Tinitial[j]=c.mpf(max(ep(Tinitial[j])[1],t))
        self.energy=pulse_energy_caps(c,self.mu,self.end.logP+self.end.logUp,logE,self.reference,B,Tinitial)
        self.errors={key:self.energy[key] for key in ('m13','m15')}
        qmax=max(end.abs_upper(self._terms(list(cell['ends']),cell)['qinit'][0]) for cell in self._cells)
        self.viscous_error=end.scalar_bound(c,self.invRv*c.mpf(qmax)/self.reference)
        self.full_errors={key:self.errors[key]+(self.viscous_error if key=='m15' else 0) for key in self.errors}
        self.units=copy.deepcopy(self.end.units)
        self.units.update(radial_interval='Rp<=R<=Rv',all_original_pulse_charts_combined=True,
            main_weighted_kernel='exp(-logE-39/(2mu)+13)*integral_0^(11/mu) exp((1.5-mu)t)*gp(mu*t)dt',
            initial_Rp_over_Rv='exp(-13/mu)',terminal_M_at_Rv_exact_zero=True)
        self.bounds=dict(whole_pulse_energy=self.energy,weighted_main_source_log=self.Kmain_log,
            positive_weighted_main_source=self.Kmain,lower_radius_ratio_log=self.Rp_over_Rv_log,
            positive_lower_radius_ratio=self.Rp_over_Rv,viscous_error=self.viscous_error,
            reference=self.reference,qinit_bound=c.mpf(qmax))
        self.source_bindings={name:end.coefficients.function_binding(fn) for name,fn in (
            ('whole_pulse_linear',pulse_linear_jets),('whole_pulse_energy',pulse_energy_caps),
            ('original_main_profile',end.original_end.CompliantAxialPulseField.main))}
        self.hashes=dict(self.end.hashes)
        for name,digest in self.index.hashes.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting source '+name)
            self.hashes[name]=digest
        for module in (end,index):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Checked combined pressure receipt required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed combined pulse input '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._state=record((self.family,self.units,self.bounds,self._cells,self.Fprime))
        self._hashes=copy.deepcopy(self.hashes)

    def _terms(self,Z,cell):
        return pulse_linear_jets(self.c,Z,self.delta,self.mu,cell['ap'],cell['C'],cell['rows'],
            self.Pi,self.Fprime,self.Wplus,self.rho,self.Kmain,self.Rp_over_Rv)

    @source_precision
    def assert_graph(self):
        self.end.assert_graph()
        if self._state!=record((self.family,self.units,self.bounds,self._cells,self.Fprime)) or self.hashes!=self._hashes:
            raise ValueError('Combined original pulse source state changed')

    def _make_view(self,Z):
        c=self.c;cell=self._cell(Z);terms=self._terms(Z,cell)
        linear={key:terms[key+'_linear']/self.reference for key in self.errors}
        full={key:linear[key][0]+c.mpf([-ep(self.errors[key])[1],ep(self.errors[key])[1]]) for key in self.errors}
        viscous=-terms['qinit'][0]*self.invRv/self.reference;full['m15']+=viscous
        return dict(Z=c.mpf(Z),source_family=self.family,source_coefficient_cell=list(cell['ends']),
            full_pressure_values=full,linear_pressure_Taylor=linear,source_linear_inputs=terms,
            original_physical_units=self.units,quadratic_source_errors=self.errors,
            weighted_viscous_endpoint=viscous,unweighted_viscous_boundary_exact_zero=True,
            all_original_pulse_charts_combined=True,terminal_own_mean_identities_retained=True,
            positive_weighted_main_source_retained=True,full_pressure_derivative_accuracy_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))


@source_precision
def run(end_owner=None,index_owner=None):
    began=time.monotonic();owner=CurrentOriginalN1PulsePressureSource(end_owner,index_owner,False)
    samples={z:owner.report(owner.evaluate(z)) for z in ('-1','-.639','0','.211','.371','.639','1')}
    models=owner.affine_models()
    raw=dict(source_family=owner.family,samples=record(samples),models=record(models),
        original_physical_units=owner.units,source_bounds=record(owner.bounds),
        source_bindings=owner.source_bindings,input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('COMBINED_ORIGINAL_PULSE_PRESSURE_COMPUTED',len(models),flush=True)
    return owner,samples,models
