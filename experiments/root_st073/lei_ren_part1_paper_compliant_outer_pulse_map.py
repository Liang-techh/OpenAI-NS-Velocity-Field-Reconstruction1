"""Actual same-family Section7.31 trial-amplitude axial correction map.

Near-equal rows are replaced by their exact analytic divided difference.
Super-exponentially small coefficients retain a formal logarithmic factor;
their finite factors and true incoming moment functions are enclosed.
The energy-selected ap still requires the later corrected swirl/heat tail.
"""
# Recomputed for the distinct compliant pressure/moment family.
# Formula origin: lei_ren_part1_paper_shared_outer_pulse_map.py; old .01 receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_outer_buffer import SharedOuterBuffer
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def raw_beta(c,x):
    lo,hi=endpoints(x);near=mp.mpf(0) if lo<=0<=hi else min(abs(lo),abs(hi));far=max(abs(lo),abs(hi))
    lower=c.exp(-1/(1-c.mpf(far)**2)) if far<1 else c.mpf(0)
    upper=c.exp(-1/(1-c.mpf(near)**2)) if near<1 else c.mpf(0)
    return c.mpf([endpoints(lower)[0],endpoints(upper)[1]])


def correction_basis(c,mu,normalization,cells=2048):
    """Enclose row1 and (row2-row1)/mu directly from positive integrals."""
    ell=c.mpf('.15');A=[c.mpf(0),c.mpf(0)];D=[c.mpf(0),c.mpf(0)];gram=c.mpf(0)
    for i in range(cells):
        r=c.mpf([-1+mp.mpf(2)*i/cells,-1+mp.mpf(2)*(i+1)/cells]);s=ell*r
        beta=raw_beta(c,r);mass=beta*(c.mpf(2)/cells)/normalization
        for n,d in enumerate((3,1)):
            x=d-s;base=c.exp(-x/2+mu*x)
            # (exp(mu*x)-1)/mu = int_0^x exp(mu*v)dv.
            quotient=c.mpf([endpoints(x)[0],endpoints(x*c.exp(mu*x))[1]])
            A[n]+=mass*base;D[n]+=mass*base*quotient
        gram+=(c.mpf(2)/cells)*beta**2*c.exp(-2*mu*s)/(ell*normalization**2)
    det=A[0]*D[1]-A[1]*D[0]
    if endpoints(det)[1]>=0:raise ArithmeticError('Actual divided-difference matrix not uniformly invertible')
    return dict(first_row=A,exact_divided_difference_row=D,divided_difference_determinant=det,
        energy_gram=gram,end_energy_weights=[c.exp(-26+6*mu)*gram,c.exp(-26+2*mu)*gram],
        beta_radius=ell,cells=cells,
        identity='row2=row1+mu*D; D integrates exp(-x/2+mu*x)*(exp(mu*x)-1)/mu*beta, x=d-s')


def pulse_rows(c,mu,panels=2048,band=48):
    """Common-saddle scaled full pulse integrals, including positive tails.

    P_i=e^(-13 lambda_i/mu) int exp(lambda_i*t) gp(mu*t)dt
       =exp(logpref)*F_i. Use k0=1/(2mu) for both rows. The row
    factor exp(i*(2+u)) is finite and never subtracts lambda1/lambda2.
    """
    k0=1/(2*mu);u0=c.exp(c.ln(4*mu)/3);L=1/u0**2;w=u0/c.sqrt(6)
    if endpoints(u0*(1+band*w))[1]>=mp.mpf('.5'):raise ArithmeticError('Saddle band outside known pulse branch')
    logpref=-2*k0-3*L+2*c.ln(u0)-c.ln(6)/2-c.ln(mu)
    cap=c.exp(-1000);core=[c.mpf(0),c.mpf(0)]
    for i in range(panels):
        x=c.mpf([-band+mp.mpf(2)*band*i/panels,-band+mp.mpf(2)*band*(i+1)/panels])
        v=1+x*w;u=u0*v;phi=x**2*(2*v+1)/(6*v**2);smooth=1/(1-u)**2
        if endpoints(-1/u**2+smooth)[1]>=-1000:raise ArithmeticError('Unresolved flat denominator cap')
        denominator=c.mpf([1,endpoints(1+cap)[1]])
        base=(c.mpf('10.99')-u)*c.exp(smooth-phi)/denominator
        for row in (1,2):core[row-1]+=(c.mpf(2)*band/panels)*base*c.exp(row*(2+u))
    # phi is strictly convex. On each far side its tangent gives a
    # Gaussian tail integral bound, without a spurious large L factor.
    gaussian=c.mpf(0)
    for sign in (-1,1):
        x=c.mpf(sign*band);v=1+x*w
        phi=x**2*(2*v+1)/(6*v**2);derivative=x*(v*v+v+1)/(3*v**3)
        absolute=c.mpf([min(abs(a) for a in endpoints(derivative)),max(abs(a) for a in endpoints(derivative))])
        gaussian+=c.exp(-phi)/absolute
    tails=[];full=[]
    for row in (1,2):
        tail=11*c.exp(4+c.mpf('2.5')*row)*gaussian
        # All remaining support u>=.5, including the startup, is positive.
        late_log=c.ln(121)-k0/2+3*L+13*row-2*c.ln(u0)+c.ln(6)/2
        if endpoints(late_log)[1]>=-1000:raise ArithmeticError('Unscaled earlier pulse support cap failed')
        tail+=cap;tails.append(tail)
        full.append(c.mpf([endpoints(core[row-1])[0],endpoints(core[row-1]+tail)[1]]))
    return dict(common_logpref=logpref,scaled_core_rows=core,scaled_full_rows=full,
        positive_scaled_tail_bounds=tails,saddle_u0=u0,saddle_L=L,panels=panels,band=band,
        flat_denominator_positive_cap=cap,
        proof='phi=L*(2v+v^-2-3)=x^2*(2v+1)/(6v^2), phi_second=1/v^4>0; tangent tails bound the entire u in(0,.5) complement; remaining u>=.5 is separately positive and capped',
        full_pulse_integrals_not_materialized=True)


def pulse_energy(c,cells=4096):
    """Directed entire gp^2 integral with exact exponential cell masses."""
    def value(x):
        if x<=0 or x>=11:return c.mpf(0)
        xx=c.mpf(x)
        if x>=mp.mpf('.02'):primitive=xx-c.mpf('.01')
        else:primitive=c.mpf([0,endpoints(xx*stable_sigma(c,50*xx)[0])[1]])
        return primitive*stable_sigma(c,11-xx)[0]
    total=c.mpf(0)
    for i in range(cells):
        a=mp.mpf(11)*i/cells;b=mp.mpf(11)*(i+1)/cells
        # Primitive is monotone; the cutoff is monotone decreasing.
        pa=c.mpf(a)-c.mpf('.01') if a>=mp.mpf('.02') else c.mpf(0)
        pb=c.mpf(b)-c.mpf('.01') if b>=mp.mpf('.02') else c.mpf(b)
        low=pa*stable_sigma(c,11-c.mpf(b))[0]
        high=pb*stable_sigma(c,11-c.mpf(a))[0]
        gp=c.mpf([max(mp.mpf(0),endpoints(low)[0]),max(mp.mpf(0),endpoints(high)[1])])
        total+=(c.exp(-2*c.mpf(a))-c.exp(-2*c.mpf(b)))*gp**2/2
    if endpoints(total)[0]<=mp.mpf('.24') or endpoints(total)[1]>=mp.mpf('.246'):
        raise ArithmeticError('Directed fixed pulse energy outside paper bounds')
    return total


class SharedOuterPulseMap:
    def __init__(self,basis_cells=2048,pulse_panels=2048):
        self.buffer=SharedOuterBuffer();self.initial=self.buffer.initial;self.ctx=c=self.buffer.ctx
        self.mu=self.buffer.params.mu;self.hashes=dict(self.buffer.hashes)
        for n in ('compliant_outer_buffer','compliant_outer_buffer_check'):
            name=PREFIX+n+'.json';record=json.loads((HERE/name).read_bytes())
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Pulse inlet changed: '+source)
                self.hashes[source]=digest
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            self.basis=correction_basis(c,self.mu,self.initial.repair.normalization,basis_cells)
            self.rows=pulse_rows(c,self.mu,pulse_panels)
            self.Kpulse=pulse_energy(c)
            self.logscale=self.rows['common_logpref']-c.ln(self.mu)
            self.incoming_cap=c.exp(-1000)
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def coefficients(self,Z,a,cells=128):
        """Actual smooth trial-amplitude c_j, represented as exp(logscale)*Cj."""
        c=self.ctx;a=c.mpf(a)
        if endpoints(a)[0]<endpoints(c.mpf('.9'))[0] or endpoints(a)[1]>endpoints(c.mpf('1.2'))[1]:
            raise ValueError('Trial a in [.9,1.2] required')
        with mp.workdps(210):
            inlet=self.buffer.power(Z,1,cells);get=lambda key:IntervalTaylor(c,inlet[key])
            u=get('Utheta_over_Pstar');invP=c.exp(-self.buffer.params.logPstar)
            incoming=[get('Mz_over_R')/u*invP,get('Mtheta_z_over_sqrt2_R_3half_Pstar')/(u*u)*invP]
            logs=[]
            for row,jet in enumerate(incoming,1):
                norm=sum(max(abs(v) for v in endpoints(jet[k])) for k in (0,1))
                if norm<=0:lognorm=None
                else:
                    lognorm=c.ln(c.mpf(norm));logbound=lognorm-13/(2*self.mu)+13*row-self.rows['common_logpref']
                    if endpoints(logbound)[1]>=-1000:raise ArithmeticError('Actual scaled incoming moment cap failed')
                logs.append(lognorm)
            cap=endpoints(self.incoming_cap)[1]
            incoming_scaled=[IntervalTaylor(c,[c.mpf([-cap,cap]),c.mpf([-cap,cap])]) for _ in (1,2)]
            # These exact sources are enclosed, not set to zero. The same
            # source functions satisfy their summed C1 bound; boxes enlarge it.
            rhs=[incoming_scaled[i]+self.rows['scaled_full_rows'][i]*a for i in (0,1)]
            r1=rhs[0]*(-self.mu);r2=-(rhs[1]-rhs[0])
            A=self.basis['first_row'];D=self.basis['exact_divided_difference_row'];det=self.basis['divided_difference_determinant']
            controls=[(r1*D[1]-r2*A[1])/det,(r2*A[0]-r1*D[0])/det]
            if not endpoints(controls[0][0])[1]<0<endpoints(controls[1][0])[0]:
                raise ArithmeticError('Trial end-bump signs unresolved')
            e=get('Mztheta_over_R_Pstar_squared')/(u*u)
            return dict(trial_a=a,common_log_coefficient_scale=self.logscale,
                scaled_coefficients=controls,exact_form='c_j=exp(common_log_coefficient_scale)*scaled_coefficients_j',
                actual_incoming_moment_Taylor_enclosures=incoming,incoming_C1_log_upper_bounds=logs,
                scaled_incoming_C1_positive_cap=self.incoming_cap,
                actual_normalized_incoming_energy=list(e.coefficients),
                actual_source_definition='rows: sum_j Aij cj=-mi exp(-13 lambda_i/mu)-a Pi; mi supplied by F32; Pi full true pulse integrals; lambda_i=.5-i*mu',
                exact_trial_linear_moment_closure=True,
                true_source_coefficients_implicitly_defined=True,
                coefficients_or_incoming_tails_not_replaced_by_zero=True,
                energy_equation='Kpulse*a^2+mu*sum(Kj*cj(a,Z)^2)=(1-exp(-26))/4-mu*incoming_energy+future_corrected_swirl_tail',
                future_corrected_swirl_tail_definition='mu/(2 Rp Utheta(Rp,Z)^2)*integral_Rv^infinity Utheta_corrected^2 dR',
                actual_ap_selected=False,heat_exterior_matched=False,temporal_recursion=False)

    def report(self):
        with mp.workdps(210):
            samples=[self.coefficients('.5',a) for a in ('.9','1','1.2')]
            whole=self.coefficients([-1,1],1)
            return dict(actual_five_defect_family_sha256=self.initial.family,
                implicit_source_sha256=self.initial.datum.source_sha,datum_enclosure_sha256=self.initial.datum.datum_sha,
                mu=self.mu,bump_basis=self.basis,pulse_rows=self.rows,fixed_pulse_energy_Kpulse=self.Kpulse,
                samples=samples,whole_axis_trial_map=whole,
                actual_trial_amplitude_linear_corrections_callable=True,
                exact_smooth_trial_linear_moment_closure_specified=True,
                near_equal_rows_recovered_by_analytic_divided_difference=True,
                actual_ap_selected=False,corrected_swirl_energy_tail_available=False,
                complete_O4_corrected_field_built=False,heat_exterior_matched=False,
                whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    field=SharedOuterPulseMap();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual same-family O.4 trial linear map: exact divided-difference inverse; positive pulse tails retained; ap pending',flush=True)
    return result


if __name__=='__main__':run()
