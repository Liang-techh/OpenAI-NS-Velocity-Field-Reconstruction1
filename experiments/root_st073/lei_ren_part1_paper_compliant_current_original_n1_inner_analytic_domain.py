"""Same-source core/actual-first-collar bounds for the genuine n1 solve.

The analytic core extension past Ra is a comparison field. The actual
prescribed-shear field is used in the collar. Inverse microscopic width
is kept formal, then cancelled by the real integration measure. These
are bounds and a Neumann-series admission, not computed solution terms.
"""
import ast
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_regular_system as parent
import lei_ren_part1_paper_compliant_actual_bridge_integrals as bridge
import lei_ren_part1_paper_compliant_current_core_common_fixed_point as common
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=parent.HERE,parent.PREFIX,parent.sha
NAME=PREFIX+'current_original_n1_inner_analytic_domain.json.gz'
RECEIPT=PREFIX+'current_original_n1_inner_analytic_domain_check.json'
GATE='current_original_n1_core_and_actual_first_collar_analytic_bounds_installed'
OPEN=parent.OPEN+('completed_leading_supports_outside_Rin_certified',
                 'n1_Volterra_iteration_prefix_computed')
ep=parent.parent.ends


def upper(c,value):return c.mpf(ep(c.mpf(value))[1])
def lower(c,value):return c.mpf(ep(c.mpf(value))[0])


class WidthBound:
    """Positive Laurent majorant sum c_p hb^p of the ORIGINAL width.

    Numeric enclosure is permitted only after all negative powers cancel.
    No representative width, midpoint or capped field is selected.
    """
    def __init__(self,c,rows):
        self.c=c
        self.rows={int(p):upper(c,v) for p,v in rows.items() if ep(v)[1]!=0}
        if any(ep(v)[0]<0 for v in self.rows.values()):raise ValueError('Nonnegative bounds required')

    @classmethod
    def scalar(cls,c,value):return cls(c,{0:c.mpf(value)})

    def _cast(self,value):
        if isinstance(value,WidthBound):
            if value.c is not self.c:raise ValueError('Same source interval context required')
            return value
        return WidthBound.scalar(self.c,value)

    def __add__(self,value):
        value=self._cast(value);zero=self.c.mpf(0)
        return WidthBound(self.c,{p:self.rows.get(p,zero)+value.rows.get(p,zero)
                                 for p in self.rows.keys()|value.rows.keys()})
    __radd__=__add__

    def __mul__(self,value):
        value=self._cast(value);rows={};zero=self.c.mpf(0)
        for p,a in self.rows.items():
            for q,b in value.rows.items():rows[p+q]=rows.get(p+q,zero)+a*b
        return WidthBound(self.c,rows)
    __rmul__=__mul__

    def __truediv__(self,value):
        if isinstance(value,WidthBound):raise ValueError('Formal width bounds may not be inverted')
        return self*(self.c.mpf(1)/value)

    def shift(self,power):return WidthBound(self.c,{p+power:v for p,v in self.rows.items()})

    def numeric_upper(self,cap):
        if any(p<0 for p in self.rows):raise ValueError('Uncancelled inverse source width')
        return upper(self.c,sum((v*cap**p for p,v in self.rows.items()),self.c.mpf(0)))

    def log_guarded_upper(self,logh,fallback_log):
        """Enclose positive width products only after an exact log proof.

        A fixed 1e-50000 width cap can swamp the cancelled source. Instead,
        bound each remaining positive-width term below exp(-100) times
        the width-free term. This changes only an outward error bound.
        """
        if any(p<0 for p in self.rows):raise ValueError('Uncancelled inverse source width')
        zero=self.rows.get(0,self.c.mpf(0))
        target=self.c.ln(zero)-100 if ep(zero)[0]>0 else self.c.mpf(fallback_log)
        total=zero;proofs={}
        for p,value in self.rows.items():
            if p==0:continue
            source_log=p*logh+self.c.ln(value)
            if ep(source_log)[1]>ep(target)[0]:raise ValueError('Positive width product needs a sharper log bound')
            cap=upper(self.c,self.c.exp(target));total+=cap
            proofs[str(p)]=dict(exact_positive_source_log_upper=source_log,
                               numerical_product_log_cap=target,product_cap=cap,
                               cap_is_error_upper_only=True,passed=True)
        return upper(self.c,total),proofs

    def export(self,logh):
        return {str(p):dict(coefficient=v,log_width_scale=p*logh,
                           exact_width_power=p) for p,v in self.rows.items()}


def unsigned_stirling1(order):
    rows=[[1]]
    for k in range(1,order+1):
        previous=rows[-1]
        rows.append([(previous[j-1] if j else 0)+(k-1)*(previous[j] if j<len(previous) else 0)
                     for j in range(k+1)])
    return rows


def source_bindings():
    controls=bridge.bridge_control_bindings()
    pressure_kernel=common.raw_flatten_pressure_function()
    interfaces=parent.parent.first.defining_boundary_bindings()
    # The flatten power comes from the actual outer definition, rather than
    # an older candidate helper's explanatory kernel description.
    tree=ast.parse((HERE/'lei_ren_part1_paper_outer.py').read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef)
            and n.name=='_evaluate_log_radius')
    wanted='log_a-Decimal(str(z_factor_log))+Decimal(str(flat_sigma*(z_factor_log-math.log(2.0))))'
    assignments=[n for n in ast.walk(fn) if isinstance(n,ast.Assign)
                 and any(ast.unparse(t)=='log_u' for t in n.targets)]
    if sum(ast.dump(n.value)==ast.dump(ast.parse(wanted,mode='eval').body) for n in assignments)!=1:
        raise ValueError('Actual flatten angular source changed')
    return dict(original_bridge_controls=controls,
        original_core_actual_first_interface_bindings=interfaces,
        original_all_order_P0_source_bindings=common.all_order_pressure_source_bindings(),
        exact_raw_flatten_function={key:str(value) for key,value in pressure_kernel.items()},
        actual_flatten_log_angular_AST_bound=True,
        squared_flatten_kernel='2**(-2*sigma)*(1+Z**2)**(-2*(1-sigma)); 0<=sigma<=1',
        analytic_log_branch='Re(1+Z**2)>=1-|Im Z|**2>0',
        comparison_core_interval='0<=y<=hb (alpha=1)',
        actual_core_collar_is_not_analytic_core_extension=True,
        exact_actual_collar_identities=dict(
            Phi_rho='chi*Phi*A',
            Phi_rhorho='Phi*(chi_rho*A+chi**2*A**2+chi*A_rho)',
            Uz_rho='chi*(Phi/Phi_core)*Uz_core_rho',
            Uz_rhorho='(Phi/Phi_core)*(chi_rho*Uz_core_rho+chi*(chi-1)*A*Uz_core_rho+chi*Uz_core_rhorho)',
            A='Phi_core_rho/Phi_core',chi='1-(1-hb)*sigma(y/hb)',
            chi_rho='-(1-hb)*sigma_prime(y/hb)/(hb*rho)'),
        identities_use_exact_reduced_core_equations=True)


def collar_radial_bounds(c,phi,uz,phi_floor,phi_value,uz_value,order=4):
    """Differentiate the defining ODE; keep every inverse hb power."""
    # rho>=4, |chi|<=1, |1-chi|<=1. Fixed real radial cutoff has no Z poles.
    sigma=sigma_jets(c,c.mpf([0,'.5']))
    sig=[upper(c,abs(v)*math.factorial(k)) for k,v in enumerate(sigma.coefficients)]
    stirling=unsigned_stirling1(order)
    scalar=lambda v:WidthBound.scalar(c,v)
    chi=[scalar(1)]+[WidthBound(c,{-j:c.mpf(stirling[k][j])*sig[j]/4**k
                                  for j in range(1,k+1)}) for k in range(1,order+1)]
    # phi_core*A=phi_core_rho determines derivatives of A.
    A=[]
    for k in range(order):
        A.append(upper(c,(phi[k+1]+sum((math.comb(k,j)*phi[j]*A[k-j]
                                       for j in range(1,k+1)),c.mpf(0)))/phi_floor))
    C=[sum((chi[j]*math.comb(k,j)*A[k-j] for j in range(k+1)),scalar(0))
       for k in range(order)]
    # ratio_rho=(chi-1)*A*ratio has the same safe positive majorants as C.
    Phi=[scalar(phi_value)];ratio=[scalar(phi_value/phi_floor)]
    for k in range(order):
        Phi.append(sum((C[j]*Phi[k-j]*math.comb(k,j) for j in range(k+1)),scalar(0)))
        ratio.append(sum((C[j]*ratio[k-j]*math.comb(k,j) for j in range(k+1)),scalar(0)))
    W=[scalar(uz_value)]
    for k in range(1,order+1):
        W.append(sum((chi[i]*ratio[j]*uz[k-i-j]*math.factorial(k-1)
                      / (math.factorial(i)*math.factorial(j)*math.factorial(k-1-i-j))
                      for i in range(k) for j in range(k-i)),scalar(0)))
    M=[scalar(uz_value)]
    for k in range(1,order+1):M.append((W[k-1]+k*M[k-1])/4)
    return dict(Phi=Phi,Uz=W,mean=M,chi=chi,A=A,sigma_derivatives=sig)


def hierarchy_bounds(c,phi,w,mean,Fmax,delta,eta,rho_max,Lmin):
    """Absolute majorants for actual F, Uz, P increment and Q, then n1."""
    scalar=lambda value:WidthBound.scalar(c,value)
    z=1+eta;d=1+z*z;dz_radius=eta/2
    f=[v*Fmax for v in phi]
    # Mean derivatives are bounded on the outer tube. A Cauchy derivative
    # in Q is taken on the middle tube; all later jets use the inner tube.
    q=[(w[i]*(2*z)+mean[i]*((1-delta)*z)+mean[i]*(d/dz_radius))/Lmin for i in range(5)]
    jet_radius=eta/4
    at=lambda rows,i,k:rows[i]*(math.factorial(k)/jet_radius**k)
    def second(rows,weight):
        f0,fz,fr=at(rows,0,0),at(rows,0,1),at(rows,1,0)
        A=f0*(abs(weight)*z)+fz*d+fr*(2*z*rho_max)
        Az=f0*abs(weight)+fz*(abs(weight-2)*z)+at(rows,0,2)*d+fr*(2*rho_max)+at(rows,1,1)*(2*z*rho_max)
        Ar=fr*(abs(weight-2)*z)+at(rows,1,1)*d+at(rows,2,0)*(2*z*rho_max)
        return (A*(abs(weight-1+delta)*z/Lmin)
                +(Az/Lmin+A*(2*delta*z/(Lmin*Lmin)))*d+Ar*(2*z*rho_max/Lmin))/Lmin
    nt=second(f,-2-delta);nz=second(w,-1-delta)
    q0,qz,qr,qrr=at(q,0,0),at(q,0,1),at(q,1,0),at(q,2,0)
    # V0=R Q has no inverse-R term; Lambda is supplied by the owner below.
    np_inertial=((q0+qz*((1-delta)*z/2)+qr*rho_max)/Lmin
                 +q0*(q0/2+qr*rho_max)+w[0]*(q0*(2*z)+qz*d+qr*(2*z*rho_max))/Lmin)/2
    np_viscous=(qr*4+qrr*(2*rho_max))/2
    return dict(F=f,Uz=w,mean=mean,Q=q,N1theta=nt,N1z=nz,
                N1p_inertial=np_inertial,N1p_viscous_per_Lambda=np_viscous,
                Z_jet_radius=jet_radius)


def regular_matrix_norms(c,b,epsilon,delta,eta,rho_max,Lmin):
    """Infinity row norms on the common middle tube; no width inverse."""
    z=1+eta;d=1+z*z;x=c.sqrt(rho_max);rZ=eta/2
    value=lambda rows,i:rows[i].numeric_upper(c.mpf('1e-50000'))
    f,fr,w,wr,q=(value(b[key],i) for key,i in (('F',0),('F',1),('Uz',0),('Uz',1),('Q',0)))
    fz,wz=f/rZ,w/rZ
    Cf=f+rho_max*fr;Zf=((2+delta)*z*f+d*fz+2*z*rho_max*fr)/Lmin
    Zw=((1+delta)*z*w+d*wz+2*z*rho_max*wr)/Lmin
    H=(1-delta)*z/2+d*w;B=q+(1+2*z*w)/Lmin;E=z*w+c.mpf('.5')
    row0=[c.mpf(1),c.mpf(1),c.mpf(1),4*epsilon*x*f,
        epsilon*(2*(q+(2-delta)*E/Lmin)+2*Zf+2*(1-delta)*z*Cf/Lmin
                 +2*(1+delta)*z*Cf/Lmin+x*B),
        epsilon*(8*epsilon*z*rho_max*f/Lmin+2*((1-delta)*E/Lmin+Zw+(1-delta)*z*rho_max*wr/Lmin)
                 +2*(1+delta)*z*rho_max*wr/Lmin+4*z/Lmin+x*B)]
    row1=[epsilon*(2*H/Lmin+4*d*Cf/Lmin),
          epsilon*(2*(H+d*rho_max*wr)/Lmin+2*d*rho_max*wr/Lmin+2*d/Lmin)]
    return upper(c,max((ep(v)[1] for v in row0))),upper(c,max(ep(v)[1] for v in row1))


def neumann_tail(c,M0,a,C0,C1,loss,N):
    """Tail AFTER terms 0..N; initial Gg uses a weighted L1 source norm.

    At most ceil(k/2) derivative blocks occur in a nonzero word. The
    initial source is integrated first: k further simplex factors give k!.
    """
    if not isinstance(N,int) or isinstance(N,bool) or N<1:raise ValueError('N>=1 required')
    alpha=C1/(C0*loss)
    def term(k):
        m=(k+1)//2
        factor=max(c.mpf(1),upper(c,alpha*m),key=lambda v:ep(v)[1])
        return upper(c,M0*(2*a*C0)**k*factor**m/math.factorial(k))
    q=upper(c,4*c.e*a*a*C0*C0*(1+alpha*(N+4)/2)/((N+2)*(N+3)))
    if ep(q)[1]>=1:raise ValueError('Two-step ratio does not certify this tail')
    return dict(retained_iteration_indices=[0,N],first_omitted_term_bound=term(N+1),
        second_omitted_term_bound=term(N+2),two_step_tail_ratio_upper=q,
        omitted_sum_bound=upper(c,(term(N+1)+term(N+2))/(1-q)),
        parameter_derivative_loss=loss,alpha=alpha,
        formula='M0*(2*a*C0)^k/k!*max(1,alpha*ceil(k/2))^ceil(k/2)',
        source_measure_is_integrated_before_simplex_bound=True,
        actual_iteration_terms_computed=False)


@dataclass(frozen=True,eq=False)
class OriginalN1InnerDomainPacket:
    source_kind:str='core_and_actual_first_collar'


class CurrentOriginalN1InnerAnalyticDomain:
    @source_precision
    def __init__(self,n1,require_checked=True):
        if type(n1) is not parent.CurrentOriginalN1RegularSystem or not n1.acceptance_loaded:
            raise ValueError('Accepted genuine n1 owner required')
        n1.assert_graph();self.n1=n1;self.c=c=n1.c;self.family=copy.deepcopy(n1.family)
        self.source=n1.source.source;self.core=n1.source.core
        self.hashes=dict(n1.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name,
                     'lei_ren_part1_paper_outer.py',PREFIX+'flat_pulse_derivatives.py',
                     PREFIX+'current_core_common_fixed_point.py',PREFIX+'K1_ledger.json'):
            self.hashes[name]=sha(name)
        self.bindings=source_bindings();self.acceptance_loaded=False;self._packets={}
        s=self.source;major=s.records['core_transfer'];tube=s.records['shared_analytic_tube'];norm=s.records['physical_norm_family']
        analytic=major['analytic_core_family_sha256']
        if not (major['contraction_proved'] and norm['uniform_analytic_fixed_point_admission_extended']
                and tube['common_complex_axis_poles_excluded']
                and norm['base_analytic_core_family_sha256']==analytic==tube['analytic_core_family_sha256']):
            raise ValueError('Same admitted analytic family required')
        ledger=json.loads((HERE/(PREFIX+'K1_ledger.json')).read_bytes())
        if not (ledger['h_b_equals_epsilon_b_by_definition'] and ledger['positive_width_not_materialized']
                and ledger['uniform_Cstar_family_sha256']==norm['uniform_Cstar_family_sha256']
                and ledger['base_analytic_core_family_sha256']==analytic):
            raise ValueError('Same exact hb=epsilon_b definition required')
        if ep(s.logh)[1]>=ep(c.ln(s.cap))[0]:raise ValueError('Original positive width cap required')
        # The live pressure owner hydrates logPstar only. Inspect the checked
        # original parameter enclosures instead of inventing missing fields.
        ps=s.records['pressure_source']['compliant_source'];p=ps['parameter_bounds']
        read=parent.parent.original.read
        if not (ep(read(c,p['log_mu']))[1]<0 and ep(read(c,p['Tw']))[0]>=0
                and ep(read(c,p['Ts']))[0]>0 and not p['coordinates']['d']['exponential_terms']
                and ep(read(c,p['coordinates']['d']['finite_offset']))[0]>=1
                and ps['stage_count_total']==len(ps['stages'])==14
                and all(ep(read(c,v['mass']))[0]>=0 for v in ps['stages'].values())):
            raise ValueError('Original nonnegative pressure measure and flatten branch premises required')
        self._owner=n1;self._context=c;self._family=copy.deepcopy(self.family)
        self._source=s;self._width=(s.logh._mpi_,s.cap._mpi_,s.logRa._mpi_)
        self._norm=(s.rebuild.phi_norm._mpi_,s.rebuild.psi_norm._mpi_)
        self._records={key:copy.deepcopy(s.records[key]) for key in
                       ('core_transfer','shared_analytic_tube','physical_norm_family','pressure_source')}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Inner domain receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Inner analytic input changed '+name)
                parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.n1 is not self._owner or self.c is not self._context or self.family!=self._family
                or self.n1.family!=self.family or self.source is not self._source
                or self.source is not self.n1.source.source or not self.n1.acceptance_loaded
                or self._width!=(self.source.logh._mpi_,self.source.cap._mpi_,self.source.logRa._mpi_)
                or self._norm!=(self.source.rebuild.phi_norm._mpi_,self.source.rebuild.psi_norm._mpi_)
                or any(self.source.records[key]!=value for key,value in self._records.items())
                or self.source.actual.core is not self.source.core
                or self.source.actual.comparison is not self.source.comparison
                or self.source.actual.logh._mpi_!=self.source.logh._mpi_
                or self.source.actual.bridge.logh._mpi_!=self.source.logh._mpi_
                or self.hashes!=self._hashes):raise ValueError('Inner analytic source graph changed')
        return self.n1.assert_graph()

    @source_precision
    def domain(self):
        self.assert_graph();c=self.c;s=self.source;core=self.core
        Bphi=upper(c,s.rebuild.phi_norm);Bpsi=upper(c,s.rebuild.psi_norm)
        h=lower(c,core.h);floor=lower(c,core.phi_floor)
        eta=lower(c,h*floor/(128*(1+Bphi)))
        joint_sum=c.mpf(5)/20+c.mpf('.5')
        variation=upper(c,16*Bphi*eta/h)
        complex_floor=lower(c,floor-variation)
        z=1+eta;Lmin=lower(c,1-core.delta*z*z);rho_max=c.mpf('4.1')
        if not (ep(eta)[0]>0 and ep(eta)[1]<ep(h/8)[0] and ep(joint_sum)[1]<1
                and ep(complex_floor)[0]>=ep(7*floor/8)[0] and ep(Lmin)[0]>0):
            raise ValueError('Same-source common nozero tube failed')
        # Cauchy on |rho|<5, Z disks h/2; rho radius .5 and Z radius h/4.
        core_phi=[upper(c,4*Bphi*math.factorial(i)*2**i) for i in range(5)]
        core_w=[upper(c,(4*z+core.j if i==0 else 0)+core.epsilon*4*Bpsi*math.factorial(i)*2**i)
                for i in range(5)]
        A0=core_phi[1]/complex_floor
        integral=upper(c,s.cap*rho_max*A0/2)
        if ep(4*c.exp(s.cap/2))[1]>=ep(rho_max)[0]:raise ValueError('Actual collar must stay inside comparison analytic rectangle')
        if ep(integral)[1]>=ep(c.mpf('1/16'))[0]:raise ValueError('Actual collar exponential bound unresolved')
        phi_b=upper(c,core_phi[0]*c.exp(integral))
        ratio_b=phi_b/complex_floor
        w_b=upper(c,core_w[0]+s.cap*rho_max*ratio_b*core_w[1]/2)
        collar=collar_radial_bounds(c,core_phi,core_w,complex_floor,phi_b,w_b)
        scalar=lambda value:WidthBound.scalar(c,value)
        core_rows=dict(Phi=[scalar(v) for v in core_phi],Uz=[scalar(v) for v in core_w],
                       mean=[scalar(v) for v in core_w])
        Fmax=upper(c,c.exp(-2*core.logLambda-1000))
        if ep(core.logC-(core.Lambda*core.Gbar+2*core.logLambda+1000))[0]<=0:
            raise ValueError('Actual complex amplitude guard required')
        bound_args=(Fmax,core.delta,eta,rho_max,Lmin)
        cb=hierarchy_bounds(c,core_rows['Phi'],core_rows['Uz'],core_rows['mean'],*bound_args)
        bb=hierarchy_bounds(c,collar['Phi'],collar['Uz'],collar['mean'],*bound_args)
        def forcing(b):
            np=b['N1p_inertial']+b['N1p_viscous_per_Lambda']*core.Lambda
            return {'4':np*(2*core.epsilon*c.sqrt(rho_max)),
                    '5':b['N1theta']*(2*core.epsilon),
                    '6':b['N1z']*(2*core.epsilon)+np*(4*core.epsilon**2*z*rho_max/Lmin)}
        gc,gb=forcing(cb),forcing(bb)
        # dx=(sqrt(rho)/2)*hb ds, and 0<=s<=1/2. All inverse hb
        # powers in g are at most one. Integration cancels them exactly.
        integrated={key:row.shift(1)*(c.sqrt(rho_max)/4) for key,row in gb.items()}
        guarded={key:row.log_guarded_upper(s.logh,-10*core.logLambda-1000) for key,row in integrated.items()}
        integrated_upper={key:row[0] for key,row in guarded.items()}
        core_upper={key:row.numeric_upper(s.cap)*2 for key,row in gc.items()}
        M0=upper(c,max(ep(core_upper[key]+integrated_upper[key])[1] for key in gc))
        # Sum, rather than select, core and bridge bounds for a uniform
        # matrix majorant. No inverse width appears in the matrix itself.
        together={key:[a+b for a,b in zip(cb[key],bb[key])] for key in ('F','Uz','Q')}
        C0,C1=regular_matrix_norms(c,together,core.epsilon,core.delta,eta,rho_max,Lmin)
        tail=neumann_tail(c,M0,c.sqrt(rho_max),C0,C1,eta/8,64)
        datum=s.core.datum;pressure=upper(c,(datum.m2+datum.stages['z_flatten']['mass'])/(1-eta*eta)**2+datum.m0)
        phi_all=[a+b for a,b in zip(core_rows['Phi'],collar['Phi'])]
        pressure_increment=[phi_all[0]*phi_all[0]*rho_max]+[
            sum((phi_all[j]*phi_all[i-1-j]*math.comb(i-1,j) for j in range(i)),scalar(0))
            for i in range(1,5)]
        def export(rows):return {key:[v.export(s.logh) for v in value] for key,value in rows.items()}
        raw=dict(source_family=self.family,source_bindings=self.bindings,
            joint_core_domain=dict(rho_complex_modulus_less_than='5',Z_real_centers=['-1','1'],
                Z_disk_radius=h/2,joint_normal_convergence_ratio=joint_sum,
                Xh_norm_majorants=dict(Phi=Bphi,Psi=Bpsi),core_value_factor=4,
                Cauchy_rho_radius='.5',Cauchy_Z_radius=h/4),
            common_tubes=dict(outer_Z_radius=eta,matrix_Z_radius=eta/2,Gg_Z_radius=eta/4,solution_Z_radius=eta/8),
            denominator_bounds=dict(Phi_core_modulus_lower=complex_floor,Phi_variation_upper=variation,L_modulus_lower=Lmin),
            original_pressure=dict(normalized_modulus_upper=pressure,physical_log_scale=2*core.logP,
                exact_function=s.core.datum.normalized_jets.__name__,all_fourteen_atoms_retained=True,
                inherited_exact_P0_binding=self.n1.source.P0_binding),
            original_flatten_parameter_premises=dict(zero_lt_mu_lt_one=True,Tw_nonnegative=True,Ts_positive=True,yd_ge_one=True),
            original_axis_amplitude=dict(complex_modulus_upper=Fmax,source_formula='exp(-selected_logCstar-Lambda*G(Z))',
                no_source_amplitude_materialized=True),
            actual_inner_geometry=dict(core='0<=rho<=4',collar='rho=4*exp(hb*s), 0<=s<=1/2',
                Rin='Ra*exp(hb/2)',Rkeep='Ra*exp(hb/4)',source_log_hb=s.logh,
                exact_positive_width_preserved=True,numerical_cap_bounds_only=s.cap,
                real_rho_upper_bound=rho_max,strict_Ra_Rkeep_Rin_order_from_positive_hb=True,
                actual_bridge_not_replaced_by_comparison=True),
            actual_collar_exponential_integral_upper=integral,
            actual_collar_radial_bounds=export({key:collar[key] for key in ('Phi','Uz','mean')}),
            actual_core_radial_bounds=export(core_rows),
            whole_inner_pressure_increment_rho_derivative_bounds=[v.export(s.logh) for v in pressure_increment],
            pressure_increment_prefactor='epsilon_core*F0(Z)**2; complex upper uses epsilon_core*Fmax**2',
            physical_R_derivative_rule='multiply rho derivative order i by exp(i*logLambda); retain formal hb powers',
            ordinary_Z_derivative_rule='for fields on the Gg tube use k!/(eta/4)^k times corresponding middle-tube bound; solution tube reserves another eta/8',
            actual_collar_regular_forcing_bounds={key:value.export(s.logh) for key,value in gb.items()},
            actual_collar_source_integral_bounds={key:value.export(s.logh) for key,value in integrated.items()},
            source_measure_identity='dx=(sqrt(rho)/2)*hb*ds; first-collar phase length 1/2',
            actual_collar_integrated_source_upper=integrated_upper,core_Gg_upper=core_upper,
            actual_collar_integrated_source_log_product_guards={key:row[1] for key,row in guarded.items()},
            initial_Gg_common_tube_bound=M0,regular_B0_infinity_norm_upper=C0,regular_B1_infinity_norm_upper=C1,
            nonadjacent_B1_support_inherited=True,regular_neumann_tail_after_degree64=tail,
            source_inner_common_holomorphic_bounds_constructed=True,
            all_fixed_radial_derivatives_holomorphic_on_same_domain_by_smooth_ODE=True,
            explicit_radial_bound_orders=[0,4],source_radial_moment_recovery_preserved=True,
            actual_assembled_outer_supports_not_inferred_from_core_or_collar=True,
            no_actual_full_n1_solution_from_bounds_only=True,
            **dict.fromkeys(OPEN,False))
        packet=OriginalN1InnerDomainPacket()
        self._packets[id(packet)]=(packet,raw,copy.deepcopy(parent.parent.first.encode(parent.parent.original.serialized(raw))))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1InnerDomainPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live inner domain packet issued by this owner required')
        if parent.parent.first.encode(parent.parent.original.serialized(entry[1]))!=entry[2]:raise ValueError('Domain packet changed')
        return dict(**copy.deepcopy(entry[1]),**{GATE:self.acceptance_loaded})


@source_precision
def run(n1):
    began=time.monotonic();owner=CurrentOriginalN1InnerAnalyticDomain(n1,require_checked=False)
    packet=owner.domain();raw=owner.report(packet)
    raw.update(input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(parent.parent.first.encode(parent.parent.original.serialized(raw)),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('ACTUAL_FIRST_COLLAR_N1_ANALYTIC_BOUNDS_READY',flush=True)
    return owner,packet
