"""Actual-source Section 11 inverse/primitive majorants at fixed fast phase.

This defines bounds of the analytic loop, not point values of its signed
inputs. All source amplitudes, eta and inverse positive scales remain logs.
Only y, Z, yy and yZ are admitted; higher mixed orders are not inferred.
"""
import ast
import json
import math
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_uniform_inputs as source

packets,bounds=source.packets,source.bounds
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_generic_shear_loop_jet_bounds.json'
RECEIPT=PREFIX+'current_generic_shear_loop_jet_bounds_check.json'
GATE='current_actual_phase_held_inverse_and_loop_y2_yZ_log_bounds_certified'
OPEN=source.OPEN
ORDERS=((0,0),(1,0),(0,1),(2,0),(1,1))
FIRST=((1,0),(0,1))
SECOND=(((1,0),(1,0)),((1,0),(0,1)))
ZERO=(0,0)


def plus(i,j):return (i[0]+j[0],i[1]+j[1])


class BoundJet:
    """Absolute caps of ordinary derivatives, with an exact-zero sentinel."""
    def __init__(self,c,rows):
        if set(rows)!=set(ORDERS):raise ValueError('Exactly y,Z,yy,yZ derivative rows required')
        if any(v.ctx is not c for v in rows.values()):raise ValueError('Same directed context required')
        self.ctx,self.rows=c,rows
    def __getitem__(self,key):return self.rows[key]
    @classmethod
    def constant(cls,c,value):
        return cls(c,{k:bounds.LogUpper.constant(c,value if k==ZERO else 0) for k in ORDERS})
    def __add__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Same bound-jet context required')
        return BoundJet(self.ctx,{k:bounds.LogUpper.add(self.ctx,[self[k],other[k]]) for k in ORDERS})
    def __mul__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Same bound-jet context required')
        c=self.ctx;rows={}
        for j,k in ORDERS:
            rows[(j,k)]=bounds.LogUpper.add(c,[
                bounds.LogUpper.constant(c,math.comb(j,i)*math.comb(k,ell))*self[(i,ell)]*other[(j-i,k-ell)]
                for i in range(j+1) for ell in range(k+1)])
        return BoundJet(c,rows)
    def scale(self,cap):return BoundJet(self.ctx,{k:v*cap for k,v in self.rows.items()})
    def inverse_positive(self,log_lower):
        one=BoundJet.constant(self.ctx,1)
        return BoundJet(self.ctx,bounds.quotient_table(self.ctx,one.rows,self.rows,log_lower,ORDERS))
    def record(self):return bounds.encode_table(self.rows)


def decode_jet(c,rows):
    def read(row):
        zero=row['exact_zero'];value=row['log_absolute_upper']
        if zero!=(value is None):raise ValueError('Exact-zero source bound encoding differs')
        return bounds.LogUpper(c,None if zero else packets.interval(c,value))
    return BoundJet(c,{k:read(rows['y%d_Z%d'%k]) for k in ORDERS})


def loop_bounds(c,inputs,*,log_a_min,log_d,log_eta,log_q_star,sigma_caps):
    """Uniform analytic bounds; the formal q-flat branch has all jets zero."""
    ep=packets.recovery.endpoints
    if ep(log_eta+c.ln(2))[1]>0:raise ValueError('Actual eta<=1/2 required')
    if set(inputs)!=set(('a','b','p2','t0','E')):raise ValueError('Actual full source derivative caps required')
    C=lambda v:bounds.LogUpper.constant(c,v)
    L=lambda value:bounds.LogUpper(c,value)
    add=lambda *values:bounds.LogUpper.add(c,values)
    one,two,half=C(1),C(2),C('.5')
    a,b,p2,t0,E=(inputs[k] for k in ('a','b','p2','t0','E'))
    inva=L(-log_a_min);inveta=L(-log_eta)
    kappa=a+(b*b)*a.inverse_positive(log_a_min)
    # On active support gamma=2eta-(kappa-2)>=eta and gamma<3.
    # The same global cap also covers the flat, exactly-zero q branch.
    q0=L(log_q_star);sig={ZERO:one};g={ZERO:q0};ell={}
    for i in FIRST:
        sig[i]=sigma_caps[1]*kappa[i]*inveta
        ell[i]=half*add(kappa[i]*inveta,a[i]*inva)
        g[i]=q0*ell[i]
    for i,j in SECOND:
        k=plus(i,j)
        sig[k]=add(sigma_caps[2]*kappa[i]*kappa[j]*inveta*inveta,
                   sigma_caps[1]*kappa[k]*inveta)
        lij=half*add(kappa[k]*inveta,kappa[i]*kappa[j]*inveta*inveta,
                    a[k]*inva,a[i]*a[j]*inva*inva)
        g[k]=q0*add(lij,ell[i]*ell[j])
    q=BoundJet(c,sig)*BoundJet(c,g)
    u=(p2*q).scale(L(-log_d));u0=u[ZERO]
    # h>=1. beta bounds h^-1 derivatives; r derivatives keep the correlation.
    beta={ZERO:one};r={ZERO:one}
    for i in FIRST:
        beta[i]=u0*u[i];r[i]=u[i]
    for i,j in SECOND:
        k=plus(i,j)
        beta[k]=add(u[i]*u[j],u0*u[k],C(3)*u0*u0*u[i]*u[j])
        r[k]=add(u[k],C(3)*u0*u[i]*u[j])
    alpha=(q*BoundJet(c,beta)).scale(two)
    # 1-|r|=1/[h(h+|u|)]>=1/(2h^2), so D>=1/(4h^4).
    h2=add(one,u0*u0);invD=L(c.ln(4)+2*h2.log)
    w0=two*invD
    wr=add(invD,C(8)*invD*invD)
    wrr=add(C(12)*invD*invD,C(64)*invD*invD*invD)
    wpsi=invD*invD
    wpsir=add(two*invD*invD,C(8)*invD*invD*invD)
    w={ZERO:w0}
    for i in FIRST:w[i]=wr*r[i]
    for i,j in SECOND:w[plus(i,j)]=add(wrr*r[i]*r[j],wr*r[plus(i,j)])
    t=t0+alpha*BoundJet(c,w)
    tpsi=alpha[ZERO]*wpsi
    tpsii={i:add(alpha[i]*wpsi,alpha[ZERO]*wpsir*r[i]) for i in FIRST}
    # v/a=1+t0^2+2q^2 is correlated and >=1 on both branches.
    ratio=BoundJet.constant(c,1)+t0*t0+(q*q).scale(two)
    K=ratio.inverse_positive(c.mpf(0)).scale(C(1/(2*c.pi)))
    invlambda=ratio[ZERO]*C(2*c.pi)
    T1=t.scale(C(2*c.pi));T2=(t*t).scale(C(2*c.pi))
    P=T2+BoundJet.constant(c,2*c.pi);Phi=K*P
    lambdajet=K*(BoundJet.constant(c,1)+t*t)
    phipsipsi=two*K[ZERO]*t[ZERO]*tpsi
    psi={ZERO:C(2*c.pi)}
    for i in FIRST:psi[i]=Phi[i]*invlambda
    for i,j in SECOND:
        k=plus(i,j)
        psi[k]=add(Phi[k],lambdajet[i]*psi[j],lambdajet[j]*psi[i],
                   phipsipsi*psi[i]*psi[j])*invlambda
    chi={ZERO:one,**{k:psi[k]*C(1/(2*c.pi)) for k in ORDERS if k!=ZERO}}
    A=(a*BoundJet(c,chi)).scale(half)
    hat={ZERO:T1[ZERO]}
    for i in FIRST:hat[i]=add(T1[i],t[ZERO]*psi[i])
    for i,j in SECOND:
        k=plus(i,j)
        hat[k]=add(T1[k],t[i]*psi[j],t[j]*psi[i],tpsi*psi[i]*psi[j],t[ZERO]*psi[k])
    # |phi|<=1 on a period; derivatives are at fixed phi.
    M=(a*BoundJet(c,hat)).scale(C(1/(2*c.pi)))+b
    B=(E*M).scale(half)
    psiphi=invlambda;psiphi2=phipsipsi*invlambda*invlambda*invlambda
    psii_phi={i:add(lambdajet[i],phipsipsi*psi[i])*invlambda*invlambda for i in FIRST}
    Aphi=half*add(a[ZERO],a[ZERO]*psiphi*C(1/(2*c.pi)))
    Aphi2=a[ZERO]*psiphi2*C(1/(4*c.pi))
    Mphi=add(a[ZERO]*t[ZERO]*psiphi*C(1/(2*c.pi)),b[ZERO])
    Mphi2=a[ZERO]*add(tpsi*psiphi*psiphi,t[ZERO]*psiphi2)*C(1/(2*c.pi))
    Bphi=half*E[ZERO]*Mphi;Bphi2=half*E[ZERO]*Mphi2
    fastA={'phi':Aphi,'phi_phi':Aphi2};fastB={'phi':Bphi,'phi_phi':Bphi2}
    for i in FIRST:
        label=('y' if i==(1,0) else 'Z')+'_phi'
        fastA[label]=half*add(a[i],add(a[i]*psiphi,a[ZERO]*psii_phi[i])*C(1/(2*c.pi)))
        Mi_phi=add(add(a[i]*t[ZERO]*psiphi,
            a[ZERO]*add(t[i],tpsi*psi[i])*psiphi,
            a[ZERO]*t[ZERO]*psii_phi[i])*C(1/(2*c.pi)),b[i])
        fastB[label]=half*add(E[i]*Mphi,E[ZERO]*Mi_phi)
    return dict(slow_phase_held_log_bounds={key:value.record() for key,value in
        dict(q=q,kappa=kappa,u=u,t=t,T1_fixed_psi=T1,T2_fixed_psi=T2,
             Phi_fixed_psi=Phi,psi=BoundJet(c,psi),A=A,B=B).items()},
        phase_and_first_slow_phase_log_bounds={
            'psi_phi':psiphi.record(),'psi_phi_phi':psiphi2.record(),
            **{'psi_'+('y' if i==(1,0) else 'Z')+'_phi':v.record() for i,v in psii_phi.items()},
            **{'A_'+key:v.record() for key,v in fastA.items()},
            **{'B_'+key:v.record() for key,v in fastB.items()}},
        positive_conditioning=dict(log_lambda_positive_lower=-invlambda.log,
            log_Poisson_D_positive_lower=-invD.log,
            correlated_v_over_a='1+t0^2+2q^2>=1',
            correlated_one_minus_abs_r='1/[h*(h+abs(u))]>=1/(2*h^2)',
            correlated_one_minus_r_squared='1/h^2',
            eta_and_Delta_not_added_to_or_subtracted_from_rounded2=True),
        derivative_orders=[list(k) for k in ORDERS],
        q_flat_branch=dict(condition='Delta=kappa-2>=eta',q_and_all_jets_exact_zero=True,
            phase='Phi=psi/(2pi)',inverse='psi=2pi*phi',A_and_B_all_jets_exact_zero=True),
        primitive_units=dict(A='dimensionless original A',B='original B/Pstar; E=Utheta/Pstar',
            Pstar_is_global_constant_under_slow_derivatives=True),
        current_source_signed_jets_materialized=False,
        input_caps_are_not_used_as_field_values=True)


def exact_theorem():
    """Bind the defining cutoff and check implicit/primitive chain rules."""
    tree=ast.parse((HERE/(PREFIX+'current_generic_shear_loop.py')).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='flat_step')
    odds=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='odds' for t in n.targets))
    if ast.dump(odds)!=ast.dump(ast.parse('1/(1-x)**2-1/x**2',mode='eval').body):
        raise ValueError('Same original flat step odds required')
    y,Z=s.symbols('y Z');psi=s.Function('psi')(y,Z)
    T=s.Function('T1');aa=s.Function('a')(y,Z);bb=s.Function('b')(y,Z);EE=s.Function('E')(y,Z)
    w=s.Symbol('w');cc=s.symbols('c0:10')
    # A fully general quadratic jet represents every partial derivative
    # entering this second-order chain rule, without Subs dummy ambiguity.
    poly=cc[0]+cc[1]*w+cc[2]*y+cc[3]*Z+cc[4]*w*w+cc[5]*w*y+cc[6]*w*Z+cc[7]*y*y+cc[8]*y*Z+cc[9]*Z*Z
    phi=s.Symbol('fixed_phi');comp=poly.subs(w,psi);that=T(psi,y,Z)
    checks={}
    fp=s.diff(poly,w).subs(w,psi)
    # Independent total derivatives of the defining composite equations.
    for i,j in ((y,y),(y,Z)):
        si,sj=s.diff(psi,i),s.diff(psi,j)
        sii=s.diff(psi,i,j)
        f_i=s.diff(poly,i).subs(w,psi)
        f_ij=s.diff(poly,i,j).subs(w,psi)
        f_pi=s.diff(poly,w,i).subs(w,psi)
        f_pj=s.diff(poly,w,j).subs(w,psi)
        f_pp=s.diff(poly,w,2).subs(w,psi)
        first=f_i+fp*si
        second=f_ij+f_pi*sj+f_pj*si+f_pp*si*sj+fp*sii
        if s.simplify(s.diff(comp,i)-first)!=0 or s.simplify(s.diff(comp,i,j)-second)!=0:
            raise ArithmeticError('Implicit inverse phase chain rule failed')
        checks['inverse_'+str(i)+str(j)]=True
    chi=phi-psi/(2*s.pi);A=aa*chi/2;M=-aa*that/(2*s.pi)-bb*phi;B=EE*M/2
    for i in (y,Z):
        if s.expand(s.diff(A,i)-(s.diff(aa,i)*chi/2-aa*s.diff(psi,i)/(4*s.pi)))!=0:
            raise ArithmeticError('Primitive first derivative failed')
        checks['A_'+str(i)]=True
    for i,j in ((y,y),(y,Z)):
        expected=s.diff(aa,i,j)*chi/2-(s.diff(aa,i)*s.diff(psi,j)+s.diff(aa,j)*s.diff(psi,i)+aa*s.diff(psi,i,j))/(4*s.pi)
        expectedB=(s.diff(EE,i,j)*M+s.diff(EE,i)*s.diff(M,j)+s.diff(EE,j)*s.diff(M,i)+EE*s.diff(M,i,j))/2
        if s.expand(s.diff(A,i,j)-expected)!=0 or s.expand(s.diff(B,i,j)-expectedB)!=0:
            raise ArithmeticError('Primitive second derivative failed')
        checks['A_B_'+str(i)+str(j)]=True
    a,b,q=s.symbols('a b q',positive=True)
    if s.cancel((a+b*b/a+2*a*q*q)/a-(1+(b/a)**2+2*q*q))!=0:
        raise ArithmeticError('Correlated positive phase coefficient failed')
    checks['positive_v_over_a_identity']=True
    return dict(passed=True,exact_chain_rule_identities=checks,original_flat_step_odds_AST_bound=True,
        fixed_global_d_star_eta_and_N_under_slow_derivatives=True,
        fast_phase_total_derivatives=dict(y='f_y_at_phi+N*f_phi',yy='f_yy_at_phi+2N*f_y_phi+N^2*f_phi_phi',
            yZ='f_yZ_at_phi+N*f_Z_phi'),
        modulate_consumes_slow_A_y_B_y_not_total_y=True)


class CurrentLoopJetBounds:
    def __init__(self):
        self.ctx=c=MPIntervalContext();c.dps=110
        self.uniform=json.loads((HERE/source.NAME).read_bytes())
        receipt=json.loads((HERE/source.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(source.GATE):
            raise ValueError('Checked actual whole generic input/log scales required')
        self.family=self.uniform['source_family']
        if receipt['source_family']!=self.family:raise ValueError('Same actual original source family required')
        self.hashes=dict(receipt['input_hashes']);self.hashes[source.RECEIPT]=sha(source.RECEIPT)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prior=json.loads((HERE/bounds.NAME).read_bytes())
        self.outer=json.loads((HERE/source.source.NAME).read_bytes())
        for record in (self.prior,self.outer):
            if record['source_family']!=self.family:raise ValueError('Source quotient bound family differs')
        flat=PREFIX+'flat_pulse_derivatives';flatcheck=json.loads((HERE/(flat+'_check.json')).read_bytes())
        if not flatcheck['all_passed'] or not flatcheck['original_radial_shape_derivatives_C4_available']:
            raise ValueError('Checked original global sigma derivative caps required')
        if flatcheck['actual_five_defect_family_sha256']!=self.family['actual_five_defect_family_sha256']:
            raise ValueError('Same original cutoff/source family required')
        for name,digest in flatcheck['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Cutoff/source hash family differs')
            self.hashes[name]=digest
        self.hashes[flat+'_check.json']=sha(flat+'_check.json')
        row=json.loads((HERE/(flat+'.json')).read_bytes())
        self.sigma_caps={n:bounds.LogUpper.constant(c,packets.interval(c,row['sigma_global_derivative_bounds'][n])) for n in (1,2)}
        self.hashes[PREFIX+'current_generic_shear_loop.py']=sha(PREFIX+'current_generic_shear_loop.py')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed actual loop derivative prerequisite: '+name)
        self.scales=self.uniform['current_actual_logarithmic_loop_scales'];self.theorem=exact_theorem()

    def chart(self,chart):
        c=self.ctx
        prior=self.prior['current_original_source_log_bound_charts'];outer=self.outer['original_O3_quotient_log_norms']
        row=prior[chart] if chart in prior else outer[chart]
        q=row['admitted_original_quotient_log_norms']
        inputs={key:decode_jet(c,q[key]) for key in ('a','b','p2','t0')}
        inputs['E']=decode_jet(c,row['ordinary_mixed_source_log_norms']['E'])
        pos=row['actual_denominator_source_certificate'] if chart in prior else row['actual_positive_denominator_theorem']
        scales=self.scales;read=lambda v:packets.interval(c,v)
        result=loop_bounds(c,inputs,log_a_min=read(pos['log_actual_a_positive_lower']),
            log_d=read(scales['logarithmic_selected_positive_lower_constants']['d_star']),
            log_eta=read(scales['selected_positive_eta_log']),
            log_q_star=read(scales['logarithmic_conservative_upper_constants']['q_star']),sigma_caps=self.sigma_caps)
        result.update(chart=chart,source_family=self.family,actual_positive_denominator_theorem=pos,
            original_signed_source_definition_receipt=PREFIX+'current_generic_shear_inputs_check.json',
            original_quotient_caps_and_full_inertial_pressure_energy_retained=True)
        return result

    def run(self):
        charts={chart:self.chart(chart) for chart in self.scales['whole_source_chart_inventory']}
        result=dict(source_family=self.family,current_actual_loop_jet_log_bounds_by_chart=charts,
            source_chart_count=len(charts),exact_inverse_primitive_theorem=self.theorem,
            sigma_global_derivative_log_caps={str(n):v.record() for n,v in self.sigma_caps.items()},
            current_actual_logarithmic_scales=self.scales,
            **{GATE:True},**dict.fromkeys(OPEN,False),
            phase_held_loop_primitive_derivative_bounds_certified=True,
            certified_phase_held_orders=[list(k) for k in ORDERS],
            Z2_y2Z_or_full_mixed4_loop_bounds_certified=False,
            signed_current_point_loop_or_inverse_jets_installed=False,
            source_graph_ancestor_constructors_called=False,
            scope='Actual analytic Section11 loop uniform log majorants on17 original source covers, fixed-phase inverse/A/B y,Z,yy,yZ and first slow-phase derivatives. No cap used as a field value, no materialized microscopic source, installed point loop, changed moments/new repair/N, global cone or recursion.',
            input_hashes=self.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        print('Actual-source phase-held inverse/A/B logarithmic derivative bounds PASS:17 charts',flush=True)
        return result


def run():return CurrentLoopJetBounds().run()


if __name__=='__main__':run()
