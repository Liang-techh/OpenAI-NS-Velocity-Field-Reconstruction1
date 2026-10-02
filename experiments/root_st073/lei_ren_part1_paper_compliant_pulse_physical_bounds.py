"""Lemma2.1 physical r/z derivative map and complete O.4 log-bound ledger.

Profiles G(y,Z), y=log R, are mapped to lambda^beta G in physical space.
All mixed r/z derivatives through total order4 use actual O.4 interval
jets. Huge absolute radii/norms remain formal finite logarithms. This
supremum ledger is not an energy integral or an outer C4/interface proof.
"""
import functools
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
UZ='Uz_over_Pstar_without_common_theta_radial_factor'
UT='Utheta_over_Pstar_without_common_theta_radial_factor'
UR='Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor'
P='P_over_Pstar_squared'
ZSYM,DSYM,BSYM=s.symbols('Z delta beta',real=True)


@functools.lru_cache(maxsize=1)
def physical_operators():
    """Coefficients of H_ab: true mixed profile derivatives, not Taylor data.

    dr^a dz^b(lambda^beta G)=lambda^(beta-a+b(delta-1))
          *(2/R)^(a/2)*H_ab.
    H_a0=prod_j=0..a-1(Dy-j/2)G.
    H_a,b+1=((beta+b(delta-1))*Z*H+d*DZ H-2Z*Dy H)/L.
    The radial -a shift cancels when DZ differentiates R^(-a/2).
    """
    result={}; z,delta,beta=ZSYM,DSYM,BSYM; L=1-delta*z*z; d=1-z*z
    radial={(0,0):s.Integer(1)}
    for a in range(5):
        row=radial
        for b in range(5-a):
            result[(a,b)]=row
            out={}
            def add(index,value):out[index]=out.get(index,s.Integer(0))+value
            for (k,n),coefficient in row.items():
                add((k,n),((beta+b*(delta-1))*z*coefficient+d*s.diff(coefficient,z))/L)
                add((k,n+1),d*coefficient/L)
                add((k+1,n),-2*z*coefficient/L)
            row={index:s.cancel(value) for index,value in out.items() if value!=0}
        out={}
        for (k,n),coefficient in radial.items():
            out[(k+1,n)]=out.get((k+1,n),s.Integer(0))+coefficient
            out[(k,n)]=out.get((k,n),s.Integer(0))-s.Rational(a,2)*coefficient
        radial={index:value for index,value in out.items() if value!=0}
    return result


def interval_expression(c,expression,z,delta,beta):
    if expression==ZSYM:return z
    if expression==DSYM:return delta
    if expression==BSYM:return beta
    if expression.is_Rational:return c.mpf(int(expression.p))/int(expression.q)
    if expression.is_Add:return sum((interval_expression(c,v,z,delta,beta) for v in expression.args),c.mpf(0))
    if expression.is_Mul:
        value=c.mpf(1)
        for v in expression.args:value*=interval_expression(c,v,z,delta,beta)
        return value
    if expression.is_Pow and expression.args[1].is_Integer:
        return interval_expression(c,expression.args[0],z,delta,beta)**int(expression.args[1])
    raise ValueError('Unsupported exact operator expression: '+str(expression))


def physical_bracket(c,grid,a,b,z,delta,beta):
    result=c.mpf(0)
    for (k,n),expression in physical_operators()[(a,b)].items():
        result+=interval_expression(c,expression,z,delta,beta)*grid['y'+str(k)+'_Z'+str(n)]
    return result


def abs_upper(c,value):return c.mpf(max(abs(v) for v in endpoints(value)))


class PulsePhysicalBounds:
    def __init__(self):
        self.ctx=c=MPIntervalContext(); c.dps=240; self.hashes={}
        name=PREFIX+'compliant_pulse_mixed_C4_check.json'; r=json.loads((HERE/name).read_bytes())
        if not r['all_passed'] or not r['factored_all_chart_mixed_derivative_boxes_available']:
            raise ValueError('Actual all-chart mixed pulse derivatives required')
        self.hashes.update(r['input_hashes'])
        for source,digest in self.hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Physical map source changed: '+source)
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.receipt=json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())
        self.family=r['actual_five_defect_family_sha256']; self.source=r['implicit_source_sha256']
        self.mu=read_interval(c,self.receipt['selected_mu']); self.delta=read_interval(c,self.receipt['selected_delta'])
        norm_name=PREFIX+'compliant_physical_norm_family.json'; norm=json.loads((HERE/norm_name).read_bytes())
        if norm['implicit_source_sha256']!=self.source or not norm['full_physical_C3_K_norms_certified_for_uniform_analytic_family']:
            raise ValueError('Admitted same-source physical radius family required')
        self.hashes[norm_name]=hashlib.sha256((HERE/norm_name).read_bytes()).hexdigest()
        self.logC=read_interval(c,norm['selected_logCstar'])
        with mp.workdps(270):
            # Exact Md and defining parameter formula are part of the
            # source-bound compliant pressure receipt's immutable definition.
            parameter_name=PREFIX+'compliant_pressure_source.json'
            parameter_record=json.loads((HERE/parameter_name).read_bytes())['compliant_source']
            if parameter_record['implicit_source_sha256']!=self.source:raise ValueError('Physical parameter source mismatch')
            parameters=parameter_record['implicit_source_definition']
            if parameters['Md']!='40' or parameters['logPstar']!='exp(Md)+11':raise ValueError('Unexpected actual Md/logP source')
            self.hashes[parameter_name]=hashlib.sha256((HERE/parameter_name).read_bytes()).hexdigest()
            self.logP=c.exp(c.mpf(parameters['Md']))+11
            self.logmu=c.ln(c.mpf(parameters['c_mu']))-4*self.logP
            self.Tw=-60*self.logmu
            # Rref=110*(Cstar Pstar)^10; Rp/Rref=exp(yd+1+Tw).
            self.logRp_parts=dict(logCstar=10*self.logC,logPstar=10*self.logP,
                finite_outer_offset=c.ln(110)+self.logP+1+self.Tw)
            self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def component_scale(self,label,a,b,chi,offset,log_tau):
        """Conservative time/radial prefactors with exact inverse-mu parts.

        On a chart t=log(R/Rp)>=chi/mu+offset. All radial factor rates
        here are nonpositive. lambda>=sqrt(tau), tau=1-physical_time>0.
        Finite offsets are stored separately from huge radial origins.
        """
        c=self.ctx; delta=self.delta; mu=self.mu
        beta=-1 if label==UR else (-2-2*delta if label==P else -1-delta)
        gamma=beta-a+b*(delta-1)
        if endpoints(gamma)[1]>=0:raise ArithmeticError('Time-sector exponent must be negative')
        radial_power=c.mpf(1-a)/2 if label==UR else -c.mpf(a)/2
        two_power=c.mpf(a-1)/2 if label==UR else c.mpf(a)/2
        p_power=2 if label==P else 1
        decay_half=c.mpf(a)/2 if label in (UR,P) else c.mpf(a+1)/2
        decay_mu=0 if label==P else 1
        chi=c.mpf(chi); offset=c.mpf(offset)
        return dict(logCstar_term=radial_power*self.logRp_parts['logCstar'],
            logPstar_term=p_power*self.logP+radial_power*self.logRp_parts['logPstar'],
            finite_radius_term=radial_power*self.logRp_parts['finite_outer_offset']+two_power*c.ln(2),
            inverse_mu_decay_term=-chi*decay_half/mu,
            finite_decay_term=-chi*decay_mu-offset*(decay_half+decay_mu*mu),
            physical_time_term=gamma*c.mpf(log_tau)/2,
            physical_lambda_exponent=gamma,original_profile_lambda_exponent=beta,
            positive_finite_prefactor_retained_formally=True)

    def report(self):
        c=self.ctx; z=c.mpf([-1,1]); charts=(
            ('entrance','whole_Z_entrance_box','0','0'),
            ('main','whole_Z_main_box','.02','0'),
            ('exit','whole_Z_exit_box','10','0'),
            ('gap_main','whole_Z_gap_main_box','11','0'),
            ('gap_end','whole_Z_gap_end_box','12','0'),
            ('end','whole_Z_end_box','13','-4'))
        rows={}; templates={}
        with mp.workdps(270):
            for (a,b),operator in physical_operators().items():
                templates['r'+str(a)+'_z'+str(b)]={str(k)+','+str(n):str(v) for (k,n),v in operator.items()}
            for chart,key,chi,offset in charts:
                grids={label:{index:read_interval(c,v) for index,v in values.items()}
                    for label,values in self.receipt[key]['physical_mixed_derivatives_total_order_le4'].items()}
                rows[chart]={}
                for label in (UZ,UT,UR,P):
                    beta=-1 if label==UR else (-2-2*self.delta if label==P else -1-self.delta)
                    rows[chart][label]={}
                    for (a,b) in physical_operators():
                        bracket=physical_bracket(c,grids[label],a,b,z,self.delta,c.mpf(beta))
                        norm=abs_upper(c,bracket); zero=endpoints(norm)[1]==0
                        lognorm=None if zero else c.ln(norm)
                        sectors={}
                        for logtau in ('-1','-10','-100'):
                            scale=self.component_scale(label,a,b,chi,offset,logtau)
                            logparts={name:value for name,value in scale.items() if name.endswith('_term')}
                            sectors[logtau]=dict(**scale,factored_bracket_log_upper=lognorm,
                                physical_derivative_log_upper=None if zero else sum(logparts.values(),c.mpf(0))+lognorm,
                                derivative_exactly_zero=zero,
                                scope='sup over this complete O4 chart and all Z in[-1,1] at tau=exp(log_tau); cylindrical component derivative')
                        rows[chart][label]['r'+str(a)+'_z'+str(b)]=dict(factored_bracket=bracket,
                            factored_bracket_absolute_upper=norm,sectors=sectors)
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                actual_delta=self.delta,actual_mu=self.mu,logPstar=self.logP,logRp_parts=self.logRp_parts,
                physical_operator_templates=templates,all_chart_physical_derivative_log_bounds=rows,
                mapping=dict(R='r²/(2lambda²)',Z='z/lambda^(1-delta)',lambda_relation='lambda²-lambda^(2delta)z²=tau=1-physical_time',
                    velocity=dict(angular_axial='lambda^(-1-delta)*Utheta/Uz',radial='lambda^(-1)*Ur'),
                    pressure='lambda^(-2-2delta)*P',chain_rule='paper Lemma2.1 (2.8)-(2.9)',
                    actual_source_radial_origin='Rref=110*(Cstar*Pstar)^10; Rp=Rref*exp(yd+1+Tw)',
                    physical_time_domain='0<tau=1-t; log-sector receipts at -1,-10,-100',
                    all_positive_radii_and_scales_remain_exact_formal_sources=True),
                all_cylindrical_r_z_derivatives_total_order_le4_mapped=True,
                complete_pulse_physical_spatial_supremum_log_ledger_available=True,
                full_pulse_C4_installed=False,full_outer_C4_certified=False,
                full_cartesian_vector_derivatives_certified=False,physical_energy_integral_certified=False,
                whole_outer_cone_certified=False,temporal_recursion=False,
                next_dependency='Exact functional pulse interfaces and post-pulse high derivatives; vector Cartesian map, physical energy and stress cone',
                input_hashes=self.hashes)


def run():
    result=PulsePhysicalBounds().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual O4 physical r/z derivatives through total order4 and all-chart finite log-bound ledger generated',flush=True)
    return result


if __name__=='__main__':run()
