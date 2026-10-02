"""Physical Cartesian vector derivatives of the accepted leading outer chain.

The moving cylindrical basis is differentiated exactly. All spatial multiindices
through order four and the first physical-time derivative are mapped on r>0.
Actual profile boxes and their original positive amplitude/radius scales are
consumed separately. This does not install the missing core/axis or time recursion.
"""
import functools
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    PulsePhysicalBounds, physical_bracket, physical_operators, abs_upper, UZ, UT, UR, P)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
CS, SN=s.symbols('cos_theta sin_theta',real=True)
COMPONENTS=('ux','uy','uz','p')
INDICES=tuple((i,j,b) for degree in range(5) for b in range(degree+1)
              for j in range(degree-b+1) for i in (degree-b-j,))


def transverse_step(row,direction):
    """Dx=c Dr-s/r Dtheta; Dy=s Dr+c/r Dtheta, including r^-q."""
    if direction not in ('x','y'):raise ValueError('Transverse direction required')
    out={}
    def add(index,value):out[index]=out.get(index,s.Integer(0))+value
    for (label,a,q),coefficient in row.items():
        angular=-SN*s.diff(coefficient,CS)+CS*s.diff(coefficient,SN)
        radial=CS if direction=='x' else SN
        angular_factor=-SN if direction=='x' else CS
        add((label,a+1,q),radial*coefficient)
        add((label,a,q+1),-q*radial*coefficient+angular_factor*angular)
    return {key:s.expand(value) for key,value in out.items() if s.expand(value)!=0}


@functools.lru_cache(maxsize=1)
def cartesian_templates():
    seeds={'ux':{(UR,0,0):CS,(UT,0,0):-SN},
           'uy':{(UR,0,0):SN,(UT,0,0):CS},
           'uz':{(UZ,0,0):s.Integer(1)},'p':{(P,0,0):s.Integer(1)}}
    out={}
    for component,seed in seeds.items():
        for i,j,b in INDICES:
            row=seed
            for _ in range(i):row=transverse_step(row,'x')
            for _ in range(j):row=transverse_step(row,'y')
            if any(a+q!=i+j for _,a,q in row):raise ArithmeticError('Transverse degree lost')
            out[(component,i,j,b)]=row
    return out


def angular_polynomial(c,expression,cosine,sine):
    value=c.mpf(0)
    for (i,j),coefficient in s.Poly(expression,CS,SN).terms():
        value+=(c.mpf(int(coefficient.p))/int(coefficient.q))*cosine**i*sine**j
    return value


def cartesian_brackets(c,grids,i,j,b,z,delta,cosine,sine,precomputed=None):
    """Separate input-component brackets; their lambda exponents differ.

    A term r^-q dr^a dz^b(lambda^beta G), a+q=N, equals
    lambda^(beta-N+b(delta-1))*R^(-N/2)*2^((a-q)/2)*H_ab.
    Original amplitude normalization is external to this identity.
    """
    if (i,j,b) not in INDICES:raise ValueError('Spatial total order <=4 required')
    beta={UR:c.mpf(-1),UT:-1-delta,UZ:-1-delta,P:-2-2*delta}
    cache={} if precomputed is None else precomputed
    out={}
    for component in COMPONENTS:
        contributions={}
        for (label,a,q),coefficient in cartesian_templates()[(component,i,j,b)].items():
            key=(label,a,b)
            if key not in cache:cache[key]=physical_bracket(c,grids[label],a,b,z,delta,beta[label])
            factor=angular_polynomial(c,coefficient,cosine,sine)*c.sqrt(2)**(a-q)
            contributions[label]=contributions.get(label,c.mpf(0))+factor*cache[key]
        out[component]=contributions
    return out


def physical_time_bracket(c,grid,z,delta,beta):
    """Lemma2.1 at fixed physical x, not derivative of a stage parameter."""
    return (-beta*grid['y0_Z0']/2+(1-delta)*z*grid['y0_Z1']/2
            +grid['y1_Z0'])/(1-delta*z*z)


class CompliantCartesianField:
    def __init__(self):
        self.pulse=PulsePhysicalBounds(); self.ctx=c=self.pulse.ctx
        self.mu=self.pulse.mu; self.delta=self.pulse.delta
        self.family=self.pulse.family; self.source=self.pulse.source
        self.hashes=dict(self.pulse.hashes); self.records={}
        required={'power_inlet_C4':'two_sided_O3_pulse_join_certified',
                  'flatten_mixed_C4':'entire_original_100_unit_flatten_mixed_C4_available',
                  'power_angular_C4':'flatten_power_and_power_angular_joins_certified',
                  'steep_waiting_C4':'angular_steep_and_internal_joins_certified',
                  'collar_Gamma_C4':'waiting_collar_and_collar_Gamma_joins_certified'}
        for name,flag in required.items():
            checkname=PREFIX+'compliant_'+name+'_check.json'
            check=json.loads((HERE/checkname).read_bytes())
            if not check['all_passed'] or not check[flag]:raise ValueError('Unaccepted Cartesian prerequisite: '+name)
            if check['actual_five_defect_family_sha256']!=self.family or check['implicit_source_sha256']!=self.source:
                raise ValueError('Cartesian family mismatch: '+name)
            for path,digest in check['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Cartesian dependency changed: '+path)
            self.hashes.update(check['input_hashes'])
            self.hashes[checkname]=hashlib.sha256((HERE/checkname).read_bytes()).hexdigest()
            filename=PREFIX+'compliant_'+name+'.json'
            self.records[name]=json.loads((HERE/filename).read_bytes())
            self.hashes[filename]=hashlib.sha256((HERE/filename).read_bytes()).hexdigest()
        self.Evparts={key:read_interval(c,value) for key,value in self.records['flatten_mixed_C4']['Ev0_squared_scale_log_parts'].items()}
        self.Lrel=read_interval(c,self.records['power_angular_C4']['postflatten_length'])
        steep=self.records['steep_waiting_C4']
        self.Ts=read_interval(c,steep['original_steep_power_length'])
        self.wait=read_interval(c,steep['original_refined_waiting_root'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def charts(self):
        c=self.ctx; rec=self.records
        yield 'local_O3_power',rec['power_inlet_C4']['whole_Z_left_box'],'pulse',c.mpf(0),c.mpf(-1),'log(R/Rp) in[-1,0]; preceding O3 not covered'
        for name,key,chi,offset in (
            ('pulse_entrance','whole_Z_entrance_box','0','0'),('pulse_main','whole_Z_main_box','.02','0'),
            ('pulse_exit','whole_Z_exit_box','10','0'),('pulse_gap_main','whole_Z_gap_main_box','11','0'),
            ('pulse_gap_end','whole_Z_gap_end_box','12','0'),('pulse_end','whole_Z_end_box','13','-4')):
            yield name,self.pulse.receipt[key],'pulse',c.mpf(chi),c.mpf(offset),'same complete O4 chart as accepted pulse receipt'
        specs=(('flatten',rec['flatten_mixed_C4']['whole_Z_flatten_box'],c.mpf(0)),
            ('following_power',rec['power_angular_C4']['whole_Z_power_box'],c.mpf(100)),
            ('angular',rec['power_angular_C4']['whole_Z_angular_box'],100+self.Lrel-4),
            ('steep_entry',rec['steep_waiting_C4']['whole_Z']['entry']['entire'],100+self.Lrel),
            ('steep_power',rec['steep_waiting_C4']['whole_Z']['power']['entire'],101+self.Lrel),
            ('steep_exit',rec['steep_waiting_C4']['whole_Z']['exit']['entire'],101+self.Lrel+self.Ts),
            ('waiting',rec['steep_waiting_C4']['whole_Z']['waiting']['entire'],102+self.Lrel+self.Ts),
            ('collar',rec['collar_Gamma_C4']['whole_Z_collar_box'],102+self.Lrel+self.Ts+self.wait),
            ('unbounded_Gamma',rec['collar_Gamma_C4']['whole_Z_unbounded_exterior'],105+self.Lrel+self.Ts+self.wait))
        for name,packet,offset in specs:
            yield name,packet,'Ev0',c.mpf(13),offset,'whole accepted stage; Rv=Rp*exp(13/mu); final Gamma offset unbounded'

    def scale(self,label,N,b,normalization,chi,offset,time=False):
        c=self.ctx; beta=c.mpf(-1) if label==UR else (-2-2*self.delta if label==P else -1-self.delta)
        gamma=beta-2 if time else beta-N+b*(self.delta-1)
        if endpoints(gamma)[1]>=0:raise ArithmeticError('Incorrect lambda bound direction')
        if normalization=='pulse':
            row=self.pulse.component_scale(label,N,b,chi,offset,0)
            parts={key:value for key,value in row.items() if key.endswith('_term') and key!='physical_time_term'}
            # The a/q-dependent factor of two is INSIDE Cartesian brackets.
            parts['finite_radius_term']-=c.mpf(N)*c.ln(2)/2
        elif normalization=='Ev0':
            radius_power=-c.mpf(N)/2; pressure=label==P
            parts=dict(logCstar_term=radius_power*self.pulse.logRp_parts['logCstar'],
                logPstar_term=(2 if pressure else 1)*self.pulse.logP+radius_power*self.pulse.logRp_parts['logPstar'],
                inverse_mu_decay_term=radius_power*chi/self.mu+(0 if pressure else self.Evparts['inverse_mu_term']/2),
                finite_radius_term=radius_power*(self.pulse.logRp_parts['finite_outer_offset']+offset),
                inlet_amplitude_term=c.mpf(0) if pressure else (self.Evparts['inlet_log']+self.Evparts['finite_offset'])/2)
        else:raise ValueError('Original scale normalization required')
        return dict(log_prefactor_parts=parts,physical_lambda_exponent=gamma,
                    log_tau_sector_terms={value:gamma*c.mpf(value)/2 for value in ('-1','-10','-100')},
                    positive_exact_scales_retained_formally=True,
                    prefactor_bound_uses_lambda_ge_sqrt_tau=True)

    def report(self):
        c=self.ctx; z=c.mpf([-1,1]); cosine=sine=c.mpf([-1,1]); charts={}
        for name,packet,normalization,chi,offset,scope in self.charts():
            grids={label:{index:read_interval(c,value) for index,value in values.items()}
                   for label,values in packet['physical_mixed_derivatives_total_order_le4'].items()}
            cache={}; rows={}; scales={}
            for i,j,b in INDICES:
                key='x'+str(i)+'_y'+str(j)+'_z'+str(b); N=i+j
                brackets=cartesian_brackets(c,grids,i,j,b,z,self.delta,cosine,sine,cache)
                rows[key]={}
                for component,contributions in brackets.items():
                    rows[key][component]={}
                    for label,value in contributions.items():
                        scale_key=label+'_transverse'+str(N)+'_z'+str(b)
                        if scale_key not in scales:scales[scale_key]=self.scale(label,N,b,normalization,chi,offset)
                        norm=abs_upper(c,value); zero=endpoints(norm)[1]==0
                        rows[key][component][label]=dict(bracket=value,absolute_upper=norm,
                            log_absolute_upper=None if zero else c.ln(norm),exactly_zero=zero,scale_key=scale_key)
            time_rows={}
            for label,grid in grids.items():
                beta=-1 if label==UR else (-2-2*self.delta if label==P else -1-self.delta)
                value=physical_time_bracket(c,grid,z,self.delta,beta); norm=abs_upper(c,value)
                time_rows[label]=dict(bracket=value,absolute_upper=norm,
                    exactly_zero=endpoints(norm)[1]==0,scale=self.scale(label,0,0,normalization,chi,offset,time=True))
            charts[name]=dict(normalization=normalization,minimum_log_radius_parts=dict(inverse_mu_coefficient=chi,finite_offset=offset),
                coordinate_scope=scope,spatial_multiindices=rows,shared_prefactor_bounds=scales,
                first_physical_time_derivative_cylindrical=time_rows,
                cartesian_time_basis={'ux':'cos(theta)*dt_ur-sin(theta)*dt_utheta',
                                      'uy':'sin(theta)*dt_ur+cos(theta)*dt_utheta','uz':'dt_uz','p':'dt_p'},
                scope='r>0, fixed tau>0, all physical Z in(-1,1), all theta; interval envelopes extend to Z=+-1')
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            selected_delta=self.delta,selected_mu=self.mu,charts=charts,
            templates={component:{'x'+str(i)+'_y'+str(j)+'_z'+str(b):
                [{'input_component':label,'r_derivative':a,'inverse_r_power':q,'angular_coefficient':str(coefficient)}
                 for (label,a,q),coefficient in cartesian_templates()[(component,i,j,b)].items()]
                for i,j,b in INDICES} for component in COMPONENTS},
            map_identity='r^-q dr^a dz^b(lambda^beta G)=lambda^(beta-N+b(delta-1))*R^(-N/2)*2^((a-q)/2)*H_ab; N=a+q',
            time_identity='dt(lambda^beta G)=lambda^(beta-2)/L*(-beta*G/2+(1-delta)*Z*G_Z/2+G_y); fixed physical x',
            pulse_radial_normalization='Ur=sqrt(Rp/2)*Pstar*exp(-mu*log(R/Rp))*grid; sqrt(R) already differentiated in source grid',
            Ev0_normalization='Ev0=Pstar*U*exp(-13/(2mu)-13); Rv=Rp*exp(13/mu), NOT Rp*exp(13/mu+26)',
            physical_bound_evaluation='Sum separate positive prefactor bounds times bracket norms; keep separated log parts until safe logsumexp evaluation',
            chart_count=len(charts),spatial_multiindex_count=len(INDICES),
            accepted_outer_chart_cartesian_spatial_C4_mapped=True,accepted_outer_chart_first_physical_time_derivative_mapped=True,
            full_cartesian_vector_derivatives_certified=False,core_axis_interfaces_certified=False,
            whole_background_installed=False,physical_energy_integral_certified=False,
            whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(270):result=CompliantCartesianField().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual local O3, O4 and full postpulse Cartesian vector spatial4 and fixed-x time derivative maps generated',flush=True)
    return result


if __name__=='__main__':run()
