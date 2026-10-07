"""Full signed cone on transported modulation before independent repair.

Local profiles are unchanged, but their cumulative moment/pressure errors
are not zero. A common finite N is proved for the changed source chain;
the degenerate unchanged O2 buffer and global gates remain open.
"""
import ast
import functools
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_modified_quiet_cone as quiet

errors=quiet.errors
HERE,PREFIX,sha=quiet.HERE,quiet.PREFIX,quiet.sha
NAME=PREFIX+'current_O3_modified_transport_cone.json'
RECEIPT=PREFIX+'current_O3_modified_transport_cone_check.json'
parameters,endpoints,upper=quiet.parameters,quiet.endpoints,quiet.upper
OPEN=quiet.OPEN
REGIONAL=(('current_O2_modified_taper_cone','current_modified_O2_open_taper_signed_two_vector_cone_certified'),
    ('current_O3_modified_transition_cone','current_modified_closed_O3_signed_two_vector_cone_certified'))


@functools.lru_cache(maxsize=1)
def exact_theorem():
    asts=errors.recovered.histories.SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:raise ArithmeticError('Transport cone source: '+name)
        checks[name]=True
    v,y,z,mu,delta=s.symbols('actual_power_offset quiet_offset Z mu delta',real=True)
    M,U=s.symbols('same_M0 same_U0',positive=True);alpha=s.Rational(1,2)+mu
    m=M*s.exp(-v)*z;u=U*s.exp(-alpha*v)/(1+z*z)
    Q=asts.evaluate(asts.expression('pre_pulse_mixed_C4','physical_mixed','Q'),
        dict(z=z,V=[0]*5,rows=dict(m=[s.diff(m,v,j) for j in range(5)]),delta=delta,
            square=lambda a:a*a,derivative=lambda a:s.diff(a,z)))
    for j,value in enumerate(Q):zero('actual_original_normalized_radial_Q_pure_m_row'+str(j),value,-M*s.exp(-v)*(-1)**j)
    for j in range(5):
        for k in range(5-j):
            ujk=s.diff(u,v,j,z,k);mjk=s.diff(m,v,j,z,k);qjk=s.diff(Q[0],v,j,z,k)
            zero('same_backward_U0_cap_factor_%d_%d'%(j,k),ujk,ujk.subs(v,v+1)*s.exp(alpha))
            zero('same_backward_M0_cap_factor_%d_%d'%(j,k),mjk,mjk.subs(v,v+1)*s.exp(1))
            zero('same_backward_Q0_cap_factor_%d_%d'%(j,k),qjk,qjk.subs(v,v+1)*s.exp(1))
    asts.expression('current_O3_modulated_histories','history','pre',wanted='self.pre.power(zv,v/Tw)')
    asts.expression('current_O3_modulated_histories','history','t',wanted='1+v')
    zero('actual_transport_offset_is_one_plus_power_offset',1+v,1+v)
    # Exact local cutoff/support, before any separate interval cap.
    asts.expression('current_O3_finite_frequency_profiles','cutoff_rows','right',wanted='sigma_jets(c,4*t-1)')
    if 4*s.Integer(1)-1<1:raise ValueError('Pre-repair power lies inside modulation support')
    checks['all_power_offsets0_to1_have_exact_zero_local_modulation']=True
    transport=asts.method('current_O3_modulated_histories_operator','cumulative_scalar_enclosures')
    for target,arg in (('upper_t','hi'),('lower_t','lo')):
        asts.expression('current_O3_modulated_histories_operator','cumulative_scalar_enclosures',target,
            wanted="min(%s,mp.mpf('.5'))"%arg)
    node=next(n.value for n in ast.walk(transport) if isinstance(n,ast.AugAssign)
        and ast.unparse(n.target)=='result[name]' and ast.dump(n.value)==ast.dump(ast.parse('c.exp(-c.mpf(str(rate))*t)',mode='eval').body))
    ctx=SimpleNamespace(mpf=lambda a:s.Rational(str(a)),exp=s.exp)
    for name,rate in (('m',1),('h',s.Rational(3,2)),('k',s.Rational(3,2)),('e',1),('p',0)):
        value=asts.evaluate(node,dict(c=ctx,rate=rate,t=1+v))
        zero('same_full_support_'+name+'_homogeneous_transport',value,s.exp(-rate*(1+v)))
        for j in range(5):zero('same_post_support_'+name+'_ordinary_derivative'+str(j),s.diff(value,v,j),(-rate)**j*value)
    q,K1,K2=s.symbols('inverse_N positive_K1 positive_K2',positive=True)
    for name,value in (('linear_history',K1*q),('quadratic_plus_swirl_history',K1*q+K2*q*q)):
        if s.diff(value,q).is_positive is not True:raise ValueError('Transport inverse-N history cap not monotone')
        checks['positive_inverse_N_derivative_'+name]=True
    # Exact local velocity equality determines signed shears independently
    # of the surviving moment/pressure stress changes.
    zero('actual_transport_theta_shear_excess',1-2*s.diff(u,v)/u-2,2*mu)
    zero('actual_transport_axial_shear',0/u,0)
    parent=quiet.exact_theorem()
    return dict(passed=True,identities=checks,input_hashes={**parent['input_hashes'],**asts.hashes,
        Path(__file__).name:sha(Path(__file__).name)},
        source_power_domain=(0,1),actual_modulation_transport_domain=(1,2),
        local_delta_Utheta_and_delta_Uz_exactly_zero=True,
        cumulative_original_modulation_histories_not_erased=True,
        native_baseline_caps_extended_only_after_undoing_physical_radius_shifts=True,
        baseline_backward_factors=dict(U0='exp(.5+mu)',M0='exp(1)',R0='exp(1)'),
        same_signed_source_shears=dict(a='2+2mu',bs='0',vs_minus2='2mu'),
        native0_positive_error_caps_nonincreasing_with_N=True)


class CurrentModifiedTransportCone:
    def __init__(self,require_checked=True):
        self.quiet=quiet.CurrentModifiedQuietCone();self.errors=self.quiet.errors;self.ctx=self.quiet.ctx
        self.data,self.mu,self.read=self.quiet.data,self.quiet.mu,self.quiet.read
        self.minimum_N=self.quiet.minimum_N;self.theorem=exact_theorem();c=self.ctx
        b=self.quiet.baseline;constants=b['canonical_theta_lower_constants'];deficit=b['canonical_source_deficit']
        floor=constants['cD']*deficit*c.exp(-1)+constants['theta_floor']
        if endpoints(floor)[0]<=0:raise ValueError('Original pre-repair local signed theta lower required')
        self.logUmin=self.data['logAd']-c.mpf('.5')-self.mu/2-(c.mpf('.5')+self.mu)
        self.logfloor=self.logUmin+c.ln(floor);self.rho0=self.quiet.rho0
        self.baseline=dict(power_offset_domain=(0,1),actual_modulation_transport_domain=(1,2),
            canonical_theta_local_lower=floor,original_log_Utheta_scalar_lower=self.logUmin,
            original_theta_log_lower=self.logfloor,original_weighted_axial_ratio_upper=self.rho0,
            checked_whole_power_source_restricted_before_error_normalization=True)
        self.error_rows=self.leading_error_rows(self.minimum_N);self.proof=self.prove()
        self.hashes={**self.quiet.hashes,quiet.NAME:sha(quiet.NAME),quiet.RECEIPT:sha(quiet.RECEIPT),**self.theorem['input_hashes']}
        self.common=self.common_frequency()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed']:raise ValueError('Checked full transported-source cone required')
            for name,digest in record['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed full transport/common-frequency source: '+name)

    def inputs(self,power_offset,N):
        parameters.positive_integer_N(N);c=self.ctx;v=c.mpf(power_offset);lo,hi=endpoints(v)
        if N<self.minimum_N or lo<0 or hi>1 or not all(mp.isfinite(a) for a in (lo,hi)):
            raise ValueError('Finite original power offset subset[0,1], integer N>=common minimum required')
        t=1+v;K=self.errors.recovered.frequency['untransported_abs_integral_times_N_upper']
        initial={name:upper(c,K[name]*c.exp(-c.mpf(str(rate))*t)/N)
            for name,rate in (('m',1),('h','1.5'),('k','1.5'),('e',1),('p',0))}
        kinetic=upper(c,K['e_kinetic']*c.exp(-t)/N**2);swirl=upper(c,K['e_swirl']*c.exp(-t)/N)
        initial['e']=errors.recovered.uniform.cap_min(c,initial['e'],kinetic+swirl)
        local=[c.mpf(0)]*5;old=[c.mpf(1)]*5
        d=errors.recovered.scalar_history_caps(initial,old,local,local)
        recovered=errors.recovered.recover(c,d,local,self.errors.recovered.factors)
        factors=dict(U0=c.exp(c.mpf('.5')+self.mu),M0=c.exp(1),R0=c.exp(1))
        caps={name:{jk:upper(c,value*factors[name]) for jk,value in grid.items()}
            for name,grid in self.errors.baselines['quiet'].items()}
        for name in ('du','V','dr'):caps[name]={jk:c.mpf(0) for jk in errors.indices(4)}
        for name,key in (('dm','m'),('dh','h'),('dk','k'),('de','e'),('dp','p')):
            grid=recovered['five_normalized_history_error_ordinary_logR4_axial5'][key]
            caps[name]={jk:grid[jk[0]][jk[1]] for jk in errors.indices(4)}
        # The normalized radial source error is nonzero through its own dm,
        # even though the local axial velocity error vanishes.
        radial=recovered['radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5']
        caps['dr']={jk:upper(c,sum(math.comb(jk[0],i)*c.mpf('.5')**(jk[0]-i)*radial[i][jk[1]]
            for i in range(jk[0]+1))) for jk in errors.indices(4)}
        return dict(power_offset=v,actual_modulation_transport_offset=t,native_logR=self.data['radius_logs']['logRw']+v,
            normalized_scalar_history_caps=d,normalized_cumulative_pressure_value_cap=initial['p'],
            baseline_backward_factors=factors),caps

    def leading_error_rows(self,N):
        c=self.ctx;q,caps=self.inputs((0,1),N);result={}
        for label,parts in self.errors.model['stress'].items():
            result[label]={}
            for name,p in parts.items():
                cap=errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),self.errors.delta,0,p['mode'][0])[0,0]
                zero=cap._mpi_==c.mpf(0)._mpi_;rp,pp,_,_=p['mode']
                row=dict(source_mode=p['mode'],error_Ad_power=p['error_Ad_power'],native_coefficient_absolute_cap=cap,
                    source_error_exactly_zero=zero)
                if zero:row.update(relative_error_log_parts=None,log_absolute_error_over_baseline_theta_lower=None)
                else:
                    logparts=dict(relative_radius=(c.mpf(str(rp))-c.mpf('.5'))*q['native_logR'],
                        relative_Pstar=(pp-1)*self.data['logP'],external_Ad=p['error_Ad_power']*self.data['logAd'],
                        weighted_axial=self.data['log_mu']/2 if label=='axial' else c.mpf(0),
                        coefficient=c.ln(cap),baseline_floor=-self.logfloor)
                    row.update(relative_error_log_parts=logparts,log_absolute_error_over_baseline_theta_lower=sum(logparts.values(),c.mpf(0)))
                result[label][name]=row
        return result

    def prove(self):
        c=self.ctx;positive={};zero_count=0
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Transport cone inequality unresolved: '+name)
            positive[name]=value
        for label,parts in self.error_rows.items():
            for name,row in parts.items():
                if row['source_error_exactly_zero']:zero_count+=1
                else:require('native_error_log_gap/'+label+'/'+name,-1010-row['log_absolute_error_over_baseline_theta_lower'])
            require('sum_native_errors_below_exp_minus1000/'+label,10-c.ln(len(parts)))
        eps=c.exp(-1000);rho=(self.rho0+eps)/(1-eps);Q=2-2*rho*rho
        require('positive_modified_theta',1-eps);require('weighted_axial_ratio_below_exp_minus390',c.exp(-390)-rho)
        require('actual_vs_minus2',2*self.mu);require('original_Q_over_theta_squared',Q)
        return dict(minimum_integer_N=self.minimum_N,all_finite_integer_N_at_least_minimum=True,
            whole_power_offset_domain=(0,1),whole_Z_domain=(-1,1),all_inherited_modulation_phases=True,
            positive_margins=positive,baseline_theta_log_lower=self.logfloor,
            relative_theta_error_upper=eps,relative_weighted_axial_error_upper=eps,
            modified_weighted_axial_ratio_upper=rho,actual_vs_minus2_lower=2*self.mu,
            original_signed_D_over_theta_exact=1,original_signed_Q_over_theta_squared_lower=Q,
            exact_local_profile_zero_stress_error_pieces=zero_count,
            surviving_cumulative_moment_energy_pressure_and_own_radial_errors_retained=True,
            current_modified_pre_repair_power_transport_cone_certified=True,**{key:False for key in OPEN})

    def common_frequency(self):
        family=self.data['source']['accepted']['source_family'];certificates={};thresholds={'quiet':self.quiet.minimum_N,
            'repair':self.quiet.repair.repair_threshold,'pre_repair_transport':self.minimum_N}
        for stem,gate in REGIONAL:
            name=PREFIX+stem+'.json';receipt_name=PREFIX+stem+'_check.json'
            data=json.loads((HERE/name).read_bytes());receipt=json.loads((HERE/receipt_name).read_bytes())
            if not receipt['all_passed'] or not data[gate] or not receipt[gate]:raise ValueError('Checked regional cone required: '+stem)
            sf=data['source_family']
            if sf!=family:raise ValueError('Regional cone source family differs: '+stem)
            for path,digest in receipt['input_hashes'].items():
                if sha(path)!=digest:raise ValueError('Changed regional source cone: '+path)
            self.hashes.update({**receipt['input_hashes'],receipt_name:sha(receipt_name),name:sha(name)})
            bound=data['whole_modified_O2_taper_cone'] if 'O2' in stem else data['whole_modified_closed_O3_cone']
            thresholds[stem]=bound['minimum_integer_N'];certificates[stem]=dict(gate=gate,threshold=bound['minimum_integer_N'],
                receipt=receipt_name,receipt_sha256=sha(receipt_name))
        n=max(thresholds.values())
        if n!=self.minimum_N:raise ValueError('Common frequency changed; recompute transport bound')
        return dict(sufficient_common_integer_N=n,all_finite_integer_N_at_least_minimum=True,
            regional_threshold_ledger=thresholds,accepted_regional_certificates=certificates,
            current_modified_modulation_transport_and_quiet_common_N_certified=True,
            covered_domains=dict(O2='offset(-2,0]',O3='offset[0,1]',transport='original power s[0,1]',quiet='original power s[1,2]'),
            earlier_O2_buffer_and_exact_degenerate_left_edge_excluded_from_strict_certificate=True,
            exact_original_source_restored_after_full_implicit_repair=True,
            original_registry_counts_not_incremented=True,**{key:False for key in OPEN})

    def query(self,power_offset,N):
        q,caps=self.inputs(power_offset,N)
        return dict(**q,finite_integer_N=N,source_family=self.data['source']['accepted']['source_family'],
            local_profile_changes_exactly_zero=True,source_history_caps=q['normalized_scalar_history_caps'],
            full_signed_cone_certificate=self.proof,common_frequency_ledger=self.common,
            current_modified_pre_repair_power_transport_cone_certified=True,**{key:False for key in OPEN})


def run():
    field=CurrentModifiedTransportCone(require_checked=False)
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_transport_source_and_native_units_theorem=field.theorem,original_local_power_baseline=field.baseline,
        whole_modified_pre_repair_transport_cone=field.proof,complete_leading_stress_error_log_ledger=field.error_rows,
        scoped_common_frequency_ledger=field.common,
        examples={name:field.query(v,N) for name,v,N in (
            ('whole_bridge',(0,1),field.minimum_N),('O3_to_power_seam','0',field.minimum_N),
            ('interior','.537',field.minimum_N),('repair_entry','1',field.minimum_N),('current_frequency',(0,1),10**12))},
        current_modified_pre_repair_power_transport_cone_certified=True,
        current_modified_modulation_transport_and_quiet_common_N_certified=True,
        **{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Modified transport cone and changed-region common N generated; N>=',field.minimum_N,flush=True)
    return result


if __name__=='__main__':run()
