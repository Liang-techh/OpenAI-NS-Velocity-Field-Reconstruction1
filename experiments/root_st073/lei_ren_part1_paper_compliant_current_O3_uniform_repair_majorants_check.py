"""Focused acceptance of uniform repair estimates and support-strip bounds."""
import json
import math
from pathlib import Path
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_uniform_repair_majorants as source


def point(value):
    lo,hi=source.endpoints(value)
    return (lo+hi)/2


def numerical_checks(field):
    counts=dict(implicit_Jacobian_ball_fixtures=0,independent_raw_beta_jet_samples=0,
                independent_log_bump_and_quiet_shear_fixtures=0,independent_partial_weight_integrals=0)
    with mp.workdps(90):
        mu=point(field.mu);R=point(field.R);c=field.ctx
        B=mp.matrix([[point(v) for v in row] for row in field.matrix['linear']])
        cross=[point(v) for v in field.matrix['cross_weights']]
        energy=[point(v) for v in field.matrix['energy_weights']]
        pressure=[point(v) for v in field.matrix['pressure_weights']]
        for N in (field.repair_threshold,field.quiet_threshold,10**12):
            query=field.query(N)
            for coefficients in ((0,0,0,0,0),(1,-1,1,-1,1),(-1,1,-1,1,-1),('.3','-.7','.8','-.2','.4')):
                a0,a2,e0,e1,e2=[mp.mpf(v)*R for v in coefficients]
                J=mp.matrix([[0]*5 for _ in range(5)])
                J[1,0]=cross[0]*e0/N;J[1,1]=cross[1]*e2/N
                J[1,2]=cross[0]*a0/N;J[1,4]=cross[1]*a2/N
                J[3,0]=2*energy[0]*a0/N;J[3,1]=2*energy[2]*a2/N
                for i,e in enumerate((e0,e1,e2)):
                    J[3,i+2]=-mu*energy[i]*e/N;J[4,i+2]=mu*pressure[i]*e/N
                norm=max(sum(abs(J[i,j]) for j in range(5)) for i in range(5))
                if norm>source.endpoints(query['quadratic_Jacobian_infinity_norm_upper'])[1]:
                    raise ArithmeticError('Independent quadratic Jacobian exceeds analytic cap')
                inverse=(B+J)**-1
                norm=max(sum(abs(inverse[i,j]) for j in range(5)) for i in range(5))
                if norm>source.endpoints(query['inverse_nonlinear_Jacobian_infinity_norm_upper'])[1]:
                    raise ArithmeticError('Independent nonlinear inverse exceeds Neumann bound')
                counts['implicit_Jacobian_ball_fixtures']+=1
        r=s.Symbol('r');raw=s.exp(-1/(1-r*r))
        derivatives=[s.lambdify(r,s.diff(raw,r,k),'mpmath') for k in range(5)]
        jets=field.jets
        # The value at r=0 attains exp(-1); evaluate it more accurately than
        # the 240-digit outward cap rather than adding a comparison tolerance.
        with mp.workdps(280):
            for rv in ('-.999','-.9','-.5','-.1','0','.1','.5','.9','.999'):
                for k,fn in enumerate(derivatives):
                    if abs(fn(mp.mpf(rv)))>source.endpoints(jets['raw_ordinary_derivative_caps'][k])[1]:
                        raise ArithmeticError('Independent raw beta derivative exceeds cap')
                    counts['independent_raw_beta_jet_samples']+=1
        norm=mp.quad(lambda v:mp.exp(-1/(1-v*v)),[-1,0,1])
        if not source.endpoints(field.normalization)[0]<=norm<=source.endpoints(field.normalization)[1]:
            raise ArithmeticError('Independent raw integral misses accepted normalization')
        ell=mp.mpf(1)/40;centers=[mp.mpf(1)/5,mp.mpf(1)/2,mp.mpf(4)/5]
        def bump(y,center):
            v=(y-center)/ell
            return mp.exp(-y-1/(1-v*v))/(ell*norm) if abs(v)<1 else mp.mpf(0)
        G=[source.endpoints(v)[1] for v in jets['log_bump_ordinary_derivative_caps']]
        N=field.quiet_threshold
        for center in centers:
            for fraction in ('-.95','-.3','0','.4','.95'):
                y=center+ell*mp.mpf(fraction)
                rows=[mp.diff(lambda v:bump(v,center),y,j) for j in range(5)]
                if any(abs(rows[j])>G[j] for j in range(5)):
                    raise ArithmeticError('Independent ordinary log bump exceeds source cap')
                for muv in ('0.001','1e-20'):
                    m=mp.mpf(muv);alpha=mp.mpf('.5')+m
                    e=-mp.mpf('.9')*R*m/N*rows[0];ey=-mp.mpf('.9')*R*m/N*rows[1]
                    u_y=mp.mpf('.8')*R*mp.sqrt(m)/N*rows[1]
                    f=mp.exp(-alpha*y);denom=mp.exp(-alpha)-m*R*G[0]/N
                    if denom<=0:raise ArithmeticError('Fixture denominator must be positive')
                    a=2*m-2*(ey+alpha*e)/(f+e);b=-2*u_y/(f+e)
                    error=2*m*R*(G[1]+alpha*G[0])/(N*denom)
                    if a<2*m-error or abs(b)>2*mp.sqrt(m)*R*G[1]/(N*denom):
                        raise ArithmeticError('Independent quiet quotient exceeds cap')
                    counts['independent_log_bump_and_quiet_shear_fixtures']+=1
        # Direct continuous integrals test new prefix and shrinking complement
        # caps. These supplement the accepted source FTC/terminal identities.
        for yv in ('0','.19','.51','.8','.824','.825','1'):
            y=mp.mpf(yv)
            for remaining in (False,True):
                W=field._weights_cap(c.mpf(yv),remaining)
                for i,center in enumerate(centers):
                    left,right=center-ell,center+ell
                    a=max(left,y) if remaining else left
                    b=right if remaining else min(right,y)
                    if b<=a:continue
                    def beta(z):
                        v=(z-center)/ell
                        return mp.exp(-1/(1-v*v))/(ell*norm) if abs(v)<1 else mp.mpf(0)
                    for name,p,multiplicity,index in (('mass',0,1,i),('I',mp.mpf('.5'),1,i),
                        ('S',-mp.mpf('.5')-mu,1,i),('Cp',-mp.mpf('1.5')-mu,1,i),
                        ('energy',-1,2,i),('pressure',-2,2,i)):
                        integral=mp.quad(lambda z:beta(z)**multiplicity*mp.exp(p*z),[a,(a+b)/2,b])
                        if integral>source.endpoints(W[name][index])[1]:
                            raise ArithmeticError('Independent partial integral exceeds '+name+' strip cap')
                        counts['independent_partial_weight_integrals']+=1
                    if i in (0,2):
                        j=0 if i==0 else 1
                        value=mp.quad(lambda z:beta(z)*z,[a,(a+b)/2,b])
                        # 0<exp(-mu*r*z)<=1; this larger integral bounds |D|.
                        if value>source.endpoints(W['D'][j])[1]+mp.mpf('1e-80'):
                            raise ArithmeticError('Independent divided-row strip exceeds cap')
                        cross_integral=mp.quad(lambda z:beta(z)**2*mp.exp(-z/2),[a,(a+b)/2,b])
                        if cross_integral>source.endpoints(W['cross'][j])[1]:
                            raise ArithmeticError('Independent cross strip exceeds cap')
                        counts['independent_partial_weight_integrals']+=2
    return counts


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed uniform repair dependency: '+name)
    field=source.CurrentUniformRepairMajorants(require_checked=False)
    if source.parameters.encoded(field.theorem)!=data['exact_uniform_repair_theorem']:
        raise ValueError('Uniform exact-source theorem changed')
    specs=(('repair_threshold',field.repair_threshold,(0,1)),('quiet_threshold',field.quiet_threshold,(0,1)),
        ('original_frequency',10**12,(0,1)),('before_first_bump',field.quiet_threshold,'0'),
        ('inside_last_strip',field.quiet_threshold,('.8','.824')),
        ('near_last_exit',field.quiet_threshold,('.824','.825')),('terminal',field.quiet_threshold,'1'))
    for name,N,y in specs:
        fresh=field.query(N,y)
        if source.parameters.encoded(fresh)!=data['examples'][name]:raise ValueError('Uniform query differs: '+name)
        if fresh['saved_frequency_controls_or_signed_defects_reused']:raise ValueError('Fixed-N coefficients relabelled')
        if N>=field.quiet_threshold and not fresh['quiet_repaired_shear']['at_least_three_halves_mu_certified']:
            raise ArithmeticError('Sufficient quiet frequency loses the source reserve')
    terminal=field.query(field.quiet_threshold,'1')
    if any(v._mpi_!=field.ctx.mpf(0)._mpi_ for v in terminal['same_source_cumulative_defect_absolute_caps'].values()):
        raise ArithmeticError('Terminal implicit closure does not yield zero defect bound')
    whole=field.query(field.quiet_threshold)
    near=field.query(field.quiet_threshold,('.824','.825'))
    if any(source.endpoints(near['remaining_repair_primitive_absolute_caps'][k])[1]>=source.endpoints(whole['remaining_repair_primitive_absolute_caps'][k])[1] for k in whole['remaining_repair_primitive_absolute_caps']):
        raise ArithmeticError('Remaining support strip does not shrink')
    for k in source.OPEN:
        if data.get(k) is not False:raise ValueError('Repair estimate promoted unresolved global gate')
    invalid=((0,(0,1)),(True,(0,1)),(field.repair_threshold-1,(0,1)),(field.quiet_threshold,(-1,0)),
        (field.quiet_threshold,(0,2)),(field.quiet_threshold,mp.inf))
    for N,y in invalid:
        try:field.query(N,y)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid uniform repair query accepted')
    counts=numerical_checks(field)
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_source_Jacobian_beta_jet_and_terminal_identities=len(field.theorem['identities']),
        variable_N_queries_recomputed=len(specs),invalid_queries_rejected=len(invalid),**counts,
        repair_only_sufficient_integer_N=field.repair_threshold,quiet_primitive_sufficient_integer_N=field.quiet_threshold,
        whole_quiet_primitive_margin_at_least_three_halves_mu=True,
        variable_N_implicit_controls_and_partial_terminal_defect_bounds_certified=True,
        fixed_frequency_coefficients_not_reused=True,**{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Uniform repair PASS:',len(field.theorem['identities']),'identities;',counts,flush=True)
    return result


if __name__=='__main__':run()
