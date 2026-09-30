"""Physical-chart and stress receipt for the actual initial exit field.

This is an unlocalized collar, not a global finite-energy NS solution.
"""
import json
from pathlib import Path
import mpmath as mp


class LocalCoreExitField:
    """One physical (x,y,z,t) callable on the core and initial exit.

    Outside its radial construction domain it raises rather than extending
    by zero. The source pressure approximation remains explicit.
    """
    def __init__(self,tangent,*,nu='.01',T=0):
        self.tangent=tangent; self.core=tangent.core
        self.precision=tangent.precision
        with mp.workdps(self.precision):
            self.nu=mp.mpf(nu); self.T=mp.mpf(T)
            if self.nu<=0:raise ValueError('nu must be positive')

    def evaluate(self,x,y,z,t):
        with mp.workdps(self.precision):
            x,y,z,t=map(mp.mpf,(x,y,z,t)); tau=self.T-t
            if tau<=0:raise ValueError('Require t<T')
            dt=self.core.delta; exponent=(1-dt)/2
            scaled_z=z/(mp.sqrt(self.nu)*tau**exponent)
            if scaled_z==0:
                Z=mp.mpf(0)
            else:
                # Delta=0 gives this exact inverse; Newton resolves the
                # source's nonzero delta without dropping it from the map.
                initial=scaled_z/mp.sqrt(1+scaled_z*scaled_z)
                equation=lambda w:w/(1-w*w)**exponent-scaled_z
                derivative=lambda w:(1-dt*w*w)/(1-w*w)**(exponent+1)
                Z=mp.findroot(equation,initial,df=derivative,solver='newton')
            if not abs(Z)<1:raise ValueError('Physical chart inversion reached axial boundary')
            q=tau/(1-Z*Z); R=(x*x+y*y)/(2*self.nu*q)
            if R<=self.tangent.r:
                sample={**self.core.evaluate(R,Z),'R':R}; region='core'
            else:
                sample=self.tangent.evaluate(mp.log(R/self.tangent.r),Z)
                region=sample.get('region','initial_exit')
            chart=physical_chart(sample,Z,mp.log(q),delta=dt,nu=self.nu,
                                 phi=mp.atan2(y,x),precision=self.precision)
            return {**chart,'region':region,'R':R,'Z':Z,
                    'scope':'Unlocalized core and constructed exit only; no complete outer matching'}

    def velocity(self,x,y,z,t):
        return self.evaluate(x,y,z,t)['uvw']


def physical_chart(sample,Z,logq,*,delta='1e-32',nu='.01',phi=0,precision=160):
    with mp.workdps(precision):
        z=mp.mpf(Z); q=mp.exp(mp.mpf(logq)); dt=mp.mpf(delta)
        viscosity=mp.mpf(nu); angle=mp.mpf(phi); R=sample['R']
        radius=mp.sqrt(2*viscosity*q*R)
        axial=mp.sqrt(viscosity)*q**((1-dt)/2)*z
        ur=mp.sqrt(viscosity/q)*sample['Ur']
        scale=mp.sqrt(viscosity)*q**(-(1+dt)/2)
        ut=scale*mp.sqrt(2*R)*sample['F']; uz=scale*sample['Uz']
        return {'xyz':(radius*mp.cos(angle),radius*mp.sin(angle),axial),
                'cylindrical_velocity':(ur,ut,uz),
                'uvw':(ur*mp.cos(angle)-ut*mp.sin(angle),
                       ur*mp.sin(angle)+ut*mp.cos(angle),uz),
                'tau':q*(1-z*z),'pressure':viscosity*q**(-1-dt)*sample['P']}


def cone_receipt(sample):
    """Section 3.21/3.23 pointwise tests; no uniform certification."""
    F=sample['F']; s=sample['stress']
    st=s['S_theta']/F; sz=s['S_z']/F
    tt=s['T_theta']/F; tz=s['T_z']/F
    if F<=0 or st>=0:
        return {'prerequisites_passed':False,'relaxed_passed':False,'admissible_passed':False}
    kappa=-(st*st+sz*sz)/st
    dot=tt*st+tz*sz; cross=-tt*sz+tz*st
    direction=dot<0
    if kappa>2:
        margin=2*dot*dot-(kappa-2)*cross*cross
        relaxed=direction and margin>0
        branch='kappa>2'
    else:
        # Divide by positive -S_theta/F, retaining tiny collar factors.
        margin=-dot/(-st)-(2-kappa)
        relaxed=direction and margin>0
        branch='kappa<=2'
    def signed_log(x):
        return {'sign':int(mp.sign(x)),'log_abs':mp.nstr(mp.log(abs(x)),35) if x else None}
    return {'prerequisites_passed':True,'kappa':signed_log(kappa),
            'T_dot_S_over_F_squared':signed_log(dot),'branch':branch,
            'branch_margin':signed_log(margin),'relaxed_passed':bool(relaxed),
            'admissible_passed':bool(relaxed and kappa>2)}


def run():
    from lei_ren_part1_paper_exit_comparison import build_comparison
    from lei_ren_part1_paper_exit_bridge import ExitBridge
    from lei_ren_part1_paper_exit_tangents import ExitTangents
    comparison=build_comparison()
    tangent=ExitTangents(comparison,steps=16,derivative_step='1e-45')
    fine=ExitTangents(comparison,steps=32,derivative_step='1e-45')
    half_step=ExitTangents(comparison,steps=32,derivative_step='5e-46')
    independent_velocity=ExitBridge(comparison,steps=32)
    rows=[]
    with mp.workdps(comparison.precision):
        def num(x):return mp.nstr(x,35)
        def signed_log(x):return {'sign':int(mp.sign(x)),'log_abs':num(mp.log(abs(x))) if x else None}
        Z=mp.mpf('.3')
        start=fine.evaluate(0,Z)
        initial_core=comparison.core.evaluate(comparison.R_a,Z)
        initial_errors={}
        for family,target in (('moments','moments'),('momentsZ','moments_Z')):
            initial_errors[family]={}
            for key,value in start[family].items():
                reference=initial_core[target][key]
                scale=max(abs(reference),abs(value))
                error=abs(value-reference)/scale if scale else mp.mpf(0)
                initial_errors[family][key]=num(error)
                if error>mp.mpf('1e-100'):
                    raise ArithmeticError('Initial core moment/jet normalization mismatch: '+family+'/'+key)
        for y in (mp.mpf(0),mp.mpf('.0025'),mp.mpf('.01')):
            a=tangent.evaluate(y,Z); b=fine.evaluate(y,Z)
            chart=physical_chart(b,Z,'-4',precision=comparison.precision)
            rows.append({'y':num(y),'Z':num(Z),'Ur':signed_log(b['Ur']),
                'strict_annulus_scope':'Ra boundary excluded from open-annulus strict condition' if y==0 else 'Interior sample',
                'Uz_Z':num(b['UZ']),'log_F_Z_over_F':num(b['FZ']/b['F']),
                'Ur_relative_refinement':num(abs(a['Ur']-b['Ur'])/abs(b['Ur'])),
                'Uz_Z_difference':num(abs(a['UZ']-b['UZ'])),
                'stress':{k:signed_log(v) for k,v in b['stress'].items()},
                'cone':cone_receipt(b),
                'physical_chart':{k:[signed_log(v) for v in chart[k]] for k in ('xyz','uvw','cylindrical_velocity')}})
        midpoint=fine.evaluate(mp.mpf('.0025'),Z)
        midpoint_half=half_step.evaluate(mp.mpf('.0025'),Z)
        driver_refinement={'y':'.0025','steps':32,'Z_steps':['1e-45','5e-46'],
            'Ur_relative_difference':num(abs(midpoint['Ur']-midpoint_half['Ur'])/abs(midpoint_half['Ur'])),
            'Uz_Z_difference':num(abs(midpoint['UZ']-midpoint_half['UZ'])),
            'F_Z_over_F_difference':num(abs(midpoint['FZ']/midpoint['F']-midpoint_half['FZ']/midpoint_half['F']))}
        # Independently differentiate Ur radially and Uz axially; this does
        # not substitute the analytic radial continuity identity for a check.
        y=mp.mpf('.0025'); h=mp.mpf('1e-5'); hz=mp.mpf('1e-44')
        radial=[fine.evaluate(y+k*h,Z)['Ur'] for k in (-2,-1,1,2)]
        ur_y=(radial[0]-8*radial[1]+8*radial[2]-radial[3])/(12*h)
        axial=[independent_velocity.evaluate(y,Z+k*hz)['Uz'] for k in (-2,-1,1,2)]
        uz_Z=(axial[0]-8*axial[1]+8*axial[2]-axial[3])/(12*hz)
        a=fine.evaluate(y,Z); R=a['R']; dt=comparison.delta
        radial_div=mp.sqrt(2/R)*ur_y+a['Ur']/mp.sqrt(2*R)
        axial_div=((1-Z*Z)*uz_Z-2*Z*a['u_y']-(1+dt)*Z*a['Uz'])/(1-dt*Z*Z)
        divergence={'y':num(y),'radial_difference_step':num(h),'Z_difference_step':num(hz),
            'q_times_divergence':num(radial_div+axial_div),
            'relative_cancellation_error':num(abs(radial_div+axial_div)/(abs(radial_div)+abs(axial_div))),
            'scope':'One finite-difference sample of actual collar; not whole-domain validation'}
    receipt={'precision':comparison.precision,'initial_moment_relative_errors':initial_errors,
             'rows':rows,'divergence':divergence,
             'driver_step_refinement':driver_refinement,
             'refinement_scope':'16/32 ODE steps at fixed Z step; driver Z step separately halved at midpoint',
             'driver_Z_derivatives':'MP fourth-order centered differences, not analytic bounds',
             'global_cone_certified':False,'global_finite_energy_certified':False,
             'outer_matching_complete':False,'scale_recursion_established':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'rows':len(rows),'divergence':divergence,'cones':[r['cone'] for r in rows]},indent=2))
    return receipt


def run_callable_checks():
    from lei_ren_part1_paper_exit_comparison import build_comparison
    from lei_ren_part1_paper_exit_tangents import ExitTangents
    comparison=build_comparison()
    tangent=ExitTangents(comparison,steps=8)
    field=LocalCoreExitField(tangent)
    rows=[]
    with mp.workdps(comparison.precision):
        Z=mp.mpf('.3')
        core_R=2/comparison.Lambda
        samples=[('core',{**comparison.core.evaluate(core_R,Z),'R':core_R}),
                 ('initial_exit',tangent.evaluate(mp.mpf('.0025'),Z))]
        for region,sample in samples:
            for logq in ('-4','-8'):
                expected=physical_chart(sample,Z,logq,delta=comparison.delta,phi='.4')
                actual=field.evaluate(*expected['xyz'],-expected['tau'])
                error=max(abs(a-b)/max(abs(a),abs(b)) if a or b else mp.mpf(0)
                          for a,b in zip(expected['cylindrical_velocity'],actual['cylindrical_velocity']))
                if actual['region']!=region or error>mp.mpf('1e-100'):
                    raise ArithmeticError('Physical callable roundtrip failed')
                rows.append({'region':region,'logq':logq,
                    'R_relative_error':mp.nstr(abs(actual['R']/sample['R']-1),35),
                    'Z_error':mp.nstr(abs(actual['Z']-Z),35),
                    'cylindrical_velocity_relative_error':mp.nstr(error,35)})
    receipt={'rows':rows,'scope':'Physical-map roundtrip at two times in local core/exit, not dynamic recursion',
             'outer_matching_complete':False,'global_finite_energy_certified':False}
    Path(__file__).with_name('lei_ren_part1_paper_exit_callable.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))
    return receipt


if __name__=='__main__':
    import sys
    run_callable_checks() if '--callable' in sys.argv else run()
