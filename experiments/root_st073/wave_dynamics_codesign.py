"""Optimize wave shape against projected complete instantaneous momentum.

Keep the original five flux equalities and growth/norm constraints while
minimizing the physical-volume residual after eliminating the fixed linear
time-derivative/pressure directions. All quadratic wave terms are retained.
The training patch is not a trajectory or an independent PDE certificate.
"""
import argparse
import hashlib
import json
import os
import time
from pathlib import Path

for _thread_variable in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_thread_variable] = '1'

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent


def decode(value):
    value = np.asarray(value, float)
    return value[...,0]+1j*value[...,1]


def pack(value):
    value = np.asarray(value)
    return np.stack((value.real,value.imag),axis=-1).tolist()


def real_matrix(matrix):
    return np.block([[matrix.real,-matrix.imag],[matrix.imag,matrix.real]])


class ProjectedMomentum:
    def __init__(self, cache, transform, rcond):
        self.A = cache['A']@transform
        self.B = cache['B']@transform
        self.L = cache['L']@transform
        self.R0 = cache['R0']
        self.weights = cache['weights']
        self.sqrtw = np.repeat(np.sqrt(self.weights),3)
        keep = cache['tangent_s'] > rcond*cache['tangent_s'][0]
        self.Q = cache['tangent_U'][:,keep]
        self.s = cache['tangent_s'][keep]
        self.Vh = cache['tangent_Vh'][keep]
        self.scales = cache['tangent_scales']
        self.rank = int(keep.sum())
        self.normalization = 1.0
        self.last_x = None
        self.last_result = None

    def evaluate(self, x):
        if self.last_x is not None and np.array_equal(x,self.last_x):
            return self.last_result
        wave = self.A@x
        gradient = self.B@x
        residual = self.R0+self.L@x+np.einsum('nij,nj->ni',gradient,wave)
        weighted = residual.reshape(-1)*self.sqrtw
        projection = self.Q.T@weighted
        remainder = weighted-self.Q@projection
        derivative = self.L+np.einsum('nijq,nj->niq',self.B,wave)
        derivative += np.einsum('nij,njq->niq',gradient,self.A)
        jac = derivative.reshape(-1,len(x))*self.sqrtw[:,None]
        value = float(remainder@remainder)/self.normalization
        # The projector is fixed and self-adjoint; projected remainder suffices.
        grad = 2*(jac.T@remainder)/self.normalization
        control = -(self.Vh.T@(projection/self.s))/self.scales
        corrected = (remainder/self.sqrtw).reshape(-1,3)
        self.last_result = (value,grad,dict(
            frozen_volume_L2=float(np.linalg.norm(weighted)),
            corrected_volume_L2=float(np.linalg.norm(remainder)),
            corrected_max=float(np.linalg.norm(corrected,axis=1).max()),
            tangent_coefficients=control.tolist()))
        self.last_x = np.array(x,copy=True)
        return self.last_result


def run(maxiter=80,rcond=1e-6,method='retracted'):
    started = time.perf_counter()
    source_path = ROOT/'wave_stress_growth_codesign.json'
    raw = source_path.read_bytes()
    source = json.loads(raw)
    problem = source['problem']
    initial_white = decode(source['selected']['coefficients_whitened'])
    variable_scale = float(np.linalg.norm(initial_white))
    transform = real_matrix(decode(problem['whitening']))*variable_scale
    x0 = np.r_[initial_white.real,initial_white.imag]/variable_scale
    cache_path = ROOT/'wave_momentum_projection.npz'
    with np.load(cache_path,allow_pickle=False) as loaded:
        cache = {key:loaded[key] for key in loaded.files}
    objective = ProjectedMomentum(cache,transform,rcond)
    objective.normalization = objective.evaluate(x0)[0]
    objective.last_x = None
    target = np.asarray(problem['target']).reshape(-1)
    flux_scales = np.maximum(abs(target),1.0)
    forms = decode(problem['flux_forms_whitened']).reshape(-1,len(initial_white),len(initial_white))
    forms = np.array([real_matrix(h)*variable_scale**2 for h in forms])
    mass = real_matrix(decode(problem['mass_whitened']))*variable_scale**2
    growth = real_matrix(decode(problem['growth_whitened']))*variable_scale**2
    floor = float(source['energy_threshold_lambda'])
    margin_matrix = growth-floor*mass
    margin_scale = max(abs(float(x0@growth@x0)),1.0)
    bound = float(problem['coefficient_bound'])/variable_scale

    def equalities(x):
        return (np.einsum('i,kij,j->k',x,forms,x)-target)/flux_scales

    def equalities_jac(x):
        return 2*np.einsum('kij,j->ki',forms,x)/flux_scales[:,None]

    def inequalities(x):
        return np.array([x@margin_matrix@x/margin_scale,1-x@x/bound**2])

    def inequalities_jac(x):
        return np.array([2*margin_matrix@x/margin_scale,-2*x/bound**2])

    def summary(x):
        value,_,metrics = objective.evaluate(x)
        y = variable_scale*(x[:len(initial_white)]+1j*x[len(initial_white):])
        wave_velocity = objective.A@x
        wave_gradient = objective.B@x
        energy = .5*float(np.sum(cache['weights'][:,None]*wave_velocity**2))
        production = -float(np.einsum('n,ni,nij,nj->',cache['weights'],
            wave_velocity,cache['mean_gradient'],wave_velocity))
        dissipation = float(source['viscosity'])*float(np.sum(
            cache['weights'][:,None,None]*wave_gradient**2))
        control = np.array(metrics['tangent_coefficients'])
        derivative = np.r_[control[36:90:2],control[37:90:2]]
        wave_time_derivative = cache['A']@derivative
        fitted_energy_rate = float(np.sum(cache['weights'][:,None]*
                                         wave_velocity*wave_time_derivative))
        return dict(objective=value,**metrics,
            flux_max_normalized=float(max(abs(equalities(x)))),
            growth_lambda=float(x@growth@x/(x@mass@x)),
            inequalities=inequalities(x).tolist(),
            training_wave_energy=energy,training_shear_production=production,
            training_viscous_dissipation=dissipation,
            training_growth_lambda=(production-dissipation)/(2*energy),
            fitted_wave_energy_rate=fitted_energy_rate,
            coefficients_whitened=pack(y),
            coefficients_original=pack(decode(problem['whitening'])@y))

    report = dict(status='optimizing',accepted=False,pde_validated=False,
        scale_recursion_established=False,constraints_maintained=False,
        source_sha256=hashlib.sha256(raw).hexdigest(),source=source_path.name,
        projection_source='wave_momentum_projection.json',
        projection_cache_sha256=hashlib.sha256(cache_path.read_bytes()).hexdigest(),
        rcond=rcond,tangent_rank=objective.rank,variable_scale=variable_scale,
        growth_floor=floor,maxiter=maxiter,method=method,initial=summary(x0),iterations=[],
        scope='Training-only complete momentum shape fit with five-node flux equalities and original growth/norm bounds; tangent mode0 changes compatibility. No independent PDE or recursive acceptance.')
    best = x0.copy()
    best_value = report['initial']['objective']
    output = ROOT/'wave_dynamics_codesign.json'

    def save():
        report['selected'] = summary(best)
        report['elapsed_seconds'] = time.perf_counter()-started
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

    def callback(x):
        nonlocal best,best_value
        value = objective.evaluate(x)[0]
        eq = float(max(abs(equalities(x))))
        ineq = float(min(inequalities(x)))
        if eq <= 1e-5 and ineq >= -1e-7 and value < best_value:
            best,best_value = x.copy(),value
        report['iterations'].append(dict(iteration=len(report['iterations'])+1,
            objective=value,flux_max_normalized=eq,inequality_min=ineq))
        if len(report['iterations']) % 5 == 0:
            save()
            print(json.dumps(report['iterations'][-1]),flush=True)

    # Directional finite-difference check before optimization, not a PDE test.
    direction = np.cos(np.arange(len(x0))+0.4)
    direction /= np.linalg.norm(direction)
    step = 1e-5
    analytic = objective.evaluate(x0)[1]@direction
    finite = (objective.evaluate(x0+step*direction)[0]-objective.evaluate(x0-step*direction)[0])/(2*step)
    error = abs(analytic-finite)/max(abs(analytic),abs(finite),1e-12)
    report['objective_gradient_check'] = dict(analytic=float(analytic),finite_difference=float(finite),relative_error=float(error))
    if error > 1e-4:
        raise ValueError('Projected objective gradient check failed')
    save()
    if method == 'slsqp':
        fit = minimize(lambda x:objective.evaluate(x)[:2],x0,jac=True,method='SLSQP',
            constraints=[dict(type='eq',fun=equalities,jac=equalities_jac),
                         dict(type='ineq',fun=inequalities,jac=inequalities_jac)],
            callback=callback,options=dict(maxiter=maxiter,ftol=1e-9,disp=False))
        callback(fit.x)
        report.update(optimizer_success=bool(fit.success),optimizer_message=str(fit.message),
                      optimizer_status=int(fit.status),final_trial=summary(fit.x))
    else:
        # Follow the feasible growth-floor surface. This deliberately searches
        # a subset of the inequality-feasible set; no global optimum is claimed.
        def surface(x):
            return np.r_[equalities(x),inequalities(x)[0]]

        def surface_jac(x):
            return np.vstack((equalities_jac(x),inequalities_jac(x)[0]))

        def retract(x):
            for _ in range(20):
                residual = surface(x)
                if max(abs(residual)) < 1e-8:
                    return x if np.linalg.norm(x) <= bound*(1+1e-9) else None
                step = np.linalg.lstsq(surface_jac(x),-residual,rcond=1e-10)[0]
                old_norm = np.linalg.norm(residual)
                for damping in (1.,.5,.25,.125,.0625):
                    trial = x+damping*step
                    if np.linalg.norm(surface(trial)) < old_norm:
                        x = trial
                        break
                else:
                    return None
            return None

        x = x0.copy()
        message = 'iteration limit; feasible descent points retained'
        report['constraint_surface_singular_values'] = np.linalg.svd(surface_jac(x),compute_uv=False).tolist()
        for _ in range(maxiter):
            value,gradient,_ = objective.evaluate(x)
            jac = surface_jac(x)
            tangent_gradient = gradient-jac.T@np.linalg.lstsq(jac.T,gradient,rcond=1e-10)[0]
            norm = np.linalg.norm(tangent_gradient)
            if norm < 1e-10:
                message = 'small projected gradient on restricted constraint surface'
                break
            direction = -tangent_gradient/norm
            for length in (.05,.025,.0125,.00625,.003125,.0015625,.00078125,.000390625):
                candidate = retract(x+length*direction)
                if candidate is not None and objective.evaluate(candidate)[0] < value-1e-8*length*norm:
                    x = candidate
                    callback(x)
                    break
            else:
                message = 'no accepted feasible descent step at the tested step lengths'
                break
        report.update(optimizer_success=False,optimizer_message=message,
                      optimizer_status=None,final_trial=summary(x),
                      restricted_search='Retraction onto ten flux equalities plus growth lambda=1000 surface; norm bound checked at every accepted step.')
    report['status'] = 'completed'
    save()
    print(json.dumps({key:report['selected'][key] for key in
          ('objective','corrected_volume_L2','flux_max_normalized','growth_lambda')}),flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--maxiter',type=int,default=80)
    parser.add_argument('--rcond',type=float,default=1e-6)
    parser.add_argument('--method',choices=('retracted','slsqp'),default='retracted')
    args = parser.parse_args()
    run(args.maxiter,args.rcond,args.method)
