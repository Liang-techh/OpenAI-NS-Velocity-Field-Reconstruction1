"""Analytic radial velocity jets from the installed continuous definitions.

Complements the existing analytic Z jets for Section 3 stress evaluation.
No coefficients, moments, pressure datum or terminal residual are changed.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import _mp,build_joined_field
from lei_ren_part1_paper_axial_correction import signed_log


def velocity_radial_jets(profile,logR,Z):
    with mp.workdps(max(profile.precision,profile.jet_precision)):
        s=profile.schedule;z=_mp(Z)
        values=profile.values_with_jets(logR,z)
        theta_y=values['Utheta']*values['logarithmic_slope']
        start=_mp(str(profile.offset(logR,s.logR_p)))
        end=_mp(str(profile.offset(logR,s.logR_v)))
        if start<0:
            axial_y=profile.incoming_provider.Uz_y(profile._incoming_y(logR),z)
            method='shared_incoming_analytic_cutoff'
        elif end>0:
            axial_y=mp.mpf(0);method='post_axial_support'
        else:
            runtime=profile.runtime(float(z));base=s.at_log_radius(logR,z)
            amplitude=mp.exp(_mp(base['log_angular_amplitude']))
            slope=_mp(base['logarithmic_slope'])
            if mp.mpf('-3.15')<=end<=mp.mpf('-.85'):
                beta=beta_y=mp.mpf(0)
                for c,center in zip(runtime.c,(-3,-1)):
                    jet=runtime.basis.values(end-center)
                    beta+=c*jet['beta'];beta_y+=c*jet['beta_s']
                axial_y=amplitude*(slope*beta+beta_y)
                method='shared_end_bump_product_jet'
            else:
                pulse=runtime.pulse.value_jet(runtime.mu*start)
                axial_y=amplitude*runtime.a*(slope*pulse['value']+runtime.mu*pulse['derivative'])
                method='shared_pulse_product_jet'
        return dict(values,Utheta_y=theta_y,Uz_y=axial_y,
            radial_velocity_derivatives_analytic=True,axial_radial_method=method,
            radial_coordinate='logR',input_uncertainty_enclosed=False,
            stress_cone_certified=False,finite_energy_certified=False,scale_recursion_certified=False)


def run():
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    folder=Path(__file__).parent
    print('building shared field for analytic radial velocity jets',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;s=p.schedule;samples=[]
    with mp.workdps(field.precision):
        z=mp.mpf('.3');step=mp.mpf('1e-5')
        points=[('incoming',p.log_at(s.logRref,'1.3')),
            ('pulse',p.log_at(s.logR_p,mp.nstr(5/_mp(s.mu),field.precision))),
            ('end_bump',p.log_at(s.logR_v,'-2.96')),
            ('heat',p.log_at(s.logR_tail,'4'))]
        for label,logR in points:
            row=velocity_radial_jets(p,logR,z)
            neighbors={i:p.values(p.log_at(logR,i*step),z) for i in (-2,-1,1,2)}
            errors=[]
            for name in ('Utheta','Uz'):
                derivative=sum(w*neighbors[i][name] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
                target=row[name+'_y']
                error=abs(derivative/target-1) if target else abs(derivative)
                errors.append(error)
            if max(errors)>mp.mpf('1e-8'):raise ArithmeticError('Analytic radial jet differs from actual velocity')
            samples.append(dict(label=label,method=row['axial_radial_method'],
                radial_derivative_relative_errors=[mp.nstr(e,40) for e in errors],
                Uz_y=signed_log(row['Uz_y'],50)))
        report=dict(samples=samples,analytic_logR_velocity_derivatives=True,
            shared_runtime_atoms=True,input_uncertainty_enclosed=False,
            stress_cone_certified=False,complete_five_moments=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'samples':[{k:v for k,v in r.items() if k!='Uz_y'} for r in samples]}),flush=True)
    return report


if __name__=='__main__':run()
