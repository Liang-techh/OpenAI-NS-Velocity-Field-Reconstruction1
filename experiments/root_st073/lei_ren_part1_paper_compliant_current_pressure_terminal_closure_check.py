"""Independent original-generator fixture and current pressure admission."""
import ast,json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_pressure_terminal_closure import (
    CurrentPressureTerminalClosure,original_preheat_partition,ExactOriginalPreheatPressureOperator,
    verify_exact_integral_witness,live_native_density_translations,live_native_raw_angular_translations,
    original_pressure_function_identification,HERE,NAME,RECEIPT,GATES,OPEN,VIEWS,
    sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def independent_original_master_fixture():
    """Original global formula versus every segmented density, with finite logs.

    Moderate arbitrary positive lengths make all identities resolvable.
    These values test source algebra, not final paper parameter admissibility.
    """
    with mp.workdps(85):
        mu=mp.mpf('.08');delta=mp.mpf('.01');eps=mp.mpf('.003')
        yd=mp.mpf(4);Tw=mp.mpf(3);L=mp.mpf('2.5');Ts=mp.mpf('2.2');W=mp.mpf('1.7')
        bp=mp.mpf('.5')+mu;r=1-mu;k=1-delta/2;bh=(1+delta)/2
        yp=yd+1+Tw;yv=yp+13/mu;yf=yv+100;yr=yf+L;ys=yr+1;yq=ys+Ts;yt=yq+1;tail=yt+W
        def sig(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/(x*x)-1/((1-x)*(1-x))))
        cache={}
        def J(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return x-mp.mpf('.5')
            if x not in cache:cache[x]=mp.quad(sig,[0,x/2,x])
            return cache[x]
        def phi(t):return mp.exp(-4/(3-t)**2) if t<3 else mp.mpf(0)
        def logA(y):return y/10-mp.mpf('.6')*J(y)-mu*J(y-yd)-r*J(y-yr)+k*J(y-yq)
        def original(y,z):
            if y<tail:
                q=1+z*z
                return mp.exp(logA(y))/q*(q/2)**sig((y-yv)/100)
            t=y-tail;bracket=(1-sig(t))*(1-eps)+sig(t)*(1-eps*phi(t))
            return mp.exp(logA(tail))/(2*(1-eps))*mp.exp(-bh*t)*bracket
        starts=dict(reference_extension=0,slope_transition_ref=0,axial_turnoff=1,
            slope_transition_mu=yd,power_buffer=yd+1,pulse_reserved=yp,z_flatten=yv,
            power_buffer_rel=yf,steep_transition_in=yr,steep_power=ys,
            steep_transition_out=yq,waiting=yt,heat_collar=tail,exterior_power_tail=tail)
        p=original_preheat_partition();sub={p['mu']:mu,p['delta']:delta,p['epsilon']:eps,
            p['yd']:yd,p['Tw']:Tw,p['L']:L,p['Ts']:Ts,p['W']:W}
        rows=0;nonzero_tail=0
        for name,shape in p['shapes'].items():
            fn=s.lambdify((p['z'],p['t']),shape.subs(sub),modules=[{'J_sigma':J,'sigma':sig,'phi':phi},'mpmath'])
            if name=='reference_extension':points=map(mp.mpf,('-.23','-2.7'))
            elif name=='exterior_power_tail':points=map(mp.mpf,('3.17','4.8'))
            else:
                length=mp.mpf(str(p['domains'][name][1].subs(sub) if isinstance(p['domains'][name][1],s.Basic) else p['domains'][name][1]))
                points=(length*mp.mpf('.23'),length*mp.mpf('.77'))
            for t in points:
                for z in map(mp.mpf,('0','.43','.89')):
                    a=original(starts[name]+t,z);b=fn(z,t)
                    if a<=0 or b<=0 or abs((a/b)**2-1)>mp.mpf('1e-55'):
                        raise ArithmeticError('Original global and native segmented pressure density differ: '+name)
                    rows+=1
                    if name=='exterior_power_tail':nonzero_tail+=1
        # A wrong c_inf factor would fail even where H=1 and sigma=1.
        if nonzero_tail!=6:raise ValueError('Infinite raw exterior fixture missing')
        return dict(original_global_to_segmented_pressure_density_rows=rows,
            fourteen_original_stages_checked=True,original_nonzero_exterior_rows=nonzero_tail,
            source_algebra_fixture_not_parameter_admission=True,passed=True)


def independent_integral_dependency_checks():
    """Unilateral mutations fail; neither integral is inferred from the other."""
    def reject(label,mutate):
        witness=ExactOriginalPreheatPressureOperator();mutate(witness)
        try:verify_exact_integral_witness(witness)
        except (ArithmeticError,ValueError):return label
        raise ArithmeticError('Unilateral exact source mutation accepted: '+label)
    rejected=[]
    rejected.append(reject('native_flatten_density_factor_two',
        lambda w:w.native_densities.__setitem__('z_flatten',2*w.native_densities['z_flatten'])))
    rejected.append(reject('flatten_variable_exponent_removed',
        lambda w:w.native_densities.__setitem__('z_flatten',w.native_densities['z_flatten'].subs(w.partition['sigma'](w.partition['t']/100),0))))
    rejected.append(reject('independent_original_P0_sign_flipped',
        lambda w:setattr(w,'datum',lambda:w.integral(False))))
    eta=s.Symbol('independent_P0_perturbation',positive=True)
    rejected.append(reject('independent_original_P0_perturbed',
        lambda w:setattr(w,'datum',lambda:-w.integral(False)+eta)))
    rejected.append(reject('axis_flatten_exit_factor_two_removed',
        lambda w:w.raw_angular_chain.__setitem__('Xf',w.raw_angular_chain['Xf']/2)))
    rejected.append(reject('raw_native_Xv_perturbed',
        lambda w:w.raw_angular_chain.__setitem__('Xv',w.raw_angular_chain['Xv']+eta)))
    rejected.append(reject('raw_waiting_root_perturbed',
        lambda w:setattr(w,'raw_waiting_root',w.raw_waiting_root+eta)))
    # Enforce dependency direction in the two exact source functions.
    source=ast.parse((HERE/(NAME[:-5]+'.py')).read_text(encoding='utf8'))
    cls=next(n for n in source.body if isinstance(n,ast.ClassDef) and n.name=='ExactOriginalPreheatPressureOperator')
    datum=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='datum')
    value=next(n.value for n in datum.body if isinstance(n,ast.Return))
    if ast.dump(value)!=ast.dump(ast.parse('-self.integral(native=False)',mode='eval').body):
        raise ValueError('Original P0 depends on native/raw total instead of original generator')
    for method in ('stage_integral','integral'):
        fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
        names={ast.unparse(n) for n in ast.walk(fn) if isinstance(n,(ast.Name,ast.Attribute))}
        if names & {'P0','self.datum','normalized_jets','m2','m0','self.m2','self.m0'}:
            raise ValueError('Raw exact integral depends on original pressure or interval mass values')
    return dict(unilateral_mutations_rejected=rejected,mutation_count=len(rejected),
        original_P0_uses_only_original_generator=True,native_integral_does_not_use_P0_or_interval_masses=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPressureTerminalClosure(require_checked=False)
    proof=original_pressure_function_identification(field.balance,field.raw)
    if encode(pack(proof))!=raw['original_pressure_function_identification'] or not all(proof['identities'].values()):
        raise ValueError('Original/native pressure function proof changed')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Pressure terminal producer claims admission')
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current pressure source/family/datum differs')
    fixture=independent_original_master_fixture();dependencies=independent_integral_dependency_checks()
    pressure_rows=0;zero_rows=0
    for name,(chart,Z,t) in VIEWS.items():
        packet=field.evaluate(chart,Z,t)
        if encode(pack(packet))!=raw['current_pressure_terminal_views'][name]:raise ValueError('Current pressure packet changed')
        if any(packet[k] for k in OPEN) or not packet['original_axis_pressure_not_replaced_or_tail_patched']:
            raise ValueError('Pressure terminal global scope or datum changed')
        for key in ('actual_Dtheta_Taylor','actual_Cp_Taylor','source_proved_original_P0_plus_native_raw_integral_Taylor'):
            if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in packet[key].coefficients):
                raise ValueError('Current pressure/axis function identity is not zero through5')
            zero_rows+=6
        original=field.angular.terminal_constants(Z)['pressure_infinity']
        if encode(pack(original))!=encode(pack(packet['original_forward_Cp_enclosure_retained'])):
            raise ValueError('Original propagated Cp diagnostic changed')
        # Zero is already proved from source functions; containment only
        # checks consistency of the original directed forward enclosures.
        for key in ('original_forward_Cp_enclosure_retained','original_P0_plus_native_raw_integral_enclosure'):
            if any(not endpoints(v)[0]<=0<=endpoints(v)[1] for v in packet[key].coefficients):
                raise ArithmeticError('Source identity contradicts original pressure enclosures')
        pressure_rows+=len(packet['stable_current_absolute_pressure_mixed4'])
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,original_pressure_function_identification=proof,
        independent_original_master_fixture=fixture,pressure_mixed4_rows=pressure_rows,
        independent_integral_dependency_checks=dependencies,
        source_zero_axial5_coefficients=zero_rows,original_forward_pressure_diagnostic_retained=True,
        all_passed=True,input_hashes=hashes,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original preheat datum/native pressure function identification and current Cp=0 PASS',flush=True)
    return result


if __name__=='__main__':run()
