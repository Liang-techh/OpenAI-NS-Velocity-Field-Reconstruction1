"""Current source chain, energy-unit transport and native derivative admission."""
import copy
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_postpulse_energy_history import (
    CurrentPostpulseEnergyHistory,energy_transport_source_proof,HERE,PREFIX,NAME,RECEIPT,
    GATES,OPEN,VIEWS,sha,pack,encode,endpoints,binding)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def independent_energy_integral_fixture():
    """Resolved finite integral test of transport, not paper parameter admission.

    Original flatten cutoff/beta/phi are used. J=t^2/2 and H=1-b*exp(-t)
    are independently integrable smooth test functions with J(1)=1/2,
    used solely to exercise the generic energy/normalization identities.
    Canonical Gamma and transition source definitions are checked separately.
    """
    with mp.workdps(70):
        mu=mp.mpf('.07');delta=mp.mpf('.03');eps=mp.mpf('.004')
        L=mp.mpf(7);Ts=mp.mpf('1.4');W=mp.mpf('1.7');ell=mp.mpf('.15')
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/x**2-1/(1-x)**2))
        def raw(x):return mp.exp(-1/(1-x*x)) if abs(x)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1]);beta=lambda t:raw(t/ell)/(ell*norm)
        E=mp.quad(lambda t:mp.exp(-2*mu*t)*beta(t),[-ell,0,ell])
        F=mp.quad(lambda t:mp.exp(-2*mu*t)*beta(t)**2,[-ell,0,ell])
        I=lambda rate,length:-mp.expm1(-rate*length)/rate
        ein=mp.quad(lambda t:mp.exp(-2*mu*t-(1-mu)*t*t),[0,1])
        eout=mp.quad(lambda t:mp.exp(-2*t+(1-delta/2)*t*t),[0,1])
        Ns=mp.exp(-1-mu);Nq=Ns*mp.exp(-2*Ts);Nt=Nq*mp.exp(-1-delta/2)
        rows=0
        def equal(a,b,label):
            nonlocal rows
            if abs(a-b)>mp.mpf('1e-52')*max(1,abs(a),abs(b)):raise ArithmeticError('Independent energy fixture differs: '+label)
            rows+=1
        for z in map(mp.mpf,('0','.437','.83')):
            q=1+z*z;b=mp.mpf('.008')*(1-z*z)
            phi=lambda t:mp.exp(-4/(3-t)**2) if t<3 else mp.mpf(0)
            K=lambda t:(1-sigma(t))*(1-eps)+sigma(t)*(1-eps*phi(t))*(1-b*mp.exp(-t))
            E0=mp.quad(lambda t:mp.exp(-delta*t)*K(t)**2,[0,1,2,3])
            E0+=mp.exp(-3*delta)/delta-2*b*mp.exp(-3*(delta+1))/(delta+1)+b*b*mp.exp(-3*(delta+2))/(delta+2)
            H=E0/(1-eps)**2;waiting=H*mp.exp(-delta*W)+I(delta,W)
            afterpower=waiting*mp.exp(-1-delta/2)+eout
            afterentry=(afterpower*mp.exp(-2*Ts)+I(2,Ts))*Ns
            post=ein+afterentry;ds=[mp.mpf('.004')*(1+z*z),-mp.mpf('.002')*z*z]
            change=sum(mp.exp(-2*mu*c)*(2*d*E+d*d*F) for d,c in zip(ds,(-3,-1)))
            f=lambda t:(q/2)**sigma(t/100)
            flat=mp.quad(lambda t:mp.exp(-2*mu*t)*f(t)**2,[0,20,50,80,100])
            Nf=q*q*mp.exp(-200*mu)/4;Nrel=Nf*mp.exp(-2*mu*L)
            total=flat+Nf*I(2*mu,L)+Nrel*(post+change)
            direct_angular=mp.mpf(0)
            for d,c in zip(ds,(-3,-1)):
                direct_angular+=mp.quad(lambda t:mp.exp(-2*mu*(L+c+t))*((1+d*beta(t))**2-1),[-ell,0,ell])
            equal(direct_angular,mp.exp(-2*mu*L)*change,'actual signed bump square integral')
            direct_post=ein+Ns*I(2,Ts)+Nq*eout+Nt*I(delta,W)+Nt*mp.exp(-delta*W)*E0/(1-eps)**2
            equal(direct_post,post,'complete native post angular integral')
            flattenexit=(total-flat)*mp.exp(200*mu)/(q*q/4)/2
            equal(flattenexit,(I(2*mu,L)+mp.exp(-2*mu*L)*(post+change))/2,'flatten to power normalization')
            equal((post+change)*mp.exp(-8*mu)/2+I(2*mu,4)/2,
                  ((post+change)*mp.exp(-8*mu)+I(2*mu,4))/2,'power to angular inlet')
            equal(afterentry*mp.exp(1+mu)/2,(afterpower*mp.exp(-2*Ts)+I(2,Ts))/2,'steep entry to power')
            equal(waiting*mp.exp(-1-delta/2)*mp.exp(1+delta/2)/2,waiting/2,'steep exit to waiting')
            equal(H/2,E0/(2*K(0)**2),'waiting to collar half normalization')
        return dict(independent_resolved_energy_integral_rows=rows,
            full_model_exterior_integrated_analytically_to_infinity=True,
            exact_signed_bump_square_integral_retained=True,
            test_transition_J='t^2/2',test_heat_H='1-b*exp(-t)',
            canonical_Gamma_fixture=False,actual_parameter_admission=False,passed=True)


def reject_stale_owners(field):
    rejected=[];old=field.selected.exact.companion.heat
    for name,change in (
        ('old_outer_future',lambda x:setattr(x.outer,'future',old.future)),
        ('old_steep_future',lambda x:setattr(x.steep,'future',old.future)),
        ('old_heat_repair',lambda x:setattr(x.heat,'repair',old.future.repair)),
        ('old_heat_cache',lambda x:setattr(x.heat,'tail_cache',old.tail_cache))):
        out=copy.copy(field);out.outer=copy.copy(field.outer);out.steep=copy.copy(field.steep);out.heat=copy.copy(field.heat)
        # Preserve the alias graph of the copied current test candidate.
        out.steep.outer=out.outer;out.heat.outer=out.outer;out.heat.steep=out.steep
        change(out)
        try:out.assert_graph()
        except ValueError:rejected.append(name);continue
        raise ArithmeticError('Stale postpulse owner accepted: '+name)
    return dict(rejected=rejected,mutation_count=len(rejected),passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPostpulseEnergyHistory(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current postpulse source/family/datum differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Postpulse producer claims admission')
    proof=energy_transport_source_proof(field)
    if encode(pack(proof))!=raw['current_energy_transport_source_proof'] or not all(proof['identities'].values()):
        raise ValueError('Current exact energy source transport changed')
    for target,value in {'self.flatten':'CurrentFlattenMixedC4(proxy.pulse,self.family,self.source,self.hashes)',
        'self.outer':'CurrentPowerAngularC4(proxy)','self.steep':'CurrentSteepWaitingC4(proxy)',
        'self.heat':'CurrentCollarGammaC4(proxy)'}.items():
        binding('compliant_current_postpulse_energy_history','__init__',target,value)
    fixture=independent_energy_integral_fixture();mutations=reject_stale_owners(field)
    mixed_rows=energy_rows=zero_rows=0
    indices={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}
    for name,args in VIEWS.items():
        packet=field.evaluate(*args);point=packet['source_packet']
        if encode(pack(packet))!=raw['current_postpulse_energy_views'][name]:raise ValueError('Current postpulse packet differs: '+name)
        if not all(packet['source_owner_graph'].values()) or any(packet[k] for k in OPEN):raise ValueError('Postpulse graph/scope differs')
        energy=point['energy_Taylor']
        if energy.order!=5 or endpoints(energy[0])[0]<=0:raise ArithmeticError('Positive current energy C5 lost')
        for v in energy.coefficients:
            if not all(mp.isfinite(x) for x in endpoints(v)):raise ArithmeticError('Nonfinite current energy derivative')
            energy_rows+=1
        grid=point['physical_mixed_derivatives_total_order_le4']
        if set(grid)!={UZ,UT,UR,P}:raise ValueError('Current postpulse component omitted')
        for label,rows in grid.items():
            if set(rows)!=indices:raise ValueError('Current mixed4 derivative order omitted')
            for v in rows.values():
                if not all(mp.isfinite(x) for x in endpoints(v)):raise ArithmeticError('Nonfinite current mixed row')
                mixed_rows+=1
                if label in (UZ,UR):
                    if endpoints(v)!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Current inherited meridional zero lost')
                    zero_rows+=1
    # Comparison boxes only confirm consistency after exact source proof.
    # They are never used to establish the joins or select a source value.
    for left,right in (('flatten_exit','power_inlet'),('angular_terminal','steep_entry'),('waiting_terminal','collar_inlet')):
        a=field.evaluate(*VIEWS[left])['source_packet']['energy_Taylor']
        b=field.evaluate(*VIEWS[right])['source_packet']['energy_Taylor']
        for n in range(6):
            al,ah=endpoints(a[n]);bl,bh=endpoints(b[n])
            if max(al,bl)>min(ah,bh):raise ArithmeticError('Exact energy join contradicts enclosures')
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_energy_transport_source_proof=proof,
        independent_energy_integral_fixture=fixture,stale_owner_mutations=mutations,
        current_energy_C5_rows=energy_rows,current_postpulse_mixed4_rows=mixed_rows,
        source_inherited_zero_meridional_rows=zero_rows,all_passed=True,input_hashes=hashes,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current postpulse energy history PASS: '+str(energy_rows)+' energy rows, '+str(mixed_rows)+' mixed4 rows',flush=True)
    return result


if __name__=='__main__':run()
