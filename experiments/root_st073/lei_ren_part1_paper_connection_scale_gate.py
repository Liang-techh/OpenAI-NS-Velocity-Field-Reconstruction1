"""Necessary Section 9 input tests from axis values, not a certificate.

A_Omega=max(-Re G) on the complex neighborhood is not the real maximum
of G. This module does not substitute one for the other.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets


def run(precision=160,*,Lambda='1e36',j='.02',logC=None,
        logPstar='14',delta='1e-32',h_b='.005'):
    with mp.workdps(precision):
        lam=mp.mpf(str(Lambda))
        logC=2*mp.log(lam) if logC is None else mp.mpf(str(logC))
        logP=mp.mpf(str(logPstar)); hb=mp.mpf(str(h_b))
        axis=RegularCoreAxisJets(j=j,Lambda=lam,logC=logC,delta=delta,precision=precision)
        endpointG=max(abs(axis.G(-1)),abs(axis.G(1)))
        # Axis values belong to the C3 norm defining A in Eq9.16.
        A_lower=10+lam*endpointG
        logK_lower=max(mp.log(10**6),mp.log(lam/4),logC,logP,
                       mp.log(A_lower),logC+lam*endpointG)
        logR=mp.log(110)+10*(logC+logP)
        logc_upper=mp.log(1/(mp.mpf(10**6)*2001))-mp.log(10**10)
        j_upper=1/(mp.mpf(32)*10**6*2001)
        required_loghb_upper=logc_upper-100*logK_lower
        required_logC_shape=40*A_lower+1+mp.log(1+A_lower)-logP
        margins={
            'Cstar_ge_exp4A':logC-4*A_lower,
            'Rref_ge_long_shape_bound':logR-mp.log(110)-400*A_lower-10-10*mp.log(1+A_lower),
            'Rz_ge_axial_bound':logR-8-mp.log(110)-2*logK_lower+2*logP,
            'hb_small_enough_for_K_lower':required_loghb_upper-mp.log(hb),
            'j_below_universal_necessary_upper':j_upper-axis.j}
        n=lambda x:mp.nstr(x,45)
        receipt={'Lambda':n(lam),'j':n(axis.j),'logCstar':n(logC),'logPstar':n(logP),
            'logRref':n(logR),'axis_endpoint_G_max':n(endpointG),
            'A_necessary_lower_bound':n(A_lower),'logK_necessary_lower_bound':n(logK_lower),
            'j_necessary_upper_bound_Kp_ge_1':n(j_upper),
            'loghb_necessary_upper_bound':n(required_loghb_upper),
            'shape_only_necessary_minimum_logCstar':n(required_logC_shape),
            'gate_margins':{k:n(v) for k,v in margins.items()},
            'necessary_gates_passed':{k:bool(v>=0) for k,v in margins.items()},
            'A_Omega_certified':False,'A_C3_norm_bounded':False,'K_C3_norm_bounded':False,
            'source_parameter_regime_certified':False,
            'scope':'Rejection from necessary lower bounds; larger Cstar alone does not certify a new shared candidate'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))
    return receipt


if __name__=='__main__':run()
