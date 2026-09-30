"""Necessary source core/reference scale relations; no theorem certification.

Sections 8--10 require Rref=110(Cstar Pstar)^10. Independent demo
outer parameters cannot be silently reused as a matched regular-core field.
"""
import json
from pathlib import Path
import mpmath as mp


def core_scale_requirements(*,j='.02',logPstar='14',logRref='10',
                            Lambda=None,A_Omega_lower_bound='0',precision=100):
    with mp.workdps(precision):
        j=mp.mpf(str(j)); logP=mp.mpf(str(logPstar)); logR=mp.mpf(str(logRref))
        A=mp.mpf(str(A_Omega_lower_bound))
        if not 0<j<=mp.mpf('.05') or A<0:
            raise ValueError('Require 0<j<=1/20 and a nonnegative domain bound')
        minimum_Lambda=max(mp.mpf(500),j**-2)
        lam=minimum_Lambda if Lambda is None else mp.mpf(str(Lambda))
        if lam<minimum_Lambda:
            raise ValueError('Lambda must be at least max(500,j^-2)')
        minimum_logC=2*mp.log(lam)+lam*A
        minimum_logR=mp.log(110)+10*(minimum_logC+logP)
        implied_logC=(logR-mp.log(110))/10-logP
        return {'j':mp.nstr(j,precision),'Lambda':mp.nstr(lam,precision),
            'logPstar':mp.nstr(logP,precision),'candidate_logRref':mp.nstr(logR,precision),
            'A_Omega_lower_bound':mp.nstr(A,precision),
            'minimum_logCstar':mp.nstr(minimum_logC,precision),
            'necessary_minimum_logRref':mp.nstr(minimum_logR,precision),
            'implied_logCstar':mp.nstr(implied_logC,precision),
            'necessary_scale_test_passed':bool(logR>=minimum_logR),
            'reference_relation':'logRref=log(110)+10*(logCstar+logPstar)',
            'source_core_certified':False,
            'scope':'Necessary lower bound only; actual complex-domain A_Omega, Section 9/10 bounds and PDE core remain required.'}


if __name__=='__main__':
    receipt=core_scale_requirements()
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))
