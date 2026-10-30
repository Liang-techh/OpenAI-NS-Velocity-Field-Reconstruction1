"""Independent checks for log-safe original flat-cutoff factorization."""
from __future__ import annotations


def run(leaf,c):
    import mpmath as mp
    sc=c.mpf(1)/2048;h=c.mpf(1)/3
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    contract={'exact_original_flat_sigma':True,'same_hb_log_and_sc':True,'hb_and_phase_Z_independent':True}
    out=leaf.build(c,sc=sc,log_hb=c.ln(h),source_family='synthetic',source_atoms=atoms,source_contract=contract)
    checks={'flat_log_proof':out['source_proof']['sigma_and_derivative_le_hb_squared'],
      'chi_hb_powers':[x['powers'][0] for x in out['chi'][0]]==[0,2,3],
      'chi_r_hb_powers':[x['powers'][0] for x in out['chi_r'][0]]==[1,2],
      'Z_independent':all(out['chi'][k][0]['coefficient'].a==0 and out['chi'][k][0]['coefficient'].b==0 and out['chi_r'][k][0]['coefficient'].a==0 and out['chi_r'][k][0]['coefficient'].b==0 for k in range(1,6)),
      'not_self_certified':out['source_proof']['caller_source_contract_independently_proved'] is False}
    def eval_terms(terms):
        x=c.mpf(0)
        for t in terms:x+=t['coefficient']*h**t['powers'][0]
        return x
    mp.mp.dps=80;mh=mp.mpf(1)/3;W=mp.mpf(1)/8192
    for q in (mp.mpf('0.5'),mp.mpf(1)):
        s=W*q;rho=4*mp.exp(mh*s);odds=1/(1-s)**2-1/s**2;e=mp.exp(odds);sig=e/(1+e)
        sigp=sig*(1-sig)*(2/(1-s)**3+2/s**3)
        actual=(1+(mh-1)*sig,(mh-1)*sigp/(mh*rho))
        for name,expected in zip(('chi','chi_r'),actual):
            v=eval_terms(out[name][0]);ref=c.mpf(mp.nstr(expected,65))
            checks[name+'_sample_'+str(q)]=bool(v.a<=ref.a and v.b>=ref.b)
    try:leaf.build(c,sc=sc,log_hb=c.mpf(-1000000000),source_family='synthetic',source_atoms=atoms,source_contract=contract)
    except ValueError:checks['reject_unproved_hb_square']=True
    else:checks['reject_unproved_hb_square']=False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,
       'synthetic_only':True,'reference':'independent high-precision original sigmoid and derivative'}
