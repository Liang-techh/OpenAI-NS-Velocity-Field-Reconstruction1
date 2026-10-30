"""Log-safe factorization of the actual first-collar flat cutoff.

For the original sigma(s) with odds exp((1-s)^-2-s^-2), this leaf proves
sigma and sigma' are bounded by the same source hb**2 before emitting
interval sector factors. It does not evaluate or select sigma, sigma', or hb.
"""
from __future__ import annotations


def build(c, *, sc, log_hb, source_family, source_atoms, source_contract):
    sc,log_hb=c.mpf(sc),c.mpf(log_hb)
    if sc.a<=0 or sc.b>=c.mpf(1)/1024:raise ValueError('original microscopic sc required')
    if log_hb.b>=c.ln(c.mpf('0.5')).a:raise ValueError('0<hb<1/2 required')
    if source_family is None or set(source_atoms)!=set(('hb','F0_anchor','Pstar_squared','Lambda')):
        raise ValueError('four same-source atoms required')
    keys=('exact_original_flat_sigma','same_hb_log_and_sc','hb_and_phase_Z_independent')
    if any(source_contract.get(k) is not True for k in keys):raise ValueError('source contracts required')
    W=sc/4
    # On 0<s<=W<=1/2: sigma<=exp(4-1/W^2), and
    # sigma'<=18 W^-3 exp(4-1/W^2). The latter uses monotonicity of
    # exp(-1/s^2)*s^-3 on (0,1/2].  Zero is the flat extension.
    value_log=4-1/(W**2)
    derivative_log=c.ln(18)+4-1/(W**2)-3*c.ln(W)
    target=2*log_hb
    if not (value_log.b<=target.a and derivative_log.b<=target.a):
        raise ValueError('source-specific sigma/hb^2 flatness not proved')
    rho_upper=(4/(1-W/2)).b
    rho=c.make_mpf((c.mpf(4)._mpi_[0],rho_upper._mpi_[1]))
    invrho=c.mpf(1)/rho
    zero=c.mpf(0);one=c.mpf(1);positive=c.make_mpf((zero._mpi_[0],one._mpi_[1]))
    def term(power,coefficient,origin):return {'powers':(power,0,0,0),'coefficient':coefficient,'origin':origin}
    # chi=1-hb^2*(sigma/hb^2)+hb^3*(sigma/hb^2).
    # chi_r=-hb*(sigma'/hb^2)/rho+hb^2*(sigma'/hb^2)/rho.
    chi0=[term(0,one,'unit'),term(2,-positive,'minus_sigma_factor'),term(3,positive,'hb_sigma_factor')]
    chir0=[term(1,-positive*invrho,'minus_sigma_prime_factor'),term(2,positive*invrho,'hb_sigma_prime_factor')]
    chi={k:([*chi0] if k==0 else [term(0,zero,('exact_Z_independent',k))]) for k in range(6)}
    chir={k:([*chir0] if k==0 else [term(0,zero,('exact_Z_independent',k))]) for k in range(6)}
    return {'source_family':source_family,'source_atoms':source_atoms,'rho':rho,'x':c.sqrt(rho),
      'chi':chi,'chi_r':chir,'source_proof':{'sigma_log_upper':value_log,'sigma_prime_log_upper':derivative_log,
        'twice_same_hb_log':target,'sigma_and_derivative_le_hb_squared':True,
        'proof_uses_original_flat_extension':True,'caller_source_contract_independently_proved':False},
      'formal_factor_semantics':'sigma/hb^2 and sigma_prime/hb^2 are unknown actual functions in [0,1]',
      'midpoint_or_endpoint_used_as_sigma':False,'shared_cutoff_factor_names_preserved':True,'interval_coefficients_not_jointly_correlated':True}
