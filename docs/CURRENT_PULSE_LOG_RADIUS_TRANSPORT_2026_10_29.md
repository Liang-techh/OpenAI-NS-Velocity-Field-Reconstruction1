# Pulse log-radius ordinary-Z transport — 2026-10-29

Goal ACTIVE / INCOMPLETE. A new source-free, caller-owned point interface carries the pulse radius in logarithmic form and includes moving-coordinate ordinary-Z derivatives without materializing the physical radius. It addresses a numerical representation needed by the F endpoint work; it does not establish an actual source-owner identity or evaluate a physical pulse. No old coordinate fixture, graph traversal, quadrature or ancestor computation was run.

For `ell=log(R)`, main/gap coordinates and the end coordinate are transported as

```
ell = logRp + xi/mu
d_Z ell = logRp_Z + xi_Z/mu - xi*mu_Z/mu^2
ell_end = logRp + 13/mu + s
d_Z ell_end = logRp_Z - 13*mu_Z/mu^2 + s_Z.
```

The returned `R_Z=exp(ell)*ell_Z` remains a formal factor. The callable requires the caller to provide the ordinary-Z input derivatives and accepts finite point numeric inputs only; interval data are rejected and no error enclosure is claimed. Rational inputs retain exact arithmetic. It carries each additive term separately: a rounded sum is not an identity or an accurate representation of a smaller term.

At the gap upper endpoint, `xi=13-4mu` and `xi_Z=-4mu_Z` give the direct structured form `ell=logRp+13/mu-4` and `ell_Z=logRp_Z-13*mu_Z/mu^2`. The end-side formula at `s=-4`, `s_Z=0` is the same **when both sides use the same input functions and coordinate interpretation**, even if mu_Z and logRp_Z are nonzero. The implementation keeps `13/mu` and `-4` separately, rather than relying on a finite-precision `13-4mu` before division. It also exposes generic main/gap and end inputs; it never sets missing `s_Z` to zero.

The single frozen, newly scoped fixture returned **23/23** controls. It checked synthetic exact Fraction inputs with nonzero mu_Z and logRp_Z, the moving upper endpoint, main/gap matching at xi=11 under shared inputs, formal radius scaling, rejection of zero/nonfinite mu and interval data, and a synthetic small-mu case where a rounded endpoint loses `4mu` but the separated `-4` survives. Those checks are algebraic/point-input tests, not a source-family identity, certified bound, quadrature or physical field value.

Static producer reading locates Gamma's `mu` and `logRp` through the historical native frame. The documented mu formula is constant-only; logRp is built from `phase_rc + log(r_minus) - 2 - 60 log(mu)`. An equivalent accepted Rc recipe is `log(110)+10 log(Cstar)+11 exp(40)+124`; its producer-level expression does not alone establish that Gamma's archived `phase_rc` handle and the current source are semantically identical or Z-independent. Those source-family/handle and derivative bindings, capture-free substitution, actual F finite endpoints, common measure, orientation and Pstar remain open. Do not promote matching templates or synthetic shared inputs to proof of F/P/G equality.

Retain the 212/213/183 inherited counts, six signed channels, one reserved contribution, four physical pulse corrections, independent reference offset and all prior errors. P1 and M remain distinct from m15. The previously documented overly broad JSON search incident remains preserved and must not be repeated. Next: establish same-source Gamma/frame parameter ownership and ordinary Z-independence, then bind historical F endpoint substitution and common physical scale before using this transport with actual source jets and certified error propagation. Global five-moment closure, heat/stress, true order-n recursion and Cartesian NS validation remain incomplete. Runtime model metadata was unavailable (`null`).
