# Fixed-order core horizon under dyadic contraction

`experiments/root_st073/core_critical_horizon.py` extends the ST073-V
local radial field beyond its registered `k<=6` time window as a
diagnostic, with `tau=.5*2**(-k)`. It samples 25 fixed similarity
locations in `0.001<=X<=1/64`, `|eta|<=.3`, keeping all parameters
except radial truncation order fixed.

| Radial order | Momentum sample max, k=14 | k=16 | k=18 | First sampled `1e-3` failure in tested k |
| ---: | ---: | ---: | ---: | ---: |
| 8 | `5.60e-4` | `4.06e-3` | `2.95e-2` | 15 |
| 10 | `2.01e-6` | `1.45e-5` | `1.04e-4` | none through 18 |
| 12 | `3.80e-9` | `3.60e-8` | `2.15e-7` | unresolved by floating-point cancellation |
| 14 | `2.24e-9` | `1.80e-8` | `8.09e-8` | unresolved by floating-point cancellation |

The eighth-order failure at `k=15` is a **fixed truncation failure**,
not evidence that the underlying recurrence fails there: the tenth-order
field has the same axis data and remains below the gate at the sampled
points through `k=18`. Both resolved fixed-order residuals increase as
physical scales contract. Near the higher orders, the four PDE terms
almost cancel. On this Windows host `np.longdouble` is only float64;
at `k=18`, the twelfth- and fourteenth-order residuals are only about
`2.2` and `0.9` times `eps * sum(term norms)`. Their reported small
values cannot establish further convergence or a critical-time bound.

The result separates two requirements: normalized profile transfer
can look stable while an absolute physical momentum error grows; and
a scale-recursive construction needs an order-uniform or otherwise
controlled remainder as `tau` approaches zero. The original eighth-
order local field is sampled below `1e-3` on its registered window,
but no fixed finite-order extrapolation here proves the requested
critical-time field. A sampled failure rules out that particular
unchanged field at that point; sampled success does not prove a domain
maximum, spatial-volume L2, finite energy, or inner/outer matching.
