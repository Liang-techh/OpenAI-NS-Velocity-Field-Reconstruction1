# Complete continuous pulse enclosure

The complete weighted pulse rows now have outward interval bounds for the
declared decimal `mu`. This is an integral certificate, not a certificate of
the inherited input parameter or the global velocity field.

For `lambda = .5 - row*mu`, `k = lambda/mu`, let `u0=(2/k)^(1/3)`,
`L=1/u0^2`, `w=u0/sqrt(6)`, `v=1+w*x`, and `u=u0*v`. The centered integrand is

```
(10.99-u)*exp(1/(1-u)^2 - x^2*(2*v+1)/(6*v^2))
 / (1+exp(-1/u^2+1/(1-u)^2)).
```

Its prefactor is `exp(-2*k-3*L)*u0^2/(mu*sqrt(6))`. Interval rectangles
enclose the entire integrand range on every panel of `[-band,band]`.
There is no unbounded quadrature remainder. The code checks that both band
edges remain inside `0<u<.5` and `band*w<.25` using outward endpoints.

The positive omitted pieces are bounded separately: the small-u cutoff outside
the band by phase monotonicity, `.5<=u<=1` by `5.5/mu*exp(-2.5*k)`, and
`0<=xi<=10` by `11/lambda*exp(-3*k)`. The startup primitive is nonnegative and
at most `xi`, including below `.02`; it is not replaced by `xi-.01` there.
The lower full bound is the lower centered contribution; the upper includes
all three positive omitted bounds. Every operation uses a local interval
context. Exact dyadic endpoints are saved; decimal displays are approximate.

Both rows' nominal centers lie inside the independent bounds at 1024 and
4096 panels. Relative centered widths are approximately 7.77% and 1.89%.
These are conservative first bounds, not 80-digit integral accuracy.
The live runtime exposes cached `pulse_enclosure(row, panels=..., precision=...)`
and checks its owned pulse provider against these bounds.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_pulse_enclosure.py`
and `python experiments/root_st073/lei_ren_part1_paper_continuous_axial_runtime.py`.
Both checks pass; terminal mass and Z residuals remain nonzero.

Next: tighter correlated pulse/basis bounds and propagation through the
ill-conditioned coefficient solve. Incoming integrals, angular input, pulse
energy target, full mean closure, finite energy, stress, and scale recursion
remain unproved. No terminal mean has been overwritten by zero.
