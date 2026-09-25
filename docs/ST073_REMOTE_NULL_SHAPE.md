# Moment-null remote swirl shape screen

The remote `[4,16]` patch now admits an extra normalized swirl mode
`N(X)=B(y)(2y-1-alpha)`, where `alpha` is its `sqrt(2X)`-weighted projection
onto `B`. Thus `integral sqrt(2X) N(X) dX = 0`. Adding `a_N N` leaves the
angular moment equation unchanged. For each `a_N`, the meridional coefficient
is re-solved from the kinetic moment and fit smoothly across the registered
eta/time slices. The physical correction is still built using the same
axisymmetric streamfunction, so the new swirl mode adds no divergence.

At `tau=.0084`, amplitudes `a_N=-1,-.5,0,.5,1` were screened at
`eta=0,.2,.35` and `X=6,10,14`. Direct slice moment defects at these nine
points remained below `9.90e-4`; the fitted physical field preserves the two
radial identities to this sampled accuracy. **None of the 45 physical stress
analog samples passed the strict local cone.** The decisive `X=10, eta=0`
stress projection onto `N` moved only from `3.67` at `a_N=-1` to `3.48` at
`a_N=.5`, remaining on the wrong side of zero. Strongly varying this single
null swirl shape is therefore not sufficient to repair the stress direction.

The tested cone is a physical-field analog, not the paper's normalized
leading-order cone. This is an exclusion of the sampled one-mode family,
not a proof that no moment-matched mean field admits a suitable wave. The
next constructive search should vary a **meridional streamfunction shape**
or several independent moment-null radial modes, including eta dependence,
and evaluate the normalized cone before constructing the nonaxisymmetric
phase/amplitude waves. Full momentum, finite energy and critical-time forcing
remain open.

Run `python experiments/root_st073/remote_null_shape_screen.py`; the report
is `experiments/root_st073/remote_null_shape_screen.json` with
`accepted: false`.

The follow-up meridional-shape screen in `ST073_REMOTE_MERIDIONAL_SHAPE.md`
can reverse the stress projection but still misses the full cone.
