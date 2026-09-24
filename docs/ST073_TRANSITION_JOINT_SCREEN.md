# Joint solenoidal transition and pressure screen

The pressure-only correction widened the local wave cone but did not reduce the momentum residual. This experiment adds a compact solenoidal poloidal mode in the axial collar:

`psi = 15 R^2 y^2 (1-y^2)^5 B(s)(1+4s)`, with `y=r/R`, `s=(|z|-zflat)/(zsupp-zflat)`, and `B(s)=1024 s^5(1-s)^5`. The code evaluates `u_r=-psi_z/r`, `u_z=psi_r/r`; its support is the existing axial collar and its velocity is divergence-free by construction. The mode is in `joint_collar_fit.py`. Four compact pressure shapes from the preceding screen are fitted independently.

The local objective in `transition_poloidal_pressure_screen.py` selects poloidal amplitude `0.008`. At `tau=0.5/64`, fixed `r=0.0056890761915166545`, all 13 constrained axial nodes pass with cone ratio below `0.8`; all 12 unfitted midpoints also pass, with largest ratio `0.7994`. The sampled radial-line residual RMS falls to `79,008`, compared with `183,876` for the original field. The sampled center maximum falls from `737,608` to `439,964`.

The full-support check **rejects promotion of this local fit**. On a Gauss6 physical-volume grid at the same time, its momentum L2 is `370.7` versus `190.1` for the original field, and its sampled maximum is `434,828` versus `91,544`. At `tau=0.5*2^-5.5`, L2 is `247.0` versus `158.0`. See `transition_poloidal_pressure_audit.py` and its JSON report.

The follow-up `transition_global_constrained_screen.py` minimizes the Gauss6 physical-volume residual while preserving the 13 local cone constraints. Across poloidal amplitudes from `-0.010` to `0.010`, only zero and positive amplitudes are feasible in this screen; the best is zero poloidal amplitude with L2 `210.6`, still above the original `190.1`. This is evidence of a conflict **within these five modes and this finite sample**, not an impossibility result for a different transition construction.

Do not use either candidate as an accepted Navier–Stokes field. The next construction should introduce a solenoidal radial/axial transition basis that targets the full-support residual while keeping the local cone open, and should fit on a volume objective from the start. After that, test the resulting exact-curl nonaxisymmetric wave over a radial-axial-time patch. Full momentum max and physical-volume L2 below `1e-3` remain unmet.

Reproduce with:

```text
python experiments/root_st073/transition_poloidal_pressure_screen.py
python experiments/root_st073/transition_poloidal_pressure_audit.py
python experiments/root_st073/transition_global_constrained_screen.py
```
