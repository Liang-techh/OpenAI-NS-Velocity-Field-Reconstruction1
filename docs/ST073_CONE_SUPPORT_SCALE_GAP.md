# The current cone band is too narrow for the selected pulse duration

`moment_bridge_cone_support.py` maps the physical full-residual stress-cone
analogue of the three-knot moment-closed bridge at `eta=-0.2`. It samples
13 radial similarity positions from `y=0.035` to `0.065` at both `k=11`
and `k=19`, with the same Gauss-12 stress primitive used to select the
late-bridge wave source. This is a numerical cone diagnostic, not the
OpenAI paper's normalized leading-profile cone theorem.

| Scale | Passing sampled `y` band | Current radial half-width | Best-centered sampled band half-width | Half-width needed for `width²/nu = 0.04 tau` |
| ---: | --- | ---: | ---: | ---: |
| `k=11` | `0.0375–0.0575` | `1.99e-5` | `4.23e-5` | `3.13e-4` |
| `k=19` | `0.0375–0.0550` | `1.24e-6` | `2.31e-6` | `1.95e-5` |

At the existing center `y=0.05`, the contiguous passing samples give
symmetric half-widths `3.17e-5` and `1.32e-6`. Their pulse half-widths
are still about **97** and **218** radial diffusion times respectively.
The first failing sample on either side bounds a connected, centered
all-pass interval; even the more generous failure-bracket widths for a
*relocated* center are only `5.29e-5` and `2.97e-6`. The diffusion-balanced
widths exceed these bounds by factors **5.91** and **6.57**.

This is a bracket from finite samples and finite-difference cone data.
It rules out enlarging a connected wave support around this sampled
passing band to the diffusion-balanced width without crossing a sampled
cone failure. It does not rule out a disjoint admissible region, a new
mean background, a different time pulse, or the paper's Appendix-C
relaxed-cone/shear-loop construction. The axial and time extent of the
cone and the complete wave residual remain unchecked.

Together with `ST073_LOCAL_PULSE_INVERSE.md`, this explains why a local
amplitude inverse on the current patch asks for a velocity correction
larger than the background. The next constructive search should target
an **alternative radial–axial region or a redesigned moment-matched mean
profile** with a quantitatively wider cone margin. It must jointly choose
the pulse duration and carrier scale, then run the paper's supported
wave–mean–moment residual cycle on more than one dyadic band. A fixed
shape carried toward `k=19` does not cure this gap.
