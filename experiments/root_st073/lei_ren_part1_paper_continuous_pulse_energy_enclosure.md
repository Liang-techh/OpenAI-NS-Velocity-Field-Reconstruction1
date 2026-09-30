# Pulse energy with a bounded startup primitive

`ContinuousPulseEnergyEnclosure` encloses the complete energy integral
`Kp = integral_0^11 exp(-2*xi)*gp(xi)^2 dxi` of the installed pulse definition.
It uses a local outward interval context and exports exact dyadic endpoints.

For `0<=xi<=.02`, write `q=50*xi` and
`gp(xi)=integral_0^q sigma(s)ds/50`. The switch sigma is increasing.
On every rational grid panel its left and right node values therefore give
lower and upper rectangle bounds. Accumulating these bounds encloses the
primitive at every node; its entire range on a panel is enclosed between the
left lower primitive and right upper primitive. Squaring this interval and
multiplying by the outward exponential range bounds the startup energy.
The endpoints sigma(0)=0 and sigma(1)=1 are explicit. No singular endpoint
division or unbounded nested quadrature is used.

On `[.02,10]`, `gp=xi-.01` and the energy antiderivative is analytic.
On `[10,11]`, `gp=(xi-.01)*sigma(11-xi)`; monotonic switch node bounds and
interval rectangles enclose the cutoff energy. All support is included.

The saved nominal Kp lies inside both the 1024- and 4096-panel bounds.
The total relative widths are approximately `1.64e-8` and `4.10e-9`.
These are proven conservative widths from range integration, not extrapolated
quadrature refinement errors. The coefficient encloser now consumes the whole
Kp interval, producing a 4096-panel relative amplitude width around `2.05e-9`.
Incoming row and energy-target uncertainty are still unbounded; no full-field
finite-energy or mean-closure certificate follows.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_pulse_energy_enclosure.py`.
The diagnostic compares the installed receipt value, rather than rerunning
the slow nested quadrature as a proof substitute. The installed materialized
pulse amplitude and coefficients remain unchanged.
