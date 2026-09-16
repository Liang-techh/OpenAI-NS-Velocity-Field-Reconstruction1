# VIS001: OpenAI Navier–Stokes visualization target

Observation date: 2026-09-16.  The public source that can be identified from
the repository and the supplied shared goal is the figure in OpenAI's
announcement:

- [OpenAI announcement: On the Navier–Stokes Millennium Prize Problem](https://openai.com/index/navier-stokes-solution/)
- [OpenAI paper: Finite Time Blowup for Navier–Stokes](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
- [Shared goal discussion](https://chatgpt.com/share/6aaab4a8-759c-83e9-b3e8-e45eb954fd16)

## Identified public media

The announcement page places the target under “The result”.  Its DOM exposes
a vortex illustration. Its caption describes local incompressible motion,
with warm and cool colors distinguishing faster and slower angular rotation,
and trajectories depicting inward motion and axial stretching.

The responsive media URLs observed on that page are:

- [light asset](https://images.ctfassets.net/kftzwdyauwt9/75EbpsuBOy5LbUgCppWXD1/88da8c19dcf76d6f4f8fd7485dcb6346/navier-stokes-light-master.png)
- [dark asset](https://images.ctfassets.net/kftzwdyauwt9/38fsDDtQmNLBW3F74od13W/7f0213895c06edc97709dc82d7905bc3/navier-stokes-dark-master.png)

The page selected the `.png` asset through a Contentful image URL with
`w=3840&q=90&fm=webp` parameters during observation.  The target element was a
static image; no frame counter, time slider, playback control, or target video
element was exposed.  The browser reported the loaded image as 557 by 557
intrinsic pixels at the observed responsive size.  The figure has a dark
rectangular canvas, a vertically oriented slender central vortex, many curved
and helical trajectories, and the visible labels “inward spiral” and “axial
stretching”.  Teal/blue trajectories occupy much of the outer and middle
region, while orange trajectories are concentrated around the faster rotating
central portions.  These are visual observations only; the figure has no
numeric axes, coordinate ticks, units, legend scale, or time annotation.

The announcement calls this a “snapshot”, so it does not establish a frame
range or a numerical value of `t`.  It also does not publish the trajectory
seeds, camera projection, color normalization, or the sampled vector data.

## Coordinate and time information from the paper

The paper defines cylindrical coordinates about the vertical axis by

```text
x = (r cos(theta), r sin(theta), z)
e_r = (cos(theta), sin(theta), 0)
e_theta = (-sin(theta), cos(theta), 0)
e_z = (0, 0, 1).
```

Here `r` is distance from the axis, `theta` is the angle around it, and `z` is
height.  The leading field is axisymmetric, so its cylindrical components do
not depend on `theta`; `u_r`, `u_theta`, and `u_z` mean radial, azimuthal, and
axial velocity.  The paper's physical description says that the central flow
spirals inward toward the axis, flows upward and downward on opposite sides of
a dividing layer near `z = 0`, and has radial outflow farther above and below.
Outside the annulus, the heat exterior is purely azimuthal and independent of
height, with speed decreasing as `r` increases.

For the leading self-similar description (paper Eq. (3.2), restated in the
leading-flow discussion), `tau = 1 - t`,

```text
tau = q (1 - eta^2),   z = q^D eta,   X = r^2 / (2 q),
A = 1/2 + h,            D = 1/2 - h,
q > 0,                   -1 < eta < 1,   0 < h < 1/100,
```

and `q` is the unique positive solution of
`q - z^2 q^(-2h) = tau`.  The leading profiles are described by

```text
u_theta^(0) = q^(-A) E(X, eta),
u_z^(0)     = q^(-A) U(X, eta),
r u_r^(0)  = V_0(X, eta),
p^(0)      = q^(-2 A) Pi(X, eta).
```

The stated core scales are `ell_r ~ tau^(1/2)`,
`ell_z ~ tau^(1/2 - h)`, `|u_theta|, |u_z| ~ tau^(-1/2 - h)`, and
`|u_r| = O(tau^(-1/2))` as `t` approaches the singular time `1`.  These are
paper-level qualitative/time-scaling constraints, not a numerical extraction
of the article image.

The paper's Figure 1 (printed page 4, PDF page 3) is a separate schematic with
three labels `t1 < t2 < t3 < 1`; it is not the color figure above and does not
provide samples of `u`, `v`, or `w`.  The announcement figure likewise gives no
numerical frame time.

## Resolution and unresolved ambiguity

For VIS004, the identified target is therefore the static announcement figure
and its qualitative properties: a slender, axis-centered, inward-spiraling,
axially stretched vortex with radius-dependent rotational coloring.  A field
constructed from the paper coordinates and scales can be compared against
those properties, but it must not be labeled an exact frame match.

If “the OpenAI visualization” means a time-resolved animation or a different
image that was shown in another conversation, that asset is not identifiable
from the cited shared goal or the public announcement.  The missing reference
needed for an exact animation comparison is the original media URL or an
attached frame/screenshot together with its frame/time index.  Until that is
provided, the target's exact `t`, domain/units, camera, trajectory seeds, color
map, and numerical `u(x,y,z,t), v(x,y,z,t), w(x,y,z,t)` remain unresolved.
