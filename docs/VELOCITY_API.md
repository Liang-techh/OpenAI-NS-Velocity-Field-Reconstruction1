# Three-dimensional velocity API

This exports the saved, nonzero coupled candidate as Cartesian components. It is an executable velocity field, not a verified reproduction of the OpenAI visualization and not an accepted Navier–Stokes solution. Target provenance/mapping is tracked in VIS001.

## Install and evaluate

From the repository root:

```sh
python -m pip install -e .
ns-velocity --point .1 0 .1 .5
ns-velocity --export artifacts/visual/velocity_api --grid-size 17
```

Without installing, set PYTHONPATH=src and run `python -m openai_ns_reconstruction.velocity_components` with the same arguments. PowerShell: `$env:PYTHONPATH="src"`.

```python
from openai_ns_reconstruction.velocity_components import velocity, u, v, w, VelocityField

ux, uy, uz = velocity(0.1, 0.0, 0.1, 0.5)
# (-0.017550887073245407, 0.07337970151893687, 0.034940231516750846)

field = VelocityField()
vectors = field.at_points([[0.1, 0, 0.1], [0, 0, 0.2]], 0.5)  # (2, 3)
values = field.grid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], [.25, .5, .75])
# values.shape == (3, 3, 3, 3, 3): time, x, y, z, component
```

The separate u/v/w functions return the x/y/z velocity components. NumPy broadcasting is supported; scalar input returns Python floats. For many points, use one `velocity` or `at_points` call instead of evaluating each component separately. Evaluation is chunked to limit temporary memory.

## Coordinate and model contract

- Right-handed Cartesian x,y,z; z is the axial coordinate. u,v,w are velocity components along those axes.
- Dimensionless model variables; there is no established SI calibration or OpenAI camera/time mapping yet.
- Valid time interval [.25,.75]; outside this interval raises ValueError rather than extrapolating.
- Exact compact zero exterior when sqrt(x²+y²)>=2 or |z|>=2. The interior is evaluated from the actual saved candidate, including the axis.
- Default coefficients are packaged at `src/openai_ns_reconstruction/data/velocity_candidate.json`, copied verbatim from `artifacts/constrained/coupled_joint/candidate.json`. `VelocityField(path)` allows an explicit compatible coupled candidate.
- The metadata includes the SHA256 of the selected candidate. Bundled data makes the installed package independent of the repository's artifacts directory.

## Mathematical structure

The implemented field has the form U=U_parent+sum(a_i S_i)+sum(b_j P_j), with90 swirl modes and27 poloidal modes. The full117 velocity coefficients, pressure data and parent field parameters are serialized in the bundled JSON. Swirl modes have Cartesian form(-y F,x F,0) times their temporal factors. Poloidal modes derive from psi=r² z F(r²,z²,t), giving(-x(F+2z² F_v),-y(F+2z² F_v),2z(F+r² F_s)), where s=r² and v=z². The code supplies compact envelopes, core constraints, coordinate scaling and temporal factors; no velocity values are invented by the API.

The formulas and coefficients are implemented by constrained_coupled.py, constrained_poloidal.py, constrained_inner_swirl.py and their parent modules. The public API is a stable entry point into these implementations.

## Saved sample

`artifacts/visual/velocity_api/grid.npz` contains x,y,z coordinate axes, times=[.25,.5,.75], and u,v,w arrays, each shape(3,17,17,17). `metadata.json` records coordinate conventions and hashes. Load with:

```python
import numpy as np
sample = np.load("artifacts/visual/velocity_api/grid.npz")
U = np.stack([sample["u"], sample["v"], sample["w"]], axis=-1)
```

Two focused tests passed for direct-candidate agreement, scalar/array behavior, axis/exterior cases, time-domain errors and chunk ordering. This does not establish visual correspondence or PDE acceptance. The last standard-step PDE validation remains failed, with sampled momentum maximum1.00221830.
