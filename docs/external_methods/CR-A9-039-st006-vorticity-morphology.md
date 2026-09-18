# CR-A9-039 — retained ST006 3-D vorticity morphology time series

## Purpose

This increment is candidate-side visualization evidence for the already published retained ST006 field. It does not add or fit a velocity family. It samples the public repository API, reconstructs Cartesian vorticity on fixed grids, and summarizes the three-dimensional enstrophy distribution by a weighted covariance eigensystem.

The resulting principal/transverse RMS extents, aspect ratio, covariance-volume proxy, vorticity strength and registered-coordinate principal axis can be used to review whether a candidate becomes more elongated or spatially concentrated across the fixed project window. They are descriptive candidate measurements only. No OpenAI numerical target, camera registration, hidden frame time or visual pass threshold is introduced.

## External result screened

- source repository: `numpy/numpy`
- screened commit: `38a7122a3ab009df1062c6a7d1139f9449865173`
- public API used: `numpy.linalg.eigh`
- upstream implementation/documentation location: `numpy/linalg/_linalg.py`
- license: NumPy BSD 3-clause terms
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none; NumPy is already a repository dependency

NumPy documents `eigh` for real symmetric/Hermitian matrices and returns the associated eigenvalues/eigenvectors. This increment uses that public API only for the local symmetric 3x3 enstrophy covariance tensor.

## Local implementation and differences

The vorticity diagnostic, weighting, grid contract, candidate identity binding, covariance construction, principal-axis sign canonicalization, multi-resolution drift report, hashing and truth semantics are local repository code. They are not NumPy algorithms or claims.

The covariance is built from weights `|curl(u)|^2` on centered-difference interior nodes. Its eigenvectors are treated as geometric axes of a candidate-side weighted point cloud, not as vortex-core lines, stability modes or PDE eigenfunctions. The principal axis is unsigned; its sign is canonicalized only to make JSON receipts deterministic.

The fixed contract is autonomous engineering choice:

- physical Cartesian box `[-2,2]^3`, already the registered evaluation box;
- times `0.25, 0.50, 0.75` inside the registered window;
- odd resolution ladder `17, 25, 33`;
- centered second-order Cartesian differences;
- no rotation, camera registration or hidden-time alignment.

The resolution ladder reports sensitivity only. It defines no visual acceptance threshold and is not the registered independent PDE derivative ladder. Derivatives computed here are explicitly forbidden as PDE-acceptance evidence.

## Truth boundary

This method does not change `[u,v,w]`, pressure, forcing, support, energy normalization, optimization/validation samples, seeds, residual definitions or scientific thresholds. It cannot establish `visualization_ready`, OpenAI visual correspondence, PDE validity, paper exactness, OpenAI-field identity or a blow-up result. All of those states remain false in the generated receipt.
