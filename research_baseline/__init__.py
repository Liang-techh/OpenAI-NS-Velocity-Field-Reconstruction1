"""Repository-local API for the retained, NOT PDE-validated ST006 research field."""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
from ._hybrid import HybridFamily
from ._spacetime import force

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / 'artifacts/research/ST006'
CANDIDATE_SHA256 = '6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3'


def verify_integrity(root: Path = ROOT) -> dict:
    """Check the immutable field and recorded evidence, not scientific acceptance."""
    root = Path(root)
    manifest = json.loads((root/'artifacts/research/ST006/manifest.json').read_text())
    if manifest.get('candidate_sha256') != CANDIDATE_SHA256:
        raise ValueError('The retained candidate identity was changed')
    for flag in ('pde_validated', 'paper_exact', 'blowup_proved', 'openai_field_identified'):
        if manifest.get(flag) is not False:
            raise ValueError(f'Unsupported scientific claim: {flag}')
    checks = {manifest['candidate_file']: CANDIDATE_SHA256, **manifest['evidence_sha256'], **manifest['runtime_sha256']}
    for name, expected in checks.items():
        path = (root/name).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError('Manifest path exits the repository')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Checksum mismatch: {name}')
    return {'candidate': 'ST006', 'files_checked': len(checks), 'artifact_integrity': True, 'pde_validated': False}


@dataclass(frozen=True)
class PublishedField:
    """Immutable parameters. Coordinates and time are dimensionless."""
    _family: HybridFamily
    _raw: np.ndarray
    sha256: str = CANDIDATE_SHA256
    pde_validated: bool = False

    def fields(self, points, time):
        """Return (velocity[...,3], pressure[...]) with broadcast time and batching."""
        points = np.asarray(points, dtype=float)
        if points.ndim < 1 or points.shape[-1] != 3 or not np.isfinite(points).all():
            raise ValueError('points must be finite with shape (...,3)')
        shape = np.broadcast_shapes(points.shape[:-1], np.shape(time))
        x = np.broadcast_to(points, shape+(3,)).reshape(-1,3)
        t = np.broadcast_to(np.asarray(time, dtype=float), shape).reshape(-1)
        if not np.isfinite(t).all() or np.any((t < .25) | (t > .75)):
            raise ValueError('time must lie in the registered [.25,.75] interval')
        u = np.empty_like(x); p = np.empty(len(x))
        for start in range(0,len(x),512):
            sl = slice(start,start+512)
            u[sl],p[sl] = self._family.fields(self._raw,x[sl],t[sl])
        return u.reshape(shape+(3,)),p.reshape(shape)

    def at_points(self, points, time):
        return self.fields(points,time)[0]

    def pressure(self, points, time):
        return self.fields(points,time)[1]

    def forcing(self, points, time):
        # Validate the same registered coordinates/time before force evaluation.
        points = np.asarray(points,dtype=float)
        if points.ndim < 1 or points.shape[-1] != 3 or not np.isfinite(points).all():
            raise ValueError('points must be finite with shape (...,3)')
        time = np.asarray(time,dtype=float)
        if not np.isfinite(time).all() or np.any((time < .25) | (time > .75)):
            raise ValueError('time must lie in [.25,.75]')
        return force(points,time,*self._raw[-2:])

    def velocity(self,x,y,z,t):
        x,y,z,t = np.broadcast_arrays(x,y,z,t)
        return self.at_points(np.stack((x,y,z),axis=-1),t)


@lru_cache(maxsize=1)
def load_best() -> PublishedField:
    """Load ST006; 'best' means retained baseline in the documented root lane."""
    verify_integrity()
    data = json.loads((ARTIFACT/'candidate.json').read_bytes())
    if data.get('schema') != 'root_st001_compact_spacetime_v1' or data.get('basis_kind') != 'hybrid_legendre7_inverse_even_v1':
        raise ValueError('Unexpected candidate schema or basis')
    if data.get('pde_validated') is not False:
        raise ValueError('Candidate claims unsupported PDE acceptance')
    raw = np.array(data['coefficients'],dtype=float)
    family = HybridFamily(data['nr'],data['nz'],data['nt'])
    family.coefficients(raw)
    if np.max(np.abs(raw[:2*family.n])) > 4 or np.max(np.abs(raw[2*family.n:-2])) > 100 or np.any((raw[-2:] < 0) | (raw[-2:] > 10)):
        raise ValueError('Candidate parameter bounds violated')
    raw.setflags(write=False)
    return PublishedField(family,raw)
