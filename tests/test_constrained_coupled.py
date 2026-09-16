from pathlib import Path
from dataclasses import replace
import numpy as np
from openai_ns_reconstruction.constrained_poloidal import PoloidalCandidate
from openai_ns_reconstruction.constrained_coupled import CoupledCandidate


def test_coupled_embedding_and_roundtrip(tmp_path):
    root=Path(__file__).resolve().parents[1]
    old=PoloidalCandidate.load(root/'artifacts/constrained/poloidal_anchor/candidate.json')
    c=CoupledCandidate(replace(old.base,coefficients=(0.,)*90),old.pressure_coefficients,old.base.coefficients+old.coefficients,old.anchor_core)
    rng=np.random.default_rng(18);x=rng.uniform(-1.7,1.7,(32,3));t=rng.uniform(.25,.75,32)
    np.testing.assert_allclose(c.velocity(x,t),old.velocity(x,t),atol=1e-14)
    np.testing.assert_allclose(c.pressure(x,t),old.pressure(x,t),atol=1e-14)
    path=tmp_path/'candidate.json';c.save(path)
    np.testing.assert_allclose(CoupledCandidate.load(path).velocity(x,t),c.velocity(x,t),atol=1e-14)
