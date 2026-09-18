"""Concrete autonomous realization of the Kokuno slow squared-partition contract.

The corrected 2026-09-09 reconstruction specifies the *structure* of the slow
partition used by the oscillatory field,

    sum_ell chi_ell(q)^2 = 1,
    eta_beta = chi_ell(q) chi_{ell,a},
    sum_beta eta_beta^2 = 1,

with beta=(ell,a), three slow chart coordinates, and per-band product-grid
mesh ``S_*^-3 = ell^-6``.  The rectangle sign sigma is attached after this
beta partition.  The public reader does not publish a unique numerical smooth
chi-family or the hidden grid origin.

This module therefore supplies an explicitly repository-autonomous C-infinity
realization of that source contract.  It produces concrete partition values
and source-normalized D_r/D_z derivatives for the existing Agent-2 localized
complete-curl family, but it is NOT a recovered Kokuno cutoff family.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_partition_family import KokunoSourcePartitionFamilyContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-compatible-autonomous-partition-v1"

_SOURCE_FORMULAS = {
    "dyadic_squared_partition": "sum_ell chi_ell(q)^2=1",
    "product_grid_squared_partition": "chi_{ell,a} is a product-grid squared partition in three slow chart coordinates",
    "band_mesh": "S_*^(-3)=ell^(-6) on the active bands",
    "slow_label": "beta=(ell,a)",
    "rectangle_label": "gamma=(ell,a,sigma), with sigma assigned after the beta partition",
    "combined_weight": "eta_beta=chi_ell(q)*chi_{ell,a}",
    "combined_squared_partition": "sum_beta eta_beta^2=1",
    "dyadic_overlap": "overlapping dyadic supports imply |ell-ell'|<=2",
    "band_threshold": "the interaction count is stated for ell>=ell_0>4",
}

_AUTONOMOUS_REALIZATION = {
    "bump": "rho(s)=exp(-1/(1-s^2)) for |s|<1 and 0 otherwise",
    "dyadic_coordinate": "y=-log_2(q), with raw band bump rho(y-ell) normalized in l2",
    "grid_origin": "one fixed three-coordinate chart origin supplied by this repository realization",
    "grid_bumps": "tensor-product rho((slow_j-origin_j)/ell^(-6)-a_j), normalized in l2 separately for each ell",
    "active_labels": "finite union of floor/floor+1 supports; beta keeps each band-specific grid label once",
    "derivatives": "analytic chain/product rule followed by exact derivative of l2 normalization",
}

_TRUTH_BOUNDARY = {
    "source_squared_partition_structure_preserved": True,
    "source_band_mesh_formula_used": True,
    "autonomous_concrete_partition_realization_executable": True,
    "source_normalized_partition_derivatives_executable": True,
    "existing_partition_family_guard_reused": True,
    "concrete_source_partition_bumps_reconstructed": False,
    "source_hidden_grid_origin_recovered": False,
    "source_actual_partition_labels_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "positive_order_background_corrections_included": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_real(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _rho_and_prime(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s = np.asarray(s, dtype=float)
    rho = np.zeros_like(s)
    prime = np.zeros_like(s)
    mask = np.abs(s) < 1.0
    if np.any(mask):
        sm = s[mask]
        den = 1.0 - sm * sm
        rm = np.exp(-1.0 / den)
        rho[mask] = rm
        prime[mask] = rm * (-2.0 * sm) / (den * den)
    return rho, prime


def _normalize_squared(
    raw: np.ndarray,
    D_r_raw: np.ndarray,
    D_z_raw: np.ndarray,
    *,
    name: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    norm2 = np.sum(raw * raw, axis=-1)
    if np.any(norm2 <= 0.0) or not np.all(np.isfinite(norm2)):
        raise ValueError(f"{name} autonomous raw partition has no active support")
    norm = np.sqrt(norm2)
    dot_r = np.sum(raw * D_r_raw, axis=-1)
    dot_z = np.sum(raw * D_z_raw, axis=-1)
    chi = raw / norm[..., None]
    D_r_chi = D_r_raw / norm[..., None] - raw * (dot_r / (norm2 * norm))[..., None]
    D_z_chi = D_z_raw / norm[..., None] - raw * (dot_z / (norm2 * norm))[..., None]
    return chi, D_r_chi, D_z_chi


def _active_integers(values: np.ndarray) -> tuple[int, ...]:
    active: set[int] = set()
    for value in np.ravel(np.asarray(values, dtype=float)):
        base = math.floor(float(value))
        active.add(base)
        active.add(base + 1)
    return tuple(sorted(active))


def _active_grid_labels(scaled: np.ndarray) -> tuple[tuple[int, int, int], ...]:
    rows = np.asarray(scaled, dtype=float).reshape(-1, 3)
    active: set[tuple[int, int, int]] = set()
    for row in rows:
        options = []
        for value in row:
            base = math.floor(float(value))
            options.append((base, base + 1))
        active.update(tuple(int(x) for x in item) for item in itertools.product(*options))
    return tuple(sorted(active))


@dataclass(frozen=True)
class KokunoSourceCompatiblePartitionRealization:
    """Autonomous smooth partition satisfying the corrected-reader structure.

    ``D_r_q``/``D_z_q`` and slow-coordinate derivatives must use the same
    source-normalized derivative convention as the complete-curl lane.  The
    realization then applies the chain rule through only its autonomous cutoffs.
    """

    ell_min: int = 5
    grid_origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    partition_atol: float = 2.0e-12
    derivative_rtol: float = 2.0e-10

    def __post_init__(self) -> None:
        if isinstance(self.ell_min, bool) or not isinstance(self.ell_min, (int, np.integer)):
            raise ValueError("ell_min must be an integer autonomous band threshold")
        ell_min = int(self.ell_min)
        if not (5 <= ell_min <= 10_000):
            raise ValueError("ell_min must satisfy the source interaction range ell_min>=5")
        try:
            origin = tuple(float(x) for x in self.grid_origin)
        except TypeError as exc:
            raise ValueError("grid_origin must be a finite three-coordinate tuple") from exc
        if len(origin) != 3 or not np.all(np.isfinite(origin)):
            raise ValueError("grid_origin must contain exactly three finite coordinates")
        partition_atol = float(self.partition_atol)
        derivative_rtol = float(self.derivative_rtol)
        if not np.isfinite(partition_atol) or not (0.0 < partition_atol <= 1.0e-8):
            raise ValueError("partition_atol must lie in (0,1e-8]")
        if not np.isfinite(derivative_rtol) or not (0.0 < derivative_rtol <= 1.0e-6):
            raise ValueError("derivative_rtol must lie in (0,1e-6]")
        object.__setattr__(self, "ell_min", ell_min)
        object.__setattr__(self, "grid_origin", origin)
        object.__setattr__(self, "partition_atol", partition_atol)
        object.__setattr__(self, "derivative_rtol", derivative_rtol)

    @property
    def partition_contract(self) -> KokunoSourcePartitionFamilyContract:
        return KokunoSourcePartitionFamilyContract(
            partition_atol=self.partition_atol,
            derivative_rtol=self.derivative_rtol,
        )

    @staticmethod
    def mesh_for_ell(ell: int) -> float:
        ell = int(ell)
        if ell <= 0:
            raise ValueError("ell must be positive for the source band mesh ell^-6")
        return float(ell) ** -6

    def _validate_inputs(
        self,
        q: Any,
        D_r_q: Any,
        D_z_q: Any,
        slow_coordinates: Any,
        D_r_slow_coordinates: Any,
        D_z_slow_coordinates: Any,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        q = _finite_real(q, "q")
        D_r_q = _finite_real(D_r_q, "D_r_q")
        D_z_q = _finite_real(D_z_q, "D_z_q")
        slow = _finite_real(slow_coordinates, "slow_coordinates")
        D_r_slow = _finite_real(D_r_slow_coordinates, "D_r_slow_coordinates")
        D_z_slow = _finite_real(D_z_slow_coordinates, "D_z_slow_coordinates")
        if D_r_q.shape != q.shape or D_z_q.shape != q.shape:
            raise ValueError("D_r_q and D_z_q must have exactly the q shape")
        expected_slow_shape = q.shape + (3,)
        if slow.shape != expected_slow_shape:
            raise ValueError("slow_coordinates must have q.shape+(3,)")
        if D_r_slow.shape != expected_slow_shape or D_z_slow.shape != expected_slow_shape:
            raise ValueError("slow-coordinate derivatives must have the slow_coordinates shape")
        if np.any(q <= 0.0) or np.any(q > 1.0):
            raise ValueError("q must lie in the source construction range (0,1]")
        y = -np.log2(q)
        if np.any(np.floor(y) < self.ell_min):
            raise ValueError("q reaches bands below autonomous ell_min; this realization fails closed")
        return q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow

    def evaluate(
        self,
        q: Any,
        D_r_q: Any,
        D_z_q: Any,
        slow_coordinates: Any,
        D_r_slow_coordinates: Any,
        D_z_slow_coordinates: Any,
    ) -> dict[str, Any]:
        """Return concrete beta weights and source-normalized first derivatives."""
        q, D_r_q, D_z_q, slow, D_r_slow, D_z_slow = self._validate_inputs(
            q, D_r_q, D_z_q, slow_coordinates, D_r_slow_coordinates, D_z_slow_coordinates
        )
        sample_shape = q.shape
        y = -np.log2(q)
        D_r_y = -D_r_q / (q * math.log(2.0))
        D_z_y = -D_z_q / (q * math.log(2.0))
        ell_labels = _active_integers(y)
        if min(ell_labels) < self.ell_min:
            raise ValueError("active dyadic support crosses below autonomous ell_min")
        ell_array = np.asarray(ell_labels, dtype=float)
        s_ell = y[..., None] - ell_array
        raw_ell, raw_ell_prime = _rho_and_prime(s_ell)
        D_r_raw_ell = raw_ell_prime * D_r_y[..., None]
        D_z_raw_ell = raw_ell_prime * D_z_y[..., None]
        chi_ell, D_r_chi_ell, D_z_chi_ell = _normalize_squared(
            raw_ell, D_r_raw_ell, D_z_raw_ell, name="dyadic"
        )

        origin = np.asarray(self.grid_origin, dtype=float)
        eta_parts: list[np.ndarray] = []
        D_r_eta_parts: list[np.ndarray] = []
        D_z_eta_parts: list[np.ndarray] = []
        beta_labels: list[tuple[int, tuple[int, int, int]]] = []
        grid_labels_by_ell: dict[int, tuple[tuple[int, int, int], ...]] = {}
        mesh_by_ell: dict[int, float] = {}

        for ell_index, ell in enumerate(ell_labels):
            mesh = self.mesh_for_ell(ell)
            mesh_by_ell[int(ell)] = mesh
            scaled = (slow - origin) / mesh
            D_r_scaled = D_r_slow / mesh
            D_z_scaled = D_z_slow / mesh
            a_labels = _active_grid_labels(scaled)
            grid_labels_by_ell[int(ell)] = a_labels
            a_array = np.asarray(a_labels, dtype=float)
            n_a = len(a_labels)
            raw_a = np.ones(sample_shape + (n_a,), dtype=float)
            D_r_raw_a = np.zeros_like(raw_a)
            D_z_raw_a = np.zeros_like(raw_a)
            rhos: list[np.ndarray] = []
            primes: list[np.ndarray] = []
            for j in range(3):
                sj = scaled[..., j, None] - a_array[:, j]
                rj, pj = _rho_and_prime(sj)
                rhos.append(rj)
                primes.append(pj)
                raw_a *= rj
            for j in range(3):
                other = np.ones_like(raw_a)
                for k in range(3):
                    if k != j:
                        other *= rhos[k]
                D_r_raw_a += primes[j] * D_r_scaled[..., j, None] * other
                D_z_raw_a += primes[j] * D_z_scaled[..., j, None] * other
            chi_a, D_r_chi_a, D_z_chi_a = _normalize_squared(
                raw_a, D_r_raw_a, D_z_raw_a, name=f"product-grid ell={ell}"
            )

            band = chi_ell[..., ell_index, None]
            D_r_band = D_r_chi_ell[..., ell_index, None]
            D_z_band = D_z_chi_ell[..., ell_index, None]
            eta_parts.append(band * chi_a)
            D_r_eta_parts.append(D_r_band * chi_a + band * D_r_chi_a)
            D_z_eta_parts.append(D_z_band * chi_a + band * D_z_chi_a)
            beta_labels.extend((int(ell), tuple(int(x) for x in a)) for a in a_labels)

        eta = np.concatenate(eta_parts, axis=-1)
        D_r_eta = np.concatenate(D_r_eta_parts, axis=-1)
        D_z_eta = np.concatenate(D_z_eta_parts, axis=-1)
        labels = tuple(beta_labels)
        checked = self.partition_contract.validate_family(eta, D_r_eta, D_z_eta, labels)
        return {
            **checked,
            "ell_labels": ell_labels,
            "chi_ell": chi_ell,
            "D_r_chi_ell": D_r_chi_ell,
            "D_z_chi_ell": D_z_chi_ell,
            "grid_labels_by_ell": grid_labels_by_ell,
            "mesh_by_ell": mesh_by_ell,
            "autonomous_realization": True,
            "source_actual_partition_recovered": False,
        }

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "slow dyadic/product-grid squared partition before localized oscillatory curls",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "autonomous_realization": dict(_AUTONOMOUS_REALIZATION),
            "parameters": {
                "ell_min": self.ell_min,
                "grid_origin": list(self.grid_origin),
                "partition_atol": self.partition_atol,
                "derivative_rtol": self.derivative_rtol,
                "origin": (
                    "repository-autonomous diagnostic realization; the smooth bump, dyadic placement, "
                    "first admissible band, grid origin and selected finite labels are not recovered source data"
                ),
            },
            "caller_contract": {
                "required": [
                    "q whose active dyadic bands stay at ell>=ell_min",
                    "source-normalized D_r q and D_z q",
                    "three slow chart coordinates",
                    "source-normalized D_r/D_z derivatives of those coordinates",
                ],
                "output": "eta_beta, D_r eta_beta, D_z eta_beta, unique beta=(ell,a) labels",
                "important": (
                    "these concrete values preserve the public source partition identities and ell^-6 mesh, "
                    "but they are not the unreleased Kokuno chi-family, hidden grid origin, or actual labels"
                ),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceCompatiblePartitionRealization":
        if not isinstance(payload, dict):
            raise ValueError("source-compatible partition payload must be an object")
        required = {
            "schema", "source", "autonomous_realization", "parameters", "caller_contract", "truth_boundary"
        }
        if set(payload) - {"sha256"} != required:
            raise ValueError("source-compatible partition payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source-compatible partition schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {
            "ell_min", "grid_origin", "partition_atol", "derivative_rtol", "origin"
        }:
            raise ValueError("source-compatible partition parameters changed")
        obj = cls(
            ell_min=params["ell_min"],
            grid_origin=tuple(params["grid_origin"]),
            partition_atol=params["partition_atol"],
            derivative_rtol=params["derivative_rtol"],
        )
        expected = obj.to_payload()
        for key in required:
            if payload[key] != expected[key]:
                raise ValueError(f"source-compatible partition {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("source-compatible partition sha256 mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceCompatiblePartitionRealization":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
