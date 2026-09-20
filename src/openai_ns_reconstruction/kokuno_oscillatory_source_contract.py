"""Typed provenance contract for the Kokuno oscillatory/complete-curl lane.

This module deliberately contains no velocity fitting, pressure, forcing, residual
optimization, or derivative tuning.  Its purpose is narrower: make the boundary
between the corrected 2026-09-09 Kokuno source statements and this repository's
independent executable realization machine-checkable.

The corrected reconstruction motivates the localized complete-curl structure.
The repository-specific public-z pullback, autonomous parameter choices,
numerical derivative diagnostics, and artifact hashes are *not* asserted to be
paper-exact hidden data.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SOURCE_DATE = "2026-09-09"
SOURCE_PATH = "KokunoYumeto updated/output.pdf"
SOURCE_VERSION = "corrected-208-page-reconstruction"
PARENT_EXACT_HEAD = "946b4d13f63ae6761f02e9e66d4c3c7ed485ffa6"


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class KokunoOscillatorySourceContract:
    """Immutable source/realization boundary for Agent 2 oscillatory fields."""

    source_date: str = SOURCE_DATE
    source_path: str = SOURCE_PATH
    source_version: str = SOURCE_VERSION
    parent_exact_head: str = PARENT_EXACT_HEAD

    # Source-bound mathematical organization recorded from the corrected reader.
    source_sections: tuple[str, ...] = ("PA.2 complete curls", "PA.3 localized waves")
    complete_curl_identity: str = (
        "curl(chi_k B^(k)) = chi_k curl(B^(k)) + grad(chi_k) x B^(k)"
    )
    complete_curl_terms: tuple[str, ...] = (
        "chi_k curl(B^(k))",
        "grad(chi_k) x B^(k)",
    )
    cutoff_gradient_term_required: bool = True
    localized_base_field_components: tuple[str, ...] = (
        "V1^beta",
        "V2^beta",
        "eps_beta^-1 V3^beta",
    )
    source_phase: str = (
        "Theta_beta = beta z / eps^2 + (beta/r) sin(beta alpha - phi_beta)"
    )
    source_normalized_axial_derivative: str = "D_z = eps_beta partial_z"
    source_structures: tuple[str, ...] = (
        "Varkappa radial envelope",
        "sigma_tilde chart",
        "Theta_beta phase",
        "k_beta phase vector",
        "Phi_beta transverse one-form/frame",
    )

    # Repository realization facts.  These are intentionally segregated from source facts.
    repository_public_z_pullback: str = "Z_beta = eps_beta z"
    repository_parameter_status: str = "bounded autonomous realization; not recovered hidden source data"
    repository_derivative_status: str = "FD diagnostics are repository numerical approximations"
    repository_artifact_status: str = "checksum/manifest semantics are repository infrastructure"
    paper_exact: bool = False
    pde_validated: bool = False

    def source_bound_payload(self) -> dict[str, Any]:
        return {
            "source_date": self.source_date,
            "source_path": self.source_path,
            "source_version": self.source_version,
            "source_sections": list(self.source_sections),
            "complete_curl_identity": self.complete_curl_identity,
            "complete_curl_terms": list(self.complete_curl_terms),
            "cutoff_gradient_term_required": self.cutoff_gradient_term_required,
            "localized_base_field_components": list(self.localized_base_field_components),
            "source_phase": self.source_phase,
            "source_normalized_axial_derivative": self.source_normalized_axial_derivative,
            "source_structures": list(self.source_structures),
        }

    def repository_realization_payload(self) -> dict[str, Any]:
        return {
            "parent_exact_head": self.parent_exact_head,
            "repository_public_z_pullback": self.repository_public_z_pullback,
            "repository_parameter_status": self.repository_parameter_status,
            "repository_derivative_status": self.repository_derivative_status,
            "repository_artifact_status": self.repository_artifact_status,
            "paper_exact": self.paper_exact,
            "pde_validated": self.pde_validated,
        }

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": "kokuno-agent2-oscillatory-source-contract-v1",
            "source_bound": self.source_bound_payload(),
            "repository_realization": self.repository_realization_payload(),
            "scientific_boundary": {
                "claim": "Kokuno-inspired localized complete-curl repository realization",
                "not_claimed": [
                    "paper-exact hidden oscillatory data",
                    "global corrected leading field",
                    "matched pressure",
                    "restricted forcing closure",
                    "held-out Navier-Stokes validation",
                ],
            },
        }

    @property
    def contract_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def validate(self) -> None:
        if self.source_date != SOURCE_DATE:
            raise ValueError("source date does not match the corrected 2026-09-09 reader")
        if self.source_path != SOURCE_PATH:
            raise ValueError("source path does not identify the corrected Kokuno reconstruction")
        if self.source_version != SOURCE_VERSION:
            raise ValueError("source version is not the corrected 208-page reconstruction")
        if not self.cutoff_gradient_term_required:
            raise ValueError("complete curl must retain grad(chi_k) x B^(k)")
        if self.complete_curl_terms != (
            "chi_k curl(B^(k))",
            "grad(chi_k) x B^(k)",
        ):
            raise ValueError("complete-curl term contract changed or lost the cutoff-gradient term")
        if self.paper_exact:
            raise ValueError("repository realization must not be labeled paper-exact")
        if self.pde_validated:
            raise ValueError("source/provenance contract cannot certify PDE validation")

    def manifest(self) -> dict[str, Any]:
        self.validate()
        payload = self.semantic_payload()
        return {
            "payload": payload,
            "contract_sha256": _sha256(payload),
        }

    def save_manifest(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.manifest(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_manifest(cls, path: str | Path) -> "KokunoOscillatorySourceContract":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "contract_sha256"}:
            raise ValueError("invalid source-contract manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict):
            raise ValueError("invalid source-contract payload")
        if _sha256(payload) != raw["contract_sha256"]:
            raise ValueError("source-contract manifest checksum mismatch")

        source = payload.get("source_bound")
        realization = payload.get("repository_realization")
        if not isinstance(source, dict) or not isinstance(realization, dict):
            raise ValueError("source/realization partitions are missing")

        contract = cls(
            source_date=source.get("source_date", ""),
            source_path=source.get("source_path", ""),
            source_version=source.get("source_version", ""),
            parent_exact_head=realization.get("parent_exact_head", ""),
            source_sections=tuple(source.get("source_sections", ())),
            complete_curl_identity=source.get("complete_curl_identity", ""),
            complete_curl_terms=tuple(source.get("complete_curl_terms", ())),
            cutoff_gradient_term_required=source.get("cutoff_gradient_term_required", False),
            localized_base_field_components=tuple(source.get("localized_base_field_components", ())),
            source_phase=source.get("source_phase", ""),
            source_normalized_axial_derivative=source.get("source_normalized_axial_derivative", ""),
            source_structures=tuple(source.get("source_structures", ())),
            repository_public_z_pullback=realization.get("repository_public_z_pullback", ""),
            repository_parameter_status=realization.get("repository_parameter_status", ""),
            repository_derivative_status=realization.get("repository_derivative_status", ""),
            repository_artifact_status=realization.get("repository_artifact_status", ""),
            paper_exact=realization.get("paper_exact", True),
            pde_validated=realization.get("pde_validated", True),
        )
        contract.validate()
        if contract.semantic_payload() != payload:
            raise ValueError("manifest semantics do not match the typed contract")
        return contract


def default_kokuno_oscillatory_source_contract() -> KokunoOscillatorySourceContract:
    contract = KokunoOscillatorySourceContract()
    contract.validate()
    return contract
