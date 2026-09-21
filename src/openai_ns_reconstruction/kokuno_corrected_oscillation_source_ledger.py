"""Exact-source ledger for the corrected Kokuno oscillation / complete-curl lane.

This module is deliberately provenance-only.  It records formulas and source
identities that are explicit in the corrected 2026-09-09, 208-page public
reconstruction, while keeping the current repository oscillatory runtime on the
other side of a hard source-to-realization boundary.

Nothing here changes a velocity coefficient, support, pressure, forcing,
viscosity, residual, optimizer, or scientific threshold.  In particular, this
ledger does not assert that the repository's autonomous oscillatory phase,
support, or parameter values are paper-exact.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


PARENT_A2_HEAD = "e36d9da4b7f037e998f5b1658f8c0ea291a76b80"
PARENT_A2_SOURCE_BLOB = "4f1dc0566c041b9de20f18b4ae9c49d9fd93c58d"

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_DATE = "2026-09-09"
SOURCE_VERSION = "corrected-208-page-reconstruction"
SOURCE_WORKBENCH_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_WORKBENCH_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_CHECKS_PATH = "navier-stokes/navier_stokes_checks.json"
SOURCE_CHECKS_BLOB = "03c5063968de0bbc3b9a34462494859df677fff7"
SOURCE_BUNDLE_PATH = "navier-stokes/navier_stokes_source_bundle.zip"
SOURCE_BUNDLE_BLOB = "b3cf96edff06e39e066a4a6b9522362975d3fbb1"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
SOURCE_PDF_SHA256 = "242fe7f83feab43b91915595230b5e93e9f344eb654e1edcb65e9669a9cd5878"
SOURCE_CHECKS_RECORDED_JSON_SHA256 = (
    "386d3fef8edcf2665ecd5daa2944269f414ced12c1afe89d9d7fc50973fbd76a"
)
SOURCE_READER_PAGES = 208

OSCILLATION_COMPONENT_SHA256: dict[str, str] = {
    "proof_sources/oscillations/01_auxiliary_geometry.md": (
        "c3304b492b9fc08f905dc7a09ce145d7cd06062df4d3d718629b9f547b5ed0df"
    ),
    "proof_sources/oscillations/02_pulse_equation.md": (
        "d4d5d37041ec829a45d28b350f14c4a0793393095f9acbd2b09642228a88c2ec"
    ),
    "proof_sources/oscillations/03_covariance_curl_and_tails.md": (
        "71432f2fffa382cfe4b8c325c51ee4a8caa3c2578cef4fc308f0db2497ca3474"
    ),
    "proof_sources/oscillations/04_admissible_cone.md": (
        "264806a5428ea6bafa113970eea4161427e1ed2df7839999905c89336bba2f6d"
    ),
    "proof_sources/oscillations/READ_SCOPE.md": (
        "e56d320f5d11fc6ccadbc40cb8eb44d9e9cc7de8e8e96a1c2a9ccbcf822a7149"
    ),
    "proof_sources/oscillations/STATUS.md": (
        "7a44930041b87f8f9a66bead2606611abcc9238d42a4316c014fb200feda6ea3"
    ),
    "proof_sources/oscillations/complete_derivations.md": (
        "759ea68a3774832763561f98b579208859c45320c68dbcd1da74f7563bb1203d"
    ),
    "proof_sources/oscillations/exact_checks.py": (
        "a44ba49990c3b604c883d5b4e1726cb33068867cd85b23f3a2c22847e4ce8184"
    ),
    "proof_sources/oscillations/exact_checks_results.json": (
        "30d7a4d934e984a77cfcb04188d72c35391c46e38bcbf25a2024e38bcb783a12"
    ),
    "proof_sources/oscillations/oscillations_body.tex": (
        "3fd5c61de39b6f65a581ead5b29c741c2e0f46fb30f62e582dc24d93bbe83615"
    ),
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class KokunoCorrectedOscillationSourceLedger:
    """Machine-checkable corrected-source facts and a fail-closed truth boundary."""

    source_repository: str = SOURCE_REPOSITORY
    source_commit: str = SOURCE_COMMIT
    source_date: str = SOURCE_DATE
    source_version: str = SOURCE_VERSION
    source_workbench_path: str = SOURCE_WORKBENCH_PATH
    source_workbench_blob: str = SOURCE_WORKBENCH_BLOB
    source_checks_path: str = SOURCE_CHECKS_PATH
    source_checks_blob: str = SOURCE_CHECKS_BLOB
    source_bundle_path: str = SOURCE_BUNDLE_PATH
    source_bundle_blob: str = SOURCE_BUNDLE_BLOB
    source_bundle_sha256: str = SOURCE_BUNDLE_SHA256
    source_pdf_sha256: str = SOURCE_PDF_SHA256
    source_checks_recorded_json_sha256: str = SOURCE_CHECKS_RECORDED_JSON_SHA256
    source_reader_pages: int = SOURCE_READER_PAGES
    parent_a2_head: str = PARENT_A2_HEAD
    parent_a2_source_blob: str = PARENT_A2_SOURCE_BLOB

    # Exact public-source chart/scaling structure used by the oscillation lane.
    chart_scaling: tuple[str, ...] = (
        "A = 1/2 + h",
        "D = 1/2 - h",
        "0 < h < 1/100",
        "Q = 2^(-ell)",
        "epsilon = Q^h",
        "S_* = ell^2",
        "R = r / Q^(1/2)",
        "Z = z / Q^D",
        "T = (1-t) / Q",
    )

    # Required background inputs.  Their construction belongs to the background lane.
    imported_background_requirement: str = (
        "On the fixed enlarged annulus, the imported background has chart radial "
        "component b=O(epsilon), while tangential V,G differ from their leading "
        "profiles by O(epsilon^2), with the stated fixed slow-derivative bounds."
    )
    required_background_inputs: tuple[str, ...] = (
        "V",
        "G",
        "R",
        "F = V/R",
        "F_R",
        "G_R",
        "g = (R F_R, G_R)",
        "representative leading shear g_0",
        "representative frame vectors N and K",
        "pulse coordinate v and length L_s",
        "B_s, sigma, u_*, x_0",
    )

    # Corrected full phase/frame formulas.
    carrier: str = "k = ceil(epsilon^(-1/2))"
    tangential_carrier_selection: str = (
        "(tilde p/R_0, p_z) = B_s (K - sigma u_* g_0/(L_s |g_0|^2)); "
        "k p is the specified nearest nonzero integer to k tilde p"
    )
    phase_hamiltonian: str = "H_Phi = p F + p_z G"
    phase: str = "Phi = p theta + p_z Z/epsilon + x_0 R - v H_Phi"
    phase_vector: tuple[str, str, str] = (
        "(n_Phi)_r = x_0 - v (H_Phi)_R",
        "(n_Phi)_theta = p/R",
        "(n_Phi)_z = p_z - epsilon v (H_Phi)_Z",
    )

    # Constrained oscillatory amplitude and its source pressure relation.
    amplitude_matrix: tuple[tuple[str, str, str], ...] = (
        ("0", "-2F", "0"),
        ("2F + R F_R", "0", "0"),
        ("G_R", "0", "0"),
    )
    damping: str = "d = epsilon k^2 |n_Phi|^2"
    constrained_amplitude_equation: str = (
        "t_m' + K t_m + m^2 d t_m + i k m n_Phi pi_m = -f_m"
    )
    constrained_amplitude_transversality: str = "n_Phi dot t_m = 0"
    source_pressure_relation: str = (
        "pi_m = (i/(k m)) (n_Phi dot K t_m - n_Phi' dot t_m + "
        "n_Phi dot f_m) / |n_Phi|^2"
    )

    # Exact normalized cylindrical complete curl.
    normalized_curl: tuple[str, str, str] = (
        "(curl_* A)_r = R^(-1) partial_theta A_z - D_z A_theta",
        "(curl_* A)_theta = D_z A_r - D_r A_z",
        "(curl_* A)_z = (D_r + R^(-1)) A_theta - R^(-1) partial_theta A_r",
    )
    vector_potential_coefficient: str = (
        "C_m = i (n_Phi cross t_m) / (k m |n_Phi|^2)"
    )
    harmonic_vector_potential: str = "A_m = C_m exp(i k m Phi)"
    exponential_curl_contribution: str = (
        "i k m (n_Phi cross C_m) = t_m when n_Phi dot t_m = 0"
    )
    exact_curl_remainder: tuple[str, str, str] = (
        "(r_m)_r = -D_z (C_m)_theta",
        "(r_m)_theta = D_z (C_m)_r - D_r (C_m)_z",
        "(r_m)_z = (D_r + R^(-1)) (C_m)_theta",
    )
    complete_amplitude: str = "a_m = t_m + r_m"
    normalized_divergence: str = (
        "div_* a = (D_r + R^(-1)) a_r + R^(-1) partial_theta a_theta + D_z a_z"
    )
    exact_curl_divergence_identity: str = "div_*(curl_* A) = 0"
    complete_amplitude_asserted_transverse: bool = False
    retained_longitudinal_identity: str = (
        "i k m (n_Phi dot a_m) = -((D_r + R^(-1))(a_m)_r + D_z(a_m)_z)"
    )
    retained_longitudinal_source: str = "n_Phi dot a_m = n_Phi dot r_m"
    physical_vector_potential_scaling: str = "A_phys = Q^(1/2-A) A_*"
    physical_pressure_scaling: str = "p_phys = Q^(-2A) p_*"
    real_wave_pairing: str = (
        "Real waves use conjugate harmonic pairs; for t cos(k Phi), velocity coefficients "
        "are t/2 for m=+1 and m=-1."
    )

    # Source support/cutoff structure.  This is not the repository cos^8 support.
    source_support_envelope: str = (
        "zeta = exp(-a_a/log^2(X/X_a) - a_b/log^2(X_b/X)), extended by zero "
        "outside X_a < X < X_b"
    )
    source_wave_support: str = (
        "Wave coefficients have prescribed closed slow/transverse/shell supports, "
        "sqrt(zeta) edge weight, pulse envelope P(v), and smooth zero extension."
    )
    cutoff_definition: str = "hat t_m = psi t_m; hat pi_m = psi pi_m"
    exact_cutoff_residual: str = (
        "hat t_m' + K hat t_m + m^2 d hat t_m + i k m n_Phi hat pi_m + f_m "
        "= (1-psi) f_m + psi' t_m"
    )
    cutoff_tail_location: str = "|v - L_s/2| >= L_s/5"
    cutoff_tail_status: str = (
        "At each fixed stage, cutoff-tail terms obey the proved Gaussian/exponential "
        "smallness and retain smooth support extension under fixed derivatives."
    )

    # Current repository/source mapping status.  These must remain hard false here.
    current_runtime_phase_source_exact: bool = False
    current_runtime_support_source_exact: bool = False
    current_runtime_complete_curl_source_equivalence_verified: bool = False
    source_to_runtime_parameter_map_complete: bool = False
    current_runtime_paper_exact: bool = False
    matched_pressure_materialized: bool = False
    restricted_forcing_materialized: bool = False
    complete_ns_residual_assessed: bool = False
    pde_validated: bool = False

    def upstream_identity_payload(self) -> dict[str, Any]:
        return {
            "repository": self.source_repository,
            "commit": self.source_commit,
            "date": self.source_date,
            "version": self.source_version,
            "reader_pages": self.source_reader_pages,
            "workbench_path": self.source_workbench_path,
            "workbench_blob": self.source_workbench_blob,
            "checks_path": self.source_checks_path,
            "checks_blob": self.source_checks_blob,
            "source_bundle_path": self.source_bundle_path,
            "source_bundle_blob": self.source_bundle_blob,
            "source_bundle_sha256": self.source_bundle_sha256,
            "pdf_sha256": self.source_pdf_sha256,
            "checks_recorded_json_sha256": self.source_checks_recorded_json_sha256,
            "oscillation_component_sha256": dict(OSCILLATION_COMPONENT_SHA256),
        }

    def formula_payload(self) -> dict[str, Any]:
        return {
            "chart_scaling": list(self.chart_scaling),
            "background": {
                "imported_requirement": self.imported_background_requirement,
                "required_inputs": list(self.required_background_inputs),
            },
            "phase": {
                "carrier": self.carrier,
                "tangential_carrier_selection": self.tangential_carrier_selection,
                "H_Phi": self.phase_hamiltonian,
                "Phi": self.phase,
                "n_Phi": list(self.phase_vector),
            },
            "constrained_amplitude": {
                "K_matrix": [list(row) for row in self.amplitude_matrix],
                "damping": self.damping,
                "equation": self.constrained_amplitude_equation,
                "transversality": self.constrained_amplitude_transversality,
                "source_pressure_relation": self.source_pressure_relation,
            },
            "complete_curl": {
                "normalized_curl": list(self.normalized_curl),
                "C_m": self.vector_potential_coefficient,
                "A_m": self.harmonic_vector_potential,
                "exponential_contribution": self.exponential_curl_contribution,
                "exact_remainder": list(self.exact_curl_remainder),
                "complete_amplitude": self.complete_amplitude,
                "normalized_divergence": self.normalized_divergence,
                "exact_divergence_identity": self.exact_curl_divergence_identity,
                "complete_amplitude_asserted_transverse": self.complete_amplitude_asserted_transverse,
                "retained_longitudinal_identity": self.retained_longitudinal_identity,
                "retained_longitudinal_source": self.retained_longitudinal_source,
                "physical_vector_potential_scaling": self.physical_vector_potential_scaling,
                "physical_pressure_scaling": self.physical_pressure_scaling,
                "real_wave_pairing": self.real_wave_pairing,
            },
            "support_and_cutoff": {
                "source_support_envelope": self.source_support_envelope,
                "source_wave_support": self.source_wave_support,
                "cutoff_definition": self.cutoff_definition,
                "exact_cutoff_residual": self.exact_cutoff_residual,
                "cutoff_tail_location": self.cutoff_tail_location,
                "cutoff_tail_status": self.cutoff_tail_status,
            },
        }

    def truth_boundary_payload(self) -> dict[str, Any]:
        return {
            "source_formula_status": "confirmed_in_corrected_208p_workbench_and_source_bundle",
            "current_runtime_mapping_status": "not_verified",
            "current_runtime_phase_source_exact": self.current_runtime_phase_source_exact,
            "current_runtime_support_source_exact": self.current_runtime_support_source_exact,
            "current_runtime_complete_curl_source_equivalence_verified": (
                self.current_runtime_complete_curl_source_equivalence_verified
            ),
            "source_to_runtime_parameter_map_complete": self.source_to_runtime_parameter_map_complete,
            "current_runtime_paper_exact": self.current_runtime_paper_exact,
            "matched_pressure_materialized": self.matched_pressure_materialized,
            "restricted_forcing_materialized": self.restricted_forcing_materialized,
            "complete_ns_residual_assessed": self.complete_ns_residual_assessed,
            "pde_validated": self.pde_validated,
            "scope": (
                "provenance/formula ledger only; no velocity, parameter, support, pressure, "
                "forcing, residual, threshold, or sibling-lane implementation change"
            ),
        }

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": "kokuno-agent2-corrected-oscillation-source-ledger-v1",
            "parent_a2": {
                "head": self.parent_a2_head,
                "source_blob": self.parent_a2_source_blob,
            },
            "upstream_source": self.upstream_identity_payload(),
            "source_formulas": self.formula_payload(),
            "truth_boundary": self.truth_boundary_payload(),
        }

    @property
    def semantic_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def validate(self) -> None:
        exact_identity = {
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_date": SOURCE_DATE,
            "source_version": SOURCE_VERSION,
            "source_workbench_path": SOURCE_WORKBENCH_PATH,
            "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
            "source_checks_path": SOURCE_CHECKS_PATH,
            "source_checks_blob": SOURCE_CHECKS_BLOB,
            "source_bundle_path": SOURCE_BUNDLE_PATH,
            "source_bundle_blob": SOURCE_BUNDLE_BLOB,
            "source_bundle_sha256": SOURCE_BUNDLE_SHA256,
            "source_pdf_sha256": SOURCE_PDF_SHA256,
            "source_checks_recorded_json_sha256": SOURCE_CHECKS_RECORDED_JSON_SHA256,
            "source_reader_pages": SOURCE_READER_PAGES,
            "parent_a2_head": PARENT_A2_HEAD,
            "parent_a2_source_blob": PARENT_A2_SOURCE_BLOB,
        }
        for attr, expected in exact_identity.items():
            if getattr(self, attr) != expected:
                raise ValueError(f"source-ledger identity drift: {attr}")

        if self.phase != "Phi = p theta + p_z Z/epsilon + x_0 R - v H_Phi":
            raise ValueError("corrected full phase formula drift")
        if self.phase_vector != (
            "(n_Phi)_r = x_0 - v (H_Phi)_R",
            "(n_Phi)_theta = p/R",
            "(n_Phi)_z = p_z - epsilon v (H_Phi)_Z",
        ):
            raise ValueError("corrected full phase-vector formula drift")
        if self.constrained_amplitude_transversality != "n_Phi dot t_m = 0":
            raise ValueError("source transverse-amplitude constraint drift")
        if self.normalized_curl != (
            "(curl_* A)_r = R^(-1) partial_theta A_z - D_z A_theta",
            "(curl_* A)_theta = D_z A_r - D_r A_z",
            "(curl_* A)_z = (D_r + R^(-1)) A_theta - R^(-1) partial_theta A_r",
        ):
            raise ValueError("complete cylindrical curl formula drift")
        if self.exact_curl_remainder != (
            "(r_m)_r = -D_z (C_m)_theta",
            "(r_m)_theta = D_z (C_m)_r - D_r (C_m)_z",
            "(r_m)_z = (D_r + R^(-1)) (C_m)_theta",
        ):
            raise ValueError("complete-curl remainder formula drift")
        if self.complete_amplitude != "a_m = t_m + r_m":
            raise ValueError("complete-amplitude formula drift")
        if self.complete_amplitude_asserted_transverse:
            raise ValueError("corrected reader does not assert complete amplitude transversality")
        if self.retained_longitudinal_source != "n_Phi dot a_m = n_Phi dot r_m":
            raise ValueError("retained longitudinal component was lost")
        if set(self.required_background_inputs) != {
            "V",
            "G",
            "R",
            "F = V/R",
            "F_R",
            "G_R",
            "g = (R F_R, G_R)",
            "representative leading shear g_0",
            "representative frame vectors N and K",
            "pulse coordinate v and length L_s",
            "B_s, sigma, u_*, x_0",
        }:
            raise ValueError("required source background-input set drift")
        if dict(OSCILLATION_COMPONENT_SHA256) != self.upstream_identity_payload()[
            "oscillation_component_sha256"
        ]:
            raise ValueError("source-bundle oscillation component hash drift")

        forbidden_promotions = {
            "current_runtime_phase_source_exact": self.current_runtime_phase_source_exact,
            "current_runtime_support_source_exact": self.current_runtime_support_source_exact,
            "current_runtime_complete_curl_source_equivalence_verified": (
                self.current_runtime_complete_curl_source_equivalence_verified
            ),
            "source_to_runtime_parameter_map_complete": self.source_to_runtime_parameter_map_complete,
            "current_runtime_paper_exact": self.current_runtime_paper_exact,
            "matched_pressure_materialized": self.matched_pressure_materialized,
            "restricted_forcing_materialized": self.restricted_forcing_materialized,
            "complete_ns_residual_assessed": self.complete_ns_residual_assessed,
            "pde_validated": self.pde_validated,
        }
        promoted = [name for name, value in forbidden_promotions.items() if value]
        if promoted:
            raise ValueError("source ledger cannot promote runtime/scientific state: " + ", ".join(promoted))

    def manifest(self) -> dict[str, Any]:
        self.validate()
        payload = self.semantic_payload()
        return {"payload": payload, "semantic_sha256": _sha256(payload)}

    def save_manifest(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.manifest(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_manifest(cls, path: str | Path) -> "KokunoCorrectedOscillationSourceLedger":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "semantic_sha256"}:
            raise ValueError("invalid corrected-source-ledger manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict):
            raise ValueError("invalid corrected-source-ledger payload")
        if _sha256(payload) != raw["semantic_sha256"]:
            raise ValueError("corrected-source-ledger checksum mismatch")

        parent = payload.get("parent_a2")
        source = payload.get("upstream_source")
        formulas = payload.get("source_formulas")
        truth = payload.get("truth_boundary")
        if not all(isinstance(item, dict) for item in (parent, source, formulas, truth)):
            raise ValueError("corrected-source-ledger partitions are missing")

        phase = formulas.get("phase", {})
        amplitude = formulas.get("constrained_amplitude", {})
        curl = formulas.get("complete_curl", {})
        support = formulas.get("support_and_cutoff", {})
        background = formulas.get("background", {})

        ledger = cls(
            source_repository=source.get("repository", ""),
            source_commit=source.get("commit", ""),
            source_date=source.get("date", ""),
            source_version=source.get("version", ""),
            source_reader_pages=source.get("reader_pages", -1),
            source_workbench_path=source.get("workbench_path", ""),
            source_workbench_blob=source.get("workbench_blob", ""),
            source_checks_path=source.get("checks_path", ""),
            source_checks_blob=source.get("checks_blob", ""),
            source_bundle_path=source.get("source_bundle_path", ""),
            source_bundle_blob=source.get("source_bundle_blob", ""),
            source_bundle_sha256=source.get("source_bundle_sha256", ""),
            source_pdf_sha256=source.get("pdf_sha256", ""),
            source_checks_recorded_json_sha256=source.get("checks_recorded_json_sha256", ""),
            parent_a2_head=parent.get("head", ""),
            parent_a2_source_blob=parent.get("source_blob", ""),
            chart_scaling=tuple(formulas.get("chart_scaling", ())),
            imported_background_requirement=background.get("imported_requirement", ""),
            required_background_inputs=tuple(background.get("required_inputs", ())),
            carrier=phase.get("carrier", ""),
            tangential_carrier_selection=phase.get("tangential_carrier_selection", ""),
            phase_hamiltonian=phase.get("H_Phi", ""),
            phase=phase.get("Phi", ""),
            phase_vector=tuple(phase.get("n_Phi", ())),
            amplitude_matrix=tuple(tuple(row) for row in amplitude.get("K_matrix", ())),
            damping=amplitude.get("damping", ""),
            constrained_amplitude_equation=amplitude.get("equation", ""),
            constrained_amplitude_transversality=amplitude.get("transversality", ""),
            source_pressure_relation=amplitude.get("source_pressure_relation", ""),
            normalized_curl=tuple(curl.get("normalized_curl", ())),
            vector_potential_coefficient=curl.get("C_m", ""),
            harmonic_vector_potential=curl.get("A_m", ""),
            exponential_curl_contribution=curl.get("exponential_contribution", ""),
            exact_curl_remainder=tuple(curl.get("exact_remainder", ())),
            complete_amplitude=curl.get("complete_amplitude", ""),
            normalized_divergence=curl.get("normalized_divergence", ""),
            exact_curl_divergence_identity=curl.get("exact_divergence_identity", ""),
            complete_amplitude_asserted_transverse=curl.get(
                "complete_amplitude_asserted_transverse", True
            ),
            retained_longitudinal_identity=curl.get("retained_longitudinal_identity", ""),
            retained_longitudinal_source=curl.get("retained_longitudinal_source", ""),
            physical_vector_potential_scaling=curl.get("physical_vector_potential_scaling", ""),
            physical_pressure_scaling=curl.get("physical_pressure_scaling", ""),
            real_wave_pairing=curl.get("real_wave_pairing", ""),
            source_support_envelope=support.get("source_support_envelope", ""),
            source_wave_support=support.get("source_wave_support", ""),
            cutoff_definition=support.get("cutoff_definition", ""),
            exact_cutoff_residual=support.get("exact_cutoff_residual", ""),
            cutoff_tail_location=support.get("cutoff_tail_location", ""),
            cutoff_tail_status=support.get("cutoff_tail_status", ""),
            current_runtime_phase_source_exact=truth.get("current_runtime_phase_source_exact", True),
            current_runtime_support_source_exact=truth.get("current_runtime_support_source_exact", True),
            current_runtime_complete_curl_source_equivalence_verified=truth.get(
                "current_runtime_complete_curl_source_equivalence_verified", True
            ),
            source_to_runtime_parameter_map_complete=truth.get(
                "source_to_runtime_parameter_map_complete", True
            ),
            current_runtime_paper_exact=truth.get("current_runtime_paper_exact", True),
            matched_pressure_materialized=truth.get("matched_pressure_materialized", True),
            restricted_forcing_materialized=truth.get("restricted_forcing_materialized", True),
            complete_ns_residual_assessed=truth.get("complete_ns_residual_assessed", True),
            pde_validated=truth.get("pde_validated", True),
        )
        ledger.validate()
        if ledger.semantic_payload() != payload:
            raise ValueError("manifest semantics do not match corrected source ledger")
        return ledger


def default_kokuno_corrected_oscillation_source_ledger() -> KokunoCorrectedOscillationSourceLedger:
    ledger = KokunoCorrectedOscillationSourceLedger()
    ledger.validate()
    return ledger
