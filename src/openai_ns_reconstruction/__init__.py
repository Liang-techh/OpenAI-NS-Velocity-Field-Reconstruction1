"""Public API for the auditable, still-partial OpenAI NS reconstruction."""
from .coordinates import (
    SimilarityPoint,
    similarity_coordinates,
    similarity_coordinates_from_tau,
    solve_q,
    solve_q_from_tau,
    coordinate_derivatives,
    weighted_profile_derivatives,
)
from .profiles import LeadingProfile, toy_gaussian_profile
from .velocity import (
    leading_velocity_cartesian,
    leading_velocity_cylindrical,
    leading_velocity_cartesian_from_tau,
    leading_velocity_cylindrical_from_tau,
    leading_velocity_from_tau,
    leading_pressure,
    leading_pressure_cartesian,
    blowup_probe,
)
from .status import construction_status
from .st054_snapshot import ST054Snapshot, available_st054_models, load_st054

__version__ = "0.2.0"

__all__ = [
    "SimilarityPoint", "similarity_coordinates", "similarity_coordinates_from_tau",
    "solve_q", "solve_q_from_tau", "coordinate_derivatives", "weighted_profile_derivatives",
    "LeadingProfile", "toy_gaussian_profile", "leading_velocity_cartesian",
    "leading_velocity_cylindrical", "leading_velocity_cartesian_from_tau",
    "leading_velocity_cylindrical_from_tau", "leading_velocity_from_tau",
    "leading_pressure", "leading_pressure_cartesian", "blowup_probe", "construction_status",
    "ST054Snapshot", "available_st054_models", "load_st054",
]
