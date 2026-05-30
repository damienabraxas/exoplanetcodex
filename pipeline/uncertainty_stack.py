"""
uncertainty_stack.py
====================
Formal uncertainty budget for all derived abundances.

Type A (random) uncertainties:
  - EW measurement uncertainty (σ_fit from profile fitting)
  - Line scatter: std(A(X)) / sqrt(N_lines)

Type B (systematic) uncertainties from stellar parameter sensitivity:
  - ΔA(X) / ΔTeff  (typically ±100 K)
  - ΔA(X) / Δlogg  (typically ±0.1 dex)
  - ΔA(X) / Δvturb (typically ±0.1 km/s)
  - ΔA(X) / Δ[Fe/H] (model grid interpolation error)

Total uncertainty: σ_tot = sqrt(σ_A² + σ_B1² + σ_B2² + σ_B3² + ...)

Linear issue: RYA-TBD
"""

# TODO: Implement uncertainty propagation
# Inputs:  abundance table from abundances_derive.py
#          stellar parameter uncertainties from params_stellar.py
# Outputs: per-element Type A + Type B uncertainty table
# Key decisions:
#   - Parameter step sizes for numerical derivatives
#   - Whether to propagate covariances between stellar parameters

from config.constants import PIPELINE, PATHS


def run(star_id: str) -> None:
    """Called by run_pipeline.py. TODO: implement uncertainty budget."""
    raise NotImplementedError("uncertainty_stack.run() is not yet implemented.")
