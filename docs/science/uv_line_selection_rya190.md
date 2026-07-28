# UV line selection + NLTE policy (RYA-190)

The canonical record of **which UV (FUV/NUV) transitions are usable for abundance work, by
what method, with what NLTE handling.** Implemented (not merely scoped) as the single source
of truth `pipeline/uv_line_selection.py`, which the UV pipeline consumes:

- **RYA-426** (UV conditioning) maps spectral windows to these diagnostics (`anchors_for_window_map`).
- **RYA-471** (HST UV loader/arm) wires `usable_diagnostics()` as the UV arm's diagnostic set
  and pairs each UV element with its **optical cross-check** (`optical_cross_check`) — the
  UV-C I vs optical-C I leg validation.

## The policy (cited)

| Species | λ (Å) | Frame | Regime | Verdict | Method | NLTE | Reason |
|--|--|--|--|--|--|--|--|
| C I | 1657.38 | vacuum | FUV | **USE** | synthesis | GRID_OWED (~+0.10) | C I UV multiplet, R~114k (E140H/E140M); strong but correctable NLTE |
| C I | 1930.90 | vacuum | FUV | **USE** | synthesis | GRID_OWED (~+0.10) | secondary C I UV multiplet |
| C I | 5052/5380 | air | optical | **USE** | EW/synth | GRID_OWED (~+0.04) | preferred OPTICAL C (HARPS); the cross-check leg for UV C I |
| O I | 1355.60 | vacuum | FUV | **USE** | synthesis | GRID_OWED (~+0.05) | semi-forbidden; cleaner than the resonance triplet |
| O I | 1302–1306 | vacuum | FUV | **DO NOT USE** | — | — | resonance triplet saturated in all FGK; no abundance info |
| N I | 1199.55 | vacuum | FUV | **DO NOT USE** | — | — | resonance, strong NLTE ~0.3–0.5 dex; unreliable |
| NH | ~3360 | air | NUV | **USE** | synthesis (molecular) | n/a | NH A–X band; molecular N. *Coverage gap*: just beyond Procyon STIS NUV max (~3160) |
| S I | 1473.99 | vacuum | FUV | **USE** | synthesis | GRID_OWED | UV S multiplet (~1425–1479, E140M); UV S access the optical barely touches |
| CH | ~4300 | air | optical | **PREFER ALT** | synthesis (molecular) | n/a | CH G-band; the standard optical C cross-check |

References: RYA-190 policy table; Amarsi et al. 2020 (A&A 642, A62) for C/N/O UV NLTE;
Lind et al. 2011 (A&A 528, A103) for Na optical NLTE; standard atomic positions.

## NLTE policy — grids, not scalars (the discipline)

The project applies NLTE/3D via **published grids** (`NLTE_CORRECTION_ELEMENTS` +
`pipeline/threed_corrections.py`), never hardcoded scalar approximations. **No UV NLTE grid
exists yet**, so every usable UV line is `nlte_status = GRID_OWED`:

- it carries a **cited expected magnitude** (informational, provenance-tagged — e.g. C I 1657
  ~+0.10, Amarsi 2020), and
- it is **LTE-flagged LOUDLY** downstream (RYA-426 gate 7) until a UV grid lands (RYA-165 path).

This module records the policy; it **does not apply** an approximate scalar into `[X/H]`.
Doing so would change solar/target science un-ratified and risk double-counting the existing
grid NLTE. The "complete implementation" (Amarsi-2020 UV NLTE grid + Teff/logg/[Fe/H]
interpolation) is the deferred sub-task (RYA-165).

## CORRECTIONS_3D — resolved (was an RYA-190 question)

RYA-190 flagged the legacy `CORRECTIONS_3D` dict (O 6300/6363 −0.07 scalars) as
"defined but never applied." **This is now superseded:** the project moved to the grid-based
`THREED_CORRECTION_ELEMENTS` + `pipeline/threed_corrections.py` (used by the Phase-C verdict).
The legacy `CORRECTIONS_3D` scalar dict is **dead (imported nowhere)** — recommend removal as
a tidy-up (left in place here to keep this PR scoped to selection/policy). UV 3D/NLTE follows
the same grid-based path.

## Multi-indicator C/O strategy

The Codex's strongest differentiator is the multi-indicator ratio with **explicit per-indicator
NLTE/3D corrections and per-arm zero-points** (never blind-averaged):

- **C** from: UV C I 1657 (synthesis, GRID_OWED) · optical C I 5052/5380 (preferred, less NLTE)
  · CO overtone IR (LTE-friendly). CH G-band 4300 is the optical molecular cross-check.
- **O** from: UV O I 1355 (semi-forbidden) · optical [O I] 6300 (Ni-blended, RYA-104) · OH IR (preferred).
- **N** from: far-red N I 7442/7468 (preferred) · NH 3360 · CN 3880 (blend-heavy). UV N I 1200 is a trap.

RYA-471 reports the UV-vs-optical agreement + zero-point per element; **disagreement is a
finding, never averaged away.**

## Open questions (carried)

1. **UV NLTE grid interpolation** — implement full Amarsi-2020 UV NLTE grids with
   Teff/logg/[Fe/H] interpolation (RYA-165), vs the current GRID_OWED + cited-magnitude policy.
2. **ISM Lyman-α** — for targets > 50 pc, use Lyman-α as a chromospheric-activity indicator
   rather than for abundances (separate thread). Nearby targets (Procyon 3.5 pc) are clear.
3. **VALD3 UV line-list extraction** — the α Cen A E140H (R~114k, ~1140–1710 Å) landmark set
   needs a dedicated VALD3 UV extraction before its pipeline run.
