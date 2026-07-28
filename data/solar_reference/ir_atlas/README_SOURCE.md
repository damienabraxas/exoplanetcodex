# IRTF Solar IR Atlas — Tier-2, DEFERRED (RYA-459)

**Status: documented + deferred — NOT staged.** This is a placeholder, not data. No
file in this directory is a solar spectrum; do not treat it as one.

## What it is

Tier-2 reference in RYA-459: a disk-integrated solar near-IR atlas, **0.94–2.5 µm**,
to feed the **CH / CO / OH** molecular diagnostics.

- Reference: **Rayner et al. 2009** (IRTF SpeX spectral library, ApJS 185, 289) per
  the RYA-459 ticket.
- Resolution: R ~ 2000 (SpeX), low for atomic lines but adequate for molecular bands.

## Why deferred (not blocking)

1. The 0.94–2.5 µm range primarily feeds CH/CO/OH, which ties to **Phase B**
   (CRIRES+ CO overtone, RYA-457) — and Phase B is **externally blocked** (STAGGER
   grid, collaborator-gated). No consumer is waiting on this atlas today.
2. The **CO Δv=2 overtone (2.29–2.35 µm)** is **already staged** at high resolution
   by RYA-390 — ACE-FTS (telluric-free, space) + NSO photatl + Wallace telluric,
   in `data/solar_reference/ir_atlases/`. That is the measured IR-CO anchor.
3. The Kitt Peak flux atlas (this ticket) already reaches **1300 nm**, covering the
   near-IR atomic lines (N I red, K I, P I) the program needs now.

So the IRTF atlas adds CH/CO/OH band coverage in 1.3–2.5 µm that only matters once
Phase B unblocks. Defer until then.

## Manual-action item for Ryan (when Phase B reopens)

- IRTF Spectral Library: http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/
  (NB: the SpeX library is STELLAR; a disk-integrated *solar* SpeX product may need a
  direct request to the IRTF / the Rayner group, or substitution by the
  Livingston & Wallace 1991 NSO IR photospheric atlas, parts of which RYA-390 already
  pulled — `https://nispdata.nso.edu/ftp/pub/atlas/photatl/`).
- Alternative measured solar IR (1.1–5.4 µm, disk-center): NSO photatl
  (Livingston & Wallace 1991, NSO TR 91-001) — same source family as RYA-390.

On stage: drop the atlas in the external store under
`Solar Calibration/IRTF Solar Atlas/`, extract CH/CO/OH segments here, and flip
`SOLAR_REFERENCE_SPECTRA['ir_atlas_irtf']['status']` off `deferred` with provenance,
coverage, units, and resolution filled from the actual product.
