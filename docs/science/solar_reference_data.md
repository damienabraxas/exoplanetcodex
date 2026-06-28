# Solar reference data — library, provenance, and coverage (RYA-459 / RYA-162)

The multi-wavelength solar reference library. HARPS-VIS (380–690 nm) cannot reach the
lines that matter for C/N/O/P/S; this library stages the atlases that can, each wired
with an explicit **measured vs cited** provenance flag.

Registry: `config.constants.SOLAR_REFERENCE_SPECTRA`.
Audit / coverage gate: `python -m pipeline.audit_solar_reference --verify`.

## Cardinal honesty rule (RYA-455 discipline)

Hubble cannot observe the Sun directly → there is **no direct solar UV spectrum**. The
UV composite is `provenance=cited-composite` and is **never** presented as a
measurement. The audit gate fails loud if any UV/composite source is tagged
`measured`. Any UV-derived abundance inherits the cited flag downstream.

## Sources

| source | provenance | coverage | resolution | flux units | citation |
|--------|-----------|----------|-----------|-----------|----------|
| `kpno_flux_atlas` | **measured** | 296–1300 nm | FTS high-res (Δλ ≈ 0.004–0.013 Å, λ/Δλ ~ 4e5) | normalized residual + irradiance µW/cm²/nm | Kurucz, Furenlid, Brault & Testerman 1984 (NSO Atlas No. 1) |
| `uv_composite` | **cited-composite** | 119.5–2695.7 nm | low (Δλ ≈ 20 Å, R ~ 150–300) | FLAM erg/s/cm²/Å | Colina, Bohlin & Castelli 1996 / Bohlin+2001 (CALSPEC) |
| `ir_atlases` (RYA-390) | measured | 2289–2350 nm (CO Δv=2) | FTS (ACE 0.02 cm⁻¹) | normalized intensity | ACE-FTS Hase+2010 / NSO photatl / Wallace telluric |
| `ir_atlas_irtf` | (measured) | 940–2500 nm | R ~ 2000 | TBD | Rayner+2009 IRTF — **DEFERRED** |

## Coverage matrix — diagnostic → reference

Kitt Peak (measured, resolved lines) is the working anchor ≥296 nm; CALSPEC (cited,
low-res) provides only the deep-UV (<296 nm) and the absolute-flux scale and cannot
resolve lines.

| diagnostic | line (Å) | element | Kitt Peak (MEASURED) | CALSPEC UV (CITED) | consumer |
|-----------|---------:|---------|----------------------|--------------------|----------|
| NH 3360 | 3360.0 | N | YES — n=6265, Δλ 0.004 Å | covered (cited) | RYA-369 N (UV band head) |
| CN violet 3883 | 3883.0 | N | YES — n=5565, Δλ 0.005 Å | covered (cited) | RYA-369 N cross-check |
| N I 7442/7468 | 7455 | N | YES — n=3369, Δλ 0.010 Å | covered (cited) | RYA-369 N I red |
| N I 8216/8223 | 8219 | N | YES — n=1874, Δλ 0.010 Å | covered (cited) | RYA-369 N I red |
| N I 8680–8718 | 8699 | N | YES — n=4275, Δλ 0.011 Å | covered (cited) | RYA-369 N I red multiplet |
| [O I] 6300 | 6300.3 | O | YES — n=1388, Δλ 0.007 Å | covered (cited) | O forbidden cross-check (RYA-455) |
| O I 777 triplet | 7773 | O | YES — n=1164, Δλ 0.009 Å | covered (cited) | O I 777 PRIMARY O (RYA-455) |
| P I 10581/10596 | 10589 | P | YES — n=2085, Δλ 0.013 Å | covered (cited) | P near-IR (alt to FUV/HST) |
| K I 7665/7699 | 7682 | K | YES — n=4767, Δλ 0.009 Å | covered (cited) | K DATA-GAP resonance doublet |
| Co I 3845 | 3845.5 | Co | YES — n=1134, Δλ 0.004 Å | covered (cited) | Co DATA-GAP probe |
| Sc II 4246 | 4246.8 | Sc | YES — n=503, Δλ 0.010 Å | covered (cited) | Sc DATA-GAP probe |

**Headline: 11/11 diagnostics now have a MEASURED (Kitt Peak) reference**, including
all **5 solar-N channels** (NH 3360, CN violet, N I 7442/7468, 8216/8223, 8680–8718)
— the RYA-369 solar-N unblock. The Kitt Peak near-IR reach (to 1300 nm) also lands
measured references for the K/Co/Sc/P DATA-GAP elements (RYA-371 Phase C), and a
**P I near-IR multiplet** (10581/10596 Å) as an alternative to the FUV P lines that
otherwise need HST/STIS (RYA-119).

## Which reference for which job

- **Resolved line work (EW, profiles, synthesis-fit) at ≥296 nm** → Kitt Peak flux
  atlas (measured, high-res). This is the differential anchor for N / CNO / K-Co-Sc-P.
- **Deep-UV (<296 nm) and absolute flux calibration** → CALSPEC composite
  (cited-composite). UV-derived numbers carry the cited flag.
- **K-band CO overtone (2.3 µm), ¹³C** → `ir_atlases` (RYA-390, measured), gated by the
  parked Phase B (RYA-457).
- **CH/CO/OH 1.3–2.5 µm** → IRTF (Tier-2, deferred until Phase B unblocks; see
  `data/solar_reference/ir_atlas/README_SOURCE.md`).

## Still missing / gated

- **IRTF 0.94–2.5 µm** — documented + deferred (no consumer until Phase B / RYA-457).
- **Deep-UV resolved lines (<296 nm)** — only the low-res cited composite exists; truly
  resolved solar UV is physically unavailable (no direct solar UV spectrometer).
