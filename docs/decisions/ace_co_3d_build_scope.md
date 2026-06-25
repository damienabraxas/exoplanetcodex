# Scoping spike — full-3D-cube CO synthesis build feasibility

**Ticket:** RYA-444 (shallow scoping spike; desk research only, half-day box).
**Date:** 2026-06-25. **Branch:** ryandamienschmitt/rya-444-scoping-spike-3d-cube-co.
**Build-feasibility verdict: AMBER** — obtainable-on-request + Orion-runnable for the
LTE case, but gated on a human inquiry to a 3D-RT code author (no free public
download). Not a clone-and-go 6-week project; not blocked-forever either.

This scopes the ENGINEERING of the full-3D-cube solar CO synthesis whose SCIENCE was
greenlit by RYA-443 (CORROBORATED vs Amarsi et al. 2021: 3D CO = 8.48, same Li2015
gf, same ACE atlas, same disc-centre intensity geometry). No code obtained, no cube
downloaded, no synthesis run, no pipeline touched.

**Critical-path statement (holds regardless of verdict):** the 3D-cube CO arm is an
ENHANCEMENT to one IR arm of the solar anchor. It is NOT on the RYA-371 Phase A/C
critical path and NOT on the alpha Cen optical critical path. It stays a scoped
side-mission; this verdict does not gate any milestone work.

**What we actually need (narrows the build):** the CO 3D correction we are chasing is
3D **LTE** molecular line formation (Amarsi 2021 did CO in 3D LTE, not NLTE), on a
PUBLISHED STAGGER solar cube (we do NOT generate models), for the single 12CO (2-0)
overtone band, with the Li2015 CO line list we already use. That is the cheap corner
of the 3D problem — the expensive corners (RHD model generation; 3D NLTE) are out of
scope.

---

## Gate 1 — 3D RT code obtainability (the dominant risk)

Turbospectrum cannot do full-3D-cube synthesis (it does 1D and the <3D> average,
confirmed in RYA-442). None of the full-3D codes is a free public download; all are
request/collaboration. Per-candidate:

| code | what it is | access path | license / language / deps | obtainable |
|---|---|---|---|---|
| **Scate** (Hayek et al. 2011, A&A 535 A12) | 3D LTE RT post-processor; **what Amarsi 2021 used for CO**; reads STAGGER natively, takes Turbospectrum-format opacities | no public repo / no availability statement found; held by the Uppsala / Amarsi group | unstated; Fortran/MPI (HPC) | **REQUEST** (collaboration with Amarsi/Uppsala) — the exact-match path |
| **Linfor3D** (Steffen, AIP) | 3D LTE+NLTE synthesis post-processor; v6+ ported to free **GDL** (no IDL licence), runs parallel on HPC | public user manual at AIP; code on request to AIP; **reads CO5BOLD cubes, NOT STAGGER** | GDL (free) / IDL; CO5BOLD-native | **REQUEST** (AIP) — but needs STAGGER->CO5BOLD conversion or a CO5BOLD solar cube |
| **M3DIS / MULTI3D in DISPATCH** (Eitner et al. 2024, A&A 688 A52) | modern RHD + 3D LTE/NLTE synthesis; MULTI3D now ingests **Turbospectrum-format linelists** (matches our Li2015) | no code-availability statement in the paper; actively developed at MPIA (Bergemann group) | Fortran/C++ HPC framework (DISPATCH) | **REQUEST** (MPIA) — modern, TS-opacity match |
| **Balder** (Amarsi et al. 2018) / **Multi3D** (Leenaarts & Carlsson 2009) | 3D **NLTE** trace-element RT (Balder = modified Multi3D) | collaboration-only (Amarsi / Oslo) | Fortran/MPI HPC | **REQUEST**, and **overkill** — NLTE not needed for LTE CO |

**Gate-1 finding:** there is NO civilian clone-and-go 3D RT synthesis code. Every
candidate is obtainable only by a direct request / light collaboration. This is the
single biggest gating item and it requires a HUMAN INQUIRY (see flag below). It is
not RED (the codes exist and are routinely shared on request to research users), but
it is not GREEN (no `git clone` path). Best-aligned single contact: **Amarsi**
(Scate) — he performed the exact CO 3D LTE analysis we are reproducing, on the same
gf list and atlas geometry; second: **MPIA/Bergemann** (M3DIS, TS-opacity match).
Sources: Hayek et al. 2011 A&A 535 A12; Amarsi et al. 2021 A&A 656 A113 (Scate use);
Steffen Linfor3D manual (AIP); Eitner et al. 2024 A&A 688 A52 (M3DIS); Amarsi et al.
2018 (Balder); Leenaarts & Carlsson 2009 (Multi3D).

## Gate 2 — STAGGER solar 3D cube availability (the snapshot cube, not the <3D> average)

**AVAILABLE (public).** We currently hold only the 8 KB <3D> averaged solar member
(RYA-442). The full snapshot cube is publicly distributed:

- **Source:** Rodriguez Diaz et al. 2024, A&A 688, A480 (arXiv:2405.07872), "An
  extended and refined grid of 3D STAGGER model atmospheres — processed snapshots for
  stellar spectroscopy." All snapshots, mesh files, and analysis scripts are
  accessible from the STAGGER-grid website. (Magic et al. 2013 is the original grid;
  Chiavassa et al. 2018, A&A 611 A11, is the STAGGER synthetic-spectra precedent.)
- **Size:** ~316 MB per snapshot (full) or ~35 MB (reduced), each carrying T, rho,
  internal energy, momentum on the mesh; mesh + EOS in separate files. Amarsi used
  ~52 snapshots -> roughly **1.8 GB (reduced) to ~16 GB (full)** for the solar set.
  Well within Sirius/Orion storage; not a blocker.
- **Format / conversion:** STAGGER native binary (+ mesh + EOS). Scate and M3DIS
  read STAGGER directly; Linfor3D would need a STAGGER->CO5BOLD conversion (or use a
  CO5BOLD solar cube instead). Conversion is a known, bounded step, not a research
  problem.

**Gate-2 finding:** the cube is not a blocker — public, modest size, native format
for the two best-aligned codes (Scate, M3DIS). Do NOT download it during this spike.

## Gate 3 — Orion compute sizing (order-of-magnitude)

The headline cluster-class numbers in the literature are for the parts we DO NOT do:

- 3D RHD MODEL GENERATION: ~10^4 CPU-h per model up to ~10^11 CPU-h for a full grid
  (Eitner et al. 2024). We use a PUBLISHED cube -> this cost is zero for us.
- 3D **NLTE** spectrum synthesis: ~10^5 x the 1D NLTE cost; ~10^6-10^7 total for a
  full species spectrum (3D NLTE RT theory reviews). Cluster-only. **We do not need
  NLTE for the CO molecular correction.**

What we DO: 3D **LTE** single-band post-processing synthesis of 12CO (2-0) over one
solar cube, plus a small A(C) grid (~5 points). LTE has no NLTE iteration — it is one
formal solution per (snapshot x column x ray x wavelength). Order-of-magnitude:
~52 snapshots x ~80x80 columns x a few rays x a few hundred wavelengths x ~5 A(C)
points, each a sub-second 1D-like formal solve, fully parallel ->

- **CPU:** order 10^1-10^3 core-hours (a single band, embarrassingly parallel).
- **RAM:** a snapshot is 35-316 MB; synthesis holds a snapshot + opacity/EOS tables
  -> a few to a few tens of GB peak (process one snapshot at a time).
- **Translation:** a single fat node — order **16-64 cores, 64-128 GB RAM,
  wall-time hours-to-a-few-days** for the band + A(C) grid.

**Gate-3 finding:** Orion-runnable for the LTE case on one fat node; **NOT
cluster-only** (that verdict applies to 3D NLTE, which we do not need). Sirius (16 GB
ProBook) cannot host it (RAM + cores), so this is an Orion-class single-node job.
Numbers are order-of-magnitude; confirm against the chosen code's own benchmark
before committing wall-time.

---

## Verdict: AMBER — feasible on a code request, Orion-runnable for LTE

| gate | finding | colour |
|---|---|---|
| 1 code | no free download; obtainable on request (Scate/Amarsi exact-match, or M3DIS/MPIA, or Linfor3D/AIP) | amber (dominant risk) |
| 2 cube | public (Rodriguez Diaz 2024), ~2-16 GB, native for Scate/M3DIS | green |
| 3 compute | 3D LTE single-band: ~10^1-10^3 core-h, ~64-128 GB, one fat node -> Orion-runnable | green (amber if forced to NLTE) |

**AMBER, not GREEN:** the only thing standing between us and a schedulable build is
obtaining a 3D RT code, and none come as a public download — it takes a direct
inquiry / light collaboration. **AMBER, not RED:** the cube is public, the LTE
compute is Orion-class, and the codes are routinely shared with research users on
request, so the build is not blocked-forever.

**Rough effort once a code is in hand:** ~3-6 weeks of integration — cube ingestion
+ EOS/opacity wiring, 12CO (2-0) band synthesis over the cube, the A(C) chi2 fit
reusing the 441 harness conventions, and validation that our 3D CO lands on Amarsi's
8.48 (the built-in cross-check). Single-node Orion compute. This is an engineering
project, not a research project — PROVIDED the code request succeeds.

## Needs a direct human inquiry (flagged)

The gating action is a HUMAN EMAIL to obtain a 3D RT synthesis code. Recommended
first contact: **A. M. Amarsi** (Uppsala) for **Scate** — he did the exact CO 3D LTE
analysis we reproduce, on the same Li2015 gf and ACE-class geometry, so his pipeline
is the lowest-friction match and a natural light collaboration. Fallbacks: MPIA /
Bergemann group for **M3DIS** (modern, Turbospectrum-format opacities), or AIP /
Steffen for **Linfor3D** (free GDL, but CO5BOLD-native -> cube conversion). Until one
of these is obtained, the build cannot be scheduled — but nothing on the 371 or
alpha Cen critical paths waits on it.

## Sources

- Hayek, Asplund, Collet et al. 2011, A&A 535, A12 — Scate (3D LTE RT).
- Amarsi, Grevesse, Asplund & Collet 2021, A&A 656, A113 (arXiv:2109.04752) — CO 3D
  LTE = 8.48 with Scate + STAGGER + Li2015; the validation target.
- Eitner, Bergemann et al. 2024, A&A 688, A52 (arXiv:2405.06338) — M3DIS/DISPATCH;
  compute scales.
- Steffen, Linfor3D user manual (AIP) — GDL port, CO5BOLD-native.
- Amarsi et al. 2018 — Balder; Leenaarts & Carlsson 2009 — Multi3D (3D NLTE).
- Rodriguez Diaz et al. 2024, A&A 688, A480 (arXiv:2405.07872) — public STAGGER
  snapshots for spectroscopy (cube availability + size + format).
- Magic et al. 2013, A&A 557, A26 — original STAGGER grid; Chiavassa et al. 2018,
  A&A 611, A11 — STAGGER synthetic spectra precedent.
- 3D NLTE vs 1D cost scaling — 3D NLTE RT reviews (e.g. arXiv:2511.04254).
