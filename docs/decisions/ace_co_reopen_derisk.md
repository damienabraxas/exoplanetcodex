# Decision record — de-risk the ACE-FTS solar-CO REOPEN before any full-3D build

**Ticket:** RYA-443 (de-risk gate before authorizing a full-3D-cube synthesis build).
**Date:** 2026-06-24. **Branch:** ryandamienschmitt/rya-443-ace-co-reopen-derisk.
**Verdict: CORROBORATED — authorize scoping the full-3D-cube build (separate ticket).**

Literature de-risk only; no new synthesis. All three checks rest on PRIMARY-SOURCE
values (cited), with the gf cross-check actually performed and the CH/CO sign
genuinely resolved (not assumed).

---

## Why this gate exists

RYA-442 reopened ACE on an implied 1D->3D correction of -0.13..-0.26 dex. But that
band was computed as `(our 1D 8.646) - (published 3D ABSOLUTES 8.39/8.47/8.52)` =
our total offset from the community solar C = `3D effect + any error in our own 1D
pipeline`. Attributing all 0.186 dex to 3D is the optimistic read; whatever is
pipeline, 3D will not fix. So REOPEN rests on one untested assumption: that our
8.646 is a FAITHFUL 1D result. This record tests that with a CO-specific,
same-gf-scale 1D anchor and the published CO-specific differential.

**Reproducibility:** the RYA-441 quarantined harness re-runs bit-identically to
A(C)_mu1 = 8.646 (ATLAS9, mu=1 intensity). No production path touched.

---

## Step 1 — faithful-1D anchor + gf cross-check  → PASS

**Primary-source correction first:** the ticket proposed Grevesse, Lambert & Sauval
1991 (A&A 242, 488) as a 1D *CO* anchor (8.60 +/- 0.05). At primary-source
precision that paper is the **CH** vibration-rotation analysis (104 CH vib-rot
lines, A(C) = 8.60 +/- 0.05), **not** CO-specific. It is therefore the wrong anchor
for a CO faithfulness test (CH and CO behave oppositely in 3D — see Step 3). It
remains a useful sanity point (a 1D molecular A(C) ~8.6) but is not load-bearing.

**The CO-specific faithful-1D anchor used instead** — Amarsi, Grevesse, Asplund &
Collet 2021 (A&A 656, A113), 3D LTE molecular CNO. Their **1D** solar A(C) from the
12C16O `Delta-nu = 2` overtone (2295-2591 nm = exactly our 2.3 um band), on the
**Li et al. 2015** CO line list:

| source | model | A(C) from CO overtone (1D) |
|---|---|---|
| Amarsi 2021 | 1D MARCS | 8.608 |
| Amarsi 2021 | 1D ATMO  | 8.621 |
| **ours (RYA-441)** | **1D ATLAS9, mu=1** | **8.646** |

Delta(ours - Amarsi 1D) = **+0.038 (vs MARCS), +0.025 (vs ATMO)** -> |Delta| ~0.03,
well within the 0.07 faithfulness bar. The residual is consistent with an
ATLAS9-vs-MARCS 1D-model difference plus our exact line subset; it is NOT a gross
pipeline error.

**gf cross-check (the convenient-agreement guard):** our CO list is
`CO_IR_Li2015.dat` (header literally `ExoMol Li2015`, 12C16O species 0608.012016).
Amarsi 2021 adopt **the same** Li et al. 2015 transition probabilities for 12C16O.
The near-agreement of the two 1D absolutes is therefore on an IDENTICAL gf scale by
provenance — it is not a coincidence masking a compensating gf offset. **gf
cross-check PASSES.**

-> **Pipeline faithful.** Our 8.646 is a genuine 1D CO result, ~0.03 dex above a
modern same-gf 1D CO analysis; it is not inflated by a pipeline defect.

## Step 2 — CO-specific 1D->3D DIFFERENTIAL (published, not absolute subtraction) → toward reference

From Amarsi 2021 (A&A 656, A113), 12C16O `Delta-nu = 2` overtone:

- A(C)_3D = **8.479**, A(C)_1D-MARCS = 8.608  ->  **CO 1D->3D differential = -0.129 dex**
- (vs 1D ATMO 8.621: -0.142 dex)

Corroborated by Asplund, Amarsi & Grevesse 2021 (A&A 653, A141), the 2020-vision
review, CO `Delta-nu = 1`: A(C)_3D = 8.487, A(C)_1D-MARCS = 8.606 -> -0.119 dex.
(Two independent reductions agree: CO 1D->3D ~ -0.12 .. -0.14.)

Applying the published CO-specific differential to OUR faithful 1D value:

> 8.646 + (-0.129) = **8.517**  ->  **+0.057 dex vs the adopted reference 8.46**

Within +/-0.10. **REOPEN holds on a proper differential, not just an absolute
coincidence.** (Note this is much smaller than RYA-442's -0.13..-0.26 absolute-
subtraction band; it lands within band only BECAUSE our 1D is faithful — exactly
the assumption this gate confirmed.)

## Step 3 — CH/CO sign resolution (the single biggest risk) → CO is canonical-negative

Modern updated-STAGGER work finds **CH** A-X 3D LTE runs **HIGHER** than 1D
(positive correction), opposite the canonical negative molecular correction
(Asplund et al. 2021; Popa et al. 2025, arXiv:2511.14289). The threat: if **CO**
shared that positive shift, REOPEN would collapse. Resolved from the primary CO
data:

| molecule | 3D-1D-MARCS correction | sign |
|---|---|---|
| 12C16O (overtone) | -0.13 (Amarsi 2021) / -0.12 (Asplund 2021) | **NEGATIVE** (toward ref) |
| CH (A-X)          | +0.045 (Amarsi 2021)                       | POSITIVE (away) |

**CO and CH have OPPOSITE-sign 3D-1D corrections.** The modern-STAGGER CH positive
shift does NOT apply to CO; CO vibration-rotation behaves as the canonical negative
molecular correction (toward the reference). The biggest threat to REOPEN is
resolved in REOPEN's favour — explicitly, from the CO numbers, not assumed.

## Bonus — reconciles the RYA-442 <3D> probe sign

RYA-442's <3D> mean-model probe gave a POSITIVE correction (+0.147) and was flagged
"sign-wrong" against a tripwire that assumed <3D> must be negative. The literature
decomposition shows that assumption was wrong FOR CO: Asplund 2021 gives, for CO,
A(C)_<3D> = 8.653 vs A(C)_1D-MARCS = 8.606, i.e. **<3D> - 1D = +0.047 (POSITIVE)**.
The full negative correction is almost entirely the `3D - <3D>` inhomogeneity term
(8.487 - 8.653 = -0.166). So the 442 probe sign was NOT a bug — the <3D> mean-T term
for CO is genuinely positive; only full 3D (with granulation inhomogeneity) is
negative. This strengthens the case that a FULL-3D cube (not a mean model) is the
right and necessary tool, and confirms 442's "partial / lower-bound" framing.

---

## Verdict: CORROBORATED

1. Pipeline faithful: our 1D CO 8.646 matches a modern same-gf (Li2015) 1D CO
   analysis (8.608/8.621) to ~0.03 dex; gf cross-check passes by provenance.
2. CO-specific published 1D->3D differential (-0.13) applied to 8.646 lands at 8.52,
   within 0.10 of reference.
3. CO 3D-1D sign is NEGATIVE (toward reference), opposite CH's modern positive shift.

-> **Authorize scoping the full-3D-cube synthesis build (Linfor3D-class) as a
separate ticket.** Expectation set by the literature: a correct full-3D disk-center
12CO (2-0) treatment should move our 8.646 by about -0.12..-0.14 dex to ~8.50-8.52,
i.e. to within ~0.05-0.06 of 8.46 — a viable solar CO source. Not a turnkey
guarantee (our +0.03 1D residual and model-generation effects remain), but the 1D
slab is sound and the sign is settled, so the 3D build is not being raised on a
crack.

## Anti-motivated-reasoning audit (this search corroborates what we want)

- The first automated pull of Amarsi 2021 returned CO 1D = 8.50 / differential
  -0.02 (which would have FLIPPED the verdict to PIPELINE-OFFSET / collapse). It was
  cross-checked against the ar5iv full text and the Asplund 2021 review and found
  WRONG; the two independent primary sources agree on 1D CO ~8.61 and differential
  ~-0.13. Held the line — did not take the first convenient (or inconvenient) number.
- gf cross-check actually performed (line-list provenance identity, both Li2015).
- CH/CO sign genuinely resolved from the CO numbers, not assumed negative.

## Sources

- Amarsi, Grevesse, Asplund & Collet 2021, A&A 656, A113 (arXiv:2109.04752) — 3D LTE
  molecular CNO; 12C16O overtone 2295-2591 nm, Li et al. 2015 gf; CO 1D MARCS 8.608 /
  3D 8.479 (-0.129); CH +0.045.
- Asplund, Amarsi & Grevesse 2021, A&A 653, A141 (arXiv:2105.01661) — 2020 vision;
  CO 1D MARCS 8.606 / <3D> 8.653 / 3D 8.487; adopted A(C) = 8.46 +/- 0.04.
- Popa, Hoppe, Bergemann et al. 2025, MNRAS (arXiv:2511.14289) — 3D NLTE CH, modern
  STAGGER (CH positive shift).
- Grevesse, Lambert, Sauval, van Dishoeck, Farmer & Norton 1991, A&A 242, 488 —
  CH (not CO) vib-rot, A(C) = 8.60 +/- 0.05, 104 CH lines.
- Li, Gordon, Rothman et al. 2015, ApJS 216, 15 — the CO line list / gf used by both
  our pipeline and Amarsi 2021.
- Reference A(C) = 8.46: constants.py SOLAR_ASPLUND2021['C'].
