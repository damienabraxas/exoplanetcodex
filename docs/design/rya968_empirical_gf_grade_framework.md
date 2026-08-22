# RYA-968 — Empirical gf-grade framework (v2, literature-grounded)

**DESIGN SPEC — for sign-off. No implementation.** Child of RYA-958.
**Scope: the ABSOLUTE layer.** RYA-969 is the differential half.
Firewall: RYA-161. Extends: RYA-855/850. Renders through: RYA-707/224/851/809. Doc sync: RYA-179.

> **v2 supersedes v1.** v1 proposed measuring each ungraded line's excess scatter over a graded
> floor and carrying it as an inflated per-line σ. The literature grounding (2026-08-22) replaces
> that with the field's actual method — **selection**, not inflation. What survives from v1 is the
> per-line architecture, the firewall structure, and the measurements; what changes is the
> mechanism. §3 records what v1 got wrong.

---

## 1. This is published practice, not our invention

The GBS third version (Soubiran+2024) and OCCASO (Carbajo-Hijarrubia+2024) use the
laboratory-graded lines as a **reference distribution** and **admit** ungraded lines whose
abundance falls within tolerance of that anchor. Rejected lines are documented, not carried with
a large σ. That is the method; we are replicating a benchmark methodology to our provenance
standard, not inventing one.

Two flag systems to mirror rather than invent:

* **Gaia-ESO / Heiter+2021** — two independent axes, `gf_flag` and `synflag`, each `Y/U/N`.
  Reliability and blending are *separate* questions and a single tier cannot carry both.
* **Elgueta+2026** (RYA-787; already vendored as `elgueta2026`) — a four-stage physical decision
  tree, **Depth → Saturation → Purity → goodness-of-fit**, each `Y/N/U`, anchored to laboratory
  data and *explicitly rejecting APOGEE-style astrophysical fine-tuning*.

🔴 **Elgueta is the RYA-161 firewall stated by an external group.** Its stage structure is the
tier logic; adopting it means our firewall is the field's, not a house rule.

**Why the problem is unavoidable:** roughly half of optical lines lack a reliable laboratory gf
(Jofré, Heiter & Soubiran 2019, ARA&A 57). There is no version of this project that avoids
deciding what to do with ungraded lines.

### Reference check
Every reference above was verified against arXiv/Crossref before being written down, and **two
attributions in the discussion were off**:

| as discussed | actually |
|---|---|
| "Heiter 2019, arXiv 1811.08041" | **Jofré**, Heiter & Soubiran 2019, ARA&A 57 — Heiter is second author (`10.1146/annurev-astro-091918-104509`) |
| "NGC6352, arXiv 0810.4832" | **Feltzing**, Primas & Johnson 2009, A&A 493, 913 (`10.1051/0004-6361:200810137`) |
| "GBS third version (Jofre, A&A aa55211-25)" | ⚠️ **that manuscript number does not resolve.** The third-version papers that do are **Soubiran+2024** (fundamental Teff/log g, `10.1051/0004-6361/202347136`) and Vitali+2026 (n-capture, `10.1051/0004-6361/202661004`). If a Jofré-led third-version *metallicity* paper is in press, **we need its identifier before citing it.** |

`elgueta2026` was already vendored with the correct DOI. `soubiran2024`,
`carbajohijarrubia2024`, `jofre2019` and `feltzing2009` added to `data/refs/bibliography.csv`
(RYA-854).

---

## 2. Scope — the absolute layer, and only that

| layer | ticket | gf | precision |
|---|---|---|---|
| **ABSOLUTE** | **RYA-968 (this)** | does **not** cancel; per-line gf error is real and irreducible | this is where the honest ~0.1–0.5 dex lives |
| DIFFERENTIAL | RYA-969 | **cancels** star-to-star | 0.01–0.02 dex (Meléndez / Bedell / Nissen) |

**The 0.17 dies in both layers, by different mechanisms**: cancelled in the differential,
measured-and-replaced here. The absolute number being large is not a failure of this framework —
it is the correct statement of where the Sun's absolute A(X) actually sits.

**The astrophysical-gf firewall follows from the split** (Feltzing+2009): solar-fitted gf is
*legitimate* for differential work, where it cancels, and *forbidden* for absolute work, where it
is circular. RYA-161 is the house name for a distinction the literature already draws.

---

## 3. 🔴 What v1 got wrong, and the two measurements that show it

v1 proposed: measure ungraded scatter, quadrature-subtract the graded floor, carry the excess
(~0.49 dex for Kurucz) as a per-line σ. **The field does not do this, and it is the wrong shape.**
It keeps every line and inflates the error; the field *selects* and keeps the error small on the
survivors. v1's ~0.49 dex is best understood as **the cost of not selecting** — a property of the
unselected Kurucz population that nobody carries.

But adopting the published method wholesale has two failure modes, and **both are measured on our
own data, not hypothesised**:

### 3.1 The published REW cut does not bite on our data

GBS admits lines in `−6.7 < REW < −4.5`. Our VIS Fe I pool spans **−6.09 to −4.90** — entirely
inside it. **249 of 249 lines survive the cut (100 %).**

⚠️ **So the published confound control removes nothing here.** And the confound is still present
*inside* the window: LAB median REW **−5.03** against Kurucz **−5.28**, LAB median EW **55.5 mÅ**
against Kurucz **31.0**. Importing the GBS window and declaring the confound handled would be
exactly the lazy move the ticket warns about. **We need a window derived from our own
distribution, or matched admission, not a borrowed constant.**

### 3.2 🔴 Admission manufactures precision the anchor cannot support

Applying reference-matched admission against our 7-line anchor (mean 7.498, sd 0.157):

| tolerance | admitted | rejected | admitted sd | admitted **sem** |
|---|---|---|---|---|
| ±0.05 dex (GBS) | 23 of 249 | 226 | **0.024** | **0.0049** |
| ±1σ of anchor (0.157) | 77 of 249 | 172 | 0.084 | 0.0096 |

The admitted pool's scatter (0.024) is **smaller than the anchor's own** (0.157) — *by
construction*, because the lines were selected for agreeing with the anchor mean. Quoting
sem = 0.005 dex as the absolute uncertainty would be nonsense:

> **The absolute zero-point is set entirely by the laboratory lines.** With n = 7 and sd = 0.157,
> the anchor supports σ_mean ≈ **0.059 dex**, and *no number of admitted ungraded lines can beat
> that*, because they carry no independent information about the scale.

**This is the central design constraint of v2**, and it is where the method would quietly become
circular if built naively. Admission buys *statistical* precision around a zero-point it cannot
improve.

⚠️ Also note **only 1 of the 7 anchor lines survives its own ±0.05 dex test** — six of seven lab
lines disagree with the lab mean by more than the GBS tolerance. GBS applies that number to a
far better-behaved pool. **Adopting ±0.05 unexamined would reject our own anchor.**

---

## 4. Method (v2)

### 4.1 Admission, in stages — Elgueta's tree, our data

Per line, in order, each `Y/N/U`, mirroring Elgueta+2026 and carrying Gaia-ESO's two axes:

1. **Depth** — measurable against the noise.
2. **Saturation** — on the linear curve of growth; a *derived* REW window (§3.1), not a borrowed one.
3. **Purity** — unblended. This is Gaia-ESO's `synflag` axis.
4. **Goodness of fit** — the profile actually describes the feature; RYA-959's width ceilings live here.
5. **gf reliability** — Gaia-ESO's `gf_flag` axis: laboratory / NIST-graded / ungraded.
6. **Anchor consistency** — does the line behave like the laboratory-anchored distribution?

Stages 1–4 are **physical and self-contained**: they never look at any abundance. Stage 6 is the
only one that compares abundances, and §5 constrains exactly what it may compare against.

### 4.2 Per-line σ, precedence unchanged from v1

**cited → self-reported (cross-engine/band spread) → inferred → fallback.** Self-reported needs
no anchor and no confound model, which is why it outranks anything inferred. The fallback is our
own pooled measurement, not Kurucz's 0.17.

### 4.3 The zero-point cap — new in v2, and non-negotiable

Any absolute product reports **two** uncertainties:

* **statistical** — from the admitted pool, shrinks with n;
* **zero-point** — from the laboratory anchor alone, `σ_anchor/√n_anchor` ⊕ the anchor's own
  cited σ. **It does not shrink when ungraded lines are admitted.**

A budget that reports only the first has manufactured precision (§3.2). Today that floor is
**≈0.059 dex** and it is set by seven lines.

---

## 5. The firewall — reconciled with admission

Admission-on-agreement looks like the accuracy-grading RYA-161 forbids. It is not, and the
distinction must be written into the code, not the intent:

| | |
|---|---|
| **FORBIDDEN** | agreement with an **external expected abundance** — a literature A(X), a solar-fitted astrophysical gf, an APOGEE-style calibration. Circular: you recover what you assumed. |
| **PERMITTED** | agreement with the **laboratory-anchored distribution measured in our own spectrum**. The lab lines carry an *independent physical calibration* — an apparatus measured their gf. "Does this line behave like lines whose gf we know?" is a **gf-quality test**, not an answer-tuning test. |

Elgueta+2026 draws exactly this line by anchoring to laboratory data while rejecting
astrophysical fine-tuning.

**Structural controls (v1's, still binding):**

* **F1 — no external reference is an input.** Grading functions receive per-line measurements and
  the *anchor built from our own lab-graded lines*; never a literature value, gold value, or
  target. Tested by a grep-level import check.
* **F2 — stages 1–5 are invariant under a constant offset.** Add δ to every abundance and the
  physical stages and every σ come out bit-identical. ⚠️ **Stage 6 is deliberately not
  offset-invariant** — that is what it is for — so F2 applies to 1–5 and **stage 6 must be
  separately auditable**: its admitted/rejected list is a required output, per line, with the
  distance from the anchor recorded.
* **F3 — thresholds declared before the data is seen**, in config with a ticket reference.
  §3.1/§3.2 show why a borrowed threshold is not exempt from this.
* **F4 — the Cr canary (+0.402) is a blocking test.** If grading shrinks it, the build stops.

**The mean offset by tier** (LAB 7.498 → VALD3 7.848) is reported as a **finding, never an
input** — it may be the Kurucz zero-point RYA-819/831 chased.

---

## 6. Failed lines are DOCUMENTED, not quarantined — and the renderer already exists

A rejected line is a **first-class deliverable**: retained in the record, **excluded from the
value**, and given an appendix entry that *shows* why it failed, with a per-line diagnostic plot.
Same mechanism the graded lines already use, one tier down: measured, documented, plotted, not
added to the number.

🔴 **968 does not build any of that.** It already exists:

| ticket | provides |
|---|---|
| **RYA-707** (done) | SPP Appendix A + the unmeasured-line proof generator |
| **RYA-224** (standing standard) | per-line diagnostic plots — generation, naming, Linear attachment |
| **RYA-851** (in progress) | the live Solar page + Fe appendix where this renders |
| **RYA-809** (done) | the reason-code taxonomy, already applied to a real Fe set |

**968's output is exactly two fields per line: a TIER and a REASON-CODE**, drawn from the
existing vocabulary (`BAD_GF`, `ATOMIC_BLEND`, `MOLECULAR_BLEND`, `TELLURIC_ADJACENT`,
`SATURATION_COG`, `CONTINUUM_LIMITED`, `NON_MINIMUM`, `ABUNDANCE_OUTLIER`, …) — **no parallel
vocabulary, no appendix schema, no plot generator.**

Deletion is reserved for nothing. Even a non-feature — zero-flux telluric core, ghost — is
documented *as such*.

---

## 7. Gates

| status | gate |
|---|---|
| 🔴 **RED** | **the anchor is 7 lines.** It no longer gates a floor estimate, but it now caps the absolute zero-point at ≈0.059 dex (§4.3) and it is too thin to define an admission distribution. Growing it toward RYA-945's 199 VIS lab lines is the highest-value unblock in this whole area. |
| 🔴 **RED** | **the REW window must be derived from our data** (§3.1). The GBS constant admits 100 % of our pool and controls nothing. |
| open | admission tolerance chosen from our anchor's behaviour, not borrowed — ±0.05 rejects 6 of our own 7 anchor lines (§3.2). |
| open | F2 passes on synthetic data for stages 1–5; stage 6 emits its per-line audit list. |
| open | the Cr canary is unmoved. |

---

## 8. Beyond Fe

Fe I is the pilot only because it *has* a laboratory anchor. Stages 1–4 are element-agnostic and
apply anywhere. Stage 6 needs an anchor and therefore **does not exist for un-anchored elements**
(Al and the rest of the backbone) — those get the self-reported and fallback routes only.

Whether an anchor-derived criterion transfers between elements is an open question and the honest
default is **that it does not**: a floor is a property of a pipeline *and* a line list.
**Recommend Fe-only in v1 of the implementation.**

---

## 9. Open for sign-off

1. **§3.2 — do you accept the zero-point cap?** The absolute A(Fe) cannot be quoted tighter than
   the lab anchor supports (≈0.059 dex on 7 lines), however many ungraded lines are admitted.
   This is the v2 equivalent of v1's "~3× wider" question, and it is the one that matters.
2. **§7 — growing the 7-line anchor.** Now gates two things, not one.
3. **§3.1 — who derives the REW window**, given the published one does not bite?
4. **Fe-only in v1?** Recommended.
5. **The unresolved GBS third-version metallicity reference** (§1) — do you have its identifier?

---

## 10. Evidence index

All §3 numbers reproduce from committed data:
`data/results/band_products/FeI_3780_6910_harps_solar_harps_molecfit_corrected_PROFILEFIT_1D-LTE_lines.csv`
(RYA-959; 249 in-aggregate lines with abundance + REW) joined ±0.02 Å to
`data/linelists/canonical_gf.csv` `gf_tier` (RYA-945). Anchor statistics from the 7 `LAB`-tier
lines in that pool. References verified against arXiv and Crossref, 2026-08-22.
