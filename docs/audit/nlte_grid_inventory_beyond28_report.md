# NLTE grid inventory beyond the 28 — RYA-758

**Audit date:** 2026-08-09 · **Scope:** every naturally occurring element Be (Z=4)
to U (Z=92) that is *not* in the canonical set · **Artifacts:**
`data/audit/nlte_grid_inventory_beyond28.csv` (54 rows),
`data/audit/nlte_grid_inventory_beyond28_arxiv_scan.json` (raw post-Lodders sweep) ·
**Generators:** `scripts/build_nlte_inventory_beyond28_rya758.py`,
`scripts/scan_arxiv_post_lodders_rya758.py`,
`scripts/spotcheck_nlte_grid_loader_rya758.py`

Discovery only. No canonical-set change, no element wiring, no edit to
`pipeline/constants.py`, `pipeline/abundances_derive.py`, `data/reference/solar/`,
or `CODEX_STATE_REGISTER.md`. The `SCIENCE_STANDARDS.md` amendment is drafted as a
proposal for ratification, not applied:
`docs/proposals/science_standards_grid_availability_amendment.md`.

---

## Executive summary

Two findings change what we thought we knew. **CRITICAL:** a full 3D NLTE
outside-28 element exists that the session's first pass missed — **silver**
(Caliskan, Amarsi, Jönsson, Grevesse & Sahoo 2026, A&A 711, A155), whose departure
grid is publicly downloadable from Zenodo in the *same PySME format we already
ingest for Cu*, making Ag the single cheapest expansion candidate in the entire
outside-28 universe. **HIGH:** the fluorine question resolves in the negative but
for the wrong reason — F is not grid-blocked, it is *diagnostic*-blocked, because no
F line is detectable in the quiet-Sun photosphere at all and the accepted solar
value comes from HF in a sunspot umbra. Beyond those, the expansion universe is
genuinely tiny: of 54 outside-28 elements, only **8** have any NLTE treatment at all
(2 in 3D, 6 in 1D), **21** have no photospheric diagnostic whatsoever, and the
remaining **25** sit on 1D/3D LTE with no model atom in existence. The canonical
inclusion criterion should therefore be sharpened around *grid tier* rather than
around science interest, because science interest is not the binding constraint —
molybdenum is the best-motivated element outside the 28 on the science and is
blocked purely by the absence of a Mo model atom, a state the most recent dedicated
Mo paper (November 2025) confirms explicitly.

---

## Findings by tier

Counts over the 54 rows of `data/audit/nlte_grid_inventory_beyond28.csv`.

| `grid_tier` | n | elements |
| --- | --: | --- |
| `3D_NLTE` | 2 | Be, **Ag** |
| `1D_NLTE` | 6 | Rb, Pd, Pr, Nd, Pb, Th |
| `LTE_only` | 25 | B, Ga, Ge, Nb, Mo, Ru, Rh, Cd, Sn, La, Ce, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, W, Os, Ir, Au |
| `not_measurable` | 21 | F, Cl, Ne, Ar, Kr, Xe, As, Se, Br, In, Sb, Te, I, Cs, Ta, Re, Pt, Hg, Tl, Bi, U |

| `codex_fit_verdict` | n | elements |
| --- | --: | --- |
| `cheap_to_adopt` | 1 | **Ag** |
| `feasible_high_effort` | 6 | Be, Rb, Pd, Nd, Pb, Th |
| `blocked_no_grid` | 25 | Ga, Ge, Nb, **Mo**, Ru, Rh, Cd, Sn, La, Ce, Pr, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, W, Os, Ir, Au |
| `blocked_spectroscopy` | 18 | B, F, Cl, As, Se, Br, In, Sb, Te, I, Cs, Ta, Re, Pt, Hg, Tl, Bi, U |
| `not_measurable` | 4 | Ne, Ar, Kr, Xe |

### `3D_NLTE` — two elements, one of them new

**Ag (Z=47)** is the CRITICAL find. Caliskan et al. 2026 (A&A 711, A155,
DOI `10.1051/0004-6361/202659578`, arXiv:2605.05356) present the *first* Ag I NLTE
model atom, with ab-initio multiconfiguration Hartree–Fock oscillator strengths, and
derive a solar 3D non-LTE A(Ag) = 1.15 ± 0.08 — **+0.19 dex above** the previous
reference, cutting the photosphere-versus-meteorite discrepancy from 0.25 to
0.06 dex. The diagnostics are the two Ag I resonance lines at 3280.68 Å and
3382.90 Å, both in the near-UV and both inside the Kitt Peak arm; Grevesse et al.
2015 measure them at 35 mÅ and 22 mÅ. See *Grid acquisition paths* for why this is
`cheap_to_adopt` rather than merely feasible.

**Be (Z=4)** is the one Lodders et al. 2025 already knew about (Amarsi et al. 2024,
A&A 690, A128). It is a cautionary case, not a green light: despite full 3D NLTE,
Lodders grade Be **C** and *inflate* the uncertainty from Amarsi's ±0.05 to ±0.14,
naming three reasons — a single Be II line at 3131.07 Å, a strong unidentified blend
at 3131.02 Å, and a UV radiation field that forms in the chromosphere, which no 1D
or 3D photospheric model carries. The spectroscopy is reachable (the line is in the
Kitt Peak near-UV arm); the modelling is what is hard, and no departure grid is
published in the Amarsi Zenodo family.

### `1D_NLTE` — six elements, three of them corrections to the first-pass table

* **Rb** — Korotin 2020 (Astron. Lett. 46, 541), 29-level Rb I model atom, solar
  NLTE A(Rb) = 2.35 ± 0.05 with a −0.12 dex correction. **Downgraded from the seed
  table's `cheap_to_adopt`** on two hard facts: no public grid exists (all three
  public NLTE services were enumerated — see below — and none serves Rb), and both
  resonance lines (7800.27, 7947.60 Å) sit past 7000 Å, outside the HARPS arm that
  carries our entire measured EW pool.
* **Pd** — **post-Lodders.** Mashonkina & Smogorzhevskii 2025 (A&A 703, A296,
  arXiv:2510.06968) build a new comprehensive Pd I model atom and report solar NLTE
  log ε = 1.61–1.70. Promotes Pd out of the seed table's LTE-only bucket. Note the
  0.09 dex internal spread across lines: that alone consumes our ±0.10 gate.
* **Pr** — **seed-table correction.** The first pass cited "1D NLTE, Lodders 2025
  p703". No such citation exists in the paper; Lodders assert only that "the only
  elements with NLTE estimates of abundances in 1D available are Sr, Y, Zr, Ba, Pr,
  Eu, Gd, and Nd" with no reference attached. The only Pr NLTE work this audit could
  locate is Mashonkina, Ryabchikova, Ryabtsev & Kildiyarova 2009 (A&A 495, 297) —
  **for A and Ap stars over Teff 7250–9500 K**, a hull that excludes the Sun
  (5772 K) and every codex target. A solar Pr NLTE number falls out of that paper as
  a by-product, so the tier is `1D_NLTE`; the verdict is `blocked_no_grid` because
  there is no grid we could interpolate for an FGK star.
* **Nd** — **post-Lodders and the largest rare-earth advance since 2015.** Dixon,
  Ezzeddine, Li, Merle, Bautista & Guo 2025 (ApJ 994, 44,
  DOI `10.3847/1538-4357/ae0b59`) build a Nd I/Nd II model atom and a grid of NLTE
  corrections for **122 Nd II lines from the UV to the H band**, with solar
  A(Nd) = 1.44 ± 0.05 at S_H = 0.1. Corrections run −0.3 to +0.3 dex and **change
  sign between bands** (positive optical/UV, negative H-band) — precisely the case
  the RYA-708 per-instrument-per-band law exists for.
* **Pb** — **seed-table correction.** The first pass swept "Hf–Bi" into LTE-only.
  Pb is the exception: Mashonkina, Ryabtsev & Frebel 2012 (A&A 540, A98, *"Non-LTE
  effects on the lead and thorium abundance determinations for cool stars"*) give a
  **+0.12 dex** correction at S_H = 0.1 on MARCS, which Grevesse et al. 2015 apply.
  The blocker is the line, not the physics — Pb I 3683.48 Å sits in the outer red
  wing of an extremely strong Co I + Fe I + V I blend and in the blue wing of
  another Fe I line.
* **Th** — same paper, +0.01 dex. Th is *not* LTE-blocked; it is blend-blocked. Of
  the 0.56 pm total feature at 4019.13 Å, Grevesse et al. attribute 0.208 pm to a
  Co I blend and 0.038 pm to a V I blend, leaving 0.314 pm for Th II with ~20%
  uncertainty on that share. We already have measured solar Co and V, which is
  exactly what a three-component deblend needs.

### `LTE_only` — 25 elements with a line and no model atom

The rare earths dominate (La, Ce, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu) and share one
root cause, stated by Grevesse et al. 2015 §4.21: *"We are not aware of any other
NLTE study of rare Earth elements in the Sun"* — the sole exception being Eu, which
is already canonical. Their abundances are Wisconsin-group 1D LTE analyses
(Lawler/Den Hartog/Sneden series) corrected only for the 3D−HM model difference.
**Gd is a seed-table correction**: it was listed as 1D NLTE on the strength of the
same uncited Lodders sentence as Pr, and this audit could locate no Gd II model atom
at all.

**Mo (Z=42) is the headline of this tier** and the only outside-28 element clearing
two science gates. The most recent dedicated study — Mishenina, Kurtukian-Nieto,
Gorbaneva, Amarsi, Psaltis & Pignatari 2026 (A&A 705, A38, arXiv:2511.21190),
submitted November 2025 — still works in LTE and closes with: *"we emphasize that
defining NLTE corrections for Mo and Ru abundance measurements is essential, as
these remain currently unknown."* There is a trap here worth recording: that paper
quotes a **+0.15 dex MPIA correction for the analogous Cr I 5208.41/5206.02 lines**
as a qualitative proxy for Mo I 5506.49/5533.03, and secondary sources misread this
as an MPIA Mo correction. MPIA Spectrum Tools was enumerated live on 2026-08-09 and
serves no Mo. Both Mo I lines are optical and inside the HARPS arm, so the day a Mo
model atom appears, Mo is immediately workable.

### `not_measurable` — 21 elements with no photospheric diagnostic

Three distinct mechanisms, and the distinction matters for how they are parked:

1. **Noble gases (Ne, Ar, Kr, Xe).** No photospheric lines at solar Teff at all.
   Lodders derive Ne from the Ne/O ratio of the quiet transition region, Ar from
   element ratios, Kr and Xe by interpolating s-process nuclide systematics. These
   are not measurements and never will be from a photospheric spectrum. Ne and Ar
   are α-elements, so they clear a science gate — and are still permanently out of
   reach. Kr/Xe carry a second problem: their "values" are s-process model
   interpolations, so using them to test s-process predictions would be circular.
2. **Sunspot-only species (F, Cl, In, Tl).** A solar value exists but comes from
   umbral spectra, which no FGK dwarf offers. In is worse than the label suggests:
   Grevesse et al. show the 4511.3 Å photospheric feature historically attributed to
   In I is contributed **<20%** by In I (Vitas et al. 2008).
3. **No identified solar line at all (As, Se, Br, Sb, Te, I, Cs, Ta, Re, Pt, Hg,
   Bi, U).** Lodders Table 2 carries "…" for each; the solar-system value is
   meteoritic. Lodders devote a whole section to Hg concluding the Sun would be the
   best source *because* meteoritic Hg is unreliable — and then record that no
   photospheric determination exists. **U is the one that stings**: there is no
   photospheric U, so any codex Th/U cosmochronometer would be half-photospheric,
   half-meteoritic. Say so or do not report it.

---

## Recommendations

**1. The canonical inclusion criterion should gate on grid tier, not on science
interest.** The evidence is one-sided: the four science gates (α / CNO /
bio-significant / n-capture) are cleared by **52 of 54** outside-28 elements at
level ≥1, by only 5 at level 2, and by **none at all** at level ≥3. They do not
discriminate. Grid tier does:
2 elements in 3D NLTE, 6 in 1D NLTE, 46 with neither. Concrete replacement language
for clause (c) is in
`docs/proposals/science_standards_grid_availability_amendment.md`, together with a
table classifying all 54 elements under it.

**2. Add an explicit `not_measurable` disposition, distinct from "owed".** Twenty-one
elements have no photospheric diagnostic. Parking them as refinement debt
(RYA-676) implies work that will never be done. They need a terminal state that says
*closed, and here is the physics that closes it* — otherwise every future session
re-litigates F and the noble gases exactly the way RYA-376 re-litigated Zn.

**3. Ag is the one expansion candidate worth opening a ticket for now.** Not on
science priority — Ag clears one gate, the same as 40 other elements — but because
it is the only element in the universe surveyed here where the physics, the grid
file, the file format, and the wavelength coverage are all simultaneously in hand.
See *Grid acquisition paths*.

**4. Do not treat "Lodders says an NLTE estimate exists" as a citation.** Three of
the seed table's rows (Pr, Nd, Gd) rested on one unreferenced sentence in Lodders
et al. 2025. Chasing it down produced one genuine grid (Nd, and only because a
paper appeared *after* Lodders), one out-of-hull A/Ap-star study (Pr), and nothing
at all (Gd). The audit-level rule is the pipeline-level rule: a claim is a
conclusion, not evidence.

**5. Two elements sit at gates_passed = 0 (Be, F).** Under any purely
science-gated criterion both are excluded — including Be, which has the best
line-formation physics of anything outside the 28. That is the clearest argument
that the criterion needs a grid-tier axis rather than more science gates.

---

## Grid acquisition paths

Only elements at `cheap_to_adopt` or `feasible_high_effort` have a path at all.
Everything else is blocked upstream of acquisition.

**All three public NLTE grid services were enumerated live on 2026-08-09.** Not
sampled, not remembered — the served-species lists were read out of the live HTML
forms and the Zenodo API:

| service | species served | any outside the 28? |
| --- | --- | --- |
| MPIA Spectrum Tools, NLTE abundance correction (`nlte.mpia.de/gui-siuAC_secE.php`) | H I, O I, Mg I, Si I, Ca I, Ca II, Ti I, Ti II, Cr I, Mn I, Fe I, Fe II, Co I | **none** |
| MPIA Spectrum Tools, Payne NLTE from SAPP (`gui-payne_fit.php`, Storm et al. 2026) | Li, C, O, Na, Mg, Al, Si, Ca, Ti, Cr, Mn, Fe, Co, Ni, Sr, Y, Ba, Eu | **none** |
| INASAN nLTE2 (`spectrum.inasan.ru/nLTE2/`, Mashonkina et al. 2023, MNRAS 524, 3526) | Ba II, Ca I, Ca II, Eu II, Fe I, Mg I, Na I, Sr II, Ti II, Zn I, Zn II | **none** |
| Amarsi group Zenodo `Grid/NLTE` (concept DOI `10.5281/zenodo.3888393`), all 8 versions | H, Li, C, N, O, Na, Mg, Al, Si, K, Ca, Mn, Ba, Fe, Ti, Cu, S, **Ag** | **Ag only** |

That table is the audit's most useful single result: **exactly one outside-28
element is served by any public NLTE grid infrastructure anywhere, and it is Ag.**

* **Ag — `cheap_to_adopt`.** `https://doi.org/10.5281/zenodo.20037437` (Grid/NLTE
  v8, published 2026-05-05, CC-BY-4.0, 4.4 GB total) ships
  `nlte_Ag_caliskan_Jan2026_pysme.grd` (960 MB), the IDLSME twin, a 2.5 GB text
  form, plus `atmos_Ag.txt` and `label_Ag.txt`. This is the **same PySME `.grd`
  family and the same first author** as our already-registered
  `Cu_Caliskan2024_PySME.csv` (Zenodo 15062813) and the Amarsi-2020 grids behind
  Na/Mg/Si/Al/K/N/S — so the derivation path is the one RYA-540 already executed.
  Route: stage on Sirius (`/mnt/codex-data/grids/nlte/`, **never the Mac**), derive
  per-line deltas in the PySME venv, register in `NLTE_CORRECTION_ELEMENTS`.
  Wavelengths 3280.68 / 3382.90 Å need the Kitt Peak near-UV arm, which we hold.
* **Be — `feasible_high_effort`.** No departure grid is published; Amarsi et al.
  2024 release Table 2 through CDS VizieR (`J/A+A/690/A128`). Acquisition is a
  direct request to the authors, and even with the grid the chromospheric radiation
  field and the 3131.02 Å blend remain unsolved by anyone.
* **Nd — `feasible_high_effort`, verification owed.** Dixon et al. 2025 state the
  122-line correction grid exists; a public download location is **not** stated in
  the abstract. Before any adoption, read the ApJ data-availability statement and
  check Zenodo. Recording this as owed rather than assuming a link exists.
* **Pd — `feasible_high_effort`.** Mashonkina & Smogorzhevskii 2025; INASAN nLTE2
  was enumerated and Pd is not among its 11 species. A request to the authors is
  owed — the same channel that RYA-433 found WAF-blocked for a *different* INASAN
  endpoint, so expect a manual pull.
* **Rb — `feasible_high_effort`.** No public grid; corrections tabulated in Korotin
  2020. The harder half is instrumental: 7800/7948 Å needs the Kitt Peak red arm for
  the Sun and a red-capable instrument for every other target.
* **Pb, Th — `feasible_high_effort`.** No grid to acquire at all: Mashonkina et al.
  2012 publish per-line corrections (+0.12 and +0.01 dex), not an interpolable grid.
  What these two need is not acquisition but the in-window blend-fit machinery built
  for RYA-581/585.

### Grid-loader compatibility spot-check (integration path proven)

`scripts/spotcheck_nlte_grid_loader_rya758.py`. The ticket names
`pipeline.nlte_grid_loader`; **no such module exists**. The live entry point is
`pipeline.nlte_corrections._load_mpia_element_grid`, registry-driven off
`config.constants.NLTE_CORRECTION_ELEMENTS` and guarded by the RYA-413
placeholder-zero refusal, reached in production through
`apply_element_nlte_corrections` / `element_grid_in_bounds`. Two already-registered
controls were loaded and interpolated at the solar node — no new grid was downloaded
(grids install on Sirius only):

```
[OK] Ca  file=Ca_Mashonkina2017.csv
     raw table shape      : (576, 6)  columns=['element','wave_A','teff_K','logg','feh','delta_nlte']
     interpolated lines   : 7  [4578.56, 5867.57, 6122.22, 6424.585, 6439.084, 6449.82, 6745.23]
     (teff,logg,feh) hull : {'teff': (5000.0, 6500.0), 'logg': (3.5, 4.5), 'feh': (-0.5, 0.3)}
     solar inside hull    : True
     delta_nlte @  6122.220 A, solar = -0.0162 dex

[OK] Ti  file=Ti_Mallinson2024_PySME.csv
     raw table shape      : (33, 7)  columns=['element','ion','wave_A','teff_K','logg','feh','delta_nlte']
     interpolated lines   : 3  [5648.565, 5662.15, 5689.46]
     (teff,logg,feh) hull : {'teff': (5100.0, 6200.0), 'logg': (4.0, 4.7), 'feh': (-0.3, 0.6)}
     solar inside hull    : True
     delta_nlte @  5689.460 A, solar = +0.0504 dex
```

The Ti control returns +0.0504 dex at 5689.46 Å, reproducing the registered
RYA-545 solar median of +0.0506 — so the loader is not merely opening files, it is
returning the physics we banked.

---

## Post-Lodders-2025 developments

Lodders, Bergemann & Palme 2025 was accepted 2025-02-07, so anything submitted from
2025-03-01 is outside its evidence base. `scripts/scan_arxiv_post_lodders_rya758.py`
issued five broad arXiv API queries (`abs:"non-LTE" AND cat:astro-ph.SR`,
`abs:"NLTE" AND cat:astro-ph.SR`, `abs:"model atom"`, `abs:"departure coefficients"`,
`abs:"solar abundance"`), kept the **173** unique submissions on or after the cutoff,
and filtered them for whole-word outside-28 element names; the raw hit list is
`data/audit/nlte_grid_inventory_beyond28_arxiv_scan.json`. Targeted searches for F,
Mo, Rb, Be, Pr, Gd and the rare earths were run separately.

Found:

| element | development | why it matters |
| --- | --- | --- |
| **Ag** | Caliskan, Amarsi, Jönsson, Grevesse & Sahoo 2026, A&A 711, A155 (arXiv:2605.05356) — first Ag I NLTE model atom, 3D NLTE solar A(Ag) = 1.15 ± 0.08; grid on Zenodo 2026-05-05. Companion atomic data: Jönsson, Sahoo, Caliskan & Amarsi 2026, A&A 709, A31. | **CRITICAL** — a 3D NLTE element outside the 28 that post-dates Lodders. |
| **Nd** | Dixon et al. 2025, ApJ 994, 44 — Nd I/Nd II model atom + 122-line correction grid, solar A(Nd) = 1.44 ± 0.05. | Promotes Nd to `1D_NLTE`; first real rare-earth NLTE grid for FGK stars. |
| **Pd** | Mashonkina & Smogorzhevskii 2025, A&A 703, A296 (arXiv:2510.06968) — new Pd I model atom, solar NLTE log ε = 1.61–1.70. | Promotes Pd to `1D_NLTE`. |
| **Mo, Ru** | Mishenina et al. 2026, A&A 705, A38 (arXiv:2511.21190) — 154 giants, Mo and Ru still derived in **LTE**; explicit statement that NLTE corrections "remain currently unknown". | Confirms `blocked_no_grid` as the *current* state, with Amarsi as a co-author. |
| **B** | Spite, Barbuy & Tan 2025 (arXiv:2510.11594) — a measurable B I line at 2089.6 Å in metal-poor stars, HST/STIS. | A second B diagnostic, still far-UV; confirms `blocked_spectroscopy`. |
| **Rb** | arXiv:2506.21332 — Rb abundances in cool giants from H-band spectra. | A new IR Rb diagnostic, but for giants; our targets are dwarfs. |

Adjacent, inside the 28 but worth flagging for the sweep: **arXiv:2511.14289**,
*"Solar carbon abundance from 3D non-LTE modelling of the diagnostic lines of the CH
molecule"* — Lodders explicitly declined molecular C indicators *because* no NLTE
modelling of C₂/CH/CO existed. That gap is now partly closed, which bears on the
codex C leg (RYA-237/359). Also **Storm et al. 2026**, the Payne NLTE model behind
the new MPIA fitting service, covering 18 elements all inside our 28.

**Method limits, stated rather than hidden.** No ADS API token is available in this
environment, so the ADS half of the ticket's source #2 was substituted with the
arXiv API sweep above plus targeted web searches; a token-backed ADS query would
cover non-arXiv venues (Astronomy Letters, ApJS) that this sweep can miss. Five
publisher sites (aanda.org full text, journals.uchicago.edu) return HTTP 403 to
programmatic clients — every affected claim was therefore verified against the arXiv
copy or the Crossref record instead, never assumed. All 22 DOIs in the CSV were
resolved through the Crossref API; four initially-recorded identifiers were **wrong
and were corrected** (Ag `…202559578`→`…202659578`; Mashonkina 2012
`…201118337`→`…201218790`; Tb `10.1086/322540`→`10.1086/323001`; Gd
`10.1086/507064`→`10.1086/508262` — the last two had resolved to entirely unrelated
papers).

**Not found, and the negative is the result.** No NLTE model atom or departure grid
for: Mo, Ru, Rh, Ga, Ge, Nb, Cd, Sn, W, Os, Ir, Au, La, Ce, Sm, Gd, Tb, Dy, Ho, Er,
Tm, Yb, Lu, Hf, or any halogen.

---

## Column conventions

`grid_tier` records the best *published* line-formation treatment for the element's
solar diagnostic, independent of whether we could obtain the file;
`codex_fit_verdict` records what it would take *this pipeline* to derive the element.
The two differ on purpose — Pr is `1D_NLTE` / `blocked_no_grid` because a solar NLTE
number exists but no FGK-applicable grid does.

`bio_significant` is a **codex-internal classification**, not a quantity lifted from
a paper: it marks an established essential role in terrestrial biochemistry, and the
role is named per row in `notes` (B — borate cross-linking of plant cell-wall
rhamnogalacturonan-II; Cl — chloride electrolyte and the photosystem-II
oxygen-evolving complex; Se — selenocysteine and glutathione peroxidase; Br —
peroxidasin sulfilimine cross-links in collagen IV; Mo — the nitrogenase
FeMo-cofactor and molybdopterin enzymes; I — thyroid hormones T3/T4; W —
tungstoenzymes in anaerobic archaea and bacteria). It carries no DOI because it is
a judgement this project is making, and it is flagged as such so nobody later mines
it as if it were sourced.

`canonical_criterion_gates_passed` is the four-gate sum (α / CNO / bio-significant /
n-capture). `cno` is `False` on every row by construction — C, N and O are all
canonical already.

---

## References

Sources are listed in the order the audit relied on them. Every DOI below was
resolved through the Crossref API on 2026-08-09.

1. **Lodders, K., Bergemann, M. & Palme, H. 2025**, *Solar System Elemental
   Abundances from the Solar Photosphere and CI-Chondrites*, Space Sci. Rev. 221, 23
   (accepted 2025-02-07). DOI `10.1007/s11214-025-01146-w`; arXiv:2502.10575. Local
   copy: `Reference documents/AbundancesLodders.pdf`. Table 2 (quality index) and the
   per-element text are the backbone of the tier assignments.
2. **Grevesse, N., Scott, P., Asplund, M. & Sauval, A. J. 2015**, *The elemental
   composition of the Sun III. The heavy elements Cu to Th*, A&A 573, A27.
   DOI `10.1051/0004-6361/201424111`; arXiv:1405.0288. Source of the diagnostic-line
   wavelengths, equivalent widths and the heavy-element NLTE corrections.
3. **Caliskan, S., Amarsi, A. M., Jönsson, P., Grevesse, N. & Sahoo, B. K. 2026**,
   *Ag I model atom and the 3D non-LTE solar silver abundance*, A&A 711, A155.
   DOI `10.1051/0004-6361/202659578`; arXiv:2605.05356.
4. **Amarsi, A. M., Ogneva, D., Buldgen, G., Grevesse, N., Zhou, Y. & Barklem,
   P. S. 2024**, *The solar beryllium abundance revisited with 3D non-LTE models*,
   A&A 690, A128. DOI `10.1051/0004-6361/202451778`; arXiv:2408.13105.
5. **Dixon, J. D., Ezzeddine, R., Li, Y., Merle, T., Bautista, M. & Guo, Y. 2025**,
   *Investigating non-LTE abundances of Neodymium (Nd) in metal-poor FGK stars*,
   ApJ 994, 44. DOI `10.3847/1538-4357/ae0b59`; arXiv:2509.22811.
6. **Mashonkina, L. & Smogorzhevskii, A. 2025**, *Understanding an origin of
   palladium in metal-poor stars based on the non-LTE analysis of Pd I lines*,
   A&A 703, A296. arXiv:2510.06968 (DOI `10.48550/arXiv.2510.06968`).
7. **Korotin, S. A. 2020**, *Non-LTE Effects in Rubidium Lines in Cool Stars*,
   Astron. Lett. 46, 541. DOI `10.1134/S1063773720080022`.
8. **Mashonkina, L., Ryabtsev, A. & Frebel, A. 2012**, *Non-LTE effects on the lead
   and thorium abundance determinations for cool stars*, A&A 540, A98.
   DOI `10.1051/0004-6361/201218790`.
9. **Mashonkina, L., Ryabchikova, T., Ryabtsev, A. & Kildiyarova, R. 2009**,
   *Non-LTE line formation for Pr II and Pr III in A and Ap stars*, A&A 495, 297.
   DOI `10.1051/0004-6361:200810258`.
10. **Mishenina, T., Kurtukian-Nieto, T., Gorbaneva, T., Amarsi, A. M., Psaltis, A.
    & Pignatari, M. 2026**, *Molybdenum and ruthenium in the Galactic disk*,
    A&A 705, A38. arXiv:2511.21190 (DOI `10.48550/arXiv.2511.21190`).
11. **Maiorca, E., Uitenbroek, H., Uttenthaler, S. et al. 2014**, *A new solar
    fluorine abundance and a fluorine determination in the two open clusters M 67
    and NGC 6404*, ApJ 788, 149. DOI `10.1088/0004-637X/788/2/149`.
12. **Lodders, K. & Fegley, B. 2023**, *Solar system abundances and condensation
    temperatures of the halogens fluorine, chlorine, bromine, and iodine*,
    Geochemistry 83, 125957. (The halogen source Lodders et al. 2025 defer to.)
13. **Cunha, K. & Smith, V. V. 1999**, *A determination of the solar photospheric
    boron abundance*, ApJ 512, 1006. DOI `10.1086/306796`.
14. **Spite, M., Barbuy, B. & Tan, K. 2025**, *Precise boron abundance in a sample of
    metal-poor stars from far-ultraviolet lines*, arXiv:2510.11594.
15. Wisconsin-group rare-earth 1D LTE analyses adopted by Grevesse et al. 2015:
    La — Lawler, Bonvallet & Sneden 2001a, ApJ 556, 452, DOI `10.1086/321549`;
    Ce — Lawler, Sneden, Cowan, Ivans & Den Hartog 2009, ApJS 182, 51,
    DOI `10.1088/0067-0049/182/1/51`;
    Sm — Lawler, Den Hartog, Sneden & Cowan 2006, ApJS 162, 227, DOI `10.1086/498213`;
    Gd — Den Hartog, Lawler, Sneden & Cowan 2006, ApJS 167, 292, DOI `10.1086/508262`;
    Tb — Lawler, Wickliffe, Cowley & Sneden 2001c, ApJS 137, 341, DOI `10.1086/323001`;
    Ho — Lawler, Sneden & Cowan 2004, ApJ 604, 850, DOI `10.1086/382068`;
    Er — Lawler, Sneden, Cowan et al. 2008, ApJS 178, 71, DOI `10.1086/589834`;
    Hf — Lawler, Den Hartog, Labby et al. 2007, ApJS 169, 120, DOI `10.1086/510368`;
    Pr, Dy, Tm, Yb, Lu — Sneden, Lawler, Cowan, Ivans & Den Hartog 2009, ApJS 182,
    80, DOI `10.1088/0067-0049/182/1/80`.
16. **Amarsi, A. M. et al. 2020**, *Non-LTE departure coefficients for large
    spectroscopic surveys*, A&A 642, A62 — the method paper behind the Zenodo
    `Grid/NLTE` family (concept DOI `10.5281/zenodo.3888393`; Ag version
    `10.5281/zenodo.20037437`, Cu `10.5281/zenodo.15062813`, Ti
    `10.5281/zenodo.10753497`, S `10.5281/zenodo.17064337`, Fe
    `10.5281/zenodo.7088951`, 2020 base `10.5281/zenodo.3982506`).
17. **Mashonkina, L., Pakhomov, Yu., Sitnova, T. et al. 2023**, MNRAS 524, 3526 —
    the INASAN nLTE2 database citation, per its own References page.
18. **Kovalev, M., Brinkmann, S., Bergemann, M. & MPIA IT 2018**, *NLTE MPIA web
    server*, `http://nlte.mpia.de` — the MPIA Spectrum Tools citation, per its own
    References page.
19. **Storm, N. et al. 2026** — the Payne NLTE model library behind MPIA's Spectrum
    Fitting service (cited by the tool itself). Earlier: **Storm, N., Barklem, P. S.,
    Yakovleva, S. A. et al. 2024**, *3D NLTE modelling of Y and Eu*, A&A 683, A200.
