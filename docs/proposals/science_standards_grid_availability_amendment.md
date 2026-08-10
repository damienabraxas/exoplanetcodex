# PROPOSAL — canonical inclusion criterion, grid-availability clause

**Status: DRAFT, awaiting Ryan's ratification. Not applied.** This is a proposed
addition to `docs/SCIENCE_STANDARDS.md`; that file is untouched by RYA-758.

**Provenance:** RYA-758 (audit), amending the draft criterion sketched in RYA-757
§"Companion — canonical inclusion criterion". Evidence:
`docs/audit/nlte_grid_inventory_beyond28_report.md` and
`data/audit/nlte_grid_inventory_beyond28.csv` (54 elements, audited 2026-08-09).

---

## 1. Why the criterion needs sharpening at all

RYA-376 excluded Zn on "not bio-significant" grounds, which was wrong, and RYA-757
reversed it. The reversal was correct but it exposed the real defect: **there was no
written criterion**, so the question got re-argued from scratch. A criterion that
only says "a registered NLTE grid or a ratified LTE-acceptable disposition exists"
does not prevent the re-argument, because it does not say what counts as a grid,
what counts as a disposition, or who decides.

The RYA-758 audit also shows the *science* half of the criterion cannot carry the
decision. Across the 54 elements outside the canonical 28:

* **52 of 54** clear at least one of the four science gates (α / CNO /
  bio-significant / n-capture).
* Only **5** clear two (Se, Br, Mo, I, W).
* **None** clear three or more.

The gates do not discriminate. Grid tier does — 2 elements in 3D NLTE, 6 in 1D NLTE,
25 LTE-only, 21 with no photospheric diagnostic at all. And two elements
(**Be** and **F**) clear *zero* science gates, one of which — Be — has the best
line-formation physics of anything outside the 28. A purely science-gated criterion
would exclude the element with the best physics and admit forty-seven that have no
grid.

## 2. Proposed replacement text for clause (c)

> **(c) Grid availability.** An element enters the canonical set only if a published
> **1D NLTE model atom together with a departure-coefficient grid or per-line
> correction table** exists in the peer-reviewed literature, **and that grid's
> (Teff, log g, [Fe/H]) hull contains the star being analysed**. Where 3D NLTE is
> absent, the element enters with 3D flagged as a registered future refinement and
> the 1D-versus-3D uncertainty carried explicitly in the reported budget, never
> silently absorbed.
>
> **LTE-only admission** is permitted in exactly one case: the element is
> scientifically necessary and has no NLTE route in existence (not merely
> unacquired), **and** the LTE-versus-NLTE systematic is estimated from the nearest
> analogous modelled species and propagated into the reported uncertainty as a named
> term. An LTE-only admission is a ratified exception with its reasoning on record,
> not a default.
>
> **Terminal exclusion.** An element with no photospheric diagnostic in an FGK
> spectrum — no identified line, or a value obtainable only from sunspot umbrae,
> solar wind, or nucleosynthesis interpolation — is recorded `not_measurable` and is
> **closed, not owed**. It does not enter the refinement-debt registry, because there
> is no work that would ever discharge it.

### The three changes from the RYA-757 draft, and why each

1. **"a registered NLTE grid" → "a model atom *and* a grid or correction table,
   whose parameter hull contains the star".** The hull clause is not pedantry; the
   audit found a live instance. Praseodymium has a published Pr II/Pr III NLTE model
   atom (Mashonkina et al. 2009, A&A 495, 297) — for **A and Ap stars over
   Teff 7250–9500 K**. The Sun is 5772 K. Without the hull clause, Pr reads as
   "NLTE available" and would pass.
2. **"a ratified LTE-acceptable disposition" → a two-part test with propagation.**
   The old wording licenses "we decided LTE is fine". The replacement requires that
   the LTE-versus-NLTE gap be *estimated and carried*. The audit shows how that is
   done in practice: Mishenina et al. 2026 estimate Mo's likely NLTE sensitivity from
   the analogous Cr I 5208.41/5206.02 transitions (≈+0.15 dex) precisely because Mo
   has no model atom. An LTE admission that cannot even do that has no error budget.
3. **A new terminal `not_measurable` state.** Twenty-one of 54 elements have no
   photospheric diagnostic. Parking them as "owed" implies work that will never
   happen and guarantees they are re-litigated every time someone reads the tracker.
   They need a closed state with the physics on record.

## 3. How the amendment classifies all 54 outside-28 elements

Rows are the audit's `codex_fit_verdict`; the amendment's disposition follows
mechanically from `grid_tier` plus the hull and diagnostic tests.

| Amendment disposition | n | elements | basis |
| --- | --: | --- | --- |
| **Admissible now** (grid exists, hull contains our stars, acquisition path known) | 1 | **Ag** | 3D NLTE (Caliskan et al. 2026, A&A 711, A155); departure grid public at `10.5281/zenodo.20037437` in the PySME format we already ingest for Cu; both resonance lines 3280.68/3382.90 Å inside the Kitt Peak arm |
| **Admissible on acquisition** (physics done; a grid needs obtaining, or a line needs harder measurement) | 6 | Be, Rb, Pd, Nd *(grid acquisition)*; Pb, Th *(in-window blend-fit)* | see report §Grid acquisition paths |
| **Ineligible — no NLTE route in existence** | 26 | B, Ga, Ge, Nb, **Mo**, Ru, Rh, Cd, Sn, La, Ce, Pr, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, W, Os, Ir, Au | no model atom located; Pr fails on the hull clause specifically; B additionally far-UV/STIS-only |
| **Terminally excluded — `not_measurable`** | 21 | F, Cl, Ne, Ar, Kr, Xe, As, Se, Br, In, Sb, Te, I, Cs, Ta, Re, Pt, Hg, Tl, Bi, U | no FGK photospheric diagnostic: noble gases (no lines), sunspot-only (F, Cl, In, Tl), or no identified solar line at all |

1 + 6 + 26 + 21 = 54. B sits in "ineligible" rather than "terminally excluded"
because its lines do exist and are reachable with HST/STIS on our non-solar targets —
it is grid-blocked *and* instrument-expensive, not undiagnosable.

### The cases that test the wording

* **Mo** — the strongest science case outside the 28 (nitrogenase FeMo-cofactor;
  the only element clearing two gates that also has usable optical lines inside the
  HARPS arm) and it is still **ineligible**, because Mishenina et al. 2026 state
  flatly that Mo NLTE corrections "remain currently unknown". This is the amendment
  doing its job: science interest does not buy an exemption, and the LTE-only escape
  hatch would require a propagated Cr-analogue systematic that no one has published
  for the Sun. Mo is the element to watch (see `docs/OPEN_QUESTIONS.md`).
* **Th** — the cosmochronology anchor the RYA-757 draft named as the archetypal
  LTE-tolerant admission. The audit says the draft had the physics backwards: Th II
  4019.13 Å *does* carry a published NLTE correction (+0.01 dex, Mashonkina et al.
  2012), so it never needed an LTE exemption. Its real blockers are a three-way
  blend (of the 0.56 pm feature, 0.208 pm is Co I and 0.038 pm is V I) and the fact
  that Lodders et al. 2025 bracket the value as "listed for reference" only. Worse,
  **U has no photospheric value at all**, so a Th/U ratio would be half
  photospheric and half meteoritic — which the amendment's propagation requirement
  would force us to state out loud.
* **Be** — passes clause (c) on the best physics available anywhere outside the 28
  and passes **zero** science gates. Whatever the final criterion looks like, it must
  be able to say yes to Be, or it is measuring the wrong thing.
* **Ne, Ar** — α-elements, so they clear a science gate, and they are permanently
  unmeasurable from a photospheric spectrum. Terminal exclusion has to be able to
  override a science gate or it is not terminal.

## 4. What ratifying this would oblige

1. Insert clause (c) into `docs/SCIENCE_STANDARDS.md` under a named canonical
   inclusion criterion heading (RYA-179 doc-drift bundle owns the surrounding
   revision).
2. Add `not_measurable` as a terminal disposition alongside the existing tiers in
   `data/audit/element_refinement_registry.csv` (RYA-676 owns that schema), so the
   21 terminally excluded elements can be recorded closed rather than owed.
3. Nothing else. No element enters or leaves the canonical set on this proposal;
   Ag would need its own intake ticket, on the RYA-757 pattern.
