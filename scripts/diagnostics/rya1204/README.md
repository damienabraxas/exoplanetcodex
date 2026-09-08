# RYA-1204 — the two near-UV opacity levers, each measured alone

Diagnostic. Nothing here is a product; nothing was calibrated to a target (RYA-161).

RYA-1202 established the near-UV Fe deficit is **under-tabulated opacity, not missing
physics**, and named two independent candidates. This measures both.

## Step 0 — Fe I bf against Bautista/NORAD: **STOP, it is not the lever**

`parse_norad.py` reads NORAD's level-resolved Fe I photoionization
(`fe1.px.txt`, Bautista & Pradhan — the Iron Project lineage the ticket names) and
`step0_final.py` builds the LTE-weighted total to compare against `jonabs` component 17,
extracted by `jonabs_fe1.py`.

**🔴 The normalisation had to be measured before anything could be compared.** jonabs
stores `sum_i g_i exp(-E_i/kT) sigma_i` — *without* dividing by the partition function.
Comparing before establishing that showed a spurious factor ~40 and would have "proved" a
lever that does not exist. It was caught because the ratio jonabs/(Bautista/U) tracks U
across the whole temperature grid (25.06 at 1000 K → 47.14 at 8000 K), which is a property
no physical discrepancy would have.

**🔴 The ground state cannot contribute at all.** Its binding energy is 7.755 eV —
threshold 1599 Å. Every photon in 3000–3780 Å is far too soft. The band's Fe I bf opacity
comes entirely from levels **3.93–7.72 eV above ground**, which is why the component is
small and temperature-steep.

Result at 5900 K over 3000–3780 Å:

| comparison | ratio |
| --- | --- |
| under-sampling alone (Bautista full vs Bautista at jonabs's 4 nodes, linear) | **1.0006x** |
| absolute scale at the nodes (jonabs / Bautista) | 1.19x median, 0.73–1.77x |
| what the code uses vs full-resolution Bautista, band mean | **1.44x** |

The 4-node linear interpolation loses **0.06%** of the band mean. The node placement is
not the problem, and our table already carries *more* Fe I bf opacity than a re-tabulation
from Bautista would give. Bell 2001 needed 2x **more**; re-tabulating would take us ~30%
the **wrong way**. Per the ticket's own Step 0 gate: **stop, do not build Lever A.**

⚠️ Caveat, stated because it bounds the claim: NORAD's file carries 421 LS states where
Bautista 1997 reports 1,117. The missing ones are high-lying — exactly those that ionize
in the near-UV — so the Bautista total here may be an underestimate and the precise 1.44x
is soft. It cannot reverse the verdict: more states move Bautista *up* toward jonabs, and
the finding is that jonabs is not 2x **low**.

Plot: `data/results/rya1204/rya1204_step0_fe1bf.png`.

## Lever B — near-UV molecular opacity: **this is the lever**

`build_mol.py` writes the 6,223 molecular lines the near-UV list excludes (NH 2271,
OH 2094, CN 1205, CH 653) as Turbospectrum `.bsyn` lists; `leverB_product.py` re-runs the
RYA-759 near-UV Fe product with them on and off, patching only `use_molecules` and the
molecules directory.

On the isotopologue — the reason they were excluded — the assignment is the dominant one
and the error it can introduce is **bounded and tiny**: 16O is 99.76% of solar O, 14N
99.64%, 12C 98.93%. If VALD's gf is the all-isotopologue total, assigning it to the
dominant one over-counts by ≤1.1%; if it is already per-dominant-isotopologue, it is
exact. The status quo is a 100% error — the band carries no molecular opacity at all.
The isotopologue cannot move a line here: these wavelengths are VALD's, not computed
term values, so it enters only through the abundance.

The measurement validates itself three ways beyond the abundance shift: every window the
script flags `[continuum off]` comes back into range, reduced chi2 falls (528 → 63 on the
worst line), and line-to-line scatter drops. A window containing no molecular lines
(3756.937 Å) returns **bit-identical** results, which controls the patch.

## Reproduce

    python scripts/diagnostics/rya1204/parse_norad.py      # needs fe1.px.txt from NORAD
    python scripts/diagnostics/rya1204/step0_final.py
    python scripts/diagnostics/rya1204/build_mol.py
    USE_MOL=1 python scripts/diagnostics/rya1204/leverB_product.py

`fe1.px.txt` / `fe1.px.gd.txt` are fetched from
`https://norad.astronomy.osu.edu/fe1/` and are not vendored here.

## Results — both levers measured alone, on the production route

RYA-759 near-UV Fe I product, `--limit 40`, ATLAS9.Castelli against the Kitt Peak atlas.
All three arms are the SAME 40 lines, so the comparison is paired.

| arm | A(Fe I) | scatter | continuum-clean | shift |
| --- | --- | --- | --- | --- |
| baseline (production today) | 7.502 | 0.374 | 36 / 40 | — |
| **Lever A** — Fe I bf re-tabulated from Bautista | 7.515 | 0.394 | 35 / 40 | **+0.013** |
| **Lever B** — near-UV molecular opacity | **7.415** | 0.337 | **40 / 40** | **−0.087** |

**Lever A does nothing detectable.** Paired: mean +0.076, SEM 0.095, **t = +0.80**, median
−0.004, and it moves 23 of 40 lines down — a coin flip. It is not an adverse lever, it is
an absent one, which is what Step 0 predicted from the cross-sections alone.

**Lever B is real.** Paired: mean −0.156, SEM 0.036, **t = −4.30**, 31 of 40 lines move
down and 6 are bit-identical (no molecular line in the window — the control).

It also validates itself independently of the abundance:

* all **4** windows flagged `[continuum off]` are repaired; **0** remain
* median reduced chi2 **118 → 81**; the worst line 3175.310 goes **529 → 63**
* line-to-line scatter **0.374 → 0.337**
* mean |synth/obs − 1| **0.108 → 0.080**

⚠️ **And it UNDER-corrects.** The signed continuum ratio is still above 1 on **24 of 40**
windows after the fix (mean 1.067 → 1.023). The synthesis still sits above the observed
pseudo-continuum, so opacity is still owed and −0.087 is a **lower bound**, not a
calibration that was stopped when it reached a target.
