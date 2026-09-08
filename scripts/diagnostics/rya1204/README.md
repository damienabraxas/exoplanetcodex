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
