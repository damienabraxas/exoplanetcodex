# RYA-1061 — UV completeness audit

This is a read-only evidence audit for **1150 <= wavelength < 3780 A**. It changes no
canonical gf row, line pool, holding, product, abundance, or astrophysical-gf authority.

## Verdict

The committed canonical list contains **21,279** rows in the audit window, all in
3000–3780 A. FUV and NUV each contain **zero** canonical rows. The near-UV list is
20,568/21,279 VALD3 rows; only
71 rows carry the primary-laboratory tier. This is a
coverage discontinuity, not evidence that no transitions exist below 3000 A.

RYA-836 remains the verified Fe I anchor: 105 primary-lab Fe I transitions below 3780 A
in its source census, of which 61 lie in 3000–3780 A. The current canonical inventory is
not equivalent to that source census (it carries 71 total LAB rows across species), so
source availability, canonical adoption, and target reachability must remain separate.

The highest-leverage acquisition is Ward et al. (2023) Cr II (268 published transitions
over 208–414 nm, clipped here at 378 nm), followed by transition-table verification for
Mn II, Co II, Sc I/II, and V II. These are **source candidates**, not approved canonical
rows. No line-level adoption verdict is possible until tables are acquired, wavelength
frames normalized, ions matched, HFS/isotope structure retained, and ambiguous matches
reported.

## Target and literature routing

The holdings registry has four UV holdings: Alpha Cen A/STIS, Alpha Cen B/STIS,
Procyon/STIS, and Procyon/COS. ASTRAL is supporting spectral/identification evidence.
Pagano et al. is routed to the chromosphere/TR/corona lane. The Alpha Cen Be II 3130 work
is a candidate photospheric-abundance source requiring exact bibliographic and gf/blend
verification. The Procyon GHRS atlas is supporting chromospheric/identification evidence;
this census found no Procyon photospheric-abundance paper to pair with its UV spectra.

## Files

- `current_uv_inventory.csv`: per-species, per-regime canonical counts by gf tier.
- `laboratory_source_census.csv`: source-level evidence and acquisition routing.
- `stellar_uv_paper_census.csv`: strict science-lane and authority classification.
- `target_uv_holdings.csv`: UV rows derived from the holdings registry.
- `target_literature_overlap.csv`: target-level spectra/literature gap.
- `candidate_gap_analysis.csv`: prioritized gaps and next actions.
- `per_species_verdict.csv`: current-state canonical verdicts only.
- `summary.json`: machine-readable headline and audit limits.

Reproduce with `python3 scripts/rya1061_uv_completeness_audit.py`.
