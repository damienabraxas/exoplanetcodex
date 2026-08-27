# RYA-1062 — UV ion-roadmap cross-scan

Read-only reconciliation of the committed ion tracker, frozen solar gold, RYA-1061 UV
inventory, and source-level laboratory leads. No canonical gf, line pool, ion label,
product, abundance, or roadmap state is changed.

## Verdict

The live tracker has 27 ion/element rows. Frozen gold disagrees with the
measured ion for **5** elements; these are routed to RYA-683 and remain
unmodified here. Five scientifically plausible ion stages (P II, S II, Ni II, Cu II,
Zn II) are absent from the roadmap and require a photospheric-line/data-reach screen
before tickets can promote them.

No surveyed ion is declared `complete`. The modern primary-lab tables are generally not
represented as LAB-tier canonical UV rows. Cr II is the clearest material opportunity:
Ward et al. (2023) publishes 268 transitions over 208–414 nm, while Codex has 826 Cr II
near-UV rows and zero primary-lab rows. Sc II, Ti II, V II, Mn II, Co II, and the
rare-earth majority ions also retain source-level reconciliation work.

`UNVERIFIED_SEED` and `UNRESOLVED_PRIMARY_FAMILY` are deliberately not evidence grades.
They are stop signs: exact bibliography, tables, wavelength frame, physical transition,
uncertainty, and HFS/isotope provenance must be resolved before adoption.

Reproduce with `python3 scripts/rya1062_uv_ion_roadmap_audit.py`.
