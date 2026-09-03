# RYA-1169 Si intake checkpoint

This is an evidence-first, reproducible checkpoint for the Silicon intake. It
reconstructs the AGSS21 reference membership, ingests the verified modern
primary-laboratory UV/blue measurements, and joins both to `canonical_gf` by
species + wavelength + excitation potential.

The output is deliberately **BLOCKED**, not `FROZEN_READY_FOR_MEASUREMENT`.
No laboratory row is called `graded` or `deepgraded` until its feature depth is
measured from an approved Solar spectrum. The broader historical primary-lab
census and row-level NIST comparison also remain open.

Run:

```bash
python3 scripts/build_si_intake_rya1169.py
python3 -m pytest tests/test_si_intake_rya1169.py
```

`literature_sources.csv` contains clickable DOI/article/arXiv links for every
paper used in this checkpoint, including direct download landing pages.

The follow-on Asplund membership/depth test is built with
`scripts/grade_si_asplund_rya1169.py` and written under
`data/results/rya1169/`. Its all-band matrix keeps UV and IR empty: the adopted
source set contains no UV rows, and Scott et al. explicitly discarded the
available near-IR Si I population because only low-accuracy theoretical gf and
large/uncertain NLTE effects were available. Those gaps must not be filled with
generic Si candidates and still called Asplund Grade.

Scientific lineage: AGSS21 adopts 7.51 dex. Amarsi & Asplund (2017) derive
7.51 dex by applying line-by-line 3D non-LTE corrections to the nine Si I
Solar indicators of Scott et al. (2015), then combining them with Scott's one
Si II indicator (whose correction is assumed zero). Two violet Si I test lines
in Amarsi & Asplund are retained here as explicitly rejected from the final
Solar mean. Scott's gf values use Garz (1973) relative measurements shifted by
+0.097 dex using O'Brian & Lawler (1991a,b) lifetimes.
