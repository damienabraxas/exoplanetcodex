# RYA-1134 source adjudication

This review consumes the frozen RYA-1132 assets and the later RYA-1173 primary-source
extractions. It does not use an abundance to choose a gf. The local document inventory
under `data/reference/asplund2021_al/` records the inspected PDF identities.

## Burheim, Hartman & Nilsson 2023

Source: A&A 672 A197, DOI **10.1051/0004-6361/202245394**, Tables 1 and 3, page 4
and page 5 of the seven-page PDF `aa45394-22.pdf`.

`burheim_table1_quantities.csv` transcribes the independent Einstein A values,
branching fractions, lifetimes and lifetime references for the twelve already-held
gf rows. The generator reproduces log(gf) from A, upper J, and vacuum wavelength;
the finite tolerance reflects the printed A and log(gf) precision. Adopted values
remain the published log(gf), not extra digits inferred from rounded A values.

Eight transitions use experimentally measured lifetimes from Buurman et al. 1986
or Buurman & Dönszelmann 1990. Four transitions ending in 4f use Papoulia 2019
**theoretical** lifetimes (Table 1 note d). Those are mixed measured-branching /
theoretical-lifetime evidence, not primary-laboratory-only gf. Three of these four
were labeled primary laboratory in the frozen manifest; the fourth is the held
11254.9 feature. RYA-946's September 3 policy allows sourced mixed/theoretical
reference inputs without changing the laboratory-only Codex/Deep rule.

Table 2 contains branching fractions without derived gf and supplies no adopted
gf here. The 41920.7 vacuum-A line has the paper's omitted weak fine-structure
blend caveat. The 11258.008 vacuum-A gf is a fine-structure transition, not the
frozen unresolved feature total at 11254.9239 air A; that feature remains HOLD.

The positive log uncertainty is `log10(1 + fractional uncertainty)`. Its confidence
level is not silently asserted to be Gaussian one sigma.

## Vujnović et al. 2002

Source: A&A 388 704–711, DOI **10.1051/0004-6361:20020560**, CDS J/A+A/388/704,
Tables 2–5 and their byte-layout ReadMe; PDF `aa7151.pdf`, sections 3–4.

All 106 raw rows are parsed independently. A-value limits, uncertainty limits,
uncertainty parentheses, approximate intensity/branching ratios, and the Table 5
theory marker survive. Other-author A values and reference strings are retained.
Table 3/4 intensity-only rows are not absolute gf measurements. Table 5's starred
3900.675 value comes from Tayal & Hibbert theory, not this laboratory measurement.

The six original UV/blue Al I promotions are Table 2 rows 1–4, 9 and 11. Their
lower levels are independently constrained by the ground-term energies 0 and
112.061 cm-1. Species, wavelength and lower energy all enter the join. Their
source upper-J identities enter the A-to-gf conversion. The source describes its
errors as summed relative errors of branching ratio and lifetime, so the ledger
retains an error bound rather than inventing a statistical confidence level.

Rows 5/6 (5s–4p) use theoretical fine-structure intensity ratios. Row 10,
3092.839 A, was not observed: its transition probability follows an LS ratio and
has no independent uncertainty. It cannot inherit row 11's 2% uncertainty.
Rows 7/8 (13150/13123 A) use indirectly measured branching loss plus lifetimes;
the paper does **not** say those two use the same assumed LS split as rows 5/6.
This refines RYA-1141's `THEORETICAL_LS_RATIO` label on rows 7/8; Burheim remains
the selected source for those transitions, and both older values remain evidence.

Direct inspection of Buurman & Dönszelmann 1990, A&A 227, 289–290, Table 1
shows 60.5(9) ns for 4p J=3/2 and 65(2) ns for J=1/2. The Vujnović CDS
ReadMe transposes those J labels. Burheim Table 1 agrees with the primary source.
The raw CDS file is preserved; its lifetime labels must not override the primary
identity. Buurman et al. 1986, A&A 164, 224–226, Table 1 independently verifies
the shared 5p lifetime 275(8) ns. Both ADS scans were visually inspected and
preserved in Reference Documents with canonical bibliography records.

The old intake used air wavelengths in its A-to-gf conversion. This verification
uses vacuum wavelength and retains the old value separately, exposing the small
deterministic correction instead of changing the frozen input.

## Johnson / Träbert / NIST Al II 2669

Johnson, Smith & Parkinson 1986, ApJ 308 1013–1017, abstract and experimental
section (`1986ApJ...308.1013J`): direct ion-storage measurement of the 3s2 1S0 –
3s3p 3P1o intersystem line, A = 3330 ± 230 s-1 at **90% confidence**. The ledger
stores that confidence bound; it is not a one-sigma uncertainty.

Vujnović's CDS ReadMe Table 1 cites the later Träbert et al. 1999 lifetime of
302(2) microseconds. Table 6 cites both Träbert and Johnson for its rounded
3.3e3 s-1 value. These are a shared source chain, not three independent
measurements to average. The later lifetime is a stronger comparison, but the
article is unavailable for direct review and its primary uncertainty/confidence
analysis is therefore not independently adjudicated. It is not silently
substituted for Johnson. Per the RYA-1134 execution decision, Johnson's measured
value is retained with an explicit `JOHNSON1986_RETAINED_TRABERT1999_UNAVAILABLE`
limitation rather than blocking the atomic handoff.

## Evaluated and theoretical sources

The held NIST ASD snapshot supplies `gi * fik`, source codes and per-line accuracy.
The u36/LS rows trace to the Opacity Project/Mendoza chain, as verified against
Kelleher & Podobedova 2008 Table 4 (`jpcrd372008911p.pdf`). An evaluated accuracy
grade is preserved as a separate field; it does not turn that calculation into
a laboratory measurement. Other source codes remain unresolved rather than
being inferred from the word NIST.

For unresolved fine-structure groups, a single component cannot replace the sum.
NIST grade E means greater than 50% error: it supplies no finite upper bound.
Neither 50% nor 100% is substituted. Such candidates remain visible and held.

The later RYA-1173 extraction of Nordlander & Lind Table A.1 supplies explicit
level identities and inherited source uncertainties beyond the narrow original
Al census matches. Its TOPbase entries may enter evaluated Reference membership;
they never gain primary-laboratory status or count as independent NIST evidence.
Its exact six-line Solar selection remains the distinct `asplund-al` replication
set, with the seventh telluric-rejected line retained.

## Component proof and its limits

The six preserved raw VALD deliveries retain lower/upper term, J and energy fields
that the merged schema omitted. The component ledger recovers these by exact
wavelength, lower energy and gf correspondence and retains source file and line
number. A group containing multiple fine-structure transitions remains explicit;
it cannot receive a single transition's total gf. For unique physical groups,
the raw fractions sum to unity and the proposed rescaling closes on the selected
source total. Source J must agree for every component; independently tabulated
upper energies must agree within printed precision (0.0006 eV).

For the frozen NIST snapshot, lower J comes from gi; upper statistical weight is
reconstructed from A, f and vacuum wavelength and accepted only within 3% of an
integer, allowing for printed rounding. The raw VALD J is checked independently.
This reconstruction is an explicit limitation of the snapshot, not a new primary
measurement. Atomic membership and permission to synthesize an exact holding
remain separate decisions. No production component deck is installed here.
