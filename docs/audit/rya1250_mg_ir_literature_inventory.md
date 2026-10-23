# Mg IR literature inventory

This audit addresses the Mg cells held because the Asplund set does not cover
the full IR range. It is a literature inventory, not yet a production line-set
admission.

## Local evidence

The repository currently registers only `MG_AGSS21` for Mg. The local Mg
litscan names Asplund et al. (2021), Scott et al. (2015), and Lodders,
Bergemann & Palme (2025) as abundance references, but it contains no separate
Mg IR line set. The existing NLTE assets include `Mg_Bergemann_MPIA.csv`, but
that is a correction/grid asset, not a reference-line registry.

## Primary IR line source

Bergemann et al. (2017), ApJ 847, 15, Table 2, is the best immediate source
for a cited Mg I IR set. The paper uses Pehlivan Rhodin et al. (2017) gf values,
reports line-level gf errors, and gives the damping inputs used in its
1D/<3D>-NLTE analysis. Its solar line table includes:

* 8806.756 Å;
* the fine-structure groups at 10811.05–10811.16 Å and 10965.386–10965.450 Å;
* 11828.185 Å and 12083.278–12083.662 Å;
* 15024.992 Å, 15748.988 Å, 15765.645–15765.842 Å, and 15886.183–15886.261 Å.

Pehlivan Rhodin et al. (2017), A&A 598, A102, is the underlying atomic-data
source. It reports 34 experimental Mg I oscillator strengths, theoretical
complements, and uncertainties as small as about 5%, with particular emphasis
on the near-IR.

Civiš et al. (2013), A&A 554, A24, is a second laboratory IR source covering
Mg I transitions in the 1300–4500 cm⁻¹ range (approximately 2.2–7.7 μm). It is
useful for validating longer-wavelength transitions, but it does not by itself
provide coverage for the whole 0.9–1.7 μm production range.

## Admission constraints

The Bergemann table must be converted into a registered `set:MG_*` with
component-level wavelengths, excitation energies, gf provenance, uncertainties,
and blend/damping metadata. Fine-structure groups must remain component sets;
they must not be collapsed to one VALD component or assigned a made-up total
uncertainty. Solar usability also needs a blend check: Bergemann explicitly
excluded 15748 Å from its solar mean because blends contribute about half of
the equivalent width.

No independent, production-ready Mg II IR reference set was found in the
reviewed material. The current Mg II IR set remains the Asplund set; the
non-Asplund IR expansion is therefore Mg I first, with Mg II kept separate
until a cited line table and component treatment are verified.

## Near-UV and K-band grade audit

The orchestrator's active grade definitions were applied to the held Mg cells.
Codex Grade requires a published per-line gf uncertainty, a laboratory gf, or a
NIST accuracy class A--C. Deep Grade is a depth-selected subset of that same
gf-graded pool. The 3000--3780 Å Mg I selection contains 11 lines and the
CRIRES+ K selection contains 78 lines, but all are currently VALD3/Kurucz or
other ungraded provenance. Consequently both pools correctly return zero
eligible lines; the bands are data-present but grade-unresolved, not empty.

Civiš et al. (2013) reports measured Mg I K-band transitions, including lines
near 2386.6, 2385.0, 2384.3, 2281.8, and 2106.9 nm (their Table 3), and
provides QDT oscillator strengths in its later transition tables. The paper
does not provide the per-line gf uncertainty needed by the current Codex/Deep
definition, so those transitions are recorded as literature candidates and are
not emitted as a falsely graded product. A future admission requires recovering
an uncertainty-bearing gf source (or an applicable NIST A--C class) for the
same component identities.
