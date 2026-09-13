# The Al lineage, in the sources' own words

Every quotation below was read in the cited copy; nothing here is paraphrase. This file exists so
the lineage can be audited without re-fetching four papers, and so that the ONE place the published
record disagrees with itself is visible rather than smoothed away.

## 1. AGSS21 — the adopted solar value

Asplund, Amarsi & Grevesse 2021, A&A 653, A141 (DOI 10.1051/0004-6361/202140445),
Sect. 3 *Aluminium (Z = 13)*, read in `Reference documents/aa40445-21.pdf` pp. 9–10:

> Nordlander & Lind (2017) recently performed full 3D non-LTE spectral line formation calculations
> for Al i lines in the Sun, based on their new model atom that includes realistic cross-sections for
> inelastic collisions with neutral hydrogen (Belyaev 2013). They employed the same 3D hydrodynamical
> Stagger model solar atmosphere and 3D non-LTE radiative transfer code Balder as we employed in this
> study. Impressively, their solar modelling also agrees well with the observed IR emission lines of
> Al i. **They adopted the same lines and line data as in Scott et al. (2015b), except that they
> excluded the 1089.1 nm line due to telluric contamination**, finding a solar photospheric Al
> abundance of log ϵ_Al = 6.43 ± 0.03 (statistical + systematic). Here, we adopted the Nordlander &
> Lind (2017) solar abundance, which is identical to the 3D LTE results of Scott et al. (2015b); for
> the Sun, departures from LTE are in fact largely unimportant for **these five Al i lines**.

**AGSS21 PUBLISHES NO Al LINE LIST.** Its Table A.2 is Fe only (that is what RYA-1109 ingested). The
Al set has to be reconstructed from the papers AGSS21 points at — which is what this directory is.

`Scott et al. (2015b)` is disambiguated from AGSS21's own reference list (p. 26):

> Scott, P., Asplund, M., Grevesse, N., Bergemann, M., & Sauval, A. J. 2015a, A&A, 573, A26
> Scott, P., Grevesse, N., Asplund, M., et al. 2015b, A&A, 573, A25

so **2015b = A&A 573, A25 = "The intermediate mass elements Na to Ca"** (the paper that contains Al),
and 2015a is the iron-group paper. The `a`/`b` ordering is by author list, not by volume order, and
getting it backwards would send this census to the wrong paper.

## 2. Nordlander & Lind 2017 — the analysis AGSS21 adopts

Nordlander & Lind 2017, A&A 607, A75 (DOI 10.1051/0004-6361/201730427), Sect. 3.1.5
*Abundance analysis*, read in the published full text at aanda.org and in arXiv:1708.01949v2:

> Results of our disk-center solar abundance analysis are illustrated in Fig. 8, where the average
> abundances are given as unweighted arithmetic mean values with uncertainties representing the
> line-to-line dispersion. **Adopting a line selection and weights from Scott et al. (2015), but
> disregarding the line at 10 891 Å due to telluric contamination**, yields our recommended solar
> abundance, A(Al) = 6.43 ± 0.03, where the error takes into account systematic errors.

Fig. 8 plots the per-line abundances, and its x axis names the lines individually:

> 6696  6698  7835  8912  10768  10872

— **six lines**, and its caption repeats the weighting:

> Adopting the weights from Scott et al. (2015) and including estimates of systematic uncertainties,
> the final 3D NLTE abundance is A(Al) = 6.43 ± 0.03.

## 3. Scott et al. 2015b — where the line identities live

A&A 573, A25, Sect. 5.3 *Aluminium*:

> **We retained seven quite weak Al i lines (Table 2).** We note that Al is essentially all Al ii in
> the solar photosphere. […] Al i transition probabilities have been discussed by Kelleher &
> Podobedova (2008b). **The data for our adopted lines come from theoretical calculations by the OP
> (Mendoza et al. 1995), under the assumption of LS-coupling.**

and Sect. 2 on the weight column:

> We gave each line a weight from 1 to 3, depending on the estimated uncertainty on our measured
> equivalent width, which was further modified in some cases to account for uncertainties in the
> atomic data (Sect. 5).

with the direction of the scale fixed in Sect. 6:

> Trendlines give equal weight to each line (unlike our mean abundances, **where we give larger
> weights to higher quality lines**).

## 4. 🔴 THE ONE PLACE THE PUBLISHED RECORD DISAGREES WITH ITSELF

    Scott et al. 2015b   SEVEN lines retained          (Sect. 5.3 + Table 2, seven rows)
    Nordlander & Lind    SIX  lines in the abundance   (Sect. 3.1.5, minus 10891; Fig. 8 names six)
    AGSS21               "these five Al i lines"       (Sect. Aluminium, prose only)

7 − 1 = 6, and Fig. 8 of the primary names all six by wavelength. **AGSS21's "five" is not
reproducible from either source it cites.** This census therefore carries the set as SIX used + ONE
published exclusion, on the authority of the primary, and records the conflict rather than resolving
it silently — the count is a published fact about a published selection, and RYA-946 says a negative
selection must be preserved, not dropped.

Two candidate explanations were considered and NEITHER is adopted, because nothing in any of the
three papers states either one:

* AGSS21 miscounted (the simplest reading — nothing else in the sentence depends on the number);
* AGSS21 counted the 6696/6698 pair as a single doublet feature, giving five *features* over six
  *lines*. Fig. 8 of the primary plots them as two separate points with different abundances, so
  this would be AGSS21's own bookkeeping, not the primary's.

`published_line_count_conflict` in the provenance file carries this forward so no downstream reader
can inherit "five" as if it were sourced.

## 5. What the six lines are NOT

Nordlander & Lind analyse a great many more Al i transitions than the six — the resonance lines at
3944/3961 Å, centre-to-limb variation in 7835 Å, the 12.33 µm emission line, the HFS-sensitive IR
lines at 13123 and 16750 Å, HST/STIS UV lines in metal-poor stars, and J-band lines in HD 122563.
Table A.1 lists all 55 rows of line data used anywhere in the paper. **None of that makes them part
of the solar abundance set**, and RYA-946 is explicit that an abundance value is not evidence a line
was used. `nordlander_lind_2017_analysis_lines.csv` carries all 55 with a `role` column, and the
role is `ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE` wherever the paper does not say per line.
