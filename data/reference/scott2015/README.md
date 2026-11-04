# Scott et al. 2015 I and II: per-line LTE abundances in five model atmospheres

Transcribed from the papers' line tables:
- Scott et al. 2015a, A&A 573, A25, Table 2 (Na to Ca). Local file `aa24109-14.pdf`.
- Scott et al. 2015b, A&A 573, A26, online Table 1 (Sc to Ni). Local file `Scott_2015_AA_573_A26_FePeak.pdf`.

Use `parse_scott2015.py` with the pypdf text dump of the table pages. There are 287 lines, and every λ row in both tables is matched.

Columns:
- `A_3D`, `A_mean3D`, `A_HM`, `A_MARCS`, `A_MISS`: the LTE abundance from the 3D, ⟨3D⟩, Holweger & Müller 1974, MARCS and MISS models.
- `wavelength_nm`: air wavelength, as printed in the paper.

Paper I's element headers lost the ion numeral in the PDF text. The first block of an element is the neutral (I) and the second is the singly ionised (II). This was checked by excitation energy: Mg II, Si II and Ca II all sit above 8 eV, or are the known Ca II lines.

What uses this table: the two atmosphere terms of Asplund, Amarsi & Grevesse 2021, Sect. 2.1, as set out in Scott 2015b Sect. 4:
- temperature structure = ½ the mean (⟨3D⟩ − HM)
- inhomogeneities = ½ the mean (3D − ⟨3D⟩)

They are averaged over the product's pool lines that appear in this table. RYA-1315, Ryan 2026-11-03.
