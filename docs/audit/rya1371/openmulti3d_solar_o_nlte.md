# OpenMULTI3D 1D NLTE corrections for solar oxygen (RYA-1371, Phase 0)

2026-11-18. First OpenMULTI3D product for an element the ORCH has already run. It follows
Ryan's ruling: OpenMULTI3D is a publication-pipeline engine wherever it has a model.

## What was done

1. **Engine:** OpenMULTI3D pinned at `32085d26`, 1D, run on the Mac.
2. **Atmosphere:** the MARCS solar model shipped with the engine,
   `p5777_g+4.4_..._t01` (Teff 5777 K, log g 4.4, [Fe/H] 0, ξ 1 km/s).
3. **Atoms:** three Keeper O atoms (see `keeper_atom_inventory.md`).
   - `o41f_ref`: the atom we already held.
   - `o41f`: the 2026-11-17 Keeper update, "Barklem + Kaulakys, closed by Drawin".
   - `o41j`: "Belyaev, Kaulakys, closed by Drawin".

   All three carry triplet log gf 0.369 / 0.223 / 0.002, the same values as the ORCH line
   list. The output pass uses the atom's own transitions (`output_from_line_list=F`). The
   engine's bundled line list gives 0.350 for 7772 and is not used for O.
4. **Grid:** 13 abundances, A(O) = 8.50–9.10 in steps of 0.05. Each run converges to a
   relative change below 1e-4 in 10 iterations and takes about 1.5 s.
   - Each run gives the NLTE, LTE and continuum **flux** (4-point Lobatto disk integration).
   - Equivalent width = ∫(1 − F/F_c) over a fixed window per line, identical for LTE and
     NLTE. Integration windows: 7772: 7768.5–7773.055 Å; 7774: 7773.055–7774.777 Å;
     7775: 7774.777–7778.5 Å; 8446 (one feature): 8443.5–8449.5 Å.
5. **NLTE correction per line:** Δ = A_NLTE(EW\*) − A_LTE, where EW\* is the engine's LTE
   equivalent width at the ORCH 1D LTE abundance for that line. This is the same differential
   definition Engine A uses (Amarsi+2019 grid). The ORCH LTE abundances are the SET-AGSS21
   1D-LTE products from the RYA-1230 nominal legs.

Driver: `scripts/run_openmulti3d_cog.py`, spec `config/experiments/rya1371_openmulti3d_solar_o_cog.json`.
Analysis: `scripts/analyze_openmulti3d_nlte_corrections.py`. Both are written per element:
atom, abundance grid and line windows come from the spec. Compact outputs are in
`data/results/rya1371/openmulti3d_solar_o_v1/`.

## Result (IAG; Kitt Peak agrees within 0.005 dex)

| Line (Å) | ORCH 1D LTE | Engine A Δ | OM Δ `o41f` | OM Δ `o41j` | OM 1D NLTE `o41f` | OM 1D NLTE `o41j` |
|---|---|---|---|---|---|---|
| 7771.944 | 8.957 | −0.156 | −0.147 | −0.215 | 8.810 | 8.742 |
| 7774.166 | 8.927 | −0.154 | −0.135 | −0.200 | 8.792 | 8.727 |
| 7775.388 | 8.918 | −0.153 | −0.111 | −0.170 | 8.807 | 8.749 |
| 8446 (3 comps) | 8.70–8.72 | −0.107 | −0.097 | −0.093 | 8.60–8.62 | 8.60–8.63 |

**Triplet means, by holding:**

| Holding | 1D LTE | Engine A | OM `o41f` | OM `o41j` |
|---|---|---|---|---|
| IAG | 8.934 | 8.780 | **8.803** | **8.739** |
| KP Kurucz2005 | 8.936 | 8.781 | 8.805 | 8.741 |
| KP molecfit | 8.938 | 8.784 | 8.807 | 8.743 |

**Asplund σ for the triplet (`o41f`).** The terms are:
- statistical error (SEM of 3 lines): 0.005
- NLTE term, ½|Δ|: 0.065
- 3D term, ½(Engine A 3D NLTE − 1D NLTE): 0.008

Combined, σ = **0.066**. For `o41j`, the NLTE term is 0.097.

## What it shows

1. **The 2026-11-17 Keeper change does not affect the triplet.** `o41f_ref` and the new
   `o41f` give identical corrections to 4 decimals. That update changed one [O I] 6300
   damping value and the header.
2. **H collisions set the 1D NLTE answer.** Belyaev (`o41j`) deepens the triplet correction
   by 0.06–0.07 dex relative to Barklem+Kaulakys (`o41f`). That is the size of the spread
   in the 3D NLTE literature: A21 8.69 vs Steffen 2015 8.76, Bergemann 2021 8.75 and Magg
   2022 8.77, which Amarsi & Grevesse 2026 attribute to the H-collision treatment. Engine A
   (Amarsi+2019) falls between the two atoms. The 8446 feature is insensitive to the
   choice (0.004 dex).
3. **The 777–8446 disagreement is already present in LTE.**
   - The ORCH 1D LTE triplet sits about 0.21 dex above 8446 (8.93 vs 8.71).
   - Every NLTE treatment narrows the gap only a little, to about 0.18 dex.
   - The fits have reduced χ² of 39–68 on every line.
   - So the discrepancy comes from the LTE fits (continuum or profile shape), not from NLTE.
4. **1D NLTE triplet 8.74–8.80 lies above the 3D NLTE literature (8.69–8.77).** This is the
   expected direction for 1D. Engine A's 3D NLTE runs 0.016 dex below its 1D NLTE.

## Limits

- This is a 1D NLTE correction product applied to the ORCH's existing LTE fits. It is not
  a new fit of the observed spectrum.
- The engine LTE equivalent widths at the ORCH abundance are 77 / 65 / 52 mÅ. The ORCH
  route stores no observed EW for this band ("no EW exists in this band"), so a direct
  observed-EW inversion and a profile fit with the engine are the next check.
- The atmosphere is MARCS 5777/4.4, while the ORCH LTE uses its interpolated solar
  atmosphere. Δ is differential, and the Bride 2x2 milestone found this mismatch moves
  areas by 1–3 mÅ, so it is second order for Δ.
- The 7774 window is only 1.7 Å wide. Wing truncation is identical for LTE and NLTE, so
  it cancels in Δ to first order.
- No 3D. A solar 3D OpenMULTI3D run is Phase 2.
