# RYA-1354 v70: beyond gf (lifetimes, radiative damping, van der Waals broadening) and the inverted-core-FS rule

Status: **CANDIDATE / HOLD** (RYA-1352). Declared on RYA-1354 (2026-11-20 18:31 UTC) before any run. No published data touched. Nothing here enters a synthesis.

Outputs:
- `data/audit/rya1354/mg1_model_v70/`
- `al2_model_v70/`
- `si3_model_v70/`

Scripts:
- `scripts/rya1354_twovalence_model.py` (rule C)
- `scripts/rya1354_lifetimes.py` (A)
- `scripts/rya1354_vdw.py` (B)

Photoionization is split out as v70b and will be declared separately.

## C. Inverted-core-FS rule (declared): PASS

An l whose core-ion NIST fine structure is inverted in the majority of levels stays unscaled (s_l = 1). The d shell is inverted in all three cores: Mg II 5/5, Al III 4/4, Si IV 2/2 levels. So s_d goes from 1.10 / 0.39 / 0.02 to 1.0.

| Ion | Blind intercomb. v69 → v70 (Δ vs lab) | |Δlog gf| p90 / p99 | Lab median / abs68 | Coverage 1σ / 2σ |
|---|---|---|---|---|
| Mg I 4571 | +0.036 → **+0.036** | 0.001 / 0.011 | +0.020 / 0.022 (n = 24) | 83 / 96 % |
| Al II 2669 | −0.009 → **−0.009** | 0.017 / 0.33 | −0.001 / 0.021 (246) | 90 / 97 % |
| Si III 1892 | +0.002 → **+0.002** | 0.057 / 0.82 | −0.004 / 0.022 (107) | 87 / 96 % |

The only lines that move are weak, d-involving, spin-changing lines (log gf ≲ −8). Those lines were where the arbitrary d scale had been acting.

Mg I variants were re-run under the rule. The same two compact-box variants are excluded by the v69 energy rule (rcomp35 0.87 %, rcomp20 2.13 %). Classes: C 1,359, D 96, E 2,762.

## A. Level lifetimes from our own A-values (declared): **PASS 6/6**

Method: τ = 1/ΣA over every model E1 decay. For ³D the comparison uses the g-weighted mean decay rate, because the lab fine structure is unresolved.

| Term | Lab τ (ns), Aldenius 2007 LIF | Model τ (ns) | Δ | z (quadrature / correlated) |
|---|---|---|---|---|
| 3s4s ³S | 11.5(1.0) | 9.12 ± 0.85 | −21 % | −1.82 / −1.46 |
| 3s5s ³S | 29(3) | 25.8 ± 1.9 | −11 % | −0.89 / −0.67 |
| 3s3d ³D | 5.9(4) | 5.66 ± 0.62 | −4 % | −0.32 / −0.29 |
| 3s4d ³D | 17.6(1.2) | 15.7 ± 1.5 | −11 % | −1.00 / −0.82 |
| 3s5d ³D | 33(3) | 33.1 ± 2.9 | 0 % | +0.02 / +0.02 |
| **3p² ³P** (doubly excited) | 2.0(2) | **1.89 ± 0.18** | −6 % | −0.42 / −0.36 |

**The b triplet has a lab tension.** For 3s4s ³S (a single decay route, the b triplet), our model's b-line gf, the Pehlivan Rhodin 2017 lab gf in our referee file, and the older lifetimes all point the same way:
- the older lifetimes are 9.7(6) ns (Kwiatkowski 1980) and 10.1(8) ns (Andersen 1972), both as tabulated by Aldenius;
- the Pehlivan Rhodin gf imply τ = 9.6 ns, and our model gives 9.1 ns.

Aldenius's 11.5 ns, and the gf derived from it (5183 Å −0.239, used by GES), sit about 0.07 dex lower. We report this; we don't adjudicate it.

**Γ_rad vs VALD/Kurucz** (61 matched lines, report only): median −0.026, abs68 0.16 dex.
- For 3p² ³P, VALD's 9.02 implies τ ≈ 1 ns. Our 8.72 (τ 1.89 ns) agrees with the measured 2.0 ns, so the model beats VALD/Kurucz where a lab lifetime exists.
- The largest differences are on autoionizing 3p nl upper levels.

Files: `levels_lifetimes.csv`, `lines_gamma_rad.csv`, `lifetimes_vs_aldenius2007.csv`, `gamma_rad_vs_vald.csv`, `lifetimes_summary.json`.

## B. van der Waals (H-collision) broadening from our wavefunctions (report only)

1. **ABO implementation control (declared ≤ 2 %): PASS.**
   - The tables are re-implemented from Turbospectrum `anstee.f` / `Anstee-1802-2014.dat` (Anstee & O'Mara 1995; Barklem & O'Mara 1997; Barklem, O'Mara & Ross 1998).
   - b triplet: σ = 733.9 a₀², α = 0.240, against GES 729 / 0.238 (0.7 %).
   - 5528 Å: 1465 / 0.312, against 1461 / 0.312 (0.3 %).
   - The σ, α → log γ₆ conversion at 10⁴ K gives −7.264, against GES −7.267.
   - 4703 Å is outside the tables (5d n* = 4.52 > 4.0); GES carries the Barklem et al. 2000 direct value.
   - Against VALD log γ₆ over all in-table lines (n = 28): median +0.0015, abs68 0.0018.
2. **n*_r from the CI ⟨r²⟩ (new).** The active electron's ⟨r²⟩ = ⟨Ψ|r₁²+r₂²|Ψ⟩ − ⟨r²⟩(Mg⁺ 3s).
   - **The hydrogenic assumption ABO relies on holds for singly-excited Mg I levels.**
     - |n*_r − n*| ≤ 0.05 from 3s3p ¹P upward.
     - Δlog γ₆ p50 0.003, p90 0.027, max 0.14 dex (n = 104).
   - It fails for the lowest levels: 3s² (Δn* +0.19) and 3s3p ³P (+0.21). There, subtracting the ion 3s ⟨r²⟩ is itself approximate.
     - Effect on 4571 Å: +0.14 dex.
     - Effect on the b triplet: +0.01 dex.
3. **Doubly-excited levels break the standard recipe.** 3p² ³P sits 5.36 by the energy formula, so standard practice treats it as a Rydberg state. Its ⟨r²⟩ is compact (n*_r ≈ 2.0).
   - Unsöld with the model ⟨r²⟩ is **~1 dex (10×) narrower** than Unsöld from the energy-based n* for the 2776–2783 Å lines.
   - 18 lines are affected (shift > 0.3 dex).
   - Elsewhere, Unsöld with the model ⟨r²⟩ equals the hydrogenic value (p90 0.009 dex, n = 2,826).
4. **Gap in our line lists.**
   - Of 508 Mg I lines matched in `linelist_solar.csv`, **470 carry vdW = 0**.
   - The GES/iSpec VIS list has zeros for 5711, 6318/6319, 7657 and others. The IR lists have no Mg vdW at all (v69).
   - 24 lines get an in-table ABO fill candidate (`vdw_fill_candidates.csv`), e.g. 15024/15040/15047 (−7.148) and 15740/15748/15765 (−6.923).
   - Most out-of-table lines (5711 has upper n* 3.47 > 3.0) need ABO beyond the tables. That is a v70c item.

Files: `vdw_levels.csv`, `vdw_lines.csv`, `vdw_vs_vald.csv`, `vdw_fill_candidates.csv`, `vdw_summary.json`.

## Sources

All of these are in Reference documents with rows in `data/refs/bibliography.csv`:
- Aldenius et al. 2007
- Anstee & O'Mara 1995
- Barklem & O'Mara 1997
- Barklem, O'Mara & Ross 1998
- Barklem, Anstee & O'Mara 1998 (PASA)
- Barklem 2016

Not held, with rows marked "NOT HELD": Barklem, Piskunov & O'Mara 2000 and Butler, Mendoza & Zeippen 1993.

## Tests

19 pass. New tests cover the ABO control, the hydrogenic ⟨r²⟩ inverse, the Unsöld scaling and the rule flag.

## Next (to be declared)

- **v70b:** photoionization from the box pseudo-continuum. Referee: TOPbase / Butler 1993 (needs the PDF).
- **v70c:**
  - a solar strong-line wing test of the vdW candidates (b triplet, 8806, 15024–15047), including ABO vs n*_r;
  - vdW for out-of-table lines;
  - the 3p² ³P lines (ours vs standard).
- **Line-list ticket:** Mg I vdW/Γ_rad zeros in the GES and IR lists; the b-triplet lab tension (Aldenius vs Pehlivan).
