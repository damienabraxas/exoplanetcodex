# RYA-1315 Ti line screen (Ca's RYA-1311 rule) and one reasoned quarantine

Ti had not been through the line screen that Ca got in RYA-1311 (34eabffe). This applies the same rule to every Ti band, on the 1D-LTE line artifacts of the live products: a line is flagged on a spectrum if robust |z| > 3 or red_chi2 > 5x the pool median, and listed if it is flagged on at least 2 spectra of its band (or on the only spectrum of a single-holding band).

- 69 lines are flagged (`flagged_lines.csv`). Several are clearly bad. Ti I 8734.71 reads 10-21 sigma low on all three spectra. Ti I 4888.71, Ti II 4529.48 and Ti II 5183.71 read about 11-13 sigma high.
- Every product's A is the median of its lines, so the screen barely moves anything (`impact_1d_lte.txt`). With all flagged lines out, every A moves by at most 0.021, all inside sigma. The one exception is Kitt Peak molecfit Ti I near-IR: -0.071 against sigma 0.176, from 10034.49 at 5.88. The HARPS VIS headline products move by at most +0.016.
- Decision: no Ti re-run to drop lines (Ryan: no extra runs; the effect is inside sigma). The list is kept here for the next Ti run. It is NOT entered in `problem_children.csv`, because that would make a re-run diverge from what is published without a decision.

Quarantined with reason (feed v1.343): Kitt Peak Kurucz-corrected Ti II near-IR, both tiers and all six treatments (12 rows). This follows the Ca near-UV precedent. The pool is 2 lines (9432.4, 9748.3) in the water band. Ti II 9432 reads 5.68 at red_chi2 66, and the product reads 5.30 +/- 0.38. On IAG and on molecfit-corrected Kitt Peak the same two lines read 4.81 and 4.84, and those products stand.

No ENGINE-A-3DNLTE for Ti: no per-line solar 3D-NLTE Ti table is published. Asplund 2021 adopts Scott+2015, which is 3D LTE with <3D> NLTE corrections for Ti I. A21 says outright that a consistent 3D non-LTE analysis is still owed. `data/reference/scott2015` holds LTE abundances only. Scott's own approximation, <3D> NLTE, is the synth-mean3D-NLTE-gerber-stagger treatment, and every Ti cell has it.

Empty cells, all declared gaps by line pool, not by engine: Ti I near-IR on IAG and Kitt Peak-Kurucz (17 lines, none with a lab gf); Ti I CRIRES Y (rya794: 4 lines, no lab gf; rya1054: 1 lab line below the depth gate); Ti II CRIRES Y/J/H/K (no lab-gf Ti II line, or no line in canonical_gf). No new Deep products (Ryan 2026-11-02).
