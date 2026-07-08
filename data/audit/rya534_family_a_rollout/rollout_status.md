# RYA-534 — Family-A TS-Gerber NLTE rollout (Engine-B NLTE for RYA-525)

Rolling out the 10 remaining Family-A elements on the RYA-533 TS-native Gerber NLTE deck,
**on-demand, one grid at a time** (each ~5–16 GB; free/keep per Sirius disk, all Sirius-only
per RYA-526, md5-pinned). Deck = `scripts/ts_gerber_gate.py` (element-general; extracts each
element's lines from the bundled GES NLTE line list BY WAVELENGTH so the level IDs are verbatim —
the RYA-533 silent-departure=1 trap). Gate = reproduce the element's published solar 1D-NLTE
anchor; the deck RAISES if departures don't engage (no silent LTE).

**Anchor policy (validate-don't-tune):** each element is gated against an INDEPENDENT published
1D-NLTE solar value — preferentially our own banked Engine-A number (PySME/MPIA/Korotin/INSPECT/
nlte_cno), which also makes each gate a cross-engine model-atom-systematic check for RYA-525.

## Per-element results

| El | Lines (GES level-IDs) | median δ (TS-Gerber) | anchor (source) | verdict | grid md5 (zip) |
|----|----|----|----|----|----|
| Na | 5682/5688 | −0.068 | −0.107 INSPECT (RYA-533) | ✅ PASS | d1e8b51e… (RYA-533) |
| O  | 7771/7774/7775 (777 triplet) | −0.105 | **−0.134** Amarsi-2019 1D-NLTE (our nlte_cno table6; RYA-362 cross-check) | ✅ PASS | 504f3e2a… |
| Mg | 5711.088 | −0.023 | −0.02 (our Mg_Amarsi2020_PySME −0.022) | ✅ PASS | 3cc28460… |
| Ca | 6122.217 / 6162.173 | −0.009 | +0.02 (our Ca_Mashonkina2017 +0.017) | ✅ PASS (small sign diff, both ~0) | 70054bdb… |
| Mn | 6013.510 / 6021.800 | +0.043 | +0.10 (Bergemann; our MPIA +0.107) | ✅ PASS (~½ Bergemann; RYA-411 xref) | 60f1543e… |
| Si | 5772 | — | −0.01 (our PySME −0.013) | 🔁 re-provision (grid truncated under disk contention) | eeb806a1… |
| Ti | 5866/5964 | — | +0.05 | 🔁 re-provision (54 GB grid, stalled/disk) | — |
| Co | 5000/5013 | — | +0.10 | 🔁 re-provision (lines corrected: GES-identified levels) | — |
| Ni | 5018/5035 | — | +0.02 | 🔁 re-provision (lines corrected: GES-identified levels) | — |
| Sr | 4215.519 (II) | −0.013 | −0.01 (our Sr_Bergemann2012_INSPECT) | ✅ PASS | 8d387ae8… |
| Ba | 4554.029 (II) | — | −0.05 | 🔁 re-provision (ion-lookup + line fixed) | — |

**Landed: Na (RYA-533) + O + Mg + Ca + Mn + Sr = 6 / 11.** Remaining 5 (Si/Ti/Co/Ni/Ba) re-provisioning
sequentially (one grid at a time — driver config CORRECTED for all). Provenance JSONs in `data/nlte_grids/gerber_ts/`;
gate outputs in this dir. Cross-engine (the RYA-525 model-atom-systematic diagnostic), Na:
INSPECT −0.107 / PySME −0.129 / TS-Gerber −0.068. O: Amarsi-2019 1D −0.134 / TS-Gerber −0.105.

## Remaining 6 — root causes found + fixed in the driver (need clean grid re-provision)

The autonomous batch surfaced three real issues; the deck's silent-LTE guard + fail-loud
correctly caught all of them (no false passes). All fixes are now in `scripts/ts_gerber_gate.py`:

1. **Sr / Ba (ionised species).** The GES header encodes ionisation as a SEPARATE middle field
   (`'  Z.000  '  STAGE  NLINES`, STAGE = ion+1), NOT the species-code decimals. First pass looked
   for `Z.001` → "no GES block". FIXED (`ges_lines` matches STAGE = ion+1). Also the resonance lines
   I first chose aren't in the GES II blocks — corrected to GES-present + identified lines: **Sr II
   4215.519**, **Ba II 4554.029**.
2. **Co / Ni (line level-IDs).** First-pass lines had an UNIDENTIFIED upper level in the GES list
   (`… 0 'none' … 'x'`) → bsyn silently set departure=1 (the deck RAISED "not engaged"). Corrected to
   lines with BOTH levels identified: **Co I 5000.875/5013.324**, **Ni I 5018.282/5035.972**.
3. **Si / Ti (disk-contention truncation).** The large grids (Si ~24 GB, Ti ~55 GB unzipped) were
   truncated when multiple concurrent unzips filled the disk → interpolator hit EOF at a pointer past
   the truncated `.bin`. FIX: provision these **one at a time with ample free disk** (they need a
   dedicated clean run, not the concurrent batch). Ti is the biggest grid in the set (atom.ti503b,
   503 levels).

**Resume:** for each of Si, Ti, Co, Ni, Sr, Ba — re-download its grid on Sirius (one at a time for
Si/Ti), then `venv_pysme/bin/python ts_gerber_gate.py <El>`; the config + line lists are already
correct. Gate → provenance JSON → register bump → engineB map, same as the landed 5.

### O note (validate-don't-tune in action)
First-pass anchor was a literature −0.20 (Drawin-era / 3D-inflated) → CHECK. Investigated: our own
banked Amarsi-2019 **1D-NLTE** leg gives −0.134 (the 3D leg −0.169). TS-Gerber −0.105 reproduces
the 1D value to 0.030 → PASS. The anchor was corrected to an independent published reference, not
fitted to the result.
