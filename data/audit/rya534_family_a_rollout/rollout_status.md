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
| Si | — | — | — | ⏳ pending | — |
| Ca | 6122.217 / 6162.173 | −0.009 | +0.02 (our Ca_Mashonkina2017 +0.017) | ✅ PASS (small sign diff, both ~0) | 70054bdb… |
| Ti | — | — | — | ⏳ pending | — |
| Mn | — | — | — | ⏳ pending | — |
| Co | — | — | — | ⏳ pending | — |
| Ni | — | — | — | ⏳ pending | — |
| Sr | — | — | — | ⏳ pending | — |
| Ba | — | — | — | ⏳ pending | — |

**Landed: Na (RYA-533) + O = 2 / 11.** Provenance JSONs in `data/nlte_grids/gerber_ts/`;
gate outputs in this dir. Cross-engine (the RYA-525 model-atom-systematic diagnostic), Na:
INSPECT −0.107 / PySME −0.129 / TS-Gerber −0.068. O: Amarsi-2019 1D −0.134 / TS-Gerber −0.105.

### O note (validate-don't-tune in action)
First-pass anchor was a literature −0.20 (Drawin-era / 3D-inflated) → CHECK. Investigated: our own
banked Amarsi-2019 **1D-NLTE** leg gives −0.134 (the 3D leg −0.169). TS-Gerber −0.105 reproduces
the 1D value to 0.030 → PASS. The anchor was corrected to an independent published reference, not
fitted to the result.
