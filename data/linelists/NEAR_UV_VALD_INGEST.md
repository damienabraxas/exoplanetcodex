# Near-UV VALD re-ingest — what it found — RYA-713

Ryan: *"what about the VALD references in the IR? Not everything needs to be NIST graded"*
→ *"So hold on check UV then — did you reject lines that we have in VALD?"*
→ *"re-ingest the near-UV VALD with references."*

## Answer to the original question: no lines were rejected on provenance

The near-UV rejection was on **measurability**, never on gf. But the check was worth making,
because it exposed that we could not have judged provenance either way: every one of the
5135 near-UV Fe lines in `linelist_solar.csv` carried the bare tag `loggf_source = VALD3`,
and `canonical_gf` spans **3780.3–9199.9 Å** — it *starts above* the near-UV.

That is why the IR worked and the UV did not. The IR lines came through `canonical_gf`
carrying real source tags (BWL, GESHRL14, BK); the near-UV was never adjudicated at all.

## The information was already on disk

`vald_solar_nearuv_2000_3780_hfson_raw.txt` (68.5 MB) has carried per-line sources the whole
time, inline in VALD long format:

```
'_   Kurucz TiII 2016   1 wl:K16   1 gf:K16   1 K16 ...   Ti+'
                                      ^^^^^^
```

`scripts/ingest_vald_references.py` reads that tag and records it verbatim. It does **not**
assign accuracy — mapping source → uncertainty is a separate citable judgement kept in
`GF_SOURCE_ACCURACY`, so a wrong estimate can be corrected without re-parsing 68 MB.

## Fe I 3000–3780 Å — 4461 lines

| source | n | dex | |
|---|---|---|---|
| K14 / K13 | **3743** | 0.200 | Kurucz semi-empirical |
| **RU** | **591** | **0.030** | **Ruffoni et al. — FTS laboratory** |
| FBHM, KOH, BRW, MPV, JBL, KCN | 127 | — | unrecognised; **accuracy NOT assumed** |

## The finding that redirects the work: the laboratory lines are Fe II

**All 591 Ruffoni FTS lines are Fe II, not Fe I.** The near-UV Fe I set is essentially all
Kurucz semi-empirical. That also explains a number noted earlier and not followed up:
Engine A serves **828 Fe II** lines in 2960–4200 Å against only **203 Fe I**.

The near-UV is Fe II territory — and Fe II is the dominant ionisation stage in the
photosphere, so it is the better abundance indicator regardless.

**105 of the 106** near-UV Fe II lines at usable depth carry laboratory gf, and 12 are also
isolated (gap ≥ 0.30 Å).

## And they are still unmeasurable

| run | candidates | measured | usable | EW range |
|---|---|---|---|---|
| Fe I 3000–3800 | 901 | 901 | **0** | 107–922 mÅ |
| **Fe II 3000–3780** | **106** | **106** | **0** | **262–728 mÅ** |

**The near-UV is unmeasurable by profile fitting regardless of gf quality or ionisation
stage.** The blocker is the spectrum — median catalogued line gap 0.086 Å, continuum median
0.283–0.805, no unabsorbed wavelength in any window. That is now established twice, on two
ions, one of which has laboratory oscillator strengths.

## The larger find: this extract IS a near-UV line list

Synthesis was blocked because the GES synthesis list spans **4200–9200 Å** and nothing
covers below it. But the same file holds **55,849 records in 3000–3780 Å across 103
species** — Co I 15779, V I 4813, Fe I 4364, Mn I 3277, Nb II 2898, V II 2340 …

That is precisely what near-UV synthesis needs: every contributor in the window, not just
the target element. The blocker was never *"no near-UV line list exists"* — it is *"the
near-UV line list we hold is not in Turbospectrum format."*

Which is a conversion task with a known precedent (RYA-503 built the ExoMol→Turbospectrum
converter), not a research problem or an acquisition.

## Owed

1. **Convert the near-UV VALD extract to Turbospectrum format** — the concrete unblock.
2. **Identify FBHM / KOH / BRW / MPV / JBL / KCN** against the VALD source bibliography; 127
   lines currently carry no assignable accuracy and are neither trusted nor dismissed.
3. Re-ingest the *other* bands' VALD extracts the same way — the red-optical worked only
   because `canonical_gf` happened to cover it.
