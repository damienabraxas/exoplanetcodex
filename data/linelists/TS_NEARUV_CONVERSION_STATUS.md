# Near-UV VALD → Turbospectrum conversion — RYA-713

Ryan: *"convert the near-UV VALD extract to Turbospectrum format."*

## Delivered

`scripts/vald_to_turbospectrum.py` → `data/linelists/ts_nearuv_3000_3780.lte`

| | |
|---|---|
| transitions | **55,849** |
| species blocks | **103** |
| span | 3000–3780 Å |
| size | 7.2 MB |
| parses with the project's own header regex (`ts_gerber_gate`) | **103 blocks, 55,849 rows, 0 structural problems** |

Top species: Co I 15779 · V I 4813 · Fe I 4364 · Mn I 3277 · Nb II 2898 · V II 2340 ·
Ce II 2085 · Cr I 1703.

This is what near-UV synthesis needs — every contributor in the window, not just a target
element. **The blocker was never acquisition; the file was on disk. It was format.**

## Status: STRUCTURALLY COMPLETE, NOT ENGINE-VALIDATED

`bsyn` accepts the file, reports it (`1 line lists: …ts_nearuv_3000_3780.lte`), runs without
error — **and emits an empty spectrum.**

**But so does the reference GES list.** Running the identical harness against
`nlte_ges_linelist_jmg17feb2022_I_II`, a known-good file, produces an empty spectrum too.

**So the fault is in the test harness, not the converted list.** The conversion is unproven
either way — it may be correct.

## What I got wrong along the way, and it is the same mistake twice in one session

I made three format corrections attributing the empty spectrum to my output, **before ever
running the reference list through the same harness**:

1. **Header alignment** — the species code is LEFT-aligned in the 20-char field
   (`'  26.000            '`), and I had right-aligned it.
2. **LTE vs NLTE row layout** — NLTE rows carry a `gamma_stark` field between `gamma_rad`
   and the parities plus trailing level-index/label/flag fields; LTE rows do not. I copied
   the NLTE layout.
3. **`MARCS-FILE: .false.`** and the **`DATA/` symlink** in the working directory — both
   real setup requirements, both taken from the working gate.

Items 1–3 are almost certainly genuine improvements: matching the reference byte-for-byte is
right regardless. But I *attributed* the failure to them without evidence, and the control
that would have redirected me took one command.

This is the same shape as the synthesis control earlier: eleven runs debugging a
reimplementation before checking whether the thing it was meant to reproduce still worked.
**Run the control on a known-good input before diagnosing your own output.**

## Owed

1. **Fix the load-test harness.** It is a bsyn input-deck problem, not a line-list problem —
   the reference file proves that. Likely candidates: the `'NLTEINFOFILE:' ''` empty value,
   or a missing key the gate supplies that my deck omits.
2. **Then** re-run the near-UV list through it; the conversion may already be correct.
3. **Two data caveats to resolve before science use**, both currently carried honestly:
   * **13,029 of 55,849 lines have `vdW = 0`** — VALD supplied no van der Waals damping.
     Turbospectrum will apply its own default; that should be a stated choice, not a silent one.
   * **51 lines carry positive vdW values** (49.0, 12.0, 11.0). VALD encodes ABO cross-sections
     as packed positives, but these are small integers and do not look like ABO. Not guessed at.
   * **VALD's Stark parameter is dropped** (the LTE row has no such field). Recorded here
     rather than forced into a column whose convention is unconfirmed.
