# `atomic_lines.tsv` — the near-UV synthesis line list (GENERATED, not tracked)

RYA-759 Move 2. iSpec-format atomic linelist for **3000–3780 Å**, built from our own
VALD holdings so the near-UV can be synthesised on the production path.

    ISPEC_DIR=/mnt/codex-data/engines/ispec_src \
      venv312/bin/python scripts/rya759_nearuv_synth.py --step linelist

~10 s, 12 MB. Deliberately **not committed**: it is fully regenerable from inputs that
*are* tracked (`data/linelists/vald_solar_nearuv_2000_3780_hfson_raw.txt` +
`pipeline/nearuv_linelist.py`), and `write()` re-reads it through
`ispec.read_atomic_linelist` on every build, so a corrupt copy cannot survive silently.

## What the build produced (2026-08-10, Sirius)

| | |
| -- | -- |
| lines / species | **55,798 / 102**, 3000.00–3780.00 Å |
| VALD records parsed | 161,526, **0 malformed** |
| vdW = 0 (VALD gave none; TS default applies) | 13,029 |
| γ_rad = 0 (VALD gave none) | 13,065 |
| per-line VALD gf source tags carried | 4,462 (Fe) |

## Excluded, on purpose, and counted

**Molecules — 6,223 lines** (NH 2271, OH 2094, CN 1205, CH 653). Turbospectrum's
molecular row needs an isotopologue code; VALD's extract does not say which
isotopologue a line belongs to, and guessing is the RYA-684 failure shape. Note the
consequence: iSpec's vendored molecular lists start at **400 nm**, so this band has no
molecular opacity from either source.

**H I — 51 lines, 3652.06–3770.66 Å.** With them present, bsyn prints `wrong H line
data file!` and aborts EVERY window in 3700–3780 Å; removing H I alone fixes it
(one-variable test). All 51 are high-order Balmer members crowding the 3646 Å limit,
which Stehle's tables do not cover. Hydrogen is not generally excluded — Hα/Hβ/Hγ come
through the optical GES list and synthesise normally. **Consequence:** between the
Balmer limit and ~3771 Å the merging Balmer series is a real opacity source this
synthesis does not reproduce.
