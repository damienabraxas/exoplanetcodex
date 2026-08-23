# RYA-1006 — what the conditioning-axis collision actually overwrote

The **summary** products displaced from `data/results/band_products/` on 2026-08-23, kept
here under names that say what they are. The committed anchor names were restored from git
and read 7.461 / 7.417 / 7.535 again.

⚠️ **The matching `_1D-LTE_lines.csv` files are NOT here.** They are gitignored (RYA-469),
they existed in no other worktree, backup or tmp copy (`find` over the whole box), and the
conditioned versions that replaced them are kept **outside git** at
`data/audit/rya1006_preserved/` on the machine that ran them. The anchor's per-line data is
**gone and must be re-derived**; `pipeline.anchor_pools.load` loud-fails until it is.

| file here | what it is | A | n |
|---|---|---|---|
| `*_SYNTH_DEEPGRADED_LOCALRENORM_*` | RYA-1000's Kitt Peak deep run, `--local-renorm` | 7.337 | 108 |
| `*_SYNTH_DEEPGRADED_LOCALRENORM_*` (harps) | RYA-1000's HARPS deep run, `--local-renorm` | 7.339 | 109 |
| `*_SYNTH_ONDISK-n38-UNIDENTIFIED_*` (kpno) | a 07:59 shallow run, flags unrecorded | 7.28 | 38 (2 excl) |
| `*_SYNTH_ONDISK-UNIDENTIFIED_*` (harps) | its HARPS twin | 7.318 | 38 (2 excl) |
| `*_SYNTH_FROMEW_*` | the RYA-986 pools, **unconditioned**, copied before the live 08:07 renorm runs could land on them | — | 1483 |

🔴 **The provenance files are byte-identical to the ones they displaced.** `diff` is clean
against the moved-aside baselines in `/private/tmp/rya1000_baselines/`. That is the finding:
before this ticket, a conditioned product was indistinguishable from an unconditioned one by
filename AND by provenance — only the number moved. Both are fixed
(`derive_band_products.conditioning_tag` and `_conditioning_note`).

The n=38 pair is recorded as UNIDENTIFIED rather than guessed at: no run of that shape is in
any log I can reach, and a label invented here would be a provenance claim nobody measured.
