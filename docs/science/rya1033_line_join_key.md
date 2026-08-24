# RYA-1033 — a rounded wavelength is not a line identity

**Status:** fixed on branch `ryandamienschmitt/rya-1033-canonical-gf-join-key`.
**Evidence:** `data/audit/rya1033_rounded_key_join/` (regenerate with
`scripts/rya1033_rounded_key_audit.py`). **Guard:** `tests/test_line_match_rya1033.py`.

## The defect

Several joins tied a measured line to its atomic data through a 2-decimal-rounded air
wavelength. On the committed solar pool this drops **17 of 421** Fe I lines: each one has a
`canonical_gf` row within **1.17 mA**, and **none** is genuinely absent. The measured pool
stores Fe I 4787.49462 and `canonical_gf` stores 4787.495 — 0.38 mA apart, but they round to
4787.49 and 4787.50 and land on different keys.

An unmatched line did not fail loudly. It became `NaN` `gf_tier` with no `loggf_reference`
and travelled on as "ungraded", which is indistinguishable in every downstream product from
a genuine Kurucz-tier line — a wrong answer wearing the shape of a normal one (RYA-833).

## The sharper half: the key is not a function of the value

    round(6136.615, 2)      -> 6136.61     Python, correctly-rounded decimal
    np.round(6136.615, 2)   -> 6136.62     numpy/pandas, scale-multiply-round

`scripts/promote_solar_ew.py` rounded with pandas; `pipeline/abundances_derive.py` rounded
the same quantity with Python. One wavelength, two keys. On this pool the choice of library
alone moves the casualty count from **17** to **18**. A tolerance is a declared
approximation; a key that disagrees with itself is a coin flip, and no amount of care at the
call sites fixes it.

## Wavelength alone cannot identify an Fe line

`canonical_gf` holds **360** Fe I clusters whose members sit within 5 mA of each other, and
in **all 360** the members disagree on gf — they are different transitions that coincide in
wavelength. **7** of the 421 measured Fe I lines land on such a cluster:

    6065.48200 -> 6065.4820  EP 2.609  log gf -1.530  NIST-C+   <- the real line
                  6065.4850  EP 4.956  log gf -3.471  KURUCZ    <- 3 mA away, 1.9 dex off

The rounded key silently returned one of them. `pipeline/line_match` discards candidates
whose EP disagrees **before** choosing the nearest wavelength, and where no EP is available
it records the line as AMBIGUOUS and refuses. This is the rule `perline_product` and
`gf_grades` already applied after RYA-780/852; it is now shared rather than repeated.

## Why not join on `line_id` / `key_z`

The ticket's first option. It is not available: no measured artifact carries a stable line
id — every EW and per-line product is keyed `(element, ion, wavelength_air_A)`. Adopting one
is a schema migration across every committed and gitignored product, and it would still not
repair the pools that already exist. The declared-tolerance nearest match is the fix that
works on the data as it is, and the ticket names it as the sanctioned alternative.

The tolerance is **derived, not chosen**: the worst table-to-table disagreement over the
whole measured Fe I pool is under 2 mA, and 5 mA clears that ~4x while sitting far below the
~180 mA VIS Fe I line spacing. It is the value `anchor_pools` and `derive_band_products`
already use.

## The graded pool count

| keying | graded tiers | n |
|---|---|---|
| 2-dp rounded key (the defect) | LAB + NIST-C+ | 63 |
| tolerance match (the fix) | LAB + NIST-C+ | **64** |
| tolerance match (the fix) | LAB only | 23 |

The rounding defect costs **one** graded line, and the ticket's reconstructed **63** is
reproduced exactly: it is the LAB+NIST-C+ pool counted through the rounded key.

⚠️ **The 67 it was compared against is a different pool, not the same one keyed
differently.** 67 is the depth-gated LAB-only synthesis selection
(`FeI_4200_6910_*_SYNTH_GRADED_1D-LTE_lines.csv`): 176 in-band LAB lines minus the 109 above
the EW depth gate. Both numbers are correct for what they count, so the rounding bug is not
the 63-vs-67 gap — the two counts were never measuring the same thing.

## Found en route, NOT fixed here

**Fe I 6705.1169 is a Ruffoni-2014 LAB line published as ungraded.** The measured band_ew
artifacts carry synthesis-list wavelengths, which for some lines sit up to **24 mA** from
`canonical_gf` — well outside the 5 mA `_graded_mask` window — while the excitation
potential matches to **0.000 eV**, proving it is the same line. 22 of the 257 in-aggregate
HARPS VIS lines are unresolvable at 5 mA on wavelength alone; 20 resolve uniquely at 30 mA
once EP is required, and one of those 20 is graded.

This is a second, independent mechanism for the pool-count instability RYA-1033 describes,
and it is left alone deliberately: widening `_graded_mask` changes which lines are graded and
therefore moves published values, which is a re-measurement decision with its own control
(the same call RYA-959 made about `sigma_max`). `_graded_mask` also takes only wavelengths,
so doing it safely means threading EP through the product path. Flagged for Ryan.
