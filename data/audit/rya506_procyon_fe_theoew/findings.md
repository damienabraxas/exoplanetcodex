# RYA-506 — Procyon Fe theo-EW collapse: root cause + FIX

Branch `ryandamienschmitt/rya-506-procyon-fe-theoew-repro`. Repro harness:
`scratch/rya506_theoew_repro.py`. Fix APPLIED to `pipeline/abundances_derive.py` (no
STAR_PARAMS / line-list edits, no filter loosening, no merge).

## Verdict — DEFINITIVE cause is the macOS `spawn` multiprocessing default
**ROOT CAUSE: iSpec's theoretical-EW step runs its SPECTRUM synthesis in a child process
that imports the C `synthesizer` extension and relies on FORK semantics. On macOS, Python
3.8+ defaults the multiprocessing start method to `spawn`; under `spawn` that child dies
(re-importing the C extension) WITHOUT enqueueing a result, and iSpec silently returns a
zero-initialised theoretical_ew array. The pipeline trusts + caches the zeros → 100% of the
Fe pool fails the `theo < 5 mÅ` quarantine → "MOOG: No abundances".** Not physics: the EWs,
params, and atmosphere are all fine; the Fe abundance reproduces once the child survives.

Causality was established in the mandated order (RYA-506 refinement 1): a fresh-EMPTY
`/tmp/ispec_codex` still fails → the missing/stale tmp dir is NOT the trigger. Forcing the
start method to `fork` then fixes it **deterministically** (3/3 clean, byte-identical Fe
pool + A(Fe)) — so the trigger is `spawn`, and the fix is elimination (restore `fork`), not
retry. (An earlier framing called it a parent/child "result race"; the deterministic
spawn-fails / fork-succeeds result supersedes that — it is the start method, not timing.)

## Evidence chain
1. Clean-from-raw stage is fine: `procyon_normalized.csv` (SNR 2272, in-gate) and
   `procyon_ew.csv` (1820 Fe EWs, 5.3–299.5 mÅ, median 116.7) — sane.
2. The **plain** `ad.run('procyon', engine='spectrum')` fails **deterministically** (5/5):
   `Fe I triage: 0 clean, 96 dropped/quarantined (theo<5 mÅ)`, `EW sanity filter: 106 → 0`,
   `MOOG: No abundances`.
3. The **same run under a transparent logging wrapper** on `_fe2_theoretical_ew` (which only
   calls the original and prints) **succeeds reliably** (3/3): theo-EW sane (Fe I 6.25–144.8
   mÅ, 0 dropped; 106 → 77), and it converges to **A(Fe I; NLTE) = 7.571, A(Fe II) = 7.535**.
   Same inputs, opposite outcome → the theo-EW output is unstable w.r.t. execution timing.
4. **A(Fe I; NLTE) = 7.571 vs RYA-322 stored 7.593** (Δ 0.022 dex) → the number reproduces
   from raw once the theo-EW succeeds. Confirms the ticket framing: the value is likely
   correct; what was undefensible is *how it was produced*.

## Mechanism (the smoking gun)
`_fe2_theoretical_ew` (`pipeline/abundances_derive.py:305`) →
`ispec.calculate_theoretical_ew_and_depth` (`ispec/synth/spectrum.py`):
- runs the SPECTRUM synthesis in a **child `multiprocessing.Process`** with a `JoinableQueue`;
- **initialises `output_ew = np.zeros(num_lines)`**;
- `while p.is_alive() and num_seconds < timeout:` polls the queue. **If the child dies without
  enqueueing a result, the loop exits with `output_ew` still all zeros** — and only the
  *timeout* branch logs; the **dead-child branch returns the zeros with NO error and NO log**.
- `_fe2_theoretical_ew` accepts these zeros and **caches** them (`_fe2_theo_cache`), so every
  later convergence iteration reuses the poisoned ~0 → permanent 100% quarantine.

The logging wrapper perturbs parent-process scheduling enough that the child reliably finishes
before the poll loop gives up — hence instrumented=success, plain=failure.

## Loud-fail-rule violations found
1. iSpec silently returns zero-initialised `output_ew` when the synthesis child dies (no raise).
2. The pipeline's theo-EW wrapper does not validate that a pool of known-strong Fe lines
   produced non-trivial theoretical EWs — an all-~0 return for solar-strength Fe at 6554 K is
   physically impossible, yet accepted as truth and cached.
3. Compounding: `write_atmosphere`/`calculate_theoretical_ew_and_depth` write into
   `tmp_dir='/tmp/ispec_codex'` **without ensuring the directory exists** → `FileNotFoundError`
   if the dir is absent (surfaced when the stale dir was cleared).

## Persistent-state dependency (the reproducibility defect)
RYA-322's A(Fe)=7.593 rode on a processed Fe pool that had already passed this step on a
machine where the child happened to survive; a clean checkout + from-raw run hits the race and
gets zero lines. No un-versioned artifact is *needed* once the theo-EW step is made reliable —
the pool regenerates correctly when the child completes.

## Applied fix (Ryan-approved approach; does NOT loosen the theo<5 filter)
All in `pipeline/abundances_derive.py`:
1. **Eliminate the root cause:** at module import, force the multiprocessing start method to
   **`fork`** (`_mp.set_start_method('fork', force=True)`), restoring the fork semantics the
   iSpec SPECTRUM child relies on. Verified: 3/3 clean, deterministic, byte-identical result.
   This is "eliminate > mitigate" — no retry (retry was ineffective: under `spawn` the child
   dies every attempt; under `fork` it never dies).
2. **Loud-fail guard (backstop):** `_fe2_theoretical_ew` now RAISES and does NOT cache if the
   synthesis returns an all-zero **batch** (`np.any(theo > 0)` is False) — the clean signature
   of a failed child. Guards on the batch signature, not a per-line "above-floor" threshold, so
   a single weak line legitimately near zero (inside a batch with non-zero neighbours) is still
   quarantined normally; only a wholesale synthesis failure trips it.
3. **Robustness:** `os.makedirs('/tmp/ispec_codex', exist_ok=True)` before the iSpec call, so a
   missing tmp dir can never crash `write_atmosphere`.
4. **Upstream (logged in `docs/OPEN_QUESTIONS.md`):** iSpec's dead-child path should raise, not
   silently return zero-init — report to the iSpec project.

## Reproduced provenance (for the RYA-322 ξ-pin / reduced-EW-slope review)
Clean-from-raw, post-fix, deterministic across runs: **A(Fe I) = 7.557 (1D-LTE) / 7.571 (NLTE,
MPIA-Bergemann), A(Fe II) = 7.534 / 7.535**, n_lines 50 / 9. A(Fe I; NLTE) 7.571 vs the stored
RYA-322 7.593 → Δ 0.022 dex; Fe I − Fe II (NLTE) ionization gap ≈ 0.036 dex. The number is now
reproducible from raw; RYA-322's ξ-pin + reduced-EW-slope should be re-run on this pool.

## Blast radius
`calculate_theoretical_ew_and_depth` is the universal theo-EW gate (every star + element,
incl. the solar anchor). The `fork` fix + guard apply pipeline-wide. The solar-anchor
clean-from-raw reproduction is tracked separately as **RYA-509** (runs after this lands).
