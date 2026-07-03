# RYA-506 — Procyon Fe theo-EW collapse: root cause (diagnose-only)

Branch `ryandamienschmitt/rya-506-procyon-fe-theoew-repro`. Analysis-only; no production
edits, no line-list / STAR_PARAMS edits, nothing banked. Repro harness:
`scratch/rya506_theoew_repro.py`.

## Verdict
**ROOT CAUSE: a silent multiprocessing failure in iSpec's theoretical-EW step returns an
all-zero theoretical_ew array, which the pipeline trusts and caches — so 100% of the Procyon
Fe pool fails the `theo < 5 mÅ` quarantine → empty EW pool → "MOOG: No abundances".** It is a
timing-sensitive race, not a physics or data problem: the measured EWs, stellar params, and
interpolated model atmosphere are all fine, and the Fe abundance IS recoverable when the race
does not fire.

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

## Proposed fix (STOP here per autonomy — non-trivial, bring to Ryan)
Do NOT loosen the theo<5 quarantine (validate-don't-tune — a wrong physics tolerance is a
separate ticket). Instead:
1. **Guard (pipeline side, safe):** in `_fe2_theoretical_ew`, if the returned theoretical_ew is
   all-≈0 for lines whose observed EW is well above the floor, **RAISE with diagnostics and do
   NOT cache** — never let a failed synthesis masquerade as truth. This refuses a bad synthesis;
   it does not admit lines.
2. **Reliability:** detect the zeros-return (dead child) and **retry** the synthesis (bounded),
   or run the theo-EW synthesis in-process/serially to remove the child-Process race.
3. **Robustness:** ensure `tmp_dir` exists before the iSpec calls.
4. **Upstream note:** iSpec's dead-child path should raise, not silently return zeros; guard on
   our side regardless.
Then re-run clean-from-raw and record the reproduced A(Fe) provenance for the RYA-322 ξ-pin /
reduced-EW-slope review.
