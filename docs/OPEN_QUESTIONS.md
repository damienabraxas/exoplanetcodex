# Open questions / upstream issues

Running log of open questions and upstream (dependency) bugs to report or track. New file
(RYA-506); add entries as they surface.

---

## iSpec: theoretical-EW / synthesis child silently returns zero-init on child death — RYA-506

**Where:** `ispec/synth/spectrum.py`, `calculate_theoretical_ew_and_depth` (and the sibling
`generate_spectrum`). The SPECTRUM work runs in a child `multiprocessing.Process`; the parent
initialises `output_ew = np.zeros(num_lines)` and polls a `JoinableQueue` with
`while p.is_alive() and num_seconds < timeout:`. **If the child dies without enqueueing a
result, the loop exits with `output_ew` still all zeros and iSpec returns them with NO
exception and NO log** — only the separate *timeout* branch logs; the dead-child branch is
silent.

**Impact for us (RYA-506):** on macOS the Python 3.8+ default start method is `spawn`, under
which the iSpec child (which imports the C `synthesizer` extension and expects fork semantics)
dies on spawn → silent all-zero theoretical EWs → 100% of the Procyon Fe pool failed the
`theo < 5 mÅ` quarantine → "MOOG: No abundances." A defensible calibration number
(A(Fe I; NLTE) ≈ 7.571) was made to look unreproducible by an unreported subprocess death.

**Our mitigation (already applied):** force `multiprocessing.set_start_method('fork')` at
pipeline import (`pipeline/abundances_derive.py`), plus a loud all-zero-batch guard in
`_fe2_theoretical_ew` that raises instead of caching/quarantining on a failed synthesis.

**Upstream ask (report to the iSpec project):** the dead-child branch of the poll loop should
**raise** (or at minimum log LOUDLY) rather than return the zero-initialised array — a silent
all-zero synthesis result is indistinguishable from "every line is a null line" and corrupts
any downstream quarantine/abundance step. Ideally iSpec should also let the caller select the
start method / use a `fork` context explicitly for the synthesizer child on macOS.
