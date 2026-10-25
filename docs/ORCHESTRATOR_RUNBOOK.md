# Orchestrator Runbook -- one solar element, A to Z

**Audience:** any agent (Claude, Codex, ...) or person taking a new element through the
pipeline. **This file owns the operating procedure and changes with the code.**
Start at `LEDGERS.md`: code establishes executable behavior but cannot override
ratified science. Report a conflict with both sources; stop the affected scientific
action until the owning authority resolves it. Follow `DEV_CYCLE.md` for authorization
and evidence handoff.

Worked examples: Si and Al (RYA-1233, science PR #583, site PRs #53/#54). Copy what they did.

---

## Architecture and invocation boundaries

The standalone front door is `run_pipeline.py`; no agent service is needed.
Resolve host paths through `config/path_register.yaml` and the committed host profile.
The source-level map below is a routing guide; inspect the current implementation and
issue comments before assuming an example still describes the available behavior.

| Layer | Owner | Contract |
|---|---|---|
| CLI | `run_pipeline.py` | `--element` selects the matrix; `--all-elements` selects the sweep; omitting both selects the announced legacy whole-star chain. |
| Prepare | `pipeline/element_prepare.py` | Reviews literature/gf, lines, culls, ion-aware NLTE and holding frames. Writes reports; `--apply` changes inputs and needs approval. |
| Sweep | `pipeline/run_sweep.py` | Derives the element axis, runs Fe first, records every element including Fe-gate omissions. It invokes the matrix, not the element publisher. |
| Identity and routing | `pipeline/run_descriptor.py`, `pipeline/band_policy.py` | Resolve holding/window/species/deck/route/pool into preconditions and stage commands; keep per-deck outputs separate. |
| Matrix execution | `pipeline/run_matrix.py` | Plans serially, dispatches permitted cells with bounded concurrency, records reasons, progress and fingerprints; unchanged recorded work can SKIP. |
| Budget and admission | `pipeline/element_publish.py` | Derives units from built synthesis cells, runs paired legs and assembly, then calls publication gates. Preserves historical artifacts and withdraws superseded synthesis products. |
| Numerical budget | `scripts/rya1230_cno_budget_legs.py`, `scripts/rya1230_cno_budget.py` | Reuse the scientific recipe; pair each response to its own pool and carry missing/invalid components as explicit holds. |
| Evidence | `pipeline/element_verdict.py`, `pipeline/coverage_report.py` | Literature diagnostics and per-cell publication/hold/staleness/work-owed evidence. A coverage DONE can include HELD; it is not scientific sign-off. |
| Human handoff | `scripts/post_run_report.py`, `DEV_CYCLE.md` | Preview reporting before an authorized Linear write; PR/review and website publication remain distinct actions. |

`--dry-run` prevents stage dispatch, but still writes progress/reports; use an isolated
fixture or the Python API's `report_dir` for a disposable destination. `--prepare`
reads scientific inputs and writes a preparation report; it is not filesystem read-only.
`--coverage` reports existing evidence. `--publish` is a mutating scientific action:
never append it to a discovery command. In the current CLI, `--all-elements` does not
forward `--publish` or `--jobs` to the sweep; use the reviewed single-element path for
budget/publication. This limitation is not permission to hand-publish around gates.

Inspect current RYA-1232/1233 and the relevant element test ticket (including RYA-1250
for Mg) chronologically. A merged development PR, a successful measurement, a SKIP,
a feed entry and a complete budget are different claims. Check the actual branch,
artifacts and treatment provenance; do not infer in-transfer 3D physics from a cited
solar correction table or a mean-3D label. The readiness audit at
`docs/audit/rya1274_codex_readiness.md` records the inspected baseline and open gaps.

---

## 0. Rules (read once, apply always)

| Rule | Why |
|---|---|
| **The orchestrator does the work; the agent starts runs and reports.** Change code only after the issue is reported to Ryan and agreed. | The goal is a pipeline that runs A to Z without an agent steering it. |
| **Data lives on Sirius.** Paths come from `config/path_register.yaml` (Sirius defaults). Never hardcode a path; never keep an env file outside the repo. | Another agent must reproduce your run. |
| **Raw spectra are never used.** Spectra are telluric-corrected and normalised ONCE per holding and then frozen. Do not re-prepare or re-gate them. | Settled for every element. |
| **Check the literature first** (Asplund, Amarsi & Grevesse 2021 = `aa40445-21.pdf` in the local reference folder; then the paper it cites for the element). We are repeating published work, not inventing it. | Weeks were lost choosing methods the literature already settled. |
| **Kill only your own processes, by PID.** Never `pkill -f`. Other sessions run on the same hosts. | A pattern kill destroyed another session's run. |
| **Never worktree under `~/actions-runner` on Sirius.** CI wipes `.git/worktrees`. | |
| **Ryan merges PRs.** Open them; do not merge. One element per PR. | |
| **Never commit files you did not create** (untracked `solar_*_latest.json` etc. belong to other work). | |
| **Linear comments ASCII only. Never log `LINEAR_API_KEY`.** | |
| **Website text is plain language** -- no ticket numbers, no internal IDs such as `SET-MG_AGSS21`. | |

---

## 1. Host check (1 minute)

```bash
# Sirius: nothing to source.            Mac (a mirror of Sirius):
source config/hosts/mac.env            # only on the Mac
python run_pipeline.py --check-paths   # must print 'ok' for every line, rc 0
```
Synthesis (iSpec/Turbospectrum) on Sirius uses `~/venv_rya1103/bin/python` (resolved
automatically); `venv_ci` cannot run iSpec. PySME NLTE work is Sirius-only.

## 2. Prepare -- inspect inputs and write the readiness report

```bash
python run_pipeline.py --star solar --element <El> --prepare
```
Writes `data/results/orchestrator/prepare/solar_<El>.{json,md}` and prints READY or
NOT READY with the blockers. Clear each blocker below, then re-run prepare until READY.

| Prepare check | Meaning | How to clear it |
|---|---|---|
| **A published sets: 0 lines** | No published solar line set is registered for the element. | Step 3. |
| **F literature: litscan MISSING** | No literature anchor file. | Step 4. |
| **D NLTE: N lines not in `<El>_Amarsi2020_PySME.csv`** | The NLTE table lacks graded lines. | Step 5 (Sirius). |
| **C cull candidates** | Saturated lines not used by the solar literature. | Step 6 (Ryan's yes/no). |
| **B gf source ... unclassified** | A gf source has no correlation entry. | Resolve the source under `data/reference/gf_error_model.csv` and its owning procedure; add only a supported, cited classification, never a guessed correlation. |
| **F2 literature need** | Bands with no error-stated gf: no Codex/Deep pool there. | Not blocking. Record it in the ticket; do not invent gf errors. |
| **E holdings not serving** | A spectrum does not cover a window at rest. | Report it; do not re-prepare spectra. |

## 3. Register the published line set (the Asplund recipe)

1. Find the element in Asplund+2021 Sect. 2 (Sect. 2.4 for Mg) and its line table (Table A.1 for Na/Mg/K/Ca,
   A.2 for Fe; otherwise the paper Asplund cites, e.g. Amarsi & Asplund 2017 for Si,
   Nordlander & Lind 2017 for Al). Note: which lines, which gf, which ions, and HOW the
   final value is combined (unweighted mean? both ions? one ion only?).
2. Transcribe the table to `data/reference/<el>_agss21/<el>_agss21_lines.csv`
   (wavelength_air_A, species, ep_eV, published_loggf, band, plus any per-line gf sigma and
   gf source the paper gives). Add a short `README.md` with the source page/table. Copy the
   shape of `data/reference/asplund2021_al/nordlander_lind_2017_analysis_lines.csv`.
3. Add one entry to `data/reference/line_sets/SOURCES.yaml` (name `<EL>_AGSS21`,
   `solar_literature: true`, bib, doi, citation, gf_basis). Copy the AL_AGSS21 entry.
   If the authors use an arithmetic mean, declare `aggregation: mean`; this follows
   the selection through nominal products and paired uncertainty responses. Existing
   selections without a declaration retain their median recipe.
4. Build it: `python -m pipeline.line_sets --element <El>` -- writes the per-species csv,
   graded csv and the REGISTRY.csv row.
5. Re-run prepare. Check A should now list the lines; lines missing from `canonical_gf`
   show an add plan. After reviewing the report: `--prepare --apply` writes the gf
   adoptions and culls. **That `--apply` is the human checkpoint -- report first.**

Check the full and graded counts separately. A registered set with zero graded lines
does not establish that the measurement uses the authors' gf: the current adoption
plan requires both a published gf and a per-line sigma. If the source gives no sigma,
leave it unknown and report the hold; do not transfer an experimental uncertainty to
a different theoretical or combined-feature gf. A collective accuracy bound is not
a published per-line standard deviation. Table wavelengths printed in nm can also
need an explicit canonical-ID join checked on species, wavelength and lower energy;
do not silently increase a wavelength-only matching tolerance.
Declared solar reference sets remain in coverage when their graded selection is empty,
with an explicit input hold rather than disappearing from the matrix.

## 4. Literature anchor (litscan)

Create `data/reference/litscan/<El>.yaml` by copying `Al.yaml` / `Si.yaml`: element, ions,
scale (3D-NLTE), and `entries` with Asplund+2021's value and sigma plus the other
published solar values (Lodders+2025, Scott+2015b, ...), each with citation and the
method. Every `bib` key must resolve. Values come from the papers, never from memory.

## 5. NLTE table (Sirius only)

Here `python` must be the PySME interpreter resolved by the path-register key
`engines.venv_pysme_python`, not the iSpec synthesis interpreter. Verify with
`python -c 'import pysme'` before the check; the synthesis environment can pass
`--check-paths` while lacking PySME (RYA-1250).

```bash
ssh sirius
cd <a checkout of this branch, NOT under ~/actions-runner>
python scripts/extend_nlte_grid.py --element <El> --ion <I-or-II> --solar-only   # non-regression check
python scripts/extend_nlte_grid.py --element <El> --ion <I-or-II> --write --jobs 3
```
Commit the updated `data/nlte_grids/<El>_Amarsi2020_PySME.csv` (+ `.prov.json`), then
re-run prepare: check D must read 0.

Coverage and failure provenance are keyed by ion. A neutral correction does not
cover an ionized line at the same wavelength. Ionized diagnostics require explicit
radiative/Stark broadening and component J matches; the synthesis verifies native
NLTE engagement and the LTE-reset environment. A non-monotonic or unbracketed
curve has no correction. Check D also accounts for documented failures, so zero
unaccounted lines does not by itself establish complete usable NLTE coverage.

### Published solar correction tables

A paper may publish LTE and NLTE abundances from the same equivalent widths
(e.g. AGSS21 Table A.1). Declare the source columns and explicit canonical joins
in `data/reference/solar_nlte_tables.yaml`. The shared adapter is solar-only and
emits cited additive corrections with the baseline atmosphere stated; it does
not represent an independent in-transfer 3D calculation. These sources are
fingerprinted with the run inputs and use the existing uncertainty recipe.

Structured prepare reports scope missing/adopted lines to each species/window.
Unknown/global blockers still hold execution. A missing infrared transition
therefore does not stop a prepared optical cell. Native Gerber NLTE routes also
consult the model atom's measured ion reach; a continuum-only ion is not an
NLTE-capable species.

## 6. Cull candidates

Prepare lists saturated lines absent from registered solar-literature selections. Post the list to the ticket;
Ryan says yes or no. Approved culls are written by `--prepare --apply`.

More precisely, the code checks membership in the **registered** solar-literature
sets. Absence from those sets is not evidence of absence from all solar literature.
Report contrary literature before approval (RYA-1250: Mg I 5711 A is absent from
AGSS21 Table A.1, but Lodders+2025 explicitly discusses its solar use).

Check the line disposition registry first. A synthesis-required disposition keeps
the line in synthesis and holds unsuitable EW/profile routes; saturation is not
a blanket exclusion. Dispositions must identify exact wavelengths that the shared
matcher can resolve. Prepare reports these as `route_requirements`.

## 7. Run and publish (hands-off)

```bash
python run_pipeline.py --star solar --element <El> --jobs 8 --publish > publish_<El>.log 2>&1
```
One command: measures every cell, builds the uncertainty budget legs, assembles, publishes
to `data/products/solar/<El>.json`, withdraws superseded products, writes the verdict and
coverage. Resolve permitted compute from current science authority and the ticket; historical
Mac runs are not general authorization. Expect a long run. Run it in the background and
watch the log; **do not edit code while it runs.** If it fails mid-run: stop, find the cause,
report to Ryan, fix only after he agrees.

Done when the coverage table (`python run_pipeline.py --star solar --element <El> --coverage`)
reads **STALE 0, NOT_RUN 0 or only host-routed cells, DONE**.

## 8. Check before pushing (local evidence and actual CI status)

```bash
python -m pytest -q -p no:cacheprovider tests/ > branch.txt   # ~36 min on the Mac
```
Run the same on `origin/main` when a baseline comparison is needed. Report both
existing failures and regressions; a matching baseline is not a passing suite. Inspect
the actual CI status rather than assuming a committed workflow is running or disabled. Coverage DONE is NOT proof the feed is clean -- the
feed tests (rya1096, rya1080, rya1034) are. Some tests rewrite tracked files
(`rya1080_guard.json`, a rya1169 csv): restore them before committing.

Apply the update triggers in `skills/codex-state-register/SKILL.md`; a state change
needs its register evidence and SEQUENCE line. Run the freshness checks after committing,
push and open the PR. Report any failing guard rather than bypassing it silently.

## 9. Website (separate repo: exoplanetcodex-site)

1. In `scripts/generate_cno_publication.py` add the element to `NAMES`, `ASPLUND2021`,
   `ASPLUND2021_SIGMA`, and a reader-facing name for its line set in `INDICATOR_NAMES`.
2. In `scripts/cno_selection.py` add it to `LINE_SET_SELECTORS` / `LINE_SET_SOURCE` if
   Asplund's value is one analysis of a line set (Si, Al, Na). If Asplund combines groups
   (e.g. Mg I and Mg II reported separately, then combined), say so in the ticket and
   implement that combination the way the C/N/O families are done. No recipe -> the page
   falls back to the best single product automatically.
3. In `scripts/generate_solar_report.py` add it to the appendix-link list.
4. Regenerate against the science checkout:
   `python scripts/generate_cno_publication.py --science-root <science> --elements <El>`
   then `python scripts/generate_solar_report.py --science-root <science>`.
5. Read the generated page: top value present, matches the Sun table, compares to
   Asplund, no internal IDs. Run `pytest tests/` (compare with main), open the PR.

## 10. Report (end of session)

First the zero-token summary: `python scripts/post_run_report.py --star solar --element <El>`
(`--dry-run` to preview; needs `LINEAR_API_KEY` in the environment; posts the verdict, cell
counts and problem lines to the element's literature dossier from `data/catalog/linear_dossiers.csv`,
e.g. Na -> RYA-735, Mg -> RYA-723; re-runs with no change post nothing). Then, by hand:
One ASCII comment on the element's ticket: the value vs Asplund, coverage numbers, what
was HELD and why, any orchestrator problem found (what an agent had to do by hand), the
PR links. Orchestrator problems also go on RYA-1233 so the pipeline gets better.

---

## Known gaps (an agent still has to do these by hand)

- Line-set transcription and litscan (steps 3-4) -- literature reading.
- NLTE table extension is run by hand on Sirius (step 5).
- Prepare and publish are two commands.
- One changed input re-runs every budget leg of the element.
- No host-wide CPU budget: do not start a second heavy run on a busy host (check `uptime`).
- Website per-element lists (step 9) are hand-edited.
