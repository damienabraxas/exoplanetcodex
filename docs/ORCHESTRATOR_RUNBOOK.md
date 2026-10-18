# Orchestrator Runbook -- one solar element, A to Z

**Audience:** any agent (Claude, Codex, ...) or person taking a new element through the
pipeline. **This file is the source of truth for the procedure; it lives in git and changes
with the code.** If a step here disagrees with the code, the code wins -- fix this file in
the same PR.

Worked examples: Si and Al (RYA-1233, science PR #583, site PRs #53/#54). Copy what they did.

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

## 2. Prepare -- what the element still needs (read-only)

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
| **B gf source ... unclassified** | A gf source has no correlation entry. | Add a row to `data/reference/gf_error_model.csv` (copy the nearest existing source; cite why). |
| **F2 literature need** | Bands with no error-stated gf: no Codex/Deep pool there. | Not blocking. Record it in the ticket; do not invent gf errors. |
| **E holdings not serving** | A spectrum does not cover a window at rest. | Report it; do not re-prepare spectra. |

## 3. Register the published line set (the Asplund recipe)

1. Find the element in Asplund+2021 Sect. 3 and its line table (Table A.1 for Na/Mg/K/Ca,
   A.2 for Fe; otherwise the paper Asplund cites, e.g. Amarsi & Asplund 2017 for Si,
   Nordlander & Lind 2017 for Al). Note: which lines, which gf, which ions, and HOW the
   final value is combined (unweighted mean? both ions? one ion only?).
2. Transcribe the table to `data/reference/<el>_agss21/<el>_agss21_lines.csv`
   (wavelength_air_A, species, ep_eV, published_loggf, band, plus any per-line gf sigma and
   gf source the paper gives). Add a short `README.md` with the source page/table. Copy the
   shape of `data/reference/asplund2021_al/nordlander_lind_2017_analysis_lines.csv`.
3. Add one entry to `data/reference/line_sets/SOURCES.yaml` (name `<EL>_AGSS21`,
   `solar_literature: true`, bib, doi, citation, gf_basis). Copy the AL_AGSS21 entry.
4. Build it: `python -m pipeline.line_sets --element <El>` -- writes the per-species csv,
   graded csv and the REGISTRY.csv row.
5. Re-run prepare. Check A should now list the lines; lines missing from `canonical_gf`
   show an add plan. After reviewing the report: `--prepare --apply` writes the gf
   adoptions and culls. **That `--apply` is the human checkpoint -- report first.**

## 4. Literature anchor (litscan)

Create `data/reference/litscan/<El>.yaml` by copying `Al.yaml` / `Si.yaml`: element, ions,
scale (3D-NLTE), and `entries` with Asplund+2021's value and sigma plus the other
published solar values (Lodders+2025, Scott+2015b, ...), each with citation and the
method. Every `bib` key must resolve. Values come from the papers, never from memory.

## 5. NLTE table (Sirius only)

```bash
ssh sirius
cd <a checkout of this branch, NOT under ~/actions-runner>
python scripts/extend_nlte_grid.py --element <El> --solar-only   # non-regression check
python scripts/extend_nlte_grid.py --element <El> --write --jobs 3
```
Commit the updated `data/nlte_grids/<El>_Amarsi2020_PySME.csv` (+ `.prov.json`), then
re-run prepare: check D must read 0.

## 6. Cull candidates

Prepare lists saturated lines that no solar analysis used. Post the list to the ticket;
Ryan says yes or no. Approved culls are written by `--prepare --apply`.

## 7. Run and publish (hands-off)

```bash
python run_pipeline.py --star solar --element <El> --jobs 8 --publish > publish_<El>.log 2>&1
```
One command: measures every cell, builds the uncertainty budget legs, assembles, publishes
to `data/products/solar/<El>.json`, withdraws superseded products, writes the verdict and
coverage. Expect hours (Si: ~9 h of budget legs on the Mac). Run it in the background and
watch the log; **do not edit code while it runs.** If it fails mid-run: stop, find the cause,
report to Ryan, fix only after he agrees.

Done when the coverage table (`python run_pipeline.py --star solar --element <El> --coverage`)
reads **STALE 0, NOT_RUN 0 or only host-routed cells, DONE**.

## 8. Check before pushing (CI is not reliable; run it yourself)

```bash
python -m pytest -q -p no:cacheprovider tests/ > branch.txt   # ~36 min on the Mac
```
Run the same on `origin/main` and compare failure sets. **Only branch-only failures
matter**; main has known failures. Coverage DONE is NOT proof the feed is clean -- the
feed tests (rya1096, rya1080, rya1034) are. Some tests rewrite tracked files
(`rya1080_guard.json`, a rya1169 csv): restore them before committing.

Then bump `CODEX_STATE_REGISTER.md` (new v-next row; the pre-push hook enforces it) and a
landing line in `SEQUENCE.md` for the ticket, push, open the PR.

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
