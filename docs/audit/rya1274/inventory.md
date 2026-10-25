# RYA-1274 baseline inventory

Baseline: `183310d8851629f43549e93269989032aa422bae`. Paths were enumerated from
`git ls-files`, with an additional `rg --files -uuu` search for nested instruction files.
No tracked or untracked root/nested AGENTS/CLAUDE file existed before this change.
The new `AGENTS.md` is the sole added instruction adapter.

## Project skills

- `skills/codex-data-audit/SKILL.md:1`
- `skills/codex-decision-log/SKILL.md:1`
- `skills/codex-merge-readiness/SKILL.md:1`
- `skills/codex-model-applicability/SKILL.md:1`
- `skills/codex-mr-code-brief/SKILL.md:1`
- `skills/codex-product-audit/SKILL.md:1`
- `skills/codex-scientific-code-review/SKILL.md:1`
- `skills/codex-state-register/SKILL.md:1`
- `skills/codex-task-intake/SKILL.md:1`
- `skills/codex-vald-extraction/SKILL.md:1`

## Agent instruction entry files

None in the baseline tracked tree.


## Hooks

- `hooks/pre-commit:1`
- `hooks/pre-push:1`

## CI

- `.github/workflows/ci.yml:1`

## Runbooks

- `docs/ORCHESTRATOR_RUNBOOK.md:1`
- `docs/runbooks/preflight_check_rya905.md:1`
- `docs/runbooks/sirius_data_drive_rya419.md:1`

## Development and authority

- `LEDGERS.md:1`
- `DEV_CYCLE.md:1`
- `DATA_SAFETY.md:1`
- `docs/SCIENCE_STANDARDS.md:1`
- `docs/ELEMENT_PROTOCOL.md:1`
- `docs/SCIENCE_PRODUCT_PACKAGE.md:1`
- `docs/CONVENTIONS.md:1`
- `docs/audit/rya926_skills_governance.md:1`

## Resource and orchestration owners

- `config/path_register.yaml:1`
- `config/hosts/mac.env:1`
- `config/stars.yaml:1`
- `data/refs/README.md:1`
- `data/refs/bibliography.csv:1`
- `data/reference/line_sets/SOURCES.yaml:1`
- `data/reference/line_sets/REGISTRY.csv:1`
- `data/reference/gf_error_model.csv:1`
- `data/catalog/model_registry.csv:1`
- `data/catalog/holdings_manifest_registry.csv:1`
- `run_pipeline.py:1`
- `pipeline/run_matrix.py:1`
- `pipeline/run_descriptor.py:1`
- `pipeline/run_sweep.py:1`
- `pipeline/element_prepare.py:1`
- `pipeline/element_publish.py:1`
- `pipeline/coverage_report.py:1`
- `pipeline/element_verdict.py:1`

## Skill review dispositions

All listed SKILL.md bodies were inspected. No new skill or duplicate procedure was needed.

| Existing procedure | Disposition |
|---|---|
| task-intake / implementation brief | Reused for issue history, current Git, authority and debt; the issue explicitly authorizes its Linear handoff and bounded follow-up filing. |
| scientific-code-review | Reused for wrong-data, silent-failure and provenance priorities. |
| model-applicability | Reused for deck/ion/domain disposition; no science run. |
| product-audit | Read; its product-signoff procedure is conditional and no product was signed off. Historical numeric examples are not canonical inputs. |
| state-register | Applied for the validated readiness findings; no scientific verdict changed. |
| data-audit / VALD extraction | Read for held inputs, flux/frame/provenance and manual VALD intake; no acquisition or intake performed. VALD raw-location prose predates RYA-1228 stewardship: resolve current path register and data policy before reuse. |
| decision-log | Read; memory is a convenience and cannot replace committed authority. No new science decision was ratified. |
| merge-readiness | Read; no merge. Sirius production validation remains a separate requirement, not established by Mac fixture tests. |

## Integration and machine-local limits

There is no committed `.codex/`, `.claude/` or `.agents/` app configuration in this tree.
The project skills are ordinary tracked files under `skills/`; discover/read them explicitly,
without assuming automatic installation. `git config --show-origin --get core.hooksPath`
resolved `hooks` from the shared repository-local Git config on this workstation. This
installation evidence does not establish installation in a fresh clone.

Machine-local presence checked: `~/.codex/config.toml` and the external project parent's
.claude/settings.local.json exist. Their secret-bearing contents and account/app UI state
were not copied into this audit or modified. No claim of Codex/Claude tool/config equivalence
is made. No global AGENTS file was found in `~/.codex`. Session-provided developer/user
instructions still govern this run.

Official instruction-discovery reference inspected on 2026-10-25:
[Custom instructions with AGENTS.md](https://learn.chatgpt.com/docs/agent-configuration/agents-md).
This adapter remains small; nested/session instructions and actual tool limitations apply.

GitHub workflow API read: `.github/workflows/ci.yml` is `active`. This contradicts historical
"CI is off" wording in `hooks/pre-push:5` and is not proof that an available runner passes
the suite. The committed workflow still contains LFS materialisation of raw VALD inputs
purged in RYA-1228; no workflow execution/download was triggered.
