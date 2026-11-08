# RYA-1352 validation and handoff

Base SHA: `3baabba95f728342f069f708dc815e729a1d0447`.
Branch/worktree: `ryandamienschmitt/rya-1352-grade-governance-audit-codexdeep`,
`codex/rya1352`. Host: Darwin; audit/test interpreter: existing
`codex/venv/bin/python`, Python 3.12.13 / pandas 3.0.3.

| Check | Exit/result |
| --- | --- |
| `git fetch origin` | Initial sandbox DNS failure; authorized retry succeeded, base pinned above. |
| `python scripts/audit_grade_governance_rya1352.py --base 3baabba95f728342f069f708dc815e729a1d0447 --output docs/audit/rya1352/baseline.json` | 0; 1,388 protected files unchanged, all-row classifier unchanged, 10 feeds with product arrays counted. Full report includes hashes for other protected feed files too. |
| `python -m pytest tests/test_grade_governance_rya1352.py tests/test_rya1233_gf_graded_definition.py tests/test_rya1233_grade_every_element.py -q` | 0; 11 passed. The audit covers all canonical accepted IDs and protected hashes; tests prove mutation detection on a disposable Git fixture, legacy grade behavior and synthetic invalid-sigma admission. |
| `python -m pipeline.refinement_debt_join --report` | 0 with project interpreter; 21 registry rows, 10 not-fired ticket rows, 5 no-ticket rows. No debt cleared. Initial system Python attempt failed missing pandas. |
| `python3.12 scripts/generate_sources_page.py --audit-library "<resolved refs.library>"` | 1; pre-existing bibliography vocabulary failures on civis2013 and mallinson2024; library traversal not certified. |
| `python3.12 scripts/check_register_freshness.py` | 0; register structure/freshness passes, existing stale generated tracker notice retained. |
| `git diff --check` | 0 after normalizing only the two added ledger row line endings. |
| Final audit to `/tmp/rya1352-final-audit.json` | 0; same canonical accepted-ID hash and zero protected changes after final edits. |
| `python3.12 scripts/check_register_freshness.py --since-main` after commits | 0; zero state surfaces changed; no register bump required. |
| `git push -u origin <task branch>` | 1; pre-push RYA-1184 hook finds existing RYA-1349 merge missing from SEQUENCE on origin/main. Version pointer and freshness pass. Used the hook's explicit `--no-verify` option to push this bounded review branch; unrelated ledgers remain untouched. |

Self-review applied the scientific-code-review skill: no scientific policy was
ratified, no classifier modified, no production evidence record adopted, no engine
or spectrum touched. The JSON contract is a field-design proposal, not a deployed
validator. Future validation/migration gates are listed as owed in the dossier.
Two characterized run-bug-ledger rows retain the conflict and latent invalid-sigma
admission. RYA-1174's existing Al provenance defect is referenced without duplicating
its ownership. One verified methodological bibliography row was added; existing
NIST ASD source key is `nist_asd`. No source page was regenerated or deployed.

Review gate: settle the August primary-only versus October published-uncertainty
policy and the new RYA-1351 contract, then review the proposed evidence contract
and migration. This PR supplies no scientific adoption or publication approval.
No merge is authorized. RYA-1352 remains In Review pending Ryan's sign-off.
Exact final SHA, PR URL, post-commit freshness and actual CI state are in the Linear
handoff; no clean CI claim is inferred from these local tests.
