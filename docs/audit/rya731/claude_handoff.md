# Claude handoff for Solar Ni

## Start here: Ni II near-UV atmosphere uncertainty

The run finished. The publisher holds eight Clear 2026 near-UV LTE candidates (six Sirius and two Mac) because `model_atmosphere` is unresolved. The three visible holdings retain no lines. Neither is a synthesis process failure.

1. Read the [RYA-1226 close-out](https://linear.app/ryans-adventure-zone/issue/RYA-1226), science PR #559 at `150bf396`, and the RYA-1032 `atmosphere_nuisance` evidence. Fe's near-UV precedent resolved this component using a paired ATLAS9 to MARCS.GES comparison on the same accepted pool, with dimensionality and NLTE held fixed. The Fe value belongs to that Fe pool alone. Check the current Asplund budget contract before selecting how to record the Ni II term.
2. Trace `scripts/rya1230_cno_budget.py` at its `scott_atmosphere()` and `marcs` branch. `pipeline/asplund_atmosphere.py` finds Scott rows for Ni I only, then falls through to a CNO-only table. With nominal-only publication, there is no explicit `marcs` leg to fill the term. That's the observed gap. It is independent of the Clear oscillator strengths and of NLTE coverage.
3. Inspect the archived [Ni II nominal comparison files](../../../data/audit/rya731_ni2_clear2026/nominal/) from Sirius and Mac. Compare the 1D LTE ATLAS9 and Gerber MARCS.GES fits for each Kitt Peak holding on matched accepted indicators. Match by wavelength and excitation. Verify the comparison isolates atmosphere, then reuse the applicable paired-budget route. Keep the mean-3D LTE contrast separate from an atmosphere-grid term. Run only the scoped comparison if the existing artifacts cannot supply it.
4. Reassemble the held Ni II units and check the publication dispositions. The pinned atom has no Ni II bound-bound transitions. Keep synthetic stage products, feed rows, and the publication report in sync. Regenerate the site from the resulting committed feed and preserve the 223 admitted Ni I rows.

## Work locations and current release

- [Science PR #633](https://github.com/damienabraxas/exoplanetcodex/pull/633), branch `ryandamienschmitt/rya-731-ni-release`, Mac worktree `~/codex/rya731-release`. The product feed is v1.161 at `d7169665`, with 223 admitted Ni I rows and one quarantined row.
- [Website PR #100](https://github.com/damienabraxas/exoplanetcodex-site/pull/100), branch `ryandamienschmitt/rya-731-solar-nickel`, Mac worktree `~/codex/rya731-site-release`, head `aa5e63f`. The page, 29-page PDF and full Sun table are pinned to science feed commit `d7169665`.
- Sirius run checkout `~/scratch/rya731`. Logs and the combined Sirius scope are under `~/scratch/rya731_run`. The Ni I VIS and near-UV label patches completed, as did both theoretical passes. Mac OpenMULTI3D completed after restoring the pinned engine loader. Neither host has an active Ni task for this release.
- Publication dispositions are in `data/results/orchestrator/budget/solar_Ni/publish_summary_sirius.json` and `publish_summary.json`. The combined run record is `data/results/orchestrator/solar_Ni_latest.json`. The Mac coverage and unit bookkeeping are host-scoped.
- OpenMULTI3D is a 1D comparison here. RYA-1392's separate Fe and full-3D oxygen work remain outside this Ni release.

No abundance values were adjusted to match a reference.
