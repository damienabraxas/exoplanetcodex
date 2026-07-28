# RYA-464 — Per-star multi-instrument arm wiring for benchmark C/O

## The wall (RYA-348 Step-0 finding)

`pipeline/cno_synthesis.py` registered a **hardcoded single region** —
`REGIONS = {'vis': HARPS_VIS}` — and the multi-arm orchestrator `run_phase_a` loaded the
ESPRESSO/UVES arms through `_load_reflected_solar_arm` **unconditionally**. That loader
reads `vesta_<inst>_manifest.csv` (Vesta reflected sunlight, RYA-370/372). So for any
**non-solar** star, the "primary O via UVES O I 777" arm would have silently synthesized
against **sunlight**, not the star. Procyon's runnable C/O scope was HARPS-VIS only.

## The fix — a per-star arm registry (the unlock)

Arms are now declared **per star**, each carrying its `RegionConfig`, diagnostics, a
spectrum-loader kind, and a **readiness flag**:

```
ArmWiring(name, region, diagnostics, loader, ready, defer_reason, provenance)
STAR_ARMS = {'solar': {...3 arms...}, 'procyon': {...4 arms...}}
```

- `star_arm_registry(star_id)` — the star's arms (substring-matched; **loud-fails** for an
  undeclared star — no silent solar-geometry default).
- `available_arms(star_id)` → `(ready_names, {deferred_name: reason})`.
- `resolve_arm_spectrum(star_id, arm)` — dispatches the loader by kind and **loud-fails
  (`ArmNotWired`)** when an arm is unready *or* when a non-solar star is routed at the
  `reflected_solar` loader. This is the line that closes the latent bug: a non-solar star
  can never be synthesized against Vesta.
- `run_phase_a(star_id, arms=None)` builds its arm set from the registry (not a hardcoded
  list), runs the ready arms, and **defers** the rest with their reasons (printed at
  runtime + recorded in `<star>_phase_a_cross_arm.json`). `arms='all'` / CLI `--arms all`
  uses the star's full declared set; `--list-arms` prints the registry.

The solar/Vesta path is **one case** of this mechanism — solar's three arms are all
`ready` via the existing loaders, so `run_phase_a('solar')` dispatches to the identical
`run_cno`/`_load_reflected_solar_arm` calls (regression-preserved; see below).

## Arm status (this build)

| star | arm | instrument | range (Å) | status | reason / loader |
|---|---|---|---|---|---|
| solar | harps | HARPS | 3780–6910 | READY | `harps_normalized` |
| solar | espresso | ESPRESSO | 3780–7890 | READY | `reflected_solar` (Vesta) |
| solar | uves | UVES | 3760–9470 | READY | `reflected_solar` (Vesta) |
| procyon | harps | HARPS | 3780–6910 | READY | `harps_normalized` (RYA-273) |
| procyon | uves | UVES | 3760–9470 | DEFERRED | RYA-272 loader not built + UVES spectra not staged (gate RYA-271). **O I 777 = primary O** |
| procyon | uv | STIS | 1150–3200 | DEFERRED | HST STIS/COS loader not built; only RYA-222/262 audit inventories |
| procyon | ir | CRIRES+ | 15000–24000 | DEFERRED | telluric-gated (RYA-373); no 2.3 µm CO overtone staged (RYA-351) |

The deferred arms are **declared, not fabricated**: the UVES diagnostic set reuses the
existing `Diagnostic` objects (O I 777 primary + [O I] 6300 + C I 5052/5380 + N I 8216 +
CN red — real line data), and the UV/IR regions carry only factual instrument spans with
**no invented line data** (their diagnostics populate from the RYA-262/425 audits).

## Why this is the right scope

The arm **loaders** (RYA-272 UVES, an HST STIS loader) are not built and the per-arm
**spectra** are not staged (only HST audit inventories exist; no `reflected_solar` Vesta
products on disk either). So no arm beyond HARPS can be *executed* honestly today. This
ticket delivers the **integration layer** — the per-star registry + loader dispatch +
loud deferral — so each arm plugs in the moment its loader+data+audit land (RYA-272/271,
HST, RYA-373/425), for **every** benchmark (Procyon, α Cen A/B, τ Boö, 55 Cnc), not just
Procyon. It unblocks RYA-348 (multi-instrument C/O) and RYA-404 (per-region 27-element).

## Regression — solar bit-identity

The Vesta `reflected_solar` products are not staged in this environment, so solar Phase A
cannot be re-executed here. Bit-identity is therefore established **structurally**: solar's
registry arms (`harps`→`VIS_DIAGNOSTICS`/`harps_normalized`; `espresso`→`ESPRESSO_OPT`/
`ESPRESSO_DIAGNOSTICS`/`reflected_solar`; `uves`→`UVES_OPT`/`UVES_DIAGNOSTICS`/
`reflected_solar`) dispatch to the **same** loaders and diagnostics in the same order as the
pre-RYA-464 code (guarded by `test_arm_registry_rya464.py`). The science payload of
`solar_phase_a_cross_arm.json` (`cross_arm`, `per_arm`) is unchanged; only descriptive
wrapper keys (`arms_run`, `arms_deferred`) were added. The committed solar baseline stands.

## Leg validation (when arms light up)

The keystone cross-check **O I 777 (UVES) vs [O I] 6300 (HARPS)** — the analog of the
solar Kitt-Peak +0.036 zero-point catch (RYA-459/460) — runs through the existing
`_print_cross_arm_table` (per-arm value, spread, FLAGGED-DISAGREEMENT never averaged) the
moment the UVES arm is `ready`. No new cross-check code is needed.
