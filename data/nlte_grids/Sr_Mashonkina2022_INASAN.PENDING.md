# Sr II Mashonkina/INASAN NLTE grid — MANUAL PULL OWED (RYA-433)

**Status: BLOCKED on a manual (browser) pull.** RYA-433 set out to vendor the named-PRIMARY
Mashonkina 2022 / INASAN Sr II 4077/4215 grid so the cited primary == the file on disk ==
the applied grid (closing the RYA-421 provenance defect: a Mashonkina-primary flag over a
Bergemann-only disk). The grid could **not** be obtained programmatically — so per the brief's
CRITICAL condition this STOPs and flags a manual pull, and does **not** silently relabel
Bergemann as primary.

## Why it is blocked (evidence)
INASAN exposes the grid only through a JS/AJAX form whose data endpoint is WAF-blocked:
- Form page: `https://spectrum.inasan.ru/nLTE2` — loads (HTTP 200).
- Data endpoint (from the page's `SubmitRequest()`):
  `GET https://spectrum.inasan.ru/nLTE.cgi?elem=Sr2&lam=<line>&Teff=<T>&logg=<g>&fe=<[Fe/H]>&abn=<[Sr/Fe]>`
- **Every programmatic attempt returns HTTP 403 Forbidden**: urllib and curl; GET and POST;
  with full browser headers (User-Agent, Referer `https://spectrum.inasan.ru/nLTE2`,
  Origin, `X-Requested-With: XMLHttpRequest`, `Sec-Fetch-*`, cookie jar). All three CGIs
  (`nLTE.cgi`, `nLTE2.cgi`, `nLTE3.cgi`) are 403. The interactive browser form works, so the
  block is on non-interactive clients (TLS-fingerprint / JS-token level we cannot replicate).

## What WAS captured (from the page JS, no 403) — for the manual pull + the audit
- **Lines (`lams_Sr2`):** `4077.71`, `4215.54` (the resonance doublet — matches our linelist).
- **Grid axes (`min_Sr2` / `max_Sr2` = [Teff, logg, [Fe/H], [Sr/Fe]]):**
  - Teff: **4000 – 6500 K**
  - logg: **0.5 – 5.0**
  - **[Fe/H]: −5 to −2  → CEILING = −2.0** (a metal-poor VMP MARCS grid; even lower than
    Bergemann/INSPECT's 0.0). There is **no solar node** — so the brief's "reproduce the solar
    anchor" cannot apply to INASAN; the in-hull anchor to reproduce on a future pull is at
    [Fe/H] ≤ −2. solar..55 Cnc(+0.32) are out-of-hull → `NLTE_unavailable`, loud (RYA-409).
  - [Sr/Fe] (`abn`): −1.5 to 1.0
- **Citation (verbatim from the INASAN landing page `https://spectrum.inasan.ru/nLTE/`):**
  > "Mashonkina, Sitnova, Pakhomov, 2016, Astronomy Letters, 42, 606" (database);
  Sr II model atom: Mashonkina et al. 2022; distribution: Mashonkina et al. 2023.
- **Retrieval date of this metadata:** 2026-06-24.

## Manual-pull recipe (run in a browser session at `https://spectrum.inasan.ru/nLTE2`)
1. Select `elem = Sr2`; for each line `lam ∈ {4077.71, 4215.54}`.
2. Sweep the in-grid nodes: Teff ∈ {4000..6500}, logg ∈ {0.5..5}, [Fe/H] ∈ {−5..−2},
   [Sr/Fe] (`abn`) ∈ {−1.5..1}. Use the batch form (`form2`: upload a parameter file,
   result-as-file) to capture ALL nodes at once — do NOT hand-transcribe a subset.
3. Save as `data/nlte_grids/Sr_Mashonkina2022_INASAN.csv` with columns
   `element, ion, wave_A, teff_K, logg, feh, sr_fe, delta_nlte`, plus a `.prov.json`
   recording URL + retrieval date + the citation above.
4. Flip `NLTE_CORRECTION_ELEMENTS['Sr']['grid']` to the INASAN file (the APPLIED primary),
   keep `Sr_Bergemann2012_INSPECT.csv` as the cross-check; reproduce the grid's in-hull
   anchor (validate-don't-tune); cross-check INASAN vs INSPECT at the **[Fe/H] −3..−2 overlap**
   (within ~0.15 dex). Sr stays GET-DATA-pending (RYA-428 — the Sr II measurement is separate).

## Interim state (this branch)
- Applied/on-disk grid stays `Sr_Bergemann2012_INSPECT.csv` (Bergemann) — the only Sr II grid
  we could vendor — but it is **labelled cross-check/working, NOT silently "primary"**;
  the registry `ref` + the regime-map note now say loudly that the Mashonkina/INASAN primary
  is owed via this manual pull (no silent flag/disk mismatch).
- The RYA-433 accuracy correction to the regime note IS applied: the Sr II resonance NLTE is
  small-NEGATIVE near-solar, sign-changes, and grows POSITIVE toward metal-poor (overionization)
  — not a monotonic shrink to zero.
