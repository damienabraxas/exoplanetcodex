# RYA-944 — Delbouille/Liège disk-center intensity atlas: intake audit

**Verdict: GO WITH CAVEATS.** The atlas is acquired, checksummed, staged on Sirius and
registered. Two of the ticket's three stated "CRITICAL intake properties" are confirmed;
**the wavelength-medium property is REFUTED** — see §2, which is the finding that matters.

| axis | ticket said | MEASURED | verdict |
|---|---|---|---|
| product type | disk-center intensity atlas | disk-center intensity, Jungfraujoch | ✅ confirmed |
| geometry | mu = 1.0 | mu = 1.0 (centre of disk) | ✅ confirmed |
| normalization | continuum = 100% | median 0.988 in a clean window | ✅ confirmed |
| **wavelength medium** | **VACUUM, convert vac→air** | **AIR as delivered** | 🔴 **REFUTED** |
| coverage | 3000–10000 Å | 3000.0000–10000.0000 Å, no gap > 0.004 Å | ✅ confirmed |
| telluric | "low but nonzero — record, do not assume zero" | present; H2O 7–8× cleaner than Kitt Peak, O2 comparable | ✅ confirmed + quantified |

SHA-256 `72a1d97b4201872425b37d2f15982357a37fa52647ad83b469f2bb28748ef794`
(10,067,835 B, 2,380,014 points).

## 1. What was acquired

BASS2000 serves **three different atlases behind one endpoint**, spliced by wavelength:

| segment | range | source | resolution |
|---|---|---|---|
| UV | 670–1609 Å | SOHO/SUMER, Curdt+2001 | 0.04 Å |
| **Visible** | **3000–10000 Å** | **Jungfraujoch — Delbouille, Neven & Roland (1972)** | 0.002 Å |
| IR | 10000–54000 Å | Kitt Peak — Delbouille, Roland, Brault & Testerman (1981) | 0.004 cm⁻¹ |

Only the **visible arm** is this holding. The splice is silent — the endpoint returns data
across the seam with no marker — and the sampling step changes at 10000 Å (0.0020/0.0040 →
0.0048 Å), which is the only outward tell. Requesting 9900–10100 Å in one call therefore
returns **two instruments concatenated**. The loader is bounded at 10000 Å for that reason.

Note the citation variance: the ticket cites *Delbouille, Roland & Neven 1973*; BASS2000
itself attributes the visible arm to *Delbouille L., Neven L., Roland G. (1972)*. The
provenance string records BASS2000's own attribution.

## 2. 🔴 The wavelength medium is AIR, not vacuum

The ticket states as a CRITICAL property, "verified this session": *"**VACUUM wavelengths.**
Pipeline is air-Å throughout. Convert vac→air at the loader boundary."* Applying that
conversion would shift every wavelength by **−1.1 to −2.4 Å** — roughly 200 solar line
widths — and silently corrupt every abundance derived from this atlas.

Measured line centroids against air rest wavelengths, on the delivered data with **no
conversion applied**:

| line | air rest | measured | Δ vs air | Δ vs vacuum |
|---|---:|---:|---:|---:|
| Fe I 4045 | 4045.812 | 4045.8203 | +8.3 mÅ | −1134.7 mÅ |
| Fe I 4383 | 4383.545 | 4383.5532 | +8.2 mÅ | −1223.5 mÅ |
| Mg I b2 | 5172.684 | 5172.6907 | +6.7 mÅ | −1434.1 mÅ |
| Mg I b1 | 5183.604 | 5183.6071 | +3.1 mÅ | −1440.6 mÅ |
| Na D2 | 5889.951 | 5889.9577 | +6.7 mÅ | −1625.7 mÅ |
| Na D1 | 5895.924 | 5895.9270 | +3.0 mÅ | −1631.0 mÅ |
| Fe I 6430 | 6430.846 | 6430.8495 | +3.5 mÅ | −1774.0 mÅ |
| Hα | 6562.797 | 6562.8107 | +13.7 mÅ | −1799.3 mÅ |
| Fe I 8688 | 8688.624 | 8688.6330 | +9.0 mÅ | −2377.6 mÅ |

Nine lines spanning 4045–8688 Å, all AIR, by a margin of ~170×. A further two Fe lines
(6430.846 / 6432.680) and six IR lines beyond this holding's range gave the same verdict.

**The residual is solar physics, not calibration error.** Expressed as a velocity the mean
is **+371 m/s**, inside the expected window for solar gravitational redshift (+633 m/s)
minus convective blueshift (−200 to −400 m/s). The decisive case is Hα: it forms high in
the chromosphere where convective blueshift nearly vanishes, and it returns **+626 m/s**
against a textbook +633 m/s. The wavelength scale is sound at the mÅ level.

**Why the ticket's own smoke test could not have caught this.** The spec's check is

```python
idx = np.argmin(np.abs(wav-6430.846)); assert abs(wav[idx]-6430.846) < 0.002
```

`argmin` returns the nearest **grid point**, and on any continuous grid with step ≤ 0.004 Å
the nearest grid point is within 0.002 Å of *any* requested value — in air, in vacuum, or
shifted by an arbitrary constant. The assertion is a grid-spacing test wearing a medium
test's label; it passes unconditionally. The medium can only be measured from a line
**centroid**, which is what the ticket's prose demanded and what was done here.

Registered as `wavelength_medium=air`, `air_vacuum_conversion_required=no`. **No conversion
is applied at the loader boundary, and none must be added.**

## 3. Tellurics: present, and weak in water only

An absence is a hypothesis, never a conclusion (RYA-833), so the atlas was measured
*between* two atlases whose telluric state is already established, with a clean-continuum
window as the control that proves the metric discriminates at all.

Percent of points below 0.5 (full table in `telluric_control.txt`):

| window | **Delbouille** | Kitt Peak (retains) | IAG-Baker2020 (corrected) |
|---|---:|---:|---:|
| O2 B-band 6867–6884 | 12.61 | 14.93 | 0.00 |
| O2 A-band 7594–7685 | 26.32 | 33.54 | 0.18 |
| H2O 8100–8400 | **0.44** | 3.53 | 0.29 |
| H2O 9280–9600 | **3.99** | 29.27 | 0.00 |
| clean continuum 6425–6525 | 1.26 | 1.07 | 1.10 |

All four atlases agree to ~0.2 pp in the clean window, so the metric is measuring the
atmosphere and not the atlas. Delbouille's minimum reaches 0.000 inside both O2 bands:
**saturated telluric absorption is present** and lines there are not measurable.

The pattern is exactly what Jungfraujoch's 3580 m altitude predicts. Water vapour scale
height is ~2 km, so H2O drops steeply with altitude — 7–8× cleaner than Kitt Peak. O2 is
well mixed, so it falls only with the pressure ratio — essentially unchanged. The atlas is
a *dry*-site atlas, not a corrected one.

⇒ `telluric_basis = line_selection`, `telluric_applied = not-applied`. Treating this atlas
as telluric-free because Jungfraujoch is high would measure the O2 A-band as solar.

## 4. Normalization

Continuum sits at **1.0** as documented ("continuum estimated locally, normalized to
100%"). Median in clean windows: 0.988 (6425–6428), 0.982 (7480–7483), 0.979 (5240–5243),
0.967 (8710–8713). The ticket's normalization assertion (0.97 < median(6425–6428) < 1.03)
**passes at 0.9876**.

Two caveats for anyone fitting against it:

- The maximum over all 2.38 M points is **0.9959** — no pixel reaches exactly 1.0. A fit
  that pins the continuum to exactly 1.0 will carry a ~0.4% offset.
- The continuum is estimated **locally in wavelength**, so intensities are in *local*
  arbitrary units. This atlas carries **no absolute flux calibration and no usable
  broad-band continuum shape**. It is the opposite of the Kurucz `irradthu` case
  (RYA-929, `PRE_NORMALISED=False`) and must be registered as pre-normalised.

## 5. Incidental findings — NOT fixed here (this ticket is data intake only)

**5a. 🔴 `iag_fts_solar_atlas` is catalogued `telluric_basis=corrected` but points at the
telluric-RETAINING file.** We hold two IAG atlases with opposite telluric states:

- `iag_reiners2016/spvis.dat.gz` — the **base** atlas, tellurics **present**
- `iag_baker2020/iag_telfree_solaratlas.fits` — the **telluric-free** derivative

`solar_reference_holdings_rya708.csv` routes `iag_fts_solar_atlas` to the **Reiners** file,
while `instrument_catalog.csv` marks that instrument `telluric_basis=corrected`. Per
`telluric_policy.exclusion()`, basis `corrected` returns `""` — **no line inside any O2 or
H2O band is excluded**. Measured on the file actually routed: 46.25% of the O2 A-band and
51.63% of H2O 9280–9600 lie below 0.5, reaching 0.001. Any line drawn from that holding
inside a telluric band is atmospheric absorption measured as solar.

The comparison table in `pipeline/telluric_policy.py` quotes IAG at 0.1% / 0.0% in those
bands. Baker2020 measures 0.18% / 0.00%; Reiners2016 measures 46.25% / 51.63%. **The
docstring's numbers are Baker's**, so the table and the routing describe different files.
Logged to the run bug ledger; needs its own ticket.

**5b. `kpno_solar_atlas` manifest row is stale.** `solar_reference_holdings_rya708.csv`
gives `host=mac`, path `/Users/ryanschmitt/Documents/Exoplanet Codex/.../Kitt Peak Flux
Atlas`. That directory is empty on the Mac; the atlas lives on Sirius at
`/mnt/codex-data/spectra/Solar Calibration/Kitt Peak Flux Atlas` after the 2026-08-16
migration.

## 6. Reproduce

```
python data/audit/rya944_delbouille_intake/fetch_bass2000.py    # 28 chunks, 3000-10000 A
python data/audit/rya944_delbouille_intake/assemble_audit.py    # dedupe seams, audit, SHA
python data/audit/rya944_delbouille_intake/telluric_control.py  # Sirius: 4-atlas control
```
