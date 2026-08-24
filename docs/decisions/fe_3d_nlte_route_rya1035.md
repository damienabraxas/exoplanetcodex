# Fe 3D-NLTE route-finding — RYA-1035

**Status: Step 0 answered. Verdict HAVE. The build branch is not scoped.**
Date 2026-08-24. Artifacts: `data/results/rya1035/fe_3d_route_census.{csv,json}`,
generator `scripts/rya1035_fe_3d_route_census.py`.

---

## The one-line answer

RYA-1035 said: *"Step 0 (cheap, do FIRST — do not build what we can consume) … A HAVE here
collapses the whole ticket to a wiring job."*

**It is a HAVE.** `NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin` has been on the same MPG
Keeper share we already fetch every Gerber deck from **since 2021**, and TSFitPy's own
canonical downloader lists it as `[Fe] 3d_bin_link`. We fetched Al, Cr, Eu and Y from that
share and never opened the Fe folder.

---

## 🔴 Three claims in the ticket and the register are wrong, and one is wrong in our favour

### 1. "Fe has NO deck (`T2_BUILD_OWED`, build-our-own only)" — FALSE

`CODEX_STATE_REGISTER.md` v117 ends with that sentence. It was read off **our disk** and
never checked against **the source**. The Fe ⟨3D⟩ deck is 72 MB zipped / 92.9 MB unpacked —
by far the *smallest* ⟨3D⟩ deck on the share, because Fe carries one A(X) per node where Al
carries 31.

This is [[feedback_absence_is_a_hypothesis]] and
[[feedback_status_surfaces_must_read_the_code]] in one: an absence measured on the wrong
side of the fetch boundary.

The same sweep, run live, says the gap is much wider than Fe:

| ships a `STAGGERmean3D` deck | we hold it |
| --- | --- |
| Al · Ba · Ca · Co · Cr · Eu · Fe · H · Mg · Mn · Na · Ni · O · Si · Sr · Ti · Y — **17 of 18** | Al · Cr · Eu · Y — **4** |

(`CH` is the 18th folder and has no mean-3D deck.) **Every "no ⟨3D⟩ deck for X" claim in
the matrix is suspect on the same grounds** — that belongs to RYA-1015/RYA-1029, not here,
but it should not wait for each element's turn to be discovered one ticket at a time.

### 2. "Fe cannot take Al's mean-deck route (τ500-vs-Rosseland, 127/6400 columns)" — a category error

The 127/6400 measurement is real (RYA-1013) and it is about **raw STAGGER cube columns**,
which are trimmed on the Rosseland scale. **The mean-deck route never touches a cube
column.** It uses `sun_avg3d_stagger.mod` — the spatially and temporally averaged ⟨3D⟩
model, *already averaged on surfaces of constant τ500*, header `TAU5000 SCALE`, 101 depths,
**log τ₅₀₀ −5.000 … +5.000**. That is four decades below the τ ≤ 1e-4 threshold only 127
raw columns reached, and Fe I optical lines form near log τ₅₀₀ ≈ −0.5 … −3, well inside it.

The atmosphere is also **element-agnostic** — it is the same file Al uses, and nothing in it
is Fe-specific. So the stated blocker does not apply to this route at all. It *would* apply
to a full-3D column-by-column build, which is the route Step 0 says we do not need to scope.

Measured end-to-end against the real deck: the solar node returns `ndep=101`, `nlev=607`,
log τ −5.000…+5.000, 100% finite, zero all-zero rows, median *b* → 1.0000 at the deepest
layer (the LTE limit). **The route runs.**

### 3. "The register currently carries Fe FULL_3D = BROKEN" — STALE

The ticket asked for this to be re-checked. Re-measured on `main` today, over the committed
RYA-817 in-domain pool:

| path | Fe I finite | Fe II finite | median Δ |
| --- | --- | --- | --- |
| **pinned axis** (what the production script uses) | **114 / 114** | **7 / 7** | +0.0389 / +0.0534 |
| unpinned per-line axis (the archived RYA-207 reading) | 33 / 114 | 1 / 7 | +0.0406 / +0.0539 |

`model_availability_matrix` says *"the MLP returns NaN for EVERY in-domain line … 114
in-domain → n=0"*. That is **overstated in both halves**: the loss is 71% on Fe I, not 100%,
and it belongs to a path `scripts/rya817_run_3dnlte_bands.py` does not take — it passes
`afe3n_axis=afe_star`, and its reactivation control passes 4/4 against Amarsi+2022 Table 6.

**Correct state: `HAVE_RUNS`, optical band only, axis railed at the grid's A(Fe)=7.5
ceiling with the rail recorded and its sensitivity measured at 0.0066 dex.** The cell is
updated; the RYA-923 fix it was waiting for is already in the code.

---

## 🔴 The real Fe-specific blocker, and it is not the one the ticket names

The Fe ⟨3D⟩ aux table has a **zeroed metallicity column on exactly the seven rows a solar
run selects.** The deck's seven Teff = 5777 members span the full metallicity axis —
`p5777g44m00 / m05 / m10 / m20 / m30 / m40 / p05` — and the shipped `[Fe/H]` column reads
**+0.00 for all seven**. The other 182 rows agree with their own names exactly, 182/182, so
the **column** is wrong and the **name** is right. The Al ⟨3D⟩ aux is the positive control:
**0 disagreements in 6345 rows**, so the referee cannot move Al.

Left alone the consequence is not a crash, it is a **wrong star**: all seven tie at the
solar node, the tie-break is on A(X) — identical at 7.50 across all seven — and the first
wins, which is `p5777g44m10`, **[Fe/H] = −1.0**. The true solar record is sixth in file
order.

**And `[Fe/H]` is not the only field that was written wrong.** The deck's own relation
A(X) = 7.50 + [Fe/H] is *exact* on all 183 clean rows and violated on exactly those six,
which ship A(X) = 7.50 at [Fe/H] = −4.0 … +0.5. The name encodes the metallicity so it can
referee it; **nothing referees A(X)**. So an overridden row is marked SUSPECT and
**refused** — the override exists to get the bad rows out of the candidate set, not to
repair them. Reconstructing A(X) from the relation would be inventing vendor data.

The Sun is unaffected: `p5777g44m00` is the one row of the seven the column happens to get
right, so it is never overridden and never suspect.

### 🔴 …and the canonical aux is the one we must NOT use

TSFitPy's downloader points at `auxData_Fe_STAGGERmean3D_May-21-2021_marcs_names.txt`, and
`Al@mean3D` is registered against its `_marcs_names` sibling. **For Fe that file is
unrecoverable.** The vendor's `convert_3d_grid_to_marcs_names.py` builds the new name *from
the `[Fe/H]` column*, so it propagated the zeroing into the name: all seven come out as the
byte-identical string `p5777_g+4.4_m0.0_t02_st_z+0.00_a+0.00_…`. Name and column now agree
and both are wrong for six of them — there is no signal left to referee with.

Measured end-to-end: through `_marcs_names` the solar node is **refused**; through the plain
aux it returns `p5777g44m00`.

> **⇒ `Fe@mean3D` must be registered against the PLAIN aux — the opposite of `Al@mean3D`.**

---

## A latent defect this turned up in the existing reader

`_node_from_model_name` read the MARCS-style alias's `_m0.0_` field as the metallicity. In
that convention **`m` is the stellar MASS**; the metallicity is `z-4.00` further along. The
two readings agree **only at the solar node** (mass 0.0, z+0.00), which is why the existing
test passed — it pinned the solar *example*, not the invariant
([[feedback_pin_the_invariant_not_the_example]]). On the Al ⟨3D⟩ aux the two readings
disagree on **123 node rows**.

It was harmless on `main` because the function was only ever applied to the *record's* name,
which is always STAGGER-form. It stops being harmless the moment it is applied to an aux id,
which is exactly what the Fe referee needs to do. Fixed, with the invariant pinned.

Also fixed: the parser did not know STAGGER's **short-Teff** form (`p50g25m40` = Teff 5000).
**182 of the Fe deck's 189 rows use it**, and `read_deck_node` refuses a record it cannot
identify — so Fe ⟨3D⟩ would have been unusable at every node except the solar one.

---

## The route, specified

| | |
| --- | --- |
| deck | `NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin` (92,945,908 B, md5 `d3fa0517…`) |
| aux | `auxData_Fe_STAGGERmean3D_May-21-2021.txt` — **plain, not `_marcs_names`** |
| atom | `atom.fe607a` (already on Sirius for the 1D deck) |
| atmosphere | `data/atmospheres/stagger_avg3d_rya442/sun_avg3d_stagger.mod` (in-repo) |
| reader | `gerber_nlte.read_deck_node` — direct, no vendor binary (RYA-821) |
| layout | 500 + 4 + 4 + 101×8 + 101×607×8 = **491,772 B**; aux pointers differ by exactly that (1001 → 492773), so 1-based after a 1000-byte header |
| abundance axis | **none** — A(X) = 7.50 + [Fe/H], one value per node |
| node coords | STAGGER 5777 / 4.44, *not* MARCS 5750 / 4.5 |
| coverage | Teff 4000–7000, logg 1.5–5.0, [Fe/H] −4 … +0.5, 189 nodes |

Because Fe has no abundance axis, **the v118 machinery Al needed does not apply**: the
pre-v118 hoisting of departures out of the χ² loop was always correct for Fe and stays
correct. Fe ⟨3D⟩ is the *simpler* of the two routes, not the harder one.

### What is still owed

1. ~~**Stage the deck on Sirius.**~~ **DONE 2026-08-24.** Both files are in
   `/mnt/codex-ext/codex-grids/nlte/gerber_ts/` (the real path, resolved with `readlink`
   before measuring free space — RYA-800), md5-verified in place against the Mac
   measurement, 109 GB free after. `atom.fe607a` was already there and its md5
   `d08dc8232ed68eec65f9bb6631e82ea8` matches, with `7.50  55.85` on line 2.
   **The route was then verified END-TO-END on the staged bytes**: six aux overrides
   detected, no abundance axis, `deck_abundance` = 7.50 with the record cross-examined by
   the aux, solar node `p5777g44m00` at ndep 101 / nlev 607 / log τ −5…+5 / A(X) = 7.50,
   all finite, median *b* → 1.0000, departures sha256 `62556a5664d0b8f7` — **identical to
   the Mac read** — and the suspect rows refused. ⚠️ Only the **PLAIN** aux was staged.
2. **Register `Fe@mean3D`** in `gerber_nlte.DECKS`, against the plain aux. *Not done here on
   purpose*: the matrix derives "WIRED" from `DECKS` (v117), so registering a deck whose
   bytes are not on Sirius would make the status surface claim a capability we cannot run —
   the exact lie v117 removed. There is a test pinning that it stays unregistered until then.
3. ~~**Resolve the deck abundance.**~~ **RESOLVED — A(Fe) = 7.50. See below.**
4. **The band product.** Same step still owed for Al (`derive_band_products` has no
   `gerber-mean3d` deck choice and loads a MARCS.GES atmosphere). Fe inherits that work; it
   is not Fe-specific and should not be re-derived per element.

## 🔴 The deck abundance, resolved: **A(Fe) = 7.50** — and 7.46 was our own input

The record said 7.46 and the grid said 7.50. Resolved **by provenance**, as the re-scope
required — not by picking the more convenient number.

### Where 7.46 came from: a closed loop with no external referee

The evidence for 7.46 was one bsyn message, read as the grid speaking:

> `Bsyn: NLTE departure coeff calculated for abundance = 7.46 while it is 7.50 here`

Traced through the vendor source, it is **us** speaking:

| step | file | what happens |
| --- | --- | --- |
| 1 | `gerber_nlte.for_node` | `a_deck = deck_abundance(element)` → reads **7.46 from our own prov.json** |
| 2 | `_interpolate` stdin | passed as the interpolator's `abu_ref` |
| 3 | `interpol_modeles_nlte.f:206` | `read(*,*) abu_ref` — **from stdin** |
| 4 | `interpol_modeles_nlte.f:761` | `write(27,1971) abu_ref` — **verbatim into the departure file** |
| 5 | `read_departure.f` | loaded back as `abundance_nlte` |
| 6 | `bsyn.f:988` | printed at us |
| 7 | `Fe_gerber2023.prov.json` | the printout recorded as *"a property of the deck, not a fitted value"* |

`abu_ref` is a **label the caller supplies**. The deck never asserts it. And this module's
own docstring already carried the disproof, measured under an earlier ticket: running the
interpolator at A = 7.36 / 7.46 / 7.56 gives departure files **byte-identical except for
the stamp itself** (444836 bytes each, 1 differing line of 132).

> ⚠️ **7.46 is the Asplund solar A(Fe).** It is exactly the number that looks right on
> arrival, which is why nobody asked where it came from. Ryan flagged this as a
> reference-proximity smell before the trace confirmed it.

### Where 7.50 comes from: the grid, three independent ways

1. **The model atom declares it.** `atom.fe607a` line 2 is `7.50  55.85` — A(Fe) and the
   atomic mass. Fetched and read; its md5 `d08dc8232ed68eec65f9bb6631e82ea8` **matches the
   one already pinned in our own prov.json**, so this is byte-identical to the copy staged
   on Sirius. This is the abundance the departure solve itself used.
2. **Both aux tables encode it.** A(X) = 7.50 + [Fe/H] **exactly** — 15,229 MARCS rows and
   183 clean ⟨3D⟩ rows. (The six exceptions are the zeroed rows above, which are refused.)
3. **Turbospectrum itself hardcodes it** as solar iron: `metal = abund(15) - 7.50`,
   `interpol_modeles_nlte.f:1177`.

### Why the two paths disagreed at all

`read_deck_node` (the ⟨3D⟩ direct route) returns the **aux's own** A(X), so it self-declares
**7.50 whatever abundance it is asked for** — measured across requests at 7.46 / 7.50 / 7.60
/ 6.90, all returning `p5777g44m00` with identical departures (sha256 `62556a5664d0b8f7`).
`read_departure_file` (the interpolator route) returns the **stamp**, i.e. the caller's own
input. One path reports the grid; the other reports us. That difference *was* the ambiguity.

### What changes, and what deliberately does not

**Nothing measured moves.** `deck_abundance()` feeds only (a) the interpolator's `abu_ref`,
a stamp the departures do not follow, and (b) a printed `A_deck=` provenance string in
`derive_band_products` and the RYA-798 control. What bsyn finally compares is written by
`as_ispec_tuple` from the **trial** abundance, so the fit is unaffected. This is a **record
correction**, which is precisely why it was safe to make without a physics review.

**bsyn's STOP stays exactly as it is.** Refusing an ambiguous abundance is the correct
failure mode. The ambiguity was resolved; the guard was not downgraded to a warning.

**And the loop is now uncloseable.** `deck_abundance()` cross-examines the provenance record
against the deck's own aux and **raises** on disagreement. A value that came from us can no
longer be handed back to the vendor binary as though it had come from the vendor. It falls
back to the record only when the aux is genuinely unreachable (Sirius-only), where the
record is all there is.

---

### Not scoped, deliberately

The ticket's build branch — full-3D via our own solver, or a τ_Rosseland↔τ500 corrected
route — was conditioned on *"only if Step 0 comes back empty"*. It did not. Note also that
RYA-1008 already established that **no full-3D NLTE stellar RT code is downloadable
anywhere**, so that branch was never a fetch-and-build; it is an access ask. Nothing here
changes that, and nothing here needs it.

---

## Two more routes worth a ticket, not worth this one

- **`iron_abundancecorr.tar.gz`** (Amarsi et al. 2016, 182 MB, live, HTTP 200) — 1D **and
  ⟨3D⟩** NLTE corrections for Fe I and Fe II, from a different group than Gerber. Its value
  is as an **independent referee** on the Gerber ⟨3D⟩ route, which is the validation-triangle
  third leg the umbrella wants. Unfetched.
- **`cofe_tools.tar.gz`** (Amarsi et al. 2019, 114 MB, live) — the C/O release we already
  hold as `data/nlte_grids/amarsi2019_cno/` also carries a **3D-NLTE Fe II** leg.
  `THREED_HOLDINGS` maps that release to C/N/O and not to Fe, so we may already hold Fe II
  3D-NLTE on disk and not know it.

## Firewall (RYA-161)

Every route above is described by **physics and access** — what structure the departures
were solved on, what depths it reaches, what box it covers, whether we can get it. No route
was preferred, ranked or rejected on where its answer sits relative to 7.46, 7.50 or any
other reference value, and no abundance was derived in this ticket.
