# RYA-969 — The differential bridge

**DESIGN SPEC — for sign-off. No implementation.**
**Scope: the DIFFERENTIAL layer.** RYA-968 is the absolute half.
Firewall: RYA-161. Feeds: RYA-116 (α Cen A). Blocked by: RYA-975. Doc sync: RYA-179.

---

## 0. Summary — the method is sound, the ladder is short, and two things block it

**The method is the field's and it is right.** Line-by-line differential analysis against a
reference star cancels the gf uncertainty — same line, same unknown gf, subtracts out — and
reaches 0.01–0.02 dex (Meléndez+2009; Bedell+2014; Nissen+2015). For non-twins, Jofré+2015 bridges
in two steps: target → same-type benchmark → Sun.

**Our ladder is shorter than the ticket assumes, and that is good news.** Measured from
`stars.yaml`, **no hop we currently need comes near Jofré's failure thresholds**, so most stars
can bridge *directly* to the Sun and the two-step machinery is reserved for the rung we do not
yet have (the M dwarf, RYA-970).

**But two blockers are structural, and one is in our own catalogue:**

* 🔴 **There is no second rung.** Every measured line product in the project is **solar**. There is
  not one measured pool for α Cen A/B, τ Ceti, ε Eri or 55 Cnc, so the line-sharing gate cannot be
  computed for a single hop — no hop has data at both ends.
* 🔴 **The conditioning is asymmetric.** `solar_harps_molecfit_corrected` is `telluric_applied =
  applied`; `alpha_cen_a_harps` and `alpha_cen_b_harps` are **`unknown`**. A differential across
  asymmetric conditioning does not cancel systematics — **it injects them**. RYA-975 is therefore a
  hard prerequisite for hop 1, not an adjacent improvement.

---

## 1. Why differential is legitimate, stated precisely

The gf term cancels because the *same physical line* is measured in both stars and the same
unknown `log gf` enters both abundances with opposite sign in the difference. This is not a
statistical trick and not tuning: it is an algebraic cancellation of a shared unknown.

**What it does NOT do is give you a scale.** A differential yields `[X/H]_target−ref`. Converting
that to an absolute `A(X)` requires the reference star's own absolute abundance, which comes from
the absolute layer (RYA-968) or a GBS reference value — **never from fitting the target**.

🔴 **The two layers compose in exactly one direction:**

> **Precision comes from the differential. Scale comes from the absolute. The scale error never
> shrinks, however precise the differential is.**

RYA-968 §3.2 puts a hard number on the second half: our absolute zero-point is capped at
**≈0.059 dex** by a seven-line laboratory anchor. So a differential quoted at 0.01 dex is a real
*relative* precision, and the *absolute* statement built on it still carries 0.059. Reporting one
without the other is the circularity RYA-968 was written to prevent, arriving from the other side.

---

## 2. The ladder, measured rather than assumed

Stellar parameters from `config/stars.yaml` (τ Ceti and ε Eri added by RYA-957):

| star | Teff | log g | [Fe/H] |
|---|---|---|---|
| Procyon | 6554 | 4.00 | +0.03 |
| α Cen A | 5792 | 4.30 | +0.20 |
| **Sun** | **5772** | **4.44** | **0.00** |
| τ Ceti | 5414 | 4.49 | −0.49 |
| α Cen B | 5231 | 4.53 | +0.20 |
| 55 Cnc A | 5196 | 4.45 | +0.31 |
| ε Eri | 5076 | 4.61 | −0.09 |

### 2.1 Every hop from the Sun is short

| hop | ΔTeff | Δlog g | Δ[Fe/H] |
|---|---|---|---|
| Sun → α Cen A | **+20 K** | −0.14 | +0.20 |
| Sun → τ Ceti | −358 K | +0.05 | −0.49 |
| Sun → α Cen B | −541 K | +0.09 | +0.20 |
| Sun → 55 Cnc A | −576 K | +0.01 | +0.31 |
| Sun → ε Eri | −696 K | +0.17 | −0.09 |
| Sun → Procyon | +782 K | −0.44 | +0.03 |

Jofré's documented failures are **Δ[Fe/H] ≈ 2.5** (Sun → HD 140283) and **ΔTeff ≈ 2000 K**
(Sun → M giant). Our largest are **782 K** and **0.49 dex** — factors of ~2.6 and ~5 inside the
known failure points.

⇒ **For the stars we hold, the two-step bridge is not required.** Sun → α Cen A at ΔTeff = 20 K is
very nearly a twin pairing, and the rest are ordinary short hops. The two-step machinery earns its
place at the **M-dwarf rung** (RYA-970), where the gap really is large.

### 2.2 🔴 Nearest neighbour is NOT the right bridge reference

Scaled nearest neighbours (Teff/100, log g/0.1, [Fe/H]/0.1):

| star | nearest | second |
|---|---|---|
| α Cen A | **Sun** (2.4) | α Cen B (6.1) |
| α Cen B | **55 Cnc A** (1.4) | ε Eri (3.4) |
| 55 Cnc A | **α Cen B** (1.4) | ε Eri (4.5) |
| ε Eri | α Cen B (3.4) | 55 Cnc A (4.5) |
| τ Ceti | ε Eri (5.4) | Sun (6.1) |
| Procyon | α Cen A (8.4) | Sun (9.0) |

α Cen B and 55 Cnc A are each other's nearest neighbours, and neither's nearest is α Cen A. A
naive nearest-neighbour chain would therefore route α Cen B through 55 Cnc A — **adding a hop for
no gain**, because α Cen B can reach the Sun directly inside the failure limits.

🔴 **Each hop adds its own uncertainty in quadrature, so the objective is the shortest PATH TO THE
SUN, not the nearest neighbour.** The config must express a *chain to the Sun*, and the selection
rule is: **take the direct hop whenever it passes the line-sharing gate; insert an intermediate
rung only when it does not.** Nearest-neighbour is the fallback for choosing that intermediate,
not the primary rule.

⚠️ **Procyon is the genuine outlier** — an F dwarf, ΔTeff +782 K and Δlog g −0.44 from the Sun, a
different atmospheric regime. It is the one current star where the direct hop is worth testing
rather than assumed.

---

## 3. The line-sharing gate

**Per hop**, compute the set of lines that are measurable and unsaturated **in both stars**:

* present in both measured pools, matched on wavelength **and** excitation potential (RYA-780);
* unsaturated in both — REW inside the linear-COG window, derived per star (RYA-968 §3.1 showed a
  borrowed REW window controls nothing);
* passing each star's own physical checks (RYA-959 width ceilings; the RYA-968 stage 1–4 tree).

**Below a declared threshold of shared lines, the hop is FLAGGED TOO LARGE and fails loudly.** It
is never silently bridged. The remedy is an intermediate rung, which is exactly how the M-dwarf
gap (RYA-970) will present itself.

🔴 **This gate is currently uncomputable for every hop** (§0), because it needs measured pools at
both ends and only the solar end exists. **Standing up one non-solar measured pool is the single
unblock that makes this design testable** — and α Cen A is both the shortest hop and the JWST
deadline target (RYA-116).

---

## 4. Products — two per target, always both

| product | source | carries |
|---|---|---|
| **Differential `[X/H]`** | this layer, line-by-line vs the bridge reference | high relative precision; gf cancelled |
| **Absolute `A(X)`** | RYA-968 | the scale, its provenance, and the zero-point cap |

Both are reported. The differential is chained to the Sun through the config chain, accumulating
per-hop uncertainty in quadrature; the absolute is not chained at all.

**Validation:** at each benchmark rung, compare against the GBS reference abundances (Jofré+2015).
⚠️ That comparison is a **check reported as a finding, never an input** — the moment a bridge is
adjusted to match GBS, the bridge is transporting an assumption instead of a difference.

---

## 5. Firewall (RYA-161)

Differential analysis is legitimate **because gf cancels**, not because it agrees with anything.
Feltzing+2009 is the literature basis for the asymmetry RYA-161 encodes: astrophysical
(solar-fitted) gf is *legitimate here*, where it cancels, and *forbidden* in the absolute layer,
where it is circular.

**Structural requirements:**

* **The reference star's abundance is an INPUT, never a fitted parameter.** It comes from the
  absolute layer or a GBS reference value. Enforced by signature: the bridge function receives the
  reference abundance and the two line sets, and has no path to alter the reference.
* **The bridge transports a difference.** A hop's output is `Δ[X/H]` plus its shared-line set; the
  target's absolute value is only ever assembled at the end, from the chain plus the anchor.
* **No hop may be re-pointed to improve agreement.** The chain is chosen by the line-sharing gate
  and parameter distance, both declared before any abundance is computed — the RYA-968 F3 rule.
* **GBS comparison is a report, not a term.**

---

## 6. 🔴 Blockers, in the order they bite

| | blocker |
|---|---|
| 1 | **No second rung exists.** Every measured product is solar; the gate has nothing to compare. **Unblock: one non-solar measured pool — α Cen A first** (shortest hop, and RYA-116's deadline). |
| 2 | **Telluric asymmetry — RYA-975.** `solar_harps_molecfit_corrected` = `applied` vs `alpha_cen_a_harps` = **`unknown`**. Differencing across that does not cancel systematics, it injects them, and `unknown` is worse than `not-applied` because the sign is not even predictable. **Hard prerequisite for hop 1.** |
| 3 | **α Cen A/B identity — RYA-971 / RYA-423.** The first bridge reference must be the right star. RYA-952 established that astrometry alone cannot separate a 4–8″ pair, and the branch assignment is still unadjudicated. A bridge built on a misassigned reference is wrong in a way no downstream check would catch. |
| 4 | **Shared-line threshold undeclared.** Must be fixed before the first run (RYA-968 F3). |

---

## 7. Element-agnostic, and what RYA-967 buys

Nothing here is Fe-specific: the cancellation is a property of measuring the same line in two
stars. **RYA-967 compounds it** — every additional line the synth route recovers is another
potentially shared line, which lengthens the shared set, shortens the *effective* gap, and lets
more elements reach the bottom of the ladder. The same graded backbone does triple duty: absolute
anchor (RYA-968), empirical grade (RYA-968), differential bridge (here).

---

## 8. Open for sign-off

1. **Accept that the two-step bridge is not needed for our current stars?** The measured geometry
   says direct-to-Sun passes for all six; the two-step exists for the M rung.
2. **Confirm shortest-path-to-Sun over nearest-neighbour** as the chain rule (§2.2).
3. **RYA-975 as a hard blocker on hop 1** — does telluric symmetry gate the bridge, or do we
   proceed and carry the asymmetry as a stated systematic?
4. **Which non-solar pool is stood up first?** α Cen A is the shortest hop and the JWST target.
5. **The shared-line threshold** — who declares it, and where?

---

## 9. References — verified 2026-08-22

| ref | verified as |
|---|---|
| Jofré+2015 | *Gaia FGK benchmark stars: abundances of alpha and iron-peak elements*, A&A 582, A81 — `10.1051/0004-6361/201526604` (arXiv:1507.00027). The two-step bridge. |
| Meléndez+2009 | *The peculiar solar composition and its possible relation to planet formation*, **ApJ 704, L66** — `10.1088/0004-637X/704/1/L66` (arXiv:0909.2299). 🔴 **Vendored as `melendez2009twins`, NOT `melendez2009`** — that key was already taken by **Meléndez & Barbuy 2009** (A&A 497, 611), the *firewalled* Fe II gf paper (RYA-161/852). Different paper, different authors, opposite role: one is the differential precision anchor, the other is a source RYA-161 exists to exclude. ⚠️ The ticket cites "aa13038-09", an A&A manuscript number matching **neither**; if the A&A companion (Ramírez, Meléndez & Asplund 2009, A&A 508, L17) was intended, confirm before relying on it. |
| Bedell+2014 | *Stellar chemical abundances: in pursuit of the highest achievable precision*, ApJ 795, 23 — `10.1088/0004-637X/795/1/23`. |
| Nissen 2015 | *High-precision abundances of elements in solar twin stars*, A&A 579, A52 — `10.1051/0004-6361/201526269` (arXiv:1504.07598). |
| Feltzing+2009 | already vendered by RYA-968 — the differential/absolute firewall basis. |

## 10. Evidence index

§2 geometry from `config/stars.yaml` via `config.constants.STAR_PARAMS`. §0/§6 blocker 1 from a
sweep of every `data/results/band_products/*_lines.csv` and `data/measured/band_ew/*` across all
Sirius worktrees — all solar. §6 blocker 2 from
`data/catalog/holdings_manifest_registry.csv` `telluric_applied`.
