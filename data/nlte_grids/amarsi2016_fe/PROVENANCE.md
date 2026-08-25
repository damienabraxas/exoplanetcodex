# Amarsi+2016 ⟨3D⟩ Fe — the scale-matched external anchor (RYA-1042)

**Citation, pulled at acquisition from the archive's own `readme.txt` — not from memory:**

> Amarsi, Lind, Asplund, Barklem & Collet 2016, **MNRAS 463, 1518**;
> DOI `10.1093/mnras/stw2077`
> *"Non-LTE line formation of Fe in late-type stars – III. 3D non-LTE analysis of
> metal-poor stars"*

DOI **verified at Crossref** on acquisition (2026-08-25): title, author list, volume and
pages all resolve as above.

- Source: `https://www.astro.uu.se/~amarsi/iron_abundancecorr.tar.gz`
- md5 `64c7f482864851f14f1d7dc1ff881d07`, 182,158,708 B
- Staged: `/srv/codex/grids-overflow/nlte/amarsi2016_fe/` (Sirius only — `nmtd_lmtd.txt`
  alone is **748 MB**). Only the solar-node slice is committed here.

## 🔴 Which of the three grids, and why it matters

The archive ships three files and its readme says exactly what each is:

| file | quantity | verdict |
| --- | --- | --- |
| `nmarcs_lmarcs.txt` | 1D non-LTE − 1D LTE, on 1D MARCS | not our scale |
| `nmtd_lmarcs.txt` | ⟨3D⟩ non-LTE − **1D** LTE | ❌ **the trap** |
| **`nmtd_lmtd.txt`** | **⟨3D⟩ non-LTE − ⟨3D⟩ LTE**, on averaged-3D STAGGER | ✅ **this one** |

Our RYA-1040 product is precisely that difference — ⟨3D⟩-NLTE minus ⟨3D⟩-LTE on the STAGGER
averaged atmosphere — so `nmtd_lmtd` is the **same physical object computed by an
independent group**. The readme says so itself: *"Only use this if you are correcting LTE
results based on ⟨3D⟩ or 3D model atmospheres."*

⚠️ **`nmtd_lmarcs` would have been the trap.** It is also ⟨3D⟩ non-LTE and also
Amarsi+2016, but its comparand is **1D LTE**, so it carries the 1D→mean-3D **atmosphere**
shift inside it. Gating our differential against it would fold an atmosphere difference
into an NLTE test — the RYA-542 confound the paired-product design exists to avoid.

⚠️ **And it is not Amarsi+2022.** RYA-817's MLP is **full 3D**; gating a mean-⟨3D⟩ deck
against it measures mean-vs-full, not deck-vs-reference.

## The audit (the file is self-documenting; nothing below is inferred)

Header, verbatim:

```
T_eff/K, log10(g/cm s^-2), log10(eps_Fe), v_turb/km s^-1, Species,
lambda_{Air,centre}/nm, E_low, log10(gf), lambda_{Air,min}/nm,
lambda_{Air,max}/nm, Clean, Abundance correction
```

| item | finding |
| --- | --- |
| **wavelength medium** | **AIR**, stated in the column name. Our line lists are air too — checked, not assumed: a medium mix-up is ~1.4 Å at 5000 Å and would mismatch every line silently. |
| **units** | wavelength **nm** (our band arguments are Å); correction in **dex** |
| **atmosphere / geometry** | averaged-3D STAGGER, per the readme; this file's LTE comparand is on the *same* atmosphere |
| **abundance axis** | **absolute `log₁₀(ε_Fe)`**, not [Fe/H] — 2.50…8.00 in 0.25 steps |
| **grid** | Teff 4000–7000 / 250 K · logg 1.50–5.00 / 0.5 · vturb 0.75, 1.50, 3.00 · Fe1 + Fe2 |
| **`Clean`** | the authors' unblended flag — everything downstream filters on it, and the counts are reported so the filter is visible |

⚠️ **The title says "metal-poor stars" and the grid still covers solar.** Paper III is a
metal-poor analysis; the released grid is not restricted to it. Checked rather than
inferred from the title — the solar node exists at **(5750, 4.50, 7.50)**, and **7.50 is
the Gerber deck's own A(Fe)** (RYA-1035), so the anchor is *read* at our abundance rather
than interpolated to it.

⚠️ **THERE IS NO `vturb = 1.0` NODE.** Solar vturb is 1.0; the grid brackets it at 0.75 and
1.50. The anchor is therefore a **bracket, not a point**, and both are emitted — an
interpolated single number would hide that the grid was never asked.

⚠️ **A SENTINEL LIVES IN THE FILE.** At vturb 3.00 one Fe I row is exactly `-4.0000`, a
floor value rather than a correction. The extractor flags and excludes rows at exactly
−4.0 and reports the count; averaging it in would drag a median.

## The anchor, at the solar node (Clean = yes, sentinels excluded)

| species · vturb | n | median | range |
| --- | --- | --- | --- |
| Fe I · 0.75 | 540 | **−0.0299** | −0.4367 … +0.1834 |
| Fe I · 1.50 | 540 | **−0.0411** | −0.3488 … +0.2598 |
| Fe I · 3.00 | 539 | −0.0597 | −0.3025 … +0.2822 (1 sentinel excluded) |
| Fe II · 0.75 | 75 | −0.0010 | −0.0176 … +0.0014 |
| Fe II · 1.50 | 75 | −0.0008 | −0.0245 … +0.0026 |
| Fe II · 3.00 | 75 | +0.0007 | −0.0427 … +0.0042 |

⇒ **Solar Fe I at vturb 1.0 brackets to ≈ −0.030 … −0.041 dex; Fe II is ≈ 0.000.** Fe II
sitting at zero while Fe I does not is itself a sanity check — Fe II is the majority
ionisation stage and close to LTE, which is what the literature expects.

Over our VIS band specifically (420–691 nm, Clean = yes, vturb 1.50): **n = 258,
median −0.0376**, range −0.0928 … +0.0024.

## Firewall (RYA-161 / RYA-1035)

This anchor is **external** and was computed by another group with a different code, a
different model atom and a different solver. It is never derived from the machinery it
validates — which is the closed-loop trap RYA-1035 found in the deck's own abundance
record (7.46 was our stdin round-tripping back through bsyn's error message).
