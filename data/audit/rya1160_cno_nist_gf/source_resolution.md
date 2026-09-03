# CNO gf source resolution — what each indicator's gf actually traces to

Read from the papers themselves, not inferred.

| indicator | gf source | status |
| --- | --- | --- |
| C i (14 of 16 permitted) | **Li, Amarsi, Papoulia, Ekman & Jönsson 2021**, MNRAS 502, 3780 — AGSS21: *"new g f-values from large-scale atomic structure calculations by Li et al. (2021)"* | ✅ acquired, C I–IV, 4,321 transitions |
| [C i] 872.7 nm | Amarsi et al. 2019 Table 1 | ✅ acquired |
| N i (5 lines) | **Tachiev & Froese Fischer 2002**, A&A 385, 716 — AGSS21: *"taken from Tachiev & Froese Fischer (2002) … reliable at the 0.03 dex level"* | ✅ acquired + A→gf converted |
| NH, CN | Brooke et al. 2014a,b, 2015 | ✅ acquired |
| CN/NH dissociation energies | Barklem & Collet 2016 | ✅ acquired (De only; table6/7 absent — RYA-1147) |
| [O i] 630.0 / 636.4 blend gf (Ni I) | Johansson et al. 2003 + log ε(Ni)=6.20 | ✅ in bibliography |
| **O i (777, 616, 844, 926)** | **NOT NAMED** | 🔴 **OPEN** |
| **[O i] 630.0 / 636.4 (the O transition itself)** | **NOT NAMED**; NIST cites TP `T4539,T5081`, repo bibliography points at `storey_zeippen2000` | 🔴 **OPEN** |

## The oxygen gap, checked in both places

**AGSS21** names Tachiev & Froese Fischer for **N i only**. It states no gf source for O i.
**Amarsi et al. 2019** (A&A 630, A104), which AGSS21 cites for its O i lines, defines the
*line set* — *"the permitted O I 616 and 777 nm multiplets and the forbidden, low-excitation
[O I] 630.0 and 636.4 nm lines … extended to include the permitted O I 844 and 926 nm
multiplets"* — but names no gf source either.

So the oxygen gf provenance is genuinely unstated in both papers. It is not something this
audit failed to find; it is not published in the chain we can see.

**Why it matters now.** Tachiev's O-like calculation covers O I and disagrees with NIST
systematically: **−0.016 dex across the whole 777 triplet** (Tachiev A = 3.556e7 vs NIST
3.69e7) and −0.005 dex across the 8446 triplet. The 777 triplet is AGSS21's dominant oxygen
indicator, so which source it used moves A(O) by ~0.016 dex. Until that is resolved, the
O i rung is **not** Asplund-grade — it is "one of two candidate sources, undetermined".

**[O I] 6300/6363 are not in Tachiev at all**: table5 contains zero 2p⁴→2p⁴ rows, so
intra-ground-configuration forbidden transitions are outside its scope. Their gf traces
elsewhere; Storey & Zeippen 2000 is the standing candidate and is not acquired.

## C I 17 vs 14 — settled from the source's own prose

Amarsi et al. 2019 §2.4: *"The C I lines are listed in the first 17 rows of Table 1 …
**The 16 permitted lines** are all of high excitation potential, χexc ≳ 7.5 eV … **the
forbidden, weak (log gf = −8.165), low-excitation (χexc = 1.264 eV) [C I] 872.7 nm line is
also included**."*

So **17 = 16 permitted + 1 forbidden**, confirming the split derived independently from the
census. AGSS21 then uses **14** of the 16 permitted for its `C i` indicator and the 872.7 nm
line for `[C i]`. Ryan's call: keep all 17, flag the 2 permitted lines outside AGSS21's
adopted set. Their identity is not determinable from either paper.

⚠️ One further nuance for RYA-1159: Amarsi 2019 states its corrections *"do not take into
account the Ni I line that blends the [O I] 630.0 nm line"*, whereas AGSS21 **does** model
that blend with the Johansson 2003 gf. The two are not interchangeable on 6300.
