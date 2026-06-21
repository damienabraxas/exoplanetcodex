# RYA-405 — Solar Fe baseline reconcile: Step-0 git-state audit → STOP (no integration gap)

**Verdict: the integration-gap hypothesis is REFUTED.** The entire solved Fe II stack —
RYA-290, **305**, **336**, 341→342, 344, 346, 347, **352** — is **already on `main`** and
**actively firing** on the solar gate run. Main's gate Fe II = 7.700 is the *fully-integrated*
305+352 result, not a fix sitting behind on a branch. Per the ticket's CRITICAL branch point
("if 305/352 ARE on main and the gate still fails → STOP; fresh RCA, do not re-land"), this
session **re-lands nothing**. The "failing gate" is a **path/arbiter mismatch**, not a missing
fix — diagnosed below and routed for a fresh RCA.

## Step 0 — git-state audit table (CONFIRMED on `main` @ `166cb4b`)

| Fix | Ticket | On main? | Commit(s) | Documented solved-value |
|-----|--------|----------|-----------|-------------------------|
| Solar anchor Fe I/II imbalance audit | RYA-290 | **YES** | `a5692ad`, `1633326` | diagnostic only; closed by 305 |
| Fe II EW-path repair (triage + synth anchor) | RYA-305 | **YES** | merge `31fb531`, `aaa6b51`, `ad8498e` | **synth** Fe I−Fe II = −0.007 @ Fe II **7.486**; EW triage (clean/recover/drop) |
| Scale-aware solar Fe gate | RYA-336 | **YES** (ticket said Backlog — STALE) | `9ccf064` | scale-robust verdict + abs-A(Fe) diagnostic window [7.44,7.58] |
| Synth-ceiling decouple + measured line set | RYA-342 | **YES** | `afe858f`, `a02dc2a` | unblocks the synth arbiter |
| Reconcile synth Fe II (crowned) | RYA-341 | (analysis-only, not merged; its FIX = 342) | — | **synth** Fe I−Fe II = **−0.015** (Fe II **7.500**); DECISION 2 ratified |
| Fe II blend-recover → STOP | RYA-344 | **YES** | `418e021` | 0/5 recoverable |
| Fe II core-residual diagnostic | RYA-346 | **YES** | `97e727d` | 5/5 no-absorber |
| Fe II atomic-data audit | RYA-347 | **YES** | `3706a66` | 5/5 RETIRE (gf + ABO-damping both ruled out) |
| Fe II EW-quality cull (supersedes 347-basis) | RYA-352 | **YES** | PR #21 `1558c19`, `74afeb4` | **clean EW pool Fe II = 7.696 (HIGH)**; 7.466 was an ARTIFACT; *EW path is not the arbiter* |

**Runtime proof the stack is live (verbatim, solar run on `main`):**
```
[Fe II EW-quality cull RYA-352] dropped 156 of 227 Fe II lines (ceiling 100mÅ, err/EW>0.5)
Fe II triage (RYA-305): 7 clean (EW kept), 1 blended→synthesis-recover (EW-excluded), 6 dropped/quarantined
```
→ 305 + 352 are not merely present, they are *executing* and produce the 7-line clean Fe II
pool whose EW median is 7.700. Nothing is un-landed.

## The premise is refuted — and why

The ticket hypothesised "main reads +0.21 above 305's solved 7.486 → main running BEHIND the
solved Fe II stack." The +0.21 is real (7.700 − 7.486 = 0.214) but it is **RYA-352's deliberate
re-basing, not a lost fix**:

- **7.486 is the SYNTHESIS Fe II** (RYA-305/341) — produced by flux-space synthesis, the
  *ratified ionization arbiter* (RYA-341 DECISION 2). `validate_fe_rya238` does **not** run
  synthesis; it runs `abundances_derive.run()` = the **EW→MOOG path**.
- **7.700 is the EW-path Fe II** = RYA-352's clean-EW-quality pool value (documented **7.696**,
  HIGH). RYA-352 explicitly: *"a88ef0f's 7.466 'Asplund landing' is an artifact … the clean
  EW-quality pool reads 7.696 (high) … the EW path is **not** forced onto Asplund; synthesis
  stays the ionization arbiter."*

So the EW gate was **never expected** to read 7.486. Main correctly reflects 352 (EW high) and
341/342 (synth balanced). The two numbers come from two paths; the gate scores the EW one.

## Step 2/3 — current-main gate (verbatim) and the residual

```
── Solar Fe gate — PRIMARY (scale-robust, RYA-336) ──
Fe I reduced-EW slope = -0.011  -> PASS
Fe I-Fe II ionization = -0.184  -> FAIL (|ΔFe| < 0.05)
Fe I scatter          =  0.138  -> FAIL (< 0.1)
── Absolute A(Fe): scale-aware DIAGNOSTIC (window [7.44,7.58]) ──
A(Fe I)  NLTE abs = 7.516  -> PASS
A(Fe II) NLTE abs = 7.700  -> FAIL   (7 clean EW lines)
►► Solar Fe VERDICT (scale-robust primary): FAIL [slope ✓ · ionization ✗ · scatter ✗]
```

The gate's two real failures, named and routed:

1. **Ionization (−0.184) + Fe II-absolute (7.700) FAIL = the EW-vs-synth path mismatch.**
   The gate computes ΔFe(I−II) and the Fe II absolute on the **EW-path Fe II (7.700)**, which
   RYA-352 documents as *high by design*. The project's **ratified ionization arbiter is
   synthesis** (RYA-341 DECISION 2 → Fe II 7.500, ΔFe −0.015 — balanced). **The gate is reading
   the wrong Fe II path for the ionization verdict.** This is a *gate-architecture* question
   (should `validate_fe_rya238`'s ionization metric source Fe II from the synthesis arbiter, or
   treat the EW-path Fe II as a known-high diagnostic like the absolute scale?), **not** a
   re-landing of 305/352. → **fresh RCA / new ticket.**

2. **Fe I scatter (0.138 vs <0.10)** — a separate, smaller item on the *Fe I* pool (62 lines,
   1D scatter 0.139), unrelated to the Fe II stack. → its own look (pool curation / per-type
   scatter floor, cf. RYA-277/395).

## RYA-336 residual verdict

**RYA-336 is already on `main`** (`9ccf064`) — the ticket's "still Backlog/never-started" is
**stale**. Its scale-aware diagnostic ([7.44,7.58], centre 7.51 = 7.46 + 0.05 1D-3D offset) is
active and **Fe I 7.516 PASSES** it. So the original 336 absolute-scale residual (MPIA 1D-NLTE
7.516 vs [7.41,7.51], the +0.006) is **resolved** for Fe I. The Fe II-absolute FAIL (7.700) is
the **EW-path bias**, not a scale-offset — so it is **not** a 336 item; it routes with failure
(1) above.

## Net

- **No integration gap; nothing re-landed** (ticket STOP honoured).
- Solar Fe baseline on main = the documented solved stack, *correctly integrated*. The
  "balanced −0.015 / 7.500" lives in the **synthesis** arbiter (RYA-341, on main via 342); the
  **EW gate** reads the **high-by-design** 7.700 (RYA-352).
- **Fresh RCA owed (new ticket):** reconcile the solar Fe **gate's** Fe II ionization metric
  with the ratified synthesis arbiter — i.e. the EW gate should not score a verdict-level
  ionization FAIL on a Fe II path the project decided is *not* the arbiter. Plus the secondary
  Fe I scatter 0.138 item.
- 281's Procyon floor stays gated on this: it is differential against the solar baseline, and
  the baseline's gate verdict is an EW-vs-synth-arbiter question, not a moving number to re-land.
