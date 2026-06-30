# RYA-486 — Scope decision: the CNO redo after RYA-485

**Status:** decided (scope-only ticket; no runs triggered here).
**Authority:** this document fixes the scope of the CNO re-runs that follow the
RYA-485 step-back, so they do not sprawl. It does **not** run solar or Procyon CNO
(those are RYA-485-loop / RYA-348 / RYA-371, each gated below).

## Why a redo is owed at all

RYA-485 moved the solar CNO foundation: it verified solar O is genuinely 3D-NLTE
(not silent-LTE), triangulated three independent non-Kitt solar 777 references
(IAG Reiners-2016, IAG-telluric Baker+2020, Wallace-2011 KPNO), and surfaced a
differential regime-mismatch (Procyon 1D vs Sun 3D). Once the iterate-solar loop
converges with those references + the IR channels opened, the solar CNO
denominator will have moved and/or gained indicators. Every downstream CNO number
is differential against that denominator, so it must be re-derived against the new
one. This ticket says *how much* gets re-derived.

## Decision A — Solar CNO redo is FULL C+N+O, not O-alone

**Grounded in code, not asserted.** `pipeline/cno_synthesis.py` (module docstring,
lines 7-9): "C, N and O are coupled through molecular equilibrium (CO / CN / CH),
so we synthesize rather than invert EWs." The diagnostic registry comment (line
242): "CN needs A(C); [O I] is tied to A(C) via CO." Moving the solar O
denominator therefore shifts the molecular carbon indicators (CH / C2) and drags N
along through CN. You cannot cleanly redo O without C, and N rides the coupling.
The RYA-371 / RYA-369 CO-equilibrium non-convergence was this coupling already
biting. **Solar O redo = solar CNO redo. One system, done properly.**

**The IR is a genuine unlock, not a re-reference.** The new atlases reach the
near-IR, adding O I 844 nm and 926 nm — both already in the Amarsi NLTE grid
(`pipeline/nlte_cno.py:87-88`: 844nm = 8445.5-8447.5 A, 926nm = 9259.5-9267.0 A,
multiplet-averaged corrections) — plus CO / CN bands. These are independent-
continuum O and C channels the optical cannot give, so the redo is *more
constrained* than the current solar CNO, capable of shrinking the spread rather
than merely re-confirming it.

**Method:** runs through the RYA-485 iterate-til-right loop — one variable per
iteration, the solar control must hold, never tune toward a literature/Asplund
target; chase internal consistency and correct RT regime. Re-bank solar C, N, O.

## Decision B — Procyon CNO redo is MANDATORY and SURGICAL (CNO only)

Procyon CNO is differential against solar CNO by definition. Solar denominator
moves to every Procyon [C/H], [N/H], [O/H] moves. Non-optional. Re-run the
RYA-348 Procyon arms (VIS + UVES O I 777 + IR if newly available) against the
**new** solar CNO. **Species scope = C, N, O only.** Carry forward and re-evaluate
against the new solar values:

- the carbon spread 0.253-dex finding (vs new solar C),
- the O I 777 provisional [O/H] +0.085 +/- 0.186 (vs new solar O + the continuum-
  lever triangulation from the three RYA-485 references),
- the UV arm stays PARKED — RYA-348 Phase 3 found three hard blockers (FUV
  pseudo-continuum synthesis not science-grade, solar UV reference unresolvable,
  FUV C I not in the Amarsi grid); the unblock work is tracked in RYA-487.

## Decision C — Full Procyon 27-element run (RYA-404) is NOT a CNO side-effect

The solar CNO work changes the C / N / O denominators and the molecular/NLTE
method for **those three species only**. It does not touch Fe, alpha-elements, or
iron-peak metals: their denominators and methods are unchanged and already
validated against the Gaia Benchmark Stars.

**Grounded in code:** the metal abundances come from the EW-inversion path
(`pipeline/abundances_derive.py`), which does **not** import the CNO synthesis
engine (zero references to `cno_synthesis` / `run_cno`). The two paths are
independent, so re-deriving CNO cannot perturb a banked metal. And RYA-404 fires
on its own gates — the F-star acceptance floor (`config/constants.py` F profile:
`nlte_available = False`, `fe1_scatter_max = 0.222`) plus the NLTE engine — not as
a CNO consequence. **Redo CNO surgically; leave the validated metals banked.**

**The one caveat that would flip this, and its current determination.** IF the
iterate-solar loop surfaces a *pipeline-wide* bug — e.g. a continuum-normalization
issue that also touches metal lines, not just CNO — THEN a full re-run is
warranted, because the metals genuinely would be affected. That is "if we find a
general bug, diagnosed," NOT "redo everything to be safe."

*As of RYA-485 (merge d27d260): the caveat does NOT fire.* RYA-485 changed only
data manifests, docs, scripts and tests — no `pipeline/` or `config/` production
code — and its findings were CNO-method-specific (the 3D-vs-1D differential
regime-mismatch, the solar-O-is-3D-NLTE verification, the three 777 references).
None of these is a continuum-normalization defect on the shared path that feeds
metal EWs. So 404 stays separately gated; metals stay banked. Re-evaluate this
determination if a later loop iteration touches shared normalization code.

## Sequencing

1. **Solar CNO redo** (full C+N+O, all three references, IR-inclusive) — via the
   RYA-485 loop. The foundation.
2. **Procyon CNO redo** — differential vs the new solar CNO. Surgical, CNO only.
3. **Full Procyon 404** — only on its own gates; only if a loop iteration surfaces
   a diagnosed pipeline-wide bug.

## Gate

Do not start the Procyon CNO redo until solar CNO is re-banked and stable across
the three references and the verified RT regime (the RYA-485 loop converged). Do
not trigger 404 unless a pipeline-wide finding justifies it; let the finding set
the scope, not a blanket "to be safe."

## Drift guard

`tests/test_scope_cno_redo_rya486.py` pins the four facts this scope rests on, so
the scope is flagged if the code later drifts out from under it: (1) the CNO
synthesis engine declares the molecular coupling; (2) O I 844/926 remain in the
Amarsi grid; (3) the metals EW path stays independent of CNO synthesis; (4) the
F-star (404) acceptance profile stays its own NLTE-unavailable gate.
