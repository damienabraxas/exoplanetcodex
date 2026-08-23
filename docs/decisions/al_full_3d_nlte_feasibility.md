# Feasibility spike — full-3D-NLTE Al

**Ticket:** RYA-1008 (Phase-0 feasibility spike; desk research only, half-day box —
modelled on RYA-444). **Date:** 2026-08-23.
**Branch:** ryandamienschmitt/rya-1008-frontier-feasibility-spike-full-3d-nlte-al-none-exists-in.

**Verdict: AMBER — COLLABORATE for the grid, BUILD-CAPABLE for one star.**
The frontier gap is real but it is *narrower and better defined* than the ticket states,
and the one piece that gates everything (the 3D NLTE RT solver) is still request-only.

No code obtained, no cube downloaded, no synthesis run, no pipeline touched.

**Critical-path statement (holds regardless of verdict):** full-3D-NLTE Al is a frontier
side-mission. It is NOT on any milestone path. Per the ticket, real investment stays
gated on the M-dwarf tier and a 5-10 star portfolio. Nothing here asks for resources.

---

## 0. PREMISE CORRECTION (read this first)

The ticket's framing is wrong in one specific, decision-relevant way. Verified at the
primary source, not the citing paper.

**The ticket says:** *"There is NO full-3D-NLTE Al grid published anywhere -- not solar,
not off-solar"* and *"Al <3D>-NLTE solar anchor: NL17's 6.43 (mean-3D, not full 3D)."*

**Nordlander & Lind 2017 (A&A 607, A75) says, in Sect. 2.1, verbatim:**

> "For the Sun, we perform **full 3D NLTE** calculations using multi3d (Leenaarts &
> Carlsson 2009). The multi3d code is used as described by Amarsi et al. (2016), with
> background opacities computed for this work. We solve the statistical equilibrium
> using 26 short characteristic rays [...] for a range of abundances using **five
> snapshots** taken from an updated version [...] of the solar radiation hydrodynamical
> simulation used by Scott et al. (2015) but resampled from the original resolution of
> 240^2 x 230 to **60^2 x 101** (horizontal x vertical)."

So **solar full-3D-NLTE Al already exists**: A(Al) = 6.43 +- 0.03 is a *full* 3D NLTE
number on STAGGER snapshots, not a mean-3D one. What is <3D> in NL17 is the **released
correction grid** (1D MARCS + <3D> STAGGER, computed with the 1D `multi` code), which is
a different object in the same paper.

**Our own merged survey already had this right.** `docs/science/rya817_3dnlte_frontier.md`
lists Al under `FULL_3D_NLTE` in the *solar treatment* table (12 elements) and separately
under `GRID_MEAN3D` in the *off-solar grid* table. The ticket collapsed those two rows
into one claim. No repo change is needed — the correction is to the ticket premise.

**The gap, stated correctly:**

> Full-3D-NLTE Al exists for exactly **one star, the Sun** (NL17). What does not exist
> anywhere is a full-3D-NLTE Al treatment **off-solar** — no grid, and not even a single
> published metal-poor or giant star.

**This does not kill the prize; it sharpens it, and it makes the spike easier:**

* We would NOT be "first to do full-3D-NLTE Al" — that claim would be false and must not
  be made. We would be first to do it **off-solar**.
* It hands us a **known-answer checkpoint for Al itself** (6.43 +- 0.03, solar), so the
  rung-2 "reproduce a known Amarsi-group result" gate can be run on the *actual target
  element*, exercising the exact atom and the exact lines — not only on O or Fe as a
  proxy. That is strictly a better test, and it is new information for the ladder in the
  ticket comments.

---

## Gate 1 — 3D NLTE RT code obtainability (still the dominant risk)

RYA-444 (2026-06-25) found no clone-and-go 3D RT code. **Re-checked today against the
claim in the RYA-1008 comment that the public path is now "fully open." That claim is
half right — and the half that is wrong is the half that gates the build.**

| piece | what it is | access verified 2026-08-23 | obtainable |
|---|---|---|---|
| **DISPATCH** (Nordlund et al. 2018, MNRAS 477, 624) | the RHD framework M3DIS is built in | `bitbucket.org/aanordlund/dispatch2` — `is_private: false`; anonymous `git ls-remote` **succeeds**; **BSD 3-Clause**; Fortran; ~617 MB; last push **2026-08-21** | **PUBLIC** (new finding) |
| DISPATCH `experiments/stellar_atmospheres` | the 3D RHD *model-generation* experiment | present and complete in the public tree | **PUBLIC** |
| **MULTI3D** (Leenaarts & Carlsson 2009) — *the 3D NLTE RT solver* | the piece that actually does 3D NLTE | `experiments/Multi3D/` in the public tree contains **one 1080-byte file**, `mpi_debug_mod.f90`. The solver is **not there**. No public repo found anywhere. | **REQUEST** |
| **TSO.jl** (`github.com/pe1995/TSO.jl`) | tabulated opacities + EoS; wraps Turbospectrum to build the tables | public, **MIT**, last push 2026-08-21 | **PUBLIC** |
| **MUST.jl** (`github.com/pe1995/MUST.jl`) | DISPATCH/MULTI3D driver + post-processing | public listing, MIT, but its own README: *"after permissions to the repository have been granted [...] To get access please contact eitner@mpia.de"* | **PUBLIC-ish / ASK** |
| **Turbospectrum_NLTE** (`github.com/bertrandplez/Turbospectrum_NLTE`) | our existing 1D/<3D> synthesis code | public, GPL-3.0 | **PUBLIC (we already run it)** |
| **Balder** (Amarsi et al. 2018) | Amarsi's 3D NLTE solver | Caliskan et al. 2026 (arXiv 2605.05356) describes it as *"Balder, a custom version of Multi3D"* | **COLLABORATION-ONLY** |
| **Scate** (Hayek et al. 2011) | 3D **LTE** post-processor | no public repo (RYA-444) | REQUEST — and **LTE only, cannot do the NLTE step** |
| **RH 1.5D** (Pereira & Uitenbroek 2015, A&A 574 A3; Uitenbroek 2001) | **NLTE** statistical-equilibrium + synthesis, massively parallel, C/MPI/HDF5, ships example model atoms | `github.com/ITA-Solar/rh` — **public**, Zenodo DOI, pushed **2026-08-14** | **PUBLIC** — but **1.5D, NOT full 3D** (see below) |

### Assessing the "the public path is fully open" claim

The RYA-1008 comment reports a find: *"MULTI3D was updated (Bergemann group, M3DIS,
Hoppe & Bergemann 2024, arXiv 2405.06338) so its opacity sources MATCH TURBOSPECTRUM's
[...] and it is PUBLIC (github.com/bertrandplez/Turbospectrum_NLTE)."*

Three corrections, each checked at source:

1. **arXiv 2405.06338 is Eitner, Bergemann, Hoppe, Nordlund, Plez & Klevas 2024**
   (A&A 688, A52) — **the same M3DIS paper RYA-444 already assessed** and marked
   REQUEST. It is not a newer, separate "Hoppe & Bergemann 2024" find. (Hoppe is the
   third author; a genuine `Hoppe & Bergemann` MULTI3D paper is cited *in* that paper as
   **"in prep."**, i.e. it was unpublished at the time and is not an availability
   statement.)
2. **The "PUBLICLY AVAILABLE" footnote in that paper attaches to Turbospectrum, not to
   MULTI3D.** The paper's own words for the URL cited are the Turbospectrum
   documentation. Turbospectrum being public is true, already known, and already ours —
   it is not a route to a 3D NLTE solver.
3. **The paper contains no code-availability statement for M3DIS, MULTI3D, or the grid.**
   Verified by reading the full text.

**What the comment gets genuinely right, and it matters:** the *opacity* piece really is
routed around. `TSO.jl` is public and MIT, and its documented workflow builds the
DISPATCH/MULTI3D opacity + EoS tables **by driving Turbospectrum** — the `opacity_tables`
branch of the very repo we already run. So the "BLUE is the collaboration gap" worry is
correctly retired: **we would not need BLUE, and we would not need to ask anyone for an
opacity package.** That is a real and useful de-risking.

**Gate-1 finding:** the gap moved, it did not close. It is no longer *opacity* — it is
now precisely and only **the MULTI3D NLTE solver itself**. Everything around it
(framework, opacity, EoS, driver, cubes, atomic data) is public. That is a much better
position than RYA-444 recorded, and a much worse one than "fully open." A 3D **LTE**
capability is buildable from public parts today; **3D NLTE is not**, because the one
non-public part is the NLTE step, and NLTE is the entire point for Al (the corrections
we are chasing reach +1 dex).

Corroborating evidence that the solver circulates person-to-person rather than by
download — NL17's own acknowledgements:

> "The authors wish to thank Anish Amarsi for making his version of the multi3d code
> available to us along with tailored background line opacities"

...and Caliskan et al. **2026** still run Balder, "a custom version of Multi3D." Nine
years apart, same pattern.

### "Can we just find a download of MULTI3D via Balder?" — no, and the direction is backwards

Asked and checked, 2026-08-23. **Balder is *downstream* of MULTI3D, not a container for
it** — it is Amarsi's privately modified branch. Obtaining Balder would mean obtaining an
unreleased personal fork: the same access problem, plus the YELLOW line (unreleased code
is asked for, never dissected). No public Balder release exists; searches return only
papers that *use* it.

**The absence has a positive control.** MULTI3D is not merely "not found" — it is absent
from **`ITA-Solar`**, Carlsson & Leenaarts' own Oslo GitHub org, which actively publishes
`rh`, `helita`, `BifrostTools.jl`, `TraceParticles.jl` and more. The institution that owns
MULTI3D publishes freely and does not publish this. That is a decision, not an oversight,
and it is not one to route around.

**What the search did turn up: RH 1.5D is public and is a real NLTE solver.** This is a
genuine gain and it was missed on the first pass.

* **It is not full 3D.** The repo's geometry directories are `rh15d`, `rhf1d`, `rhsphere`
  — there is **no `rh2d` / `rh3d`**. 1.5D solves each atmospheric column independently and
  **neglects horizontal radiative transfer**, which is precisely the term that matters for
  Al's strong, scattering-dominated 3961 A resonance line. It **cannot** produce a
  full-3D Al result and must never be labelled as one.
* Uitenbroek's original RH *does* carry 2D/3D geometry, but its NSO page is dead
  (`www.nso.edu/staff/uitenbr/rh.html` -> HTTP 404), so that is not a clone-and-go route
  either.
* **Licensing caveat:** the repo has **no LICENSE file** despite "open source" in its
  README. The README explicitly invites use subject to citing Pereira & Uitenbroek 2015
  and Uitenbroek 2001, so intent is clear — but under the "was this made available for
  this use" rule, the absence of a formal license is worth noting before we depend on it.

**Why it still matters for the ladder:** RH 1.5D unblocks the **departure solve** on
public parts. Rung 3 ("3D-NLTE for one atomic line — add the departure solve") can be
built and validated as **model atom + statistical equilibrium** with RH 1.5D *now*,
swapping in the full-3D solver when access is granted. The 2025 review names 1.5D as the
recognised cheaper intermediate, so this is a legitimate rung rather than a fudge — it
just does not, and cannot, close the full-3D gap.

## Gate 2 — 3D cubes, including off-solar (NOT a blocker)

RYA-444 cleared the *solar* cube. Al's prize is off-solar, so this re-checks coverage.

* **Source:** Rodriguez Diaz et al. 2024, A&A 688, A480 (arXiv 2405.07872), "An extended
  and refined grid of 3D STAGGER model atmospheres — processed snapshots for stellar
  spectroscopy." **243 models** released, covering FGK (only the [Fe/H] = -4.00 models
  are excluded), so the metal-poor and subgiant nodes Al needs **are in the release**.
* **Site liveness checked today:** `https://3dsim.oca.eu/en/the-stagger-grid-2-0`
  returns HTTP 200 (and the `/fr/` variant likewise).
* **We hold essentially nothing 3D today:** `data/threed_grids/` is a single 4 KB CSV
  (`solar3d_metals_rya399.csv`). Every cube is a fresh fetch. Not downloaded in this spike.
* **A published cost lever, and it is large.** Rodriguez Diaz et al. find that
  **as few as two carefully selected snapshots** suffice for accurate equivalent widths
  (max abundance error **0.01 dex**), with the horizontal mesh downsampled 240 -> 80.
  NL17 used five snapshots; Caliskan 2026 used ten. Adopting the 2-snapshot recipe cuts
  the dominant cost term by ~2.5-5x **with a citable error budget**, rather than by
  guessing. This is the single biggest thing that makes a one-star attempt tractable.

**Gate-2 finding: GREEN.** Public, live, covers the off-solar regime, and comes with a
published recipe for spending less.

## Gate 3 — compute cost (order-of-magnitude)

Scaling for an ALI 3D NLTE solve:

```
cost ~ N_snapshots x N_cells x N_rays x N_frequencies x N_iterations x c
```

Anchored on NL17's actual solar setup (5 snapshots, 60^2 x 101 = 3.6e5 cells, 26 rays,
42-level Al I atom) and the Rodriguez Diaz 2-snapshot / 80-cell recipe
(80^2 x 101 = 6.5e5 cells); N_freq ~ 1e3 (frequency quadrature over ~1e2 bound-bound
transitions), N_iter ~ 50-100, c ~ 1e-8 - 1e-7 s per cell-ray-frequency short-characteristic
step:

| quantity | estimate |
|---|---|
| per snapshot, per abundance point | ~1e0 - 1e2 core-hours |
| **one star** (2 snapshots x ~5 abundance points) | **~1e2 - 1e3 core-hours** |
| peak working set (J_nu dominates: N_cells x N_freq x 8 B) | **~5 - 20 GB** |
| **the off-solar GRID** (~1e2 - 2.4e2 nodes) | **~1e4 - 1e5+ core-hours** |

Against the hardware we actually have:

* **Mac (Apple M4, 10 cores, 16 GiB, 69 GiB free)** — the main science box. One star is
  ~10-100 wall-hours at 10 cores: *arithmetically* in reach. But the ~5-20 GB working set
  against 16 GiB total is the binding constraint; it would need frequency-batching or MPI
  domain decomposition (which is exactly why MULTI3D is domain-decomposed) and would
  otherwise thrash. **Marginal, not comfortable.**
* **Sirius (2 physical cores, 16 GB)** — not viable; ~1e2-1e3 core-hours on 2 cores is
  weeks, and per [[feedback_grids_sirius_only]] it is the grid box, not a compute box.
* **The grid** at ~1e4-1e5+ core-hours is **cluster-only** on any reading. This is
  consistent with Caliskan et al. 2026 running a *smaller-atom* (57-level Ag), solar-only
  3D NLTE calculation on **Tetralith, a national supercomputer** — not on a workstation.

**Gate-3 finding: SPLIT.** One star is workstation-marginal / fat-node-comfortable. The
grid — the actual frontier prize — is cluster-class and always was. Note this is *worse
than RYA-444's* Gate 3, and legitimately so: RYA-444 costed 3D **LTE** CO (no NLTE
iteration) and concluded Orion-runnable. Al needs the NLTE solve, which is the expensive
term RYA-444 explicitly excluded. **Do not transfer RYA-444's GREEN compute verdict to
Al.**

These are order-of-magnitude figures from a scaling argument, not a benchmark. Confirm
against the chosen code's own timings before committing any wall-time.

---

## The easy-lift test target

**Recommendation: HD 140283** (the "Methuselah" metal-poor subgiant), Al I resonance
line 3961 A, one star, 2 snapshots.

Why this one:

* **It is genuinely undone.** NL17 did *full* 3D NLTE for the Sun only; HD 140283 it did
  in 1D and <3D> NLTE. A full-3D Al calculation there is unclaimed work, which is what
  the ticket asks for — *"a star/line regime NOT already done by Amarsi's group."*
* **It comes with a published bracket to land inside**, from NL17's own benchmark table
  (Teff 5777, log g 3.67, [Fe/H] -2.40): **LTE 3.14 +- 0.14 -> <3D> NLTE 3.57 +- 0.12,
  1D NLTE 3.52 +- 0.12**. The LTE->NLTE swing is **+0.4 dex** — large and unambiguous, so
  the test cannot pass for the wrong reason (this is the same criterion the ticket
  comments used to prefer O over Fe, and it is *better* satisfied here: Fe's solar 3D-NLTE
  net is -0.002, which RYA-981/817 already flagged as a weak test).
* **Its Teff is 5777 K** — the solar value. The atmosphere differs from the Sun almost
  entirely in metallicity and gravity, which isolates the variable of interest instead of
  confounding it.
* **A cube exists for it** in the public STAGGER 2.0 release (metal-poor FGK nodes are in
  the 243).
* **It is where our own machinery already gives up.** `pipeline/differential_bridge.py`
  names Sun -> HD 140283 (~2.5 dex) as the *failure* case that exceeds Jofre's
  differential limits. A direct 3D-NLTE calculation is precisely the tool for the regime
  the differential bridge cannot reach — so the capability would plug a known hole, not
  just produce a trophy.

**Run the solar case first anyway, as the gate.** Because of the premise correction
above, the Sun is now a *known-answer* rung for Al itself: reproduce 6.43 +- 0.03 with our
own atom and pipeline before touching HD 140283. Ordering: solar Al (known answer,
validates) -> HD 140283 (unclaimed, publishable). This is the ticket comments' rung-2
discipline applied to the target element directly.

**Cheaper rung if even one star is too much:** the **1.5D approximation** (independent 1D
NLTE solve per 3D column) is a recognised intermediate that, per the 2025 3D-NLTE review
(arXiv 2511.04254), "yields a major speed up." It is not full 3D and must never be
labelled as such, but it exercises the atom, the cube handling and the opacity chain at a
fraction of the cost.

---

## Build-vs-collaborate verdict

The ticket asks whether Al lands the same way as RYA-444/445. **Answer: yes for the
grid, and no for the capability — and that split is the finding.**

| | verdict |
|---|---|
| **The off-solar Al GRID** (the frontier prize) | **COLLABORATE.** ~1e4-1e5+ core-hours is cluster-class, and it needs the one non-public component (MULTI3D/Balder). Same landing as RYA-444/445. |
| **A one-star full-3D-NLTE Al capability** | **BUILD-CAPABLE, gated on one ask.** Every part is public *except* the NLTE solver. |
| **A 3D LTE capability** | **BUILD NOW if ever wanted** — DISPATCH (BSD-3) + TSO.jl (MIT) + Turbospectrum (ours) + public cubes is a complete public toolchain. Insufficient for Al (NLTE is the whole effect), but it is the honest smoke-test rung. |

**The collaboration ask is now much smaller and much more askable than RYA-445 assumed.**
RYA-445 was scoped as "we need your opacity package." That ask is dead — `TSO.jl` +
Turbospectrum replaces it. The remaining ask is one specific thing: **access to MULTI3D**
(Bergemann/MPIA, via `eitner@mpia.de`, which MUST.jl's own README already names as the
contact route) or Balder (Amarsi/Uppsala). Asking for one named code is a far lighter
request than asking for a proprietary opacity package.

**Timing is unchanged and should stay unchanged.** Per Ryan's standing rule the knock
happens with work in hand — Solar done and presented, alpha Cen minimum, 5-10 star
portfolio. Nothing here argues for accelerating that. The value of this spike is that
**when that knock happens, the ask is one code, not a shopping list** — and we can say
exactly what we would do with it and on which star.

**Recommended near-term posture (start at the bottom of the ladder, climb to Al).**
Ryan's direction, confirmed 2026-08-23: *"we start with the basics and move up to Al."*
That is now executable rather than blocked, because rung 1 needs only the public parts:

1. **Rung 1 — H Balmer in 3D LTE, the toolchain smoke test.** Buildable today from
   DISPATCH (BSD-3) + TSO.jl (MIT) + Turbospectrum (ours) + a public STAGGER snapshot.
   No ask, no collaboration, no cluster. This is the honest "does the 3D toolchain run
   and make a spectrum" rung, and it is the cheapest real intelligence we can hold.
2. **In parallel, the one email.** Rungs 2, 3 and 5 all wait on the same MULTI3D answer,
   so asking early costs nothing and unblocks three rungs at once. This is a *code-access*
   question, separate from the RYA-445 science collaboration — which stays gated on
   Solar-done-and-presented + alpha Cen + the 5-10 star portfolio, per the standing rule
   that a peer knocks with work in hand.
3. **Rung 4 (CO in 3D LTE) is also unblocked** and can be climbed without the answer if
   rung 2 stalls — useful, because it keeps the ladder moving while the ask is pending.
4. **Build the Al model atom + departure solve on public RH 1.5D** while the ask is
   pending. This is the most Al-specific work that can be done with no access at all: the
   42-level atom is assembled from published data (GREEN), and RH 1.5D validates the
   statistical equilibrium against NL17's 1D/<3D> numbers. When MULTI3D lands, the atom is
   already built and tested and only the RT geometry changes. **It does not and cannot
   substitute for full 3D** — 1.5D drops the horizontal RT that Al's 3961 A resonance line
   depends on.

Al itself (rung 5) stays unscheduled and unresourced, exactly as the ticket requires.
Nothing here asks to bring it forward; the point of the spike is that when the ladder
does reach Al, the path and its one blocker are already known.

## Where this sits on the phased ladder — and the one thing that changes

The ladder from the ticket comments, with today's Gate-1 evidence applied to each rung.
RYA-1008 is **rung 5**. Nothing below is scheduled; this is the map, not a plan.

| rung | what it is | needs | status after this spike |
|---|---|---|---|
| 1 | H Balmer, **3D LTE** — toolchain smoke test | 3D LTE only | **UNBLOCKED — fully public.** DISPATCH (BSD-3) + TSO.jl (MIT) + Turbospectrum (ours, GPL-3) + STAGGER 2.0 cubes is a complete public toolchain for 3D LTE. |
| 2 | reproduce a known Amarsi element on the Sun (O I 777, then Fe) | **3D NLTE** | **BLOCKED on the MULTI3D ask.** O I 777's published result is 3D *NLTE*; reproducing it needs the solver. |
| 3 | 3D-NLTE for one atomic line — add the departure solve | **3D NLTE** | **PARTLY UNBLOCKED** — the *departure solve itself* (model atom + statistical equilibrium) can be built and validated on public **RH 1.5D**; the *full-3D* version still needs the ask. |
| 4 | CO in 3D (vs Amarsi 2021 solar CO 8.48) | 3D **LTE** molecular (RYA-444) | **UNBLOCKED for the RT** — 444's own scoping notes CO is 3D LTE, not NLTE. Remaining work is molecular EOS in 3D, not code access. |
| 5 | **Al 3D-NLTE off-solar — this ticket** | **3D NLTE** | **BLOCKED on the same ask**, plus cluster-class compute for the grid. |

**The finding for the phased plan: the public/private line does not run along the
ladder's difficulty axis — it runs along the LTE/NLTE axis.** Rungs 1 and 4 are LTE and
are now buildable entirely from public parts. Rungs 2, 3 and 5 are NLTE and all sit
behind **one** ask, the same ask.

This qualifies the near-term plan in the newest ticket comment — *"Rungs 1-2 (smoke test,
reproduce known O then Fe) run on public opacity now."* **Rung 1 does. Rung 2 does not**,
and not for the reason the comment retires: the blocker was never the opacity (that part
is correctly resolved — TSO.jl + Turbospectrum genuinely replaces BLUE), it is the NLTE
solver itself. Reconnaissance can begin, but it stops at the end of rung 1 until MULTI3D
access is answered.

**So the cheapest informative next step is not rung 2 — it is the email.** One question
to `eitner@mpia.de` (MUST.jl's own stated contact route) resolves rungs 2, 3 and 5
simultaneously. Everything else on the NLTE half of the ladder is unschedulable until it
is answered, and no amount of public-parts work substitutes for it.

## Provenance check (the GREEN / YELLOW lines from the comments)

Every component recommended here is on the GREEN side — *available for this use*, not
merely citable:

* **DISPATCH** — BSD 3-Clause. Permits use, modification, redistribution. GREEN.
* **TSO.jl**, **MUST.jl** — MIT. GREEN, with one caveat below.
* **Turbospectrum_NLTE** — GPL-3.0, already in our stack. GREEN.
* **STAGGER 2.0 cubes** — published data release. GREEN (the "we build the synthesis,
  not the convection sim" scope from the comments).
* **Al model atom** — NL17's is 42 Al I levels built from **NIST energy levels, Kurucz
  transition probabilities, and Barklem broadening**: all published atomic data. Building
  our own from the same published sources is the FoRMATo-style GREEN path, i.e. normal
  science, exactly as the comments describe. It is also *small* — 42 levels vs Ag's 57 —
  which is a real cost advantage.
* **MULTI3D / Balder** — **ASK, never dissect.** This is the YELLOW boundary in the
  comments and it is exactly where we are. The correct move is the request; there is no
  version of reverse-engineering the solver that is acceptable, and the note recommends
  none.

**One caveat worth naming:** `MUST.jl` is a public MIT listing whose own README says
*"after permissions to the repository have been granted [...] To get access please
contact eitner@mpia.de."* A public listing plus an author who says permission is required
should be treated as **ask-first**, not as a silent green light — the comments' bright
line is "was this made available for this use," and the author has told us how he wants
that question asked.

## Firewall (RYA-161)

Any Codex 3D Al result validates against NL17's solar 6.43 as a **floor, never a tuning
target**. The solar rung is a pass/fail gate on a published number; the off-solar work is
unvalidated frontier and must be cross-checked against 1D/<3D> NLTE *trends* (NL17's
benchmark table gives five stars' worth), not against any single reference value.
Per [[feedback_a_borrowed_threshold_is_not_a_control]], the +0.4 dex bracket for
HD 140283 is a comparison, not an acceptance window — do not convert it into one.

## Sources

Nordlander & Lind 2017, A&A 607, A75 (arXiv 1708.01949) — Al model atom (42 Al I levels),
full 3D NLTE solar Al, 1D/<3D> grids, benchmark stars ·
Eitner, Bergemann, Hoppe, Nordlund, Plez & Klevas 2024, A&A 688, A52 (arXiv 2405.06338) — M3DIS ·
Nordlund, Ramsey, Popovas & Kuffmeier 2018, MNRAS 477, 624 — DISPATCH ·
Leenaarts & Carlsson 2009, ASPC 415, 87 — MULTI3D ·
Amarsi et al. 2018 — Balder ·
Caliskan, Amarsi, Jonsson, Grevesse & Sahoo 2026 (arXiv 2605.05356) — Ag I 3D NLTE, Balder on Tetralith ·
Rodriguez Diaz et al. 2024, A&A 688, A480 (arXiv 2405.07872) — STAGGER 2.0, 243 models, 2-snapshot recipe ·
Magic et al. 2013 — STAGGER grid · Hayek et al. 2011, A&A 535, A12 — Scate ·
arXiv 2511.04254 — 3D NLTE review, 1.5D approximation ·
Gustafsson et al. 2008 — MARCS ·
repo: `docs/science/rya817_3dnlte_frontier.md`, `pipeline/differential_bridge.py`,
`data/threed_grids/`.
