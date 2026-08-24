"""RYA-1015 — the authoritative element x model-type availability matrix.

Ryan's call: STOP re-deriving grid/model availability in conversation. This module
builds ONE reconciled matrix from the sources that already exist, so "what atoms have
what models" is READ OFF A GRID, never re-litigated.

It PROMOTES rather than rebuilds. Four independent sources are reconciled:

  1. CSV claim   -- data/curation/nlte_grid_availability.csv        (RYA-462)
  2. 3D claim    -- data/curation/threednlte_availability.csv       (RYA-817)
  3. CODE truth  -- config.constants NLTE/THREED_CORRECTION_ELEMENTS (what actually runs)
  4. DISK truth  -- a Sirius `find -L` snapshot                      (RYA-1015 scan)

**Every disagreement becomes a loud PROBLEM cell carrying the triggering fact.** A
matrix built from any ONE source reproduces the loop this ticket exists to end — the
RYA-597 Ti drift (CSV said Bergemann2011, code ran Mallinson2024) is exactly what a
single-source matrix cannot see.

THE DISK HALF IS A SNAPSHOT, NOT A LIVE READ. Only Mr. Code can reach Sirius, and CI
cannot, so the scan is committed as a dated artifact. The snapshot records the
`find -L` control result; a snapshot whose control FAILED is refused (see
`load_disk_snapshot`) because a blind scan reports absences that are pure artifact.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The canonical 28 = 27 (incl. Fe II counted separately, RYA-109) + Zn (RYA-757).
#: See pipeline/element_freeze.py — the record key is (element, ion), not element.
CANONICAL_28: tuple[str, ...] = (
    "Li", "C", "N", "O", "Na", "Mg", "Al", "Si", "P", "S", "K", "Ca", "Sc", "Ti",
    "V", "Cr", "Mn", "Fe", "Fe II", "Co", "Ni", "Cu", "Zn", "Sr", "Y", "Zr", "Ba", "Eu",
)

#: Model-type axis (RYA-400 "The Beast": 1D-vs-3D x LTE-vs-NLTE).
MODEL_TYPES: tuple[str, ...] = ("1D_LTE", "1D_NLTE", "MEAN3D_NLTE", "FULL_3D_NLTE")

#: WHAT A CELL MEANS -- fixed here so it is never re-argued.
#:
#: A cell states the availability of a RUNNABLE OR APPLICABLE GRID/DECK of that model
#: type for that element. It is a CAPABILITY statement, not a literature statement.
#:
#: The distinction is load-bearing and it is the thing that keeps getting re-litigated:
#: a VENDORED POST-HOC CORRECTION TABLE (e.g. data/nlte_grids/amarsi2019_cno/ for C/O,
#: vendor/1L-3NErrors/ for Fe) lets us APPLY somebody's 3D result. It does NOT give us
#: a full-3D model we can run. Per RYA-1008 no public full-3D NLTE stellar RT code
#: exists at all, so FULL_3D_NLTE is NONE for every element -- while the vendored
#: correction that IS applied in production is recorded in the cell's `facts` so the
#: capability gap and the production reality are both visible.
#:
#: The RYA-817 CSV (threednlte_availability.csv) encodes LITERATURE availability, which
#: is a different question and therefore reports different values for the same element.
#: Both are kept; neither is silently resolved into the other.

#: <3D> decks now consumed by an Engine-B route (RYA-821). Keyed element -> deck id.
_MEAN3D_WIRED: dict[str, str] = {"Al": "Al@mean3D"}

#: <3D> deck PRESENT on Sirius (element -> deck filename stem). Wiring one is the Al
#: pattern: a registry entry plus its own atmosphere. Not a research task.
_MEAN3D_DECKS: frozenset = frozenset({"Al", "Cr", "Eu", "Y"})

#: Model atoms held on Sirius (gerber_ts/atom.*). Half of a departure solve: with the
#: <3D> STAGGER atmosphere these become reachable the moment the tier-2 solver works.
_MODEL_ATOMS: frozenset = frozenset({
    "Al", "Ba", "Ca", "Co", "Cr", "Eu", "Fe", "H", "Mg", "Mn", "Na", "Ni", "O",
    "Si", "Sr", "Ti", "Y"})

#: 🔴 THE WHOLE MEAN3D_NLTE COLUMN *IS* TIER 2. RYA-1013 defines tier-2 as
#: "<3D>/1.5D NLTE", so every cell here is a tier-2 capability -- the question is never
#: "which tier", it is WHICH ROUTE delivers it and how far along that route we are.
#: RYA-1013 names the routes: CONSUME a published deck, BUILD-OUR-OWN, or the full-3D
#: check. These states track route + readiness, and deliberately do NOT reuse the word
#: "tier" for anything else.
T2_CONSUME_VALIDATED = "T2_CONSUME_VALIDATED"  # published deck wired AND proven end-to-end
T2_CONSUME_WIRED = "T2_CONSUME_WIRED"          # published deck wired; run still owed
T2_CONSUME_READY = "T2_CONSUME_READY"          # deck + atom on disk, nothing consumes them
T2_BUILD_OWED = "T2_BUILD_OWED"                # no deck; atom + <3D> atmosphere held ->
                                               # reachable only by the build-our-own route

#: Cell states.
HAVE = "HAVE"                  # CSV says + disk confirms + code uses
CODE_USES = "CODE_USES"        # the pipeline applies it now (code is ground truth)
CSV_ONLY = "CSV_ONLY"          # claimed but not on disk -> PROBLEM
DISK_ONLY = "NEEDS_WIRING"     # the grid IS on Sirius; no code path consumes it yet
REQUEST_ONLY = "REQUEST_ONLY"  # exists in literature, no public download
NONE = "NONE"                  # genuinely absent -> acquisition task
PROBLEM = "BROKEN"             # we hold it and it FAILS; `error` + `fix` say why/how

#: RYA-1015 disk-parse traps. `CO` in a Gerber filename is the CO MOLECULE, not
#: cobalt, and `MN` is Mn shouted. Mapping CO -> Co would invent a cobalt NLTE grid
#: we do not have; this is a real defect the naive parse produces.
_FILENAME_ELEMENT_FIXUPS = {"MN": "Mn"}
_NOT_ELEMENTS = {"CO"}  # molecular decks, never an atomic element cell



#: REPO-SIDE 3D HOLDINGS. The Sirius scan covers the departure decks only; the actual
#: 3D capability lives IN THE REPO and was invisible to a Sirius-only scan. This is
#: what we CAN RUN, keyed by the element it can correct.
#:
#: A matrix that reported FULL_3D_NLTE = NONE everywhere was WRONG: we hold a 3D-NLTE
#: Fe engine, 3D-NLTE C/O tables, and a 3D metals increment.
THREED_HOLDINGS: dict[str, dict] = {
    "Fe": {
        "path": "vendor/1L-3NErrors/",
        "kind": "FULL_3D_NLTE",
        "engine": "ENGINE-A-3DNLTE",
        "what": "Amarsi, Liljegren & Nissen 2022 (A&A 668 A68) 3D-NLTE Fe MLP "
                "(fe1_model_gt02.p / fe1_model_lt02.p / fe2_model.p)",
        "blocked_by": "RYA-923",
        "blocker": "URGENT/OPEN: the MLP returns NaN for EVERY in-domain line on main "
                   "(114 in-domain -> n=0). 1D-LTE legs still PASS, so only the "
                   "correction path regressed. Committed cells carry values from when "
                   "it worked (Fe I 7.604 n=114, Fe II 7.642 n=7) -- so the capability "
                   "is REAL but currently UNRUNNABLE.",
    },
    "C": {"path": "data/nlte_grids/amarsi2019_cno/", "kind": "FULL_3D_NLTE",
          "engine": "cno-3dnlte",
          "what": "Amarsi, Nissen & Skuladottir 2019 (A&A 630 A104) line-by-line "
                  "3D-NLTE / 1D-NLTE tables; 3D leg below Teff 6500 K"},
    "O": {"path": "data/nlte_grids/amarsi2019_cno/", "kind": "FULL_3D_NLTE",
          "engine": "cno-3dnlte",
          "what": "Amarsi 2019 3D-NLTE O I 777; [O I] 6300 is forbidden-LTE by "
                  "construction (RYA-447)"},
    "N": {"path": "data/nlte_grids/amarsi2019_cno/", "kind": "FULL_3D_NLTE",
          "engine": "cno-3dnlte",
          "what": "Amarsi 2019 CNO synthesis leg (N atomic departures are the separate "
                  "1D registry grid)"},
    # RYA-820, fetched 2026-08-23 into the Sirius grids-overflow spill area on the
    # root drive (codex-ext is 89% full). All three verified full 3D-NLTE at the
    # PRIMARY source before download, per RYA-820's do-not-mislabel rule.
    "Li": {"path": "grids-overflow/nlte/threed_offsolar/Li_wang2021_breidablik/",
           "kind": "FULL_3D_NLTE", "engine": "ENGINE-A-3DNLTE (Breidablik)",
           "what": "Wang et al. 2021, MNRAS 500, 2159 -- 3D-NLTE Li over the full "
                   "STAGGER grid (Teff 4000-7000, logg 1.5-5.0, [Fe/H] -4..0.5, "
                   "A(Li) -0.5..4.0; 610.4/670.8/812.6 nm). Zenodo 10.5281/zenodo.13829605"},
    "Mg": {"path": "grids-overflow/nlte/threed_offsolar/Mg_matsuno2024/",
           "kind": "FULL_3D_NLTE", "engine": "ENGINE-A-3DNLTE",
           "what": "Matsuno et al. 2024, A&A 688, A72 -- 3D-NLTE Mg corrections, "
                   "Balder on Stagger, 2646 rows. VizieR J/A+A/688/A72"},
    "Na": {"path": "grids-overflow/nlte/threed_offsolar/Na_canocchi2026/",
           "kind": "FULL_3D_NLTE", "engine": "ENGINE-A-3DNLTE",
           "what": "Canocchi et al. 2026, A&A 709, A90 -- 3D-NLTE Na, nine Na I lines, "
                   "Teff 4000-6500, logg 1.5-5.0, [Fe/H] -4..+0.5, ships RBF + FFNN "
                   "interpolators. Zenodo 10.5281/zenodo.19396611"},
    "Si": {"path": "data/threed_grids/solar3d_metals_rya399.csv", "kind": "MEAN3D_NLTE",
           "engine": "THREED_CORRECTION_ELEMENTS",
           "what": "Amarsi & Asplund 2017 (MNRAS 464, 264) solar 3D increment (RYA-399)"},
    "Ti": {"path": "data/threed_grids/solar3d_metals_rya399.csv", "kind": "MEAN3D_NLTE",
           "engine": "THREED_CORRECTION_ELEMENTS",
           "what": "Scott et al. 2015 Paper II (A&A 573, A26) solar 3D Ti (RYA-399)"},
    "Cr": {"path": "data/threed_grids/solar3d_metals_rya399.csv", "kind": "MEAN3D_NLTE",
           "engine": "THREED_CORRECTION_ELEMENTS",
           "what": "Scott et al. 2015 Paper II (A&A 573, A26) solar 3D Cr (RYA-399)"},
}

#: The <3D> STAGGER solar atmosphere we hold -- the model any <3D> route needs.
STAGGER_MEAN3D_ATMOSPHERE = "data/atmospheres/stagger_avg3d_rya442/sun_avg3d_stagger.mod"



#: GAP -> LINEAR TICKET, from the full RYA-1015 ticket sweep (all 1015 issues, 2026-08-23).
#: Every non-HAVE cell should name the ticket that owns it, so the matrix routes work
#: instead of only reporting absence.
GAP_TICKETS: dict[tuple[str, str], str] = {
    # Gerber TS-native decks: ON SIRIUS, unwired. RYA-710 is the umbrella.
    ("Al", "1D_NLTE"): "RYA-1005 (deck fully staged 74 GB md5-pinned; `gerber_nlte` "
                       "registers only Fe -- Engine-B-NLTE refuses on a REGISTRY LINE, "
                       "not on missing data) / RYA-801 / RYA-710",
    ("Cr", "MEAN3D_NLTE"): "RYA-800 (fetch+wire Cr Gerber ~25 GB) / RYA-710",
    ("Y", "1D_NLTE"): "RYA-802 (acquire+hold Y Gerber ~42 GB; corrects nothing today, "
                      "Y II unmeasured) / RYA-710",
    ("Y", "MEAN3D_NLTE"): "RYA-802 / RYA-710",
    ("Eu", "1D_NLTE"): "RYA-803 (Eu Gerber ~47 GB; model_atom source returns 0 bytes -- "
                       "BLOCKED, needs a live mirror) / RYA-710",
    ("Eu", "MEAN3D_NLTE"): "RYA-803 / RYA-710",
    ("Al", "MEAN3D_NLTE"): "RYA-821 (Al <3D>-NLTE Nordlander & Lind 2017 -- mean-3D, "
                           "wire honestly) / RYA-801",
    ("Li", "1D_NLTE"): "RYA-540 (grid STAGED on Sirius; PySME derivation solar delta "
                       "-0.030 does NOT reproduce the Lind-2009 small-positive anchor "
                       "-> NOT wired, validate-don't-tune STOP) / RYA-103",
    ("Ni", "1D_NLTE"): "RYA-710 (Gerber Ni MARCS deck on disk, unwired). Note RYA-731: "
                       "Ni's blocker is gf-scale (BAD_GF), not NLTE.",
    # Genuinely absent -> acquisition tickets
    ("Zn", "1D_NLTE"): "RYA-757 (Zn intake: 27->28 canonical + Sitnova+2022 NLTE grid -- "
                       "grid EXISTS in the literature, never acquired)",
    ("V", "1D_NLTE"): "RYA-470 / RYA-363 (NLTE_VOID: no Amarsi/GALAH grid and no usable "
                      "neutral V I model atom -- the only genuine void, RYA-404)",
    ("P", "1D_NLTE"): "RYA-717 (no NLTE grid; DATA_GAP, FUV needs HST/STIS)",
    ("Sc", "1D_NLTE"): "RYA-732 (no NLTE grid; HFS single blue line 4246)",
    ("Co", "1D_NLTE"): "RYA-727 (no NLTE grid; continuum-limited blue-edge lines)",
    ("Zr", "1D_NLTE"): "RYA-739 (no NLTE grid NEEDED -- Zr II is the majority ion, "
                       "LTE-robust; route is synthesis, RYA-560)",
    # 3D
    ("Fe", "FULL_3D_NLTE"): "RYA-923 (URGENT/OPEN: MLP returns NaN for every in-domain "
                            "line) / RYA-817 / RYA-924 (cells were computed on the Mac "
                            "from a temp scratchpad; inputs unrecoverable)",
    ("Fe II", "FULL_3D_NLTE"): "RYA-923 / RYA-817 / RYA-924",
    ("Li", "FULL_3D_NLTE"): "RYA-820 (off-solar 3D-NLTE Li, Wang 2021 -- PUBLIC, not fetched)",
    ("Mg", "FULL_3D_NLTE"): "RYA-820 (Matsuno 2024 Mg -- PUBLIC, not fetched)",
    ("Na", "FULL_3D_NLTE"): "RYA-820 (Canocchi 2024 Na -- PUBLIC, not fetched)",
}

#: Molecule gaps -> tickets.
MOLECULE_TICKETS: dict[str, str] = {
    "TiO": "RYA-751 (COOL-STAR MOLECULES CaH/MgH/TiO/SiO sweep -- deferred to the "
           "M-dwarf phase). Blocks the M-dwarf tier.",
    "VO": "RYA-751 (cool-star molecule, not carried)",
    "ZrO": "RYA-751 (cool-star molecule, not carried)",
    "MgH": "RYA-751 (COOL-STAR MOLECULES sweep -- deferred)",
    "FeH": "RYA-751 (cool-star molecule, not carried)",
    "SiH": "RYA-751 / RYA-189 (IR molecular line handling)",
    "CaH": "RYA-751 (COOL-STAR MOLECULES sweep -- deferred)",
    "H2O": "RYA-751 / RYA-503 (ExoMol/HITRAN->Turbospectrum converter EXISTS, CO-validated)",
    "C2": "RYA-742 (C2 Swan solar C indicator)",
    "CH": "RYA-743 (CH G-band -- 'THE BITCH')",
    "CN": "RYA-746 (CN solar N indicator, needs C prior)",
    "CO": "RYA-744 (CO dv=1/dv=2, CRIRES+ target)",
    "NH": "RYA-745 (NH solar N indicator, IR only)",
    "OH": "RYA-747 (OH solar O indicator -- NOVELTY CANDIDATE)",
}



#: KNOWN ANOMALIES / TRAPS per (element, model_type) -- defects, quirks and gotchas that
#: are true of the DATA, not of our wiring. Carried on the row so a reader of the CSV
#: meets the caveat at the same moment they meet the number.
ANOMALIES: dict[tuple[str, str], str] = {
    ("Ti", "1D_NLTE"): "RYA-597 DRIFT (FIXED by RYA-1015): the availability CSV cited "
        "Ti_Bergemann2011_MPIA.csv while the code has run Ti_Mallinson2024_PySME.csv "
        "since RYA-545. Bergemann-2011 used scaled-Drawin H collisions and was INFLATED "
        "~2x (+0.108 vs Mallinson's +0.052, RYA-546). The old file is retained as an "
        "immutable v1 reference, flagged superseded -- not deleted.",
    ("Ti", "MEAN3D_NLTE"): "Engine-B still runs the SCALED-DRAWIN atom.ti503b; swapping "
        "it to ab-initio is owed (RYA-548). Largest grid in the Gerber set "
        "(~26 GB zip / ~55 GB bin, 503 levels); Keeper throttled the fetch heavily.",
    ("Fe", "1D_NLTE"): "Fe runs through the fe-nlte subsystem, NOT "
        "NLTE_CORRECTION_ELEMENTS -- a registry-only reading reports it unwired. "
        "RYA-549: code=register=matrix reconciled; the NLTE delta is TINY "
        "(+0.010 solar / +0.004 55 Cnc), NOT the inflated Ti/Mn/Cr class.",
    ("Fe", "FULL_3D_NLTE"): "RYA-924: the committed 3D-NLTE cells were computed ON THE "
        "MAC from a temp scratchpad -- off-Sirius, inputs now unrecoverable, and "
        "assert_on_sirius was never wired to that route. RYA-922: the 3D-NLTE route was "
        "HARDCODED to Kitt Peak (module-level INSTRUMENT, no --instrument flag) so it "
        "could not produce 3D-NLTE for any other arm.",
    ("Fe II", "FULL_3D_NLTE"): "Same MLP and the same RYA-923 ceiling failure as Fe I; "
        "the Fe II leg used only 7 lines (n=7) against Fe I's 114.",
    ("Li", "1D_NLTE"): "ANOMALY (RYA-540): grid staged, derivation ATTEMPTED, and the "
        "solar delta came out -0.030 -- the WRONG SIGN against the Lind-2009 "
        "small-positive anchor, with the grid's J-label resolving oddly. Stopped under "
        "validate-don't-tune rather than shipped. Li I 6707 is a resonance line with a "
        "CN molecular blend (RYA-103) and needs the dedicated derivation, not EW.",
    ("Al", "1D_NLTE"): "ANOMALY (RYA-1005): the Gerber Al deck is FULLY STAGED (74 GB, "
        "md5-pinned) and Engine-B-NLTE still refuses Al -- `gerber_nlte` registers only "
        "Fe. The blocker is a REGISTRY LINE, not missing data. Also RYA-773: the "
        "Amarsi-2020 Al departure grid does NOT cover the clean 7835/36 + 8772/73 "
        "doublet -- the best lines are uncovered.",
    ("Al", "MEAN3D_NLTE"): "Nordlander & Lind 2017 is <3D>, NOT full 3D -- and its "
        "released GRID is <3D>/1D while its SOLAR value (6.43) is full 3D NLTE. "
        "Conflating the two is the RYA-1008 premise error.",
    ("Cu", "1D_NLTE"): "PENDING-OK: registered grid with a non-LOCKED verdict. Cu NLTE "
        "is small (+0.001 solar, reproduces Shi-2014). Cu production stays GET-DATA on "
        "measured-line QUALITY (RYA-395: 5105/5218/5782 re-measure), NOT on the grid.",
    ("Mn", "1D_NLTE"): "HFS-split lines: EW cannot measure them: the unlock was "
        "HFS-resolved SYNTHESIS on the Den Hartog e6S->z6P triplet (RYA-473). "
        "Registry hygiene open (RYA-566): problem_children says A(Mn) 5.554, the "
        "committed verdict is 5.466.",
    ("Ca", "1D_NLTE"): "RYA-413: the MPIA Ca I 6166 NLTE node was a PLACEHOLDER ZERO "
        "across all 72 nodes -- a registered-grid data defect. Ca NLTE primary "
        "(Amarsi vs defect-fixed MPIA) is still unadjudicated (RYA-414).",
    ("Sr", "1D_NLTE"): "RYA-433: the CITED primary (Mashonkina/INASAN) is NOT the "
        "vendored file -- Bergemann INSPECT was vendored and is demoted to cross-check. "
        "Sr II 4077/4215 are SATURATED resonance lines (RYA-430).",
    ("V", "1D_NLTE"): "THE ONLY GENUINE VOID: no Amarsi/GALAH grid AND no usable neutral "
        "V I model atom (verified absent from Zenodo 3982506). Interim route is the "
        "V II ionization anchor (RYA-470) -- measure around the unsolved atom.",
    ("Ni", "1D_NLTE"): "A Gerber Ni MARCS deck IS on Sirius, but Ni's blocker is "
        "gf-scale (BAD_GF), not NLTE -- wiring the grid does not unblock the element.",
    ("Zn", "1D_NLTE"): "Zn is the 28th canonical (RYA-757) and is NOT in TARGET_ELEMENTS "
        "(26) or the RYA-400 regime map -- it is canonical on paper and absent from the "
        "code's element roster. Sitnova+2022 grid exists and was never acquired.",
    ("Zr", "1D_NLTE"): "NO NLTE GRID NEEDED -- Zr II is the majority ion and LTE-robust "
        "(the Sr II / V II logic). Route is synthesis of the strong Zr II lines. "
        "Recorded so 'no grid' is never read as a gap.",
    ("Eu", "1D_NLTE"): "BLOCKED (RYA-803): the Gerber Eu model_atom source returns "
        "0 BYTES -- needs a live mirror before the ~47 GB fetch is worth starting. "
        "Eu II 6645 is also at the noise floor (6.8 mA, RYA-565/102).",
    ("Y", "1D_NLTE"): "Corrects NOTHING today: Y II (the dominant ion) is absent from "
        "the measured pool, so the grid would sit unused (RYA-802 is fetch-and-hold).",
    ("O", "1D_NLTE"): "O runs via the cno-3dnlte subsystem, not the registry. "
        "[O I] 6300 is forbidden-LTE by construction (RYA-447) and is blended with "
        "Ni I 6300.34, whose gf was a stale duplicate until RYA-543/365.",
    ("C", "1D_NLTE"): "3D leg applies BELOW Teff 6500 K, 1D-NLTE above (RYA-359/237). "
        "C I 5380 is excluded by Codex charter while Amarsi 2019 includes it (RYA-748).",
    ("N", "1D_NLTE"): "4 cool-metal-rich weak-line rails are excluded -> out-of-hull "
        "loud-flag (RYA-526). N I atomic departures are a SEPARATE registry grid from "
        "the CNO subsystem leg.",
}

#: Cross-cutting findings that belong to no single element. Emitted as their own rows so
#: the CSV carries the audit, not just the inventory.
AUDIT_ANOMALIES: list[tuple[str, str, str]] = [
    ("scan-method", "find without -L returns NOTHING through the grids symlink",
     "The Sirius grids directory is a SYMLINK onto the ntfs3 external drive. A plain find "
     "silently reports zero grids. Positive control MANDATORY before trusting any "
     "absence (RYA-1013 trap). Verified: find -L = 1 hit, plain find = 0 hits."),
    ("scan-coverage", "all THREE Sirius drives verified, not just the grids subtree",
     "sdb2 ext4 468G root / sda1 ext4 458G codex-data (291G) / sdc1 ntfs3 932G "
     "codex-ext (827G, 89% FULL). Full sweep set-diffed against the snapshot: 36 grid "
     "files both sides, SETS IDENTICAL. The narrow path was complete, now verified."),
    ("false-duplicate", "the '119.9 GB of duplicate grids' does NOT exist",
     "The .grd files inside worktrees are SYMLINKS (52-74 bytes) to the canonical "
     "copies. `find -L -printf %s` follows them and counts the target size repeatedly. "
     "Real usage of those worktrees is ~2 GB. Corrected in RYA-1015."),
    ("filename-trap", "NLTEgrid4TS_CO_MARCS is carbon monoxide, NOT cobalt",
     "The naive element parse maps CO -> Co and INVENTS a cobalt NLTE grid we do not "
     "have. MN is Mn shouted and must still resolve. Pinned by test."),
    ("provenance", "the Gerber decks have NO DOI and live on a MUTABLE Seafile share",
     "MPG Keeper, unversioned -- md5-pinned on our side; a re-fetch verifies against "
     "those hashes. The Zenodo/Amarsi decks are DOI-backed and versioned. Different "
     "provenance grades sitting in the same directory."),
    ("bookkeeping", "the Gerber H deck is staged on Sirius and unrecorded",
     "NLTEgrid_H_MARCS_May-10-2021.bin is present; RYA-804 is bookkeeping only -- commit "
     "its prov.json, do NOT re-download. H is the reference element (A(H)=12 by "
     "definition) so it is not a CANONICAL_28 row."),
    ("roster", "the canonical 28 is 27 (incl. Fe II counted separately) + Zn",
     "TARGET_ELEMENTS in the code is 26 and the RYA-400 regime map is 26; both omit "
     "Fe II as a distinct species and Zn entirely. The record key must be (element, "
     "ion), not element -- keying by element alone collides Fe I and Fe II."),
    ("capability", "we cannot COMPUTE full 3D for any element",
     "RYA-1008: no public full-3D NLTE stellar photosphere RT code exists. MULTI3D, "
     "Balder, M3DIS, Linfor3D and PORTA are all private; the public codes (RH 1.5D, "
     "Lightweaver) stop at 1.5D/2D. Every 3D capability we have is a PUBLISHED GRID WE "
     "APPLY, never a calculation we run."),
    ("unreconciled", "RYA-776 engine_coverage.csv is a 6th source NOT reconciled here",
     "data/catalog/engine_coverage.csv (225 rows) carries engine x WAVELENGTH coverage. "
     "This matrix has no wavelength axis, so band-level contradictions between the two "
     "are possible and unchecked. Follow-up owed."),
]


@dataclass
class Cell:
    element: str
    model_type: str
    state: str
    csv_claim: str | None = None
    code_grid: str | None = None
    disk_paths: list[str] = field(default_factory=list)
    facts: list[str] = field(default_factory=list)
    error: str = ""      # what is actually wrong (BROKEN cells)
    fix: str = ""        # the concrete next action (BROKEN / NEEDS_WIRING cells)
    routes: dict = field(default_factory=dict)   # route -> status (tier-2 cells)

    @property
    def is_problem(self) -> bool:
        return self.state in (PROBLEM, CSV_ONLY, DISK_ONLY)

    def as_dict(self) -> dict:
        return {
            "element": self.element,
            "model_type": self.model_type,
            "state": self.state,
            "csv_claim": self.csv_claim,
            "code_grid": self.code_grid,
            "disk_paths": self.disk_paths,
            "facts": self.facts,
            "error": self.error,
            "fix": self.fix,
            "routes": self.routes,
            "problem": self.is_problem,
        }


class DiskSnapshotError(RuntimeError):
    """The Sirius snapshot is missing, or its find -L control did not pass."""


#: Sources for the REPO-SIDE holdings (the Sirius prov.json files cover only the
#: departure decks). Cited so the matrix is usable as a reference, not just a checklist.
REPO_SOURCES: dict[str, dict] = {
    "vendor/1L-3NErrors/": {
        "citation": "Amarsi, Liljegren & Nissen 2022, A&A 668, A68 (3D-NLTE Fe MLP)",
        "source_url": "https://github.com/AlexanderLiljegren/1L-3NErrors (vendored)",
        "caveat": "Training domain from Jofre et al. 2014, A&A 564, A133 Tables 4/5 "
                  "('golden' Fe I/II lines). A(Fe) axis ceiling 7.5 -- see the BROKEN "
                  "root cause.",
    },
    "data/nlte_grids/amarsi2019_cno/": {
        "citation": "Amarsi, Nissen & Skuladottir 2019, A&A 630, A104",
        "source_url": "CDS VizieR J/A+A/630/A104",
        "caveat": "3D leg below Teff 6500 K, 1D-NLTE leg above (RYA-359/237).",
    },
    "data/threed_grids/solar3d_metals_rya399.csv": {
        "citation": "Amarsi & Asplund 2017 (MNRAS 464, 264) for Si; "
                    "Scott et al. 2015 Paper II (A&A 573, A26) for Ti/Cr",
        "source_url": "published tables, transcribed (RYA-399)",
        "caveat": "SOLAR increment only -- no off-solar parameter axis.",
    },
}


def load_grid_provenance() -> dict:
    """Where each Sirius deck came from: URL, citation, md5, caveat (RYA-1015 capture)."""
    import json
    f = ROOT / "data" / "audit" / "rya1015" / "grid_provenance.json"
    if not f.exists():
        return {}
    try:
        return json.loads(f.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def load_disk_snapshot(path: Path) -> dict[tuple[str, str], list[str]]:
    """Parse the committed Sirius `find -L` snapshot into {(element, model_type): paths}.

    LOUD-FAILS if the snapshot's positive control did not pass. A `find` that does not
    follow symlinks returns NOTHING for /srv/codex/grids (it is a symlink to the ntfs3
    drive), so an uncontrolled scan manufactures absences -- the RYA-1013 trap. We
    refuse to build a matrix on a blind scan rather than silently report NONE.
    """
    if not path.exists():
        raise DiskSnapshotError(
            f"Sirius disk snapshot missing: {path}. Regenerate with the RYA-1015 "
            f"scan (find -L + positive control) -- do NOT build the matrix without it."
        )
    text = path.read_text()
    if "CONTROL=PASS" not in text:
        raise DiskSnapshotError(
            f"{path}: positive control did not PASS. The scan was blind (find without "
            f"-L returns nothing through the /srv/codex/grids symlink), so every "
            f"absence in it is an artifact. Refusing to build the matrix."
        )
    out: dict[tuple[str, str], list[str]] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("|")
        if len(parts) < 2:
            continue
        directory, fname = parts[0], parts[1]
        parsed = classify_disk_file(fname)
        if parsed is None:
            continue
        element, model_type = parsed
        out.setdefault((element, model_type), []).append(f"{directory}/{fname}")
    return out


def classify_disk_file(fname: str) -> tuple[str, str] | None:
    """Map a grid filename to (element, model_type), or None if it is not an atomic grid.

    Model atoms (atom.*) are supporting inputs, not availability cells, so they map to
    None -- an atom without a departure grid does not make an element NLTE-capable.
    """
    # Gerber TS-native <3D> deck: NLTEgrid[4TS]_<El>_STAGGERmean3D_<date>.bin
    m = re.match(r"NLTEgrid(?:4TS)?_([A-Za-z]{1,2})_STAGGERmean3D_", fname)
    if m:
        el = _normalise_element(m.group(1))
        return (el, "MEAN3D_NLTE") if el else None
    # Gerber TS-native 1D deck: NLTEgrid4TS_<El>_MARCS_<date>.bin
    m = re.match(r"NLTEgrid(?:4TS)?_([A-Za-z]{1,2})_MARCS_", fname)
    if m:
        el = _normalise_element(m.group(1))
        return (el, "1D_NLTE") if el else None
    # Amarsi GALAH PySME departure grid: nlte_<El>_*_pysme.grd
    m = re.match(r"nlte_([A-Za-z]{1,2})_.*pysme\.grd$", fname)
    if m:
        el = _normalise_element(m.group(1))
        return (el, "1D_NLTE") if el else None
    return None


def _normalise_element(token: str) -> str | None:
    if token in _NOT_ELEMENTS:
        return None
    if token in _FILENAME_ELEMENT_FIXUPS:
        return _FILENAME_ELEMENT_FIXUPS[token]
    return token[0].upper() + token[1:].lower() if len(token) == 2 else token.upper()


def load_csv_claims(path: Path) -> dict[str, list[dict]]:
    """Rows of the RYA-462 availability CSV, grouped by element."""
    out: dict[str, list[dict]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            out.setdefault(row["element"].strip(), []).append(row)
    return out


def load_threed_claims(path: Path) -> dict[str, dict]:
    """Rows of the RYA-817 3D availability CSV, keyed by element."""
    out: dict[str, dict] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            out[row["element"].strip()] = row
    return out


def _code_truth() -> tuple[dict[str, dict], dict[str, dict]]:
    """What the pipeline ACTUALLY applies right now -- the ground truth."""
    import sys
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from config import constants as C
    return dict(C.NLTE_CORRECTION_ELEMENTS), dict(C.THREED_CORRECTION_ELEMENTS)


def build_matrix(snapshot_path: Path | None = None) -> dict:
    """Reconcile CSV claim vs disk reality vs code usage into the full matrix."""
    curation = ROOT / "data" / "curation"
    snapshot_path = snapshot_path or (ROOT / "data" / "audit" / "rya1015"
                                      / "sirius_scan_raw.txt")
    disk = load_disk_snapshot(snapshot_path)
    csv_claims = load_csv_claims(curation / "nlte_grid_availability.csv")
    threed = load_threed_claims(curation / "threednlte_availability.csv")
    nlte_code, threed_code = _code_truth()

    cells: list[Cell] = []
    for element in CANONICAL_28:
        base = element.split()[0] if element.startswith("Fe") else element
        for mt in MODEL_TYPES:
            cells.append(_reconcile(element, base, mt, disk, csv_claims,
                                    threed, nlte_code, threed_code))

    for c in cells:
        t = GAP_TICKETS.get((c.element, c.model_type))
        if t and c.state != HAVE:
            c.facts.append(f"OWNED BY: {t}")

    problems = [c for c in cells if c.is_problem]
    return {
        "generated": date.today().isoformat(),
        "generator": "pipeline/model_availability_matrix.py (RYA-1015)",
        "snapshot": str(snapshot_path.relative_to(ROOT)),
        "sources": {
            "csv_claim": "data/curation/nlte_grid_availability.csv (RYA-462)",
            "threed_claim": "data/curation/threednlte_availability.csv (RYA-817)",
            "code_truth": "config.constants NLTE/THREED_CORRECTION_ELEMENTS",
            "disk_truth": f"{snapshot_path.name} (Sirius find -L, control PASS)",
        },
        "elements": list(CANONICAL_28),
        "model_types": list(MODEL_TYPES),
        "cells": [c.as_dict() for c in cells],
        "problem_count": len(problems),
        "problems": [c.as_dict() for c in problems],
    }


def _reconcile(element: str, base: str, mt: str, disk, csv_claims, threed,
               nlte_code, threed_code) -> Cell:
    cell = Cell(element=element, model_type=mt, state=NONE)
    disk_paths = disk.get((base, mt), [])
    cell.disk_paths = disk_paths

    if mt == "1D_LTE":
        # Every element is reachable by 1D-LTE synthesis; this is the baseline route,
        # not a grid. Stated explicitly so the column is not mistaken for a gap.
        cell.state = HAVE
        cell.facts.append("1D-LTE synthesis is the universal baseline route (no grid).")
        return cell

    if mt == "1D_NLTE":
        code_entry = nlte_code.get(base)
        # cno-3dnlte is a WIRED PRODUCTION ROUTE for C/N/O, not a registry entry --
        # omitting it made O read DISK_ONLY when O is in fact live in production.
        rows = [r for r in csv_claims.get(base, [])
                if r["subsystem"] in ("registry-nlte", "fe-nlte", "cno-3dnlte")]
        claim = rows[0]["grid_file"] if rows else None
        # Not every wired element goes through NLTE_CORRECTION_ELEMENTS: C/N/O run via
        # the cno-3dnlte subsystem and Fe via fe-nlte. Treating the registry as the only
        # wiring route reports live production elements as unwired DISK_ONLY.
        wired_elsewhere = [r for r in rows
                           if r["subsystem"] in ("cno-3dnlte", "fe-nlte")
                           and r["wired"].strip().lower() == "true"]
        if wired_elsewhere and not code_entry:
            cell.csv_claim = claim
            cell.state = HAVE
            cell.facts.append(
                f"Wired via the {wired_elsewhere[0]['subsystem']} subsystem "
                f"({claim}), not NLTE_CORRECTION_ELEMENTS. "
                f"{len(disk_paths)} departure grid(s) on Sirius.")
            return cell
        cell.csv_claim = claim
        cell.code_grid = code_entry.get("grid") if code_entry else None

        if code_entry and claim and claim != code_entry["grid"]:
            cell.state = PROBLEM
            cell.facts.append(
                f"DRIFT: CSV claims '{claim}' but code applies "
                f"'{code_entry['grid']}'. Code is ground truth.")
        elif code_entry:
            cell.state = HAVE if (claim or disk_paths) else CODE_USES
            if not disk_paths:
                cell.facts.append(
                    "Applied from the vendored CSV in data/nlte_grids/; the "
                    "departure grid itself lives on Sirius only.")
        elif disk_paths:
            cell.state = DISK_ONLY
            cell.error = (
                f"{len(disk_paths)} departure grid(s) are ON SIRIUS but the element has "
                f"NO entry in NLTE_CORRECTION_ELEMENTS and no subsystem route. The data "
                f"was fetched and cannot be used.")
            cell.fix = (
                "Register the element in config/constants.py NLTE_CORRECTION_ELEMENTS "
                "(Engine-A) and/or the gerber_nlte registry (Engine-B), then derive the "
                "vendored correction CSV into data/nlte_grids/. See the OWNED BY ticket.")
        elif claim:
            cell.state = CSV_ONLY
            cell.facts.append(f"CSV claims '{claim}' but neither disk nor code has it.")
        else:
            cell.state = NONE
        return cell

    # --- 3D axis, from the RYA-817 CSV + disk ---
    row = threed.get(base, {})
    if mt == "MEAN3D_NLTE":
        offsolar = (row.get("offsolar_3d_nlte") or "").strip()
        solar = (row.get("solar_3d_nlte") or "").strip()
        # ORDER MATTERS. A held deck outranks a solar-only scalar increment: marking Cr
        # HAVE on solar3d_metals (a SOLAR-ONLY number) hid its full parameter-space <3D>
        # deck sitting unwired. A HAVE that masks a NEEDS_WIRING is the worst kind of
        # wrong -- it looks finished.
        # RYA-1013 names TWO independent routes to the same tier-2 number. An element
        # can sit at different readiness on each, and collapsing them to one state hides
        # exactly the case worth having: BOTH available => a cross-check, because the
        # routes share no machinery. Al is the first element in that position.
        cell.routes = {
            "consume": ("WIRED" if base in _MEAN3D_WIRED else
                        "READY" if base in _MEAN3D_DECKS else "no deck"),
            "build_our_own": ("OWED — atom + <3D> atmosphere held"
                              if base in _MODEL_ATOMS else "blocked — no model atom"),
        }
        if base in _MEAN3D_WIRED and base in _MODEL_ATOMS:
            cell.routes["cross_check"] = (
                "BOTH routes available — consume is wired and build-our-own is "
                "reachable, and they share no machinery. Running both gives an "
                "independent check of the same number (the RYA-1013 triangle logic, "
                "applied to this element).")

        if base in _MEAN3D_WIRED:
            cell.state = T2_CONSUME_WIRED
            cell.code_grid = _MEAN3D_WIRED[base]
            cell.facts.append(
                f"WIRED (RYA-821) via gerber_nlte deck '{_MEAN3D_WIRED[base]}'. Keyed at "
                f"STAGGER coordinates (Teff 5777 / logg 4.44), NOT MARCS 5750/4.5 -- "
                f"measured 0 rows at the MARCS node, 31 at the STAGGER node -- so the "
                f"atmosphere is carried per-deck.")
            cell.state = PROBLEM
            cell.error = (
                "THE DECK IS FINE; THE VENDOR INTERPOLATOR IS NOT (verified 2026-08-23). "
                "Reading the deck DIRECTLY in Python at the solar node gives n_dep=101, "
                "n_lev=354, log tau -5.000..+5.000 -- EXACTLY our <3D> atmosphere's "
                "range -- with physical departures (b~1.0 in deep layers, deviating to "
                "0.41 at the surface), 100% finite, zero all-zero rows. Record layout "
                "verified: 500(id) + 4(n_dep) + 4(n_lev) + 101*8(tau) + 101*354*8(bvals) "
                "= 287,348 bytes = the observed stride across all 6,345 records. "
                "interpol_multi_nlte produced ALL-ZERO b-values and then corrupted the "
                "heap (glibc malloc assertion) on the SAME record this reader parses "
                "correctly -- so the fault is that binary, not the data. "
                "EARLIER ATTEMPT DETAIL: The registry wiring is "
                "correct; interpol_modeles_nlte cannot use this deck. TWO faults, both "
                "measured, with the MARCS Al deck passing as a control in the same "
                "harness: (1) fed the <3D> atmosphere it CRASHES -- 'Bad real number in "
                "item 1 of list input' at interpol_modeles_nlte.f:1284. babsma reads "
                "TAU5000 SCALE fine, this binary does NOT; validating one says nothing "
                "about the other. (2) fed a MARCS corner it exits rc=0 and writes NO "
                "departure file -- 'ERROR: no match found for model 1', because the deck "
                "is keyed at Teff 5777 / logg 4.44 and MARCS models exist at 5750 / 4.50. "
                "It also prints 'Min abund is ********', a Fortran field overflow "
                "consistent with the nan entries in this aux table's mass/Vturb columns.")
            cell.fix = (
                "BLOCKED ON AN INPUT THAT IS NOT PUBLICLY DISTRIBUTED. Traced to the "
                "source: interpol_modeles_nlte.f reads ONLY native MARCS (its two flags "
                "are `test` and binary-vs-ascii MARCS -- there is no averaged-format "
                "mode), and native MARCS requires tauR and Pg per depth. The public <3D> "
                "STAGGER models carry neither: BOTH archives on the MPG Keeper share "
                "(average_stagger_grid_forTSv20.zip and ..._forMULTI1D.zip) hold the same "
                "5-column mul23 TAU5000 form -- log tau500, T, n_e, V, v_mic -- and the "
                "MULTI1D copy is byte-identical in content to the one we already hold. "
                "So the deck's aux names MARCS files (p5777_g+4.4...) that the "
                "distribution does not ship; Gerber's group evidently converted <3D> "
                "models to MARCS internally. "
                "BEST OPTION -- READ THE DECK DIRECTLY. The layout is fully determined "
                "and verified, so the vendor binary is not needed at all, and this also "
                "dissolves the corner-model problem: a direct reader takes the node from "
                "the aux instead of from 8 MARCS files. "
                "FALLBACKS: (a) ask the Gerber/Bergemann group for the MARCS-format <3D> "
                "models the aux indexes -- one email, same ask-shape as RYA-1008; "
                "(b) derive Pg and tauR ourselves (Pg from the EOS tables we hold; tauR "
                "needs Rosseland opacities, the harder half); (c) reach Al <3D> by the "
                "BUILD-OUR-OWN route instead, which does not use this binary. "
                "DO NOT feed a 5750/4.50 MARCS corner to force a run -- it pairs <3D> "
                "departures with a different atmosphere at a node the deck does not "
                "contain, and rc=0 makes that look like success.")
        elif base in _MEAN3D_DECKS:
            cell.state = T2_CONSUME_READY
            cell.disk_paths = disk_paths or cell.disk_paths
            cell.error = ("A <3D> STAGGERmean3D deck AND its model atom are BOTH on "
                          "Sirius, and nothing consumes them.")
            cell.fix = ("Wire on the Al pattern (RYA-821): add a `<El>@mean3D` entry to "
                        "gerber_nlte.DECKS carrying its OWN atmosphere. Mechanical, not "
                        "research.")
            if h := THREED_HOLDINGS.get(base):
                cell.facts.append(
                    f"NOTE: a solar-only increment also exists ({h['path']}) -- that is "
                    f"a SCALAR at the solar node, not a parameter-space grid. It must "
                    f"not be reported as though it were this deck.")
        elif base in _MODEL_ATOMS:
            cell.state = T2_BUILD_OWED
            cell.error = ("No <3D> deck. We DO hold the model atom and the <3D> STAGGER "
                          "atmosphere -- the missing piece is the departure solve.")
            cell.fix = ("Reachable only by the RYA-1013 BUILD-OUR-OWN route -- this IS tier-2 "
                        "work, not a different tier. NOT an acquisition. "
                        "Nothing to fetch and nothing to ask for.")
        elif solar in ("FULL_3D_NLTE", "MEAN3D_NLTE") or offsolar == "GRID_MEAN3D":
            cell.state = REQUEST_ONLY
            cell.facts.append(
                f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}. "
                f"Published, but we hold NEITHER a deck NOR a model atom -- this one "
                f"genuinely has to be asked for.")
        else:
            cell.state = NONE
            if solar or offsolar:
                cell.facts.append(
                    f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}.")
        return cell

    # FULL_3D_NLTE -- capability, not literature. See the MODEL_TYPES note above.
    solar = (row.get("solar_3d_nlte") or "").strip()
    offsolar = (row.get("offsolar_3d_nlte") or "").strip()
    holding = (row.get("our_holding") or "").strip()
    # NOTE: THREED_CORRECTION_ELEMENTS holds the solar 3D INCREMENT (solar3d_metals),
    # which is a <3D> product and belongs to the MEAN3D_NLTE column. Populating
    # code_grid here made the full-3D cell name a grid while reading NONE -- a cell that
    # contradicted itself. The grid is reported once, in the column that owns it.
    code_entry = threed_code.get(base)

    # What we CAN RUN, from the repo-side holdings -- NOT a blanket NONE. We hold a
    # 3D-NLTE Fe MLP and 3D-NLTE C/N/O tables; a Sirius-only scan could not see them.
    h = THREED_HOLDINGS.get(base)
    if h and h["kind"] == "FULL_3D_NLTE":
        cell.code_grid = h["path"]
        if h.get("blocked_by"):
            cell.state = PROBLEM
            cell.error = (
                "ROOT CAUSE FOUND (RYA-1015, 2026-08-23): NOT a code regression -- a "
                "BOUNDARY condition on the A(Fe;3N) axis. _apply_aberr_to_line lets the "
                "axis drift to (line's own 1D-LTE abundance + its own correction). For "
                "solar Fe that is ~7.45 + ~0.053 = 7.503, which is 0.003 dex ABOVE the "
                "MLP's published ceiling of A(Fe)=7.5, so amarsi3d.aberr_for_line "
                "refuses it with stellar_ok=False and returns NaN. Measured cliff: "
                "a_1dlte 7.447 -> +0.0528 (in domain); 7.450 -> NaN. The solar Fe I pool "
                "centres at 7.45-7.47, so nearly every line falls just outside. The "
                "LINE-level checks (feature/delta_E/level) all PASS, which is why 114 "
                "lines look in-domain and then yield n=0. "
                "Original symptom: the MLP returns NaN for EVERY in-domain line. "
                "114 Fe I lines PASS the domain check and then produce n=0. The 1D-LTE "
                "legs still PASS to 3 decimals, so the line list, atmosphere, star "
                "params and EW inversion are all fine -- only the correction path "
                "regressed. Committed cells still carry Fe I 7.604 / Fe II 7.642 from "
                "when it worked, so the products LOOK healthy while the engine is dead.")
            cell.fix = (
                "PIN THE AXIS, do not widen the domain. Pass afe3n_axis = the star's "
                "converged A(Fe;3N) (7.46 solar) to _apply_aberr_to_line / "
                "amarsi3d.aberr_for_line instead of letting each line's own value drift. "
                "Verified: pin=7.46 -> +0.0513 in-domain; pin=7.51 -> NaN. This is what "
                "the RYA-817 afe3n_axis parameter exists for, and it matches the vendor "
                "README (a SINGLE stellar A(Fe;3N) iterated to convergence, not a "
                "per-line value). It does NOT relax the domain check, which RYA-923 "
                "forbids. Then: make the reactivation control BLOCK instead of emitting "
                "n=0 cells, and re-derive both cells (1D-LTE base moved n=152->322 under "
                "the PR #315 width fix). PHYSICS CHANGE -- needs Ryan's say-so before "
                "any cell is published (RYA-923). "
                "Superseded hypothesis: 1) Bisect nlte_corrections.py 5278efb..a7ff4e0 -- "
                "suspects are _apply_aberr_to_line returning NaN, or _in_grid "
                "disagreeing with amarsi3d.domains() so a line passes the MLP domain "
                "check and is then rejected by the grid check. "
                "2) Make the reactivation control BLOCK: it currently FAILs and lets the "
                "run emit n=0 cells that read as 'no lines in domain'. "
                "3) Re-derive the two committed cells (their 1D-LTE base moved n=152->322 "
                "under the PR #315 width fix). "
                "DO NOT relax the domain check to get numbers out (RYA-923).")
            cell.facts.append(f"BROKEN -- owned by {h['blocked_by']}.")
            cell.facts.append(f"Engine {h['engine']}: {h['what']} ({h['path']}).")
        else:
            cell.state = HAVE
            cell.facts.append(f"CAN RUN via {h['engine']}: {h['what']} ({h['path']}).")
            if h["path"].startswith("grids-overflow/"):
                cell.facts.append(
                    "FETCHED 2026-08-23 (RYA-820) and NOT YET WIRED to a driver -- the "
                    "grid is on disk and verified; the Engine-A-3DNLTE wiring + domain "
                    "check (RYA-817 discipline) is still owed.")
                cell.fix = (
                    "Wire as an Engine-A-3DNLTE treatment with a transition-energy / "
                    "parameter domain check; flag out-of-domain lines, never "
                    "extrapolate. Mind the A(X) ceiling trap that broke Fe (RYA-923).")
    else:
        cell.state = NONE
        cell.facts.append(
            "No 3D-NLTE correction or engine we can run for this element. "
            "(We cannot COMPUTE full 3D from scratch for anything -- RYA-1008: no "
            "public full-3D NLTE RT code -- so every 3D capability here is a "
            "published grid/model we apply.)")
    if solar or offsolar:
        cell.facts.append(
            f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}.")
    return cell


#: The two engines (pipeline/engine_selection.py). There is no "primary" — both are
#: products that get presented, and higher reach is BROADER, not better.
ENGINE_A = "Engine-A (EW + grid delta)"
ENGINE_B = "Engine-B (synthesis)"


def build_engine_matrix(matrix: dict) -> list[dict]:
    """Per element: what each ENGINE can actually run, and where its grid lives.

    Engine-A runs 1D-NLTE from a VENDORED departure CSV in data/nlte_grids/ (in-repo).
    Engine-B runs synthesis, and goes NLTE only when a TS-native Gerber deck is on
    Sirius. Either engine falls back to LTE, so "no grid" never means "no engine" --
    it means that engine is LTE-only for that element, which is the distinction that
    keeps getting lost.
    """
    import sys
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from config import constants as C
    nlte_code = dict(C.NLTE_CORRECTION_ELEMENTS)

    by = {(c["element"], c["model_type"]): c for c in matrix["cells"]}
    grids_dir = ROOT / "data" / "nlte_grids"
    rows = []
    for element in CANONICAL_28:
        base = element.split()[0] if element.startswith("Fe") else element
        code_entry = nlte_code.get(base)
        onedee = by[(element, "1D_NLTE")]
        mean3d = by[(element, "MEAN3D_NLTE")]

        # --- Engine-A: needs a vendored CSV it can read from the repo ---
        vendored = code_entry.get("grid") if code_entry else onedee.get("csv_claim")
        a_present = bool(vendored) and (grids_dir / vendored).exists() \
            if vendored and vendored.endswith(".csv") else bool(vendored)
        if code_entry:
            a_mode, a_where = "1D-NLTE", f"data/nlte_grids/{code_entry['grid']}"
        elif onedee["state"] == HAVE and vendored:
            a_mode, a_where = "1D-NLTE", f"wired via subsystem ({vendored})"
        else:
            a_mode, a_where = "LTE only", "no departure grid"

        # --- Engine-B: NLTE only with a TS-native deck on Sirius ---
        ts_1d = [p for p in onedee["disk_paths"] if "gerber_ts" in p]
        ts_3d = [p for p in mean3d["disk_paths"] if "gerber_ts" in p]
        if ts_3d:
            b_mode = "<3D>-NLTE deck on disk (UNWIRED)"
            b_where = "; ".join(Path(p).name for p in ts_3d)
        elif ts_1d:
            b_mode = "1D-NLTE (TS-native)"
            b_where = "; ".join(Path(p).name for p in ts_1d)
        else:
            b_mode, b_where = "LTE only", "no TS-native deck"

        # Engine-A-3DNLTE / the 3D route: a published 3D model we APPLY.
        h3 = THREED_HOLDINGS.get(base)
        if h3 and h3.get("blocked_by"):
            c_mode = f"3D-NLTE BROKEN ({h3['blocked_by']})"
            c_where = h3["path"]
        elif h3:
            c_mode = ("3D-NLTE" if h3["kind"] == "FULL_3D_NLTE" else "<3D> increment")
            c_where = h3["path"]
        else:
            c_mode, c_where = "none", "no 3D model held"

        rows.append({
            "engine_c_mode": c_mode,
            "engine_c_where": c_where,
            "element": element,
            "engine_a_mode": a_mode,
            "engine_a_where": a_where,
            "engine_a_present": a_present,
            "engine_b_mode": b_mode,
            "engine_b_where": b_where,
            "two_engine": a_mode != "LTE only" and not b_mode.startswith("LTE"),
        })
    return rows


def render_markdown(matrix: dict) -> str:
    """Compact element x model-type grid for humans."""
    by = {(c["element"], c["model_type"]): c for c in matrix["cells"]}
    mts = matrix["model_types"]
    out = [f"# Element x model availability (RYA-1015)", "",
           f"Generated {matrix['generated']} by `{matrix['generator']}`.",
           f"Disk half: `{matrix['snapshot']}` (Sirius `find -L`, control PASS).", "",
           "| element | " + " | ".join(mts) + " |",
           "|---|" + "---|" * len(mts)]
    for el in matrix["elements"]:
        row = [el]
        for mt in mts:
            c = by[(el, mt)]
            mark = "**" + c["state"] + "**" if c["problem"] else c["state"]
            row.append(mark)
        out.append("| " + " | ".join(row) + " |")
    out += ["", f"**PROBLEM cells: {matrix['problem_count']}**", ""]
    for c in matrix["problems"]:
        out.append(f"- `{c['element']}` / `{c['model_type']}` -> **{c['state']}**: "
                   + " ".join(c["facts"]))
    return "\n".join(out) + "\n"


#: Molecules matter for the M-dwarf tier and for C/N/O, and they live in a DIFFERENT
#: place from the atomic grids, which is exactly why they get forgotten: the LTE
#: linelists are vendored in-repo while the only molecular NLTE deck is on Sirius.
MOLECULES_EXPECTED: tuple[str, ...] = (
    # carried today (C/N/O coupling, RYA-360)
    "C2", "CH", "CN", "CO", "NH", "OH",
    # NOT carried -- the M-dwarf / cool-star gap, listed so the hole is visible
    "TiO", "VO", "ZrO", "MgH", "FeH", "SiH", "CaH", "H2O",
)


def build_molecule_matrix(snapshot_path: Path | None = None) -> list[dict]:
    """What we hold per MOLECULE: vendored LTE linelist, and any NLTE deck.

    Absent molecules are listed explicitly rather than omitted -- an inventory that
    only shows what you have cannot show you a gap.
    """
    snapshot_path = snapshot_path or (ROOT / "data" / "audit" / "rya1015"
                                      / "sirius_scan_raw.txt")
    lists_dir = ROOT / "data" / "linelists" / "molecular" / "turbospectrum"
    have_lists = {p.name for p in lists_dir.iterdir() if p.is_dir()} \
        if lists_dir.exists() else set()

    decks: dict[str, list[str]] = {}
    for line in snapshot_path.read_text().splitlines():
        if line.startswith("#") or "|" not in line:
            continue
        _, fname = line.split("|")[0], line.split("|")[1]
        m = re.match(r"NLTEgrid(?:4TS)?_([A-Za-z0-9]{2,4})_MARCS_", fname)
        if m and m.group(1).upper() in {x.upper() for x in MOLECULES_EXPECTED}:
            decks.setdefault(m.group(1).upper(), []).append(fname)

    rows = []
    for mol in MOLECULES_EXPECTED:
        lte = mol in have_lists
        nlte = decks.get(mol.upper(), [])
        rows.append({
            "molecule": mol,
            "lte_linelist": ("data/linelists/molecular/turbospectrum/" + mol)
                            if lte else None,
            "nlte_deck": nlte[0] if nlte else None,
            "state": "HAVE" if lte else "NONE",
            "ticket": MOLECULE_TICKETS.get(mol, ""),
        })
    return rows


def _molecule_table(molecules: list[dict] | None) -> str:
    if not molecules:
        return ""
    rows = []
    for m in molecules:
        lte = (f'<span class="st" style="color:#5fd38d">HAVE</span>'
               f'<div class=g>{m["lte_linelist"]}</div>') if m["lte_linelist"] \
              else '<span class="st" style="color:#6a7690">NONE</span>'
        nlte = (f'<span class="st" style="color:#5fd38d">HAVE</span>'
                f'<div class=g>{m["nlte_deck"]}</div>') if m["nlte_deck"] \
               else '<span class="st" style="color:#6a7690">NONE</span>'
        cls = "" if m["lte_linelist"] else ' class="prob"'
        rows.append(f'<tr><th>{m["molecule"]}</th><td{cls}>{lte}</td><td>{nlte}</td></tr>')
    return ('<h2 style="font-size:1.15rem;margin:2rem 0 .5rem">Molecules</h2>'
            '<p class="sub">Molecular data lives apart from the atomic grids, which is '
            'why it goes missing. LTE linelists are vendored in-repo; the only molecular '
            'NLTE deck is on Sirius. Absent molecules are listed, not omitted &mdash; an '
            'inventory that shows only what you have cannot show a gap.</p>'
            '<table><thead><tr><th>molecule</th><th>LTE linelist</th>'
            '<th>NLTE deck</th></tr></thead><tbody>'
            + "".join(rows) + '</tbody></table>')


def render_html(matrix: dict, engine_rows: list[dict],
                molecules: list[dict] | None = None) -> str:
    """Self-contained page: every canonical element x model type, with REAL grid names.

    Static by design -- the live site is GitHub Pages, so a fetch()-driven page would
    show nothing. What you can see is the point of a tracker.
    """
    by = {(c["element"], c["model_type"]): c for c in matrix["cells"]}
    eng = {r["element"]: r for r in engine_rows}

    def td(el, mt):
        c = by[(el, mt)]
        cls = {"HAVE": "have", "REQUEST_ONLY": "req", "NONE": "none",
               "DISK_ONLY": "prob", "CSV_ONLY": "prob", "PROBLEM": "prob",
               "CODE_USES": "have"}.get(c["state"], "none")
        names = "<br>".join(Path(p).name for p in c["disk_paths"]) or ""
        claim = c["code_grid"] or c["csv_claim"] or ""
        detail = "<br>".join(x for x in (claim, names) if x)
        return (f'<td class="{cls}"><span class="st">{c["state"]}</span>'
                f'{"<div class=g>" + detail + "</div>" if detail else ""}</td>')

    rows = []
    for el in matrix["elements"]:
        e = eng[el]
        rows.append(
            f'<tr><th>{el}</th>'
            + "".join(td(el, mt) for mt in matrix["model_types"])
            + f'<td class="eng">A: {e["engine_a_mode"]}<br>B: {e["engine_b_mode"]}</td></tr>')

    heads = "".join(f"<th>{mt.replace('_', '-')}</th>" for mt in matrix["model_types"])
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Element x model availability - The Exoplanet Codex</title>
<style>
:root{{color-scheme:dark}}
body{{background:#0a0e1a;color:#d8e0f0;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:0;padding:2rem}}
h1{{font-size:1.5rem;margin:0 0 .25rem}} .sub{{color:#8fa0c0;margin:0 0 1.5rem}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border:1px solid #223;padding:.5rem .6rem;vertical-align:top;text-align:left}}
thead th{{background:#141c30;position:sticky;top:0}}
tbody th{{background:#101830;font-weight:700;width:4rem}}
.st{{font-weight:700;font-size:11px;letter-spacing:.04em}}
.g{{color:#8fa0c0;font-size:11px;font-family:ui-monospace,Menlo,monospace;margin-top:.3rem;word-break:break-all}}
.have .st{{color:#5fd38d}} .req .st{{color:#e8c060}} .none .st{{color:#6a7690}}
.prob{{background:#2a1420}} .prob .st{{color:#ff7b8a}}
.eng{{color:#a8b6d0;font-size:11px}}
.legend{{margin:1.25rem 0;color:#8fa0c0;font-size:12px}}
.legend b{{color:#d8e0f0}}
</style></head><body>
<h1>Element &times; model availability</h1>
<p class="sub">All {len(matrix['elements'])} canonical species &times; {len(matrix['model_types'])} model types, with the actual grid on disk.
Generated {matrix['generated']} by <code>{matrix['generator']}</code> &mdash; reconciled across
the RYA-462 CSV, the RYA-817 3D CSV, the code, and a Sirius <code>find -L</code> scan (control PASS).</p>
<table><thead><tr><th>species</th>{heads}<th>engines</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table>
{_molecule_table(molecules)}
<p class="legend">
<b>HAVE</b> wired and usable &nbsp;|&nbsp; <b>REQUEST_ONLY</b> published, no public deck &nbsp;|&nbsp;
<b>NONE</b> genuinely absent &nbsp;|&nbsp; <b class="pr" style="color:#ff7b8a">DISK_ONLY</b> on Sirius but nothing consumes it
({matrix['problem_count']} such cells).<br>
<b>Engine-A</b> = EW + grid delta (vendored CSV, in repo). <b>Engine-B</b> = synthesis
(Turbospectrum; TS-native NLTE from a Gerber deck on Sirius). Neither is primary.<br>
<b>1D-LTE</b> is available for every species by synthesis and needs no grid.
<b>full-3D-NLTE</b> is NONE everywhere: no public full-3D NLTE stellar RT code exists (RYA-1008).
</p></body></html>
"""


def write_findings_csv(matrix: dict, engine_rows: list[dict],
                       molecules: list[dict], out: Path) -> Path:
    """One flat CSV of EVERY finding: species x model type, engines, and molecules.

    Long format, one row per (subject, model_type), so it sorts and filters in a
    spreadsheet without unpacking anything.
    """
    eng = {r["element"]: r for r in engine_rows}
    prov = load_grid_provenance()

    # Carry the NOTE the source CSVs already wrote -- RYA-462 and RYA-817 both hold
    # per-row prose that never reached the matrix.
    csv_notes: dict[tuple[str, str], str] = {}
    curation = ROOT / "data" / "curation"
    try:
        for row in csv.DictReader((curation / "nlte_grid_availability.csv").open()):
            sub = row["subsystem"]
            mt = ("MEAN3D_NLTE" if sub == "metals-3d"
                  else "FULL_3D_NLTE" if sub == "cno-3dnlte" else "1D_NLTE")
            key = (row["element"].strip(), mt)
            if row.get("note"):
                csv_notes[key] = (csv_notes.get(key, "") + " " + row["note"]).strip()
        for row in csv.DictReader((curation / "threednlte_availability.csv").open()):
            if row.get("note"):
                k = (row["element"].strip(), "FULL_3D_NLTE")
                csv_notes[k] = (csv_notes.get(k, "") + " " + row["note"]).strip()
    except (OSError, csv.Error):
        pass

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "subject", "model_type", "state", "can_run",
                    "engine", "grid_or_path", "disk_decks",
                    "source_url", "citation", "md5", "provenance_caveat",
                    "route_consume", "route_build_our_own", "cross_check",
                    "error", "fix", "owned_by_ticket", "anomalies",
                    "source_csv_note", "notes"])
        for c in matrix["cells"]:
            e = eng.get(c["element"], {})
            if c["model_type"] == "1D_LTE":
                engine, can = "Engine-B (synthesis)", "yes"
            elif c["model_type"] == "1D_NLTE":
                engine = f'Engine-A: {e.get("engine_a_mode","?")} | ' \
                         f'Engine-B: {e.get("engine_b_mode","?")}'
                can = "yes" if c["state"] in ("HAVE", "CODE_USES") else "no"
            else:
                engine = e.get("engine_c_mode", "none")
                can = "yes" if c["state"] in ("HAVE", "CODE_USES") else "no"
            owner = next((f[10:] for f in c["facts"] if f.startswith("OWNED BY:")), "")
            gp = (c["code_grid"] or c["csv_claim"] or "")
            src = REPO_SOURCES.get(gp, {})
            # Only cite a source when the cell actually HOLDS something. A NONE cell
            # naming a URL is the same self-contradiction as the solar3d_metals leak.
            if not src and not c["disk_paths"] and not c["code_grid"]:
                src = {}
            elif not src:
                # An element can hold decks from BOTH families (e.g. Al and O carry an
                # Amarsi GALAH .grd AND a Gerber TS deck). Cite every source we hold,
                # not just the first -- the point of the column is reference.
                base = c["element"].split()[0]
                fams = sorted({("gerber_ts" if "gerber" in p else "amarsi_galah")
                               for p in c["disk_paths"]})
                parts = [prov.get(base, {}).get(f, {}) for f in fams]
                parts = [x for x in parts if x]
                src = {
                    "source_url": " | ".join(f"[{f}] {x.get('source_url','')}"
                                             for f, x in zip(fams, parts)
                                             if x.get("source_url")),
                    "citation": " | ".join(dict.fromkeys(
                        x.get("citation") or "Amarsi et al. 2020, A&A 642, A62 "
                        "(GALAH DR3 departure grids)" for x in parts)),
                    "md5": " | ".join(x.get("md5", "") for x in parts if x.get("md5")),
                    "caveat": " | ".join(dict.fromkeys(
                        x.get("caveat", "") for x in parts if x.get("caveat"))),
                } if parts else {}
            w.writerow([
                "element", c["element"], c["model_type"], c["state"], can, engine,
                gp, "; ".join(Path(p).name for p in c["disk_paths"]),
                src.get("source_url", ""), src.get("citation", ""),
                src.get("md5", ""), src.get("caveat", ""),
                (c.get("routes") or {}).get("consume", ""),
                (c.get("routes") or {}).get("build_our_own", ""),
                (c.get("routes") or {}).get("cross_check", ""),
                c.get("error", ""), c.get("fix", ""), owner,
                ANOMALIES.get((c["element"], c["model_type"]), ""),
                csv_notes.get((c["element"], c["model_type"]), ""),
                " ".join(f for f in c["facts"] if not f.startswith("OWNED BY:")),
            ])
        for m in molecules:
            w.writerow(["molecule", m["molecule"], "LTE_linelist",
                        "HAVE" if m["lte_linelist"] else "NONE",
                        "yes" if m["lte_linelist"] else "no",
                        "Engine-B (synthesis)", m["lte_linelist"] or "", "",
                        "iSpec input/linelists/turbospectrum/molecules (RYA-360)"
                        if m["lte_linelist"] else "", "", "", "", "", "",
                        "", "acquire the linelist" if not m["lte_linelist"] else "",
                        m.get("ticket", ""),
                        "" if m["lte_linelist"] else
                        "NOT CARRIED -- cool-star/M-dwarf molecule. An inventory of only "
                        "what you have cannot show a gap.", "", ""])
            w.writerow(["molecule", m["molecule"], "NLTE_deck",
                        "HAVE" if m["nlte_deck"] else "NONE",
                        "yes" if m["nlte_deck"] else "no",
                        "Engine-B (TS-native)", "", m["nlte_deck"] or "",
                        "MPG Keeper (Seafile)" if m["nlte_deck"] else "",
                        "Gerber, Bergemann et al. 2023, A&A 669, A43"
                        if m["nlte_deck"] else "", "",
                        "NO DOI, mutable Seafile -> md5-pinned"
                        if m["nlte_deck"] else "", "", "", "",
                        "", "", m.get("ticket", ""),
                        "Only ONE molecular NLTE deck exists in our holdings (CO). "
                        "Every other molecule is LTE-only." if not m["nlte_deck"] else "",
                        "", ""])
        # Cross-cutting audit findings: they belong to no element, so they get their own
        # rows. Without these the CSV is an inventory; with them it is the audit.
        for kind, title, detail in AUDIT_ANOMALIES:
            w.writerow(["anomaly", kind, "", "FINDING", "", "", "", "", "", "", "", "",
                        "", "", "", "", "", title, "", detail])
    return out
