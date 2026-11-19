"""The published PROVENANCE CATEGORY (`grade`) of a solar C/N/O product. RYA-1230.

Ryan, 2026-09-26: *"grade was missing from all the previous products so make sure we
account for that as you rerun."* RYA-1220 stamped it once by script (22ea6567), and that
commit never reached main, so every C/N/O row in the store published without one -- and
the site's "Reference Grade preferred" rule silently had nothing to prefer. It is now
stamped by the PUBLISHER, so a re-run cannot lose it again.

The rule is RYA-1220's, on Fe's own axis (a product measured on an external PUBLISHED line
set is Reference Grade; the Codex's own selection is Codex Grade), with one addition that
the same axis already implies:

  selector SET-<name>        -> Theoretical Codex  when line_sets.SOURCES.yaml adopts <name>
                                                  as THEORETICAL_CODEX (RYA-1351): theoretical
                                                  gf, no published uncertainty, Ryan-adopted
  selector SET-<name>        -> Reference Grade  any other published reference line set
                                                  (SET-AGSS21, SET-LBP25)
  selector MOL-CN_AX_IR      -> Reference Grade  its fit windows ARE AGSS21's own published
                                                  CN A-X (0-0) line positions
                                                  (cno_synthesis._cn_ir_windows)
  selector REFERENCE         -> Reference Grade  the full laboratory-gf pool (Fe's axis)
  anything else              -> Codex Grade      the Codex's own indicator selection

No CNO pool is depth-split, so Deep Grade is never assigned.
"""
from __future__ import annotations

CNO_ELEMENTS = ("C", "N", "O")
#: Molecular diagnostics whose windows are a published reference set's own line positions.
REFERENCE_MOLECULAR_SELECTORS = ("MOL-CN_AX_IR",)


#: RYA-1233: the rule is Fe's axis, not a C/N/O one, so it applies to every element that
#: does not stamp its own -- Si published 26 products with NO grade, and the site's
#: "Reference Grade preferred" headline had nothing to prefer (it picked a 2-line row).
SELF_STAMPING = ("Fe",)


def grade_for(row: dict) -> str | None:
    """`grade` for a published row, or None for an element that stamps its own (Fe)."""
    if str(row.get("element")) in SELF_STAMPING:
        return None
    # RYA-470: the IR hyperfine special case grades as the pool it measured.
    from pipeline.free_continuum import base_selector
    sel = str(base_selector(row.get("selector")) or "")
    # RYA-1351: an adopted Theoretical Codex set is its own grade. Checked before the SET-
    # rule, which would otherwise stamp theoretical gf "Reference Grade".
    if sel.startswith("SET-"):
        from pipeline.line_sets import theoretical_codex
        if theoretical_codex(sel[len("SET-"):]):
            return "Theoretical Codex"
    if sel.startswith("SET-") or sel in REFERENCE_MOLECULAR_SELECTORS:
        return "Reference Grade"
    # RYA-1311: the full laboratory-gf pool is Reference Grade on Fe's own axis
    # (rya1178_emit_fe_schema.GRADE_FOR_LINE_SET["reference"]); it fell through to Codex Grade.
    if sel == "REFERENCE" or str(row.get("tier") or "") == "REFERENCE":
        return "Reference Grade"
    # RYA-1233: the Codex's DEEP pool (gf-graded lines above the 0.60 depth gate) is Deep
    # Grade, Fe's own axis -- it used to fall through to Codex Grade.
    if sel.startswith("DEEPGRADED") or str(row.get("tier") or "").startswith("DEEPGRADED"):
        return "Deep Grade"
    return "Codex Grade"
