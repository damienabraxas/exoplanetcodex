"""
RYA-1080 — the feed and the repo must not be able to disagree.
==============================================================
`data/products/solar/<El>.json` is the surface that PUBLISHES. `data/results/band_products/`
is the surface people GREP. When a product is in the first and not the second, every check
that reads the tree answers a question about the feed with an answer about the disk. That
produced two wrong diagnoses in one session — "HARPS has no graded product" and "Fe II has
never been re-run" — the second of which nearly triggered a from-scratch re-run of finished
work.

🔴 THE TRAP THIS GUARD EXISTS TO CLOSE IS NOT "copied_to IS None".
Eight rows had a non-null `copied_to` pointing into `/private/tmp/g3d/` and
`/private/tmp/sirius_orphans/` — outside the repo, untracked, gone at the next reboot. A
non-null field READS as reconciled. So the question asked here is always

    is there a file AT A REPO-RELATIVE PATH, GIT-TRACKED, whose sha256 is the one the
    feed recorded?

and never "is the field set?". `data/results/band_products/` is gitignored (.gitignore:87),
so presence on disk is not tracking either — both are checked.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Findings that mean the repo cannot back what the feed publishes. Any of these is a
#: merge blocker; REGENERABILITY_GAP is recorded and reported but is not one.
BLOCKING = ("MISSING_COPIED_TO", "COPIED_TO_OUTSIDE_REPO", "ARTIFACT_ABSENT",
            "ARTIFACT_UNTRACKED", "SHA_MISMATCH", "NO_RECORDED_SHA")


@dataclass(frozen=True)
class Finding:
    kind: str
    feed: str
    element: str
    ion: str
    band: str
    treatment: str
    detail: str

    @property
    def blocking(self) -> bool:
        return self.kind in BLOCKING

    def key(self) -> tuple:
        """Identity for pinning a known set, stable against wording changes in `detail`."""
        return (self.kind, self.feed, self.element, self.ion, self.band, self.treatment)


def _tracked(root: Path) -> set[str]:
    out = subprocess.run(["git", "ls-files"], cwd=root, capture_output=True, text=True)
    return set(out.stdout.split("\n")) if out.returncode == 0 else set()


def check(root: Path | None = None, feed_dir: Path | None = None) -> list[Finding]:
    root = Path(root) if root else ROOT
    feed_dir = Path(feed_dir) if feed_dir else root / "data" / "products" / "solar"
    tracked = _tracked(root)
    findings: list[Finding] = []

    for feed in sorted(feed_dir.glob("*.json")):
        doc = json.loads(feed.read_text())
        for p in doc.get("products", []):
            prov = p.get("provenance", {})
            f = lambda kind, detail: Finding(  # noqa: E731
                kind, feed.name, str(p.get("element")), str(p.get("ion")),
                str(p.get("band")), str(p.get("treatment")), detail)

            ct = prov.get("copied_to")
            if ct is None:
                findings.append(f("MISSING_COPIED_TO",
                                  f"published from {prov.get('host')}:{prov.get('path')} "
                                  f"with no committed artifact"))
                continue

            ctp = Path(ct)
            if ctp.is_absolute():
                # The /private/tmp class. Resolve it: an absolute path INSIDE the repo is
                # merely unportable, one outside it is not committed at all.
                try:
                    rel = str(ctp.resolve().relative_to(root.resolve()))
                except ValueError:
                    findings.append(f("COPIED_TO_OUTSIDE_REPO",
                                      f"{ct} is not inside the repository — a non-null "
                                      f"copied_to is not evidence of being committed"))
                    continue
            else:
                rel = str(ctp)

            target = root / rel
            if not target.exists():
                findings.append(f("ARTIFACT_ABSENT", f"{rel} does not exist"))
                continue
            if rel not in tracked:
                findings.append(f("ARTIFACT_UNTRACKED",
                                  f"{rel} exists but is not git-tracked "
                                  f"(band_products is gitignored — use git add -f)"))
                continue

            expect = prov.get("sha256")
            if not expect:
                findings.append(f("NO_RECORDED_SHA", f"{rel} committed with no feed sha256"))
                continue
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            if actual != expect:
                findings.append(f("SHA_MISMATCH",
                                  f"{rel}: feed recorded {expect[:12]}, committed file is "
                                  f"{actual[:12]} — not the file the feed measured"))
                continue

            if str(prov.get("host", "")).lower() == "mac":
                findings.append(f("REGENERABILITY_GAP",
                                  f"{rel} is committed but was produced on the Mac; its "
                                  f"holding is not on the committed runner (RYA-1011)"))
    return findings


def report(findings: list[Finding]) -> str:
    if not findings:
        return "OK: every live feed product has a committed, checksum-matching artifact."
    blocking = [x for x in findings if x.blocking]
    gaps = [x for x in findings if not x.blocking]
    out = []
    if blocking:
        out.append(f"🔴 {len(blocking)} BLOCKING finding(s) — the feed publishes what the "
                   f"repo cannot back:")
        for x in blocking:
            out.append(f"   {x.kind:<24} {x.element} {x.ion} {x.band} {x.treatment}"
                       f"\n       {x.detail}")
    if gaps:
        out.append(f"\n⚠️ {len(gaps)} regenerability gap(s) — committed, but not "
                   f"reproducible on the committed runner (recorded, not silent).")
    return "\n".join(out)


def main() -> int:
    findings = check()
    print(report(findings))
    Path(ROOT / "data" / "results" / "rya1080").mkdir(parents=True, exist_ok=True)
    (ROOT / "data" / "results" / "rya1080" / "rya1080_guard.json").write_text(
        json.dumps([asdict(x) for x in findings], indent=2) + "\n")
    return 1 if any(x.blocking for x in findings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
