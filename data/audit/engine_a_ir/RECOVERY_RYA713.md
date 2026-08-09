# Engine-A label recovery — 14 elements unblocked — RYA-713

## Before

| verdict | n |
|---|---|
| IR-CAPABLE-NOW | **1** — Ti only |
| IR-CAPABLE-BLOCKED-ON-LABEL | **14** |

Fourteen elements had a multi-GB departure-coefficient grid on disk that could not be
pointed at any line, because `label_{El}.txt` — the 9–64 kB file mapping level index to
species/configuration/term/J/energy — was never extracted. RYA-477's fetch (2026-06-29)
kept only the `.grd` and freed the archive; RYA-545's later Ti fetch happened to keep the
text files, which is the only reason Ti worked.

## After

| verdict | n | elements |
|---|---|---|
| **IR-CAPABLE-NOW** | **13** | Al, Ba, C, Ca, K, Li, Mg, Mn, N, Na, O, Si, Ti |
| IR-CAPABLE-BLOCKED-ON-LABEL | 1 | **Cu** — different upstream (Caliskan 2024), not in Zenodo 3982506 |
| IR-CAPABLE-BLOCKED-ON-LINELIST | 1 | **H** — grid and label staged, but the GES level-identified list carries no H |
| PER-LINE-TABLE-ONLY | 3 | Fe, Cr, Sr — Bergemann MPIA serves per-*line* corrections upstream |
| NO-ENGINE-A | 10 | Co, Eu, Ni, P, S, Sc, V, Y, Zn, Zr |

**13× improvement for the cost of a download that was deleted afterwards.**

## Method

`_label_recovery_rya713.log`. Per element: fetch the `pysme` tarball from Zenodo 3982506,
**verify md5 against the record**, extract only `label_*.txt` and `atmos_*.txt`, delete the
archive. Every `.grd` was already held, so nothing large was re-downloaded to keep.

* Transient traffic: ~5.3 GB. Peak transient disk: one tarball (1.4 GB, Mn).
* Retained: 14 × ~9–64 kB.
* Disk after: `/mnt/codex-data` 130 GB free (was 134) — the difference is elsewhere, not this.
* Every element md5-verified before extraction; a mismatch would have skipped that element
  rather than staging an unverified file.

## The two that remain, and why they are different

**Cu** is a genuine gap: its grid came from Caliskan 2024, a different upstream, so the
Zenodo 3982506 recovery could not reach it. Same fix, different source.

**H** is not a label problem at all — grid and label are both staged. The GES
level-identified linelist simply carries no hydrogen, so there is no line to attach the
departure coefficients to. That is a *linelist* gap and it is worth noting against the
standing Balmer-Teff question: we use H lines for normalisation and RV but have never
derived an independent Teff from them.

## What this does and does not unlock

**Does:** 13 elements can now have Engine-A departure corrections computed for *any*
transition whose two levels are in the model atom — including the IR and the near-UV. The
grid never had a wavelength restriction; we simply could not address it.

**Does not:** it does not supply oscillator strengths. The Fe I IR pass found 123 of 174
failures were atomic-data faults and **0 of 271 passing lines carry a NIST grade**. A
departure coefficient corrects a line's *physics*; it cannot rescue a line whose *gf* is
unknown to 0.1–0.3 dex. Engine-A reach and line-list quality are independent problems and
only the first just moved.
