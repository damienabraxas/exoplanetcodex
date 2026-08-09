# The optical control FAILED — the EW method does not reproduce a known answer

Ryan: *"we want what we know to be close right? Fe in Visible should be pretty dang close to
Asplund"* · *"in IR and UV it is the wild wild west"* · *"Is the IR Uncertainty a systematic
in that region? or the instrument? etc"*

**The control was run and it failed. That is the control working.** No IR abundance was
reported on the back of this method, because the optical leg was tested first.

## The test

Same harness, same code path, same continuum policy, run on Fe I 3924–6905 Å — the band
where the answer is known (Asplund 7.46; our banked 1D-NLTE anchor 7.466, RYA-553, on 808
HARPS lines). 146 lines measured by **both** the HARPS pool and this harness on Kitt Peak.

| | |
|---|---|
| median Kitt Peak / HARPS EW ratio | **0.773** |
| implied abundance shift | **−0.112 dex** |
| MAD of the ratio | **0.490** |
| range | **0.03 – 5.04** |
| within ±10 % | **22 of 146 (15 %)** |

A −0.11 dex offset against a known anchor, with a 5× spread, is a failed control. Nothing
built on these EWs is reportable.

## Answering the question directly: region, instrument, or method?

**Not the region.** The ratio does trend with wavelength (0.52 → 0.70 → 0.88, blue to red),
but that is a *consequence*, not a cause — see below.

**Not the instrument.** An instrumental or calibration difference is a scale factor: one
ratio, tight scatter. Observed is a **5× spread** (0.03–5.04). No calibration does that.

**It is the METHOD, and specifically my window rule.**

| window half-width | median KP/HARPS | n | mechanism |
|---|---|---|---|
| crowded, < 0.20 Å | **0.294** | 69 | line **wings clipped** — most of the EW is outside the window |
| medium, 0.20–0.35 Å | 0.809 | 28 | |
| isolated, > 0.35 Å | **1.088** | 49 | window wider than the line — extra flux included |

`corr(half-width, ratio) = +0.491`.

* **93 under-measured** lines: median window **0.165 Å**, median HARPS EW **87 mÅ**. Strong
  lines, narrow windows. A strong line has broad damping wings; a 0.165 Å half-width cuts
  them off and the EW is lost.
* **31 over-measured** lines: median window **0.450 Å**, median HARPS EW **15.6 mÅ**. Weak
  lines, wide windows.

And the wavelength trend resolves: **the blue is more crowded**, so windows there are
narrower, so more wing is clipped. Region is the proxy; crowding is the variable; the window
rule is the fault.

## The root error

`window_half_width()` sets the window to half the distance to the nearest catalogued line.
That is a **blend-avoidance** rule, and it is the correct instinct for *deciding whether a
line is measurable*. It is the wrong rule for *measuring* one: an equivalent width needs the
whole profile, wings included, and a strong line's wings routinely extend past the midpoint
to its neighbour.

The two requirements are in direct conflict — you cannot both exclude the neighbour and
include the wings by choosing an interval. **That is why the reference method fits profiles
rather than integrating intervals** (RYA-429's EW fitter): a fit models the neighbour and
the wings simultaneously instead of drawing a line between them.

## What survives and what does not

**Does NOT survive — the EW values.** The 271 IR "PASS" numbers are *window-integrated
absorption*, which is a different quantity from an equivalent width. They must not be
inverted to abundances, and the earlier IR EW figures in this directory carry the same
defect.

**Survives — everything that is about the spectrum rather than the integral:**

* the **feature verification** (is the line present? at the catalogued position? is the depth
  consistent with the parameters?) — these read the profile, not the integral
* the **root-cause attribution** — 123 of 174 atomic-data faults stands; a ghost is a ghost
  regardless of how its EW was computed
* the **appendix panels** — they plot the spectrum, and now also document this failure
* the **continuum policy** finding (Kitt Peak is pre-normalised; a second normalisation
  biased EWs by a median −11.7 %) — independent and still correct
* the **Engine A/B comparison** — it never used these EWs; both engines ran synthesis on the
  same lines

## What this changes for the frontier bands

This is the concrete answer to *"is the IR uncertainty systematic in that region?"* — **no.**
It is a method systematic that operates in every band, and it happens to be *worse* in
crowded spectrum. The IR is less crowded than the blue, so the IR would have looked
*better* while being wrong in exactly the same way — and with no reference value to catch it.

**That is precisely the wild-west risk, and it is why the optical control is a precondition
rather than a courtesy.** A frontier number carrying an honest wide error bar is good
science. A frontier number carrying a −0.11 dex method bias and an honest wide error bar is
just a wrong answer with good manners.

## Owed

1. **Replace window integration with profile fitting** — Gaussian/Voigt with a blend-aware
   model, reusing the RYA-429 fitter rather than a second implementation.
2. **Re-run the optical control** and require it to reproduce the anchor before any frontier
   band is measured again.
3. **Then** re-measure the IR, and carry the *measured optical residual* into the frontier
   error budget as the harness's own systematic term rather than assuming it is zero.
