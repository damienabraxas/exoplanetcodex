# ACE-FTS non-CO molecular reach (recon only)

Atlas axis: wavenumber_cm; coverage 700.0-4430.0 cm^-1 (2.26-14.29 um).

_Windows are approximate; confirm band-heads vs HITRAN/ExoMol. Abundances inherit the CO-leg disk-center mu~1 + 1D->3D framework (RYA-444); NOT measured here._

| species | system | window cm^-1 | in_cov | npts | depth | n_lines | ~SNR |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C | CO 1-0 fundamental | 1900-2250 | yes | 70001 | 0.34 | 51359 | 2092 |
| C | CO 2-0 overtone (anchor) | 4150-4360 | yes | 42001 | 0.34 | 17250 | 272 |
| O | OH 1-0 fundamental | 2600-3600 | yes | 200001 | 0.29 | 131406 | 4449 |
| C | CH 1-0 fundamental | 2650-3100 | yes | 90001 | 0.29 | 58147 | 3932 |
| N | NH 1-0 fundamental | 3000-3500 | yes | 100001 | 0.28 | 68426 | 5343 |

**Blend flag:** OH + NH + CH all populate ~3000-3500 cm^-1 -> O/N/C disentangling there is a simultaneous-synthesis problem (Turbospectrum), never isolated EW.