# RYA-1219 CRIRES+ Vesta J/K telluric correction

This audit records the completed CRIRES+ Vesta J and K telluric-correction pass. The four J settings (J1226, J1228, two J1232 exposures) and four K settings (K2148, K2166, K2192, K2217) were processed as full-arm, per-chip products with `molecfit` and the matching GDAS profile for each observing night. Each product contains a non-unity `MTRANS` extension and `TELLAPP=True`; the product flux was continuum-normalized per segment after correction.

All eight D1 residual gates pass. The before/after residuals are recorded in `corrected_products_manifest.csv`; after values are 0.01394--0.01914. This is a telluric-correction result and does not require CNO line validation or reflected-solar radial-velocity conditioning. Rest-frame conditioning remains a separate downstream step, so consumers must not interpret these products as abundance-ready solely from this ticket.

K2166 required a controlled retry because the optional well-mixed CO refit was numerically unstable. The retained product uses the successful broad-band molecfit model while bypassing only that failed refit; this is recorded in the manifest and should be inspected when selecting CO lines. K2148 and K2166 do not place the CO bandhead on a detector chip, so a CO-specific bandhead diagnostic is unavailable for those settings; this does not block correction of the telluric absorption present in their recorded orders.

The raw inventory remains in `idp_inventory.csv`. Corrected product paths, SHA256 checksums, MTRANS evidence, GDAS profiles, and gate values are in `corrected_products_manifest.csv`. The products themselves are staged under `data/results/rya1219_crires_products/{J,K}/` (FITS are ignored by the repository; the manifest is the durable record).
