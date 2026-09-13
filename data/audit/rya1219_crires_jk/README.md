# RYA-1219 processing record

The Vesta CRIRES+ inventory was rechecked from the staged IDPs. It contains
three Y, four J, seven H, and four K exposures. This ticket processes only J/K;
the pre-existing corrected Y/H products are separate holdings.

The validated K pass used the existing `pipeline.crires_telluric.condition_co_arm`
path with ESO molecfit and the observation-night GDAS profiles. K2192 and K2217
contained the 2.3-µm CO order and produced corrected, continuum-normalized
per-frame products. The D1 residuals were 0.0066 and 0.0101. The reflected-solar
RV gate did not reach the required clean-line count, so the products remain
topocentric/provisional and no rest-frame coadd is claimed.

The J full-frame trial reached molecfit calctrans but did not complete within the
available local run window; it was stopped without registering a product. No
J/K holding was changed to `telluric_applied=applied`. A complete J/K run must
finish all orders, pass the per-line CNO checks, and resolve the reflected-solar
rest-frame gate before registration.

The raw inventory and K convergence evidence are committed in this directory.
