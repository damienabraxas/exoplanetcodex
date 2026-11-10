# Separate CO source-coordinate check

Predeclared supplement to the CH/CN/C2/NH/OH automatic-bsyn census, RYA1281/1278,
2026-11-07. The held canonical CO manifest declares its conversion from ExoMol
vacuum wavenumber to vacuum A. The legacy file is CO_IR_Li2015.dat; inspected
iSpec wrapper automatic routing globs *.bsyn. Neither fact establishes the
behavior of dedicated production CO code. No global ORCH defect claimed.

Before evaluating: verify installed/vendor equality and all declared117783rows,
preserve canonical IDs/status, compare raw-coordinate window membership with
explicit canonical vacuum-to-air membership in7763-7783A. Do not append CO to
the previous651AIR-consumer rows without conversion, promote its metadata or
compute absorption/abundances. Same scientific HOLD and source protection rules.
Command: project Mac interpreter,1thread,
`scripts/audit_bride_co_source.py --output <fresh-co-source-directory>`.

## Results

Initial supplement code/plan10204c87, source-ID membership extension31d4b33f
precommitted before fresh v2. Both runs exit0. All117783declared source rows
parse; installed and vendor files byte-identical.192rows lie in the numeric
7763-7783interval before conversion,191in that AIR interval after explicit
conversion.171source IDs are shared, symmetric difference41. Similar counts
therefore do not establish the same selected rows. This compares membership of
the same source IDs, not a wavelength-only transition crossmatch.

CO.dat does not match the inspected automatic *.bsyn glob. This identifies only
that wrapper route, not a production CO failure or missing CO chemical depletion.
No dedicated CO processing route audited or modified here. Canonical intake
status/probability/source metadata preserved; no strength/EOS/abundance admission.

Original v1 preserved in `experiments.bride_co_source_first`; accepted extended
v2at `experiments.bride_co_source`. All191data rows.csv bytes are identical v1/v2;
v2adds membership summary, not altered transitions/coordinates. Report/rows copied
byte-exact to `data/results/rya1281/missing_opacity_v1/co_source`; actual generator
registered. Separate191COrows must not be described as part of the651automatic
CH/CN/C2/NH/OHconsumer-AIR rows without this explicit source/frame distinction.
No new RT/source overwrite/production or adopted result.
