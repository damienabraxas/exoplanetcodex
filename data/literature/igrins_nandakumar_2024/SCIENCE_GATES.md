# RYA-1056 science gates and phase status

## Acquisition result

The complete CDS catalogue is preserved: `ReadMe`, `stars.dat`, and 21 element
tables (23 files). The normalized product contains 76 wavelength diagnostics,
3,800 star/line cells, and 3,426 published measurements. There are 59 nearby
canonical wavelength candidates and 17 without one. All 76 remain discovery
records, not confirmed physical-transition matches.

## Fe gate

The CDS catalogue contains no Fe table and therefore provides **zero direct Fe
transitions** for the requested Fe band-gap experiment. The paper's associated
synthesis line list is required before any Fe baseline/candidate/accepted-set
or influence analysis is scientifically possible. No Fe line or gf threshold
was changed.

## Al gate

All seven published Al wavelengths are in
`igrins_al_completeness_audit.csv`. A nearby RYA-1001 census wavelength is only
a candidate identity: CDS supplies no EP, level pair, log(gf), gf reference,
or HFS components. Rows therefore remain either
`WAVELENGTH_CANDIDATE_ALREADY_IN_CENSUS_HOLD_IDENTITY` or
`NEW_DISCOVERY_HOLD_IDENTITY_GF_HFS`. The Burheim hierarchy is unchanged.

## CNO gate

The 21 CDS atomic-element tables contain no C, N, or O table. Molecular OH/CO
usage in the stellar-parameter method cannot be reconstructed from this
catalogue. Existing Codex C/O synthesis products, the N curation debt, the
vendored CH/CO/CN/NH/OH lists, and the Amarsi 2019 CNO NLTE grids are preserved
unchanged. A modern transition-level CNO audit needs separate primary atomic
and molecular source acquisition; atomic and molecular results must remain
separate.

## Promotion rule

No row produced here can enter GF-LAB or an evaluated rung. Promotion requires
species/ion plus lower/upper transition identity, authoritative gf provenance,
published uncertainty where available, and HFS/isotopic adjudication, followed
by CRIRES+ coverage, blend, and telluric gates. IR abundance reporting remains
blocked until telluric correction is verified.
