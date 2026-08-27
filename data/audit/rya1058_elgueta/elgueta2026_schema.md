# Elgueta et al. 2026 field dictionary (RYA-1058)

Authority: committed CDS `ReadMe` for `J/A+A/710/A111`. The three `atomic*.dat`
files are fixed-width 805-byte records. `Wave` (bytes 1–12) is the transition
wavelength in Å despite the CDS `0.1nm` label; `El` (14–17) is compact species
and ion (for example `FeI`); `LEP` (19–27) is lower excitation potential in eV;
`loggf` (29–37) is oscillator strength; `C6` (39–47) is van der Waals broadening.
Blank numeric and flag fields mean unavailable/not measured.

Six stellar-type blocks contain synthetic/observed depth and wavelength offset,
measurement uncertainty/warnings, then the paper's selection flags: `Depth`,
`Sat` (`Y` means unsaturated), `PurEW`, `Pur`, fit residual statistics, `GoF`,
and final `Rob`. The blocks map exactly to Procyon/F dwarf, Sun/G dwarf,
epsilon Eri/K dwarf, beta Hyi/FGK subgiant, Arcturus/FGK giant, and gamma
Sge/M giant. `Rob=Y` is the paper's final conjunction; blanks are preserved and
never recoded as passes. The raw tables provide no per-line gf reference, so
empirical robustness never establishes laboratory provenance.
