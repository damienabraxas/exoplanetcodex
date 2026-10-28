# Reverse charge-transfer contract: first verification step

Continuation of the native-charge review at f3aa992c. This is a Mac LTE-reservoir
inspection, not a new RT/SE solve, rate adoption or producer identification.
Ryan's sequence is oxygen, then carbon and nitrogen (CNO milestone), then Fe,
then another star. The ordered work list is recorded on RYA-1281 and RYA-1278.

## What the held files establish

The unchanged atom.o41f has 24 CH0 records: first indices 6 through 29, second
index 40. These connect neutral oxygen excited states to the O II ground state.
The parser preserves raw first/second indices as MultiCollision.j/i; the native
adapter translates them to source/destination explicitly. The held GENCL source
at pinned pymulti 4d6aa3c, mul23_opacu.f lines 4721-4738, explicitly describes
input-order direction and multiplies CH0 by nh(1); CH+ uses nh(6). Neither branch
constructs a reverse rate. The atom has no CH+ table and no identified reverse
CH0 rows. These are source-code facts, not identification of the modified grid
producer or the underlying quantum collision calculation.

The reversible adapter currently adds the per-oxygen reverse rate

    C_reverse = q_forward * n(H I ground) * nstar(O I i)/nstar(O II k).

Its units are s^-1 after q is converted from cm^3/s to m^3/s. It is an effective
rate which closes local LTE balance by construction. It is not itself the
microscopic reverse coefficient (m^3/s). A specific reaction requires a specified
reverse collider. For the ion-pair hypothesis

    O I(i) + H(1s) <-> O II(k) + H-,

one would infer, under LTE reservoirs,

    q_reverse = q_forward * nstar(O I i)*nstar(H 1s)
                            / [nstar(O II k)*nstar(H-)].

This algebra does not establish that the held coefficients describe that
reaction, that their source has the same level/statistical-weight conventions,
or that the producer applied this reconstruction. A proton collider cannot be
substituted for H- merely because both are called charge transfer.

For independent physical context, Belyaev et al. 2019, Sections 1-2 and Table 1,
explicitly distinguish O+ + H- and O- + H+ channels and discuss inverse mutual
neutralization: https://doi.org/10.1093/mnras/stz1511 . This is corroboration of
the reaction distinction, NOT provenance or numerical validation of atom.o41f.
No coefficients from that paper were ingested. Barklem 2018 publisher full text
was not accessible in this check; search metadata is not used to assign rates.

## Bounded reservoir probe

scripts/audit_bride_charge_reservoir.py reconstructs the same hash-pinned native
5750/4.5 LTE context and inspects H ground, H+, H- and the oxygen LTE ratio.
It independently evaluates ground-H Saha attachment using the installed engine's
constants, separately from Lightweaver's existing Hminus field (which uses summed
H populations in the inspected implementation). It does not change either field.
The reversible adapter's LTE-H requirement remains in force.

For the lower triplet state (one-based level 6):

| logtau500 | T (K) | H- LW / ground-Saha H- | inferred qreverse/qforward if H- |
|---|---:|---:|---:|
| -4.9177 | 4045.9 | 1.00000018 | 106926 |
| -3.0030 | 4483.7 | 1.00000083 | 37727 |
| -1.2621 | 4964.0 | 1.00000334 | 14852 |
| -0.0102 | 6423.3 | 1.00038068 | 2054 |

The two H- evaluations differ negligibly at these sampled outer/photospheric
points. This probe does not support blaming that particular distinction for the
large outer-population mismatch. It also makes clear why a per-particle reverse
rate and a microscopic reverse coefficient must not be confused. Ratios are
conditional calculations, not measured or independently sourced collision data.

Original exploratory probe/report remain immutable at
CODEX_WORK_ROOT/bride_mac_checks/rya1281-charge-reservoir-v1. The committed CLI's
separate run is rya1281-charge-reservoir-v2/report.json; it includes source hashes,
probe hash, all 24 raw endpoint indices and numerical collider densities.

    source config/hosts/mac.env
    "$CODEX_BRIDE_EVAL/.venv/bin/python" scripts/audit_bride_charge_reservoir.py \
      --reference "$CODEX_WORK_ROOT/bride_inputs/rya1281/reference-5750-4.5-A8.70-v1.json" \
      --output "$CODEX_WORK_ROOT/bride_mac_checks/<new-immutable-directory>"

## Exact next verification

Locate the original coefficient provenance and reverse-process prescription;
compare a microscopic reverse table using matched level/statistical weights and
hydrogen attachment convention. Until then preserve the effective LTE-restoring
pair as a hypothesis. Next outer-layer controls must retain that label and vary
radiation/background/boundary assumptions separately. No abundance fitting,
source changes, default-converter change, production emission or admission.
