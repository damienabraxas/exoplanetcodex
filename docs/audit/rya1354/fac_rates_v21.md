# FAC minimal-model rate experiment v21

Scientific HOLD; EXPOSED_NONBLIND. No accepted lifetime, uncertainty, covariance, calibrated sigma or grade.

The first native run v20 stopped after E1 because our driver requested separate database files per operator. FAC permits one filename per database type/session. The failed attempt is retained in `fac_rates_failed_v20` (99 native paths, 531740 compressed bytes). The correction uses a single database with six explicit operator blocks and one final print. No physical or numerical settings changed; v21 ran cold.

All four v21 cases completed on Sirius: coherent10/coherent100, Babushkin/Coulomb, fully relativistic mode0, zero explicit table cutoffs. Each queried all13 total-J-allowed channels through E1/M2/E3 to the five lower even states and M1/E2/M3 to the lower quintet. Both numerical guards passed and regenerated C/CM bytes exactly match the previously qualified states. The archive retains136 paths in one770433-byte part, including original/derived frontend, exact diff, build logs, source documents, native raw queries and standard tables. Original FAC source/library/executable remain unchanged; the isolated diagnostic frontend calls the same native library before writer pruning.

Independent replay verifies exact declared input, mapped-state bytes, numerical diagnostics and raw-query completeness. Each case has11 finite native rates plus two E3 native refusals (status-1, angular_nz0, NaN amplitudes/rates) to levels0 and3. These remain unresolved, not zero-filled. Totals and lifetimes remain null. Four targeted tests pass, including missing channels, wrong gauge and zero-filled refusal rejection.

Resume: audit pinned TRMultipole/AngularZMix selection rules against the exact minimal basis to disposition E3 refusals; independently check standard-table float32 amplitude quantization and printed rates against raw library outputs; verify full one-photon channel accounting before deriving any exploratory model total. Keep both gauges and the minimal-model versus correlated-GRASP distinction. This checkpoint does not establish physical accuracy or permit line grading/adoption.
