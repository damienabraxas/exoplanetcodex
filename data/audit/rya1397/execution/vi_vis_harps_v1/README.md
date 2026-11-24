# Solar V I HARPS VIS component insertion v1

Executed on Sirius from pipeline revision `925caac8` using the scoped [runner](run_hfs_insertion.py). The previous 5 mÅ selector missed three published V I transitions with physically corresponding rows already in GES. The runner joins the laboratory transition to the GES group by lower and upper energy and J, checks the canonical-scaled sum of gf, and uses the strongest existing component as the fit center through a versioned copy of the input set. The atomic list and the published set stayed byte-for-byte intact. Each synthesis window contains all of its existing components. GES represents 5731 Å as a single row, so this result is conditional on that representation.

| Laboratory air Å | Fit center air Å | GES rows | 1D-LTE line abundance |
|---:|---:|---:|---:|
| 4875.486 | 4875.508 | 5 | 3.915 |
| 5731.249 | 5731.241 | 1 | 3.755 |
| 6111.651 | 6111.673 | 2 | 3.860 |

The ledger confirms exactly three added and fitted features, 27 stored fits reused in each fitting treatment, and zero fallback refits or unmatched requests. All three passed the synthesis aggregation gates and all 27 old abundances and acceptance flags match. Both 1D LTE and Engine B moved from 3.855 to 3.857 dex. The stage statistical term changed from 0.0444 to 0.0401 dex, with the stage systematic term at 0.17 dex. The mixed gf pool is ungraded. Engine A supplied zero finite departures and emitted no product. These are review outputs and no publication budget was assembled.

The [raw products, mapping and input manifest](products/), [before-and-after products](product_comparison.csv), [line comparison](line_comparison.csv), and [run log](rya1397-v-vis-hfs-v1.log) document the test. Sirius verification found no changes to the committed Solar feed or V stage. The data mount retained 73 GiB free.

Additional input checksums captured on Sirius:

- GES VIS deck `9231e3109703afa9abcabf24609ff7a669c7b60731d85d5dafdd86123b263f17`
- Canonical gf `c90e12a7397a549e2895cebe2398ec3c44d2f197bb9eab9aa48d446b8f6097de`
- HARPS spectrum `698ced10e784ec49a0fe88764496af7b3e10f4ca89456249d5cbbcf39422aefe`
