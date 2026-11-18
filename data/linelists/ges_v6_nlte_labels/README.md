# `nlte_labels.tsv`: NLTE level labels for the GES v6 synthesis list

The GES list lives in the iSpec install, so its sidecar is committed here. Written by
`run_matrix._ensure_ges_labels` when GES labels fewer than half of an element's graded lines and
the element has no published feed. Only that element's cells read it, and it replaces that
element's own GES labels. See `pipeline/nlte_label_sidecar.py`.

Cr (RYA-728), measured: 0 of 27 graded Cr I and 1 of 8 graded Cr II had both levels labelled in
GES, and that one (4588.199) disagreed with VALD's terms. Sidecar: 27 of 27 and 8 of 8; Cr I 4310
of 8152 rows and Cr II 3197 of 6648 labelled, with 17 Cr II rows left unlabelled because the routes
named different levels. Routes: VALD line record, VALD level table, and RYA-818's term join. Energy
matching is never used (RYA-818: wrong level 18% of the time for Cr I, 35% for Cr II).
