# Contributing to the Spectroscopic Frontier atlas

Contribute a bounded source/diagnostic entry or a correction through a reviewed PR.
A proposed source is not adopted science. No direct production publication, protected
FITS/gold changes or recomputation of abundances is part of atlas maintenance.

1. Search existing atlas IDs and canonical source registries; reuse them. Name the
   species, ion/isotopologue and component type. Never infer identity from wavelength
   alone. For molecular transitions follow the existing molecular-transition schema;
   for atomic lines preserve source-qualified levels, components and origin.
2. Inspect the primary release, article and exact code/data version. Add a source
   pointer and bibliography key, not a second author/DOI table. Record attempted and
   failed access separately; `not_yet_searched` is preferable to a fabricated negative.
3. State what was searched, when, which versions/diagnostics/parameter domain were
   covered, what remains unknown and why the resource matters. A database landing
   page is discovery evidence, not a completed release audit. Do not extrapolate a
   result to an entire element, band or stellar class.
4. Document the physics separately: model atom/rates; atmosphere/EOS/geometry;
   solver versus formal synthesis; LTE/non-LTE and dimensionality; partition and
   dissociation data; line strength/grade origin; magnetic/broadening applicability;
   and observational continuum/telluric/LSF conditioning. Nulls and unassessed
   constraints must remain explicit. Do not create default cross-resource mixing.
5. Check license/attribution against the exact upstream version. Pending means
   unresolved; public access is not a permissive license. Preserve recommended
   software and original data citations. Link issues for gaps with a next gate;
   classify priority from scientific impact, not a completion percentage.
6. Edit `config/frontier/atlas.v1.json`. When advancing the science base, inspect the
   corresponding main revision and all changed source receipts first. Regenerate
   JSON and Markdown, run the checks, and review the diff for unintended metadata
   changes. Do not hand-edit generated counts, axis values, transitions or hashes.
7. Obtain independent scientific review of scope, identity, domain, uncertainties,
   provenance and claims, separately from code/schema review. This draft does not
   support scientific adoption; any validated release needs a new reviewed contract
   and explicit sign-off evidence. A passing schema or numerical control is not
   scientific certification. Ryan authorizes merge and publication.

Corrections retain the previous version in Git history and identify affected entry
IDs, source versions, cause and impact. If a scientific product relied on a wrong
entry, route its disposition to that product's owning workflow; never silently
rewrite a released result from an atlas update. Freeze/hash a reviewed atlas release
before public consumption. Update the bibliography and canonical owner when the
source error belongs there, rather than using this atlas as a permanent override.

Future website search should filter component, species/isotopologue, band,
parameter scope, treatment, availability and validation independently, with visible
provenance, audit dates, gaps and accessible keyboard controls. HTML implementation
waits for RYA-1370 review/approval; no private API/runtime secret is necessary.
