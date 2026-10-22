CERTIFIED GLOBAL ANNIHILATION BARRIERS IN AN INFINITE XY LATTICE
English manuscript and verification snapshot, revision 7, 22 October 2026
Prepared for author review before a possible SIADS submission.
No journal submission, repository publication, or contact with an editor was performed by this revision.
Repository configured: https://github.com/rsumoza/global_annihilation_barriers_infinite_-XY_lattice
Planned tag: siads-v7-2026-10-22 (publication unverified).
Current distribution instructions: REPOSITORY_RELEASE_ES.txt and index.txt.
Later revision notes below are historical; the current links supersede their pending-URL statements.

READ FIRST
main.tex and sections/*.tex contain the complete article and its proofs.
SIADS_XY_Manuscript_v7_2026-10-22.pdf is the compiled review copy.
supplement/Supplementary_Reproducibility.tex and .pdf identify the frozen
sources, fresh execution, scope, and arithmetic assumptions.
supplement/INDEX.txt and INDEX.html list every proposed supplement item.
AUTHOR_NEXT_STEPS_ES.txt identifies the remaining author metadata/review.

The central result is a matched global lower and upper annihilation
barrier for a dipole of separation 12 on the infinite homogeneous square
XY lattice, with exact cosine and all strong-energy-space escape paths.
The new analytic extension treats the complete ell1 coupling ball
||b||1<=1/500. It proves critical-state continuation and a positive global
barrier with |B_b-B_0|<=2||b||1. It does not identify a perturbed minimax
with a particular continued saddle or prove a mechanism switch.
The Green proof gives effective inverse-radius tails in the same class.
The conservative finite-box budgets remain impractically large.

BUILD THE ARTICLE
From the extracted directory, run twice (or until references stabilize):

    pdflatex -interaction=nonstopmode -halt-on-error main.tex

Build the companion from its directory, again twice:

    cd supplement
    pdflatex -interaction=nonstopmode -halt-on-error Supplementary_Reproducibility.tex

The delivered PDFs use the SIADS-permitted alternative format: 10-point
type, 6in by 8in text block, mathematical line numbering. Standard TeX
Live packages are used (amsmath/amssymb/amsthm, lineno, cleveref, hyperref,
lmodern, graphicx, booktabs, array, xcolor; xurl and microtype for the
companion). Main references are in references.tex; no BibTeX run is needed.

OPTIONAL OFFICIAL SIAM TEMPLATE
main.tex automatically uses siamonline250211.cls when the class is on
the TeX input path. Install the COMPLETE official distribution from:

https://epubs.siam.org/pb-assets/macros/online/siamonline_250211.zip
https://epubs.siam.org/pb-assets/macros/online/docsiamonline.pdf

The class and all its dependencies must be installed; algorithm.sty is
one such dependency (usually supplied by a TeX Live science collection).
Do not use a stub to prepare a submission PDF. No isolated SIAM class is
redistributed here, because the class license requires its complete
distribution. The conditional source makes no changes to the SIAM class
font sizes or margins. Recheck pagination/references after switching.

REPRODUCE THE NEW ARITHMETIC
Only Python's standard library is required for the rational checks. The
environment preflight additionally invokes g++. From this directory, run:

    python3 checks/verify_supported_environment.py --output checks/supported_environment.json
    python3 checks/verify_robustness_constants.py
    python3 checks/verify_green_tail_budget.py
    python3 checks/verify_infinite_class.py --self-test --output checks/infinite_class_report.json

Use paths relative to this directory, or absolute paths. The robustness
checker defaults to this snapshot's proof_record directory and checks
subsidiary hashes against the global certificate, and its global hash
against the saved reproduction record. To check records from a different
fresh run, pass --certificate-root /absolute/path/results_fresh and inspect
that run's own integrity record. A different platform can produce harmless
format/build metadata differences; a status is not a proof by itself.

The Green checker verifies the exact rational constants and dyadic budget
for the original path, plus the uniform perturbed-state constants. It
does not prove the analytic kernel estimate, solve new equilibria,
materialize a huge box, or rerun the old certificate. The full analytic
proofs are in the article and appendices.

PAYLOAD/NORM CHECK
The additional verify_path_inputs.py audits the original binary arrays
and rational norm estimates. With the two frozen archives accessible:

    python3 checks/verify_path_inputs.py \
      --connection-zip /path/XY_Conexion_Global_Reproducible_2026-10-19.zip \
      --optimality-zip /path/XY_Optimalidad_Extension_Reproducible_2026-10-19.zip \
      --output proof_record/path_audit.json

This input audit consumes authenticated geometry/correction enclosures;
it does not replace their separate full-space interval verification.

REPRODUCE THE GLOBAL CERTIFICATE
S1 (this small source bundle) is NOT the complete global input package.
The three original archives S2-S4 in supplement/INDEX.txt are required.
Verify their exact sizes and SHA-256, then extract alongside each other.
Run checks/verify_supported_environment.py from S1 first. This reviewed
reproduction target requires a 64-bit long-double significand, binary64
little-endian input and round-to-nearest; wider formats are not validated
by the frozen source guards that merely accept at least 64 bits.
From xy_optimality_extension, with absent/empty results_fresh, run:

    python3 verify_optimality.py \
      --separator-package ../xy_global_separator \
      --connection-package ../xy_global_connection \
      --output results_fresh --workers 2

Do not choose options that reuse validated outputs for this reproduction.
The new run on 21 October took about 12m17s. Its summary, environment,
stdout/stderr, full final certificate and all 548 fresh output files
are under proof_record. The preceding 20 October record remains in its
history directory. The complete original source/proposal
archives remain separate indexed supplement items, rather than being
silently omitted or replaced by result JSON files. Their names/hashes are
also recorded in proof_record/reproduction_summary.json.

This is a fresh reproduction of the SAME software with fixed proposals.
It is not an independent interval implementation, proof-assistant proof,
fresh candidate search, complete historical attack run, or fresh execution
of the separate heteroclinic extension. The analytic proofs and interval
arithmetic assumptions remain part of the trust boundary.

CONTENTS AND INTEGRITY
MANIFEST.sha256 lists every included source and record, other than the
manifest itself. The ZIP's own hash is provided externally. No third-party
class, downloaded helper, hidden credential, cached binary, or unrelated
historical model is included. Original package licenses and author rights
are not changed by this snapshot.

The supplement index and author note are editable. The supplied affiliation and
contact information have been inserted. Funding, author responsibility and
final review must
be confirmed by the human author before submission. SIAM's disclosure
policy is reflected in the Declarations; it is not a statement that
the author's personal review has already occurred.

REVISION 2: RESPONSE TO SUPPLIED EXTERNAL COMMENTS
The revised article adds the Dirichlet/Rayleigh dissipation derivation,
a two-panel conceptual geometry figure, shorter abstract sentences,
the radius-versus-search-complexity distinction, explicit future kernel
cancellations, and pointwise bounds for dependent bond angles. No
certified barrier interval, critical-point claim, or coupling budget was
changed. No public archival deposit or journal submission was made.

Rebuild the vector/preview schematic with:
    python3 figures/make_certificate_geometry.py
This uses NumPy and Matplotlib and draws conceptual geometry only; it
does not read or display calibrated certificate coordinates. Caption
and alt text are in figures/certificate_geometry_caption.txt; alt text
is also embedded as comments in the article's figure source.

The separately delivered SIADS_XY_Response_to_Review_2026-10-20.pdf
explains every comment and policy clarification. Its editable source
Response_to_Review.tex accompanies this source bundle. It is a working
response to the supplied pre-submission review, not a journal decision.
The proposed permanent public deposit remains an author action; a
separate DOI-before-acceptance mandate was not found in the official
SIADS/SIAM instructions checked. Verify and insert the actual deposit
identifier only when one exists.

REVISION 3: TARGETED INSPECTION AND MATHEMATICAL PRECISIONS
The source now defines the energy space by the closed gradient subspace,
formalizes closed-cut boundary retention and tangencies, and displays a
rational cap coefficient of at least 19/40000. The closed Newton ball was
already present and is retained explicitly. The abstract identifies the
fixed separation and saddle-minus-minimum energy; the bibliography has
three additions and two corrected records. Actual author metadata and
typed cross-references are included.

A targeted source audit and a conditional geometry replay are recorded in
verification_notes/current_targeted. These DO NOT constitute a new full
global execution: separation and connection raw inputs were unavailable
for this inspection. Original optimality files and scientific pins were
not modified. General NaN/serialization limitations of helper/historical
paths are documented rather than silently characterized as resolved.

The Green rational checker now uses explicit checks under python -O and
rejects nonpositive tolerances. The new environment preflight restricts
usage to the reviewed long-double64 target; it is not a CAP verifier.
The gate tests, reports and limitations are in CURRENT_AUDIT_ES.txt.

REVISION 4: SEPARATOR AUDIT AND FORCE-BALANCE CONTEXT
The unchanged standalone separator was rebuilt from scratch on 21 October
2026 and passed. Independent integer/rational code verifies all 715
geometric inclusions and three maxima. The new records and diagnostic
guarded replay are under verification_notes/separator_2026_10_21. This
additional execution covers the separator, not a new complete optimality
chain. The connection attachment was unavailable at that stage; revision 5
below resolves this historical input obligation.

All 17 prior bibliography entries were checked against primary metadata;
Braun--Hudson arXiv:2609.39642v1 is added as reference 18. Section 5 makes
the exact source cancellation explicit. The barrier, coupling budget and
critical-point statements have not changed. CURRENT_AUDIT_ES.txt is the
authoritative current assessment; older notes remain dated historical
records. The updated main PDF and supplementary PDF were rebuilt and
visually reviewed. No original scientific archive or pin was changed.

REVISION 5: COMPLETE REPRODUCTION AND CONNECTION AUDIT
All five connection volumes were reconstructed and matched the pinned
441730856-byte S3 original. checks/reconstruct_connection.py validates
each volume, ZIP CRC and final archive identity before atomic replacement.
The complete unchanged three-package chain was rebuilt in a new directory
and passed with FRESH_OPTIMAL_ANNIHILATION_BARRIER_PASS in 737.332688 seconds.
All 548 global output hashes and nested output hashes were checked.
The barrier and scientific scope are unchanged.

Appendix D proves exact finite-load mean zero, the Neumann comparison,
and both explicit exterior residual terms without a symmetry assumption
on the correction. The actual path payload audit reproduces M=6007.
Independent exact normal/adjoint geometry supplements the earlier exact
atlas geometry. Instrumented fixed-data executions, exact scalar export
checks, negative controls and source bindings are documented under
verification_notes/connection_2026_10_21. They are diagnostic copies, not
silent patches to the frozen originals or a second interval engine.
The generic congruence helper's injected-NaN limitation remains disclosed.

Use CURRENT_AUDIT_ES.txt for the final current balance. Earlier audit
notes describe the evidence available at their dates, including temporary
missing-input states that no longer apply. The bibliography retains the
18 primary-source-checked entries from revision 4; no universal novelty
or exhaustive-literature claim is added by this computational audit.

REVISION 6: EXPLICIT INFINITE CLASS AND GITHUB RELEASE PREPARATION
The main fixed-instance global certificate and its original sources are
unchanged. A new section proves inheritance of its communication height
for a full X-ball of radius 1/20 around the certified minimum. A normalized
block format defines a noncompact infinite class with an effective tail.
The arbitrary-accuracy barrier enclosure terminates for those inputs.
This is reuse of the existing global certificate, not discovery of new
mechanisms or a complete algorithm for all configurations. The root-name
and general enclosure construction is proved; its arbitrary-input
software interface is not implemented. The new exact rational script
checks constants, finite prefixes, interface bonds and negative controls.

README.md and REPOSITORY_RELEASE_ES.txt describe the ready-to-publish
repository layout. Original archives S2-S4 must be GitHub Release assets,
not ordinary Git objects. repository_state.json deliberately has null
URL/tag/commit until the author publishes the actual release. The
supplement and index still accompany the initial SIADS submission.
External access is not assumed to waive the journal's review requirements.

REVISION 7: AUTHOR-SUPPLIED REPOSITORY
https://github.com/rsumoza/global_annihilation_barriers_infinite_-XY_lattice
Planned tagged release: https://github.com/rsumoza/global_annihilation_barriers_infinite_-XY_lattice/releases/tag/siads-v7-2026-10-22
Documentation and availability only; mathematics and frozen inputs unchanged.
Remote publication and source commit remain unverified.
