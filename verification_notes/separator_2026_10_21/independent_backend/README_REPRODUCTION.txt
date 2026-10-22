Separator backend audit, 21 October 2026

Purpose and limits
------------------
These are diagnostic copies, not replacements for the frozen certified package.
The original interval arithmetic formulas and fixed inputs are retained. Finite
endpoint/order checks are inserted in interval constructors; finite checks are
inserted before and after dn/up; congruence input arrays are checked. Exact-hex
sidecars record the internal bounds before their decimal serialization.
All scientific stdout and coefficient exports are compared byte for byte with
an untouched, fresh run of the original supervisor. Exact Fractions verify the
one-sided direction of every consumed exported endpoint, the baseline PAD,
and the binary chart residual dependency. This is not a second arithmetic
implementation, proof-assistant formalization, or proof for arbitrary inputs.

Prerequisites and directory layout
----------------------------------
Linux x86_64, g++ with 64-bit long-double significand, Python >=3.11.
Use the same compiler flags recorded in instrumentation.json. The manuscript's
environment preflight should be run first. In an otherwise fresh directory:

  audit_root/
    xy_global_separator/      untouched, extracted original package
    independent_backend/     this audit source directory
    fresh_results/           output of untouched original supervisor

Original scientific payload SHA256 pins are recorded by verify_separator.py.
Frozen original source files must not be patched.

Reproduction from the audit_root directory
-----------------------------------------
1. Run the original supervisor, freshly rebuilding every numerical premise:

   python3 xy_global_separator/verify_separator.py \
       xy_global_separator/distribution --output fresh_results --workers 2

2. Construct instrumented copies and fixed case inputs, then install exporters:

   python3 independent_backend/prepare_instrumented.py
   python3 independent_backend/add_exact_sidecars.py

3. Run the guarded baseline (first create results/ if absent):

   mkdir -p independent_backend/results
   cd independent_backend/baseline
   ./verify_guarded min > ../results/baseline.json 2> ../results/baseline.log
   cd ../..

4. Run all 109 guarded charts and the complete neighbor geometry:

   python3 independent_backend/run_instrumented.py

5. Build and run the two small negative-control probes using the SAME flags:

   g++ -std=c++17 -O2 -DNDEBUG -fno-fast-math -ffp-contract=off \
       -frounding-math independent_backend/probe_original.cpp \
       -o independent_backend/probe_original
   g++ -std=c++17 -O2 -DNDEBUG -fno-fast-math -ffp-contract=off \
       -frounding-math independent_backend/probe_guarded.cpp \
       -o independent_backend/probe_guarded
   independent_backend/probe_original > independent_backend/probe_original.txt
   independent_backend/probe_guarded > independent_backend/probe_guarded.txt

6. Check all exact exports and compare scientific bytes with the fresh originals:

   python3 independent_backend/check_exact_exports.py

Every exact-hex sidecar needed by the checker is copied to results/ (including
results/cNNN/coefficients_exact.hex and results/geometry_exact.hex). Raw proposal
cases are not needed when inspecting the completed reports.

The expected diagnostic status is FIXED_EXPORTS_EXACT_COMPARISON_PASS. It must
not be relabeled as a new global optimality certificate. The original complete
separator supervisor status is a separate result.

Files not needed for source distribution
---------------------------------------
The reconstructed cases/ contain raw proposals hardlinked to originals. The
built verify_guarded and probe binaries are regenerated. Do not redistribute
them as replacements for the original archive. Patches, source scripts, exact
sidecars, stdout/log files, and the compact reports are sufficient audit records.
complete_sidecars.py only completed a historical run begun before sidecars
were installed; it is not needed in the fresh sequence above.

Known generic counterexamples
----------------------------
probe_original.cpp shows that positive_by_congruence accepts a NaN matrix and
that bare precision22 output can round an upper bound downward. The first is
a defensive-programming defect; it does not establish that fixed certified
inputs produce a NaN. Guarded complete fixed-input replay addresses that
reachability question only for this execution. The second is repaired for the
separator's four consumed baseline bounds by explicit PAD=1e-18 in its exact
postprocessor; every actual PAD comparison is checked from exact binary hex.
Chart and geometry endpoints use precision30 plus an extra dn/up step. Their
actual exported directions are checked individually, without a float parser.
