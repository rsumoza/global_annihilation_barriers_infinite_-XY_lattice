GUARDED REPLAY OF THE EXACT-CAP EXTENSION — 21 OCTOBER 2026

Scope
-----
Four C++ executables are reconstructed: projection, center_distances,
adjoint and chart. The first three run for both reflections; the chart
executable runs for the five enlarged charts c042,c043,c044,c053,c055.
The eleven executions use the original interval formulas and fixed inputs,
with added finite/order checks and exact hexadecimal export sidecars.
This is not a second interval implementation and not a new global theorem.

Frozen inputs
-------------
The three complete extracted source packages are required. Originals are
never modified. prepare_guarded.py verifies their pinned scientific
manifest identities and each payload hash before creating an independent
work/ directory. Data reconstruction follows verify_optimality.py exactly.
Compiled binaries and raw reconstructed work/ arrays are not distributed
with this compact audit record.

Reproduce
---------
First use the manuscript's supported-environment preflight: IEEE binary64,
little-endian, nearest rounding and a 64-bit long-double significand.
Pass the actual package paths explicitly:

 python3 prepare_guarded.py \
   --extension /path/to/xy_optimality_extension \
   --separator /path/to/xy_global_separator \
   --connection /path/to/xy_global_connection

 python3 run_guarded.py

After a successful untouched fresh full supervisor run:

 python3 check_exports.py --clean /path/to/full_fresh

The separate independent normal/adjoint-coordinate check can be run with:

 python3 check_normal_exact.py

It uses exact dyadics and an integer degree-four Laplacian, not the C++
interval engine. It proves G*c=n from the rational numerators and verifies
all eight reported normal/adjoint coordinate enclosures by rational square
comparisons. Only that fixed geometric identity is independently checked;
cap coercivity, energy and residual bounds remain separate obligations.

Small export-interface negative controls (also valid under Python -O):

 python3 negative_controls.py

probe_guarded.cpp additionally tests NaN matrix and reversed interval
rejection by the instrumented C++ constructors. Its build flags must be
the same as the other guarded programs.

The full supervisor remains a separate required execution. This component
runner never calls it, never reruns the 4096-side saddle proof and never
accepts archived JSON as a premise for a new full-space certificate.

Evidence and serialization
--------------------------
Every jout/ji interval export is recorded in order before decimal-specific
padding, as exact hexadecimal long-double endpoints. The checker recursively
binds these records to their printed JSON intervals and compares the endpoint
directions with Fraction. Every chart coefficient interval and four one-sided
header bounds are recorded and checked separately. Scientific JSON and
coefficient files must match the untouched comparison run byte for byte.

The compiled rho, candidate-complement floor, compressed floor m and chart
residual constants are checked as their actual binary values, not as decimal
spellings. The tight outer LEVEL versus exact rational saddle-upper comparison
is retained without floating-point conversion; its float display is explicitly
informational only.

The checker can also perform an explicitly labeled archived comparison:

 python3 check_exports.py \
   --clean /path/to/xy_optimality_extension/results_provided \
   --separator-results /path/to/xy_optimality_extension/antecedent_results/separator \
   --connection-results /path/to/xy_optimality_extension/antecedent_results/connection \
   --label archived --output results/exact_exports_archived_report.json

This mode is a preliminary regression check. It is not a fresh global run.

Files
-----
instrumented/: guarded sources plus complete diffs against original sources.
input_bindings.json: frozen payload hashes and package identities.
instrumentation.json: original/guarded source hashes and compiler flags.
results/: eleven scientific outputs, sidecars, logs and summary records.
prepare_guarded.py,run_guarded.py,check_exports.py: reproducible audit workflow.

Limits
------
The original generic NaN-congruence rejection deficiency is not silently
repaired in the frozen source. Guarded replay addresses finite/order behavior
of these fixed executions. The platform/compiler/basic correctly rounded
operations and the analytic core/exterior inequalities remain assumptions
or separate audit obligations. A matching byte comparison does not establish
all possible inputs, all architectures, an independent backend, or formal
proof-assistant verification. The final global barrier still requires the
source-level mathematical reductions and the complete supervisor chain.
