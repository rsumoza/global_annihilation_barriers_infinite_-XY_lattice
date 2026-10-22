Connection backend audit, 21 October 2026

This diagnostic retains the frozen interval arithmetic formulas and proposal
bytes. It inserts finite/order guards and writes separate exact hexadecimal
long-double sidecars. Scientific stdout and path tables are compared byte for
byte with a fresh original execution. Decimal summaries are parsed as Fractions;
their actual PAD corrections and exact compiled-literal dependencies are checked.
This is not a second interval engine or a formal verification of the compiler.

Expected directory layout
-------------------------
workspace/
  audit_separator_2026_10_21/independent_backend/
    ... completed guarded baseline, exact stderr log and source bindings ...
  audit_connection_2026_10_21/
    xy_global_connection/       untouched extracted original
    full_fresh/fresh_connection/  freshly recomputed original results
    independent_backend/       these source files and completed compact results

The original complete connection supervisor is:

 python3 xy_global_connection/verify_global.py \
    xy_global_connection/distribution --output full_fresh/fresh_connection

The coordinated three-package original replay may create the same fresh output
via verify_optimality.py. Supplied historical JSON is not accepted as the clean
reference by check_exact_exports.py; the fresh supervisor must finish first.

Instrumented reproduction
-------------------------
From the independent_backend/ directory:

 python3 prepare_instrumented.py
 python3 prepare_baseline_reuse.py
 python3 run_small.py
 python3 run_saddle_monitored.py
 python3 select_saddle_guard.py --variant original
 python3 check_exact_exports.py

prepare_instrumented.py builds guarded baseline/saddle/barrier/path/join programs
with exactly the recorded flags. It hardlinks immutable raw proposals into
separate working directories; none of the frozen source or raw inputs is written.
Join arrays are extracted from the actual first/last fixed leg nodes.

The minimum has exactly the same original engine and all six fixed input files
as the completed separator baseline. baseline_reuse_binding.json records every
identity. The replay reuses that guarded computation; a bounded scalar probe
reconstructs mu=dn(c-rho) and the branch literal from its exact c/rho sidecar.
The probe's actual c/rho arguments are also exported and checked. A standalone
alternative is to run the built stage2/baseline/verify_guarded executable from
that directory and retain its fresh exact sidecar; this audit's delivered
checker explicitly documents and checks the baseline reuse route instead.

run_small.py runs the marked sublevel at mesh 640, paths with 8/19/34 nodes, and
the exact stored join comparison against the 4096 candidate. No approximate
optimizer or previously supplied result is trusted. run_saddle_monitored.py
runs the complete 4096 saddle with finite/order guards. Its allocation monitor
uses aggregate anon+kernel+shmem, excluding reclaimable disk page cache. It
requires that budget below 4 GiB before starting, samples every 2 seconds, and
terminates only its own guarded child if the budget exceeds 6.8 GiB or an OOM
event increments. Serial execution after the original is the simplest route.

Exact checks and scope
----------------------
results/exact_export_report.json records:
- PAD direction for every actually consumed summary scalar and interval;
- exact baseline dependencies against barrier's compiled rho/delta/mustar;
- equality of the baseline and barrier 1.12L branch literal (same source/compiler);
- the optimality overlay's additional unpadded exact rational thresholds;
- optional PAD enclosure of all path node/segment diagnostic tables;
- original scientific manifest, working proposal identities, sidecar/export
  hashes, audit source hashes and compiler flags.

The bare node/segment and barrier polynomial dump decimals are approximate
records, not universally outward endpoints. They are not numerical premises of
the original connection proof; only padded summaries and internal checks are.

Package only scripts, C++ sources/patches, source/input binding records and
results/. Working raw inputs under stage2/ and connections/, compiled binaries
and __pycache__ are regenerated and should not be included in the compact audit.
The original frozen package remains the authoritative source/proposal bundle.

Lazy diagnostic messages (same checks)
-------------------------------------
The initial guarded saddle constructed std::string messages on successful hot
guard calls. A bounded 2-million-addition benchmark measured 0.4335 s versus
0.1495 s when diagnostic strings are constructed only after a predicate fails.
The exact final interval and every predicate/counter were unchanged. The
lazy-message variant is prepared and run with:

 python3 prepare_fast_guard.py
 python3 run_saddle_fast_monitored.py

It has a separate source directory, executable, stdout, exact sidecar and memory
monitor record. It never overwrites the original frozen package. Only a
completed successful run with scientific stdout byte-identical to the fresh
original may be selected. If both variants run, the incomplete extra is stopped
after the completed variant is accepted and is excluded from proof evidence.
selected_saddle_guard_binding.json identifies the actual chosen executable,
source hashes, stdout and sidecar (select_saddle_guard.py --variant lazy is
used after a completed lazy run); root-level canonical saddle files are copied
from that completed variant. Benchmarks and partial output are diagnostic only.

Integer-string budget
---------------------
Exact extended64 subnormal endpoints can yield rational denominators with about
5000 decimal digits. check_exact_exports.py sets a 30000-digit serialization
budget explicitly; this changes report string conversion only, not arithmetic
or the numerical bounds. The initial 4300-digit default could not stringify a
first-node subnormal diagnostic endpoint; all bound checks passed after the
reporting budget was adjusted.
