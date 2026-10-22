RECORD OF THE SEPARATOR AUDIT — 21 OCTOBER 2026

Read INFORME_AUDITORIA_ES.txt for the current Spanish assessment.
This directory records a new execution of the standalone separator and
several independent or instrumented component checks. It does not record
a new complete connection or optimality run.

The large frozen input arrays are not duplicated. Extract the original
separator archive into this directory, yielding xy_global_separator/.
Its original ZIP identity is stated in the report. Do not modify it.

Run the narrow environment preflight from the manuscript directory:
  python3 checks/verify_supported_environment.py

For a new untouched separator execution, from this directory:
  python3 xy_global_separator/verify_separator.py \
    xy_global_separator/distribution --output fresh_results_NEW --workers 2

The included fresh_results/ is this audit's executed output. The complete
byte/hash comparison with the original can be repeated by:
  python3 compare_fresh.py
To compare your additional execution, keep the included record and use:
  python3 compare_fresh.py --fresh fresh_results_NEW \
    --output fresh_comparison_NEW.json

For the independent exact geometry audit, see independent_geometry/README_ES.txt.
Its script accepts explicit package and geometry paths, so no source copying
is necessary when the original package is already elsewhere.

For source binding and the conditional rational boundary replay:
  python3 independent_supervisor/audit_supervisor.py
That particular replay intentionally uses archived local enclosures. The
complete fresh supervisor above generates those enclosures anew.

For full finite/order instrumentation and exact exported-endpoint checks,
follow independent_backend/README_REPRODUCTION.txt. It recreates large raw
case inputs from the frozen package and compiles diagnostic copies. The
compact results and exact hexadecimal sidecars are included; compiled
binaries, raw reconstructions and Python caches are excluded.

The original generic NaN acceptance and non-outward bare-decimal examples
remain documented. Their relevance is assessed for each actual consumer
and fixed execution; neither is silently rebranded as a refutation of the
fixed separator theorem or a proof of arbitrary-input correctness.

bibliography_spot/ records primary-source metadata checks and the newly
identified preprint. The manuscript's actual final wording and bibliography
are in the parent source bundle, not the proposed-insertion note.

FILES.sha256 inventories this audit directory. Paths and timing fields in
historical outputs describe the actual run and are not installation paths.
