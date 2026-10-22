# Independent audit of the separator supervisor and Python consumers

Date: 2026-10-21. Frozen package: `xy_global_separator`.

## Scope and conclusion

This audit checked source-level acceptance contracts, the exact scientific
ledger, the C++/Python decimal interface, and the rational boundary checker.
It did not recompute the C++ minimum or local charts; that fresh run is a
separate task. The included replay uses the archived local enclosures and is
explicitly conditional on their validity. No fixed-input counterexample to
the separator was found in this part of the review.

The baseline's bare decimal output is not, by itself, an outward interval.
For this separator, however, every baseline scalar subsequently used in the
proof is padded by the existing exact rational `PAD=1e-18`. The fixed outputs
are small enough that this pad dominates the 22-significant-digit printing
error. Consequently the generic printing counterexample does not establish
a failure of this separator's acceptance chain.

## Reproduced checks

Run `python3 audit_supervisor.py` from any directory. It leaves the frozen
package unchanged. Results are recorded in `audit_summary.json`.

- All 564 scientific manifest members match their recorded SHA-256 hashes.
- Before those checks, the scientific manifest and supervisor themselves are
  checked against their frozen external pins, respectively
  `991c744eab7e812d629662a0c0bc5cbc177e9dd2e9daa67caeefab7a814f941b`
  and `d4646d988e538075b51d0c0bdd2d970976d1ccb53156ad52d414ec87f5c0095c`.
- All 109 case IDs, grid indices, center filenames, and geometry records agree.
- The 109 chart configurations use the expected dimension 128; their radius
  strings are canonical exact binary64 hexadecimal values.
- The Green archive has no duplicate edge identity and its size agrees with
  the column ledger.
- A new rational boundary replay agrees with the archived boundary output,
  after removing paths, run time, and path-dependent input-hash fields. It
  checks 109 charts and performs 2150 coordinate-face cell evaluations, with
  zero unresolved radial, circular, or coordinate-face cells.
- Four deliberately invalid baseline records are rejected with exit code 2:
  NaN residual, residual at the admissibility threshold, radius at the
  admissibility threshold, and nonpositive coercivity. No pass output is
  produced by these negative controls.

## Acceptance and source binding

`verify_separator.py` pins the scientific manifest by its SHA-256, copies each
scientific file to a private temporary directory, verifies the copied bytes,
compiles all three C++ verifiers, reconstructs every Green proposal and checks
its per-case hash, and generates each consumed numerical result afresh. It
does not read `results_provided` as a premise. A successful child return code
alone is insufficient: JSON statuses, positivity and interior conditions,
the boundary checker and final separation status are also required.

The output directory need not be empty. This is weaker operational hygiene
than a mandatory fresh directory, but no stale numerical premise was found:
all files consumed along this source path are freshly overwritten, constructed
from the current ledger, or copied only after the corresponding child returns
successfully. Extra unrelated old files may appear in the final output hash
inventory. A clean output directory is preferable for an unambiguous record.

The Python geometric checker uses exact `Fraction` conversions and dyadic
outward operations on the 2^160 grid. Domain radii are read from their exact
hexadecimal binary64 representations, not informational decimal radii. The
closed coordinate-face intervals include both endpoints and corners; only
cells proved strictly outside the disk are discarded. A low face cell is
accepted only if it is assigned to the strict interior of another chart,
using its transverse distance bound, or has a certified energy floor at least
the target level. Failure to prove coverage yields `UNRESOLVED` and a nonzero
return code, not an escape or instability claim.

Geometry/chart interval intersections check compatibility, but intersection
alone would not prove that two independent quantities are identical. Here the
needed identity is supplied by the pinned common `q.raw`, the case center
filename binding, and the source using those same bytes. The independent
ledger check confirms all these identities on the present frozen payload.

## Baseline decimal interface

The baseline calls `engine.cpp::run`, whose output has 22 significant digits
without a print-specific `nextafter` expansion. The proof consumer
`postprocess_separator.py` reads the decimals exactly with `parse_float=Fraction`
and uses only these four baseline numerical quantities:

| Quantity | Direction applied by the consumer |
| --- | --- |
| `delta` | `delta + PAD` |
| `rho` | `rho + PAD` |
| `c` | `c - PAD` |
| `energy_candidate[0]` | `energy_candidate[0] - PAD` |

All four printed magnitudes are below 100. A conservative absolute formatting
error bound of 5e-20 for nearest 22-significant-digit decimal output is below
`PAD=1e-18`. This supplies the required outward direction at this interface.
The remaining residual slack after padding is exactly

`542000672032347 / 1000000000000000000000000`,

approximately 5.42000672032347e-10. It also dominates the extended64 rounding
error of the local chart's compiled decimal residual threshold
`.001166548L`. The radius and convex-bridge inequalities have still larger
slack. The geometry and polynomial endpoints used by the Python checker are
instead printed with 30 significant digits after a further outward
`nextafter` step; their validity assumes the supported extended64 environment.

The coarse chart JSON values use ordinary decimals but are not transferred
as quantitative boundary-proof bounds: the supervisor uses their positivity
and an internally calculated `interior_pass` Boolean, while the rational
boundary calculation consumes the separately exported outward coefficients.

## Limits and recommended precision

- This audit does not turn the conditional checker into a universal validator
  of fabricated JSON. The full supervisor supplies its source and raw-data
  prerequisites. Standalone postprocessing intentionally trusts its declared
  prerequisites.
- The source guards permit `LDBL_MANT_DIG >= 64`, but the 30-digit padded
  decimal proof is justified here for `LDBL_MANT_DIG == 64`. The manuscript's
  restricted-environment preflight is necessary; a broad statement covering
  arbitrary higher-precision long double implementations would be unsupported.
- Nonfinite arithmetic inside the C++ congruence routine remains a separate
  defensive issue. Rejection of a NaN final baseline residual does not prove
  that every intermediate nonfinite quantity is rejected. No claim of that
  stronger contract is made by this audit.
- These checks establish neither the connection upper path nor the final
  optimal-saddle theorem. They concern the standalone global separator only.
- An independent conditional replay and source audit are different from an
  independently implemented interval backend or a proof-assistant
  formalization. Neither stronger description is warranted here.
