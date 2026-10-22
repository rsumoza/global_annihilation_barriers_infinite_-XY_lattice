# Independent audit of connection supervision and scalar consumers

Date: 2026-10-21. Scope: the frozen connection package, its scalar
postprocessors, and its consumption by the optimality overlay. The C++ backend
and full fresh execution are audited separately. No original file was changed.

## Result and reproducible evidence

`python3 audit_supervisor.py` verifies the external manifest/supervisor pins,
checks all 44 scientific members and the nested stage-2 manifest identity,
checks the shared minimum/saddle/mode raw identities, and replays both exact
scalar postprocessors using archived local enclosures. Both outputs match
their archived versions as parsed exact data.

Seven deliberately invalid records are rejected with nonzero exit status and
no pass output: wrong Morse index, missing path-time coverage, nonzero final
integer lift, charged endpoint, leg energy above the saddle, wrong join
direction, and failed Hermite energy gap. The script and `audit_summary.json`
record the exact values, raw endpoint hashes and rejection messages.

No fixed-data counterexample was found to the Python/source-binding portion
of this connection chain. This conclusion is narrower than a fresh interval
proof and does not establish arbitrary-input robustness of the backend.

## Source identity and endpoint binding

`verify_global.py` checks the externally pinned distribution manifest, copies
and hashes every member into a private snapshot, verifies the fixed reference,
fine saddle and negative mode, then invokes the pinned stage-2 supervisor.
That supervisor separately checks its own nested scientific manifest and
copies the very same base field to the minimum, saddle-reference and barrier
inputs. Thus the energies are in one normalization and refer to the same
candidate minimum, not merely similarly named files.

The return leg starts at the exact zero perturbation; the supervisor tests its
first node byte-for-byte against zero. It extracts the *last* return-leg node
as the plus join, and the *first* annihilation-leg node as the minus join. The
join-point C++ program receives these extracted bytes. They are not a pair of
independently supplied endpoint proposals. The raw leg configurations and
payload lengths are checked, and have respectively 8 and 19 nodes on side 256.
The separately certified literal full-annihilation path is supplementary; it
does not replace the Hermite attachments used for exact saddle-height equality.

The connection supervisor regenerates all consumed results. It does not use
`results_provided` as premises. A pre-existing output directory is permitted,
but all numerical prerequisites and consumed path/join records are freshly
overwritten; no stale-result premise was identified on this source path.

## Scalar inequalities and serialization

The baseline has bare 22-significant-digit decimal output; saddle, path and
join summaries have bare 24-significant-digit output. These strings are not
automatically outward enclosures. Their consumers account for this explicitly:

1. The stage-2 supervisor applies `1e-18` in the conservative direction to its
   minimum and marked-sublevel dependency thresholds.
2. `postprocess_certificate.py` reads exact fractions and pads saddle residual,
   signed Hessian bounds, correction radius, branch margin and relative
   energy endpoints. It explicitly checks the declared subunit domain for
   the saddle quantities. The correction energy error is computed from the
   padded values as `min(rho^2/2, delta*rho/2+rho^3/12)`.
3. The minimum-reference energy error is recomputed from the bound minimum
   certificate as `delta_min^2/(2*mu_min)`, with residual padded upward and
   segment coercivity padded downward. It is checked against the stated
   `0.00001168099` allowance.
4. `postprocess_global.py` pads saddle relative-energy endpoints, each join
   distance and signed projection, join energy uppers, path-height uppers,
   and the marked first-exit bound. Its inequalities then use exact fractions.

For the fixed outputs these quantities have magnitude below one, except the
baseline reference energy, whose magnitude is below 100. `1e-18` dominates
the decimal formatting errors for those ranges. Ordinary values in the large
path tables are diagnostic: no proof consumer treats an unpadded table decimal
as its final global upper bound.

The optimality overlay makes some additional unpadded comparisons of printed
saddle/minimum scalars with safe constants. Independently inserting the same
pad preserves all of them on these fixed inputs. In particular:

- `0.012825 - (rho_saddle + pad)` is approximately `6.163439498475e-8`.
- `(m_candidate - pad) - 0.0719549999` is approximately `9.99999697452829e-11`.

Its outer-level comparison has a much smaller positive margin, approximately
`1.2192816430179687e-21`. This comparison is nevertheless between the exact
rational stored saddle upper bound, which already contains serialization pads,
and an exact decimal target. It does not compare rounded display decimals.
The script records the full rational difference. No additional unaccounted
display-rounding error is introduced by that comparison.

## Global geometric connection

For each affine path segment, the C++ path checker sums the complete compact
perturbation support and bounds the interior energy by the maximum endpoint
upper bound plus `||D(z_next-z_previous)||_2^2/8`. This follows from the
universal bound `H >= -I` (equivalently here the established `||H|| <= 1`
suffices) in the homogeneous energy norm; it covers the
whole parameter interval, not only sampled nodes.

The join checker measures distances in the full compact difference between
the side-4096 saddle candidate, side-2048 reference and side-256 perturbation,
including the zero extensions and their boundary edges. The postprocessor
adds the full saddle correction radius before estimating the exact-saddle
distance and projected sign.

For a join with corrected distance bound `r`, signed projection lower bound
`p0`, and local passage dimensions `S,R`, the exact checks require:

- `p0 > S`;
- `R*p0 - S*r > 0`, so the ray reaches its longitudinal face before a lateral
  boundary;
- `E_saddle_lower - E_join_upper - r^3/6 > 0`, so the Hermite bridge from the
  exact saddle descends strictly away from it.

The two endpoint signs are required to be opposite. The fixed Hermite gaps
are approximately `1.12927448117e-4` and `8.06002874120e-5`; face margins are
approximately `6.62195266345e-5` and `6.58979413250e-5`. Both complete global
legs lie below the certified saddle lower level, with margins approximately
`4.23213565178e-5` and `4.82741535919e-5`. These large strict margins provide
the intended attachments to the two local descending components.

For the terminal leg, the checker proves no remaining plaquette charge,
consistent compact integer lifting and strict principal-branch bounds inside
and outside the support. The postprocessor additionally requires that the
compact lift be identically zero. The resulting global real lift has all
bond differences strictly between `-pi` and `pi`, so scalar contraction to
zero decreases every cosine bond energy. This supplies the continuation to
a uniform state in the homogeneous energy space. It is a continuous
variational connection, not a heteroclinic orbit assertion.

The short segment from the exact minimum to its common reference is controlled
by the validated convex ball. The finite-path/Hermite pieces therefore attach
to the exact minimum and exact saddle rather than merely numerical candidates.

## Remaining limits

- The generic NaN acceptance issue in an internal congruence routine is not
  resolved by these scalar postprocessors. Backend instrumentation is needed
  to determine whether fixed intermediate inputs remain finite.
- The original platform guards permit all mantissas at least 64 bits. The
  justified submission environment should be restricted to the tested
  extended64 format, as in the revised manuscript preflight. Broad platform
  portability has not been established.
- Exact source identity and successful fresh execution are separate facts.
  This report's conditional scalar replays cannot replace the parent's full
  fresh minimum/saddle/path/oracle execution.
- The standalone connection establishes the specified upper path and local
  descending attachments. The assertion that this height is the *global*
  annihilation barrier additionally requires the complete separator and
  exact-saddle-cap overlay; the original connection certificate correctly
  labels global optimality and heteroclinic existence as unproved there.
- No universal XY configuration theorem or proof-assistant formalization
  follows from this case-specific audit.
