# Certified global annihilation barriers in an infinite XY lattice

Author: Rodolfo Leonardo Sumoza Matos. Prepared for SIAM Journal on Applied
Dynamical Systems (SIADS), documentation revision 7, 22 October 2026.

Repository: [https://github.com/rsumoza/global_annihilation_barriers_infinite_-XY_lattice](https://github.com/rsumoza/global_annihilation_barriers_infinite_-XY_lattice)

Planned frozen Release: [siads-v7-2026-10-22](https://github.com/rsumoza/global_annihilation_barriers_infinite_-XY_lattice/releases/tag/siads-v7-2026-10-22). This is a configured destination;
publication and remote download have not been verified. See `repository_state.json`.
No DOI or journal submission is claimed.

## Result and scope

The fixed global certificate concerns a neutral dipole at separation 12
on the infinite exact-cosine square XY lattice. The certified barrier is
`0.00567867 < B_ann < 0.00570116`. An explicit infinite class of initial
perturbations with controlled block tails reuses the same communication height.
The manuscript proves arbitrary-precision enclosure termination for that class;
the general refinement interface is not implemented. Completeness for arbitrary
XY configurations remains open. Revision 7 changes distribution and links,
not mathematical statements or the frozen scientific inputs.

## Start here

`index.txt` is the detailed English supplement and reproduction index.
`REPOSITORY_RELEASE_ES.txt` tells the author exactly what to upload.
`README.txt` preserves the detailed scientific and historical execution record.
`CURRENT_AUDIT_ES.txt` explains limitations.

## Small arithmetic checks

```sh
python3 checks/verify_robustness_constants.py
python3 checks/verify_green_tail_budget.py
python3 checks/verify_infinite_class.py --self-test --output /tmp/xy-infinite-class.json
```

These commands use Python's standard library. They check rational consequences,
representative finite prefixes and negative controls, rather than rerunning the
global certificate or computing arbitrary nontrivial fibers.

## Full reproduction

The Git tree and GitHub's automatic source-code ZIP do not include S2-S4.
After the Release has been published, download its three frozen inputs:

```sh
python3 checks/download_release_assets.py --assets-dir ../release_assets
python3 checks/verify_release_assets.py --assets-dir ../release_assets
python3 checks/verify_supported_environment.py --output /tmp/xy-supported-environment.json
```

`checks/download_release_assets.py --list` prints the intended URLs without
accessing the network. The downloader checks byte counts and SHA-256 before
installing any file. See `index.txt` for extraction and the exact full-replay
command. Keep new output separate from archived `proof_record/`.

The reviewed platform is x86_64 Linux with a 64-bit long-double significand,
little-endian binary64 and round-to-nearest. The frozen backend is reused;
this is not a second interval implementation or proof-assistant verification.

## Build

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
cd supplement
pdflatex -interaction=nonstopmode -halt-on-error Supplementary_Reproducibility.tex
pdflatex -interaction=nonstopmode -halt-on-error Supplementary_Reproducibility.tex
```

Optional figure regeneration requires NumPy and Matplotlib. The arithmetic
checks do not require them. See `README.txt` for TeX packages and optional SIAM macros.

## Release and integrity

Upload the contents of this directory to Git, retaining `checks/`, `figures/`,
`sections/`, `supplement/`, `proof_record/`, and `verification_notes/`.
Attach the three original ZIPs to the planned Release, together with the
separately supplied `SIADS_XY_Repository_v7_2026-10-22.zip` and `SOURCE_SNAPSHOT.sha256`. Preserve their exact
bytes and canonical filenames. `MANIFEST.sha256` covers this source tree;
`release-assets.sha256` covers the three original ZIPs. Archived audit manifests
remain historical records and are not rewritten.

The initial journal submission also needs the supplement and index, plus its
required frozen material through an accepted submission route. GitHub does not
automatically waive that requirement. Original licenses and author rights are
unchanged; this package does not grant a new license to third-party material.
