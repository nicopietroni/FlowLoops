# FlowLoops: Robustly Extracting Descriptive Curve Networks from Standard Triangle Meshes
Edward Su, Nico Pietroni, Nicholas Vining, Alla Sheffer <br/>
*SIGGRAPH ASIA 2026* <br/>

![alt text](teaser.jpg)

## Abstract
Sparse, descriptive networks of 3D curves can efficiently represent complex free-form shapes, and serve as both a useful design tool and a powerful compact shape representation. Despite their utility, extracting these descriptive curve networks from shapes described using standard representations such as triangular meshes remains a challenging and open problem. Previous approach for generating such networks required well-behaved, spiral-free, curvature-aligned, quad-dominant meshes as inputs, but these meshes are themselves challenging to generate reliably. We present FlowLoops, an algorithm for extracting descriptive, compact 3D curve networks from general triangular meshes.

We observe that descriptive curves are expected to align with the curvature field in anisotropic surface areas and smoothly extend into isotropic areas following approximately geodesic directions, and that the descriptiveness of a given curve network can be explicitly evaluated by using it to reconstruct a surface and measuring the difference between the reconstructed and input surfaces. We consequently cast the problem of extracting a descriptive curve network as one of tracing a minimal set of field-aligned curves across a smooth, curvature-aligned, branched covering of the surface. When placing the curves we achieve sparsity by maximizing the spacing between them while ensuring accurate reconstruction. We validate our approach by extracting descriptive curve networks from a wide range of challenging inputs, and demonstrate that our method produces curve networks that accurately reproduce the original input surface. We further compare our results against potential alternative demonstrating their generality and superior descriptive power.

**BibTex**
```
To be added
```

### Summary
This repository is the official repository of the paper *FlowLoops: Robustly Extracting Descriptive Curve Networks from Standard Triangle Meshes*. <br>
Currently this page is under development and changes addressing missing areas will be added sortly.

### Data
Complete output for four models is tracked in `data/sample/`, and two input meshes
for a first run are in `test_data/`. The complete dataset — results for all 125
models plus the input meshes that produced them — is published as release
assets; see Download.

See [`data/README.md`](data/README.md) for the per-model file layout.

### Download

| Asset | Size | Contents |
|---|---|---|
| [`flowloops-inputs.tar.gz`](https://github.com/nicopietroni/LoopWeaver/releases/latest/download/flowloops-inputs.tar.gz) | 50 MB | the 125 input triangle meshes |
| [`flowloops-results.tar.gz`](https://github.com/nicopietroni/LoopWeaver/releases/latest/download/flowloops-results.tar.gz) | 598 MB | curve networks, surfaced meshes, cross fields, timings, screenshots |

```bash
curl -LO https://github.com/nicopietroni/LoopWeaver/releases/latest/download/flowloops-inputs.tar.gz
curl -LO https://github.com/nicopietroni/LoopWeaver/releases/latest/download/flowloops-results.tar.gz
```

> **Cloning:** earlier revisions of this repository tracked the full 1.5 GB
> dataset, so those blobs remain in the git history and a plain `git clone`
> still downloads them. To skip them:
>
> ```bash
> git clone --depth 1 --filter=blob:none https://github.com/nicopietroni/LoopWeaver.git
> ```

### Build

**Requirements**

- CMake ≥ 3.16
- A C++20 compiler
- OpenGL, GLEW, GLFW 3
- [TracingFieldLib](https://github.com/nicopietroni/TracingFieldLib) (header
  only), downloaded separately from its library page — see below

All the other third-party dependencies (Eigen, libigl, Dear ImGui, Triangle,
Curve-Surfacing-Core) are vendored under `external/`. See
[`external/README.md`](external/README.md) for their upstream sources, pinned
versions and licenses.

Verified on macOS 26 (Apple Silicon) with Apple clang 21 and CMake 4.3:

```bash
brew install cmake glew glfw
```

The source carries a non-Apple OpenGL path, so Linux should work, but it has
not been tested. The dependencies there would be:

```bash
sudo apt install cmake build-essential libglew-dev libglfw3-dev libglu1-mesa-dev
```

**TracingFieldLib**

TracingFieldLib is not included in this repository. Download it from
[github.com/nicopietroni/TracingFieldLib](https://github.com/nicopietroni/TracingFieldLib)
and place it next to LoopWeaver:

```
<parent folder>/
├── TracingFieldLib/       <- downloaded from the library page
│   └── include/
└── LoopWeaver/            <- this repository
    ├── external/
    └── loopweaver/
```

CMake looks for `../TracingFieldLib` by default; if it is somewhere else, pass
`-DTRACINGFIELDLIB_DIR=<path>` when configuring. Only its `include/` folder is
used (its bundled `external_lib/` is not needed). The lookup table
`include/Lut.txt` is read at run time from its absolute path, so do not move
TracingFieldLib after building.

**Build**

```bash
git clone https://github.com/nicopietroni/TracingFieldLib.git
git clone https://github.com/nicopietroni/LoopWeaver.git
cd LoopWeaver
cmake -S . -B build
cmake --build build -j8
```

This defaults to a Release build. Pass `-DCMAKE_BUILD_TYPE=Debug` for
debugging — expect an order-of-magnitude slowdown in the solver loops.
Other options: `-DENABLE_SANITIZER=ON` for ASan/UBSan, and
`-DENABLE_NATIVE_ARCH=ON` to tune Release for the local CPU (faster, but the
resulting binary is not portable to other machines).

Do **not** substitute a system Eigen. The vendored Eigen 3.4.0 is the version
this code and the vendored libigl were developed against; Eigen 5.x will not
compile.

### Usage

```bash
./build/loop_weaver test_data/ellipse.obj            # GUI, mesh loaded
./build/loop_weaver test_data/ellipse.obj -batch     # whole pipeline, save, quit
./build/loop_weaver                                  # print all options
```

`loop_weaver` is an interactive OpenGL application. With `-batch` it runs the
full pipeline (symmetry cut, cleaning, cross field, decomposition, relaxation,
final reconstruction), writes its results and exits; the window is kept hidden
in that mode, but it still needs a display, so it will not run over a plain
SSH session or in headless CI.

**Parameters** are set in three layers, each overriding the previous one:

1. the compiled defaults;
2. `loopweaver_settings.txt` next to the executable (`build/` for a normal
   build), one parameter per line. It is written with the compiled defaults
   the first time the program runs; delete it to write it again after an
   update that changed them. `-settings <file>` reads another file instead;
3. the command-line options.

Every run also saves its parameters as `<name>_settings.txt`, so it can be
repeated with `-settings <name>_settings.txt`.

**Main options** — the full list, with defaults, is in
[`README_commandline.md`](README_commandline.md):

| Option | Meaning |
|---|---|
| `<mesh>.obj` | Input triangle mesh. Required, and must come before any field file. |
| `<field>.field` | Optional pre-computed cross field; without it the field is computed from the curvature. |
| `-batch` | Run the whole pipeline without interaction, save everything, then exit. |
| `-settings <file>` | Settings file to use instead of `loopweaver_settings.txt`. |
| `-symm` / `-no_symm` | Cut the mesh with the symmetry plane x=0, decompose one half and mirror the result. |
| `-angle <deg>` | Sharp feature angle. `180` means no sharp features. |
| `-error_setup <0-3>` | Reconstruction error preset (see below). A negative value disables the reconstruction condition. |
| `-error <ratio>`, `-normal_error <deg>` | Set the reconstruction distance and normal-angle limits directly. |
| `-side_angle <deg>`, `-bezier <pct>`, `-gauss <deg>` | Patch conditions: side turning angle, Bézier fit error, total Gaussian curvature. Each has a `-no_` form to disable it. |
| `-subsampl <n>` | Subsampling of the paths in the reconstruction condition. |
| `-thinplate` | Thin-plate-only surfaces, in the condition and the final reconstruction. |
| `-dynamic` / `-no_dynamic`, `-relax` / `-no_relax` | Retrace candidates after each insertion; relax and prune the paths at the end. Both on by default. |

`-error_setup` presets, from strictest to loosest:

| Value | Max error ratio | Max normal angle | Angle percentile |
|---|---|---|---|
| `0` | 0.02 | 30° | 0.2 |
| `1` | 0.03 | 40° | 0.3 |
| `2` | 0.04 | 50° | 0.5 |
| `3` | 0.05 | 60° | 0.6 |

Without `-error_setup` the defaults are 0.02, 30° and 0.5, which match none of
the presets exactly.

The options `-curvchange`, `-create_mesh`, `-minsides` and `-curvangle` of
earlier versions are no longer accepted; unknown options are reported and
ignored.

**Outputs** are written to the **current working directory** (not next to the
input mesh), named after the input mesh `<name>.obj`:

| File | Content |
|---|---|
| `<name>.path` | The loops, as vertex indices of `<name>_remeshed.obj`. With symmetry they cover the whole mirrored mesh. |
| `<name>_remeshed.obj` | The processed mesh the paths refer to. |
| `<name>_edge.obj` | The loops as an edge mesh. |
| `<name>_solved.obj` | The final reconstructed surface. |
| `<name>.field` | The cross field used. |
| `<name>_timings.txt` | Field, tracing and surfacing times, and the total loop length. |
| `<name>_settings.txt` | All the parameters of the run. |

The surfacing solver's exchange files (`.curve`, `.cycle`, `.normalcurve`,
`.feat`, `.curve_types`) are deleted at the end of the run. Earlier versions
kept them and also saved `_config.txt` and `_original.png` /
`_boundaries.png` screenshots, which is why `data/` contains those files.

**Batch runs over a folder**

`run_scripts/` contains wrappers that run every `.obj` in a folder with a fixed
option set, move each mesh's outputs into `build/done/`, and log every run to
`build/timeouts.txt`:

```bash
./run_scripts/run_weaver_features_symm.sh /path/to/meshes [extra options]
```

| Script | Options |
|---|---|
| `run_weaver_features_symm.sh` | `-symm -angle 45 -error_setup 1 -subsampl 2` |
| `run_weaver_features_no_symm.sh` | `-no_symm -angle 45 -error_setup 1 -subsampl 2` |
| `run_weaver_no_features_symm.sh` | `-symm -angle 180 -error_setup 1 -subsampl 2` |
| `run_weaver_no_features_no_symm.sh` | `-no_symm -angle 180 -error_setup 1 -subsampl 2` |

Options after the folder are added to every run (for example
`-no_dynamic -thinplate`). `run_all_weaver.sh [dataset folder] [extra options]`
runs all four variants over a dataset folder containing `no_symm_features/`,
`symm_features/`, `no_symm_no_features/` and `symm_no_features/`, writing to
`build/output/<variant>/es<preset>/`.

The scripts expect the binary at `build/loop_weaver` and require GNU coreutils
for `gtimeout` (`brew install coreutils` on macOS). Each mesh is given a 1800 s
budget.

In case you have technical issues with the project, please raise an issue on GitHub
