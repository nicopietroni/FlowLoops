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

All third-party dependencies are vendored under `external/` — there are no git
submodules, and the build makes no network requests. See
[`external/README.md`](external/README.md) for upstream sources, pinned
versions and licenses.

**Requirements**

- CMake ≥ 3.16
- A C++20 compiler
- OpenGL, GLEW, GLFW 3

Verified on macOS 26 (Apple Silicon) with Apple clang 21 and CMake 4.3:

```bash
brew install cmake glew glfw
```

The source carries a non-Apple OpenGL path, so Linux should work, but it has
not been tested. The dependencies there would be:

```bash
sudo apt install cmake build-essential libglew-dev libglfw3-dev libglu1-mesa-dev
```

**Build**

```bash
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
./build/loop_weaver test_data/ellipse.obj
```

`loop_weaver` is an interactive OpenGL application and **always opens a window,
including in `-batch` mode** — it needs a display, so it will not run over a
plain SSH session or in headless CI. In `-batch` mode it runs the full pipeline,
writes its results and screenshots, and then exits on its own.

Outputs are written to the **current working directory**, named after the input
mesh: `.curve`, `.cycle`, `.feat`, `.field`, `.normalcurve`, `.path`,
`_edge.obj`, `_remeshed.obj`, `_solved.obj`, `_config.txt`, `_timings.txt`, and
`_original.png` / `_boundaries.png` screenshots.

**Arguments**

| Argument | Meaning |
|---|---|
| `<mesh>.obj` | Input triangle mesh. Required, and must come before any field file. |
| `<field>.field`, `<field>.ffield` | Optional pre-computed cross field. |
| `-batch` | Run the whole pipeline without interaction, then exit. |
| `-symm` | Enable symmetry detection and symmetric tracing. |
| `-split` | Enable split removal in patch decomposition. |
| `-curvchange` | Use curvature-change nodes. |
| `-create_mesh` | Remesh instead of tracing on the original mesh connectivity. |
| `-angle <deg>` | Feature (sharp edge) detection angle. Default `50`. |
| `-minsides <n>` | Minimum patch sides, clamped to 2–3. Default `2`. |
| `-curvangle <deg>` | Max summed side angle, clamped to ≤ 180. Default `180`; `≤ 0` disables the side-curvature condition. |
| `-bezier <pct>` | Max Bezier fitting error %, clamped to ≤ 5. Default `2`; `≤ 0` disables the condition. |
| `-subsampl <f>` | Remesh / subsample factor. |
| `-error_setup <0-3>` | Reconstruction error preset (see below). A negative value disables the reconstruction condition. |

`-error_setup` presets, from strictest to loosest:

| Value | Max error ratio | Max normal angle | Angle percentile |
|---|---|---|---|
| `0` | 0.02 | 30° | 0.2 |
| `1` | 0.03 | 40° | 0.3 |
| `2` | 0.04 | 50° | 0.5 |
| `3` | 0.05 | 60° | 0.6 |

Preset `0` is also the built-in default when `-error_setup` is not passed.

**Batch runs over a folder**

`run_scripts/` contains wrappers that run every `.obj` in a folder with a fixed
flag combination, move each result into `build/done/`, and log timeouts:

```bash
./run_scripts/run_weaver_features_symm.sh /path/to/meshes
```

These expect the binary at `build/loop_weaver` and require GNU coreutils for
`gtimeout` (`brew install coreutils` on macOS). Each mesh is given a 1200 s
budget.

In case you have technical issues with the project, please raise an issue on GitHub
