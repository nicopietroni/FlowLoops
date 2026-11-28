# Vendored dependencies

Everything LoopWeaver needs to build, other than CMake, a C++20 compiler, the
three system packages (OpenGL, GLEW, GLFW) and
[TracingFieldLib](https://github.com/nicopietroni/TracingFieldLib), is copied
into this directory. TracingFieldLib is downloaded separately from its library
page; see the Build section of the main [README](../README.md). There are no
git submodules and the build makes no network calls.

Each tree is a **pruned** copy — only the parts the build actually compiles or
includes were taken. Upstream licenses are kept alongside the code they cover.

| Directory | Upstream | Pinned at | License |
|---|---|---|---|
| `curve-surfacing-core/` | [Edisuism/CurveSurfacingLib-V2](https://github.com/Edisuism/CurveSurfacingLib-V2), `Curve-Surfacing-Core/` | `7604729` | Research / non-commercial — see upstream [`NOTICE.md`](https://github.com/Edisuism/CurveSurfacingLib-V2/blob/master/NOTICE.md) |
| `curve-surfacing-core/external/Directional/` | [avaxman/Directional](https://github.com/avaxman/Directional) | vendored via Curve-Surfacing-Core | MPL-2.0 |
| `eigen-3.4.0/` | [Eigen 3.4.0](https://eigen.tuxfamily.org) | 3.4.0 | MPL-2.0 |
| `libigl/` | [libigl/libigl](https://github.com/libigl/libigl) | 2.4.0 | MPL-2.0 |
| `triangle/` | [libigl/triangle](https://github.com/libigl/triangle) (Shewchuk's Triangle) | `3ee6cac2230f0fe1413879574f741c7b6da11221` | Non-commercial, notices must stay intact |
| `imgui/` | [ocornut/imgui](https://github.com/ocornut/imgui) | `9f7d8086` (v1.92.9 WIP) | MIT |

`curve-surfacing-core/include/curve_surfacing/CycleDiscovery/tetgen/` is TetGen
1.4.x (Hang Si, WIAS), non-commercial per its source headers.
