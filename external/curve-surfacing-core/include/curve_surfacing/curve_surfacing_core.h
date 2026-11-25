#pragma once

#include <vector>
#include <string>
#include <optional>
#include <Eigen/Dense>
#include "app_loader.h"
#include "Mesh3D.h"
#include "CurvNet.h"

namespace CurveSurfacing {

// Solver type for least-squares problems in curvature smoothing
enum class LeastSquaresSolverType {
    CHOLESKY,  // Fast Cholesky decomposition (requires positive definite matrix, may fail on dense/irregular meshes)
    LSCG       // Robust iterative Least-Squares Conjugate Gradient solver (handles rank-deficient systems, slower but more stable)
};

// Vertex-relaxation variant used by the per-iteration remeshing step
// (paper Sec. 5.4): redistributes vertices tangentially before the
// Delaunay pass so triangle quality survives repeated iterations.
enum class RelaxationType {
    Tangential,           // In-house damped Laplacian with first-order tangent-plane
                          // projection back onto the local surface. The default.
    TflSmooth,            // TracingFieldLib's SmoothLaplacian: damped uniform Laplacian
                          // with EXACT reprojection onto the pre-smoothing surface
                          // (spatial-grid closest point). Experimental alternative.
    TflSmoothSlideCurves  // TflSmooth + interior curve vertices additionally slide
                          // ALONG their curve (1D feature smoothing, reprojected onto
                          // the original polyline) — redistributes dense boundary
                          // sampling without losing curve interpolation.
};

struct CurveSurfacingResult {
    bool        success = false;
    std::string error_message;

    // Comparison metrics (populated only when params.ref_mesh is set)
    double hausdorff_distance      = -1.0;
    double mean_distance           = -1.0;
    double symmetric_mean_distance = -1.0;
    double area_ratio              = -1.0;
    double curvature_ratio         = -1.0;

    // Mesh information
    int output_vertices = 0;
    int output_faces    = 0;

    // How the iteration ended: passes actually run, and the average vertex
    // movement of the last one as a fraction of the bounding-box diagonal
    // (the convergence measure; below 0.001 counts as converged).
    int    iterations_run = 0;
    double final_movement = 0.0;

    // Why the iteration stopped: "converged", "max iterations",
    // "solver not converged", "diverged: ..." (see stop_on_divergence).
    // For the last two the mesh is the one of the previous iteration.
    std::string stop_reason;
    bool   diverged          = false;
    double max_flipped_ratio = 0.0; // worst fraction of flipped free-vertex normals
    double time_ms           = 0.0;

    // Per cycle, per curve in the cycle's own order: 0 = CRVT_LINE (the
    // normals along it transported consistently and constrain the surface),
    // 1 = TRIM_LINE (they did not; the curve is a boundary only). The same
    // information the <output_file_prefix>.curve_types side file summarises
    // per curve, but per use and in memory.
    std::vector<std::vector<int>> cycle_curve_types;

    // Final mesh result
    Geex::MeshData output_mesh;
};

// ── Optional parameters for curve_surfacing_core and load_mesh_from_structures ──
//
// All fields have sensible defaults so callers only need to set what differs.
// Use std::optional for inputs whose *absence* is meaningful.
//
// Example — minimal call:
//   curve_surfacing_core(curves, cycles);
//
// Example — with patches and a reference mesh:
//   CurveSurfacingParams p;
//   p.patches_data = my_patches;
//   p.ref_mesh     = my_ref_mesh;
//   p.iterations   = 20;
//   curve_surfacing_core(curves, cycles, p);
//
struct CurveSurfacingParams {
    // ── Core behaviour ───────────────────────────────────────────────────────
    int  iterations     = 15;
    bool enable_logging = false;

    // Stop iterating once the average vertex movement drops below 0.1% of the
    // bounding-box diagonal (paper Sec. 6). Set false to always run the full
    // iteration count (e.g. for convergence experiments).
    bool stop_at_convergence = true;

    // Number of consecutive iterations the movement must stay below that
    // threshold to count as converged. The movement is not monotone: starting
    // from a flat initial surface it can dip below 0.1% for one iteration and
    // grow again while the patch is still inflating (1 = previous behaviour).
    int convergence_iterations = 2;

    // Curvature-variation control weight t ∈ [0,1) (paper Sec. 5.2, Fig. 4).
    // Governs how strongly the flow-field magnitudes are smoothed across each
    // patch — i.e. how far curvature prescribed at the boundary curves
    // propagates into patch interiors. 0.3 is the paper's default; raise it
    // (e.g. 0.9) when large patches should inflate to full roundness.
    double curvature_smooth_weight = 0.3;

    // Final fairing pass (paper Sec. 5.4, optional post-process): one global
    // Laplacian-regularization solve with surface-deviation penalization.
    // Smooths the ridge bands the local optimization leaves along network
    // curves. Curve and boundary vertices stay fixed, so curve interpolation
    // is preserved. fairing_deviation_weight balances shape preservation
    // (higher = closer to the unfaired surface, less smoothing).
    bool   final_fairing           = false;
    double fairing_deviation_weight = 0.3;

    LeastSquaresSolverType solver_type = LeastSquaresSolverType::LSCG;

    // Which vertex-relaxation variant the remeshing step uses (see enum).
    RelaxationType relaxation = RelaxationType::Tangential;

    // Preservation weight anchoring the fitted curvature of faces adjacent to
    // a network curve (paper Sec. 5.2 "known curvature field constraints").
    // These faces carry the RMF boundary conditions into the flow-field
    // smoothing — anchored too weakly, the smoothness term erases the boundary
    // seed and patch interiors never inflate. 50 is a robust default; lower it
    // (e.g. 10-25) to let the smoothness term win on noisy boundary fits.
    double boundary_anchor_weight = 50.0;

    // ── Initial surface ──────────────────────────────────────────────────────
    // Target edge length for the per-patch remesh of the auto-generated
    // initial surface (TriRemesh, borders fixed). <= 0 keeps the default:
    // 1.5 x the mean segment length of the cycle's curves.
    double remesh_target_length = 0.0;

    // Multiplier of the mean curve segment length used as remesh target when
    // remesh_target_length <= 0 (per patch, from the patch's own cycle).
    // Applies to auto-generated patches and, with remesh_patches, to
    // user-supplied ones. 1.5 is the previous hard-coded value.
    double remesh_edge_multiplier = 1.5;

    // Also remesh user-supplied patches (patches_data) with the same remesher
    // after loading them. Off by default: a user mesh is inserted verbatim,
    // which preserves its vertex order.
    bool remesh_patches = false;

    // Seed for the C library's rand(), which CycleDiscovery draws on; set
    // once at the start of every call so the same input gives the same
    // initial surface.
    unsigned random_seed = 1;

    // ── Normal classification ────────────────────────────────────────────────
    // A curve's supplied normals are trusted (CRVT_LINE) when the normal at
    // one end, transported along the curve, meets the other end's within
    // this angle in degrees (or within it of the opposite direction);
    // otherwise the curve is a TRIM_LINE and its normals are ignored.
    double rmf_consistency_angle = 10.0;

    // Skip that check: when normal_data is given, every curve is CRVT_LINE
    // and the supplied normals constrain the surface everywhere. For a
    // caller whose normals are a genuine surface normal field that may
    // legitimately twist along a curve.
    bool trust_normal_data = false;

    // ── Patch selection ──────────────────────────────────────────────────────
    // Empty vector = surface all patches
    std::vector<int> patch_ids;

    // ── Optional inputs (nullopt = feature not used) ─────────────────────────
    // Sharp edges to preserve during surfacing
    std::optional<Geex::FeatureEdgeData> feature_edges;

    // Reference mesh for hausdorff/mean-distance comparison metrics
    std::optional<Geex::MeshData> ref_mesh;

    // User-supplied patch meshes (skips SurfaceMaker auto-generation)
    std::optional<std::vector<Geex::MeshData>> patches_data;

    // Pre-computed RMF normals (skips internal RMF computation)
    std::optional<Geex::NormalCurveData> normal_data;

    // ── Initial surface ─────────────────────────────────────────────────────
    // Replace the interior of each (active) patch of the initial surface with
    // the thin-plate fill tangent to the border normals (TracingFieldLib
    // Geo::BiharmonicFillTangent; border normals = the curve normals, from
    // normal_data or the RMF): starts the iteration from an inflated surface
    // instead of the flat one of the triangulation. biharmonic_tension > 0
    // moves it towards the harmonic fill (less overshoot).
    bool   biharmonic_initial = false;
    double biharmonic_tension = 0.0;

    // ── Progress and safety ──────────────────────────────────────────────────
    // 0 = silent, 1 = one summary line per call, 2 = also one line per
    // iteration (independent of enable_logging, which is the full debug log)
    int progress_level = 0;

    // Stop when the iteration starts deforming the mesh instead of
    // converging, and keep the previous iterate:
    //  - more than divergence_flip_ratio of the free vertices have a normal
    //    flipped (> 90 deg) with respect to the initial surface, or
    //  - the movement grows by more than divergence_growth times the previous
    //    one (and is above 0.5% of the bounding-box diagonal)
    bool   stop_on_divergence    = true;
    double divergence_flip_ratio = 0.02;
    double divergence_growth     = 2.0;

    // ── Output ───────────────────────────────────────────────────────────────
    // Base path/prefix for side-car output files written by the core.
    // When non-empty, a <output_file_prefix>.curve_types file is written after
    // surfacing, listing each curve as CRVT_LINE or TRIM_LINE.
    // Example: "/path/to/project/my_model"  →  my_model.curve_types
    std::string output_file_prefix;
};

// ── Primary entry points ─────────────────────────────────────────────────────

CurveSurfacingResult curve_surfacing_core(
    const std::vector<Geex::CurveData>& curv_data,
    const std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params = {}
);

Geex::Mesh3D* load_mesh_from_structures(
    const std::vector<Geex::CurveData>& curv_data,
    const std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params = {}
);

// ── Utility ──────────────────────────────────────────────────────────────────

std::vector<std::vector<Eigen::Vector3d>> loadCurvesFromCycleFormat(const std::string& filename);
std::vector<std::vector<Eigen::Vector3d>> loadCurves(const std::string& filename);

// Extract RMF normals for each cycle's boundary vertices.
// Returns: vector[cycle_id][vertex_idx] = {position, RMF_normal}
std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>
extract_cycle_rmf_data(
    std::vector<Geex::CurveData>& curv_data,
    std::vector<Geex::CycleData>& cycle_data,
    std::vector<Geex::MeshData>&  patches_data,
    const Geex::FeatureEdgeData&  feature_edges = Geex::FeatureEdgeData()
);

} // namespace CurveSurfacing

#ifndef CURVE_SURFACING_STATIC
#  include "curve_surfacing_core.cpp"
#endif