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

    // ── Output ───────────────────────────────────────────────────────────────
    // Base path/prefix for side-car output files written by the core.
    // When non-empty, a <output_file_prefix>.curve_types file is written after
    // surfacing, listing each curve as CRVT_LINE or TRIM_LINE.
    // Example: "/path/to/project/my_model"  →  my_model.curve_types
    std::string output_file_prefix;
};

// ── Primary entry points ─────────────────────────────────────────────────────

CurveSurfacingResult curve_surfacing_core(
    std::vector<Geex::CurveData>& curv_data,
    std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params = {}
);

Geex::Mesh3D* load_mesh_from_structures(
    std::vector<Geex::CurveData>& curv_data,
    std::vector<Geex::CycleData>& cycle_data,
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