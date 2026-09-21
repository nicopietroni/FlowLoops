#pragma once

#include "app_data.h"
#include "least_squares_solver.h"
#include "constraint.h"

#include <vector>

// Forward-declare Directional types so consumers don't need directional/ on their include path.
// The full definitions are only required inside app_function.cpp.
namespace directional {
    struct TriMesh;
    class IntrinsicFaceTangentBundle;
}

// Crossfield generation
directional::IntrinsicFaceTangentBundle generateDirectionalCrossfield(
    Geex::Mesh3D* mesh,
    directional::IntrinsicFaceTangentBundle& tb,
    directional::TriMesh& triMesh);

// Cache for optimized (repeated) crossfield generation
struct CrossfieldCache {
    directional::TriMesh* triMesh;
    directional::IntrinsicFaceTangentBundle* tb;
    Eigen::VectorXi constFaces;
    Eigen::MatrixXd constVectors;
    bool initialized = false;

    CrossfieldCache();
    ~CrossfieldCache();
};

directional::IntrinsicFaceTangentBundle generateDirectionalCrossfieldOptimized(
    Geex::Mesh3D* mesh, CrossfieldCache& cache);

void facet_differential(const Geex::Mesh3D* pmesh, const Geex::HE_face* face,
                        const Eigen::Vector3d& u_axis, const Eigen::Vector3d& v_axis,
                        std::vector<Eigen::Vector2d>& dx_list, std::vector<Eigen::Vector2d>& dn_list,
                        const Geex::HE_face* prev_face = NULL);

void fit_curvature_tensor(Geex::Mesh3D* p_mesh, bool use_cross_field,
                          directional::IntrinsicFaceTangentBundle* tb);

// boundary_anchor_weight: preservation weight anchoring the fitted curvature
// of faces adjacent to a network curve (they carry the RMF boundary
// conditions, paper Sec. 5.2). See CurveSurfacingParams::boundary_anchor_weight.
void smooth_curvature_variation(Geex::Mesh3D* p_mesh, double curvature_smooth_weight_,
                                Geex::LeastSquaresSolverType solver_type = Geex::LeastSquaresSolverType::LSCG,
                                double boundary_anchor_weight = 50.0);

std::vector<Geex::LinearConstraint<double>> edgecurv_fillmatrix(
    Geex::HE_edge* edge, Eigen::Vector3d target_curvnorm, double coeff = 1.0);

std::vector<Geex::LinearConstraint<double>> fillmatrix_normalConstraint(
    Geex::HE_edge* edge, std::vector<Geex::HE_vert*>& selected_vertices);

Eigen::Vector3d compute_target_curvature(Geex::HE_edge* eh);

// Returns the average vertex displacement of the update.
double update_vertex_positions(Geex::Mesh3D* p_mesh,
                               std::vector<Geex::HE_vert*> map_var_vtx,
                               const std::vector<double>& solution);

// Tangential Laplacian relaxation (paper Sec. 5.4): redistribute free vertices
// toward their one-ring centroid, projected onto the tangent plane so the
// surface shape is preserved. Keeps the mesh a good discretization between
// adaptation steps; without it triangle quality degrades each iteration and
// the curvature fit / targets get progressively noisier.
// selected: if non-empty, only these vertices are relaxed (patch-subset mode).
void tangential_relaxation(Geex::Mesh3D* p_mesh,
                           const std::vector<Geex::HE_vert*>& selected,
                           int passes = 3, double step = 0.5);

// EXPERIMENTAL alternative to tangential_relaxation backed by TracingFieldLib's
// Geo::SmoothLaplacian: damped uniform Laplacian with EXACT reprojection onto
// the pre-smoothing surface (spatial-grid closest point, vs. our first-order
// tangent-plane projection). When slide_curves is true, curve vertices are
// additionally smoothed ALONG their curve (1D feature smoothing) and
// reprojected onto the original curve polyline — redistributing dense curve
// sampling without losing curve interpolation. Corner vertices stay fixed.
void tfl_surface_relaxation(Geex::Mesh3D* p_mesh,
                            const std::vector<Geex::HE_vert*>& selected,
                            bool slide_curves = false);

// Final fairing pass (paper Sec. 5.4): one global least-squares solve of
// Laplacian regularization with surface-deviation penalization. Smooths the
// ridge/fillet bands the local curvature optimization leaves along the
// network curves, while the deviation term (and the fixed curve vertices)
// preserve the overall shape and curve interpolation.
// selected: if non-empty, only these vertices are faired (patch-subset mode).
void laplacian_fairing(Geex::Mesh3D* p_mesh,
                       double deviation_weight,
                       Geex::LeastSquaresSolverType solver_type,
                       const std::vector<Geex::HE_vert*>& selected);

bool match_target_curvature(
    Geex::Mesh3D* p_mesh,
    std::vector<double>& v_residuals_out,
    std::vector<double>& m_residual_out,
    std::vector<double>& c_residual_out,
    std::vector<Geex::HE_vert*>& selected_vertices,
    double* avg_movement_out = nullptr);

#ifndef CURVE_SURFACING_STATIC
#  include "app_function.cpp"
#endif