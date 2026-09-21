#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "curve_surfacing_core.h"
#include "cs_logging.h"
#include "app_data.h"
#include "app_function.h"
#include "app_loader.h"
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/TriMesh.h>
#include <unordered_set>

namespace CurveSurfacing {

// ==============================================================================
// INTERNAL HELPERS
// ==============================================================================

// Resolve the feature-edges value from params (empty struct when absent).
static const Geex::FeatureEdgeData& resolve_feature_edges(const CurveSurfacingParams& params)
{
    static const Geex::FeatureEdgeData s_empty;
    return params.feature_edges.has_value() ? *params.feature_edges : s_empty;
}

// Classify which vertices belong to the active patch set and which sit on the
// shared boundary between active and inactive patches.
static void classify_vertices(
    Geex::Mesh3D* mesh,
    const std::unordered_set<Geex::PolyGraph_Cycle*>& active_patches,
    bool surface_all_patches,
    std::vector<bool>& is_active_vertex,
    std::vector<bool>& is_boundary_vertex,
    Geex::CsLog& log)
{
    const ptrdiff_t nv = mesh->get_num_of_vertices();
    is_active_vertex.assign(nv, false);
    is_boundary_vertex.assign(nv, false);

    // First pass: mark vertices that belong to at least one active patch face.
    for (ptrdiff_t f = 0; f < mesh->get_num_of_faces(); ++f) {
        Geex::HE_face* face = mesh->get_face(f);
        if (face && face->curve_cycle && active_patches.count(face->curve_cycle)) {
            for (Geex::HE_vert* v : face->get_vertices()) {
                if (v) is_active_vertex[v->id] = true;
            }
        }
    }

    if (surface_all_patches) return;

    // Second pass: a vertex is a boundary vertex if it touches both an active
    // and an inactive patch face.
    for (ptrdiff_t v = 0; v < nv; ++v) {
        if (!is_active_vertex[v]) continue;

        Geex::HE_vert* vert = mesh->get_vertex(v);
        if (!vert) continue;

        Geex::HE_edge* start_edge = vert->edge;
        if (!start_edge) continue;

        bool touches_active   = false;
        bool touches_inactive = false;

        Geex::HE_edge* cur = start_edge;
        do {
            if (cur->face && cur->face->curve_cycle) {
                if (active_patches.count(cur->face->curve_cycle))
                    touches_active   = true;
                else
                    touches_inactive = true;
            }
            cur = cur->pair;
            if (cur) cur = cur->next;
        } while (cur && cur != start_edge);

        if (touches_active && touches_inactive)
            is_boundary_vertex[v] = true;
    }

    // Logging stats.
    int num_active = 0, num_boundary = 0, num_optimizable = 0;
    for (ptrdiff_t v = 0; v < nv; ++v) {
        if (is_active_vertex[v]) {
            ++num_active;
            if (is_boundary_vertex[v]) ++num_boundary;
            else                       ++num_optimizable;
        }
    }
    log << "Vertex classification: " << num_active << " active ("
        << num_optimizable << " optimizable, " << num_boundary << " boundary fixed), "
        << (nv - num_active) << " inactive (fixed)" << std::endl;
}

// Extract the final mesh state from Mesh3D into a MeshData result.
static void extract_mesh_data(Geex::Mesh3D* mesh, CurveSurfacingResult& result)
{
    const ptrdiff_t nv = mesh->get_num_of_vertices();
    const ptrdiff_t nf = mesh->get_num_of_faces();

    result.output_mesh.num_vertices = nv;
    result.output_mesh.num_faces    = nf;

    result.output_mesh.vertices.clear();
    result.output_mesh.vertices.reserve(nv);
    for (ptrdiff_t i = 0; i < nv; ++i) {
        Geex::HE_vert* v = mesh->get_vertex(i);
        if (v) result.output_mesh.vertices.push_back(v->pos);
    }

    result.output_mesh.faces.clear();
    result.output_mesh.faces.reserve(nf);
    result.output_mesh.face_cycle_ids.clear();
    result.output_mesh.face_cycle_ids.reserve(nf);

    for (ptrdiff_t i = 0; i < nf; ++i) {
        Geex::HE_face* face = mesh->get_face(i);
        if (!face || !face->edge) continue;

        std::vector<int> face_indices;
        Geex::HE_edge* he    = face->edge;
        Geex::HE_edge* start = he;
        do {
            int vertex_index = -1;
            for (ptrdiff_t j = 0; j < nv; ++j) {
                if (mesh->get_vertex(j) == he->vert) { vertex_index = j; break; }
            }
            if (vertex_index >= 0) face_indices.push_back(vertex_index);
            he = he->next;
        } while (he && he != start);

        if (!face_indices.empty()) {
            result.output_mesh.faces.push_back(face_indices);
            result.output_mesh.face_cycle_ids.push_back(
                face->curve_cycle ? face->curve_cycle->id : -1);
        }
    }

    result.output_vertices = nv;
    result.output_faces    = nf;
}

// ==============================================================================
// CORE IMPLEMENTATION
// ==============================================================================

CurveSurfacingResult curve_surfacing_core(
    std::vector<Geex::CurveData>& curv_data,
    std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params)
{
    CurveSurfacingResult result;
    Geex::cs_logging_enabled() = params.enable_logging;
    Geex::CsLog log;  // gated by the flag above

    const Geex::FeatureEdgeData& feature_edges = resolve_feature_edges(params);
    if (!feature_edges.empty())
        log << "Feature edges loaded: " << feature_edges.edges.size()
            << " sharp edges to preserve" << std::endl;

    if (params.normal_data.has_value() && !params.normal_data->empty())
        log << "Normal data loaded: " << params.normal_data->num_cycles
            << " cycles with normal data for RMF" << std::endl;

    try {
        Geex::CurvNet* curvNet = new Geex::CurvNet();

        // Load mesh — four possible combinations of patches_data × normal_data.
        Geex::Mesh3D* mesh = nullptr;
        const bool has_patches = params.patches_data.has_value();
        const bool has_normals = params.normal_data.has_value();

        if (!has_patches && !has_normals) {
            mesh = curvNet->load_cycle_discovery_data_from_structures(
                curv_data, cycle_data, feature_edges);
        } else if (has_patches && !has_normals) {
            mesh = curvNet->load_cycle_discovery_data_from_structures(
                curv_data, cycle_data, *params.patches_data, feature_edges);
        } else if (!has_patches && has_normals) {
            mesh = curvNet->load_cycle_discovery_data_from_structures(
                curv_data, cycle_data, feature_edges, *params.normal_data);
        } else {
            // has_patches && has_normals
            mesh = curvNet->load_cycle_discovery_data_from_structures(
                curv_data, cycle_data, *params.patches_data, feature_edges, *params.normal_data);
        }

        if (!mesh) {
            result.error_message = "Failed to create mesh from structures";
            log << result.error_message << std::endl;
            delete curvNet;
            return result;
        }

        if (mesh->get_num_of_vertices() == 0 || mesh->get_num_of_faces() == 0) {
            result.error_message = cycle_data.empty()
                ? "No cycles provided: cycle data (.cycle) is required to build initial patches"
                : "Initial mesh is empty: no patch could be triangulated from the given cycles";
            Geex::cs_err() << "Error: " << result.error_message << std::endl;
            delete curvNet;
            return result;
        }

        log << "Successfully created mesh with " << mesh->get_num_of_vertices()
            << " vertices and " << mesh->get_num_of_faces() << " faces" << std::endl;

        // Build the active-patch set.
        std::unordered_set<Geex::PolyGraph_Cycle*> active_patches;
        const bool surface_all_patches = params.patch_ids.empty();

        if (surface_all_patches) {
            log << "Surfacing ALL patches" << std::endl;
            for (int i = 0; i < curvNet->m_p_polygraph_->get_num_of_cycles(); ++i)
                active_patches.insert(curvNet->m_p_polygraph_->get_cycle(i));
        } else {
            log << "Surfacing " << params.patch_ids.size() << " selected patches: ";
            for (int pid : params.patch_ids) {
                if (pid >= 0 && pid < curvNet->m_p_polygraph_->get_num_of_cycles()) {
                    active_patches.insert(curvNet->m_p_polygraph_->get_cycle(pid));
                    log << pid << " ";
                } else {
                    log << std::endl << "Warning: Invalid patch_id " << pid
                        << " (valid range: 0-"
                        << curvNet->m_p_polygraph_->get_num_of_cycles() - 1 << ")" << std::endl;
                }
            }
            log << std::endl;
        }

        // Vertex classification.
        std::vector<bool> is_active_vertex, is_boundary_vertex;
        classify_vertices(mesh, active_patches, surface_all_patches,
                          is_active_vertex, is_boundary_vertex, log);

        // Set up processing structures.
        Geex::APP_DataStructs* delaunay_ = new Geex::APP_DataStructs();
        delaunay_->m_init_mesh_   = mesh;
        delaunay_->m_p_curve_net_ = curvNet;

        // Optionally load reference mesh for comparison.
        if (params.ref_mesh.has_value() &&
            !params.ref_mesh->vertices.empty() && !params.ref_mesh->faces.empty())
        {
            Geex::MeshData ref_copy = *params.ref_mesh;
            delaunay_->load_reference_mesh_mod(ref_copy);
        }

        // Cross-field and initial curvature tensor.
        directional::IntrinsicFaceTangentBundle tb;
        directional::TriMesh triMesh;
        tb = generateDirectionalCrossfield(delaunay_->m_init_mesh_, tb, triMesh);
        log << "Fitting curvature tensor..." << std::endl;
        fit_curvature_tensor(delaunay_->m_init_mesh_, true, &tb);

        // Iteration loop.
        {
            std::vector<double> v_residuals, m_residuals, c_residuals;
            const Geex::LeastSquaresSolverType geo_solver =
                (params.solver_type == LeastSquaresSolverType::CHOLESKY)
                ? Geex::LeastSquaresSolverType::CHOLESKY
                : Geex::LeastSquaresSolverType::LSCG;

            log << "Running " << params.iterations << " surfacing iterations..." << std::endl;

            for (int itr = 0; itr < params.iterations; ++itr) {
                // Build the vertex selection for this iteration.
                std::vector<Geex::HE_vert*> selected_vertices;
                if (!surface_all_patches) {
                    for (ptrdiff_t v = 0; v < mesh->get_num_of_vertices(); ++v) {
                        if (is_active_vertex[v] && !is_boundary_vertex[v]) {
                            Geex::HE_vert* vert = mesh->get_vertex(v);
                            if (vert && !mesh->is_on_boundary(vert) && !mesh->vertex_on_curve(vert))
                                selected_vertices.push_back(vert);
                        }
                    }
                }

                // Remeshing step (paper Sec. 5.4): redistribute vertices
                // tangentially, then restore the Delaunay property. Keeps the
                // mesh a good discretization so the flow-field fit stays
                // reliable across iterations — without it, triangle quality
                // degrades every iteration, the curvature fit gets noisier,
                // and the process diverges on feature-rich models.
                // The relaxation variant is selected by params.relaxation
                // (see RelaxationType).
                delaunay_->m_init_mesh_->update_normal();
                if (params.relaxation == RelaxationType::Tangential)
                    tangential_relaxation(delaunay_->m_init_mesh_, selected_vertices);
                else
                    tfl_surface_relaxation(delaunay_->m_init_mesh_, selected_vertices,
                                           params.relaxation == RelaxationType::TflSmoothSlideCurves);

                bool topo = false;
                delaunay_->m_init_mesh_->make_delaunay(topo);
                delaunay_->m_init_mesh_->update_normal();
                tb = generateDirectionalCrossfield(delaunay_->m_init_mesh_, tb, triMesh);
                fit_curvature_tensor(delaunay_->m_init_mesh_, true, &tb);
                smooth_curvature_variation(delaunay_->m_init_mesh_,
                                           params.curvature_smooth_weight, geo_solver,
                                           params.boundary_anchor_weight);

                double avg_movement = 0.0;
                const bool converged = match_target_curvature(
                    delaunay_->m_init_mesh_, v_residuals, m_residuals, c_residuals,
                    selected_vertices, &avg_movement);
                if (!converged) {
                    log << "  Solver did not converge at iteration " << (itr + 1)
                        << ". Stopping and using previous iteration result." << std::endl;
                    break;
                }
                delaunay_->m_init_mesh_->update_normal();

                // Termination criterion from the paper (Sec. 6): stop when the
                // average vertex movement falls below 0.1% of the bounding-box
                // diagonal. The iteration has no convergence guarantee, so
                // running past this point degrades the surface instead of
                // improving it.
                Geex::Mesh3D* m = delaunay_->m_init_mesh_;
                const double bbox_diag =
                    (Eigen::Vector3d(m->xmax, m->ymax, m->zmax) -
                     Eigen::Vector3d(m->xmin, m->ymin, m->zmin)).norm();
                const double rel_movement =
                    (bbox_diag > 1e-12) ? avg_movement / bbox_diag : 0.0;
                log << "  Iteration " << (itr + 1) << " complete (avg movement "
                    << rel_movement * 100.0 << "% of bbox)" << std::endl;
                if (rel_movement < 0.001) {
                    if (params.stop_at_convergence) {
                        log << "  Average movement below 0.1% of bbox — converged, stopping."
                            << std::endl;
                        break;
                    }
                    log << "  Average movement below 0.1% of bbox (converged; continuing "
                           "because stop_at_convergence is off)." << std::endl;
                }
            }

            // Final fairing pass (paper Sec. 5.4, optional post-process):
            // smooths the ridge bands the local curvature optimization leaves
            // along the network curves. Curve vertices are fixed, so curve
            // interpolation is untouched.
            if (params.final_fairing) {
                log << "Final Laplacian fairing pass (deviation weight "
                    << params.fairing_deviation_weight << ")..." << std::endl;

                std::vector<Geex::HE_vert*> selected_vertices;
                if (!surface_all_patches) {
                    for (ptrdiff_t v = 0; v < mesh->get_num_of_vertices(); ++v) {
                        if (is_active_vertex[v] && !is_boundary_vertex[v]) {
                            Geex::HE_vert* vert = mesh->get_vertex(v);
                            if (vert && !mesh->is_on_boundary(vert) && !mesh->vertex_on_curve(vert))
                                selected_vertices.push_back(vert);
                        }
                    }
                }
                laplacian_fairing(delaunay_->m_init_mesh_,
                                  params.fairing_deviation_weight,
                                  geo_solver, selected_vertices);
                delaunay_->m_init_mesh_->update_normal();
            }
        }

        extract_mesh_data(delaunay_->m_init_mesh_, result);
        result.success = true;
        log << "Core curve surfacing completed successfully!" << std::endl;

        // Write per-curve type listing if a file prefix was supplied.
        if (!params.output_file_prefix.empty()) {
            delaunay_->save_curve_types(params.output_file_prefix);
        }

        // Cleanup — null out curvNet pointer first to prevent double-deletion.
        delaunay_->m_p_curve_net_ = nullptr;
        delete delaunay_;
        delete curvNet;

    } catch (const std::exception& e) {
        result.success       = false;
        result.error_message = std::string("Exception in curve_surfacing_core: ") + e.what();
        log << result.error_message << std::endl;
    }

    return result;
}

// ==============================================================================
// LOAD MESH
// ==============================================================================

Geex::Mesh3D* load_mesh_from_structures(
    std::vector<Geex::CurveData>& curv_data,
    std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params)
{
    Geex::cs_logging_enabled() = params.enable_logging;

    const Geex::FeatureEdgeData& feature_edges = resolve_feature_edges(params);
    const bool has_patches = params.patches_data.has_value();
    const bool has_normals = params.normal_data.has_value();

    Geex::CurvNet* curvNet = new Geex::CurvNet();

    Geex::Mesh3D* mesh = nullptr;
    if (!has_patches && !has_normals) {
        mesh = curvNet->load_cycle_discovery_data_from_structures(
            curv_data, cycle_data, feature_edges);
    } else if (has_patches && !has_normals) {
        mesh = curvNet->load_cycle_discovery_data_from_structures(
            curv_data, cycle_data, *params.patches_data, feature_edges);
    } else if (!has_patches && has_normals) {
        mesh = curvNet->load_cycle_discovery_data_from_structures(
            curv_data, cycle_data, feature_edges, *params.normal_data);
    } else {
        mesh = curvNet->load_cycle_discovery_data_from_structures(
            curv_data, cycle_data, *params.patches_data, feature_edges, *params.normal_data);
    }

    // Note: the caller is responsible for managing both the CurvNet and Mesh3D
    // objects. The mesh holds references into the CurvNet's internal structures,
    // so the CurvNet must remain alive for as long as the mesh is in use.
    return mesh;
}

// ==============================================================================
// EXTRACT CYCLE RMF DATA  (signature unchanged — not part of the new params API)
// ==============================================================================

std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>
extract_cycle_rmf_data(
    std::vector<Geex::CurveData>& curv_data,
    std::vector<Geex::CycleData>& cycle_data,
    std::vector<Geex::MeshData>&  patches_data,
    const Geex::FeatureEdgeData&  feature_edges)
{
    std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> result;
    result.resize(cycle_data.size());

    try {
        Geex::CurvNet* curvNet = new Geex::CurvNet();
        Geex::Mesh3D*  mesh    = nullptr;

        if (patches_data.empty()) {
            mesh = curvNet->load_cycle_discovery_data_from_structures(
                curv_data, cycle_data, feature_edges);
        } else {
            mesh = curvNet->load_cycle_discovery_data_from_structures(
                curv_data, cycle_data, patches_data, feature_edges);
        }

        if (!mesh) {
            std::cerr << "Failed to create mesh for RMF extraction" << std::endl;
            delete curvNet;
            return result;
        }

        for (size_t cycle_id = 0; cycle_id < cycle_data.size(); ++cycle_id) {
            Geex::PolyGraph_Cycle* pg_cycle = curvNet->m_p_polygraph_->get_cycle(cycle_id);
            if (!pg_cycle || !pg_cycle->edge) continue;

            std::map<Geex::HE_vert*, std::vector<Eigen::Vector3d>> vertex_normals;

            Geex::PolyGraph_Edge* pg_edge = pg_cycle->edge;
            do {
                for (int e = 0; e < mesh->get_num_of_edges(); ++e) {
                    Geex::HE_edge* edge = mesh->get_edge(e);
                    if (edge && edge->cycle_curv == pg_edge && edge->vert) {
                        Eigen::Vector3d normal = edge->tangent;
                        if (normal.norm() > 1e-10)
                            vertex_normals[edge->vert].push_back(normal);
                    }
                }
                pg_edge = pg_edge->next;
            } while (pg_edge != pg_cycle->edge);

            for (const auto& entry : vertex_normals) {
                Geex::HE_vert* vert = entry.first;
                Eigen::Vector3d avg = Eigen::Vector3d::Zero();
                for (const auto& n : entry.second) avg += n.normalized();
                if (!entry.second.empty() && avg.norm() > 1e-10) {
                    avg.normalize();
                    result[cycle_id].push_back({vert->pos, avg});
                } else {
                    result[cycle_id].push_back({vert->pos, Eigen::Vector3d(0, 0, 0)});
                }
            }
        }

        delete mesh;
        delete curvNet;

    } catch (const std::exception& e) {
        std::cerr << "Exception in extract_cycle_rmf_data: " << e.what() << std::endl;
    }

    return result;
}

} // namespace CurveSurfacing