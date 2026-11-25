#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "curve_surfacing_core.h"
#include "cs_logging.h"
#include "app_data.h"
#include "app_function.h"
#include "app_loader.h"
#include <cstdlib>
#ifdef CS_USE_DIRECTIONAL
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/TriMesh.h>
#endif
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <biharmonic_fill.h>   // TracingFieldLib: Geo::BiharmonicFillTangent
#include <chrono>
#include <iostream>

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

// Thin-plate fill of the interior of each active patch, tangent to the
// border normals stored on the curve half-edges (edge->tangent). Border =
// half-edges of the patch without a face of the same patch on the other side.
static void biharmonic_fill_patches(
    Geex::Mesh3D* mesh,
    const std::unordered_set<Geex::PolyGraph_Cycle*>& active_patches,
    double tension, Geex::CsLog& log)
{
    std::map<Geex::PolyGraph_Cycle*, std::vector<Geex::HE_face*>> patch_faces;
    for (ptrdiff_t f = 0; f < mesh->get_num_of_faces(); ++f) {
        Geex::HE_face* face = mesh->get_face(f);
        if (face && face->curve_cycle && active_patches.count(face->curve_cycle))
            patch_faces[face->curve_cycle].push_back(face);
    }

    int num_filled = 0;
    for (auto& entry : patch_faces) {
        Geex::PolyGraph_Cycle* cycle = entry.first;
        std::unordered_map<Geex::HE_vert*, int> local;
        std::vector<Geex::HE_vert*> verts;
        std::vector<Geo::Point3<double>> pos;
        std::vector<std::vector<int>> faces;
        for (Geex::HE_face* face : entry.second) {
            std::vector<int> fv;
            for (Geex::HE_vert* v : face->get_vertices()) {
                auto it = local.find(v);
                if (it == local.end()) {
                    it = local.insert({v, (int)verts.size()}).first;
                    verts.push_back(v);
                    pos.push_back(Geo::Point3<double>(v->pos.x(), v->pos.y(), v->pos.z()));
                }
                fv.push_back(it->second);
            }
            faces.push_back(fv);
        }

        // border vertices (fixed) with the sum of their border edge normals
        std::map<int, Geo::Point3<double>> border_normals;
        for (Geex::HE_face* face : entry.second) {
            // half-edges of the face in order (start vertex of hes[k] is the
            // end vertex of hes[k-1])
            std::vector<Geex::HE_edge*> hes;
            Geex::HE_edge* he = face->edge;
            Geex::HE_edge* start = he;
            do { hes.push_back(he); he = he->next; } while (he && he != start);
            for (size_t k = 0; k < hes.size(); ++k) {
                Geex::HE_edge* e = hes[k];
                const bool border = !e->pair || !e->pair->face ||
                                    e->pair->face->curve_cycle != cycle;
                if (!border)
                    continue;
                Geo::Point3<double> n(0, 0, 0);
                if (e->tangent.norm() > 1e-10)
                    n = Geo::Point3<double>(e->tangent.x(), e->tangent.y(), e->tangent.z());
                Geex::HE_vert* v0 = hes[(k + hes.size() - 1) % hes.size()]->vert;
                for (Geex::HE_vert* v : {v0, e->vert}) {
                    int id = local[v];
                    auto bn = border_normals.find(id);
                    if (bn == border_normals.end())
                        border_normals[id] = n;
                    else
                        bn->second = bn->second + n;
                }
            }
        }
        // curve vertices stay fixed even inside the patch
        for (size_t i = 0; i < verts.size(); ++i)
            if (mesh->vertex_on_curve(verts[i]) && !border_normals.count((int)i))
                border_normals[(int)i] = Geo::Point3<double>(0, 0, 0);

        if (!Geo::BiharmonicFillTangent<double>(pos, faces, border_normals, 10.0, tension))
            continue;
        for (size_t i = 0; i < verts.size(); ++i)
            if (!border_normals.count((int)i))
                verts[i]->pos = Eigen::Vector3d(pos[i].X, pos[i].Y, pos[i].Z);
        ++num_filled;
    }
    mesh->update_normal();
    log << "Thin-plate initial surface on " << num_filled << " patches" << std::endl;
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
    // vertex pointer -> index in output_mesh.vertices (replaces the linear
    // search per face corner below, which made this O(nv * nf))
    std::unordered_map<const Geex::HE_vert*, int> vert_index;
    vert_index.reserve(nv);
    for (ptrdiff_t i = 0; i < nv; ++i) {
        Geex::HE_vert* v = mesh->get_vertex(i);
        //if (v) result.output_mesh.vertices.push_back(v->pos);
        if (v) {
            vert_index[v] = (int)result.output_mesh.vertices.size();
            result.output_mesh.vertices.push_back(v->pos);
        }
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
            //for (ptrdiff_t j = 0; j < nv; ++j) {
            //    if (mesh->get_vertex(j) == he->vert) { vertex_index = j; break; }
            //}
            auto vit = vert_index.find(he->vert);
            if (vit != vert_index.end()) vertex_index = vit->second;
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

CS_INLINE CurveSurfacingResult curve_surfacing_core(
    const std::vector<Geex::CurveData>& curv_data,
    const std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params)
{
    CurveSurfacingResult result;
    const auto t_start = std::chrono::steady_clock::now();
    auto elapsed_ms = [&t_start]() {
        return std::chrono::duration<double, std::milli>(
                   std::chrono::steady_clock::now() - t_start).count();
    };
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
        std::srand(params.random_seed);
        Geex::CurvNet* curvNet = new Geex::CurvNet();
        curvNet->remesh_target_length_ = params.remesh_target_length;
        curvNet->remesh_patches_       = params.remesh_patches;
        curvNet->remesh_edge_multiplier_ = params.remesh_edge_multiplier;
        curvNet->rmf_consistency_angle_deg_ = params.rmf_consistency_angle;
        curvNet->trust_normal_data_    = params.trust_normal_data;

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
        if (params.progress_level >= 2)
            std::cout << "[CurveSurfacing] initial surface: " << cycle_data.size()
                      << " cycles, " << mesh->get_num_of_vertices() << " vertices, "
                      << mesh->get_num_of_faces() << " faces (" << (int)elapsed_ms()
                      << " ms)" << std::endl;

        // The classification the load step decided on, per cycle edge.
        for (int c = 0; c < curvNet->m_p_polygraph_->get_num_of_cycles(); ++c) {
            std::vector<int> types;
            Geex::PolyGraph_Cycle* cycle = curvNet->m_p_polygraph_->get_cycle(c);
            if (cycle && cycle->edge) {
                Geex::PolyGraph_Edge* e = cycle->edge;
                do {
                    types.push_back(e->curve_type == Geex::CURVE_TYPE::CRVT_LINE ? 0 : 1);
                    e = e->next;
                } while (e && e != cycle->edge);
            }
            result.cycle_curve_types.push_back(std::move(types));
        }

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

        // Optional inflated initial surface (thin-plate fill).
        if (params.biharmonic_initial)
            biharmonic_fill_patches(mesh, active_patches, params.biharmonic_tension, log);

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
#ifdef CS_USE_DIRECTIONAL
        directional::IntrinsicFaceTangentBundle tb;
        directional::TriMesh triMesh;
        tb = generateDirectionalCrossfield(delaunay_->m_init_mesh_, tb, triMesh);
        log << "Fitting curvature tensor..." << std::endl;
        fit_curvature_tensor(delaunay_->m_init_mesh_, true, &tb);
#else
        generateCrossfield(delaunay_->m_init_mesh_);
        log << "Fitting curvature tensor..." << std::endl;
        fit_curvature_tensor(delaunay_->m_init_mesh_, true);
#endif

        // Iteration loop.
        {
            std::vector<double> v_residuals, m_residuals, c_residuals;
            const Geex::LeastSquaresSolverType geo_solver =
                (params.solver_type == LeastSquaresSolverType::CHOLESKY)
                ? Geex::LeastSquaresSolverType::CHOLESKY
                : Geex::LeastSquaresSolverType::LSCG;

            log << "Running " << params.iterations << " surfacing iterations..." << std::endl;

            // Divergence check data: normals of the initial surface (vertex
            // set does not change, only connectivity via Delaunay flips) and
            // the positions at the start of each iteration, to go back to.
            mesh->update_normal();
            const ptrdiff_t nv_all = mesh->get_num_of_vertices();
            std::vector<Eigen::Vector3d> initial_normals(nv_all);
            std::vector<char> is_free(nv_all, 0);
            int num_free = 0;
            for (ptrdiff_t v = 0; v < nv_all; ++v) {
                Geex::HE_vert* vert = mesh->get_vertex(v);
                if (!vert) continue;
                initial_normals[v] = vert->normal;
                bool free_v = !mesh->is_on_boundary(vert) && !mesh->vertex_on_curve(vert);
                if (!surface_all_patches)
                    free_v = free_v && is_active_vertex[v] && !is_boundary_vertex[v];
                if (free_v) { is_free[v] = 1; ++num_free; }
            }
            std::vector<Eigen::Vector3d> iterate_start_pos(nv_all);
            auto save_positions = [&]() {
                for (ptrdiff_t v = 0; v < nv_all; ++v)
                    if (Geex::HE_vert* vert = mesh->get_vertex(v)) iterate_start_pos[v] = vert->pos;
            };
            auto restore_positions = [&]() {
                for (ptrdiff_t v = 0; v < nv_all; ++v)
                    if (Geex::HE_vert* vert = mesh->get_vertex(v)) vert->pos = iterate_start_pos[v];
                mesh->update_normal();
            };
            auto flipped_ratio = [&]() {
                if (num_free == 0) return 0.0;
                int flipped = 0;
                for (ptrdiff_t v = 0; v < nv_all; ++v) {
                    if (!is_free[v]) continue;
                    Geex::HE_vert* vert = mesh->get_vertex(v);
                    if (vert && vert->normal.dot(initial_normals[v]) < 0) ++flipped;
                }
                return double(flipped) / double(num_free);
            };
            double prev_movement = -1.0;
            int below_threshold = 0; // consecutive iterations under 0.1% bbox
            result.stop_reason = "max iterations";

            for (int itr = 0; itr < params.iterations; ++itr) {
                save_positions();
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
#ifdef CS_USE_DIRECTIONAL
                tb = generateDirectionalCrossfield(delaunay_->m_init_mesh_, tb, triMesh);
                fit_curvature_tensor(delaunay_->m_init_mesh_, true, &tb);
#else
                generateCrossfield(delaunay_->m_init_mesh_);
                fit_curvature_tensor(delaunay_->m_init_mesh_, true);
#endif
                smooth_curvature_variation(delaunay_->m_init_mesh_,
                                           params.curvature_smooth_weight, geo_solver,
                                           params.boundary_anchor_weight);

                double avg_movement = 0.0;
                int solver_its = 0;
                double solver_err = 0.0;
                const bool converged = match_target_curvature(
                    delaunay_->m_init_mesh_, v_residuals, m_residuals, c_residuals,
                    selected_vertices, &avg_movement, &solver_its, &solver_err);
                if (!converged) {
                    log << "  Solver did not converge at iteration " << (itr + 1)
                        << ". Stopping and using previous iteration result." << std::endl;
                    // the relaxation of this iteration already moved the
                    // vertices: really go back to the previous iterate
                    restore_positions();
                    result.stop_reason = "solver not converged";
                    if (params.progress_level >= 1)
                        std::cout << "[CurveSurfacing] WARNING: least-squares solver not converged at iteration "
                                  << (itr + 1) << " (" << solver_its << " its, error " << solver_err
                                  << "), keeping iteration " << itr << std::endl;
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
                // Divergence check: flipped normals or growing movement
                const double flip = flipped_ratio();
                result.max_flipped_ratio = std::max(result.max_flipped_ratio, flip);
                if (params.progress_level >= 2)
                    std::cout << "[CurveSurfacing]   it " << (itr + 1) << "/" << params.iterations
                              << "  move " << rel_movement * 100.0 << "% bbox"
                              << "  flipped " << flip * 100.0 << "%"
                              << "  solver " << solver_its << " its"
                              << "  (" << (int)elapsed_ms() << " ms)" << std::endl;
                if (params.stop_on_divergence) {
                    std::string why;
                    if (flip > params.divergence_flip_ratio)
                        why = "normals flipped on " + std::to_string(flip * 100.0) + "% of free vertices";
                    else if (prev_movement > 0 && rel_movement > 0.005 &&
                             rel_movement > params.divergence_growth * prev_movement)
                        why = "movement grew from " + std::to_string(prev_movement * 100.0) +
                              "% to " + std::to_string(rel_movement * 100.0) + "% of bbox";
                    if (!why.empty()) {
                        restore_positions();
                        result.diverged = true;
                        result.stop_reason = "diverged: " + why;
                        if (params.progress_level >= 1)
                            std::cout << "[CurveSurfacing] WARNING: iteration " << (itr + 1)
                                      << " diverging (" << why << "), keeping iteration "
                                      << itr << std::endl;
                        break;
                    }
                }
                prev_movement = rel_movement;

                result.iterations_run = itr + 1;
                result.final_movement = rel_movement;
                log << "  Iteration " << (itr + 1) << " complete (avg movement "
                    << rel_movement * 100.0 << "% of bbox)" << std::endl;
                // if (rel_movement < 0.001) {
                //     result.stop_reason = "converged";
                //     if (params.stop_at_convergence) {
                if (rel_movement < 0.001) ++below_threshold;
                else below_threshold = 0;
                if (below_threshold >= std::max(1, params.convergence_iterations)) {
                    result.stop_reason = "converged";
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
        result.time_ms = elapsed_ms();
        log << "Core curve surfacing completed successfully!" << std::endl;
        if (params.progress_level >= 1)
            std::cout << "[CurveSurfacing] " << cycle_data.size() << " cycles, "
                      << result.output_vertices << " vertices: " << result.iterations_run
                      << "/" << params.iterations << " iterations, last move "
                      << result.final_movement * 100.0 << "% bbox, max flipped "
                      << result.max_flipped_ratio * 100.0 << "%, " << result.stop_reason
                      << " (" << (int)result.time_ms << " ms)" << std::endl;

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

CS_INLINE Geex::Mesh3D* load_mesh_from_structures(
    const std::vector<Geex::CurveData>& curv_data,
    const std::vector<Geex::CycleData>& cycle_data,
    const CurveSurfacingParams&   params)
{
    Geex::cs_logging_enabled() = params.enable_logging;

    const Geex::FeatureEdgeData& feature_edges = resolve_feature_edges(params);
    const bool has_patches = params.patches_data.has_value();
    const bool has_normals = params.normal_data.has_value();

    //Geex::CurvNet* curvNet = new Geex::CurvNet();
    // Same seed and loader knobs as curve_surfacing_core, otherwise the
    // initial surface built here differs from the one the core starts from
    std::srand(params.random_seed);
    Geex::CurvNet* curvNet = new Geex::CurvNet();
    curvNet->remesh_target_length_      = params.remesh_target_length;
    curvNet->remesh_patches_            = params.remesh_patches;
    curvNet->remesh_edge_multiplier_    = params.remesh_edge_multiplier;
    curvNet->rmf_consistency_angle_deg_ = params.rmf_consistency_angle;
    curvNet->trust_normal_data_         = params.trust_normal_data;

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

CS_INLINE std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>>
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