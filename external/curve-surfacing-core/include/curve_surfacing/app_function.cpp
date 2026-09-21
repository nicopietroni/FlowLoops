#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "app_function.h"
#include "cs_logging.h"

#include <unordered_map>
#include <smooth_mesh.h>   // TracingFieldLib: Geo::SmoothLaplacian
#include <directional/principal_matching.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/power_field.h>
#include <directional/CartesianField.h>
#include <directional/power_to_raw.h>

// CrossfieldCache owns its Directional members on the heap so the header
// can forward-declare those types without pulling in directional/ headers.
inline CrossfieldCache::CrossfieldCache()
    : triMesh(new directional::TriMesh())
    , tb(new directional::IntrinsicFaceTangentBundle())
    , initialized(false)
{}

inline CrossfieldCache::~CrossfieldCache()
{
    delete triMesh;
    delete tb;
}


// Function to generate a proper 4-RoSy crossfield using Directional
directional::IntrinsicFaceTangentBundle generateDirectionalCrossfield(Geex::Mesh3D* mesh, directional::IntrinsicFaceTangentBundle& tb, directional::TriMesh& triMesh) {
    // Convert Mesh3D to Eigen matrices
    auto* verts = mesh->get_vertices_list();
    auto* faces_list = mesh->get_faces_list();
    int num_vertices = verts->size();
    int num_faces = faces_list->size();
    Eigen::MatrixXd V(num_vertices, 3);
    Eigen::MatrixXi F(num_faces, 3);
    for (int i = 0; i < num_vertices; ++i) {
        V.row(i) = Eigen::Vector3d((*verts)[i]->pos.x(), (*verts)[i]->pos.y(), (*verts)[i]->pos.z());
    }
    for (int i = 0; i < num_faces; ++i) {
        auto face_verts = (*faces_list)[i]->get_vertices();
        for (int j = 0; j < 3; ++j) {
            F(i, j) = face_verts[j]->id;
        }
    }
    triMesh.set_mesh(V, F);
    tb.init(triMesh);
    std::vector<Eigen::Vector3d> boundaryDirections;

    std::vector<int> constFaceslist;
    std::vector<Eigen::Vector3d> constVectorslist;
    for (int i=0;i<triMesh.EF.rows();i++){
        if ((triMesh.EF(i,1)==-1)||(triMesh.EF(i,0)==-1)){  //boundary
        int face = (triMesh.EF(i,1)==-1 ? triMesh.EF(i,0) : triMesh.EF(i,1));
        constFaceslist.push_back(face);
        constVectorslist.push_back((triMesh.V.row(triMesh.EV(i,0))-triMesh.V.row(triMesh.EV(i,1))).normalized());
        continue;
        }
    }

    Eigen::VectorXi constFaces;
    Eigen::MatrixXd constVectors;
    constFaces.resize(constFaceslist.size());
    constVectors.resize(constVectorslist.size(),3);
    for (int i=0;i<constFaces.size();i++){
    constFaces(i)=constFaceslist[i];
    constVectors.row(i)=constVectorslist[i];
    }

    //generating the viewing fields
    int N = 4;
    Eigen::MatrixXd rawFieldConstraints=Eigen::MatrixXd::Zero(triMesh.F.rows(),N*3);
    Eigen::VectorXi posInFace=Eigen::VectorXi::Zero(triMesh.F.rows());
    for (int i=0;i<constFaces.size();i++){
        rawFieldConstraints.block(constFaces(i),3*posInFace(constFaces(i)), 1,3)=constVectors.row(i);
        posInFace(constFaces(i))++;
    }
  
    //Just to show the other direction if N is even, since we are by default constraining it
    if (N%2==0)
        rawFieldConstraints.middleCols(rawFieldConstraints.cols()/2, rawFieldConstraints.cols()/2)=-rawFieldConstraints.middleCols(0, rawFieldConstraints.cols()/2);
    
    directional::CartesianField constraintsField;
    constraintsField.init(tb, directional::fieldTypeEnum::RAW_FIELD, N);
    constraintsField.set_extrinsic_field(rawFieldConstraints);

    // Generate 4-RoSy field using power_field
    directional::CartesianField powerField;
    directional::power_field(tb, constFaces, constVectors, Eigen::VectorXd::Constant(constFaces.size(),1.0), 4, powerField);

    // Convert to raw representation (actual 3D vectors)
    directional::CartesianField rawField;
    directional::power_to_raw(powerField, 4, rawField, true);
    assert(rawField.extField.cols() == 3 * 4); // 12 for 4-RoSy -> true
    directional::principal_matching(rawField); // fills rawField.matching

    const directional::TangentBundle* ftb = rawField.tb;
    const Eigen::VectorXi &matching = rawField.matching;
    const Eigen::MatrixXi &adj = ftb->adjSpaces; // adj.rows() == matching.size()
    // For each tangent bundle adjacency (edge between faces)
    for (int a = 0; a < matching.size(); ++a) {
        int f0 = adj(a, 0);
        int f1 = adj(a, 1);
        int rot = matching(a); // Number of 90-degree rotations needed

        // Find the corresponding half-edge in your mesh
        Geex::HE_face* face0 = mesh->get_face(f0);
        Geex::HE_face* face1 = mesh->get_face(f1);
        if (!face0 || !face1) continue;

        // For each half-edge of face0, check if its pair's face is face1
        Geex::HE_edge* start = face0->edge;
        Geex::HE_edge* eh = start;
        do {
            if (eh->pair && eh->pair->face == face1) {
                eh->pij = rot;
                // Also set the opposite half-edge's pij (should be the inverse rotation)
                // Use proper modulo to avoid overflow: ((4 - (rot % 4)) + 4) % 4
                eh->pair->pij = ((4 - (rot % 4)) + 4) % 4;
                break; // Found the matching half-edge
            }
            eh = eh->next;
        } while (eh != start);
    }

    // Store the crossfield angles in face->theta using Directional's canonical tangent basis
    for (int fi = 0; fi < num_faces; ++fi) {
        auto f = (*faces_list)[fi];
        // Get canonical tangent basis from Directional's tangent bundle
        Eigen::Vector3d b1 = tb.mesh->FBx.row(fi);
        Eigen::Vector3d b2 = tb.mesh->FBy.row(fi);
        for (int k = 0; k < 4; ++k) {
            Eigen::Vector3d dir_world = rawField.extField.row(fi).segment(3*k, 3).transpose();
            if (dir_world.squaredNorm() < 1e-18) {
                f->theta[k] = 0.0;
                continue;
            }
            dir_world.normalize();
            double x = dir_world.dot(b1);
            double y = dir_world.dot(b2);
            double angle = std::atan2(y, x);
            f->theta[k] = angle;
        }
    }

    return tb; // Return the tangent bundle containing the crossfield
}

// Optimized version that caches mesh structures between iterations
directional::IntrinsicFaceTangentBundle generateDirectionalCrossfieldOptimized(
    Geex::Mesh3D* mesh, CrossfieldCache& cache) {

    // Initialize cache members explicitly
    if (!cache.initialized) {
        cache.constFaces = Eigen::VectorXi();
        cache.constVectors = Eigen::MatrixXd();
    }
    
    auto* verts = mesh->get_vertices_list();
    auto* faces_list = mesh->get_faces_list();
    int num_vertices = verts->size();
    int num_faces = faces_list->size();

    // Convert Mesh3D to Eigen matrices (always needed for position updates)
    Eigen::MatrixXd V(num_vertices, 3);
    Eigen::MatrixXi F(num_faces, 3);
    
    for (int i = 0; i < num_vertices; ++i) {
        V.row(i) = Eigen::Vector3d((*verts)[i]->pos.x(), (*verts)[i]->pos.y(), (*verts)[i]->pos.z());
    }
    for (int i = 0; i < num_faces; ++i) {
        auto face_verts = (*faces_list)[i]->get_vertices();
        for (int j = 0; j < 3; ++j) {
            F(i, j) = face_verts[j]->id;
        }
    }
    
    if (!cache.initialized) {
        // First time initialization
        cache.triMesh->set_mesh(V, F);
        cache.tb->init(*cache.triMesh);
        
        // Cache boundary constraints (topology-dependent, doesn't change between iterations)
        std::vector<int> constFaceslist;
        std::vector<Eigen::Vector3d> constVectorslist;
        
        for (int i = 0; i < cache.triMesh->EF.rows(); i++) {
            if ((cache.triMesh->EF(i,1) == -1) || (cache.triMesh->EF(i,0) == -1)) {  // boundary
                int face = (cache.triMesh->EF(i,1) == -1 ? cache.triMesh->EF(i,0) : cache.triMesh->EF(i,1));
                constFaceslist.push_back(face);
                constVectorslist.push_back((cache.triMesh->V.row(cache.triMesh->EV(i,0)) - 
                                           cache.triMesh->V.row(cache.triMesh->EV(i,1))).normalized());
            }
        }
        
        cache.constFaces.resize(constFaceslist.size());
        cache.constVectors.resize(constVectorslist.size(), 3);
        for (int i = 0; i < cache.constFaces.size(); i++) {
            cache.constFaces(i) = constFaceslist[i];
            cache.constVectors.row(i) = constVectorslist[i];
        }
        
        cache.initialized = true;
    } else {
        // Update vertex positions in the TriMesh (much faster than reconstruction)
        cache.triMesh->V = V;
        // Reinitialize the tangent bundle with updated vertex positions (still faster than full rebuild)
        cache.tb->init(*cache.triMesh);
        
        // Update boundary constraint vectors (directions may change due to vertex movement)
        for (int i = 0; i < cache.constFaces.size(); i++) {
            int face_id = cache.constFaces(i);
            // Find the boundary edge for this face
            for (int ei = 0; ei < cache.triMesh->EF.rows(); ei++) {
                if ((cache.triMesh->EF(ei,1) == -1 && cache.triMesh->EF(ei,0) == face_id) ||
                    (cache.triMesh->EF(ei,0) == -1 && cache.triMesh->EF(ei,1) == face_id)) {
                    cache.constVectors.row(i) = (cache.triMesh->V.row(cache.triMesh->EV(ei,0)) - 
                                                 cache.triMesh->V.row(cache.triMesh->EV(ei,1))).normalized();
                    break;
                }
            }
        }
    }
    
    // Generate crossfield using cached structures
    directional::CartesianField powerField;
    directional::power_field(*cache.tb, cache.constFaces, cache.constVectors, 
                           Eigen::VectorXd::Constant(cache.constFaces.size(), 1.0), 4, powerField);
    
    // Convert to raw representation
    directional::CartesianField rawField;
    directional::power_to_raw(powerField, 4, rawField, true);
    directional::principal_matching(rawField);
    
    // Store matching information in mesh half-edges
    const directional::TangentBundle* ftb = rawField.tb;
    const Eigen::VectorXi& matching = rawField.matching;
    const Eigen::MatrixXi& adj = ftb->adjSpaces;

    for (int a = 0; a < matching.size(); ++a) {
        int f0 = adj(a, 0);
        int f1 = adj(a, 1);
        int rot = matching(a);
        
        Geex::HE_face* face0 = mesh->get_face(f0);
        Geex::HE_face* face1 = mesh->get_face(f1);
        if (!face0 || !face1) continue;
        
        Geex::HE_edge* start = face0->edge;
        Geex::HE_edge* eh = start;
        do {
            if (eh->pair && eh->pair->face == face1) {
                eh->pij = rot;
                eh->pair->pij = (4 - rot) % 4;
                break;
            }
            eh = eh->next;
        } while (eh != start);
    }
    
    // Store crossfield angles in faces
    for (int fi = 0; fi < num_faces; ++fi) {
    auto f = (*faces_list)[fi];

    for (int k = 0; k < 4; ++k) {
        f->theta[k] = 0.0;  // Initialize to safe default
    }
    
    // Add bounds checking for tangent bundle access
    if (fi >= cache.tb->mesh->FBx.rows() || fi >= cache.tb->mesh->FBy.rows()) {
        continue;
    }
    
    Eigen::Vector3d b1 = cache.tb->mesh->FBx.row(fi);
    Eigen::Vector3d b2 = cache.tb->mesh->FBy.row(fi);
    
    for (int k = 0; k < 4; ++k) {
        // Add bounds checking for rawField access
        if (fi >= rawField.extField.rows() || (3*k + 2) >= rawField.extField.cols()) {
            continue;
        }
        
        Eigen::Vector3d dir_world = rawField.extField.row(fi).segment(3*k, 3).transpose();
        if (dir_world.squaredNorm() < 1e-18) {
            f->theta[k] = 0.0;
            continue;
        }
        dir_world.normalize();
        double x = dir_world.dot(b1);
        double y = dir_world.dot(b2);
        double angle = std::atan2(y, x);
        
        // Validate angle before assignment
        if (std::isfinite(angle)) {
            f->theta[k] = angle;
        }
    }
}
    return *cache.tb;
}

static double clamped_cot(const Eigen::Vector3d& a, const Eigen::Vector3d& b);

void facet_differential(const Geex::Mesh3D* pmesh, const Geex::HE_face* face,
    const Eigen::Vector3d& u_axis, const Eigen::Vector3d& v_axis, std::vector<Eigen::Vector2d>& dx_list, std::vector<Eigen::Vector2d>& dn_list, const Geex::HE_face* prev_face) {

    Geex::HE_edge* eh = face->edge;
    do
    {
        // Boundary curve edge: no neighbouring face, but the RMF frame gives
        // the KNOWN surface normal at the curve (paper Sec. 5.2 "known
        // curvature field constraints"). Feed it into the fit as a sample —
        // without this the boundary tilt prescribed by the RMF normals never
        // enters the flow-field magnitudes, the field stays ~0 and patch
        // interiors never inflate.
        if (!eh->pair->face && eh->cycle_curv && eh->tangent.squaredNorm() > 1e-16)
        {
            // The variation of the RMF normals ALONG the curve is the known
            // normal curvature of the surface there (curves are curvature
            // lines), and unlike a normal-vs-face comparison it is mesh-
            // resolution independent. Sample it as dN/ds between this curve
            // edge and the next one sharing the end vertex.
            Eigen::Vector3d fn(face->normal.x(), face->normal.y(), face->normal.z());
            if (fn.squaredNorm() > 1e-16) {
                fn.normalize();

                Eigen::Vector3d n_this = eh->tangent.normalized();
                if (n_this.dot(fn) < 0) n_this = -n_this;

                // Find the adjacent boundary curve edge at the end vertex by
                // circulating its one-ring.
                Geex::HE_vert* v_end = eh->vert;
                Geex::HE_edge* next_curve = nullptr;
                Geex::HE_edge* c = v_end->edge;
                if (c) {
                    Geex::HE_edge* start_c = c;
                    do {
                        if (c != eh->pair && c->pair && !c->pair->face &&
                            c->cycle_curv && c->tangent.squaredNorm() > 1e-16) {
                            next_curve = c;
                            break;
                        }
                        c = c->pair;
                        if (!c) break;
                        c = c->next;
                    } while (c && c != start_c);
                }

                if (next_curve) {
                    Eigen::Vector3d n_next = next_curve->tangent.normalized();
                    if (n_next.dot(fn) < 0) n_next = -n_next;

                    Eigen::Vector3d mid_this = 0.5 * (eh->vert->pos + eh->prev->vert->pos);
                    Eigen::Vector3d mid_next = 0.5 * (next_curve->vert->pos + next_curve->prev->vert->pos);

                    Eigen::Vector3d dx_vec = mid_next - mid_this;   // along the curve
                    Eigen::Vector3d dn_vec = n_next - n_this;       // dN over that arc

                    Eigen::Vector2d dx(dx_vec.dot(u_axis), dx_vec.dot(v_axis));
                    Eigen::Vector2d dn(dn_vec.dot(u_axis), dn_vec.dot(v_axis));
                    if (dx.squaredNorm() > 1e-20) {
                        dx_list.push_back(dx);
                        dn_list.push_back(dn);
                    }
                }
            }
        }

        if (eh->pair->face && eh->pair->face != prev_face)
        {
            Eigen::Vector3d ddir;
            Eigen::Vector3d ndiff;

            double cot_weight = pmesh->ortho_dual_edge_data(eh, ddir, ndiff, true);
            eh->ortho_dual = ddir;

            Eigen::Vector2d dx(ddir.dot(u_axis), ddir.dot(v_axis));
            Eigen::Vector2d dn(ndiff.dot(u_axis), ndiff.dot(v_axis));

            dx_list.push_back(dx);
            dn_list.push_back(dn);
        }

        eh = eh->next;
    } while (eh != face->edge);
}


void fit_curvature_tensor(Geex::Mesh3D* p_mesh, bool use_cross_field, directional::IntrinsicFaceTangentBundle* tb) {
    //compute best fit curvature tensor per face.
    int processed_count = 0;
    int skipped_count = 0;

    for (int fitr = 0; fitr < p_mesh->get_num_of_faces(); fitr++)
    {
        Geex::HE_face* fh = p_mesh->get_face(fitr);

        if (!fh->tag) {
            skipped_count++;
            continue;
        }
        processed_count++;

        //initialize.
        fh->tensor_axis[0] = Eigen::Vector3d(0,0,0);
        fh->tensor_axis[1] = Eigen::Vector3d(0,0,0);
        fh->tensor_diag[0] = 0;
        fh->tensor_diag[1] = 0;

        // old edge-based fallback basis (kept for fallback)
        Eigen::Vector3d faxis = (fh->edge->vert->pos - fh->edge->prev->vert->pos).normalized();
        Eigen::Vector3d fn = Eigen::Vector3d(fh->normal.x(), fh->normal.y(), fh->normal.z());
        if (fn.squaredNorm() > 1e-16) fn.normalize();
        Eigen::Vector3d binormal = fn.cross(faxis).normalized();

        Eigen::Vector3d u_axis(0,0,0), v_axis(0,0,0);

        if (use_cross_field && tb != nullptr) {
            // Use TB basis (FBx, FBy) and face theta to construct u/v in world coords
            Eigen::Vector3d b1 = tb->mesh->FBx.row(fitr).transpose();
            Eigen::Vector3d b2 = tb->mesh->FBy.row(fitr).transpose();

            // defensive normalization
            if (b1.squaredNorm() > 1e-16) b1.normalize();
            if (b2.squaredNorm() > 1e-16) b2.normalize();

            // enforce right-handedness: (b1 x b2) should agree with face normal
            if (b1.cross(b2).dot(fn) < 0.0) {
                b2 = -b2;
            }

            // build u/v from theta (same basis used to compute theta)
            double phi = fh->theta[0];
            u_axis = (std::cos(phi) * b1 + std::sin(phi) * b2);
            v_axis = (std::cos(phi + 0.5 * M_PI) * b1 + std::sin(phi + 0.5 * M_PI) * b2);

            // project to tangent plane & normalize
            u_axis -= fn.dot(u_axis) * fn;
            v_axis -= fn.dot(v_axis) * fn;
            if (u_axis.squaredNorm() > 1e-16) u_axis.normalize();
            if (v_axis.squaredNorm() > 1e-16) v_axis.normalize();
        }
        else
        {
            // fallback: original behavior (edge-based axes)
            if (faxis.squaredNorm() < 1e-12) { // safe fallback if degenerate
                faxis = fn.unitOrthogonal();
            }
            faxis.normalize();
            binormal = fn.cross(faxis).normalized();

            if (use_cross_field) {
                double cos_theta = std::cos(fh->theta[0]), sin_theta = std::sin(fh->theta[0]);
                u_axis = cos_theta * faxis + sin_theta * binormal;
                v_axis = -sin_theta * faxis + cos_theta * binormal;
            } else {
                u_axis = faxis;
                v_axis = binormal;
            }

            // project to tangent and normalize
            u_axis -= fn.dot(u_axis) * fn; if (u_axis.squaredNorm() > 1e-16) u_axis.normalize();
            v_axis -= fn.dot(v_axis) * fn; if (v_axis.squaredNorm() > 1e-16) v_axis.normalize();
        }

        std::vector<Eigen::Vector2d> dx_list, dn_list;

        double lambda_1 = 0, lambda_2 = 0;

        facet_differential(p_mesh, fh, u_axis, v_axis, dx_list, dn_list);

        if (dx_list.size() >= 2)
        {
            if (use_cross_field)
            {
                // given cross field, the unkowns are the diagonal entries of tensor in u/v basis.
                Eigen::MatrixXd mat_A((int)dx_list.size(), 2);
                Eigen::VectorXd rhs_b((int)dx_list.size());

                for (int dx_itr = 0; dx_itr < (int)dx_list.size(); dx_itr++)
                {
                    Eigen::Vector2d u_dx = dx_list[dx_itr].normalized();
                    mat_A(dx_itr, 0) = u_dx.x() * dx_list[dx_itr].x();
                    mat_A(dx_itr, 1) = u_dx.y() * dx_list[dx_itr].y();
                    rhs_b(dx_itr) = u_dx.dot(dn_list[dx_itr]);
                }

                Eigen::Matrix2d mat_AtA = mat_A.transpose() * mat_A;
                Eigen::Vector2d vec_Atb = mat_A.transpose() * rhs_b;

                // Regularize for stability
                double reg = 1e-12;
                mat_AtA += reg * Eigen::Matrix2d::Identity();

                Eigen::Vector2d lambdas = mat_AtA.ldlt().solve(vec_Atb);
                lambda_1 = lambdas(0);
                lambda_2 = lambdas(1);

                // store world-space tensor axes & diag
                fh->tensor_axis[0] = Eigen::Vector3d(u_axis.x(), u_axis.y(), u_axis.z());
                fh->tensor_axis[1] = Eigen::Vector3d(v_axis.x(), v_axis.y(), v_axis.z());
                fh->tensor_diag[0] = lambda_1;
                fh->tensor_diag[1] = lambda_2;
            }
        }

        else {
            fh->tensor_axis[0] = u_axis;
            fh->tensor_axis[1] = v_axis;
            fh->tensor_diag[0] = 0.00001;
            fh->tensor_diag[1] = 0.00001;
        }
    }
}

void smooth_curvature_variation(Geex::Mesh3D* p_mesh, double curvature_smooth_weight_, Geex::LeastSquaresSolverType solver_type, double boundary_anchor_weight) {
    Geex::cs_log() << "  Starting curvature smoothing with weight " << curvature_smooth_weight_;
    if (solver_type == Geex::LeastSquaresSolverType::CHOLESKY) {
        Geex::cs_log() << " (Cholesky solver)..." << std::endl;
    } else {
        Geex::cs_log() << " (LSCG solver)..." << std::endl;
    }
    
    // Smooth the curvature tensor for reducing curvature variation
    if (curvature_smooth_weight_ <= 0.0) {
        Geex::cs_log() << "  Skipping smoothing (weight <= 0)" << std::endl;
        return;  // Nothing to do if weight is zero or negative
    }
        
    double smooth_weight = curvature_smooth_weight_;
    double smooth_coeff = sqrt(smooth_weight);
    double preserve_coeff = (1.0 - smooth_weight);
        
    const int face_num = p_mesh->get_num_of_faces();
    Geex::cs_log() << "  Processing " << face_num << " faces..." << std::endl;
        
    // Step 1: Set up variables for each face's principal curvatures
    int var_num = 0;
    std::vector<MyPair<int>> map_face_varid(face_num);
        
    for (int fitr = 0; fitr < face_num; fitr++) {
        Geex::HE_face* fh = p_mesh->get_face(fitr);
        fh->id = fitr;
            
        // Skip untagged faces
        if (!fh->tag) {
            map_face_varid[fitr] = MyPair<int>(-1, -1);
        } else {
            // Each face has two principal curvatures (2 variables)
            var_num += 2;
            map_face_varid[fitr] = MyPair<int>(var_num - 2, var_num - 1);
        }
    }
        
    // Used to track weights for preservation terms
    std::vector<double> lambda_weights(var_num, 0);
        
    // Container for all constraints (smoothness + preservation)
    std::vector<Geex::LinearConstraint<double>> all_constraints;
        
    // Step 2: Create smoothness constraints between adjacent faces
    p_mesh->reset_edges_tag(false);
    for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++) {
        Geex::HE_edge* eh = p_mesh->get_edge(eitr);
            
        // Only process each edge once and only if both faces are valid
        if (!eh->tag && eh->face && eh->pair->face && eh->face->tag && eh->pair->face->tag) {
            eh->tag = eh->pair->tag = true;
                
            double wij[2] = {0.0, 0.0};
                
            Geex::HE_face* fh[2];
            fh[0] = eh->face;
            fh[1] = eh->pair->face;

            int adj_axis_map[2];
            adj_axis_map[0] = mod(eh->pij,2);
            adj_axis_map[1] = (adj_axis_map[0]+1)%2;
                
            // Get dual edge information
            Eigen::Vector3d dual_edge, ndiff;
            p_mesh->ortho_dual_edge_data(eh, dual_edge, ndiff, true);
            double dual_len = dual_edge.norm();
                
            // Compute weights based on alignment of principal directions with dual edge
            if (dual_len < 1e-5) {  // cocircular case
                wij[0] = wij[1] = 1.0;
            } 
            else {
                    dual_edge = dual_edge.normalized();
                    Eigen::Vector3d pair_dual_edge, pair_ndiff;
                    p_mesh->ortho_dual_edge_data(eh->pair, pair_dual_edge, pair_ndiff, true);
                    pair_dual_edge = pair_dual_edge.normalized();
                    
                    for (int axis_itr = 0; axis_itr < 2; axis_itr++)
                    {
                        wij[axis_itr] = 0.5*( abs(dual_edge.dot(fh[0]->tensor_axis[axis_itr])) + abs(pair_dual_edge.dot(fh[1]->tensor_axis[adj_axis_map[axis_itr]])) );
                    }
                
            }
                
            // For each principal direction, create smoothness constraint between faces
            for (int axis_itr = 0; axis_itr < 2; axis_itr++) {
                double coeff = smooth_coeff*wij[axis_itr];
                double diff[2] = {coeff, -coeff};

                int varid[2] = {map_face_varid[fh[0]->id].id[axis_itr], map_face_varid[fh[1]->id].id[adj_axis_map[axis_itr]]};
                
                double cur_val[2] = { fh[0]->tensor_diag[axis_itr], fh[1]->tensor_diag[adj_axis_map[axis_itr]] };

                    
                // Only create constraints if both faces have valid variables
                if (varid[0] > -1 && varid[1] > -1) {
                    // Accumulate weights for preservation term
                    lambda_weights[varid[0]] += wij[axis_itr];
                    lambda_weights[varid[1]] += wij[axis_itr];
                        
                    // Create smoothness constraint: κ_i(f1) - κ_j(f2) = 0
                    // with weight = smooth_coeff * wij[axis_itr]
                    Geex::LinearConstraint<double> smoothConstraint;
                    smoothConstraint.Weight = smooth_coeff * wij[axis_itr];
                    smoothConstraint.RHSide = 0.0;
                        
                    // Left hand side: κ_i(f1) - κ_j(f2)
                    smoothConstraint.LHSide.push_back(Geex::VarTerm<double>(varid[0], 0, 1.0));
                    smoothConstraint.LHSide.push_back(Geex::VarTerm<double>(varid[1], 0, -1.0));
                        
                    all_constraints.push_back(smoothConstraint);
                } 
                // Handle case where one face has fixed value
                else if (varid[0] > -1 && varid[1] == -1) {
                    lambda_weights[varid[0]] += wij[axis_itr];
                        
                    // Constraint: κ_i(f1) = κ_j(f2) (fixed value from second face)
                    Geex::LinearConstraint<double> smoothConstraint;
                    smoothConstraint.Weight = smooth_coeff * wij[axis_itr];
                    smoothConstraint.RHSide = cur_val[1];  // Fixed value from face 2
                        
                    smoothConstraint.LHSide.push_back(Geex::VarTerm<double>(varid[0], 0, 1.0));
                    all_constraints.push_back(smoothConstraint);
                }
                else if (varid[0] == -1 && varid[1] > -1) {
                    lambda_weights[varid[1]] += wij[axis_itr];
                        
                    // Constraint: κ_j(f2) = κ_i(f1) (fixed value from first face)
                    Geex::LinearConstraint<double> smoothConstraint;
                    smoothConstraint.Weight = smooth_coeff * wij[axis_itr];
                    smoothConstraint.RHSide = cur_val[0];  // Fixed value from face 1
                        
                    smoothConstraint.LHSide.push_back(Geex::VarTerm<double>(varid[1], 0, 1.0));
                    all_constraints.push_back(smoothConstraint);
                }
            }
        }
    }
        
    // Step 3: Create preservation constraints for each face
    // Anchor strength for faces adjacent to a network curve. Their fitted
    // curvature carries the RMF boundary conditions (the paper's "known
    // curvature field constraints", Sec. 5.2) — they must be anchored hard,
    // not skipped, or the smoothness term erases the boundary seed and patch
    // interiors never receive any curvature to propagate.
    const double boundary_anchor = boundary_anchor_weight;

    for (int fitr = 0; fitr < face_num; fitr++) {
        MyPair<int>& var_ids = map_face_varid[fitr];
        Geex::HE_face* fh = p_mesh->get_face(fitr);

        // Does this face sit on a network curve with a valid RMF frame?
        bool on_curve_boundary = false;
        if (fh->edge) {
            Geex::HE_edge* he = fh->edge;
            do {
                if (he->pair && !he->pair->face && he->cycle_curv &&
                    he->tangent.squaredNorm() > 1e-16) {
                    on_curve_boundary = true;
                    break;
                }
                he = he->next;
            } while (he && he != fh->edge);
        }

        // Skip non-curve boundary faces (open boundaries without RMF data)
        if (!on_curve_boundary && p_mesh->is_on_boundary(fh)) {
            continue;
        }

        // For each principal direction
        for (int axis_itr = 0; axis_itr < 2; axis_itr++) {
            if (var_ids.id[axis_itr] > -1) {
                // Compute preservation weight based on accumulated smoothness weights
                double lweight = lambda_weights[var_ids.id[axis_itr]] * 0.5;
                lweight *= lweight;

                // Preservation constraint: κ_i(f) = κ_i_original(f)
                Geex::LinearConstraint<double> preserveConstraint;
                preserveConstraint.Weight = preserve_coeff * lweight;
                if (on_curve_boundary)
                    preserveConstraint.Weight = std::max(preserveConstraint.Weight,
                                                         preserve_coeff) * boundary_anchor;
                preserveConstraint.RHSide = fh->tensor_diag[axis_itr];

                preserveConstraint.LHSide.push_back(Geex::VarTerm<double>(var_ids.id[axis_itr], 0, 1.0));
                all_constraints.push_back(preserveConstraint);
            }
        }
    }
    
    Geex::cs_log() << "  Setting up solver with " << var_num << " variables and " 
              << all_constraints.size() << " constraints..." << std::endl;
    
    // Step 4: Set up and solve the linear system using Eigen
    Geex::LinearConstraintLeastSQSolver<double> solver(solver_type);
    Geex::cs_log() << "  Initializing solver..." << std::flush;
    solver.initialize_system(all_constraints);
    Geex::cs_log() << " done" << std::endl;
    
    // Initial solution: current curvature values
    std::vector<double> solution(var_num);
    for (int fitr = 0; fitr < face_num; fitr++) {
        MyPair<int>& var_ids = map_face_varid[fitr];
        for (int axis_itr = 0; axis_itr < 2; axis_itr++) {
            if (var_ids.id[axis_itr] > -1) {
                solution[var_ids.id[axis_itr]] = p_mesh->get_face(fitr)->tensor_diag[axis_itr];
            }
        }
    }
    
    // Solve the system
    Geex::cs_log() << "  Solving linear system..." << std::flush;
    if (!solver.solve_system(solution)) {
        Geex::cs_err() << "  Curvature smoothing solve failed; keeping unsmoothed curvatures" << std::endl;
        return;
    }
    Geex::cs_log() << " done" << std::endl;
    
    // Step 5: Update face curvature tensors with the solution
    for (int fitr = 0; fitr < face_num; fitr++) {
        MyPair<int>& var_ids = map_face_varid[fitr];
        
        for (int axis_itr = 0; axis_itr < 2; axis_itr++) {
            if (var_ids.id[axis_itr] > -1) {
                double smoothed_curv = solution[var_ids.id[axis_itr]];
                if (!std::isnan(smoothed_curv)) {
                    p_mesh->get_face(fitr)->tensor_diag[axis_itr] = smoothed_curv;
                }
            }
        }
    }
    
    Geex::cs_log() << "  Curvature smoothing complete" << std::endl;
}

// Cotangent of the angle between vectors a and b, clamped so that
// near-degenerate (needle) triangles cannot inject unbounded weights into the
// least-squares system. Unclamped, |a×b| → 0 sends the cot to Inf/NaN and a
// handful of such rows dominate the solve, exploding the whole surface
// (observed on espresso/toilet/beetle). Clamp ≈ cot(1°).
static double clamped_cot(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
    const double kCotClamp = 57.29;  // cot(1 degree)
    const double cross_norm = a.cross(b).norm();
    const double dot = a.dot(b);
    if (!(cross_norm > 0.0))  // degenerate or NaN
        return (dot >= 0.0) ? kCotClamp : -kCotClamp;
    const double cot = dot / cross_norm;
    return std::max(-kCotClamp, std::min(kCotClamp, cot));
}

std::vector<Geex::LinearConstraint<double>> edgecurv_fillmatrix(Geex::HE_edge* edge, Eigen::Vector3d target_curvnorm, double coeff){

    Geex::HE_vert* vts_list[4];

    vts_list[0] = edge->next->vert;
    vts_list[1] = edge->prev->vert;
    vts_list[2] = edge->pair->next->vert;
    vts_list[3] = edge->vert;

    // compute the cotangent of the angle between two vectors using the dot product divided by the norm of the cross product.
    double cot013 = clamped_cot(vts_list[3]->pos - vts_list[1]->pos,
                                vts_list[0]->pos - vts_list[1]->pos);

    double cot031 = clamped_cot(vts_list[0]->pos - vts_list[3]->pos,
                                vts_list[1]->pos - vts_list[3]->pos);

    double cot213 = clamped_cot(vts_list[2]->pos - vts_list[1]->pos,
                                vts_list[3]->pos - vts_list[1]->pos);

    double cot231 = clamped_cot(vts_list[2]->pos - vts_list[3]->pos,
                                vts_list[1]->pos - vts_list[3]->pos);
    
    double diff[4];
    diff[0] = -cot013 - cot031;
    diff[1] = cot031 + cot231;
    diff[2] = -cot213 - cot231;
    diff[3] = cot213 + cot013;

    std::vector<Geex::LinearConstraint<double>> constraints;

    // Equation 5. minimize squared differences between each edge curvature normal and target curvature normal
    // (diff_i*diff_0*coeff) x_{0,c} + (diff_i*diff_1*coeff) x_{1,c} + (diff_i*diff_2*coeff) x_{2,c} + (diff_i*diff_3*coeff) x_{3,c}
    // divide that whole equation by diff_i * coeff (when diff_i*coeff != 0)
    // where p0_c is the c-coordinate of vertex p0 and target_c is the c_component of target_curvnorm
    // diff0 * p0_c  +  diff1 * p1_c  +  diff2 * p2_c  +  diff3 * p3_c  =  target_c
    for (int coord_itr = 0; coord_itr < 3; coord_itr++) {
        for (int v_itr0 = 0; v_itr0 < 4; v_itr0++) { // 4 vertex frame
            if (vts_list[v_itr0]->var_id > -1) { // free variable
                Geex::LinearConstraint<double> constr;
                constr.RHSide = diff[v_itr0] * target_curvnorm[coord_itr] * coeff;
                constr.Weight = 1.0;

                for (int v_itr1 = 0; v_itr1 < 4; v_itr1++) {
                    if (vts_list[v_itr1]->var_id > -1) {
                        // Add matrix entry for free variable
                        Geex::VarTerm<double> term(vts_list[v_itr1]->var_id, coord_itr, diff[v_itr0] * diff[v_itr1] * coeff);
                        constr.LHSide.push_back(term);
                    } else {
                        // Add RHS contribution from fixed variable
                        constr.RHSide -= diff[v_itr0] * diff[v_itr1] * vts_list[v_itr1]->pos[coord_itr] * coeff;
                    }
                }
                constraints.push_back(constr);
            }
        }
    }

    return constraints;
}

//debugging by visualising varids
std::vector<Geex::LinearConstraint<double>> fillmatrix_normalConstraint(Geex::HE_edge* edge, std::vector<Geex::HE_vert*>& selected_vertices){
    std::vector<Geex::LinearConstraint<double>> constraints;
    if (edge && edge->face && edge->cycle_curv && (edge->cycle_curv->curve_type == Geex::CURVE_TYPE::CRVT_LINE || edge->cycle_curv->curve_type == Geex::CURVE_TYPE::TRIM_LINE))
    {
        Geex::HE_vert* vh = edge->next->vert;
        double coeff = 1.0;

        Eigen::Vector3d normal = edge->tangent;
        Eigen::Vector3d rhs = 2.0 * normal * (normal.dot(edge->vert->pos + edge->pair->vert->pos));
        Eigen::Vector3d lhs[3] = {
            4.0 * normal.x() * normal,
            4.0 * normal.y() * normal,
            4.0 * normal.z() * normal
        };

        for (int row = 0; row < 3; ++row) {
            Geex::LinearConstraint<double> constr;
            constr.Weight = 1.0; // tried 100, no difference
            constr.RHSide = rhs[row] * coeff;
            if (vh->var_id >= 0) {
                // Each constraint couples x, y, z of the vertex
                selected_vertices.push_back(vh);
                for (int col = 0; col < 3; ++col) {
                    Geex::VarTerm<double> var(vh->var_id, col, lhs[row][col]);
                    constr.LHSide.push_back(var);
                }
            }
            constraints.push_back(constr);
        }
    }

    return constraints;
}


Eigen::Vector3d compute_target_curvature(Geex::HE_edge* eh)
{
    Eigen::Vector2d dx0(eh->ortho_dual.dot(eh->face->tensor_axis[0]), eh->ortho_dual.dot(eh->face->tensor_axis[1]));
    Eigen::Vector2d dn0(eh->face->tensor_diag[0] * dx0[0], eh->face->tensor_diag[1] * dx0[1]);
    Eigen::Vector2d dx0_dir(0, 0);
    if (!std::isnan(dx0.norm()) && dx0.norm() > 1e-10) {
        dx0_dir = dx0.normalized();
    }
    double dn0_len = dn0.dot(dx0_dir);

    Eigen::Vector2d dx1(eh->pair->ortho_dual.dot(eh->pair->face->tensor_axis[0]), eh->pair->ortho_dual.dot(eh->pair->face->tensor_axis[1]));
    Eigen::Vector2d dn1(eh->pair->face->tensor_diag[0] * dx1[0], eh->pair->face->tensor_diag[1] * dx1[1]);
    Eigen::Vector2d dx1_dir(0, 0);
    if (!std::isnan(dx1.norm()) && dx1.norm() > 1e-10) {
        dx1_dir = dx1.normalized();
    }
    double dn1_len = dn1.dot(dx1_dir);

    double sin_dihedral = (dn0_len + dn1_len) / 2.0;

    if (std::isnan(sin_dihedral) || abs(sin_dihedral)>1)
    {
        sin_dihedral = 0;
    }

    double sin_half_dihedral = sign(sin_dihedral)*sqrt( abs(0.5*( 1 - sqrt(abs(1-sin_dihedral*sin_dihedral)) )) );
    eh->target_sin_hdi = eh->pair->target_sin_hdi = sin_half_dihedral;

    //compute target mean curvature normal at edges
    // Guard the normalization: for folded geometry the two face normals can be
    // near-opposite, and normalizing their (near-zero) sum yields NaN.
    Eigen::Vector3d normal_sum = eh->face->normal + eh->pair->face->normal;
    Eigen::Vector3d curv_norm = Eigen::Vector3d::Zero();
    if (normal_sum.squaredNorm() > 1e-16) {
        curv_norm = normal_sum.normalized();
        curv_norm *= sin_half_dihedral*2*(eh->vert->pos-eh->pair->vert->pos).norm();
    }

    return curv_norm;
}

// Curve-edge policy for the curvature-matching solve (paper Sec. 4/5: smooth
// transitions vs. sharp features). Returns true if the edge lies on a network
// curve classified as a sharp feature (TRIM_LINE) — such hinges must stay
// unconstrained so the crease can keep its dihedral angle. Zero-targeting them
// instead (the old behaviour) demands coplanarity across every crease and
// pushes whole panels away on feature-rich models (espresso, toilet, beetle).
inline bool edge_on_sharp_feature_curve(const Geex::HE_edge* eh)
{
    const Geex::PolyGraph_Edge* pge =
        eh->cycle_curv ? eh->cycle_curv : eh->pair->cycle_curv;
    return pge && pge->curve_type == Geex::CURVE_TYPE::TRIM_LINE;
}

// Returns the average vertex displacement of this update (used by callers to
// implement the paper's termination criterion: stop when the average movement
// is below 0.1% of the bounding-box diagonal).
double update_vertex_positions(Geex::Mesh3D* p_mesh, std::vector<Geex::HE_vert*> map_var_vtx, const std::vector<double>& solution)
{
    double avg_delta = 0.0;
    int valid_count = 0;

    // Compute bounding box diagonal
    Eigen::Vector3d bb_ll(p_mesh->xmin, p_mesh->ymin, p_mesh->zmin);
    Eigen::Vector3d bb_ur(p_mesh->xmax, p_mesh->ymax, p_mesh->zmax);
    double diagonal_len = (bb_ll - bb_ur).norm();

    // Per-iteration step cap: a well-conditioned solve moves vertices a tiny
    // fraction of the bbox per iteration, so anything beyond 10% of the
    // diagonal is a numerical outlier (e.g. a vertex ringed by needle
    // triangles), not signal. Clamp the step length, keep the direction.
    const double max_step = 0.1 * diagonal_len;

    for (int var_itr = 0; var_itr < map_var_vtx.size(); var_itr++) {
        Geex::HE_vert* vh = map_var_vtx[var_itr];
        Eigen::Vector3d new_pos(solution[3*var_itr], solution[3*var_itr+1], solution[3*var_itr+2]);

        if (!new_pos.allFinite()) {
            continue;
        }

        Eigen::Vector3d delta = new_pos - vh->pos;
        double step = delta.norm();
        if (diagonal_len > 1e-12 && step > max_step) {
            new_pos = vh->pos + delta * (max_step / step);
            step = max_step;
        }

        avg_delta += step;
        vh->pos = new_pos;
        valid_count++;
    }

    if (valid_count > 0) {
        avg_delta /= valid_count;
    }

    return avg_delta;
}

void tfl_surface_relaxation(Geex::Mesh3D* p_mesh,
                            const std::vector<Geex::HE_vert*>& selected,
                            bool slide_curves)
{
    const ptrdiff_t nv = p_mesh->get_num_of_vertices();
    const ptrdiff_t nf = p_mesh->get_num_of_faces();

    std::vector<bool> allowed(nv, selected.empty());
    for (Geex::HE_vert* vh : selected)
        if (vh && vh->id >= 0 && vh->id < nv) allowed[vh->id] = true;

    // Extract vertex/face arrays (and normalize vertex ids to list order).
    std::vector<Geo::Point3<double>> V(nv);
    for (ptrdiff_t i = 0; i < nv; ++i) {
        Geex::HE_vert* vh = p_mesh->get_vertex(i);
        vh->id = i;
        V[i] = Geo::Point3<double>(vh->pos.x(), vh->pos.y(), vh->pos.z());
    }
    std::vector<std::vector<int>> F;
    F.reserve(nf);
    for (ptrdiff_t f = 0; f < nf; ++f) {
        Geex::HE_face* fh = p_mesh->get_face(f);
        if (!fh || !fh->edge) continue;
        std::vector<int> tri;
        Geex::HE_edge* he = fh->edge;
        do { tri.push_back((int)he->vert->id); he = he->next; } while (he != fh->edge);
        if (tri.size() == 3) F.push_back(tri);
    }

    // Curve/boundary handling. Curve vertices with exactly two same-curve
    // boundary edges may slide along the curve (feature smoothing +
    // reprojection); everything else on a curve/boundary stays fixed.
    Geo::SmoothParam<double> SP;
    std::vector<std::pair<int,int>> feature_edges;                 // boundary curve edges
    std::vector<std::vector<Geex::HE_edge*>> vert_curve_edges(nv); // face-side boundary curve edges per vertex

    for (int e = 0; e < p_mesh->get_num_of_edges(); ++e) {
        Geex::HE_edge* eh = p_mesh->get_edge(e);
        if (!eh || !eh->face || !eh->pair || eh->pair->face) continue;  // face-side boundary only
        if (!eh->cycle_curv) continue;
        int a = (int)eh->prev->vert->id, b = (int)eh->vert->id;
        feature_edges.push_back({a, b});
        vert_curve_edges[a].push_back(eh);
        vert_curve_edges[b].push_back(eh);
    }

    for (ptrdiff_t i = 0; i < nv; ++i) {
        Geex::HE_vert* vh = p_mesh->get_vertex(i);
        const bool on_b = p_mesh->is_on_boundary(vh) || p_mesh->vertex_on_curve(vh);
        if (!allowed[i]) { SP.FixedVert.push_back((int)i); continue; }
        if (!on_b) continue;

        bool can_slide = false;
        if (slide_curves && vert_curve_edges[i].size() == 2) {
            Geex::HE_edge* e0 = vert_curve_edges[i][0];
            Geex::HE_edge* e1 = vert_curve_edges[i][1];
            // Same underlying curve on both sides → interior curve vertex,
            // not a corner.
            if (e0->cycle_curv && e1->cycle_curv &&
                e0->cycle_curv->curve_ == e1->cycle_curv->curve_)
                can_slide = true;
        }
        if (!can_slide) SP.FixedVert.push_back((int)i);
    }

    // Reference copies for exact reprojection (surface + feature polylines).
    std::vector<Geo::Point3<double>> RefV = V;
    std::vector<std::vector<int>>    RefF = F;

    SP.NumIte = slide_curves ? 10 : 3;
    SP.Damp   = 0.5;
    SP.SType  = Geo::SMLaplacian;
    SP.ReprojVert  = &RefV;
    SP.ReprojFaces = &RefF;
    if (slide_curves) {
        SP.Features       = feature_edges;
        SP.ReprojFeatures = feature_edges;
    }

    Geo::SmoothLaplacian(V, F, SP);

    std::vector<bool> is_fixed(nv, false);
    for (int idx : SP.FixedVert) is_fixed[idx] = true;
    for (ptrdiff_t i = 0; i < nv; ++i) {
        if (is_fixed[i]) continue;
        Eigen::Vector3d p(V[i].X, V[i].Y, V[i].Z);
        if (p.allFinite()) p_mesh->get_vertex(i)->pos = p;
    }
}

void tangential_relaxation(Geex::Mesh3D* p_mesh,
                           const std::vector<Geex::HE_vert*>& selected,
                           int passes, double step)
{
    const ptrdiff_t nv = p_mesh->get_num_of_vertices();

    // Restrict to the caller's selection when given (patch-subset surfacing).
    std::vector<bool> allowed(nv, selected.empty());
    for (Geex::HE_vert* vh : selected)
        if (vh && vh->id >= 0 && vh->id < nv) allowed[vh->id] = true;

    std::vector<Eigen::Vector3d> new_pos(nv);

    for (int pass = 0; pass < passes; ++pass) {
        std::vector<bool> moved(nv, false);

        for (ptrdiff_t vi = 0; vi < nv; ++vi) {
            Geex::HE_vert* vh = p_mesh->get_vertex(vi);
            if (!vh || !allowed[vi]) continue;
            if (p_mesh->is_on_boundary(vh) || p_mesh->vertex_on_curve(vh)) continue;

            Geex::HE_edge* start = vh->edge;
            if (!start) continue;

            // One-ring centroid.
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            int cnt = 0;
            Geex::HE_edge* cur = start;
            bool ring_ok = true;
            do {
                if (!cur->vert) { ring_ok = false; break; }
                centroid += cur->vert->pos;
                ++cnt;
                cur = cur->pair;
                if (!cur) { ring_ok = false; break; }
                cur = cur->next;
            } while (cur && cur != start);
            if (!ring_ok || cnt == 0) continue;
            centroid /= cnt;

            // Move toward the centroid, tangentially only: project out the
            // normal component so the shape is untouched.
            Eigen::Vector3d lap = centroid - vh->pos;
            Eigen::Vector3d n = vh->normal;
            if (n.squaredNorm() > 1e-16) {
                n.normalize();
                lap -= n.dot(lap) * n;
            }
            if (!lap.allFinite()) continue;

            new_pos[vi] = vh->pos + step * lap;
            moved[vi] = true;
        }

        for (ptrdiff_t vi = 0; vi < nv; ++vi)
            if (moved[vi]) p_mesh->get_vertex(vi)->pos = new_pos[vi];
    }
}

void laplacian_fairing(Geex::Mesh3D* p_mesh,
                       double deviation_weight,
                       Geex::LeastSquaresSolverType solver_type,
                       const std::vector<Geex::HE_vert*>& selected)
{
    const ptrdiff_t nv = p_mesh->get_num_of_vertices();

    // Free-vertex selection: same rule as match_target_curvature.
    std::vector<bool> allowed(nv, selected.empty());
    for (Geex::HE_vert* vh : selected)
        if (vh && vh->id >= 0 && vh->id < nv) allowed[vh->id] = true;

    std::vector<Geex::HE_vert*> map_var_vtx;
    for (ptrdiff_t vi = 0; vi < nv; ++vi) {
        Geex::HE_vert* vh = p_mesh->get_vertex(vi);
        if (!vh) continue;
        vh->var_id = -1;
        if (!allowed[vi]) continue;
        if (p_mesh->is_on_boundary(vh) || p_mesh->vertex_on_curve(vh)) continue;
        vh->var_id = (int)map_var_vtx.size();
        map_var_vtx.push_back(vh);
    }
    if (map_var_vtx.empty()) return;

    std::vector<Geex::LinearConstraint<double>> constraints;
    constraints.reserve(map_var_vtx.size() * 6);

    for (Geex::HE_vert* vh : map_var_vtx) {
        // Collect the one-ring.
        std::vector<Geex::HE_vert*> ring;
        Geex::HE_edge* start = vh->edge;
        if (!start) continue;
        Geex::HE_edge* cur = start;
        bool ok = true;
        do {
            if (!cur->vert) { ok = false; break; }
            ring.push_back(cur->vert);
            cur = cur->pair;
            if (!cur) { ok = false; break; }
            cur = cur->next;
        } while (cur && cur != start);
        if (!ok || ring.empty()) continue;

        const double inv_k = 1.0 / (double)ring.size();

        for (int c = 0; c < 3; ++c) {
            // Fairness: x_i - mean(one-ring) = 0
            Geex::LinearConstraint<double> fair;
            fair.Weight = 1.0;
            fair.RHSide = 0.0;
            fair.LHSide.push_back(Geex::VarTerm<double>(vh->var_id, c, 1.0));
            for (Geex::HE_vert* nb : ring) {
                if (nb->var_id > -1)
                    fair.LHSide.push_back(Geex::VarTerm<double>(nb->var_id, c, -inv_k));
                else
                    fair.RHSide += inv_k * nb->pos[c];  // fixed neighbour
            }
            constraints.push_back(fair);

            // Surface-deviation penalization: x_i = x_i_current
            Geex::LinearConstraint<double> dev;
            dev.Weight = deviation_weight;
            dev.RHSide = vh->pos[c];
            dev.LHSide.push_back(Geex::VarTerm<double>(vh->var_id, c, 1.0));
            constraints.push_back(dev);
        }
    }

    Geex::LinearConstraintLeastSQSolver<double> solver(solver_type);
    solver.initialize_system(constraints);

    // Warm start from current positions.
    std::vector<double> solution(map_var_vtx.size() * 3);
    for (size_t i = 0; i < map_var_vtx.size(); ++i)
        for (int c = 0; c < 3; ++c)
            solution[3 * i + c] = map_var_vtx[i]->pos[c];

    if (!solver.solve_system(solution)) {
        Geex::cs_err() << "  Fairing solve failed; keeping unfaired surface" << std::endl;
        return;
    }
    update_vertex_positions(p_mesh, map_var_vtx, solution);
}

// This version automatically finds free vertices according to logic by Hao
bool match_target_curvature(
    Geex::Mesh3D* p_mesh,
    std::vector<double>& v_residuals_out,
    std::vector<double>& m_residual_out,
    std::vector<double>& c_residual_out,
    std::vector<Geex::HE_vert*>& selected_vertices,
    double* avg_movement_out
){

    // Determine free vertices
    std::vector<Geex::HE_vert*> map_var_vtx;
    
    // If selected_vertices is specified, use only those vertices
    // Otherwise, use all vertices that are not on boundary or curve
    if (!selected_vertices.empty()) {
        // Use only the provided selected vertices
        for (Geex::HE_vert* vh : selected_vertices) {
            if (vh) {
                vh->var_id = map_var_vtx.size();
                map_var_vtx.push_back(vh);
            }
        }
        // Mark all other vertices as non-variable
        for (int vitr = 0; vitr < p_mesh->get_num_of_vertices(); vitr++) {
            Geex::HE_vert* vh = p_mesh->get_vertex(vitr);
            if (std::find(selected_vertices.begin(), selected_vertices.end(), vh) == selected_vertices.end()) {
                vh->var_id = -1;
            }
        }
    } else {
        // Original behavior: use all non-boundary, non-curve vertices
        for (int vitr = 0; vitr < p_mesh->get_num_of_vertices(); vitr++)
        {
            Geex::HE_vert* vh = p_mesh->get_vertex(vitr);
            vh->var_id = -1;
            if (!p_mesh->is_on_boundary(vh) && !p_mesh->vertex_on_curve(vh))
            {
                map_var_vtx.push_back(vh);
                vh->var_id = map_var_vtx.size()-1;
            }
        }
    }
    // Aggregate container for individual constraints
    std::vector<Geex::LinearConstraint<double>> all_constraints;

    //compute target normal differences per edge, according to curvature tensor, and build linear system for edge-based curvatures.
    p_mesh->reset_edges_tag(false);
    bool complained = false;
    size_t num_normal_constraints = 0;

    for (int eitr = 0; eitr < p_mesh->get_num_of_edges(); eitr++)
    {
        // Skip if already been tagged true to avoid duplicate work
        Geex::HE_edge* eh = p_mesh->get_edge(eitr);
        if (eh->tag || eh->pair->pair != eh) {
            continue;
        }
        eh->tag = eh->pair->tag = true;

        // Compute normal curvature contributions and directions projected onto principal curvature directions for each face.
        if (eh->face && eh->pair->face && eh->face->tag && eh->pair->face->tag)
        {
            // Sharp feature curves stay unconstrained (free crease); smooth
            // transition curves are treated like interior edges, using the
            // tensor-derived target so the surface stays smooth across them.
            if ((eh->cycle_curv || eh->pair->cycle_curv) && edge_on_sharp_feature_curve(eh))
            {
                continue;
            }

            Eigen::Vector3d curv_norm = compute_target_curvature(eh);

            std::vector<Geex::LinearConstraint<double>> edgecurv_constraints = edgecurv_fillmatrix(eh, curv_norm);
            all_constraints.insert(all_constraints.end(), edgecurv_constraints.begin(), edgecurv_constraints.end());
        }

        // RMF normal boundary conditions on curve edges (paper Sec. 5.3)
        if ( (!eh->face || !eh->pair->face))
        {
            if ((eh->cycle_curv && (eh->cycle_curv->curve_type==Geex::CURVE_TYPE::CRVT_LINE || eh->cycle_curv->curve_type==Geex::CURVE_TYPE::TRIM_LINE)) ||
                (eh->pair->cycle_curv && (eh->pair->cycle_curv->curve_type==Geex::CURVE_TYPE::CRVT_LINE || eh->pair->cycle_curv->curve_type==Geex::CURVE_TYPE::TRIM_LINE)) )
            {
                std::vector<Geex::LinearConstraint<double>> normal_constraints = fillmatrix_normalConstraint(eh, selected_vertices);
                std::vector<Geex::LinearConstraint<double>> normal_constraints_pair = fillmatrix_normalConstraint(eh->pair, selected_vertices);

                num_normal_constraints += normal_constraints.size() + normal_constraints_pair.size();

                all_constraints.insert(all_constraints.end(), normal_constraints.begin(), normal_constraints.end());
                all_constraints.insert(all_constraints.end(), normal_constraints_pair.begin(), normal_constraints_pair.end());
            }
        }
    }
    Geex::cs_log() << "  Normal boundary constraints: " << num_normal_constraints
              << " of " << all_constraints.size() << " total" << std::endl;

    Geex::LinearConstraintLeastSQSolver<double> solver;
    solver.initialize_system(all_constraints);
    // Seed the iterative solver with current vertex positions so near-null
    // (tangential) components stay put instead of collapsing toward the origin.
    std::vector<double> solution(map_var_vtx.size() * 3);
    for (size_t si = 0; si < map_var_vtx.size(); si++)
        for (int sc = 0; sc < 3; sc++)
            solution[3 * si + sc] = map_var_vtx[si]->pos[sc];

    double energy = 0.0;
    bool converged = solver.solve_system(solution);
    if (!converged) {
        return false;  // Solver failed to converge
    }
    double avg_movement = update_vertex_positions(p_mesh, map_var_vtx, solution);
    if (avg_movement_out) *avg_movement_out = avg_movement;

    v_residuals_out = solver.calculate_per_vertex_constraint_residuals(solution, map_var_vtx, all_constraints);
    solver.compute_residual(solution, m_residual_out);
    solver.compute_constraint_residuals(solution, c_residual_out);

    return true;  // Success
}
