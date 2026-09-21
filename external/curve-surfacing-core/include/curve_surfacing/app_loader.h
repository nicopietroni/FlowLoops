#ifndef APP_LOADER_H
#define APP_LOADER_H

#include <queue>
#include <unordered_set>
#include <chrono>
#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <Eigen/Dense>

namespace Geex {

// Data structures for loading cycle discovery data
struct CurveData {
    int num_pts;
    int closed;
    int capacity;
    std::vector<Eigen::Vector3d> points;
};

struct CycleData {
    int num_curves;
    std::vector<std::array<int, 3>> curve_info; // [curve_id, direction, instance]
};

struct MeshData {
    int num_vertices;
    int num_faces;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces; // Each face is a vector of vertex indices
    std::vector<int> face_cycle_ids; // Cycle ID for each face 
};

// Feature edge data structures for sharp edge preservation
struct FeatureEdge {
    Eigen::Vector3d start;  // First endpoint (x0, y0, z0)
    Eigen::Vector3d end;    // Second endpoint (x1, y1, z1)
    
    // Helper: compute edge direction
    Eigen::Vector3d direction() const {
        return (end - start).normalized();
    }
    
    // Helper: compute edge length
    double length() const {
        return (end - start).norm();
    }
    
    // Helper: compute midpoint
    Eigen::Vector3d midpoint() const {
        return (start + end) / 2.0;
    }
};

struct FeatureEdgeData {
    int num_edges = 0;
    std::vector<FeatureEdge> edges;
    
    // Check if data is empty
    bool empty() const {
        return edges.empty();
    }
};

// Data structures for normal curve files (.normalcurve)
// These store edge normals for computing RMF (Rotation Minimizing Frames)
struct EdgeNormalData {
    Eigen::Vector3d vertex1;     // First vertex of the edge
    Eigen::Vector3d vertex2;     // Second vertex of the edge
    Eigen::Vector3d edge_normal; // Normal vector associated with the edge
    
    // Helper: compute edge direction
    Eigen::Vector3d direction() const {
        return (vertex2 - vertex1).normalized();
    }
    
    // Helper: compute edge length
    double length() const {
        return (vertex2 - vertex1).norm();
    }
    
    // Helper: compute edge midpoint
    Eigen::Vector3d midpoint() const {
        return (vertex1 + vertex2) / 2.0;
    }
};

struct CycleNormalData {
    int num_edges = 0;
    std::vector<EdgeNormalData> edges;
    
    // Check if data is empty
    bool empty() const {
        return edges.empty();
    }
    
    // Get total number of vertices in the cycle
    int num_vertices() const {
        return num_edges; // For a closed cycle, num_vertices == num_edges
    }
};

struct NormalCurveData {
    int num_cycles = 0;
    std::vector<CycleNormalData> cycles;
    
    // Check if data is empty
    bool empty() const {
        return cycles.empty();
    }
};

// New modular loading functions
std::vector<CurveData> load_curves_from_file(const std::string& filename);
std::vector<CycleData> load_cycles_from_file(const std::string& filename);
MeshData loadMeshDataFromOBJ(const std::string& filename);
std::vector<MeshData> load_patch_meshes_from_files(const std::string& base_filename, int num_patches);
FeatureEdgeData load_feature_edges_from_file(const std::string& filename);
NormalCurveData load_normal_curves_from_file(const std::string& filename);

} // namespace Geex

#ifndef CURVE_SURFACING_STATIC
#  include "app_loader.cpp"
#endif

#endif // APP_LOADER_H