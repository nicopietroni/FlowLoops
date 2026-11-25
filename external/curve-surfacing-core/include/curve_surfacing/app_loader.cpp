#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "app_loader.h"
#include "cs_logging.h"
#include <fstream>
#include <sstream>
#include <string>

#include <fstream>
#include <sstream>
#include <string>

namespace Geex {

	// Helper function to check curve point ordering
	// Prints a warning if curve appears to be in reversed order
	// Returns true if curve appears reversed, false if it appears correct
	CS_INLINE bool ensure_correct_curve_ordering(std::vector<Eigen::Vector3d>& points, int curve_index) {
		if (points.size() < 3) {
			return false; // Not enough points to check ordering
		}
		
		int negative_dot_products = 0;
		int total_checks = 0;
		
		// Check consecutive point pairs
		for (size_t i = 0; i < points.size() - 2; i++) {
			Eigen::Vector3d vec1 = points[i + 1] - points[i];
			Eigen::Vector3d vec2 = points[i + 2] - points[i + 1];
			
			double dot = vec1.dot(vec2);
			total_checks++;
			
			if (dot < 0) {
				negative_dot_products++;
			}
		}
		
		// If more than half of the dot products are negative, curve appears reversed
		if (negative_dot_products > total_checks / 2) {
			cs_log() << "  Warning: Curve " << curve_index << " appears to have reversed point ordering "
			          << "(" << negative_dot_products << "/" << total_checks << " negative dot products)" << std::endl;
			return true; // Curve appears reversed
		}
		
		return false; // Curve appears in correct order
	}

    CS_INLINE Geex::MeshData loadMeshDataFromOBJ(const std::string& filename) {
    Geex::MeshData mesh_data;
    mesh_data.num_vertices = 0;
    mesh_data.num_faces = 0;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        cs_err() << "Error: Could not open OBJ file: " << filename << std::endl;
        return mesh_data;
    }
    
    std::string line;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        
        if (token == "v") {
            // Parse vertex
            double x, y, z;
            if (iss >> x >> y >> z) {
                vertices.push_back(Eigen::Vector3d(x, y, z));
            }
        }
        else if (token == "f") {
            // Parse face
            std::vector<int> face_indices;
            std::string vertex_data;
            
            while (iss >> vertex_data) {
                // Handle vertex/texture/normal format (e.g., "1/1/1" or just "1")
                size_t first_slash = vertex_data.find('/');
                std::string vertex_index_str = vertex_data.substr(0, first_slash);
                
                int vertex_index = std::stoi(vertex_index_str) - 1; // OBJ indices are 1-based
                face_indices.push_back(vertex_index);
            }
            
            if (face_indices.size() >= 3) {
                faces.push_back(face_indices);
            }
        }
    }
    
    file.close();
    
    // Populate MeshData structure
    mesh_data.num_vertices = static_cast<int>(vertices.size());
    mesh_data.num_faces = static_cast<int>(faces.size());
    mesh_data.vertices = vertices;
    mesh_data.faces = faces;
    
    cs_log() << "Loaded OBJ file: " << filename << " with " << mesh_data.num_vertices 
              << " vertices and " << mesh_data.num_faces << " faces" << std::endl;
    
    return mesh_data;
}


	// New modular loading functions
	CS_INLINE std::vector<Geex::CurveData> load_curves_from_file(const std::string& filename) {
		std::vector<Geex::CurveData> curves_data;
		std::ifstream curv_file(filename);
		
		if (!curv_file.is_open()) {
			cs_err() << ">>> Could not load " << filename << std::endl;
			return curves_data;
		}
		
		int curv_num = 0;
		curv_file >> curv_num;
		curves_data.reserve(curv_num);
		
		int reversed_count = 0;
		
		for (int c_itr = 0; c_itr < curv_num; c_itr++) {
			Geex::CurveData curve;
			curv_file >> curve.num_pts >> curve.closed >> curve.capacity;
			curve.points.reserve(curve.num_pts);
			
			for (int pt_itr = 0; pt_itr < curve.num_pts; pt_itr++) {
				double x, y, z;
				curv_file >> x >> y >> z;
				curve.points.push_back(Eigen::Vector3d(x, y, z));
			}
			
			// Check curve point ordering and warn if reversed
			if (ensure_correct_curve_ordering(curve.points, c_itr)) {
				reversed_count++;
			}
			
			curves_data.push_back(curve);
		}
		
		curv_file.close();
		cs_log() << "Loaded " << curv_num << " curves from " << filename;
		if (reversed_count > 0) {
			cs_log() << " (Warning: " << reversed_count << " curves appear to have reversed ordering)";
		}
		cs_log() << std::endl;
		return curves_data;
	}
	
	CS_INLINE std::vector<Geex::CycleData> load_cycles_from_file(const std::string& filename) {
		std::vector<Geex::CycleData> cycles_data;
		
        std::string cycle_filename = filename + ".cycle";
        std::ifstream cycle_file(cycle_filename);

		if (!cycle_file.is_open()) {
			cs_err() << ">>> Could not load cycle file" << std::endl;
			return cycles_data;
		}
		
		int num_cycles = 0;
		cycle_file >> num_cycles;
		cycles_data.reserve(num_cycles);
		
		for (int cy_itr = 0; cy_itr < num_cycles; cy_itr++) {
			Geex::CycleData cycle;
			cycle_file >> cycle.num_curves;
			cycle.curve_info.reserve(cycle.num_curves);
			
			for (int curv_itr = 0; curv_itr < cycle.num_curves; curv_itr++) {
				int curv_id, curv_dir, instance;
				cycle_file >> curv_id >> curv_dir >> instance;
				cycle.curve_info.push_back({curv_id, curv_dir, instance});
			}
			
			cycles_data.push_back(cycle);
		}
		
		cycle_file.close();
		cs_log() << "Loaded " << num_cycles << " cycles from " << cycle_filename << std::endl;
		return cycles_data;
	}
	
	CS_INLINE std::vector<Geex::MeshData> load_patch_meshes_from_files(const std::string& base_filename, int num_patches) {
		std::vector<Geex::MeshData> patches_data;
		patches_data.reserve(num_patches);
		
		for (int cy_itr = 0; cy_itr < num_patches; cy_itr++) {
			std::string mesh_filename = base_filename + std::to_string(cy_itr) + ".off";
			
			std::ifstream mesh_file(mesh_filename);
			
			if (!mesh_file.is_open()) {
				cs_err() << ">>> Failed loading patch " << cy_itr << std::endl;
				patches_data.push_back(Geex::MeshData()); // Add empty patch
				continue;
			}
			
			Geex::MeshData patch;
			
			// Read OFF header
			std::string offmark;
			mesh_file >> offmark;
			
			int unused;
			mesh_file >> patch.num_vertices >> patch.num_faces >> unused;
			
			// Read vertices
			patch.vertices.reserve(patch.num_vertices);
			for (int vitr = 0; vitr < patch.num_vertices; vitr++) {
				double x, y, z;
				mesh_file >> x >> y >> z;
				patch.vertices.push_back(Eigen::Vector3d(x, y, z));
			}
			
			// Read faces
			patch.faces.reserve(patch.num_faces);
			for (int fitr = 0; fitr < patch.num_faces; fitr++) {
				int face_degree;
				mesh_file >> face_degree;
				
				std::vector<int> face_vertices(face_degree);
				for (int fv_itr = 0; fv_itr < face_degree; fv_itr++) {
					mesh_file >> face_vertices[fv_itr];
				}
				
				patch.faces.push_back(face_vertices);
			}
			
			mesh_file.close();
			patches_data.push_back(patch);
			cs_log() << "Loaded patch " << cy_itr << std::endl;
		}
		
		return patches_data;
	}

CS_INLINE FeatureEdgeData load_feature_edges_from_file(const std::string& filename) {
    FeatureEdgeData data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        cs_err() << "Warning: Could not open feature edge file: " << filename << std::endl;
        cs_err() << "Continuing without feature edge constraints." << std::endl;
        return data;
    }
    
    file >> data.num_edges;
    
    if (data.num_edges < 0) {
        cs_err() << "Error: Invalid number of feature edges: " << data.num_edges << std::endl;
        data.num_edges = 0;
        return data;
    }
    
    data.edges.reserve(data.num_edges);
    
    for (int i = 0; i < data.num_edges; ++i) {
        FeatureEdge edge;
        if (!(file >> edge.start.x() >> edge.start.y() >> edge.start.z()
                   >> edge.end.x() >> edge.end.y() >> edge.end.z())) {
            cs_err() << "Error: Failed to read feature edge " << i << std::endl;
            break;
        }
        data.edges.push_back(edge);
    }
    
    file.close();
    
    cs_log() << "Loaded " << data.edges.size() << " feature edges from " << filename << std::endl;
    
    return data;
}

CS_INLINE NormalCurveData load_normal_curves_from_file(const std::string& filename) {
    NormalCurveData data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        cs_err() << "Error: Could not open normal curve file: " << filename << std::endl;
        return data;
    }
    
    // Read number of cycles
    file >> data.num_cycles;
    
    if (data.num_cycles < 0) {
        cs_err() << "Error: Invalid number of cycles: " << data.num_cycles << std::endl;
        data.num_cycles = 0;
        return data;
    }
    
    data.cycles.reserve(data.num_cycles);
    
    // Read each cycle
    for (int cycle_idx = 0; cycle_idx < data.num_cycles; ++cycle_idx) {
        CycleNormalData cycle;
        
        // Read number of edges in this cycle
        if (!(file >> cycle.num_edges)) {
            cs_err() << "Error: Failed to read number of edges for cycle " << cycle_idx << std::endl;
            break;
        }
        
        if (cycle.num_edges < 0) {
            cs_err() << "Error: Invalid number of edges in cycle " << cycle_idx << ": " << cycle.num_edges << std::endl;
            //continue;
            // keep an empty entry: cycles[i] must stay aligned with cycle i
            cycle.num_edges = 0;
            data.cycles.push_back(cycle);
            continue;
        }
        
        cycle.edges.reserve(cycle.num_edges);
        
        // Read each edge in the cycle
        for (int edge_idx = 0; edge_idx < cycle.num_edges; ++edge_idx) {
            EdgeNormalData edge;
            
            // Read: x1 y1 z1 x2 y2 z2 nx ny nz
            if (!(file >> edge.vertex1.x() >> edge.vertex1.y() >> edge.vertex1.z()
                       >> edge.vertex2.x() >> edge.vertex2.y() >> edge.vertex2.z()
                       >> edge.edge_normal.x() >> edge.edge_normal.y() >> edge.edge_normal.z())) {
                cs_err() << "Error: Failed to read edge " << edge_idx << " in cycle " << cycle_idx << std::endl;
                break;
            }
            
            cycle.edges.push_back(edge);
        }
        
        data.cycles.push_back(cycle);
    }
    
    file.close();

    // A truncated/corrupted file stops early: report it and keep num_cycles
    // consistent with what was actually read
    if ((int)data.cycles.size() != data.num_cycles) {
        cs_err() << "Error: normal curve file " << filename << " declares " << data.num_cycles
                 << " cycles but only " << data.cycles.size() << " could be read" << std::endl;
        data.num_cycles = (int)data.cycles.size();
    }
    
    // Report statistics
    int total_edges = 0;
    for (const auto& cycle : data.cycles) {
        total_edges += cycle.edges.size();
    }
    
    cs_log() << "Loaded normal curve file: " << filename << std::endl;
    cs_log() << "  " << data.cycles.size() << " cycles, " << total_edges << " total edges" << std::endl;
    
    return data;
}

}