#ifndef _CURVE_NETWORK_
#define _CURVE_NETWORK_


#include <vector>

#include "bspline.h"

#include "util_functions.h"

// Include for FeatureEdgeData support
#include "app_loader.h"

namespace Geex {

	class CurvNet;

	// Forward declarations — Mesh3D.h is NOT included here to break the
	// circular dependency. Pointers to these types are sufficient.
	class HE_edge;
	class HE_face;
	class Mesh3D;

	class PolyGraph_Cycle;
	class PolyGraph_Edge;
	class PolyGraph;

	// Forward declarations for data structures
	struct CurveData;
	struct CycleData;
	struct MeshData;

	class ParaEdge;


	enum SPLINE_TYPE {BSPLINE, RBSPLINE};
	enum CURVE_TYPE {CRVT_LINE, TRIM_LINE};

	class SplineCurv{
	public:
		int id;
		CurvNet* net;
		enum SPLINE_TYPE type;//spline type: BSPLINE, RBSPLINE
		int degree;
		float para_range[2];
		std::vector<int> m_cpoints;
		std::vector<float> knots;

		// Method to sample points along the spline curve
    	std::vector<Eigen::Vector3d> sample_points(int num_samples);


		//samples inside the curve.
		std::vector<Eigen::Vector3d> m_samples; 
		//parameter values of samples.
		std::vector<double> m_sample_para;
		//interpolating tangents at samples.
		std::vector<Eigen::Vector3d> m_sample_interp_tgt;

		//regenerate num samples from the curve.
		void resample(int num=2);
		//regenerate samples close to the given points.
		void resample(std::vector<Eigen::Vector3d>& pt_list);
		//build tangents at sample points, given two normals at end points.
		void buildInterpTangent(double u0, double u1, Eigen::Vector3d n0, Eigen::Vector3d n1);

		// Method to get the position of a control point !CHANGE!
		Eigen::Vector3d get_control_point_position(int index);

		std::vector<PolyGraph_Edge*> pg_edges; //one of the corresponding mesh edges.
		bool selected;

		SplineCurv(CurvNet* pnet);
		~SplineCurv();
		float* get_gl_ctlarray();

		void UpdateCurve(){
			if(gl_ctlarray_)
				delete gl_ctlarray_;
			gl_ctlarray_=NULL;

			if (arc_)
				delete arc_;
			arc_ = NULL;
		}

		Eigen::Vector3d project_pt(Eigen::Vector3d pt, double* para=NULL);
		Eigen::Vector3d get_value_at(double para);
		void get_drev(double para, Eigen::Vector3d& order1);

		int get_endpt_idx(double para){
			if (abs(para-para_range[0]) < abs(para-para_range[1]))
			{
				return m_cpoints.front();
			}
			else
			{
				return m_cpoints.back();
			}
		}
	protected:
		float* gl_ctlarray_;
		BSplineCurve* arc_;

		void buildCurve();
	};



	class PolyGraph_Cycle{
	public:
		PolyGraph_Edge * edge;
		int degree;
		int id;
		bool selected;

		bool is_planar;

		PolyGraph_Cycle(){
			edge = NULL;
			degree = 0;
			id = -1;
			selected = false;
			is_planar=false;
		}
	};

	class PolyGraph_Edge{
	public:
		int id;

		SplineCurv* curve_;
		bool order_;
		
		PolyGraph_Cycle* cycle_;
		enum CURVE_TYPE curve_type;//curve type: curvature line (CRVT_LINE) or trim line (TRIM_LINE).

		PolyGraph_Edge* prev;
		PolyGraph_Edge* next;

		void compute_endpts(Eigen::Vector3d& start, Eigen::Vector3d& end);

		PolyGraph_Edge(){
			id = -1;
			curve_ = NULL;
			order_ = true;
			cycle_ = NULL;
			curve_type = CURVE_TYPE::CRVT_LINE;
			prev = NULL;
			next = NULL;
		}

		PolyGraph_Edge(SplineCurv* curve_ptr, bool c_order){

			id = -1;

			curve_ = curve_ptr;
			order_ = c_order;

			cycle_ = NULL;
			curve_type = CURVE_TYPE::CRVT_LINE;
			prev = NULL;
			next = NULL;
		}
	};

	class PolyGraph_Vert{
	public:
		Eigen::Vector3d pos;
		PolyGraph_Edge* edge;

		PolyGraph_Vert(){
			pos = Eigen::Vector3d(0,0,0);
			edge = NULL;
		}
		PolyGraph_Vert(Eigen::Vector3d pt){
			pos = pt;
			edge = NULL;
		}
	};

	//The class is datastructure for cycles of network.
	class PolyGraph
	{
	public:
		PolyGraph();
		~PolyGraph();

		PolyGraph_Edge* insert_edge(SplineCurv* curv_ptr, bool order, PolyGraph_Vert* v0 = NULL, PolyGraph_Vert* v1 = NULL){
			PolyGraph_Edge* nedge = new PolyGraph_Edge(curv_ptr, order);
			edges.push_back(nedge);
			nedge->id = edges.size()-1;

			curv_ptr->pg_edges.push_back(nedge);

			return nedge;
		}

		PolyGraph_Cycle* insert_cycle(const std::vector<PolyGraph_Edge*> cedges);

		PolyGraph_Cycle* get_cycle(int c_id) {
			if (c_id > -1 && c_id < cycles.size())
			{
				return cycles[c_id];
			}
			return NULL;
		}

		PolyGraph_Edge* get_edge(int e_id) {
			if(e_id > -1 && e_id < edges.size()) 
				return edges[e_id];
			return NULL;
		}

		// Method to get all edges
		std::vector<PolyGraph_Edge*> get_all_edges() const {
			return edges;
		}

		

		int get_num_of_cycles(){return cycles.size(); }

	private:
		std::vector<PolyGraph_Vert*> vertices;
		std::vector<PolyGraph_Edge*> edges;
		std::vector<PolyGraph_Cycle*> cycles;
	};

	

	class CurvNet{
	public:
		std::vector<Vector4d> m_control_points;
		std::vector<SplineCurv*> m_curves;
		std::set<small_set<PolyGraph_Cycle*, 2>> m_g1_cyclepairs;

		PolyGraph* m_p_polygraph_;

	public:
		CurvNet();
		~CurvNet();

		void clear();

		int load_obj(std::string filename);
		Mesh3D* load_ma(std::string filename);

		Mesh3D* combinePatchMeshes(const std::vector<Mesh3D*>& patch_meshes);
		double calculateAverageBoundaryEdgeLength(Mesh3D* mesh3d);
		void split_edge(Mesh3D* mesh3d, int max_iterations = 10000);
		void split_faces_optimized(Mesh3D* mesh3d, int max_iterations = 10000, int batch_size = 1000);

		//load data from those exported by Ju/Carr's paper.
		Mesh3D* load_cycle_discovery_data(std::string filename);


	// Load cycle discovery data with SurfaceMaker-generated patches
	Mesh3D* load_cycle_discovery_data_from_structures(
		const std::vector<CurveData>& curves_data,
		const std::vector<CycleData>& cycles_data,
		const FeatureEdgeData& feature_edges = FeatureEdgeData());
	
	// Load cycle discovery data with user-supplied patches
	Mesh3D* load_cycle_discovery_data_from_structures(
		const std::vector<CurveData>& curves_data,
		const std::vector<CycleData>& cycles_data,
		const std::vector<MeshData>& patches_data,
		const FeatureEdgeData& feature_edges = FeatureEdgeData());
	
	// Overloaded version that accepts normal curve data for direct RMF population
	Mesh3D* load_cycle_discovery_data_from_structures(
		const std::vector<CurveData>& curves_data,
		const std::vector<CycleData>& cycles_data,
		const FeatureEdgeData& feature_edges,
		const NormalCurveData& normal_data);
	
	// Overloaded version with user-supplied patches AND normal curve data
	Mesh3D* load_cycle_discovery_data_from_structures(
		const std::vector<CurveData>& curves_data,
		const std::vector<CycleData>& cycles_data,
		const std::vector<MeshData>& patches_data,
		const FeatureEdgeData& feature_edges,
		const NormalCurveData& normal_data);

		//save cycle data with marked trimlines using the format of .cycle file.
		void save_cycles();

		//save a per-curve type listing to data_file_name.curve_types, recording
		//which curves are curvature lines (CRVT_LINE) and which are trim lines (TRIM_LINE).
		//Call this after computeRMF / populateRMFFromNormalData so that curve_type is set.
		void save_curve_types();

		//Convenience overload: sets data_file_name to file_prefix then calls save_curve_types().
		void save_curve_types(const std::string& file_prefix);

		//save g1 smoothness marking data to file, data_file_name.g1_marking .
		void save_g1_smooth_marking();

		//collect all pairs of curve cycles that have shared boundary of g1 smoothness.
		void collect_g1_cycle_pairs();

		//resample the curves, and smooth the curves.
		void smooth_curve_network();

		//merge close endpts in the curve network, to make sure the close end points are numerically identical.
		void merge_close_points_in_curve_network();

		// Method to get the position of a control point
    	Eigen::Vector3d get_control_point(int index) const {
        // Implement this method to return the actual position of the control point
        // For example, if control points are stored in a vector of Eigen::Vector3d:
        return m_control_points[index].head<3>();

    }

		void match_patchboundary_curves(Mesh3D* patches);

		void check_cycle_planarity(PolyGraph_Cycle* cycle);

		// Helper function to export a mesh patch to OFF format

		//This function computes RMF vectors on cycle curves, and marks if the cycle curves are curvature lines.
		//It assumes the cycles and the surface meshes are initialized and associated by match_patchboundary_curves().
		void computeCurveRMF(Mesh3D* patches);

		void computeCornerNormals(PolyGraph_Edge* pg_edge, Eigen::Vector3d end_normals[2], double end_para[2], const std::vector<ParaEdge>* mesh_edges = NULL);
	// Populate RMF frames from normal curve data
	void populateRMFFromNormalData(Mesh3D* patches, const NormalCurveData& normal_data);
		// Fix RMF frames that point inward by analyzing cycle geometry
		void averageOverlappingEdgeFrames(Mesh3D* mesh, double tolerance);
		void fixInwardRMFFrames(Mesh3D* patches);
		void averageOverlappingRMFFrames(Mesh3D* mesh, const FeatureEdgeData& feature_edges = FeatureEdgeData());
		bool areEdgesOverlapping(HE_edge* edge1, HE_edge* edge2);
		void averageRMFFramesInGroup(const std::vector<HE_edge*>& edges);
		Eigen::Vector3d computeFaceNormal(HE_face* face);
		bool shouldFlipRMFFrame(PolyGraph_Edge* pg_edge, HE_edge* mesh_edge, bool is_clockwise);
		bool isCycleClockwise(PolyGraph_Cycle* cycle);
		std::vector<HE_edge*> findAssociatedMeshEdges(Mesh3D* patches, PolyGraph_Edge* pg_edge);
		Eigen::Vector3d computeExpectedOutwardDirection(PolyGraph_Edge* pg_edge, HE_edge* mesh_edge, bool is_clockwise);
		Eigen::Vector3d computeCycleCenterDirection(PolyGraph_Edge* pg_edge, double param);

	protected:


		//The function computes the RMF normals along given sample points (v_pts), given the boundary conditions that the curve 
		// have normals (n) at parameter values (u). The RMF normals are assigned to the corresponding mesh edges (and their pairing half edges).
		//The function also determines if the PolyGraph_Edge is a curvature line or trim line, based on the consistency of RMF normal vectors.
		void computeRMFonCurve(PolyGraph_Edge* pg_edge, SplineCurv* curve, double  u [2], Eigen::Vector3d  n [2], std::vector<ParaEdge>& v_pts);

		std::string data_file_name;
	};


	class Para_Point{
	public:
		double para;
		Eigen::Vector3d point;
		Para_Point(double u, Eigen::Vector3d pt){
			para = u;
			point = pt;
		}

		bool operator< (const Para_Point& rhs) const{
			return para < rhs.para;
		}
	};

	class ParaEdge{
	public:
		double para;
		HE_edge* edge;
		Eigen::Vector3d mid_pt;
		ParaEdge(){
			para=0;
			edge=NULL;
			mid_pt.x() = mid_pt.y() = mid_pt.z() = 0.0;
		}
		ParaEdge(double p, HE_edge* eh, Eigen::Vector3d mpt){
			para=p;
			edge=eh;
			mid_pt=mpt;
		}

		bool operator<(const ParaEdge& rhs) const{
			return para < rhs.para;
		}
	};

}


// ── Directional-style header-only include ────────────────────────────────
#ifndef CURVE_SURFACING_STATIC
#  include "CurvNet.cpp"
#endif

#endif