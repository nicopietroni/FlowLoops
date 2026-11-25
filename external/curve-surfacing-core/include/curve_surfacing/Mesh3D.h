#ifndef MESH3D_H
#define MESH3D_H

#include <vector>
#include <string>
#include <cstring>
#include <cmath> 
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <queue>
#include <cstdlib>
#include <atomic>

#include <Eigen/Dense>


inline std::atomic<int> g_he_edge_alive{0};
inline std::atomic<int> g_he_face_alive{0};
inline std::atomic<int> g_he_vert_alive{0};


/** \defgroup MeshCore Mesh data structure */

namespace Geex
{

	class HE_vert;
	class HE_face;
	class HE_edge;
	class Mesh3D;

	class SplineCurv;
	class PolyGraph_Edge;
	class PolyGraph_Cycle;

	enum PatchType {REGULAR_PATCH, SPHERICAL_PATCH, UNDEFINED_PATCH} ;

	enum VertexType {FREE_VTX, FIXED_VTX, APPROX_BOUNDARY_VTX};

	class BD_Corner{
	public:
		HE_vert* p_corner_pos;
		HE_edge* p_corner_drv;
		int step;
		BD_Corner(){
			p_corner_pos = NULL;
			p_corner_drv = NULL;
			step = 0;
		}
		BD_Corner(HE_vert* cptr, HE_edge* eptr, int u);
	};

	class BD_Seg{
	public:
		BD_Corner two_ends[2];

		HE_edge* p_bdseg_drv;
		HE_vert* p_bdseg_pos;
		int step;

	private:
		int counter;
	public:
		BD_Seg(){
			counter = 0;
			p_bdseg_drv = NULL;
			p_bdseg_pos = NULL;
			step = 0;
		}

		BD_Seg(HE_vert* bd_v, HE_edge* bd_dir){
			counter = 0;
			p_bdseg_pos = bd_v;
			p_bdseg_drv = bd_dir;
			step = 0;
		}

		void push_corner(BD_Corner& cn){
			assert(counter<2);
			two_ends[counter++] = cn;
		}

		bool intersect_at(BD_Seg& other, bool isFront){
			if (other.two_ends[0].p_corner_pos == two_ends[isFront?0:1].p_corner_pos)
			{
				return true;
			}
			else if (other.two_ends[1].p_corner_pos == two_ends[isFront?0:1].p_corner_pos)
			{
				std::swap(other.two_ends[0], other.two_ends[1]);
				return true;
			}
			return false;
		}
	};

	/*!
	*	The basic vertex class for half-edge structure.
	*/

	class HE_vert
	{
	public:
		Eigen::Vector3d pos;		//!< 3D coordinate
		HE_edge* edge;			//!< one of the half-edges_list emanating from the vertex
		Eigen::Vector3d normal;	    //!< vertex normal
		ptrdiff_t id;					//!< index 
		unsigned int degree;			//!< the degree of vertex 
		bool tag;				//!< tag for programming easily

		Eigen::Vector3d pdir_min;		//principal direction, the minor axis.
		Eigen::Vector3d pdir_max;		//principal direction, the major axis.

		int singular_idx;

		bool  frame_fixed;// !< indicate whether the frame of cyclide patches are fixed , bob

		bool fixed;	

		bool fairing;

		bool selected;

		int var_id; //store variable id, hpan.

		double quality; //holder for various quality measures.


		std::vector<HE_edge*> pg_edges; //used for polygraph construction only.


		double	m_fCurvePara; //parameter of corresponding point on curve. hpan
		SplineCurv*	m_pParaCurve; //pointer to curve with corresponding point. hpan
		int	m_dOnCurveParaID; //id of parameter variable, used for indexing the array of variables. hpan

		enum VertexType m_vtype; //type of the vertex. hpan

		std::vector<BD_Seg> hermite_data;	//to keep the boundary Hermite data. hpan

	public:
		HE_vert(const Eigen::Vector3d &v)
			:pos(v), edge(0), id(-1), degree(0), tag(false), frame_fixed(false), fixed(false), selected(false), fairing(true), pdir_min(Eigen::Vector3d(0,0,0)),
			pdir_max(Eigen::Vector3d(0,0,0)), m_fCurvePara(0), m_pParaCurve(NULL), m_dOnCurveParaID(-1), m_vtype(VertexType::FREE_VTX), singular_idx(0), var_id(-1)
		{
			++g_he_vert_alive;  // Track instance creation
			quality =  0.0;


		}
		~HE_vert()
		{
			--g_he_vert_alive;  // Track instance destruction
		}
	};

	/*!
	*	The basic edge class for half-edge structure.
	*/

	class HE_edge
	{ 
	public:
		HE_vert* vert;	//!< vertex at the end of the half-edge
		HE_edge* pair;	//!< oppositely oriented adjacent half-edge 
		HE_face* face;	//!< face the half-edge borders
		HE_edge* next;	//!< next half-edge around the face
		HE_edge* prev;	//!< prev half-edge around the face
		ptrdiff_t id;			//!< index
		bool tag;				//!< tag for programming easily

		bool selected;

		std::vector<HE_vert*> sample_vh; //array of sample vertices. for merging subdivided polygraph patches.

		int pij; //period jump of cross directions.
		double kij; //chain curvature along dual edge.

		int var_id; //variable id.
		
		Eigen::Vector3d ortho_dual; //the orthogonal dual. hpan

		double target_sin_hdi;
		double dual_len;

		PolyGraph_Edge* cycle_curv;

		bool fixed;
		Eigen::Vector3d tangent;

	public:

		HE_edge()
			:vert(0), pair(0), face(0), next(0), prev(0), id(-1), tag(false), selected(false), fixed(false), tangent(Eigen::Vector3d(0,0,0)),
			pij(0), kij(0), var_id(-1), target_sin_hdi(0), dual_len(0.0)
		{
			++g_he_edge_alive;  // Track instance creation
			ortho_dual = Eigen::Vector3d(0,0,0);

			cycle_curv = NULL;


		}
		~HE_edge()
		{
			--g_he_edge_alive;  // Track instance destruction
		}

		inline Eigen::Vector3d GetMidPoint()
		{
			return 0.5*(vert->pos+pair->vert->pos);
		}
	};


	/*!
	*	The basic face class for half-edge structure.
	*/

	class HE_face
	{ 
	public:
		HE_edge* edge;		//!< one of the half-edges_list bordering the face	
		unsigned int valence;		//!< the number of edges_list 
		Eigen::Vector3d normal;	//!< face normal
		int id;					    //!< index
		bool tag;				    //!< tag for programming easily
		std::vector<ptrdiff_t>  texture_indices; //! texture indices
		std::vector<ptrdiff_t>  normal_indices; //! texture indices
		bool fixed;

		bool selected;

		double theta[4]; //cross field angle. use two for conjugate field, four for 4-RoSy field. hpan.
		int var_id[2]; //var id. hpan.

		Eigen::Vector3d tensor_axis [2]; //axis of estimated curvature tensor. i.e. the principal directions.
		double tensor_diag [2]; //the estimated principal curvature values.

		Eigen::Vector3d unrestricted_tensor [2]; //principal curvature directions.
		double unrestricted_pcv [2]; //principal curvature values.

		PolyGraph_Cycle* curve_cycle; //curve cycle this facet belongs to.

		double quality; //holder for quality measurements.
		bool is_umbilic;


	private:


	public:
		HE_face()
			:edge(0), id(-1), tag(false), fixed(false), selected(false)
		{
			++g_he_face_alive;  // Track instance creation
			var_id[0] = var_id[1] = -1;
			theta[0] = 0.0;
			theta[1] = M_PI*0.5;

			tensor_axis[0] = tensor_axis[1] = Eigen::Vector3d(0,0,0);
			tensor_diag[0] = tensor_diag[1] = 0;

			unrestricted_tensor[0] = unrestricted_tensor[1] = Eigen::Vector3d(0,0,0);

			curve_cycle = NULL;

			quality = 0.0;
			is_umbilic = false;

			unrestricted_pcv[0] = unrestricted_pcv[1] = 0.0;


		}
		~HE_face()
		{
			--g_he_face_alive;  // Track instance destruction
		}


		std::vector<HE_vert*> get_vertices() const {
			std::vector<HE_vert*> vertices;
			HE_edge* start_edge = edge;
			HE_edge* current_edge = edge;
			do {
				vertices.push_back(current_edge->vert);
				current_edge = current_edge->next;
			} while (current_edge != start_edge);
			return vertices;
		}

		inline Eigen::Vector3d GetCentroid()
		{
			Eigen::Vector3d V(0, 0, 0);
			HE_edge* he = edge;
			int i = 0;
			do 
			{
				V+= he->vert->pos;
				he = he->next;
				i++;
			} while(he!=edge);
			return V/double(i);
		}

		double get_area(){

			 Eigen::Vector3d x[3];
			 HE_edge* cur_edge = edge;
			 for (int i = 0; i < 3; i ++)
			 {
				 x[i] = cur_edge->vert->pos;
				 cur_edge = cur_edge->next;
			 }

			 return 0.5*(x[2]-x[0]).cross(x[1]-x[0]).norm();
		}

		bool is_degenerate(){
			 Eigen::Vector3d x[3];
			 HE_edge* cur_edge = edge;
			 for (int i = 0; i < 3; i ++)
			 {
				 x[i] = cur_edge->vert->pos;
				 cur_edge = cur_edge->next;
			 }

			 double edge_len = ((x[1]-x[0]).norm() + (x[2]-x[1]).norm() + (x[0]-x[2]).norm())/3.0;

			 double ratio = (x[2]-x[0]).cross(x[1]-x[0]).norm()/(edge_len*edge_len);
			 // Degenerate-triangle test: area/edge² below sin(5°) — i.e. the
			 // triangle's smallest angle is under ~5 degrees.
			 const double kDegenerateAngleDeg = 5.0;
			 const double threshold = std::sin(kDegenerateAngleDeg/180.0*M_PI);
			 if (ratio < threshold)
			 {
				 return true;
			 }

			 return false;
		}

		void get_cross_axis(Eigen::Vector3d& u_axis, Eigen::Vector3d& v_axis){
			Eigen::Vector3d faxis = (edge->vert->pos - edge->prev->vert->pos).normalized();
			Eigen::Vector3d binormal = normal.cross(faxis);

			u_axis = Eigen::Vector3d(0,0,0);
			v_axis = Eigen::Vector3d(0,0,0);

			double cos_theta = std::cos(theta[0]), sin_theta = std::sin(theta[0]);

			u_axis = cos_theta*faxis + sin_theta*binormal;
			v_axis = -sin_theta*faxis + cos_theta*binormal;
		}

		int index(const HE_vert* vh){
			int idx = 0;
			HE_edge* eh = edge;
			do 
			{
				if (eh->vert == vh)
				{
					return idx;
				}
				idx++;
				eh = eh->next;
			} while (eh != edge);
			return -1;
		}
		bool has_texture_map()
		{
			return (!texture_indices.empty()) && (texture_indices.size() == (size_t)valence);
		}
		bool has_normal_map()
		{
			return (!normal_indices.empty()) && (normal_indices.size() == (size_t)valence);
		}

		void copy_regular_tri_data(Mesh3D* p_mesh, HE_face* rhs);

		HE_edge* has_edge(int v0, int v1) {
			HE_edge* eh = edge;
			do
			{
				if (eh->prev->vert->id == v0 && eh->vert->id == v1)
				{
					return eh;
				}
				eh = eh->next;
			} while (eh != edge);

			return NULL;
		}

	};

	
	bool inline CompareEdgeID (HE_edge* he1, HE_edge* he2 )
	{
		return he1->id < he2->id;
	}

	bool inline CompareVertexID (HE_vert* hv1, HE_vert* hv2 )
	{
		return hv1->id < hv2->id;
	}

	bool inline CompareFaceID (HE_face* hf1, HE_face* hf2 )
	{
		return hf1->id < hf2->id;
	}


	/*!
	* a half-edge based mesh data structure
	* For understanding half-edge structure, 
	* please read the article in http://www.flipcode.com/articles/article_halfedge.shtml
	*/

	class Mesh3D
	{

	public:

		typedef std::vector<HE_vert* > VERTEX_LIST;
		typedef std::vector<HE_face* > FACE_LIST;
		typedef std::vector<HE_edge* > EDGE_LIST;

		typedef VERTEX_LIST* PTR_VERTEX_LIST;
		typedef FACE_LIST* PTR_FACE_LIST;
		typedef EDGE_LIST* PTR_EDGE_LIST;

		typedef  VERTEX_LIST::iterator VERTEX_ITER;
		typedef  FACE_LIST::iterator FACE_ITER;
		typedef  EDGE_LIST::iterator EDGE_ITER;

		typedef  VERTEX_LIST::reverse_iterator VERTEX_RITER;
		typedef  FACE_LIST::reverse_iterator FACE_RITER;
		typedef  EDGE_LIST::reverse_iterator EDGE_RITER;
		typedef std::pair<HE_vert*, HE_vert* > PAIR_VERTEX;

	protected:


		PTR_VERTEX_LIST vertices_list;		//!< store vertices
		PTR_EDGE_LIST edges_list;			//!< store edges
		PTR_FACE_LIST faces_list;			//!< store faces


		bool m_closed;						//!< indicate whether the mesh is closed
		bool m_quad;						//!< indicate whether the mesh is quadrilateral
		bool m_tri;							//!< indicate whether the mesh is triangular
		bool m_hex;							//!< indicate whether the mesh is hexagonal
		bool m_pentagon;					//!< indicate whether the mesh is pentagonal

		std::map<PAIR_VERTEX, HE_edge*> m_edgemap;	


		int m_num_components;				//!< number of components
		int m_num_boundaries;				//!< number of boundaries
		int m_genus;						//!< the genus value

		bool m_encounter_non_manifold;


	public:

		double xmax, xmin, ymax, ymin, zmax, zmin;

		std::vector<std::vector<HE_vert*> > boundaryvertices;

	public:
		Mesh3D(void);

		~Mesh3D(void);

		Mesh3D(const Mesh3D&) = delete;
		Mesh3D& operator=(const Mesh3D&) = delete;

		inline PTR_VERTEX_LIST get_vertices_list()
		{
			return vertices_list;
		}

		inline PTR_EDGE_LIST get_edges_list()
		{
			return edges_list;
		}

		inline PTR_FACE_LIST get_faces_list()
		{
			return faces_list;
		}

		inline ptrdiff_t get_num_of_vertices()
		{	
			return vertices_list?(ptrdiff_t)vertices_list->size():0;
		}

		inline ptrdiff_t get_num_of_faces()
		{
			return faces_list?(ptrdiff_t)faces_list->size():0;
		}

		inline ptrdiff_t get_num_of_edges()
		{
			return edges_list?(ptrdiff_t)edges_list->size():0;
		}

		inline HE_vert* get_vertex(ptrdiff_t id) 
		{
			return id >= get_num_of_vertices()||id<0? NULL:(*vertices_list)[id];
		}

		inline HE_edge* get_edge(ptrdiff_t id) 
		{
			return id >= get_num_of_edges()||id<0? NULL: (*edges_list)[id];
		}

		inline HE_face* get_face(ptrdiff_t id) 
		{
			return id >= get_num_of_faces()||id<0? NULL: (*faces_list)[id];
		}

		inline int get_num_of_components()
		{
			return m_num_components;
		}

		inline int get_num_of_boundaries()
		{
			return m_num_boundaries;
		}

		inline int genus()
		{
			return m_genus;
		}

		inline bool is_valid() 
		{
			if( get_num_of_vertices()==0 || get_num_of_faces() == 0 )
				return false;
			return true;
		}

		inline bool is_closed() 
		{
			return m_closed;
		}

		inline bool is_tri()
		{
			return m_tri;
		}

		inline bool is_quad()
		{
			return m_quad;
		}

		inline bool is_hex()
		{
			return m_hex;
		}

		inline bool is_pentagon()
		{
			return m_pentagon;
		}

		/*!
		*	\param v a 3d point
		*	\return a pointer to the created vertex
		*/
		HE_vert* insert_vertex(Eigen::Vector3d &v);

		/*!
		*	\param vec_hv the vertices list of a face
		*	\param texture the pointer of texture vector
		*	\param normal the normal of texture vector
		*	\return a pointer to the created face
		*/
		HE_face* insert_face(VERTEX_LIST& vec_hv, std::vector<ptrdiff_t> * texture = 0, std::vector<ptrdiff_t> * normal = 0);

		bool is_on_boundary(HE_vert* hv); 
		bool is_on_boundary(HE_face* hf); 
		bool is_on_boundary(HE_edge* he); 


		bool load_off(const char* fins);
		void write_off(const char* fouts);
		bool load_obj(const char* fins);
		bool load_from_mesh_data(const std::vector<Eigen::Vector3d>& vertices, 
		                         const std::vector<std::vector<int>>& faces);
		void write_obj(const char* fouts);
		void write_vtk(const char* fouts);

		/*! 
		*	call it when you have created the mesh
		*/
		void update_mesh(); 

		/*!
		*	compute all the normals of vertices and faces
		*/
		void update_normal(bool onlyupdate_facenormal = false);

		void compute_boundingbox();

		Mesh3D* make_copy();

		Mesh3D* make_copy_simple();

		Mesh3D* reverse_orientation();

		/*!
		*	for a pair of edges, only one of them is tagged to be true.
		*/
		void init_edge_tag();

		void reset_vertices_tag(bool tag_status);
		void reset_faces_tag(bool tag_status);
		void reset_edges_tag(bool tag_status);
		void reset_all_tag(bool tag_status);

		void translate(const Eigen::Vector3d& tran_V);
		void scale(double factorx, double factory, double factorz);

		inline bool is_encounter_nonmanifold()
		{
			return m_encounter_non_manifold;
		}

		void make_delaunay( bool& topo_updated);

		bool delaunay_test( HE_edge* triedge );

		void swap_edge( HE_edge* triedge );

		void split_face(HE_face* face, HE_vert* vtx);

		void split_edge(HE_edge* edge, HE_vert* vtx);

		/*!
		*	\param iterations number of smoothing iterations
		*	\param damping_factor smoothing strength (0.0 = no smoothing, 1.0 = full smoothing)
		*	\param preserve_boundary if true, boundary vertices are not moved
		*/
		void laplacian_smooth(int iterations = 1, double damping_factor = 0.5, bool preserve_boundary = true);

		/**
		* @brief Laplacian smoothing using linear least-squares constraint solver
		* 
		* This version formulates Laplacian smoothing as a linear constraint problem:
		* For each free vertex v_i, create constraint: v_i = (1/|N_i|) * sum(v_j for j in N_i)
		* where N_i is the set of neighboring vertices.
		* 
		* @param lambda smoothing weight (higher = more smoothing towards Laplacian)
		* @param preserve_boundary if true, boundary vertices are fixed as constraints
		*/
		void laplacian_smooth_constraints(double lambda = 1.0, bool preserve_boundary = true);


		HE_edge* pg_insert_edge(HE_vert* vstart, HE_vert* vend);
		HE_face* pg_insert_face(std::vector<HE_edge*>& edges);

		PolyGraph_Edge* vertex_on_curve(HE_vert* vh);
		int vertex_curves(HE_vert* vh);

		Mesh3D* make_copy_regular();

		double edge_chain_crvt(HE_edge* eij);

		double bc_dual_edge_length(HE_edge* edge);

		double ortho_dual_edge_data(HE_edge* edge, Eigen::Vector3d& dir, Eigen::Vector3d& ndiff, bool use_cot_w = false) const;

		void get_cvt_fg( double& f, std::vector<double>& g );

		Eigen::Vector3d get_cvt_centroid(HE_vert* vh);

		void init_vef_list();


	public:
		// Delete every vertex, edge and face and drop the lists; the next
		// insert recreates them. What remesh_patch uses to rebuild in place.
		void clear_data();

	private:

		HE_edge* insert_edge(HE_vert* vstart, HE_vert* vend); 

		void clear_vertices();
		void clear_edges();
		void clear_faces();

		void check_closed();
		void check_meshtype();

		void compute_faces_list_normal(); 
		void compute_perface_normal(HE_face* hf);
		void compute_vertices_list_normal();
		void compute_pervertex_normal(HE_vert* hv); 

		void compute_num_components();
		void compute_num_boundaries();
		void compute_genus();

		void set_nextedge_for_border_vertices(); 

		void remove_hanged_vertices();

		/*!
		set mesh's edge(vertex(startid), vertex(endid)->id = edgeid.
		only used in make_copy
		*/
		void copy_edge_id(ptrdiff_t edgeid, ptrdiff_t startid, ptrdiff_t endid, Mesh3D* mesh);


		static const int facet_sample_num_ = 30;
	};


} //end of namespace


#ifndef CURVE_SURFACING_STATIC
#  include "Mesh3D.cpp"
#endif

#endif //MESH3D_H