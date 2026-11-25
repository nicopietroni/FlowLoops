#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "app_data.h"
#include "cs_logging.h"

#include <igl/hausdorff.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/volume.h>
#include <igl/doublearea.h>
#include <igl/principal_curvature.h>
#include <igl/gaussian_curvature.h>
#include <igl/edge_lengths.h>
#include <igl/massmatrix.h>
#include <igl/random_points_on_mesh.h>
#include <fstream>
#include <queue>
#include <unordered_map>




namespace Geex {

    
    CS_INLINE APP_DataStructs::APP_DataStructs() {
		m_init_mesh_ = NULL;
		m_reference_mesh_ = NULL;
		m_b_bdry_fixed_ = false;

		m_p_curve_net_ = NULL;

		m_int_curve_samples_=0;

		m_p_current_cycle_ = NULL;
		m_p_current_cycle_curve_ = NULL;

		m_p_current_vtx_ = NULL;
		m_vertex_type_ = VertexType::FREE_VTX;

		m_cyclepair_[0] = m_cyclepair_[1] = NULL;

		b_cross_conjugate_switch_ = true;

		m_b_align_bdry_curv_ = false;
    }

    CS_INLINE APP_DataStructs::~APP_DataStructs() {
		if(m_init_mesh_) delete m_init_mesh_;
		if(m_reference_mesh_) delete m_reference_mesh_;
		if(m_p_curve_net_) delete m_p_curve_net_;
    }

    CS_INLINE void APP_DataStructs::save(const std::string& filename) {
        std::ofstream out(filename.c_str()) ;
    }
    
    CS_INLINE void APP_DataStructs::load(const std::string& filename) {
        cs_log() << "loading " << filename << std::endl ;

		if (m_init_mesh_) delete m_init_mesh_;
		m_init_mesh_ = new Mesh3D;
        m_init_mesh_->load_obj(filename.c_str());

		double bbox[6];
		get_bbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);

		cs_log() << "end." << std::endl;
    }

	CS_INLINE void APP_DataStructs::convert_to_objfile(const std::string& filename){
		if (m_p_curve_net_)
		{
			std::ofstream file (filename);

			if (file.is_open())
			{
				for (int cpt_itr = 0; cpt_itr < m_p_curve_net_->m_control_points.size(); cpt_itr++)
				{
					Eigen::Vector4d cpt = m_p_curve_net_->m_control_points[cpt_itr];
					file << "v " << cpt.x() << " " << cpt.y() << " " << cpt.z() << std::endl;
				}

				for (int citr = 0; citr < m_p_curve_net_->m_curves.size(); citr++)
				{
					SplineCurv* curve = m_p_curve_net_->m_curves[citr];

					file << "cstype " << ((curve->type==BSPLINE)?"bspline":"rat bspline") << std::endl;
					file << "deg " << curve->degree << std::endl;
					file << "curv " << curve->para_range[0] << " " << curve->para_range[1];

					for (int cpt_itr = 0; cpt_itr < curve->m_cpoints.size(); cpt_itr++)
					{
						
						file << " " << curve->m_cpoints[cpt_itr] + 1;
					}
					file << std::endl;

					file << "parm u";
					for (int knot_itr = 0; knot_itr < curve->knots.size(); knot_itr++)
					{
						file << " " << curve->knots[knot_itr];
					}

					file << std::endl;

					file << "end" << std::endl;
				}
			}

			file.close();
		}
	}

	CS_INLINE void APP_DataStructs::convert_to_curvefile(const std::string& filename){
		if (m_p_curve_net_)
		{
			std::ofstream file (filename);

			if (file.is_open())
			{
				file << m_p_curve_net_->m_curves.size() << std::endl;

				for (int citr = 0; citr < m_p_curve_net_->m_curves.size(); citr++)
				{
					SplineCurv* curve = m_p_curve_net_->m_curves[citr];

					file << curve->m_cpoints.size() << " 0 2" << std::endl;

					for (int cpt_itr = 0; cpt_itr < curve->m_cpoints.size(); cpt_itr++)
					{
						Eigen::Vector4d cpt = m_p_curve_net_->m_control_points[curve->m_cpoints[cpt_itr]];
						file << cpt.x() << " " << cpt.y() << " " << cpt.z() << std::endl;
					}
				}
			}

			file.close();
		}
	}

	CS_INLINE void APP_DataStructs::save_cycles(){

		if (m_p_curve_net_)
		{
			m_p_curve_net_->save_cycles();
		}

	}

	CS_INLINE void APP_DataStructs::save_curve_types(){

		if (m_p_curve_net_)
		{
			m_p_curve_net_->save_curve_types();
		}

	}

	CS_INLINE void APP_DataStructs::save_curve_types(const std::string& file_prefix){

		if (m_p_curve_net_)
		{
			m_p_curve_net_->save_curve_types(file_prefix);
		}

	}

	CS_INLINE void APP_DataStructs::reverse_component_orientation(Mesh3D* mesh, HE_face* seed_fh){
		if (!mesh || !seed_fh)
		{
			return;
		}

		
		mesh->reset_faces_tag(false);
		mesh->reset_edges_tag(false);
		std::queue<HE_face*> q_faces;
		q_faces.push(seed_fh);
		q_faces.front()->tag = true;

		while(!q_faces.empty())
		{
			HE_face* cur_fh = q_faces.front();
						
			q_faces.pop();

			std::vector<HE_edge*> edges_to_invert;
						
			HE_edge* cur_fedge = cur_fh->edge;
			do
			{
				if (cur_fedge->pair->face && !cur_fedge->pair->face->tag)
				{
					q_faces.push(cur_fedge->pair->face);
					cur_fedge->pair->face->tag = true;
				}

				edges_to_invert.push_back(cur_fedge);

				if (!cur_fedge->pair->face)
				{
					edges_to_invert.push_back(cur_fedge->pair);
				}

				cur_fedge = cur_fedge->next;
			} while (cur_fedge != cur_fh->edge);

			for (int invt_eitr = 0; invt_eitr < edges_to_invert.size(); invt_eitr++)
			{
				HE_edge* eh_invt = edges_to_invert[invt_eitr];
				if (!eh_invt->tag && !eh_invt->pair->tag)
				{
					std::swap(eh_invt->vert, eh_invt->pair->vert);
					eh_invt->tag = eh_invt->pair->tag = true;

					if (eh_invt->vert->edge == eh_invt)
					{
						eh_invt->vert->edge = eh_invt->pair;
					}
					if (eh_invt->pair->vert->edge == eh_invt->pair)
					{
						eh_invt->pair->vert->edge = eh_invt;
					}
				}
							
				std::swap(eh_invt->prev, eh_invt->next);
			}
		}
	}

	CS_INLINE void APP_DataStructs::pairup_g1(PolyGraph_Cycle* cycle_0, PolyGraph_Cycle* cycle_1){
		if (cycle_0 && cycle_1 && !(cycle_0->is_planar ^ cycle_1->is_planar))
		{
			PolyGraph_Cycle* cycle[2] = {cycle_0, cycle_1};

			//find shared edge(s).
			PolyGraph_Edge* cy0_eh = cycle[0]->edge;
			std::vector<std::pair<PolyGraph_Edge*, PolyGraph_Edge*>> shared_eh_pairs;
			do
			{
				std::vector<PolyGraph_Edge*>& adj_edges = cy0_eh->curve_->pg_edges;
				for (int eitr = 0; eitr < adj_edges.size(); eitr++)
				{
					if (adj_edges[eitr]->cycle_ == cycle[1])
					{
						shared_eh_pairs.push_back( std::pair<PolyGraph_Edge*, PolyGraph_Edge*>(cy0_eh, adj_edges[eitr]) );
						break;
					}
				}

				cy0_eh = cy0_eh->next;
			} while (cy0_eh != cycle[0]->edge);

			if (shared_eh_pairs.size() == 0)
			{
				return;
			}

			for (int share_eh_itr = 0; share_eh_itr < shared_eh_pairs.size(); share_eh_itr++)
			{
				PolyGraph_Edge* shared_eh[2] = {shared_eh_pairs[share_eh_itr].first, shared_eh_pairs[share_eh_itr].second};
				//zip up
				std::vector<std::pair<HE_edge*, double>> edge_para_map[2];
				for (int eitr = 0; eitr < m_init_mesh_->get_num_of_edges(); eitr ++)
				{
					HE_edge* eh = m_init_mesh_->get_edge(eitr);
					if (!eh->face && eh->cycle_curv)
					{
						if (eh->cycle_curv == shared_eh[0])
						{
							Eigen::Vector3d midpt = eh->GetMidPoint();
							double para;
							eh->cycle_curv->curve_->project_pt(midpt, &para);

							edge_para_map[0].push_back(std::pair<HE_edge*, double>(eh,para));
						}
						else if(eh->cycle_curv == shared_eh[1]){
							Eigen::Vector3d midpt = eh->GetMidPoint();
							double para;
							eh->cycle_curv->curve_->project_pt(midpt, &para);

							edge_para_map[1].push_back(std::pair<HE_edge*, double>(eh,para));
						}
					}
				}

				struct {
					bool operator() (const std::pair<HE_edge*, double>& lhs, const std::pair<HE_edge*, double>& rhs){
						if (lhs.second < rhs.second)
						{
							return true;
						}
						else{
							return false;
						}
					}
				}edge_para_comp;

				std::sort(edge_para_map[0].begin(), edge_para_map[0].end(), edge_para_comp);
				std::sort(edge_para_map[1].begin(), edge_para_map[1].end(), edge_para_comp);

				//note the orientation of the two patches must match.
				//check orientation, if mismatch, invert one.
				{
					if ( (edge_para_map[0][0].first->next == edge_para_map[0][1].first && 
						edge_para_map[1][0].first->next == edge_para_map[1][1].first) 
						||
						(edge_para_map[0][0].first->prev == edge_para_map[0][1].first && 
						edge_para_map[1][0].first->prev == edge_para_map[1][1].first) 
						)
					{
						//invert the first patch.
						reverse_component_orientation(m_init_mesh_, edge_para_map[0][0].first->pair->face);
					}
				}

				//reset pointers.
				
				for (int eitr = 0; eitr < edge_para_map[0].size(); eitr++)
				{
					HE_edge* eh[2] = {edge_para_map[0][eitr].first, edge_para_map[1][eitr].first};

					std::pair<HE_vert*, HE_vert*> vertices_to_merge [2];
					vertices_to_merge[0] = (std::pair<HE_vert*, HE_vert*>(eh[0]->vert, eh[1]->pair->vert));
					vertices_to_merge[1] = (std::pair<HE_vert*, HE_vert*>(eh[0]->pair->vert, eh[1]->vert));
					for (int vp_itr = 0; vp_itr < 2; vp_itr++)
					{
						std::pair<HE_vert*, HE_vert*>& vtx_pair = vertices_to_merge[vp_itr];

						HE_vert* v0 = vtx_pair.first, *v1 = vtx_pair.second;

						if (!v1->edge || v0 == v1)
						{
							continue;
						}

						HE_edge* v1_eh = v1->edge;
						do
						{
							if (v1_eh->pair->vert == v1)
							{
								v1_eh->pair->vert = v0;
							}
							v1_eh = v1_eh->pair->next;
						} while (v1_eh != v1->edge);

						v1->degree = 0;
						v1->edge = NULL;

						v0->degree = 0;
						HE_edge* v0_eh = v0->edge;
						do
						{
							v0->degree ++;

							v0_eh = v0_eh->pair->next;
						} while (v0_eh != v0->edge);
					
					}

					eh[0]->pair->pair = eh[1]->pair;
					eh[1]->pair->pair = eh[0]->pair;

					eh[0]->prev->next = eh[1]->next;
					eh[0]->next->prev = eh[1]->prev;

					eh[1]->prev->next = eh[0]->next;
					eh[1]->next->prev = eh[0]->prev;

					for (int pair_itr = 0; pair_itr < 2; pair_itr++)
					{
						if (eh[pair_itr]->pair->vert->edge == eh[pair_itr])
						{
							eh[pair_itr]->pair->vert->edge = eh[pair_itr]->pair->next;
						}
					}
				}

				
			}
			

			m_init_mesh_->update_normal();
		}
	}

	CS_INLINE void APP_DataStructs::save_g1_smooth_marking(){
		if (m_p_curve_net_)
		{
			m_p_curve_net_->save_g1_smooth_marking();
		}
	}

	CS_INLINE void APP_DataStructs::apply_g1_smooth_marking(){
		if (m_init_mesh_ && m_p_curve_net_)
		{
			for (std::set<small_set<PolyGraph_Cycle*, 2>>::iterator cy_pair_itr = m_p_curve_net_->m_g1_cyclepairs.begin(); 
				cy_pair_itr != m_p_curve_net_->m_g1_cyclepairs.end(); cy_pair_itr++)
			{
				const small_set<PolyGraph_Cycle*, 2>& cy_pair = (*cy_pair_itr);
				pairup_g1(cy_pair[0], cy_pair[1]);
			}
		}
	}

	CS_INLINE void APP_DataStructs::auto_mark_g1_smoothness(){
		if (m_init_mesh_ && m_p_curve_net_)
		{
			m_p_curve_net_->collect_g1_cycle_pairs();
			apply_g1_smooth_marking();
		}
	}

	CS_INLINE void APP_DataStructs::load_reference_mesh(const std::string& filename){
		if (m_reference_mesh_ != NULL)
		{
			delete m_reference_mesh_;
		}

		m_reference_mesh_ = new Mesh3D;
		m_reference_mesh_->load_obj(filename.c_str());

		m_init_mesh_->compute_boundingbox();
		m_reference_mesh_->compute_boundingbox();

		//match scale by scaling the refmesh.
		Eigen::Vector3d cur_bb_ll (m_init_mesh_->xmin, m_init_mesh_->ymin, m_init_mesh_->zmin);
		Eigen::Vector3d cur_bb_ur (m_init_mesh_->xmax, m_init_mesh_->ymax, m_init_mesh_->zmax);

		Eigen::Vector3d ref_bb_ll (m_reference_mesh_->xmin, m_reference_mesh_->ymin, m_reference_mesh_->zmin);
		Eigen::Vector3d ref_bb_ur (m_reference_mesh_->xmax, m_reference_mesh_->ymax, m_reference_mesh_->zmax);

		Eigen::Vector3d cur_center = 0.5*(cur_bb_ll+cur_bb_ur);
		Eigen::Vector3d ref_center = 0.5*(ref_bb_ll+ref_bb_ur);

		double scale_ratio = (cur_bb_ll-cur_bb_ur).norm()/(ref_bb_ll-ref_bb_ur).norm();

		for (auto vh : *(m_reference_mesh_->get_vertices_list()))
		{
			vh->pos = (vh->pos - ref_center)*scale_ratio + cur_center;
		}

		m_reference_mesh_->update_mesh();

		cs_log() << "Ref mesh loaded."<<std::endl;
	}

		CS_INLINE void APP_DataStructs::load_reference_mesh_mod(MeshData ref_mesh_data){
		// Check for empty or invalid mesh data to prevent crashes
		if (ref_mesh_data.vertices.empty() || ref_mesh_data.faces.empty()) {
			cs_err() << "Warning: Empty reference mesh data provided. Skipping load." << std::endl;
			return;
		}

		if (m_reference_mesh_ != NULL)
		{
			delete m_reference_mesh_;
		}

		m_reference_mesh_ = new Mesh3D;
		m_reference_mesh_->load_from_mesh_data(ref_mesh_data.vertices, ref_mesh_data.faces);

		// Additional check after loading to ensure mesh has valid data
		if (m_reference_mesh_->get_num_of_vertices() == 0) {
			cs_err() << "Warning: Loaded mesh has no vertices. Skipping KD-tree build." << std::endl;
			return;
		}

		m_init_mesh_->compute_boundingbox();
		m_reference_mesh_->compute_boundingbox();

		//match scale by scaling the refmesh.
		Eigen::Vector3d cur_bb_ll (m_init_mesh_->xmin, m_init_mesh_->ymin, m_init_mesh_->zmin);
		Eigen::Vector3d cur_bb_ur (m_init_mesh_->xmax, m_init_mesh_->ymax, m_init_mesh_->zmax);

		Eigen::Vector3d ref_bb_ll (m_reference_mesh_->xmin, m_reference_mesh_->ymin, m_reference_mesh_->zmin);
		Eigen::Vector3d ref_bb_ur (m_reference_mesh_->xmax, m_reference_mesh_->ymax, m_reference_mesh_->zmax);

		Eigen::Vector3d cur_center = 0.5*(cur_bb_ll+cur_bb_ur);
		Eigen::Vector3d ref_center = 0.5*(ref_bb_ll+ref_bb_ur);

		double scale_ratio = (cur_bb_ll-cur_bb_ur).norm()/(ref_bb_ll-ref_bb_ur).norm();

		for (auto vh : *(m_reference_mesh_->get_vertices_list()))
		{
			vh->pos = (vh->pos - ref_center)*scale_ratio + cur_center;
		}

		m_reference_mesh_->update_mesh();

		cs_log() << "Ref mesh loaded."<<std::endl;
	}


	CS_INLINE void APP_DataStructs::load_curv_mesh(const std::string& filename){
		cs_log() << "loading curve network and surface from " << filename << std::endl;

		if(!m_p_curve_net_) 
			m_p_curve_net_ = new CurvNet;

		m_p_curve_net_->clear();


		if (m_init_mesh_)
		{
			delete m_init_mesh_;
		}
		m_init_mesh_ = m_p_curve_net_->load_cycle_discovery_data(filename);
		apply_g1_smooth_marking();
		

		double bbox[6];
		get_bbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);

		cs_log() << "Data bounding box size: " << (bbox[3] - bbox[0]) << " " << (bbox[4] - bbox[1]) << " " << (bbox[5] - bbox[2]) << std::endl;

		cs_log() << "Mesh size: " << m_init_mesh_->get_num_of_vertices() << " vertices, " << m_init_mesh_->get_num_of_faces() << " facets, " <<
			m_init_mesh_->get_num_of_edges()/2 << " edges. " << std::endl;
	}

	CS_INLINE void APP_DataStructs::load_curv_mesh_mod(std::vector<CurveData> curv_data, std::vector<CycleData> cycle_data){

		if(!m_p_curve_net_) 
			m_p_curve_net_ = new CurvNet;

		m_p_curve_net_->clear();


		if (m_init_mesh_)
		{
			delete m_init_mesh_;
		}
		m_init_mesh_ = m_p_curve_net_->load_cycle_discovery_data_from_structures(curv_data, cycle_data);
		apply_g1_smooth_marking();
		

		double bbox[6];
		get_bbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);

		cs_log() << "Data bounding box size: " << (bbox[3] - bbox[0]) << " " << (bbox[4] - bbox[1]) << " " << (bbox[5] - bbox[2]) << std::endl;

		cs_log() << "Mesh size: " << m_init_mesh_->get_num_of_vertices() << " vertices, " << m_init_mesh_->get_num_of_faces() << " facets, " <<
			m_init_mesh_->get_num_of_edges()/2 << " edges. " << std::endl;
	}


	CS_INLINE void APP_DataStructs::load_curvnet(const std::string& filename){
		cs_log() << "loading curve network from "<<filename<<std::endl;
		if(filename.substr(filename.length()-4,4) == ".obj"){
			if(!m_p_curve_net_) 
				m_p_curve_net_ = new CurvNet;

			m_p_curve_net_->clear();

			m_p_curve_net_->load_obj(filename);

			load_polygraph(filename.substr(0,filename.length()-4)+".plyg");
		}
		else if (filename.substr(filename.length()-3,3) == ".ma"){

			if(!m_p_curve_net_) 
				m_p_curve_net_ = new CurvNet;

			m_p_curve_net_->clear();

			if (m_init_mesh_)
			{
				delete m_init_mesh_;
			}
			m_init_mesh_ = m_p_curve_net_->load_ma(filename);

			apply_g1_smooth_marking();
		}
		else{
			cs_err() << "Unknown file format." << std::endl;
		}

		double bbox[6];
		get_bbox(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
	}

	CS_INLINE void APP_DataStructs::load_polygraph(const std::string& filename){

	}

	CS_INLINE void APP_DataStructs::save_polygraph(const std::string& filename){

	}


    CS_INLINE void APP_DataStructs::get_bbox(
        double& x_min, double& y_min, double& z_min,
        double& x_max, double& y_max, double& z_max
	) {
        x_min = y_min = z_min =  1e30 ;
        x_max = y_max = z_max = -1e30 ;

		if (m_init_mesh_ && m_init_mesh_->is_valid())
		{
			x_min = m_init_mesh_->xmin;
			x_max = m_init_mesh_->xmax;
			y_min = m_init_mesh_->ymin;
			y_max = m_init_mesh_->ymax;
			z_min = m_init_mesh_->zmin;
			z_max = m_init_mesh_->zmax;
		}

		if(m_p_curve_net_){
			for(int i=0; i<m_p_curve_net_->m_control_points.size(); i++){
				Eigen::Vector3d pt (m_p_curve_net_->m_control_points[i].x(),m_p_curve_net_->m_control_points[i].y(),m_p_curve_net_->m_control_points[i].z());
				x_min = std::min(pt.x(),x_min);
				x_max = std::max(pt.x(),x_max);
				y_min = std::min(pt.y(),y_min);
				y_max = std::max(pt.y(),y_max);
				z_min = std::min(pt.z(),z_min);
				z_max = std::max(pt.z(),z_max);
			}
		}
    }

	CS_INLINE void APP_DataStructs::savea( const std::string& filename )
	{

		double xmin,ymin,zmin,xmax,ymax,zmax;
		get_bbox(xmin,ymin,zmin,xmax,ymax,zmax);

		//save curve.
		std::string curve_file_name = filename;
		curve_file_name.append("_curve.obj");
		if (m_p_curve_net_)
		{
			std::ofstream ofile (curve_file_name);
			if (ofile.is_open())
			{
				for (auto vtx : m_p_curve_net_->m_control_points)
				{
					ofile << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << std::endl;
				}

				for (auto curve : m_p_curve_net_->m_curves)
				{
					ofile << "l";
					for (auto cpt : curve->m_cpoints)
					{
						ofile << " " <<(cpt+1) ;
					}
					ofile << std::endl;
				}

				int vidx_offset = m_p_curve_net_->m_control_points.size();
				ofile <<"v " << xmin << " " << ymin << " "<< zmin << std::endl;
				ofile <<"v " << xmin << " " << ymax << " "<< zmin << std::endl;
				ofile <<"v " << xmax << " " << ymax << " "<< zmin << std::endl;
				ofile <<"f " << vidx_offset+1 << " " << vidx_offset+2 << " " << vidx_offset+3 <<std::endl;
			}

			ofile.close();
		}

		//save mesh.
		std::string mesh_file_name = filename;
		mesh_file_name.append("_mesh.obj");
		if (m_init_mesh_)
		{
			Mesh3D* pmesh = m_init_mesh_;
			std::ofstream ofile (mesh_file_name);
			if (ofile.is_open())
			{
				ofile << "mtllib mat.mtl" << std::endl;

				for (auto vh : *pmesh->get_vertices_list())
				{
					ofile << "v " << vh->pos << std::endl;
				}

				std::vector<std::vector<HE_face*>> map_cycle_facets (m_p_curve_net_->m_p_polygraph_->get_num_of_cycles(), std::vector<HE_face*>());
				for (auto fh : *pmesh->get_faces_list())
				{
					map_cycle_facets[fh->curve_cycle->id].push_back(fh);
				}

				int cyc_id = 0;
				for (auto& mitr : map_cycle_facets)
				{
					ofile << "usemtl material_"<<(cyc_id++)<<std::endl;
					for (auto fh : mitr)
					{
						ofile << "f";
						HE_edge* eh = fh->edge;
						do
						{
							ofile <<" " <<(eh->pair->vert->id+1);
							eh = eh->next;
						} while (eh != fh->edge);
						ofile << std::endl;
					}
				}
			}

			ofile.close();
		}

		//save reference mesh.
		std::string ref_mesh_file_name = filename;
		ref_mesh_file_name.append("_reference_mesh.obj");
		if (m_reference_mesh_)
		{
			m_reference_mesh_->write_obj(ref_mesh_file_name.c_str());
		}

		cs_log() << "end saving meshes."<<std::endl;
	}

	CS_INLINE void APP_DataStructs::save_mesh(const std::string& filename) {
		std::string mesh_name = filename;
		mesh_name.append(".obj");
		if (m_init_mesh_)
		{
			m_init_mesh_->write_obj(mesh_name.c_str());
			cs_log() << "end saving " << mesh_name << std::endl;
		}
	}


	// Helper function to convert Mesh3D to Eigen matrices for libigl
	CS_INLINE void mesh3d_to_eigen(Mesh3D* mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
		if (!mesh || !mesh->is_valid()) {
			V.resize(0, 3);
			F.resize(0, 3);
			return;
		}
		
		auto* vertices = mesh->get_vertices_list();
		auto* faces = mesh->get_faces_list();
		
		// Convert vertices
		V.resize(vertices->size(), 3);
		for (size_t i = 0; i < vertices->size(); ++i) {
			HE_vert* v = (*vertices)[i];
			V(i, 0) = v->pos.x();
			V(i, 1) = v->pos.y();
			V(i, 2) = v->pos.z();
		}
		
		// Convert faces
		F.resize(faces->size(), 3);
		for (size_t i = 0; i < faces->size(); ++i) {
			HE_face* f = (*faces)[i];
			HE_edge* edge = f->edge;
			
			// Get the three vertices of the triangle
			int v0 = edge->pair->vert->id;
			edge = edge->next;
			int v1 = edge->pair->vert->id;
			edge = edge->next;
			int v2 = edge->pair->vert->id;
			
			F(i, 0) = v0;
			F(i, 1) = v1;
			F(i, 2) = v2;
		}
	}

}