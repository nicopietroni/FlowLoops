#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes (full types needed by method bodies) ──
#include "CurvNet.h"
#include "Mesh3D.h"
#include "bspline.h"
#include "SurfaceMaker.h"
#include "app_loader.h"
#include "cs_logging.h"

#include <sstream>
#include <iostream>
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>
#include <eigen_interface.h>
#include <box3.h>
#include <mesh_functions.h>
#include <triangular_remesh.h>
#include <mesh_clean.h>


namespace Geex {

	// Target-edge-length multiplier for per-patch remeshing of the initial
	// triangulation. 1.0 keeps the raw cycle-discovery density; 1.5 coarsens
	// slightly, trading detail for triangle quality ("balanced mesh quality").
	const double kRemeshEdgeLenMultiplier = 1.5;

	// Curve-network fairing solve (control-point Laplacian smoothing).
	const int    kFairingCGMaxIterations = 1000;
	const double kFairingCGTolerance     = 1e-6;


	CS_INLINE SplineCurv::SplineCurv(CurvNet* pnet){
		id = -1;
		net=pnet;
		type=SPLINE_TYPE::BSPLINE;
		degree=0;
		gl_ctlarray_ = NULL;
		
		selected = false;
		arc_ = NULL;
		
	}

	CS_INLINE SplineCurv::~SplineCurv(){
		if(gl_ctlarray_) delete gl_ctlarray_;
		if (arc_)	delete arc_;
	}

	CS_INLINE void SplineCurv::buildCurve(){
		if (!arc_)
		{
			arc_ = new BSplineCurve;
			arc_->degree=degree;
			arc_->dense=arc_->DEFAULT_DENSE;
			arc_->knot=new double[knots.size()];
			for(int i=0;i<knots.size();i++)
				arc_->knot[i]=knots[i];
			arc_->l=knots.size()-1-2*degree;
			arc_->n_points=m_cpoints.size();
			arc_->weights=new double[arc_->n_points];
			arc_->Xcoeff =new double[arc_->n_points];
			arc_->Ycoeff =new double[arc_->n_points];
			arc_->Zcoeff =new double[arc_->n_points];
		
			if(type == SPLINE_TYPE::RBSPLINE){
				for(int i=0; i<m_cpoints.size(); i++){
					Vector4d cpt = net->m_control_points[m_cpoints[i]];
					arc_->Xcoeff[i] = cpt.x()*cpt.w();
					arc_->Ycoeff[i] = cpt.y()*cpt.w();
					arc_->Zcoeff[i] = cpt.z()*cpt.w();
					arc_->weights[i] = cpt.w();
				}
			}
			else if (type == SPLINE_TYPE::BSPLINE)
			{
				for(int i=0; i<m_cpoints.size(); i++){
					Vector4d cpt = net->m_control_points[m_cpoints[i]];
					arc_->Xcoeff[i] = cpt.x();
					arc_->Ycoeff[i] = cpt.y();
					arc_->Zcoeff[i] = cpt.z();
					arc_->weights[i] = 1;
				}
			}
		}
	}

	CS_INLINE Eigen::Vector3d SplineCurv::project_pt(Eigen::Vector3d pt, double* para){

		buildCurve();

		double ftx,fty,ftz;
		double rpara = arc_->PointProjection(pt.x(),pt.y(),pt.z(),&ftx,&fty,&ftz);

		if (para)
		{
			*para = rpara;
		}

		return Eigen::Vector3d(ftx,fty,ftz);
	}

	CS_INLINE void SplineCurv::get_drev(double para, Eigen::Vector3d& order1){
		
		buildCurve();

		BSplineCurve *first_order=arc_->DerivativeOrderOne();

		double px,py,pz,pw;

		px = py = pz = 0;
		pw = 1;

		if (first_order->degree >= 1)
		{
			first_order->GetValueAt(para,&px,&py,&pz,&pw);
		}
		else {
			arc_->PolylineTangent(para,&px,&py,&pz);
		}
		
		if (type==SPLINE_TYPE::BSPLINE)
		{
			order1.x() = px; order1.y() = py; order1.z() = pz;
		}
		else if (type==SPLINE_TYPE::RBSPLINE)
		{
			double cx, cy, cz, cw;
			arc_->GetValueAt(para, &cx, &cy, &cz, &cw);

			order1.x() = (px*cw - cx*pw)/(cw*cw);
			order1.y() = (py*cw - cy*pw)/(cw*cw);
			order1.z() = (pz*cw - cz*pw)/(cw*cw);
		}

		delete first_order;
	}

	CS_INLINE Eigen::Vector3d SplineCurv::get_value_at(double para){
		buildCurve();

		double px,py,pz,pw;
		arc_->GetValueAt(para,&px,&py,&pz,&pw);

		return Eigen::Vector3d(px/pw,py/pw,pz/pw);
	}

	CS_INLINE std::vector<Eigen::Vector3d> SplineCurv::sample_points(int num_samples){
		std::vector<Eigen::Vector3d> sampled_points;
        double step = 1.0 / (num_samples - 1);
        for (int i = 0; i < num_samples; ++i) {
            double t = i * step;
            sampled_points.push_back(get_value_at(t));
        }
        return sampled_points;
	}

	CS_INLINE void SplineCurv::resample(int num){

		m_samples.clear();
		m_sample_para.clear();

		//equi-curve-dist sampling

		std::vector<Eigen::Vector3d> all_cpts;
		for (int citr = 0; citr < m_cpoints.size(); citr++)
		{
			Vector4d cpt = net->m_control_points[m_cpoints[citr]];
			all_cpts.push_back(Eigen::Vector3d(cpt.x(), cpt.y(), cpt.z()));
		}

		double totalDist=0;
		for (int citr = 0; citr < m_cpoints.size()-1; citr++)
		{
			totalDist += (all_cpts[citr]-all_cpts[citr+1]).norm();
		}

		double seg_len = totalDist/(double)(num-1);

		double seg_para_len = (para_range[1]-para_range[0])/(double)(num-1);

		m_samples.clear();
		m_sample_para.clear();

		m_samples.push_back(all_cpts[0]);
		m_sample_para.push_back(para_range[0]);

		int left_cidx=0;
		double acc_len=0.0;
		for (int i = 0; i < num-2; i++)
		{

			while (acc_len<seg_len)
			{
				double cur_seg_len = (all_cpts[left_cidx+1] - all_cpts[left_cidx]).norm();
				acc_len += cur_seg_len;
				left_cidx++;
			}
			left_cidx--;
			double cur_seg_len = (all_cpts[left_cidx+1] - all_cpts[left_cidx]).norm();
			double para_t = (acc_len-seg_len)/cur_seg_len;
			Eigen::Vector3d npt = para_t*all_cpts[left_cidx] + (1-para_t)*all_cpts[left_cidx+1];

			m_samples.push_back(npt);
			m_sample_para.push_back(para_range[0]+(i+1)*seg_para_len);

			acc_len -= seg_len;
			left_cidx++;
		}
		m_samples.push_back(all_cpts.back());
		m_sample_para.push_back(para_range[1]);
	}

	CS_INLINE void SplineCurv::resample(std::vector<Eigen::Vector3d>& pt_list){
		buildCurve();

		m_samples.clear();
		m_sample_para.clear();
		double px,py,pz,pw;

		arc_->GetValueAt(knots[0], &px, &py, &pz, &pw);
		m_samples.push_back(Eigen::Vector3d(px/pw,py/pw,pz/pw));
		m_sample_para.push_back(knots[0]);

		std::vector<Para_Point> uplist;
		for (int i = 0; i < pt_list.size(); i++)
		{
			Eigen::Vector3d pt = pt_list[i];
			double fpx, fpy, fpz, para;
			para = arc_->PointProjection(pt.x(), pt.y(), pt.z(), &fpx, &fpy, &fpz);
			uplist.push_back(Para_Point(para, Eigen::Vector3d(fpx,fpy,fpz)));
		}

		std::sort(uplist.begin(), uplist.end());

		for (int i = 0; i < uplist.size(); i++)
		{
			m_samples.push_back(uplist[i].point);
			m_sample_para.push_back(uplist[i].para);
		}
		
		arc_->GetValueAt(knots.back(), &px, &py, &pz, &pw);
		m_samples.push_back(Eigen::Vector3d(px/pw,py/pw,pz/pw));
		m_sample_para.push_back(knots.back());
	}

	CS_INLINE void SplineCurv::buildInterpTangent(double u0, double u1, Eigen::Vector3d n0, Eigen::Vector3d n1){
		buildCurve();

		m_sample_interp_tgt.clear();

		if (u0>u1)
		{
			std::swap(u0,u1);
			std::swap(n0,n1);
		}

		//linear interpolation.
		for (int i = 0; i < m_sample_para.size(); i++)
		{
			double upara = m_sample_para[i];
			double t = (upara-u0)/(u1-u0);
			Eigen::Vector3d lerped_n = (n0*(1-t) + n1*t).normalized();
			Eigen::Vector3d drev1;
			this->get_drev(upara, drev1);
			drev1 = drev1.normalized();
			m_sample_interp_tgt.push_back(drev1.cross(lerped_n).normalized());
		}

		//reorient coherently.
		Eigen::Vector3d ref_dir = m_sample_interp_tgt[0];
		for (int i = 1; i < m_sample_interp_tgt.size(); i++)
		{
			if (ref_dir.dot(m_sample_interp_tgt[i]) < 0)
			{
				m_sample_interp_tgt[i] *= -1;
			}

			ref_dir = m_sample_interp_tgt[i];
		}
	}

	// Method to get the position of a control point
	CS_INLINE Eigen::Vector3d SplineCurv::get_control_point_position(int index) {
        return net->get_control_point(index);
    }

	CS_INLINE float* SplineCurv::get_gl_ctlarray(){
		if(!gl_ctlarray_){
			if(type==SPLINE_TYPE::BSPLINE){
				gl_ctlarray_ = new float[m_cpoints.size()*3];

				for(int j=0; j<m_cpoints.size(); j++){
					int cpidx = m_cpoints[j];
					Vector4d cpt = net->m_control_points[cpidx];
					gl_ctlarray_[3*j  ] = cpt.x();
					gl_ctlarray_[3*j+1] = cpt.y();
					gl_ctlarray_[3*j+2] = cpt.z();
				}
			}
			else if(type==SPLINE_TYPE::RBSPLINE){
				gl_ctlarray_ = new float[m_cpoints.size()*4];

				for(int j=0; j<m_cpoints.size(); j++){
					int cpidx = m_cpoints[j];
					Vector4d cpt = net->m_control_points[cpidx];
					gl_ctlarray_[4*j  ] = cpt.x()*cpt.w();
					gl_ctlarray_[4*j+1] = cpt.y()*cpt.w();
					gl_ctlarray_[4*j+2] = cpt.z()*cpt.w();
					gl_ctlarray_[4*j+3] = cpt.w();
				}
			}
		}

		return gl_ctlarray_;
	}

	CS_INLINE Eigen::Vector3d IntersectionPoint(SplineCurv* curve0, SplineCurv* curve1){
		Eigen::Vector3d endpts0[2] = {curve0->get_value_at(curve0->para_range[0]), curve0->get_value_at(curve0->para_range[1])};
		Eigen::Vector3d endpts1[2] = {curve1->get_value_at(curve1->para_range[0]), curve1->get_value_at(curve1->para_range[1])};

		double min_dist = 1e10;
		Eigen::Vector3d min_dist_pt (0,0,0);
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if ((endpts0[i]-endpts1[j]).norm() < min_dist)
				{
					min_dist = (endpts0[i]-endpts1[j]).norm();
					min_dist_pt = 0.5*(endpts0[i]+endpts1[j]);
				}
			}
		}

		return min_dist_pt;
	}

	CS_INLINE void PolyGraph_Edge::compute_endpts(Eigen::Vector3d& start, Eigen::Vector3d& end){
		SplineCurv* curve = curve_;
		
		SplineCurv* adj_curve[2];
		adj_curve[0] = prev->curve_;
		adj_curve[1] = next->curve_;

		start = IntersectionPoint(curve, adj_curve[0]);
		end = IntersectionPoint(curve, adj_curve[1]);
	}

	CS_INLINE CurvNet::CurvNet(){
		m_p_polygraph_ = NULL;
	}

	CS_INLINE CurvNet::~CurvNet(){	
		for(int i=0; i<m_curves.size(); i++)
			delete m_curves[i];

		if (m_p_polygraph_)
		{
			delete m_p_polygraph_;
		}
	}

	CS_INLINE void CurvNet::clear(){
		m_control_points.clear();
		for(int i=0; i<m_curves.size(); i++)
			delete m_curves[i];
		m_curves.clear();

		if (m_p_polygraph_)
		{
			delete m_p_polygraph_;
			m_p_polygraph_ = NULL;
		}
	}

	CS_INLINE int CurvNet::load_obj(std::string filename){

		FILE *m_pFile = fopen(filename.c_str(), "r");

		if(!m_pFile){
			cs_err()<<"Error: failed opening file "<<filename<<std::endl;
			return 0;
		}

		char *tok;

#define LINE_BUF_SIZE 1024
		{
			char pLine[LINE_BUF_SIZE];
			while(fgets(pLine, LINE_BUF_SIZE, m_pFile))
			{
				//read vertices
				if(pLine[0] == 'v' && pLine[1] == ' ')
				{
					std::string vcoord (pLine);
					vcoord.pop_back();
					while (vcoord.back() == '\\')
					{
						vcoord.back() = ' ';
						fgets(pLine, LINE_BUF_SIZE, m_pFile);
						vcoord += std::string(pLine);
						vcoord.pop_back();
					}

					std::istringstream iss_vcoord (vcoord);
					std::string vlabel;
					iss_vcoord >> vlabel;

					Vector4d nvv(0,0,0,1.0);
					
					for (int i=0; i<4 && iss_vcoord; i++) 
					{
						std::string coord;
						iss_vcoord >> coord;
						nvv[i] = atof(coord.c_str());
					}

					m_control_points.push_back(nvv);
				}
				else if(pLine[0] == 'v' && pLine[1] == 't')
				{
				}
				else if(pLine[0] == 'v' && pLine[1] == 'n')
				{
				}
				//read curves
				else if(std::string(pLine,pLine+7) == "cstype "){

					SplineCurv* spln_curv = new SplineCurv(this);

					if(std::string(pLine+7,pLine+14) == "bspline"){
						spln_curv->type = SPLINE_TYPE::BSPLINE;
					}
					else if(std::string(pLine+7,pLine+18) == "rat bspline"){
						spln_curv->type = SPLINE_TYPE::RBSPLINE;
					}

					fgets(pLine, LINE_BUF_SIZE, m_pFile);
					spln_curv->degree = atoi(pLine+3);

					fgets(pLine, LINE_BUF_SIZE, m_pFile);
					std::string cpts (pLine);
					cpts.pop_back();
					while (cpts.back() == '\\')
					{
						cpts.back() = ' ';
						fgets(pLine, LINE_BUF_SIZE, m_pFile);
						cpts += std::string(pLine);
						cpts.pop_back();
					}

					std::istringstream iss_cv (cpts);
					std::string	iss_cv_sub;
					iss_cv >> iss_cv_sub;
					iss_cv >> iss_cv_sub;
					spln_curv->para_range[0] = atof(iss_cv_sub.c_str());
					iss_cv >> iss_cv_sub;
					spln_curv->para_range[1] = atof(iss_cv_sub.c_str());
					iss_cv >> iss_cv_sub;

					do
					{
						int vtx_id = atoi(iss_cv_sub.c_str());
						spln_curv->m_cpoints.push_back(vtx_id-1);
						iss_cv >> iss_cv_sub;
					} while (iss_cv);

					fgets(pLine, LINE_BUF_SIZE, m_pFile);
					std::string param (pLine);
					param.pop_back();
					while(param.back() == '\\'){
						param.back() = ' ';
						fgets(pLine, LINE_BUF_SIZE, m_pFile);
						param += std::string(pLine);
						param.pop_back();
					}

					std::istringstream iss (param);
					std::string sub;
					iss >> sub;
					iss >> sub;
					iss >> sub;
					do
					{
						spln_curv->knots.push_back(atof(sub.c_str()));	
						iss >> sub;
					} while (iss);

					m_curves.push_back(spln_curv);

				}
			}

		}

		fclose(m_pFile);

		//assign id to curves.
		for (int i = 0; i < m_curves.size(); i++)
		{
			m_curves[i]->id = i;
		}

		return 1;
	}

	CS_INLINE void load_cycle_patch_ma(std::string base_name, Mesh3D* nmesh, int cyc_count, PolyGraph_Cycle* cyc_ptr){
		//load cycle patch here.

		std::stringstream mesh_filename ;
		mesh_filename << base_name << "_o" << (cyc_count) << ".obj";
	
		std::ifstream mesh_file (mesh_filename.str());

		if (mesh_file.is_open())
		{
						
			const unsigned int vid_offset = nmesh->get_num_of_vertices();

			while (!mesh_file.eof())
			{
				char line_head;
				mesh_file >> line_head;

				if (line_head == 'v')
				{
					double px, py, pz;
					mesh_file >> px >> py >> pz;
					Eigen::Vector3d temp_vertex(px,py,pz);
					nmesh->insert_vertex(temp_vertex);
				}
				else if (line_head == 'f')
				{
					int vid[4];
					mesh_file >> vid[0] >> vid[1] >> vid[2] >> vid[3];

					for (int tri_itr = 0; tri_itr < 2; tri_itr++)
					{
						std::vector<HE_vert*> face_vts (3, NULL);

						for (int fv_itr = 0; fv_itr < 3; fv_itr++)
						{
							face_vts[fv_itr] = nmesh->get_vertex(vid[(fv_itr + tri_itr*2)%4] - 1 +  vid_offset);
						}

						HE_face* fh = nmesh->insert_face(face_vts);
						if (!fh)
						{
							cs_err() << "Failed inserting facet" << std::endl;
						}
						else
						{
							fh->curve_cycle = cyc_ptr;
						}
					}

				}
			}
		}
		else
		{
			cs_err() << ">>> failed loading patch "<<cyc_count<<std::endl;
		}

		mesh_file.close();
	}


	CS_INLINE Mesh3D* CurvNet::load_ma(std::string filename){

		std::map<std::string, SplineCurv*> map_name_curve;

		std::vector<std::vector<SplineCurv*>> list_cycles;

		data_file_name = filename.substr(0,filename.length()-3);

		std::ifstream file (filename);

		if (!file.is_open())
		{
			cs_err() << "Failed opening " << filename <<". No such file."<< std::endl;

			return NULL;
		}

#undef LINE_BUF_SIZE
#define LINE_BUF_SIZE 512

		char line[LINE_BUF_SIZE];
		file.getline(line, LINE_BUF_SIZE);

		while (line[0] != '\0')
		{
			std::istringstream iss (line);
			std::string command;
			iss >> command;

			if (command == "createNode")
			{
				std::string nodetype;
				iss >> nodetype;
				if (nodetype.find("nurbsCurve")!=std::string::npos)
				{
					//get curve parent name.
					std::string curve_shape_name_opt, curve_shape_name, curve_name_opt, curve_name;
					iss >> curve_shape_name_opt;
					iss >> curve_shape_name;
					iss >> curve_name_opt;
					iss >> curve_name;

					if (curve_name_opt.length() ==0 && curve_name.length()==0)
					{
						//read next line.
						file.getline(line, LINE_BUF_SIZE);
						iss = std::istringstream (line);
						iss >> curve_name_opt >> curve_name;
					}

					curve_name.erase(curve_name.begin());
					curve_name.erase(curve_name.end()-1);
					curve_name.erase(curve_name.end()-1);

					//look forward the next few lines, and search for attributes.
					do
					{
						file.getline(line, LINE_BUF_SIZE);
						std::istringstream inn_iss (line);
						std::string attribute, att_name;
						inn_iss >> attribute >> att_name;

						if (attribute == "setAttr" && att_name == "\".cc\"")
						{
							file.getline(line, LINE_BUF_SIZE);
							std::istringstream l1 (line);
							//first line: degree, number of spans, form (0=open, 1=closed, 2=periodic), rational (yes/no), dimension
							SplineCurv* nspln_crv = new SplineCurv(this);
							l1 >> nspln_crv->degree;
							int num_span;
							l1 >> num_span;
							int form;
							l1 >> form;
							std::string rational;
							l1 >> rational;
							nspln_crv->type = (rational=="yes")?SPLINE_TYPE::RBSPLINE : SPLINE_TYPE::BSPLINE;

							file.getline(line, LINE_BUF_SIZE);
							std::istringstream l2 (line);
							//second line: number of knots, list of knot values
							int num_knots;
							l2 >> num_knots;
							for (int knot_itr = 0; knot_itr < num_knots; knot_itr++)
							{
								double knot_val;
								if(l2 >> knot_val)
									nspln_crv->knots.push_back(knot_val);
								else{
									file.getline(line, LINE_BUF_SIZE);
									l2 = std::istringstream(line);
									knot_itr --;
								}

								if (knot_itr == 0 || knot_itr == num_knots-1)
								{
									nspln_crv->knots.push_back(knot_val);
								}
							}

							nspln_crv->para_range[0] = nspln_crv->knots.front();
							nspln_crv->para_range[1] = nspln_crv->knots.back();

							//third line: number of CVs.
							file.getline(line, LINE_BUF_SIZE);
							int num_cvs = atoi(line);
							while(num_cvs == 0){
								file.getline(line, LINE_BUF_SIZE);
								num_cvs = atoi(line);
							}

							//fourth and later lines: CV positions in x y z w(if rational)
							for (int cv_itr = 0; cv_itr < num_cvs; cv_itr++)
							{
								file.getline(line, LINE_BUF_SIZE);
								std::istringstream llater (line);

								double x,y,z,w=0;
								llater >> x >> y >> z;
								if (nspln_crv->type == SPLINE_TYPE::RBSPLINE)
								{
									llater >> w;
								}
								this->m_control_points.push_back(Vector4d(x,y,z,w));

								nspln_crv->m_cpoints.push_back(this->m_control_points.size()-1);
							}

							this->m_curves.push_back(nspln_crv);

							map_name_curve[curve_name] = nspln_crv;

							break;
						}
					} while (true);
					

				}
				else if(nodetype.find("transform")!=std::string::npos){
					//get transform type.
					std::string node_name_opt, node_name;
					iss >> node_name_opt >> node_name;
					if (node_name.find("polySurface") != std::string::npos) //this is tailored to read data from QuadDesign paper.
					{
						//skip the addAttr line.
						file.getline(line, LINE_BUF_SIZE);

						//get the number of curves in this loop.
						file.getline(line, LINE_BUF_SIZE);
						std::istringstream set_att_curvenum_iss (line);

						std::string set_att_com, set_att_opt;
						int cycle_curve_num;
						set_att_curvenum_iss >> set_att_com >> set_att_opt >> cycle_curve_num;

						std::vector<SplineCurv*> cycle_curves;

						do
						{
							file.getline(line, LINE_BUF_SIZE);
							std::istringstream inn_iss (line);
							std::string attribute, att_name, type_opt, type_val, cyc_curve_name;
							//setAttr ".curveLoop[0]" -type "string" "curve100";
							inn_iss >> attribute >> att_name >> type_opt >> type_val >> cyc_curve_name;

							if (attribute == "setAttr"){
								
								cs_log() << cyc_curve_name.substr(1,cyc_curve_name.size()-3) << std::endl;

								cycle_curves.push_back(map_name_curve[cyc_curve_name.substr(1,cyc_curve_name.size()-3)]);
							}
							else
							{
								break;
							}
						}while(true);

						list_cycles.push_back(cycle_curves);
					}
				}
			}
			
			file.getline(line, LINE_BUF_SIZE);
			
		}

		file.close();

		//assign id to curves.
		for (int i = 0; i < m_curves.size(); i++)
		{
			m_curves[i]->id = i;
		}

		//build polygraph.

		Mesh3D* nmesh = new Mesh3D;

		if (!m_p_polygraph_)
		{
			m_p_polygraph_ = new PolyGraph;
		}
		else
		{
			delete m_p_polygraph_;
			m_p_polygraph_ = new PolyGraph;
		}

		//check if .cycle file exists.
		std::string cycle_filename = data_file_name + "_marked.cycle";
		std::ifstream cycle_file (cycle_filename);
		if (cycle_file.is_open())
		{
			int num_cycles = 0;
			cycle_file >> num_cycles;
			for (int cy_itr = 0; cy_itr < num_cycles; cy_itr++)
			{
				int ccurv_num;
				cycle_file >> ccurv_num;

				std::vector<PolyGraph_Edge*> cycle_edges;

				for (int curv_itr = 0; curv_itr < ccurv_num; curv_itr++)
				{
					int curv_id, curv_dir, instance;
					cycle_file >> curv_id >> curv_dir >> instance;

					cycle_edges.push_back(m_p_polygraph_->insert_edge(m_curves[curv_id], (curv_dir==1)));

					{
						cycle_edges.back()->curve_type = (CURVE_TYPE)instance;
					}
				}

				PolyGraph_Cycle* cur_cycle = m_p_polygraph_->insert_cycle(cycle_edges);
				

				//load cycle patch here.
				load_cycle_patch_ma(data_file_name, nmesh, cy_itr, cur_cycle);
			}

			cs_log() << ".cycle file loaded." << std::endl;

			nmesh->update_mesh();

			match_patchboundary_curves(nmesh);
		}
		else if (list_cycles.size() > 0)
		{

			for (int cy_itr = 0; cy_itr < list_cycles.size(); cy_itr++)
			{
				std::vector<SplineCurv*>& cycle_curves = list_cycles[cy_itr];
				std::vector<PolyGraph_Edge*> cycle_edges;

				for (int curv_itr = 0; curv_itr < cycle_curves.size(); curv_itr++)
				{
					cycle_edges.push_back(m_p_polygraph_->insert_edge(cycle_curves[curv_itr],true));
				}

				PolyGraph_Cycle* cur_cycle = m_p_polygraph_->insert_cycle(cycle_edges);
				
				//load cycle patch here.
				load_cycle_patch_ma(data_file_name, nmesh, cy_itr, cur_cycle);
				
			}

			nmesh->update_mesh();
			cs_log() << "cycle data loaded from ma file." << std::endl;

			match_patchboundary_curves(nmesh);
		}

		//load g1 smoothness markings.
		std::string g1smooth_filename = data_file_name + ".g1_marking";

		std::ifstream g1smooth_file (g1smooth_filename);

		if (g1smooth_file.is_open()){
			int cy_pair_num;
			g1smooth_file >> cy_pair_num;

			m_g1_cyclepairs.clear();

			for (int cy_pair_itr = 0; cy_pair_itr < cy_pair_num; cy_pair_itr++)
			{
				int cy0, cy1;
				g1smooth_file >> cy0 >> cy1;

				small_set<PolyGraph_Cycle*, 2> cypair;
				cypair.insert(m_p_polygraph_->get_cycle(cy0));
				cypair.insert(m_p_polygraph_->get_cycle(cy1));

				m_g1_cyclepairs.insert(cypair);
			}

			g1smooth_file.close();
			cs_log() << "g1 smoothness marking file loaded." << std::endl;
		}

		nmesh->update_normal();

		return nmesh;
	}

	CS_INLINE void CurvNet::save_g1_smooth_marking(){
		std::string filename = data_file_name + ".g1_marking";

		cs_log() << "saving g1 markings to " << filename << " ... ";

		std::ofstream ofile (filename);

		if (ofile.is_open())
		{
			ofile << m_g1_cyclepairs.size() << std::endl;
			for (std::set<small_set<PolyGraph_Cycle*, 2>>::iterator cy_pair_itr = m_g1_cyclepairs.begin(); cy_pair_itr != m_g1_cyclepairs.end(); cy_pair_itr++)
			{
				ofile << (*cy_pair_itr)[0]->id << " " << (*cy_pair_itr)[1]->id << std::endl;
			}

			ofile.close();
		}

		cs_log() << "end." << std::endl;
	}

	CS_INLINE void CurvNet::merge_close_points_in_curve_network(){
	}

	CS_INLINE void CurvNet::smooth_curve_network(){
		//resample the curves.
		std::vector<Vector4d> new_control_pts;
		for (int citr = 0; citr < m_curves.size(); citr++)
		{
			SplineCurv* curve = m_curves[citr];
			curve->resample(curve->m_cpoints.size());

			curve->m_cpoints.clear();
			curve->knots.clear();
			curve->knots.push_back(curve->para_range[0]);
			for (int sitr = 0; sitr < curve->m_samples.size(); sitr++)
			{
				Eigen::Vector3d sample = curve->m_samples[sitr];
				new_control_pts.push_back(Vector4d(sample.x(), sample.y(), sample.z(), 1.0));
				curve->knots.push_back(curve->m_sample_para[sitr]);
				curve->m_cpoints.push_back(new_control_pts.size()-1);
			}
			curve->knots.push_back(curve->para_range[1]);

			curve->UpdateCurve();
		}
		m_control_points = new_control_pts;

		//find coincident points and tangent lines.
		class PtBin{
		public:
			Eigen::Vector3d pt;
			std::vector<std::pair<SplineCurv*,double>> v_curve_para;
			PtBin(Eigen::Vector3d point) : pt(point){}
		};
		class TgtBin{
		public:
			Eigen::Vector3d pt;
			Eigen::Vector3d tgt;
			std::vector<std::pair<SplineCurv*,double>> v_curve_para;
			TgtBin(Eigen::Vector3d point, Eigen::Vector3d tangent) : pt(point), tgt (tangent) {}
		};

		const double endpt_closeness_threshold = 1e-4;

		std::vector<PtBin> pt_bins;
		for (int crv_itr_0 = 0; crv_itr_0 < m_curves.size(); crv_itr_0++)
		{
			double endparas_0[2] = {m_curves[crv_itr_0]->para_range[0], m_curves[crv_itr_0]->para_range[1]};
			Eigen::Vector3d endpts_0[2] = {m_curves[crv_itr_0]->get_value_at(endparas_0[0]), m_curves[crv_itr_0]->get_value_at(endparas_0[1])};
			
			for (int crv_itr_1 = crv_itr_0+1; crv_itr_1 < m_curves.size(); crv_itr_1++)
			{
				double endparas_1[2] = {m_curves[crv_itr_1]->para_range[0], m_curves[crv_itr_1]->para_range[1]};
				Eigen::Vector3d endpts_1[2] = {m_curves[crv_itr_1]->get_value_at(endparas_1[0]), m_curves[crv_itr_1]->get_value_at(endparas_1[1])};

				for (int pitr0 = 0; pitr0 < 2; pitr0++)
				{
					for (int pitr1 = 0; pitr1 < 2; pitr1++)
					{
						if ((endpts_0[pitr0]-endpts_1[pitr1]).norm() < endpt_closeness_threshold)
						{
							PtBin npt_bin (endpts_0[pitr0]);
							npt_bin.v_curve_para.push_back(std::pair<SplineCurv*, double>(m_curves[crv_itr_0], endparas_0[pitr0]));
							npt_bin.v_curve_para.push_back(std::pair<SplineCurv*, double>(m_curves[crv_itr_1], endparas_1[pitr1]));

							pt_bins.push_back(npt_bin);
							break;
						}
					}
				}
			}
		}


		std::vector<std::vector<bool>> curve_endpt_added (m_curves.size(), std::vector<bool>(2,false));
		
		for (auto& ptbin_itr : pt_bins)
		{
			for (int crv_itr_int = 0; crv_itr_int < m_curves.size(); crv_itr_int ++)
			{
				SplineCurv* crv_itr = m_curves[crv_itr_int];
				double endparas[2] = {crv_itr->para_range[0], crv_itr->para_range[1]};
				Eigen::Vector3d endpts[2] = {crv_itr->get_value_at(endparas[0]), crv_itr->get_value_at(endparas[1])};
				for (auto eptr = 0; eptr < 2; eptr++)
				{
					if ( !curve_endpt_added[crv_itr_int][eptr] && (endpts[eptr]-ptbin_itr.pt).norm() < endpt_closeness_threshold)
					{
						curve_endpt_added[crv_itr_int][eptr] = true;

						if (std::find_if(ptbin_itr.v_curve_para.begin(), ptbin_itr.v_curve_para.end(), 
							[&](const std::pair<SplineCurv*,double>& a) {return a.first == crv_itr;}) == ptbin_itr.v_curve_para.end()
							)
						{
							ptbin_itr.v_curve_para.push_back(std::pair<SplineCurv*,double>(crv_itr, endparas[eptr]));
							break;
						}
					}
				}
			}
		}

		const double tgt_parallel_threshold = std::sin(40.0/180.0*M_PI);

		std::vector<TgtBin> tgt_bins;
		for (auto& ptbin_itr : pt_bins)
		{
			if (ptbin_itr.v_curve_para.size()<3)
			{
				continue;
			}
			for (int crv0 =0; crv0 < ptbin_itr.v_curve_para.size(); crv0 ++)
			{
				Eigen::Vector3d tgt0;
				ptbin_itr.v_curve_para[crv0].first->get_drev(ptbin_itr.v_curve_para[crv0].second, tgt0);
				for (int crv1 = crv0+1; crv1 < ptbin_itr.v_curve_para.size(); crv1++)
				{
					Eigen::Vector3d tgt1;
					ptbin_itr.v_curve_para[crv1].first->get_drev(ptbin_itr.v_curve_para[crv1].second, tgt1);

					if ((tgt0.normalized().cross(tgt1.normalized())).norm() < tgt_parallel_threshold)
					{
						SplineCurv* curve0 = ptbin_itr.v_curve_para[crv0].first;
						SplineCurv* curve1 = ptbin_itr.v_curve_para[crv1].first;

						int endpt_idx0 = curve0->get_endpt_idx(ptbin_itr.v_curve_para[crv0].second);
						int endpt_idx1 = curve1->get_endpt_idx(ptbin_itr.v_curve_para[crv1].second);

						int nbpt_idx0 = (endpt_idx0 == curve0->m_cpoints.front()) ? 
										curve0->m_cpoints[1] : 
										curve0->m_cpoints[curve0->m_cpoints.size() - 2];
						int nbpt_idx1 = (endpt_idx1 == curve1->m_cpoints.front()) ? 
										curve1->m_cpoints[1] : 
										curve1->m_cpoints[curve1->m_cpoints.size() - 2];

						Eigen::Vector3d nbpt0(
							m_control_points[nbpt_idx0].x(), 
							m_control_points[nbpt_idx0].y(), 
							m_control_points[nbpt_idx0].z()
						);
						Eigen::Vector3d nbpt1(
							m_control_points[nbpt_idx1].x(), 
							m_control_points[nbpt_idx1].y(), 
							m_control_points[nbpt_idx1].z()
						);

						if ((nbpt0 - ptbin_itr.pt).dot(nbpt1 - ptbin_itr.pt) < 0)
						{
							TgtBin ntgt_bin(ptbin_itr.pt, tgt0.normalized());
							ntgt_bin.v_curve_para.push_back(ptbin_itr.v_curve_para[crv0]);
							ntgt_bin.v_curve_para.push_back(ptbin_itr.v_curve_para[crv1]);

							tgt_bins.push_back(ntgt_bin);
						}
					}
				}
			}
		}

		//build index map.
		std::vector<int> map_vtx_var (m_control_points.size(), -1);
		int var_count=0;
		for (auto curve : m_curves)
		{
			for (int ccpt_itr = 1; ccpt_itr < curve->m_cpoints.size()-1; ccpt_itr++)
			{
				map_vtx_var[curve->m_cpoints[ccpt_itr]] = var_count ++;
			}
		}
		for (auto& pt_bin : pt_bins)
		{
			if (pt_bin.v_curve_para.size()<3)
			{
				continue;
			}
			for (auto curve_para : pt_bin.v_curve_para)
			{
				int cpt_idx = curve_para.first->get_endpt_idx(curve_para.second);
				map_vtx_var[cpt_idx] = var_count;
			}
			var_count++;
		}

		//build matrix.
		Eigen::SparseMatrix<double> fairing_mat(var_count * 3, var_count * 3);
		Eigen::VectorXd rhs(var_count * 3);  // RHS vector

		//closeness term.
		double closeness_coeff = 1.0;

		for (auto curve : m_curves)
		{
			for (int ccpt_itr = 1; ccpt_itr < curve->m_cpoints.size() - 1; ccpt_itr++)
			{
				int vtxitr = curve->m_cpoints[ccpt_itr];
				if (map_vtx_var[vtxitr] > -1)
				{
					for (int coord_itr = 0; coord_itr < 3; coord_itr++)
					{
						int index_ij = 3 * map_vtx_var[vtxitr] + coord_itr;
						
						// Insert closeness_coeff into the sparse matrix at (index_ij, index_ij)
						fairing_mat.insert(index_ij, index_ij) = closeness_coeff;

						// Update the RHS vector at the corresponding index
						rhs[index_ij] = closeness_coeff * m_control_points[vtxitr][coord_itr];
					}
				}
			}
		}

		for (auto& pt_bin : pt_bins)
		{
			if (pt_bin.v_curve_para.size() < 3)
			{
				continue;
			}

			int cpt_idx = pt_bin.v_curve_para[0].first->get_endpt_idx(pt_bin.v_curve_para[0].second);

			if (map_vtx_var[cpt_idx] > -1)
			{
				for (int coord_itr = 0; coord_itr < 3; coord_itr++)
				{
					int index_ij = 3 * map_vtx_var[cpt_idx] + coord_itr;
					
					// Insert closeness_coeff into the sparse matrix at (index_ij, index_ij)
					fairing_mat.insert(index_ij, index_ij) = closeness_coeff;

					// Update the RHS vector at the corresponding index
					rhs[index_ij] = closeness_coeff * m_control_points[cpt_idx][coord_itr];
				}
			}
		}

		//fairing term.
		double intra_curve_fairing_coeff = 10.0;
		double inter_curve_fairing_coeff = 10.0;
		for (auto curve : m_curves)
		{
			for (int ccpt_itr = 0; ccpt_itr < curve->m_cpoints.size()-3; ccpt_itr++)
			{
				double diff [4] = {-1.0, 3.0, -3.0, 1.0};
				int varid [4] = {map_vtx_var[curve->m_cpoints[ccpt_itr]], map_vtx_var[curve->m_cpoints[ccpt_itr+1]],
									map_vtx_var[curve->m_cpoints[ccpt_itr+2]], map_vtx_var[curve->m_cpoints[ccpt_itr+3]]};
				Vector4d vtx_pos [4] = {m_control_points[curve->m_cpoints[ccpt_itr]], m_control_points[curve->m_cpoints[ccpt_itr+1]],
										m_control_points[curve->m_cpoints[ccpt_itr+2]], m_control_points[curve->m_cpoints[ccpt_itr+3]]};

				for (int coord_itr = 0; coord_itr < 3; coord_itr++)
				{
					for (int vitr0 = 0; vitr0 < 4; vitr0++)
					{
						if (varid[vitr0]<0)
						{
							continue;
						}
						int rid = 3*varid[vitr0] + coord_itr;
						for (int vitr1 = 0; vitr1 < 4; vitr1++)
						{
							if (varid[vitr1] < 0)
							{
								rhs[rid] -= intra_curve_fairing_coeff * diff[vitr0] * diff[vitr1] * vtx_pos[vitr1][coord_itr];
							}
							else
							{
								int cid = 3 * varid[vitr1] + coord_itr;
								if (rid >= cid)
								{
									fairing_mat.insert(rid, cid) = intra_curve_fairing_coeff * diff[vitr0] * diff[vitr1];
								}
							}
						}
					}
				}
			}
		}

		//inter curve
		for (auto tgt_bin : tgt_bins)
		{
			if (tgt_bin.v_curve_para.size()<2)
			{
				continue;
			}
			SplineCurv* crv0 = tgt_bin.v_curve_para[0].first, *crv1 = tgt_bin.v_curve_para[1].first;
			double para0 = tgt_bin.v_curve_para[0].second, para1 = tgt_bin.v_curve_para[1].second;

			int cpt_idx_0 [3] = {-1,-1,-1};
			int cpt_idx_1 [3] = {-1,-1,-1};

			int idx = crv0->get_endpt_idx(para0);
			if (idx==crv0->m_cpoints.front())
			{
				cpt_idx_0[0] = crv0->m_cpoints[0];
				cpt_idx_0[1] = crv0->m_cpoints[1];
				cpt_idx_0[2] = crv0->m_cpoints[2];
			}
			else
			{
				int back_idx = crv0->m_cpoints.size()-1;
				cpt_idx_0[0] = crv0->m_cpoints[back_idx];
				cpt_idx_0[1] = crv0->m_cpoints[back_idx-1];
				cpt_idx_0[2] = crv0->m_cpoints[back_idx-2];
			}

			idx = crv1->get_endpt_idx(para1);
			if (idx == crv1->m_cpoints.front())
			{
				cpt_idx_1[0] = crv1->m_cpoints[0];
				cpt_idx_1[1] = crv1->m_cpoints[1];
				cpt_idx_1[2] = crv1->m_cpoints[2];
			}
			else
			{
				int back_idx = crv1->m_cpoints.size()-1;
				cpt_idx_1[0] = crv1->m_cpoints[back_idx];
				cpt_idx_1[1] = crv1->m_cpoints[back_idx-1];
				cpt_idx_1[2] = crv1->m_cpoints[back_idx-2];
			}

			if (map_vtx_var[cpt_idx_0[0]] != map_vtx_var[cpt_idx_1[0]])
			{
				cs_err() << "Error: endpoint index does not match." << std::endl;
			}

			int tuple_vtx_idx [5] = {cpt_idx_0[2], cpt_idx_0[1], cpt_idx_0[0], cpt_idx_1[1], cpt_idx_1[2]};
			for (int ccpt_itr = 0; ccpt_itr < 2; ccpt_itr++)
			{
				double diff [4] = {-1.0, 3.0, -3.0, 1.0};
				int varid [4] = {map_vtx_var[tuple_vtx_idx[ccpt_itr]], map_vtx_var[tuple_vtx_idx[ccpt_itr+1]],
									map_vtx_var[tuple_vtx_idx[ccpt_itr+2]], map_vtx_var[tuple_vtx_idx[ccpt_itr+3]]};
				Vector4d vtx_pos [4] = {m_control_points[tuple_vtx_idx[ccpt_itr]], m_control_points[tuple_vtx_idx[ccpt_itr+1]],
										m_control_points[tuple_vtx_idx[ccpt_itr+2]], m_control_points[tuple_vtx_idx[ccpt_itr+3]]};

				for (int coord_itr = 0; coord_itr < 3; coord_itr++)
				{
					for (int vitr0 = 0; vitr0 < 4; vitr0++)
					{
						if (varid[vitr0]<0)
						{
							continue;
						}
						int rid = 3*varid[vitr0] + coord_itr;
						for (int vitr1 = 0; vitr1 < 4; vitr1++)
						{
							if (varid[vitr1]<0)
							{
								rhs[rid] -= inter_curve_fairing_coeff * diff[vitr0] * diff[vitr1] * vtx_pos[vitr1][coord_itr];
							}
							else
							{
								int cid = 3*varid[vitr1] + coord_itr;
								if (rid >= cid)
								{
									fairing_mat.insert(rid, cid) = inter_curve_fairing_coeff * diff[vitr0] * diff[vitr1];
								}
							}
						}
					}
				}
			}
		}

		//solve the optimization.
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
		solver.setMaxIterations(kFairingCGMaxIterations);
		solver.setTolerance(kFairingCGTolerance);
		Eigen::VectorXd solution = solver.compute(fairing_mat).solve(rhs);
		//read solution.
		// Read the solution and update control points
		for (int cpitr = 0; cpitr < m_control_points.size(); cpitr++)
		{
			if (map_vtx_var[cpitr] < 0)
			{
				continue;
			}
			Eigen::Vector3d npt(
				solution[map_vtx_var[cpitr] * 3],
				solution[map_vtx_var[cpitr] * 3 + 1],
				solution[map_vtx_var[cpitr] * 3 + 2]
			);

			m_control_points[cpitr].x() = npt.x();
			m_control_points[cpitr].y() = npt.y();
			m_control_points[cpitr].z() = npt.z();
		}

		for (auto crv : m_curves)
		{
			crv->UpdateCurve();
		}
	}

	CS_INLINE void CurvNet::collect_g1_cycle_pairs(){
		if (m_p_polygraph_)
		{
			m_g1_cyclepairs.clear();
			for (int curv_itr = 0; curv_itr < m_curves.size(); curv_itr++)
			{
				SplineCurv* curve = m_curves[curv_itr];
				for (int ccyc_itr_0 = 0; ccyc_itr_0 < curve->pg_edges.size(); ccyc_itr_0++)
				{
					PolyGraph_Edge* pgedge_0 = curve->pg_edges[ccyc_itr_0];
					double para[2];
					Eigen::Vector3d normal[2];
					computeCornerNormals(pgedge_0, normal, para);
					for (int ccyc_itr_1 = ccyc_itr_0+1; ccyc_itr_1 < curve->pg_edges.size(); ccyc_itr_1++)
					{
						PolyGraph_Edge* pgedge_1 = curve->pg_edges[ccyc_itr_1];
						double para_inner[2];
						Eigen::Vector3d normal_inner[2];
						computeCornerNormals(pgedge_1, normal_inner, para_inner);

						const double SINE_G1_SMOOTH_ANGLE = std::sin(10.0/180.0*M_PI);
						if (normal[0].cross(normal_inner[0]).norm() < SINE_G1_SMOOTH_ANGLE && 
							normal[1].cross(normal_inner[1]).norm() < SINE_G1_SMOOTH_ANGLE)
						{
							small_set<PolyGraph_Cycle*, 2> pair;
							pair.insert(pgedge_0->cycle_);
							pair.insert(pgedge_1->cycle_);
							m_g1_cyclepairs.insert(pair);
						}
					}
				}
			}
		}
	}

	CS_INLINE void CurvNet::save_cycles(){

		std::string filename = data_file_name + "_marked.cycle";

		cs_log() << "saving curve cycle data to " << filename <<" ... ";

		int cycle_num = m_p_polygraph_->get_num_of_cycles();
		std::ofstream ofile (filename);

		if (ofile.is_open())
		{
			ofile << cycle_num << std::endl;

			for (int cyc_itr = 0; cyc_itr < cycle_num; cyc_itr++)
			{
				PolyGraph_Cycle* cycle = m_p_polygraph_->get_cycle(cyc_itr);
				ofile << cycle->degree << std::endl;

				PolyGraph_Edge* cycle_edge = cycle->edge;
				do
				{
					ofile << cycle_edge->curve_->id << " " << (cycle_edge->order_?1:2) << " " << cycle_edge->curve_type << std::endl;

					cycle_edge = cycle_edge->next;
				} while (cycle_edge != cycle->edge);
			}

			ofile.close();
		}

		cs_log() << "end." << std::endl;
	}

	// Write a human-readable .curve_types file listing, for every SplineCurv, whether it
	// is a curvature line (CRVT_LINE) or a trim line (TRIM_LINE).
	//
	// The classification is determined by the PolyGraph_Edge::curve_type values that are
	// assigned during computeRMF / populateRMFFromNormalData.  A curve whose pg_edges
	// disagree (rare edge case) is reported as MIXED with a per-edge breakdown.
	//
	// File format:
	//   # curve_id  curve_type  (CRVT_LINE | TRIM_LINE | MIXED)
	//   <curve_id> <type_string>
	//   ...
	CS_INLINE void CurvNet::save_curve_types() {

		std::string filename = data_file_name + ".curve_types";

		cs_log() << "saving curve type listing to " << filename << " ... ";

		std::ofstream ofile(filename);
		if (!ofile.is_open()) {
			cs_err() << "ERROR: could not open " << filename << std::endl;
			return;
		}

		ofile << "# curve_id  curve_type  (CRVT_LINE | TRIM_LINE | MIXED)\n";

		int total_crvt  = 0;
		int total_trim  = 0;
		int total_mixed = 0;

		for (SplineCurv* curve : m_curves) {
			if (!curve) continue;

			int n_crvt = 0;
			int n_trim = 0;
			for (PolyGraph_Edge* pge : curve->pg_edges) {
				if (!pge) continue;
				if (pge->curve_type == CURVE_TYPE::CRVT_LINE) ++n_crvt;
				else                                           ++n_trim;
			}

			std::string type_str;
			if (n_crvt > 0 && n_trim == 0) {
				type_str = "CRVT_LINE";
				++total_crvt;
			} else if (n_trim > 0 && n_crvt == 0) {
				type_str = "TRIM_LINE";
				++total_trim;
			} else if (n_crvt == 0 && n_trim == 0) {
				// No pg_edges yet — type not yet determined; default is CRVT_LINE (see ctor)
				type_str = "CRVT_LINE";
				++total_crvt;
			} else {
				// Mixed: majority vote, report detail
				type_str = (n_crvt >= n_trim) ? "MIXED(majority:CRVT_LINE)" : "MIXED(majority:TRIM_LINE)";
				++total_mixed;
			}

			ofile << curve->id << " " << type_str << "\n";
		}

		ofile << "# summary: " << total_crvt  << " CRVT_LINE, "
		                       << total_trim  << " TRIM_LINE, "
		                       << total_mixed << " MIXED\n";
		ofile.close();

		cs_log() << "end.  ("
		          << total_crvt  << " curvature, "
		          << total_trim  << " trim, "
		          << total_mixed << " mixed)\n";
	}

	CS_INLINE void CurvNet::save_curve_types(const std::string& file_prefix) {
		data_file_name = file_prefix;
		save_curve_types();
	}

	CS_INLINE double CurvNet::calculateAverageBoundaryEdgeLength(Mesh3D* mesh3d) {
    double total_length = 0.0;
    int boundary_edge_count = 0;

    for (int e = 0; e < mesh3d->get_num_of_edges(); ++e) {
        HE_edge* edge = mesh3d->get_edge(e);
        if (!edge || !edge->vert || !edge->pair || !edge->pair->vert) continue;

        // Check if edge is on boundary
        if (mesh3d->is_on_boundary(edge)) {
            Eigen::Vector3d v0 = edge->vert->pos;
            Eigen::Vector3d v1 = edge->pair->vert->pos;
            double length = (v1 - v0).norm();
            total_length += length;
            boundary_edge_count++;
        }
    }

    return boundary_edge_count > 0 ? total_length / boundary_edge_count : 0.0;
}

// Face splitting function without visualization - takes a mesh and performs optimized face splitting
CS_INLINE void CurvNet::split_faces_optimized(Mesh3D* mesh3d, int max_iterations, int batch_size)
{
    if (!mesh3d) {
        cs_err() << "Error: mesh3d is null" << std::endl;
        return;
    }

    // 1. OPTIMIZED AREA CALCULATION (assumes triangular mesh)
    auto calculateTriangleArea = [](HE_face* face) -> double {
        if (!face || !face->edge) return 0.0;
        
        HE_edge* e1 = face->edge;
        HE_edge* e2 = e1->next;
        HE_edge* e3 = e2->next;
        
        if (!e1 || !e2 || !e3 || !e1->vert || !e2->vert || !e3->vert) return 0.0;
        
        const Eigen::Vector3d& v0 = e1->vert->pos;
        const Eigen::Vector3d& v1 = e2->vert->pos;
        const Eigen::Vector3d& v2 = e3->vert->pos;
        
        // Direct cross product calculation (more efficient)
        Eigen::Vector3d cross = (v1 - v0).cross(v2 - v0);
        return 0.5 * cross.norm();
    };

    // 2. PRIORITY QUEUE for efficient largest face selection
    struct FaceAreaPair {
        HE_face* face;
        double area;
        int face_id;
        
        bool operator<(const FaceAreaPair& other) const {
            return area < other.area; // Max-heap
        }
    };

    std::priority_queue<FaceAreaPair> face_queue;
    std::unordered_set<HE_face*> valid_faces;
    
    for (int f = 0; f < mesh3d->get_num_of_faces(); ++f) {
        HE_face* face = mesh3d->get_face(f);
        if (face) {
            double area = calculateTriangleArea(face);
            face_queue.push({face, area, f});
            valid_faces.insert(face);
        }
    }
    
    // === MAIN SPLITTING LOOP ===
    int iteration = 0;
    int splits_since_update = 0;
    int total_batch_updates = 0;
    
    while (iteration < max_iterations && !face_queue.empty()) {
        // Get largest face
        FaceAreaPair largest = face_queue.top();
        face_queue.pop();
        
        // Validate face (skip if invalidated by previous splits)
        if (valid_faces.find(largest.face) == valid_faces.end()) {
            continue;
        }
        
        // 4. SHAPE-PRESERVING SPLIT POINT CALCULATION
        HE_edge* e1 = largest.face->edge;
        HE_edge* e2 = e1->next;
        HE_edge* e3 = e2->next;
        
        if (!e1 || !e2 || !e3 || !e1->vert || !e2->vert || !e3->vert) {
            continue;
        }
        
        const Eigen::Vector3d& v0 = e1->vert->pos;
        const Eigen::Vector3d& v1 = e2->vert->pos;
        const Eigen::Vector3d& v2 = e3->vert->pos;
        
        // SHAPE-PRESERVING SPLIT POINT OPTIONS:
        
        // Option 1: Planar centroid (stays in triangle plane)
        Eigen::Vector3d planar_centroid = (v0 + v1 + v2) / 3.0;
        
        // Option 3: Edge-midpoint average (better for thin triangles)
        Eigen::Vector3d edge_midpoint_avg = ((v0 + v1) + (v1 + v2) + (v2 + v0)) / 6.0;
        
        // Choose split point based on triangle shape
        double edge1_len = (v1 - v0).norm();
        double edge2_len = (v2 - v1).norm();
        double edge3_len = (v0 - v2).norm();
        double max_edge = std::max({edge1_len, edge2_len, edge3_len});
        double min_edge = std::min({edge1_len, edge2_len, edge3_len});
        double aspect_ratio = max_edge / min_edge;
        
        Eigen::Vector3d split_point;
        if (aspect_ratio > 3.0) {
            // For thin triangles: use edge midpoint average for better shape
            split_point = edge_midpoint_avg;
        } else {
            // For well-shaped triangles: use planar centroid
            split_point = planar_centroid;
        }
        
        // Create vertex and split face
        HE_vert* new_vert = mesh3d->insert_vertex(split_point);
        valid_faces.erase(largest.face);
        mesh3d->split_face(largest.face, new_vert);
        
        iteration++;
        splits_since_update++;
        
        // 5. BATCH MESH UPDATES 
        if (splits_since_update >= batch_size || iteration >= max_iterations) {
            auto batch_start = std::chrono::high_resolution_clock::now();

            bool topo = false;
            mesh3d->make_delaunay(topo);

            // Rebuild priority queue
            while (!face_queue.empty()) face_queue.pop();
            valid_faces.clear();
            
            for (int f = 0; f < mesh3d->get_num_of_faces(); ++f) {
                HE_face* face = mesh3d->get_face(f);
                if (face) {
                    double area = calculateTriangleArea(face);
                    face_queue.push({face, area, f});
                    valid_faces.insert(face);
                }
            }
            
            auto batch_end = std::chrono::high_resolution_clock::now();
            auto batch_time = std::chrono::duration_cast<std::chrono::milliseconds>(batch_end - batch_start);
            
            total_batch_updates++;
            splits_since_update = 0;
        
        }
    }
    
    // Final mesh update test_group_patches
    mesh3d->update_mesh();
    mesh3d->update_normal();
    bool topo = false;
    mesh3d->make_delaunay(topo); 
}




CS_INLINE void CurvNet::split_edge(Mesh3D* mesh3d, int max_iterations)
{
    // Iteratively split longest interior edges
    double avg_boundary_length = calculateAverageBoundaryEdgeLength(mesh3d);

    int iteration = 0;
    while (iteration < max_iterations) {
        double max_len = -1.0;
        HE_edge* longest_edge = nullptr;

        // Find the longest interior edge (avoid boundary edges)
        for (int e = 0; e < mesh3d->get_num_of_edges(); ++e) {
            HE_edge* edge = mesh3d->get_edge(e);
            if (!edge || !edge->vert || !edge->pair || !edge->pair->vert) continue;
            
            // Skip boundary edges - only split interior edges
            if (mesh3d->is_on_boundary(edge)) continue;
            
            // Make sure edge has two adjacent faces (required for split_edge)
            if (!edge->face || !edge->pair->face) continue;
            
            Eigen::Vector3d v0 = edge->vert->pos;
            Eigen::Vector3d v1 = edge->pair->vert->pos;
            double len = (v1 - v0).norm();
            if (len > max_len) {
                max_len = len;
                longest_edge = edge;
            }
        }

        // Check if longest interior edge is significantly longer than average boundary length
        if (longest_edge && max_len > avg_boundary_length * 1.5) { // tolerance
            // Get vertices of the edge to split
            HE_vert* v0 = longest_edge->vert;
            HE_vert* v1 = longest_edge->pair->vert;

            // Create new vertex at midpoint
            Eigen::Vector3d midpoint = 0.5 * (v0->pos + v1->pos);
            HE_vert* new_vert = mesh3d->insert_vertex(midpoint);
           
            mesh3d->split_edge(longest_edge, new_vert);
            iteration++;
        } 
        else {
            break;
        }

        if (iteration % 100 == 0) {
            bool topo;
            mesh3d->make_delaunay(topo);
        }
    }

    // Apply constraint-based Laplacian smoothing to improve mesh quality after edge splitting
    auto constraint_smooth_start = std::chrono::high_resolution_clock::now();

    mesh3d->laplacian_smooth_constraints(1.0, true);

    auto constraint_smooth_end = std::chrono::high_resolution_clock::now();
    auto constraint_smooth_time = std::chrono::duration_cast<std::chrono::milliseconds>(constraint_smooth_end - constraint_smooth_start);
}



CS_INLINE Mesh3D* CurvNet::combinePatchMeshes(const std::vector<Mesh3D*>& patch_meshes)
{
    Mesh3D* combined_mesh = new Mesh3D();
    std::vector<int> patch_vertex_offsets; // Track vertex offset for each patch
    int total_vertices = 0;
    int total_faces = 0;
    
    // First pass: calculate total counts and vertex offsets
    for (unsigned patch_id = 0; patch_id < patch_meshes.size(); patch_id++) {
        Mesh3D* patch_mesh = patch_meshes[patch_id];
        if (patch_mesh) {
            patch_vertex_offsets.push_back(total_vertices);
            total_vertices += patch_mesh->get_num_of_vertices();
            total_faces += patch_mesh->get_num_of_faces();
        } else {
            patch_vertex_offsets.push_back(-1); // Mark null patches
        }
    }
    
    // Second pass: insert all vertices from all patches (no merging)
    std::vector<std::vector<HE_vert*>> patch_vertex_maps; // Map original patch vertices to combined mesh vertices
    patch_vertex_maps.resize(patch_meshes.size());
    
    for (unsigned patch_id = 0; patch_id < patch_meshes.size(); patch_id++) {
        Mesh3D* patch_mesh = patch_meshes[patch_id];
        if (!patch_mesh) continue;
        
        auto* patch_verts = patch_mesh->get_vertices_list();
        patch_vertex_maps[patch_id].resize(patch_verts->size());
        
        // Insert each vertex from this patch into the combined mesh
        for (unsigned v_id = 0; v_id < patch_verts->size(); v_id++) {
            HE_vert* original_vert = (*patch_verts)[v_id];
            HE_vert* new_vert = combined_mesh->insert_vertex(original_vert->pos);
            patch_vertex_maps[patch_id][v_id] = new_vert;
        }
    }
    
    // Third pass: insert all faces using the vertex mapping
    for (unsigned patch_id = 0; patch_id < patch_meshes.size(); patch_id++) {
        Mesh3D* patch_mesh = patch_meshes[patch_id];
        if (!patch_mesh) continue;
        
        auto* patch_faces = patch_mesh->get_faces_list();
        
        for (unsigned f_id = 0; f_id < patch_faces->size(); f_id++) {
            HE_face* original_face = (*patch_faces)[f_id];
            auto face_vertices = original_face->get_vertices();
            
            if (face_vertices.size() >= 3) {
                std::vector<HE_vert*> new_face_verts;
                
                // Map the face vertices to the combined mesh vertices
                for (unsigned i = 0; i < 3; i++) { // Only use first 3 vertices for triangles
                    int original_v_id = face_vertices[i]->id;
                    if (original_v_id < patch_vertex_maps[patch_id].size()) {
                        new_face_verts.push_back(patch_vertex_maps[patch_id][original_v_id]);
                    }
                }
                
                if (new_face_verts.size() == 3) {
                    HE_face* new_face = combined_mesh->insert_face(new_face_verts);
                    
                    // Preserve the curve cycle association from the original face
                    if (new_face && original_face->curve_cycle) {
                        new_face->curve_cycle = original_face->curve_cycle;
                    }
                }
            }
        }
    }

    return combined_mesh;
}

	// Helper function for Laplacian smoothing with fixed boundaries
	CS_INLINE void laplacian_smooth_patch(Mesh3D* patch_mesh, int iterations = 5, double lambda = 0.5) {
		if (!patch_mesh || patch_mesh->get_num_of_vertices() == 0) {
			return;
		}
		
		auto* vert_list = patch_mesh->get_vertices_list();
		
		// Identify boundary vertices (keep them fixed)
		std::unordered_set<HE_vert*> boundary_verts;
		for (HE_vert* v : *vert_list) {
			if (patch_mesh->is_on_boundary(v)) {
				boundary_verts.insert(v);
			}
		}
		
		// Perform Laplacian smoothing iterations
		for (int iter = 0; iter < iterations; ++iter) {
			// Store new positions (don't update in-place to avoid order dependency)
			std::unordered_map<HE_vert*, Eigen::Vector3d> new_positions;
			
			for (HE_vert* v : *vert_list) {
				// Skip boundary vertices
				if (boundary_verts.count(v)) {
					continue;
				}
				
				// Calculate centroid of neighboring vertices
				Eigen::Vector3d centroid(0, 0, 0);
				int neighbor_count = 0;
				
				// Traverse one-ring neighborhood using outgoing edges
				HE_edge* start_edge = v->edge;
				if (!start_edge) continue;
				
				HE_edge* current_edge = start_edge;
				do {
					// Get the vertex at the end of this edge
					if (current_edge->next && current_edge->next->vert) {
						centroid += current_edge->next->vert->pos;
						neighbor_count++;
					}
					// Move to next outgoing edge from this vertex
					if (current_edge->pair) {
						current_edge = current_edge->pair->next;
					} else {
						break;  // Hit boundary
					}
				} while (current_edge && current_edge != start_edge);
				
				if (neighbor_count > 0) {
					centroid /= neighbor_count;
					// Weighted average between current position and centroid
					new_positions[v] = (1.0 - lambda) * v->pos + lambda * centroid;
				}
			}
			
			// Apply new positions
			for (auto& pair : new_positions) {
				pair.first->pos = pair.second;
			}
		}
		
		// Update mesh normals after smoothing
		patch_mesh->update_normal();
	}
	
	// Helper function to remesh a single patch using the remesher library
	CS_INLINE void remesh_patch(Mesh3D* patch_mesh, const std::vector<PolyGraph_Edge*>& cycle_edges, double resolution_multiplier = 1.5) {
		if (!patch_mesh || patch_mesh->get_num_of_faces() == 0) {
			return;
		}
		
		// Calculate average boundary edge length from the cycle curves
		double total_boundary_length = 0.0;
		int boundary_sample_count = 0;
		
		for (PolyGraph_Edge* pg_edge : cycle_edges) {
			SplineCurv* curve = pg_edge->curve_;
			// Sum up distances between consecutive control points
			for (size_t i = 0; i < curve->m_cpoints.size() - 1; i++) {
				int idx0 = curve->m_cpoints[i];
				int idx1 = curve->m_cpoints[i + 1];
				Vector4d& cp0 = curve->net->m_control_points[idx0];
				Vector4d& cp1 = curve->net->m_control_points[idx1];
				
				Eigen::Vector3d p0(cp0.x(), cp0.y(), cp0.z());
				Eigen::Vector3d p1(cp1.x(), cp1.y(), cp1.z());
				double edge_len = (p1 - p0).norm();
				total_boundary_length += edge_len;
				boundary_sample_count++;
			}
		}
		
		double avg_boundary_length = (boundary_sample_count > 0) 
			? total_boundary_length / boundary_sample_count 
			: 0.01;  // Fallback
		
		// Use boundary edge length as target (with optional multiplier)
		double target_length = resolution_multiplier * avg_boundary_length;
		
		// Convert Mesh3D to Eigen matrices
		int num_verts = patch_mesh->get_num_of_vertices();
		int num_faces = patch_mesh->get_num_of_faces();
		
		Eigen::MatrixXd V(num_verts, 3);
		Eigen::MatrixXi F(num_faces, 3);
		
		// Fill vertices
		auto* vert_list = patch_mesh->get_vertices_list();
		int v_idx = 0;
		for (HE_vert* v : *vert_list) {
			V(v_idx, 0) = v->pos.x();
			V(v_idx, 1) = v->pos.y();
			V(v_idx, 2) = v->pos.z();
			v->id = v_idx;  // Update ID for face indexing
			v_idx++;
		}
		
		// Fill faces
		auto* face_list = patch_mesh->get_faces_list();
		int f_idx = 0;
		for (HE_face* f : *face_list) {
			HE_edge* e = f->edge;
			F(f_idx, 0) = e->vert->id;
			e = e->next;
			F(f_idx, 1) = e->vert->id;
			e = e->next;
			F(f_idx, 2) = e->vert->id;
			f_idx++;
		}
		
		// Setup remeshing parameters
		Geo::TriRemParam<double> RemP;
		RemP.steps = 6;  
		RemP.fix_borders = true;
		RemP.iterative = false;
		RemP.adaptive = false;
		RemP.feature_angle = -1;  // No feature detection
		RemP.collapse_small = true;
		RemP.refine_big = true;
		RemP.writeMSG = false;  // Suppress messages for individual patches
		
		RemP.TargetL = target_length;

		// Initialize feature arrays
		RemP.IsVertCorner.resize(V.rows(), false);
		RemP.IsFaceEFeature.resize(F.rows(), std::vector<bool>(3, false));
		
		// Perform remeshing
Geo::TriRemesh<double>(V, F, RemP);
		
		// Rebuild Mesh3D from remeshed data
		// Delete all existing mesh elements
		auto* old_verts = patch_mesh->get_vertices_list();
		auto* old_edges = patch_mesh->get_edges_list();
		auto* old_faces = patch_mesh->get_faces_list();
		
		// Clear lists (note: this doesn't delete the objects, just removes them from lists)
		old_verts->clear();
		old_edges->clear();
		old_faces->clear();
		
		// Insert new vertices
		std::vector<HE_vert*> new_verts;
		for (int i = 0; i < V.rows(); i++) {
			Eigen::Vector3d pos(V(i, 0), V(i, 1), V(i, 2));
			HE_vert* v = patch_mesh->insert_vertex(pos);
			new_verts.push_back(v);
		}
		
		// Insert new faces
		for (int i = 0; i < F.rows(); i++) {
			std::vector<HE_vert*> face_verts = {
				new_verts[F(i, 0)],
				new_verts[F(i, 1)],
				new_verts[F(i, 2)]
			};
			patch_mesh->insert_face(face_verts);
		}
		
		// Update mesh
		patch_mesh->update_mesh();
		patch_mesh->update_normal();
	}
	
	CS_INLINE Mesh3D* CurvNet::load_cycle_discovery_data_from_structures(
		const std::vector<CurveData>& curves_data,
		const std::vector<CycleData>& cycles_data,
		const FeatureEdgeData& feature_edges) {
		
		// Load curves into CurvNet
		for (size_t c_itr = 0; c_itr < curves_data.size(); c_itr++) {
			const CurveData& curve_data = curves_data[c_itr];
			
			SplineCurv* new_curv = new SplineCurv(this);
			new_curv->id = c_itr;
			new_curv->degree = 1;
			new_curv->type = SPLINE_TYPE::BSPLINE;
			
			new_curv->para_range[0] = 0;
			new_curv->para_range[1] = 1;
			new_curv->knots.push_back(0);
			
			for (size_t pt_itr = 0; pt_itr < curve_data.points.size(); pt_itr++) {
				const Eigen::Vector3d& pt = curve_data.points[pt_itr];
				this->m_control_points.push_back(Vector4d(pt.x(), pt.y(), pt.z(), 0));
				new_curv->m_cpoints.push_back(this->m_control_points.size() - 1);
				new_curv->knots.push_back((double)pt_itr / (double)(curve_data.num_pts - 1));
			}
			new_curv->knots.push_back(1);
			
			this->m_curves.push_back(new_curv);
		}
		
		cs_log() << "Curves loaded into CurvNet." << std::endl;
		
		// Create combined mesh
		Mesh3D* nmesh = new Mesh3D;
		
		// Initialize polygraph
		if (!m_p_polygraph_) {
			m_p_polygraph_ = new PolyGraph;
		} else {
			delete m_p_polygraph_;
			m_p_polygraph_ = new PolyGraph;
		}
		
		// Load cycles and create patches using SurfaceMaker
		std::vector<Mesh3D*> temp_patch_meshes; // Store temporary meshes to delete later
		
		for (size_t cy_itr = 0; cy_itr < cycles_data.size(); cy_itr++) {
			const CycleData& cycle_data = cycles_data[cy_itr];
			
			// Create cycle edges
			std::vector<PolyGraph_Edge*> cycle_edges;
			for (const auto& curve_info : cycle_data.curve_info) {
				int curv_id = curve_info[0];
				int curv_dir = curve_info[1];
				int instance = curve_info[2];
				
				cycle_edges.push_back(m_p_polygraph_->insert_edge(m_curves[curv_id], (curv_dir == 1)));
			}
			
			PolyGraph_Cycle* cur_cycle = m_p_polygraph_->insert_cycle(cycle_edges);
			
			// Create patch mesh using SurfaceMaker
			cs_log() << "Creating surface patch for cycle " << cy_itr << "..." << std::endl;
			
			// Create a new mesh for this patch
			Mesh3D* patch_mesh = new Mesh3D;
			
			// Convert curve data to AML::double3 format expected by SurfaceMaker
			std::vector<std::vector<AML::double3>> cycle_curve_network;
			
			for (PolyGraph_Edge* edge : cycle_edges) {
				SplineCurv* curve = edge->curve_;
				std::vector<AML::double3> curve_points;
				
				// Sample points along the curve
				for (int cpoint_idx : curve->m_cpoints) {
					Vector4d& cp = this->m_control_points[cpoint_idx];
					AML::double3 aml_point;
					aml_point[0] = cp.x();
					aml_point[1] = cp.y(); 
					aml_point[2] = cp.z();
					curve_points.push_back(aml_point);
				}
				
				// Reverse the curve if needed based on order
				if (edge->order_) {
					std::reverse(curve_points.begin(), curve_points.end());
				}
				
				cycle_curve_network.push_back(curve_points);
			}
			
			// Create SurfaceMaker instance and generate surface
			SurfaceMaker surface_maker;
			surface_maker.OnOpenCurveNetwork(cycle_curve_network);
			
			// Get the generated triangle surface
			const GraphStructure::TriangleSurface& triangle_surface = surface_maker.GetTriangleSurface();
			
			if (!triangle_surface.empty()) {
				cs_log() << "Generated surface with " << triangle_surface.size() << " triangle cycles" << std::endl;
				
				// Extract points and faces for this specific patch (with proper vertex deduplication)
				std::vector<AML::double3> patch_points;
				std::vector<std::vector<unsigned>> patch_faces;
				
				// Process triangles
				for (const GraphStructure::TriangleCycle& triangle_cycle : triangle_surface) {
					for (unsigned j = 0; j < triangle_cycle.size(); j++) {
						const std::vector<AML::double3>& triangle = triangle_cycle[j];
						std::vector<unsigned> tf;
						
						for (unsigned k = 0; k < triangle.size(); k++) {
							const AML::double3& tp = triangle[k];
							unsigned pID;
							
							// Find if point already exists in this patch's point list
							auto it = std::find(patch_points.begin(), patch_points.end(), tp);
							if (it == patch_points.end()) {
								patch_points.push_back(tp);
								pID = patch_points.size() - 1;
							} else {
								pID = it - patch_points.begin();
							}
							tf.push_back(pID);
						}
						patch_faces.push_back(tf);
					}
				}
				
				cs_log() << "Patch " << cy_itr << ": " << patch_points.size() << " points, " << patch_faces.size() << " faces" << std::endl;
				
				// Insert vertices into the patch mesh
				std::vector<HE_vert*> patch_vertices;
				for (unsigned i = 0; i < patch_points.size(); i++) {
					const AML::double3& coord = patch_points[i];
					Eigen::Vector3d v(coord[0], coord[1], coord[2]);
					HE_vert* vert = patch_mesh->insert_vertex(v);
					patch_vertices.push_back(vert);
				}
				
				// Insert faces into the patch mesh
				for (unsigned i = 0; i < patch_faces.size(); i++) {
					const std::vector<unsigned>& tri = patch_faces[i];
					if (tri.size() >= 3) {
						std::vector<HE_vert*> face_verts;
						for (unsigned j = 0; j < 3; j++) {  // Only use first 3 vertices for triangles
							if (tri[j] < patch_vertices.size()) {
								face_verts.push_back(patch_vertices[tri[j]]);
							}
						}
						if (face_verts.size() == 3) {
							HE_face* fh = patch_mesh->insert_face(face_verts);
							// Note: Some face insertions may fail due to mesh topology, this is normal
						}
					}
				}

				// Check if we have any faces in the patch mesh
				if (patch_mesh->get_num_of_faces() == 0) {
					cs_err() << ">>> Error: No faces were successfully inserted for cycle " << cy_itr << std::endl;
					delete patch_mesh;
					continue;
				}
				
				// Update the patch mesh after construction
				patch_mesh->update_mesh();
				patch_mesh->update_normal();
				
				// Clean up non-manifold geometry before remeshing
				{
					std::vector<Geo::Point3<double>> mesh_pos;
					std::vector<std::vector<int>> mesh_faces;
					
					// Convert Mesh3D to Geo format
					for (int i = 0; i < patch_mesh->get_num_of_vertices(); i++) {
						HE_vert* v = patch_mesh->get_vertex(i);
						mesh_pos.push_back(Geo::Point3<double>(v->pos.x(), v->pos.y(), v->pos.z()));
					}
					
					for (int i = 0; i < patch_mesh->get_num_of_faces(); i++) {
						HE_face* f = patch_mesh->get_face(i);
						std::vector<int> face_verts;
						HE_edge* edge = f->edge;
						HE_edge* start_edge = edge;
						do {
							face_verts.push_back(edge->vert->id);
							edge = edge->next;
						} while (edge != start_edge);
						mesh_faces.push_back(face_verts);
					}
					
					int initial_verts = mesh_pos.size();
					int initial_faces = mesh_faces.size();
					
					// Clean non-manifold vertices by deleting isolated faces
Geo::CleanMesh<double>::SolveNonManifoldByDeletingFaces(mesh_pos, mesh_faces);
					
					int cleaned_verts = mesh_pos.size();
					int cleaned_faces = mesh_faces.size();
					
					if (initial_verts != cleaned_verts || initial_faces != cleaned_faces) {
						cs_log() << "Cleaned non-manifold geometry: " 
						          << (initial_verts - cleaned_verts) << " vertices removed, "
						          << (initial_faces - cleaned_faces) << " faces removed" << std::endl;
						
						// Rebuild the patch mesh from cleaned data
						delete patch_mesh;
						patch_mesh = new Mesh3D;
						
						// Add cleaned vertices
						for (size_t i = 0; i < mesh_pos.size(); i++) {
							Eigen::Vector3d v(mesh_pos[i].X, mesh_pos[i].Y, mesh_pos[i].Z);
							patch_mesh->insert_vertex(v);
						}
						
						// Add cleaned faces
						for (size_t i = 0; i < mesh_faces.size(); i++) {
							std::vector<HE_vert*> face_verts;
							for (int v_idx : mesh_faces[i]) {
								face_verts.push_back(patch_mesh->get_vertex(v_idx));
							}
							patch_mesh->insert_face(face_verts);
						}
						
						patch_mesh->update_mesh();
						patch_mesh->update_normal();
					}
				}
				
				cs_log() << "Remeshing patch " << cy_itr << " (" << patch_mesh->get_num_of_vertices() 
				          << " verts, " << patch_mesh->get_num_of_faces() << " faces)..." << std::flush;
				
				// Apply remeshing to the patch (replaces split_edge)
				remesh_patch(patch_mesh, cycle_edges, kRemeshEdgeLenMultiplier);

				cs_log() << " done (" << patch_mesh->get_num_of_vertices()
				          << " verts, " << patch_mesh->get_num_of_faces() << " faces)" << std::endl;

				// Insert vertices directly from patch_mesh using iterators
				Geex::Mesh3D::PTR_VERTEX_LIST patch_vertex_list = patch_mesh->get_vertices_list();
				std::vector<HE_vert*> vertex_mapping(patch_vertex_list->size());
				int vertex_index = 0;
				
				for (HE_vert* patch_vert : *patch_vertex_list) {
					HE_vert* new_vert = nmesh->insert_vertex(patch_vert->pos);
					vertex_mapping[vertex_index] = new_vert;
					vertex_index++;
				}
				
				// Insert faces directly from patch_mesh
				Geex::Mesh3D::PTR_FACE_LIST patch_face_list = patch_mesh->get_faces_list();
				for (HE_face* patch_face : *patch_face_list) {
					auto face_vertices = patch_face->get_vertices();
					std::vector<HE_vert*> face_vts;
					face_vts.reserve(face_vertices.size());
					
					for (HE_vert* patch_vert : face_vertices) {
						face_vts.push_back(vertex_mapping[patch_vert->id]);
					}
					
					HE_face* fh = nmesh->insert_face(face_vts);
					if (fh) {
						fh->curve_cycle = cur_cycle;
					}
				}
				
				// Store the patch mesh for cleanup later
				temp_patch_meshes.push_back(patch_mesh);
				
			} else {
				cs_err() << "Failed to generate surface for cycle " << cy_itr << std::endl;
				if (patch_mesh) {
					delete patch_mesh;
					patch_mesh = nullptr;
				}
			}
		}
		
		cs_log() << "Cycles and patches loaded." << std::endl;
		
		nmesh->update_mesh();
		cs_log() << "Finished loading data from structures." << std::endl;
		
		// Clean up temporary patch meshes safely
		// The mesh data was copied (not moved) into nmesh via insert_vertex/insert_face
		// We need to clear the patch mesh lists to avoid double-deletion
		for (Mesh3D* temp_mesh : temp_patch_meshes) {
			// Manually clear the lists without deleting objects
			// This prevents the destructor from trying to delete them
			if (temp_mesh->get_vertices_list()) {
				temp_mesh->get_vertices_list()->clear();
			}
			if (temp_mesh->get_edges_list()) {
				temp_mesh->get_edges_list()->clear();
			}
			if (temp_mesh->get_faces_list()) {
				temp_mesh->get_faces_list()->clear();
			}
			// Now safe to delete - destructor will see empty lists
			delete temp_mesh;
		}
		temp_patch_meshes.clear();
		
		match_patchboundary_curves(nmesh);
		computeCurveRMF(nmesh);
		
		// Fix inward-pointing RMF frames using geometric analysis
		fixInwardRMFFrames(nmesh);
		
		// Average overlapping edge RMF frames for consistency (skip feature edges)
		averageOverlappingRMFFrames(nmesh, feature_edges);
		
		// Check cycle planarity
		for (int cyitr = 0; cyitr < m_p_polygraph_->get_num_of_cycles(); cyitr++) {
			check_cycle_planarity(m_p_polygraph_->get_cycle(cyitr));
		}

	cs_log() << "All cycle planarity checks complete." << std::endl;
		
	return nmesh;
}

// Overloaded version that accepts user-supplied patch meshes
CS_INLINE Mesh3D* CurvNet::load_cycle_discovery_data_from_structures(
	const std::vector<CurveData>& curves_data,
	const std::vector<CycleData>& cycles_data,
	const std::vector<MeshData>& patches_data,
	const FeatureEdgeData& feature_edges) {
	
	// Validate that we have the same number of cycles and patches
	if (cycles_data.size() != patches_data.size()) {
		cs_err() << ">>> Error: Number of cycles (" << cycles_data.size() 
		          << ") does not match number of patches (" << patches_data.size() << ")" << std::endl;
		return nullptr;
	}
	
	// Load curves into CurvNet
	for (size_t c_itr = 0; c_itr < curves_data.size(); c_itr++) {
		const CurveData& curve_data = curves_data[c_itr];
		
		SplineCurv* new_curv = new SplineCurv(this);
		new_curv->id = c_itr;
		new_curv->degree = 1;
		new_curv->type = SPLINE_TYPE::BSPLINE;
		
		new_curv->para_range[0] = 0;
		new_curv->para_range[1] = 1;
		new_curv->knots.push_back(0);
		
		for (size_t pt_itr = 0; pt_itr < curve_data.points.size(); pt_itr++) {
			const Eigen::Vector3d& pt = curve_data.points[pt_itr];
			this->m_control_points.push_back(Vector4d(pt.x(), pt.y(), pt.z(), 0));
			new_curv->m_cpoints.push_back(this->m_control_points.size() - 1);
			new_curv->knots.push_back((double)pt_itr / (double)(curve_data.num_pts - 1));
		}
		new_curv->knots.push_back(1);
		
		this->m_curves.push_back(new_curv);
	}
	
	cs_log() << "Curves loaded into CurvNet." << std::endl;
	
	// Create combined mesh
	Mesh3D* nmesh = new Mesh3D;
	
	// Initialize polygraph
	if (!m_p_polygraph_) {
		m_p_polygraph_ = new PolyGraph;
	} else {
		delete m_p_polygraph_;
		m_p_polygraph_ = new PolyGraph;
	}
	
	// Load cycles and user-supplied patches
	for (size_t cy_itr = 0; cy_itr < cycles_data.size(); cy_itr++) {
		const CycleData& cycle_data = cycles_data[cy_itr];
		const MeshData& patch_data = patches_data[cy_itr];
		
		// Create cycle edges
		std::vector<PolyGraph_Edge*> cycle_edges;
		for (const auto& curve_info : cycle_data.curve_info) {
			int curv_id = curve_info[0];
			int curv_dir = curve_info[1];
			int instance = curve_info[2];
			
			cycle_edges.push_back(m_p_polygraph_->insert_edge(m_curves[curv_id], (curv_dir == 1)));
		}
		
		PolyGraph_Cycle* cur_cycle = m_p_polygraph_->insert_cycle(cycle_edges);
		
		// Load the user-supplied patch into the combined mesh
		cs_log() << "Loading user-supplied patch " << cy_itr << " with " 
		          << patch_data.num_vertices << " vertices and " 
		          << patch_data.num_faces << " faces" << std::endl;
		
		const unsigned int vid_offset = nmesh->get_num_of_vertices();
		
		// Insert vertices
		for (const auto& vertex : patch_data.vertices) {
			nmesh->insert_vertex(const_cast<Eigen::Vector3d&>(vertex));
		}
		
		// Insert faces
		for (const auto& face_vertices : patch_data.faces) {
			std::vector<HE_vert*> face_vts;
			face_vts.reserve(face_vertices.size());
			
			for (int vid : face_vertices) {
				face_vts.push_back(nmesh->get_vertex(vid + vid_offset));
			}
			
			HE_face* fh = nmesh->insert_face(face_vts);
			if (fh) {
				fh->curve_cycle = cur_cycle;
			}
		}
	}
	
	cs_log() << "Cycles and user-supplied patches loaded." << std::endl;
	
	nmesh->update_mesh();
	cs_log() << "Finished loading data from structures." << std::endl;
	
	match_patchboundary_curves(nmesh);
	computeCurveRMF(nmesh);
	
	// Fix inward-pointing RMF frames using geometric analysis
	fixInwardRMFFrames(nmesh);
	
	// Average overlapping edge RMF frames for consistency (skip feature edges)
	averageOverlappingRMFFrames(nmesh, feature_edges);
	
	// Check cycle planarity
	for (int cyitr = 0; cyitr < m_p_polygraph_->get_num_of_cycles(); cyitr++) {
		check_cycle_planarity(m_p_polygraph_->get_cycle(cyitr));
	}

	cs_log() << "All cycle planarity checks complete." << std::endl;
	
	return nmesh;
}

// Helper function to populate RMF frames from normal curve data
CS_INLINE void CurvNet::populateRMFFromNormalData(Mesh3D* patches, const NormalCurveData& normal_data) {
	cs_log() << "Populating RMF frames from normal curve data..." << std::endl;
	
	if (!patches || !m_p_polygraph_ || normal_data.empty()) {
		cs_err() << "Error: Invalid input to populateRMFFromNormalData" << std::endl;
		return;
	}
	
	// Validate that we have matching number of cycles
	if (normal_data.num_cycles != m_p_polygraph_->get_num_of_cycles()) {
		cs_err() << "Warning: Number of cycles in normal data (" << normal_data.num_cycles 
		          << ") doesn't match PolyGraph cycles (" << m_p_polygraph_->get_num_of_cycles() << ")" << std::endl;
	}
	
	int edges_populated = 0;
	int total_mesh_edges = 0;
	int curves_marked_crvt = 0;
	int curves_marked_trim = 0;
	
	// The normal curve file stores normals for directed edges in cycle traversal order
	// We walk through curves in cycle order, respecting direction, and sequentially consume normals
	
	// Iterate through each cycle
	for (int cy_idx = 0; cy_idx < std::min((int)normal_data.cycles.size(), m_p_polygraph_->get_num_of_cycles()); cy_idx++) {
		const CycleNormalData& cycle_normal_data = normal_data.cycles[cy_idx];
		PolyGraph_Cycle* cur_cycle = m_p_polygraph_->get_cycle(cy_idx);
		
		if (!cur_cycle || cycle_normal_data.edges.empty()) continue;
		
		// Walk through the cycle's PolyGraph_Edges in order
		std::vector<PolyGraph_Edge*> cycle_pg_edges;
		PolyGraph_Edge* pg_edge = cur_cycle->edge;
		do {
			cycle_pg_edges.push_back(pg_edge);
			pg_edge = pg_edge->next;
		} while (pg_edge != cur_cycle->edge);
		
		// Sequential consumption of normals matching the directed edge traversal
		int normal_edge_idx = 0;
		
		for (PolyGraph_Edge* pg_edge : cycle_pg_edges) {
			SplineCurv* curve = pg_edge->curve_;
			if (!curve) continue;
			
			// Get curve vertices - these are the control points that define the polyline
			int num_curve_vertices = curve->m_cpoints.size();
			int num_curve_edges = std::max(0, num_curve_vertices - 1);
			
			if (num_curve_edges == 0) continue;
			
			// Collect normals for this curve's directed edges from the normal data
			std::vector<Eigen::Vector3d> curve_normals;
			for (int i = 0; i < num_curve_edges && normal_edge_idx < (int)cycle_normal_data.edges.size(); i++) {
				curve_normals.push_back(cycle_normal_data.edges[normal_edge_idx].edge_normal.normalized());
				normal_edge_idx++;
			}
			
			if (curve_normals.empty()) continue;

			// Determine start/end normals in increasing curve parameter order.
			// If the cycle traverses the curve in reverse, the normal sequence is reversed
			// with respect to curve parameterization.
			Eigen::Vector3d n0 = pg_edge->order_ ? curve_normals.front() : curve_normals.back();
			Eigen::Vector3d n1 = pg_edge->order_ ? curve_normals.back() : curve_normals.front();
			
			// Find all mesh edges that belong to this PolyGraph_Edge (curve)
			std::vector<ParaEdge> mesh_edges_on_curve;
			
			for (int fitr = 0; fitr < patches->get_num_of_faces(); fitr++) {
				HE_face* fh = patches->get_face(fitr);
				if (fh->curve_cycle != cur_cycle) continue;
				
				HE_edge* eh = fh->edge;
				do {
					if (eh->cycle_curv == pg_edge) {
						total_mesh_edges++;
						Eigen::Vector3d midpt = eh->GetMidPoint();
						double para;
						curve->project_pt(midpt, &para);
						mesh_edges_on_curve.push_back(ParaEdge(para, eh, midpt));
					}
					eh = eh->next;
				} while (eh != fh->edge);
			}
			
			// Sort mesh edges by parametric position along the curve
			std::sort(mesh_edges_on_curve.begin(), mesh_edges_on_curve.end());

			// Use the same end-frame consistency idea as computeRMFonCurve to classify curve type.
			// Transport n0 along the curve using double-reflection and compare with n1.
			if (n0.norm() > 1e-12 && n1.norm() > 1e-12) {
				n0.normalize();
				n1.normalize();

				Eigen::Vector3d transported_r = n0;

				double u0 = 0.0;
				double u1 = 1.0;
				if (!mesh_edges_on_curve.empty()) {
					u0 = std::max(0.0, std::min(1.0, mesh_edges_on_curve.front().para));
					u1 = std::max(0.0, std::min(1.0, mesh_edges_on_curve.back().para));
				}
				if (u0 > u1) {
					std::swap(u0, u1);
					std::swap(n0, n1);
				}

				Eigen::Vector3d prev_pt = curve->get_value_at(u0);
				double prev_para = u0;

				for (const ParaEdge& pe : mesh_edges_on_curve) {
					double cur_para = std::max(u0, std::min(u1, pe.para));
					Eigen::Vector3d cur_pt = pe.mid_pt;

					Eigen::Vector3d v1 = cur_pt - prev_pt;
					double c1 = v1.dot(v1);
					if (c1 < 1e-16) {
						prev_pt = cur_pt;
						prev_para = cur_para;
						continue;
					}

					Eigen::Vector3d rl = transported_r - (2.0 / c1) * v1.dot(transported_r) * v1;

					Eigen::Vector3d fo, so;
					curve->get_drev(prev_para, fo);
					curve->get_drev(cur_para, so);
					Eigen::Vector3d tl = fo - (2.0 / c1) * v1.dot(fo) * v1;
					Eigen::Vector3d v2 = so - tl;
					double c2 = v2.dot(v2);
					if (c2 > 1e-16) {
						Eigen::Vector3d rll = rl - (2.0 / c2) * v2.dot(rl) * v2;
						if (rll.norm() > 1e-12) {
							transported_r = rll.normalized();
						}
					}

					prev_pt = cur_pt;
					prev_para = cur_para;
				}

				// Final transport to u1 endpoint
				Eigen::Vector3d end_pt = curve->get_value_at(u1);
				Eigen::Vector3d v1_end = end_pt - prev_pt;
				double c1_end = v1_end.dot(v1_end);
				if (c1_end > 1e-16) {
					Eigen::Vector3d rl = transported_r - (2.0 / c1_end) * v1_end.dot(transported_r) * v1_end;

					Eigen::Vector3d fo, so;
					curve->get_drev(prev_para, fo);
					curve->get_drev(u1, so);
					Eigen::Vector3d tl = fo - (2.0 / c1_end) * v1_end.dot(fo) * v1_end;
					Eigen::Vector3d v2 = so - tl;
					double c2 = v2.dot(v2);
					if (c2 > 1e-16) {
						Eigen::Vector3d rll = rl - (2.0 / c2) * v2.dot(rl) * v2;
						if (rll.norm() > 1e-12) {
							transported_r = rll.normalized();
						}
					}
				}

				double dotv = transported_r.dot(n1);
				dotv = std::max(-1.0, std::min(1.0, dotv));
				double theta = std::acos(dotv);
				const double RMF_CONSISTENCY_ANGLE = 10.0 / 180.0 * M_PI;

				if (theta < RMF_CONSISTENCY_ANGLE || theta > M_PI - RMF_CONSISTENCY_ANGLE) {
					pg_edge->curve_type = CURVE_TYPE::CRVT_LINE;
					curves_marked_crvt++;
				} else {
					pg_edge->curve_type = CURVE_TYPE::TRIM_LINE;
					curves_marked_trim++;
				}
			}

			// Assign normals to mesh edges using nearest-neighbor (no interpolation)
			// curve_normals[i] is the normal for the edge from vertex i to vertex i+1
			// The normals are already in the correct order from the normal file (matching cycle traversal)
			for (const ParaEdge& pe : mesh_edges_on_curve) {
				double t = pe.para; // Parametric position in [0, 1]
				
				// If curve is reversed (order_ == false), flip the parametric coordinate
				// This ensures we map to the correct normal in the sequence
				if (!pg_edge->order_) {
					t = 1.0 - t;
				}
				
				// Map parametric position to nearest normal index (round instead of interpolate)
				double normal_idx_float = t * (num_curve_edges - 1);
				int nearest_idx = std::max(0, std::min((int)std::round(normal_idx_float), (int)curve_normals.size() - 1));
				
				// Use the nearest normal directly - no interpolation
				Eigen::Vector3d normal = curve_normals[nearest_idx];
				
				pe.edge->tangent = normal;
				if (pe.edge->pair) {
					pe.edge->pair->tangent = normal;
				}
				edges_populated++;
			}
		}
	}
	
	cs_log() << "RMF population complete: " << edges_populated << " of " << total_mesh_edges 
	          << " mesh edges populated from normal data" << std::endl;
	cs_log() << "Curve type classification from supplied normals: "
	          << curves_marked_crvt << " CRVT_LINE, "
	          << curves_marked_trim << " TRIM_LINE" << std::endl;
}

// Overloaded version that accepts normal curve data for direct RMF population
CS_INLINE Mesh3D* CurvNet::load_cycle_discovery_data_from_structures(
	const std::vector<CurveData>& curves_data,
	const std::vector<CycleData>& cycles_data,
	const FeatureEdgeData& feature_edges,
	const NormalCurveData& normal_data) {
	
	// Load curves into CurvNet
	for (size_t c_itr = 0; c_itr < curves_data.size(); c_itr++) {
		const CurveData& curve_data = curves_data[c_itr];
		
		SplineCurv* new_curv = new SplineCurv(this);
		new_curv->id = c_itr;
		new_curv->degree = 1;
		new_curv->type = SPLINE_TYPE::BSPLINE;
		
		new_curv->para_range[0] = 0;
		new_curv->para_range[1] = 1;
		new_curv->knots.push_back(0);
		
		for (size_t pt_itr = 0; pt_itr < curve_data.points.size(); pt_itr++) {
			const Eigen::Vector3d& pt = curve_data.points[pt_itr];
			this->m_control_points.push_back(Vector4d(pt.x(), pt.y(), pt.z(), 0));
			new_curv->m_cpoints.push_back(this->m_control_points.size() - 1);
			new_curv->knots.push_back((double)pt_itr / (double)(curve_data.num_pts - 1));
		}
		new_curv->knots.push_back(1);
		
		this->m_curves.push_back(new_curv);
	}
	
	cs_log() << "Curves loaded into CurvNet." << std::endl;
	
	// Create combined mesh
	Mesh3D* nmesh = new Mesh3D;
	
	// Initialize polygraph
	if (!m_p_polygraph_) {
		m_p_polygraph_ = new PolyGraph;
	} else {
		delete m_p_polygraph_;
		m_p_polygraph_ = new PolyGraph;
	}
	
	// Load cycles and create patches using SurfaceMaker
	std::vector<Mesh3D*> temp_patch_meshes; // Store temporary meshes to delete later
	
	for (size_t cy_itr = 0; cy_itr < cycles_data.size(); cy_itr++) {
		const CycleData& cycle_data = cycles_data[cy_itr];
		
		// Create cycle edges
		std::vector<PolyGraph_Edge*> cycle_edges;
		for (const auto& curve_info : cycle_data.curve_info) {
			int curv_id = curve_info[0];
			int curv_dir = curve_info[1];
			int instance = curve_info[2];
			
			cycle_edges.push_back(m_p_polygraph_->insert_edge(m_curves[curv_id], (curv_dir == 1)));
		}
		
		PolyGraph_Cycle* cur_cycle = m_p_polygraph_->insert_cycle(cycle_edges);
		
		// Create patch mesh using SurfaceMaker
		cs_log() << "Creating surface patch for cycle " << cy_itr << "..." << std::endl;
		
		// Create a new mesh for this patch
		Mesh3D* patch_mesh = new Mesh3D;
		
		// Convert curve data to AML::double3 format expected by SurfaceMaker
		std::vector<std::vector<AML::double3>> cycle_curve_network;
		
		for (PolyGraph_Edge* edge : cycle_edges) {
			SplineCurv* curve = edge->curve_;
			std::vector<AML::double3> curve_points;
			
			// Sample points along the curve
			for (int cpoint_idx : curve->m_cpoints) {
				Vector4d& cp = this->m_control_points[cpoint_idx];
				AML::double3 aml_point;
				aml_point[0] = cp.x();
				aml_point[1] = cp.y(); 
				aml_point[2] = cp.z();
				curve_points.push_back(aml_point);
			}
			
			// Reverse the curve if needed based on order
			if (edge->order_) {
				std::reverse(curve_points.begin(), curve_points.end());
			}
			
			cycle_curve_network.push_back(curve_points);
		}
		
		// Create SurfaceMaker instance and generate surface
		SurfaceMaker surface_maker;
		surface_maker.OnOpenCurveNetwork(cycle_curve_network);
		
		// Get the generated triangle surface
		const GraphStructure::TriangleSurface& triangle_surface = surface_maker.GetTriangleSurface();
		
		if (!triangle_surface.empty()) {
			cs_log() << "Generated surface with " << triangle_surface.size() << " triangle cycles" << std::endl;
			
			// Extract points and faces for this specific patch (with proper vertex deduplication)
			std::vector<AML::double3> patch_points;
			std::vector<std::vector<unsigned>> patch_faces;
			
			// Process triangles
			for (const GraphStructure::TriangleCycle& triangle_cycle : triangle_surface) {
				for (unsigned j = 0; j < triangle_cycle.size(); j++) {
					const std::vector<AML::double3>& triangle = triangle_cycle[j];
					std::vector<unsigned> tf;
					
					for (unsigned k = 0; k < triangle.size(); k++) {
						const AML::double3& tp = triangle[k];
						unsigned pID;
						
						// Find if point already exists in this patch's point list
						auto it = std::find(patch_points.begin(), patch_points.end(), tp);
						if (it == patch_points.end()) {
							patch_points.push_back(tp);
							pID = patch_points.size() - 1;
						} else {
							pID = it - patch_points.begin();
						}
						tf.push_back(pID);
					}
					patch_faces.push_back(tf);
				}
			}
			
			cs_log() << "Patch " << cy_itr << ": " << patch_points.size() << " points, " << patch_faces.size() << " faces" << std::endl;
			
			// Insert vertices into the patch mesh
			std::vector<HE_vert*> patch_vertices;
			for (unsigned i = 0; i < patch_points.size(); i++) {
				const AML::double3& coord = patch_points[i];
				Eigen::Vector3d v(coord[0], coord[1], coord[2]);
				HE_vert* vert = patch_mesh->insert_vertex(v);
				patch_vertices.push_back(vert);
			}
			
			// Insert faces into the patch mesh
			for (unsigned i = 0; i < patch_faces.size(); i++) {
				const std::vector<unsigned>& tri = patch_faces[i];
				if (tri.size() >= 3) {
					std::vector<HE_vert*> face_verts;
					for (unsigned j = 0; j < 3; j++) {
						if (tri[j] < patch_vertices.size()) {
							face_verts.push_back(patch_vertices[tri[j]]);
						}
					}
					if (face_verts.size() == 3) {
						HE_face* fh = patch_mesh->insert_face(face_verts);
					}
				}
			}
			
			// Check if we have any faces in the patch mesh
			if (patch_mesh->get_num_of_faces() == 0) {
				cs_err() << ">>> Error: No faces were successfully inserted for cycle " << cy_itr << std::endl;
				delete patch_mesh;
				continue;
			}
			
			// Update the patch mesh after construction
			patch_mesh->update_mesh();
			patch_mesh->update_normal();
			
			// Clean up non-manifold geometry before remeshing
			{
				std::vector<Geo::Point3<double>> mesh_pos;
				std::vector<std::vector<int>> mesh_faces;
				
				// Convert Mesh3D to Geo format
				for (int i = 0; i < patch_mesh->get_num_of_vertices(); i++) {
					HE_vert* v = patch_mesh->get_vertex(i);
					mesh_pos.push_back(Geo::Point3<double>(v->pos.x(), v->pos.y(), v->pos.z()));
				}
				
				for (int i = 0; i < patch_mesh->get_num_of_faces(); i++) {
					HE_face* f = patch_mesh->get_face(i);
					std::vector<int> face_verts;
					HE_edge* edge = f->edge;
					HE_edge* start_edge = edge;
					do {
						face_verts.push_back(edge->vert->id);
						edge = edge->next;
					} while (edge != start_edge);
					mesh_faces.push_back(face_verts);
				}
				
				int initial_verts = mesh_pos.size();
				int initial_faces = mesh_faces.size();
				
				// Clean non-manifold vertices by deleting isolated faces
				Geo::CleanMesh<double>::SolveNonManifoldByDeletingFaces(mesh_pos, mesh_faces);
				
				int cleaned_verts = mesh_pos.size();
				int cleaned_faces = mesh_faces.size();
				
				if (initial_verts != cleaned_verts || initial_faces != cleaned_faces) {
					cs_log() << "Cleaned non-manifold geometry: " 
					          << (initial_verts - cleaned_verts) << " vertices removed, "
					          << (initial_faces - cleaned_faces) << " faces removed" << std::endl;
					
					// Rebuild the patch mesh from cleaned data
					delete patch_mesh;
					patch_mesh = new Mesh3D;
					
					// Add cleaned vertices
					for (size_t i = 0; i < mesh_pos.size(); i++) {
						Eigen::Vector3d v(mesh_pos[i].X, mesh_pos[i].Y, mesh_pos[i].Z);
						patch_mesh->insert_vertex(v);
					}
					
					// Add cleaned faces
					for (size_t i = 0; i < mesh_faces.size(); i++) {
						std::vector<HE_vert*> face_verts;
						for (int v_idx : mesh_faces[i]) {
							face_verts.push_back(patch_mesh->get_vertex(v_idx));
						}
						patch_mesh->insert_face(face_verts);
					}
					
					patch_mesh->update_mesh();
					patch_mesh->update_normal();
				}
			}
			
			cs_log() << "Remeshing patch " << cy_itr << " (" << patch_mesh->get_num_of_vertices() 
			          << " verts, " << patch_mesh->get_num_of_faces() << " faces)..." << std::flush;
			
			// Apply remeshing to the patch
			remesh_patch(patch_mesh, cycle_edges, kRemeshEdgeLenMultiplier);

			// ALTERNATIVE: Use Laplacian smoothing instead of remeshing
			laplacian_smooth_patch(patch_mesh, 20, 0.8);  // 10 iterations, lambda=0.5
			
			cs_log() << " done (" << patch_mesh->get_num_of_vertices()
			          << " verts, " << patch_mesh->get_num_of_faces() << " faces)" << std::endl;

			// Insert vertices directly from patch_mesh using iterators
			Geex::Mesh3D::PTR_VERTEX_LIST patch_vertex_list = patch_mesh->get_vertices_list();
			std::vector<HE_vert*> vertex_mapping(patch_vertex_list->size());
			int vertex_index = 0;
			
			for (HE_vert* patch_vert : *patch_vertex_list) {
				HE_vert* new_vert = nmesh->insert_vertex(patch_vert->pos);
				vertex_mapping[vertex_index] = new_vert;
				vertex_index++;
			}
			
			// Insert faces directly from patch_mesh
			Geex::Mesh3D::PTR_FACE_LIST patch_face_list = patch_mesh->get_faces_list();
			for (HE_face* patch_face : *patch_face_list) {
				auto face_vertices = patch_face->get_vertices();
				std::vector<HE_vert*> face_vts;
				face_vts.reserve(face_vertices.size());
				
				for (HE_vert* patch_vert : face_vertices) {
					face_vts.push_back(vertex_mapping[patch_vert->id]);
				}
				
				HE_face* fh = nmesh->insert_face(face_vts);
				if (fh) {
					fh->curve_cycle = cur_cycle;
				}
			}
			
			// Store the patch mesh for cleanup later
			temp_patch_meshes.push_back(patch_mesh);
			
		} else {
			cs_err() << "Failed to generate surface for cycle " << cy_itr << std::endl;
			if (patch_mesh) {
				delete patch_mesh;
				patch_mesh = nullptr;
			}
		}
	}
	
	cs_log() << "Cycles and patches loaded." << std::endl;
	
	nmesh->update_mesh();
	cs_log() << "Finished loading data from structures." << std::endl;
	
	// Clean up temporary patch meshes safely
	for (Mesh3D* temp_mesh : temp_patch_meshes) {
		if (temp_mesh->get_vertices_list()) {
			temp_mesh->get_vertices_list()->clear();
		}
		if (temp_mesh->get_edges_list()) {
			temp_mesh->get_edges_list()->clear();
		}
		if (temp_mesh->get_faces_list()) {
			temp_mesh->get_faces_list()->clear();
		}
		delete temp_mesh;
	}
	temp_patch_meshes.clear();
	
	match_patchboundary_curves(nmesh);
	
	// Check if normal data is actually provided
	if (!normal_data.empty()) {
		// Populate RMF frames from the normal data instead of computing them
		cs_log() << "Using provided normal data to populate RMF frames..." << std::endl;
		populateRMFFromNormalData(nmesh, normal_data);
		cs_log() << "RMF frames populated from normal data (skipped computeCurveRMF, fixInwardRMFFrames, averageOverlappingRMFFrames)" << std::endl;
	} else {
		// Normal data is empty, fallback to computing RMF frames
		cs_log() << "Normal data is empty, computing RMF frames automatically..." << std::endl;
		computeCurveRMF(nmesh);
		fixInwardRMFFrames(nmesh);
		averageOverlappingRMFFrames(nmesh, feature_edges);
		cs_log() << "RMF frames computed and fixed." << std::endl;
	}
	
	// Check cycle planarity
	for (int cyitr = 0; cyitr < m_p_polygraph_->get_num_of_cycles(); cyitr++) {
		check_cycle_planarity(m_p_polygraph_->get_cycle(cyitr));
	}

	cs_log() << "All cycle planarity checks complete." << std::endl;
		
	return nmesh;
}

// Overload with user-supplied patches AND normal curve data
CS_INLINE Mesh3D* CurvNet::load_cycle_discovery_data_from_structures(
	const std::vector<CurveData>& curves_data,
	const std::vector<CycleData>& cycles_data,
	const std::vector<MeshData>& patches_data,
	const FeatureEdgeData& feature_edges,
	const NormalCurveData& normal_data) {
	
	// Call the version without normal_data first to load the mesh with user patches
	Mesh3D* nmesh = load_cycle_discovery_data_from_structures(curves_data, cycles_data, patches_data, feature_edges);
	
	if (!nmesh) {
		return nullptr;
	}
	
	// Now populate RMF frames from the normal data if provided
	if (!normal_data.empty()) {
		cs_log() << "Using provided normal data to populate RMF frames on user-supplied patches..." << std::endl;
		populateRMFFromNormalData(nmesh, normal_data);
		cs_log() << "RMF frames populated from normal data." << std::endl;
	} else {
		cs_log() << "Warning: Normal data is empty, RMF frames already computed by previous call." << std::endl;
	}
	
	return nmesh;
}

CS_INLINE Mesh3D* CurvNet::load_cycle_discovery_data(std::string filename){		data_file_name = filename;

		//load curves.
		std::string curv_filename = filename + "_scaled.curve";
		std::ifstream curv_file (curv_filename);

		if (curv_file.is_open())
		{
			int curv_num = 0;
			// reads the first line which is the number of curves and assigns it to curv_num
			curv_file >> curv_num;
			for (int c_itr = 0; c_itr < curv_num; c_itr++)
			{
				// the next 3 lines are the number of points, whether the curve is closed, and the capacity (not used here)
				int num_pts, closed, capacity;
				curv_file >> num_pts >> closed >> capacity;

				// default to linear B-spline curve
				SplineCurv* new_curv = new SplineCurv(this);
				new_curv->id = c_itr;
				new_curv->degree = 1;
				new_curv->type = SPLINE_TYPE::BSPLINE;
				
				new_curv->para_range[0] = 0;
				new_curv->para_range[1] = 1;
				// knots determine which control points influence which part of the curve
				new_curv->knots.push_back(0);
				for (int pt_itr = 0; pt_itr < num_pts; pt_itr++)
				{
					// next three numbers are the x,y,z coordinates of each control point
					double x,y,z;
					curv_file >> x >> y >> z;

					this->m_control_points.push_back(Vector4d(x, y, z, 0));
					new_curv->m_cpoints.push_back(this->m_control_points.size()-1);

					new_curv->knots.push_back((double)pt_itr/(double)(num_pts-1));
				}
				new_curv->knots.push_back(1);

				this->m_curves.push_back(new_curv);
			}
			cs_log() << ".curve file loaded." << std::endl;
		}
		else
		{
			cs_err() << ">>> Could not load " << curv_filename << std::endl;
			return NULL;
		}

		curv_file.close();

		//load cycles.
		Mesh3D* nmesh = new Mesh3D;

		bool found_marked_cycle_file = true;
		std::string cycle_filename = filename + "_marked.cycle";

		std::ifstream cycle_file (cycle_filename);

		if (!cycle_file.is_open())
		{
			found_marked_cycle_file = false;

			cycle_filename = filename + ".cycle";
			cycle_file = std::ifstream(cycle_filename);
			cs_log() << ">>> Did not find marked cycle file." << std::endl;
		}

		if (cycle_file.is_open())
		{
			if (!m_p_polygraph_)
			{
				m_p_polygraph_ = new PolyGraph;
			}
			else
			{
				delete m_p_polygraph_;
				m_p_polygraph_ = new PolyGraph;
			}

			// first line is the number of cycles
			int num_cycles = 0;
			cycle_file >> num_cycles;
			for (int cy_itr = 0; cy_itr < num_cycles; cy_itr++)
			{
				// second line is the number of curves in the current cycle
				int ccurv_num;
				cycle_file >> ccurv_num;

				std::vector<PolyGraph_Edge*> cycle_edges;

				for (int curv_itr = 0; curv_itr < ccurv_num; curv_itr++)
				{
					// next three numbers are curve ID, curve direction (1 or 2), and instance (used only if marked cycle file)
					int curv_id, curv_dir, instance;
					cycle_file >> curv_id >> curv_dir >> instance;

					cycle_edges.push_back(m_p_polygraph_->insert_edge(m_curves[curv_id], (curv_dir==1)));

					if (found_marked_cycle_file)
					{
						cycle_edges.back()->curve_type = (CURVE_TYPE)instance;
					}
				}

				PolyGraph_Cycle* cur_cycle = m_p_polygraph_->insert_cycle(cycle_edges);
				

				//load cycle patch here.
				{
					std::string mesh_filename = filename + std::to_string(cy_itr) + ".off";
	
					std::ifstream mesh_file (mesh_filename);

					if (mesh_file.is_open())
					{
						// first line is "OFF"
						std::string offmark;
						mesh_file >> offmark;

						// second line is number of vertices, number of faces, and an unused number
						int vt_num, face_num, unused;
						mesh_file >> vt_num >> face_num >> unused;

						const unsigned int vid_offset = nmesh->get_num_of_vertices();

						for (int vitr = 0; vitr < vt_num; vitr++)
						{
							// next three numbers are the x,y,z coordinates of each vertex
							double x,y,z;
							mesh_file >> x >> y >> z;
							Eigen::Vector3d temp_vertex(x,y,z);
							nmesh->insert_vertex(temp_vertex);
						}

						for (int fitr = 0; fitr < face_num; fitr++)
						{
							// first number is the degree of the face, followed by the vertex IDs that make up the face
							int face_degree;
							mesh_file >> face_degree;

							std::vector<HE_vert*> face_vts (face_degree, NULL);

							for (int fv_itr = 0; fv_itr < face_degree; fv_itr++)
							{
								int vid;
								mesh_file >> vid;
								face_vts[fv_itr] = nmesh->get_vertex(vid +  vid_offset);
							}

							HE_face* fh = nmesh->insert_face(face_vts);
							if (!fh)
							{
								cs_err() << "Failed inserting facet" << std::endl;
							}
							else
							{
								// associates the current cycle with the face
								fh->curve_cycle = cur_cycle;
							}
						}

						cs_log() << "loaded patch "<<cy_itr<<std::endl;
					}
					else
					{
						cs_err() << ">>> failed loading patch "<<cy_itr<<std::endl;
					}

					mesh_file.close();
				}
			}

			cs_log() << ".cycle file loaded." << std::endl;
		}

		cycle_file.close();

		nmesh->update_mesh();
		cs_log() << "Finished loading data from cyclediscovery." << std::endl;

		match_patchboundary_curves(nmesh);
		computeCurveRMF(nmesh);
		
		// Fix inward-pointing RMF frames using geometric analysis
		fixInwardRMFFrames(nmesh);
		
		// Average overlapping edge RMF frames for consistency
		averageOverlappingRMFFrames(nmesh);

		//check cycle planarity
		for (int cyitr = 0; cyitr < m_p_polygraph_->get_num_of_cycles(); cyitr ++)
		{
			check_cycle_planarity(m_p_polygraph_->get_cycle(cyitr));
		}

		return nmesh;

	}

	CS_INLINE void CurvNet::computeRMFonCurve(PolyGraph_Edge* pg_edge, SplineCurv* curve, double  u [], Eigen::Vector3d  n [], std::vector<ParaEdge>& v_pts){

		std::vector<Eigen::Vector3d> m_sample_interp_rgt;
		std::vector<Eigen::Vector3d> m_sample_interp_sgt;
		std::vector<double> m_sample_dist;

		double u0=u[0], u1=u[1];
		Eigen::Vector3d n0=n[0], n1=n[1];
	
		if (u0>u1){
			std::swap(u0,u1);
			std::swap(n0,n1);
		}

		std::vector<Eigen::Vector3d> t;
		Eigen::Vector3d tan;
		m_sample_interp_rgt.push_back(n0);
		curve->get_drev(u0, tan);
		m_sample_interp_sgt.push_back(n0.cross(tan.normalized()));
		t.push_back(tan);
		m_sample_dist.push_back(0);

		//RMF without condition
		Eigen::Vector3d fo, so;
		Eigen::Vector3d prev_pt = curve->get_value_at(u0);
		double prev_para = u0;
		for(int i=0;i<v_pts.size();i++){
			Eigen::Vector3d cur_pt = v_pts[i].mid_pt;
			double cur_para = v_pts[i].para;

			Eigen::Vector3d v1 = cur_pt - prev_pt;
			double c1 = v1.dot(v1);
			Eigen::Vector3d rl = m_sample_interp_rgt[i] - (2/c1)*v1.dot(m_sample_interp_rgt[i])*v1;
			curve->get_drev(prev_para, fo);
			curve->get_drev(cur_para, so);
			t.push_back(so);
			Eigen::Vector3d tl = fo - (2/c1)*v1.dot(fo)*v1;
			Eigen::Vector3d v2 = so - tl;
			double c2 = v2.dot(v2);
			Eigen::Vector3d rll = (rl - (2/c2) * v2.dot(rl) * v2).normalized();
			Eigen::Vector3d sll = rll.cross(so).normalized();

			m_sample_interp_rgt.push_back(rll);
			m_sample_interp_sgt.push_back(sll);

			m_sample_dist.push_back((cur_pt-prev_pt).norm()+m_sample_dist.back());
			prev_pt = cur_pt;
			prev_para = cur_para;
		}
		{
			Eigen::Vector3d cur_pt = curve->get_value_at(u1);
			double cur_para = u1;

			Eigen::Vector3d v1 = cur_pt - prev_pt;
			double c1 = v1.dot(v1); 
			Eigen::Vector3d rl = m_sample_interp_rgt.back() - (2/c1) * v1.dot(m_sample_interp_rgt.back()) * v1;
			curve->get_drev(prev_para, fo);
			curve->get_drev(cur_para, so);
			t.push_back(so);
			Eigen::Vector3d tl = fo - (2/c1) * v1.dot(fo) * v1;
			Eigen::Vector3d v2 = so - tl;
			double c2 = v2.dot(v2);
			Eigen::Vector3d rll = (rl - (2/c2) * v2.dot(rl) * v2).normalized();  // Eigen's normalize method
    		Eigen::Vector3d sll = rll.cross(so).normalized();
			m_sample_interp_rgt.push_back(rll);
			m_sample_interp_sgt.push_back(sll);

			m_sample_dist.push_back((cur_pt-prev_pt).norm()+m_sample_dist.back());
		}

		//Adjust the direction
		double theta = std::acos(m_sample_interp_rgt.back().dot(n1));

		const double RMF_CONSISTENCY_ANGLE = 10.0/180.0*M_PI;
		if (theta < RMF_CONSISTENCY_ANGLE || theta > M_PI - RMF_CONSISTENCY_ANGLE)
		{
			pg_edge->curve_type = CURVE_TYPE::CRVT_LINE;
		}
		else
		{
			pg_edge->curve_type = CURVE_TYPE::TRIM_LINE;
		}
		
		bool mark = false;
		if(theta > M_PI/2.0){
			theta = M_PI - theta;
			Eigen::Vector3d tmp = n1*(-1);
			if(m_sample_interp_sgt.back().dot(tmp) < 0)
				mark = true;
		}
		else{
			if(m_sample_interp_sgt.back().dot(n1) < 0)
				mark = true;
		}
		theta /= m_sample_dist.back();

		if(theta > 0.0)
		{
			for(int i=0;i<m_sample_interp_rgt.size();i++){
				Eigen::Vector3d tmp1;
				if(!mark)
					tmp1 = m_sample_interp_rgt[i]*cos(m_sample_dist[i]*theta) + m_sample_interp_sgt[i]*sin(m_sample_dist[i]*theta);
				else
					tmp1 = m_sample_interp_rgt[i]*cos(m_sample_dist[i]*theta) - m_sample_interp_sgt[i]*sin(m_sample_dist[i]*theta);
				Eigen::Vector3d tmp2 = tmp1.cross(t[i]);
				m_sample_interp_rgt[i] = tmp1.normalized();
				m_sample_interp_sgt[i] = tmp2.normalized();
			}
		}

		for (int eitr = 0; eitr < v_pts.size(); eitr++)
		{
			v_pts[eitr].edge->tangent = v_pts[eitr].edge->pair->tangent = m_sample_interp_rgt[eitr+1];
		}
	}



	CS_INLINE void IntersectionNormal(SplineCurv* curve, SplineCurv* adj_curve, double& para, Eigen::Vector3d& normal, 
							const std::vector<ParaEdge>* mesh_edges = nullptr){
		double c_para[2], adjc_para[2];
		Eigen::Vector3d c_ept[2], adjc_ept[2];
		Eigen::Vector3d c_drv[2], adjc_drv[2];

		for (int i = 0; i < 2; i++)
		{
			c_para[i] = curve->para_range[i];
			c_ept[i] = curve->get_value_at(c_para[i]);
			curve->get_drev(c_para[i], c_drv[i]);

			adjc_para[i] = adj_curve->para_range[i];
			adjc_ept[i] = adj_curve->get_value_at(adjc_para[i]);
			adj_curve->get_drev(adjc_para[i], adjc_drv[i]);
		}

		int c_idx=-1, adjc_idx=-1;
		double min_dist=1e10;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if ((c_ept[i]-adjc_ept[j]).norm() < min_dist)
				{
					min_dist = (c_ept[i]-adjc_ept[j]).norm();
					c_idx=i;
					adjc_idx=j;
				}
			}
		}

		para = c_para[c_idx];
		normal = c_drv[c_idx].cross(adjc_drv[adjc_idx]);

		if (normal.norm() < 1.e-7)
		{
			normal = Eigen::Vector3d(0.,0.,0.);
		}
		else
		{
			normal = normal.normalized();
			
				// Validate against mesh face normals if available
				if (mesh_edges && !mesh_edges->empty())
				{
					// Find mesh edges near this corner (parameter value)
					Eigen::Vector3d corner_pt = c_ept[c_idx];
					Eigen::Vector3d avg_face_normal(0, 0, 0);
					int valid_face_count = 0;
					
					// Collect face normals from nearby mesh edges (within small parameter range)
					const double param_tolerance = 0.1; // Check edges within 10% of parameter range
					for (const auto& pedge : *mesh_edges)
					{
						if (std::abs(pedge.para - para) < param_tolerance)
						{
							// Get face normal from this edge
							HE_edge* eh = pedge.edge;
							if (eh->face && eh->face->normal.norm() > 1e-7)
							{
								avg_face_normal += eh->face->normal.normalized();
								valid_face_count++;
							}
							if (eh->pair && eh->pair->face && eh->pair->face->normal.norm() > 1e-7)
							{
								avg_face_normal += eh->pair->face->normal.normalized();
								valid_face_count++;
							}
						}
					}
					
					// If we found nearby faces, validate and correct the normal direction
					if (valid_face_count > 0)
					{
						avg_face_normal.normalize();
						
						// If corner normal points opposite to average face normal, flip it
						if (normal.dot(avg_face_normal) < 0)
						{
							normal = -normal;
						}
					}
				}
			}
		}

	CS_INLINE void CurvNet::computeCornerNormals(PolyGraph_Edge* pg_edge, Eigen::Vector3d end_normals[2], double end_para[2],
							   const std::vector<ParaEdge>* mesh_edges){
		SplineCurv* curve = pg_edge->curve_;
		
		SplineCurv* adj_curve[2];
		adj_curve[0] = pg_edge->prev->curve_;
		adj_curve[1] = pg_edge->next->curve_;

		IntersectionNormal(curve, adj_curve[0], end_para[0], end_normals[0], mesh_edges);
		IntersectionNormal(curve, adj_curve[1], end_para[1], end_normals[1], mesh_edges);

		if (end_para[0] > end_para[1])
		{
			std::swap(end_para[0], end_para[1]);
			std::swap(end_normals[0], end_normals[1]);
		}
	}

	CS_INLINE void CurvNet::computeCurveRMF(Mesh3D* patches){
		cs_log() << "Compute RMF normals on curves..." << std::endl;
		
		int total_faces = 0;
		int faces_with_curve_cycle = 0;
		int edges_with_cycle_curv = 0;
		int rmf_computations = 0;
		
		if (m_p_polygraph_&&patches)
		{
			for (int cy_itr = 0; cy_itr < m_p_polygraph_->get_num_of_cycles(); cy_itr++)
			{
				PolyGraph_Cycle* cur_cy = m_p_polygraph_->get_cycle(cy_itr);

				

				std::unordered_map<PolyGraph_Edge*, std::vector<ParaEdge>> map_pgedge_meshedge;

				for (int fitr = 0; fitr < patches->get_num_of_faces(); fitr++)
				{
					total_faces++;
					
					HE_face* fh=patches->get_face(fitr);
					if (fh->curve_cycle != cur_cy)
					{
						continue;
					}
					faces_with_curve_cycle++;
					
					HE_edge* eh=fh->edge;
					do
					{
						if (eh->cycle_curv)
						{
							edges_with_cycle_curv++;
							PolyGraph_Edge* pg_edge = eh->cycle_curv;
							Eigen::Vector3d midpt = eh->GetMidPoint();
							double para;
							pg_edge->curve_->project_pt(midpt,&para);
							map_pgedge_meshedge[pg_edge].push_back(ParaEdge(para,eh,midpt));
						}
						eh=eh->next;
					} while (eh!=fh->edge);
				}

				for (std::unordered_map<PolyGraph_Edge*, std::vector<ParaEdge>>::iterator cyc_pg_iter = map_pgedge_meshedge.begin();
					cyc_pg_iter != map_pgedge_meshedge.end(); cyc_pg_iter++)
				{
					std::sort(cyc_pg_iter->second.begin(), cyc_pg_iter->second.end());
					PolyGraph_Edge* pg_edge = cyc_pg_iter->first;

					Eigen::Vector3d end_normals[2];
					double end_para[2];
					// Pass mesh edge information for corner normal validation
					computeCornerNormals(pg_edge, end_normals, end_para, &(cyc_pg_iter->second));

					if (end_normals[0].norm()>1e-5 && end_normals[1].norm()>1e-5)
					{
						rmf_computations++;
						computeRMFonCurve(pg_edge, pg_edge->curve_, end_para, end_normals, map_pgedge_meshedge[pg_edge]);
					}
					
				}
			}
		}
		
		cs_log() << "RMF Debug - RMF computations performed: " << rmf_computations << std::endl;
	}


	CS_INLINE void CurvNet::fixInwardRMFFrames(Mesh3D* mesh) {
    for (int e = 0; e < mesh->get_num_of_edges(); ++e) {
        HE_edge* eh = mesh->get_edge(e);
        if (!eh || !eh->face) continue;
        if (eh->tangent.norm() < 1e-8) continue; // RMF undefined

        Eigen::Vector3d n_face = computeFaceNormal(eh->face);
        Eigen::Vector3d n_rmf  = eh->tangent.normalized();

        double d = n_face.dot(n_rmf);
        if (d < 0.0) {
            eh->tangent = -eh->tangent;
            if (eh->pair)
                eh->pair->tangent = -eh->pair->tangent;
        }
    }
}

	CS_INLINE void CurvNet::averageOverlappingRMFFrames(Mesh3D* mesh, const FeatureEdgeData& feature_edges) {
		if (!mesh) return;
		
		int overlapping_groups_processed = 0;
		int total_frames_averaged = 0;
		
		// Helper function to check if an edge is a feature edge
		auto is_feature_edge = [&feature_edges](HE_edge* edge, double tolerance = 1e-3) -> bool {
			if (feature_edges.empty()) return false;
			
			Eigen::Vector3d v0 = edge->vert->pos;
			Eigen::Vector3d v1 = edge->pair->vert->pos;
			
			for (const auto& fe : feature_edges.edges) {
				// Check if edge matches feature edge in either direction
				double dist_forward = (v0 - fe.start).norm() + (v1 - fe.end).norm();
				double dist_backward = (v0 - fe.end).norm() + (v1 - fe.start).norm();
				
				if (dist_forward < tolerance || dist_backward < tolerance) {
					return true;
				}
			}
			return false;
		};
		
		// Create spatial index for quick proximity searches
		std::vector<HE_edge*> boundary_edges_with_rmf;
		
		// Collect all boundary edges with valid RMF data (excluding feature edges)
		for (int eitr = 0; eitr < mesh->get_num_of_edges(); eitr++) {
			HE_edge* eh = mesh->get_edge(eitr);
			if (!eh || !eh->pair) continue;
			
			// Skip feature edges - they should not have their RMF frames averaged
			if (is_feature_edge(eh)) continue;
			
			// Check if it's a boundary edge with RMF data
			bool is_boundary = (!eh->face || !eh->pair->face);
			if (is_boundary && eh->tangent.norm() > 1e-12) {
				boundary_edges_with_rmf.push_back(eh);
			}
		}

		// Process each edge to find overlapping groups
		std::vector<bool> processed(boundary_edges_with_rmf.size(), false);
		
		for (size_t i = 0; i < boundary_edges_with_rmf.size(); i++) {
			if (processed[i]) continue;
			
			HE_edge* primary_edge = boundary_edges_with_rmf[i];
			std::vector<HE_edge*> overlapping_group;
			std::vector<size_t> group_indices;
			
			overlapping_group.push_back(primary_edge);
			group_indices.push_back(i);
			processed[i] = true;
			
			// Find all edges that overlap spatially with the primary edge
			for (size_t j = i + 1; j < boundary_edges_with_rmf.size(); j++) {
				if (processed[j]) continue;
				
				HE_edge* candidate = boundary_edges_with_rmf[j];
				
				if (areEdgesOverlapping(primary_edge, candidate)) {
					overlapping_group.push_back(candidate);
					group_indices.push_back(j);
					processed[j] = true;
				}
			}
			
			// If we found overlapping edges, average their RMF frames
			if (overlapping_group.size() > 1) {
				averageRMFFramesInGroup(overlapping_group);

				// Verify all directions are now within threshold
				Eigen::Vector3d reference_dir = overlapping_group[0]->tangent.normalized();
				bool all_consistent = true;
				for (size_t k = 1; k < overlapping_group.size(); k++) {
					Eigen::Vector3d current_dir = overlapping_group[k]->tangent.normalized();
					double alignment = std::abs(current_dir.dot(reference_dir));
					if (alignment < 0.999) { // Allow small numerical errors
						all_consistent = false;
					}
				}

				overlapping_groups_processed++;
				total_frames_averaged += overlapping_group.size();
			}
		}
	}
	
	CS_INLINE bool CurvNet::areEdgesOverlapping(HE_edge* edge1, HE_edge* edge2) {
		if (!edge1 || !edge2 || !edge1->pair || !edge2->pair) return false;
		
		// Get edge endpoints
		Eigen::Vector3d e1_p1 = edge1->vert->pos;
		Eigen::Vector3d e1_p2 = edge1->pair->vert->pos;
		Eigen::Vector3d e2_p1 = edge2->vert->pos;
		Eigen::Vector3d e2_p2 = edge2->pair->vert->pos;
		
		const double position_tolerance = 1e-6; // Spatial tolerance for "same position"
		
		// Check if endpoints are very close (indicating overlapping edges)
		bool p1_matches_p1 = (e1_p1 - e2_p1).norm() < position_tolerance;
		bool p1_matches_p2 = (e1_p1 - e2_p2).norm() < position_tolerance;
		bool p2_matches_p1 = (e1_p2 - e2_p1).norm() < position_tolerance;
		bool p2_matches_p2 = (e1_p2 - e2_p2).norm() < position_tolerance;
		
		// Edges overlap if both endpoints match (in either order)
		bool same_direction = p1_matches_p1 && p2_matches_p2;
		bool opposite_direction = p1_matches_p2 && p2_matches_p1;
		
		if (same_direction || opposite_direction) {
			return true;
		}
		
		// Alternative check: edges with midpoints very close and similar directions
		Eigen::Vector3d e1_midpoint = 0.5 * (e1_p1 + e1_p2);
		Eigen::Vector3d e2_midpoint = 0.5 * (e2_p1 + e2_p2);
		Eigen::Vector3d e1_direction = (e1_p2 - e1_p1).normalized();
		Eigen::Vector3d e2_direction = (e2_p2 - e2_p1).normalized();
		
		double midpoint_distance = (e1_midpoint - e2_midpoint).norm();
		double direction_alignment = std::abs(e1_direction.dot(e2_direction));
		
		// Edges overlap if midpoints are very close and directions are aligned
		return (midpoint_distance < position_tolerance && direction_alignment > 0.9);
	}
	
	CS_INLINE void CurvNet::averageRMFFramesInGroup(const std::vector<HE_edge*>& edges) {
		if (edges.empty()) return;
		
		// Compute average RMF direction ensuring consistent orientation
		Eigen::Vector3d average_direction(0, 0, 0);
		
		// Use first edge as reference direction
		Eigen::Vector3d reference_direction = edges[0]->tangent.normalized();
		average_direction += reference_direction;
		
		// Add other edges, flipping if necessary to maintain consistent orientation
		for (size_t i = 1; i < edges.size(); i++) {
			Eigen::Vector3d current_direction = edges[i]->tangent.normalized();
			
			// Check alignment with reference direction
			double alignment = current_direction.dot(reference_direction);
			
			if (alignment < 0) {
				// Flip direction to align with reference
				current_direction = -current_direction;
			}
			
			average_direction += current_direction;
		}
		
		// Normalize the average
		if (average_direction.norm() > 1e-12) {
			average_direction.normalize();
			
			// Apply averaged direction to all edges in the group
			for (HE_edge* edge : edges) {
				// Determine if we need to flip the edge to match the average
				Eigen::Vector3d current_normalized = edge->tangent.normalized();
				double alignment = current_normalized.dot(average_direction);
				
				if (alignment < 0) {
					// Flip to match average direction
					edge->tangent = -edge->tangent.norm() * average_direction;
				} else {
					// Keep same magnitude but use average direction
					edge->tangent = edge->tangent.norm() * average_direction;
				}
				
				// Update pair edge if it exists
				if (edge->pair) {
					edge->pair->tangent = edge->tangent;
				}
			}
		}
	}
	
	CS_INLINE Eigen::Vector3d CurvNet::computeFaceNormal(HE_face* face) {
		if (!face || !face->edge) return Eigen::Vector3d(0, 0, 0);
		
		// Get three consecutive vertices of the face
		HE_edge* e1 = face->edge;
		HE_edge* e2 = e1->next;
		HE_edge* e3 = e2->next;
		
		if (!e1 || !e2 || !e3) return Eigen::Vector3d(0, 0, 0);
		
		Eigen::Vector3d v1 = e1->vert->pos;
		Eigen::Vector3d v2 = e2->vert->pos;
		Eigen::Vector3d v3 = e3->vert->pos;
		
		// Compute face normal using cross product
		Eigen::Vector3d edge1 = v2 - v1;
		Eigen::Vector3d edge2 = v3 - v1;
		Eigen::Vector3d normal = edge1.cross(edge2);
		
		if (normal.norm() > 1e-12) {
			return normal.normalized();
		}
		
		return Eigen::Vector3d(0, 0, 0);
	}
	

	CS_INLINE void CurvNet::check_cycle_planarity(PolyGraph_Cycle* cycle){
		std::vector<Eigen::Vector3d> cycle_points;

		Eigen::Vector3d centroid (0.,0.,0.);

		PolyGraph_Edge* pg_edge = cycle->edge;
		do
		{
			SplineCurv* curve = pg_edge->curve_;
			for (auto cpt_idx : curve->m_cpoints)
			{
				Vector4d cpt = m_control_points[cpt_idx];
				cycle_points.push_back(Eigen::Vector3d(cpt.x(), cpt.y(), cpt.z()));
				centroid += cycle_points.back();
			}

			pg_edge = pg_edge->next;

		} while (pg_edge != cycle->edge);

		centroid /= cycle_points.size();

		Eigen::Matrix3d normal_mat = Eigen::Matrix3d::Zero();
		for (auto cpt : cycle_points)
		{
			Eigen::Vector3d diff = cpt - centroid;
			normal_mat += diff * diff.transpose();
		}

		normal_mat /= (double)cycle_points.size();

		// normal_mat is symmetric positive-semidefinite (it's a covariance matrix),
		// so its eigenvalues are equivalent to the squared singular values used by
		// the original TNT/JAMA SVD-based check. SelfAdjointEigenSolver returns
		// eigenvalues in ascending order.
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(normal_mat);
		Eigen::Vector3d sing_vals = eigensolver.eigenvalues();

		if (std::abs(sing_vals[0]/sing_vals[2]) < 1e-2)
		{
			cycle->is_planar = true;
		}
		else
		{
			cycle->is_planar = false;
		}
	}

	CS_INLINE void CurvNet::match_patchboundary_curves(Mesh3D* patches){
		int boundary_edges_processed = 0;
		int edges_matched = 0;
		
		if (patches)
		{

			for (int eitr = 0; eitr < patches->get_num_of_edges(); eitr++)
			{
				HE_edge* eh = patches->get_edge(eitr);
				if (!eh || !eh->face)
				{
					boundary_edges_processed++;
					
					// Check if we have valid pair and face
					if (!eh || !eh->pair || !eh->pair->face) continue;
					
					HE_face* facet = eh->pair->face;
					PolyGraph_Cycle* patch_cycle = facet->curve_cycle;

					// Check if patch_cycle is valid
					if (!patch_cycle || !patch_cycle->edge) continue;

					PolyGraph_Edge* cycle_edge = patch_cycle->edge;

					//project eh->midpoint to cycle_edge's curve.
					//find the overlapping one.
					const Eigen::Vector3d mid_pt = eh->GetMidPoint();
					PolyGraph_Edge* min_dist_edge = NULL;
					double min_dist = 1e10;
					PolyGraph_Edge* start_edge = cycle_edge;
					int safety_counter = 0;
					const int max_iterations = 1000; // Safety limit
					
					do
					{
						// Safety checks for cycle_edge and its curve
						if (!cycle_edge || !cycle_edge->curve_) break;
						
						Eigen::Vector3d prj_pt = cycle_edge->curve_->project_pt(mid_pt);

						double dist = (mid_pt-prj_pt).norm();

						if (dist < min_dist)
						{
							min_dist = dist;
							min_dist_edge = cycle_edge;
						}

						cycle_edge = cycle_edge->next;
						safety_counter++;
						
						// Safety check to prevent infinite loops
						if (safety_counter > max_iterations) {
							break;
						}
					} while (cycle_edge && cycle_edge != start_edge);

					if (min_dist_edge) {
						eh->cycle_curv = eh->pair->cycle_curv = min_dist_edge;
						edges_matched++;
					}
				}
			}
		}
	}

	CS_INLINE PolyGraph::PolyGraph()
	{
	}

	CS_INLINE PolyGraph::~PolyGraph()
	{
		for (int vitr = 0; vitr < vertices.size(); vitr++)
		{
			delete vertices[vitr];
		}

		for (int eitr = 0; eitr < edges.size(); eitr++)
		{
			delete edges[eitr];
		}

		for (int citr = 0; citr < cycles.size(); citr++)
		{
			delete cycles[citr];
		}
	}

	CS_INLINE PolyGraph_Cycle* PolyGraph::insert_cycle(const std::vector<PolyGraph_Edge*> cedges){
		if (cedges.size() > 0)
		{
			PolyGraph_Cycle* ncycle = new PolyGraph_Cycle;
			this->cycles.push_back(ncycle);
			ncycle->id = this->cycles.size() - 1;

			for (int eitr = 0; eitr < cedges.size(); eitr++)
			{
				PolyGraph_Edge* eh = cedges[eitr];
				eh->cycle_ = ncycle;
				eh->next = cedges[(eitr+1)%cedges.size()];
				eh->next->prev = eh;
			}

			ncycle->degree = cedges.size();
			ncycle->edge = cedges.front();

			return ncycle;
		}
		return NULL;
	}

}