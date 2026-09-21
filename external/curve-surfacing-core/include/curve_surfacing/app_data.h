/*
*  _____ _____ ________  _
* /  __//  __//  __/\  \//
* | |  _|  \  |  \   \  / 
* | |_//|  /_ |  /_  /  \ 
* \____\\____\\____\/__/\\
*
* Graphics Environment for EXperimentations.
*  Copyright (C) 2006 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: 
*
*     ALICE Project - INRIA
*     INRIA Lorraine, 
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX 
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs. 
*/

#ifndef __APP_DATASTRUCTURE__
#define __APP_DATASTRUCTURE__

#include "Mesh3D.h"
#include "CurvNet.h"
#include "app_loader.h"

namespace Geex {

    //______________________________________________________________________________________
    class APP_DataStructs {

    public:
        APP_DataStructs() ;
        virtual ~APP_DataStructs() ;

        void save(const std::string& filename) ;
        void load(const std::string& filename) ;
		void savea(const std::string& filename) ;
		void save_mesh(const std::string& filename) ;

		void load_curvnet(const std::string& filename);
		void save_polygraph(const std::string& filename);
		void load_polygraph(const std::string& filename);

		void load_curv_mesh(const std::string& filename);
		void load_curv_mesh_mod(std::vector<CurveData> curv_data, std::vector<CycleData> cycle_data);
		void convert_to_curvefile(const std::string& filename);
		void convert_to_objfile(const std::string& filename);
		void save_cycles();
		void save_curve_types();
		void save_curve_types(const std::string& file_prefix);
		void save_g1_smooth_marking();
		void apply_g1_smooth_marking();
		void pairup_g1(PolyGraph_Cycle* cycle_0, PolyGraph_Cycle* cycle_1);
		void auto_mark_g1_smoothness();

		void load_reference_mesh(const std::string& filename);
		void load_reference_mesh_mod(MeshData ref_mesh_data);

		//reverse the orientation of facets connected to seed_fh.
		void reverse_component_orientation(Mesh3D* mesh, HE_face* seed_fh);

        void get_bbox(
            double& x_min, double& y_min, double& z_min,
            double& x_max, double& y_max, double& z_max
            ) ;

		int& get_curve_sample_num() {return m_int_curve_samples_;}
		int get_curve_type() {
			if (m_p_current_cycle_curve_)
			{
				return m_p_current_cycle_curve_->curve_type;
			}
			return CURVE_TYPE::CRVT_LINE;
		}
		int& get_vert_type() {return m_vertex_type_;}

		bool& switch_align_bdry_curv() {return m_b_align_bdry_curv_;}

    // protected:

		bool m_b_bdry_fixed_;

		bool m_b_align_bdry_curv_;

		//switch for cross field or conjugate field computation.
		bool b_cross_conjugate_switch_; //true for cross, false for conjugate.

		Mesh3D* m_init_mesh_;
		Mesh3D* m_reference_mesh_;

		CurvNet* m_p_curve_net_;

		int m_int_curve_samples_;
		
		PolyGraph_Cycle* m_p_current_cycle_;
		PolyGraph_Edge* m_p_current_cycle_curve_;
		int m_vertex_type_;

		PolyGraph_Cycle* m_cyclepair_[2];

		HE_vert* m_p_current_vtx_;
    } ;

    //______________________________________________________________________________________

}

#ifndef CURVE_SURFACING_STATIC
#  include "app_data.cpp"
#endif

#endif