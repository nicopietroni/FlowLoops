#pragma once
#include "curve_surfacing_inline.h"

#include "Mesh3D.h"
#include <set>
#include "constraint.h"
#include "least_squares_solver.h"

namespace Geex
{


	CS_INLINE Mesh3D::Mesh3D(void)
	{
		vertices_list = NULL;
		faces_list = NULL;
		edges_list = NULL;

		xmax = ymax = zmax = (double)1.0;
		xmin = ymin = zmin = (double)-1.0;

		m_closed = true;
		m_quad = false;
		m_tri = false;
		m_hex = false;
		m_pentagon = false;

		m_num_components = 0;
		m_num_boundaries = 0;
		m_genus = 0;
		m_encounter_non_manifold = false;


	}

	CS_INLINE Mesh3D::~Mesh3D(void)
	{
		clear_data();

		
	}
	CS_INLINE void Mesh3D::clear_data()
	{
		clear_faces();
		clear_edges();
		clear_vertices();
	}
	CS_INLINE void Mesh3D::clear_vertices()
	{
		if (vertices_list==NULL) 
			return;

		size_t deleted = 0;
		for (VERTEX_ITER viter = vertices_list->begin(); viter != vertices_list->end(); viter++) {
			if (*viter) {  // Safety check
				delete *viter;
				++deleted;
			}
		}

		delete vertices_list;
		vertices_list = NULL;
	}
	CS_INLINE void Mesh3D::clear_edges()
	{
		m_edgemap.clear();
		
		if (edges_list==NULL) 
			return;

		size_t deleted = 0;
		for (EDGE_ITER eiter = edges_list->begin(); eiter != edges_list->end(); eiter++) {
			if (*eiter) {  // Safety check
				delete *eiter;
				++deleted;
			}
		}

		delete edges_list;
		edges_list = NULL;
	}
	CS_INLINE void Mesh3D::clear_faces()
	{
		if (faces_list==NULL) 
			return;

		size_t deleted = 0;
		for (FACE_ITER fiter = faces_list->begin(); fiter!=faces_list->end(); fiter++) {
			if (*fiter) {  // Safety check
				delete *fiter;
				++deleted;
			}
		}

		delete faces_list;
		faces_list = NULL;
	}
	CS_INLINE HE_vert* Mesh3D::insert_vertex(Eigen::Vector3d &v)
	{
		HE_vert* hv = new HE_vert(v);
		if (vertices_list==NULL) 
			vertices_list = new VERTEX_LIST;

		hv->id = (ptrdiff_t)vertices_list->size();
		vertices_list->push_back(hv);

		return hv;
	}
	CS_INLINE HE_face* Mesh3D::insert_face(VERTEX_LIST & vec_hv, std::vector<ptrdiff_t> * texture, std::vector<ptrdiff_t> * normal)
	{
		int vsize = (int)vec_hv.size();
		if (vsize < 3) 
			return NULL;

		bool b_find = false;
		for (int i = 0; i < vsize; i++)
		{
			HE_edge* c_he = m_edgemap[PAIR_VERTEX(vec_hv[i], vec_hv[(i+1)%vsize])];
			if (c_he && c_he->face)
			{
				
				b_find = true;
				break;
			}
		}

		if (b_find)
		{
			m_encounter_non_manifold = true;

			for (int i = 0; i < vsize; i++)
			{
				HE_edge* c_he = m_edgemap[PAIR_VERTEX(vec_hv[(i+1)%vsize], vec_hv[i])];
				if (c_he && c_he->face)
				{
					return NULL; 
				}
			}
			std::reverse(vec_hv.begin(), vec_hv.end());
			if (texture)
			{
				std::reverse(texture->begin(), texture->end());
			}
			if (normal)
			{
				std::reverse(normal->begin(), normal->end());
			}
		}


		if (faces_list==NULL) 
			faces_list = new FACE_LIST;

		HE_face* hf = new HE_face;
		hf->valence = vsize;
		VERTEX_ITER viter = vec_hv.begin();
		VERTEX_ITER nviter = vec_hv.begin();
		nviter++;

		HE_edge *he1, *he2;
		std::vector<HE_edge* > v_he;
		int i;
		for (i=0; i<vsize-1; i++) 
		{
			he1 = insert_edge( *viter, *nviter);
			he2 = insert_edge( *nviter, *viter);

			if (hf->edge==NULL) 
				hf->edge = he1;

			he1->face = hf;
			he1->pair = he2;
			he2->pair = he1;
			v_he.push_back(he1);
			viter++, nviter++;
		}

		nviter = vec_hv.begin();

		he1 = insert_edge(*viter, *nviter);
		he2 = insert_edge(*nviter , *viter);
		he1->face = hf;
		if (hf->edge==NULL) 
			hf->edge = he1;

		he1->pair = he2;
		he2->pair = he1;
		v_he.push_back(he1);

		for (i=0; i<vsize-1; i++) 
		{
			v_he[i]->next = v_he[i+1];
			v_he[i+1]->prev = v_he[i];
		}
		v_he[i]->next = v_he[0];
		v_he[0]->prev = v_he[i];

		hf->id = (int)faces_list->size();
		faces_list->push_back(hf);


		if (texture)
		{
			hf->texture_indices.resize(texture->size());
			std::copy(texture->begin(), texture->end(), hf->texture_indices.begin());
		}
		if (normal)
		{
			hf->normal_indices.resize(normal->size());
			std::copy(normal->begin(), normal->end(), hf->normal_indices.begin());
		}

		return hf;
	}
	
	CS_INLINE HE_edge* Mesh3D::insert_edge(HE_vert* vstart, HE_vert* vend)
	{
		if (vstart == NULL || vend == NULL)
		{
			return NULL;
		}

		if (edges_list==NULL) 
			edges_list = new EDGE_LIST;

		if( m_edgemap[PAIR_VERTEX(vstart, vend)] != NULL )
			return m_edgemap[PAIR_VERTEX(vstart, vend)];

		HE_edge* he = new HE_edge;

		he->vert = vend;
		he->vert->degree++;
		vstart->edge = he;
		m_edgemap[PAIR_VERTEX(vstart, vend)] = he;

		he->id = (int)edges_list->size();
		edges_list->push_back(he);

		return he;
	}

	CS_INLINE HE_edge* Mesh3D::pg_insert_edge(HE_vert* vstart, HE_vert* vend){
		if (vstart == NULL || vend == NULL)
		{
			return NULL;
		}

		if (edges_list==NULL) 
			edges_list = new EDGE_LIST;


		HE_edge* he = new HE_edge;

		he->vert = vend;
		he->vert->degree++;
		vstart->pg_edges.push_back(he); //hpan.
		vstart->edge = he;

		he->id = (int)edges_list->size();
		edges_list->push_back(he);

		return he;
	}

	CS_INLINE HE_face* Mesh3D::pg_insert_face(std::vector<HE_edge*>& edges){
		if (faces_list==NULL) 
			faces_list = new FACE_LIST;

		HE_face* hf = new HE_face;
		hf->valence = edges.size();
	
		for(int i=0; i<edges.size(); i++){
			edges[i]->face = hf;

			HE_edge* eh0 = edges[i], *eh1 = edges[(i+1)%edges.size()];
			if(eh0->vert == eh1->pair->vert)
			{
				eh0->next = eh1;
				eh1->prev = eh0;
			}
			else if (eh0->pair->vert == eh1->vert)
			{
				eh0->prev = eh1;
				eh1->next = eh0;
			}
			else{
			}
		}

		hf->edge = edges[0];
		hf->id = (int)faces_list->size();
		faces_list->push_back(hf);
		compute_perface_normal(hf);

		return hf;
	}

	CS_INLINE void Mesh3D::init_vef_list(){
		if(!vertices_list) vertices_list = new VERTEX_LIST;
		if(!edges_list) edges_list = new EDGE_LIST;
		if(!faces_list) faces_list = new FACE_LIST;
	}

	
	CS_INLINE void Mesh3D::set_nextedge_for_border_vertices()
	{
		if (!is_valid()) 
			return;

		EDGE_ITER eiter = edges_list->begin();

		for (; eiter!=edges_list->end(); eiter++) 
		{
			if ((*eiter)->next==NULL&& (*eiter)->face==NULL)
				(*eiter)->pair->vert->edge = *eiter;

		}

		for (eiter = edges_list->begin(); eiter!=edges_list->end(); eiter++) 
		{
			if ( (*eiter)->next ==NULL ) 
			{
				HE_vert* hv = (*eiter)->vert;
				if (hv->edge != (*eiter)->pair)
				{
					(*eiter)->next = hv->edge;
					hv->edge->prev = *eiter;
				}
			}
		}
	}
	
	CS_INLINE bool Mesh3D::is_on_boundary(HE_vert* hv)
	{
		HE_edge* edge = hv->edge; 
		do {
			if (edge==NULL||edge->pair->face==NULL||edge->face==NULL) 
				return true;

			edge = edge->pair->next; 

		} while (edge != hv->edge);

		return false;
	}
	
	CS_INLINE bool Mesh3D::is_on_boundary(HE_face* hf)
	{
		HE_edge* edge = hf->edge; 

		do {
			if (is_on_boundary(edge)) 
				return true;
			edge = edge->next;
		} while (edge != hf->edge);

		return false;
	}
	
	CS_INLINE bool Mesh3D::is_on_boundary(HE_edge* he)
	{
		if(he->face==NULL||he->pair->face==NULL)
			return true;
		return false;
	}
	
	CS_INLINE void Mesh3D::check_closed()
	{
		m_closed = true; 
		VERTEX_ITER viter = vertices_list->begin();
		for (;viter!=vertices_list->end(); viter++) 
		{
			if (is_on_boundary(*viter)) {
				m_closed = false;
				return;
			}
		}
	}
	
	CS_INLINE void Mesh3D::check_meshtype()
	{
		if (get_num_of_faces() == 0)
		{
			m_tri = false;
			m_quad = false;
			m_hex = false;
			return;
		}

		FACE_ITER fiter = faces_list->begin();
		m_tri = true;
		m_quad = true;
		m_hex = true;
		m_pentagon = true;
		for (;fiter!=faces_list->end(); fiter++)
		{
			int d = (*fiter)->valence;
			if (d != 3)
			{
				m_tri = false;
			}
			if (d != 4)
			{
				m_quad = false;
			}
			if (d != 5)
			{
				m_pentagon = false;
			}
			if (d !=6)
			{
				m_hex = false;
			}
		}
	}
	
	CS_INLINE void Mesh3D::update_mesh()
	{
		if (!is_valid())
			return;

		remove_hanged_vertices();
		set_nextedge_for_border_vertices();
		check_closed();
		check_meshtype();
		m_edgemap.clear();
		update_normal();
		compute_boundingbox();
		compute_genus();

	}
	
	CS_INLINE bool Mesh3D::load_off(const char* fins)
	{
		std::ifstream fin(fins);

		try{

			clear_data();

			int vsize, fsize, esize;

			std::string head;

			fin>>head;

			if (head == "OFF") {
			}
			else
				return false;

			fin>>vsize>>fsize>>esize;

			double x, y, z;
			for (int i = 0; i<vsize; i++) {
				fin>>x>>y>>z;
				Eigen::Vector3d nvv(x, y, z) ;
				insert_vertex(nvv);
			}

			for (int i = 0; i<fsize; i++) {

				VERTEX_LIST v_list;
				int valence;
				fin>>valence;

				for (int j=0; j<valence; j++) {
					int id;
					fin>>id;
					HE_vert* hv = get_vertex(id);

					bool findit = false;
					for (int i = 0; i <(int) v_list.size(); i++)
					{
						if (hv == v_list[i])
						{
							findit = true;
							break;
						}
					}
					if (findit == false && hv != NULL)
					{
						v_list.push_back(hv);
					}

				}

				if ((int)v_list.size() >= 3)
				{
					insert_face(v_list);
				}
			}

			update_mesh();
		}
		catch (...) {
			clear_data();
			xmax = ymax = zmax = (double)1.0;
			xmin = ymin = zmin = (double)-1.0;
			fin.close();
			return false;
		}
		fin.close();
		return is_valid();
	}
	
	CS_INLINE void Mesh3D::write_off(const char* fouts)
	{

		std::ofstream fout(fouts);
		fout.precision(16);
		fout<<"OFF\n";
		fout<<(int)vertices_list->size()<<" "<<(int)faces_list->size()<<" "<<(int)edges_list->size()/2<<"\n";

		VERTEX_ITER viter = vertices_list->begin();
		for (;viter!=vertices_list->end(); viter++) {

			fout<<std::scientific<<(*viter)->pos<<"\n";
		}


		FACE_ITER fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) 
		{
			fout<<(*fiter)->valence;

			HE_edge* edge = (*fiter)->edge; 

			do {
				fout<<" "<<edge->pair->vert->id;
				edge = edge->next;

			} while (edge != (*fiter)->edge);
			fout<<"\n";
		}

		fout.close();
	}
	
	CS_INLINE bool Mesh3D::load_obj(const char* fins)
	{
		FILE *m_pFile = fopen(fins, "r");

		if (!m_pFile)
		{
			return false;
		}

		char *tok;

		try
		{
			clear_data();
			fseek(m_pFile, 0, SEEK_SET);
			char pLine[512];
			while(fgets(pLine, 512, m_pFile))
			{
				if(pLine[0] == 'v' && pLine[1] == ' ')
				{
					Eigen::Vector3d nvv;
					tok = strtok(pLine," ");
					for (int i=0; i<3; i++)
					{
						tok = strtok(NULL," ");
						// Malformed line (fewer than 3 coordinates): treat the
						// missing ones as 0 instead of dereferencing NULL.
						nvv[i] = tok ? atof(tok) : 0.0;
					}
					Eigen::Vector3d temp_vertex(nvv.x(),nvv.y(),nvv.z());
					insert_vertex(temp_vertex);
				}
				else if(pLine[0] == 'v' && pLine[1] == 't')
				{
				}
				else if(pLine[0] == 'v' && pLine[1] == 'n')
				{
				}
			}


			fseek(m_pFile, 0, SEEK_SET);


			while(fgets(pLine, 512, m_pFile))
			{
				char *pTmp = pLine;
				if(pTmp[0] == 'f')
				{
					VERTEX_LIST s_faceid;
					std::vector<ptrdiff_t> normal_ind, texture_ind;
					tok = strtok(pLine," ");
					while ((tok = strtok(NULL," ")) != NULL)
					{
						const std::string token(tok);   // v, v/t, v/t/n, or v//n
						size_t len = token.size();
						int start_pos = 0;
						int pos[2] = {-1,-1};
						for (int k = 0; k < (int)len && start_pos < 2; k++)
						{
							if (token[k] == '/')
							{
								pos[start_pos] = k;
								start_pos++;
							}
						}
						size_t end_pos = len;
						if (start_pos == 1)
						{
							end_pos = pos[0];
							int id = (int)strtol(token.c_str() + pos[0] + 1, NULL, 10) - 1;
							texture_ind.push_back(id);
						}
						else if (start_pos == 2)
						{
							end_pos = pos[0];
							if (pos[0]+1!=pos[1])
							{
								std::string mstr = token.substr(pos[0]+1, pos[1]-pos[0]-1);
								int id = (int)strtol(mstr.c_str(), NULL, 10) - 1;
								texture_ind.push_back(id);
							}
							int id2 = (int)strtol(token.c_str() + pos[1] + 1, NULL, 10) - 1;
							normal_ind.push_back(id2);

						}
						std::string mstr = token.substr(0, end_pos);
						int id = (int)strtol(mstr.c_str(), NULL, 10) - 1;
						// Out-of-range index in a malformed file must not read
						// past the vertex array.
						HE_vert* hv = (id >= 0 && id < get_num_of_vertices())
						              ? get_vertex(id) : NULL;
						bool findit = false;
						for (int i = 0; i <(int) s_faceid.size(); i++)
						{
							if (hv == s_faceid[i])	//remove the redundant vertex id if it exists
							{
								findit = true;
								break;
							}
						}
						if (findit == false && hv != NULL)
						{
							s_faceid.push_back(hv);
						}
					}
					if (s_faceid.size() >= 3)
					{
						insert_face(s_faceid, texture_ind.size()==s_faceid.size()?&texture_ind:0, normal_ind.size() == s_faceid.size()?&normal_ind:0);
					}
				}
			}
			update_mesh();
		}
		catch (...)
		{
			clear_data();
			xmax = ymax = zmax = (double)1.0;
			xmin = ymin = zmin = (double)-1.0;

			fclose(m_pFile);
			return false;
		}

		fclose(m_pFile);
		return is_valid();
	}
	
	CS_INLINE bool Mesh3D::load_from_mesh_data(const std::vector<Eigen::Vector3d>& vertices, 
	                                  const std::vector<std::vector<int>>& faces)
	{
		try
		{
			clear_data();
			
			for (size_t i = 0; i < vertices.size(); i++) {
				Eigen::Vector3d vertex_pos = vertices[i];
				insert_vertex(vertex_pos);
			}
			
			for (size_t i = 0; i < faces.size(); i++) {
				const std::vector<int>& face_indices = faces[i];
				VERTEX_LIST face_vertices;
				face_vertices.reserve(face_indices.size());
				
				for (int vertex_idx : face_indices) {
					if (vertex_idx >= 0 && vertex_idx < (int)vertices.size()) {
						HE_vert* hv = get_vertex(vertex_idx);
						if (hv != NULL) {
							face_vertices.push_back(hv);
						}
					}
				}
				
				if (face_vertices.size() >= 3) {  // Only add valid faces
					insert_face(face_vertices);
				}
			}
			
			update_mesh();
		}
		catch (...)
		{
			clear_data();
			xmax = ymax = zmax = (double)1.0;
			xmin = ymin = zmin = (double)-1.0;
			return false;
		}

		return is_valid();
	}
	
	
	CS_INLINE void Mesh3D::write_obj(const char* fouts)
	{
		std::ofstream fout(fouts);

		fout<<"g object\n";
		fout.precision(16);
		VERTEX_ITER viter = vertices_list->begin();
		for (;viter!=vertices_list->end(); viter++) {

			fout<<"v "<< std::scientific <<(*viter)->pos<<"\n";
		}


		FACE_ITER fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) {

			fout<<"f";

			HE_edge* edge = (*fiter)->edge; 
			unsigned int count = 0;
			do {
				fout<<" "<<edge->pair->vert->id+1;
				if ((*fiter)->has_texture_map())
				{
					fout<<"/" <<(*fiter)->texture_indices[count]+1;
				}
				else
					fout <<"/";

				if ((*fiter)->has_normal_map())
				{
					fout<<"/" <<(*fiter)->normal_indices[count]+1;
				}
				else
					fout<<"/" <<edge->pair->vert->id+1;


				edge = edge->next;

			} while (edge != (*fiter)->edge);
			fout<<"\n";
		}
		fout.close();
	}
	
	CS_INLINE void Mesh3D::write_vtk(const char* fouts)
	{
		std::ofstream fout(fouts);
		fout.precision(16);
		fout<<"# vtk DataFile Version 3.0\n"
			<<"Polygonal Mesh\n"
			<<"ASCII\n"
			<<"DATASET POLYDATA\n"
			<<"POINTS " << vertices_list->size() <<" double\n";

		VERTEX_ITER viter = vertices_list->begin();
		for (;viter!=vertices_list->end(); viter++) {
			fout<<std::scientific<<(*viter)->pos<<"\n";
		}

		FACE_ITER fiter = faces_list->begin();
		size_t count = 0;
		for (;fiter!=faces_list->end(); fiter++) 
		{
			count += (*fiter)->valence + 1;
		}

		fout <<"POLYGONS " << faces_list->size() << " " << count << "\n";


		fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) 
		{
			fout<<(*fiter)->valence;

			HE_edge* edge = (*fiter)->edge; 

			do {
				fout<<" "<<edge->pair->vert->id;
				edge = edge->next;

			} while (edge != (*fiter)->edge);
			fout<<"\n";
		}

		fout.close();
	}
	
	CS_INLINE void Mesh3D::compute_boundingbox()
	{
		if (vertices_list == NULL || vertices_list->size() < 3) 
			return;

		xmax = ymax = zmax = (double)-10e10;
		xmin = ymin = zmin = (double)10e10;


		VERTEX_ITER viter = vertices_list->begin();
		for (; viter!=vertices_list->end(); viter++) 
		{
			xmin = (*viter)->pos[0]<xmin? (*viter)->pos[0]:xmin;
			ymin = (*viter)->pos[1]<ymin? (*viter)->pos[1]:ymin;
			zmin = (*viter)->pos[2]<zmin? (*viter)->pos[2]:zmin;
			xmax = (*viter)->pos[0]>xmax? (*viter)->pos[0]:xmax;
			ymax = (*viter)->pos[1]>ymax? (*viter)->pos[1]:ymax;
			zmax = (*viter)->pos[2]>zmax? (*viter)->pos[2]:zmax;
		}
	}
	
	CS_INLINE void Mesh3D::compute_faces_list_normal()
	{
		for (FACE_ITER fiter = faces_list->begin(); fiter!=faces_list->end(); fiter++) 
		{
			compute_perface_normal(*fiter);
		}
	}
	
	CS_INLINE void Mesh3D::compute_perface_normal(HE_face* hf)
	{
		size_t i = 0;
		HE_edge* pedge = hf->edge;
		HE_edge* nedge = hf->edge->next;

		hf->normal = Eigen::Vector3d(0, 0, 0);
		for (i=0; i<hf->valence; i++) 
		{
			HE_vert* p = pedge->vert;
			HE_vert* c = pedge->next->vert;
			HE_vert* n = nedge->next->vert;
			Eigen::Vector3d pc, cn;
			pc = c->pos - p->pos; 
			cn = n->pos - c->pos; 

			hf->normal += pc.cross(cn);
			pedge = nedge;
			nedge = nedge->next;
			if (hf->valence==3) 
			{
				break;
			}
		}

		{
			hf->normal = (hf->normal).normalized();
		}
		
	}
	
	CS_INLINE void Mesh3D::compute_vertices_list_normal()
	{
		VERTEX_ITER viter = vertices_list->begin();

		for (; viter!=vertices_list->end(); viter++) 
		{
			compute_pervertex_normal(*viter);
		}
	}
	
	CS_INLINE void Mesh3D::compute_pervertex_normal(HE_vert* hv)
	{

		HE_edge* edge = hv->edge; 

		if (edge==NULL) {
			hv->normal = Eigen::Vector3d(0, 0, 0);
			return;
		}
		hv->normal = Eigen::Vector3d(0, 0, 0);

		do {
			if (edge->face!=NULL)
			{
				if (edge->face->get_area() > 1e-10) //check degenerate facets. hpan.
				{
					hv->normal += edge->face->normal;
				}
				
			}
			edge = edge->pair->next; 

		} while (edge && edge != hv->edge);

		hv->normal = (hv->normal).normalized();
	}
	
	CS_INLINE void Mesh3D::update_normal(bool onlyupdate_facenormal)
	{
		compute_faces_list_normal();
		if (onlyupdate_facenormal == false)
		{
			compute_vertices_list_normal();
		}

	}
	
	CS_INLINE void Mesh3D::copy_edge_id(ptrdiff_t edgeid, ptrdiff_t startid, ptrdiff_t endid, Mesh3D* mesh)
	{
		mesh->m_edgemap[PAIR_VERTEX( mesh->get_vertex(startid), mesh->get_vertex(endid))]->id = edgeid;
	}


	CS_INLINE Mesh3D* Mesh3D::make_copy_simple(){
		Mesh3D* new_mesh = new Mesh3D;

		VERTEX_ITER viter = vertices_list->begin();
		for (; viter!=vertices_list->end(); viter++)
		{
			HE_vert* vh = new_mesh->insert_vertex((*viter)->pos);

			if ((*viter)->fixed)
			{
				vh->fixed = true;
			}
			vh->m_vtype = (*viter)->m_vtype;
		}

		FACE_ITER fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) 
		{
			HE_face* hf = *fiter;
			HE_edge* edge = hf->edge;
			VERTEX_LIST mvlist;
			do {
				mvlist.push_back(new_mesh->get_vertex(edge->pair->vert->id));
				edge = edge->next;
			} while(edge!=hf->edge);

			HE_face* rfh = new_mesh->insert_face(mvlist, hf->texture_indices.empty()?0:&hf->texture_indices, hf->normal_indices.empty()?0:&hf->normal_indices);
		}

		new_mesh->update_mesh();

		PTR_EDGE_LIST edges_list = new_mesh->get_edges_list();
		for (int i = 0; i < edges_list->size(); i++)
		{
			(*edges_list)[i]->cycle_curv = (*edges_list)[i]->cycle_curv;
		}


		return new_mesh;
	}

	
	CS_INLINE Mesh3D* Mesh3D::make_copy()
	{
		Mesh3D* new_mesh = new Mesh3D;

		VERTEX_ITER viter = vertices_list->begin();
		for (; viter!=vertices_list->end(); viter++)
		{
			HE_vert* vh = new_mesh->insert_vertex((*viter)->pos);

			if ((*viter)->fixed)
			{
				vh->fixed = true;
			}
			vh->m_vtype = (*viter)->m_vtype;
		}

		FACE_ITER fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) 
		{
			HE_face* hf = *fiter;
			HE_edge* edge = hf->edge;
			VERTEX_LIST mvlist;
			do {
				mvlist.push_back(new_mesh->get_vertex(edge->pair->vert->id));
				edge = edge->next;
			} while(edge!=hf->edge);

			HE_face* rfh = new_mesh->insert_face(mvlist, hf->texture_indices.empty()?0:&hf->texture_indices, hf->normal_indices.empty()?0:&hf->normal_indices);

			rfh->curve_cycle = hf->curve_cycle;
		}

		for (EDGE_ITER eiter = edges_list->begin(); eiter!=edges_list->end(); eiter++)
		{
			HE_edge* he = *eiter;
			copy_edge_id(he->id, he->pair->vert->id, he->vert->id, new_mesh);
		}

		PTR_EDGE_LIST edgelist = new_mesh->get_edges_list();

		std::sort(edgelist->begin(), edgelist->end(), CompareEdgeID);

		new_mesh->update_mesh();

		VERTEX_ITER cviter = new_mesh->get_vertices_list()->begin();
		for (viter = vertices_list->begin(); viter!=vertices_list->end(); viter++, cviter++)
		{
			if ((*viter)->edge != NULL)
			{
				(*cviter)->edge = new_mesh->get_edge((*viter)->edge->id);
			}
		}

		FACE_ITER cfiter = new_mesh->get_faces_list()->begin();
		for (fiter = faces_list->begin(); fiter!=faces_list->end(); fiter++ , cfiter++)
		{
			(*cfiter)->edge = new_mesh->get_edge((*fiter)->edge->id);
		}

		for (int i = 0; i < edgelist->size(); i++)
		{
			(*edgelist)[i]->cycle_curv = (*edges_list)[i]->cycle_curv;
		}


		return new_mesh;
	}
	
	CS_INLINE Mesh3D* Mesh3D::reverse_orientation()
	{
		Mesh3D* new_mesh = new Mesh3D;

		VERTEX_ITER viter = vertices_list->begin();
		for (; viter!=vertices_list->end(); viter++)
		{
			new_mesh->insert_vertex((*viter)->pos);
		}

		FACE_ITER fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) 
		{
			HE_face* hf = *fiter;
			HE_edge* edge = hf->edge;
			VERTEX_LIST mvlist;
			do {
				mvlist.push_back(new_mesh->get_vertex(edge->pair->vert->id));
				edge = edge->next;
			} while(edge!=hf->edge);
			std::reverse(mvlist.begin(), mvlist.end());
			new_mesh->insert_face(mvlist);
		}

		new_mesh->update_mesh();
		return new_mesh;
	}

	
	CS_INLINE void Mesh3D::compute_num_components()
	{
		FACE_ITER fiter = faces_list->begin();
		for (; fiter != faces_list->end(); fiter++ )
		{
			(*fiter)->tag = false;
		}

		m_num_components = 0;

		for (fiter = faces_list->begin(); fiter != faces_list->end(); fiter++ )
		{
			HE_face* hf = *fiter;
			if (hf->tag == false)
			{
				m_num_components++;

				std::queue<HE_face* > facets;
				facets.push(hf);
				hf->tag = true;

				while (!facets.empty())
				{
					HE_face* pFacet = facets.front();
					facets.pop();
					pFacet->tag = true;

					HE_edge* he = pFacet->edge;
					do 
					{
						if(he->pair->face!=NULL && he->pair->face->tag==false)
						{
							facets.push(he->pair->face);
							he->pair->face->tag = true;

							HE_edge* vhe = he->vert->edge;

							do 
							{
								if (vhe->face && vhe->face->tag == false)
								{
									facets.push(vhe->face);
									vhe->face->tag = true;

									HE_edge* mvhe = vhe->vert->edge;
									do 
									{
										if (mvhe->face && mvhe->face->tag == false)
										{
											facets.push(mvhe->face);
											mvhe->face->tag = true;
										}
										mvhe = mvhe->pair->next;
									} while(mvhe != vhe->vert->edge);

								}
								vhe = vhe->pair->next;
							} while(vhe != he->vert->edge);

						}
						he = he->next;
					} while(he != pFacet->edge);

					he = pFacet->edge;
				}
			}
		}
	}
	
	CS_INLINE void Mesh3D::compute_num_boundaries()
	{

		VERTEX_ITER viter = vertices_list->begin();
		for (; viter != vertices_list->end(); viter++ )
		{
			(*viter)->tag = false;
		}

		boundaryvertices.clear();

		m_num_boundaries = 0;

		for(viter = vertices_list->begin(); viter != vertices_list->end(); viter++)
		{

			HE_vert* hv = *viter;

			if (is_on_boundary(hv) && hv->tag == false)
			{

				std::vector<HE_vert*> vedges;
				vedges.push_back(hv);

				m_num_boundaries++;
				std::queue<HE_vert* > vertices;
				hv->tag = true;
				vertices.push(hv);

				while (!vertices.empty())
				{
					HE_vert* pVertex = vertices.front();
					pVertex->tag = true;
					vertices.pop();

					HE_edge* he = pVertex->edge;

					do 
					{
						if (is_on_boundary(he) && is_on_boundary(he->vert) && he->vert->tag==false)
						{
							he->vert->tag = true;
							vertices.push(he->vert);
							vedges.push_back(he->vert);

						}
						he = he->pair->next;
					} while(he!=pVertex->edge);

				}

				HE_vert* hvtmp = vedges[0];
				vedges[0] = vedges[1];
				vedges[1] = hvtmp;
				boundaryvertices.push_back(vedges);
			}
		}
	}
	
	CS_INLINE void Mesh3D::compute_genus()
	{
		compute_num_components();
		compute_num_boundaries();
		int c = m_num_components;
		int b = m_num_boundaries;
		ptrdiff_t v = get_num_of_vertices();
		ptrdiff_t e = get_num_of_edges()/2;
		ptrdiff_t f = get_num_of_faces();

		m_genus = (int)(2*c+e-b-f-v)/2;
	}
	
	CS_INLINE void Mesh3D::remove_hanged_vertices()
	{

		bool find = false;
		VERTEX_ITER viter = vertices_list->begin();
		for (;viter!=vertices_list->end(); viter++)
		{
			if ((*viter)->edge==NULL)
			{
				find = true;
				break;
			}
		}

		if (find)
		{
			VERTEX_LIST* new_vertices_list = new VERTEX_LIST;
			int i = 0;
			for (viter = vertices_list->begin();viter!=vertices_list->end(); viter++)
			{
				if ((*viter)->edge!=NULL)
				{
					new_vertices_list->push_back(*viter);
					(*viter)->id = i;
					i++;
				}
				else
					delete *viter;
			}
			delete vertices_list;
			vertices_list = new_vertices_list;

			EDGE_LIST* new_edges_list = new EDGE_LIST;
			i = 0;
			for (EDGE_ITER eiter = edges_list->begin();eiter!=edges_list->end(); eiter++)
			{
				if ((*eiter)->face!=NULL||(*eiter)->pair->face!=NULL)
				{
					new_edges_list->push_back(*eiter);
					(*eiter)->id = i;
					i++;
				}
				else
					delete *eiter;
			}
			delete edges_list;
			edges_list = new_edges_list;
		}
	}
	
	CS_INLINE void Mesh3D::reset_vertices_tag(bool tag_status)
	{
		VERTEX_ITER viter = vertices_list->begin();
		for (; viter != vertices_list->end(); viter++)
			(*viter)->tag = tag_status;
	}
	
	CS_INLINE void Mesh3D::reset_faces_tag(bool tag_status)
	{
		FACE_ITER fiter = faces_list->begin();
		for (; fiter !=faces_list->end(); fiter++)
			(*fiter)->tag = tag_status;
	}
	
	CS_INLINE void Mesh3D::reset_edges_tag(bool tag_status)
	{
		EDGE_ITER eiter = edges_list->begin();
		for (; eiter != edges_list->end(); eiter++)
			(*eiter)->tag = tag_status;
	}
	
	CS_INLINE void Mesh3D::reset_all_tag(bool tag_status)
	{
		reset_edges_tag(tag_status);
		reset_faces_tag(tag_status);
		reset_vertices_tag(tag_status);
	}
	
	CS_INLINE void Mesh3D::translate(const Eigen::Vector3d& tran_V)
	{
		VERTEX_ITER viter = vertices_list->begin();
		for (; viter != vertices_list->end(); viter++)
		{
			HE_vert* hv = *viter;
			hv->pos += tran_V;
		}
	}
	
	CS_INLINE void Mesh3D::scale(double factorx, double factory, double factorz)
	{
		VERTEX_ITER viter = vertices_list->begin();
		for (; viter != vertices_list->end(); viter++)
		{
			HE_vert* hv = *viter;
			hv->pos[0] *= factorx;
			hv->pos[1] *= factory;
		}
	}
	
	CS_INLINE void Mesh3D::init_edge_tag()
	{
		reset_edges_tag(false);
		EDGE_ITER eiter = edges_list->begin();
		for (; eiter != edges_list->end(); eiter++)
		{
			HE_edge* he = *eiter;
			if (he->tag==false && he->pair->tag == false)
			{
				he->tag = true;
			}
		}
	}


	CS_INLINE void Mesh3D::swap_edge( HE_edge* triedge )
	{
		if (triedge==NULL || triedge->face==NULL || triedge->pair->face==NULL || triedge->face->valence!=3 || triedge->pair->face->valence != 3)
		{
			return;
		}

		HE_vert* hv2 = triedge->vert;
		HE_vert* hv1 = triedge->pair->vert;
		HE_vert* hv3 = triedge->next->vert;
		HE_vert* hv4 = triedge->pair->next->vert;


		HE_edge* che = hv3->edge;
		do 
		{
			if (che->vert == hv4)
			{
				return;
			}
			che = che->pair->next;
		} while (che != hv3->edge);


		HE_face* hf1 = triedge->face;
		HE_face* hf2 = triedge->pair->face;


		HE_edge* he1 = triedge->next->next;
		HE_edge* he2 = triedge->next;
		HE_edge* he3 = triedge->pair->next->next;
		HE_edge* he4 = triedge->pair->next;


		triedge->vert = hv3;
		triedge->next = he1;
		triedge->prev = he4;

		triedge->pair->vert = hv4;
		triedge->pair->next = he3;
		triedge->pair->prev = he2;

		hf1->edge = triedge;
		hf2->edge = triedge->pair;

		he1->face = hf1;
		he2->face = hf2;
		he3->face = hf2;
		he4->face = hf1;

		he1->next = he4;
		he1->prev = triedge;
		he2->next = triedge->pair;
		he2->prev = he3;
		he3->next = he2;
		he3->prev = triedge->pair;
		he4->next = triedge;
		he4->prev = he1;

		if (hv1->edge == triedge)
		{
			hv1->edge = he1->pair;
		}
		if (hv2->edge == triedge->pair)
		{
			hv2->edge = he3->pair;
		}

		hv1->degree--;
		hv2->degree--;
		hv3->degree++;
		hv4->degree++;

	}

		
	CS_INLINE Eigen::Vector3d project_to_segment(Eigen::Vector3d pt, Eigen::Vector3d s0, Eigen::Vector3d s1){
		Eigen::Vector3d dir = (s1-s0).normalized();
		double t = (pt-s0).dot(dir)/(s1-s0).norm();

		if (t<0)
		{
			t=0;
		}
		else if (t>1)
		{
			t=1;
		}

		return (1-t)*s0 + t*s1;
	}


	CS_INLINE void Mesh3D::split_face(HE_face* face, HE_vert* vtx)
	{
		int nv = (int)vertices_list->size(); 
		int ne = (int)edges_list->size();
		int nf = (int)faces_list->size();
		int face_valence = face->valence;

		assert(face_valence == 3);

		HE_vert* new_vertex = vtx;
		new_vertex->degree = face_valence;

		std::vector<HE_face* > new_faces(face_valence);
		std::vector<HE_edge* > new_edges(2*face_valence);

		new_edges.resize(2*face_valence);
		new_faces.resize(face_valence);

		std::vector<HE_edge*> face_edges(face_valence);
		std::vector<HE_vert*> face_vertices(face_valence);
		HE_edge* he = face->edge;
		for (int i = 0; i < face_valence; i++)
		{
			face_vertices[i] = he->vert;
			face_vertices[i]->degree++;
			face_edges[i] = he;
			he = he->next;
			if (i == 0)
			{
				new_faces[i] = face;
			}
			else
			{
				new_faces[i] = new HE_face;
				new_faces[i]->id = nf;
				nf++;
			}
			new_faces[i]->valence = 3;
			new_faces[i]->edge = face_edges[i];
			face_edges[i]->face = new_faces[i];
		}

		for (int i = 0; i < face_valence; i++)
		{
			new_edges[2*i] = new HE_edge;
			new_edges[2*i+1] = new HE_edge;
			new_edges[2*i]->pair = new_edges[2*i+1];
			new_edges[2*i+1]->pair = new_edges[2*i];
			new_edges[2*i]->face = new_faces[i];
			new_edges[2*i+1]->face = new_faces[(i+1)%face_valence];
			new_edges[2*i]->vert = new_vertex;
			new_edges[2*i+1]->vert = face_vertices[i];
			new_edges[2*i]->id = ne;
			new_edges[2*i+1]->id = ne+1;
			ne += 2;
		}

		new_vertex->edge = new_edges[1];

		for (int i = 0; i < face_valence; i++)
		{
			new_edges[2*i]->prev = face_edges[i];
			new_edges[2*i]->next = i==0? new_edges[2*face_valence-1] : new_edges[2*i-1];
			new_edges[2*i+1]->prev = i==face_valence-1? new_edges[0] : new_edges[2*i+2];
			new_edges[2*i+1]->next = face_edges[(i+1)%face_valence];

			new_edges[2*i]->prev->next = new_edges[2*i];
			new_edges[2*i]->next->prev = new_edges[2*i];
			new_edges[2*i+1]->prev->next = new_edges[2*i+1];
			new_edges[2*i+1]->next->prev = new_edges[2*i+1];
		}
		edges_list->insert(edges_list->end(), new_edges.begin(), new_edges.end());
		std::vector<HE_face*>::iterator fiter = new_faces.begin();
		fiter++;
		faces_list->insert(faces_list->end(), fiter, new_faces.end());

	}

CS_INLINE void Mesh3D::split_edge(HE_edge* edge, HE_vert* vtx) {
    if (!edge || !vtx || !edge->face || !edge->pair || !edge->pair->face) {
        return;
    }

    HE_edge* e0 = edge;             // half-edge (v0 -> v1)
    HE_edge* ep0 = edge->pair;      // opposite half-edge (v1 -> v0)
    HE_face* f0 = e0->face;         // triangle (v0, v1, v2)
    HE_face* f1 = ep0->face;        // triangle (v1, v0, v3)

    HE_edge* e1 = e0->next;   // (v1 -> v2)
    HE_edge* e2 = e1->next;   // (v2 -> v0)
    HE_edge* ep1 = ep0->next; // (v0 -> v3)
    HE_edge* ep2 = ep1->next; // (v3 -> v1)

    HE_vert* v0 = e2->vert;   // head of e2 (v0)
    HE_vert* v1 = e0->vert;   // head of e0 (v1)
    HE_vert* v2 = e1->vert;   // head of e1 (v2)
    HE_vert* v3 = ep1->vert;  // head of ep1 (v3)

    HE_vert* vm = vtx;
    vm->degree = 4;

    int base_id = (int)edges_list->size();

    HE_edge* e0a = new HE_edge; HE_edge* e0b = new HE_edge;
    HE_edge* ep0a = new HE_edge; HE_edge* ep0b = new HE_edge;
    HE_edge* e_vm_v2 = new HE_edge; HE_edge* e_v2_vm = new HE_edge;
    HE_edge* ep_vm_v3 = new HE_edge; HE_edge* ep_v3_vm = new HE_edge;

    e0a->id = base_id + 0; e0b->id = base_id + 1;
    ep0a->id = base_id + 2; ep0b->id = base_id + 3;
    e_vm_v2->id = base_id + 4; e_v2_vm->id = base_id + 5;
    ep_vm_v3->id = base_id + 6; ep_v3_vm->id = base_id + 7;

    e0a->pair = ep0b; ep0b->pair = e0a;
    e0b->pair = ep0a; ep0a->pair = e0b;
    e_vm_v2->pair = e_v2_vm; e_v2_vm->pair = e_vm_v2;
    ep_vm_v3->pair = ep_v3_vm; ep_v3_vm->pair = ep_vm_v3;

    e0a->vert = vm;    // (v0 -> vm)
    e0b->vert = v1;    // (vm -> v1)
    ep0a->vert = vm;   // (v1 -> vm)
    ep0b->vert = v0;   // (vm -> v0)

    e_vm_v2->vert = v2; // (vm -> v2)  <-- important: vm->v2 has head v2
    e_v2_vm->vert = vm; // (v2 -> vm)

    ep_vm_v3->vert = v3; // (vm -> v3)  <-- important: vm->v3 has head v3
    ep_v3_vm->vert = vm; // (v3 -> vm)

    HE_face* f0_new = new HE_face; // (vm,v1,v2)
    HE_face* f1_new = new HE_face; // (vm,v0,v3)

    f0_new->id = (int)faces_list->size();
    f1_new->id = (int)faces_list->size() + 1;
    f0_new->valence = 3;
    f1_new->valence = 3;


    e0a->next = e_vm_v2;
    e_vm_v2->next = e2;
    e2->next = e0a;

    e0a->prev = e2; e_vm_v2->prev = e0a; e2->prev = e_vm_v2;

    e0a->face = f0;
    e_vm_v2->face = f0;
    e2->face = f0;
    f0->edge = e0a;
    f0->valence = 3;

    e0b->next = e1;
    e1->next = e_v2_vm;
    e_v2_vm->next = e0b;

    e0b->prev = e_v2_vm; e1->prev = e0b; e_v2_vm->prev = e1;

    e0b->face = f0_new;
    e1->face = f0_new;
    e_v2_vm->face = f0_new;
    f0_new->edge = e0b;

    ep0a->next = ep_vm_v3;
    ep_vm_v3->next = ep2;
    ep2->next = ep0a;

    ep0a->prev = ep2; ep_vm_v3->prev = ep0a; ep2->prev = ep_vm_v3;

    ep0a->face = f1;
    ep_vm_v3->face = f1;
    ep2->face = f1;
    f1->edge = ep0a;
    f1->valence = 3;

    ep0b->next = ep1;
    ep1->next = ep_v3_vm;
    ep_v3_vm->next = ep0b;

    ep0b->prev = ep_v3_vm; ep1->prev = ep0b; ep_v3_vm->prev = ep1;

    ep0b->face = f1_new;
    ep1->face = f1_new;
    ep_v3_vm->face = f1_new;
    f1_new->edge = ep0b;

    vm->edge = e0b;

    if (v0->edge == e0 || v0->edge == ep0) {
        if (e2->pair) v0->edge = e2->pair;
        else v0->edge = ep0b;
    }
    if (v1->edge == e0 || v1->edge == ep0) {
        if (e1->pair) v1->edge = e1->pair;
        else v1->edge = e0b;
    }
    if (v2->edge == e1) v2->edge = e_v2_vm;
    if (v3->edge == ep1) v3->edge = ep_v3_vm;

    for (auto* vertex : {v0, v1, v2, v3}) {
        if (vertex->edge == e0 || vertex->edge == ep0) {
            bool found_replacement = false;
            for (auto* edge : *edges_list) {
                if (edge->pair && edge->pair->vert == vertex && edge != e0 && edge != ep0) {
                    vertex->edge = edge;
                    found_replacement = true;
                    break;
                }
            }
            if (!found_replacement) {
                if (vertex == v0) vertex->edge = ep0b;
                else if (vertex == v1) vertex->edge = e0b;
                else if (vertex == v2) vertex->edge = e_v2_vm;
                else if (vertex == v3) vertex->edge = ep_v3_vm;
            }
        }
    }

    v0->degree += 1;
    v1->degree += 1;
    v2->degree += 1;
    v3->degree += 1;

    auto erase_edge_from_list = [&](HE_edge* e) {
        for (auto it = edges_list->begin(); it != edges_list->end(); ++it) {
            if ((*it) == e) {
                edges_list->erase(it);
                delete e;
                return;
            }
        }
        delete e;
    };

    for (auto* vertex : *vertices_list) {
        if (vertex->edge == e0 || vertex->edge == ep0) {
            bool found_replacement = false;
            for (auto* candidate_edge : *edges_list) {
                if (candidate_edge->pair && 
                    candidate_edge->pair->vert == vertex && 
                    candidate_edge != e0 && 
                    candidate_edge != ep0) {
                    vertex->edge = candidate_edge;
                    found_replacement = true;
                    break;
                }
            }
            
            if (!found_replacement) {
                if (vertex == v0) vertex->edge = ep0b;
                else if (vertex == v1) vertex->edge = e0b; 
                else if (vertex == v2) vertex->edge = e_v2_vm;
                else if (vertex == v3) vertex->edge = ep_v3_vm;
                else {
                    if (vm->edge && vm->edge->pair && vm->edge->pair->vert == vertex) {
                        vertex->edge = vm->edge;
                    } else {
                        vertex->edge = e0b; // fallback
                    }
                }
            }
        }
    }

    e0->pair = nullptr; e0->next = nullptr; e0->prev = nullptr; e0->face = nullptr;
    ep0->pair = nullptr; ep0->next = nullptr; ep0->prev = nullptr; ep0->face = nullptr;

    erase_edge_from_list(e0);
    erase_edge_from_list(ep0);

    vertices_list->push_back(vm);

    edges_list->push_back(e0a);
    edges_list->push_back(e0b);
    edges_list->push_back(ep0a);
    edges_list->push_back(ep0b);
    edges_list->push_back(e_vm_v2);
    edges_list->push_back(e_v2_vm);
    edges_list->push_back(ep_vm_v3);
    edges_list->push_back(ep_v3_vm);

    faces_list->push_back(f0_new);
    faces_list->push_back(f1_new);
}

CS_INLINE void Mesh3D::laplacian_smooth(int iterations, double damping_factor, bool preserve_boundary)
{
    if (!vertices_list || vertices_list->empty()) {
        return;
    }

    damping_factor = std::max(0.0, std::min(1.0, damping_factor));

    int vertices_moved = 0;
    double total_movement = 0.0;

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<Eigen::Vector3d> new_positions(vertices_list->size());
        std::vector<bool> should_update(vertices_list->size(), true);
        
        for (size_t i = 0; i < vertices_list->size(); ++i) {
            new_positions[i] = (*vertices_list)[i]->pos;
        }

        for (size_t i = 0; i < vertices_list->size(); ++i) {
            HE_vert* vertex = (*vertices_list)[i];
            
            if (preserve_boundary && is_on_boundary(vertex)) {
                should_update[i] = false;
                continue;
            }

            if (!vertex->edge) {
                should_update[i] = false;
                continue;
            }

            std::vector<HE_vert*> neighbors;
            std::set<HE_vert*> neighbor_set; // Use set to avoid duplicates
            
            for (int e = 0; e < edges_list->size(); ++e) {
                HE_edge* edge = (*edges_list)[e];
                if (!edge) continue;
                
                if (edge->vert == vertex && edge->pair && edge->pair->vert) {
                    neighbor_set.insert(edge->pair->vert);
                } else if (edge->pair && edge->pair->vert == vertex && edge->vert) {
                    neighbor_set.insert(edge->vert);
                }
            }
            
            for (HE_vert* neighbor : neighbor_set) {
                neighbors.push_back(neighbor);
            }

            if (!neighbors.empty()) {
                Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
                for (HE_vert* neighbor : neighbors) {
                    centroid += neighbor->pos;
                }
                centroid /= static_cast<double>(neighbors.size());
                
                Eigen::Vector3d old_pos = vertex->pos;
                new_positions[i] = vertex->pos + damping_factor * (centroid - vertex->pos);
                
                double movement = (new_positions[i] - old_pos).norm();
                total_movement += movement;
                vertices_moved++;
                
            } else {
                should_update[i] = false;
            }
        }

        for (size_t i = 0; i < vertices_list->size(); ++i) {
            if (should_update[i]) {
                (*vertices_list)[i]->pos = new_positions[i];
            }
        }
    }

    update_normal();

    double avg_movement = vertices_moved > 0 ? total_movement / vertices_moved : 0.0;
}

CS_INLINE void Mesh3D::laplacian_smooth_constraints(double lambda, bool preserve_boundary)
{
    if (!vertices_list || vertices_list->empty()) {
        return;
    }

    std::vector<std::vector<int>> vertex_neighbors(vertices_list->size());
    std::vector<bool> is_boundary_vertex(vertices_list->size(), false);
    
    for (size_t e = 0; e < edges_list->size(); ++e) {
        HE_edge* edge = (*edges_list)[e];
        if (!edge || !edge->vert || !edge->pair || !edge->pair->vert) continue;
        
        int v0_id = edge->vert->id;
        int v1_id = edge->pair->vert->id;
        
        if (v0_id < vertices_list->size() && v1_id < vertices_list->size()) {
            vertex_neighbors[v0_id].push_back(v1_id);
            vertex_neighbors[v1_id].push_back(v0_id);
            
            if (is_on_boundary(edge)) {
                is_boundary_vertex[v0_id] = true;
                is_boundary_vertex[v1_id] = true;
            }
        }
    }
    
    for (auto& neighbors : vertex_neighbors) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    std::vector<int> free_vertex_ids;
    for (size_t i = 0; i < vertices_list->size(); ++i) {
        if (!preserve_boundary || !is_boundary_vertex[i]) {
            if (!vertex_neighbors[i].empty()) { // Only include vertices with neighbors
                free_vertex_ids.push_back(static_cast<int>(i));
            }
        }
    }
    
    if (free_vertex_ids.empty()) {
        return;
    }

    std::map<int, int> vertex_to_var;
    for (size_t i = 0; i < free_vertex_ids.size(); ++i) {
        vertex_to_var[free_vertex_ids[i]] = static_cast<int>(i);
    }

    std::vector<Geex::LinearConstraint<double>> constraints;

    for (size_t i = 0; i < free_vertex_ids.size(); ++i) {
        int vertex_id = free_vertex_ids[i];
        const auto& neighbors = vertex_neighbors[vertex_id];
        
        if (neighbors.empty()) continue;
        
        for (int coord = 0; coord < 3; ++coord) {
            Geex::LinearConstraint<double> constraint;
            constraint.Weight = lambda; // Smoothing weight
            constraint.RHSide = 0.0; // Will be computed as weighted sum
            
            Geex::VarTerm<double> vertex_term(i, coord, 1.0);
            constraint.LHSide.push_back(vertex_term);
            
            double neighbor_weight = -1.0 / static_cast<double>(neighbors.size());
            double rhs_contribution = 0.0;
            
            for (int neighbor_id : neighbors) {
                auto it = vertex_to_var.find(neighbor_id);
                if (it != vertex_to_var.end()) {
                    Geex::VarTerm<double> neighbor_term(it->second, coord, neighbor_weight);
                    constraint.LHSide.push_back(neighbor_term);
                } else {
                    HE_vert* neighbor_vert = (*vertices_list)[neighbor_id];
                    double neighbor_coord_val = (coord == 0) ? neighbor_vert->pos.x() : 
                                                (coord == 1) ? neighbor_vert->pos.y() : 
                                                               neighbor_vert->pos.z();
                    rhs_contribution += neighbor_weight * neighbor_coord_val;
                }
            }
            
            constraint.RHSide = -rhs_contribution; // Note: constraint is LHS = RHS
            constraints.push_back(constraint);
        }
    }

    double position_weight = 0.1; // Small weight to prefer staying close to current position
    for (size_t i = 0; i < free_vertex_ids.size(); ++i) {
        int vertex_id = free_vertex_ids[i];
        HE_vert* vertex = (*vertices_list)[vertex_id];
        
        for (int coord = 0; coord < 3; ++coord) {
            Geex::LinearConstraint<double> constraint;
            constraint.Weight = position_weight;
            
            double current_coord = (coord == 0) ? vertex->pos.x() : 
                                   (coord == 1) ? vertex->pos.y() : 
                                                  vertex->pos.z();
            constraint.RHSide = current_coord;
            
            Geex::VarTerm<double> var_term(i, coord, 1.0);
            constraint.LHSide.push_back(var_term);
            
            constraints.push_back(constraint);
        }
    }

    Geex::LinearConstraintLeastSQSolver<double> solver;
    solver.initialize_system(constraints);
    
    std::vector<double> solution;
    solver.solve_system(solution);
    
    if (solution.size() != free_vertex_ids.size() * 3) {
        return;
    }

    double total_movement = 0.0;
    int vertices_moved = 0;
    
    for (size_t i = 0; i < free_vertex_ids.size(); ++i) {
        int vertex_id = free_vertex_ids[i];
        HE_vert* vertex = (*vertices_list)[vertex_id];
        
        Eigen::Vector3d old_pos = vertex->pos;
        Eigen::Vector3d new_pos;
        
        new_pos.x() = solution[i * 3 + 0];
        new_pos.y() = solution[i * 3 + 1];
        new_pos.z() = solution[i * 3 + 2];
        
        vertex->pos = new_pos;
        
        double movement = (new_pos - old_pos).norm();
        total_movement += movement;
        vertices_moved++;
    }

    update_normal();
    
    double avg_movement = vertices_moved > 0 ? total_movement / vertices_moved : 0.0;

    std::vector<double> residuals;
    solver.compute_constraint_residuals(solution, residuals);
    double total_residual_energy = 0.0;
    for (double r : residuals) {
        total_residual_energy += r * r;
    }
}


	CS_INLINE void Mesh3D::make_delaunay( bool& topo_updated)
	{
		topo_updated = false;

		std::vector<HE_edge*> edge_v;

		reset_edges_tag(false);
		int edge_counter = 0;
		for (EDGE_ITER eitr = edges_list->begin(); eitr != edges_list->end(); ++eitr)
		{
			HE_edge* eh = *eitr;
			if (eh->tag || eh->pair->tag)
			{
				continue;
			}
			edge_v.push_back(eh);
			eh->tag = true;
			eh->pair->tag = true;
		}


		const int max_itr = edges_list->size()*5;
		int itr_count = 0;

		while (edge_v.size()>0 && itr_count < max_itr)
		{
			itr_count ++;

			HE_edge* eh = edge_v.back();
			edge_v.pop_back();
			eh->tag = false;
			eh->pair->tag = false;

			if (eh->face == NULL || eh->pair->face == NULL)
			{
				continue;
			}

			if (eh->cycle_curv)
			{
				continue;
			}

			if (!delaunay_test(eh))
			{

				if (0)
				{//unflippable


				} else{
					topo_updated = true;

					swap_edge(eh);
					topo_updated = true;

					if (!eh->next->tag && !eh->next->pair->tag)
					{
						eh->next->tag = eh->next->pair->tag = true;
						edge_v.push_back(eh->next);
					}
					if (!eh->prev->tag && !eh->prev->pair->tag)
					{
						eh->prev->tag = eh->prev->pair->tag = true;
						edge_v.push_back(eh->prev);
					}
					if (!eh->pair->next->tag && !eh->pair->next->pair->tag)
					{
						eh->pair->next->tag = eh->pair->next->pair->tag = true;
						edge_v.push_back(eh->pair->next);
					}
					if (!eh->pair->prev->tag && !eh->pair->prev->pair->tag)
					{
						eh->pair->prev->tag = eh->pair->prev->pair->tag = true;
						edge_v.push_back(eh->pair->prev);
					}
				}
			}
		}


	}


	CS_INLINE bool Mesh3D::delaunay_test( HE_edge* triedge )
	{
		Eigen::Vector3d v0, v1, e0, e1;
		v0 = triedge->next->vert->pos;
		v1 = triedge->pair->next->vert->pos;
		e0 = triedge->vert->pos;
		e1 = triedge->pair->vert->pos;

		double angle1 = acos((e0-v0).normalized().dot((e1-v0).normalized()));
		double angle2 = acos((e0 - v1).normalized().dot((e1 - v1).normalized()));

		if (angle1 + angle2 < M_PI)
		{
			return true;
		}

		return false;
	}

	CS_INLINE PolyGraph_Edge* Mesh3D::vertex_on_curve(HE_vert* vh){
		HE_edge* eh = vh->edge;
		do
		{


			if(eh->cycle_curv)
				return eh->cycle_curv;

			if (eh->pair->cycle_curv)
				return eh->pair->cycle_curv;

			eh=eh->pair->next;
		} while (eh!=vh->edge);

		return NULL;
	}

	CS_INLINE int Mesh3D::vertex_curves(HE_vert* vh){
		int num = 0;
		HE_edge* eh = vh->edge;
		do
		{

			if(eh->cycle_curv || eh->pair->cycle_curv)
				num++;

			eh=eh->pair->next;
		} while (eh!=vh->edge);

		return num;
	}

	CS_INLINE double Mesh3D::edge_chain_crvt(HE_edge* eij){
		HE_face* f0 = eij->face, *f1 = eij->pair->face;
		if (f0 && f1)
		{
			Eigen::Vector3d f0_axis = (f0->edge->vert->pos - f0->edge->prev->vert->pos).normalized();
			Eigen::Vector3d f1_axis = (f1->edge->vert->pos - f1->edge->prev->vert->pos).normalized();
			Eigen::Vector3d vec_ij = (eij->vert->pos - eij->pair->vert->pos).normalized();
			Eigen::Vector3d vec_ji = -vec_ij;

			Eigen::Vector3d rdir0 = f0_axis.cross(vec_ij);
			Eigen::Vector3d rdir1 = vec_ji.cross(f1_axis);

			double angle0 , angle1;

			if (f0->edge == eij)
			{
				angle0 = 0;
			}
			else
			{
				angle0 = acos(f0_axis.dot(vec_ij));
				if (rdir0.dot(f0->normal) < 0)
				{
					angle0 *= -1;
				}
			}

			if (f1->edge == eij->pair)
			{
				angle1 = 0;
			}
			else
			{
				angle1 = acos(f1_axis.dot(vec_ji));
				if (rdir1.dot(f1->normal) < 0)
				{
					angle1 *= -1;
				}
			}

			double kij = angle0 + angle1 - M_PI;
			if (kij > M_PI)
			{
				kij -= 2*M_PI;
			}
			else if (kij <= - M_PI)
			{
				kij += 2*M_PI;
			}

			if (!std::isnan(kij))
			{
				return -kij;
			}
			else
			{
				return 0;
			}

			
		}

		return 0;
	}

	CS_INLINE double Mesh3D::bc_dual_edge_length(HE_edge* edge){
		if (edge->face && edge->pair->face)
		{
			double l0 = (edge->face->GetCentroid() - edge->GetMidPoint()).norm();
			double l1 = (edge->pair->face->GetCentroid() - edge->GetMidPoint()).norm();

			return l0+l1;
		}
		return (edge->face->GetCentroid() - edge->GetMidPoint()).norm();
	}

	CS_INLINE double Mesh3D::ortho_dual_edge_data(HE_edge* edge, Eigen::Vector3d& dir, Eigen::Vector3d& ndiff, bool use_cot_w /*= false*/) const{

		double edge_len = (edge->vert->pos - edge->prev->vert->pos).norm();
		

		dir = -(edge->face->normal).cross(edge->vert->pos - edge->prev->vert->pos); //normalize
		if (edge->face && edge->pair->face)
		{
			ndiff = (edge->pair->face->normal - edge->face->normal);//*edge_len

			if (use_cot_w)
			{
				Eigen::Vector3d pt[4];
				pt[0] = edge->prev->vert->pos;
				pt[1] = edge->vert->pos;
				pt[2] = edge->next->vert->pos;
				pt[3] = edge->pair->next->vert->pos;

				Eigen::Vector3d e20 = (pt[0] - pt[2]).normalized();
				Eigen::Vector3d e21 = (pt[1] - pt[2]).normalized();
				Eigen::Vector3d e30 = (pt[0] - pt[3]).normalized();
				Eigen::Vector3d e31 = (pt[1] - pt[3]).normalized();

				double cos_2 = e20.dot(e21);
				double cos_3 = e30.dot(e31);

				double sin_2 = e20.cross(e21).norm();
				double sin_3 = e30.cross(e31).norm();

				double cot_23 = cos_2/sin_2 + cos_3/sin_3;

				double dual_ratio = cot_23*0.5;

				if (std::isnan(dual_ratio)) // || dual_len > 1e6
				{
					dual_ratio = 0;
				}
				else
				{
					// Clamp: degenerate (needle) triangles otherwise produce
					// unbounded cot weights that destabilize the solves
					// downstream. Bound ≈ cot(1°).
					const double kCotClamp = 57.29;
					dual_ratio = std::max(-kCotClamp, std::min(kCotClamp, dual_ratio));
				}

				dir *= dual_ratio ;

				return dual_ratio;
			}
		}

		return 1;
	}


	CS_INLINE Mesh3D* Mesh3D::make_copy_regular(){
		Mesh3D* new_mesh = new Mesh3D;

		VERTEX_ITER viter = vertices_list->begin();
		for (; viter!=vertices_list->end(); viter++)
		{
			HE_vert* vh = new_mesh->insert_vertex((*viter)->pos);

			if ((*viter)->fixed)
			{
				vh->fixed = true;
			}
		}

		new_mesh->edges_list = new EDGE_LIST;

		for (EDGE_ITER eiter = edges_list->begin(); eiter!=edges_list->end(); eiter++)
		{
			HE_edge* he = *eiter;
			HE_edge* nedge = new HE_edge;
			nedge->vert = new_mesh->get_vertex(he->vert->id);
			nedge->id = he->id;
			new_mesh->edges_list->push_back(nedge);
		}

		PTR_EDGE_LIST edgelist = new_mesh->get_edges_list();

		std::sort(edgelist->begin(), edgelist->end(), CompareEdgeID);

		
		VERTEX_ITER cviter = new_mesh->get_vertices_list()->begin();
		for (viter = vertices_list->begin(); viter!=vertices_list->end(); viter++, cviter++)
		{
			if ((*viter)->edge != NULL)
			{
				(*cviter)->edge = new_mesh->get_edge((*viter)->edge->id);
			}
		}

		for (EDGE_ITER ceitr = new_mesh->get_edges_list()->begin(), eitr = edges_list->begin(); ceitr != new_mesh->get_edges_list()->end(); ceitr++, eitr++)
		{
			HE_edge* neh = *ceitr, *eh = *eitr;
			neh->next = new_mesh->get_edge(eh->next->id);
			neh->prev = new_mesh->get_edge(eh->prev->id);
			neh->pair = new_mesh->get_edge(eh->pair->id);
		}

		FACE_ITER fiter = faces_list->begin();
		new_mesh->faces_list = new FACE_LIST;

		for (;fiter!=faces_list->end(); fiter++) 
		{
			HE_face* hf = *fiter;
			HE_edge* edge = hf->edge;

			HE_face* nface = new HE_face;
			nface->edge = new_mesh->get_edge(edge->id);
			nface->valence = hf->valence;

			nface->id = hf->id;
			new_mesh->faces_list->push_back(nface);
			
			do {
				
				new_mesh->get_edge(edge->id)->face = nface;

				edge = edge->next;
			} while(edge!=hf->edge);

		}


		new_mesh->update_mesh();


		for (int i = 0; i < edgelist->size(); i++)
		{
			(*edgelist)[i]->cycle_curv = (*edges_list)[i]->cycle_curv;
		}


		return new_mesh;
	}


}
#undef  _CRT_SECURE_NO_WARNINGS