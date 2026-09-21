/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

#ifndef UV_POISSON_SOLVER
#define UV_POISSON_SOLVER

#include <point2.h>
#include <point3.h>
#include <eigen_interface.h>
#include <box3.h>
#include <box2.h>
#include <mesh_functions.h>
#include <connectivity.h>
#include <UV_utils.h>
#include <UV_distortion.h>
#include "normal.h"

namespace Geo
{

    template <class ScalarType>
    class UVPoissonSolver
    {

        const std::vector<Geo::Point3<ScalarType>> &VertPos;
        const std::vector<std::vector<int>> &Faces;
        const std::vector<Field::CrossF<ScalarType>> &FaceCurv;

        // ///solver data
        // std::map<VertexType*,int> VertexToInd;
        // std::map<int, VertexType*> IndToVertex;

       
        Eigen::SparseMatrix<ScalarType> A; // A
        Eigen::VectorX<ScalarType> b, x;   // x and b

        // number of variables
        unsigned int n_vert_vars;
        /// total system size
        unsigned int total_size;
        /// number of fixed variables
        unsigned int n_fixed_vars;

        /// if you intend to follow the cross field
        bool use_direction_field, correct_fixed;

        /// size of the scalar field
        ScalarType fieldScale;

        // int VertexIndex(VertexType* v)
        // {
        //     typename std::map<VertexType*,int>::iterator iteMap=VertexToInd.find(v);
        //     assert(iteMap!=VertexToInd.end());
        //     return ((*iteMap).second);
        // }

        // VertexType* IndexVertex(int index)
        // {
        //     typename std::map<int,VertexType*>::iterator iteMap=IndToVertex.find(index);
        //     assert(iteMap!=IndToVertex.end());
        //     return ((*iteMap).second);
        // }

        // void AddVertexIndex(VertexType* v,int index)
        // {
        //     VertexToInd.insert(std::pair<VertexType*,int>(v,index));
        //     IndToVertex.insert(std::pair<int,VertexType*>(index,v));
        // }

        /// set the value of A of the system Ax=b
        void SetValA(int Xindex, int Yindex, ScalarType val)
        {
            assert(0 <= Xindex && Xindex < int(total_size));
            assert(0 <= Yindex && Yindex < int(total_size));
            A.coeffRef(Xindex, Yindex) += val;
        }

        // void FindFarthestVert(VertexType* &v0, VertexType* &v1)
        // {
        //   v0=v1=NULL;
        //   const int bestAxis = mesh.bbox.MaxDim();
        //   for(VertexType &vv : mesh.vert) {
        //     if(vv.P()[bestAxis] <= mesh.bbox.min[bestAxis]) v0 = &vv;
        //     if(vv.P()[bestAxis] >= mesh.bbox.max[bestAxis]) v1 = &vv;
        //   }
        //   assert( (v0!=v1) && v0 && v1);
        // }

        /// set the value of b of the system Ax=b
        void SetValB(int Xindex,
                     ScalarType val)
        {
            /*S.b()[Xindex] += val;*/
            b[Xindex] += val;
        }

        /// add the area term, scalefactor is used to sum up
        /// and normalize on the overlap zones
        void AddAreaTerm(int index[3][3][2], ScalarType ScaleFactor)
        {
            const ScalarType entry = 0.5 * ScaleFactor;
            ScalarType val[3][3] = {{0, entry, -entry},
                                    {-entry, 0, entry},
                                    {entry, -entry, 0}};

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    /// add for both u and v
                    int Xindex = index[i][j][0] * 2;
                    int Yindex = index[i][j][1] * 2;

                    SetValA(Xindex + 1, Yindex, -val[i][j]);
                    SetValA(Xindex, Yindex + 1, val[i][j]);
                }
        }

        /// set the diagonal of the matrix (which is zero at the beginning)
        /// as the sum of the other element inverted by sign
        void SetDiagonal(ScalarType val[3][3])
        {
            for (int i = 0; i < 3; i++)
            {
                ScalarType sum = 0;
                for (int j = 0; j < 3; j++)
                    sum += val[i][j];
                val[i][i] = -sum;
            }
        }

        /// add this values to the right hand side
        void AddRHS(ScalarType b[6],
                    int index[3])
        {
            for (int i = 0; i < 3; i++)
            {
                ScalarType valU = b[i * 2];
                ScalarType valV = b[(i * 2) + 1];
                SetValB((index[i] * 2), valU);
                SetValB((index[i] * 2) + 1, valV);
            }
        }

        /// add a 3x3 block matrix to the system matrix...
        /// indexes are specified in the 3x3 matrix of x,y pairs
        /// indexes must be multiplied by 2 cause u and v
        void Add33Block(ScalarType val[3][3], int index[3][3][2])
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    /// add for both u and v
                    int Xindex = index[i][j][0] * 2;
                    int Yindex = index[i][j][1] * 2;
                    assert(Xindex < int(n_vert_vars * 2));
                    assert(Yindex < int(n_vert_vars * 2));
                    SetValA(Xindex, Yindex, val[i][j]);
                    SetValA(Xindex + 1, Yindex + 1, val[i][j]);
                }
        }

        /// add a 3x3 block matrix to the system matrix...
        /// indexes are specified in the 3x3 matrix of x,y pairs
        /// indexes must be multiplied by 2 cause u and v
        void Add44Block(ScalarType val[4][4], int index[4][4][2])
        {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                {
                    /// add for both u and v
                    int Xindex = index[i][j][0] * 2;
                    int Yindex = index[i][j][1] * 2;
                    assert(Xindex < (n_vert_vars * 2));
                    assert(Yindex < (n_vert_vars * 2));
                    SetValA(Xindex, Yindex, val[i][j]);
                    SetValA(Xindex + 1, Yindex + 1, val[i][j]);
                }
        }

        // /// return the LHS for a given face
        // void perElementLHS(FaceType *f,
        //                    ScalarType val[3][3],
        //                    int index[3][3][2])
        // {
        //     /// initialize to zero
        //     for (int x = 0; x < 3; x++)
        //         for (int y = 0; y < 3; y++)
        //             val[x][y] = 0;

        //     /// get the vertices
        //     VertexType *v[3];
        //     v[0] = f->V(0);
        //     v[1] = f->V(1);
        //     v[2] = f->V(2);

        //     /// get the indexes of vertex instance (to consider cuts)
        //     /// for the current face
        //     int Vindexes[3];
        //     Vindexes[0] = VertexIndex(f->V(0));
        //     Vindexes[1] = VertexIndex(f->V(1));
        //     Vindexes[2] = VertexIndex(f->V(2));

        //     /// initialize the indexes for the block
        //     for (int x = 0; x < 3; x++)
        //         for (int y = 0; y < 3; y++)
        //         {
        //             index[x][y][0] = Vindexes[x];
        //             index[x][y][1] = Vindexes[y];
        //         }

        //     /// initialize edges
        //     CoordType e[3];
        //     for (int k = 0; k < 3; k++)
        //         e[k] = v[(k + 2) % 3]->P() - v[(k + 1) % 3]->P();

        //     /// then consider area but also considering scale factor dur to overlaps
        //     ScalarType areaT = ((f->P(1) - f->P(0)) ^ (f->P(2) - f->P(0))).Norm() / 2.0;
        //     for (int x = 0; x < 3; x++)
        //         for (int y = 0; y < 3; y++)
        //             if (x != y)
        //             {
        //                 ScalarType num = (e[x] * e[y]);
        //                 val[x][y] = num / (4.0 * areaT);
        //             }

        //     /// set the matrix as diagonal
        //     SetDiagonal(val);
        // }

        /// return the LHS for a given face
        void perElementLHS(const int IndexF,
                           ScalarType val[3][3],
                           int index[3][3][2])
        {
            /// initialize to zero
            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                    val[x][y] = 0;

            // /// get the vertices
            // VertexType *v[3];
            // v[0] = f->V(0);
            // v[1] = f->V(1);
            // v[2] = f->V(2);

            /// get the indexes of vertex instance (to consider cuts)
            /// for the current face
            int Vindexes[3];
            Vindexes[0] = Faces[IndexF][0]; // VertexIndex(f->V(0));
            Vindexes[1] = Faces[IndexF][1]; // VertexIndex(f->V(1));
            Vindexes[2] = Faces[IndexF][2]; // VertexIndex(f->V(2));

            /// initialize the indexes for the block
            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                {
                    index[x][y][0] = Vindexes[x];
                    index[x][y][1] = Vindexes[y];
                }

            // /// initialize edges
            // CoordType e[3];
            // for (int k = 0; k < 3; k++)
            //     e[k] = v[(k + 2) % 3]->P() - v[(k + 1) % 3]->P();

            /// initialize edges
            Geo::Point3<ScalarType> e[3];
            for (int k = 0; k < 3; k++)
            {
                // e[k] = v[(k + 2) % 3]->P() - v[(k + 1) % 3]->P();
                int V0 = Vindexes[(k + 2) % 3];
                int V1 = Vindexes[(k + 1) % 3];
                Geo::Point3<ScalarType> P0 = VertPos[V0];
                Geo::Point3<ScalarType> P1 = VertPos[V1];
                e[k] = P0 - P1;
            }

            /// then consider area but also considering scale factor dur to overlaps
            // ScalarType areaT = ((f->P(1) - f->P(0)) ^ (f->P(2) - f->P(0))).Norm() / 2.0;
            ScalarType areaT = Area(VertPos, Faces, IndexF);

            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                    if (x != y)
                    {
                        //ScalarType num = (e[x] * e[y]);
                        ScalarType num = Geo::Dot(e[x],e[y]);
                        val[x][y] = num / (4.0 * areaT);
                    }

            /// set the matrix as diagonal
            SetDiagonal(val);
        }

        // /// return the RHS for a given face
        // void perElementRHS(FaceType *f,
        //                    ScalarType b[6],
        //                    ScalarType vector_field_scale = 1)
        // {

        //     /// then set the rhs
        //     CoordType scaled_Kreal;
        //     CoordType scaled_Kimag;
        //     CoordType fNorm = f->N();
        //     fNorm.Normalize();
        //     CoordType p[3];
        //     p[0] = f->P0(0);
        //     p[1] = f->P0(1);
        //     p[2] = f->P0(2);

        //     CoordType neg_t[3];
        //     neg_t[0] = fNorm ^ (p[2] - p[1]);
        //     neg_t[1] = fNorm ^ (p[0] - p[2]);
        //     neg_t[2] = fNorm ^ (p[1] - p[0]);

        //     CoordType K1 = CoordType::Construct(f->PD1());
        //     CoordType K2 = CoordType::Construct(f->PD2());
        //     K1.Normalize();
        //     K2.Normalize();

        //     scaled_Kreal = K1 * (vector_field_scale); /// 2);
        //     scaled_Kimag = K2 * (vector_field_scale); /// 2);

        //     b[0] = scaled_Kreal * neg_t[0];
        //     b[1] = scaled_Kimag * neg_t[0];
        //     b[2] = scaled_Kreal * neg_t[1];
        //     b[3] = scaled_Kimag * neg_t[1];
        //     b[4] = scaled_Kreal * neg_t[2];
        //     b[5] = scaled_Kimag * neg_t[2];
        // }

        /// return the RHS for a given face
        void perElementRHS(const int IndexF,
                           ScalarType b[6],
                           ScalarType vector_field_scale = 1)
        {

            /// then set the rhs
            Geo::Point3<ScalarType> scaled_Kreal;
            Geo::Point3<ScalarType> scaled_Kimag;
            // CoordType fNorm = f->N();
            Geo::Point3<ScalarType> fNorm = Geo::Normal<ScalarType>(VertPos, Faces, IndexF);
            fNorm.Normalize();

            // CoordType p[3];
            // p[0] = f->P0(0);
            // p[1] = f->P0(1);
            // p[2] = f->P0(2);
            Geo::Point3<ScalarType> p[3];
            int IndexV0 = Faces[IndexF][0];
            int IndexV1 = Faces[IndexF][1];
            int IndexV2 = Faces[IndexF][2];
            p[0] = VertPos[IndexV0];
            p[1] = VertPos[IndexV1];
            p[2] = VertPos[IndexV2];

            // CoordType neg_t[3];
            // neg_t[0] = fNorm ^ (p[2] - p[1]);
            // neg_t[1] = fNorm ^ (p[0] - p[2]);
            // neg_t[2] = fNorm ^ (p[1] - p[0]);
            Geo::Point3<ScalarType> neg_t[3];
            neg_t[0] = Geo::Cross(fNorm, (p[2] - p[1]));
            neg_t[1] = Geo::Cross(fNorm, (p[0] - p[2]));
            neg_t[2] = Geo::Cross(fNorm, (p[1] - p[0]));

            // CoordType K1 = CoordType::Construct(f->PD1());
            // CoordType K2 = CoordType::Construct(f->PD2());
            // K1.Normalize();
            // K2.Normalize();
            Geo::Point3<ScalarType> K1 = FaceCurv[IndexF].D[0];
            Geo::Point3<ScalarType> K2 = FaceCurv[IndexF].D[1];
            K1.Normalize();
            K2.Normalize();

            scaled_Kreal = K1 * (vector_field_scale); /// 2);
            scaled_Kimag = K2 * (vector_field_scale); /// 2);

            // b[0] = scaled_Kreal * neg_t[0];
            // b[1] = scaled_Kimag * neg_t[0];
            // b[2] = scaled_Kreal * neg_t[1];
            // b[3] = scaled_Kimag * neg_t[1];
            // b[4] = scaled_Kreal * neg_t[2];
            // b[5] = scaled_Kimag * neg_t[2];
            b[0] = Geo::Dot(scaled_Kreal, neg_t[0]);
            b[1] = Geo::Dot(scaled_Kimag, neg_t[0]);
            b[2] = Geo::Dot(scaled_Kreal, neg_t[1]);
            b[3] = Geo::Dot(scaled_Kimag, neg_t[1]);
            b[4] = Geo::Dot(scaled_Kreal, neg_t[2]);
            b[5] = Geo::Dot(scaled_Kimag, neg_t[2]);
        }

        // /// return the LHS and RHS for a given face
        // void PerElementSystemReal(FaceType *f,
        //                           ScalarType val[3][3],
        //                           int index[3][3][2],
        //                           ScalarType b[6],
        //                           ScalarType vector_field_scale = 1.0)
        // {
        //     perElementLHS(f, val, index);

        //     if (use_direction_field)
        //         perElementRHS(f, b, vector_field_scale);
        // }

        /// return the LHS and RHS for a given face
        void PerElementSystemReal(const int IndexF,
                                  ScalarType val[3][3],
                                  int index[3][3][2],
                                  ScalarType b[6],
                                  ScalarType vector_field_scale = 1.0)
        {
            // perElementLHS(f, val, index);
            perElementLHS(IndexF, val, index);

            if (use_direction_field)
                // perElementRHS(f, b, vector_field_scale);
                perElementRHS(IndexF, b, vector_field_scale);
        }

        // void FixPointLSquares()
        // {
        //     ScalarType penalization = 1000000;
        //     int offset_row = n_vert_vars;
        //     assert(to_fix.size() > 0);
        //     for (size_t i = 0; i < to_fix.size(); i++)
        //     {
        //         /// take a vertex
        //         VertexType *v = to_fix[i];
        //         assert(!v->IsD());
        //         int index = VertexIndex(v);
        //         // v->vertex_index[0];
        //         int indexvert = index * 2;
        //         int indexRow = (offset_row + i) * 2;

        //         SetValA(indexRow, indexRow, penalization);
        //         SetValA(indexRow + 1, indexRow + 1, penalization);

        //         /// add values to the B vector
        //         ScalarType U = v->T().U() * penalization;
        //         ScalarType V = v->T().V() * penalization;
        //         SetValB(indexRow, U);
        //         SetValB(indexRow + 1, V);

        //         /*///set upper right part
        //         SetValA(indexvert,indexCol,penalization);
        //         SetValA(indexvert+1,indexCol+1,penalization);*/

        //         SetValA(indexvert, indexvert, penalization);
        //         SetValA(indexvert + 1, indexvert + 1, penalization);
        //         SetValA(indexRow, indexRow, penalization);
        //         SetValA(indexRow + 1, indexRow + 1, penalization);
        //         SetValA(indexvert, indexRow, -penalization);
        //         SetValA(indexvert + 1, indexRow + 1, -penalization);
        //         SetValA(indexRow, indexvert, -penalization);
        //         SetValA(indexRow + 1, indexvert + 1, -penalization);
        //         // SetValA(indexCol+1,indexCol+1,-1);
        //     }
        // }

        void FixPointLSquares()
        {
            ScalarType penalization = 1000000;
            int offset_row = n_vert_vars;
            assert(to_fix.size() > 0);
            for (size_t i = 0; i < to_fix.size(); i++)
            {
                /// take a vertex
                // VertexType *v = to_fix[i];
                // assert(!v->IsD());
                // int index = VertexIndex(v);

                int index = to_fix[i];

                // v->vertex_index[0];
                int indexvert = index * 2;
                int indexRow = (offset_row + i) * 2;

                SetValA(indexRow, indexRow, penalization);
                SetValA(indexRow + 1, indexRow + 1, penalization);

                /// add values to the B vector
                // ScalarType U = v->T().U() * penalization;
                // ScalarType V = v->T().V() * penalization;
                ScalarType U = to_fixUV[i].X * penalization;
                ScalarType V = to_fixUV[i].Y * penalization;

                SetValB(indexRow, U);
                SetValB(indexRow + 1, V);

                /*///set upper right part
                SetValA(indexvert,indexCol,penalization);
                SetValA(indexvert+1,indexCol+1,penalization);*/

                SetValA(indexvert, indexvert, penalization);
                SetValA(indexvert + 1, indexvert + 1, penalization);
                SetValA(indexRow, indexRow, penalization);
                SetValA(indexRow + 1, indexRow + 1, penalization);
                SetValA(indexvert, indexRow, -penalization);
                SetValA(indexvert + 1, indexRow + 1, -penalization);
                SetValA(indexRow, indexvert, -penalization);
                SetValA(indexRow + 1, indexvert + 1, -penalization);
                // SetValA(indexCol+1,indexCol+1,-1);
            }
        }

        // // build the laplacian matrix cyclyng over all rangemaps
        // // and over all faces
        // void BuildLaplacianMatrix(double vfscale = 1)
        // {

        //     /// then for each face
        //     for (unsigned int j = 0; j < mesh.face.size(); j++)
        //     {

        //         FaceType *f = &mesh.face[j];
        //         if (f->IsD())
        //             continue;

        //         int var_idx[3]; // vertex variable indices
        //         for (int k = 0; k < 3; ++k)
        //         {
        //             VertexType *v = f->V(k);
        //             var_idx[k] = VertexIndex(v);
        //         }
        //         ScalarType val[3][3];
        //         int index[3][3][2];
        //         ScalarType b[6];
        //         PerElementSystemReal(f, val, index, b, vfscale);

        //         // Add the element to the matrix
        //         Add33Block(val, index);

        //         /////add area term.. to test if needed
        //         /*if (!use_direction_field)
        //             AddAreaTerm(index,1.0);//f->area);*/
        //         /*ScalarType area=((f->P(1)-f->P(0))^(f->P(2)-f->P(0))).Norm();
        //         if (!use_direction_field)
        //             AddAreaTerm(index,area);*/

        //         // ScalarType area=((f->P(1)-f->P(0))^(f->P(2)-f->P(0))).Norm();
        //         if (!use_direction_field)
        //             AddAreaTerm(index, 1);

        //         /// add right hand side
        //         if (use_direction_field)
        //             AddRHS(b, var_idx);
        //     }
        // }

        // build the laplacian matrix cyclyng over all rangemaps
        // and over all faces
        void BuildLaplacianMatrix(double vfscale = 1)
        {

            /// then for each face
            // for (unsigned int j = 0; j < mesh.face.size(); j++)
            for (unsigned int j = 0; j < Faces.size(); j++)
            {

                // FaceType *f = &mesh.face[j];
                // if (f->IsD())
                //     continue;

                int var_idx[3]; // vertex variable indices
                for (int k = 0; k < 3; ++k)
                {
                    // VertexType *v = f->V(k);
                    // var_idx[k] = VertexIndex(v);
                    var_idx[k] = Faces[j][k];
                }
                ScalarType val[3][3];
                int index[3][3][2];
                ScalarType b[6];
                // PerElementSystemReal(f, val, index, b, vfscale);
                PerElementSystemReal(j, val, index, b, vfscale);

                // Add the element to the matrix
                Add33Block(val, index);

                /////add area term.. to test if needed
                /*if (!use_direction_field)
                    AddAreaTerm(index,1.0);//f->area);*/
                /*ScalarType area=((f->P(1)-f->P(0))^(f->P(2)-f->P(0))).Norm();
                if (!use_direction_field)
                    AddAreaTerm(index,area);*/

                // ScalarType area=((f->P(1)-f->P(0))^(f->P(2)-f->P(0))).Norm();
                if (!use_direction_field)
                    AddAreaTerm(index, 1);

                /// add right hand side
                if (use_direction_field)
                    AddRHS(b, var_idx);
            }
        }

        void FindSizes()
        {
            // tag vertices and compute numbers of equations to determine the number of rows in the matrix
            // TagVertices_Constrained(n_vert_vars, n_transition_eqs, n_align_sharp_eqs);
            n_vert_vars = VertPos.size(); // mesh.vn;

            /// initialize matrix size
            total_size = (n_fixed_vars + n_vert_vars) * 2; /// must be multiplied by 2 becasue of u and v
        }

        void AllocateSystem()
        {
            //--- Allocates the data for Ax=b
            A = Eigen::SparseMatrix<double>(total_size, total_size); // A
            A.reserve(Eigen::VectorXi::Constant(total_size, 32));    // This prealloaction trick greatly speed up the acc
            b = Eigen::VectorXd::Zero(total_size);                   // x and b
        }

        /// intitialize the whole matrix
        void InitMatrix()
        {
            FindSizes();
            AllocateSystem();
        }

        bool Solve()
        {
            // return true;
            A.finalize();
            Eigen::SparseMatrix<double> As = Eigen::SparseMatrix<double>(A);
            As.finalize();

            Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(As);
            x = solver.solve(b);
            return (solver.info() == Eigen::Success);
        }

        // void InitIndex()
        // {
        //     for (size_t i = 0; i < mesh.vert.size(); i++)
        //         if (!mesh.vert[i].IsD())
        //             AddVertexIndex(&mesh.vert[i], i);
        // }

        /// map back values to vertex
        /// if normalize==true then set the
        /// coordinates between 0 and 1
        void MapCoords(std::vector<Geo::Point2<ScalarType>> &UVPos,
                       bool normalize = true,
                       ScalarType fieldScale = 1.0)
        {
            // /// clear Visited Flag
            // if (correct_fixed)
            //     tri::UpdateFlags<MeshType>::VertexClearV(mesh);
            // // set fixed to V
            // for (size_t i = 0; i < to_fix.size(); i++)
            //     to_fix[i]->SetV();

            Box2<ScalarType> bbox;
            if (normalize)
            {
                for (size_t i = 0; i < n_vert_vars; i++)
                {
                    ScalarType U = x[i * 2];
                    ScalarType V = x[(i * 2) + 1];
                    bbox.Add(Point2<ScalarType>(U, V));
                }
            }

            // for each vertex
            for (size_t i = 0; i < n_vert_vars; i++)
            {
                // VertexType *v = IndexVertex(i);
                //  take U and V
                ScalarType U = x[i * 2];
                ScalarType V = x[(i * 2) + 1];
                Point2<ScalarType> p;
                // if (!v->IsV())
                p = Point2<ScalarType>(U, V);
                // else
                //     p = v->T().P();

                // p/=fieldScale;
                if (normalize)
                {
                    //p -= bbox.min;
                    p -= bbox.Min;
                    p = p * 1 / bbox.Diag();
                }

                // v->T().P() = p;
                UVPos[i] = p;
            }

            // /// then copy to faces
            // for (size_t i = 0; i < mesh.face.size(); i++)
            // {
            //     FaceType *f = &mesh.face[i];
            //     for (int j = 0; j < 3; j++)
            //     {
            //         VertexType *v = f->V(j);
            //         Point2<ScalarType> p = v->T().P();
            //         f->WT(j).P() = p;
            //     }
            // }
        }

    public:
     /// vertices to fix
        // std::vector<VertexType *> to_fix;
        std::vector<size_t> to_fix;
        std::vector<Geo::Point2<ScalarType>> to_fixUV;

        /// return true if is possible to
        bool IsFeasible()
        {
            return true;
            //assert(0);
            // tri::UpdateTopology<MeshType>::FaceFace(mesh);
            // int NNmanifoldE = tri::Clean<MeshType>::CountNonManifoldEdgeFF(mesh);
            // if (NNmanifoldE != 0)
            // {
            //     printf("Non Manifold Edges \n");
            //     return false;
            // }
            // int NNmanifoldV = tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
            // if (NNmanifoldV != 0)
            // {
            //     printf("Non Manifold Vertices \n");
            //     return false;
            // }
            // int H = tri::Clean<MeshType>::CountHoles(mesh);
            // if (H == 0)
            //     return false;

            // int G = tri::Clean<MeshType>::MeshGenus(mesh);
            // if (G != 0)
            // {
            //     printf("Genus %d\n", G);
            //     return false;
            // }

            // return (true);
        }

        // /// set the border as fixed
        // void SetBorderAsFixed()
        // {
        //     for (size_t i = 0; i < mesh.vert.size(); i++)
        //     {
        //         VertexType *v = &mesh.vert[i];
        //         if (v->IsD())
        //             continue;
        //         if (v->IsB())
        //             to_fix.push_back(v);
        //     }
        //     std::sort(to_fix.begin(), to_fix.end());
        //     auto new_end = std::unique(to_fix.begin(), to_fix.end());
        //     to_fix.resize(distance(to_fix.begin(), new_end));
        // }

        // /// set selected vertices as fixed
        // void SetSelectedAsFixed()
        // {
        //     for (int i = 0; i < mesh.vert.size(); i++)
        //     {
        //         VertexType *v = &mesh.vert[i];
        //         if (v->IsD())
        //             continue;
        //         if (v->IsS())
        //             to_fix.push_back(v);
        //     }
        //     std::sort(to_fix.begin(), to_fix.end());
        //     typename std::vector<VertexType *>::iterator new_end = std::unique(to_fix.begin(), to_fix.end());
        //     int dist = distance(to_fix.begin(), new_end);
        //     to_fix.resize(dist);
        // }

        /// intialize parameters and setup fixed vertices vector
        void Init(bool _use_direction_field = false,
                  bool _correct_fixed = true,
                  ScalarType _fieldScale = 1.0)
        {
            use_direction_field = _use_direction_field;
            correct_fixed = _correct_fixed;
            fieldScale = _fieldScale;
            to_fix.clear();
        }

        /// solve the system, it return false if the matrix is singular
        bool SolvePoisson(std::vector<Geo::Point2<ScalarType>> &UVPos,
                          bool _write_messages = false,
                          ScalarType fieldScale = 1.0)
        {
            int t0, t1, t2, t3;

            /// Initializing Matrix
            if (_write_messages)
            {
                printf("\n INITIALIZING THE MATRIX \n");
                t0 = clock();
            }

            /// set vertex indexes
            // InitIndex();

            if (use_direction_field)
            {
                assert(to_fix.size() > 0);
            }
            else
            {
                assert(to_fix.size() > 1);
            }

            n_fixed_vars = to_fix.size();
            /// initialize the matrix ALLOCATING SPACE
            InitMatrix();

            /// build the laplacian system
            BuildLaplacianMatrix(fieldScale);

            ////add the lagrange multiplier
            FixPointLSquares();

            if (_write_messages)
            {
                t1 = clock();
                printf("\n time:%d \n", t1 - t0);
                printf("\n SOLVING \n");
            }

            // int n_vars=(n_vert_vars)*2;
            // int integer_constr_size=(n_transition_vars+n_fixed_vars+n_bary_transition_vars)*2;
            // X=std::vector< double >(n_vars+n_fixed_vars*2);
            bool done = Solve();
            if (!done)
                return false;
            if (_write_messages)
            {
                t2 = clock();
                printf("\n time:%d \n", t2 - t1);
                printf("\n ASSIGNING COORDS \n");
            }

            UVPos.clear();
            UVPos.resize(VertPos.size());
            MapCoords(UVPos, false, fieldScale);
            if (_write_messages)
            {
                t3 = clock();
                printf("\n time:%d \n", t3 - t2);
            }

            /// then check if majority of faces are folded
            // if (!solve_global_fold)
            //     return true;
            // if (tri::Distortion<MeshType, false>::GloballyUnFolded(mesh))
            // {
            //     tri::UV_Utils<MeshType>::GloballyMirrorX(mesh);
            //     bool isUnfolded = tri::Distortion<MeshType, false>::GloballyUnFolded(mesh);
            //     assert(!isUnfolded);
            // }
            return true;
        }

        // PoissonSolver(MeshType &_mesh) : mesh(_mesh)
        // {
        //     assert(mesh.vert.size() >= 3);
        //     assert(mesh.face.size() >= 1);
        // }
        UVPoissonSolver(const std::vector<Geo::Point3<ScalarType>> &_VertPos,
                        const std::vector<std::vector<int>> &_Faces,
                        const std::vector<Field::CrossF<ScalarType>> &_FaceCurv) : VertPos(_VertPos), Faces(_Faces), FaceCurv(_FaceCurv)
        {
        }
    };

}

#endif
