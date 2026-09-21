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

#ifndef LAPLACIAN_SOLVE
#define LAPLACIAN_SOLVE

#include "point3.h"
#include "connectivity.h"
#include <Eigen/Sparse>
#include <set>
#include "solve/linear_constraint.h"
#include "solve/least_squares_lin_constraint_solver.h"
#include <cotangent.h>
#include <mesh_functions.h>

#define PENALTY_FIXED 1000000

namespace Geo
{

template <class ScalarType>
void FindLaplacianForMesh(const std::vector<Point3<ScalarType> > &pos,
                          const std::vector<std::vector<int> > &Kernel,
                          std::vector< std::vector<Point3<ScalarType> > > &Laplacian)
{
    typedef Point3<ScalarType> CoordType;

    //check symmetry
    for (size_t i=0;i<Kernel.size();i++)
    {
        int IndexV0=i;
        for (size_t j=0;j<Kernel[i].size();j++)
        {
            int IndexV1=Kernel[i][j];
            assert(std::find(Kernel[IndexV1].begin(), Kernel[IndexV1].end(), IndexV0) != Kernel[IndexV1].end());
        }
    }

    Laplacian.clear();
    Laplacian.resize(pos.size());
    for (size_t i=0;i<Kernel.size();i++)
    {
        CoordType PVert=pos[i];
        for (size_t j=0;j<Kernel[i].size();j++)
        {
            CoordType PNeigh=pos[Kernel[i][j]];
            Laplacian[i].push_back(PVert-PNeigh);
        }
    }
}

// template <class ScalarType>
// void SetUniformLaplacianWeightForMesh(const std::vector<std::vector<int> > &Kernel,
//                                       std::vector<std::vector<ScalarType> > &Weight)
// {
//     Weight.clear();
//     Weight.resize(Kernel.size());
//     for (size_t i=0;i<Weight.size();i++)
//     {
//         Weight[i]=std::vector<ScalarType>(Kernel[i].size(),0.5);//((ScalarType)1)/((ScalarType)Kernel[i].size()));
//     }
// }

template <class ScalarType>
void SetUniformLaplacianWeightForMesh(const std::vector<std::vector<int> > &Kernel,
                                      std::vector<std::vector<ScalarType> > &Weight)
{
    Weight.clear();
    Weight.resize(Kernel.size());
    for (size_t i=0;i<Weight.size();i++)
    {
        if (Kernel[i].size()>0)
        {
            ScalarType Value=((ScalarType)1)/((ScalarType)Kernel[i].size());
            Weight[i]=std::vector<ScalarType>(Kernel[i].size(),Value);
        }
    }
}

template <class ScalarType>
void SetUniformLaplacianWeightForMesh(const std::vector<std::vector<int> > &Connectivity,
                                      const size_t NumVert,
                                      std::vector<std::vector<ScalarType> > &Weight)
{
    std::vector<std::vector<int> > Kernel;
    Geo::FindVertexStar(Connectivity,NumVert,Kernel);
    Geo::SetUniformLaplacianWeightForMesh(Kernel,Weight);
}

template <class ScalarType>
void FindLaplacianCotangentWeightForTriMesh(const std::vector<Point3<ScalarType> > &VertPos,
                                            const std::vector<std::vector<int> > &Faces,
                                            const std::vector<std::vector<int> > &Kernel,
                                            std::vector<std::vector<ScalarType> > &Weight,
                                            bool check_angle=true)
{
    std::vector<std::vector<int> > NextF,NextE;
    getFFAdjacency(Faces,NextF,NextE);
    std::vector<std::vector<ScalarType> > CotWeight;
    CotangentWeightTriMesh(VertPos,Faces,NextF,NextE,CotWeight,check_angle);

    std::map<std::pair<int,int>, ScalarType> EdgeWMap;
    for (size_t i=0;i<Faces.size();i++)
    {
        assert(Faces[i].size()==3);
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV0=Faces[i][j];
            int IndexV1=Faces[i][(j+1)%3];
            std::pair<int,int> Entry(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            EdgeWMap[Entry]=CotWeight[i][j];
        }
    }

    Weight.clear();
    Weight.resize(Kernel.size());
    for (size_t i=0;i<Kernel.size();i++)
        for (size_t j=0;j<Kernel[i].size();j++)
        {
            int IndexV0=i;
            int IndexV1=Kernel[i][j];
            std::pair<int,int> Entry(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            assert(EdgeWMap.count(Entry)>0);
            Weight[i].push_back(EdgeWMap[Entry]);
        }
}

template <class ScalarType>
void FindLaplacianCotangentWeightForTetMesh(const std::vector<Point3<ScalarType> > &VertPos,
                                            const std::vector<std::vector<int> > &Tetra,
                                            const std::vector<std::vector<int> > &Kernel,
                                            std::vector<std::vector<ScalarType> > &Weight)
{
    std::vector<std::vector<int> > NextT,NextF;
    Tetra::GetPPAdjacency(Tetra,NextT,NextF);

    std::map<std::pair<int,int>, ScalarType > EdgeWMap;
    CotangentWeightTetraMesh(VertPos,Tetra,NextT,NextF,EdgeWMap);
    Weight.resize(Kernel.size());
    for (size_t i=0;i<Kernel.size();i++)
        for (size_t j=0;j<Kernel[i].size();j++)
        {
            int IndexV0=i;
            int IndexV1=Kernel[i][j];
            std::pair<int,int> Entry(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            assert(EdgeWMap.count(Entry)>0);
            Weight[i].push_back(EdgeWMap[Entry]);
        }
}

template <class ScalarType>
void FindLaplacianCotangentWeight(const std::vector<Point3<ScalarType> > &VertPos,
                                  const std::vector<std::vector<int> > &Connections,
                                  const std::vector<std::vector<int> > &Kernel,
                                  std::vector<std::vector<ScalarType> > &Weights,
                                  MType TypeMesh)
{
    switch (TypeMesh)
    {
    case MT_Triangle:
        //std::cout<<"Using cotangent Tris"<<std::endl;
        FindLaplacianCotangentWeightForTriMesh(VertPos,Connections,Kernel,Weights);
        break;

    case MT_Tetra:
        //std::cout<<"Using cotangent Tetra"<<std::endl;
        FindLaplacianCotangentWeightForTetMesh(VertPos,Connections,Kernel,Weights);
        break;

    default:
        //std::cout<<"No Cotangent"<<std::endl;
        SetUniformLaplacianWeightForMesh(Kernel,Weights);
        // code to be executed if
        // expression doesn't match any constant
    }
}

template <class ScalarType>
void FindLaplaciaDataForMesh(const std::vector<Point3<ScalarType> > &pos,
                             const std::vector<std::vector<int> > &connections,
                             std::vector<std::vector<int> > &Kernel,
                             std::vector< std::vector<Point3<ScalarType> > > &Laplacian,
                             std::vector<std::vector<ScalarType> > &Weights,
                             MType TypeMesh,
                             bool UseCotangentWeights)
{
    FindVertexStar(connections,pos.size(),Kernel);

    FindLaplacianForMesh(pos,Kernel,Laplacian);

    if (!UseCotangentWeights)
        SetUniformLaplacianWeightForMesh(Kernel,Weights);
    else
    {
        FindLaplacianCotangentWeight(pos,connections,Kernel,Weights,TypeMesh);
    }
}

template <class ScalarType>
void RemoveLaplacianAlongNormal(const std::vector<Point3<ScalarType> > &pos,
                                const std::vector<std::vector<int> > &connections,
                                std::vector< std::vector<Point3<ScalarType> > > &Laplacian)
{
    std::vector<Point3<ScalarType> > perFaceNormals,perVertexNormals;
    Geo::ComputeNormals(pos,connections,perFaceNormals,perVertexNormals);
    for (size_t i=0;i<Laplacian.size();i++)
    {
        for (size_t j=0;j<Laplacian[i].size();j++)
        {
            ScalarType DotVal=Geo::Dot<ScalarType>(Laplacian[i][j],perVertexNormals[i]);
            Laplacian[i][j]-=(perFaceNormals[i]*DotVal);
        }
    }
}

template <class ScalarType>
class LaplacianSolve
{
    typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;
    typedef Point3<ScalarType> CoordType;

    Eigen::SparseMatrix<ScalarType> L,L0;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType > > *solver;

    std::vector<std::vector<int> > Kernels;
    std::vector<std::vector<ScalarType> > Weights;
    std::vector<int> FixedVert;

    std::vector<LinearConstraint<ScalarType> > LinConstr;

    void InitLaplacianMatrix()
    {
        std::vector<std::pair<int,int> > index;
        std::vector<ScalarType> entry;

        //INITIALIZE Laplacian
        for (size_t i = 0; i < Kernels.size(); i++)
            for (size_t j = 0; j < Kernels[i].size(); j++)
            {
                ScalarType wi = Weights[i][j];
                int IndexV0=(i*3);
                int IndexV1=(Kernels[i][j]*3);
                assert(IndexV0!=IndexV1);
                for (size_t k=0;k<3;k++)
                {
                    index.push_back(std::pair<int,int>(IndexV0+k,IndexV0+k));
                    entry.push_back(wi);
                    index.push_back(std::pair<int,int>(IndexV0+k,IndexV1+k));
                    entry.push_back(-wi);
                }
            }

        size_t nextRow=Kernels.size()*3;
        for (size_t i=0;i<LinConstr.size();i++)
        {
            int CurrRow=nextRow+i;
            for (size_t j=0;j<LinConstr[i].LHSide.size();j++)
            {
                size_t IndexVar=LinConstr[i].LHSide[j].IndexVar;
                size_t IndexCoord=LinConstr[i].LHSide[j].IndexCoord;
                size_t VarCoord=IndexVar*3+IndexCoord;
                ScalarType W=LinConstr[i].LHSide[j].W;
                if (W==0)continue;
                index.push_back(std::pair<int,int>(CurrRow,VarCoord));
                entry.push_back(W*LinConstr[i].Weight);
            }
        }

        //int sizeM=Kernels.size()*3+FixedVert.size()*3+LinConstr.size();
        int sizeM=Kernels.size()*3+LinConstr.size();
        InitSparse(index,entry,sizeM,Kernels.size()*3,L0);

        L=L0.transpose()*L0;
    }

    void AddFixedConstraints()
    {
        startFixedConstr=LinConstr.size();
        for (size_t i=0;i<FixedVert.size();i++)
        {
            int IndexV=FixedVert[i];
            for (size_t j=0;j<3;j++)
            {
                Geo::LinearConstraint<ScalarType> Constr;
                Constr.Weight=PENALTY_FIXED;
                Geo::VarTerm<ScalarType> Var(IndexV,j,1);
                Constr.LHSide.push_back(Var);
                LinConstr.push_back(Constr);
            }
        }
    }

    void UpdateFixedConstraints(const std::vector<CoordType> &FixPos)
    {
        for (size_t i=0;i<FixedVert.size();i++)
        {
            //int IndexV=FixedVert[i];
            for (size_t j=0;j<3;j++)
            {
                int CurrIdx=startFixedConstr+(i*3)+j;
                LinConstr[CurrIdx].RHSide=FixPos[i].cV(j);
            }
        }
    }

    int startFixedConstr;

public:

    int NumVariables()
    {
        assert(Kernels.size()==Weights.size());
        //        for (size_t i=0;i<FixedVert.size();i++)
        //        {assert(FixedVert[i]<Kernels.size());}
        return Kernels.size();
    }

    void Init(const std::vector<std::vector<int> > &_Kernels,
              const std::vector<std::vector<ScalarType> > &_Weights,
              const std::vector<int> &_FixedVert,
              const std::vector<LinearConstraint<ScalarType> > &_LinConstr=std::vector<LinearConstraint<ScalarType> >(),
              bool WriteDBG=false)
    {
        Kernels=_Kernels;
        Weights=_Weights;
        FixedVert=_FixedVert;

        LinConstr=_LinConstr;

        AddFixedConstraints();

        assert(Kernels.size()==Weights.size());

        InitLaplacianMatrix();

        solver= new Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType > > (L);
        assert(solver->info() == Eigen::Success);

        if (WriteDBG)
            std::cout<<"SUCCESS PREFACTORIZATION"<<std::endl;
    }

    void SolveLaplacian(const std::vector<std::vector<CoordType> > &LaplValues,
                        const std::vector<CoordType> &FixPos,
                        std::vector<CoordType> &VertPos)
    {

        //then get the rhs
        //int sizeM=Kernels.size()*3+FixedVert.size()*3+LinConstr.size();
        int sizeM=Kernels.size()*3+LinConstr.size();
        MatrixXm V;
        //V=MatrixXm(sizeM*3,1);
        V=MatrixXm(sizeM,1);
        V.setZero();

        for (size_t i=0;i<LaplValues.size();i++)
        {
            CoordType CurrL(0,0,0);
            for (size_t j=0;j<LaplValues[i].size();j++)
                CurrL+=LaplValues[i][j]*Weights[i][j];

            //CurrL/=LaplValues.size();
            int IndexVX=i*3;
            int IndexVY=i*3+1;
            int IndexVZ=i*3+2;

            V(IndexVX,0)+=(CurrL.X);
            V(IndexVY,0)+=(CurrL.Y);
            V(IndexVZ,0)+=(CurrL.Z);
        }

        //update RHS for fixed ones
        UpdateFixedConstraints(FixPos);

        size_t nextRow=(Kernels.size())*3;
        for (size_t i=0;i<LinConstr.size();i++)
        {
            int CurrRow=nextRow+i;
            V(CurrRow)+=LinConstr[i].RHSide*LinConstr[i].Weight;
        }

        V=L0.transpose()*V;
        //solve the system
        V = solver->solve(V).eval();

        //update the position
        for (size_t i=0;i<VertPos.size();i++)
        {
            int indexV=(i*3);

            VertPos[i].X=V(indexV,0);
            VertPos[i].Y=V(indexV+1,0);
            VertPos[i].Z=V(indexV+2,0);
        }
    }

    LaplacianSolve()
    {
        startFixedConstr=0;
        solver=NULL;
    }

    ~LaplacianSolve()
    {
        if (solver!=NULL)
            delete (solver);
    }
};
}
#endif
