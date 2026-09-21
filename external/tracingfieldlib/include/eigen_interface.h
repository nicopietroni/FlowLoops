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

#ifndef EIGEN_INTERFACE
#define EIGEN_INTERFACE

#include "point3.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <cmath>
#include <stdio.h>
#include <vector>
// #include "mesh_functions.h"

// return mesh as vector of vertices and faces
template <class ScalarType> class EigenInterface {

public:
  typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
      MatrixXm;

  static void GetVertData(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                          MatrixXm &vert) {
    // create eigen matrix of vertices
    vert = MatrixXm(VertPos.size(), 2);

    // copy vertices
    for (int i = 0; i < (int)VertPos.size(); i++)
      for (int j = 0; j < 2; j++)
        vert(i, j) = VertPos[i].cV(j);
  }

  static void GetVertData(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          MatrixXm &vert) {
    // create eigen matrix of vertices
    vert = MatrixXm(VertPos.size(), 3);

    // copy vertices
    for (int i = 0; i < (int)VertPos.size(); i++)
      for (int j = 0; j < 3; j++)
        vert(i, j) = VertPos[i].cV(j);
  }

  static void
  GetTriMeshData(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 Eigen::MatrixXi &faces, MatrixXm &vert) {
    // assert(!HasUnfererencedVert<ScalarType>(VertPos, Faces));

    // create eigen matrix of vertices
    vert = MatrixXm(VertPos.size(), 3);

    // copy vertices
    for (int i = 0; i < (int)VertPos.size(); i++)
      for (int j = 0; j < 3; j++)
        vert(i, j) = VertPos[i].cV(j);

    // create eigen matrix of faces
    faces = Eigen::MatrixXi(Faces.size(), 3);

    // copy faces
    for (int i = 0; i < (int)Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (int j = 0; j < 3; j++)
        faces(i, j) = Faces[i][j];
    }
  }

  static void
  GetTriMeshData(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 Eigen::MatrixXi &faces, MatrixXm &vert) {
    assert(!HasUnfererencedVert(VertPos, Faces));

    // create eigen matrix of vertices
    vert = MatrixXm(VertPos.size(), 2);

    // copy vertices
    for (int i = 0; i < (int)VertPos.size(); i++)
      for (int j = 0; j < 2; j++)
        vert(i, j) = VertPos[i].cV(j);

    // create eigen matrix of faces
    faces = Eigen::MatrixXi(Faces.size(), 3);

    // copy faces
    for (int i = 0; i < (int)Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (int j = 0; j < 3; j++)
        faces(i, j) = Faces[i][j];
    }
  }

  static void
  GetTriMeshData(const Eigen::MatrixXi &faces, 
                 const MatrixXm &vert,
                 std::vector<Geo::Point2<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &Faces) {
    
   VertPos.resize(vert.rows());
    // copy vertices
    for (int i = 0; i < (int)vert.rows(); i++)
      for (int j = 0; j < 2; j++)
        VertPos[i].V(j)=vert(i, j);

    Faces=std::vector<std::vector<int>>(faces.rows(),
                                     std::vector<int>(3,-1));          
    // copy faces
    for (int i = 0; i < (int)faces.rows(); i++) {
      assert(Faces[i].size() == 3);
      for (int j = 0; j < 3; j++)
        Faces[i][j]=faces(i, j);
    }
  }

  static void GetTriMeshData(const Eigen::MatrixXi &eigen_f,
                             const MatrixXm &eigen_v,
                             std::vector<Geo::Point3<ScalarType>> &VertPos,
                             std::vector<std::vector<int>> &Faces) {
    Faces.clear();
    VertPos.clear();

    for (int i = 0; i < eigen_v.rows(); ++i) {
      Geo::Point3<ScalarType> curr_V;
      for (int j = 0; j < eigen_v.cols(); ++j) {
        curr_V.V(j) = eigen_v(i, j);
      }
      VertPos.push_back(curr_V);
    }

    for (int i = 0; i < eigen_f.rows(); ++i) {
      std::vector<int> FaceIdx;
      for (int j = 0; j < eigen_f.cols(); ++j) {
        FaceIdx.push_back(eigen_f(i, j));
      }
      Faces.push_back(FaceIdx);
    }
  }

  static Eigen::Vector3<ScalarType> ToEigen(const Geo::Point3<ScalarType> &V) {
    return (Eigen::Vector3<ScalarType>(V.X, V.Y, V.Z));
  }

  static Geo::Point3<ScalarType> ToPoint3(const Eigen::Vector3<ScalarType> &V) {
    return (Geo::Point3<ScalarType>(V(0), V(1), V(2)));
  }

  static Eigen::Matrix3<ScalarType>
  FramesToEigen(const Geo::Point3<ScalarType> &V0,
                const Geo::Point3<ScalarType> &V1,
                const Geo::Point3<ScalarType> &V2) {
    Eigen::Matrix3<ScalarType> m;

    m(0, 0) = V0.X;
    m(1, 0) = V0.Y;
    m(2, 0) = V0.Z;

    m(0, 1) = V1.X;
    m(1, 1) = V1.Y;
    m(2, 1) = V1.Z;

    m(0, 2) = V2.X;
    m(1, 2) = V2.Y;
    m(2, 2) = V2.Z;

    return m;
  }

  static void EigenToFrames(const Eigen::Matrix3<ScalarType> &m,
                            Geo::Point3<ScalarType> &V0,
                            Geo::Point3<ScalarType> &V1,
                            Geo::Point3<ScalarType> &V2) {
    V0.X = m(0, 0);
    V0.Y = m(1, 0);
    V0.Z = m(2, 0);

    V1.X = m(0, 1);
    V1.Y = m(1, 1);
    V1.Z = m(2, 1);

    V2.X = m(0, 2);
    V2.Y = m(1, 2);
    V2.Z = m(2, 2);
  }

  static void NormalizeColumns(Eigen::Matrix3<ScalarType> &m) {
    Geo::Point3<ScalarType> V0, V1, V2;
    EigenToFrames(m, V0, V1, V2);
    V0.Normalize();
    V1.Normalize();
    V2.Normalize();
    m = FramesToEigen(V0, V1, V2);
  }

  // Point3<S> operator * ( const Point3<S> & v ) const
  //{
  //     Point3<S> t;

  //    t[0] = a[0]*v[0] + a[1]*v[1] + a[2]*v[2];
  //    t[1] = a[3]*v[0] + a[4]*v[1] + a[5]*v[2];
  //    t[2] = a[6]*v[0] + a[7]*v[1] + a[8]*v[2];
  //    return t;
  //}

  static Eigen::Vector3<ScalarType>
  ToEigenVect(const Geo::Point3<ScalarType> &P) {
    Eigen::Vector3<ScalarType> ret;
    ret[0] = P.X;
    ret[1] = P.Y;
    ret[2] = P.Z;
    return ret;
  }

  static std::vector<Eigen::Vector3<ScalarType>>
  ToEigenVect(const std::vector<Geo::Point3<ScalarType>> &P) {
    std::vector<Eigen::Vector3<ScalarType>> ret;
    for (size_t i = 0; i < P.size(); i++)
      ret.push_back(ToEigenVect(P[i]));
    return ret;
  }

  static bool hasNaNInf(const Eigen::SparseMatrix<ScalarType> &matrix) {
    for (int k = 0; k < matrix.outerSize(); ++k) {
      for (typename Eigen::SparseMatrix<ScalarType>::InnerIterator it(matrix,
                                                                      k);
           it; ++it) {
        if (std::isnan(it.value()) || std::isinf(it.value())) {
          return true;
        }
      }
    }
    return false;
  }

  static Eigen::VectorX<ScalarType>
  toEigenVector(const std::vector<ScalarType> &std_vec) {
    Eigen::VectorX<ScalarType> eigen_vec(std_vec.size());
    for (std::size_t i = 0; i < std_vec.size(); ++i) {
      eigen_vec(i) = std_vec[i];
    }
    return eigen_vec;
  }

};

template <class ScalarType>
Geo::Point3<ScalarType> operator*(const Eigen::Matrix4<ScalarType> &m,
                                  const Geo::Point3<ScalarType> &p) {
  ScalarType w;
  Geo::Point3<ScalarType> s;
  s.X = m(0, 0) * p.X + m(0, 1) * p.Y + m(0, 2) * p.Z + m(0, 3);
  s.Y = m(1, 0) * p.X + m(1, 1) * p.Y + m(1, 2) * p.Z + m(1, 3);
  s.Z = m(2, 0) * p.X + m(2, 1) * p.Y + m(2, 2) * p.Z + m(2, 3);
  w = m(3, 0) * p.X + m(3, 1) * p.Y + m(3, 2) * p.Z + m(3, 3);
  if (w != 0)
    s /= w;
  return s;
}

template <class ScalarType>
Geo::Point3<ScalarType> operator*(const Eigen::Matrix3<ScalarType> &m,
                                  const Geo::Point3<ScalarType> &p) {
  Eigen::Vector3<ScalarType> pEig = EigenInterface<ScalarType>::ToEigenVect(p);
  pEig = m * pEig;
  Geo::Point3<ScalarType> rot_p = EigenInterface<ScalarType>::ToPoint3(pEig);
  return rot_p;
}

//// return normals of the mesh
// static void GetNormalData(const MeshType &mesh,
//                           MatrixXm &Nvert,
//                           MatrixXm &Nface)
//{
//     // create eigen matrix of vertices
//     Nvert=MatrixXm(mesh.VN(), 3);
//     Nface=MatrixXm(mesh.FN(), 3);

//    // per vertices normals
//    for (int i = 0; i < mesh.VN(); i++)
//        for (int j = 0; j < 3; j++)
//            Nvert(i,j) = mesh.vert[i].cN()[j];

//    // per vertices normals
//    for (int i = 0; i < mesh.FN(); i++)
//        for (int j = 0; j < 3; j++)
//            Nface(i,j) = mesh.face[i].cN()[j];
//}

//// get face to face adjacency
// static void GetTriFFAdjacency(MeshType &mesh,
//                               Eigen::MatrixXi &FFp,
//                               Eigen::MatrixXi &FFi)
//{
//     tri::UpdateTopology<MeshType>::FaceFace(mesh);
//     FFp = Eigen::MatrixXi(mesh.FN(),3);
//     FFi = Eigen::MatrixXi(mesh.FN(),3);

//    for (int i = 0; i < mesh.FN(); i++)
//        for (int j = 0; j < 3; j++)
//        {
//            FaceType *AdjF=mesh.face[i].FFp(j);
//            if (AdjF==&mesh.face[i])
//            {
//                FFp(i,j)=-1;
//                FFi(i,j)=-1;
//            }
//            else
//            {
//                FFp(i,j)=tri::Index(mesh,AdjF);
//                FFi(i,j)=mesh.face[i].FFi(j);
//            }
//        }
//}

// static void GetUVData(const MeshType &mesh,
//                       MatrixXm & uv)
//{
//     tri::RequireVertexCompactness(mesh);
//     tri::RequirePerVertexTexCoord(mesh);

//    uv = MatrixXm(mesh.VN(), 2);

//    // per vertices uv
//    for (int i = 0; i < mesh.VN(); i++)
//    {
//        uv(i,0) = mesh.vert[i].cT().U();
//        uv(i,1) = mesh.vert[i].cT().V();
//    }
//}

//// get edge to face and edge to vertex adjacency
// static void GetTriEdgeAdjacency(const MeshType &mesh,
//                                 Eigen::MatrixXi& EV,
//                                 Eigen::MatrixXi& FE,
//                                 Eigen::MatrixXi& EF)
//{
//     Eigen::MatrixXi faces;
//     MatrixXm vert;
//     GetTriMeshData(mesh,faces,vert);
//     GetTriEdgeAdjacency(vert,faces,EV,FE,EF);
// }

// static Eigen::Vector3d VectorFromCoord(CoordType v)
//{
//     Eigen::Vector3d ret(v[0],v[1],v[2]);
//     return ret;
// }

// template< class VecType >
// static void PerVertexArea(const MeshType &m, VecType &h)
//{
//     tri::RequireCompactness(m);
//     h.resize(m.vn);
//     fill(h.begin(),h.end(),0);
//     for (int i = 0; i < m.FN(); ++i)
//     {
//         ScalarType a = DoubleArea(m.face[i])/6.0;
//         for(int j=0;j<m.face[i].VN();++j)
//             h[tri::Index(m,m.face[i].cV(j))] += a;
//     }
// }

// template< class VecType >
// static void PerFaceArea(const MeshType &m, VecType &h)
//{
//     tri::RequireCompactness(m);
//     h.resize(m.fn);
//     for(int i=0;i<m.fn;++i)
//         h[i] =DoubleArea(m.face[i])/2.0;
// }

// static void MassMatrixEntry(MeshType &m,
//                             std::vector<std::pair<int,int> > &index,
//                             std::vector<ScalarType> &entry,
//                             bool vertexCoord=true)
//{
//     tri::RequireCompactness(m);

//    typename MeshType::template PerVertexAttributeHandle<ScalarType> h =
//            tri::Allocator<MeshType>:: template
//            GetPerVertexAttribute<ScalarType>(m, "area");
//    for(int i=0;i<m.vn;++i) h[i]=0;

//    for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
//    {
//        ScalarType a = DoubleArea(*fi);
//        for(int j=0;j<fi->VN();++j)
//            h[tri::Index(m,fi->V(j))] += a;
//    }
//    ScalarType maxA=0;
//    for(int i=0;i<m.vn;++i)
//        maxA = std::max(maxA,h[i]);

//    //store the index and the scalar for the sparse matrix
//    for (size_t i=0;i<m.vert.size();i++)
//    {
//        if (vertexCoord)
//        {
//            for (size_t j=0;j<3;j++)
//            {
//                int currI=(i*3)+j;
//                index.push_back(std::pair<int,int>(currI,currI));
//                entry.push_back(h[i]/maxA);
//            }
//        }
//        else
//        {
//            int currI=i;
//            index.push_back(std::pair<int,int>(currI,currI));
//            entry.push_back(h[i]/maxA);
//        }
//    }
//    tri::Allocator<MeshType>::template
//    DeletePerVertexAttribute<ScalarType>(m,h);
//}

// static void GetLaplacianEntry(MeshType &mesh,
//                               FaceType &f,
//                               std::vector<std::pair<int,int> > &index,
//                               std::vector<ScalarType> &entry,
//                               bool cotangent,
//                               ScalarType weight = 1,
//                               bool vertexCoord=true)
//{
//     if (cotangent) vcg::tri::MeshAssert<MeshType>::OnlyTriFace(mesh);

//    for (int i=0;i<f.VN();i++)
//    {

//        if (cotangent)
//        {
//            weight=Harmonic<MeshType>::template
//            CotangentWeight<ScalarType>(f,i);
//        }

//        //get the index of the vertices
//        int indexV0=Index(mesh,f.V0(i));
//        int indexV1=Index(mesh,f.V1(i));

//        if (vertexCoord)
//        {
//            //then assemble the matrix
//            for (int j=0;j<3;j++)
//            {
//                //multiply by 3 and add the component
//                int currI0=(indexV0*3)+j;
//                int currI1=(indexV1*3)+j;

//                index.push_back(std::pair<int,int>(currI0,currI0));
//                entry.push_back(weight);
//                index.push_back(std::pair<int,int>(currI0,currI1));
//                entry.push_back(-weight);

//                index.push_back(std::pair<int,int>(currI1,currI1));
//                entry.push_back(weight);
//                index.push_back(std::pair<int,int>(currI1,currI0));
//                entry.push_back(-weight);
//            }
//        }
//        else
//        {
//            int currI0=(indexV0);
//            int currI1=(indexV1);

//            index.push_back(std::pair<int,int>(currI0,currI0));
//            entry.push_back(weight);
//            index.push_back(std::pair<int,int>(currI0,currI1));
//            entry.push_back(-weight);

//            index.push_back(std::pair<int,int>(currI1,currI1));
//            entry.push_back(weight);
//            index.push_back(std::pair<int,int>(currI1,currI0));
//            entry.push_back(-weight);
//        }
//    }
//}

// static void GetLaplacianMatrix(MeshType &mesh,
//                                std::vector<std::pair<int,int> > &index,
//                                std::vector<ScalarType> &entry,
//                                bool cotangent,
//                                ScalarType weight = 1,
//                                bool vertexCoord=true )
//{
//     //store the index and the scalar for the sparse matrix
//     for (size_t i=0;i<mesh.face.size();i++)
//         GetLaplacianEntry(mesh,mesh.face[i],index,entry,cotangent,weight,vertexCoord);
// }

template <class ScalarType>
inline void SetScale(Eigen::Matrix4<ScalarType> &m, const ScalarType sx,
              const ScalarType sy, const ScalarType sz) {
  m.setZero();
  m(0, 0) = sx;
  m(1, 1) = sy;
  m(2, 2) = sz;
  m(3, 3) = 1;
}

template <class ScalarType>
inline void SetTranslate(Eigen::Matrix4<ScalarType> &m, const ScalarType tx,
                  const ScalarType ty, const ScalarType tz) {
  m.setIdentity();
  m(0, 3) = tx;
  m(1, 3) = ty;
  m(2, 3) = tz;
}

inline bool canMultiply(const Eigen::MatrixXd &A, const Eigen::VectorXd &V) {
  return A.cols() == V.rows();
}

// template<typename ScalarType>
// bool hasNaNInf(const Eigen::SparseMatrix<ScalarType>& matrix) {
//     for (int k = 0; k < matrix.outerSize(); ++k) {
//         for (typename Eigen::SparseMatrix<ScalarType>::InnerIterator
//         it(matrix, k); it; ++it) {
//             if (std::isnan(it.value())||std::isinf(it.value())) {
//                 return true;
//             }
//         }
//     }
//     return false;
// }

#endif