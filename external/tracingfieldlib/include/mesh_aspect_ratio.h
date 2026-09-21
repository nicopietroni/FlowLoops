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

#ifndef MESH_ASPECT_RATIO
#define MESH_ASPECT_RATIO

#include <cassert>
#include <mesh_functions.h>
#include <triangle3.h>

namespace Geo{

enum QType{QTRadii,QTMeanRatio};

template<class ScalarType>
ScalarType TriMeshQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                          const std::vector<std::vector<int> > &Faces,
                          const int &IndexF,QType QT=QTRadii)
{
   assert(IndexF>=0);
   assert(IndexF<Faces.size());
   assert(Faces[IndexF].size()==3);
   Geo::Traingle3<ScalarType> Tri=Geo::TriangleOfFace(VertPos,Faces,IndexF);
   if (QT==QTRadii)
       return (Geo::QualityRadii(Tri));
   else
       return (Geo::QualityMeanRatio(Tri));
}

template<class ScalarType>
ScalarType TriMeshAvgQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                             const std::vector<std::vector<int> > &Faces,
                             QType QT=QTRadii)
{
    ScalarType AvgQ=0;
    for (size_t i=0;i<Faces.size();i++)
        AvgQ+=TriMeshQuality(VertPos,Faces,i,QT);
    AvgQ/=Faces.size();
    return (AvgQ);
}

template<class ScalarType>
ScalarType TriMeshAvgQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                             const std::vector<std::vector<int> > &Faces,
                             const std::vector<int> &FacesIdx,QType QT=QTRadii)
{
    ScalarType AvgQ=0;
    for (size_t i=0;i<FacesIdx.size();i++)
        AvgQ+=TriMeshQuality(VertPos,Faces,FacesIdx[i],QT);
    AvgQ/=FacesIdx.size();
    return (AvgQ);
}

template<class ScalarType>
ScalarType TriMeshWorstQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                              const std::vector<std::vector<int> > &Faces,
                              const std::vector<int> &FacesIdx,QType QT=QTRadii)
{
    assert(FacesIdx.size()>0);
    ScalarType WorstW=1;
    for (size_t i=0;i<FacesIdx.size();i++)
        WorstW=std::min(TriMeshQuality(VertPos,Faces,FacesIdx[i],QT),WorstW);
    return (WorstW);
}

template<class ScalarType>
ScalarType TriMeshWorstQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                              const std::vector<std::vector<int> > &Faces,
                              QType QT=QTRadii)
{
    ScalarType WorstW=1;
    for (size_t i=0;i<Faces.size();i++)
        WorstW=std::min(TriMeshQuality(VertPos,Faces,i,QT),WorstW);
    return (WorstW);
}
//   assert(IsAllTris(VertPos,Faces));
//   for (size)

}


#endif
