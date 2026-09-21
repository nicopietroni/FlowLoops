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

#ifndef LINEAR_CONSTRAINT
#define LINEAR_CONSTRAINT

#include <cstddef>
#include <vector>
#include <point3.h>

namespace Geo
{

template <class ScalarType>
struct VarTerm
{
    size_t IndexVar;
    size_t IndexCoord;
    ScalarType W;

    VarTerm(size_t _IndexVar,
            size_t _IndexCoord,
            ScalarType _W)
    {
        IndexVar=_IndexVar;
        IndexCoord=_IndexCoord;
        W=_W;
    }

    VarTerm()
    {
        IndexVar=0;
        IndexCoord=0;
        W=0;
    }
};

template <class ScalarType>
class LinearConstraint
{

public:

    ScalarType RHSide;
    std::vector<VarTerm<ScalarType> > LHSide;
    ScalarType Weight;

    LinearConstraint()
    {
        RHSide=0;
        Weight=1;
    }

    bool hasNaNInf() 
    {
        if (isnan(RHSide)||
            isinf(RHSide))
        {std::cout<<"Error on LinConstr RHS"<<std::endl;return true;}

        if (isnan(Weight)||
            isinf(Weight))
        {std::cout<<"Error on LinConstr Weight"<<std::endl;return true;}

        for (size_t i=0;i<LHSide.size();i++)
        {
        if (isnan(LHSide[i].W)||
            isinf(LHSide[i].W))
        {std::cout<<"Error on LinConstr Entries"<<std::endl;return true;}
        }
        return false;
    }
};

//template <class ScalarType>
//void FixedLinearConstraints(const std::vector<Geo::Point3<ScalarType> > &Pos,
//                            const std::vector<int> &FixedV,
//                            const ScalarType Weight,
//                            std::vector<LinearConstraint<ScalarType> > &FixConstr)
//{
//    FixConstr.clear();
//    for (size_t i=0;i<Pos.size();i++)
//        for (size_t j=0;j<3;j++)
//        {
//            Geo::LinearConstraint<ScalarType> Constr;
//            Constr.RHSide=Pos[i].cV(j);
//            Constr.Weight=Weight;
//            Geo::VarTerm<ScalarType> Var(i,j,1);
//            Constr.LHSide.push_back(Var);
//            FixConstr.push_back(Constr);
//        }
//}

}

#endif
