#ifndef CURVE_SURFACING_LINEAR_CONSTRAINT
#define CURVE_SURFACING_LINEAR_CONSTRAINT

#include <cstddef>
#include <vector>
#include "cs_logging.h"

namespace Geex
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
        {cs_err()<<"Error on LinConstr RHS"<<std::endl;return true;}

        if (isnan(Weight)||
            isinf(Weight))
        {cs_err()<<"Error on LinConstr Weight"<<std::endl;return true;}

        for (size_t i=0;i<LHSide.size();i++)
        {
        if (isnan(LHSide[i].W)||
            isinf(LHSide[i].W))
        {cs_err()<<"Error on LinConstr Entries"<<std::endl;return true;}
        }
        return false;
    }
};

}

#endif