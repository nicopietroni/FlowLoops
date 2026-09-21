#ifndef _TRIANGLE_INFO_H_
#define _TRIANGLE_INFO_H_
#include <limits>

class TriangleInfo {
public:
	TriangleInfo();
	~TriangleInfo();

	int edgeIndex[3];   // index of {1st,2nd,3rd} edge in edgeInfoList
	int triIndex[3];    // index of this triangle in the $side list of edge {1,2,3}
	float optCost[3];   // optimal cost per sub-edge; UNCOMPUTED_COST until filled
	int   optTile[3];   // optimal tiling index; -1 until filled

	// Sentinel used to mark optCost slots that have not yet been computed.
	// Deliberately negative so any legitimate cost (>= 0) is distinguishable.
	static constexpr float UNCOMPUTED_COST = -1.0f;

	int getSize();
	int getOptSize();
};


#ifndef CURVE_SURFACING_STATIC
#  include "TriangleInfo.cpp"
#endif

#endif
