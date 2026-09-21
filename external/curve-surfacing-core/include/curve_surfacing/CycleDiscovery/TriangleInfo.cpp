#pragma once
#include "../curve_surfacing_inline.h"

// ── implementation includes ──
#include "TriangleInfo.h"


TriangleInfo::TriangleInfo() {
	// Mark all sub-problem costs as uncomputed.
	// Previously used FLT_MIN (a tiny positive float) as sentinel, which was
	// ambiguous when compared against real costs. UNCOMPUTED_COST is negative
	// so it can never be confused with a valid cost.
	optCost[0] = optCost[1] = optCost[2] = UNCOMPUTED_COST;
	optTile[0] = optTile[1] = optTile[2] = -1;
}

CS_INLINE TriangleInfo::~TriangleInfo() {
}

int TriangleInfo::getOptSize() {
	return sizeof(optCost) + sizeof(optTile);
}

int TriangleInfo::getSize() {
	return sizeof(optCost) + sizeof(optTile)
		 + sizeof(edgeIndex) + sizeof(triIndex);
}