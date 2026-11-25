#pragma once
#include "../curve_surfacing_inline.h"

// ── implementation includes ──
#include "EdgeInfo.h"


CS_INLINE EdgeInfo::EdgeInfo()
	: v1(-1), v2(-1)
	, leftsize(-1), rightsize(-1)
	// Initialise all owning pointers to nullptr so the destructor is safe
	// even if buildList() throws before allocating them.
	, leftTris(nullptr), leftEdgeInd(nullptr)
	, rightTris(nullptr), rightEdgeInd(nullptr)
{
}

CS_INLINE EdgeInfo::~EdgeInfo() {
	delete[] leftTris;
	delete[] leftEdgeInd;
	delete[] rightTris;
	delete[] rightEdgeInd;
}

CS_INLINE int EdgeInfo::getSize() {
	int totalSize = 0;
	if (leftsize  > 0) totalSize += leftsize  * (sizeof(int) + sizeof(char));
	if (rightsize > 0) totalSize += rightsize * (sizeof(int) + sizeof(char));
	return totalSize;
}