#ifndef SURFACEMAKER_H
#define SURFACEMAKER_H

#include <vector>
#include "CycleDiscovery/amlVec.h"
#include "CycleDiscovery/GraphStructure.h"

// Forward declarations
class cycleUtils;

// Thin façade over CycleDiscovery: feeds a curve network in, runs cycle
// discovery + surface building, and exposes the resulting triangle surface.
class SurfaceMaker {
public:
    SurfaceMaker();
    ~SurfaceMaker();

    void OnOpenCurveNetwork(const std::vector<std::vector<AML::double3>>& curveNetwork);

    const GraphStructure::TriangleSurface& GetTriangleSurface() const;

private:
    cycleUtils* m_cycleUtils;
};


#ifndef CURVE_SURFACING_STATIC
#  include "SurfaceMaker.cpp"
#endif

#endif
