#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "SurfaceMaker.h"
#include "CycleDiscovery/cycleUtils.h"
#include "cs_logging.h"


CS_INLINE SurfaceMaker::SurfaceMaker() {
    m_cycleUtils = new cycleUtils();
}

CS_INLINE SurfaceMaker::~SurfaceMaker() {
    delete m_cycleUtils;
}

CS_INLINE void SurfaceMaker::OnOpenCurveNetwork(const std::vector<std::vector<AML::double3>>& curveNetwork)
{
	if (curveNetwork.empty()) {
		Geex::cs_err() << "No curve data is loaded!\n";
		m_cycleUtils->resetAll();
		return;
	}

	m_cycleUtils->resetAll();
	m_cycleUtils->m_cycleDiscovery.m_curves = curveNetwork;

	// Set default capacities (assuming a default capacity of 2 for each curve)
	m_cycleUtils->m_cycleDiscovery.m_curveCapacitys.assign(curveNetwork.size(), 2);

	m_cycleUtils->m_cycleDiscovery.constructNetwork();

	// Processing logic
	m_cycleUtils->m_cycleDiscovery.constructRotationGraphbyPoleGraphEx();
	m_cycleUtils->m_cycleDiscovery.constructCycles();
	m_cycleUtils->m_cycleDiscovery.cycleBreaking();
	m_cycleUtils->surfaceBuilding();
}

CS_INLINE const GraphStructure::TriangleSurface& SurfaceMaker::GetTriangleSurface() const
{
    return m_cycleUtils->m_triangleSurface;
}
