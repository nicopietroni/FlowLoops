#ifndef CYCLE_UTILS
#define CYCLE_UTILS

#include <vector>
#include <string>
#include "amlVec.h"
#include "GraphStructure.h"
#include <time.h>
#include "CycleDiscovery.h"

using namespace GraphStructure;

class cycleUtils {
public:
    cycleUtils();
    ~cycleUtils();

    // ── Core data ─────────────────────────────────────────────────────────────
    std::string m_filename;

    cycleDiscovery m_cycleDiscovery;

    // Surface-building weights and smoothing parameters
    double m_weightTri, m_weightEdge, m_weightBiTri, m_weightTriBd, m_weightWorsDih;
    bool   m_SurfaceSmooth;
    int    m_subdivisonSmooth, m_laplacianSmooth;
    bool   m_isUpdateSubdivision, m_isUpdateLaplacian;

    TriangleSurface m_triangleSurface;
    TriangleSurface m_triangleSurfaceNormal;

    // ── State reset ───────────────────────────────────────────────────────────
    void resetAllData();
    void resetAll();

    // ── Capacity / constraint helpers ─────────────────────────────────────────
    bool updateConstraintList();
    int  getSelectArcCapacity();
    bool updateCapacity() {
        if (!m_cycleDiscovery.m_errorNode.empty()) return false;
        for (int i = 0; i < (int)m_cycleDiscovery.m_latestUpdateArcs.size(); i++) {
            m_cycleDiscovery.m_curveNet.arcs[m_cycleDiscovery.m_latestUpdateArcs[i]].capacity =
                m_cycleDiscovery.m_latestUpdateArcsCapacity[i];
            m_cycleDiscovery.m_curveCapacitys[m_cycleDiscovery.m_latestUpdateArcs[i]] =
                m_cycleDiscovery.m_latestUpdateArcsCapacity[i];
        }
        return true;
    }
    void setCapacity(int selection);

    void resetUpdateArcs() {
        m_cycleDiscovery.m_latestUpdateArcs.clear();
        m_cycleDiscovery.m_latestUpdateArcsCapacity.clear();
        m_cycleDiscovery.m_selectArcList.back().clear();
        m_cycleDiscovery.m_errorNode.clear();
    }
    void resetUpdateNodes() {
        m_cycleDiscovery.m_latestUpdateNodes.clear();
        for (int i = 0; i < (int)m_cycleDiscovery.m_latestUpdateArcs.size(); i++) {
            m_cycleDiscovery.m_latestUpdateNodes.push_back(
                m_cycleDiscovery.m_curveNet.arcs[m_cycleDiscovery.m_latestUpdateArcs[i]].endNodesID.first);
            m_cycleDiscovery.m_latestUpdateNodes.push_back(
                m_cycleDiscovery.m_curveNet.arcs[m_cycleDiscovery.m_latestUpdateArcs[i]].endNodesID.second);
        }
        m_cycleDiscovery.m_latestUpdateArcs.clear();
        for (int i = 0; i < (int)m_cycleDiscovery.m_latestUpdateNodes.size(); i++) {
            std::vector<int> tarcs =
                m_cycleDiscovery.m_curveNet.nodes[m_cycleDiscovery.m_latestUpdateNodes[i]].arcID;
            for (int j = 0; j < (int)tarcs.size(); j++)
                m_cycleDiscovery.m_latestUpdateArcs.push_back(tarcs[j]);
        }
    }

    // ── Surface building ──────────────────────────────────────────────────────
    void surfaceBuilding();
    void subDivideCycles();
    void writeConfigurations();
};


#ifndef CURVE_SURFACING_STATIC
#  include "cycleUtils.cpp"
#endif

#endif