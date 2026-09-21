#pragma once
#include "../curve_surfacing_inline.h"

// ── implementation includes ──
#include "cycleUtils.h"

#include "DrT.h"
#include <iostream>
#include <fstream>
#include <cfloat>
#include <algorithm>

CS_INLINE cycleUtils::cycleUtils()
{
    resetAll();
}

CS_INLINE cycleUtils::~cycleUtils()
{
}

CS_INLINE void cycleUtils::resetAllData()
{
    m_cycleDiscovery.m_timePre    = 0;
    m_cycleDiscovery.m_timeRun    = 0;
    m_cycleDiscovery.m_timeNewRun = 0;

    m_weightTri = m_weightEdge = 0;
    m_weightBiTri = m_weightTriBd = 1;
    m_weightWorsDih = 0;
    m_SurfaceSmooth     = true;
    m_subdivisonSmooth  = 1;
    m_laplacianSmooth   = 3;

    m_cycleDiscovery.m_isUpdatetwistNormNum      = false;
    m_cycleDiscovery.m_isUpdateangleWeight       = false;
    m_cycleDiscovery.m_isUpdatetwistWeight       = false;
    m_cycleDiscovery.m_isUpdatecurveWeight       = false;
    m_cycleDiscovery.m_isUpdatenodeWeight        = false;
    m_cycleDiscovery.m_isUpdatebestNeightboreNum = false;
    m_cycleDiscovery.m_isUpdaterotationGraphNum  = false;
    m_cycleDiscovery.m_isUpdatedihedralWeight    = false;
    m_cycleDiscovery.m_isUpdatestateNum          = false;
    m_cycleDiscovery.m_isUpdateRoGraphConnect    = false;
    m_cycleDiscovery.m_isUpdateDoCycleBreak      = false;
    m_cycleDiscovery.m_isUpdateWeightTri =
    m_cycleDiscovery.m_isUpdateWeightEdge =
    m_cycleDiscovery.m_isUpdateWeightBiTri =
    m_cycleDiscovery.m_isUpdateWeightTriBd =
    m_cycleDiscovery.isUpdateWeightWorsDih =
    m_cycleDiscovery.m_isUpdateSurfaceSmooth = true;
    m_isUpdateSubdivision = m_isUpdateLaplacian = false;
    m_cycleDiscovery.m_isUpdateCycle = false;

    m_cycleDiscovery.m_selectArcList.clear();
    m_cycleDiscovery.m_selectArcList.resize(1);
    m_cycleDiscovery.m_userDefinedPairsInNode.clear();
    m_cycleDiscovery.m_userDefinedPairsInArc.clear();
    m_cycleDiscovery.m_latestUpdateNodes.clear();
    m_cycleDiscovery.m_latestUpdateArcs.clear();
    m_cycleDiscovery.m_latestUpdateArcsCapacity.clear();
    m_cycleDiscovery.m_errorNode.clear();
    m_cycleDiscovery.m_selectedNodeInPole.clear();

    m_triangleSurface.clear();
    m_triangleSurfaceNormal.clear();
}

CS_INLINE void cycleUtils::resetAll()
{
    resetAllData();
}

// ── Capacity helpers ──────────────────────────────────────────────────────────

CS_INLINE int cycleUtils::getSelectArcCapacity()
{
    int selectedArcID = m_cycleDiscovery.m_selectArcList.back().empty()
                        ? 0
                        : m_cycleDiscovery.m_selectArcList.back().back();

    auto it = std::find(m_cycleDiscovery.m_latestUpdateArcs.begin(),
                        m_cycleDiscovery.m_latestUpdateArcs.end(), selectedArcID);
    if (it != m_cycleDiscovery.m_latestUpdateArcs.end()) {
        int ind = it - m_cycleDiscovery.m_latestUpdateArcs.begin();
        return m_cycleDiscovery.m_latestUpdateArcsCapacity[ind];
    }
    return m_cycleDiscovery.m_curveNet.arcs[selectedArcID].capacity;
}

CS_INLINE void cycleUtils::setCapacity(int selection)
{
    int m_capacity = selection;
    m_cycleDiscovery.m_errorNode.clear();

    if (m_cycleDiscovery.m_selectArcList.back().empty()) return;
    int m_selectedArcID = m_cycleDiscovery.m_selectArcList.back().back();

    auto& latestArcs = m_cycleDiscovery.m_latestUpdateArcs;
    auto& latestCaps = m_cycleDiscovery.m_latestUpdateArcsCapacity;

    if (m_cycleDiscovery.m_curveNet.arcs[m_selectedArcID].capacity == m_capacity) {
        auto it = std::find(latestArcs.begin(), latestArcs.end(), m_selectedArcID);
        if (it != latestArcs.end()) {
            int ind = it - latestArcs.begin();
            latestArcs.erase(latestArcs.begin() + ind);
            latestCaps.erase(latestCaps.begin() + ind);
        }
    } else {
        auto it = std::find(latestArcs.begin(), latestArcs.end(), m_selectedArcID);
        if (it != latestArcs.end()) {
            latestCaps[it - latestArcs.begin()] = m_capacity;
        } else {
            latestArcs.push_back(m_selectedArcID);
            latestCaps.push_back(m_capacity);
        }
    }

    // Validate node parity
    std::vector<int> useNode;
    for (int i = 0; i < (int)latestArcs.size(); i++) {
        int arcID = latestArcs[i];
        int nodeID[] = { m_cycleDiscovery.m_curveNet.arcs[arcID].endNodesID.first,
                         m_cycleDiscovery.m_curveNet.arcs[arcID].endNodesID.second };
        for (int ni : nodeID) {
            if (std::find(useNode.begin(), useNode.end(), ni) == useNode.end())
                useNode.push_back(ni);
        }
    }
    for (int i = 0; i < (int)useNode.size(); i++) {
        std::vector<int> adjArc = m_cycleDiscovery.m_curveNet.nodes[useNode[i]].arcID;
        int cap = 0;
        for (int j = 0; j < (int)adjArc.size(); j++) {
            auto it = std::find(latestArcs.begin(), latestArcs.end(), adjArc[j]);
            if (it != latestArcs.end())
                cap += latestCaps[it - latestArcs.begin()];
            else
                cap += m_cycleDiscovery.m_curveNet.arcs[adjArc[j]].capacity;
        }
        if ((cap % 2) == 1)
            m_cycleDiscovery.m_errorNode.push_back(useNode[i]);
    }
}

CS_INLINE bool cycleUtils::updateConstraintList()
{
    auto& m_curveNet = m_cycleDiscovery.m_curveNet;

    if (m_curveNet.arcs.empty()) return false;
    if (m_cycleDiscovery.m_userDefinedPairsInNode.empty())
        m_cycleDiscovery.m_userDefinedPairsInNode.resize(m_curveNet.nodes.size());
    if (m_cycleDiscovery.m_userDefinedPairsInArc.empty())
        m_cycleDiscovery.m_userDefinedPairsInArc.resize(m_curveNet.arcs.size());

    // De-duplicate selection list
    auto& selBack = m_cycleDiscovery.m_selectArcList.back();
    for (int i = 0; i < (int)selBack.size() - 1; i++)
        for (int j = i + 1; j < (int)selBack.size(); j++)
            if (selBack[i] == selBack[j]) { selBack.erase(selBack.begin() + j); j--; }

    std::vector<std::vector<int>> jointArcs(m_curveNet.nodes.size());
    std::vector<std::vector<std::vector<int>>> adjArcs(m_curveNet.arcs.size(),
                                                        std::vector<std::vector<int>>(2));
    std::vector<bool> usedArcs(adjArcs.size(), false);

    for (int i = 0; i < (int)selBack.size(); i++) {
        int a1 = selBack[i];
        int n[] = { m_curveNet.arcs[a1].endNodesID.first, m_curveNet.arcs[a1].endNodesID.second };
        jointArcs[n[0]].push_back(a1);
        jointArcs[n[1]].push_back(a1);
        usedArcs[a1] = true;
        for (int side = 0; side < 2; side++) {
            auto& leftAdj    = m_curveNet.nodes[n[side]].arcID;
            auto& leftAdjDir = m_curveNet.nodes[n[side]].arcDirection;
            for (int j = 0; j < (int)leftAdj.size(); j++) {
                if (leftAdj[j] != a1)
                    adjArcs[leftAdj[j]][leftAdjDir[j] - 1].push_back(a1);
            }
        }
    }

    bool isRightCycle = true;
    for (int i = 0; i < (int)jointArcs.size(); i++)
        if (jointArcs[i].size() > 2) { isRightCycle = false; break; }
    if (!isRightCycle) return isRightCycle;

    isRightCycle = false;
    for (int i = 0; i < (int)jointArcs.size(); i++) {
        if (jointArcs[i].size() < 2) continue;
        isRightCycle = true;
        int a1 = jointArcs[i][0], a2 = jointArcs[i][1];
        int order1 = (m_curveNet.arcs[a1].endNodesID.first == i)
                     ? m_curveNet.arcs[a1].posInNode.first
                     : m_curveNet.arcs[a1].posInNode.second;
        int order2 = (m_curveNet.arcs[a2].endNodesID.first == i)
                     ? m_curveNet.arcs[a2].posInNode.first
                     : m_curveNet.arcs[a2].posInNode.second;
        if (order1 > order2) std::swap(order1, order2);
        m_cycleDiscovery.m_userDefinedPairsInNode[i].push_back({order1, order2});
        m_cycleDiscovery.m_latestUpdateNodes.push_back(i);
    }
    for (int i = 0; i < (int)adjArcs.size(); i++) {
        if (adjArcs[i][0].size() != 1 || adjArcs[i][1].size() != 1) continue;
        if (!usedArcs[i]) continue;
        int a1 = adjArcs[i][0][0], a2 = adjArcs[i][1][0];
        int n[] = { m_curveNet.arcs[i].endNodesID.first, m_curveNet.arcs[i].endNodesID.second };
        int order1 = (m_curveNet.arcs[a1].endNodesID.first == n[0])
                     ? m_curveNet.arcs[a1].posInNode.first
                     : m_curveNet.arcs[a1].posInNode.second;
        int order2 = (m_curveNet.arcs[a2].endNodesID.first == n[1])
                     ? m_curveNet.arcs[a2].posInNode.first
                     : m_curveNet.arcs[a2].posInNode.second;
        std::pair<int,int> temp = {order1, order2};
        bool isExist = false;
        for (auto& p : m_cycleDiscovery.m_userDefinedPairsInArc[i])
            if (p.first == temp.first && p.second == temp.second) { isExist = true; break; }
        if (!isExist) {
            m_cycleDiscovery.m_userDefinedPairsInArc[i].push_back(temp);
            m_cycleDiscovery.m_latestUpdateArcs.push_back(i);
        }
    }
    return isRightCycle;
}

// ── writeConfigurations, surfaceBuilding, subDivideCycles ─────────────────────
// (implementations unchanged from original — copied verbatim below)

CS_INLINE void cycleUtils::writeConfigurations()
{
    std::fstream File;
    char FileName[400] = {};
    for (int i = (int)strlen(FileName) - 1; i >= 0; i--) {
        if (FileName[i] == '.') { FileName[i] = '\0'; break; }
    }
    char FileNameExt[400] = {};
    int ind = 0;
    for (int i = (int)strlen(FileName) - 1; i >= 0; i--) {
        if (FileName[i] == '/' || FileName[i] == '\\') { FileName[i] = '\0'; break; }
        FileNameExt[ind] = FileName[i]; ind++;
    }
    FileNameExt[ind] = '\0';
    char ModelName[400] = {};
    for (int i = 0; i < ind; i++) ModelName[i] = FileNameExt[ind - i - 1];
    ModelName[ind] = '\0';

    strcat(FileName, "\\result.txt");
    File.open(FileName, std::ios::app);

    int len = 12;
    std::vector<int> textSize(len, 10);
    textSize[0]  = 30;
    textSize[10] = 11;

    File.seekg(0, File.end);
    int fileLength = File.tellp();
    if (fileLength == 0) {
        File << std::endl;
        for (int i = 0; i < len; i++) {
            for (int j = 0; j < textSize[i]; j++) File << '-';
            if (i != (len - 1)) { File.width(3); File << "-|-"; }
            else                { File.width(2); File << "-|"; }
        }
        File << std::endl;
    }

    if ((int)strlen(ModelName) > textSize[0]) ModelName[textSize[0]] = '\0';
    File.width(textSize[0]);  File << ModelName;     File.width(3); File << " | ";
    int curveSize = m_cycleDiscovery.m_curveNet.arcs.size();
    File.width(textSize[1]);  File << curveSize;     File.width(3); File << " | ";
    int nonManifold = 0;
    for (int i = 0; i < (int)m_cycleDiscovery.m_curveCapacitys.size(); i++)
        if (m_cycleDiscovery.m_curveCapacitys[i] != 2) nonManifold++;
    File.width(textSize[2]);  File << nonManifold;   File.width(3); File << " | ";
    int nodeSize = m_cycleDiscovery.m_curveNet.nodes.size();
    File.width(textSize[3]);  File << nodeSize;      File.width(3); File << " | ";
    int maxDegree = 2;
    for (int i = 0; i < (int)m_cycleDiscovery.m_curveNet.nodes.size(); i++) {
        int d = m_cycleDiscovery.m_curveNet.nodes[i].arcID.size();
        if (maxDegree < d) maxDegree = d;
    }
    File.width(textSize[4]);  File << maxDegree;     File.width(3); File << " | ";
    File.width(textSize[5]-2); File << m_cycleDiscovery.m_timePre + m_cycleDiscovery.m_timeRun;
    File.width(2); File << "mm"; File.width(3); File << " | ";
    File.width(textSize[6]-2); File << m_cycleDiscovery.m_timePre;
    File.width(2); File << "mm"; File.width(3); File << " | ";
    File.width(textSize[7]-2); File << m_cycleDiscovery.m_timeRun;
    File.width(2); File << "mm"; File.width(3); File << " | ";
    int totalConstraint = 0;
    for (int i = 0; i < (int)m_cycleDiscovery.m_selectArcList.size() - 1; i++)
        totalConstraint += m_cycleDiscovery.m_selectArcList[i].size();
    File.width(textSize[8]); File << totalConstraint; File.width(3); File << " | ";
    int constraintCycle = (int)m_cycleDiscovery.m_selectArcList.size() - 1;
    File.width(textSize[9]); File << constraintCycle; File.width(3); File << " | ";
    int genus = nodeSize - curveSize + (int)m_cycleDiscovery.m_cycleSetBreaked.size() - 2;
    genus /= -2;
    File.width(textSize[10]);
    nonManifold == 0 ? genus >= 0 ? File << genus : File << "Disconnect" : File << "NonManifold";
    File.width(3); File << " | ";
    File.width(textSize[11]);
    m_cycleDiscovery.m_isCycleBreak ? File << m_cycleDiscovery.m_breakNum : File << "No";
    File.width(3); File << " | ";
    File << std::endl;
    File.close();
}

CS_INLINE void cycleUtils::surfaceBuilding()
{
    auto& m_curveNet        = m_cycleDiscovery.m_curveNet;
    auto& m_cycleSetBreaked = m_cycleDiscovery.m_cycleSetBreaked;
    auto& m_cycleNormal     = m_cycleDiscovery.m_cycleNormal;
    auto& m_normalsTable    = m_cycleDiscovery.m_normalsTable;

    if (m_curveNet.arcs.empty() || m_cycleSetBreaked.empty()) return;

    if (!(m_cycleDiscovery.m_isUpdateWeightTri  || m_cycleDiscovery.m_isUpdateWeightEdge  ||
          m_cycleDiscovery.m_isUpdateWeightBiTri || m_cycleDiscovery.m_isUpdateWeightTriBd ||
          m_cycleDiscovery.isUpdateWeightWorsDih  || m_cycleDiscovery.m_isUpdateSurfaceSmooth ||
          m_isUpdateSubdivision || m_isUpdateLaplacian) && !m_cycleDiscovery.m_isUpdateCycle)
        return;

    m_cycleNormal.clear();

    Graph      net    = m_curveNet;
    CycleSet   cycles = m_cycleSetBreaked;
    TriangleSurface surface, norms;

    float weights[] = { float(m_weightTri), float(m_weightEdge),
                        float(m_weightBiTri), float(m_weightTriBd),
                        float(m_weightWorsDih) };
    bool dosmooth = m_SurfaceSmooth;
    int  subs     = m_subdivisonSmooth;
    int  laps     = m_laplacianSmooth;

    for (int c = 0; c < (int)cycles.size(); c++) {
        Cycle cycle = cycles[c];
        LinearCurveNet curves;
        for (int i = 0; i < (int)cycle.size(); i++) {
            int arcID = cycle[i].arcID;
            curves.push_back(net.arcs[arcID].vertexList);
            int strNode = m_curveNet.arcs[arcID].endNodesID.first;
            int endNode = m_curveNet.arcs[arcID].endNodesID.second;
            curves.back().insert(curves.back().begin(), m_curveNet.nodes[strNode].pos);
            curves.back().push_back(m_curveNet.nodes[endNode].pos);
            if (cycle[i].strEndID == 2)
                reverse(curves.back().begin(), curves.back().end());
        }

        Curve pointList;
        for (int i = 0; i < (int)curves.size(); i++)
            for (int j = 0; j < (int)curves[i].size() - 1; j++)
                pointList.push_back(curves[i][j]);

        int     point_num = pointList.size();
        double* points    = new double[point_num * 3];
        for (int i = 0; i < point_num; i++) {
            points[i*3+0] = pointList[i].x;
            points[i*3+1] = pointList[i].y;
            points[i*3+2] = pointList[i].z;
        }

        std::vector<double> tile_list, newPoints;
        std::vector<float>  newNormals;
        int res;

        if (m_normalsTable.empty()) {
            res = delaunayRestrictedTriangulation(points, point_num,
                                                  newPoints, tile_list,
                                                  weights, dosmooth, subs, laps);
        } else {
            // Build per-curve normals via parallel transport (unchanged from original)
            CycleSegment cyl = cycle.back();
            cycle.insert(cycle.begin(), cyl);
            cycle.push_back(cycle[1]);
            LinearCurveNet cycleNormal, endNorm;
            std::vector<double> variance;
            std::vector<int>    dirts;
            for (int j = 1; j < (int)cycle.size() - 1; j++) {
                int arcID  = cycle[j].arcID;
                int nodeID[] = { m_curveNet.arcs[arcID].endNodesID.first,
                                 m_curveNet.arcs[arcID].endNodesID.second };
                int a1 = cycle[j-1].arcID, a2 = cycle[j+1].arcID;
                auto getInd = [&](int arc, int node) {
                    return (m_curveNet.arcs[arc].endNodesID.first == node)
                           ? m_curveNet.arcs[arc].posInNode.first
                           : m_curveNet.arcs[arc].posInNode.second;
                };
                int aind1 = 0, aind2 = 0;
                if      (m_curveNet.arcs[a1].endNodesID.first  == nodeID[0]) aind1 = getInd(a1, nodeID[0]);
                else if (m_curveNet.arcs[a1].endNodesID.second == nodeID[0]) aind1 = getInd(a1, nodeID[0]);
                else if (m_curveNet.arcs[a2].endNodesID.first  == nodeID[0]) aind1 = getInd(a2, nodeID[0]);
                else if (m_curveNet.arcs[a2].endNodesID.second == nodeID[0]) aind1 = getInd(a2, nodeID[0]);
                if      (m_curveNet.arcs[a1].endNodesID.first  == nodeID[1]) aind2 = getInd(a1, nodeID[1]);
                else if (m_curveNet.arcs[a1].endNodesID.second == nodeID[1]) aind2 = getInd(a1, nodeID[1]);
                else if (m_curveNet.arcs[a2].endNodesID.first  == nodeID[1]) aind2 = getInd(a2, nodeID[1]);
                else if (m_curveNet.arcs[a2].endNodesID.second == nodeID[1]) aind2 = getInd(a2, nodeID[1]);

                cycleNormal.push_back(m_normalsTable[arcID][aind1][aind2]);
                if (cycle[j].strEndID == 1) {
                    dirts.push_back(1);
                } else {
                    dirts.push_back(2);
                    reverse(cycleNormal.back().begin(), cycleNormal.back().end());
                    for (auto& v : cycleNormal.back()) v = -v;
                }
                cycleNormal.back().pop_back();
                cycleNormal.back().erase(cycleNormal.back().begin());
                std::vector<AML::double3> temp = { cycleNormal.back().front(), cycleNormal.back().back() };
                endNorm.push_back(temp);
                variance.push_back(m_cycleDiscovery.m_twistTablesConfidence[arcID][aind1][aind2]);
            }
            // Wrap end-norm bookends
            endNorm.insert(endNorm.begin(), endNorm.back());
            endNorm.push_back(endNorm[1]);
            variance.insert(variance.begin(), variance.back());
            variance.push_back(variance[1]);

            // Transport and blend normals (unchanged logic from original)
            for (int j = 0; j < (int)cycleNormal.size(); j++) {
                AML::double3 tempL = endNorm[j][1];
                AML::double3 tempR = endNorm[j+2][0];
                Curve pts = curves[j];
                std::vector<std::vector<double>> transportMatrices;
                cycleDiscovery::computeTransportMatrixAll(pts, transportMatrices);
                auto buildTransMatrices = [&](const std::vector<std::vector<double>>& tms) {
                    std::vector<std::vector<AML::double3>> result;
                    for (auto& tm : tms) {
                        std::vector<AML::double3> mat;
                        for (int r = 0; r < 3; r++) {
                            AML::double3 row; for (int k = 0; k < 3; k++) row[k] = tm[r + k*3];
                            mat.push_back(row);
                        }
                        result.push_back(mat);
                    }
                    return result;
                };
                auto transMatrices = buildTransMatrices(transportMatrices);
                std::vector<AML::double3> nmlistL;
                for (auto& mat : transMatrices) {
                    AML::double3 tnm;
                    for (int l = 0; l < 3; l++) tnm[l] = tempL.dot(mat[l]);
                    tnm.normalize(); nmlistL.push_back(tnm);
                }
                reverse(pts.begin(), pts.end());
                transportMatrices.clear();
                cycleDiscovery::computeTransportMatrixAll(pts, transportMatrices);
                transMatrices = buildTransMatrices(transportMatrices);
                std::vector<AML::double3> nmlistR;
                for (auto& mat : transMatrices) {
                    AML::double3 tnm;
                    for (int l = 0; l < 3; l++) tnm[l] = tempR.dot(mat[l]);
                    tnm.normalize(); nmlistR.push_back(tnm);
                }
                reverse(nmlistR.begin(), nmlistR.end());
                for (int k = 0; k < (int)cycleNormal[j].size(); k++) {
                    AML::double3 tnm = variance[j]*nmlistL[k] + variance[j+2]*nmlistR[k] + cycleNormal[j][k];
                    tnm.normalize();
                    cycleNormal[j][k] = tnm;
                }
            }
            Curve normalList;
            for (auto& cn : cycleNormal) for (auto& v : cn) normalList.push_back(v);
            int    normal_num = normalList.size();
            float* normals    = new float[normal_num * 3];
            for (int i = 0; i < normal_num; i++) {
                normals[i*3+0] = normalList[i].x;
                normals[i*3+1] = normalList[i].y;
                normals[i*3+2] = normalList[i].z;
            }
            res = delaunayRestrictedTriangulation(points, normals, point_num,
                                                  newPoints, newNormals,
                                                  tile_list, weights, dosmooth, subs, laps);
            // Store for normals visualisation data (kept for surfaceBuilding contract)
            LinearCurveNet cycleNormalForVis = cycleNormal;
            for (int j = 0; j < (int)cycleNormalForVis.size(); j++) {
                cycleNormalForVis[j].push_back(AML::double3(0,0,0));
                if (dirts[j] == 2) reverse(cycleNormalForVis[j].begin(), cycleNormalForVis[j].end());
            }
            m_cycleNormal.push_back(cycleNormalForVis);
            delete[] normals;
        }

        // Collect triangles
        TriangleCycle triangleCycle, triangleCycleNorm;
        if (res == 1) {
            int tileNum     = (int)(tile_list.size() / 9);
            int newPointNum = (int)(newPoints.size() / 3);
            bool hasNormals = !newNormals.empty();
            for (int i = 0; i < tileNum; i++) {
                std::vector<AML::double3> triangle, triangleNorm;
                for (int j = 0; j < 3; j++) {
                    AML::double3 point;
                    point.x = tile_list[i*9+j*3+0];
                    point.y = tile_list[i*9+j*3+1];
                    point.z = tile_list[i*9+j*3+2];
                    triangle.push_back(point);
                    AML::double3 norm;
                    if (hasNormals) {
                        int vertexIdx = i*3+j;
                        int safeIdx   = (vertexIdx < newPointNum) ? vertexIdx : vertexIdx % newPointNum;
                        norm.x = newNormals[safeIdx*3+0];
                        norm.y = newNormals[safeIdx*3+1];
                        norm.z = newNormals[safeIdx*3+2];
                    }
                    triangleNorm.push_back(norm);
                }
                triangleCycle.push_back(triangle);
                triangleCycleNorm.push_back(triangleNorm);
            }
        }
        surface.push_back(triangleCycle);
        norms.push_back(triangleCycleNorm);

        delete[] points;
    }

    // Normal consistency pass (unchanged from original)
    const double SPATIAL_TOLERANCE = 1e-6;
    std::vector<std::vector<int>> arcToCycles(m_curveNet.arcs.size());
    for (int i = 0; i < (int)cycles.size(); i++)
        for (int j = 0; j < (int)cycles[i].size(); j++)
            arcToCycles[cycles[i][j].arcID].push_back(i);

    auto pointsAreClose = [SPATIAL_TOLERANCE](const AML::double3& p1, const AML::double3& p2) {
        return (p1 - p2).length() < SPATIAL_TOLERANCE;
    };
    auto findPt = [&](const std::vector<AML::double3>& tri, const AML::double3& p) {
        for (int i = 0; i < 3; i++) if (pointsAreClose(tri[i], p)) return i;
        return -1;
    };
    auto edgeOrient = [](int v1, int v2) {
        if ((v1==0&&v2==1)||(v1==1&&v2==2)||(v1==2&&v2==0)) return  1;
        if ((v1==1&&v2==0)||(v1==2&&v2==1)||(v1==0&&v2==2)) return -1;
        return 0;
    };
    auto flipPatch = [](TriangleCycle& sc, TriangleCycle& nc) {
        for (int t = 0; t < (int)sc.size(); t++) {
            std::swap(sc[t][0], sc[t][1]);
            for (int i = 0; i < 3; i++) nc[t][i] = -nc[t][i];
        }
    };

    std::vector<bool> processed(cycles.size(), false);
    if (!cycles.empty()) {
        processed[0] = true;
        std::vector<int> queue = {0};
        while (!queue.empty()) {
            int cur = queue.front(); queue.erase(queue.begin());
            for (int ai = 0; ai < (int)cycles[cur].size(); ai++) {
                int arcID = cycles[cur][ai].arcID;
                for (int adj : arcToCycles[arcID]) {
                    if (adj == cur || processed[adj]) continue;
                    std::vector<AML::double3> arcPts;
                    int sn = m_curveNet.arcs[arcID].endNodesID.first;
                    int en = m_curveNet.arcs[arcID].endNodesID.second;
                    arcPts.push_back(m_curveNet.nodes[sn].pos);
                    for (auto& p : m_curveNet.arcs[arcID].vertexList) arcPts.push_back(p);
                    arcPts.push_back(m_curveNet.nodes[en].pos);
                    if (arcPts.size() < 2) continue;
                    bool foundEdge = false, needsFlip = false;
                    for (int t1 = 0; t1 < (int)surface[cur].size() && !foundEdge; t1++) {
                        for (int t2 = 0; t2 < (int)surface[adj].size() && !foundEdge; t2++) {
                            for (int p = 0; p < (int)arcPts.size()-1 && !foundEdge; p++) {
                                int i1 = findPt(surface[cur][t1], arcPts[p]);
                                int i2 = findPt(surface[cur][t1], arcPts[p+1]);
                                int j1 = findPt(surface[adj][t2], arcPts[p]);
                                int j2 = findPt(surface[adj][t2], arcPts[p+1]);
                                if (i1>=0 && i2>=0 && j1>=0 && j2>=0) {
                                    foundEdge = true;
                                    int o1 = edgeOrient(i1, i2), o2 = edgeOrient(j1, j2);
                                    if (o1 != 0 && o2 != 0 && o1 == o2) needsFlip = true;
                                }
                            }
                        }
                    }
                    if (foundEdge && needsFlip) flipPatch(surface[adj], norms[adj]);
                    processed[adj] = true;
                    queue.push_back(adj);
                }
            }
        }
    }

    m_triangleSurface       = surface;
    m_triangleSurfaceNormal = norms;

    m_cycleDiscovery.m_isUpdateWeightTri   =
    m_cycleDiscovery.m_isUpdateWeightEdge  =
    m_cycleDiscovery.m_isUpdateWeightBiTri =
    m_cycleDiscovery.m_isUpdateWeightTriBd =
    m_cycleDiscovery.isUpdateWeightWorsDih =
    m_cycleDiscovery.m_isUpdateSurfaceSmooth =
    m_isUpdateSubdivision = m_isUpdateLaplacian = false;
    m_cycleDiscovery.m_isUpdateCycle = false;
}

CS_INLINE void cycleUtils::subDivideCycles()
{
    // Placeholder — subdivision logic to be implemented
    for (unsigned i = 0; i < m_cycleDiscovery.m_curveNet.arcs.size(); i++) {
        auto curve = m_cycleDiscovery.m_curveNet.arcs[i].vertexList;
        (void)curve;
    }
}