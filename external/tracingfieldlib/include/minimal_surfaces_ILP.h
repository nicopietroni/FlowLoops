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

#ifndef MINIMAL_SURFACES_ILP
#define MINIMAL_SURFACES_ILP

#include "box3.h"
#include "connectivity.h"
#include "point3.h"
#include <cassert>
#include <cotangent.h>
#include <cstddef>
#include <include/gurobi_c++.h>
#include <iostream>
#include <smooth_mesh.h>
#include <string>
#include <utility>
#include <vector>

template <class ScalarType> class MinimalSurfaceSolver {

  // enum InterpLapType { ILapUniform, ILapCot, ILapDynamic };

  std::vector<Geo::Point3<ScalarType>> &VertPos;
  std::vector<std::vector<int>> &Faces;

  std::vector<Geo::Point3<ScalarType>> SolvedPos;
  std::vector<int> Fixed;

  bool FixAnchor;
  bool UseCotangent;
  bool SecondStepWeight;
  bool ThirdStepForces;
  // bool FindForces;
  bool AnchorMaxDist;

  int MinAnchor;
  ScalarType WeightTolW;
  ScalarType DistW;
  ScalarType UseAnchorW;
  ScalarType MaxD;
  ScalarType MinSolutionGap;
  ScalarType TimeLimit;
  ScalarType MinForce;
  ScalarType MaxForce;
  ScalarType MinIntW;
  ScalarType MaxIntW;
  // InterpLapType ILapType;

  std::vector<std::vector<std::pair<int, ScalarType>>> VertLaplW;
  std::map<std::pair<int, int>, int> EdgeIdx;

  ScalarType SumWeight(const int &IndexV) {
    assert(IndexV >= 0);
    assert(IndexV < VertPos.size());
    ScalarType SumW = 0;
    for (size_t i = 0; i < VertLaplW[IndexV].size(); i++)
      SumW += VertLaplW[IndexV][i].second;
    return SumW;
  }

  int VertexVariable(const int &IndexV, const int &IndexC) {
    assert(IndexV >= 0);
    assert(IndexV < VertPos.size());
    assert(IndexC >= 0);
    assert(IndexC < 3);
    int IndexVar = IndexV * 3 + IndexC;
    return IndexVar;
  }

  void VertexComponent(const int &IndexVariable, int &IndexV, int &IndexC) {
    IndexV = IndexVariable / 3;
    IndexC = IndexVariable % 3;
  }

  int NumVertexVariables() { return (VertPos.size() * 3); }

  void InitLaplacianConstraint() {

    std::vector<std::vector<ScalarType>> CotWeight;
    Geo::CotangentWeightTriMesh(VertPos, Faces, CotWeight);
    assert(CotWeight.size() == Faces.size());

    // save cotangent per edge
    std::map<std::pair<int, int>, ScalarType> EdgeCtg;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        EdgeCtg[Key] = CotWeight[i][j];
      }

    std::vector<bool> IsB;
    Geo::getBorderVert(Faces, IsB);
    VertLaplW.clear();
    VertLaplW.resize(VertPos.size());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        ScalarType CurrW = 1;
        if (UseCotangent) {
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
          assert(EdgeCtg.count(Key) > 0);
          CurrW = EdgeCtg[Key];
        }
        if (!IsB[IndexV0])
          VertLaplW[IndexV0].push_back(
              std::pair<int, ScalarType>(IndexV1, CurrW));
        if (!IsB[IndexV1])
          VertLaplW[IndexV1].push_back(
              std::pair<int, ScalarType>(IndexV0, CurrW));
      }
    std::vector<std::pair<int, int>> BorderE;
    Geo::getBorderEdges(Faces, BorderE);
    for (size_t i = 0; i < BorderE.size(); i++) {
      int IndexV0 = BorderE[i].first;
      int IndexV1 = BorderE[i].second;
      ScalarType CurrW = 1;
      if (UseCotangent) {
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        assert(EdgeCtg.count(Key) > 0);
        CurrW = EdgeCtg[Key];
      }
      VertLaplW[IndexV0].push_back(std::pair<int, ScalarType>(IndexV1, 1));
      VertLaplW[IndexV1].push_back(std::pair<int, ScalarType>(IndexV0, 1));
    }
  }

  void InitEdgeIdx() {
    int CurrIdx = 0;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> EdgeK(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        if (EdgeIdx.count(EdgeK) == 0) {
          EdgeIdx[EdgeK] = CurrIdx;
          CurrIdx++;
        }
      }
  }

  // bounding box diagonal used to
  // tune weights and constraints
  ScalarType BBdiag;

  // vert position variables
  std::vector<GRBVar> x;
  // tolerance distance variable
  std::vector<GRBVar> d;
  // boolean variables for the anchors
  std::vector<GRBVar> b;
  // inteerpolation weights optimized
  std::vector<std::vector<GRBVar>> IWeights;
  //force variables
  std::vector<GRBVar> f;

  // Final per edge tension
  std::map<std::pair<int, int>, ScalarType> PerEdgeTens;

  // the quadratic objective
  GRBQuadExpr objective;

  void ResetGRBVariables() {
    x.clear();
    d.clear();
    b.clear();
    IWeights.clear();
    f.clear();
  }

  void InitVertPosVariables(GRBModel &model,
                            const std::set<int> &Fixed = std::set<int>()) {
    // x = std::vector<GRBVar>(NumVertexVariables());
    // for (int i = 0; i < NumVertexVariables(); i++) {
    //   int IndexV, IndexC;
    //   VertexComponent(i, IndexV, IndexC);
    //   x[i] =
    //       model.addVar(-GRB_INFINITY, GRB_INFINITY,
    //       VertPos[IndexV].V(IndexC),
    //                    GRB_CONTINUOUS, "x_" + std::to_string(i));
    // }
    x = std::vector<GRBVar>(NumVertexVariables());
    for (int i = 0; i < NumVertexVariables(); i++) {
      int IndexV, IndexC;
      VertexComponent(i, IndexV, IndexC);
      ScalarType CurrVal = SolvedPos[IndexV].V(IndexC);
      if ((Fixed.count(IndexV) > 0) && (FixAnchor))
        x[i] = model.addVar(CurrVal, CurrVal, CurrVal, GRB_CONTINUOUS,
                            "x_" + std::to_string(i));
      else
        x[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, CurrVal,
                            GRB_CONTINUOUS, "x_" + std::to_string(i));
    }
  }

  void InitInterpWVariables(GRBModel &model) {
    // interpolation weights variables
    IWeights =
        std::vector<std::vector<GRBVar>>(VertPos.size(), std::vector<GRBVar>());
    // one interp variable per vert per neigh
    for (int i = 0; i < VertPos.size(); i++) {
      int sizeNeigh = VertLaplW[i].size();
      IWeights[i].resize(sizeNeigh);
      for (int j = 0; j < IWeights[i].size(); j++)
        IWeights[i][j] =
            model.addVar(MinIntW, MaxIntW, 0.33333, GRB_CONTINUOUS, "baryW");
    }
  }

  void InitToleranceInterpolationVariables(GRBModel &model) {
    d = std::vector<GRBVar>(NumVertexVariables());
    for (int i = 0; i < NumVertexVariables(); i++) {
      d[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                          "d_" + std::to_string(i));
    }
  }

  void InitBooleanAnchorVariables(GRBModel &model) {
    b = std::vector<GRBVar>(VertPos.size());
    for (int i = 0; i < VertPos.size(); i++) {
      b[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "b_" + std::to_string(i));
    }
  }

  void InitForceVariables(GRBModel &model)
  {
    f=std::vector<GRBVar>(EdgeIdx.size());
      for (int i = 0; i < EdgeIdx.size(); i++) {
        f[i] = model.addVar(MinForce, MaxForce,(MaxForce + MinForce) / 2.0,
                            GRB_CONTINUOUS, "f_" + std::to_string(i));
      }
  }

  void SetMaxDistanceCondition(GRBModel &model, ScalarType DistVal) {
    // set max distance condition
    for (int i = 0; i < VertPos.size(); i++) {
      GRBQuadExpr QExpr = 0.0;
      int IndexV = i;
      // for each component
      for (size_t c = 0; c < 3; c++) {
        int IndexVar = VertexVariable(IndexV, c);
        ScalarType RestPos = VertPos[IndexV].V(c);
        QExpr += (x[IndexVar] - RestPos) * (x[IndexVar] - RestPos);
      }
      model.addQConstr(QExpr <= DistVal, "min_dist");
    }
  }

  void InitQuadrativObjectiveFunction(GRBModel &model, bool InterpTol,
                                      bool DistRest, bool Anchor) {
    objective = 0;
    // ADD INTERPOL TOLERANCE W
    // for (int i = 0; i < NumVertexVariables(); i++) {
    if (InterpTol)
      for (int i = 0; i < d.size(); i++) {
        objective += WeightTolW / BBdiag * (d[i] * d[i]);
      }

    // add distance term
    if ((DistRest) && (DistW > 0)) {
      for (int i = 0; i < VertPos.size(); i++) {
        int IndexV = i;
        for (int c = 0; c < 3; c++) {
          int IndexVar = VertexVariable(IndexV, c);
          ScalarType RestPos = VertPos[IndexV].V(c);
          objective += DistW / BBdiag * (x[IndexVar] - RestPos) *
                       (x[IndexVar] - RestPos);
        }
      }
    }

    // add the boolean
    if (Anchor)
      for (int i = 0; i < b.size(); i++) {
        objective += b[i] * UseAnchorW;
      }
  }

  bool SolveGurobi(GRBModel &model, bool UpdatePos, bool UpdateAnchor) {
    model.update();
    model.setObjective(objective, GRB_MINIMIZE);
    //    model.set(GRB_DoubleParam_Heuristics, 0.5);
    model.set(GRB_DoubleParam_MIPGap, MinSolutionGap); // 1% gap
    //   model.set(GRB_DoubleParam_MIPGapAbs, 10); // Absolute gap
    model.set(GRB_DoubleParam_TimeLimit, TimeLimit);
    //   model.set(GRB_IntParam_MIPFocus, 1);
    //   model.set(GRB_IntParam_Threads, 8); // Use 4 threads
    //   model.set(GRB_IntParam_Cuts, 2); // Moderate cut aggressiveness
    //   model.set(GRB_IntParam_NodeMethod, 2);
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
      std::cout << "Optimal solution found!\n";

    int sol_count = model.get(GRB_IntAttr_SolCount);
    if (sol_count > 0) {
      std::cout << " found solution!" << std::endl;

      // if updating positions
      if (UpdatePos && (x.size() > 0)) {
        SolvedPos = VertPos;
        for (int i = 0; i < x.size(); ++i) {
          for (int i = 0; i < x.size(); ++i) {
            double ValueVar = x[i].get(GRB_DoubleAttr_X);
            int IndexV, IndexC;
            VertexComponent(i, IndexV, IndexC);
            SolvedPos[IndexV].V(IndexC) = ValueVar;
          }
        }
      }

      if (UpdateAnchor && (b.size() > 0)) {
        Fixed.clear();
        for (int i = 0; i < b.size(); ++i) {

          double Valueb = b[i].get(GRB_DoubleAttr_X);
          if (Valueb > 0.5)
            Fixed.push_back(i);
        }
      }
      return true;
    } else {
      std::cout << " NOT found solution!" << std::endl;
      return false;
    }
  }

  bool OptimalAnchorStep() {
    // Create environment
    SolvedPos.clear();
    Fixed.clear();

    ResetGRBVariables();

    try {
      GRBEnv env = GRBEnv(true);
      env.set("LogFile", "minimalsurf.log");
      env.start();

      GRBModel model = GRBModel(env);

      InitVertPosVariables(model);

      InitToleranceInterpolationVariables(model);

      InitBooleanAnchorVariables(model);

      // ADD THE INTERPOLATION CONSTR INDICATOR
      // IF ANCHOR BOOLEAN ==1 THEN ACTIVATE THE
      // INTERPOLATION CONDITION

      // for each vertex
      for (int i = 0; i < VertPos.size(); i++) {
        int IndexV = i;
        // get sum of neighbours
        ScalarType SumW = SumWeight(IndexV);
        assert(SumW > 0);
        // for each component
        for (int c = 0; c < 3; c++) {

          // get the vertex variable considering the component
          int IndexVar = VertexVariable(IndexV, c);

          // sum up all neighbours by their weight
          GRBLinExpr avgExpr = 0.0;
          assert(VertLaplW[IndexV].size() > 0);
          // for each neighbour
          for (size_t n = 0; n < VertLaplW[IndexV].size(); n++) {
            // get index of vertex
            int IndexN = VertLaplW[IndexV][n].first;
            assert(IndexN != IndexV);
            // the interpolation weight
            ScalarType WeightN = VertLaplW[IndexV][n].second;
            assert(WeightN > 0);
            // and relative vertex component variable
            int VarNeigh = VertexVariable(IndexN, c);
            assert(VarNeigh >= 0);
            assert(VarNeigh < x.size());
            // sum up neigh times its weight
            avgExpr += (WeightN / SumW) * x[VarNeigh];
          }

          assert(IndexVar < d.size());
          assert(IndexVar < x.size());

          // then if not an anchor (b==0) then x should be interpolation of
          // neighbours plus a tolerance variable
          model.addGenConstrIndicator(b[IndexV], 0,
                                      x[IndexVar] - avgExpr + d[IndexVar],
                                      GRB_EQUAL, 0.0, "indicator_constraint");
        }
      }

      // ADD THE CONSTR INDICATOR
      // THAT FIX THE ANCHOR WHEN CHOSEN,
      // if the user want to fix anchor
      if (FixAnchor) {
        // for each vertex
        for (int i = 0; i < VertPos.size(); i++) {
          int IndexV = i;
          // for each component
          for (int c = 0; c < 3; c++) {
            // get the variable
            int IndexVar = VertexVariable(IndexV, c);
            // get the currnet component of its rest position
            ScalarType RestPos = VertPos[IndexV].V(c);
            // add the constraint
            model.addGenConstrIndicator(
                b[IndexV], 1, x[IndexVar], GRB_EQUAL, RestPos,
                "fixed_constraint" + std::to_string(IndexVar));
          }
        }
      }

      // SET THE MAXIMUM DISTANCE CONDITION
      // THE MAX VALUE IS PARAMETRIC WRT BBDIAGONAL
      ScalarType DistVal = pow(BBdiag * MaxD, 2);
      SetMaxDistanceCondition(model, DistVal);

      // SET THE CONSTRAINTS THAT
      // THE MINIMAL ANCHOR
      if (MinAnchor > 0) {
        GRBLinExpr fixExpr = 0.0;
        for (int i = 0; i < VertPos.size(); i++)
          fixExpr += b[i];
        model.addConstr(fixExpr >= MinAnchor, "Minimal Fix");
      }

      InitQuadrativObjectiveFunction(model, true, true, true);

      bool solved = SolveGurobi(model, true, true);
      return solved;

    } catch (GRBException &e) {
      std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (const std::exception &e) {
      std::cerr << "Standard exception: " << e.what() << "\n";
    }
    return false;
  }

  bool OptimalWeightStep() {
    std::set<int> FixedSet(Fixed.begin(), Fixed.end());
    try {
      GRBEnv env = GRBEnv(true);
      env.set("LogFile", "minimalsurf.log");
      env.start();

      GRBModel model = GRBModel(env);

      InitVertPosVariables(model, FixedSet);

      InitInterpWVariables(model);

      // InitToleranceInterpolationVariables(model);

      // ADD THE CONSTRAINT
      // THAT SUM OF WEIGHT SHOULD BE 1
      for (int i = 0; i < IWeights.size(); i++) {
        GRBLinExpr sumExpr = 0.0;
        for (size_t n = 0; n < IWeights[i].size(); n++)
          sumExpr += IWeights[i][n];
        model.addConstr(sumExpr == 1, "Sum Bary");
      }

      // ADD THE CONSTR INDICATOR
      // THAT EACH VERTEX (EXCEPT FIXED) EQUAL INTERPOLATION OF VALUES
      for (int i = 0; i < VertPos.size(); i++) {
        int IndexV = i;

        // no interpolation constraint
        if (FixedSet.count(IndexV) > 0)
          continue;

        // for each component
        for (int c = 0; c < 3; c++) {
          // int ki = (int)N[i].size();
          int IndexVar = VertexVariable(IndexV, c);

          GRBQuadExpr InterpExpr = 0.0;
          assert(VertLaplW[IndexV].size() > 0);
          for (size_t n = 0; n < VertLaplW[IndexV].size(); n++) {
            int IndexN = VertLaplW[IndexV][n].first;
            assert(IndexN != IndexV);
            int VarNeigh = VertexVariable(IndexN, c);
            assert(VarNeigh >= 0);
            assert(VarNeigh < x.size());
            assert(IndexV < IWeights.size());
            assert(n < IWeights[IndexV].size());
            InterpExpr += x[VarNeigh] * IWeights[IndexV][n];
          }
          assert(IndexVar >= 0);
          assert(IndexVar < x.size());
          // model.addQConstr(x[IndexVar] - InterpExpr + d[IndexVar] == 0,
          //                  "Interp Constraint");
          model.addQConstr(x[IndexVar] - InterpExpr == 0, "Interp Constraint");
        }
      }

      // SET THE MAXIMUM DISTANCE CONDITION
      // THE MAX VALUE IS PARAMETRIC WRT BBDIAGONAL
      ScalarType DistVal = pow(BBdiag * MaxD, 2);
      SetMaxDistanceCondition(model, DistVal);

      // add distance term
      GRBQuadExpr objective = 0;

      InitQuadrativObjectiveFunction(model, false, true, false);

      // NB: Quadratic equality constraints need this
      model.set(GRB_IntParam_NonConvex, 2);

      bool solved = SolveGurobi(model, true, false);
      return solved;
    } catch (GRBException &e) {
      std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (const std::exception &e) {
      std::cerr << "Standard exception: " << e.what() << "\n";
    }
    return false;
  }

  bool Findforces()
  {

    try {
      GRBEnv env = GRBEnv(true);
      env.set("LogFile", "minimalsurf.log");
      env.start();

      GRBModel model = GRBModel(env);

      // one variable per edge
      // std::map<std::pair<int, int>, int> EdgeIdx;

      // std::vector<GRBVar> f(EdgeIdx.size());
      // for (int i = 0; i < EdgeIdx.size(); i++) {
      //   f[i] = model.addVar(MinForce, (MaxForce + MinForce) / 2.0, MaxForce,
      //                       GRB_CONTINUOUS, "f_" + std::to_string(i));
      // }
      InitForceVariables(model);

      std::set<int> FixedSet(Fixed.begin(), Fixed.end());

      // for each component
      for (int c = 0; c < 3; c++) {
        // for each vertex
        for (int i = 0; i < VertPos.size(); i++) {
          int IndexV0 = i;
          if (FixedSet.count(IndexV0) > 0)
            continue;
          // add one linear eq constraint
          GRBLinExpr eqExpr = 0.0;

          // int IndexVar0=iVertexVariable(IndexV0, c);
          // for each neighbor
          for (size_t n = 0; n < VertLaplW[IndexV0].size(); n++) {
            int IndexVneigh = VertLaplW[IndexV0][n].first;
            // tension point outside the vertex
            Geo::Point3<ScalarType> Dir =
                (SolvedPos[IndexVneigh] - SolvedPos[IndexV0]);
            Dir.Normalize();
            // int IndexVar1=iVertexVariable(IndexVneigh, c);
            std::pair<int, int> Key(std::min(IndexV0, IndexVneigh),
                                    std::max(IndexV0, IndexVneigh));
            assert(EdgeIdx.count(Key) > 0);
            int IndexVarE = EdgeIdx[Key];
            assert(IndexVarE >= 0);
            assert(IndexVarE < f.size());
            eqExpr += Dir.V(c) * f[IndexVarE];
          }
          model.addConstr(eqExpr == 0, "Equilibrium");
        }
      }

      model.update();
      GRBLinExpr objective = 0;
      // for (int i = 0; i < EdgeIdx.size(); i++) {
      //   objective += f[i];
      // }
      model.setObjective(objective, GRB_MINIMIZE);
      //    model.set(GRB_DoubleParam_Heuristics, 0.5);
      model.set(GRB_DoubleParam_MIPGap, MinSolutionGap); // 1% gap
      model.optimize();

      if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        std::cout << "Optimal solution found!\n";

      int sol_count = model.get(GRB_IntAttr_SolCount);
      if (sol_count > 0) {
        std::cout << " found solution!" << std::endl;
        // SolvedPos = VertPos;
        // for (int i = 0; i < x.size(); ++i) {

        //   for (int i = 0; i < x.size(); ++i) {
        //     double ValueVar = x[i].get(GRB_DoubleAttr_X);
        //     int IndexV, IndexC;
        //     VertexComponent(i, IndexV, IndexC);
        //     SolvedPos[IndexV].V(IndexC) = ValueVar;
        //   }
        // }
        return true;
      } else {
        std::cout << " NOT found solution!" << std::endl;
        return false;
      }
    } catch (GRBException &e) {
      std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (const std::exception &e) {
      std::cerr << "Standard exception: " << e.what() << "\n";
    }
    return false;
  }

public:

  bool Solve(std::vector<int> &_Fixed,
             std::vector<Geo::Point3<ScalarType>> &_SolvedPos) {

    Geo::Box3<ScalarType> BBox;
    BBox.Init(VertPos);
    BBdiag = BBox.Diag();

    InitLaplacianConstraint();
    InitEdgeIdx();

    SolvedPos = VertPos;

    bool solvedAnchor = OptimalAnchorStep();
    if (!solvedAnchor)
      return false;

    if (SecondStepWeight)
      bool solvedWeight = OptimalWeightStep();

    bool ForcesOK=Findforces();
    if (ForcesOK)
      std::cout << "Solved Forces " << std::endl;
    else
      std::cout << "NOT Solved Forces " << std::endl;

    _Fixed = Fixed;
    _SolvedPos = SolvedPos;

    bool solvedOpt = false;
    bool solvedForces = false;

    ScalarType MaxDist = 0;

    std::cout << "*** RECAP ***" << std::endl;

    for (size_t i = 0; i < VertPos.size(); i++)
      MaxDist = std::max(MaxDist, (VertPos[i] - SolvedPos[i]).Norm());
    std::cout << "** MaxDist:" << MaxDist << " limit " << BBdiag * MaxD
              << std::endl;

    // if (solvedInit)
    //   std::cout << "** Solved Initial! **" << std::endl;
    // else
    //   std::cout << "** NOT Solved Initial! **" << std::endl;

    // if (solvedOpt)
    //   std::cout << "** Solved Optimal Interp W! **" << std::endl;
    // else
    //   std::cout << "** NOT Solved Optimal Interp W! **" << std::endl;

    // if (solvedForces)
    //   std::cout << "** Solved Forces! **" << std::endl;
    // else
    //   std::cout << "** NOT Solved Forces! **" << std::endl;

    return true; //(solvedInit && solvedOpt && solvedForces);
  }

  MinimalSurfaceSolver(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                       std::vector<std::vector<int>> &_Faces)
      : VertPos(_VertPos), Faces(_Faces) {
    FixAnchor = false;

    UseCotangent = true;

    SecondStepWeight = false;
    AnchorMaxDist = true;

    MinAnchor = 4;
    WeightTolW = 10000;
    DistW = 0.01;
    UseAnchorW = 1;
    MaxD = 0.03;
    MinSolutionGap = 0.1;
    TimeLimit = 60;
    MinForce = 0.1;
    MaxForce = 1;
    MinIntW = 0.2;
    MaxIntW = 0.8;
    // ILapType = ILapUniform;
    // ILapType = ILapCot;
  }
};

#endif
