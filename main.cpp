// including external gui and GL stuff

#include <GL/glew.h>
#include <random>
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <cstdio>
#define SAVE_STATUS_REMOVE
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>

// including std libs
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

// including TracineLib Stuff
#include <GL_functions/GLDrawField.h>
#include <GL_functions/GLDrawMesh.h>
#include <GL_functions/GLDrawOutline.h>

#include <IO/file_path_functions.h>
#include <IO/io_field.h>
#include <IO/io_obj.h>
#include <UI/orbitball.h>
#include <box3.h>
#include <connectivity.h>
#include <cross_field_functions.h>
#include <math_utils.h>
#include <mesh_features.h>
#include <edge_mesh_functions.h>
#include <mesh_functions.h>
#include <point3.h>
#include <smooth_cross_field.h>

// including FieldGraph Stuff
#include "loopweaver/curve_solver.h"
#include "loopweaver/curve_solver_interface.h"
#include <field_graph/graph.h>
#include <field_graph/mesh_preprocess.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_optimize.h>
#include <field_graph/path_functions.h>
#include <field_graph/path_sampling.h>
#include <field_graph/refine_for_tracing.h>

#include "loopweaver/loop_reconstruction_condition.h"
#include "loopweaver/loop_reconstruction_local_condition.h"
#include "loopweaver/reconstruction_error.h"
#include "loopweaver/mesh_preprocessing.h"
#include "loopweaver/settings.h"

#include "field_graph/normal_esteem.h"
#include "loopweaver/MechanicalIllustrativeShader.h"
#include <Space_Query/point_grid_3D.h>
#include <field_graph/normal_approx_decomposition_condition.h>
#include <field_graph/side_curvature_condition.h>
#include <field_graph/gaussian_curvature_condition.h>
#include <hausdorff.h>
#include <local_operations/tri_remesh_parameter.h>
#include <mesh_create.h>
#include <mesh_patch_decomposition.h>
#include <mesh_subdivide.h>
#include <mesh_symmetry.h>
#include <tangent_space_smooth.h>
#include <triangular_remesh.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "loopweaver/stb_image_write.h"

GLFWwindow *window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball, trackballUV;

ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

float line_thick = 1.f;
int selected_rend_item = 5;

int DrawMeshMode = 0;
int DrawColorMode = 0;
int OldDrawColorMode = -1;

bool drawOriginalBorders = true;
bool showMesh = true;
bool showFeatures = true;
bool showResult = false;
bool showCross = true;
bool showSing = true;
bool has_remeshed = false;
bool has_paths = false;
bool has_result = false;
bool showBoundaries = false;
bool ShowSymmPlane = false;
bool ShowOutline = true;
bool UseSymmetry = false;
bool UseCurvChangeNodes = false;
//bool final_extraction = false;
bool final_extraction = true; // the reconstruction is always done at the end
bool do_clean = true; // collapse bad faces and solve folds before the field (as the quadrangulator)
// graph arc costs weighted by the curvature anisotropy of the field (denser
// loops where the curvature is high and anisotropic)
bool curvature_graph = false;
// float curvature_graph_weight = 1;
float curvature_graph_weight = 2; // multiplier in [1, 1 + weight]

bool DrawCurv = false;
bool DrawColChange = false;
//bool DrawChangePos = true;
bool DrawChangePos = false; // curvature change nodes removed
typedef double ScalarType;

//bool use_toon_shader = true;
bool use_toon_shader = false; // mechanical illustrative shader not used
// used for visualization purposes
ScalarType AvEdge = 0;

// the symmetry plane used to split the mesh
Geo::Plane3<ScalarType> SymmetryPlane;

// MESH DATA
bool loadedMesh = false;

// mesh data and original mesh data
std::vector<Geo::Point3<ScalarType>> VertPos, VertPos0;
std::vector<std::vector<int>> Connectivity, Connectivity0;
std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;
std::vector<std::vector<int>> NextF, NextE;
Geo::Box3<ScalarType> MeshBox;

std::vector<Geo::Point3<ScalarType>> SolvedVertPos;
std::vector<std::vector<int>> SolvedConnectivity;
std::vector<Geo::Point3<ScalarType>> SolvedFaceNormals;
std::vector<Geo::Point3<ScalarType>> SolvedVertNormals;

// cross field data
bool processForTracing = true;
std::vector<Field::CrossF<ScalarType>> VertCurv;
std::vector<Field::CrossF<ScalarType>> FaceCurv;
bool has_cross_field = false;
// the field has been loaded from file (it refers to the input mesh: kept
// across the runs, while a computed field is recomputed)
bool has_loaded_field = false;
std::vector<Field::CrossF<ScalarType>> FaceCurvLoaded, VertCurvLoaded;
int KernelNring = 3;
ScalarType minQCrossVert, maxQCrossVert;
ScalarType minQCrossFace, maxQCrossFace;
std::vector<int> SingIndex, SingValue;
std::vector<int> FeatureSingIndex, FeatureSingValue;
std::vector<std::pair<int, int>> RemainingF;
std::vector<int> NewBorder;

FieldSmoothParam<ScalarType> SParam;

// Interface values
float GlobalSmoothVal = (float)SParam.GlobalSmoothVal;

std::vector<std::pair<int, int>> Features;
std::vector<int> Corners;

std::vector<std::vector<std::pair<int, int>>> Boundary;
std::vector<int> BoundaryVerts;

bool loop_recon_cond = true;
//bool loop_recon_local = false; // per patch version (IsCorrect only, no global update)
bool loop_recon_local = true; // per patch version (IsCorrect only, no global update)
int FinalSurfacingIterations = 10; // iterations of the final extraction (also used by the local condition)
bool final_no_smooth = false;       // final extraction without smoothing the patches (test)
// surface = thin-plate fill tangent to the border normals, no curve surfacing
bool cond_biharmonic_only = false;  // in the loop reconstruction condition
bool final_biharmonic_only = false; // in the final extraction
bool normal_approx_cond = false;
bool side_curvature_cond = true;
// developability (total absolute Gaussian curvature of a patch); default false
bool gaussian_curvature_cond = false;
bool beier_error_cond = true;
bool saved_last_screenshot = false;
int SmoothPathSteps = 20;

ScalarType maxErrRatio = 0.02;
ScalarType OldmaxErrRatio = maxErrRatio;

ScalarType maxNormAngle = 30;        // 45;
//ScalarType maxAnglePercentile = 0.2; // 0.5;
ScalarType maxAnglePercentile = 0.5; // lower values make the condition too noisy
ScalarType OldmaxNormAngle = maxNormAngle;
// the loop reconstruction condition also checks the normal error (it can
// split into too small patches, so it can be disabled)
bool use_normal_error = true;

Geo::Point3<ScalarType> mesh_color(0.9, 0.9, 0.9);
Geo::Point3<ScalarType> solved_mesh_color(0.4, 0.8, 0.66);

// error and the rest
std::vector<ScalarType> ErrorTarget;
std::vector<ScalarType> ErrorReconstructed;

std::vector<ScalarType> ErrorNormTarget;
std::vector<ScalarType> ErrorNormReconstructed;

std::vector<Geo::Point3<ScalarType>> FaceColorTarget;
std::vector<Geo::Point3<ScalarType>> FaceColorReconstructed;

//int numSamples = 10000;
int numSamples = 6000; // as the quadrangulator
int rendercycles = 0;
//bool dynamic_updates = true;
//bool dynamic_updates = false;
bool dynamic_updates = true;
// the condition used for loop reconstruction
LoopReconstructionCondition<ScalarType> LoopCond(SolvedVertPos,
                                                 SolvedConnectivity, Features);
// local version, configured from LoopCond when selected
LoopReconstructionLocalCondition<ScalarType>
    LoopCondLocal(SolvedVertPos, SolvedConnectivity, Features);

// final reconstruction: all patches together (global) or one at a time
bool final_per_patch = false;
std::vector<Geo::Point3<ScalarType>> FinalLocalUnusedV;
std::vector<std::vector<int>> FinalLocalUnusedF;
LoopReconstructionLocalCondition<ScalarType>
    FinalLocalSolver(FinalLocalUnusedV, FinalLocalUnusedF, Features);
// the displayed mesh has been mirrored back (symmetry)
bool mesh_reassembled = false;
// the paths on the mirrored mesh (both halves), saved instead of the
// decomposition ones when the mesh has been reassembled
std::vector<std::vector<int>> ReassembledVertPaths;
// total length of the layout paths (written with the timings)
ScalarType TotEdgeLen = 0;
// the default parameters (settings file), see RegisterSettings
LoopWeaverSettings Settings;

ScalarType NormalCondError = 20.0;
ScalarType NormalCondPercent = 0.2;
// Geo::NormalApproxCondition<ScalarType> NormalCond(NormalCondError,
//                                                   NormalCondPercent);
Geo::NormalLoopReconstructionCondition<ScalarType>
    NormalCond(Features, NormalCondError, NormalCondPercent);

//ScalarType MaxSideSumAngle = 180;
ScalarType MaxSideSumAngle = 120; // max turning of a patch side (degrees)
Geo::SideCurvatureCondition<ScalarType> SideAngleCond(MaxSideSumAngle);
// Taubin steps on the side polyline before measuring its turning (removes
// the zig-zag along the mesh edges); 0: not used
int SideSmoothSteps = 10;

// max total absolute Gaussian curvature of a patch (degrees); denoising:
// diffusion steps of the angle defect (Taubin removed)
ScalarType MaxGaussianCurvature = 90;
// int GaussianTaubinSteps = 0;
int GaussianDefectSmooth = 5;
Geo::GaussianCurvatureCondition<ScalarType>
    GaussCurvCond(MaxGaussianCurvature, GaussianDefectSmooth);

ScalarType MaxBezierErrorPerc = 2;
Geo::SideCurvatureBezierCondition<ScalarType> SideBezierCond;

// new version
MeshPatchDecomposition<ScalarType> MPatchDeco(VertPos, Connectivity, FaceCurv,
                                              Features);

MeshPatchDecomposition<ScalarType> MPatchDecoChangeP(VertPos, Connectivity,
                                                     FaceCurv, Features);

MechanicalIllustrativeShader *illustrative;

// CurveSolverInterface<ScalarType>::ExtractParam FinalExtrParam;

// std::vector<std::vector<std::vector<int>>> PatchNormalVerts;
std::vector<std::vector<std::vector<std::pair<int, int>>>> PatchNormalEdges;
std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> PatchNormals;

bool get_screenshot_original = false;
bool get_screenshot_boundaries = false;

bool do_Batch_compute = false;

std::vector<int> SolvedPatchIndex;

int AngleThrInt = 30;
int LernelCurvInt = 5;

std::vector<Geo::Point3<ScalarType>> CurvChangePos;
std::vector<Geo::Point3<ScalarType>> CurvChangeDir;

std::vector<Geo::Point3<ScalarType>> CurvEdgeMeshPos;
std::vector<std::pair<int, int>> CurvEdges;
std::vector<Geo::Point3<ScalarType>> CurvEdgeColors;
std::vector<Geo::Point3<ScalarType>> ClusterEdgeColors;

CurveSolver<ScalarType> CurveSolv(Features, SolvedVertPos, SolvedConnectivity,
                                  SolvedPatchIndex);

int MsecsField = 0;
int MsecsTracing = 0;
int MsecsSurfacing = 0;
int MsecsTotal = 0;

// BASE CONDITIONS FOR LOOPWEAVER
void InitDefaultParam() {

  // not used in the end using same solver at the end

  LoopCond.CurveSolv.use_original_meshing =
      true; // use_original_meshing_for_optimization;
  LoopCond.CurveSolv.resample_paths = true;

  CurveSolv.use_original_meshing =
      true; // use_original_meshing_for_final_extraction;
  CurveSolv.resample_paths = true;

  MPatchDeco.single_boundary_cond = true;
  MPatchDeco.disk_like_cond = true;
  MPatchDeco.self_connection_cond = true;
  MPatchDeco.single_sing_cond = false;
  MPatchDeco.match_sing_cond = false;
  MPatchDeco.prefer_feature_features = false;
  MPatchDeco.split_removal = false; // "Split when final remove", off by default
  // MPatchDeco.split_removal = true; // split at crossings and remove: stable, on by default
  // MPatchDeco.force_all_features = true; // as the quadrangulator
  MPatchDeco.force_all_features = false; // "Complete all features", off by default
  MPatchDeco.progress_every = 5; // progress message every 5 added/removed paths
  MPatchDeco.relax_paths = true; // spread the paths at the end
  MPatchDeco.CCAbility = -1;
  MPatchDeco.MinSides = 2; // 3;
  MPatchDeco.MaxSides = 6;
  beier_error_cond = true;
  loop_recon_cond = true;
  // dynamic_updates = true;
  // dynamic_updates = false; // "Dynamic Update", off by default
  dynamic_updates = true; // "Dynamic Update", on by default
  MPatchDeco.dynamicSmoothing = true;

  // copy all flags to the change point version
  MPatchDecoChangeP.single_boundary_cond = MPatchDeco.single_boundary_cond;
  MPatchDecoChangeP.disk_like_cond = MPatchDeco.disk_like_cond;
  MPatchDecoChangeP.self_connection_cond = MPatchDeco.self_connection_cond;
  MPatchDecoChangeP.single_sing_cond = MPatchDeco.single_sing_cond;
  MPatchDecoChangeP.match_sing_cond = MPatchDeco.match_sing_cond;
  MPatchDecoChangeP.prefer_feature_features =
      MPatchDeco.prefer_feature_features;
  MPatchDecoChangeP.split_removal = MPatchDeco.split_removal;
  MPatchDecoChangeP.CCAbility = MPatchDeco.CCAbility;
  MPatchDecoChangeP.MinSides = MPatchDeco.MinSides;
  MPatchDecoChangeP.MaxSides = MPatchDeco.MaxSides;
  MPatchDecoChangeP.dynamicSmoothing = MPatchDeco.dynamicSmoothing;
}

void UpdateFaceColor() {
  // if no error computed force use constant color
  std::cout << "Target size " << ErrorTarget.size() << " Connectivity size "
            << Connectivity.size() << std::endl;
  std::cout << "Reconstructed size " << ErrorReconstructed.size()
            << " SolvedConnectivity size " << SolvedConnectivity.size()
            << std::endl;

  if (ErrorTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorReconstructed.size() != SolvedConnectivity.size())
    DrawColorMode = 0;

  if (ErrorNormTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorNormReconstructed.size() != SolvedConnectivity.size())
    DrawColorMode = 0;

  OldDrawColorMode = DrawColorMode;

  std::cout << "Updating face color mode " << DrawColorMode << std::endl;
  if (DrawColorMode == 0) {
    FaceColorTarget =
        std::vector<Geo::Point3<ScalarType>>(Connectivity.size(), mesh_color);
    FaceColorReconstructed = std::vector<Geo::Point3<ScalarType>>(
        SolvedConnectivity.size(), solved_mesh_color);
    return;
  }

  if (DrawColorMode == 1) {
    ScalarType MaxDistance = MeshBox.Diag() * maxErrRatio;
    GetColorByScalar<ScalarType>(ErrorTarget, MaxDistance, 0, FaceColorTarget);
    GetColorByScalar<ScalarType>(ErrorReconstructed, MaxDistance, 0,
                                 FaceColorReconstructed);
    // color the one over the error in black
    for (size_t i = 0; i < ErrorTarget.size(); i++) {
      if (ErrorTarget[i] >= MaxDistance)
        FaceColorTarget[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    for (size_t i = 0; i < ErrorReconstructed.size(); i++) {
      if (ErrorReconstructed[i] >= MaxDistance)
        FaceColorReconstructed[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    return;
  }
  if (DrawColorMode == 2) {
    GetColorByScalar<ScalarType>(ErrorNormTarget, maxNormAngle, 0,
                                 FaceColorTarget);
    GetColorByScalar<ScalarType>(ErrorNormReconstructed, maxNormAngle, 0,
                                 FaceColorReconstructed);

    // color the one over the error in black
    for (size_t i = 0; i < ErrorNormTarget.size(); i++) {
      if (ErrorNormTarget[i] >= maxNormAngle)
        FaceColorTarget[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    for (size_t i = 0; i < ErrorNormReconstructed.size(); i++) {
      if (ErrorNormReconstructed[i] >= maxNormAngle)
        FaceColorReconstructed[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    return;
  }
}

void UpdateAfterRemesh() {

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert, maxQCrossVert);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace, maxQCrossFace);

  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);
  UpdateFaceColor();
}

int ero_steps = 3;
int dil_steps = 3;
ScalarType feature_angle = 50;//45; // 30;
ScalarType olfeature_angle = feature_angle;

void UpdateFeatures() {
  Features.clear();
  Corners.clear();
  Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
      VertPos, Connectivity, feature_angle, Features, Corners, ero_steps,
      dil_steps);
}

void InitSymmetryPlane() {
  Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
  Geo::Point3<ScalarType> Dir(1, 0, 0);
  SymmetryPlane.Init(Center, Dir);
}

// test switches of the LoopWeaver panel
bool show_test_normals = false;
bool show_test_bezier = false;
std::vector<Geo::Polyline3<ScalarType>> BezierPolylines;
std::vector<Geo::Point3<ScalarType>> BezierColorError;

void TestBezierPaths() {
  MPatchDeco.GetAllSidesPolyLines(BezierPolylines);
  std::vector<ScalarType> SideErrors;
  ScalarType AllMaxError = 0;
  for (size_t i = 0; i < BezierPolylines.size(); i++) {
    Geo::Polyline3<ScalarType> SmoothBezier;
    SmoothBezier =
        Geo::BezierFitting<ScalarType>::SampleFittedQuadraticBezierPolyline(
            BezierPolylines[i]);
    ScalarType MaxError = 0;
    ScalarType PolyL = BezierPolylines[i].Lenght();
    for (size_t j = 0; j < SmoothBezier.PolyPos.size(); j++) {
      ScalarType CurrError =
          (SmoothBezier.PolyPos[j] - BezierPolylines[i].PolyPos[j]).Norm();
      // CurrError=pow(CurrError,2)/pow(PolyL,2);
      MaxError = std::max(MaxError, CurrError);
    }
    AllMaxError = std::max(AllMaxError, MaxError);
    BezierPolylines[i] = SmoothBezier;
    SideErrors.push_back(MaxError);
  }
  ScalarType MaxError = MaxBezierErrorPerc * 0.01 * MeshBox.Diag();
  GetColorByScalar<ScalarType>(SideErrors, MaxError * 2, 0, BezierColorError);
  // set the one over the error in red
  for (size_t i = 0; i < SideErrors.size(); i++) {
    if (SideErrors[i] >= MaxError)
      BezierColorError[i] = Geo::Point3<ScalarType>(1, 0, 0);
  }
  std::cout << "Max Allowed Error Condition:" << MaxError << std::endl;
}

// //#define NEW_SIMMETRY_APPROACH
// #ifndef NEW_SIMMETRY_APPROACH
// void SplitMeshPlane() {
//   RefineMesh<ScalarType>::SplitMeshWithPlane(VertPos, Connectivity,
//                                              SymmetryPlane);

//   RemoveFacesOnNegativeSide<ScalarType>(VertPos, Connectivity, SymmetryPlane);

//   Geo::TriRemParam<ScalarType> RemP;
//   RemP.TargetL = AvEdge;

//   std::cout << "REMESHING" << std::endl;

//   RemP.feature_angle = feature_angle;
//   RemP.max_error = AvEdge * 0.5;
//   RemP.steps = 5;
//   Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);

//   UpdateFeatures();
// }
// #else
// void SplitMeshPlane() {
//   // Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
//   // SymmetryPlane.Init(Center, DirSymm);

//   ScalarType Tolerance = AvgEdgeLen(VertPos, Connectivity) * 0.1;
//   Geo::MeshSymmetry<ScalarType>::CutMeshWithSymmetryPlane(
//       VertPos, Connectivity, Features, Corners, SymmetryPlane, Tolerance);
//   UpdateFeatures();
//   ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
//   Geo::getFFAdjacency(Connectivity, NextF, NextE);

//   // Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
//   //     Connectivity, RemP.IsFaceEFeature, Features);
//   //
//   // Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(RemP.IsVertCorner,Corners);
// }
// #endif

// cut with the symmetry plane as the quadrangulation command line
// (MeshSymmetry::CutMeshWithSymmetryPlane, no remesh)
void SplitMeshPlane() {
  LoopWeaverPreprocess<ScalarType>::CutWithSymmetryPlane(
      VertPos, Connectivity, Features, Corners, SymmetryPlane);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
}

// mirror the reconstruction (symmetry); the errors, if any, are mirrored too
void ReassembleSolved() {
  if (!has_result)
    return;
  std::set<int> MiddleV;
  MirrorMesh<ScalarType>(SolvedVertPos, SolvedConnectivity, SymmetryPlane,
                         MiddleV);
  ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
                 SolvedVertNormals);
  // errors are recomputed on demand: back to constant color
  DrawColorMode = 0;
  UpdateFaceColor();
}

void ReassembleMesh() {
  // the paths of the half mesh, by position (mirrored below)
  std::vector<std::vector<int>> HalfVertPaths;
  MPatchDeco.GetVertexPaths(HalfVertPaths);
  std::vector<Geo::Point3<ScalarType>> HalfVertPos = VertPos;

  // mirror the mesh
  int NumV0 = VertPos.size();

  std::set<int> MiddleV;

  std::vector<std::vector<int>> Connectivity0 = Connectivity;
  MirrorMesh<ScalarType>(VertPos, Connectivity, SymmetryPlane, MiddleV);
  std::vector<Field::CrossF<ScalarType>> FaceCurvMirror = FaceCurv;
  MirrorCross<ScalarType>(FaceCurvMirror, SymmetryPlane);
  FaceCurv.insert(FaceCurv.end(), FaceCurvMirror.begin(), FaceCurvMirror.end());

  // transform boundaryes on Points pairs
  typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
  std::vector<std::vector<PosPair>> BoundaryPos;

  for (size_t i = 0; i < Boundary.size(); i++) {
    std::vector<PosPair> BPath;
    for (size_t j = 0; j < Boundary[i].size(); j++) {
      Geo::Point3<ScalarType> P0 = VertPos[Boundary[i][j].first];
      Geo::Point3<ScalarType> P1 = VertPos[Boundary[i][j].second];
      BPath.push_back(std::make_pair(P0, P1));
    }
    BoundaryPos.push_back(BPath);
  }

  // then save their position
  RemoveDuplicatedVert(VertPos, Connectivity);
  MergeCloseVert(VertPos, Connectivity, 1e-4);

  RemoveUnfererencedVert<ScalarType>(VertPos, Connectivity);

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  // then mirror the solved mesh
  ReassembleSolved();
  mesh_reassembled = true;
  // if (has_result) {
  //   std::set<int> MiddleV1;
  //   MirrorMesh<ScalarType>(SolvedVertPos, SolvedConnectivity, SymmetryPlane,
  //                          MiddleV);
  //   ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
  //                  SolvedVertNormals);
  //   // append errors to itself
  //   ErrorReconstructed.insert(ErrorReconstructed.end(),
  //                             ErrorReconstructed.begin(),
  //                             ErrorReconstructed.end());

  //   ErrorNormReconstructed.insert(ErrorNormReconstructed.end(),
  //                                 ErrorNormReconstructed.begin(),
  //                                 ErrorNormReconstructed.end());

  //   ErrorTarget.insert(ErrorTarget.end(), ErrorTarget.begin(),
  //                      ErrorTarget.end());

  //   ErrorNormTarget.insert(ErrorNormTarget.end(), ErrorNormTarget.begin(),
  //                          ErrorNormTarget.end());
  // }

  // then add the mirrored boundaries
  std::vector<std::vector<PosPair>> MirroredBoundaryPos;
  for (size_t i = 0; i < BoundaryPos.size(); i++) {
    std::vector<PosPair> MBPath;
    for (size_t j = 0; j < BoundaryPos[i].size(); j++) {
      Geo::Point3<ScalarType> P0 = BoundaryPos[i][j].first;
      Geo::Point3<ScalarType> P1 = BoundaryPos[i][j].second;

      // if too close to the plane to not add a duplicated edge
      if ((SymmetryPlane.Distance(P0) < 1e-6) &&
          (SymmetryPlane.Distance(P1) < 1e-6))
        continue;

      Geo::Point3<ScalarType> PM0 = SymmetryPlane.Mirror(P0);
      Geo::Point3<ScalarType> PM1 = SymmetryPlane.Mirror(P1);
      MBPath.push_back(std::make_pair(PM0, PM1));
    }
    MirroredBoundaryPos.push_back(MBPath);
  }

  BoundaryPos.insert(BoundaryPos.end(), MirroredBoundaryPos.begin(),
                     MirroredBoundaryPos.end());

  // the paths on the whole mesh (saved in the .path file)
  LoopWeaverPreprocess<ScalarType>::MirrorPaths(
      HalfVertPos, HalfVertPaths, SymmetryPlane, VertPos, MeshBox.Diag(),
      ReassembledVertPaths);

  Geo::PointSpatialIndex<ScalarType> VertexGrid;
  VertexGrid.Init(VertPos);
  Boundary.clear();
  Boundary.resize(BoundaryPos.size());
  for (size_t i = 0; i < BoundaryPos.size(); i++) {
    std::vector<int> BPath;
    for (size_t j = 0; j < BoundaryPos[i].size(); j++) {
      Geo::Point3<ScalarType> P0 = BoundaryPos[i][j].first;
      Geo::Point3<ScalarType> P1 = BoundaryPos[i][j].second;
      int v0 = -1;
      int v1 = -1;
      bool found0 = VertexGrid.GridClosest(VertPos, P0, MeshBox.Diag(), v0);
      if (!found0)
        continue;
      bool found1 = VertexGrid.GridClosest(VertPos, P1, MeshBox.Diag(), v1);
      if (!found1)
        continue;
      Boundary[i].push_back(std::pair<int, int>(v0, v1));
    }
  }

  //   for (size_t i = 0; i < Boundary.size(); i++) {
  //     for (size_t j = 0; j < Boundary.size(); j++) {
  //     Geo::Point3<ScalarType> P0 = VertPos[Boundary[i][j]].first];
  //     Geo::Point3<ScalarType> P1 = VertPos[Boundary[i][j]].second];

  //     PathPos.insert(std::make_pair(std::min(P0, P1), std::max(P0, P1)));

  //     // append the mirrored path
  //     Geo::Point3<ScalarType> PM0 = P0;
  //     if (SymmetryPlane.Distance(P0) > 1e-6)
  //       PM0 = SymmetryPlane.Mirror(P0);
  //     Geo::Point3<ScalarType> PM1 = P1;
  //     if (SymmetryPlane.Distance(P1) > 1e-6)
  //       PM1 = SymmetryPlane.Mirror(P1);
  //     PathPos.insert(std::make_pair(std::min(PM0, PM1), std::max(PM0,
  //     PM1)));
  //   }
  //   Boundary.clear();
  //   BoundaryVerts.clear();
  //   // then cycle over al faces to recreate the boundary
  //   for (size_t i = 0; i < Connectivity.size(); i++) {
  //     for (size_t j = 0; j < Connectivity[i].size(); j++) {
  //       int v0 = Connectivity[i][j];
  //       int v1 = Connectivity[i][(j + 1) % Connectivity[i].size()];
  //       Geo::Point3<ScalarType> P0 = VertPos[v0];
  //       Geo::Point3<ScalarType> P1 = VertPos[v1];
  //       auto search =
  //           PathPos.find(std::make_pair(std::min(P0, P1), std::max(P0,
  //           P1)));
  //       if (search != PathPos.end()) {
  //         Boundary.push_back(std::make_pair(v0, v1));
  //         BoundaryVerts.push_back(v0);
  //         BoundaryVerts.push_back(v1);
  //       }
  //     }
  //   }
  // }

  UpdateFaceColor();
  UpdateFeatures();
}

void InitFieldByCurvature() {
  SParam.FixedCross.clear();
  ComputeCurvatureField<ScalarType>(VertPos, Connectivity, FaceNormals,
                                    VertNormals, VertCurv, FaceCurv,
                                    KernelNring);
  Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert, maxQCrossVert);
  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);

  has_cross_field = true;
}

std::string folderProjectName;

void UpdateMeshFieldNormals() {
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  // update cross F poistion
  for (size_t i = 0; i < VertCurv.size(); i++)
    VertCurv[i].Pos = VertPos[i];

  AdaptFieldToMesh(VertPos, Connectivity, VertCurv, FaceCurv);
}

void SmoothPaths() {
  Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  std::vector<std::vector<int>> VertPaths;
  MPatchDeco.GetVertexPaths(VertPaths);
  // PathSampl.GetVertexSelectedPaths(VertPaths);
  //  Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

  Geo::PatchOptimize<ScalarType>::SmoothPaths(VertPos, Connectivity, VertPaths,
                                              Features, 0.5, SmoothPathSteps);

  UpdateMeshFieldNormals();
}

// void FinalExtractSurface() {
//   if (has_paths) {
//     // write last update mesh

//     // WriteOBJ("./test_solved_0.obj", LoopCond.DData.CurrSolvedVertPos,
//     //            LoopCond.DData.CurrSolvedFaces);

//     // WriteOBJ("./test_solved_1.obj", LoopCond.SolvedVertPos,
//     //            LoopCond.SolvedFaces);

//     // WriteOBJ("./test_solved_2.obj", SolvedVertPos,SolvedConnectivity);
//     // ScalarType MaxAbsError = MeshBox.Diag() * MaxBezierErrorPerc * 0.01;
//     // LoopCond.Init(MaxAbsError, maxNormAngle, maxAnglePercentile);
//     // LoopCond.CurveSolv.subsample_factor = 1;
//     LoopCond.CurveSolv.use_previous_solution_as_initial = false;
//     LoopCond.CurveSolv.FileName = GetFilanameNoExtension(PathMesh);

//     CurveSolv.use_previous_solution_as_initial = true;
//     CurveSolv.FileName = GetFilanameNoExtension(PathMesh);

//     // LoopCond.CurveSolv.smooth_pdeco_steps = 0;
//     // LoopCond.CurveSolv.save_patch_meshes = false;

//     // LoopCond.CurveSolv.only_updated_patches = false;
//     // LoopCond.CurveSolv.iteration = 10;

//     // LoopCond.CurveSolv.iteration = 10;
//     if (loop_recon_cond) {
//       // MPatchDeco.PatchManager().VertPos = LoopCond.CurrPmanVPos;
//       // MPatchDeco.PatchManager().Faces = LoopCond.CurrPmanFaces;
//       //LoopCond.CurveSolv.use_original_meshing = false;
//       // use_original_meshing_for_final_extraction;
//       LoopCond.CurveSolv.save_patch_meshes = false;
//       //LoopCond.CurveSolv.iteration = 10;
//       LoopCond.CurveSolv.iteration = FinalSurfacingIterations;
//       LoopCond.CurveSolv.stop_at_convergence = false; // final: all iterations
//       // if (final_no_smooth)
//       //   LoopCond.CurveSolv.smooth_original_meshing = false;
//       ScalarType OldRemeshFactor = LoopCond.CurveSolv.remesh_facctor;
//       LoopCond.CurveSolv.remesh_facctor = -1;
//       bool OldBiharmonicOnly = LoopCond.CurveSolv.biharmonic_only;
//       LoopCond.CurveSolv.biharmonic_only = final_biharmonic_only;
//       LoopCond.CurveSolv.UpdateSolvedMesh(MPatchDeco.PatchManager());
//       LoopCond.CurveSolv.biharmonic_only = OldBiharmonicOnly;
//       LoopCond.CurveSolv.remesh_facctor = OldRemeshFactor;
//       LoopCond.CurveSolv.stop_at_convergence = true; // back for the optimization

//       ErrorTarget = LoopCond.CurveSolv.TargetFDist;
//       // std::cout << "Target FDist size: " << ErrorTarget.size()
//       //           << std::endl;
//       // std::cout << "Connectivity size: " << Connectivity.size()
//       //           << std::endl;
//       ErrorReconstructed = LoopCond.CurveSolv.RemeshedFDist;
//       // std::cout << "Reconstructed FDist size: " << ErrorReconstructed.size()
//       //           << std::endl;
//       ErrorNormTarget = LoopCond.CurveSolv.TargetNErr;
//       // std::cout << "Target NErr size: " << ErrorNormTarget.size()
//       //           << std::endl;
//       ErrorNormReconstructed = LoopCond.CurveSolv.RemeshedNErr;
//     } else {
//       //CurveSolv.iteration = 10;
//       CurveSolv.iteration = FinalSurfacingIterations;
//       CurveSolv.stop_at_convergence = false; // final: all iterations
//       CurveSolv.use_original_meshing = true;
//       // use_original_meshing_for_final_extraction;
//       ScalarType OldRemeshFactor = CurveSolv.remesh_facctor;
//       CurveSolv.remesh_facctor = -1;
//       CurveSolv.biharmonic_only = final_biharmonic_only;
//       CurveSolv.UpdateSolvedMesh(MPatchDeco.PatchManager());
//       CurveSolv.remesh_facctor = OldRemeshFactor;
//       CurveSolv.stop_at_convergence = true;
//       ErrorTarget = CurveSolv.TargetFDist;
//       // std::cout << "Target FDist size: " << ErrorTarget.size()
//       //           << std::endl;
//       // std::cout << "Connectivity size: " << Connectivity.size()
//       //           << std::endl;
//       ErrorReconstructed = CurveSolv.RemeshedFDist;
//       // std::cout << "Reconstructed FDist size: " << ErrorReconstructed.size()
//       //           << std::endl;
//       ErrorNormTarget = CurveSolv.TargetNErr;
//       // std::cout << "Target NErr size: " << ErrorNormTarget.size()
//       //           << std::endl;
//       ErrorNormReconstructed = CurveSolv.RemeshedNErr;
//     }
//     // WriteOBJ("./test_solved_3.obj", SolvedVertPos,SolvedConnectivity);

//     // std::cout << "Reconstructed NErr size: " << ErrorNormReconstructed.size()
//     //           << std::endl;
//     ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
//                    SolvedVertNormals);

//     // UpdateMeshFieldNormals();
//     has_result = true;
//     showResult = true;
//     // DrawColorMode = 1;

//     UpdateFaceColor();
//   }
// }

// Final reconstruction from the patch layout, with the options of the
// reconstruction section (CurveSolv): all patches together or one at a time,
// curve surfacing or thin-plate fill only. Errors are not taken from here,
// they are computed on demand (Test menu).
void FinalExtractSurface() {
  if (!has_paths)
    return;
  CurveSolv.use_previous_solution_as_initial = true;
  CurveSolv.FileName = GetFilanameNoExtension(PathMesh);
  CurveSolv.iteration = FinalSurfacingIterations;
  CurveSolv.stop_at_convergence = false; // final: all iterations
  CurveSolv.biharmonic_only = final_biharmonic_only;
  if (final_per_patch) {
    FinalLocalSolver.SetSolverSettings(CurveSolv);
    FinalLocalSolver.SolveAllPatches(MPatchDeco.PatchManager(), SolvedVertPos,
                                     SolvedConnectivity, SolvedPatchIndex);
  } else
    CurveSolv.UpdateSolvedMesh(MPatchDeco.PatchManager());
  CurveSolv.stop_at_convergence = true;

  ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
                 SolvedVertNormals);
  has_result = true;
  showResult = true;
  DrawColorMode = 0;
  UpdateFaceColor();
  // the files exchanged with Curve-Surfacing-Core are not needed anymore
  CurveSolv.RemoveExchangeFiles();
}

void RefineForFieldComputation() {
  Geo::RefineForTracing<ScalarType>::RefineForFieldComputation(
      VertPos, Connectivity, Features);

  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
}
bool find_curv_change = false;

void BatchDecompose() {

  if (find_curv_change) {
    MPatchDecoChangeP.num_samples = numSamples;
    MPatchDecoChangeP.dynamicUpdates = dynamic_updates;
    MPatchDecoChangeP.InitConditions();

    MPatchDecoChangeP.ExtractPatches(Boundary, BoundaryVerts);
    return;
  }

  MPatchDeco.num_samples = numSamples;
  MPatchDeco.dynamicUpdates = dynamic_updates;
  MPatchDeco.curvature_aniso_weight =
      curvature_graph ? (ScalarType)curvature_graph_weight : 0;
  MPatchDeco.InitConditions();

  if (normal_approx_cond) {
    NormalCond.MaxErr = NormalCondError;
    NormalCond.MaxErrPercent = NormalCondPercent;
    MPatchDeco.AddExtraCondition(&NormalCond);
  }

  if (side_curvature_cond) {
    SideAngleCond.MaxSideAngle = MaxSideSumAngle;
    SideAngleCond.smooth_poly_steps = SideSmoothSteps;
    SideAngleCond.Init(VertPos, Connectivity);
    MPatchDeco.AddExtraCondition(&SideAngleCond);
  }
  if (gaussian_curvature_cond) {
    GaussCurvCond.MaxTotalCurvature = MaxGaussianCurvature;
    // GaussCurvCond.TaubinSteps = GaussianTaubinSteps;
    GaussCurvCond.DefectSmoothSteps = GaussianDefectSmooth;
    MPatchDeco.AddExtraCondition(&GaussCurvCond);
  }
  if (beier_error_cond) {
    // SideBezierCond.MaxErrorRatio = MaxBezierErrorPerc * 0.01;
    ScalarType MaxAbsError = MeshBox.Diag() * MaxBezierErrorPerc * 0.01;
    SideBezierCond.Init(MaxAbsError, VertPos, Connectivity);
    MPatchDeco.AddExtraCondition(&SideBezierCond);
  }

  if (loop_recon_cond) {
    ScalarType MaxAbsErr = MeshBox.Diag() * maxErrRatio;
    // LoopCond.Init(MaxAbsErr, maxNormAngle, maxAnglePercentile);
    // a negative angle disables the normal error check
    LoopCond.Init(MaxAbsErr, use_normal_error ? maxNormAngle : -1,
                  maxAnglePercentile);
    LoopCond.match_sing_cond = MPatchDeco.match_sing_cond;
    LoopCond.single_sing_cond = MPatchDeco.single_sing_cond;
    LoopCond.MinSides = MPatchDeco.MinSides;
    LoopCond.MaxSides = MPatchDeco.MaxSides;
    LoopCond.CurveSolv.biharmonic_only = cond_biharmonic_only;
    // MPatchDeco.AddExtraCondition(&LoopCond);
    if (loop_recon_local) {
      LoopCondLocal.CopySettingsFrom(LoopCond, FinalSurfacingIterations);
      LoopCondLocal.ClearCache();
      MPatchDeco.AddExtraCondition(&LoopCondLocal);
    } else
      MPatchDeco.AddExtraCondition(&LoopCond);
  }

  MPatchDeco.ExtractPatches(Boundary, BoundaryVerts);

  UpdateAfterRemesh();

  UpdateFaceColor();
  showBoundaries = true;
  has_paths = true;

  // MPatchDeco.RestoreOriginalPosOnSplitBorders();
}

void SmoothField() {
  auto start = std::chrono::high_resolution_clock::now();
  if (!has_cross_field)
    return;

  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

  SParam.Features = Features;

  SmoothGlobalPolyvector(VertPos, Connectivity, FaceNormals, FaceCurv, SParam);

  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);
  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);

  auto end = std::chrono::high_resolution_clock::now();
  MsecsField = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

}

void SaveTimings()
{
  //MsecsSurfacing = LoopCond.timeSolve;
  MsecsSurfacing =
      loop_recon_local ? LoopCondLocal.timeSolveLocal : LoopCond.timeSolve;
  MsecsTracing = MsecsTotal - MsecsField - MsecsSurfacing;
  std::string ProjName = GetFilanameNoExtension(PathMesh);
  std::string Timename = ProjName + std::string("_timings.txt");
  std::ofstream out(Timename);
  out << "Field: " << MsecsField << " ms" << std::endl;
  out << "Tracing: " << MsecsTracing << " ms" << std::endl;
  out << "Surfacing: " << MsecsSurfacing << " ms" << std::endl;
  out << "Total: " << MsecsTotal << " ms" << std::endl;
  out << "Edge Mesh Length: " << TotEdgeLen << std::endl;
  out.close();
}

void BatchProcessCurv() {
  if (processForTracing)
    RefineForFieldComputation();

  SParam.FixedCross.clear();
  InitFieldByCurvature();

  SmoothField();

  LoopWeaverPreprocess<ScalarType>::NormalizeFeatures(Features);

  UpdateAfterRemesh();
}

bool LoadMesh(const std::string &path) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (loadedMesh) {
    std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size()
              << " vertices" << std::endl;
    MeshBox.Init(VertPos);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);
    // LoopWeaverPreprocess<ScalarType>::CheckManifold(VertPos, Connectivity);
    // split the non-manifold elements, if any (the quadrangulator refuses
    // non-manifold meshes)
    if (LoopWeaverPreprocess<ScalarType>::MakeManifold(VertPos, Connectivity)) {
      ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
      Geo::getFFAdjacency(Connectivity, NextF, NextE);
    }
    AvEdge = AvgEdgeLen(VertPos, Connectivity);
    VertPos0 = VertPos;
    Connectivity0 = Connectivity;

    InitSymmetryPlane();
    UpdateFaceColor();
    UpdateFeatures();
    return true;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

bool LoadField(const std::string &path) {
  // std::cout<<"Loading Field"<<path.c_str()<<std::endl;
  has_cross_field = ReadField(VertPos, Connectivity, path, VertCurv, FaceCurv);
  if ((VertCurv.size() == 0) || (VertCurv.size() != VertPos.size())) {
    std::cout << "Init Vert Curv" << std::endl;
    SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv,
                         VertCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);
  }
  if (FaceCurv.size() == 0) {
    SetFaceCrossVectorFromVert(Connectivity, FaceNormals, VertCurv, FaceCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);
  }
  if (has_cross_field) {
    std::cout << "Loaded Field" << std::endl;
    showCross = true;
    showSing = true;

    Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
                                          maxQCrossVert);
    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
                         SingIndex, SingValue);
    return true;
  } else {
    std::cout << "Field not loaded" << std::endl;
    return false;
  }
}

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void cursor_position_callback(GLFWwindow *window, double xpos,
                                     double ypos) {
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  xpos *= xscale;
  ypos *= yscale;
  ypos = (display_h)-ypos;

  trackball.MouseMove((int)xpos, (int)ypos);
}

void mouse_button_callback(GLFWwindow *window, int button, int action,
                           int mods) {
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  xpos *= xscale;
  ypos *= yscale;
  ypos = (display_h)-ypos;

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    trackball.MouseDown((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    trackball.MouseUp((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
  trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {

  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
    trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
    trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);
}

void InitGLFW_Window() {
  if (!glfwInit()) {
    exit(0);
  }

  // In batch mode the window is only used as an offscreen render target for
  // the output screenshots, so keep it hidden
  if (do_Batch_compute)
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

  window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL2 example", NULL,
                            NULL);

  if (window == NULL) {
    exit(0);
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetKeyCallback(window, key_callback);
}

void InitIMGui() {
  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;

  // Setup Dear ImGui style
  ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL2_Init();
}

bool lineF = (SParam.NDir == 2);
bool crossF = (SParam.NDir == 4);

// void SetRenderBar() {
//   ImGui::Begin("Render", NULL, ImGuiWindowFlags_AlwaysAutoResize);

//   ImGui::Checkbox("Draw Shader", &use_toon_shader);
//   ImGui::Checkbox("Draw Symm Plane", &ShowSymmPlane);
//   ImGui::Checkbox("Draw Curv Color", &DrawCurv);
//   // curvature change nodes removed
//   // ImGui::Checkbox("Draw Change Color", &DrawColChange);
//   // ImGui::Checkbox("Draw Change Points", &DrawChangePos);
//   if (ImGui::CollapsingHeader("MESH")) {
//     ImGui::Checkbox("Draw Original", &showMesh);
//     ImGui::Checkbox("Draw Result", &showResult);
//     ImGui::Checkbox("Draw Original Borders", &drawOriginalBorders);
//     ImGui::Checkbox("Draw Paths", &showBoundaries);
//     ImGui::Checkbox("Draw Outline", &ShowOutline);
//     ImGui::Checkbox("Draw Features", &showFeatures);
//     static const char *itemsDrawMode[] = {"Smooth", "Flat", "Smooth Wire",
//                                           "Flat Wire"};

//     ImGui::Combo("Mesh Mode", &DrawMeshMode, itemsDrawMode,
//                  IM_ARRAYSIZE(itemsDrawMode));

//     // error coloring moved to the Test menu of the LoopWeaver panel
//     // static const char *ItemsColorMode[] = {"Constant", "Dist Error",
//     //                                        "Norm Error"};
//     // ImGui::Combo("Color Mode", &DrawColorMode, ItemsColorMode,
//     //              IM_ARRAYSIZE(ItemsColorMode));
//     //
//     // if (OldDrawColorMode != DrawColorMode) {
//     //   UpdateFaceColor();
//     // }
//     // OldDrawColorMode = DrawColorMode;
//   }

//   if (ImGui::CollapsingHeader("FIELD")) {
//     ImGui::Checkbox("Show Cross", &showCross);
//     ImGui::Checkbox("Show Singularities ", &showSing);
//   }
//   if (ImGui::CollapsingHeader("PATCHES")) {
//     ImGui::Checkbox("Show Boundaries", &showBoundaries);
//   }
//   ImGui::End();
// }

void SetRenderBar() {
  ImGui::Begin("Render", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  // (removed: shader switch, curvature color, original borders switch, the
  // PATCHES header that duplicated Draw Paths; shader and original borders
  // stay on)
  ImGui::Checkbox("Draw Symm Plane", &ShowSymmPlane);
  ImGui::Checkbox("Draw Original", &showMesh);
  ImGui::Checkbox("Draw Result", &showResult);
  ImGui::Checkbox("Draw Paths", &showBoundaries);
  ImGui::Checkbox("Draw Outline", &ShowOutline);
  ImGui::Checkbox("Draw Features", &showFeatures);
  static const char *itemsDrawMode[] = {"Smooth", "Flat", "Smooth Wire",
                                        "Flat Wire"};
  ImGui::Combo("Mesh Mode", &DrawMeshMode, itemsDrawMode,
               IM_ARRAYSIZE(itemsDrawMode));

  if (ImGui::CollapsingHeader("FIELD")) {
    ImGui::Checkbox("Show Cross", &showCross);
    ImGui::Checkbox("Show Singularities ", &showSing);
  }
  ImGui::End();
}

void GetStreenShotMesh(std::string type) {
  // write code to get a screeshot of the openGL window will call in the
  // rendering loop also possible disable the GUI for the screenshot
  std::string ProjName = GetFilanameNoExtension(PathMesh);
  std::string ScreenshotName = ProjName + "_" + type + std::string(".png");
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  std::vector<unsigned char> pixels(width * height * 3);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
  // flip the image vertically
  std::vector<unsigned char> flippedPixels(width * height * 3);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      for (int c = 0; c < 3; c++) {
        flippedPixels[((height - 1 - y) * width + x) * 3 + c] =
            pixels[(y * width + x) * 3 + c];
      }
    }
  }
  stbi_write_png(ScreenshotName.c_str(), width, height, 3, flippedPixels.data(),
                 width * 3);
  std::cout << "Screenshot saved to " << ScreenshotName << std::endl;
}
// ScalarType TotEdgeLen = 0; (declared above)

void SaveConfigFile() {
  std::string ProjName = GetFilanameNoExtension(PathMesh);
  std::string ConfigName = ProjName + std::string("_config.txt");
  FILE *fout = fopen(ConfigName.c_str(), "w");
  if (fout) {
    fprintf(fout, "Symmetry: %d\n", UseSymmetry ? 1 : 0);
    fprintf(fout, "Sharp Angle: %d\n", (int)feature_angle);
    // fprintf(fout, "Use Curvature Change Nodes: %d\n",
    //         UseCurvChangeNodes ? 1 : 0);

    fprintf(fout, "Sampling Density: %d\n", numSamples);
    fprintf(fout, "dynamic_updates: %d\n", dynamic_updates ? 1 : 0);
    fprintf(fout, "Loop Reconstruction Condition: %d\n",
            loop_recon_cond ? 1 : 0);
    if (loop_recon_cond) {
      fprintf(fout, " - maxErrRatio: %f\n", maxErrRatio);
      fprintf(fout, " - maxNormAngle: %f\n", maxNormAngle);
      fprintf(fout, " - maxAnglePercentile: %f\n", maxAnglePercentile);
      fprintf(fout, " - local: %d\n", loop_recon_local ? 1 : 0);
    }
    fprintf(fout, "Bezier Error Condition: %d\n", beier_error_cond ? 1 : 0);
    if (beier_error_cond) {
      fprintf(fout, " - MaxBezierErrorPerc: %f\n", MaxBezierErrorPerc);
    }
    fprintf(fout, "Side Curvature Condition: %d\n",
            side_curvature_cond ? 1 : 0);
    if (side_curvature_cond) {
      fprintf(fout, " - MaxSideSumAngle: %f\n", MaxSideSumAngle);
      fprintf(fout, " - SideSmoothSteps: %d\n", SideSmoothSteps);
    }
    fprintf(fout, "Gaussian Curvature Condition: %d\n",
            gaussian_curvature_cond ? 1 : 0);
    if (gaussian_curvature_cond) {
      fprintf(fout, " - MaxGaussianCurvature: %f\n", MaxGaussianCurvature);
      // fprintf(fout, " - GaussianTaubinSteps: %d\n", GaussianTaubinSteps);
      fprintf(fout, " - GaussianDefectSmooth: %d\n", GaussianDefectSmooth);
    }
    fprintf(fout, "Normal Approximation Condition: %d\n",
            normal_approx_cond ? 1 : 0);
    if (normal_approx_cond) {
      fprintf(fout, " - NormalCondError: %f\n", NormalCondError);
      fprintf(fout, " - NormalCondPercent: %f\n", NormalCondPercent);
    }
    fprintf(fout, "Edge Mesh Length: %f\n", TotEdgeLen);
  }
  fclose(fout);
}

// the saved paths lying on the sharp features (all their edges are feature
// edges of the saved mesh, detected as in UpdateFeatures)
void FindSharpPaths(const std::vector<std::vector<int>> &VertPaths,
                    std::vector<bool> &IsSharp) {
  std::vector<std::pair<int, int>> SharpE;
  std::vector<int> SharpCorners;
  Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
      VertPos, Connectivity, feature_angle, SharpE, SharpCorners, ero_steps,
      dil_steps);
  std::set<std::pair<int, int>> SharpSet;
  for (size_t i = 0; i < SharpE.size(); i++)
    SharpSet.insert(std::make_pair(std::min(SharpE[i].first, SharpE[i].second),
                                   std::max(SharpE[i].first, SharpE[i].second)));
  IsSharp.assign(VertPaths.size(), false);
  int NumSharp = 0;
  for (size_t i = 0; i < VertPaths.size(); i++) {
    int NumE = 0, NumSharpE = 0;
    for (size_t j = 0; j + 1 < VertPaths[i].size(); j++) {
      int V0 = VertPaths[i][j];
      int V1 = VertPaths[i][j + 1];
      if (V0 == V1)
        continue;
      NumE++;
      if (SharpSet.count(std::make_pair(std::min(V0, V1), std::max(V0, V1))))
        NumSharpE++;
    }
    IsSharp[i] = ((NumE > 0) && (NumSharpE == NumE));
    if (IsSharp[i])
      NumSharp++;
  }
  std::cout << "Sharp paths: " << NumSharp << " of " << VertPaths.size()
            << std::endl;
}

void SaveAll() {
  std::string ProjName = GetFilanameNoExtension(PathMesh);

  // SaveTimings(); (below, after the edge length)

  std::vector<std::vector<int>> VertPaths;
  MPatchDeco.GetVertexPaths(VertPaths);
  // with symmetry the paths of both halves on the mirrored mesh
  if (mesh_reassembled)
    VertPaths = ReassembledVertPaths;

  // save the path mesh file as obj
  std::vector<Geo::Point3<ScalarType>> EdgeMeshVertPos;
  std::vector<std::vector<int>> ConnectivityEdges;
  // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromIdxSequences(VertPos,VertPaths,EdgeMeshVertPos,ConnectivityEdges);
  Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromVertPairs(
      VertPos, Connectivity, Boundary, EdgeMeshVertPos, ConnectivityEdges);
  WriteOBJ(ProjName + std::string("_edge.obj"), EdgeMeshVertPos,
           ConnectivityEdges);
           
  TotEdgeLen = Geo::EdgeMeshFunctions<ScalarType>::Lenght(EdgeMeshVertPos, ConnectivityEdges);
  std::cout<<"*** Edge Mesh Length: "<<TotEdgeLen<<std::endl;

  SaveTimings();

  // OLD (partial list of the parameters)
  // SaveConfigFile();
  // all the parameters used, in the settings format (can be loaded back
  // with -settings <name>_settings.txt)
  Settings.Save(ProjName + std::string("_settings.txt"));

  std::string LoopsName = ProjName + std::string(".path");
  std::string FieldName = ProjName + std::string(".field");
  // WritePath(LoopsName, VertPaths);
  // the paths on the sharp features marked ("P <i> S")
  std::vector<bool> IsSharpPath;
  FindSharpPaths(VertPaths, IsSharpPath);
  WritePath(LoopsName, VertPaths, IsSharpPath);
  WriteField(FieldName, VertCurv, FaceCurv);
  WriteOBJ(ProjName + std::string("_remeshed.obj"), VertPos, Connectivity);
  WriteOBJ(ProjName + std::string("_solved.obj"), SolvedVertPos,
           SolvedConnectivity);
}

void ProcessAllLoops() {
  VertPos = VertPos0;
  Connectivity = Connectivity0;
  mesh_reassembled = false;

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  // SPLIT IF NEEDED
  if (UseSymmetry) {
    SplitMeshPlane();
  }
  // CLEAN (not with a loaded field, which refers to the input mesh)
  if (do_clean && (!has_cross_field)) {
    LoopWeaverPreprocess<ScalarType>::Clean(VertPos, Connectivity);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);
  }
  // features on the processed mesh
  if (UseSymmetry || (do_clean && (!has_cross_field)))
    UpdateFeatures();
  // FIND THE FIELD

  if (!has_cross_field)
    BatchProcessCurv();
  // PatchM.SetSingularities(SingIndex, SingValue);

  // BATCH DECOMPOSE
  // BatchDecompose();
  BatchDecompose();
  showCross = false;
  showSing = false;
  has_paths = true;

  // SmoothPaths();
  if (final_extraction)
    FinalExtractSurface();
}

void UpdateCurvatureChangeNodes() {

  bool old_loop_recon_cond = loop_recon_cond;
  bool old_normal_approx_cond = normal_approx_cond;
  bool old_side_curvature_cond = side_curvature_cond;
  bool old_gaussian_curvature_cond = gaussian_curvature_cond;
  bool old_beier_error_cond = beier_error_cond;

  loop_recon_cond = false;
  normal_approx_cond = false;
  side_curvature_cond = false;
  gaussian_curvature_cond = false;
  beier_error_cond = false;

  std::vector<Geo::Point3<ScalarType>> VertPosTest = VertPos;
  std::vector<std::vector<int>> ConnectivityTest = Connectivity;
  std::vector<Field::CrossF<ScalarType>> FaceCurvTest = FaceCurv;
  std::vector<std::pair<int, int>> FeaturesTest = Features;
  find_curv_change = true;
  ProcessAllLoops();
  find_curv_change = false;
  // GET CURVATURE CHANGE NODES ON THE EXTRACTED PATCHES
  ScalarType Ratio = (ScalarType)LernelCurvInt / (ScalarType)100;
  ScalarType KernelValue = MeshBox.Diag() * Ratio;
  ScalarType ClusterAngleDeg = static_cast<ScalarType>(AngleThrInt);
  std::vector<ScalarType> EdgeCurvDeg;
  std::vector<int> EdgeCluster;

  CurvChangePos.clear();
  CurvChangeDir.clear();
  MPatchDecoChangeP.GetCurvatureChangeClustersNodeDirection(
      KernelValue, feature_angle, ClusterAngleDeg, CurvEdgeMeshPos, CurvEdges,
      EdgeCurvDeg, EdgeCluster, CurvChangePos, CurvChangeDir);

  // VISUALIZATION
  ScalarType maxCurv =
      *std::max_element(EdgeCurvDeg.begin(), EdgeCurvDeg.end());
  ScalarType minCurv =
      *std::min_element(EdgeCurvDeg.begin(), EdgeCurvDeg.end());
  minCurv = std::min((ScalarType)0, minCurv);

  CurvEdgeColors.clear();
  for (size_t i = 0; i < EdgeCurvDeg.size(); i++) {
    Geo::Point3<ScalarType> CurrCol =
        Geo::ColorRGBRamp(minCurv, maxCurv, EdgeCurvDeg[i]);
    CurvEdgeColors.push_back(CurrCol);
  }

  int ClusterNum = *std::max_element(EdgeCluster.begin(), EdgeCluster.end());

  ClusterEdgeColors.clear();
  for (size_t i = 0; i < EdgeCluster.size(); i++) {
    Geo::Point3<ScalarType> CurrCol =
        Geo::ScatterColor<ScalarType>(ClusterNum, EdgeCluster[i]);
    ClusterEdgeColors.push_back(CurrCol);
  }
  // restore values
  VertPos = VertPosTest;
  Connectivity = ConnectivityTest;
  FaceCurv = FaceCurvTest;
  Features = FeaturesTest;

  Boundary.clear();
  BoundaryVerts.clear();
  has_paths = false;
  has_cross_field = false;
  has_remeshed = false;
  has_result = false;

  loop_recon_cond = old_loop_recon_cond;
  normal_approx_cond = old_normal_approx_cond;
  side_curvature_cond = old_side_curvature_cond;
  gaussian_curvature_cond = old_gaussian_curvature_cond;
  beier_error_cond = old_beier_error_cond;

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  UpdateFaceColor();

  MPatchDeco.PrioritySamplePos = CurvChangePos;
  MPatchDeco.PrioritySampleDir = CurvChangeDir;
}

// field, decomposition and final reconstruction (mirrored with symmetry)
// back to the loaded mesh (and loaded field, if any) before a new run: the
// previous run cut, cleaned and refined the mesh and computed a field on it
void ResetProcessing() {
  VertPos = VertPos0;
  Connectivity = Connectivity0;
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  UpdateFeatures();

  has_cross_field = has_loaded_field;
  if (has_loaded_field) {
    FaceCurv = FaceCurvLoaded;
    VertCurv = VertCurvLoaded;
    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
                         SingIndex, SingValue);
  } else {
    FaceCurv.clear();
    VertCurv.clear();
    SingIndex.clear();
    SingValue.clear();
  }

  has_paths = false;
  has_result = false;
  has_remeshed = false;
  mesh_reassembled = false;
  ReassembledVertPaths.clear();
  Boundary.clear();
  BoundaryVerts.clear();
  SolvedVertPos.clear();
  SolvedConnectivity.clear();
  DrawColorMode = 0;
  UpdateFaceColor();
  // the decomposer keeps the data of the previous decomposition
  MPatchDeco.ResetState();
}

void RunBatch() {
  // a previous run changed the mesh and the field: start again from the
  // loaded ones
  if (has_paths || has_result || mesh_reassembled ||
      (has_cross_field && (!has_loaded_field)))
    ResetProcessing();
  // the test drawings refer to the old decomposition
  show_test_normals = show_test_bezier = false;
  PatchNormalEdges.clear();
  PatchNormals.clear();
  BezierPolylines.clear();
  BezierColorError.clear();
  auto start = std::chrono::high_resolution_clock::now();
  // if (UseCurvChangeNodes)
  //   UpdateCurvatureChangeNodes();
  final_extraction = true;
  ProcessAllLoops();

  if (UseSymmetry)
    ReassembleMesh();
  auto end = std::chrono::high_resolution_clock::now();
  MsecsTotal = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

// recompute the final reconstruction only, with the current options
void UpdateReconstruction() {
  if (!has_paths)
    return;
  FinalExtractSurface();
  if (mesh_reassembled)
    ReassembleSolved();
}

void SaveReconstruction() {
  std::string ProjName = GetFilanameNoExtension(PathMesh);
  WriteOBJ(ProjName + std::string("_solved.obj"), SolvedVertPos,
           SolvedConnectivity);
  std::cout << "Saved " << ProjName << "_solved.obj" << std::endl;
}

// test: per-face reconstruction error (recomputed, not the conditions' one)
// shown as color on both meshes; mode 1 distance, 2 normal, 0 constant
void ColorByReconstructionError(int mode) {
  if ((mode != 0) && has_result) {
    ReconstructionError<ScalarType> Err;
    Err.Compute(VertPos, Connectivity, SolvedVertPos, SolvedConnectivity);
    ErrorTarget = Err.TargetDist;
    ErrorReconstructed = Err.SolvedDist;
    ErrorNormTarget = Err.TargetNormErr;
    ErrorNormReconstructed = Err.SolvedNormErr;
    std::cout << "Reconstruction error: max distance "
              << std::max(Err.MaxOf(Err.TargetDist), Err.MaxOf(Err.SolvedDist)) /
                     MeshBox.Diag()
              << " (x diag), max normal "
              << std::max(Err.MaxOf(Err.TargetNormErr), Err.MaxOf(Err.SolvedNormErr))
              << " deg" << std::endl;
  }
  DrawColorMode = has_result ? mode : 0;
  UpdateFaceColor();
}

// batch: process and save everything
void ProcessAll() {
  RunBatch();
  SaveAll();
  // ExtractDual();
  // get_screenshot_original = true;
}

// condition name: yellow, upper case
static bool ConditionCheckbox(const char *label, bool *value) {
  ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.85f, 0.1f, 1.0f));
  bool changed = ImGui::Checkbox(label, value);
  ImGui::PopStyleColor();
  return changed;
}

// remesh edge factor as an integer slider: -1 (no remesh), 1 .. 4
static void RemeshFactorSlider(const char *label, ScalarType &factor) {
  int idx = (factor <= 0) ? 0 : std::max(1, std::min(4, (int)std::round(factor)));
  ImGui::SliderInt(label, &idx, 0, 4, (idx == 0) ? "-1 (no remesh)" : "%d",
                   ImGuiSliderFlags_AlwaysClamp);
  factor = (idx == 0) ? ScalarType(-1) : ScalarType(idx);
}

// options of the surfacing tessellation, shared by condition and final
static void TessellationOptions(CurveSolver<ScalarType> &CS, const char *suffix) {
  std::string sfx(suffix);
  // ImGui::Checkbox(("Use original tessellation" + sfx).c_str(),
  //                 &CS.use_original_meshing);
  if (ImGui::Checkbox(("Use original tessellation" + sfx).c_str(),
                      &CS.use_original_meshing) &&
      !CS.use_original_meshing)
    CS.resample_paths = true; // default with the generated tessellation
  if (CS.use_original_meshing) {
    RemeshFactorSlider(("Remesh Edge Factor" + sfx).c_str(), CS.remesh_facctor);
  } else {
    // resample first: the subsampling then keeps evenly spaced points
    ImGui::Checkbox(("Resample Path" + sfx).c_str(), &CS.resample_paths);
    ImGui::SliderInt(("Subsample Factor (1 = none)" + sfx).c_str(),
                     &CS.subsample_factor, 1, 4, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    // ImGui::Checkbox(("Resample Path" + sfx).c_str(), &CS.resample_paths);
  }
  CS.MakeParametersCoherent();
}

void SetConditionsBar() {
  ImGui::Begin("ConditionS", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  // ImGui::Checkbox("Loop Reconstruction Condition", &loop_recon_cond);
  ConditionCheckbox("LOOP RECONSTRUCTION CONDITION", &loop_recon_cond);
  // global / local version, on the same line as the condition checkbox
  if (loop_recon_cond) {
    int loop_recon_mode = loop_recon_local ? 1 : 0;
    ImGui::SameLine();
    ImGui::RadioButton("Global", &loop_recon_mode, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Local", &loop_recon_mode, 1);
    loop_recon_local = (loop_recon_mode == 1);
  }
  if (loop_recon_cond) {
    //ImGui::Checkbox("Local (per patch)", &loop_recon_local);
    float maxErrRatiof = maxErrRatio;
    ImGui::SliderFloat("Reconstruction Error", &maxErrRatiof, 0.002f, 0.1f,
                       "%.3f", ImGuiSliderFlags_AlwaysClamp);
    maxErrRatio = maxErrRatiof;
    if (OldmaxErrRatio != maxErrRatio) {
      UpdateFaceColor();
      OldmaxErrRatio = maxErrRatio;
    }

    ImGui::Checkbox("Use normal error", &use_normal_error);
    ImGui::BeginDisabled(!use_normal_error);
    int maxNormAnglei = maxNormAngle;
    ImGui::SliderInt("Recon: Norm Angle Err", &maxNormAnglei, -1, 45, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    maxNormAngle = maxNormAnglei;
    if (OldmaxNormAngle != maxNormAngle) {
      UpdateFaceColor();
      OldmaxNormAngle = maxNormAngle;
    }

    float maxAnglePercentilef = maxAnglePercentile;
    ImGui::SliderFloat("Recon: Max Norm Angle Err % ", &maxAnglePercentilef,
                       0.0f, 0.5f, "%.1f", ImGuiSliderFlags_AlwaysClamp);
    maxAnglePercentile = maxAnglePercentilef;
    ImGui::EndDisabled();

    ImGui::Checkbox("Thin-plate only (condition)", &cond_biharmonic_only);
    // ImGui::Checkbox("Original Mesh 0",
    //                 &LoopCond.CurveSolv.use_original_meshing);
    // if (LoopCond.CurveSolv.use_original_meshing) {
    //   ImGui::Checkbox("Smooth Original Surface 0",
    //                   &LoopCond.CurveSolv.smooth_original_meshing);
    //   float remesh_factorf = LoopCond.CurveSolv.remesh_facctor;
    //   ImGui::InputFloat("Remesh Edge Factor", &remesh_factorf);
    //   LoopCond.CurveSolv.remesh_facctor = remesh_factorf;
    // } else {
    //   ImGui::InputInt("Subsample Factor", &LoopCond.CurveSolv.subsample_factor);
    //   ImGui::Checkbox("Resample Path", &LoopCond.CurveSolv.resample_paths);
    // }
    // LoopCond.CurveSolv.MakeParametersCoherent();
    TessellationOptions(LoopCond.CurveSolv, "");
  }
  ImGui::Separator();

  // ImGui::Checkbox("Side Curvature Condition", &side_curvature_cond);
  ConditionCheckbox("SIDE CURVATURE CONDITION", &side_curvature_cond);
  if (side_curvature_cond) {
    int MaxSideSumAnglei = MaxSideSumAngle;
    ImGui::SliderInt("Max Side Angle", &MaxSideSumAnglei, 30, 200, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    MaxSideSumAngle = MaxSideSumAnglei;
    ImGui::SliderInt("Side Smooth Steps", &SideSmoothSteps, 0, 30, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
  }
  ImGui::Separator();

  // developability: total absolute Gaussian curvature of a patch
  ConditionCheckbox("GAUSSIAN CURVATURE CONDITION", &gaussian_curvature_cond);
  if (gaussian_curvature_cond) {
    int MaxGaussi = MaxGaussianCurvature;
    ImGui::SliderInt("Max Gaussian Curvature", &MaxGaussi, 5, 360, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    MaxGaussianCurvature = MaxGaussi;
    // ImGui::SliderInt("Taubin Steps", &GaussianTaubinSteps, 0, 10, "%d",
    //                  ImGuiSliderFlags_AlwaysClamp);
    ImGui::SliderInt("Defect Smooth Steps", &GaussianDefectSmooth, 0, 20,
                     "%d", ImGuiSliderFlags_AlwaysClamp);
  }
  ImGui::Separator();

  // ImGui::Checkbox("Bezier Error Condition", &beier_error_cond);
  ConditionCheckbox("BEZIER ERROR CONDITION", &beier_error_cond);
  if (beier_error_cond) {
    float MaxBezierErrorf = MaxBezierErrorPerc;
    ImGui::SliderFloat("Max Bezier Error", &MaxBezierErrorf, 0.1f, 5.0f, "%.1f",
                       ImGuiSliderFlags_AlwaysClamp);
    MaxBezierErrorPerc = MaxBezierErrorf;
  }
  ImGui::Separator();

  // ImGui::Checkbox("Normal Approximation Condition", &normal_approx_cond);
  ConditionCheckbox("NORMAL APPROXIMATION CONDITION", &normal_approx_cond);
  if (normal_approx_cond) {

    int maxNormAngleI = NormalCondError;
    ImGui::SliderInt("Norm:Max Angle Err", &maxNormAngleI, 5, 45, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    NormalCondError = maxNormAngleI;

    float maxAnglePercentilef = NormalCondPercent;
    ImGui::SliderFloat("Norm:Max Angle Err % ", &maxAnglePercentilef, 0.05f,
                       1.f, "%.1f", ImGuiSliderFlags_AlwaysClamp);
    NormalCondPercent = maxAnglePercentilef;
  }
  ImGui::Separator();

  // number of sides per patch: one widget for the [min, max] range
  ImGui::DragIntRange2("Sides", &MPatchDeco.MinSides, &MPatchDeco.MaxSides,
                       0.1f, 1, 12, "Min: %d", "Max: %d",
                       ImGuiSliderFlags_AlwaysClamp);
  ImGui::Separator();

  float CCAbilityf = MPatchDeco.CCAbility;
  ImGui::SliderFloat("CCAbility", &CCAbilityf, -1.0f, 1.0f, "%.1f",
                     ImGuiSliderFlags_AlwaysClamp);
  MPatchDeco.CCAbility = CCAbilityf;

  ImGui::Checkbox("Single patch sing", &MPatchDeco.single_sing_cond);
  ImGui::Checkbox("Match sing values", &MPatchDeco.match_sing_cond);
  // ImGui::InputInt("Min Sides", &MPatchDeco.MinSides);
  // ImGui::InputInt("Max Sides", &MPatchDeco.MaxSides);

  ImGui::End();
}

// void SetToolBar() {
//   ImGui::Begin("Loop Weaver", NULL, ImGuiWindowFlags_AlwaysAutoResize);

//   bool OldUseSymmetry = UseSymmetry;
//   ImGui::Checkbox("Use Symmetry", &UseSymmetry);
//   ImGui::Checkbox("Use Curvature Change Nodes", &UseCurvChangeNodes);

//   int CurrA = feature_angle;
//   ImGui::SliderInt("Sharp Angle", &CurrA, 5, 180, "%d",
//                    ImGuiSliderFlags_AlwaysClamp);
//   feature_angle = (ScalarType)CurrA;
//   if (olfeature_angle != feature_angle) {
//     UpdateFeatures();
//     olfeature_angle = feature_angle;
//   }
//   ImGui::SliderFloat("Field Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
//                      ImGuiSliderFlags_AlwaysClamp);
//   SParam.GlobalSmoothVal = GlobalSmoothVal;

//   // ImGui::Checkbox("Single Patch Singularities", &single_sing_cond);
//   ImGui::Separator();
//   ImGui::InputInt("Sampling Density", &numSamples);
//   ImGui::Checkbox("Dynamic Update", &dynamic_updates);
//   ImGui::Checkbox("Split Removal", &MPatchDeco.split_removal);
//   ImGui::Checkbox("DO Final Extraction", &final_extraction);

//   if (ImGui::CollapsingHeader("Final Extraction")) {
//     ImGui::Checkbox("Thin-plate only (final)", &final_biharmonic_only);
//     // ImGui::Checkbox("Original Mesh Final  ", &CurveSolv.use_original_meshing);
//     // ImGui::Checkbox("Smooth Original Surface Final",
//     //                 &CurveSolv.smooth_original_meshing);
//     // ImGui::Checkbox("Resample Path Final", &CurveSolv.resample_paths);
//     // ImGui::InputInt("Subsample Factor Final", &CurveSolv.subsample_factor);
//     // CurveSolv.MakeParametersCoherent();
//     TessellationOptions(CurveSolv, " (final)");
//   }

//   ImGui::Separator();
//   if (ImGui::Button("Process and Save"))
//     ProcessAll();

//   if (ImGui::Button("Batch Process")) {
//     if (UseCurvChangeNodes)
//       UpdateCurvatureChangeNodes();
//     ProcessAllLoops();
//   }

//   if (ImGui::Button("Get Screenshot")) {
//     get_screenshot_original = true;
//   }

//   if (ImGui::Button("Test Normals")) {
//     bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
//         MPatchDeco.PatchManager(), Features, PatchNormalEdges, PatchNormals);

//     if (done)
//       std::cout << "Esteem Normal Correct" << std::endl;
//     else
//       std::cout << "Esteem Normal Error" << std::endl;
//   }

//   if (ImGui::Button("Test Bezier")) {
//     TestBezierPaths();
//   }

//   ImGui::SliderInt("Angle Threshold", &AngleThrInt, 0, 100);
//   ImGui::SliderInt("Kernal Curvature", &LernelCurvInt, 1, 10);

//   if (ImGui::Button("Test Curvature Change Nodes")) {
//     UpdateCurvatureChangeNodes();
//   }

//   if ((UseSymmetry) && (!OldUseSymmetry))
//     ShowSymmPlane = true;
//   if (!UseSymmetry)
//     ShowSymmPlane = false;

//   if (ImGui::CollapsingHeader("Symmetry")) {

//     if (ImGui::Button("Split Symm")) {
//       SplitMeshPlane();
//     }
//     if (ImGui::Button("Reasemble")) {
//       ReassembleMesh();
//     }
//   }

//   if (ImGui::CollapsingHeader("Field")) {

//     ImGui::Checkbox("Process mesh for Tracing", &processForTracing);

//     ImGui::SliderInt("Ring curvature", &KernelNring, 2, 5, "%d",
//                      ImGuiSliderFlags_AlwaysClamp);

//     ImGui::SliderFloat("Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
//                        ImGuiSliderFlags_AlwaysClamp);
//     SParam.GlobalSmoothVal = GlobalSmoothVal;

//     ImGui::Checkbox("Align Border", &SParam.align_borders);

//     if (ImGui::Button("FIND FIELD")) {
//       BatchProcessCurv();
//       // PatchM.SetSingularities(SingIndex, SingValue);
//     }
//   }

//   if (ImGui::CollapsingHeader("Decomposition")) {

//     if (has_cross_field) {

//       ImGui::InputInt("Smooth Path Steps", &SmoothPathSteps, 1);
//       if (SmoothPathSteps < 0)
//         SmoothPathSteps = 0;
//       if (SmoothPathSteps > 50)
//         SmoothPathSteps = 50;

//       if (ImGui::Button("Decompose ")) {
//         BatchDecompose();
//       }

//       if (ImGui::Button("Smooth Paths"))
//         SmoothPaths();
//     }
//   }
//   if (has_paths) {

//     if (ImGui::Button("Extract Surface")) {
//       FinalExtractSurface();
//     }
//   }
//   ImGui::Separator();
//   if (ImGui::CollapsingHeader("Save")) {
//     if (ImGui::Button("SAVE ALL")) {
//       SaveAll();
//     }
//   }
//   ImGui::End();
// }

void SetToolBar() {
  ImGui::Begin("Loop Weaver", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  bool OldUseSymmetry = UseSymmetry;
  ImGui::Checkbox("Use Symmetry", &UseSymmetry);
  if ((UseSymmetry) && (!OldUseSymmetry))
    ShowSymmPlane = true;
  if (!UseSymmetry)
    ShowSymmPlane = false;

  int CurrA = feature_angle;
  ImGui::SliderInt("Sharp Angle", &CurrA, 5, 180, "%d",
                   ImGuiSliderFlags_AlwaysClamp);
  feature_angle = (ScalarType)CurrA;
  if (olfeature_angle != feature_angle) {
    UpdateFeatures();
    olfeature_angle = feature_angle;
  }
  ImGui::SliderFloat("Field Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
                     ImGuiSliderFlags_AlwaysClamp);
  SParam.GlobalSmoothVal = GlobalSmoothVal;

  // decomposition
  ImGui::Separator();
  ImGui::Checkbox("Split when final remove", &MPatchDeco.split_removal);
  ImGui::Checkbox("Complete all features", &MPatchDeco.force_all_features);
  ImGui::Checkbox("Relax paths", &MPatchDeco.relax_paths);
  // (curvature priority removed: it did not help, see the commented lines)
  // // candidate priority x curvature anisotropy (loops around curved tubes)
  // ImGui::Checkbox("Curvature priority", &MPatchDeco.curvature_priority);
  // (curvature weighted graph removed from the panel: -curvgraph W only)
  // ImGui::Checkbox("Curvature weighted graph", &curvature_graph);
  // if (curvature_graph)
  //   ImGui::SliderFloat("Curvature weight", &curvature_graph_weight, 0.1f,
  //                      4.0f, "%.1f", ImGuiSliderFlags_AlwaysClamp);
  ImGui::InputInt("Sampling Density", &numSamples);
  ImGui::Checkbox("Dynamic Update", &dynamic_updates);

  // final reconstruction (always done at the end of the batch)
  ImGui::Separator();
  if (ImGui::CollapsingHeader("FINAL RECONSTRUCTION",
                              ImGuiTreeNodeFlags_DefaultOpen)) {
    int per_patch = final_per_patch ? 1 : 0;
    ImGui::RadioButton("Global", &per_patch, 0);
    ImGui::SameLine();
    ImGui::RadioButton("One patch at a time", &per_patch, 1);
    final_per_patch = (per_patch == 1);
    ImGui::Checkbox("Thin-plate only", &final_biharmonic_only);
    TessellationOptions(CurveSolv, " (final)");
    // only once the decomposition is available
    if (has_paths) {
      if (ImGui::Button("Update Reconstruction"))
        UpdateReconstruction();
    }
    // writes only the reconstructed mesh (<name>_solved.obj)
    if (has_result) {
      ImGui::SameLine();
      if (ImGui::Button("Save Reconstruction"))
        SaveReconstruction();
    }
  }

  // commands
  ImGui::Separator();
  if (ImGui::Button("Batch Process"))
    RunBatch();
  ImGui::SameLine();
  if (ImGui::Button("Save All"))
    SaveAll();

  ImGui::Separator();
  // if (ImGui::CollapsingHeader("Test")) {
  //   if (ImGui::Button("Normals")) {
  //     bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
  //         MPatchDeco.PatchManager(), Features, PatchNormalEdges, PatchNormals);
  //     if (done)
  //       std::cout << "Esteem Normal Correct" << std::endl;
  //     else
  //       std::cout << "Esteem Normal Error" << std::endl;
  //   }
  //   ImGui::SameLine();
  //   if (ImGui::Button("Bezier"))
  //     TestBezierPaths();
// 
  //   // reconstruction error on both meshes
  //   if (ImGui::Button("Distance Error"))
  //     ColorByReconstructionError(1);
  //   ImGui::SameLine();
  //   if (ImGui::Button("Normal Error"))
  //     ColorByReconstructionError(2);
  //   ImGui::SameLine();
  //   if (ImGui::Button("Reset Color"))
  //     ColorByReconstructionError(0);
  // }
  if (ImGui::CollapsingHeader("Test")) {
    // switches: recomputed every time they are turned on, hidden when off
    if (ImGui::Checkbox("Normals", &show_test_normals)) {
      PatchNormalEdges.clear();
      PatchNormals.clear();
      if (show_test_normals && has_paths) {
        bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
            MPatchDeco.PatchManager(), Features, PatchNormalEdges,
            PatchNormals);
        if (done)
          std::cout << "Esteem Normal Correct" << std::endl;
        else
          std::cout << "Esteem Normal Error" << std::endl;
      }
      show_test_normals &= has_paths;
    }
    ImGui::SameLine();
    if (ImGui::Checkbox("Bezier", &show_test_bezier)) {
      BezierPolylines.clear();
      BezierColorError.clear();
      if (show_test_bezier && has_paths)
        TestBezierPaths();
      show_test_bezier &= has_paths;
    }

    // mesh color: constant or reconstruction error (both meshes)
    int ColorMode = DrawColorMode;
    ImGui::RadioButton("Constant", &ColorMode, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Distance Error", &ColorMode, 1);
    ImGui::SameLine();
    ImGui::RadioButton("Normal Error", &ColorMode, 2);
    if (ColorMode != DrawColorMode)
      ColorByReconstructionError(ColorMode);
  }
  ImGui::End();
}

bool old_use_toon_shader = use_toon_shader;
int old_DrawMeshMode = DrawMeshMode;
bool old_showMesh = showMesh;
bool old_showFeatures = showFeatures;
bool old_showResult = showResult;
bool old_showCross = showCross;
bool old_showSing = showSing;
bool old_showBoundaries = showBoundaries;
bool old_ShowSymmPlane = ShowSymmPlane;
bool old_ShowOutline = ShowOutline;
bool old_drawOriginalBorders = drawOriginalBorders;

void SaveOldRenderMode() {
  old_use_toon_shader = use_toon_shader;
  old_DrawMeshMode = DrawMeshMode;
  old_showMesh = showMesh;
  old_showFeatures = showFeatures;
  old_showResult = showResult;
  old_showCross = showCross;
  old_showSing = showSing;
  old_showBoundaries = showBoundaries;
  old_ShowSymmPlane = ShowSymmPlane;
  old_ShowOutline = ShowOutline;
  old_drawOriginalBorders = drawOriginalBorders;
}

void RestoreOldRenderMode() {
  use_toon_shader = old_use_toon_shader;
  DrawMeshMode = old_DrawMeshMode;
  showMesh = old_showMesh;
  showFeatures = old_showFeatures;
  showResult = old_showResult;
  showCross = old_showCross;
  showSing = old_showSing;
  showBoundaries = old_showBoundaries;
  ShowSymmPlane = old_ShowSymmPlane;
  ShowOutline = old_ShowOutline;
  drawOriginalBorders = old_drawOriginalBorders;
}

void SetBaseScreenshotRenderMode() {
  use_toon_shader = false;
  DrawMeshMode = 0;
  showMesh = true;
  showFeatures = false;
  showResult = false;
  showCross = false;
  showSing = false;
  showBoundaries = false;
  ShowSymmPlane = false;
  ShowOutline = false;
  drawOriginalBorders = false;
}

void SetGLRotationUpperLeft() {
  // rotate current model viwo matrix of 45 degrrex alongy and 45 degrees along
  // x
  glRotatef(45, 0, 1, 0);
  glRotatef(45, 1, 0, 0);
}

void GLDrawMesh() {
  if (do_Batch_compute)
    rendercycles++;
  bool wire = true;
  bool smoothshade = true;
  switch (DrawMeshMode) {
  case 0:
    wire = false;
    smoothshade = true;
    break;
  case 1:
    wire = false;
    smoothshade = false;
    break;
  case 2:
    wire = true;
    smoothshade = true;
    break;
  default:
    wire = true;
    smoothshade = false;
    break;
  }

  // (screenshots removed)
  // if (get_screenshot_original) {
  //   SaveOldRenderMode();
  //   SetBaseScreenshotRenderMode();
  // }
  // if (get_screenshot_boundaries) {
  //   SaveOldRenderMode();
  //   SetBaseScreenshotRenderMode();
  //   showBoundaries = true;
  // }

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w,
               clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(40, (GLdouble)display_w / (GLdouble)display_h, 0.1, 100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 3.5f, 0, 0, 0, 0, 1, 0);

  trackball.GetView();

  trackball.Apply();

  // if (get_screenshot_original || get_screenshot_boundaries)
  //   SetGLRotationUpperLeft();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  bool DrawSurface = true;

  if (DrawChangePos) {
    GLDraw::DrawPoints<ScalarType>(CurvChangePos,
                                   Geo::Point3<ScalarType>(1, 0, 1), 40);
    // CurvChangeDir
    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        CurvChangeSegments;
    for (size_t i = 0; i < CurvChangePos.size(); i++) {
      Geo::Point3<ScalarType> P0 = CurvChangePos[i];
      Geo::Point3<ScalarType> P1 = CurvChangePos[i] + CurvChangeDir[i] * AvEdge;
      CurvChangeSegments.push_back(std::make_pair(P0, P1));
    }
    GLDraw::GLDrawSegments<ScalarType>(CurvChangeSegments, 20, 20,
                                       Geo::Point3<ScalarType>(0, 0, 0),
                                       Geo::Point3<ScalarType>(0, 0, 0));
  }

  // GLDraw::GLDrawPolylines(LoopPolylines, LoopPolylinesCurvColors,30);
  //  (const std::vector<Geo::Point3<ScalarType>> &Pos,
  //                 const std::vector<std::pair<int, int>> &Edges,
  //                 const ScalarType LineWidth,
  //                 const std::vector<Geo::Point3<ScalarType>> &colorEdges)

  //   std::vector<Geo::Point3<ScalarType>> CurvEdgeMeshPos;
  // std::vector<std::pair<int, int>> CurvEdges;

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  if (glGetError() != GL_NO_ERROR) {
    std::cout << "OpenGL Error before drawing edges: " << glGetError()
              << std::endl;
    exit(0);
  }

  if (DrawCurv)
    GLDraw::GLDrawEdges<ScalarType>(CurvEdgeMeshPos, CurvEdges, (ScalarType)30,
                                    CurvEdgeColors);
  if (DrawColChange)
    GLDraw::GLDrawEdges<ScalarType>(CurvEdgeMeshPos, CurvEdges, (ScalarType)30,
                                    ClusterEdgeColors);

  glLineWidth(line_thick);

  GLDraw::glColor(Geo::Point3<ScalarType>(1, 0.9, 0.55));

  // DRAW THE BORDERS
  int sizeBordInput = 5;
  int sizeBoundariesInput = 15;
  int sizeBoundariesVertsInput = 30;
  int sizeBordDef = sizeBordInput;
  int sizeBoundariesDef = sizeBoundariesInput;
  int sizeFiexedVDef = sizeBoundariesVertsInput;

  // DRAW PATCHES
  if (showBoundaries) {
    for (size_t i = 0; i < Boundary.size(); i++)
      GLDraw::GLDrawEdges<ScalarType>(VertPos, Boundary[i], sizeBoundariesInput,
                                      Geo::Point3<ScalarType>(0, 0, 0));
  }

  if (showResult) {
    // GLDraw::glColor(solved_mesh_color);
    assert(FaceColorReconstructed.size() == SolvedConnectivity.size());
    GLDraw::DrawSurfaceMesh<ScalarType>(
        SolvedVertPos, SolvedConnectivity, SolvedFaceNormals, SolvedVertNormals,
        DrawSurface, wire, smoothshade, false, GLDraw::TRTriangle,
        &FaceColorReconstructed);
  }

  if (showMesh) {

    if (illustrative == nullptr)
      illustrative = new MechanicalIllustrativeShader();

    if (use_toon_shader) {
      illustrative->use();
      illustrative->setBaseColor(0.96f, 0.97f, 1.00f);
      illustrative->setLineColor(0.08f, 0.08f, 0.08f);
      illustrative->setLightDirView(0.1f, 0.9f, 1.2f);
      illustrative->setRim(0.35f, 2.0f);
      illustrative->setEdgeWidth(0.20f);
    }

    glDisable(GL_LIGHTING);
    GLDraw::glColor(mesh_color);
    if ((showResult) && (has_result)) {
      glColor4f(.5f, 1.f, .8f, .3f);

      bool useBlend = showResult;
      GLDraw::DrawSurfaceMesh<ScalarType>(
          VertPos, Connectivity, FaceNormals, VertNormals, DrawSurface, wire,
          smoothshade, true, GLDraw::TRTriangle);
    } else {
      assert(FaceColorTarget.size() == Connectivity.size());
      GLDraw::DrawSurfaceMesh<ScalarType>(
          VertPos, Connectivity, FaceNormals, VertNormals, DrawSurface, wire,
          smoothshade, false, GLDraw::TRTriangle, &FaceColorTarget, NULL,
          !use_toon_shader);
    }
    if (use_toon_shader)
      glUseProgram(0);

    if (drawOriginalBorders)
      GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF, 20);
  }

  // DRAW THE FEATURES
  if (showFeatures) {
    GLDraw::GLDrawEdges<ScalarType>(VertPos, Features, sizeBoundariesInput,
                                    Geo::Point3<ScalarType>(0, 0, 0));
  }

  if (ShowSymmPlane)
    GLDraw::glDrawPlane<ScalarType>(SymmetryPlane, MeshBox.Diag() / 2);

  glPopAttrib();

  // draw the bezier paths
  for (size_t i = 0; i < BezierPolylines.size(); i++) {
    GLDraw::GLDrawPolyline<ScalarType>(BezierPolylines[i].PolyPos,
                                       BezierColorError[i], 20);
  }

  for (size_t i = 0; i < PatchNormalEdges.size(); i++)
    for (size_t j = 0; j < PatchNormalEdges[i].size(); j++)
      for (size_t k = 0; k < PatchNormalEdges[i][j].size(); k++) {
        int IndexV0 = PatchNormalEdges[i][j][k].first;
        int IndexV1 = PatchNormalEdges[i][j][k].second;
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        Geo::Point3<ScalarType> AvgSegment = (P0 + P1) * (ScalarType)0.5;
        Geo::Point3<ScalarType> CurrNormal = PatchNormals[i][j][k];
        Geo::Point3<ScalarType> P2 = AvgSegment + CurrNormal * AvEdge * 2;
        GLDraw::GLDrawSegment<ScalarType>(AvgSegment, P2, 10, 0,
                                          Geo::Point3<ScalarType>(0, 1, 1));
      }

  if (has_cross_field) {
    ScalarType scaleVal = AvEdge * (ScalarType)0.5;
    if (showCross)
      GLDraw::DrawCrossFields<ScalarType>(FaceCurv, scaleVal, maxQCrossVert,
                                          minQCrossVert);

    if (showSing)
      GLDraw::DrawSingularity<ScalarType>(VertPos, SingIndex, SingValue);
  }

  if (showMesh && ShowOutline) {
    GLDraw::DrawOutline<ScalarType>(VertPos, Connectivity, line_thick);
  }
  if (showResult) {
    GLDraw::DrawOutline<ScalarType>(SolvedVertPos, SolvedConnectivity,
                                    line_thick);
  }
//   if (get_screenshot_boundaries) {
//     glFlush();
//     glFinish();
//     GetStreenShotMesh("boundaries");
//     glFinish();
//     get_screenshot_boundaries = false;
//     saved_last_screenshot = true;
//     RestoreOldRenderMode();
//     if (do_Batch_compute) {
//       do_Batch_compute = false;
//       ImGui_ImplOpenGL2_Shutdown();
//       ImGui_ImplGlfw_Shutdown();
//       ImGui::DestroyContext();
//       glfwDestroyWindow(window);
//       glfwTerminate();
//       glFlush();
//       glFinish();
//       exit(0);
//       // glfwSetWindowShouldClose(window, GLFW_TRUE);
//     }
//   }
//   if (get_screenshot_original) {
//     glFlush();
//     glFinish();
//     GetStreenShotMesh("original");
//     glFinish();
//     get_screenshot_original = false;
//     get_screenshot_boundaries = true;
//     RestoreOldRenderMode();
//   }
//   if (do_Batch_compute && (!has_result)) {
//     ProcessAll();
//     // do_Batch_compute = false;
//   }
  // batch (-batch): process, save everything and quit
  if (do_Batch_compute) {
    do_Batch_compute = false;
    ProcessAll();
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    exit(0);
  }
}

// the default parameters that can be changed from the settings file
// (Settings is declared above)
void RegisterSettings(LoopWeaverSettings &S) {
  S.Section("decomposition");
  S.Add("use_symmetry", UseSymmetry, "cut with the symmetry plane x=0");
  S.Add("sharp_angle", feature_angle, "sharp features angle (degrees)");
  S.Add("feature_erode_steps", ero_steps, "sharp features cleaning");
  S.Add("feature_dilate_steps", dil_steps, "sharp features cleaning");
  S.Add("field_smoothness", GlobalSmoothVal, "cross field smoothness");
  S.Add("clean_mesh", do_clean, "collapse bad faces and solve folds");
  S.Add("sampling_density", numSamples, "number of tracing samples");
  S.Add("dynamic_update", dynamic_updates, "retrace the candidates");
  S.Add("split_when_final_remove", MPatchDeco.split_removal, "");
  S.Add("complete_all_features", MPatchDeco.force_all_features, "");
  S.Add("relax_paths", MPatchDeco.relax_paths, "spread the paths at the end");
  S.Add("relax_max_step", MPatchDeco.RelaxParam.max_step,
        "max seed shift of the relaxation (graph steps)");
  S.Add("relax_max_rounds", MPatchDeco.RelaxParam.max_rounds, "");
  S.Add("relax_remove_after", MPatchDeco.RelaxParam.remove_after,
        "try to remove paths after the relaxation, then relax again");
  // (curvature priority removed from the settings)
  // S.Add("curvature_priority", MPatchDeco.curvature_priority,
        // "candidate priority x (1 + weight * curvature anisotropy)");
  // S.Add("curvature_priority_weight", MPatchDeco.curvature_priority_weight,
        // "multiplier in [1, 1 + weight]");
  // S.Add("curvature_priority_min_dir", MPatchDeco.curvature_priority_min_dir,
        // "1: nodes along the min curvature, 0: along the max one");
  // S.Add("curvature_priority_both_dirs", MPatchDeco.curvature_priority_both_dirs,
        // "1: nodes along both curvature directions (min_dir ignored)");
  // (curvature seeds removed: not used, see the commented lines)
  // S.Add("curvature_seeds", MPatchDeco.curvature_seeds,
        // "first samples on the most anisotropic curvature");
  // S.Add("curvature_seeds_ratio", MPatchDeco.curvature_seeds_ratio,
        // "fraction of the vertices used");
  // S.Add("curvature_seeds_max_dir", MPatchDeco.curvature_seeds_max_dir,
        // "1: along the max curvature, 0: along the min one");
  S.Add("progress_every", MPatchDeco.progress_every,
        "progress message every N added/removed paths (-1 none)");

  S.Section("conditions");
  S.Add("min_sides", MPatchDeco.MinSides, "");
  S.Add("max_sides", MPatchDeco.MaxSides, "");
  S.Add("cc_ability", MPatchDeco.CCAbility, "-1 not used");
  S.Add("single_patch_sing", MPatchDeco.single_sing_cond, "");
  S.Add("match_sing_values", MPatchDeco.match_sing_cond, "");
  S.Add("loop_reconstruction_cond", loop_recon_cond, "");
  S.Add("loop_reconstruction_local", loop_recon_local,
        "local (per patch) version");
  S.Add("reconstruction_error", maxErrRatio, "max distance (x bbox diagonal)");
  S.Add("use_normal_error", use_normal_error, "");
  S.Add("normal_angle_error", maxNormAngle, "degrees");
  S.Add("normal_angle_error_percent", maxAnglePercentile,
        "fraction of faces allowed over the angle");
  S.Add("cond_thin_plate_only", cond_biharmonic_only, "");
  S.Add("cond_use_original_tessellation",
        LoopCond.CurveSolv.use_original_meshing, "");
  S.Add("cond_remesh_edge_factor", LoopCond.CurveSolv.remesh_facctor,
        "-1 no remesh");
  S.Add("cond_resample_paths", LoopCond.CurveSolv.resample_paths, "");
  S.Add("cond_subsample_factor", LoopCond.CurveSolv.subsample_factor, "");
  S.Add("side_curvature_cond", side_curvature_cond, "");
  S.Add("max_side_angle", MaxSideSumAngle, "degrees");
  S.Add("side_smooth_steps", SideSmoothSteps,
        "Taubin steps on the sides before measuring the angle (0: none)");
  S.Add("gaussian_curvature_cond", gaussian_curvature_cond,
        "developability of the patches");
  S.Add("max_gaussian_curvature", MaxGaussianCurvature,
        "total absolute Gaussian curvature of a patch (degrees)");
  // S.Add("gaussian_taubin_steps", GaussianTaubinSteps,
  //       "Taubin steps on the patch before measuring");
  S.Add("gaussian_defect_smooth", GaussianDefectSmooth,
        "diffusion steps of the angle defect before measuring");
  S.Add("bezier_error_cond", beier_error_cond, "");
  S.Add("max_bezier_error", MaxBezierErrorPerc, "percent of bbox diagonal");
  S.Add("normal_approximation_cond", normal_approx_cond, "");
  S.Add("normal_approx_error", NormalCondError, "degrees");
  S.Add("normal_approx_percent", NormalCondPercent, "");

  S.Section("final reconstruction");
  S.Add("final_one_patch_at_a_time", final_per_patch, "0: global");
  S.Add("final_thin_plate_only", final_biharmonic_only, "");
  S.Add("final_iterations", FinalSurfacingIterations, "");
  S.Add("final_use_original_tessellation", CurveSolv.use_original_meshing,
        "");
  S.Add("final_remesh_edge_factor", CurveSolv.remesh_facctor,
        "-1 no remesh");
  S.Add("final_resample_paths", CurveSolv.resample_paths, "");
  S.Add("final_subsample_factor", CurveSolv.subsample_factor, "");
}

// prints the command line options
void PrintUsage() {
  std::cout << "usage: loop_weaver mesh.obj [mesh.field] [options]\n"
            << " -batch                 process, save everything and quit\n"
            << " -settings <file>       default settings file\n"
            << " -symm / -no_symm       symmetry (plane x=0)\n"
            << " -split / -no_split     split when final remove\n"
            << " -complete_features / -no_complete_features\n"
            << " -dynamic / -no_dynamic dynamic update\n"
            << " -angle <deg>           sharp features angle (180: none)\n"
            << " -subsampl <n>          subsampling of the paths (condition)\n"
            << " -thinplate / -no_thinplate  thin-plate only surfaces\n"
            << "                        (condition and final reconstruction)\n"
            << " -error <ratio>         reconstruction error (x bbox diagonal)\n"
            << " -no_error              no loop reconstruction condition\n"
            << " -normal_error <deg>    normal error angle\n"
            << " -no_normal             no normal error\n"
            << " -side_angle <deg>      max side angle\n"
            << " -no_side               no side curvature condition\n"
            << " -side_smooth <n>       smoothing steps on the sides before\n"
            << "                        measuring the angle (0: none)\n"
            << " -gauss <deg>           max total Gaussian curvature of a\n"
            << "                        patch (developability)\n"
            << " -no_gauss              no Gaussian curvature condition\n"
            // << " -gauss_taubin <n>      Taubin steps before measuring\n"
            << " -gauss_smooth <n>      defect diffusion steps (default 5)\n"
            << " -bezier <percent>      max bezier error\n"
            << " -no_bezier             no bezier error condition\n"
            << " -relax / -no_relax     path relaxation at the end\n"
            // << " -curv_seeds / -curv_seeds_min  first samples on the most\n"
            // << "                        anisotropic curvature (max / min dir)\n"
            // << " -curv_priority / -no_curv_priority  candidates on the\n"
            // << "                        anisotropic curvature first\n"
            // << " -curv_priority_max / -curv_priority_both  along the max\n"
            // << "                        curvature or both directions\n"
            << " -error_setup <0..3>    preset of the reconstruction errors\n"
            << "                        (-1: no reconstruction condition)\n"
            << "Without options the defaults (settings file) are used."
            << std::endl;
}

// the command line options (they override the settings file); returns true
// if the mesh has been loaded
bool ParseCommandLine(int argc, char *argv[]) {
  bool has_load_mesh = false;
  for (int i = 1; i < argc; i++) {
    std::string Arg(argv[i]);
    std::string Ext = GetFileExtension(Arg);
    bool has_value = (i + 1 < argc);

    // mesh and field
    if (Ext == std::string("obj")) {
      PathMesh = Arg;
      has_load_mesh = LoadMesh(PathMesh);
      if (!has_load_mesh)
        exit(0);
      continue;
    }
    if ((Ext == std::string("ffield")) || (Ext == std::string("field"))) {
      if (!has_load_mesh) {
        std::cout << "PASS MESH BEFORE FIELD " << std::endl;
        exit(0);
      }
      PathField = Arg;
      has_cross_field = LoadField(PathField);
      if (!has_cross_field)
        exit(0);
      has_loaded_field = true;
      FaceCurvLoaded = FaceCurv;
      VertCurvLoaded = VertCurv;
      continue;
    }

    if (Arg == "-batch") {
      do_Batch_compute = true;
      continue;
    }
    // already read before the settings
    if ((Arg == "-settings") && has_value) {
      i++;
      continue;
    }

    // decomposition
    if (Arg == "-symm") {
      UseSymmetry = true;
      continue;
    }
    if (Arg == "-no_symm") {
      UseSymmetry = false;
      continue;
    }
    if (Arg == "-split") {
      MPatchDeco.split_removal = true;
      continue;
    }
    if (Arg == "-no_split") {
      MPatchDeco.split_removal = false;
      continue;
    }
    if (Arg == "-complete_features") {
      MPatchDeco.force_all_features = true;
      continue;
    }
    if (Arg == "-no_complete_features") {
      MPatchDeco.force_all_features = false;
      continue;
    }
    if (Arg == "-dynamic") {
      dynamic_updates = true;
      continue;
    }
    if ((Arg == "-no_dynamic") || (Arg == "-nodynamic")) {
      dynamic_updates = false;
      continue;
    }
    if (Arg == "-relax") {
      MPatchDeco.relax_paths = true;
      continue;
    }
    if ((Arg == "-no_relax") || (Arg == "-norelax")) {
      MPatchDeco.relax_paths = false;
      continue;
    }

    // // candidate priority multiplied by the curvature anisotropy (loops
    // // around the curved tubes first)
    // if (Arg == "-curv_priority") {
      // MPatchDeco.curvature_priority = true;
      // continue;
    // }
    // // multiplier on the nodes along the max curvature instead of the min one
    // if (Arg == "-curv_priority_max") {
      // MPatchDeco.curvature_priority = true;
      // MPatchDeco.curvature_priority_min_dir = false;
      // continue;
    // }
    // // multiplier on the nodes of both curvature directions
    // if (Arg == "-curv_priority_both") {
      // MPatchDeco.curvature_priority = true;
      // MPatchDeco.curvature_priority_both_dirs = true;
      // continue;
    // }
    // if (Arg == "-no_curv_priority") {
      // MPatchDeco.curvature_priority = false;
      // continue;
    // }
    // // first sampling from the most anisotropic curvature vertices, along the
    // // max (loops around tubes) or the min curvature direction
    // if (Arg == "-curv_seeds") {
      // MPatchDeco.curvature_seeds = true;
      // MPatchDeco.curvature_seeds_max_dir = true;
      // continue;
    // }
    // if (Arg == "-curv_seeds_min") {
      // MPatchDeco.curvature_seeds = true;
      // MPatchDeco.curvature_seeds_max_dir = false;
      // continue;
    // }
    // sharp features angle (180: no features)
    if ((Arg == "-angle") && has_value) {
      feature_angle = atof(argv[++i]);
      olfeature_angle = feature_angle;
      continue;
    }

    // surfaces
    // subsampling of the paths in the condition (as before: also the remesh
    // factor, made coherent by the solver)
    if ((Arg == "-subsampl") && has_value) {
      LoopCond.CurveSolv.remesh_facctor = atof(argv[i + 1]);
      LoopCond.CurveSolv.subsample_factor = atoi(argv[i + 1]);
      LoopCond.CurveSolv.MakeParametersCoherent();
      i++;
      continue;
    }
    if (Arg == "-thinplate") {
      cond_biharmonic_only = true;
      final_biharmonic_only = true;
      continue;
    }
    if (Arg == "-no_thinplate") {
      cond_biharmonic_only = false;
      final_biharmonic_only = false;
      continue;
    }

    // conditions
    if ((Arg == "-error") && has_value) {
      maxErrRatio = atof(argv[++i]);
      loop_recon_cond = true;
      continue;
    }
    if (Arg == "-no_error") {
      loop_recon_cond = false;
      continue;
    }
    if ((Arg == "-normal_error") && has_value) {
      maxNormAngle = atof(argv[++i]);
      use_normal_error = true;
      continue;
    }
    if (Arg == "-no_normal") {
      use_normal_error = false;
      continue;
    }
    if ((Arg == "-side_angle") && has_value) {
      MaxSideSumAngle = atof(argv[++i]);
      side_curvature_cond = true;
      continue;
    }
    if (Arg == "-no_side") {
      side_curvature_cond = false;
      continue;
    }
    if ((Arg == "-side_smooth") && has_value) {
      SideSmoothSteps = atoi(argv[++i]);
      continue;
    }
    // developability: max total absolute Gaussian curvature of a patch
    if ((Arg == "-gauss") && has_value) {
      MaxGaussianCurvature = atof(argv[++i]);
      gaussian_curvature_cond = true;
      continue;
    }
    if (Arg == "-no_gauss") {
      gaussian_curvature_cond = false;
      continue;
    }
    // if ((Arg == "-gauss_taubin") && has_value) {
    //   GaussianTaubinSteps = atoi(argv[++i]);
    //   continue;
    // }
    if ((Arg == "-gauss_smooth") && has_value) {
      GaussianDefectSmooth = atoi(argv[++i]);
      continue;
    }
    if ((Arg == "-bezier") && has_value) {
      MaxBezierErrorPerc = atof(argv[++i]);
      beier_error_cond = true;
      continue;
    }
    if (Arg == "-no_bezier") {
      beier_error_cond = false;
      continue;
    }

    // preset of the reconstruction errors (as before)
    if ((Arg == "-error_setup") && has_value) {
      int SetupType = atoi(argv[++i]);
      if (SetupType < 0) {
        loop_recon_cond = false;
        continue;
      }
      if (SetupType > 3) {
        std::cout << "Invalid Setup Type for -error_setup" << std::endl;
        exit(0);
      }
      loop_recon_cond = true;
      switch (SetupType) {
      case 0:
        maxErrRatio = 0.02;
        maxNormAngle = 30;
        maxAnglePercentile = 0.2;
        break;
      case 1:
        maxErrRatio = 0.03;
        maxNormAngle = 40;
        maxAnglePercentile = 0.3;
        break;
      case 2:
        maxErrRatio = 0.04;
        maxNormAngle = 50;
        maxAnglePercentile = 0.5;
        break;
      case 3:
        maxErrRatio = 0.05;
        maxNormAngle = 60;
        maxAnglePercentile = 0.6;
        break;
      }
      continue;
    }

    std::cout << "Unknown option " << Arg << " (ignored)" << std::endl;
  }
  return has_load_mesh;
}

int main(int argc, char *argv[]) {
  InitDefaultParam();

  // default settings from the text file next to the executable (written
  // with the compiled defaults if missing); -settings <file> to use another
  // one. The command line options below override them.
  {
    std::string Exe(argv[0]);
    size_t Slash = Exe.find_last_of("/\\");
    std::string SettingsPath =
        ((Slash == std::string::npos) ? std::string("")
                                      : Exe.substr(0, Slash + 1)) +
        std::string("loopweaver_settings.txt");
    for (int i = 1; i + 1 < argc; i++)
      if (std::string(argv[i]) == std::string("-settings"))
        SettingsPath = std::string(argv[i + 1]);
    RegisterSettings(Settings);
    Settings.LoadOrCreate(SettingsPath);
    // values derived from the settings
    SParam.GlobalSmoothVal = GlobalSmoothVal;
    olfeature_angle = feature_angle;
  }

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    PrintUsage();
    exit(0);
  }

  // OLD command line parsing (replaced by ParseCommandLine)
//   bool has_load_mesh = false;
//   for (int i = 1; i < argc; i++) {

//     if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
//       PathMesh = std::string(argv[i]);
//       has_load_mesh = LoadMesh(PathMesh);
//       if (!has_load_mesh)
//         exit(0);
//       continue;
//     }

//     if ((GetFileExtension(std::string(argv[i])) == std::string("ffield")) ||
//         (GetFileExtension(std::string(argv[i])) == std::string("field"))) {
//       if (!has_load_mesh) {
//         std::cout << "PASS MESH BEFORE FIELD " << std::endl;
//         exit(0);
//       }
//       PathField = std::string(argv[i]);
//       has_cross_field = LoadField(PathField);

//       if (!has_cross_field)
//         exit(0);
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-batch")) {
//       do_Batch_compute = true;
//     }
//     if (std::string(argv[i]) == std::string("-localcond")) {
//       loop_recon_local = true;
//       continue;
//     }
//     // verbosity of the surfacing messages: 0 silent, 1 steps, 2 per patch /
//     // per iteration, 3 also per iteration inside the local condition
//     // iterations of the final extraction
//     if (std::string(argv[i]) == std::string("-finalits")) {
//       FinalSurfacingIterations = atoi(argv[i + 1]);
//       i++;
//       continue;
//     }
//     // // final extraction from the original (not smoothed) patches
//     // if (std::string(argv[i]) == std::string("-finalnosmooth")) {
//     //   final_no_smooth = true;
//     //   continue;
//     // }
//     // initial surface of the patches: the thin-plate fill tangent to the
//     // border normals is the default, -harmonic goes back to the implicit
//     // Laplacian (flat) one
//     if (std::string(argv[i]) == std::string("-harmonic")) {
//       LoopCond.CurveSolv.biharmonic_initial = false;
//       CurveSolv.biharmonic_initial = false;
//       continue;
//     }
//     // surface = thin-plate fill only, in the condition / in the final extraction
//     if (std::string(argv[i]) == std::string("-thinplate_cond")) {
//       cond_biharmonic_only = true;
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-thinplate_final")) {
//       final_biharmonic_only = true;
//       continue;
//     }
//     // tension of the thin-plate fill (0 = pure thin plate)
//     if (std::string(argv[i]) == std::string("-tension")) {
//       LoopCond.CurveSolv.biharmonic_tension = atof(argv[i + 1]);
//       CurveSolv.biharmonic_tension = LoopCond.CurveSolv.biharmonic_tension;
//       i++;
//       continue;
//     }
//     // overshoot check of the thin-plate fill (fraction of the patch diagonal, <= 0 off)
//     if (std::string(argv[i]) == std::string("-maxbulge")) {
//       LoopCond.CurveSolv.biharmonic_max_bulge = atof(argv[i + 1]);
//       CurveSolv.biharmonic_max_bulge = LoopCond.CurveSolv.biharmonic_max_bulge;
//       i++;
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-progress")) {
//       LoopCond.CurveSolv.progress_level = atoi(argv[i + 1]);
//       CurveSolv.progress_level = LoopCond.CurveSolv.progress_level;
//       i++;
//       continue;
//     }

//     if (std::string(argv[i]) == std::string("-symm")) {
//       UseSymmetry = true;
//     }

//     if (std::string(argv[i]) == std::string("-nonormalerror")) {
//       use_normal_error = false;
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-noclean")) {
//       do_clean = false;
//       continue;
//     }

//     if (std::string(argv[i]) == std::string("-split")) {
//       MPatchDeco.split_removal = true;
//     }

//     if (std::string(argv[i]) == std::string("-angle")) {
//       feature_angle = atoi(argv[i + 1]);
//       olfeature_angle = feature_angle;
//       i++;
//       continue;
//     }

//     // if (std::string(argv[i]) == std::string("-curvchange")) {
//     //   UseCurvChangeNodes = true;
//     //   continue;
//     // }

//     if (std::string(argv[i]) == std::string("-create_mesh")) {
//       LoopCond.CurveSolv.use_original_meshing = false;
//       CurveSolv.use_original_meshing = false;
//       // default with the generated tessellation
//       LoopCond.CurveSolv.resample_paths = true;
//       CurveSolv.resample_paths = true;
//       continue;
//     }

//     if (std::string(argv[i]) == std::string("-subsampl")) {
//       LoopCond.CurveSolv.remesh_facctor = atof(argv[i + 1]);
//       LoopCond.CurveSolv.subsample_factor = atoi(argv[i + 1]);
//       LoopCond.CurveSolv.MakeParametersCoherent();
//       i++;
//       continue;
//     }

//     if (std::string(argv[i]) == std::string("-minsides")) {
//       MPatchDeco.MinSides = atoi(argv[i + 1]);
//       MPatchDeco.MinSides = std::min(3, MPatchDeco.MinSides);
//       MPatchDeco.MinSides = std::max(2, MPatchDeco.MinSides);
//       i++;
//       continue;
//     }

//     if (std::string(argv[i]) == std::string("-curvangle")) {
//       MaxSideSumAngle = atof(argv[i + 1]);
//       MaxSideSumAngle = std::min((ScalarType)180, MaxSideSumAngle);
//       i++;
//       if (MaxSideSumAngle <= 0) {
//         side_curvature_cond = false;
//       } else {
//         SideAngleCond.MaxSideAngle = MaxSideSumAngle;
//         side_curvature_cond = true;
//       }
//       continue;
//     }

//     if (std::string(argv[i]) == std::string("-bezier")) {
//       MaxBezierErrorPerc = atof(argv[i + 1]);
//       i++;
//       MaxBezierErrorPerc = std::min(MaxBezierErrorPerc, (ScalarType)5);
//       if (MaxBezierErrorPerc <= 0) {
//         beier_error_cond = false;
//       } else {
//         beier_error_cond = true;
//       }
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-ccability")) {
//       MPatchDeco.CCAbility = atof(argv[i + 1]);
//       i++;
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-matchsing")) {
//       MPatchDeco.match_sing_cond = true;
//       MPatchDeco.single_sing_cond = true;
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-nodynamic")) {
//       dynamic_updates = false;
//       continue;
//     }
//     // no smoothing of the mesh around each added path (as the quadrangulator)
//     if (std::string(argv[i]) == std::string("-nodynsmooth")) {
//       MPatchDeco.dynamicSmoothing = false;
//       continue;
//     }
//     // decomposition progress message every N added/removed paths (-1 none)
//     if (std::string(argv[i]) == std::string("-decoprogress")) {
//       MPatchDeco.progress_every = atoi(argv[i + 1]);
//       continue;
//     }
//     // at the end move the paths to spread them more uniformly
//     if (std::string(argv[i]) == std::string("-relax")) {
//       MPatchDeco.relax_paths = true;
//       continue;
//     }
//     // graph arc costs weighted by the curvature anisotropy, with its weight
//     if (std::string(argv[i]) == std::string("-curvgraph")) {
//       curvature_graph = true;
//       curvature_graph_weight = atof(argv[i + 1]);
//       continue;
//     }
//     if (std::string(argv[i]) == std::string("-norelax")) {
//       MPatchDeco.relax_paths = false;
//       continue;
//     }
//     // maximum shift of the seeds of the path relaxation (graph steps)
//     if (std::string(argv[i]) == std::string("-relaxstep")) {
//       MPatchDeco.RelaxParam.max_step = atoi(argv[i + 1]);
//       continue;
//     }
//     // number of tracing samples
//     if (std::string(argv[i]) == std::string("-samples")) {
//       numSamples = atoi(argv[i + 1]);
//       continue;
//     }
    
//     if (std::string(argv[i]) == std::string("-error_setup")) {
//       int SetupType = atoi(argv[i + 1]);

//       if (SetupType < 0) {
//         loop_recon_cond = false;
//         continue;
//       }
//       if ((SetupType < 0) || (SetupType > 3)) {
//         std::cout << "Invalid Setup Type for -error_setup" << std::endl;
//         exit(0);
//       }

//       loop_recon_cond = true;

//       switch (SetupType) {
//       case 0:
//         maxErrRatio = 0.02;
//         maxNormAngle = 30;
//         maxAnglePercentile = 0.2;
//         break;
//       case 1:
//         maxErrRatio = 0.03;
//         maxNormAngle = 40;
//         maxAnglePercentile = 0.3;
//         break;
//       case 2:
//         maxErrRatio = 0.04;
//         maxNormAngle = 50;
//         maxAnglePercentile = 0.5;
//         break;
//       case 3:
//         maxErrRatio = 0.05;
//         maxNormAngle = 60;
//         maxAnglePercentile = 0.6;
//         break;
//       }
//     }
//   }

  bool has_load_mesh = ParseCommandLine(argc, argv);

  if (!has_load_mesh) {
    std::cout << "ERROR MESH NOT LOADED" << std::endl;
    exit(0);
  }

  InitGLFW_Window();

  InitIMGui();

  trackball.center = Geo::Point3<float>(0, 0, 0);
  trackball.radius = 1;

  trackballUV.center = Geo::Point3<float>(0, 0, 0);
  trackballUV.radius = 1;

  glewInit();

  // Main loop
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    SetRenderBar();
    SetToolBar();
    SetConditionsBar();
    // Rendering
    ImGui::Render();

    GLDrawMesh();

    // if ((saved_last_screenshot)&&(do_Batch_compute)&&(rendercycles>10)) {
    //   do_Batch_compute=false;
    //   exit(0);
    // }

    // If you are using this code with non-legacy OpenGL header/contexts
    // (which you should not, prefer using imgui_impl_opengl3.cpp!!), you may
    // need to backup/reset/restore other state, e.g. for current shader using
    // the commented lines below.
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  // delete(illustrative);
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}