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
 */

#ifndef SKETCH_LOG
#define SKETCH_LOG

// The log of the sketch actions of an interactive session (quadrangulate,
// draw a line, remove a sketch, remove all), with the settings they used, to
// replay a session from the loaded mesh (e.g. up to a crash). Every action is
// written and flushed before it runs, so after a crash the file ends with the
// action that crashed.
// One action per line:
//   Q <settings>                 quadrangulate (not used by the GUI any more)
//   L <settings> N x y z ...     a drawn line of N points (as drawn)
//   R <settings> Index           remove the sketch Index
//   A <settings>                 remove all the sketches
// settings: view requad_all region_mode refit patch_rings

#include "point3.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType> struct SketchAction {
  char Type = 'Q';
  int View = 0, RequadAll = 0, RegionMode = 0, Refit = 0, PatchRings = 0;
  int Index = -1;
  std::vector<Geo::Point3<ScalarType>> Line;
  // snap to the constraints, in quad edges (0 = off); last in the row, an
  // old log without it reads 0
  float Snap = 0;
};

template <class ScalarType> class SketchLog {
  std::string Path;

public:
  // a new log; the log of the previous session (e.g. the one that crashed)
  // is kept as PrevPath(), the file emptied
  void Start(const std::string &_Path) {
    Path = _Path;
    // (only a log with drawn lines: one with only quadrangulations or
    // removals would overwrite the previous one, e.g. the one of a crash,
    // with nothing to replay)
    if (HasSketches(Path))
      std::rename(Path.c_str(), PrevPath().c_str());
    std::ofstream F(Path.c_str(), std::ios::trunc);
  }

  // the log of the previous session
  std::string PrevPath() const { return Path + ".prev"; }

  // the actions to replay: the ones of this session, or if no sketch yet the
  // ones of the previous session (the program started again after a crash)
  bool LoadLast(std::vector<SketchAction<ScalarType>> &Actions) const {
    const std::string &From = HasSketches(Path) ? Path : PrevPath();
    return Load(From, Actions) && !Actions.empty();
  }

  // a log with at least one drawn line (only quadrangulations or removals,
  // e.g. of the sketches of a replay, have nothing to replay)
  static bool HasSketches(const std::string &FilePath) {
    std::vector<SketchAction<ScalarType>> Actions;
    Load(FilePath, Actions);
    for (size_t i = 0; i < Actions.size(); i++)
      if (Actions[i].Type == 'L')
        return true;
    return false;
  }

  // appended and flushed at once
  void Add(const SketchAction<ScalarType> &A) const {
    if (Path.empty())
      return;
    std::ofstream F(Path.c_str(), std::ios::app);
    F.precision(17);
    F << A.Type << " " << A.View << " " << A.RequadAll << " " << A.RegionMode
      << " " << A.Refit << " " << A.PatchRings;
    if (A.Type == 'L') {
      F << " " << A.Line.size();
      for (size_t i = 0; i < A.Line.size(); i++)
        F << " " << A.Line[i].X << " " << A.Line[i].Y << " " << A.Line[i].Z;
    }
    if (A.Type == 'R')
      F << " " << A.Index;
    F << " " << A.Snap;
    F << std::endl;
  }

  static bool Load(const std::string &FilePath,
                   std::vector<SketchAction<ScalarType>> &Actions) {
    Actions.clear();
    std::ifstream F(FilePath.c_str());
    if (!F.is_open())
      return false;
    std::string Row;
    while (std::getline(F, Row)) {
      std::istringstream In(Row);
      SketchAction<ScalarType> A;
      if (!(In >> A.Type >> A.View >> A.RequadAll >> A.RegionMode >> A.Refit >>
            A.PatchRings))
        continue;
      if (A.Type == 'L') {
        size_t N = 0;
        In >> N;
        for (size_t i = 0; i < N; i++) {
          Geo::Point3<ScalarType> P;
          In >> P.X >> P.Y >> P.Z;
          A.Line.push_back(P);
        }
      }
      if (A.Type == 'R')
        In >> A.Index;
      if (!(In >> A.Snap))
        A.Snap = 0;
      Actions.push_back(A);
    }
    return true;
  }
};

} // namespace Geo

#endif
