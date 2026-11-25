#ifndef QUAD_REGION_HISTORY
#define QUAD_REGION_HISTORY

// Remembers, for every face of a quad mesh edited by region remeshes, which
// region remesh it comes from (an id, -1 = the full quadrangulation). The
// owner (e.g. the GUI) keeps one instance per quad mesh and passes it to the
// tools that change the mesh (QuadRegionRemesher::Remesh,
// QuadLayoutRefit::Refit), which keep it up to date. Then the faces of a
// region can be found again later, e.g. to remesh exactly the area of a
// sketch being removed.

#include <map>
#include <set>
#include <vector>

namespace Geo {

template <class ScalarType> class QuadRegionHistory {
  std::vector<int> RegionId; // per face
  int NextId = 0;
  int LastIdVal = -1;

public:
  // a new quad mesh (e.g. a full quadrangulation): no region yet
  void Reset(size_t NumFaces) {
    RegionId.assign(NumFaces, -1);
    LastIdVal = -1;
  }

  // the id of the last region remesh (-1: none)
  int LastId() const { return LastIdVal; }

  // the region of every face
  const std::vector<int> &Ids() const { return RegionId; }

  // the faces still coming from the region remesh Id
  std::vector<int> Faces(const int Id) const {
    std::vector<int> F;
    if (Id < 0)
      return F;
    for (size_t i = 0; i < RegionId.size(); i++)
      if (RegionId[i] == Id)
        F.push_back((int)i);
    return F;
  }

  // after a region remesh: the faces in Replaced (indices before) are gone,
  // the others kept in the same order, then NumNew new faces at the end
  // (NumFacesAfter in all); the new faces get a new id, returned
  int AfterRemesh(const std::set<int> &Replaced, size_t NumFacesAfter,
                  size_t NumNew) {
    std::vector<int> NewId;
    for (size_t i = 0; i < RegionId.size(); i++)
      if (Replaced.count((int)i) == 0)
        NewId.push_back(RegionId[i]);
    // (a history out of date with the mesh: the kept faces unknown)
    NewId.resize(NumFacesAfter - NumNew, -1);
    LastIdVal = NextId++;
    for (size_t i = 0; i < NumNew; i++)
      NewId.push_back(LastIdVal);
    RegionId = NewId;
    return LastIdVal;
  }

  // after a rebuild that keeps some faces: OldFaceOf[i] is the index before
  // of face i (-1 = a new face), Label the label of every face after (e.g.
  // its layout patch), OldLabel the label of every face before. A new face
  // takes the id its label had most before the rebuild (a patch tessellated
  // again keeps the region of the patch; before, a patch tessellated again
  // all, with no kept face, lost it). Without OldLabel the kept faces with
  // its label decide; with empty Label the new faces get -1
  void AfterRebuild(const std::vector<int> &OldFaceOf,
                    const std::vector<int> &Label,
                    const std::vector<int> *OldLabel = NULL) {
    std::map<int, std::map<int, int>> LabelIdCount;
    if ((OldLabel != NULL) && (OldLabel->size() == RegionId.size())) {
      for (size_t i = 0; i < RegionId.size(); i++)
        LabelIdCount[(*OldLabel)[i]][RegionId[i]]++;
    } else if (Label.size() == OldFaceOf.size())
      for (size_t i = 0; i < OldFaceOf.size(); i++)
        if ((OldFaceOf[i] >= 0) && (OldFaceOf[i] < (int)RegionId.size()))
          LabelIdCount[Label[i]][RegionId[OldFaceOf[i]]]++;
    std::vector<int> NewId(OldFaceOf.size(), -1);
    for (size_t i = 0; i < OldFaceOf.size(); i++) {
      int O = OldFaceOf[i];
      if ((O >= 0) && (O < (int)RegionId.size())) {
        NewId[i] = RegionId[O];
        continue;
      }
      if (Label.size() != OldFaceOf.size())
        continue;
      int Best = -1, BestN = 0;
      std::map<int, int> &C = LabelIdCount[Label[i]];
      for (std::map<int, int>::iterator It = C.begin(); It != C.end(); It++)
        if ((*It).second > BestN) {
          BestN = (*It).second;
          Best = (*It).first;
        }
      NewId[i] = Best;
    }
    RegionId = NewId;
  }
};

} // namespace Geo

#endif
