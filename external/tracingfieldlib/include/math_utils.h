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

#ifndef MATH_UTILS
#define MATH_UTILS

#include <algorithm> // for std::minmax_element
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <utility>
#include <iostream>

// #include <math.h>
#define M_PI 3.14159265358979323846264338327950288 /* pi             */

namespace Geo {

// Hash function
struct hashPairIntFunction {
  size_t operator()(const std::pair<int, int> &x) const {
    int h1 = std::hash<int>{}(x.first);
    int h2 = std::hash<int>{}(x.second);
    // return h1 ^ h2;
    return h1 ^ (h2 << 1);
  }
};

struct hashPairSizeTFunction {
  std::size_t operator()(const std::pair<size_t, size_t> &p) const {
    return std::hash<size_t>()(p.first) ^ (std::hash<size_t>()(p.second) << 1);
  }
};

template <class ScalarType> inline ScalarType ToDeg(const ScalarType &a) {
  return a * (ScalarType)180.0f / ScalarType(M_PI);
}

template <class ScalarType> inline ScalarType normalizeAngle(ScalarType angle) {
  angle = (ScalarType)fmod((float)angle, (float)2 * M_PI);
  if (angle < 0) {
    angle += 2 * M_PI;
  }
  return angle;
}

template <class ScalarType>
inline ScalarType MaxPercentile(const std::vector<ScalarType> &Values,
                                const ScalarType &Percentile) {
  assert((Percentile >= 0) && (Percentile <= 1));

  // then compute the max (considering the percentile)
  std::vector<ScalarType> ValuesSorted = Values;
  std::sort(ValuesSorted.begin(), ValuesSorted.end());
  int IndexMax = floor(
      ((ScalarType)ValuesSorted.size() * ((ScalarType)1 - Percentile)) + 0.5);
  IndexMax = std::min(IndexMax, (int)ValuesSorted.size() - 1);
  IndexMax = std::max(IndexMax, 0);
  ScalarType TestVal = ValuesSorted[IndexMax];
  return TestVal;
}

//void PrintIntStats(const std::string &label, std::vector<int> V) {
inline void PrintIntStats(const std::string &label, std::vector<int> V) { // inline: defined in a header
  if (V.empty()) {
    std::cout << label << ": (empty)" << std::endl;
    return;
  }
  std::sort(V.begin(), V.end());
  size_t n = V.size();
  long sum = 0;
  for (int x : V)
    sum += x;
  double mean = (double)sum / n;
  int median = V[n / 2];
  int p10 = V[(size_t)(0.10 * (n - 1))];
  int p90 = V[(size_t)(0.90 * (n - 1))];
  std::cout << label << " (n=" << n << "): min=" << V.front()
            << " p10=" << p10 << " median=" << median << " mean=" << mean
            << " p90=" << p90 << " max=" << V.back() << std::endl;
}

template <class ScalarType>
inline void Trim(std::vector<ScalarType> &Values, const ScalarType &MaxVal) {
  for (size_t i = 0; i < Values.size(); i++)
    Values[i] = std::min(Values[i], MaxVal);
}

template <class ScalarType>
inline ScalarType LinearMap(ScalarType &Value, const ScalarType Min0,
                            const ScalarType Max0, const ScalarType Min1,
                            const ScalarType Max1) {
  assert(Max0 > Min0);
  Value = std::max(Value, Min0);
  Value = std::min(Value, Max0);

  assert(Max1 > Min1);
  ScalarType ratio = (Value - Min0) / (Max0 - Min0);
  return (Min1 + (Max1 - Min1) * ratio);
}

template <class ScalarType>
inline void LinearMap(std::vector<ScalarType> &Values, const ScalarType Min1,
                      const ScalarType Max1) {
  assert(Max1 > Min1);

  ScalarType MinV0 = std::numeric_limits<ScalarType>::max();
  ScalarType MaxV0 = std::numeric_limits<ScalarType>::lowest();
  for (size_t i = 0; i < Values.size(); i++) {
    MinV0 = std::min(MinV0, Values[i]);
    MaxV0 = std::max(MaxV0, Values[i]);
  }

  for (size_t i = 0; i < Values.size(); i++)
    Values[i] = LinearMap(Values[i], MinV0, MaxV0, Min1, Max1);
}

template <class ScalarType>
inline void TrimPercentile(std::vector<ScalarType> &Values,
                           const ScalarType minP, const ScalarType maxP) {
  assert(minP >= 0);
  assert(minP <= 1);
  assert(maxP >= 0);
  assert(maxP <= 1);

  std::vector<ScalarType> ValuesSorted = Values;
  std::sort(ValuesSorted.begin(), ValuesSorted.end());
  int IndexMax =
      floor(((ScalarType)ValuesSorted.size() * ((ScalarType)1 - maxP)) + 0.5);
  int IndexMin =
      floor(((ScalarType)ValuesSorted.size() * ((ScalarType)minP)) + 0.5);

  IndexMax = std::min(IndexMax, (int)ValuesSorted.size() - 1);
  IndexMax = std::max(IndexMax, 0);
  IndexMin = std::min(IndexMin, (int)ValuesSorted.size() - 1);
  IndexMin = std::max(IndexMin, 0);

  ScalarType ValueMin = ValuesSorted[IndexMin];
  ScalarType ValueMax = ValuesSorted[IndexMax];

  for (size_t i = 0; i < Values.size(); i++) {
    Values[i] = std::max(Values[i], ValueMin);
    Values[i] = std::min(Values[i], ValueMax);
  }
}

template <class ScalarType>
void ClusterValues(const std::vector<ScalarType> &Values, const bool &IsLoop,
                   const ScalarType &SplitThresold,
                   std::vector<std::vector<int>> &Clusters) {

  // start with one cluster per edge
  Clusters.clear();
  std::vector<ScalarType> ClusterAngles;
  for (size_t i = 0; i < Values.size(); i++) {
    Clusters.push_back({(int)i});
    ClusterAngles.push_back(Values[i]);
  }

  // iteratively merge clusters with similar angles with a greedy approach
  // that merges first the ones that have smallest difference in angle, for
  // each cluster average the angles of its edges weighted by number of edges
  // in the cluster stop when no more merges are possible under the angle
  // threshold also consider if the mesh is a loop and allow merging first and
  // last cluster in that case
  while (true) {
    int numClusters = (int)Clusters.size();
    if (numClusters <= 1)
      break;
    // find the pair of clusters with smallest angle difference
    int mergeA = -1, mergeB = -1;
    ScalarType minAngleDiff = std::numeric_limits<ScalarType>::max();
    for (int i = 0; i < numClusters; i++) {
      int next = (i + 1) % numClusters;
      if (!IsLoop && next == 0)
        continue;
      ScalarType angleDiff = std::abs(ClusterAngles[i] - ClusterAngles[next]);
      if (angleDiff < minAngleDiff) {
        minAngleDiff = angleDiff;
        mergeA = i;
        mergeB = next;
      }
    }
    if (minAngleDiff > SplitThresold)
      break;
    // merge clusters mergeA and mergeB

    ClusterAngles[mergeA] = (ClusterAngles[mergeA] * Clusters[mergeA].size() +
                             ClusterAngles[mergeB] * Clusters[mergeB].size()) /
                            (Clusters[mergeA].size() + Clusters[mergeB].size());

    Clusters[mergeA].insert(Clusters[mergeA].end(), Clusters[mergeB].begin(),
                            Clusters[mergeB].end());

    Clusters.erase(Clusters.begin() + mergeB);
    ClusterAngles.erase(ClusterAngles.begin() + mergeB);
  }
}

template <class ScalarType>
void ClusterValuesHysteresis(const std::vector<ScalarType> &Values, const bool &IsLoop,
                              const ScalarType &SplitThreshold,
                              std::vector<std::vector<int>> &Clusters,
                              const ScalarType &ToleranceRatio = ScalarType(0.2))
{
    Clusters.clear();
    const int N = static_cast<int>(Values.size());
    if (N == 0) return;

    const ScalarType Tolerance  = SplitThreshold * ToleranceRatio;
    const ScalarType ThreshHigh = SplitThreshold + Tolerance;
    const ScalarType ThreshLow  = SplitThreshold - Tolerance;

    std::vector<bool> IsHigh(N);
    bool State = (Values[0] > SplitThreshold);
    for (int i = 0; i < N; ++i)
    {
        if      (Values[i] > ThreshHigh) State = true;
        else if (Values[i] < ThreshLow)  State = false;
        IsHigh[i] = State;
    }

    int StartIdx = 0;
    if (IsLoop)
        for (int i = 0; i < N; ++i)
            if (IsHigh[i] != IsHigh[(i + 1) % N])
            { StartIdx = (i + 1) % N; break; }

    Clusters.push_back({});
    for (int k = 0; k < N; ++k)
    {
        int Curr = (StartIdx + k) % N;
        int Next = (StartIdx + k + 1) % N;
        Clusters.back().push_back(Curr);
        if (k == N - 1) break;
        if (IsHigh[Curr] != IsHigh[Next])
            Clusters.push_back({});
    }
}


template <class ScalarType>
void RefineClusterLloyd(const std::vector<ScalarType> &Values, const bool &IsLoop,
                        std::vector<std::vector<int>> &Clusters,
                        const int MaxIter = 20)
{
    const int N = static_cast<int>(Values.size());
    const int K = static_cast<int>(Clusters.size());
    if (K <= 1 || N == 0) return;

    // Ricostruisce una mappa index → cluster per accesso rapido
    std::vector<int> Assignment(N, 0);
    for (int c = 0; c < K; ++c)
        for (int idx : Clusters[c])
            Assignment[idx] = c;

    for (int iter = 0; iter < MaxIter; ++iter)
    {
        // 1. Calcola centroide di ogni cluster
        std::vector<ScalarType> Centroids(K, ScalarType(0));
        std::vector<int>        Counts(K, 0);
        for (int i = 0; i < N; ++i)
        {
            Centroids[Assignment[i]] += Values[i];
            Counts[Assignment[i]]++;
        }
        for (int c = 0; c < K; ++c)
            if (Counts[c] > 0)
                Centroids[c] /= ScalarType(Counts[c]);

        // 2. Sposta i confini tra cluster adiacenti
        //    Solo i punti di confine possono cambiare cluster (contiguità)
        bool Changed = false;
        for (int c = 0; c < K; ++c)
        {
            int CNext = (c + 1) % K;
            if (!IsLoop && CNext == 0) break; // sequenza aperta: stop all'ultimo

            // Punto di confine: ultimo del cluster c
            int BorderIdx = Clusters[c].back();
            int NextIdx   = Clusters[CNext].front();

            // Prova a spostare BorderIdx → CNext
            if (Clusters[c].size() > 1)
            {
                ScalarType DistCurr = std::abs(Values[BorderIdx] - Centroids[c]);
                ScalarType DistNext = std::abs(Values[BorderIdx] - Centroids[CNext]);
                if (DistNext < DistCurr)
                {
                    Clusters[c].pop_back();
                    Clusters[CNext].insert(Clusters[CNext].begin(), BorderIdx);
                    Assignment[BorderIdx] = CNext;
                    Changed = true;
                    continue;
                }
            }

            // Prova a spostare NextIdx → c
            if (Clusters[CNext].size() > 1)
            {
                ScalarType DistCurr = std::abs(Values[NextIdx] - Centroids[CNext]);
                ScalarType DistNext = std::abs(Values[NextIdx] - Centroids[c]);
                if (DistNext < DistCurr)
                {
                    Clusters[CNext].erase(Clusters[CNext].begin());
                    Clusters[c].push_back(NextIdx);
                    Assignment[NextIdx] = c;
                    Changed = true;
                }
            }
        }

        if (!Changed) break; // converge
    }
}

// template <class ScalarType>
// void ClusterValuesDP(const std::vector<ScalarType> &Values, const bool &IsLoop,
//                      const int &K,
//                      std::vector<std::vector<int>> &Clusters)
// {
//     Clusters.clear();
//     const int N = static_cast<int>(Values.size());
//     if (N == 0 || K <= 0) return;
//     const int ActualK = std::min(K, N);

//     // Precalcola costo (varianza) di ogni segmento [i, j]
//     // cost[i][j] = varianza dei valori da i a j
//     auto SegmentCost = [&](int i, int j) -> ScalarType
//     {
//         ScalarType Sum = ScalarType(0), SumSq = ScalarType(0);
//         int Len = j - i + 1;
//         for (int x = i; x <= j; ++x)
//         {
//             Sum   += Values[x];
//             SumSq += Values[x] * Values[x];
//         }
//         ScalarType Mean = Sum / ScalarType(Len);
//         return SumSq / ScalarType(Len) - Mean * Mean; // E[x²] - E[x]²
//     };

//     // DP: dp[k][i] = costo minimo per partizionare Values[0..i] in k cluster
//     const ScalarType INF = std::numeric_limits<ScalarType>::max();
//     std::vector<std::vector<ScalarType>> dp(ActualK + 1, std::vector<ScalarType>(N, INF));
//     std::vector<std::vector<int>>        split(ActualK + 1, std::vector<int>(N, -1));

//     // Base: 1 cluster = costo del segmento intero [0..i]
//     for (int i = 0; i < N; ++i)
//         dp[1][i] = SegmentCost(0, i);

//     // Riempie la tabella
//     for (int k = 2; k <= ActualK; ++k)
//         for (int i = k - 1; i < N; ++i)
//             for (int m = k - 2; m < i; ++m)
//             {
//                 ScalarType Cost = dp[k-1][m] + SegmentCost(m + 1, i);
//                 if (Cost < dp[k][i])
//                 {
//                     dp[k][i]    = Cost;
//                     split[k][i] = m;
//                 }
//             }

//     // Ricostruisce i cluster dal backtracking
//     std::vector<int> Breakpoints;
//     int idx = N - 1;
//     for (int k = ActualK; k >= 2; --k)
//     {
//         Breakpoints.push_back(split[k][idx]);
//         idx = split[k][idx];
//     }
//     std::reverse(Breakpoints.begin(), Breakpoints.end());

//     // Costruisce i cluster
//     int Start = 0;
//     for (int b : Breakpoints)
//     {
//         std::vector<int> Cluster;
//         for (int i = Start; i <= b; ++i) Cluster.push_back(i);
//         Clusters.push_back(Cluster);
//         Start = b + 1;
//     }
//     std::vector<int> Last;
//     for (int i = Start; i < N; ++i) Last.push_back(i);
//     Clusters.push_back(Last);

//     // Caso loop: riordina i cluster partendo da un breakpoint reale
//     if (IsLoop && Clusters.size() > 1)
//     {
//         // Trova il breakpoint con il maggiore salto di media tra cluster adiacenti
//         int BestBreak = 0;
//         ScalarType BestGap = ScalarType(0);
//         for (int c = 0; c < (int)Clusters.size(); ++c)
//         {
//             int Prev = (c - 1 + Clusters.size()) % Clusters.size();
//             ScalarType MeanC = ScalarType(0), MeanP = ScalarType(0);
//             for (int i : Clusters[c])    MeanC += Values[i];
//             for (int i : Clusters[Prev]) MeanP += Values[i];
//             MeanC /= Clusters[c].size();
//             MeanP /= Clusters[Prev].size();
//             ScalarType Gap = std::abs(MeanC - MeanP);
//             if (Gap > BestGap) { BestGap = Gap; BestBreak = c; }
//         }
//         std::rotate(Clusters.begin(), Clusters.begin() + BestBreak, Clusters.end());
//     }
// }

// template <class ScalarType>
// void ClusterValuesMaj(const std::vector<ScalarType> &Values, const bool &IsLoop,
//                       const ScalarType &SplitThreshold,
//                       std::vector<std::vector<int>> &Clusters) {
//   Clusters.clear();

//   const int n = (int)Values.size();
//   if (n == 0)
//     return;

//   // Start with one cluster per value.
//   for (int i = 0; i < n; ++i)
//     Clusters.push_back({i});

//   // Error of a cluster: maximum pairwise value difference inside the cluster.
//   auto ClusterError = [&](const std::vector<int> &cluster) {
//     ScalarType maxDiff = ScalarType(0);

//     for (size_t i = 0; i < cluster.size(); ++i) {
//       const int vi = cluster[i];
//       assert(vi >= 0 && vi < n);

//       for (size_t j = i + 1; j < cluster.size(); ++j) {
//         const int vj = cluster[j];
//         assert(vj >= 0 && vj < n);

//         const ScalarType diff = std::abs(Values[vi] - Values[vj]);
//         if (diff > maxDiff)
//           maxDiff = diff;
//       }
//     }

//     return maxDiff;
//   };

//   // Error after merging two clusters.
//   auto MergedClusterError = [&](const std::vector<int> &a,
//                                 const std::vector<int> &b) {
//     ScalarType maxDiff = ScalarType(0);

//     // Existing internal error of cluster a
//     maxDiff = std::max(maxDiff, ClusterError(a));

//     // Existing internal error of cluster b
//     maxDiff = std::max(maxDiff, ClusterError(b));

//     // Cross-cluster error
//     for (int ia : a) {
//       assert(ia >= 0 && ia < n);

//       for (int ib : b) {
//         assert(ib >= 0 && ib < n);

//         const ScalarType diff = std::abs(Values[ia] - Values[ib]);
//         if (diff > maxDiff)
//           maxDiff = diff;
//       }
//     }

//     return maxDiff;
//   };

//   // Greedily merge adjacent clusters.
//   // At each iteration, merge the adjacent pair with the smallest resulting
//   // maximum internal error, provided that this error is below the threshold.
//   while (true) {
//     const int numClusters = (int)Clusters.size();
//     if (numClusters <= 1)
//       break;

//     int mergeA = -1;
//     int mergeB = -1;
//     ScalarType bestMergeError = std::numeric_limits<ScalarType>::max();

//     for (int i = 0; i < numClusters; ++i) {
//       const int next = (i + 1) % numClusters;

//       // If this is not a loop, do not test last-first adjacency.
//       if (!IsLoop && next == 0)
//         continue;

//       const ScalarType err = MergedClusterError(Clusters[i], Clusters[next]);

//       if (err < bestMergeError) {
//         bestMergeError = err;
//         mergeA = i;
//         mergeB = next;
//       }
//     }

//     if (mergeA < 0 || mergeB < 0)
//       break;

//     if (bestMergeError > SplitThreshold)
//       break;

//     // Merge cluster B into cluster A.
//     Clusters[mergeA].insert(Clusters[mergeA].end(), Clusters[mergeB].begin(),
//                             Clusters[mergeB].end());

//     Clusters.erase(Clusters.begin() + mergeB);
//   }
// }

} // namespace Geo

#endif
