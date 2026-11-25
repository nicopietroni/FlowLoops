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

#ifndef POLYLINE_FITTING_H
#define POLYLINE_FITTING_H

// BezierFitting (moved here from polyline3.h, unchanged) and SketchPolyline:
// a line drawn by the user made ready to be embedded, once: cleaned,
// resampled uniformly, closed if it is almost a loop, smoothed by local
// Bezier fitting. The result is the line used everywhere afterwards (the
// embedding does not smooth it again).

#include "point3.h"
#include "polyline3.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>

namespace Geo {

template <class ScalarType> class BezierFitting {

  // We assume Point3<ScalarType> supports:
  //  - +, -, *, / with scalars
  //  - Dot(p)  (or Dot(p,q))
  //  - optionally Norm(), otherwise we compute sqrt(p.Dot(p))

  // Clamp a value in [0,1]
  static inline ScalarType Clamp01(ScalarType x) {
    if (x < ScalarType(0))
      return ScalarType(0);
    if (x > ScalarType(1))
      return ScalarType(1);
    return x;
  }

  // // Compute Euclidean norm of a 3D vector
  // template <class ScalarType, class Point3>
  // static inline ScalarType Norm(const Point3& p)
  // {
  //   return std::sqrt(p.Dot(p));
  // }

  // Evaluate a quadratic Bezier curve at parameter t
  // B(t) = (1 - t)^2 P0 + 2 (1 - t) t C + t^2 P2
  static inline Point3<ScalarType>
  EvalQuadraticBezier(const Point3<ScalarType> &P0, const Point3<ScalarType> &C,
                      const Point3<ScalarType> &P2, ScalarType t) {
    t = Clamp01(t);

    ScalarType u = ScalarType(1) - t;
    ScalarType a = u * u;
    ScalarType b = ScalarType(2) * u * t;
    ScalarType c = t * t;

    return P0 * a + C * b + P2 * c;
  }

  // Compute chord-length parameterization for an open polyline
  // The parameter t_i is proportional to the accumulated arc-length
  static inline std::vector<ScalarType>
  ChordLengthParams(const std::vector<Point3<ScalarType>> &Q) {
    const size_t n = Q.size();
    std::vector<ScalarType> t(n, ScalarType(0));

    if (n <= 1)
      return t;

    ScalarType totalLength = ScalarType(0);
    for (size_t i = 1; i < n; ++i)
      totalLength += (Q[i] - Q[i - 1]).Norm();
    // totalLength += Norm<ScalarType, Point3<ScalarType>>(Q[i] - Q[i - 1]);

    // If the polyline is degenerate, fall back to uniform parameters
    if (totalLength <= ScalarType(1e-20)) {
      for (size_t i = 0; i < n; ++i)
        t[i] = (n == 1) ? ScalarType(0) : ScalarType(i) / ScalarType(n - 1);
      return t;
    }

    ScalarType acc = ScalarType(0);
    t[0] = ScalarType(0);

    for (size_t i = 1; i < n; ++i) {
      // acc += Norm<ScalarType, Point3<ScalarType>>(Q[i] - Q[i - 1]);
      acc += (Q[i] - Q[i - 1]).Norm();
      t[i] = acc / totalLength;
    }

    t.back() = ScalarType(1);
    return t;
  }

  // Fit the best quadratic Bezier curve (least squares)
  // Endpoints are fixed: P0 = Q[0], P2 = Q[last]
  // Only the control point C is solved
  static inline bool FitQuadraticBezierControlPoint_LS(
      const std::vector<Point3<ScalarType>> &Q, Point3<ScalarType> &P0_out,
      Point3<ScalarType> &C_out, Point3<ScalarType> &P2_out,
      ScalarType eps = ScalarType(1e-20)) {
    const size_t n = Q.size();
    if (n < 2)
      return false;

    const Point3<ScalarType> P0 = Q.front();
    const Point3<ScalarType> P2 = Q.back();

    // With only two points the curve is a straight line
    // Any C on the segment is valid; choose the midpoint
    if (n == 2) {
      C_out = (P0 + P2) * ScalarType(0.5);
      P0_out = P0;
      P2_out = P2;
      return true;
    }

    std::vector<ScalarType> t = ChordLengthParams(Q);

    ScalarType denom = ScalarType(0);
    Point3<ScalarType> numer = P0 * ScalarType(0); // zero vector

    // Solve least squares in closed form
    for (size_t i = 0; i < n; ++i) {
      const ScalarType ti = Clamp01(t[i]);
      const ScalarType u = ScalarType(1) - ti;

      const ScalarType a = u * u;
      const ScalarType b = ScalarType(2) * u * ti;
      const ScalarType c = ti * ti;

      // Endpoints do not constrain C because b = 0 there
      if (std::abs(b) <= ScalarType(1e-30))
        continue;

      // We want: a*P0 + b*C + c*P2 ≈ Q[i]
      // => b*C ≈ Q[i] - a*P0 - c*P2
      const Point3<ScalarType> rhs = Q[i] - P0 * a - P2 * c;

      numer = numer + rhs * b;
      denom += b * b;
    }

    if (denom <= eps) {
      // Degenerate configuration: fallback to midpoint
      C_out = (P0 + P2) * ScalarType(0.5);
    } else {
      C_out = numer / denom;
    }

    P0_out = P0;
    P2_out = P2;
    return true;
  }

public:
  static inline std::vector<Point3<ScalarType>> SmoothLocalBezierPolyline(
      const std::vector<Point3<ScalarType>> &polyIn, int iters, int winSize,
      int deg, bool sampleUniformT, bool fitUseChordLength = true) {
    std::vector<Point3<ScalarType>> poly = polyIn;

    const int N = (int)poly.size();
    if (N < 3)
      return poly;

    if (deg < 1)
      return poly;
    if (winSize < deg + 1)
      return poly;
    if (winSize < 2)
      return poly;

    const bool isLoop = (poly.front() == poly.back());

    // Work on unique points if loop
    int n = N;
    if (isLoop) {
      if (N < 4)
        return poly; // need at least 3 unique + closure
      n = N - 1;
    } else {
      if (N < 3)
        return poly;
    }

    // If winSize is even, we still handle it, but "centered" window will be
    // asymmetric by 1.
    const int halfW = winSize / 2;

    // Helper to build window indices around a center
    auto buildWindow = [&](int center, std::vector<int> &indices,
                           std::vector<Point3<ScalarType>> &window) {
      indices.clear();
      window.clear();
      indices.reserve(winSize);
      window.reserve(winSize);

      // centered range: [center-halfW, center+(winSize-halfW-1)]
      // for odd winSize, this is symmetric.
      const int left = center - halfW;
      const int right = center + (winSize - halfW - 1);

      for (int idx = left; idx <= right; ++idx) {
        int id = idx;
        if (isLoop) {
          id = (id % n + n) % n;
        } else {
          if (id < 0 || id >= n)
            continue;
        }
        indices.push_back(id);
        window.push_back(poly[id]);
      }
    };

    for (int it = 0; it < iters; ++it) {

      // Copy current unique vertices into P
      std::vector<Point3<ScalarType>> P(n);
      for (int i = 0; i < n; ++i)
        P[i] = poly[i];

      std::vector<Point3<ScalarType>> accum(n, P[0] * ScalarType(0));
      std::vector<int> counts(n, 0);

      std::vector<Point3<ScalarType>> window;
      std::vector<int> indices;

      // Slide windows centered at each vertex (excluding endpoints for open)
      for (int center = 0; center < n; ++center) {

        // For open polylines, keep endpoints fixed: do not generate windows
        // that would move them.
        if (!isLoop && (center < halfW || center >= n - halfW))
          continue;

        buildWindow(center, indices, window);

        if ((int)window.size() < deg + 1)
          continue;

        // Fit + sample on the window (fix endpoints inside the window = true)
        std::vector<Point3<ScalarType>> fitted =
            BezierFitting<ScalarType>::SampleFittedBezierPolyline(
                window, deg, sampleUniformT, fitUseChordLength,
                /*fixEndpoints=*/true);

        if (fitted.size() != window.size())
          continue;

        // Accumulate predictions for each vertex in this window
        for (size_t k = 0; k < indices.size(); ++k) {
          const int id = indices[k];
          accum[id] = accum[id] + fitted[k];
          counts[id] += 1;
        }
      }

      // Average predictions
      for (int i = 0; i < n; ++i) {
        if (counts[i] > 0) {
          P[i] = accum[i] / ScalarType(counts[i]);
        } else {
          // untouched vertex: keep original
          P[i] = poly[i];
        }
      }

      // Write back, preserving closure and open endpoints
      if (isLoop) {
        poly.resize(n + 1);
        for (int i = 0; i < n; ++i)
          poly[i] = P[i];
        poly.back() = poly.front();
      } else {
        // endpoints unchanged by construction (we skipped their centers)
        for (int i = 0; i < n; ++i)
          poly[i] = P[i];
      }
    }

    return poly;
  }

  static inline Geo::Polyline3<ScalarType>
  SmoothLocalBezierPolyline(const Geo::Polyline3<ScalarType> &polyIn, int iters,
                            int winSize, int deg, bool sampleUniformT,
                            bool fitUseChordLength = true) {
    Geo::Polyline3<ScalarType> polySmooth;
    polySmooth.PolyPos = SmoothLocalBezierPolyline(
        polyIn.PolyPos, iters, winSize, deg, sampleUniformT, fitUseChordLength);
    return polySmooth;
  }
  // ===================== ADDED: general-degree Bezier fitting & sampling
  // =====================

private:
  // Uniform params on [0,1]
  static inline std::vector<ScalarType> UniformParams(size_t n) {
    std::vector<ScalarType> t(n, ScalarType(0));
    if (n <= 1)
      return t;
    for (size_t i = 0; i < n; ++i)
      t[i] = ScalarType(i) / ScalarType(n - 1);
    return t;
  }

  // Binomial coefficient C(n,k)
  static inline ScalarType Binom(int n, int k) {
    if (k < 0 || k > n)
      return ScalarType(0);
    if (k == 0 || k == n)
      return ScalarType(1);
    k = std::min(k, n - k);
    ScalarType res = ScalarType(1);
    for (int i = 1; i <= k; ++i) {
      res *= ScalarType(n - (k - i));
      res /= ScalarType(i);
    }
    return res;
  }

  // Bernstein basis for all j=0..deg at parameter t
  static inline void BernsteinAll(int deg, ScalarType t,
                                  std::vector<ScalarType> &B) {
    B.assign(deg + 1, ScalarType(0));
    t = Clamp01(t);
    const ScalarType u = ScalarType(1) - t;

    for (int j = 0; j <= deg; ++j) {
      const ScalarType c = Binom(deg, j);
      B[j] = c * std::pow(u, ScalarType(deg - j)) * std::pow(t, ScalarType(j));
    }
  }

  // Evaluate Bezier (degree = P.size()-1)
  static inline Point3<ScalarType>
  EvalBezier(const std::vector<Point3<ScalarType>> &P, ScalarType t) {
    const int deg = int(P.size()) - 1;
    std::vector<ScalarType> B;
    BernsteinAll(deg, t, B);

    Point3<ScalarType> out = P[0] * ScalarType(0);
    for (int j = 0; j <= deg; ++j)
      out = out + P[j] * B[j];
    return out;
  }

  // Solve SPD linear system M x = b with Cholesky.
  // M is kxk (ScalarType), b and x are length k (Point3)
  static inline bool
  SolveSPD_Cholesky(const std::vector<std::vector<ScalarType>> &M,
                    const std::vector<Point3<ScalarType>> &b,
                    std::vector<Point3<ScalarType>> &x,
                    ScalarType eps = ScalarType(1e-20)) {

    const int k = int(M.size());
    x.assign(k,
             Point3<ScalarType>(ScalarType(0), ScalarType(0), ScalarType(0)));
    if (k == 0)
      return true;

    std::vector<std::vector<ScalarType>> L(
        k, std::vector<ScalarType>(k, ScalarType(0)));

    for (int i = 0; i < k; ++i) {
      for (int j = 0; j <= i; ++j) {
        ScalarType sum = M[i][j];
        for (int s = 0; s < j; ++s)
          sum -= L[i][s] * L[j][s];

        if (i == j) {
          if (sum <= eps)
            return false;
          L[i][j] = std::sqrt(sum);
        } else {
          L[i][j] = sum / L[j][j];
        }
      }
    }

    // Forward: L y = b
    std::vector<Point3<ScalarType>> y(
        k, Point3<ScalarType>(ScalarType(0), ScalarType(0), ScalarType(0)));
    for (int i = 0; i < k; ++i) {
      Point3<ScalarType> sum = b[i];
      for (int s = 0; s < i; ++s)
        sum = sum - y[s] * L[i][s];
      y[i] = sum / L[i][i];
    }

    // Backward: L^T x = y
    for (int i = k - 1; i >= 0; --i) {
      Point3<ScalarType> sum = y[i];
      for (int s = i + 1; s < k; ++s)
        sum = sum - x[s] * L[s][i];
      x[i] = sum / L[i][i];
    }

    return true;
  }

  // Fit Bezier of degree 'deg' by least squares.
  // If fixEndpoints=true, P0=Q.front and Pdeg=Q.back are fixed, we solve only
  // internal controls. fitUseChordLength decides how to parameterize the data
  // points for the LS.
  static inline bool FitBezierControlPoints_LS(
      const std::vector<Point3<ScalarType>> &Q, int deg,
      std::vector<Point3<ScalarType>> &P_out, bool fitUseChordLength = true,
      bool fixEndpoints = true, ScalarType ridge = ScalarType(1e-10),
      ScalarType eps = ScalarType(1e-20)) {

    const int m = int(Q.size());
    if (m < 2)
      return false;
    if (deg < 1)
      return false;

    // Data parameterization for fitting
    std::vector<ScalarType> t =
        fitUseChordLength ? ChordLengthParams(Q) : UniformParams(Q.size());

    P_out.assign(deg + 1, Point3<ScalarType>(ScalarType(0), ScalarType(0),
                                             ScalarType(0)));

    Point3<ScalarType> P0 = Q.front();
    Point3<ScalarType> Pn = Q.back();

    // Unknown count
    const int k = fixEndpoints ? std::max(0, deg - 1) : (deg + 1);

    if (fixEndpoints) {
      P_out[0] = P0;
      P_out[deg] = Pn;
    }

    if (k == 0) {
      // deg==1 with fixed endpoints: straight segment
      P_out[0] = P0;
      P_out[deg] = Pn;
      return true;
    }

    // Normal equations: M = A^T A (kxk), b = A^T rhs (k Point3)
    std::vector<std::vector<ScalarType>> M(
        k, std::vector<ScalarType>(k, ScalarType(0)));
    std::vector<Point3<ScalarType>> b(
        k, Point3<ScalarType>(ScalarType(0), ScalarType(0), ScalarType(0)));

    std::vector<ScalarType> B;
    for (int i = 0; i < m; ++i) {
      BernsteinAll(deg, t[i], B);

      Point3<ScalarType> rhs = Q[i];
      if (fixEndpoints) {
        rhs = rhs - P0 * B[0] - Pn * B[deg];
      }

      for (int u = 0; u < k; ++u) {
        const int j_u = fixEndpoints ? (u + 1) : u;
        const ScalarType a_u = B[j_u];

        b[u] = b[u] + rhs * a_u;

        for (int v = 0; v <= u; ++v) {
          const int j_v = fixEndpoints ? (v + 1) : v;
          const ScalarType a_v = B[j_v];
          M[u][v] += a_u * a_v;
        }
      }
    }

    // Symmetrize + ridge
    for (int u = 0; u < k; ++u) {
      for (int v = 0; v < u; ++v)
        M[v][u] = M[u][v];
      M[u][u] += ridge;
    }

    std::vector<Point3<ScalarType>> x;
    if (!SolveSPD_Cholesky(M, b, x, eps)) {
      // Fallback: linear guess along segment
      for (int j = 0; j <= deg; ++j) {
        const ScalarType a = ScalarType(j) / ScalarType(deg);
        P_out[j] = P0 * (ScalarType(1) - a) + Pn * a;
      }
      return true;
    }

    for (int u = 0; u < k; ++u) {
      const int j = fixEndpoints ? (u + 1) : u;
      P_out[j] = x[u];
    }

    if (!fixEndpoints) {
      // if endpoints are free, they are already part of x and were set above
    }

    return true;
  }

public:
  // Sample the fitted quadratic Bezier at the same number of points
  // and with the same chord-length parameterization as the input polyline
  static inline std::vector<Point3<ScalarType>>
  SampleFittedQuadraticBezierPolyline(
      const std::vector<Point3<ScalarType>> &Q) {
    Point3<ScalarType> P0, C, P2;
    std::vector<Point3<ScalarType>> out;

    if (!FitQuadraticBezierControlPoint_LS(Q, P0, C, P2))
      return out;

    std::vector<ScalarType> t = ChordLengthParams(Q);
    out.resize(Q.size());

    for (size_t i = 0; i < Q.size(); ++i)
      out[i] = EvalQuadraticBezier(P0, C, P2, t[i]);

    return out;
  }

  // Sample the fitted quadratic Bezier at the same number of points
  // and with the same chord-length parameterization as the input polyline
  static inline Polyline3<ScalarType>
  SampleFittedQuadraticBezierPolyline(const Polyline3<ScalarType> &PolyL) {
    std::vector<Point3<ScalarType>> Q = PolyL.PolyPos;
    Polyline3<ScalarType> outPoly;
    SampleFittedQuadraticBezierPolyline(Q);
    outPoly.PolyPos = SampleFittedQuadraticBezierPolyline(Q);
    return outPoly;
  }

  // Like SampleFittedQuadraticBezierPolyline, but for arbitrary degree.
  // - deg: degree of Bezier
  // - sampleUniformT:
  //     true  -> sample t uniformly in [0,1] on the final Bezier
  //     false -> sample using chord-length params of the input polyline (like
  //     your quadratic method)
  // - fitUseChordLength:
  //     if true, uses chord-length params for fitting (recommended)
  //     if false, uses uniform params for fitting
  // - fixEndpoints: if true, forces P0=Q.front and Pdeg=Q.back (like your
  // quadratic)
  static inline std::vector<Point3<ScalarType>>
  SampleFittedBezierPolyline(const std::vector<Point3<ScalarType>> &Q, int deg,
                             bool sampleUniformT, bool fitUseChordLength = true,
                             bool fixEndpoints = true) {
    std::vector<Point3<ScalarType>> out;
    if (Q.size() < 2)
      return out;
    if (deg < 1)
      return out;

    std::vector<Point3<ScalarType>> P;
    if (!FitBezierControlPoints_LS(Q, deg, P, fitUseChordLength, fixEndpoints))
      return out;

    std::vector<ScalarType> ts =
        sampleUniformT ? UniformParams(Q.size()) : ChordLengthParams(Q);

    out.resize(Q.size());
    for (size_t i = 0; i < Q.size(); ++i)
      out[i] = EvalBezier(P, ts[i]);

    return out;
  }

  static inline Polyline3<ScalarType>
  SampleFittedBezierPolyline(const Polyline3<ScalarType> &PolyL, int deg,
                             bool sampleUniformT, bool fitUseChordLength = true,
                             bool fixEndpoints = true) {
    Polyline3<ScalarType> outPoly;
    outPoly.PolyPos = SampleFittedBezierPolyline(
        PolyL.PolyPos, deg, sampleUniformT, fitUseChordLength, fixEndpoints);
    return outPoly;
  }
};
// // Local Bezier smoothing using sliding windows + averaging.
// // - iters: number of smoothing passes
// // - winSize: window size (recommended odd). Must be >= deg+1.
// // - deg: Bezier degree
// // - sampleUniformT: if true, sample fitted Bezier uniformly in t; else use
// // chord-length of the window
// // - fitUseChordLength: parameterization used for FITTING (recommended true)
// //
// // Behavior:
// // - If polyline is a loop (first==last): preserves closure.
// // - If open: endpoints are preserved (never changed).
// static inline Polyline3<ScalarType>
// SmoothLocalBezierPolyline(const Polyline3<ScalarType> &polyIn, int iters,
//                           int winSize, int deg, bool sampleUniformT = true,
//                           bool fitUseChordLength = true) {
//   Polyline3<ScalarType> out = polyIn;
//   out.PolyPos = SmoothLocalBezierPolyline(polyIn.PolyPos, iters, winSize,
//   deg,
//                                           sampleUniformT, fitUseChordLength);
//   return out;
// }

template <class ScalarType> struct SketchPolylineParam {
  // distance between the samples (absolute); at least MinSamples samples
  ScalarType Step = -1;
  int MinSamples = 10;
  // closure: the closest approach between the first and the last
  // CloseSpanRatio of the line (a hand drawn circle overshoots or stops
  // short) closer than CloseDistRatio times the length, the line going on
  // in the same direction there (tangents within CloseAngleDeg) -> the line
  // is cut there and closed
  bool close_loops = true;
  ScalarType CloseSpanRatio = 0.2;
  ScalarType CloseDistRatio = 0.1;
  ScalarType CloseAngleDeg = 60;
  // The fitting used
  enum FitMode { FitBSpline = 0, FitLocalBezier = 1 };
  FitMode Mode = FitBSpline;
  // FitBSpline: one least squares cubic B-spline through all the samples,
  // a control point every ControlSpacing (absolute) along the line; the
  // extremes of an open line are kept, a closed line is periodic
  ScalarType ControlSpacing = -1;
  // FitLocalBezier: local Bezier smoothing
  // (BezierFitting::SmoothLocalBezierPolyline);
  // 0 iterations = off
  int SmoothIters = 3;
  int SmoothWindow = 7;
  // When > 0 the window covers this length (in
  // samples: WindowLength / Step, odd, at least SmoothWindow); a window of a
  // few samples smoothed only the scale of one or two samples
  ScalarType WindowLength = -1;
  int SmoothDegree = 3;
};

template <class ScalarType> class SketchPolyline {
  typedef Geo::Point3<ScalarType> CoordType;

  static ScalarType Length(const std::vector<CoordType> &P) {
    ScalarType L = 0;
    for (size_t i = 1; i < P.size(); i++)
      L += (P[i] - P[i - 1]).Norm();
    return L;
  }

public:
  // N samples uniformly by arc length, the extremes included (a closed line,
  // first == last, stays closed)
  static void Resample(const std::vector<CoordType> &In, int N,
                       std::vector<CoordType> &Out) {
    Out.clear();
    if (In.size() < 2 || N < 2) {
      Out = In;
      return;
    }
    std::vector<ScalarType> Arc(In.size(), 0);
    for (size_t i = 1; i < In.size(); i++)
      Arc[i] = Arc[i - 1] + (In[i] - In[i - 1]).Norm();
    ScalarType L = Arc.back();
    if (!(L > 0)) {
      Out = In;
      return;
    }
    size_t Seg = 0;
    for (int k = 0; k < N; k++) {
      ScalarType S = L * (ScalarType)k / (ScalarType)(N - 1);
      while ((Seg + 2 < In.size()) && (Arc[Seg + 1] < S))
        Seg++;
      ScalarType D = Arc[Seg + 1] - Arc[Seg];
      ScalarType T = (D > 0) ? (S - Arc[Seg]) / D : 0;
      T = std::max((ScalarType)0, std::min((ScalarType)1, T));
      Out.push_back(In[Seg] * (1 - T) + In[Seg + 1] * T);
    }
    Out.front() = In.front();
    Out.back() = In.back();
  }

  // If the line is almost a loop, the
  // indices I (near the start) and J (near the end) where to cut it and
  // close it; false otherwise
  static bool LoopCut(const std::vector<CoordType> &P,
                      const SketchPolylineParam<ScalarType> &Param, int &I,
                      int &J) {
    int N = (int)P.size();
    if (N < 8)
      return false;
    if (P.front() == P.back())
      return false; // already closed
    ScalarType L = Length(P);
    if (!(L > 0))
      return false;
    int Span = std::max(2, (int)std::floor(Param.CloseSpanRatio * N));
    ScalarType BestD = std::numeric_limits<ScalarType>::max();
    I = -1;
    J = -1;
    for (int i = 0; i < Span && i < N; i++)
      for (int j = std::max(i + N / 2, N - 1 - Span); j < N; j++) {
        ScalarType D = (P[i] - P[j]).Norm();
        if (D < BestD) {
          BestD = D;
          I = i;
          J = j;
        }
      }
    if ((I < 0) || (BestD > Param.CloseDistRatio * L))
      return false;
    // the tangents there (central differences, one sided at the extremes)
    auto Tangent = [&](int k) {
      CoordType T = P[std::min(k + 1, N - 1)] - P[std::max(k - 1, 0)];
      T.Normalize();
      return T;
    };
    return (Geo::Dot(Tangent(I), Tangent(J)) >
            std::cos(Param.CloseAngleDeg * (ScalarType)M_PI / 180));
  }

  // Least squares cubic B-spline through the points P
  // (parametrized by arc length), M control points uniformly spaced
  // (M from the length and Spacing). Open: clamped knots, the first and the
  // last control points fixed on the extremes (the curve passes through
  // them). Closed (first == last): periodic. Out: the curve evaluated at the
  // parameters of the points (same number). False if too few points.
  static bool FitCubicBSpline(const std::vector<CoordType> &P, bool closed,
                              ScalarType Spacing,
                              std::vector<CoordType> &Out) {
    Out.clear();
    std::vector<CoordType> Q = P;
    if (closed && (Q.size() > 1) && (Q.front() == Q.back()))
      Q.pop_back();
    int N = (int)Q.size();
    if (N < 5)
      return false;
    // arc length parameters in [0,1] (closed: the closing segment too)
    std::vector<ScalarType> S(N, 0);
    for (int i = 1; i < N; i++)
      S[i] = S[i - 1] + (Q[i] - Q[i - 1]).Norm();
    ScalarType L = S.back() + (closed ? (Q.front() - Q.back()).Norm() : 0);
    if (!(L > 0))
      return false;
    for (int i = 0; i < N; i++)
      S[i] /= L;
    // number of control points
    int Spans = std::max(closed ? 3 : 1, (int)std::floor(0.5 + L / Spacing));
    int M = closed ? Spans : Spans + 3;
    // (no more unknowns than points)
    M = std::min(M, closed ? N - 1 : N);
    if (M < (closed ? 3 : 4))
      return false;
    // basis values of every point: (control index, weight)
    std::vector<std::vector<std::pair<int, ScalarType>>> B(N);
    auto Uniform = [](ScalarType u, ScalarType b[4]) {
      ScalarType u2 = u * u, u3 = u2 * u;
      b[0] = (1 - u) * (1 - u) * (1 - u) / 6;
      b[1] = (3 * u3 - 6 * u2 + 4) / 6;
      b[2] = (-3 * u3 + 3 * u2 + 3 * u + 1) / 6;
      b[3] = u3 / 6;
    };
    if (closed) {
      for (int i = 0; i < N; i++) {
        ScalarType T = S[i] * M;
        int Seg = std::min((int)std::floor(T), M - 1);
        ScalarType b[4];
        Uniform(T - Seg, b);
        for (int k = 0; k < 4; k++)
          B[i].push_back(std::make_pair((Seg + k) % M, b[k]));
      }
    } else {
      // clamped uniform knots: 0,0,0,0,1,...,M-4,M-3,M-3,M-3,M-3
      int D = M - 3; // domain [0, D]
      std::vector<ScalarType> Knot;
      for (int k = 0; k < 4; k++)
        Knot.push_back(0);
      for (int k = 1; k < D; k++)
        Knot.push_back((ScalarType)k);
      for (int k = 0; k < 4; k++)
        Knot.push_back((ScalarType)D);
      for (int i = 0; i < N; i++) {
        ScalarType T = S[i] * D;
        // Cox - de Boor, degree 3
        int Span = 3;
        while ((Span < M - 1) && !(T < Knot[Span + 1]))
          Span++;
        ScalarType Nb[4] = {1, 0, 0, 0};
        ScalarType Left[4], Right[4];
        for (int j = 1; j <= 3; j++) {
          Left[j] = T - Knot[Span + 1 - j];
          Right[j] = Knot[Span + j] - T;
          ScalarType Saved = 0;
          for (int r = 0; r < j; r++) {
            ScalarType Den = Right[r + 1] + Left[j - r];
            ScalarType Tmp = (Den != 0) ? Nb[r] / Den : 0;
            Nb[r] = Saved + Right[r + 1] * Tmp;
            Saved = Left[j - r] * Tmp;
          }
          Nb[j] = Saved;
        }
        for (int k = 0; k < 4; k++)
          B[i].push_back(std::make_pair(Span - 3 + k, Nb[k]));
      }
    }
    // least squares on the free control points (open: the first and the
    // last fixed on the extremes)
    int First = closed ? 0 : 1, Last = closed ? M - 1 : M - 2;
    int K = Last - First + 1;
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> AtA =
        Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Zero(K, K);
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> AtB =
        Eigen::Matrix<ScalarType, Eigen::Dynamic, 3>::Zero(K, 3);
    CoordType C0 = Q.front(), CM = Q.back();
    for (int i = 0; i < N; i++) {
      CoordType Rhs = Q[i];
      for (size_t a = 0; a < B[i].size(); a++) {
        int c = B[i][a].first;
        if (!closed && (c == 0))
          Rhs = Rhs - C0 * B[i][a].second;
        if (!closed && (c == M - 1))
          Rhs = Rhs - CM * B[i][a].second;
      }
      for (size_t a = 0; a < B[i].size(); a++) {
        int ca = B[i][a].first;
        if ((ca < First) || (ca > Last))
          continue;
        ScalarType wa = B[i][a].second;
        AtB(ca - First, 0) += wa * Rhs.X;
        AtB(ca - First, 1) += wa * Rhs.Y;
        AtB(ca - First, 2) += wa * Rhs.Z;
        for (size_t b = 0; b < B[i].size(); b++) {
          int cb = B[i][b].first;
          if ((cb < First) || (cb > Last))
            continue;
          AtA(ca - First, cb - First) += wa * B[i][b].second;
        }
      }
    }
    // (a tiny ridge: control points with few samples stay defined)
    for (int k = 0; k < K; k++)
      AtA(k, k) += (ScalarType)1e-9;
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> X = AtA.ldlt().solve(AtB);
    std::vector<CoordType> Ctrl(M);
    if (!closed) {
      Ctrl[0] = C0;
      Ctrl[M - 1] = CM;
    }
    for (int k = 0; k < K; k++)
      Ctrl[First + k] = CoordType(X(k, 0), X(k, 1), X(k, 2));
    // evaluated at the parameters of the points
    for (int i = 0; i < N; i++) {
      CoordType Pt(0, 0, 0);
      for (size_t a = 0; a < B[i].size(); a++)
        Pt = Pt + Ctrl[B[i][a].first] * B[i][a].second;
      Out.push_back(Pt);
    }
    if (closed)
      Out.push_back(Out.front());
    return true;
  }

  // the drawn line made ready: resampled, closed if almost a loop, smoothed;
  // Closed (optional) says if it came out closed (first == last)
  static void Process(const std::vector<CoordType> &In,
                      const SketchPolylineParam<ScalarType> &Param,
                      std::vector<CoordType> &Out, bool *Closed = NULL) {
    // drop repeated points
    std::vector<CoordType> P;
    for (size_t i = 0; i < In.size(); i++)
      if (P.empty() || !(In[i] == P.back()))
        P.push_back(In[i]);
    if (P.size() < 2) {
      Out = P;
      if (Closed != NULL)
        *Closed = false;
      return;
    }
    ScalarType L = Length(P);
    auto NumSamples = [&](ScalarType Len) {
      int N = Param.MinSamples;
      if (Param.Step > 0)
        N = std::max(N, (int)std::floor(0.5 + Len / Param.Step) + 1);
      return N;
    };
    std::vector<CoordType> R;
    Resample(P, NumSamples(L), R);
    // closure
    bool closed = (R.front() == R.back());
    int CutI = -1, CutJ = -1;
    if (!closed && Param.close_loops && LoopCut(R, Param, CutI, CutJ)) {
      // What sticks out cut, then closed
      R = std::vector<CoordType>(R.begin() + CutI, R.begin() + CutJ + 1);
      R.push_back(R.front());
      std::vector<CoordType> R2;
      Resample(R, NumSamples(Length(R)), R2);
      R = R2;
      R.back() = R.front();
      closed = true;
    }
    // One least squares fit
    if ((Param.Mode == SketchPolylineParam<ScalarType>::FitBSpline) &&
        (Param.ControlSpacing > 0)) {
      std::vector<CoordType> F;
      if (FitCubicBSpline(R, closed, Param.ControlSpacing, F)) {
        std::vector<CoordType> R2;
        Resample(F, NumSamples(Length(F)), R2);
        if (closed)
          R2.back() = R2.front();
        R = R2;
      }
    }
    // local Bezier smoothing (endpoints of an open line stay; a loop stays
    // closed)
    int Win = Param.SmoothWindow;
    if ((Param.WindowLength > 0) && (Param.Step > 0)) {
      Win = std::max(Win, (int)std::floor(0.5 + Param.WindowLength / Param.Step));
      if ((Win % 2) == 0)
        Win++;
    }
    // (a window longer than the line: the whole line, odd)
    Win = std::min(Win, ((int)R.size() % 2 == 1) ? (int)R.size() - 2
                                                 : (int)R.size() - 1);
    if ((Param.Mode == SketchPolylineParam<ScalarType>::FitLocalBezier) &&
        (Param.SmoothIters > 0) && (Win > Param.SmoothDegree) &&
        ((int)R.size() > Win)) {
      R = BezierFitting<ScalarType>::SmoothLocalBezierPolyline(
          R, Param.SmoothIters, Win, Param.SmoothDegree,
          /*sampleUniformT=*/false, /*fitUseChordLength=*/true);
      // The fitting shortens the zigzags, the samples
      // are no longer evenly spaced
      std::vector<CoordType> R2;
      Resample(R, NumSamples(Length(R)), R2);
      if (closed)
        R2.back() = R2.front();
      R = R2;
    }
    Out = R;
    if (Closed != NULL)
      *Closed = closed;
  }
};

} // namespace Geo

#endif
