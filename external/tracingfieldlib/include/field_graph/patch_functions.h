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

#ifndef PATCH_FUNCTIONS
#define PATCH_FUNCTIONS

#include "patch_managing.h"
#include <algorithm>
#include <cassert>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class PatchFunctions {

using Edge = std::pair<int,int>;

// Build list of vertices from a side: v0, v1, ..., vL (edges are (v0,v1),(v1,v2),...)
static std::vector<int> vertices_from_side(const std::vector<Edge>& side) {
    std::vector<int> V;
    if (side.empty()) return V;
    V.reserve(side.size()+1);
    V.push_back(side.front().first);
    for (const auto& e : side) V.push_back(e.second);
    return V;
}

// Rebuild a side from a vertex list
static std::vector<Edge> side_from_vertices(const std::vector<int>& V) {
    std::vector<Edge> S;
    if (V.size() < 2) return S;
    S.resize(V.size()-1);
    for (size_t i=0;i+1<V.size();++i) S[i] = {V[i], V[i+1]};
    return S;
}

// Split a side at an EXISTING interior vertex near the middle (by edge count).
// Returns that existing vertex id and replaces sides[idx] with two sides:
//   left:  start ... mid_vertex
//   right: mid_vertex ... end
// If the side has <2 edges (no interior vertex), returns -1 and does nothing.
static int split_side_at_mid_existing(std::vector<std::vector<Edge>>& sides, size_t idx) {
    auto& S = sides[idx];
    const int L = (int)S.size();          // number of edges
    if (L < 2) return -1;                 // no interior vertex to split on

    auto V = vertices_from_side(S);       // |V| = L+1
    int k = L / 2;                        // floor
    if (k <= 0) k = 1;
    if (k >= L) k = L-1;

    const int mid_vertex = V[k];          // EXISTING vertex id

    // left:  V[0..k]
    // right: V[k..L]
    std::vector<int> leftV(V.begin(),      V.begin()+k+1);
    std::vector<int> rightV(V.begin()+k,   V.end());

    std::vector<Edge> left  = side_from_vertices(leftV);
    std::vector<Edge> right = side_from_vertices(rightV);

    S.swap(left);
    sides.insert(sides.begin()+idx+1, std::move(right));
    return mid_vertex;
}

// Helper: rebuild corners from side starts (authoritative order)
static void rebuild_corners_from_sides(const std::vector<std::vector<Edge>>& sides,
                                       std::vector<int>& corners)
{
    corners.clear();
    corners.reserve(sides.size());
    for (const auto& s : sides) {
        if (s.empty()) throw std::runtime_error("Side has no edges.");
        corners.push_back(s.front().first);
    }
}

public:
/**
 * Ensures 3 (preferred) up to 4 corners/sides by splitting the longest sides at existing interior vertices.
 * - Uses only existing vertex IDs from the edge sequences.
 * - Preserves the starts of current sides as current corners.
 * - First targets 3 sides; then, if possible, 4. Stops early if no side has >=2 edges.
 *
 * Precondition: sides are in ring order and contiguous:
 *   sides[i].back().second == sides[(i+1)%k].front().first
 *
 * @return final number of sides/corners (2, 3, or 4)
 */
static int make_four_corners_by_splitting(std::vector<std::vector<Edge>>& sides,
                                         std::vector<int>& corners)
{
    if (sides.empty())
        throw std::runtime_error("ensure_up_to_four_corners_from_sides: sides must not be empty.");

    // Align corners with side starts if sizes differ
    if (corners.size() != sides.size())
        rebuild_corners_from_sides(sides, corners);

    // We only consider starting cases 0,1,2 per your spec; if >2, we just align and return.
    if (corners.size() > 2) {
        // Trim to first 4 if you want, but spec says we don't consider >=3 as input case.
        // Keep as-is:
        return (int)sides.size();
    }

    auto split_once_if_possible = [&](bool& did_split)->void {
        did_split = false;
        // pick the longest side with L>=2
        size_t best = (size_t)-1;
        int bestLen = 1; // strictly >1
        for (size_t i=0;i<sides.size();++i) {
            int L = (int)sides[i].size();
            if (L > bestLen) { bestLen = L; best = i; }
        }
        if (best != (size_t)-1) {
            int new_corner = split_side_at_mid_existing(sides, best);
            if (new_corner != -1) {
                // insert the new corner as start of the new (right) side
                corners.insert(corners.begin() + (best + 1), new_corner);
                did_split = true;
            }
        }
    };

    // Target 3 sides first
    while (sides.size() < 3) {
        bool did = false;
        split_once_if_possible(did);
        if (!did) break; // cannot reach 3
    }

    // Then target 4 sides
    while (sides.size() < 4) {
        bool did = false;
        split_once_if_possible(did);
        if (!did) break; // cannot reach 4 -> stay at 3 (or 2)
    }

    // Align corners with sides in case of any discrepancy
    if (corners.size() != sides.size())
        rebuild_corners_from_sides(sides, corners);

    // Optional continuity/asserts (no throw in release)
#ifndef NDEBUG
    for (size_t i=0;i<sides.size();++i) {
        if (sides[i].empty()) throw std::runtime_error("Empty side found.");
        const int last = sides[i].back().second;
        const int next0 = sides[(i+1)%sides.size()].front().first;
        if (last != next0) throw std::runtime_error("Continuity broken.");
        if (corners[i] != sides[i].front().first) throw std::runtime_error("Corner/start mismatch.");
    }
#endif

    return (int)sides.size();
}

// Remove sides (and corresponding corners) until their number <= max_corners.
// The function modifies both `sides` and `corners` in place, preserving order.
static void remove_corners_until_from_sides(
    std::vector<std::vector<Edge>>& sides,
    std::vector<int>& corners,
    int max_corners)
{
    if ((int)sides.size() <= max_corners) return; // nothing to remove

    // While there are too many sides (== corners)
    while ((int)sides.size() > max_corners) {
        // Find the shortest side (fewest edges)
        size_t shortest_idx = 0;
        size_t shortest_len = sides[0].size();
        for (size_t i = 1; i < sides.size(); ++i) {
            if (sides[i].size() < shortest_len) {
                shortest_len = sides[i].size();
                shortest_idx = i;
            }
        }

        // Merge this shortest side with the *previous* one in cyclic order
        size_t prev_idx = (shortest_idx == 0) ? sides.size() - 1 : shortest_idx - 1;

        // Get vertex lists for both sides
        auto leftV  = vertices_from_side(sides[prev_idx]);
        auto rightV = vertices_from_side(sides[shortest_idx]);

        // Merge: previous side keeps its start, shortest side keeps its end
        // Ensure shared corner vertex (rightV.front() == leftV.back())
        if (leftV.back() != rightV.front()) {
            throw std::runtime_error("Non-contiguous sides: cannot merge.");
        }

        // Concatenate (remove duplicate shared vertex)
        leftV.insert(leftV.end(), rightV.begin() + 1, rightV.end());

        // Rebuild merged side
        sides[prev_idx] = side_from_vertices(leftV);

        // Remove the merged (shortest) side
        sides.erase(sides.begin() + shortest_idx);

        // Update the corners vector coherently (remove the corner starting that side)
        rebuild_corners_from_sides(sides, corners);
    }
}

static void print_sides_and_corners(const std::vector<std::vector<Edge>>& sides,
                             const std::vector<int>& corners)
{
    std::cout << "Corners (" << corners.size() << "): ";
    for (size_t i = 0; i < corners.size(); ++i)
        std::cout << corners[i] << (i + 1 < corners.size() ? ", " : "");
    std::cout << "\n";

    std::cout << "Sides (" << sides.size() << "):\n";
    for (size_t i = 0; i < sides.size(); ++i) {
        std::cout << "  Side " << i << " (len=" << sides[i].size() << "): ";
        for (auto& e : sides[i])
            std::cout << "(" << e.first << "," << e.second << ") ";
        std::cout << "\n";
    }

    // Optional continuity check
    bool ok = true;
    for (size_t i = 0; i < sides.size(); ++i) {
        const auto last_v  = sides[i].back().second;
        const auto next_v0 = sides[(i + 1) % sides.size()].front().first;
        if (last_v != next_v0) {
            ok = false;
            std::cerr << "  !! continuity error between side " << i
                      << " and side " << (i + 1) % sides.size()
                      << " (" << last_v << " != " << next_v0 << ")\n";
        }
    }
    if (ok)
        std::cout << "Continuity OK ✅\n";
    else
        std::cout << "Continuity errors detected ⚠️\n";
}
};

} // namespace Geo
#endif
