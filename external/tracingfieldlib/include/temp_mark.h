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

#ifndef TEMP_MARK
#define TEMP_MARK

#include <vector>

namespace Geo
{
    struct TMark
    {
        std::vector<int> ElMark;
        int GlobMark;

        void Allocate(const int &numEl)
        {
            assert(numEl >= 0);
            ElMark = std::vector<int>(numEl, GlobMark);
        }

        void AllocateIfNeeded(const int &numEl)
        {
            assert(numEl >= 0);
            if (ElMark.size()!=numEl)
                ElMark = std::vector<int>(numEl, GlobMark);
        }

        void UnMarkAll()
        {
            GlobMark++;
        }

        void Mark(const int &IndexEl)
        {
            assert(IndexEl >= 0);
            assert(IndexEl < ElMark.size());
            ElMark[IndexEl]=GlobMark;
        }

        bool IsMarked(const int &IndexEl)
        {
            assert(IndexEl >= 0);
            assert(IndexEl < ElMark.size());
            return (ElMark[IndexEl]==GlobMark);
        }

        TMark()
        {
            GlobMark = 0;
        }
    };
}

#endif
