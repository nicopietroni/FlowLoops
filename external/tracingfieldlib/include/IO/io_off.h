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

#ifndef IO_OFF
#define IO_OFF

#include <cstddef>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "../point3.h"
#include "file_path_functions.h"


template <class ScalarType>
bool WriteOFF(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<std::vector<int> > &Connectivity)
{
    std::ofstream out(meshFile);
    if (!out.is_open())
        return false;

    out << "OFF\n";
    out << VertPos.size() << " " << Connectivity.size() << " 0\n";

    for (const auto& p : VertPos)
        out << p.X << " " << p.Y << " " << p.Z << "\n";

    for (const auto& face : Connectivity) {
        out << face.size();
        for (int idx : face)
            out << " " << idx;
        out << "\n";
    }

    out.close();
    return true;
}

#endif
