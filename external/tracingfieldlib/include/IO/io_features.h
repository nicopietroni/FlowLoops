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

#ifndef IO_FEATURES
#define IO_FEATURES

#include <IO/file_path_functions.h>
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <point3.h>
#include <sstream>
#include <string>
#include <vector>

inline bool LoadClist(const std::string &filename,
                      std::vector<std::pair<int, int>> &Features,
                      bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Loading:" << filename.c_str() << std::endl;
  std::ifstream file(filename.c_str());

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false;
  }

  std::string line;

  // Read each line from the file
  int CurrV = 0;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<int> numbers;
    int NextV;

    // Split the line by spaces and store each number in the vector
    while (iss >> NextV) {
      Features.push_back(std::pair<int, int>(CurrV, NextV));
      // numbers.push_back(num);
    }

    // Process or output the numbers as needed (example: print the numbers)
    // std::cout << "Read numbers: ";
    // for (const int& number : numbers) {
    //     //std::cout << number << " ";
    // }
    // std::cout << std::endl;
    CurrV++;
  }

  file.close();
  return true;
}

inline bool WriteFeatures(const std::string &filename,
                          const std::vector<std::pair<int, int>> &Features,
                          bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  for (const auto &feature : Features) {
    outfile << feature.first << " " << feature.second << "\n";
    if (!outfile) {
      return false; // Write operation failed
    }
  }

  outfile.close();
  return true;
}

inline bool ReadFeatures(const std::string &filename,
                         std::vector<std::pair<int, int>> &Features,
                         bool writeInfo = true) {

  if (writeInfo)
    std::cout << "Reading:" << filename.c_str() << std::endl;

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  Features.clear(); // Clear the vector to ensure it's empty before reading new data

  int first, second;
  while (infile >> first >> second) {
    Features.emplace_back(first, second);
  }

  std::sort(Features.begin(), Features.end());
  auto last = std::unique(Features.begin(), Features.end());
  Features.erase(last, Features.end());

  if (!infile.eof()) {
    // An error occurred before reaching the end of the file
    return false;
  }

  infile.close();
  return true;
}

inline bool ReadSharpFile(const std::string &filename,
                   const std::vector<std::vector<int>> &Faces,
                   std::vector<std::pair<int, int>> &Features,
                   bool writeInfo = true)
{
    if (writeInfo)
      std::cout << "Reading:" << filename.c_str() << std::endl;

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Cannot open file\n";
        return false;
    }

    int numEntries;
    file >> numEntries;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    //entries.reserve(numEntries);

    std::string line;
    Features.clear();
    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        int type,face,edge;
        char comma;

        ss >> type >> comma >> face >> comma >> edge;
        int v0 = Faces[face][edge];
        int v1 = Faces[face][(edge + 1) % Faces[face].size()];
        std::pair<int, int> e = std::make_pair(std::min(v0, v1), std::max(v0, v1));
        Features.push_back(e);
    }

    std::sort(Features.begin(), Features.end());
    auto last = std::unique(Features.begin(), Features.end());
    Features.erase(last, Features.end());

    return true;
}

template< class ScalarType>
inline bool
SaveFeatureCoord(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::pair<int, int>> &Features,
                 const std::string PathSave) {
  std::ofstream file(PathSave.c_str());
  if (!file.is_open()) {
    std::cerr << "Error: Unable to open file: feature_coords.txt" << std::endl;
    return false;
  }
  // write first number of features
  file << Features.size() << std::endl;
  for (size_t i = 0; i < Features.size(); i++) {
    int v0 = Features[i].first;
    int v1 = Features[i].second;
    Geo::Point3<ScalarType> P0 = VertPos[v0];
    Geo::Point3<ScalarType> P1 = VertPos[v1];
    file << P0.X << " " << P0.Y << " " << P0.Z << " " << P1.X << " " << P1.Y
         << " " << P1.Z << std::endl;
  }
  file.close();
  return true;
}

#endif
