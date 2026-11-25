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

#ifndef IO_STL
#define IO_STL

#include <cstddef>
#include <cstdint>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <point3.h>
#include <mesh_functions.h>

// Read a binary STL (one separate vertex triple per facet).
// Returns false if the file size does not match a binary STL, so the
// caller can try the ascii format (some binary files start with "solid").
template <class ScalarType>
bool LoadSTLBinary(const std::string &meshFile,
                   std::vector<Geo::Point3<ScalarType> > &VertPos,
                   std::vector<std::vector<int> > &Connectivity,
                   std::vector<Geo::Point3<ScalarType> > &FaceNormals,
                   bool strictSize=true)
{
    VertPos.clear();
    Connectivity.clear();
    FaceNormals.clear();

    std::ifstream fin(meshFile.c_str(), std::ios::binary);
    if(!fin) return false;

    fin.seekg(0, std::ios::end);
    std::streamoff fileSize = fin.tellg();
    if (fileSize < 84) return false;

    // 80 bytes header, then the number of facets
    fin.seekg(80, std::ios::beg);
    uint32_t numF = 0;
    fin.read(reinterpret_cast<char *>(&numF), 4);
    if (!fin) return false;
    // OLD: only files of the exact binary size
    // if ((std::streamoff)84 + (std::streamoff)50 * (std::streamoff)numF != fileSize)
    //     return false;
    if ((std::streamoff)84 + (std::streamoff)50 * (std::streamoff)numF != fileSize)
    {
        if (strictSize) return false;
        // TOLERANT, but not on text: a header made only of printable
        // characters is an ascii STL the ascii reader could not parse,
        // reading it as binary gives random triangles
        std::vector<char> Head(84);
        fin.seekg(0, std::ios::beg);
        fin.read(Head.data(), 84);
        bool IsText = true;
        for (size_t k = 0; k < Head.size(); k++)
            if (!(std::isprint((unsigned char)Head[k]) || std::isspace((unsigned char)Head[k])))
                IsText = false;
        if (IsText) return false;
        fin.clear();
        fin.seekg(84, std::ios::beg);
        // TOLERANT: wrong facet count in the header (trailing bytes or a
        // truncated file): read the facets that fit in the file
        uint32_t numFit = (uint32_t)((fileSize - 84) / 50);
        std::cout << "Warning: binary STL size mismatch (header " << numF
                  << " facets, file " << numFit << "), reading " << numFit
                  << std::endl;
        numF = numFit;
    }

    VertPos.reserve((size_t)numF * 3);
    Connectivity.reserve(numF);
    FaceNormals.reserve(numF);
    for (uint32_t i = 0; i < numF; i++)
    {
        // normal + 3 vertices (12 floats, little endian) + 2 bytes attribute
        float Data[12];
        uint16_t Attr;
        fin.read(reinterpret_cast<char *>(Data), 48);
        fin.read(reinterpret_cast<char *>(&Attr), 2);
        if (!fin) return false;

        FaceNormals.push_back(Geo::Point3<ScalarType>(Data[0], Data[1], Data[2]));
        std::vector<int> FIdx;
        for (size_t j = 0; j < 3; j++)
        {
            VertPos.push_back(Geo::Point3<ScalarType>(Data[3 + 3 * j],
                                                      Data[4 + 3 * j],
                                                      Data[5 + 3 * j]));
            FIdx.push_back(VertPos.size() - 1);
        }
        Connectivity.push_back(FIdx);
    }
    fin.close();
    return true;
}

// Read one number token of an ascii STL: false (value 0) when it is not a
// finite number. The token is always consumed, so the parsing goes on.
template <class ScalarType>
bool ReadSTLNumber(std::istream &fin, ScalarType &Val)
{
    Val = 0;
    std::string tok;
    if (!(fin >> tok)) return false;
    char *end = NULL;
    double d = std::strtod(tok.c_str(), &end);
    if ((end == tok.c_str()) || (*end != 0) || (!std::isfinite(d))) return false;
    Val = (ScalarType)d;
    return true;
}

// Read an ascii STL (one separate vertex triple per facet).
// Facets with a number of vertices other than 3 are skipped.
template <class ScalarType>
bool LoadSTLAscii(const std::string &meshFile,
                  std::vector<Geo::Point3<ScalarType> > &VertPos,
                  std::vector<std::vector<int> > &Connectivity,
                  std::vector<Geo::Point3<ScalarType> > &FaceNormals,
                  bool writeInfo=false)
{
    VertPos.clear();
    Connectivity.clear();
    FaceNormals.clear();

    std::ifstream fin(meshFile.c_str());
    if(!fin) return false;

    Geo::Point3<ScalarType> CurrN(0, 0, 0);
    std::vector<Geo::Point3<ScalarType> > CurrVert;
    int skipped = 0;
    // ROBUST: keywords in any case (some exporters write FACET, VERTEX),
    // numbers read as tokens (nan, inf or broken values do not stop the
    // stream as with fin >> double); a facet with a bad vertex is skipped
    bool BadFacet = false;
    std::string str;
    while (fin >> str)
    {
        for (size_t k = 0; k < str.size(); k++)
            str[k] = (char)std::tolower((unsigned char)str[k]);
        if (str == std::string("facet"))
        {
            CurrVert.clear();
            BadFacet = false;
            CurrN = Geo::Point3<ScalarType>(0, 0, 0);
            std::string key;
            fin >> key;
            for (size_t k = 0; k < key.size(); k++)
                key[k] = (char)std::tolower((unsigned char)key[k]);
            if (key == std::string("normal"))
            {
                // OLD: fin >> CurrN.X >> CurrN.Y >> CurrN.Z;
                // a bad normal is harmless, it is only kept as info
                ScalarType N[3];
                for (int k = 0; k < 3; k++)
                    ReadSTLNumber(fin, N[k]);
                CurrN = Geo::Point3<ScalarType>(N[0], N[1], N[2]);
            }
            continue;
        }
        if (str == std::string("vertex"))
        {
            // OLD: a bad number stopped the whole file
            // Geo::Point3<ScalarType> P;
            // fin >> P.X >> P.Y >> P.Z;
            // if (!fin) return false;
            ScalarType V[3];
            bool ok = true;
            for (int k = 0; k < 3; k++)
                ok = ReadSTLNumber(fin, V[k]) && ok;
            if (!fin) break;
            if (!ok) BadFacet = true;
            CurrVert.push_back(Geo::Point3<ScalarType>(V[0], V[1], V[2]));
            continue;
        }
        if (str == std::string("endfacet"))
        {
            // if (CurrVert.size() != 3)
            if ((CurrVert.size() != 3) || (BadFacet))
            {
                skipped++;
                CurrVert.clear();
                continue;
            }
            std::vector<int> FIdx;
            for (size_t j = 0; j < 3; j++)
            {
                VertPos.push_back(CurrVert[j]);
                FIdx.push_back(VertPos.size() - 1);
            }
            Connectivity.push_back(FIdx);
            FaceNormals.push_back(CurrN);
            CurrVert.clear();
        }
    }
    if ((skipped > 0) && (writeInfo))
        std::cout << "Skipped " << skipped << " non triangular or invalid facets" << std::endl;
    fin.close();
    return (Connectivity.size() > 0);
}

// Load a binary or ascii STL. With mergeVert (default) coincident vertices
// are unified by position (RemoveDuplicatedVert) and the faces that become
// degenerate are removed (RemoveDegenerateConnectivityTriangles).
template <class ScalarType>
bool LoadSTL(std::string meshFile,
             std::vector<Geo::Point3<ScalarType> > &VertPos,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false,
             bool mergeVert=true)
{
    if (writeInfo)
        std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;

    std::vector<Geo::Point3<ScalarType> > FaceNormals;
    bool has_loaded = LoadSTLBinary(meshFile, VertPos, Connectivity, FaceNormals);
    if (!has_loaded)
        has_loaded = LoadSTLAscii(meshFile, VertPos, Connectivity, FaceNormals,
                                  writeInfo);
    // TOLERANT: neither format read: binary with a wrong facet count
    if (!has_loaded)
        has_loaded = LoadSTLBinary(meshFile, VertPos, Connectivity, FaceNormals,
                                   false);
    if (!has_loaded)
    {
        std::cout << "Error: cannot read the STL file " << meshFile << std::endl;
        VertPos.clear();
        Connectivity.clear();
        return false;
    }

    size_t NumReadF = Connectivity.size();
    if (mergeVert)
    {
        RemoveDuplicatedVert(VertPos, Connectivity);
        RemoveDegenerateConnectivityTriangles(VertPos, Connectivity,
                                              FaceNormals);
    }

    if (writeInfo)
    {
        std::cout<<"Read:"<<VertPos.size()<<" vertices"<<std::endl;
        std::cout<<"Read:"<<Connectivity.size()<<" faces"<<std::endl;
    }
    if (Connectivity.empty())
        std::cout << "Error: no valid faces in the STL (" << NumReadF
                  << " facets read, all degenerate or empty file)" << std::endl;
    return (Connectivity.size() > 0);
}

#endif
