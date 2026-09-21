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

#ifndef FILE_PATH_FUNCTIONS
#define FILE_PATH_FUNCTIONS

#include <string>

inline std::string GetDirectory(const std::string &curr_path)
{
    std::string projM=curr_path;
    size_t indexExt=projM.find_last_of("/");
    projM=projM.substr(0,indexExt);
    return(projM);
}

inline std::string RemoveExtension(const std::string &curr_path)
{
    std::string projM=curr_path;
    size_t indexExt=projM.find_last_of(".");
    projM=projM.substr(0,indexExt);
    return(projM);
}

inline std::string GetFilaname(const std::string &curr_path)
{
    std::string projM=curr_path;
    size_t indexExt=projM.find_last_of("/");
    projM=projM.substr(indexExt+1);
    return(projM);
}

inline std::string GetFilanameNoExtension(const std::string &curr_path)
{
    std::string nameF=GetFilaname(curr_path);
    nameF=RemoveExtension(nameF);
    return nameF;
}

inline std::string GetFileExtension(const std::string &curr_path)
{
    size_t lastDotPos = curr_path.find_last_of(".");

    std::string fileExtension;

    // Check if a dot was found and extract the file extension
    if (lastDotPos != std::string::npos)
        fileExtension = curr_path.substr(lastDotPos + 1);

    return fileExtension;
}

#endif
