cmake_minimum_required(VERSION 3.16)
project(loop_weaver CXX)

# ==============================
# C++ Standard
# ==============================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==============================
# Build mode
# ==============================
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
option(ENABLE_SANITIZER "Enable sanitizers" OFF)
option(ENABLE_NATIVE_ARCH "Tune Release builds for this exact CPU (not portable)" OFF)

message(STATUS "Build mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Sanitizer: ${ENABLE_SANITIZER}")

# ==============================
# Base flags
# ==============================
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
if(ENABLE_NATIVE_ARCH)
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native")
endif()

# ==============================
# Sanitizer flags (can be used in BOTH Debug and Release)
# ==============================
if(ENABLE_SANITIZER)
  add_compile_options(
    -g
    -O1
    -fno-omit-frame-pointer
    -fsanitize=address,undefined
    -fsanitize-address-use-after-scope
    -fno-optimize-sibling-calls
  )
  add_link_options(
    -fsanitize=address,undefined
  )
endif()

# ==============================
# Vendored dependencies
#
# Everything lives under external/ — no submodules, no network access during
# the build. See external/README.md for upstream sources and pinned commits.
# ==============================
set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/external")
set(TRACEFIELD_DIR "${EXTERNAL_DIR}/tracingfieldlib/include")
set(EIGENPATH      "${EXTERNAL_DIR}/eigen-3.4.0")
set(IGLLIB         "${EXTERNAL_DIR}/libigl/include")
set(IMGUI_DIR      "${EXTERNAL_DIR}/imgui")
set(TRIANGLE_DIR   "${EXTERNAL_DIR}/triangle")
set(CURVE_SURFACING_CORE_DIR "${EXTERNAL_DIR}/curve-surfacing-core")
set(TETGEN_DIR "${CURVE_SURFACING_CORE_DIR}/include/curve_surfacing/CycleDiscovery/tetgen")

foreach(_dep
    "${TRACEFIELD_DIR}"
    "${EIGENPATH}/Eigen"
    "${IGLLIB}"
    "${IMGUI_DIR}"
    "${TRIANGLE_DIR}/triangle.cpp"
    "${CURVE_SURFACING_CORE_DIR}/include"
    "${TETGEN_DIR}")
  if(NOT EXISTS "${_dep}")
    message(FATAL_ERROR
      "Missing vendored dependency: ${_dep}\n"
      "The external/ directory is incomplete. See external/README.md.")
  endif()
endforeach()

# Runtime lookup table, baked in at compile time.
set(LUT_PATH "${TRACEFIELD_DIR}/Lut.txt")
add_definitions(-DLUT_PATH=\"${LUT_PATH}\")

# ==============================
# System dependencies
# ==============================
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)

# Use the vendored Eigen 3.4.0. A system Eigen must NOT be substituted:
# 5.x is incompatible with this code and with the vendored libigl.
if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3_Eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS Eigen3_Eigen)
    target_include_directories(Eigen3_Eigen INTERFACE "${EIGENPATH}")
endif()

# ==============================
# Vendored static libraries
# ==============================

# Triangle (Shewchuk), via the libigl fork. Built here rather than through its
# own CMakeLists, which declares cmake_minimum_required(VERSION 3.1) and is
# rejected outright by CMake 4.x.
add_library(triangle STATIC "${TRIANGLE_DIR}/triangle.cpp")
target_include_directories(triangle PUBLIC "${TRIANGLE_DIR}")
target_compile_definitions(triangle PRIVATE TRILIBRARY ANSI_DECLARATORS)
set_target_properties(triangle PROPERTIES POSITION_INDEPENDENT_CODE ON)

# TetGen 1.4, vendored inside Curve-Surfacing-Core.
add_library(tetgen STATIC
    "${TETGEN_DIR}/tetgen.cxx"
    "${TETGEN_DIR}/predicates.cxx"
)
target_include_directories(tetgen PUBLIC "${TETGEN_DIR}")
target_compile_definitions(tetgen PRIVATE TETLIBRARY)
set_target_properties(tetgen PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ==============================
# Executable
# ==============================
add_executable(${PROJECT_NAME}
    main.cpp
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
)

# ==============================
# Includes
# ==============================
# Eigen must come before everything else.
target_include_directories(${PROJECT_NAME} BEFORE PRIVATE ${EIGENPATH})

target_include_directories(${PROJECT_NAME} PRIVATE
    ${TRACEFIELD_DIR}
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${CURVE_SURFACING_CORE_DIR}/include
    ${CURVE_SURFACING_CORE_DIR}/include/curve_surfacing
    ${TETGEN_DIR}
    ${TRIANGLE_DIR}
    ${CURVE_SURFACING_CORE_DIR}/external/Directional/include
    ${IGLLIB}
)

# ==============================
# Link libraries
# Curve-Surfacing-Core itself is consumed header-only — its headers self-include
# their .cpp files — so there is no libcurvesurfacing-core.a to link.
# ==============================
target_link_libraries(${PROJECT_NAME}
    tetgen
    triangle
    OpenGL::GL
    OpenGL::GLU
    GLEW::GLEW
    glfw
    Eigen3::Eigen
)
