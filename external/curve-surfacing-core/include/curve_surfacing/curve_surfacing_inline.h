#pragma once

// When CURVE_SURFACING_STATIC is defined, functions are compiled into a
// static library and this macro expands to nothing — no 'inline' keyword,
// so each definition is compiled exactly once into the lib.
//
// Otherwise (header-only mode, the default), every function definition that
// appears in a .cpp file #included by its corresponding .h must be marked
// 'inline' to prevent duplicate-symbol linker errors when multiple
// translation units include the same header.
//
// Usage: prefix every free function and method definition in the .cpp files
// with CS_INLINE, e.g.:
//
//   CS_INLINE void CurvNet::clear() { ... }
//   CS_INLINE Mesh3D* CurvNet::load_ma(std::string filename) { ... }

#ifdef CURVE_SURFACING_STATIC
#  define CS_INLINE
#else
#  define CS_INLINE inline
#endif