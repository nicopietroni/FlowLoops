#pragma once

// Wrapper around tetgen.h for header-only mode.
//
// tetgen.h defines a number of short, generic macros (Orient, DOT, CROSS,
// DIST, SWAP2, ...) that it only needs for its own inline members. In
// header-only mode CycleDiscovery is compiled in the same translation unit as
// TracingFieldLib / libigl / Directional, where those names are ordinary
// identifiers (e.g. "std::vector<int> Orient(K, 0);" in mesh_clean.h), so the
// macros must not leak. tetgen.cxx includes tetgen.h directly and is not
// affected by this file.

#include "tetgen.h"

#undef fastlookup
#undef Div2
#undef Mod2
#undef Orient
#undef EdgeRing
#undef smarktest
#undef sunmarktest
#undef smarktested
#undef NORM2
#undef DIST
#undef DOT
#undef CROSS
#undef SETVECTOR3
#undef SWAP2
