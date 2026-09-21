#pragma once
#include <iostream>

// Library-wide logging for Curve-Surfacing-Core.
//
//   cs_log() — progress/info stream. Gated by the flag below, which the
//              curve_surfacing_core() entry points set from
//              CurveSurfacingParams::enable_logging. Silent by default,
//              writes to stdout when enabled.
//   cs_err() — error stream. Always on, writes to stderr. Use for genuine
//              failures the caller should see regardless of logging.
//
// Header-only; safe in both build modes (header-only and
// CURVE_SURFACING_STATIC).

namespace Geex {

inline bool& cs_logging_enabled() {
    static bool enabled = false;
    return enabled;
}

struct CsLog {
    template <typename T>
    CsLog& operator<<(const T& msg) {
        if (cs_logging_enabled()) std::cout << msg;
        return *this;
    }
    CsLog& operator<<(std::ostream& (*manip)(std::ostream&)) {
        if (cs_logging_enabled()) manip(std::cout);
        return *this;
    }
};

inline CsLog cs_log() { return CsLog{}; }

inline std::ostream& cs_err() { return std::cerr; }

} // namespace Geex
