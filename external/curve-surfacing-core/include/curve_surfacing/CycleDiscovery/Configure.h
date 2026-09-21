#ifndef _CONF_H_
#define _CONF_H_

#include "Errors.h"
#include <cmath>

#define DO_EXP false
#define GO_CMD 1

#define SAVE_FACE 0
#define SAVE_TILE 0
#define SAVE_NEWCURVE 0

// Use constexpr for type-safe, full-precision constants
constexpr double PI   = 3.14159265358979323846;
constexpr double hfPI = 1.57079632679489661923;

// Edge protection tuning parameters.
// BADEDGE_LIMIT: max unprotected edges before input is flagged as bad.
// plainEPS: coplanarity tolerance (relative to curve scale).
// plainPTB: per-step perturbation magnitude for near-planar degenerate inputs.
#define BADEDGE_LIMIT 30
#define plainEPS 0.001
#define plainPTB 0.00000001

//measurements
#define USEONLYWORSTBITRI 0

#define SAVEOBJ 0

#endif
