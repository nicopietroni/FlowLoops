#ifndef _ERRORS_H_
#define _ERRORS_H_

#include <stdexcept>
#include <string>

enum ShapeType {
	ERROR_CURVE,
	ERROR_NORM,
	ERROR_MATCH,
	ERROR_LIMIT,
	ERROR_TETGEN,
	ERROR_OTHERS
};

// Throws std::runtime_error with a descriptive message.
// Using const char* (not char*) since filename is never modified.
inline void errors(int type, const char* filename) {
	std::string msg;
	switch (type) {
	case ERROR_CURVE:
		msg = "Cannot open curve file";
		break;
	case ERROR_NORM:
		msg = "Cannot open normal file";
		break;
	case ERROR_MATCH:
		msg = "Number of points and number of normals do not match";
		break;
	case ERROR_LIMIT:
		msg = "Number of points exceeds the limit";
		break;
	case ERROR_TETGEN:
		msg = "TetGen error during edge protection";
		break;
	case ERROR_OTHERS:
	default:
		msg = "Unknown error";
		break;
	}
	if (filename && filename[0] != '\0') {
		msg += " [file: ";
		msg += filename;
		msg += "]";
	}
	throw std::runtime_error(msg);
}

#endif
