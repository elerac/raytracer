#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>

#include <functional>
#include <random>

template <typename T>
T clamp(T x, T xmin, T xmax){
	if(x<xmin)return xmin;
	else if(x>xmax)return xmax;
	else return x;
}

// Constants

const double infinity = std::numeric_limits<double>::infinity();

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180;
}



#endif
