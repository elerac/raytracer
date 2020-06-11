#ifndef RAY_H
#define RAY_H
#include "vec3.h"

class Ray{
	public:
		Vec3 origin;
		Vec3 direction;
		double tmin = 0.001;
		double tmax = 10000;
		double distance = 0.0;

	Ray(const Vec3& _origin, const Vec3& _direction) : origin(_origin), direction(_direction) {};

	Vec3 operator()(double t) const{
		return origin + t*direction;
	};

};

#endif
