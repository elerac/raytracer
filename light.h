#ifndef LIGHT_H
#define LIGHT_H
#include "vec3.h"
#include "hit.h"

class Light {
	public:
		virtual RGB Le(const Hit& res, const Vec3& wo) const = 0;
};

class Diffuse_light : public Light {
	public:
		std::shared_ptr<Texture> emit;

		Diffuse_light(const std::shared_ptr<Texture> _emit) : emit(_emit) {};

		Vec3 Le(const Hit& res, const Vec3& wo)
			const {
				return emit->value(res.u, res.v, res.hitPos);
			};
};

#endif
