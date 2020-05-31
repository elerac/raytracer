#ifndef SHAPE_H
#define SHAPE_H
#include "vec3.h"
#include "ray.h"
#include "hit.h"
#include "primitive.h"

class Shape {
	public:
		virtual bool intersect(const Ray& ray, Hit& res) const = 0;
};

class Sphere : public Shape {
	public:
		Vec3 center; //中心位置
		double radius; //半径

		Sphere(const Vec3& _center, double _radius) : center(_center), radius(_radius) {};

		//rayと交差しているか判定する．交差している場合は衝突情報をresに格納する
		bool intersect(const Ray& ray, Hit& res) const {
			//二次方程式の係数
			double b = dot(ray.direction, ray.origin-center);
			double c = (ray.origin - center).length2() - radius*radius;
			//判別式
			double D = b*b - c;

			//解の候補距離
			double t1 = -b - std::sqrt(D);
			double t2 = -b + std::sqrt(D);

			if(t1>ray.tmax | t2<ray.tmin)return false;
			double t = t1;
			if(t<ray.tmin){
				t = t2;
				if(t2>ray.tmax)return false;
			}

			//衝突情報を格納
			res.t = t;
			res.hitPos = ray(t);
			res.hitNormal = normalize(res.hitPos - center);

			double phi = std::atan2(res.hitNormal.z, res.hitNormal.x);
			if(phi < 0) phi += 2*M_PI;
			double theta = std::acos(res.hitNormal.y);
			res.u = phi/(2*M_PI);
			res.v = theta/M_PI;

			return true;
		};
};

class XY_rect : public Shape {
	public:
		double x0, x1, y0, y1, k;
		XY_rect(
            double _x0, double _x1, double _y0, double _y1, double _k
        ) : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k) {};

		bool intersect(const Ray& ray, Hit& res) const {
			double t = (k-ray.origin.z) / ray.direction.z;
			if (t < ray.tmin || t > ray.tmax)
				return false;

			double x = ray.origin.x + t*ray.direction.x;
    		double y = ray.origin.y + t*ray.direction.y;
    		if (x < x0 || x > x1 || y < y0 || y > y1)
        		return false;

			res.u = (x-x0)/(x1-x0);
			res.v = (y-y0)/(y1-y0);
			res.t = t;
			Vec3 normal(0, 0, 1);
			if(dot(ray.direction, normal) > 0)
				normal *= -1;
			res.hitNormal = normal;
			res.hitPos = ray(t);

			return true;
		};
};

class XZ_rect : public Shape {
	public:
		double x0, x1, z0, z1, k;
		XZ_rect(
            double _x0, double _x1, double _z0, double _z1, double _k
        ) : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k) {};

		bool intersect(const Ray& ray, Hit& res) const {
			double t = (k-ray.origin.y) / ray.direction.y;
			if (t < ray.tmin || t > ray.tmax)
				return false;
				
			double x = ray.origin.x + t*ray.direction.x;
    		double z = ray.origin.z + t*ray.direction.z;
    		if (x < x0 || x > x1 || z < z0 || z > z1)
        		return false;

			res.u = (x-x0)/(x1-x0);
			res.v = (z-z0)/(z1-z0);
			res.t = t;
			Vec3 normal(0, 1, 0);
			if(dot(ray.direction, normal) > 0)
				normal *= -1;
			res.hitNormal = normal;
			res.hitPos = ray(t);

			return true;
		};
};

class YZ_rect : public Shape {
	public:
		double y0, y1, z0, z1, k;
		YZ_rect(
            double _y0, double _y1, double _z0, double _z1, double _k
        ) : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k) {};

		bool intersect(const Ray& ray, Hit& res) const {
			double t = (k-ray.origin.x) / ray.direction.x;
			if (t < ray.tmin || t > ray.tmax)
				return false;
				
			double y = ray.origin.y + t*ray.direction.y;
    		double z = ray.origin.z + t*ray.direction.z;
    		if (y < y0 || y > y1 || z < z0 || z > z1)
        		return false;

			res.u = (y-y0)/(y1-y0);
			res.v = (z-z0)/(z1-z0);
			res.t = t;
			Vec3 normal(1, 0, 0);
			if(dot(ray.direction, normal) > 0)
				normal *= -1;
			res.hitNormal = normal;
			res.hitPos = ray(t);

			return true;
		};
};

class Translate : public Shape {
	public:
		std::shared_ptr<Shape> shape;
		Vec3 offset;
		Translate(std::shared_ptr<Shape> _shape, const Vec3& _offset
		) : shape(_shape), offset(_offset) {};

		bool intersect(const Ray& ray, Hit& res) const {
			Ray moved_ray(ray.origin - offset, ray.direction);
			if(!shape->intersect(moved_ray, res))
				return false;

			res.hitPos += offset;
			
			return true;
		};
};

class Rotate_y : public Shape {
	public:
		std::shared_ptr<Shape> shape;
		double angle;
		double sin_theta, cos_theta;

		Rotate_y(std::shared_ptr<Shape> _shape, double _angle
		) : shape(_shape), angle(_angle) {
			sin_theta = std::sin(angle*M_PI/180.0);
			cos_theta = std::cos(angle*M_PI/180.0);
		};

		bool intersect(const Ray& ray, Hit& res) const {
			Vec3 origin = ray.origin;
			Vec3 direction = ray.direction;

			origin.x = cos_theta*ray.origin.x - sin_theta*ray.origin.z;
			origin.z = sin_theta*ray.origin.x + cos_theta*ray.origin.z;

			direction.x = cos_theta*ray.direction.x - sin_theta*ray.direction.z;
			direction.z = sin_theta*ray.direction.x + cos_theta*ray.direction.z;

			Ray rotated_ray(origin, direction);

			if(!shape->intersect(rotated_ray, res))
				return false;
			
			Vec3 hitPos = res.hitPos;
			Vec3 hitNormal = res.hitNormal;

			hitPos.x = cos_theta*res.hitPos.x + sin_theta*res.hitPos.z;
			hitPos.z = -sin_theta*res.hitPos.x + cos_theta*res.hitPos.z;

			hitNormal.x = cos_theta*res.hitNormal.x + sin_theta*res.hitNormal.z;
			hitNormal.z = -sin_theta*res.hitNormal.x + cos_theta*res.hitNormal.z;

			res.hitPos = hitPos;
			res.hitNormal = hitNormal;

			return true;
		};

};

class Box : public Shape {
	public:
		Vec3 p0, p1;

		Box(const Vec3& _p0, const Vec3& _p1
		) : p0(_p0), p1(_p1) {
			add(std::make_shared<XY_rect>(p0.x, p1.x, p0.y, p1.y, p1.z));
			add(std::make_shared<XY_rect>(p0.x, p1.x, p0.y, p1.y, p0.z));
			add(std::make_shared<XZ_rect>(p0.x, p1.x, p0.z, p1.z, p1.y));
			add(std::make_shared<XZ_rect>(p0.x, p1.x, p0.z, p1.z, p0.y));
			add(std::make_shared<YZ_rect>(p0.y, p1.y, p0.z, p1.z, p1.x));
			add(std::make_shared<YZ_rect>(p0.y, p1.y, p0.z, p1.z, p0.x));
		};

		bool intersect(const Ray& ray, Hit& res) const {
			bool hit = false;
			for(auto side : sides) {
				Hit res_temp;
				if(side->intersect(ray, res_temp)) {
					if(res_temp.t < res.t) {
						hit = true;
						res = res_temp;
					}
				}
			}
			return hit;
		};

	private:
		std::vector<std::shared_ptr<Shape>> sides;
		void add(const std::shared_ptr<Shape>& shape) {
			sides.push_back(shape);
		};
};

#endif