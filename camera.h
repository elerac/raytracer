#ifndef CAMERA_H
#define CAMERA_H

#include "vec3.h"
#include "ray.h"

class Camera {
    public:
        Camera() {};
        
        Camera(
			Vec3 lookfrom, //カメラ位置
			Vec3 lookat, //視線方向
			Vec3 vup,  //センサーの上方向（通常(0, 1, 0)）
			double vfov //垂直方向の画角
		) {
			double theta = vfov*M_PI/180;
            dist_origin2sensor = 1.0/tan(theta/2);
			forward = normalize(lookat-lookfrom);
			right = normalize(cross(forward, vup));
            up = normalize(cross(right, forward));
            //orthonormalBasis(forward, right, up);
            //up *= -1;
            origin = lookfrom;
        }
		
        Ray get_ray(double u, double v) const {
			return Ray(origin, normalize(right*u + up*v + forward*dist_origin2sensor));
		}

    private:
        Vec3 origin;
        Vec3 forward;
        Vec3 right;
        Vec3 up;
		double dist_origin2sensor;
};

#endif
