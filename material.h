#ifndef MATERIAL_H
#define MATERIAL_H

#include "vec3.h"
#include "ray.h"
#include "texture.h"
#include "hit.h"
#include "random.h"

//反射ベクトルを返す
//v：反射させるベクトル
//u：法線ベクトル
inline Vec3 reflect(const Vec3& v, const Vec3& n){
	return -v + 2*dot(v, n)*n;
}

//ローカル座標系のベクトルからcosThetaを計算する
//戻り値：コサイン
//v：ローカル座標系のベクトル
inline double cosTheta(const Vec3& v){
	return v.y;
}

//ローカル座標系のベクトルからcosThetaの絶対値を計算する
//戻り値: コサインの絶対値
//v: ローカル座標系のベクトル
inline double absCosTheta(const Vec3& v) {
  return std::abs(v.y);
}

//フレネル反射率を返す
//v: 入射ベクトル
//n: 法線
//n1: 入射側屈折率
//n2: 出射側屈折率
inline double fresnel(const Vec3& v, const Vec3& n, double n1, double n2) {
  double f0 = std::pow((n1 - n2)/(n1 + n2), 2.0);
  double cos = absCosTheta(v);
  return f0 + (1 - f0)*std::pow(1 - cos, 5.0);
}

//屈折方向を計算する
//全反射が起きた場合はfalseを返す
//v: 入射ベクトル
//r: 屈折ベクトル
//n: 法線
//n1: 入射側屈折率
//n2: 出射側屈折率
inline bool refract(const Vec3& v, Vec3& r, const Vec3& n, double n1, double n2) {
  double cos = absCosTheta(v);
  double sin = std::sqrt(std::max(1 - cos*cos, 0.0));
  double alpha = n1/n2 * sin;

  //全反射
  if(alpha*alpha > 1.0) return false;

  //屈折ベクトル
  r = n1/n2 * (-v + dot(v, n)*n) - std::sqrt(1 - alpha*alpha)*n;

  return true;
}



class Material {
    public:
        virtual RGB sample(
            const Hit& res, const Vec3& wo, Vec3& wi
        ) const {
            return RGB(0,0,0);
        }
};

class Lambertian : public Material {
	public:
		std::shared_ptr<Texture> rho; //反射率

		Lambertian(std::shared_ptr<Texture> _rho) : rho(_rho) {};

		RGB sample(const Hit& res, const Vec3& wo, Vec3& wi)
			const{
				//一様乱数
				double u = rnd();
				double v = rnd();

				//(theta, phi)の計算
				double theta = 0.5 * std::acos(1-2*u);
				double phi = 2*M_PI * v;

				//(x, y, z)の計算
				double x = std::cos(phi)*std::sin(theta);
				double y = std::cos(theta);
				double z = std::sin(phi)*std::sin(theta);

				//サンプリングされた方向
				wi = Vec3(x, y, z);

				//確率密度関数の値
				double pdf = std::cos(theta)/M_PI;

				return rho->value(res.u, res.v, res.hitPos)/M_PI/pdf;
		};
};

class Mirror : public Material {
	public:
		std::shared_ptr<Texture> rho; //反射率
		Mirror(std::shared_ptr<Texture> _rho) : rho(_rho) {};
		
		RGB sample(const Hit& res, const Vec3& wo, Vec3& wi) const {
			//サンプリング方向は反射ベクトルのみ
			wi = reflect(wo, Vec3(0, 1, 0));

			//デルタ関数の部分は打ち消される
			double pdf = 1.0;
			return 1/cosTheta(wi)/pdf * rho->value(res.u, res.v, res.hitPos);
		};
};

class Glass : public Material {
  public:
    double ior; //屈折率

    Glass(double _ior) : ior(_ior) {};

    RGB sample(const Hit& res, const Vec3& wo, Vec3& wi) const {
      //ガラスに入射しているかのフラグ
      bool isEntering = cosTheta(wo) > 0;

      //入射と出射で屈折率、法線を入れ替える
      double n1; //入射側屈折率
      double n2; //出射側屈折率
      Vec3 normal; //法線 入射と出射で入れ替える必要がある
      if(isEntering) {
        n1 = 1.0;
        n2 = ior;
        normal = Vec3(0, 1, 0);
      }
      else {
        n1 = ior;
        n2 = 1.0;
        normal = Vec3(0, -1, 0);
      }

      //フレネル反射率
      double fr = fresnel(wo, normal, n1, n2);

      //反射
      if(rnd() < fr) {
        wi = reflect(wo, normal);
        double pdf = fr;
        return fr/absCosTheta(wi) * Vec3(1) / pdf;
      }
      //屈折
      else {
        if(refract(wo, wi, normal, n1, n2)) {
          double pdf = 1 - fr;
          return std::pow(n1/n2, 2.0) * (1 - fr)/absCosTheta(wi) * Vec3(1) / pdf;
        }
        //全反射
        else {
          wi = reflect(wo, normal);
          double pdf = 1 - fr;
          return (1 - fr)/absCosTheta(wi) * Vec3(1) / pdf;
        }
      }
    };
};

#endif
