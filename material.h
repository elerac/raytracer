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



#include "math.h"
#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

#define merl_attenution 10
#define RED_SCALE (1.0/1500.0)*merl_attenution
#define GREEN_SCALE (1.15/1500.0)*merl_attenution
#define BLUE_SCALE (1.66/1500.0)*merl_attenution

// cross product of two vectors
void cross_product (double* v1, double* v2, double* out)
{
	out[0] = v1[1]*v2[2] - v1[2]*v2[1];
	out[1] = v1[2]*v2[0] - v1[0]*v2[2];
	out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

// normalize vector
void _normalize(double* v)
{
	// normalize
	double len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] = v[0] / len;
	v[1] = v[1] / len;
	v[2] = v[2] / len;
}

// rotate vector along one axis
void rotate_vector(double* vector, double* axis, double angle, double* out)
{
	double temp;
	double cross[3];
	double cos_ang = cos(angle);
	double sin_ang = sin(angle);

	out[0] = vector[0] * cos_ang;
	out[1] = vector[1] * cos_ang;
	out[2] = vector[2] * cos_ang;

	temp = axis[0]*vector[0]+axis[1]*vector[1]+axis[2]*vector[2];
	temp = temp*(1.0-cos_ang);

	out[0] += axis[0] * temp;
	out[1] += axis[1] * temp;
	out[2] += axis[2] * temp;

	cross_product (axis,vector,cross);
	
	out[0] += cross[0] * sin_ang;
	out[1] += cross[1] * sin_ang;
	out[2] += cross[2] * sin_ang;
}


// convert standard coordinates to half vector/difference vector coordinates
void std_coords_to_half_diff_coords(double theta_in, double fi_in, double theta_out, double fi_out,
								double& theta_half,double& fi_half,double& theta_diff,double& fi_diff )
{

	// compute in vector
	double in_vec_z = cos(theta_in);
	double proj_in_vec = sin(theta_in);
	double in_vec_x = proj_in_vec*cos(fi_in);
	double in_vec_y = proj_in_vec*sin(fi_in);
	double in[3]= {in_vec_x,in_vec_y,in_vec_z};
	_normalize(in);


	// compute out vector
	double out_vec_z = cos(theta_out);
	double proj_out_vec = sin(theta_out);
	double out_vec_x = proj_out_vec*cos(fi_out);
	double out_vec_y = proj_out_vec*sin(fi_out);
	double out[3]= {out_vec_x,out_vec_y,out_vec_z};
	_normalize(out);


	// compute halfway vector
	double half_x = (in_vec_x + out_vec_x)/2.0f;
	double half_y = (in_vec_y + out_vec_y)/2.0f;
	double half_z = (in_vec_z + out_vec_z)/2.0f;
	double half[3] = {half_x,half_y,half_z};
	_normalize(half);

	// compute  theta_half, fi_half
	theta_half = acos(half[2]);
	fi_half = atan2(half[1], half[0]);


	double bi_normal[3] = {0.0, 1.0, 0.0};
	double normal[3] = { 0.0, 0.0, 1.0 };
	double temp[3];
	double diff[3];

	// compute diff vector
	rotate_vector(in, normal , -fi_half, temp);
	rotate_vector(temp, bi_normal, -theta_half, diff);
	
	// compute  theta_diff, fi_diff	
	theta_diff = acos(diff[2]);
	fi_diff = atan2(diff[1], diff[0]);

}

// Lookup theta_half index
// This is a non-linear mapping!
// In:  [0 .. pi/2]
// Out: [0 .. 89]
inline int theta_half_index(double theta_half)
{
	if (theta_half <= 0.0)
		return 0;
	double theta_half_deg = ((theta_half / (M_PI/2.0))*BRDF_SAMPLING_RES_THETA_H);
	double temp = theta_half_deg*BRDF_SAMPLING_RES_THETA_H;
	temp = sqrt(temp);
	int ret_val = (int)temp;
	if (ret_val < 0) ret_val = 0;
	if (ret_val >= BRDF_SAMPLING_RES_THETA_H)
		ret_val = BRDF_SAMPLING_RES_THETA_H-1;
	return ret_val;
}

// Lookup theta_diff index
// In:  [0 .. pi/2]
// Out: [0 .. 89]
inline int theta_diff_index(double theta_diff)
{
	int tmp = int(theta_diff / (M_PI * 0.5) * BRDF_SAMPLING_RES_THETA_D);
	if (tmp < 0)
		return 0;
	else if (tmp < BRDF_SAMPLING_RES_THETA_D - 1)
		return tmp;
	else
		return BRDF_SAMPLING_RES_THETA_D - 1;
}

// Lookup phi_diff index
inline int phi_diff_index(double phi_diff)
{
	// Because of reciprocity, the BRDF is unchanged under
	// phi_diff -> phi_diff + M_PI
	if (phi_diff < 0.0)
		phi_diff += M_PI;

	// In: phi_diff in [0 .. pi]
	// Out: tmp in [0 .. 179]
	int tmp = int(phi_diff / M_PI * BRDF_SAMPLING_RES_PHI_D / 2);
	if (tmp < 0)	
		return 0;
	else if (tmp < BRDF_SAMPLING_RES_PHI_D / 2 - 1)
		return tmp;
	else
		return BRDF_SAMPLING_RES_PHI_D / 2 - 1;
}

// Given a pair of incoming/outgoing angles, look up the BRDF.
void lookup_brdf_val(double* brdf, double theta_in, double fi_in,
			  double theta_out, double fi_out, 
			  double& red_val,double& green_val,double& blue_val)
{
	// Convert to halfangle / difference angle coordinates
	double theta_half, fi_half, theta_diff, fi_diff;
	
	std_coords_to_half_diff_coords(theta_in, fi_in, theta_out, fi_out,
		       theta_half, fi_half, theta_diff, fi_diff);


	// Find index.
	// Note that phi_half is ignored, since isotropic BRDFs are assumed
	int ind = phi_diff_index(fi_diff) +
		  theta_diff_index(theta_diff) * BRDF_SAMPLING_RES_PHI_D / 2 +
		  theta_half_index(theta_half) * BRDF_SAMPLING_RES_PHI_D / 2 *
					         BRDF_SAMPLING_RES_THETA_D;

	red_val = brdf[ind] * RED_SCALE;
	green_val = brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2] * GREEN_SCALE;
	blue_val = brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D] * BLUE_SCALE;

	
	/*if (red_val < 0.0 || green_val < 0.0 || blue_val < 0.0)
		fprintf(stderr, "Below horizon.\n");
*/
}

// Read BRDF data
bool read_brdf(const char *filename, double* &brdf)
{
	FILE *f = fopen(filename, "rb");
	if (!f)
		return false;

	int dims[3];
	fread(dims, sizeof(int), 3, f);
	int n = dims[0] * dims[1] * dims[2];
	if (n != BRDF_SAMPLING_RES_THETA_H *
		 BRDF_SAMPLING_RES_THETA_D *
		 BRDF_SAMPLING_RES_PHI_D / 2) 
	{
		fprintf(stderr, "Dimensions don't match\n");
		fclose(f);
		return false;
	}

	brdf = (double*) malloc (sizeof(double)*3*n);
	fread(brdf, sizeof(double), 3*n, f);

	fclose(f);
	return true;
}

class MERL : public Material {
	public:
		double* brdf;
		
		MERL(const char* filename) {
			read_brdf(filename, brdf);
		}

		RGB sample(const Hit& res, const Vec3& wo, Vec3& wi) const {
			//一様乱数
			double u = rnd();
			double v = rnd();
			/*
			//(theta, phi)の計算
			//double theta_in = 0.5 * std::acos(1-2*u);
			double theta_in = std::acos(std::sqrt(1-u));
			double phi_in = 2*M_PI * v;

			//確率密度関数の値
			double pdf = std::cos(theta_in)/M_PI;
			*/
			

			//(theta, phi)の計算
			double theta_in = 0.5*M_PI * u;
			double phi_in = 2*M_PI * v;

			//確率密度関数の値
			double pdf = 1/(2*M_PI);
			

			//(x, y, z)の計算
			double x = std::cos(phi_in)*std::sin(theta_in);
			double y = std::cos(theta_in);
			double z = std::sin(phi_in)*std::sin(theta_in);

			//サンプリングされた方向
			wi = Vec3(x, y, z);

        	double theta_out = std::acos(absCosTheta(wo));
        	Vec3 wi_xz(wi.x, 0, wi.z);
        	Vec3 wo_xz(wo.x, 0, wo.z);
        	double phi_out = phi_in + std::acos(dot(normalize(wi_xz), normalize(wo_xz)));

        	double r, g, b;
        	lookup_brdf_val(brdf, theta_in, phi_in, theta_out, phi_out, r, g, b);
        	Vec3 rho(r, g, b);

			return rho/M_PI/pdf;
		};
};

#endif
