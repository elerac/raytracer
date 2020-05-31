#ifndef VEC3_H
#define VEC3_H
#include <iostream>
#include <cmath>

class Vec3{
	public:
		double x;
		double y;
		double z;

	Vec3(){ x=y=z=0; };
	Vec3(double _x){ x=y=z=_x;};
	Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

	Vec3 operator-() const{
		return Vec3(-x, -y, -z);
	};

	Vec3& operator+=(const Vec3& v){
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	};
	
	Vec3& operator-=(const Vec3& v){
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	};
	
	Vec3& operator*=(const Vec3& v){
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	};
	
	Vec3& operator/=(const Vec3& v){
		x /= v.x;
		y /= v.y;
		z /= v.z;
		return *this;
	};

	//ベクトルの長さ
	double length() const{
		return std::sqrt(x*x + y*y + z*z);
	};

	//ベクトルの長さの二乗
	double length2() const{
		return x*x + y*y + z*z;
	};
};
typedef Vec3 RGB;

//ベクトル同士の演算
Vec3 operator+(const Vec3& v1, const Vec3& v2){
	return Vec3(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}

Vec3 operator-(const Vec3& v1, const Vec3& v2){
	return Vec3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}

Vec3 operator*(const Vec3& v1, const Vec3& v2){
	return Vec3(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);
}

Vec3 operator/(const Vec3& v1, const Vec3& v2){
	return Vec3(v1.x/v2.x, v1.y/v2.y, v1.z/v2.z);
}


//ベクトルとスカラーの演算
Vec3 operator+(const Vec3& v, double k){
	return Vec3(v.x+k, v.y+k, v.z+k);
}

Vec3 operator+(double k, const Vec3& v){
	return v+k;
}

Vec3 operator-(const Vec3& v, double k){
	return Vec3(v.x-k, v.y-k, v.z-k);
}

Vec3 operator-(double k, const Vec3& v){
	return Vec3(k-v.x, k-v.y, k-v.z);
}

Vec3 operator*(const Vec3& v, double k){
	return Vec3(v.x*k, v.y*k, v.z*k);
}

Vec3 operator*(double k, const Vec3& v){
	return v*k;
}

Vec3 operator/(const Vec3& v, double k){
	return Vec3(v.x/k, v.y/k, v.z/k);
}

Vec3 operator/(double k, const Vec3& v){
	return Vec3(k/v.x, k/v.y, k/v.z);
}

//内積
inline double dot(const Vec3& v1, const Vec3& v2){
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//外積
inline Vec3 cross(const Vec3& v1, const Vec3& v2){
	return Vec3(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}

//正規化
Vec3 normalize(const Vec3& v){
	return v/v.length();
}

//コンソール出力
std::ostream& operator<<(std::ostream& stream, const Vec3& v){
	stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return stream;
}

//正規直交基底を作り出す
inline void orthonormalBasis(const Vec3& v1, Vec3& v2, Vec3& v3){
	if(std::abs(v1.x) > 0.9)v2 = Vec3(0, 1, 0);
	else v2 = Vec3(1, 0, 0);

	v2 = normalize(v2 - dot(v1, v2)*v1);
	v3 = cross(v1, v2);
}

//ワールド座標系からローカル座標系に変換する
inline Vec3 worldToLocal(const Vec3& w, const Vec3& n, const Vec3& s, const Vec3& t) {
    return Vec3(dot(s, w), dot(n, w), dot(t, w));
}

//ローカル座標系からワールド座標系に変換する
inline Vec3 localToWorld(const Vec3& w, const Vec3& n, const Vec3& s, const Vec3& t) {
    return Vec3(s.x*w.x + n.x*w.y + t.x*w.z, s.y*w.x + n.y*w.y + t.y*w.z, s.z*w.x + n.z*w.y + t.z*w.z);
}


//

#endif
