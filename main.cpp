#include <memory>
#include <random>
#include <omp.h>
#include "vec3.h"
#include "ray.h"
#include "hit.h"
#include "camera.h"
#include "image.h"
#include "aggregate.h"
#include "primitive.h"
#include "shape.h"
#include "material.h"
#include "light.h"
#include "random.h"
#include "sky.h"
#include "texture.h"
#include <glob.h>

// g++ --std=c++11 -Xpreprocessor -fopenmp -O3 -mtune=native -march=native -lomp main.cpp

const int MAX_DEPTH = 15; //最大反射回数
const double ROULETTE = 0.99; //ロシアンルーレットの確率

//rayの方向から来る放射輝度の値を計算して返す
//init_ray: 最初のレイ
//aggregate: 物体集合
//sky：空の色を表現するSkyクラス
Vec3 radiance(Ray& ray, const Aggregate& aggregate, const Sky& sky) {
  Vec3 col; //最終的に返す色
  Vec3 throughput(1); //スループット
  //Ray ray = init_ray; //更新していくレイ

  //級数の評価
  for(int depth = 0; depth < MAX_DEPTH; depth++) {
    Hit res;
    //レイがシーンと衝突した場合
    if(aggregate.intersect(ray, res)) {
		ray.distance += (res.hitPos - ray.origin).length();

      //法線
      Vec3 n = res.hitNormal;
      //ローカル座標系の構築
      Vec3 s, t;
      orthonormalBasis(n, s, t);
      //出射方向をローカル座標系に変換
      Vec3 wo_local = worldToLocal(-ray.direction, n, s, t);

      //マテリアルと光源
      auto hitMaterial = res.hitPrimitive->material;
      auto hitLight = res.hitPrimitive->light;
      
	  //Leの加算
	  if(hitLight != nullptr) {
		  col += throughput*hitLight->Le(res, wo_local);
		  break;
	  }

      //方向のサンプリングとBRDFの評価
      Vec3 brdf;
      Vec3 wi_local;
      brdf = hitMaterial->sample(res, wo_local, wi_local);
      //コサイン
      double cos = absCosTheta(wi_local);
      //サンプリングされた方向をワールド座標系に変換
      Vec3 wi = localToWorld(wi_local, n, s, t);

      //スループットの更新
      throughput *= brdf*cos;

      //次のレイを生成
      //ray = Ray(res.hitPos, wi);
	  ray.origin = res.hitPos;
	  ray.direction = wi;
    }
    //空に飛んでいった場合
    else {
      col += throughput*sky.getRadiance(ray);
	  ray.distance = std::numeric_limits<float>::infinity();;
      break;
    }

    //ロシアンルーレット
    if(rnd() >= ROULETTE) break;
    else throughput /= ROULETTE;
  }
  return col;
}

Aggregate sample_scene(){
  Aggregate aggregate;
  //shape
  auto sphere1 = std::make_shared<Sphere>(Vec3(6,2,0), 2);
  auto sphere3 = std::make_shared<Sphere>(Vec3(-3,2,-3), 2.0);
  auto sphere2 = std::make_shared<Sphere>(Vec3(-10,4,3), 4.0);
  auto large_sphere = std::make_shared<Sphere>(Vec3(0,-1000,0), 1000);
  auto sphere1_trans = std::make_shared<Sphere>(Vec3(0, 2, 0), 2.0);
  auto xy_rect = std::make_shared<XY_rect>(2, 6, 0, 4, -8);
  auto wall1 = std::make_shared<Translate>(std::make_shared<Rotate_y>(std::make_shared<XY_rect>(-3, 20, 0, 15, -8), 80), Vec3(-10, 0, 10)); //奥
  auto wall2 = std::make_shared<Translate>(std::make_shared<Rotate_y>(std::make_shared<XY_rect>(-30, 15, 0, 15, -8), 5), Vec3(0, 0, 13)); //右側

  //material
  auto blue = std::make_shared<Lambertian>(std::make_shared<Solid_color>(Vec3(0.2, 0.2, 0.8)));
  auto red = std::make_shared<Lambertian>(std::make_shared<Solid_color>(Vec3(0.8, 0.2, 0.2)));
  auto green = std::make_shared<Lambertian>(std::make_shared<Solid_color>(Vec3(0.2, 0.8, 0.2)));
  auto white = std::make_shared<Lambertian>(std::make_shared<Solid_color>(Vec3(0.8, 0.8, 0.8)));
  auto light = std::make_shared<Diffuse_light>(std::make_shared<Solid_color>(Vec3(7)));

  //primitive
  auto ground = std::make_shared<Primitive>(large_sphere, white);
  auto red_sphere = std::make_shared<Primitive>(sphere1, red);
  auto blue_sphere = std::make_shared<Primitive>(sphere3, blue);
  auto green_sphere = std::make_shared<Primitive>(sphere2, green);
  //auto light_sphere = std::make_shared<Primitive>(sphere1, light1);
  auto rect = std::make_shared<Primitive>(xy_rect, light);
  auto white_wall1 = std::make_shared<Primitive>(wall1, white);
  auto white_wall2 = std::make_shared<Primitive>(wall2, white);

  //add
  aggregate.add(red_sphere);
  aggregate.add(white_wall1);
  aggregate.add(white_wall2);
  aggregate.add(blue_sphere);
  aggregate.add(ground);
  aggregate.add(rect);
  return aggregate;
}


int main() {
  const int spp = 50000; //サンプリング数
  const int width = 1280;
  const int height = 720;
  TransientImage t_img(width, height, 20, 100, 0.8);
  Image img(width, height);

  Aggregate aggregate;
  
  aggregate = sample_scene();
  Vec3 lookfrom(30, 3, -7);
  Vec3 lookat(0, 2, -2);
  Vec3 vup(0, 1, 0);
  double vfov = 18.0;
  Camera cam(lookfrom, lookat, vup, vfov);
  UniformSky sky(Vec3(0));

  
  double t_min_real = 10000000;
  double t_max_real = 0;
  #pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < img.width; i++) {
    for(int j = 0; j < img.height; j++) {
      for(int k = 0; k < spp; k++) {
        //(u, v)の計算
        double u = (2.0*(i + rnd()) - img.width)/img.height;
        double v = (2.0*(j + rnd()) - img.height)/img.height;

        //レイを生成
        Ray ray = cam.get_ray(u, -v);

        //放射輝度を計算
        Vec3 col = radiance(ray, aggregate, sky);
        
        //サンプルを加算
		img.addPixel(i, j, col);

		//transient imagintの処理
		double t = ray.distance; //sec
		//t_min_real = std::min(t, t_min_real);
		//t_max_real = std::max(t, t_max_real);
		//t = img.time_min;
		if(t < t_img.time_min or t_img.time_max<t) continue;

        //サンプルを加算
        t_img.addPixel(i, j, t, col);
      }

      //進捗状況を出力
      if(omp_get_thread_num() == 0) {
        std::cout << double(j + i*img.height)/(img.width*img.height) * 100 << "\r" << std::flush;
      }
    }
  }

  std::cout << "scene t_min: " << t_min_real <<std::endl;
  std::cout << "scene t_max: " << t_max_real <<std::endl;

  //サンプリング数で割る
  img.divide(spp);
  t_img.divide(spp);

  //ガンマ補正
  img.gamma_correction();
  t_img.gamma_correction();

  //PPM出力
  img.ppm_output("transient/main.ppm");
  t_img.ppm_output("transient/transient.ppm");

  return 0;
}
