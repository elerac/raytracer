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

//g++ --std=c++11 -Xpreprocessor -fopenmp -O3 -mtune=native -march=native -lomp main.cpp

const int MAX_DEPTH = 10; //最大反射回数
const double ROULETTE = 0.95; //ロシアンルーレットの確率

//rayの方向から来る放射輝度の値を計算して返す
//init_ray: 最初のレイ
//aggregate: 物体集合
//sky：空の色を表現するSkyクラス
Vec3 radiance(const Ray& init_ray, const Aggregate& aggregate, const Sky& sky) {
  Vec3 col; //最終的に返す色
  Vec3 throughput(1); //スループット
  Ray ray = init_ray; //更新していくレイ

  //級数の評価
  for(int depth = 0; depth < MAX_DEPTH; depth++) {
    Hit res;
    //レイがシーンと衝突した場合
    if(aggregate.intersect(ray, res)) {
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
      ray = Ray(res.hitPos, wi);
    }
    //空に飛んでいった場合
    else {
      col += throughput*sky.getRadiance(ray);
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
  auto sphere1 = std::make_shared<Sphere>(Vec3(0,0,0), 1.0);
  auto sphere2 = std::make_shared<Sphere>(Vec3(2,0,0), 1.0);
  auto sphere3 = std::make_shared<Sphere>(Vec3(-2,0,0), 1.0);
  auto large_sphere = std::make_shared<Sphere>(Vec3(0,-10001,0), 10000);
  auto sphere1_trans = std::make_shared<Translate>(sphere1, Vec3(1, 1, -1));
  auto xy_rect = std::make_shared<Rotate_y>(std::make_shared<XY_rect>(-2, 2, -2, 2, 0), 45);
  auto box = std::make_shared<Translate>(std::make_shared<Box>(Vec3(-1, -1, -1), Vec3(1, 1, 1)), Vec3(0,3,0));
  
  //material
  auto earth_surface = std::make_shared<Lambertian>(std::make_shared<Image_texture>("earthmap.jpg"));
  auto checker = std::make_shared<Lambertian>(std::make_shared<Checker_texture>(std::make_shared<Solid_color>(Vec3(1)), std::make_shared<Solid_color>(Vec3(1))));
  auto blue = std::make_shared<Lambertian>(std::make_shared<Solid_color>(Vec3(0.2, 0.2, 0.8)));
  auto red = std::make_shared<Lambertian>(std::make_shared<Solid_color>(Vec3(0.8, 0.2, 0.2)));
  auto mirror = std::make_shared<Mirror>(std::make_shared<Solid_color>(Vec3(0.9)));
  auto light1 = std::make_shared<Diffuse_light>(std::make_shared<Solid_color>(Vec3(1)));
  auto glass = std::make_shared<Glass>(1.5);
  auto merl_gold = std::make_shared<MERL>("merl_brdf/gold-metallic-paint.binary");
  auto merl_yellow = std::make_shared<MERL>("merl_brdf/yellow-plastic.binary");
  auto merl_steel = std::make_shared<MERL>("merl_brdf/steel.binary");

  //primitive
  auto ground = std::make_shared<Primitive>(large_sphere, checker);
  auto red_sphere = std::make_shared<Primitive>(sphere3, red);
  auto earth_sphere = std::make_shared<Primitive>(sphere2, earth_surface);
  auto mirror_sphere = std::make_shared<Primitive>(sphere1, blue);
  auto blue_box = std::make_shared<Primitive>(box, earth_surface);
  auto light_sphere = std::make_shared<Primitive>(sphere3, light1);
  auto rect = std::make_shared<Primitive>(xy_rect, blue);
  auto glass_sphere = std::make_shared<Primitive>(sphere1, glass);
  auto merl_gold_sphere = std::make_shared<Primitive>(sphere1, merl_gold);
  auto merl_yellow_sphere = std::make_shared<Primitive>(sphere2, merl_yellow);
  auto merl_steel_sphere = std::make_shared<Primitive>(sphere3, merl_steel);

  //add
  //aggregate.add(red_sphere);
  //aggregate.add(earth_sphere);
  //aggregate.add(mirror_sphere);
  aggregate.add(merl_gold_sphere);
  aggregate.add(merl_yellow_sphere);
  aggregate.add(merl_steel_sphere);
  //aggregate.add(glass_sphere);
  //aggregate.add(blue_box);
  //aggregate.add(light_sphere);
  aggregate.add(ground);
  //aggregate.add(rect);
  return aggregate;
}

Aggregate my_cornell_box() {
  Aggregate world;
  auto red = std::make_shared<Lambertian>(std::make_shared<Solid_color>(.65, .05, .05));
  auto white = std::make_shared<Lambertian>(std::make_shared<Solid_color>(.73, .73, .73));
  auto green = std::make_shared<Lambertian>(std::make_shared<Solid_color>(.12, .45, .15));
  auto light = std::make_shared<Diffuse_light>(std::make_shared<Solid_color>(Vec3(50)));
  auto checker = std::make_shared<Lambertian>(std::make_shared<Checker_texture>(std::make_shared<Solid_color>(Vec3(0.73)), std::make_shared<Solid_color>(Vec3(0.1))));
  auto glass = std::make_shared<Glass>(1.5);
  auto mirror = std::make_shared<Mirror>(std::make_shared<Solid_color>(Vec3(0.9)));

  world.add(std::make_shared<Primitive>(std::make_shared<YZ_rect>(0, 555, 0, 555, 555), green));
  world.add(std::make_shared<Primitive>(std::make_shared<YZ_rect>(0, 555, 0, 555, 0), red));
  world.add(std::make_shared<Primitive>(std::make_shared<XZ_rect>(213, 343, 227, 332, 554), light));
  world.add(std::make_shared<Primitive>(std::make_shared<XZ_rect>(0, 555, 0, 555, 555), white));
  world.add(std::make_shared<Primitive>(std::make_shared<XZ_rect>(0, 555, 0, 555, 0), checker));
  world.add(std::make_shared<Primitive>(std::make_shared<XY_rect>(0, 555, 0, 555, 555), white));

  std::shared_ptr<Shape> sphere1 = std::make_shared<Sphere>(Vec3(390, 120, 350), 120);
  world.add(std::make_shared<Primitive>(sphere1, mirror));

  std::shared_ptr<Shape> sphere2 = std::make_shared<Sphere>(Vec3(150, 100, 170), 100);
  world.add(std::make_shared<Primitive>(sphere2, glass));

  return world;
}

void scene_merl(Camera &cam, Aggregate &aggregate) {
  //merl brdfのファイル名を全て読み込み
  glob_t globbuf;
  std::vector<std::string> files;
  glob("merl_brdf/*.binary", 0, NULL, &globbuf);
  for (int i = 0; i < globbuf.gl_pathc; i++) {
    files.push_back(globbuf.gl_pathv[i]);
    }
  globfree(&globbuf);

  // シャッフル
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::shuffle(files.begin(), files.end(), engine);

  double space = 3;
  int width= 9;
  int height = 6;
  for(int x=0; x<width; x++) {
    for(int z=0; z<height; z++) {
      auto merl = std::make_shared<MERL>(files[x*height+z].c_str());
      auto sphere = std::make_shared<Sphere>(Vec3(x*space, 0, z*space), 1);
      aggregate.add(std::make_shared<Primitive>(sphere, merl));
    }
  }

  auto checker = std::make_shared<Lambertian>(std::make_shared<Checker_texture>(std::make_shared<Solid_color>(Vec3(0.8)), std::make_shared<Solid_color>(Vec3(0.1))));
  auto ground = std::make_shared<XZ_rect>(-space/2, (width-0.5)*space, -space/2, (height-0.5)*space, -1);
  //auto ground = std::make_shared<XZ_rect>(-space/2*10, (width-0.5)*space*10, -space/2*10, (height-0.5)*space*10, -1);
  aggregate.add(std::make_shared<Primitive>(ground, checker));

  double x_center = (double)(width-1)*space/2;
  double z_center = (double)(height-1)*space/2;
  Vec3 lookfrom(x_center, 24, -14);
  Vec3 lookat(x_center, 0, z_center*0.7);
  Vec3 vup(0, 1, 0);
  double vfov = 32.0;
  cam = {lookfrom, lookat, vup, vfov};
}


int main() {
  const int spp = 500; //サンプリング数

  Image img(1280, 720);

  Aggregate aggregate;
  Camera cam;
  scene_merl(cam, aggregate);
  //IBL sky("envmap/PaperMill_E_3k.hdr");
  IBL sky("envmap/skylit_garage_8k.hdr");
  //SimpleSky sky; 
  /*
  aggregate = sample_scene();
  Vec3 lookfrom(0, 0, 3);
  Vec3 lookat(0, 0, 0);
  Vec3 vup(0, 1, 0);
  double vfov = 80.0;
  //Camera cam(lookfrom, lookat, vup, vfov);
  IBL sky("envmap/PaperMill_E_3k.hdr");*/

  /*
  aggregate = my_cornell_box();
  Vec3 lookfrom(278, 278, -800);
  Vec3 lookat(278, 278, 0);
  Vec3 vup(0, 1, 0);
  double vfov = 40.0;
  Camera cam(lookfrom, lookat, vup, vfov);
  UniformSky sky(Vec3(0));
  */
  //IBL sky("veranda/veranda_4k.hdr");
 

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
      }

      //進捗状況を出力
      if(omp_get_thread_num() == 0) {
        std::cout << double(j + i*img.height)/(img.width*img.height) * 100 << "\r" << std::flush;
      }
    }
  }

  //サンプリング数で割る
  img.divide(spp);

  //ガンマ補正
  img.gamma_correction();

  //PPM出力
  img.ppm_output("main.ppm");

  return 0;
}
