#include <memory>
#include <random>
//#include <omp.h>
#include "vec3.h"
#include "ray.h"
#include "hit.h"
#include "sphere.h"
#include "camera.h"
#include "image.h"
#include "aggregate.h"
#include "material.h"
#include "light.h"
#include "random.h"
#include "sky.h"



const int MAX_DEPTH = 25; //最大反射回数
const double ROULETTE = 0.9; //ロシアンルーレットの確率


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
      auto hitMaterial = res.hitSphere->material;
      auto hitLight = res.hitSphere->light;

      //Leの加算
      col += throughput*hitLight->Le();

      //方向のサンプリングとBRDFの評価
      Vec3 brdf;
      Vec3 wi_local;
      double pdf;
      brdf = hitMaterial->sample(wo_local, wi_local, pdf);
      //コサイン
      double cos = cosTheta(wi_local);
      //サンプリングされた方向をワールド座標系に変換
      Vec3 wi = localToWorld(wi_local, n, s, t);

      //スループットの更新
      throughput *= brdf*cos/pdf;

      //次のレイを生成
      ray = Ray(res.hitPos + 0.001*res.hitNormal, wi);
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


int main() {
  const int spp = 20; //サンプリング数

  Image img(300, 200);

  Vec3 lookfrom(0, 0, 0);
  Vec3 lookat(0, 0, -3);
  Vec3 vup(0, 1, 0);
  double vfov = 80.0;
  Camera cam(lookfrom, lookat, vup, vfov);

  auto mat1 = std::make_shared<Diffuse>(Vec3(0.9)); //緑色
  auto mat2 = std::make_shared<Diffuse>(Vec3(0.2, 0.2, 0.8)); //青色
  auto mat3 = std::make_shared<Diffuse>(Vec3(0.8, 0.2, 0.2)); //赤色
  auto mirror = std::make_shared<Mirror>();

  auto light1 = std::make_shared<Light>(Vec3(0));

  Aggregate aggregate;
  aggregate.add(std::make_shared<Sphere>(Vec3(0, -100001, 0), 100000, mat1, light1)); //床
  aggregate.add(std::make_shared<Sphere>(Vec3(0, 0, -3), 1, mirror, light1)); //球
  aggregate.add(std::make_shared<Sphere>(Vec3(2, 0, -3), 1, mat2, light1)); //球
  aggregate.add(std::make_shared<Sphere>(Vec3(-2, 0, -3), 1, mat3, light1)); //球

  IBL sky("PaperMill_Ruins_E/PaperMill_E_3k.hdr");
  //SimpleSky sky;

//#pragma omp parallel for schedule(dynamic, 1)
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
      //if(omp_get_thread_num() == 0) {
        std::cout << double(j + i*img.height)/(img.width*img.height) * 100 << "\r" << std::flush;
      //}
    }
  }

  //サンプリング数で割る
  img.divide(spp);

  //ガンマ補正
  img.gamma_correction();

  //PPM出力
  img.ppm_output("use_ibl.ppm");

  return 0;
}
