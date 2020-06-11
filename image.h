#ifndef IMAGE_H
#define IMAGE_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "vec3.h"
#include "util.h"

class Image{
	public:
		int width; //横幅
		int height; //縦幅
		Vec3* data; //ピクセルの配列

		Image(int _width, int _height) : width(_width), height(_height){
			data = new Vec3[width*height];
		};

		~Image(){
			delete[] data;
		};

		//ピクセル(i, j)の色を取り出す
		Vec3 getPixel(int i, int j) const{
			return data[i+width*j];
		};

		//ピクセル(i, j)に色をセットする
		void setPixel(int i, int j, const Vec3& c){
			data[i+width*j] = c;
		};

		//ピクセル(i, j)に色を追加する
    	void addPixel(int i, int j, const Vec3& c) {
      		data[i+width*j] = data[i+width*j] + c;
    	};

		//全てのピクセルを一定の値で割る
		void divide(double k){
			for(int i=0; i<width; i++){
				for(int j=0; j<height; j++){
					this->setPixel(i, j, this->getPixel(i, j)/k);
				}
			}
		};

		//ガンマ補正
		void gamma_correction(){
			for(int i=0; i<width; i++){
				for(int j=0; j<height; j++){
					Vec3 c = this->getPixel(i, j);
					this->setPixel(i, j, Vec3(std::pow(c.x, 1/2.2), std::pow(c.y, 1/2.2), std::pow(c.z, 1/2.2)));
				}
			}
		};

		//PPM画像出力
		void ppm_output(const std::string& filename) const{
			std::ofstream file(filename);
			file << "P3" << std::endl;
			file << width << " " << height << std::endl;
			file << "255" << std::endl;
			for(int j=0; j<height; j++){
				for(int i=0; i<width; i++){
					Vec3 c = this->getPixel(i, j);
					int r = clamp(int(255*c.x), 0, 255);
					int g = clamp(int(255*c.y), 0, 255);
					int b = clamp(int(255*c.z), 0, 255);
					file << r << " " << g << " " << b << std::endl;
				}
			}

			file.close();
		};

};

class TransientImage{
	public:
		int width; //横幅
		int height; //縦幅
		double time_min; //最小時間
		double time_max; //最大時間
		double delta_t; //１フレームあたりの露光時間
		std::vector<std::shared_ptr<Image>> data; //ピクセルの配列
		int frame_num; //フレーム枚数

		TransientImage(int _width, int _height, double _time_min, double _time_max, double _delta_t) : width(_width), height(_height), time_min(_time_min), time_max(_time_max), delta_t(_delta_t){
			frame_num = (int)((time_max - time_min)/delta_t);
			for(int k=0; k<frame_num; k++) {
				auto img = std::make_shared<Image>(width, height);
				data.push_back(img);
			}
		};
/*
		~TransientImage(){
			delete[] data;
		};
		*/

		//ピクセル(i, j)の色を取り出す
		Vec3 getPixel(int i, int j, int k) const{
			return data[k]->getPixel(i, j);
		};

		//ピクセル(i, j)に色をセットする
		void setPixel(int i, int j, int k, const Vec3& c){
			data[k]->setPixel(i, j, c);
		};

		//ピクセル(i, j)に色を追加する
    	void addPixel(int i, int j, double t, const Vec3& c) {
			int k = time2index(t);
      		data[k]->addPixel(i, j, c);
    	};

		//全てのピクセルを一定の値で割る
		void divide(double val){
			for(int k=0; k<frame_num; k++){
				for(int i=0; i<width; i++){
					for(int j=0; j<height; j++){
						this->setPixel(i, j, k, this->getPixel(i, j, k)/val);
					}
				}
			}
		};

		//ガンマ補正
		void gamma_correction(){
			for(int k=0; k<frame_num; k++){
				for(int i=0; i<width; i++){
					for(int j=0; j<height; j++){
						Vec3 c = this->getPixel(i, j, k);
						this->setPixel(i, j, k, Vec3(std::pow(c.x, 1/2.2), std::pow(c.y, 1/2.2), std::pow(c.z, 1/2.2)));
					}
				}
			}
		};

		//PPM画像出力
		void ppm_output(const std::string& filename) const{
			int i = filename.find_last_of(".");
			std::string ext = filename.substr(i, filename.length() - i);
			std::string name = filename.substr(0, i);
			for(int k=0; k<frame_num; k++) {
				std::string filename_k = name + "-" + std::to_string(k+1) + ext;
				data[k]->ppm_output(filename_k);
			}
		};

	private:
		//時間tに対応する画像のindexを返す
		int time2index(const double t) {
			int index = (int)((t - time_min) / delta_t);
			if(index < 0) return 0;
			if(index >= frame_num) return frame_num-1;
			return index;
		};

};


#endif
