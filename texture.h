#ifndef TEXTURE_H
#define TEXTURE_H

#include "vec3.h"
//#include "perlin.h"

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif
#include <iostream>


class Texture {
	public:
		virtual RGB value(double u, double v, const Vec3& p) const = 0;
};

class Solid_color : public Texture {
    public:
        Solid_color() {}
        Solid_color(RGB c) : color_value(c) {}

        Solid_color(double red, double green, double blue)
          : Solid_color(RGB(red,green,blue)) {}

        virtual RGB value(double u, double v, const Vec3& p) const {
            return color_value;
        }

    private:
        RGB color_value;
};

class Checker_texture : public Texture {
	public:
		Checker_texture() {}
		Checker_texture(std::shared_ptr<Texture> t0, std::shared_ptr<Texture> t1): even(t0), odd(t1) {}

		virtual RGB value(double u, double v, const Vec3& p) const {
			double sines = std::sin(freq*p.x) /** std::sin(freq*p.y)*/ * std::sin(freq*p.z);
			if (sines < 0)
				return odd->value(u, v, p);
			else
				return even->value(u, v, p);
		}

	public:
        double freq = 0.0625;
		std::shared_ptr<Texture> odd;
		std::shared_ptr<Texture> even;
};
/*
class noise_texture : public Texture {
    public:
        noise_texture() {}
		noise_texture(double sc) : scale(sc) {}

        virtual RGB value(double u, double v, const Vec3& p) const {
			return color(1,1,1) * 0.5 * (1 + sin(scale*p.z() + 10*noise.turb(p)));
        }

    public:
        perlin noise;
		double scale;
};
*/

class Image_texture : public Texture {
    public:
        const static int bytes_per_pixel = 3;

        Image_texture()
          : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}

        Image_texture(const char* filename) {
            auto components_per_pixel = bytes_per_pixel;

            data = stbi_load(
                filename, &width, &height, &components_per_pixel, components_per_pixel);

            if (!data) {
                std::cerr << "ERROR: Could not load Texture image file '" << filename << "'.\n";
                width = height = 0;
            }

            bytes_per_scanline = bytes_per_pixel * width;
        }

        ~Image_texture() {
            delete data;
        }

        virtual RGB value(double u, double v, const Vec3& p) const {
            // If we have no Texture data, then return solid cyan as a debugging aid.
            if (data == nullptr)
                return RGB(0,1,1);

            // Clamp input Texture coordinates to [0,1] x [1,0]
            u = clamp(u, 0.0, 1.0);
            v = 1.0 - clamp(v, 0.0, 1.0);  // Flip V to image coordinates

            auto i = static_cast<int>(u * width);
            auto j = static_cast<int>(v * height);

            // Clamp integer mapping, since actual coordinates should be less than 1.0
            if (i >= width)  i = width-1;
            if (j >= height) j = height-1;

            const auto color_scale = 1.0 / 255.0;
            auto pixel = data + j*bytes_per_scanline + i*bytes_per_pixel;

            return RGB(color_scale*pixel[0], color_scale*pixel[1], color_scale*pixel[2]);
        }

    private:
        unsigned char *data;
        int width, height;
        int bytes_per_scanline;
};

#endif
