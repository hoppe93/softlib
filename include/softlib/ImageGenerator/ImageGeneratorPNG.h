#ifndef _IMAGE_GENERATOR_PNG_H
#define _IMAGE_GENERATOR_PNG_H

#include <string>
#include <softlib/ImageGenerator/ImageColormap.h>
#include <softlib/ImageGenerator/ImageGenerator.h>


class ImageGeneratorPNG : public ImageGenerator {
public:
    typedef struct {
        uint8_t r, g, b;
    } spng_pixel_t;
    typedef struct {
        spng_pixel_t *pixels;
        size_t width, height;
    } spng_bitmap_t;

    virtual void SaveImage(
        const std::string&, const slibreal_t*,
        const unsigned int, const unsigned int,
        ImageColormap*
    ) override;
};

#endif/*_IMAGE_GENERATOR_PNG_H*/
