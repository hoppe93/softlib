#ifndef _IMAGE_GENERATOR_H
#define _IMAGE_GENERATOR_H

#include <string>
#include <softlib/config.h>
#include <softlib/ImageGenerator/ImageColormap.h>
#include <softlib/SOFTLibException.h>


class ImageGenerator {
public:
    static constexpr ImageColormap::scolor_t icmap_gerimap[9] = {
        {0,0,0},{38,38,128},{76,38,191},{153,51,128},
        {255,64,38},{230,128,0},{230,191,26},
        {230,230,128},{255,255,255}
    };
    static constexpr ImageColormap::scolor_t icmap_gray[2] = {
        {0,0,0},{255,255,255}
    };

    ImageColormap *cmap_gerimap, *cmap_gray;

    enum colormap_type {
        COLORMAP_GERIMAP,
        COLORMAP_GRAY
    };
    enum image_type {
        IMAGE_TYPE_PNG,
        IMAGE_TYPE_PPM
    };

    ImageGenerator() { InitColormaps(); }
    virtual ~ImageGenerator() {}

    void InitColormaps() {
        this->cmap_gerimap = new ImageColormap(
            icmap_gerimap, 9, 0.0, 1.0
        );
        this->cmap_gray = new ImageColormap(
            icmap_gray, 2, 0.0, 1.0
        );
    }

    virtual void SaveImage(
        const std::string& name, const slibreal_t* img,
        const unsigned int m, const unsigned int n,
        const enum colormap_type cmap
    ) {
        switch (cmap) {
            case COLORMAP_GERIMAP:
                return SaveImage(name, img, m, n, this->cmap_gerimap);
            case COLORMAP_GRAY:
                return SaveImage(name, img, m, n, this->cmap_gray);

            default:
                throw SOFTLibException("Unrecognized colormap requested: %d.", cmap);
        }
    }

    virtual void SaveImage(
        const std::string&, const slibreal_t*,
        const unsigned int, const unsigned int,
        ImageColormap*
    ) = 0;

    static void gen(
        enum image_type, const std::string&, const slibreal_t*,
        const unsigned int, const unsigned int, enum colormap_type
    );
    static void gen(
        enum image_type, const std::string&, const slibreal_t*,
        const unsigned int, const unsigned int, ImageColormap *cmap=nullptr
    );
    static ImageGenerator *GetImageGenerator(enum image_type);
};

#endif/*_IMAGE_GENERATOR_H*/
