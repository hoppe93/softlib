#ifndef _IMAGE_GENERATOR_PPM_H
#define _IMAGE_GENERATOR_PPM_H

#include <string>
#include <softlib/ImageGenerator/ImageColormap.h>
#include <softlib/ImageGenerator/ImageGenerator.h>


class ImageGeneratorPPM : public ImageGenerator {
public:
    virtual void SaveImage(
        const std::string&, const slibreal_t*,
        const unsigned int, const unsigned int,
        ImageColormap*
    ) override;
};

#endif/*_IMAGE_GENERATOR_PPM_H*/
