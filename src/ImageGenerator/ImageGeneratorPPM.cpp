/**
 * Implementation of module for generating PPM images.
 */

#include <cstdio>
#include <softlib/ImageGenerator/ImageColormap.h>
#include <softlib/ImageGenerator/ImageGeneratorPPM.h>

void ImageGeneratorPPM::SaveImage(
    const std::string& name, const slibreal_t *img,
    const unsigned int nrows, const unsigned int ncols,
    ImageColormap *cmap
) {
    throw SOFTLibException("PPM support not yet implemented...");
}
