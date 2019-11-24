/**
 * Implementation of the common 'ImageGenerator' interface 'gen()'.
 */

#include <softlib/ImageGenerator/ImageGenerator.h>
#ifdef HAS_LIBPNG
#   include <softlib/ImageGenerator/ImageGeneratorPNG.h>
#endif
#include <softlib/ImageGenerator/ImageGeneratorPPM.h>
#include <softlib/SOFTLibException.h>


/**
 * Save the given image 'img' to a file of type 'type' with
 * name 'name'.
 *
 * type:     Type of image to generate.
 * name:     Name of file to save image to.
 * img:      Image to save (logically two-dimensional).
 * nrows:    Vertical image dimension (first).
 * ncols:    Horizontal image dimension (second).
 * colormap: Colormap to use when generating image.
 */
void ImageGenerator::gen(
    enum image_type type, const std::string& name, const slibreal_t *img,
    const unsigned int nrows, const unsigned int ncols, enum colormap_type colormap
) {
    ImageGenerator *igen = GetImageGenerator(type);
    igen->SaveImage(name, img, nrows, ncols, colormap);
    delete igen;
}
void ImageGenerator::gen(
    enum image_type type, const std::string& name, const slibreal_t *img,
    const unsigned int nrows, const unsigned int ncols, ImageColormap *colormap
) {
    ImageGenerator *igen = GetImageGenerator(type);
    igen->SaveImage(name, img, nrows, ncols, colormap);
    delete igen;
}

/**
 * Returns an 'ImageGenerator' object of the requested type.
 *
 * type: Type of ImageGenerator object to construct.
 */
ImageGenerator *ImageGenerator::GetImageGenerator(enum image_type type) {
    switch (type) {
#ifdef HAS_LIBPNG
        case ImageGenerator::IMAGE_TYPE_PNG:
            return new ImageGeneratorPNG();
#endif

        case ImageGenerator::IMAGE_TYPE_PPM:
            return new ImageGeneratorPPM();

        default:
            throw SOFTLibException(
                "Unrecognized image format specified to 'ImageGenerator': '%d'.", type
            );
    }
}
