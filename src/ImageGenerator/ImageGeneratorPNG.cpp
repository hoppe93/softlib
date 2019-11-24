/**
 * Implementation of module for generating PNG images.
 */

#include <cstdio>
#include <png.h>
#include <softlib/ImageGenerator/ImageColormap.h>
#include <softlib/ImageGenerator/ImageGeneratorPNG.h>

void ImageGeneratorPNG::SaveImage(
    const std::string& name, const slibreal_t *img,
    const unsigned int nrows, const unsigned int ncols,
    ImageColormap *cmap
) {
    spng_bitmap_t *png = new spng_bitmap_t;

    png->pixels = new spng_pixel_t[nrows*ncols];
    png->width  = ncols;
    png->height = nrows;

    for (long long int i = 0, index = 0; i < (long long int)ncols; i++) {
        for (long long int j = 0; j < (long long int)nrows; j++, index++) {
            ImageColormap::scolor_t clr = cmap->Eval(img[j*ncols+i]);

            png->pixels[index].r = clr.red;
            png->pixels[index].g = clr.green;
            png->pixels[index].b = clr.blue;
        }
    }

    FILE *f = fopen(name.c_str(), "wb");
    if (!f)
        throw SOFTLibException("Unable to create PNG file.");

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL)
        throw SOFTLibException("Unable to create PNG 'write struct'.");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL)
        throw SOFTLibException("Unable to create PNG 'info struct'.");

    if (setjmp(png_jmpbuf(png_ptr)))
        throw SOFTLibException("Unable to write PNG.");

    int depth = 8;
    png_set_IHDR(
        png_ptr,
        info_ptr,
        png->width,
        png->height,
        depth,
        PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );

    // Initialize rows of PNG
    png_byte **row_pointers = (png_byte**)png_malloc(png_ptr, png->height * sizeof(png_byte*));
    int pixel_size = 3; // Number of uint8_t's per pixel
    for (size_t y = 0; y < png->height; y++) {
        png_byte *row = (png_byte*)png_malloc(png_ptr, sizeof(uint8_t)*png->width*pixel_size);
        row_pointers[y] = row;

        for (size_t x = 0; x < png->width; x++) {
            spng_pixel_t *pixel = png->pixels+(y*png->width + x);

            *row++ = pixel->r;
            *row++ = pixel->g;
            *row++ = pixel->b;
        }
    }

    png_init_io(png_ptr, f);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    // Dealloc
    for (size_t y = 0; y < png->height; y++)
        png_free(png_ptr, row_pointers[y]);
    png_free(png_ptr, row_pointers);

    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(f);
}

